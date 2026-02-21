/*********************************************************************
 *This code does the warping path estimation and DTW calculation.
 * Weights are [1 1 1] --> Horizontal, Diagonal, Edge movement.
 ********************************************************************/
// Back tracing from minimum end point
#include <matrix.h>
#include <mex.h>

double min_fun( double x, double y, double z );
int min_fun_ind( double x, double y, double z );
int find_min_value_ind(double *P, int M, int N);
double& createMatlabScalar (mxArray*& ptr);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
//declare variables
    mxArray *D_in_m, *S_out_m, *P_out_m, *T_out_m;
    const mwSize *dims;
    double *D, *S, *P, *T;
    int M, N, numdims;
    int nxIndex;
    int m,n;
    
//associate inputs
    D_in_m = mxDuplicateArray(prhs[0]);
    
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    M = (int)dims[0]; N = (int)dims[1];
//associate outputs
    double& dist = createMatlabScalar(plhs[0]);
    double& ep = createMatlabScalar(plhs[1]);
    S_out_m = plhs[2] = mxCreateDoubleMatrix(M,N,mxREAL);
    T_out_m = plhs[3] = mxCreateDoubleMatrix(M,N,mxREAL);
    P_out_m = plhs[4] = mxCreateDoubleMatrix(M,N,mxREAL);
    
    
//associate pointers
    D = mxGetPr(D_in_m);
    S = mxGetPr(S_out_m);
    P = mxGetPr(P_out_m);
    T = mxGetPr(T_out_m);
    D_in_m = (mxArray *)mxMalloc(M*N * sizeof(double));
    S_out_m = (mxArray *)mxMalloc(M*N * sizeof(double));
    P_out_m = (mxArray *)mxMalloc(M*N * sizeof(double));
    T_out_m = (mxArray *)mxMalloc(M*N * sizeof(double));
    
//do something
    // First column initialization
    n=0;
    for(m=0;m<M;m++)
    {
        P[m+M*n] = m+1;
        S[m+M*n] = D[m+M*n];
        T[m+M*n] = 1;
    }
    
    // First one Row Initialization
    m=0;
    
    for (n=1;n<N;n++) // SKIPPING first (column) entry (which was already initialized)
    {
        P[m+M*n]= m+1;
        S[m+M*n]= S[m+M*(n-1)]+D[m+M*n];  // Accumulation
        T[m+M*n]= n+1;
    }
    
    
    // Second column Initialization
    n=1;       double S1,D1,E1;
    
    for (m=1;m<M;m++) // SKIPPING first (column) and first (row) entry (which was already initialized)
    {
        D1=(D[m+M*n]+S[(m-1)+M*(n-1)]);
        S1=(D[m+M*n]+S[(m-1)+M*(n-0)]);
        nxIndex = min_fun_ind(D1+1,D1,S1);
        
        switch (nxIndex)
        {
            case 2:
                T[m+M*n]=T[(m-1)+M*(n-1)]+2;
                P[m+M*n]=P[(m-1)+M*(n-1)];
                S[m+M*n]=D1;
                break;
            case 3:
                T[m+M*n]=T[(m-1)+M*(n-0)]+1;
                P[m+M*n]=P[(m-1)+M*(n-0)];
                S[m+M*n]=S1;
                break;
            default:
                break;
        }
        
    }
    
    
    
    for(m=1;m<M;m++)  //  N1 (SKIIPPING first two rows)
    {
        for (n=2;n<N;n++) //N2 (SKIPPING first one column)
        {
            E1=(D[m+M*n]+S[(m-1)+M*(n-2)]);
            D1=(D[m+M*n]+S[(m-1)+M*(n-1)]);
            S1=(D[m+M*n]+S[(m-1)+M*(n-0)]);
            
            nxIndex = min_fun_ind(E1,D1,S1);
            
            switch (nxIndex)
            {
                case 1:
                    T[m+M*n]=T[(m-1)+M*(n-2)]+1;
                    P[m+M*n]=P[(m-1)+M*(n-2)];
                    S[m+M*n]=E1;
                    break;
                case 2:
                    T[m+M*n]=T[(m-1)+M*(n-1)]+1;
                    P[m+M*n]=P[(m-1)+M*(n-1)];
                    S[m+M*n]=D1;
                    break;
                default:
                    T[m+M*n]=T[(m-1)+M*(n-0)]+1;
                    P[m+M*n]=P[(m-1)+M*(n-0)];
                    S[m+M*n]=S1;
                    break;
            }
        }
    }
//
// Score
    
    int i=find_min_value_ind(S,M,N);
    int EP=i+M*(N-1);
    dist=(S[EP]/T[EP]);
    ep=i+1;
    
    mxFree(D_in_m);
    mxFree(S_out_m);
    mxFree(P_out_m);
    mxFree(T_out_m);
    return;
}
double& createMatlabScalar (mxArray*& ptr) {
    ptr = mxCreateDoubleMatrix(1,1,mxREAL);
    return *mxGetPr(ptr);
}

double min_fun( double x, double y, double z )
{
    if( ( x <= y ) && ( x <= z ) ) return x;
    if( ( y <= x ) && ( y <= z ) ) return y;
    if( ( z <= x ) && ( z <= y ) ) return z;
}

int min_fun_ind( double x, double y, double z )
{
    if( ( z <= x ) && ( z <= y ) ) return 3;
    if( ( y <= x ) && ( y <= z ) ) return 2;
    if( ( x <= y ) && ( x <= z ) ) return 1;
    
}

int find_min_value_ind(double *S, int M, int N)
{
    double temp1=S[0+M*(N-1)],temp2;
    int min_ind=0,m;
    for(m=1;m<M;m++)
    {
        temp2=S[m+M*(N-1)];
        if(temp2<temp1)
        {
            min_ind=m;
            temp1=temp2;
        }
    }
    return min_ind;
}