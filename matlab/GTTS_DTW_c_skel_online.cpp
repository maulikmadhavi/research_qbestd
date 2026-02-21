/*********************************************************************
 *This code does the warping path estimation and DTW calculation.
 * Weights are [1 1 1] --> Horizontal, Diagonal, Vertical movement.
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
    mxArray *D_in_m, *A_out_m, *P_out_m, *L_out_m;
    const mwSize *dims;
    double *D, *A, *P, *L;
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
    A_out_m = plhs[2] = mxCreateDoubleMatrix(M,N,mxREAL);
    L_out_m = plhs[3] = mxCreateDoubleMatrix(M,N,mxREAL);
    P_out_m = plhs[4] = mxCreateDoubleMatrix(M,N,mxREAL);
    
//associate pointers
    D = mxGetPr(D_in_m);
    P = mxGetPr(P_out_m);
    A = mxGetPr(A_out_m);
    L = mxGetPr(L_out_m);
    D_in_m = (mxArray *)mxMalloc(M*N * sizeof(double));
    P_out_m = (mxArray *)mxMalloc(M*N * sizeof(double));
    A_out_m = (mxArray *)mxMalloc(M*N * sizeof(double));
    L_out_m = (mxArray *)mxMalloc(M*N * sizeof(double));
    
//do something
    
    // First column initialization
    n=0;
    for(m=0;m<M;m++)
    {
        P[m+M*n] = m+1;
        A[m+M*n] = D[m+M*n];
        L[m+M*n] = 1;
    }
    
    // First Row Initialization
    m=0;
    for (n=1;n<N;n++) // SKIPPING first (column) entry
    {
        P[m+M*n]= 1;
        A[m+M*n]= A[m+M*(n-1)]+D[m+M*n];  // Accumulation
        L[m+M*n]= n+1;
    }
    
    double S1,D1,V1;
    for(m=1;m<M;m++)  //  N1 For every row
        for (n=1;n<N;n++) //N2 For every column
        {
            {
//                V1=(D[m+M*n]+A[(m-1)+M*(n-0)]);
                //               D1=(D[m+M*n]+A[(m-1)+M*(n-1)]);
                //              S1=(D[m+M*n]+A[(m-0)+M*(n-1)]);
                V1=(D[m+M*n]+A[(m-1)+M*(n-0)])/(L[(m-1)+M*(n-0)]+1);
                D1=(D[m+M*n]+A[(m-1)+M*(n-1)])/(L[(m-1)+M*(n-1)]+1);
                S1=(D[m+M*n]+A[(m-0)+M*(n-1)])/(L[(m-0)+M*(n-1)]+1);
                
                nxIndex = min_fun_ind(V1,D1,S1);
                switch (nxIndex)
                {
                    case 1:
                        P[m+M*n]=P[(m-1)+M*(n-0)];
//                        A[m+M*n]=V1;
                        A[m+M*n]=A[(m-1)+M*(n-0)]+D[m+M*n];
                        L[m+M*n]=L[(m-1)+M*(n-0)]+1;
                        break;
                    case 2:
                        P[m+M*n]=P[(m-1)+M*(n-1)];
                        //                      A[m+M*n]=D1;
                        A[m+M*n]=A[(m-1)+M*(n-1)]+D[m+M*n];
                        L[m+M*n]=L[(m-1)+M*(n-1)]+1;
                        break;
                    default:
                        P[m+M*n]=P[(m-0)+M*(n-1)];
//                        A[m+M*n]=S1;
                        A[m+M*n]=A[(m-0)+M*(n-1)]+D[m+M*n];
                        L[m+M*n]=L[(m-0)+M*(n-1)]+1;
                        break;
                }
            }
        }
    
    // Score
    int i=find_min_value_ind(A,M,N);
    int EP=i+M*(N-1);
    dist=(A[EP]/L[EP]);
    ep=i+1;
    
    mxFree(D_in_m);
    mxFree(P_out_m);
    mxFree(A_out_m);
    mxFree(L_out_m);
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

int find_min_value_ind(double *A, int M, int N)
{
    double temp1=A[0+M*(N-1)],temp2;
    int min_ind=0,m;
    for(m=1;m<M;m++)
    {
        temp2=A[m+M*(N-1)];
        if(temp2<temp1)
        {
            min_ind=m;
            temp1=temp2;
        }
    }
    return min_ind;
}
