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
    mxArray *D_in_m, *P_out_m, *P1_out_m ;
    const mwSize *dims;
    double *D, *P, *P1;
    int M, N, numdims;
    int nxIndex;
    int m,n;
    int nr,nc;
    int wt[3]={1,1,1};
    
//associate inputs
    D_in_m = mxDuplicateArray(prhs[0]);
    nr=mxGetM(prhs[0]);
    nc=mxGetN(prhs[0]);
    
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    M = (int)dims[0]; N = (int)dims[1];
//associate outputs
    double& dist = createMatlabScalar(plhs[0]);
    double& ep = createMatlabScalar(plhs[1]);
    P_out_m = plhs[2] = mxCreateDoubleMatrix(M,N,mxREAL);
    P1_out_m = plhs[3] = mxCreateDoubleMatrix(M,N,mxREAL);
    
//associate pointers
    D = mxGetPr(D_in_m);
    P = mxGetPr(P_out_m);
    P1 = mxGetPr(P1_out_m);

    D_in_m = (mxArray *)mxMalloc(M*N * sizeof(double));
    P_out_m = (mxArray *)mxMalloc(M*N * sizeof(double));
    P1_out_m = (mxArray *)mxMalloc(M*N * sizeof(double));
        
//do something
    P[0]=D[0];
    P1[0]=1;
    
    m=0;
    for(n=1;n<N;n++)
    {
        P[m+M*n] = wt[2]*D[m+M*n]+P[m+M*(n-1)]; //row filling
        P1[m+M*n] = n+1; //row filling
    }
    
    n=0;
    for(m=1;m<M;m++)
    {
        P[m+M*n] = wt[1]*D[m+M*n]+P[(m-1)+M*n]; //column filling
        P1[m+M*n] = m+1; //row filling
    }
    
    double S,V,H;
    for(m=1;m<M;m++)
        for (n=1;n<N;n++)
        {
            {

                S=P[(m-1)+M*(n-1)]+wt[0]*D[m+M*n];
                V=P[(m-1)+M*n]+wt[1]*D[m+M*n];
                H=P[m+M*(n-1)]+wt[2]*D[m+M*n];
                P[m+M*n] = min_fun(S,V,H);
                nxIndex = min_fun_ind(S,V,H);
                switch (nxIndex)
                {
                    case 1:
                        P1[m+M*n]=P1[(m-1)+M*(n-1)]+1; // The subtraction plays a role when there is no change in the consecutive diagonal entries of local distance matrix,
                        break;
                    case 2:
                        P1[m+M*n]=P1[(m-1)+M*n]+1;
                        break;
                        
                    default:
                        P1[m+M*n]=P1[m+M*(n-1)]+1;
                        break;
                }

            }
        }
    
    int EP=find_min_value_ind(P,M,N)+M*(N-1);
    dist=P[EP]/P1[EP];
    ep=find_min_value_ind(P,M,N);
    
    
    /*
     * // Back tracking
     * int i=M-1,j=N-1;
     * i=find_min_value_ind(P,M,N);
     * w1[0]=find_min_value_ind(P,M,N);
     * w2[0]=j;
     *
     * int cnt_path=1,cnt=1;
     *
     * while(1)
     * {
     * S=P[(i-1)+M*(j-1)];
     * V=P[(i-1)+M*j];
     * H=P[i+M*(j-1)];
     * int nxIndex;
     * if(i>0 && j >0)
     * {
     * nxIndex = min_fun_ind(S,V,H); // Preference to diagonal
     * }
     * else
     * {
     * if(i==0)
     * {
     * nxIndex =3;// at first row, so only scope for deletion
     * }
     * else
     * {
     * nxIndex=2; // at first column, so only scope of insertion
     * }
     * }
     * switch (nxIndex)
     * {
     * case 1:
     * i= i-1;
     * j = j-1;
     * cnt_path=cnt_path+wt[0]; // The subtraction plays a role when there is no change in the consecutive diagonal entries of local distance matrix,
     * break;
     * case 2:
     * i = i-1;
     * cnt_path=cnt_path+wt[1];
     * break;
     *
     * default:
     * j = j-1;
     * cnt_path=cnt_path+wt[2];
     * break;
     * }
     *
     * // break condition
     * if(i==0 && j==0)
     * {
     * break;
     * }
     * w1[cnt]=i;
     * w2[cnt]=j;
     * cnt=cnt+1;
     * }
     * int EP=w1[0]+M*w2[0];
     * dist=P[EP]/cnt_path;*/
    
        mxFree(D_in_m);
        mxFree(P_out_m);
        mxFree(P1_out_m);
        
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
    if( ( x <= y ) && ( x <= z ) ) return 1;
    if( ( y <= x ) && ( y <= z ) ) return 2;
    if( ( z <= x ) && ( z <= y ) ) return 3;
}

int find_min_value_ind(double *P, int M, int N)
{
    double temp1=P[0+M*(N-1)],temp2;
    int min_ind=0,m;
    for(m=1;m<M;m++)
    {
        temp2=P[m+M*(N-1)];
        if(temp2<temp1)
        {
            min_ind=m;
            temp1=temp2;
        }
    }
    return min_ind;
}
/* The piece of MATLAB code used to convert the output aligned path into
 * MATLAB index format from c index format
 * /* W3=[w1;w2];
 * w11=flipud(W3(1,sum(W3)>0)'+1);
 * w12=flipud(W3(2,sum(W3)>0)'+1);
 * W=[w11,w12];
 * W1=[1 1;W];
 */