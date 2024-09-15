#include <iostream>
#include <mpi.h>

using std :: cout;
using std :: endl;

double lookup(int n, double *x, double *y, double xval)
{
  for ( int i = 0 ; i < n ; ++i)
      if ( xval >= x[i] && xval <= x[i+1] )
      return y[i] + (xval - x[i]) * (y[i+1]-y[i]) / (x[i+1]-x[i]);
  return 0;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    
    int numPE, myPE;
    
    // Initializing number of processors and the rank of the processor
    MPI_Comm_size(MPI_COMM_WORLD, &numPE);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
    
    
    //Initilazing the basic data structures that are needed for the linear interpolation task
    int n = 100;
    double x[n], y[n];
    
    for (int i = 0; i<n; ++i) {
        x[i]=i;
        y[i]=i*i;
    }
    
    
    /*
     Creating the data structures for xval and yval and yvallocal.
     yvallocal is used by each individual processor to store the output of the linear interpolation.
    */
    
    int m = 20;
    double xval[m], yvallocal[5], yval[m];
    
    for (int i = 0; i<m; ++i) {
        xval[i] = 2*i;
        yval[i]=0;
    }
    
    /*
     Broadcasting the values x, y and the array xval so that each other processor can use the values in order to compute
     the associated yval value.
     */
    
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(xval, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    MPI_Bcast(yval, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
    
    /*
     
     Each process will only compute the yval values associated with portions of xval. Since each processor should
     perform roughly the same amount of work, the total of m=20 values can be divided equally among the 4 processors.
     
     Each processor will compute 4 yval values and store them in yvallocal buffer. All the yvallocal buffers are then finally
     gathered and merged into the yval buffer.
     
     
     */
    
    if( myPE == 0) {
        for (int j=0; j<5; ++j) {
            yvallocal[j] = lookup(n, x, y, xval[j]);
        }
    }
    
    /*
     
     Indices for the yvallocal buffer need to be computed properly for processors 1,2 and 3.
     
     */
    
    else if ( myPE == 1) {
        for (int j=5; j<10; ++j) {
            yvallocal[j-5] = lookup(n, x, y, xval[j]);
//            cout << yvallocal[j] << endl;         Line used for debugging
        }
    }
    else if ( myPE == 2) {
        for (int j=10; j<15; ++j) {
            yvallocal[j-10] = lookup(n, x, y, xval[j]);
        }
    }
    else if ( myPE == 3) {
        for (int j=15; j<20; ++j) {
            yvallocal[j-15] = lookup(n, x, y, xval[j]);
        }
    }
    
    // Gather the yvallocal buffers into the yval buffer at process 0.
    MPI_Gather(yvallocal, 5, MPI_DOUBLE, yval, 5, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
    
    if(myPE==0) {
        cout << "Final yval values: " << endl;
        for(int i=0;i<m;i++) {
            cout << "yval[" << i << "]: " << yval[i] << endl;
        }
    }
    
    MPI_Finalize();
    return 0;
    
}


