//
//  lookup.c
//  
//
//  Created by Varun ShankarHoskere on 9/3/24.
//

#include <stdio.h>
#include <iostream>
#include <mpi.h>
using std :: cout;
using std :: endl;

double lookupVal(int n, double *x, double *y, double xval)
{
  for ( int i = 0 ; i < n ; ++i)
      if ( xval >= x[i] && xval <= x[i+1] )
      return y[i] + (xval - x[i]) * (y[i+1]-y[i]) / (x[i+1]-x[i]);
  return 0;
}

void lookupE(int myPE, int n, int m, double *x, double *y, double *xval, double *yval)
{
    if ( myPE == 0) {
        for ( int j = 0; j < 4;  ++j) {
            for ( int i = 0 ; i < n ; ++i)
                if ( xval[j] >= x[i] && xval[j] <= x[i+1] )
                yval[i] = y[i] + (xval[j] - x[i]) * (y[i+1]-y[i]) / (x[i+1]-x[i]);
        }
    }
    
    else if ( myPE == 1) {
        for ( int j = 4; j < 10;  ++j) {
            for ( int i = 0 ; i < n ; ++i)
                if ( xval[j] >= x[i] && xval[j] <= x[i+1] )
                yval[i] = y[i] + (xval[j] - x[i]) * (y[i+1]-y[i]) / (x[i+1]-x[i]);
        }
    }
    
    else if ( myPE == 2) {
        for ( int j = 9; j < 15;  ++j) {
            for ( int i = 0 ; i < n ; ++i)
                if ( xval[j] >= x[i] && xval[j] <= x[i+1] )
                yval[i] = y[i] + (xval[j] - x[i]) * (y[i+1]-y[i]) / (x[i+1]-x[i]);
        }
    }
    
    else if ( myPE == 3) {
        for ( int j = 15; j < 20;  ++j) {
            for ( int i = 0 ; i < n ; ++i)
                if ( xval[j] >= x[i] && xval[j] <= x[i+1] )
                yval[i] = y[i] + (xval[j] - x[i]) * (y[i+1]-y[i]) / (x[i+1]-x[i]);
        }
    }
}



int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int numPE, myPE;
    
    MPI_Comm_size(MPI_COMM_WORLD,  &numPE);     // Assign the number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE  );     // Current processor #
    
    
    int m = 20;
    int n = 100;
    double x[n], y[n];

    for ( int i = 0 ; i < n ; ++i)
    {
      x[i] = i;
      y[i] = i*i;
    }

//    double xval = 2.5;  
//    double yval = lookupVal(n,x,y,xval);
    
    double xval[m];
    double yval[m];
    
    for (int i = 0; i<m; ++i) {
        xval[i] = 2*i;
    }
    
    MPI_Bcast(x, n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(y, n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(xval, m, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(yval, m, MPI_INT, 0, MPI_COMM_WORLD);
    
    lookupE(myPE, n, m, x, y, xval, yval);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    for (int i = 0; i<m; ++i) {
        cout << yval[i] << endl;
    }

    return 0;
  
}
