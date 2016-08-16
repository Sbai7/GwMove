/*****************************************************************************
*
* This file is part of GwMove hydrogeological software developed by
* Dr. M. A. Sbai
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*
*****************************************************************************/

// Project headers
#include <csr.h>

// C++ headers
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


namespace happ
{
   void CsrMatrix::Init()
   {
      unsigned int i;

      matrix = new double[14 * dimension];
      diagonal = new double[dimension];
      pd = new unsigned int[dimension];
      pc = new unsigned int[14 * dimension];

      for (i = 0; i < 14 * dimension; i++) {
         matrix[i] = 0.;
         pc[i] = 0;
      }
      for (i = 0; i < dimension; i++) {
         diagonal[i] = 0.;
         pd[i] = 0;
      }
   }

   CsrMatrix::CsrMatrix(unsigned int dim)
   {
      dimension = dim;
      connectivity = NULL;
      Init();
   }

   CsrMatrix::CsrMatrix(unsigned int dim, unsigned int *connect, unsigned int nelem)
   {
      size = 8; // FOR 3D ONLY
      connectivity = connect;
      ne = nelem;
      dimension = dim;
      Init();
      BuildPointers();
   }

   CsrMatrix::~CsrMatrix()
   {
      if (matrix)   delete[] matrix;
      if (diagonal) delete[] diagonal;
      if (pd)       delete[] pd;
      if (pc)       delete[] pc;
   }

   void CsrMatrix::InitMatrix()
   {
      unsigned int i;

      for (i = 0; i < 14 * dimension; i++)
         matrix[i] = 0.;

      for (i = 0; i < dimension; i++)
         diagonal[i] = 0.;
   }

   bool CsrMatrix::BuildPointers()
   {
      unsigned int i, k;

      // Initialize local arrays
      unsigned int *nr = new unsigned int[dimension];
      unsigned int *nelem = new unsigned int[size*dimension];
      for (i = 0; i < dimension; i++) {
         nr[i] = 0;
         nelem[i] = 0;
      }

      // Counts and stores elements connected to each node
      for (i = 0; i < ne; i++) {
         // if ( soilId(i) <= 0 ) continue; // Must have soilId info from mesh object
         // to add excluding cells feature
         for (k = 0; k < size; k++) {
            int n = connectivity[size*i + k];
            nr[n]++;
            nelem[size*n + nr[n] - 1] = i;
         }
      }

      /************************ make pointers ***************************/

      // for first matrix line (i.e. for mesh node)
      pd[0] = 0;
      pc[0] = 0;
      unsigned int kmax = 0;

      // for the rest of the lines (or nodes)
      for (unsigned int n = 1; n < dimension; n++) {

         for (i = 0; i < nr[n]; i++) {
            int ie = nelem[size*n + i];
            for (unsigned int j = 0; j < size; j++) {
               unsigned int in = connectivity[size*ie + j];
               if (in >= n) continue;

               bool cycle = false;
               for (k = pd[n - 1] + 1; k <= kmax; k++) {
                  if (pc[k] == in) {
                     cycle = true;
                     break;
                  }
               }
               if (!cycle) {
                  kmax++;
                  pc[kmax] = in;
               }

            }
         }
         kmax++;
         pc[kmax] = n;
         pd[n] = kmax;

      }

      delete[] nr;
      delete[] nelem;

      return true;
   }

   void CsrMatrix::AdjustPositiveDefiniteness()
   {
      unsigned int i, j, k;

      double *sum = new double[dimension];
      for (i = 0; i < dimension; i++)
         sum[i] = 0.;

      // calculate the sum of off-diagonal terms, then
      for (i = 1; i < dimension; i++) {
         for (k = pd[i - 1] + 1; k <= pd[i] - 1; k++) {
            j = pc[k];
            sum[i] += matrix[k];
            sum[j] += matrix[k];
         }
      }

      // consider diagonal terms, and add negative rows on diagonals
      for (i = 0; i < dimension; i++) {
         j = pd[i];
         sum[i] += matrix[j];
         if (sum[i] < 0.) matrix[j] -= sum[i];
      }

      // free memory
      delete[] sum;
   }

   bool CsrMatrix::IncompleteFactorization()
   {
      unsigned int i = 0;
      unsigned int j, k;

      // make an M-matrix transformation
      for (k = 0; k < pd[dimension - 1] + 1; k++) {
         if (k == pd[i]) {
            // find the diagonal element
            diagonal[i] += matrix[k];
            i++;
         }
         else {
            if (matrix[k] > 0.) {
               // add to the diagonal
               diagonal[i] += matrix[k];
               j = pc[k];
               diagonal[j] += matrix[k];
            }
         }
      }

      // calculate the diagonals
      for (i = 1; i < dimension; i++) {
         for (k = pd[i - 1] + 1; k < pd[i]; k++) {
            j = pc[k];
            if (matrix[k] < 0.)
               diagonal[i] -= matrix[k] * matrix[k] / diagonal[j];
         }
         if (diagonal[i] <= 0.) return false;
      }

      return true;
   }

   void CsrMatrix::Multiply(const unsigned int N, double *x, double *y)
   {
      unsigned int i, j, k;

      // calculate contribution from matrix diagonals
      for (i = 0; i < N; i++) {
         y[i] = matrix[pd[i]] * x[i];
      }

      // add other terms
      for (i = 1; i < N; i++) {
         for (k = pd[i - 1] + 1; k < pd[i]; k++) {
            j = pc[k];

            // add lower triangular term
            y[i] += matrix[k] * x[j];

            // add upper triangular term
            y[j] += matrix[k] * x[i];
         }
      }
   }

   void CsrMatrix::BackSubSolve(const unsigned int N, double *x, double *y)
   {
      double *u = new double[N];
      double *z = new double[N];
      unsigned int i, j, k;

      // Solve L*u = x
      u[0] = x[0] / diagonal[0];
      for (i = 1; i < N; i++) {
         double sum = x[i];
         for (j = pd[i - 1] + 1; j < pd[i]; j++) {
            k = pc[j];
            if (matrix[j] < 0.)
               sum -= matrix[j] * u[k];
         }
         u[i] = sum / diagonal[i];
      }

      // solve Dz = u
      for (i = 0; i < N; i++) {
         z[i] = diagonal[i] * u[i];
      }

      // solve L'*y = z
      for (i = N - 1; i > 0; i--) {
         j = pd[i];
         y[i] = z[i] / diagonal[i];
         for (k = pd[i] - 1; k >= pd[i - 1] + 1; k--) {
            int ii = pc[k];
            if (matrix[k] < 0.)
               z[ii] -= y[i] * matrix[k];
         }
      }
      y[0] = z[0] / diagonal[0];

      delete[] u;
      delete[] z;
   }

   void CsrMatrix::PrintVector(unsigned int N, double *vec)
   {
      filebuf file;
      file.open("debug.txt", ios::out);
      ostream os(&file);

      os << setiosflags(ios::scientific) << setprecision(8);
      for (unsigned int i = 0; i < N; i++)
         os << vec[i] << endl;

      file.close();
   }
}
