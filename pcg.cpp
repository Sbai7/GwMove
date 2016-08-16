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
#include <pcg.h>

// C++ headers
#include <iomanip>


namespace happ
{
   // ======================================================
   PcgSolution::PcgSolution()
      : KrylovSolution()
   { }

   // ======================================================
   PcgSolution::PcgSolution(CsrMatrix *matrix, double *b, double *unknown)
      : KrylovSolution(matrix, b, unknown)
   { }

   // ======================================================
   PcgSolution::PcgSolution(const PcgSolution &another_pcg)
      : KrylovSolution(another_pcg)
   { }

   // ======================================================
   PcgSolution::~PcgSolution()
   { }

   // ======================================================
   int PcgSolution::Solve(unsigned int &iterations, double &error)
   {
      if (!_A || !_x || !_rhs || !_r)
         return ERROR_NULL_POINTERS;

      // work vectors
      double *p = new double[_dimension];
      double *q = new double[_dimension];
      double *w = new double[_dimension];

      unsigned int j;
      bool converged = false;

      // do incomplete factorization of system matrix
      if (!_A->IncompleteFactorization())
         return ERROR_ILU0; // error in incomplete factorization

      // calculate initial residual vector ...
      if (_nCalls == 0) _A->BackSubSolve(_dimension, _rhs, _x);
      _A->Multiply(_dimension, _x, _r);
      for (j = 0; j < _dimension; j++) {
         _r[j] = _rhs[j] - _r[j];
      }
      _A->BackSubSolve(_dimension, _r, p);
      error = dot_product(_dimension, _r, p);

      // Start iteration process ...
      for (iterations = 1; iterations <= _maxIterations; iterations++) {

         // error calculation
         if (!silent) {
            std::cout << scientific;
            std::cout << "* Iteration = " << setw(5) << iterations
               << " - Error = " << setw(16) << error
               << " - Tolerance = " << setw(16) << _tolerance << endl;
         }
         if (error < _tolerance) {
            _nCalls++; // we increase number of calls only if we succeed
            converged = true;
            break; // successful return
         }

         _A->Multiply(_dimension, p, q);
         double pq = dot_product(_dimension, p, q);

         // corrections
         double alpha = error / pq;
         for (j = 0; j < _dimension; j++) {
            _x[j] += alpha*p[j];
            _r[j] -= alpha*q[j];
         }
         _A->BackSubSolve(_dimension, _r, w);

         // w contains (LDL')^-1 * r_{i+1}
         double error_new = dot_product(_dimension, _r, w);
         double beta = error_new / error;
         for (j = 0; j < _dimension; j++) {
            p[j] = w[j] + beta*p[j];
         }
         error = error_new;
      }

      delete[] p;
      delete[] q;
      delete[] w;

      if (converged)
         return 1;

      else
         // maximum number of iterations is excedded
         return ERROR_MAX_ITERATIONS;
   }
}