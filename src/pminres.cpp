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
#include <pminres.h>

// C++ headers
#include <iomanip>


namespace happ
{
   // ======================================================
   PminresSolution::PminresSolution()
      : KrylovSolution()
   { }

   // ======================================================
   PminresSolution::PminresSolution(CsrMatrix *matrix, double *b, double *unknown)
      : KrylovSolution(matrix, b, unknown)
   { }

   // ======================================================
   PminresSolution::PminresSolution(const PminresSolution &another_pminres)
      : KrylovSolution(another_pminres)
   { }

   // ======================================================
   PminresSolution::~PminresSolution()
   { }

   // ======================================================
   int PminresSolution::Solve(unsigned int &iterations, double &error)
   {
      if (!_A || !_x || !_rhs || !_r)
         return ERROR_NULL_POINTERS;

      // work vectors
      double *u = new double[_dimension];
      double *p = new double[_dimension];

      unsigned int j;
      bool converged = false;

      // do incomplete factorization of system matrix
      if (!_A->IncompleteFactorization())
         return ERROR_ILU0; // error in incomplete factorization

      // calculate initial residual vector ...
      _A->BackSubSolve(_dimension, _rhs, _x);
      _A->Multiply(_dimension, _x, _r);
      for (j = 0; j < _dimension; j++) {
         _r[j] = _rhs[j] - _r[j];
      }

      // Start iteration process ...
      for (iterations = 1; iterations <= _maxIterations; iterations++) {

         _A->BackSubSolve(_dimension, _r, u);
         _A->Multiply(_dimension, u, p);
         double rp = dot_product(_dimension, _r, p);
         double rr = dot_product(_dimension, _r, _r);

         // corrections
         double alpha = rr / rp;
         for (j = 0; j < _dimension; j++) {
            _x[j] += alpha*u[j];
            _r[j] -= alpha*p[j];
         }

         // error calculation
         double error = dot_product(_dimension, _r, _r);
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

      }

      delete[] p;
      delete[] u;

      if (converged)
         return 1;

      else
         // maximum number of iterations is excedded
         return ERROR_MAX_ITERATIONS;
   }
}
