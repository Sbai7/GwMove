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
#include <pbicgstab.h>

// C++ headers
#include <iomanip>


namespace happ
{
   // ======================================================
   PbicgstabSolution::PbicgstabSolution()
      : KrylovSolution()
   { }

   // ======================================================
   PbicgstabSolution::PbicgstabSolution(CsrMatrix *matrix,
      double *b, double *unknown) : KrylovSolution(matrix, b, unknown)
   { }

   // ======================================================
   PbicgstabSolution::PbicgstabSolution(const PbicgstabSolution &another_pbicgstab)
      : KrylovSolution(another_pbicgstab)
   { }

   // ======================================================
   PbicgstabSolution::~PbicgstabSolution()
   { }


   // ======================================================
   int PbicgstabSolution::Solve(unsigned int &iterations, double &error)
   {
      if (!_A || _x || _rhs || _r)
         return ERROR_NULL_POINTERS;

      // work vectors
      double *r0 = new double[_dimension];
      double *r1 = new double[_dimension];
      double *r2 = new double[_dimension];
      double *u0 = new double[_dimension];
      double *u1 = new double[_dimension];
      double *u2 = new double[_dimension];
      double *t0 = new double[_dimension];
      double *t1 = new double[_dimension];
      double *v0 = new double[_dimension];
      double *v1 = new double[_dimension];
      double *v2 = new double[_dimension];

      unsigned int i, j;
      bool converged = false;

      for (i = 0; i < _dimension; i++) u0[i] = 0.;

      // do incomplete factorization of system matrix
      if (!_A->IncompleteFactorization())
         return ERROR_ILU0; // error in incomplete factorization

      // calculate initial residual vector ...
      _A->BackSubSolve(_dimension, _rhs, _x);
      _A->Multiply(_dimension, _x, _r);
      for (j = 0; j < _dimension; j++) {
         _r[j] = _rhs[j] - _r[j];
      }

      // initialize iterative parameters
      double beta, xsi, error, sigma, gamma_prime, gamma_sec, tho;
      std::vector<double> rho, gamma;
      rho.resize(2); gamma.resize(2);
      rho[0] = 1.;
      double omega = 1;
      double alpha = 0.;
      r0 = _r;


      // Start iteration process ...
      for (iterations = 1; iterations <= _maxIterations; iterations++) {

         rho[0] *= -omega;
         j = 0;

      REPEAT:
         if (j == 0) rho[1] = dot_product(_dimension, _r, r0);
         if (j == 1) rho[1] = dot_product(_dimension, r1, r0);

         // Exception ...
         if (rho[0] == 0.) {
            converged = false;
            return ERROR_BICGSTAB_BREAKDOWN;
         }
         else
            beta = (rho[1] / rho[0])*alpha;

         rho[0] = rho[1];
         if (j == 0) {
            for (i = 0; i < _dimension; i++)
               u0[i] = _r[i] - beta*u0[i];
            _A->BackSubSolve(_dimension, u0, t0);
            _A->Multiply(_dimension, t0, u1);
            xsi = dot_product(_dimension, u1, r0);
         }
         else if (j == 1) {
            for (i = 0; i < _dimension; i++) {
               u0[i] = _r[i] - beta*u0[i];
               u1[i] = r1[i] - beta*u1[i];
            }
            _A->BackSubSolve(_dimension, u1, t1);
            _A->Multiply(_dimension, t1, u2);
            xsi = dot_product(_dimension, u2, r0);
         }

         if (xsi == 0.) {
            converged = false;
            return ERROR_BICGSTAB_BREAKDOWN;
         }
         else
            alpha = rho[0] / xsi;

         if (j == 0) {
            for (i = 0; i < _dimension; i++)
               _r[i] -= alpha*t1[i];
            _A->BackSubSolve(_dimension, r0, v0);
            _A->Multiply(_dimension, v0, r1);
         }
         else {
            for (i = 0; i < _dimension; i++) {
               v0[i] -= alpha*t1[i];
               r1[i] -= alpha*u2[i];
            }
            _A->BackSubSolve(_dimension, r1, v1);
            _A->Multiply(_dimension, v1, r2);
         }

         for (i = 0; i < _dimension; i++)
            _x[i] += alpha*t0[i];

         j++;
         if (j <= 1) goto REPEAT;

         // error calculation
         error = dot_product(_dimension, _r, _r);
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

         sigma = dot_product(_dimension, v1, v1);
         if (sigma == 0.) {
            converged = false;
            return ERROR_BICGSTAB_BREAKDOWN;
         }
         else
            gamma_prime = dot_product(_dimension, v0, v1) / sigma;

         tho = dot_product(_dimension, v2, v1) / sigma;
         for (i = 0; i < _dimension; i++)
            r2[i] -= tho*v1[i];
         for (i = 0; i < gamma.size(); i++)
            gamma[i] = gamma_prime - tho*omega;
         gamma_sec = omega;

         for (i = 0; i < _dimension; i++) {
            _x[i] += gamma[0] * v0[i] + gamma_sec*v1[i];
            r0[i] = v0[i] - gamma_prime*v1[i] - omega*r2[i];
            u0[i] = t0[i] + gamma[0] * t1[i] - gamma[1] * u2[i];
         }

      }

      // free memory 
      delete[] r0;
      delete[] r1;
      delete[] r2;
      delete[] u0;
      delete[] u1;
      delete[] u2;
      delete[] t0;
      delete[] t1;
      delete[] v0;
      delete[] v1;
      delete[] v2;

      if (converged)
         return 1;

      else
         // maximum number of iterations is excedded
         return ERROR_MAX_ITERATIONS;
   }
}