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

#ifndef HAPP_KRYLOV
#define HAPP_KRYLOV

// Project headers
#include <csr.h>

// C++ headers
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector> 
using namespace std;

#define ERROR_NULL_POINTERS			-1
#define ERROR_ILU0					-2
#define ERROR_MAX_ITERATIONS		-3
#define ERROR_BICGSTAB_BREAKDOWN	-4


namespace happ
{
   /**
    * Abstract base class for an iterative subspace iterative
    * krylov solver.
    * This class must be used to subclass iterative krylov
    * solvers for solving sparse linear systems of equations
    * such as the conjugate gradient method, the biconjugate
    * gradient variants, etc.
    */
   class KrylovSolution
   {
   protected:
      /// linear equations matrix (in sparse form)
      CsrMatrix *_A;

      /// unknown solution vector
      double *_x;

      /// right-hand side vector
      double *_rhs;

      /// residual vector
      double *_r;

      /// system dimension
      unsigned int _dimension;

      /// maximum number of iterations
      unsigned int _maxIterations;

      /// tolerance criteria for stopping the iteration process
      double _tolerance;

      /// number of calls to Solve routine
      unsigned int _nCalls;

      void _PrintVector(const unsigned int N, const double *vec);

   public:
      /// indicates whether we print iteration messages or not to stdout
      bool silent;

      /// default constructor
      KrylovSolution() {
         _A = NULL;
         _x = NULL;
         _rhs = NULL;
         _r = NULL;
         _nCalls = 0;
         silent = true;
      };

      /// constructor from linear system components
      KrylovSolution(CsrMatrix *matrix, double *b, double *unknown) {
         _A = matrix; // we need to implement copy constructor for CsrMatrix class
         _x = unknown;
         _rhs = b;

         _dimension = _A->GetDim(); // set system dimension from matrix dimension

         // residual vector
         _r = new double[_dimension];
         for (unsigned int i = 0; i < _dimension; i++) _r[i] = 0.;

         // number of external solver calls
         _nCalls = 0;

         // silent iteration process (by default)
         silent = true;
      };

      /// Copy constructor 
      KrylovSolution(const KrylovSolution &another_solver) {
         _A = another_solver._A;
         _x = another_solver._x;
         _rhs = another_solver._rhs;
         _dimension = another_solver._dimension;
         _maxIterations = another_solver._maxIterations;
         _tolerance = another_solver._tolerance;
         silent = another_solver.silent;
         _nCalls = 0;

         _r = another_solver._r;
         for (unsigned int i = 0; i < _dimension; i++) _r[i] = 0.;
      };

      /// default destructor
      virtual ~KrylovSolution() {
         if (_r) delete[] _r;
      };

      /// Solves the linear system of equations
      virtual int Solve(unsigned int &iterations, double &error) = 0;

      /// Returns the new solution
      double *GetSolution();

      /// Sets matrix vector
      void SetMatrix(CsrMatrix *matrix);

      /// Sets initial unknowns from a pointer 
      void SetUnknowns(double *X);

      /// Sets initial unknowns from a std::vector 
      void SetUnknowns(std::vector<double>& vec);

      /**
       * Sets the right-hand side vector.
       * useful for problems with changing rhs and same sparse matrix.
       */
      void SetRhs(double *b);

      /// Sets rhs from a std::vector 
      void SetRhs(std::vector<double>& vec);

      /// Sets max allowed number of iterations
      void SetMaxIterations(const unsigned int maxIter);

      /// Gets max allowed number of iterations
      unsigned int GetMaxIterations();

      /// Sets convergence tolerance
      void SetTolerance(const double tol);

      /// Gets convergence tolerance
      double GetTolerance();

   };


   // ======================================================
   // PcgSolution inlined member functions

   inline
      double *KrylovSolution::GetSolution() {
      return _x;
   }

   // ======================================================
   inline
      void KrylovSolution::SetMatrix(CsrMatrix *matrix) {
      _A = matrix;
      _dimension = _A->GetDim();
   }

   // ======================================================
   inline
      void KrylovSolution::SetUnknowns(double *X) {
      _x = X;
   }

   // ======================================================
   inline
      void KrylovSolution::SetUnknowns(std::vector<double>& vec) {
      _x = vec.data();
   }

   // ======================================================
   inline
      void KrylovSolution::SetRhs(double *b) {
      _rhs = b;
   }

   // ======================================================
   inline
      void KrylovSolution::SetRhs(std::vector<double>& vec) {
      _rhs = vec.data();
   }

   // ======================================================
   inline
      void KrylovSolution::SetMaxIterations(const unsigned int maxIter) {
      _maxIterations = maxIter;
   }

   // ======================================================
   inline
      unsigned int KrylovSolution::GetMaxIterations() {
      return _maxIterations;
   }

   // ======================================================
   inline
      void KrylovSolution::SetTolerance(const double tol) {
      _tolerance = tol;
   }

   // ======================================================
   inline
      double KrylovSolution::GetTolerance() {
      return _tolerance;
   }

   inline
      void KrylovSolution::_PrintVector(const unsigned int N, const double *vec) {
      filebuf file;
      file.open("debug.txt", ios::out);
      ostream os(&file);

      os << setiosflags(ios::scientific) << setprecision(8);
      for (unsigned int i = 0; i < N; i++)
         os << vec[i] << endl;

      file.close();
   }
}

#endif // HAPP_KRYLOV
