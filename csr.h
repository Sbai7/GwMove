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

#ifndef HAPP_CSR
#define HAPP_CSR


namespace happ
{
   /**
    * Compressed sparse row square matrix class.
    */
   class CsrMatrix
   {
   private:
      /// matrix dimension
      unsigned int dimension;

      /// number of elements
      unsigned int ne;

      /// finite elements connectivity
      unsigned int *connectivity;

      /// initialize main class members.
      void Init();

      /// builds diagonal and column pointers of sparse matrix
      bool BuildPointers();

      void PrintVector(unsigned int N, double *vec);

   public:
      /// periodic size of FE connectivity (e.g. 8 for hexahedral elements-nodes connectivity)
      unsigned int size;

      /// diagonal elements
      double *diagonal;

      /// diagonal pointers
      unsigned int *pd;

      /// column pointers
      unsigned int *pc;

      /// compressed sparse matrix (non-zero entries)
      double *matrix;


      /// default constructor
      CsrMatrix(unsigned int dim);

      /// second constructor (from FEM connectivity list)
      CsrMatrix(unsigned int dim, unsigned int *connect, unsigned int nelem);

      /// default destructor
      virtual ~CsrMatrix();

      /// initialize the matrix elements
      void InitMatrix();

      /// return matrix dimension
      unsigned int GetDim() {
         return dimension;
      }

      /// sets matrix dimension
      void SetDim(unsigned int dim) {
         dimension = dim;
      }

      /// adjusts positive definitness of the sparse matrix
      void AdjustPositiveDefiniteness();

      /**
       * Incomplete factorisation of a symmetric and positive definite sparse matrix.
       * The decomposition is done such that L is the strictely lower part, P^-1 is
       * a diagonal matrix so that the resultant matrix LDLt, having the same nonzero
       * structure as the original matrix, and the same diagonal entries.
       */
      bool IncompleteFactorization();

      /**
       * Matrix vector multiplication operation. We multiply a symmetric matrix
       * given in the lower triangular compact form by vector x. The resultant
       * is vector y.
       */
      void Multiply(const unsigned int N, double *x, double *y);

      /**
       * Solves algebraic equations by back-substitution of lower triangular
       * matrix in LDL'y=x, where x is known and y unknown.
       */
      void BackSubSolve(const unsigned int N, double *x, double *y);

   };

   /// Dot-product function of two vectors.
   inline
      double dot_product(const unsigned int N, double *a, double *b)
   {
      unsigned int i;
      double sum = 0.;

      for (i = 0; i < N; i++) {
         sum += a[i] * b[i];
      }

      return sum;
   }
}

#endif // HAPP_CSR