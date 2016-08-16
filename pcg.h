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

#ifndef HAPP_PCG
#define HAPP_PCG

// Project headers
#include <krylov.h>


namespace happ
{
   /**
    * Preconditioned conjugate gradient solver class.
    * Solves the linear algebraic system of equations Ax=b by
    * the preconditioned conjugate gradient method, where A is
    * a positive definite symmetric real square matrix, x is
    * the unknown vector and b is the right hand-side vector.
    */
   class PcgSolution : public KrylovSolution
   {
   public:
      /// default constructor
      PcgSolution();

      /// constructor from linear system components
      PcgSolution(CsrMatrix *matrix, double *b, double *unknown);

      /// Copy constructor 
      PcgSolution(const PcgSolution &another_pcg);

      /// default destructor
      virtual ~PcgSolution();

      /// Solves the linear system of equations using the Pcg method
      virtual int Solve(unsigned int &iterations, double &error);
   };
}

#endif // HAPP_PCG