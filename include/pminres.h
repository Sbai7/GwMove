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

#ifndef HAPP_PMINRES
#define HAPP_PMINRES

// Project headers
#include <krylov.h>


namespace happ
{
   /**
    * Preconditioned minimum residual solver class.
    * Solves the linear algebraic system of equations Ax=b by
    * the preconditioned minimum residual method, where A is a
    * non-symmetric real square matrix, x is the unknown vector
    * and b is the right hand-side vector.
    */
   class PminresSolution : public KrylovSolution
   {
   public:

      /// default constructor
      PminresSolution();

      /// constructor from linear system components
      PminresSolution(CsrMatrix *matrix, double *b, double *unknown);

      /// Copy constructor 
      PminresSolution(const PminresSolution &another_pminres);

      /// default destructor
      virtual ~PminresSolution();

      /// Solves the linear system of equations using the PMinres method
      virtual int Solve(unsigned int &iterations, double &error);
   };
}

#endif // HAPP_PMINRES
