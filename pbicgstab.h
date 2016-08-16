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

#ifndef HAPP_PBICGSTAB
#define HAPP_PBICGSTAB

// Project headers
#include <krylov.h>


namespace happ
{
   /**
    * Preconditioned biconjugate gradient stabilised solver class.
    * Solves the linear algebraic system of equations Ax=b by
    * the preconditioned biconjugate gradient stabilized method,
    * where A is a non symmetric real square matrix, x is the
    * unknown vector and b is the right hand-side vector.
    */
   class PbicgstabSolution : public KrylovSolution
   {
   public:

      /// default constructor
      PbicgstabSolution();

      /// constructor from linear system components
      PbicgstabSolution(CsrMatrix *matrix, double *b, double *unknown);

      /// Copy constructor 
      PbicgstabSolution(const PbicgstabSolution &another_pbicgstab);

      /// default destructor
      virtual ~PbicgstabSolution();

      /// Solves the linear system of equations using the PBicgStab method
      virtual int Solve(unsigned int &iterations, double &error);
   };
}

#endif // HAPP_PBICGSTAB
