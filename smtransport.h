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

#ifndef _MASS_TRANSPORT_STEADY_FEM_H_
#define _MASS_TRANSPORT_STEADY_FEM_H_

// project headers 
#include <sflow.h>

// C++ headers 
#include <vector>
using namespace std;
using namespace happ; 


/**
 * Class for steady-state mass/solute transport modelling of confined/unconfined
 * and variably saturated multilayered aquifer systems. This class could be reused
 * for other purposes such as modelling groundwater age, etc. 
 */
class SteadyMassTransport : public SteadyFlow
{
private:

protected:
   /// nodal flow rate (previously computed from steady state flow) 
   double *Q; 

   /// preconditionned biconjugate gradient stabilised solution (BiCGStab) 
 //todo   PbicgstabSolution bicgstab; // class to build next

   /// does finite element transport matrix assembly from local matrices
   virtual int AssembleStiffnessMatrix();

   /// does assembly of rhs vector from boundary conditions
   virtual int AssembleRhs();

   /// eliminates algebraic equations rows corresponding to fixed potential nodes
   virtual void TidyEquations();

public:
   /// Initialize a solute transport FE problem from a given mesh
   SteadyMassTransport(Project*,
                       StrucMesh*,
                       double*,
                       BndCondition*,
                      _Pcg_Params_&);

   /// default destructor
   virtual ~SteadyMassTransport();

   /**
   * Prints final summary and mass balance error to standard output
   * stream at the end of the numerical simulation.
   */
   virtual void WriteFinalSummary();

};

#endif // #ifndef _MASS_TRANSPORT_STEADY_FEM_H_