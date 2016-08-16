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

#ifndef LAM_STEADY_SWI
#define LAM_STEADY_SWI

// Project headers
#include <lamsflow.h>


/**
* A Layered Adaptive Meshing class for modelling unconfined ground-water flow
* in subsurface aquifers with the sharp-interface saltwater intrusion process. 
* It has an additional support for salt-fresh water interface to which the 
* finite element mesh is automatically adjusted. 
* This class inherits all boundary conditions supported by the unconfined flow 
* process and adds support for the outflow sea face (SGD) boundary condition 
* from which fresh groundwater discharges to the sea. 
*/
class LamSteadySWI : public LamSteadyFlow
{
private:

   /// additional tolerance or stopping criterion for saltwater interface iteration
   double xiTolerance;

   /// saltwater interface position error 
   double swiError;

protected:

   /// User-specfified detailed mass balance analysis 
   virtual void UserMBA(const double *flowRate) const {};

public:

   /// Initialize a saltwater intrusion problem in an unconfined aquifer
   LamSteadySWI(Project*,
      StrucMesh*,
      BndCondition*,
      _Pcg_Params_&);

   /// default destructor
   virtual ~LamSteadySWI() {};

   /**
   * Adjusts the height of nodes lying on (i) the top water table slice,
   * and on (ii) on the lower saltwater intrusion interface slice. 
   * The closest mesh slice close to the m.s.l is the one which is 
   * taken as a reference (i.e. for which Z positions are fixed). 
   * Oher slices heights are found by linear interpolation between
   * the most upper slice and the lowest moving slice (and vice-versa) 
   * in the layered three-dimensional mesh.
   */
   virtual void SetInterfacePosition();

   /// Assembly of RHS array for specialized boundary conditions to this class.
   virtual int AssembleRhs() {return 1;};

   /// The conforming finite element solver for steady saltwater intrusion in 
   /// an unconfined aquifer process using a 'double' layared adaptive mesh technique.
   virtual int Solve() { return 1; };

   /**
   * Prints final summary and mass balance error to standard output
   * stream at the end of the numerical simulation.
   */
   virtual void WriteFinalSummary() {};
};

#endif // LAM_STEADY_SWI
