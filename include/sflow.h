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

#ifndef STEADY_FLOW
#define STEADY_FLOW

// Project headers
#include <project.h>
#include <strucmesh.h>
#include <hex8.h>
#include <bndcondition.h>
#include <pcg.h>
#include <soil.h>

// C++ headers
#include <vector>
using namespace std;
using namespace happ; 


/**
 * Conforming finite element class problem definition.
 */
class SteadyFlow
{
protected:
   /// problem dimension
   unsigned int dimension;

   /// flag indicating if the simulation is steady state
   bool isTransient;

   /// nodal unknowns (pressure heads, potential heads, etc)
   double *var;

   /// finite element matrix in CSR form
   CsrMatrix *mat;

   /// mass-lumped diagonal storage matrix
   double *mass;

   /// right-hand side of fem matrix equations arising from boundary conditions
   double *rhs;

   /// internaly used copy of conductance matrix
   double *matrixCopy;

   /// nodalwise groundwater table depth from topographic surface 
   double* gwDepth;

   /// pointer to mesh used for computational tasks
   StrucMesh *mesh;

   /// pointer to boundary conditions of the FE problem
   BndCondition *bc;

   /// pointer to project data 
   Project *project; 

   /// preconditionned conjugate gradient solver
   PcgSolution *pcg;

    /// Default constructor 
   SteadyFlow() {};

  /// does finite elements conductance matrix assembly from local matrices
   virtual int AssembleStiffnessMatrix();

   /// assembly of the mass matrix (storage term for groundwater flow)
   virtual int AssembleMassMatrix();

   /// does assembly of rhs vector from boundary conditions
   virtual int AssembleRhs();

   /// eliminates algebraic equations rows corresponding to fixed potential rows
   virtual void TidyEquations();

public:
   /// Initialize a Fem problem from a given mesh
   SteadyFlow(Project*,
              StrucMesh*,
              BndCondition*,
             _Pcg_Params_ &);

   /// default Fem destructor
   virtual ~SteadyFlow();

   /// Sets the simulation as transient or steady state
   void SetTransient(bool transient = true);

   /// returns the simulation mode (steady state or transient)
   bool IsTransient();

   /// initialize unknowns from an existing solution
   virtual void InitializeVar(double *init);

   /// returns computed variable
   double *GetVar();

   /// returns the RHS of the system (i.e. nodal flow rate array)
   inline double* GetRhs() const {
      return rhs; 
   }

   /// returns groundwater table depth 
   inline virtual double* GetGwDepth() const {
      return gwDepth; 
   }

   /// does finite element solution of the actual problem
   virtual int Solve();

   /// Abstract virtual function
   virtual void SetWaterTableTolerance(double tol) {
   };

    /**
     * Prints final summary and mass balance error to standard output
     * stream at the end of the numerical simulation.
     */
    virtual void WriteFinalSummary();

};

#endif // STEADY_FLOW
