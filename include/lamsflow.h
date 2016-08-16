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

#ifndef LAM_STEADY_FLOW
#define LAM_STEADY_FLOW

// Project headers
#include <sflow.h>
#include <function.h>


/**
 * A Layered Adaptive Meshing class for modelling unconfined ground-water flow
 * in subsurface aquifers.
 * 
 * Todo: 
 * - Add support for salt-fresh water interface (in another derived class) 
 * - Add support for outflow sea face (SGD) boundary condition 
 */
class LamSteadyFlow : public SteadyFlow
{
private:

   /// logical vector of flooded topography nodes during water table rise 
   std::vector<bool> IsFlooded; 

   /// total number of PCG iterations 
   unsigned int totalPcgIter; 

protected:

    /// tolerance or stopping criterion for water table iteration
   double tolerance;

   /// water table iteration number
   unsigned int iteration;

    /// water table error 
   double wtError; 

   /// temporary array for moving water table iterations
   double* oldZ;

   /// maximum allowed number of mesh iterations
   unsigned int maxIterations = 100;

   /// vector of Z coordinates @ nodes in the topographic surface 
   std::vector<double> zTopo; 

   /// function where the mesh iterations are stored to plot later
   TimeFunction iterPlot; 

   /// water table convergence rate array 
   vector<double> convRate;

   /**
   * Checks for top surface nodes flooded by upward water movement
   * following an eventual water table rise. Then, adjusts
   * the list of 'active' boundary conditions accordingly.
   */
   void CheckFloodingBC();

   /// Perform a single loop of the mesh / water-table iteration
   int DoSingleLoop();

   /// User-specified detailed mass balance analysis 
   virtual void UserMBA(const double *flowRate) const; 

   /// User-specified detailed statistical analysis 
   virtual void UserStata(const double *head, 
                          const double *depth, 
                          const double *flowRate) const;

public:

   /// Initialize a Fem problem from a given mesh
   LamSteadyFlow(Project*,
                 StrucMesh*,
                 BndCondition*,
                 _Pcg_Params_&);

   /// default Fem destructor
   virtual ~LamSteadyFlow();

   /// sets water table tolerance
   inline void SetWaterTableTolerance(double tol) {
      tolerance = tol;
   }

   /// returns water table tolerance
   inline double GetWaterTableTolerance() {
      return tolerance;
   }

   /**
    * Adjusts the height of nodes lying on the water table slice,
    * other slices heights are found by linear interpolation between
    * the most upper slice and the lowest moving slice in the layered
    * three-dimensional mesh.
    */
   virtual void SetInterfacePosition();

   /// returns water table error-norm between two successive iterations
   double GetConvergenceErrorNorm();

   // Overloaded methods

   /// Assembly of RHS array for specialized boundary conditions to this class.
   virtual int AssembleRhs();

   /// The conforming finite element solver for steady unconfined flow using a
   /// layared adaptive mesh technique.
   virtual int Solve();

   /**
    * Prints final summary and mass balance error to standard output
    * stream at the end of the numerical simulation.
    */
   virtual void WriteFinalSummary();
};

#endif // LAM_STEADY_FLOW
