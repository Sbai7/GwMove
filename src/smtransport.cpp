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

// project headers 
#include <smtransport.h>


// ======================================================
SteadyMassTransport::SteadyMassTransport(Project *_project,
                                               StrucMesh *nw_mesh,
                                               double *flowRate,
                                               BndCondition *BC,
                                               _Pcg_Params_ &pcgParams)
{
   mesh = nw_mesh;

   dimension = mesh->GetNnodes();

   // allocate memory for components of algebraic equations
   var = new double[dimension];
   rhs = new double[dimension];
   Q   = new double[dimension];

   // initialization
   unsigned int i;
   for (i = 0; i < dimension; i++) {
      var[i] = 0.;
      rhs[i] = 0.;
      Q[i]   = flowRate[i]; 
   }

   // build a Sparse matrix object
   mat = new CsrMatrix(dimension, mesh->GetConnectivity(), mesh->GetNelements());
//todo   mat = new CsrMatrix(dimension, mesh->connectivity, mesh->GetNelements(), false); // last argument sets symmetry flag (true=sym, false=unsym) 
//temp   matrixCopy = new double[14 * mat->GetDim()];

   // get pointer to boundary conditions
   bc = BC;
//todo   bc = <dynamic_cast> (MassTransportBC) BC;

   // build linear solver (BiCGStab)
//todo   bicgstab = new PbicgstabSolution(mat, rhs, var);
//todo   bicgstab->silent = pcgParams.silent;
//todo   bicgstab->SetMaxIterations(pcgParams.iterations);
//todo   bicgstab->SetTolerance(pcgParams.tolerance);

   // and soil types vector
//AS   Soils = project->soils;

   // by default we consider steady state 
   isTransient = false;
}


// ======================================================
SteadyMassTransport::~SteadyMassTransport()
{
//todo   if (bicgstab) delete bicgstab; 
   if (Q) delete[] Q;
}

// ======================================================
int SteadyMassTransport::AssembleStiffnessMatrix()
{
   // quit in case of invalid pointers
   if (!mat || !mesh) return 0;

   unsigned int ne = mesh->GetNelements();
   unsigned int ie;

   // loop over all elements in the mesh
   for (ie = 0; ie < ne; ie++) {

      // exclude empty elements from the finite element mesh
      if (mesh->soils[ie] == 0) continue;

      // calculate local conductance matrix
      // ...
      // ... 
      //todo   hex.ComputeMassTransportStiffness(); 

      // for every pair of nodes add contribution to global FE matrix
      unsigned int ll, kk, l, j, k, jmin;
      for (ll = 0; ll < 8; ll++) {
         l = mesh->GetNode(ie, ll);
         for (kk = 0; kk < 8; kk++) {

            k = mesh->GetNode(ie, kk);
            if (l == 1) jmin = 0;
            else jmin = mat->pd[l-1]+1;
            for (j = jmin; j < mat->pd[l]; j++) {
               if (k == mat->pc[j]) {
                  //todo   mat->lowerMatrix[j] += hex.matrix[kk][ll]; // contribution to lower matrix
                  //todo   mat->upperMatrix[j] += hex.matrix[ll][kk]; // contribution to upper matrix
               }
            }
         }
      }
   }
}

// ======================================================
int SteadyMassTransport::AssembleRhs()
{
   // quit in case of invalid pointers
   if (!bc || !rhs || !mat || !var) return 0;

   unsigned int i;
   //double sv;

   // consider all boundary conditions types
   for (i = 0; i < dimension; i++) {

//todo      bc->fixedConc[i] = false;
      int j = mat->pd[i];

      switch (bc->bcType[i]) {

      case FIXED_CONC_BC:
      //todo   rhs[i] = mat->lowerMatrix[j] * bc->GetValue(i, 1);
      //todo   bc->fixedConc[i] = true;
         break;

      case FIXED_INPUT_CONC_BC:
         if (Q[i] > 0.) rhs[i] = Q[i] * bc->GetValue(i,1);
         break;

      case FIXED_MASS_BC:
      case FIXED_MASS_FLUX_BC:
         rhs[i] = bc->GetValue(i,2);
         break; 

      case MIXED_BC:
      case MIXED_FLUX_BC:
         //todo   mat->lowerMatriw[j] += bc->GetValue(i,2);
         //todo   mat->upperMatriw[j] += bc->GetValue(i,2);
         if (Q[i] > 0.) rhs[i] = (Q[i] + bc->GetValue(i,2)) * bc->GetValue(i,1);
         else rhs[i] = bc->GetValue(i,2) * bc->GetValue(i,1); 
         break;

      case SINGLE_NODES:
         var[i] = 0.;
         //todo   mat->lowerMatrix[j] = 1.;
         //todo   mat->upperMatrix[j] = 1.;
         rhs[i] = 0.;
         break;

      default:
         rhs[i] = 0.;
         break; 

      }
   }
}

// ======================================================
void SteadyMassTransport::TidyEquations()
{

}

// ======================================================
void SteadyMassTransport::WriteFinalSummary()
{

}
