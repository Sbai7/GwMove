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

// Project headers
#include <sflow.h>

// C++ headers
#include <cmath>


SteadyFlow::SteadyFlow(Project *_project,
                       StrucMesh *_mesh,
                       BndCondition *_bc,
                      _Pcg_Params_ &pcgParams)
{
   // pointer to project
   project = _project; 

   // mesh pointer 
   mesh = _mesh;

   dimension = mesh->GetNnodes();

   // allocate memory for components of algebraic equations
   var = new double [dimension];
   rhs = new double [dimension];

   // initialization
   unsigned int i;
   for (i = 0; i < dimension; i++) {
      var[i] = 0.;
      rhs[i] = 0.;
   }

   // build a Sparse matrix object
   mat = new CsrMatrix(dimension, mesh->GetConnectivity(), mesh->GetNelements());
   matrixCopy = new double [14*mat->GetDim()];

   // get pointer to boundary conditions
   bc = _bc;

   // build linear solver (PCG)
   pcg = new PcgSolution(mat,rhs,var);
   pcg->silent = pcgParams.silent;
   pcg->SetMaxIterations( pcgParams.iterations );
   pcg->SetTolerance( pcgParams.tolerance );

   // by default we consider steady state 
   isTransient = false;

   // groundwater depth 
   gwDepth = NULL;
}

SteadyFlow::~SteadyFlow()
{
   // free memory
   if (var)  delete[] var;
   if (rhs)  delete[] rhs;
   if (gwDepth) delete[] gwDepth;
   if (isTransient) {if (mass) delete[] mass;}
   if (matrixCopy) delete[] matrixCopy;
   if (mat)  delete mat;
   if (pcg)  delete pcg;
}

void SteadyFlow::SetTransient(bool transient)
{
   isTransient = transient;
   
   unsigned int i;
   if ( isTransient && !mass) {
      // allocate memory for mass matrix if not done yet
      mass = new double [dimension];
      for (i = 0; i < dimension; mass[i++] = 0);
   }
}

bool SteadyFlow::IsTransient()
{
   return isTransient;
}

int SteadyFlow::AssembleStiffnessMatrix()
{
   // quit in case of invalid pointers
   if ( !mat || !mesh ) return 0;

   unsigned int ne = mesh->GetNelements();
   unsigned int ie;

   // loop over all elements in the mesh
   for (ie = 0; ie < ne; ie++) {

      // exclude empty elements from the finite element mesh
      if (mesh->soils[ie] == 0 ) continue;

      // calculate local conductance matrix
      unsigned int connec[8];
      double coord[3][8];
      double K[3]; // hydraulic conductivities of element ie
      mesh->GetConnectivity(ie,connec);
      mesh->GetNeighboorNodesCoord(ie,coord);
      Soil soil = project->soils[mesh->soils[ie] - 1];
      soil.GetK(K);
      Hex8 hex(K,coord,connec);
      hex.dimension = 3; // Consider Hexahedral elements
      hex.ComputeStiffness();

      // for every pair of nodes add contribution to global FE matrix
      unsigned int ll, kk, l, j, k;
      for (ll = 0; ll < 8; ll++) {
         l = mesh->GetNode(ie,ll);
         for (kk = 0; kk < 8; kk++) {

            k = mesh->GetNode(ie,kk);
            if (k > l)
               continue; // exclude strictly upper triangular elements
            else if (k < l) {
               for (j = mat->pd[l-1]+1; j < mat->pd[l]+1; j++) {
                  if (k == mat->pc[j]) break;
               }
               mat->matrix[j] += hex.matrix[kk][ll];
            }
            else {
               j = mat->pd[l];
               mat->matrix[j] += hex.matrix[kk][ll];
            }

         }
      }

   }

   return 1;
}

int SteadyFlow::AssembleMassMatrix()
{
   // quit in case of invalid pointers
   if ( !mass || !mesh ) return 0;

   unsigned int ne = mesh->GetNelements();
   unsigned int ie;

   // loop over all elements in the mesh
   for (ie = 0; ie < ne; ie++) {

      // exclude empty elements from the finite element mesh
      if (mesh->soils[ie] == 0 ) continue;

      // calculate local mass matrix
      unsigned int connec[8];
      double coord[3][8];
      double K[3]; // hydraulic conductivities of element ie
      mesh->GetConnectivity(ie,connec);
      mesh->GetNeighboorNodesCoord(ie,coord);
      Soil soil = project->soils[mesh->soils[ie] - 1];
      double s0 = soil.GetS0();
      double satWC = soil.GetSatWC();       // saturated water content
      double resWC = soil.GetResidualWC();  // residual  water content
      soil.GetK(K);
      Hex8 hex(K,coord,connec);
      hex.dimension = 3; // Consider Hexahedral elements
      hex.ComputeMass();

      // for each node add contribution to mass matrix
      unsigned int ll, l;
      for (ll = 0; ll < 8; ll++) {

         l = mesh->GetNode(ie,ll);
         
         // calculate storage of top layer nodes (or water table nodes)
         if ( l <= mesh->GetNi()*mesh->GetNj() ) {
            // add contributions from water table elastic storage term
            // and saturated zone storativity term
            mass[l] += (satWC-resWC) * mesh->surface[l] + s0 * hex.mass[ll];
         }
         else {
            // calculate storage of non water-table nodes
            mass[l] += s0 * hex.mass[ll];
         }
      }
   }

   return 1;
}

int SteadyFlow::AssembleRhs()
{
   // quit in case of invalid pointers
   if ( !bc || !rhs || !mat || !var ) return 0;

   unsigned int i;
   double sv;

   // consider all boundary conditions types
   for (i = 0; i < dimension; i++) {

      bc->fixedHead[i] = false;
      int j = mat->pd[i];

      switch (bc->bcType[i]) {

         case FIXED_HEAD_BC:
            rhs[i] = mat->matrix[j] * bc->GetValue(i,1);
            bc->fixedHead[i] = true;
            break;

         case FIXED_FLOW_BC:
            rhs[i] = bc->GetValue(i,1);
            break;

         case FIXED_FLUX_BC:
            rhs[i] = bc->GetValue(i,1);
            break;

         case RIVER_BC:
            mat->matrix[j] += bc->GetValue(i,2);
            rhs[i] = bc->GetValue(i,1) * bc->GetValue(i,2);
            break;

         case LEAKAGE_BC:
            mat->matrix[j] += bc->GetValue(i,2);
            rhs[i] = bc->GetValue(i,1) * bc->GetValue(i,2);
            break;

         case INFILTRATION_BC:
            // change boundary condition to fixed flux
            bc->SetType(i,FIXED_FLUX_BC);
            rhs[i] = bc->GetValue(i,1);
            break;

         case DRAINAGE_BC:
            if (var[i] <= bc->GetValue(i,1))
               rhs[i] = 0.;
            else {
               rhs[i] = mat->matrix[j] * bc->GetValue(i,1);
               bc->fixedHead[i] = true;
            }
            break;

         case ABSTRACTION_BC:
            if (var[i] >= bc->GetValue(i,1))
               rhs[i] = bc->GetValue(i,2);
            else
               rhs[i] = 0.;
            break;

         case SEEPAGE_BC:
            if (var[i] <= mesh->GetZ(i))
               rhs[i] = mat->matrix[j] * mesh->GetZ(i);
            else
               rhs[i] = 0.;
            break;

         case RECHARGE_BC:
            if (var[i] > bc->GetValue(i,1))
               rhs[i] = bc->GetValue(i,2);
            else if (var[i] > bc->GetValue(i,3))
               rhs[i] = (bc->GetValue(i,2)*(var[i]-bc->GetValue(i,3)) -
                         bc->GetValue(i,4)*(var[i]-bc->GetValue(i,1))) /
                         (bc->GetValue(i,1) - bc->GetValue(i,3));
            else
               rhs[i] = bc->GetValue(i,4);
            break;

         case EVAPOTRANSPIRATION_BC:
            sv = mesh->GetZ(i) - var[i];
            if (sv > bc->GetValue(i,1))
               rhs[i] = bc->GetValue(i,2);
            else if (sv > bc->GetValue(i,3))
               rhs[i] = (bc->GetValue(i,2)*(sv-bc->GetValue(i,3)) -
                         bc->GetValue(i,4)*(sv-bc->GetValue(i,1))) /
                         (bc->GetValue(i,1) - bc->GetValue(i,3));
            else
               rhs[i] = bc->GetValue(i,4);
            break;

         case SINGLE_NODES:
            var[i] = 0.;
            mat->matrix[j] = 1.;
            rhs[i] = 0.;
            break;

         default:
            rhs[i] = 0.;
            break;

      }

   }

   return 1;
}

void SteadyFlow::TidyEquations()
{
   unsigned int i,j,k;
   
   for (i = 1; i < dimension; i++) {
      for (k = mat->pd[i-1]+1; k < mat->pd[i]; k++) {
         j = mat->pc[k];
         if (bc->fixedHead[i]) {
            if (!bc->fixedHead[j]) rhs[j] -= mat->matrix[k] * bc->GetValue(i,1);
            mat->matrix[k] = 0.;
         }
         else {
            if (bc->fixedHead[j]) {
               rhs[i] -= mat->matrix[k] * bc->GetValue(j,1);
               mat->matrix[k] = 0.;
            }
         }
      }
   }
}

void SteadyFlow::InitializeVar(double *init)
{
   for (unsigned int i = 0; i < dimension; i++)
      var[i] = init[i];
}

double* SteadyFlow::GetVar()
{
   return var;
}

int SteadyFlow::Solve()
{
   int ret;
   unsigned int i;

   // set output stream format
   std::cout << setiosflags(ios::scientific);

    double error = 1.E+20;

    // Assemble global conductance matrix
    mat->InitMatrix(); // re-initialize the matrix
    ret = AssembleStiffnessMatrix();
    if (ret != 1) return ret;

    // Adjusts positive definiteness of conductance matrix
    mat->AdjustPositiveDefiniteness();

    // save a copy of the global conductance matrix
    for (i = 0; i < 14*dimension; i++)
        matrixCopy[i] = mat->matrix[i];

    // Set RHS from boundary conditions
    ret = AssembleRhs();
    if (ret != 1) return ret;

    // Tidy finite element equations
    TidyEquations();

    // Solve finite element equations
    pcg->SetMatrix(mat);
    pcg->SetRhs(rhs);
    pcg->SetUnknowns(var);

    unsigned int iterations;
    ret = pcg->Solve(iterations, error);
    switch (ret) {
        case ERROR_ILU0:
            std::cerr << " >> Error: incomplete factorization of global conductance matrix."
                         << endl;
            break;
        case ERROR_MAX_ITERATIONS:
            std::cerr << " >> Error: max. number of PCG iterations excedded."
                         << endl;
            break;
    }
    if (ret != 1) return ret;
    else
    std::cout << "*************** PCG Iterations = " << iterations << " - Error = "
                 << error << " ****************" << endl;

    // Retrieve new solution
    var = pcg->GetSolution();

    // calculate iteration error
   //error = GetPressureHeadError();

      //std::cout << "*************** Water table iterations error = " << error
      //          << " ****************" << endl << endl;

    return 1;

}

void SteadyFlow::WriteFinalSummary()
{
   unsigned int i;

   // consider boundary conditions again to get rhs values
   AssembleRhs();

   // restore global conductance matrix
   for (i = 0; i < 14*dimension; i++)
      mat->matrix[i] = matrixCopy[i];


   // calculate mass balance error and error-norm (sum of squared residuals)
   double *error = new double [dimension];
   mat->Multiply(dimension, var, error);

   double ssq = 0.;
   double total = 0.;
   double totalAbs = 0.;
   for (i = 0; i < dimension; i++) {
      if (bc->fixedHead[i])
         rhs[i] = error[i];
      else
         ssq += pow( rhs[i]-error[i], 2 );
      totalAbs += abs(error[i]);
      total += error[i];
   }

   std::cout << endl << "Mass Balance Error = " << 100.*abs(total)/(totalAbs/2.) 
             << "%" << endl;
   std::cout << "Squared Residuals Error Norm = " << sqrt(ssq/dimension) << endl;

   // free memory
   delete[] error;
}
