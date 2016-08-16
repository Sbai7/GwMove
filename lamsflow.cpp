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
#include <lamsflow.h>
#include <vector.h>

// C++ headers
#include <cmath>

using namespace happ;


LamSteadyFlow::LamSteadyFlow(Project *_project,
                            StrucMesh *_mesh,   	   // Mesh pointer
                            BndCondition *_bc,          // Boundary conditions
                           _Pcg_Params_ &pcgParams)   // Solver parameters
: SteadyFlow(_project,_mesh,_bc,pcgParams)
{
    // initialize water table iteration counter
    iteration = 0;

    // make zTopo from upper slice of the initial mesh 
    unsigned int i; 
    unsigned int nij = mesh->GetNi() * mesh->GetNj();
    zTopo.clear();
    for (i = 0; i < nij; i++) {
        zTopo.push_back( mesh->GetZ(i) );
    }

    // initialize vector of flooded nodes 
    for (i = 0; i < nij; i++) {
        IsFlooded.push_back(false);
    }

    // Initialize groundwater depth 
    gwDepth = new double[nij]; 
    for (i = 0; i < nij; i++) gwDepth[i] = 0.0;

    // ...
    iterPlot.Name("Layered adaptive mesh convergence history");
    totalPcgIter = 0;
    wtError = 1.E+99;
}

LamSteadyFlow::~LamSteadyFlow()
{
    zTopo.clear();
    IsFlooded.clear();
    convRate.clear();
} 

void LamSteadyFlow::CheckFloodingBC()
{
   // quit in case of invalid pointers
   if (!bc || !rhs || !mat || !var) return;

   unsigned int i;
   for (i = 0; i < dimension; i++) bc->fixedHead[i] = false;

   // Now, check for flooding nodes @ the topographic surface 
   unsigned int nij = mesh->GetNi() * mesh->GetNj();
   for (i = 0; i < nij; i++) {
      int j = mat->pd[i];
      if (var[i] > zTopo[i] + tolerance) {
         if (bc->bcType[i] > NO_BC && bc->bcType[i] < FLOOD_BC) {
            // this node had peviously a given BC which will be 
            // deactivated during current iteration 
            bc->Deactivate(i);
         }
         rhs[i] = mat->matrix[j] * zTopo[i];
         bc->fixedHead[i] = true;
         IsFlooded[i] = true;
      }
      else {
         // retreive back any desactivated BC if this is the case
         if (bc->bcType[i] > NO_BC && bc->bcType[i] < FLOOD_BC) {
            // activate it again
            bc->Activate(i);
         }
         bc->fixedHead[i] = false;
         IsFlooded[i] = false;
      }
   }
}

int LamSteadyFlow::AssembleRhs()
{
    // quit in case of invalid pointers
    if ( !bc || !rhs || !mat || !var ) return 0;

    unsigned int i;
    double sv;

    // consider all boundary conditions types
    for (i = 0; i < dimension; i++) {
       if (!bc->isActive(i)) continue; 

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
            mat->matrix[j] = 1.;
            var[i] = 0.;
            rhs[i] = 0.;
            break;

         default:
            rhs[i] = 0.;
            break;

      }

   }

   return 1;
}

int LamSteadyFlow::DoSingleLoop()
{
   unsigned int i;
   int ret;

   // update oldZ array
   for (i = 0; i < dimension; i++) oldZ[i] = mesh->GetZ(i);

   // increment iteration count
   iteration++;
   std::cout << "Mesh iteration # " << iteration << ":" << endl;

   // Assemble global conductance matrix
   mat->InitMatrix(); // re-initialize the matrix 
   ret = AssembleStiffnessMatrix();
   if (ret != 1) return ret;

   // Adjusts positive definiteness of conductance matrix
   mat->AdjustPositiveDefiniteness();

   // save a copy of the global conductance matrix
   for (i = 0; i < 14 * dimension; i++) matrixCopy[i] = mat->matrix[i];

   // Check for occurence of water table flooding 
   CheckFloodingBC();

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
   double pcgError;
   ret = pcg->Solve(iterations, pcgError);
   totalPcgIter += iterations;

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
      std::cout << "------------- PCG iterations = " << iterations << " - error = "
      << pcgError << " -------------" << endl;

   // Retrieve new solution
   var = pcg->GetSolution();

   // adjusts moving slices positions
   SetInterfacePosition();

   // calculate iteration error
   wtError = GetConvergenceErrorNorm();
   std::cout << "------------- Water table iteration  error = " << wtError
      << " -------------" << endl << endl;

   double itReal = (double)iteration;
   iterPlot.Add(itReal, wtError);

   return 1;
}

int LamSteadyFlow::Solve()
{
   int ret;

   // set output stream format
   std::cout << setiosflags(ios::scientific);

   // initialize oldZ array
   oldZ = new double [dimension];

   do {
      ret = DoSingleLoop(); 
      if (ret != 1) return ret; 
   } 
   while (wtError > tolerance && iteration < maxIterations);

   iterPlot.Deriv(convRate);

   // Calculate groundwater depth 
   unsigned int i, nij = mesh->GetNi() * mesh->GetNj();
   for (i = 0; i < nij; i++) gwDepth[i] = std::abs( zTopo[i]-var[i] );

   // free memory
   delete[] oldZ;

   return 1;
}

void LamSteadyFlow::SetInterfacePosition()
{
   unsigned int l,i,j,jj;
   unsigned int nij = mesh->GetNi() * mesh->GetNj();
   unsigned int nw  = mesh->GetNmovSlices();

   for (l = 1; l < nw; l++) {
      for (i = 0; i < nij; i++) {
         j  = i+(l-1)*nij;
         jj = i+(nw-1)*nij;
         if (bc->bcType[i] != SINGLE_NODES) {
            double newZ = (var[i]*(nw-l)+mesh->GetZ(jj)*(l-1))/(nw-1);
            mesh->SetZ(j,newZ);
         }
      }
   }
}

double LamSteadyFlow::GetConvergenceErrorNorm()
{
   double error = 0.;
   for (unsigned int i = 0; i < mesh->GetNnodes(); i++) {
      double diff = fabs(oldZ[i]-mesh->GetZ(i));
      if (diff > error) error = diff;
   }
   return error;
}

void LamSteadyFlow::WriteFinalSummary()
{
   unsigned int i;

   /* (1) Output of mesh / water table iteration & convergence summary */

   // mean convergence rate 
   double meanConvRate = 0;
   for (i = 0; i < convRate.size(); i++) meanConvRate += convRate[i];
   meanConvRate /= convRate.size();
   meanConvRate = abs(meanConvRate);

   // write convergence summary table 
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << " Water table convergence summary:                          " << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << " Total number of PCG iterations = " << setw(15) << totalPcgIter << endl;
   std::cout << " Number of mesh iterations      = " << setw(15) << iterPlot.Size() << endl;
   std::cout << " Water table error              = " << setw(15) << wtError << endl;
   std::cout << " Mean convergence rate          = " << setw(15) << meanConvRate << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << endl;


   /* (2) Output of global mass-balance analysis for each BC's type */

   // consider boundary conditions again to get rhs values
   AssembleRhs();

   // restore global conductance matrix
   for (i = 0; i < 14*dimension; i++)
      mat->matrix[i] = matrixCopy[i];

   // calculate nodal flow rates vector from computed heads 
   double *q = new double[dimension];
   mat->Multiply(dimension, var, q);
   Vector Q(q, dimension);

   // calculate mass balance by type of each BC
   double Q_fixedHeadBC = 0, Q_fixedFlowBC = 0, Q_fixedFluxBC   = 0, Q_riverBC   = 0,
          Q_leakageBC   = 0, Q_drainageBC  = 0, Q_abstractionBC = 0, Q_seepageBC = 0,
          Q_rechargeBC  = 0, Q_evtBC       = 0;
   double Qp_fixedHeadBC = 0, Qp_fixedFlowBC = 0, Qp_fixedFluxBC   = 0, Qp_riverBC   = 0,
          Qp_leakageBC   = 0, Qp_drainageBC  = 0, Qp_abstractionBC = 0, Qp_seepageBC = 0,
          Qp_rechargeBC  = 0, Qp_evtBC       = 0;
   double Qm_fixedHeadBC = 0, Qm_fixedFlowBC = 0, Qm_fixedFluxBC   = 0, Qm_riverBC   = 0,
          Qm_leakageBC   = 0, Qm_drainageBC  = 0, Qm_abstractionBC = 0, Qm_seepageBC = 0,
          Qm_rechargeBC  = 0, Qm_evtBC       = 0;
   for (i = 0; i < dimension; i++) {
      if (!bc->isActive(i)) continue; //TODO: printout the list of inactive BC's 

      switch (bc->bcType[i]) {

      case FIXED_HEAD_BC:
         Q_fixedHeadBC += Q(i);
         if (Q(i) > 0) Qp_fixedHeadBC += Q(i);
         if (Q(i) < 0) Qm_fixedHeadBC += Q(i);
         break;

      case FIXED_FLOW_BC:
         Q_fixedFlowBC += Q(i);
         if (Q(i) > 0) Qp_fixedFlowBC += Q(i);
         if (Q(i) < 0) Qm_fixedFlowBC += Q(i);
         break;

      case FIXED_FLUX_BC:
         Q_fixedFluxBC += Q(i);
         if (Q(i) > 0) Qp_fixedFluxBC += Q(i);
         if (Q(i) < 0) Qm_fixedFluxBC += Q(i);
         break;

      case RIVER_BC:
         Q_riverBC += Q(i);
         if (Q(i) > 0) Qp_riverBC += Q(i);
         if (Q(i) < 0) Qm_riverBC += Q(i);
         break;

      case LEAKAGE_BC:
         Q_leakageBC += Q(i);
         if (Q(i) > 0) Qp_leakageBC += Q(i);
         if (Q(i) < 0) Qm_leakageBC += Q(i);
         break;

      case DRAINAGE_BC:
         Q_drainageBC += Q(i);
         if (Q(i) > 0) Qp_drainageBC += Q(i);
         if (Q(i) < 0) Qm_drainageBC += Q(i);
         break;

      case ABSTRACTION_BC:
         Q_abstractionBC += Q(i);
         if (Q(i) > 0) Qp_abstractionBC += Q(i);
         if (Q(i) < 0) Qm_abstractionBC += Q(i);
         break;

      case SEEPAGE_BC:
         Q_seepageBC += Q(i);
         if (Q(i) > 0) Qp_seepageBC += Q(i);
         if (Q(i) < 0) Qm_seepageBC += Q(i);
         break;

      case RECHARGE_BC:
         Q_rechargeBC += Q(i);
         if (Q(i) > 0) Qp_rechargeBC += Q(i);
         if (Q(i) < 0) Qm_rechargeBC += Q(i);
         break;

      case EVAPOTRANSPIRATION_BC:
         Q_evtBC += Q(i);
         if (Q(i) > 0) Qp_evtBC += Q(i);
         if (Q(i) < 0) Qm_evtBC += Q(i);
         break;
      }
   }

   // strings for in/out (i.e. + or -) 
   String s_head,  s_flux, s_flow, s_river, s_leak, s_drain,
          s_abstr, s_seep, s_rech, s_evt;
   s_head  = (Q_fixedHeadBC   > 0) ? "in" : "out"; if (Q_fixedHeadBC == 0)   s_head  = "";
   s_flux  = (Q_fixedFlowBC   > 0) ? "in" : "out"; if (Q_fixedFlowBC == 0)   s_flux  = "";
   s_flow  = (Q_fixedFluxBC   > 0) ? "in" : "out"; if (Q_fixedFluxBC == 0)   s_flow  = "";
   s_river = (Q_riverBC       > 0) ? "in" : "out"; if (Q_riverBC == 0)       s_river = "";
   s_leak  = (Q_leakageBC     > 0) ? "in" : "out"; if (Q_leakageBC == 0)     s_leak  = "";
   s_drain = (Q_drainageBC    > 0) ? "in" : "out"; if (Q_drainageBC == 0)    s_drain = "";
   s_abstr = (Q_abstractionBC > 0) ? "in" : "out"; if (Q_abstractionBC == 0) s_abstr = "";
   s_seep  = (Q_seepageBC     > 0) ? "in" : "out"; if (Q_seepageBC == 0)     s_seep  = "";
   s_rech  = (Q_rechargeBC    > 0) ? "in" : "out"; if (Q_rechargeBC == 0)    s_rech  = "";
   s_evt   = (Q_evtBC         > 0) ? "in" : "out"; if (Q_evtBC == 0)         s_evt   = "";

   double Qp_total = Qp_fixedHeadBC + Qp_fixedFlowBC + Qp_fixedFluxBC   + Qp_riverBC   +
                     Qp_leakageBC   + Qp_drainageBC  + Qp_abstractionBC + Qp_seepageBC +
                     Qp_rechargeBC  + Qp_evtBC;

   double Qm_total = Qm_fixedHeadBC + Qm_fixedFlowBC + Qm_fixedFluxBC   + Qm_riverBC   +
                     Qm_leakageBC   + Qm_drainageBC  + Qm_abstractionBC + Qm_seepageBC +
                     Qm_rechargeBC  + Qm_evtBC;

   double Q_total = Q.Sum();
   double Q_totalAbs = Q.Norml1(); 
   // mass balance error (in %)
   double mbe = 100.*abs(Q_total) / (Q_totalAbs / 2.); 

   // calculate NormL2 of flow rate residual
   double Q_norml2 = 0.;
   for (i = 0; i < dimension; i++) {
      if (bc->fixedHead[i])
         rhs[i] = Q[i];
      else
         Q_norml2 += pow(rhs[i] - Q(i), 2);
   }
   Q_norml2 = sqrt(Q_norml2 / dimension);

   std::cout << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << " Mass balance analysis:                          " << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << " BC TYPE -------------- FLOW RATE ------------- Q_IN ------- Q_OUT -----" << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << " Fixed head nodes     " << setw(15) << Q_fixedHeadBC
                                         << setw(5)  << s_head 
                                         << setw(15) << Qp_fixedHeadBC 
                                         << setw(15) << Qm_fixedHeadBC  
                                         << endl;
   std::cout << " Fixed flow nodes     " << setw(15) << Q_fixedFlowBC   
                                         << setw(5)  << s_flux  
                                         << setw(15) << Qp_fixedFlowBC 
                                         << setw(15) << Qm_fixedFlowBC  
                                         << endl;
   std::cout << " Fixed flux nodes     " << setw(15) << Q_fixedFluxBC   
                                         << setw(5)  << s_flow  
                                         << setw(15) << Qp_fixedFluxBC 
                                         << setw(15) << Qm_fixedFluxBC  
                                         << endl;
   std::cout << " River nodes          " << setw(15) << Q_riverBC       
                                         << setw(5)  << s_river 
                                         << setw(15) << Qp_riverBC
                                         << setw(15) << Qm_riverBC
                                         << endl;
   std::cout << " Leakage nodes        " << setw(15) << Q_leakageBC     
                                         << setw(5)  << s_leak  
                                         << setw(15) << Qp_leakageBC
                                         << setw(15) << Qm_leakageBC
                                         << endl;
   std::cout << " Drainage nodes       " << setw(15) << Q_drainageBC    
                                         << setw(5)  << s_drain 
                                         << setw(15) << Qp_drainageBC
                                         << setw(15) << Qm_drainageBC
                                         << endl;
   std::cout << " Abstraction nodes    " << setw(15) << Q_abstractionBC 
                                         << setw(5)  << s_abstr 
                                         << setw(15) << Qp_abstractionBC
                                         << setw(15) << Qm_abstractionBC
                                         << endl;
   std::cout << " Seepage nodes        " << setw(15) << Q_seepageBC     
                                         << setw(5)  << s_seep  
                                         << setw(15) << Qp_seepageBC
                                         << setw(15) << Qm_seepageBC
                                         << endl;
   std::cout << " Recharge nodes       " << setw(15) << Q_rechargeBC    
                                         << setw(5)  << s_rech  
                                         << setw(15) << Qp_rechargeBC
                                         << setw(15) << Qm_rechargeBC
                                         << endl;
   std::cout << " EVT nodes            " << setw(15) << Q_evtBC         
                                         << setw(5)  << s_evt   
                                         << setw(15) << Qp_evtBC
                                         << setw(15) << Qm_evtBC
                                         << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << " TOTAL (in + out)     " << setw(15) << Q_total 
                                         << setw(5)  << "" 
                                         << setw(15) << Qp_total
                                         << setw(15) << Qm_total
                                         << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << " Residual L2-Norm   = " << setw(15) << Q_norml2 << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << " Mass Balance Error = " << setw(15) << mbe << " %" << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << endl;


   /* (3) Output of miscellaneous informations when surface flooding occurs */

   // detect if flooding occurs ... 
   unsigned int floodingOK = 0; 
   unsigned int size = IsFlooded.size();
   for (i = 0; i < size; i++) {
       if (IsFlooded[i]) floodingOK++;
   }
   if (floodingOK) {
       std::cout << "-------------------------------------------------" << endl;
       std::cout << "Flooding has occured through the following nodes:" << endl;
       std::cout << "-------------------------------------------------" << endl;
       double totalFloodQ = 0.; 
       for (i = 0; i < size; i++) {
           if (IsFlooded[i]) {
               std::cout << " + node " << setw(8) << i << "Outflow " << setw(15) <<
                   setprecision(8) << std::abs(rhs[i]) << endl; 
               totalFloodQ += std::abs(rhs[i]);
           }
       }
       std::cout << "-------------------------------------------------" << endl;
       std::cout << "Total flow rate = " << setw(31) << setprecision(8) <<
           totalFloodQ << endl; 
   }
   else {
       std::cout << "Flooding of topographic nodes has not been detected." << endl;
   }
   std::cout << endl;

    // Ouput user-selected MBA 
   UserMBA(q); 

   /* (4) Output of global groundwater variables statistics */

   // Output user-selected STATA 
 //AS-Debug  UserStata(var, gwDepth, q);


   // extract flow rates @ BC nodes to erase small values @ internal nodes
   bc->ExtractAtBC(Q);

   unsigned int nij = mesh->GetNi()*mesh->GetNj();
   Vector depth(gwDepth, nij);
   Vector head(var, dimension);
   double hPerc[11], dPerc[11], qPerc[11];
   head.Percentiles(hPerc); 
   depth.Percentiles(dPerc); 
   Q.Percentiles(qPerc);

   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << " Groundwater variables statistics:                          " << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << " Groundwater:   head --------- depth -------- Flow rate ---------------- " << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << " Max:     " << setw(15) << head.Max() 
                             << setw(15) << depth.Max() 
                             << setw(15) << Q.Max() 
                             << endl;
   std::cout << " Min:     " << setw(15) << head.Min() 
                             << setw(15) << depth.Min() 
                             << setw(15) << Q.Min() 
                             << endl;
   std::cout << " Median:  " << setw(15) << head.Mean() 
                             << setw(15) << depth.Mean() 
                             << setw(15) << Q.Mean() 
                             << endl;
   std::cout << " Variance:" << setw(15) << head.Variance() 
                             << setw(15) << depth.Variance() 
                             << setw(15) << Q.Variance() 
                             << endl;
   std::cout << " Std-Dev: " << setw(15) << head.StandardDev() 
                             << setw(15) << depth.StandardDev() 
                             << setw(15) << Q.StandardDev() 
                             << endl;
   std::cout << " Skewness:" << setw(15) << head.Skewness() 
                             << setw(15) << depth.Skewness()
                             << setw(15) << Q.Skewness()
                             << endl;
   std::cout << " Kurtosis:" << setw(15) << head.Kurtosis()
                             << setw(15) << depth.Kurtosis()
                             << setw(15) << Q.Kurtosis()
                             << endl;
   std::cout << " 10% pctl:" << setw(15) << hPerc[1]
                             << setw(15) << dPerc[1]
                             << setw(15) << qPerc[1]
                             << endl;
   std::cout << " 20% pctl:" << setw(15) << hPerc[2]
                             << setw(15) << dPerc[2]
                             << setw(15) << qPerc[2]
                             << endl;
   std::cout << " 30% pctl:" << setw(15) << hPerc[3]
                             << setw(15) << dPerc[3]
                             << setw(15) << qPerc[3]
                             << endl;
   std::cout << " 40% pctl:" << setw(15) << hPerc[4]
                             << setw(15) << dPerc[4]
                             << setw(15) << qPerc[4]
                             << endl;
   std::cout << " 50% pctl:" << setw(15) << hPerc[5]
                             << setw(15) << dPerc[5]
                             << setw(15) << qPerc[5]
                             << endl;
   std::cout << " 60% pctl:" << setw(15) << hPerc[6]
                             << setw(15) << dPerc[6]
                             << setw(15) << qPerc[6]
                             << endl;
   std::cout << " 70% pctl:" << setw(15) << hPerc[7]
                             << setw(15) << dPerc[7]
                             << setw(15) << qPerc[7]
                             << endl;
   std::cout << " 80% pctl:" << setw(15) << hPerc[8]
                             << setw(15) << dPerc[8]
                             << setw(15) << qPerc[8]
                             << endl;
   std::cout << " 90% pctl:" << setw(15) << hPerc[9]
                             << setw(15) << dPerc[9]
                             << setw(15) << qPerc[9]
                             << endl;
   std::cout << "------------------------------------------------------------------------" << endl;
   std::cout << endl;



   // free memory
   Q.Destroy();
}

void LamSteadyFlow::UserMBA(const double *flowRate) const
{
   filebuf outfile;
   String mbaFile = project->GetFullName() + ".mba";
   outfile.open(mbaFile.GetBuffer(), ios::out);
   ostream os(&outfile);
   os << setiosflags(ios::scientific);

   unsigned int i, n;
   BcZone *bcZn;
   double Q, Qplus, Qminus; 
   double Q_total = 0, Qp_total = 0, Qm_total = 0;
   
   /* (1) mass-balance analysis for each BC zone */
   if (project->mba.bcMBA > 0) {

      os << "------------------------------------------------------------------------" << endl;
      os << " Mass balance analysis by boundary conditions zone:                     " << endl;
      os << "------------------------------------------------------------------------" << endl;
      os << " BC ZONE -------------- FLOW RATE ------------- Q_IN ------- Q_OUT -----" << endl;
      os << "------------------------------------------------------------------------" << endl;

      n = bc->bcZone.size();
      for (i = 0; i < n; i++) {
         bcZn = bc->bcZone[i];
         bcZn->ExtractVar(flowRate, mesh);
         Q = bcZn->SumVar();
         Qplus = bcZn->PositiveVar();
         Qminus = bcZn->NegativeVar();
         Q_total += Q;
         Qp_total += Qplus;
         Qm_total += Qminus;

         String label = bcZn->Name();
         String s_bc = (Q > 0) ? "in" : "out"; if (Q == 0) s_bc = "";

         os << setw(22) << label << setw(15) << Q
            << setw(5) << s_bc
            << setw(15) << Qplus
            << setw(15) << Qminus
            << endl;

      }
      os << "------------------------------------------------------------------------" << endl;
      os << " TOTAL (in + out)     " << setw(15) << Q_total
         << setw(5) << ""
         << setw(15) << Qp_total
         << setw(15) << Qm_total
         << endl;

      os << "------------------------------------------------------------------------" << endl;
      os << endl;
   }

   /* (2) mass-balance analysis for each user-specified zone */  
   ijkZone *ijkZn;
   Q_total = 0, Qp_total = 0, Qm_total = 0;
   n = project->mba.userZone.size(); 
   if (n > 0) {

      os << endl;
      os << "------------------------------------------------------------------------" << endl;
      os << " Mass balance analysis by user-specified zones:                         " << endl;
      os << "------------------------------------------------------------------------" << endl;
      os << " BC ZONE -------------- FLOW RATE ------------- Q_IN ------- Q_OUT -----" << endl;
      os << "------------------------------------------------------------------------" << endl;

      for (i = 0; i < n; i++) {
         ijkZn = project->mba.userZone[i];
         ijkZn->ExtractVar(flowRate, mesh);
         Q = ijkZn->SumVar();
         Qplus = ijkZn->PositiveVar();
         Qminus = ijkZn->NegativeVar();
         Q_total += Q;
         Qp_total += Qplus;
         Qm_total += Qminus;

         String label = project->mba.label[i];
         String s_zn = (Q > 0) ? "in" : "out"; if (Q == 0) s_zn = "";

         os << setw(22) << label << setw(15) << Q
            << setw(5) << s_zn
            << setw(15) << Qplus
            << setw(15) << Qminus
            << endl;

      }
      os << "------------------------------------------------------------------------" << endl;
      os << " TOTAL (in + out)     " << setw(15) << Q_total
         << setw(5) << ""
         << setw(15) << Qp_total
         << setw(15) << Qm_total
         << endl;

      os << "------------------------------------------------------------------------" << endl;
      os << endl;
   }

   outfile.close();
}

void LamSteadyFlow::UserStata(const double *Head,
   const double *Depth,
   const double *FlowRate) const
{
   filebuf outfile;
   String stataFile = project->GetFullName() + ".sta";
   outfile.open(stataFile.GetBuffer(), ios::out);
   ostream os(&outfile);
   os << setiosflags(ios::scientific);


   /* statistical analysis for each user-specified zone */
   unsigned int i, n;
   ijkZone *ijkZn;
   double hPerc[11], dPerc[11], qPerc[11];
   unsigned int nij = mesh->GetNi()*mesh->GetNj();


   n = project->stata.userZone.size();
   if (n > 0) {

      os << endl;
      os << "------------------------------------------------------------------------" << endl;
      os << " Groundwater variables statistics by user-specified zones:               " << endl;
      os << "------------------------------------------------------------------------" << endl;
      os << " Groundwater:   head --------- depth -------- Flow rate ---------------- " << endl;
      os << "------------------------------------------------------------------------" << endl;

      for (i = 0; i < n; i++) {
         ijkZn = project->stata.userZone[i];

         // groundwater head 
         ijkZn->ExtractVar(Head, mesh);
         std::vector<double> h = ijkZn->GetBCVar();
         Vector head(h.data(), dimension);
         head.Percentiles(hPerc);

         // groundwater depth 
         ijkZn->ExtractVar(Depth, mesh, false);
         std::vector<double> d = ijkZn->GetBCVar();
         Vector depth(d.data(), nij);
         depth.Percentiles(dPerc);

         // flow rate 
         ijkZn->ExtractVar(FlowRate, mesh);
         std::vector<double> q = ijkZn->GetBCVar(); 
         Vector Q(q.data(), dimension);
         Q.Percentiles(qPerc);
         
         String label = project->stata.label[i];

         os << "Zone = " << setw(15) << label << endl; 
         os << " Max:     " << setw(15) << head.Max()
            << setw(15) << depth.Max()
            << setw(15) << Q.Max()
            << endl;
         os << " Min:     " << setw(15) << head.Min()
            << setw(15) << depth.Min()
            << setw(15) << Q.Min()
            << endl;
         os << " Median:  " << setw(15) << head.Mean()
            << setw(15) << depth.Mean()
            << setw(15) << Q.Mean()
            << endl;
         os << " Variance:" << setw(15) << head.Variance()
            << setw(15) << depth.Variance()
            << setw(15) << Q.Variance()
            << endl;
         os << " Std-Dev: " << setw(15) << head.StandardDev()
            << setw(15) << depth.StandardDev()
            << setw(15) << Q.StandardDev()
            << endl;
         os << " Skewness:" << setw(15) << head.Skewness()
            << setw(15) << depth.Skewness()
            << setw(15) << Q.Skewness()
            << endl;
         os << " Kurtosis:" << setw(15) << head.Kurtosis()
            << setw(15) << depth.Kurtosis()
            << setw(15) << Q.Kurtosis()
            << endl;
         os << " 10% pctl:" << setw(15) << hPerc[1]
            << setw(15) << dPerc[1]
            << setw(15) << qPerc[1]
            << endl;
         os << " 20% pctl:" << setw(15) << hPerc[2]
            << setw(15) << dPerc[2]
            << setw(15) << qPerc[2]
            << endl;
         os << " 30% pctl:" << setw(15) << hPerc[3]
            << setw(15) << dPerc[3]
            << setw(15) << qPerc[3]
            << endl;
         os << " 40% pctl:" << setw(15) << hPerc[4]
            << setw(15) << dPerc[4]
            << setw(15) << qPerc[4]
            << endl;
         os << " 50% pctl:" << setw(15) << hPerc[5]
            << setw(15) << dPerc[5]
            << setw(15) << qPerc[5]
            << endl;
         os << " 60% pctl:" << setw(15) << hPerc[6]
            << setw(15) << dPerc[6]
            << setw(15) << qPerc[6]
            << endl;
         os << " 70% pctl:" << setw(15) << hPerc[7]
            << setw(15) << dPerc[7]
            << setw(15) << qPerc[7]
            << endl;
         os << " 80% pctl:" << setw(15) << hPerc[8]
            << setw(15) << dPerc[8]
            << setw(15) << qPerc[8]
            << endl;
         os << " 90% pctl:" << setw(15) << hPerc[9]
            << setw(15) << dPerc[9]
            << setw(15) << qPerc[9]
            << endl;
         os << "------------------------------------------------------------------------" << endl;
      }
      os << endl;
   }

   outfile.close();
}
