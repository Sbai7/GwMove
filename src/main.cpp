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
#include <project.h>
#include <timer.h>
#include <sflow.h>
#include <lamsflow.h>
#include <periods.h>
#include <velocity.h>
#include <tecplot.h>

// C++ headers
#include <list>
using namespace std;
using namespace happ; 


// Prints the program's logo to the given output stream
void display_banner(ostream & os);

int main(int argc, char* argv[])
{
   String projectFile;
   Timer timer;

   display_banner(cout); 

   if (argc < 2) {
      std::cerr << "Error: Unspecified project filename at command line argument." 
                << endl;
      // print program usage to screen, such that the user gets informed
      return 0;
   }


   /**************************************************************************/
   /* Section 1: Project construction from input file(s) + FE mesh setup     */
   /*          + BC's object setup                                           */
   /**************************************************************************/

   std::cout << "Projet file = " << argv[1] << endl << endl;
   // take project name as program argument ...
    projectFile = String(argv[1]);
   
    // initialize a new project object
   Project project;
   
   // set project object name
   project.SetFullName(projectFile);

   // starts timing the whole program
   timer.Start();

   // read main data file where we have problem parameters
   if (!project.Read())
      // error in reading the input data file
      return 0;

   // build a new mesh object and sets its dimension
   StrucMesh *mesh = new StrucMesh();
   mesh->SetDim(3);
   
    // load mesh data structure from the file (<project>.fem)
   mesh->Read(projectFile);

   // set number of moving layers in the mesh
   _Wt_Params_ wtParams = project.GetWaterTableParams();
   mesh->SetNmovSlices( wtParams.nMoving );

   // at this point do some model consistency checks:
   // check that all soil type id's defined in mesh file exists
   // in the data file, and that soils parameters seems to be logic
   std::list <unsigned int> soilId;
   if ( !mesh->GetListOfSoilId(soilId) ) {
      // print error message on invalid soil id's in the mesh file
      return 0;
   }
   if ( !project.CheckSoilIds(soilId) ) {
      // print error message for soil type incompatibility
      return 0;
   }

   // build boundary conditions object
   BndCondition *bc = new BndCondition(mesh);
   
   // read boundary conditions file
   bc->Read(projectFile);


   /**********************************************************/
   /* Section 2: New instance of groundwater model object +  */ 
   /*          Iterative solution phase setup                */
   /**********************************************************/

   // get user choice of PCG related parameters
   _Pcg_Params_ pcgParams = project.GetPcgParams();

    // build a new groundwater flow instance
    SteadyFlow *flow;
   
    if (project.GetFlow() == UNCONFINED_FLOW)
        flow = new LamSteadyFlow(&project, mesh, bc, pcgParams);
    else
        flow = new SteadyFlow(&project, mesh, bc, pcgParams);


    // initialize potential heads (if necessary in transient mode)
    // otherwise, Fem class initialize them to zero

    // set water table tolerance, max iterations
    flow->SetWaterTableTolerance( wtParams.tolerance );

    // call Solve method to iteratively find water table position
    if ( flow->Solve() != 1) {
        // analyse error type ...
        return 0;
    }


    /*********************************************************/
    /* Section 3: Printing global and user-selected detailed */
    /*          mass balance analysis reports                */
    /*********************************************************/

    // prints final summary of the simulation, such as mass balance error
    flow->WriteFinalSummary();

    // if the model had converged ... compute velocities
    NodalVelocity *vel = new NodalVelocity(mesh, project.soils, flow->GetVar());
    bool velOk = vel->Compute();


    /*****************************************************************/
    /* Section 4: Writing Output in native Tecplot ASCII file format */
    /*****************************************************************/

    // if convergence is attained write solution output to a tecplot file
    TecplotIO *tec = new TecplotIO( mesh );
    String tecFile = projectFile + ".plt";
    tec->SetFilename(tecFile);
    tec->SetTitle("Generated by GwMove groundwater modelling software");
    tec->SetZoneTitle("Steady-state solution");

   /* add groundwater variables */
   
   // 4.1. Simulated groundwater heads & flows 
   double* h = flow->GetVar();
   tec->AddDataset(h, "head");

   double* Q = flow->GetRhs();
   tec->AddDataset(Q, "flow rate");

   ////// 4.2. Pressure heads ( = h-z) which is exactly 0 @ the water-table 
   //////      interface(s) and positive in all other nodes 
   ////unsigned int i, nn = mesh->GetNnodes();
   ////double* z; mesh->GetZ(z);
   ////double* pressure = new double[nn]; 
   ////for (i = 0; i < nn; i++) pressure[i] = h[i] - z[i];
   ////tec->AddDataset(pressure, "pressure");

   ////// 4.3. Simulated groundwater velocity components 
   ////if ( velOk ) {
   ////   tec->AddDataset(vel->Xvelocity, "Vx");
   ////   tec->AddDataset(vel->Yvelocity, "Vy");
   ////   tec->AddDataset(vel->Zvelocity, "Vz");
   ////}

   ////// 4.4. Simulated groundwater velocity norm 
   ////double* vel_norm = new double[nn];
   ////for (i = 0; i < nn; i++) {
   ////   vel_norm[i] = sqrt(vel->Xvelocity[i] * vel->Xvelocity[i] +
   ////                      vel->Yvelocity[i] * vel->Yvelocity[i] +
   ////                      vel->Zvelocity[i] * vel->Zvelocity[i]);
   ////}
   ////tec->AddDataset(vel_norm, "V_norm");

   ///tec->SetDataPacking(DATAPACKING::BLOCK);

   // does streaming output
   if ( tec->Write() )
      std::cout << endl << "Output have been written to '" << tecFile.GetBuffer() 
                << "' file." << endl;
   else
      std::cout << endl << " >> Error: Problem in writting output to '" 
                << tecFile.GetBuffer() << "' file." << endl;

   tec->Dispose();

   ////////////////////////////////////
   // Temporary code: Write water file. 
   String watFile = projectFile + ".w00";
   tec->WriteWaterFile(watFile);
   ////////////////////////////////////

   ////////////////////////////////////////////////
   // Temporary code: Write groundwater depth file. 
   std::filebuf gwDepthFile;
   String gwFile = projectFile + ".d00";
   gwDepthFile.open(gwFile.GetBuffer(), ios::out);
   ostream os(&gwDepthFile);
   // number of nodes 
   unsigned int nn = mesh->GetNnodes();
   unsigned int i, nij = mesh->GetNi()*mesh->GetNj();
   os << setw(8) << nij << "   # Number of nodes" << endl;
   os << "# Variable 1 = X-coordinate [L]" << endl;
   os << "# Variable 2 = Y-coordinate [L]" << endl;
   os << "# Variable 3 = Z-coordinate [L]" << endl;
   os << "# Variable 4 = Groundwater depth [L]" << endl;
   os << scientific;
   double* D = flow->GetGwDepth(); 
   for (i = 0; i < nij; i++) {
      Tpoint3D pt; 
      mesh->GetCoord(i, pt);
      os << setw(16) << pt.x 
         << setw(16) << pt.y 
         << setw(16) << pt.z
         << setw(16) << D[i]
         << endl;
   }
   gwDepthFile.close();
   ////////////////////////////////////////////////

   ////delete[] pressure;
   ////delete[] vel_norm; 

   /****************************************************************/
   /* Section 5: Free memory of all model objects + stop the timer */
   /****************************************************************/

   if (tec) delete tec;
   if (flow) delete flow;
   if (bc)   delete bc;
   if (mesh) delete mesh;

   // stop timing the program and report cpu
   timer.Stop();
   std::cout << "Total execution time = " << timer.ClockTime() << " (seconds)" << endl;

    return 0;
}

// Prints the Volta ASCII logo to output stream
void display_banner(ostream & os)
{
   os << "__________________________________________________________________" << endl
      << "                                                                  " << endl
      << "                                                                  " << endl
      << " .d8888b.                888b     d888                            " << endl
      << "d88P  Y88b               8888b   d8888                            " << endl
      << "888    888               88888b.d88888                            " << endl
      << "888        888  888  888 888Y88888P888  .d88b.  888  888  .d88b.  " << endl
      << "888  88888 888  888  888 888 Y888P 888 d88**88b 888  888 d8P  Y8b " << endl
      << "888    888 888  888  888 888  Y8P  888 888  888 Y88  88P 88888888 " << endl 
      << "Y88b  d88P Y88b 888 d88P 888   *   888 Y88..88P  Y8bd8P  Y8b.     " << endl
      << " *Y8888P88   *8888888P*  888       888  *Y88P*    Y88P    *Y8888  " << endl
      << "                                                                  " << endl
      << "                                                      Version 1.0 " << endl
      << "__________________________________________________________________" << endl
      << "                                                                  " << endl
      << "   Software for groundwater modelling with mutltiple moving       " << endl
      << "   interfaces in multilayer 3D subsurface aquifers.               " << endl
      << "                                                                  " << endl
      << "   This program is a module of HydroApps multiscale-multiphysics  " << endl 
      << "   modelling framework for advanced hydrogeological applications. " << endl
      << "                                                                  " << endl
      << "                   Algorithm and code developed by A. Sbai @ BRGM " << endl
      << "__________________________________________________________________" << endl
      << "                                                                  " << endl
      //<< "                                                                  " << endl
      //<< "                   888888b.   8888888b.   .d8888b.  888b     d888 " << endl
      //<< "                   888  *88b  888   Y88b d88P  Y88b 8888b   d8888 " << endl 
      //<< "                   888  .88P  888    888 888    888 88888b.d88888 " << endl
      //<< "                   8888888K.  888   d88P 888        888Y88888P888 " << endl
      //<< "                   888  *Y88b 8888888P*  888  88888 888 Y888P 888 " << endl
      //<< "                   888    888 888 T88b   888    888 888  Y8P  888 " << endl 
      //<< "                   888   d88P 888  T88b  Y88b  d88P 888   *   888 " << endl
      //<< "                   8888888P*  888   T88b  *Y8888P88 888       888 " << endl
      //<< "                                                                  " << endl
      //<< "                                                                  " << endl
      //<< "==================================================================" << endl
      << endl << flush;
}
