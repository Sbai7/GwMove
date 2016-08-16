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
#include <FileReader.h>


Project::Project()
{
   // set default PCG parameters
   pcgMaxIter   = 500;
   pcgTolerance = 1.E-12;
   pcgSilent    = true;

   // set default water table parameters
   wtMaxIter   = 50;
   wtTolerance = 0.1;
   nMoving     = 1;

   // set default timing parameters
   initialTime = 0.;
   finalTime   = 100.;
   timeStep    = 10.;

   // initialize the list of soil types
    soils.clear();

    // Groundwater flow process by default
    flow = SATURATED_FLOW; 

    // set types of user-specified type analysis 
    mba.analysisType   = _Zone_Analysis_Params_::analysis_::mass_balance;
    stata.analysisType = _Zone_Analysis_Params_::analysis_::statistics;
}

Project::~Project()
{
   // Free pointers to zones in members of 'mba' and 'stata' 
}

void Project::SetFullName(String full_name)
{
   fullName = full_name;
}

String Project::GetFullName()
{
   return fullName;
}

_Pcg_Params_ Project::GetPcgParams()
{
   _Pcg_Params_ pcg;
   pcg.iterations = pcgMaxIter;
   pcg.tolerance  = pcgTolerance;
   pcg.silent     = pcgSilent;

   return pcg;
}

_Wt_Params_ Project::GetWaterTableParams()
{
   _Wt_Params_ wt;
   wt.iterations = wtMaxIter;
   wt.tolerance  = wtTolerance;
   wt.nMoving    = nMoving;

   return wt;
}

_Time_Params_ Project::GetTimeParams()
{
   _Time_Params_ tp;
   tp.t0 = initialTime;
   tp.tf = finalTime;
   tp.dt = timeStep;

   return tp;
}

ostream& operator <<(ostream& os, Project& prj)
{
   os << "TITLE = " << prj.title << endl;
   os << "FLOW =  " << prj.flow << endl;
   os << "PCG   = " << setw(8) << prj.pcgMaxIter
      << setw(16) << prj.pcgTolerance
      << setw(8) << prj.pcgSilent << endl;
   os << "WT    = " << setw(8) << prj.wtMaxIter
      << setw(16) << prj.wtTolerance
      << setw(8) << prj.nMoving << endl;
   int nSoils = (int)prj.soils.size();
   os << "SOILS = " << nSoils << endl;
   for (int i = 0; i < nSoils; i++)
      os << setw(16) << prj.soils[i].name << setw(8) << prj.soils[i].GetId()
      << setw(16) << prj.soils[i].GetKx() << setw(16) << prj.soils[i].GetKy()
      << setw(16) << prj.soils[i].GetKz() << setw(16) << prj.soils[i].GetS0()
      << endl;
   os << "TIME  = " << setw(16) << prj.initialTime <<
      setw(16) << prj.finalTime <<
      setw(16) << prj.timeStep << endl;
   os << endl;
   return os;
}

bool Project::Read()
{
   unsigned int i, j, nSoils;
   unsigned int i1, j1, k1, i2, j2, k2;
   String label;
   ijkZone *ijkZn; 

   filebuf file;
   String dataFile = fullName + ".dat";
   
   file.open(dataFile.GetBuffer(), ios::in);
   istream is(&file);
   FileReader fr(is);

   // main loop for parsing the input project file 
   for (i = 0; i < NUM_KEYS; i++) { 

      int nTokens = fr.GetLine();
      if (nTokens < 3) {
         // print error message including the line where is occured
         return false;
      }
      
      // (1) project TITLE 
      if (String(fr.GetToken(0)) == "TITLE") {
         title = "";
         for (j = 2; j < nTokens; j++) title += String(fr.GetToken(j));
      }
      
      // (2) FLOW processes 
      else if (String(fr.GetToken(0)) == "FLOW_PROCESS" || 
               String(fr.GetToken(0)) == "FP") {
         if      (String(fr.GetToken(2)) == "saturated")   flow = SATURATED_FLOW;
         else if (String(fr.GetToken(2)) == "unconfined")  flow = UNCONFINED_FLOW;
         else if (String(fr.GetToken(2)) == "unsaturated") flow = UNSATURATED_FLOW;
         else flow = SATURATED_FLOW; 
            // add a handler for an error occurence (later) in this context
      }
      
      // (3) PCG: preconditionned conjugate gradient section 
      else if (String(fr.GetToken(0)) == "PCG") {
         pcgMaxIter   = atoi(fr.GetToken(2));
         pcgTolerance = atof(fr.GetToken(3));
         pcgSilent    = (bool) atoi(fr.GetToken(4));
      }
      
      // (4) WT: Water Table itrative parameters section 
      else if (String(fr.GetToken(0)) == "WTI" || 
               String(fr.GetToken(0)) == "WATER_TABLE") {
         wtMaxIter   = atoi(fr.GetToken(2));
         wtTolerance = atof(fr.GetToken(3));
         nMoving     = atoi(fr.GetToken(4));
      }
      
      // (5) SOILS: soil types and parameters section 
      else if (String(fr.GetToken(0)) == "SOILS" || 
               String(fr.GetToken(0)) == "MATERIALS") {
         nSoils =  atoi(fr.GetToken(2));
         for (j = 0; j < nSoils; j++) {
            nTokens = fr.GetLine();
            if (nTokens < 6) {
               // error message
               return false;
            }
            Soil soil;
            soil.name = String(fr.GetToken(0));
            soil.SetId(atoi(fr.GetToken(1)));
            soil.SetKx(atof(fr.GetToken(2)));
            soil.SetKy(atof(fr.GetToken(3)));
            soil.SetKz(atof(fr.GetToken(4)));
            soil.SetS0(atof(fr.GetToken(5)));
            soils.push_back(soil);
         }
      }

      // (6) MASS_BALANCE for choice of mass-balance analysis types 
      else if (String(fr.GetToken(0)) == "MASS_BALANCE" ||
               String(fr.GetToken(0)) == "MB") {
         mba.number = atoi(fr.GetToken(2)); 
         if (mba.number > 0) {
            for (j = 0; j < mba.number; j++) {
               nTokens = fr.GetLine(); 
               if (nTokens < 2) {
                  // error message 
                  return false; 
               }
               else if (nTokens == 2) { // we expect a BC type analysis  
                  label = String(fr.GetToken(0)); 
                  int doIt = atoi(fr.GetToken(1)); 
                  if (label == "BC" && doIt > 0) mba.bcMBA = true; 
               }
               else if (nTokens == 7) { // we expect a zone type analysis 
                  label = String(fr.GetToken(0)); 
                  mba.label.push_back(label); 
                  i1 = atoi(fr.GetToken(1));
                  j1 = atoi(fr.GetToken(2));
                  k1 = atoi(fr.GetToken(3));
                  i2 = atoi(fr.GetToken(4));
                  j2 = atoi(fr.GetToken(5));
                  k2 = atoi(fr.GetToken(6));
                  ijkZn = new ijkZone(i1, j1, k1, i2, j2, k2);
                  mba.Add(ijkZn);
               }
            }
         }
      }

      // (6) STATISTICS for choice of statistical analysis types 
      else if (String(fr.GetToken(0)) == "STATISTICS" ||
         String(fr.GetToken(0)) == "STAT") {
         stata.number = atoi(fr.GetToken(2));
         if (stata.number > 0) {
            for (j = 0; j < stata.number; j++) {
               nTokens = fr.GetLine();
               if (nTokens < 2) {
                  // error message 
                  return false;
               }
               else if (nTokens == 7) { // we expect a zone type analysis 
                  label = String(fr.GetToken(0));
                  stata.label.push_back(label);
                  i1 = atoi(fr.GetToken(1));
                  j1 = atoi(fr.GetToken(2));
                  k1 = atoi(fr.GetToken(3));
                  i2 = atoi(fr.GetToken(4));
                  j2 = atoi(fr.GetToken(5));
                  k2 = atoi(fr.GetToken(6));
                  ijkZn = new ijkZone(i1, j1, k1, i2, j2, k2);
                  stata.Add(ijkZn);
               }
            }
         }
      }

      // (8) TIME section for stress periods and time stepping schemes
      else if (String(fr.GetToken(0)) == "TIME") {
         initialTime = atof(fr.GetToken(2));
         finalTime   = atof(fr.GetToken(3));
         timeStep    = atof(fr.GetToken(4));
      }
   }

    return true;
}

bool Project::CheckSoilIds(std::list<unsigned int>& id)
{
    std::list <unsigned int>::iterator iter;

    for (iter = id.begin(); iter != id.end(); iter++) {
        bool found = false;
        for (unsigned int i = 0; i < soils.size(); i++) {
         if ( soils[i].GetId() == *iter ) {
            found = true;
            break;
         }
      }
      if ( !found ) return false;
   }

   return true;
}

void Project::SetFlow(unsigned int flo)
{
    flow = flo;
}

unsigned int Project::GetFlow()
{
   return flow;
}

