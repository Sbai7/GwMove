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

#ifndef _PROJECT_H
#define _PROJECT_H

// Project headers
#include <AdString.h>
#include <soil.h>
#include <bndCondition.h>

// C++ headers
#include <iostream>
#include <vector>
#include <list>
using namespace std;
using namespace happ; 

/// Types of groundwater flow types/solvers supported by the computer program
#define SATURATED_FLOW   		1
#define UNCONFINED_FLOW			2
#define UNSATURATED_FLOW  		3  // not implemented yet
#define DENSE_FLOW   			4  // not implemented yet
#define FORCHEIMER_FLOW   		5  // not implemented yet
#define STOCKES_FLOW   			6 	// not implemented yet
#define BOUSSINESQ_FLOW       7  // not implemented yet
#define COUPLED_FLOW          8  // not implemented yet

/// Types of solute transport variants supported by the computer program
#define SATURATED_TRANSPORT   1  // not implemented yet
#define UNCONFINED_TRANSPORT  2  // not implemented yet
#define UNSATURATED_TRANSPORT 3  // not implemented yet
#define DENSE_TRANSPORT       4  // not implemented yet
#define PARTICLES_TRANSPORT   5  // not implemented yet 

/// Number of keywords in the project file 
#define NUM_KEYS              8  // to change when a new keyword is added 

struct _Pcg_Params_ {
   unsigned int iterations;
   double       tolerance;
   bool         silent;
};

struct _Wt_Params_ {
   unsigned int iterations;
   double       tolerance;
   unsigned int nMoving;
};

struct _Time_Params_ {
   double t0, tf, dt;
};

/**
 * Basic structure for IJK & BC's zone analysis
 * such as Mass-Balance or Statistics, etc.
 */
struct _Zone_Analysis_Params_ {
   /// Possible analysis types 
   enum analysis_
   {
      mass_balance, statistics 
   };

   /// type of analysis to do 
   analysis_ analysisType;

   /// total number of mass balance analysis 
   unsigned int number; 

   /// do we perform MBA by BC object ?
   bool bcMBA = false;

   /// labels of Mass balance analysis by user-given zones  
   std::vector<String> label;

   /// user-given zones of mass balance analysis 
   std::vector<ijkZone*> userZone;

   /// default constructor 
   _Zone_Analysis_Params_() {
      number = 0;
   }

   /// destructor 
   virtual ~_Zone_Analysis_Params_() {
      unsigned int i, n = userZone.size();
      for (i = 0; i < n; i++) delete userZone[i]; 
      userZone.clear(); 
   }

   /// Insert a user zone 
   inline void Insert(const unsigned int i, ijkZone *zone) {
      delete userZone[i];
      userZone.insert(userZone.begin() + i, zone);
   }

   /// Add a user zone 
   inline void Add(ijkZone *zone) {
      userZone.push_back(zone);
   }

   /// Remove a user zone 
   inline void Remove(const unsigned int i) {
      delete userZone[i]; 
      userZone.erase(userZone.begin() + i);
   }
};


/**
 * Class for management of project input and output files.
 *
 * Todo:
 * - Add output control options support
 * - Add stress periods support 
 * - Separate entries for hydrodynamic and transport parameters of soils 
 * - Add support for redefinition of soil zones from project file 
 * - Add support for observation points where measurements are given 
 * - Add support for many scearii's ...
 * 
 * - Error control when parsing the input file
 * - Design a top class for input file management from which this project
 *   class is to be derived
 * - input/output units management 
 * - for transient simulations: time, date, management 
 * 
 */
class Project
{
private:
   /// full path to the project name
   String fullName;

   /// project title
   String title;

   /// short name of the prject (excluding the path)
   String name;

   /// max. iterations of the PCG linear solver
   unsigned int pcgMaxIter;

   /// tolerance of the PCG linear solver
   double pcgTolerance;

   /// flag for silent or verbose PCG iterations
   bool pcgSilent;

   /// max. number of water table iterations
   unsigned int wtMaxIter;

   /// tolerance of water table iteration process
   double wtTolerance;

   /// number of moving mesh slices
   unsigned int nMoving;

   /// initial simulation time
   double initialTime;

   /// final simulation time
   double finalTime;

   /// simulation time step (constant through)
   double timeStep;

   /// Flow type (i.e. Saturated, unconfined, unsaturated, etc.)
   unsigned int flow;

public:
   /// list of soil types
   vector<Soil> soils;

   /// Mass balance analysis options 
   _Zone_Analysis_Params_ mba;

   /// Statistical analysis options 
   _Zone_Analysis_Params_ stata;

   /// constructor
   Project();

   /// destructor
   ~Project();

   /// Sets full name of the project
   void SetFullName(String full_name);

   /// returns full name of the project
   String GetFullName();

   /// returns PCG parameters 
   _Pcg_Params_ GetPcgParams();

   /// returns water table related parameters
   _Wt_Params_ GetWaterTableParams();

   /// returns time stepping related parameters
   _Time_Params_ GetTimeParams();

   /// overloading of output operator
   friend ostream& operator <<(ostream& os, Project& prj);

   /// Reads project file
   bool Read();

   /// checks if all soil id's in the given list exist
   bool CheckSoilIds(std::list<unsigned int>& id);

   /// Set type of flow or solver to be used
   void SetFlow(unsigned int flo);

   /// Return the type of flow and specialized solver being used
   unsigned int GetFlow();

};

#endif // #ifndef _PROJECT_H
