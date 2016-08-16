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

#ifndef _PARTICLE_TRACKER_H
#define _PARTICLE_TRACKER_H

// Project headers
#include <hex8.h>
#include <strucmesh.h>
#include <flowline.h>
#include <soil.h>

// C++ headers 
#include <vector>
#include <list> 

// Eigen headers
//#include <Eigen/Core>
//#include <Eigen/Dense>

using namespace std;
using namespace happ; 
//using namespace Eigen;



/**
* Particle tracking class for tracing & visualization of three-dimensional 
* ground-water flow / velocity fields in subsurface aquifers.
*/
class ParticleTracker
{
private:

   /// Flowlines to be computed by particle tracking algorithm 
   std::vector<FlowLine> flowline; 

   /// initial positions of particles to be released 
   std::vector<Tpoint3D> initialCoords;

   /// is particle tracking in Forward or backwards direction?
   bool isForward;

   /// maximal time of particle tracking
   double maxTime;

   /// pointer to mesh used for computational tasks
   StrucMesh *mesh;

   /// pointer to computed nodal groundwater heads
   double *head;

   /// pointer to computed nodal groundwater flow rates 
   double *flow;

   /// current time 
   double time; 

   /// local element-wise (CFL ?) time steps
   double *dt;

   /// temporary local nodal-wise interpolation coefficients
   double weight[8][3]; 

   /// minimal element-wise coordinates 
   std::vector<double> xemin, yemin, zemin;

   /// maximal element-wise coordinates 
   std::vector<double> xemax, yemax, zemax;

   /// min and max coordinates of the computational domain 
   double xmin[3], xmax[3]; 
 
   /// soil types of the groundwater model 
   std::vector<Soil> Soils;


   /// Calculate element-wise coordinate bounds 
   void ComputeElemCoordBounds();

   /// Compute elemental time steps 
   void ComputeTimeSteps();

   /// Naive implementation of hex element search containing a given point 
   bool NaiveFindElement(const Tpoint3D pt, unsigned int &elem, Tpoint3D & pt_local);

   /// Finds in which hexahedral element a point is situated 
   bool FindElement(const Tpoint3D pt, unsigned int &elem, Tpoint3D &pt_local);

   /// Finds in which quadrilateral 2D element a 2D point in situated
   bool FindQuadElement(const Tpoint3D pt, unsigned int &quadElem);

   /** 
    * Finds in which layer (hex element) a 3D point is situated knowning in which 
    * quad element its 2D projection is located. 
    * Here, we assume that the mesh has a multilayer structure with identical nodal 
    * X,Y coordinates in all mesh slices. 
    */ 
   bool FindLayer(const Tpoint3D pt, const unsigned int quadElem, unsigned int &elem, Tpoint3D &pt_local);

   /// Checks wether a point is inside an element?
   bool IsInside(const Tpoint3D pt, const unsigned int elem, Tpoint3D &pt_local);

   /// detects a source or sink element during particle tracking 
   bool IsSourceOrSink(const unsigned int elem, Tpoint3D &pt);

   

protected:

public:
   /// Initialize a particle tracking problem 
   ParticleTracker(StrucMesh *nw_mesh,
      std::vector<Soil> soils,
      std::vector<Tpoint3D> coords, 
      double *Head, 
      double *Flow);

   /// default Fem destructor
   virtual ~ParticleTracker();


   /// Execution of the particle tracking solver  
   virtual bool Solve(); 

   // object streaming ... 
};

#endif // #ifndef _PARTICLE_TRACKER_H
