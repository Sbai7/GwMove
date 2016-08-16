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

#ifndef _VELOCITY_H
#define _VELOCITY_H

// Project headers
#include <strucmesh.h>
#include <hex8.h>
#include <soil.h>

// C++ headers
#include <vector>
using namespace std;
using namespace happ; 


/**
 * Defines post-processed velocity components at nodal points
 * of the finite element mesh
 */
class NodalVelocity
{
private:
   /// pointer to FE mesh
   StrucMesh *mesh;

   /// pointer to model scalar variable
   double *var;

   /// size of the velocity vectors 
   unsigned int dimension;

   /// Nodal control volume
   double *volume;

   /// reference to vector of soil types defining the FE problem
   vector<Soil> Soils;

public:
   /// Nodal X-velocity component
   double *Xvelocity;

   /// Nodal Y-velocity component
   double *Yvelocity;

   /// Nodal Z-velocity component
   double *Zvelocity;

   /// Default constructor
   NodalVelocity();

   /// Constructor from FE mesh and groundwater heads
   NodalVelocity(StrucMesh *pMesh, vector<Soil> soils, double *pVar);

   /// Destructor
   virtual ~NodalVelocity();

   /// Computes vector velocity components
   bool Compute();

};

#endif // #ifndef _VELOCITY_H
