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
#include <lamsswi.h>
#include <vector.h>

// C++ headers
#include <cmath>

using namespace happ;


LamSteadySWI::LamSteadySWI(Project *_project,
                           StrucMesh *_mesh,   	       // Mesh pointer
                           BndCondition *_bc,          // Boundary conditions
                           _Pcg_Params_ &pcgParams)    // Solver parameters
   : LamSteadyFlow(_project, _mesh, _bc, pcgParams)
{

}

void LamSteadySWI::SetInterfacePosition()
{
   unsigned int l, i, j, jj;
   unsigned int nij = mesh->GetNi() * mesh->GetNj();
   unsigned int nk = mesh->GetNslices();
   unsigned int nw = mesh->GetNmovSlices();
   unsigned int ns; // = mesh->GetNSwiMovSlices(); 
   double newZ; 

   for (l = 1; l < nw; l++) {
      for (i = 0; i < nij; i++) {
         j  = i + (l - 1) *nij;
         jj = i + (nw - 1)*nij;
         if (bc->bcType[i] != SINGLE_NODES) {
            newZ = (var[i] * (nw - l) + mesh->GetZ(jj)*(l - 1)) / (nw - 1);
            mesh->SetZ(j, newZ);
         }
      }
   }

   for (l = nk; l < nk - ns + 1; l--) {
      for (i = 0; i < nij; i++) {
         j  = i + (l - 1) *nij; 
         jj = i + (nk - ns)*nij;
         if (bc->bcType[i] != SINGLE_NODES) {
            newZ = (var[i] * (ns - nk - 1 + l) + mesh->GetZ(jj)*(nk - l)) / (ns - 1);
            mesh->SetZ(j, newZ);
         }
      }
   }
}

