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
#include <velocity.h>


NodalVelocity::NodalVelocity()
{
   Xvelocity = NULL;
   Yvelocity = NULL;
   Zvelocity = NULL;
   volume    = NULL;
}

NodalVelocity::NodalVelocity(StrucMesh *pMesh, vector<Soil> soils, double *pVar)
{
   mesh = pMesh;
   var  = pVar;
   dimension = mesh->GetNnodes();
   Soils = soils;

   // memory allocation
   Xvelocity = new double [dimension];
   Yvelocity = new double [dimension];
   Zvelocity = new double [dimension];
   volume    = new double [dimension];

   // initialization
   unsigned int i;
   for (i = 0; i < dimension; i++) {
      Xvelocity[i] = Yvelocity[i] = Zvelocity[i] = 0.;
      volume[i] = 0.;
   }
}

NodalVelocity::~NodalVelocity()
{
   if (Xvelocity) delete[] Xvelocity;
   if (Yvelocity) delete[] Yvelocity;
   if (Zvelocity) delete[] Zvelocity;
   if (volume)    delete[] volume;
}

bool NodalVelocity::Compute()
{
   unsigned int ie, l, k, i, j ;
   double grad[3];
   
   // loop over mesh elements
   for (ie = 0; ie < mesh->GetNelements(); ie++) {
   
      for (l = 0; l < 8; l++) {
      
         // calculate local conductance matrix
         unsigned int connec[8];
         double coord[3][8];
         double K[3]; // hydraulic conductivities of element ie
         mesh->GetConnectivity(ie,connec);
         mesh->GetNeighboorNodesCoord(ie,coord);
         Soil soil = Soils[mesh->soils[ie]-1];
         soil.GetK(K);
         Hex8 hex(K,coord,connec);
         hex.dimension = 3; // Consider Hexahedral elements
         hex.SetGaussNumber(l);
         hex.ComputeBasis();
         hex.ComputeGradients();
         double det = hex.det;

         // add nodal contribution 
         for (i = 0; i < 3; i++) {
            grad[i] = 0.;
            for (j = 0; j < 8; j++)
               grad[i] += var[connec[j]]*hex.bx[i][j];
         }

         // in case of unsaturated zone handling weighted pressure potential 
         // must be calculated in this place --- THIS IS TO ADD LATER ------

         // calculate right hand side of equations
         k = connec[l];
         double sum = 0.;
         Xvelocity[k] += -grad[0]*K[0]*det; // if ZNS multiply it by kr/wc factor
         Yvelocity[k] += -grad[1]*K[1]*det;
         Zvelocity[k] += -grad[2]*K[2]*det;

         // add to total control volume of node k
         volume[k]    += det;
      }

   }

   for (i = 0; i < mesh->GetNelements(); i++) {
      Xvelocity[i] /= volume[i];
      Yvelocity[i] /= volume[i];
      Zvelocity[i] /= volume[i];
   }

   return true;
}
