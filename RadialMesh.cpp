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

#define _USE_MATH_DEFINES

#include <RadialMesh.h>

#include <cstdlib>
#include <math.h>


namespace happ
{

   bool RadialMesh::Generate()
   { 
      if (ni == 0 && nj == 0 && nSlices == 0) return false;

      dimension = 3; // 3D mesh
      nnodes = ni*nj*nSlices;

      // allow for one more element along each i-direction 
      // which is theta-oriented in this particular context 
      if (isClosed) nelements = ni*(nj-1)*(nSlices-1);
      else nelements = (ni - 1)*(nj - 1)*(nSlices - 1);

      /* 1) Calculate coordinates of nodal points in the mesh */
      coordinates = new double[dimension*nnodes];
      connectivity = new unsigned int[8*nelements];

      unsigned int i, j, k;
      unsigned int index;
      //com unsigned int count = 0;

      double t = 0.0;    // angle in radians along i-direction 
      double r;          // radius along r-direction 
      double dr0 = InitialSpacing();         // spacing variable 
      double dr = dr0;
      double z = 0.0;    // z-coordinate
      double f, F;

      // for each k layer do:
      for (k = 0; k < nSlices; k++) {
         if (k > 0) z -= length.z/(nSlices-1);
         else z = 0.0;

         // for each j radius do:
         for (j = 0; j < nj; j++) {
            if (j > 0) {
               switch (radialSpacing) {
               case spacing_::uniform:
                  r += dr;
                  break;
               case spacing_::log:
                  f = (double) (nj-j);
                  F = (double) nj;
                  r = length.x - (rOuter*log10(f)/log10(F));
                  break;
               case spacing_::geometric:
                  r += dr; 
                  dr *= multiplier; 
               }
               
            }
            else {
               r  = rInner;
               dr = dr0;
            }
            
            // for each i angle do:
            for (i = 0; i < ni; i++) {
               if (i > 0) t += arcFraction * (2 * M_PI) / ni;
               else t = 0.0;
               index = i + j*ni + k*ni*nj;
               coordinates[dimension*index]     = r * cos(t);
               coordinates[dimension*index + 1] = r * sin(t);
               coordinates[dimension*index + 2] = z;
            }

         }

      }

      /* 2) Build element-nodes connectivity table of the mesh */
      unsigned int el_index = 0, iend;
      if (isClosed) iend = ni;
      else iend = ni-1; 

      for (k = 0; k < nSlices-1; k++) {
         for (j = 0; j < nj-1; j++) {
            for (i = 0; i < iend; i++) {
               index    = i + j*ni + k*ni*nj;
               // connectvity in the upper plane 
               connectivity[8 * el_index] = index;
               connectivity[8 * el_index + 1] = index + ni;
               if (i < ni-1 || !isClosed)
                  connectivity[8 * el_index + 3] = index + 1;
               else // i == ni -1 && isClosed == true
                  connectivity[8 * el_index + 3] = j*ni + k*ni*nj;
               connectivity[8 * el_index + 2] = connectivity[8 * el_index + 3] + ni;

               // then in the lower plane 
               connectivity[8 * el_index + 4] = connectivity[8 * el_index]     + ni*nj;
               connectivity[8 * el_index + 5] = connectivity[8 * el_index + 1] + ni*nj;
               connectivity[8 * el_index + 6] = connectivity[8 * el_index + 2] + ni*nj;
               connectivity[8 * el_index + 7] = connectivity[8 * el_index + 3] + ni*nj;
               el_index++; 
            }
         }
      }

      return true;
   }

   double RadialMesh::InitialSpacing()
   {
      unsigned int j; 
      double dr;
      double sum = 0.0;
   
      switch (radialSpacing) {

      case spacing_::uniform:
         dr = length.x / (nj - 1);
         break;

      case spacing_::log:
         multiplier = 10.0;
         for (j = 0; j < nj - 1; j++)
            sum += std::pow(multiplier, j);
         dr = length.x / sum;
         break;

      case spacing_::geometric:
         for (j = 0; j < nj - 1; j++)
            sum += std::pow(multiplier, j);
         dr = length.x / sum;
         break;

      }

      return dr;
   }

   ostream& operator <<(ostream& os, RadialMesh& mesh)
   {
      StrucMesh m = (StrucMesh) mesh;
      os << m; 
      return os;
   }
}
