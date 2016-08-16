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
#include <strucmesh.h>
#include <FileReader.h>

// C++ headers
#include <iomanip>


namespace happ
{
   StrucMesh::StrucMesh()
   {
      ni = 0; nj = 0; nSlices = 0;
      nnodes = 0; nelements = 0; 
      coordinates = NULL;
      connectivity = NULL;
      soils        = NULL;
      surface      = NULL;
      isRadial     = false;
   }

   void StrucMesh::Init(unsigned int dim, double *coord, unsigned int *connec, unsigned int *soil_id)
   {
      // sets mesh dimension
      dimension = dim;

      // initialize class pointers
      coordinates  = coord;
      connectivity = connec;
      soils        = soil_id;
      surface      = NULL;
   }

   StrucMesh::~StrucMesh()
   {
      if (coordinates)  delete[] coordinates;
      if (connectivity) delete[] connectivity;
      if (soils)        delete[] soils;
      if (surface)      delete[] surface;
   }


   bool StrucMesh::Read(String projectFile)
   {
      filebuf file;
      String meshFile = projectFile + ".fem";

      file.open(meshFile.GetBuffer(), ios::in);
      istream is(&file);
      FileReader fr(is);

      /* 1) Read file header */
      
      // number of nodes along each direction
      unsigned int nTokens = fr.GetLine();
      ni = atoi(fr.GetToken(0));
      
      nTokens = fr.GetLine();
      nj = atoi(fr.GetToken(0));

      nTokens = fr.GetLine();
      nSlices = atoi(fr.GetToken(0));

      // total number of nodes 
      nTokens = fr.GetLine();
      nnodes = atoi(fr.GetToken(0));

      // total number of elements 
      nTokens = fr.GetLine(); 
      nelements = atoi(fr.GetToken(0));

      /* 2) Read nodal coordinates section */

      // calculate number of elements from connectivity list
      unsigned int size;
      if (dimension == 2) size = 4;
      else if (dimension == 3) size = 8;

      coordinates = new double[dimension*nnodes];
      connectivity = new unsigned int[size*nelements];
      soils = new unsigned int[nelements];

      unsigned int i, j;
      for (j = 0; j < nnodes; j++) {

         nTokens = fr.GetLine();

         if (nTokens < dimension) {
            // print error message
            return false;
         }

         if (dimension == 2) {
            coordinates[dimension*j]     = atof(fr.GetToken(0));   // X-coordinate
            coordinates[dimension*j + 1] = atof(fr.GetToken(1));   // Y-coordinate
         }

         else if (dimension == 3) {
            coordinates[dimension*j]     = atof(fr.GetToken(0));  // X-coordinate
            coordinates[dimension*j + 1] = atof(fr.GetToken(1));  // Y-coordinate
            coordinates[dimension*j + 2] = atof(fr.GetToken(2));  // Z-coordinate
         }

      }

      /* 3) Read element-nodes connectivity list */
      for (j = 0; j < nelements; j++) {

         nTokens = fr.GetLine();

         if (nTokens < size + 1) {
            // print error message
            return false;
         }

         for (i = 0; i < size; i++)
            connectivity[size*j + i] = atoi(fr.GetToken(i)) - 1;

         soils[j] = atoi(fr.GetToken(size));

      }

      return true;
   }

   ostream& operator <<(ostream& os, StrucMesh& mesh)
   {
      if (!mesh.coordinates || !mesh.connectivity || !mesh.soils)
         return os;

      unsigned int i, j;

      os << setiosflags(ios::scientific);

      unsigned int size;
      if (mesh.dimension == 2) size = 4;
      else if (mesh.dimension == 3) size = 8;

      // Write file header
      os << mesh.ni        << "   # Number of nodes in x-direction" << endl; 
      os << mesh.nj        << "   # Number of nodes in y-direction" << endl;
      os << mesh.nSlices   << "   # Number of nodes in z-direction" << endl;
      os << mesh.nnodes    << "   # Number of nodes" << endl; 
      os << mesh.nelements << "   # Number of cells" << endl;

      // Write nodal-coordinates section
      for (i = 0; i < mesh.nnodes; i++) {

         if (mesh.dimension == 2)
            os << mesh.coordinates[mesh.dimension*i] <<        // X-coordinate
            mesh.coordinates[mesh.dimension*i + 1] << endl;    // Y-coordinate

         if (mesh.dimension == 3)
            os << setw(16) << mesh.coordinates[mesh.dimension*i] <<        // X-coordinate
            setw(16) << mesh.coordinates[mesh.dimension*i + 1] <<          // Y-coordinate
            setw(16) << mesh.coordinates[mesh.dimension*i + 2] << endl;    // Z-coordinate

      }

      // Write finite element connectivity section
      for (j = 0; j < mesh.nelements; j++) {

         for (i = 0; i < size; i++)
            os << setw(8) << mesh.connectivity[size*j + i] + 1;
         os << setw(5) << mesh.soils[j] << endl;

      }

      return os;
   }

   bool StrucMesh::GetListOfSoilId(std::list<unsigned int>& id)
   {
      if (!soils) return false;

      // empty the list
      id.clear();

      // fills the list with mesh elements soil id's
      for (unsigned int i = 0; i < nelements; i++) {
         id.push_back(soils[i]);
      }

      // sorts the list in ascending order
      id.sort();

      // removes adjacent duplicate elements
      id.unique();

      return true;
   }

   void StrucMesh::ComputeNodalSurfaces()
   {
      unsigned int locNodes = ni*nj;
      unsigned int locElems = (ni - 1)*(nj - 1);
      unsigned int i, j, ie, n1, n, gauss;

      if (!surface) surface = new double[ni*nj];
      for (i = 0; i < locNodes; surface[i++] = 0.);

      // transformer les coord 3D des noeuds 
      // aux coord 2D dans le repére (R,S) du plan;
      // et renumerotation dans le plan

      // define constants and local element coordinates
      const double pgauss = 1. / sqrt(3.);
      const int nxsi[4] = { -1, 1, 1, -1 };
      const int neta[4] = { -1, -1, 1, 1 };
      double bxsi[4], beta[4];
      double Jacobian[4], determinant;

      for (ie = 0; ie < locElems; ie++) {	// loop over quad elements

         Tpoint3D a[4]; // 4 quad nodes
         div_t result = div((int)ie, (int)nj - 1);
         n1 = (result.quot)*nj + result.rem;
         GetCoord(n1, a[0]);
         GetCoord(n1 + 1, a[1]);
         GetCoord(n1 + nj, a[2]);
         GetCoord(n1 + nj + 1, a[3]);

         // calcul du determinant du Jacobian
         for (gauss = 0; gauss < 4; gauss++) { // boucle sur les points de Gauss
            for (j = 0; j < 4; j++) {
               bxsi[j] = nxsi[j] * (1. + pgauss*neta[j] * neta[gauss]) / 4.;
               beta[j] = neta[j] * (1. + pgauss*nxsi[j] * nxsi[gauss]) / 4.;
               Jacobian[j] = 0.; // initialisation
            }
            for (j = 0; j < 4; j++) {
               Jacobian[0] += beta[j] * a[j].x;
               Jacobian[1] += bxsi[j] * a[j].x;
               Jacobian[2] += beta[j] * a[j].y;
               Jacobian[3] += bxsi[j] * a[j].y;
            }
            determinant = fabs(Jacobian[0] * Jacobian[3] - Jacobian[1] * Jacobian[2]);

            // add element contribution to nodal surface area
            for (j = 0; j < 4; j++) {
               switch (j) {
               case 0:
                  n = n1;
                  break;
               case 1:
                  n = n1 + 1;
                  break;
               case 2:
                  n = n1 + nj;
                  break;
               case 3:
                  n = n1 + nj + 1;
                  break;
               }
               surface[n] += determinant*(1. + pgauss*nxsi[j] * nxsi[gauss])
                  *(1. + pgauss*neta[j] * neta[gauss]) / 4.;
            }
         }
      }

   }

}
