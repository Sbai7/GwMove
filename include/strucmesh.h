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

#ifndef HAPP_STRUC_MESH
#define HAPP_STRUC_MESH

// Project headers
#include <AdString.h>
#include <point3d.h>

// C++ headers
#include <iostream>
#include <list>
using namespace std;


namespace happ
{
   /**
    * Logically rectangular IJK structured mesh class. 
    * Works on a structured IJK mesh for groundwater modelling apps. 
    *
    * \todo: 
    *
    * 1) Add a boolean variable to indicate if the mesh is being used 
    *    as a child mesh (i.e. in a multiscale framework) or not. 
    * 2) Add pointer to parent mesh in the hierarchy from which this 
    *    mesh has been built. 
    * 3) For child meshes add methods for automatic generation of the 
    *    objects from user and computational requirements. 
    * 4) Add methods to export a full hierarchy of meshes in many zones 
    *    of a single Tecplot file (useful for debugging) purposes. 
    * 5) Add methods to transfer data (i.e. computed variables) from 
    *    a parent to child mesh and vice-versa. 
    * 6) ... 
    */
   class StrucMesh
   {
   protected:
      /// mesh dimension
      unsigned int dimension;

      /// number of columns
      unsigned int ni;

      /// number of rows
      unsigned int nj;

      /// number of nodes
      unsigned int nnodes;

      /// number of elements
      unsigned int nelements;

      /// number of faces;
      unsigned int nfaces;

      /// number of slices in multilayer mesh
      unsigned int nSlices;

      /// number of top moving slices in multilayer mesh
      unsigned int nMovingSlices;

      /// nodal coordinates of finite element mesh stored in full interlace mode,
      /// that is X1, Y1, Z1, ... , Xn, Yn, Zn
      double *coordinates;

      /// connectivity list of elements-nodes given in full interlace mode
      unsigned int *connectivity;

    public:
       /// The mesh is radial/circular along I-direction 
       bool isRadial;

       /// soil id's table
       unsigned int *soils;

       /// area's of nodal control patchs on top surface
       double *surface;

       /// default constructor
       StrucMesh();

       /// initialize the mesh and allocate memory
       void Init(unsigned int dim, double *coord, unsigned int *connec, unsigned int *soil_id);

       /// default destructor
       virtual ~StrucMesh();

       /// Sets mesh dimension
       inline void SetDim(const unsigned int dim) {
          dimension = dim;
       }

       /// Gets mesh dimension
       inline unsigned int GetDim() const {
          return dimension;
       }

       /// Sets number of columns
       inline void SetNi(const unsigned int nColumns) {
          ni = nColumns;
       }

       /// Gets number of columns
       inline unsigned int GetNi() const {
          return ni;
       }

       /// Sets number of rows
       inline void SetNj(const unsigned int nRows) {
          nj = nRows;
       }

       /// Gets number of rows
       inline unsigned int GetNj() const {
          return nj;
       }

       /// Sets number of nodes
       inline void SetNnodes(const unsigned int nn) {
          nnodes = nn;
       }

       /// Gets number of elements
       inline unsigned int GetNnodes() const {
          return nnodes;
       }

       /// Sets number of elements
       inline void SetNelements(const unsigned int ne) {
          nelements = ne;
       }

       /// Gets number of elements
       inline unsigned int GetNelements() const {
          return nelements;
       }

       /// Sets number of Slices
       inline void SetNslices(const unsigned int ns) {
          nSlices = ns;
       }

       /// Gets number of slices
       inline unsigned int GetNslices() const {
          return nSlices;
       }

       /// Sets number of moving slices
       inline void SetNmovSlices(const unsigned int nms) {
          nMovingSlices = nms;
       }

       /// Gets number of moving slices
       inline unsigned int GetNmovSlices() const {
          return nMovingSlices;
       }

       /// returns local node number j in element ie
       inline unsigned int GetNode(const unsigned int ie, const unsigned int j) {
          if (dimension == 2)
             return connectivity[4 * ie + j];
          else
             return connectivity[8 * ie + j];
       }

       // returns X-ccordinates array 
       inline void GetX(double* X) const {
          if (!X) return;
          unsigned int j;
          for (j = 0; j < nnodes; j++) X[j] = coordinates[dimension*j];
       }

       // returns Y-ccordinates array 
       inline void GetY(double* Y) const {
          if (!Y) return;
          unsigned int j;
          for (j = 0; j < nnodes; j++) Y[j] = coordinates[dimension*j + 1];

       }

       // returns Z-coordinates array 
       inline void GetZ(double* Z) const {
          if (!Z) return;
          unsigned int j;
          for (j = 0; j < nnodes; j++) Z[j] = coordinates[dimension*j + 2];
       }

       // Sets z-coordinate of node i
       inline void SetZ(const unsigned int i, const double z) {
          coordinates[dimension*i + 2] = z;
       }

       // returns z-coordinate of node i
       inline double GetZ(unsigned int i) const {
          return coordinates[dimension*i + 2];
       }

       // returns coordinates of node i
       inline void GetCoord(const unsigned int i, double coord[3]) const {
          unsigned int k;

          for (k = 0; k < 3; k++)
             coord[k] = coordinates[dimension*i + k];
       }

       inline void GetCoord(const unsigned int i, Tpoint3D &pt) const {
          double coord[3];
          GetCoord(i, coord);

          pt.x = coord[0];
          pt.y = coord[1];
          pt.z = coord[2];
       }

       /// returns the connectivity table 
       inline unsigned int* GetConnectivity() const {
          return connectivity;
       }

       /// returns connectivity list of element ie
       inline void GetConnectivity(const unsigned int ie, unsigned int conn[8]) const {
          if (dimension == 3) {
             for (int i = 0; i < 8; i++)
                conn[i] = connectivity[8 * ie + i];
          }
          else {
             for (int i = 0; i < 8; i++)
                conn[i] = 0;
          }
       }

       /// returns coordinates of local nodes connected to element ie
       inline void GetNeighboorNodesCoord(const unsigned int ie, double coord[3][8]) const {
          unsigned int i, j;
          for (j = 0; j < 3; j++) {
             for (i = 0; i < 8; i++) {
                unsigned n = connectivity[8 * ie + i]; // node number
                coord[j][i] = coordinates[dimension*n + j];
             }
          }
       }

       /// overloading of output operator
       friend ostream& operator <<(ostream& os, StrucMesh& mesh);

       /// Reads mesh file
       bool Read(String projectFile);

       /**
       * Returns a sorted list of soil id's in ascending order.
       * Duplicate id's are eliminated from the list.
       */
       bool GetListOfSoilId(std::list<unsigned int>& id);

       /// Computes nodal surfaces of top most mesh slice (i.e. topographic surface)
       void ComputeNodalSurfaces();


       /* New methods which need to be implemented */

       /// Detect, extract, and sort slices of material discontinuity in the mesh
       // we need only to go through one Z-column of the mesh to do the job. 
       // Assuming a multilayer aquifer system where each layer is homogeneous 
       void ExtractMaterialInterfaces(std::list<unsigned int>& slicesId);

       /// Topological mesh transformation involving face opening
       // i.e. a desaturated node in a confining layer leading to all 
       // neigbouring elements to split along the Z-direction !
       // all such nodes (nodes) will duplicate themselves to construct 
       // a new slice of the nodes (but some are identical to those in 
       // the upper slice)
       bool InsertSlice(unsigned int aboveSliceID, double *newZ);

       /// Topological mesh transformation involving saturation of a confining 
       /// layer and closeup of an existing water table slice. 
       // This is the reciprocal action of IncertSlice method
       bool RemoveSlice(unsigned int sliceId);

       /// Laplacian mesh smooting except for nodes at solid and 
       /// water table(s) interfaces. 
       bool LaplacianSmooting();

       /// Computes mesh quality with a scaled Jacobian metric
       void ComputeScaledJacobianQuality();

   };
}

#endif // HAPP_STRUC_MESH
