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

#ifndef HAPP_RADIAL_MESH
#define HAPP_RADIAL_MESH

// Project headers 
#include <strucmesh.h>
#include <point3d.h>

using namespace std;


namespace happ
{
   /**
    * Class for internal build of a 3D radial mesh which is oriented 
    * along theta-, r-, and z- directions by mapping them to I-, J-
    * and K- directions respectively. 
    *
    *
    * \todo: 
    *
    * 1) use the second component of 'length' variable to generalize 
    *    IJ discretization for an elliptical external boundary. 
    *
    * 2) add pointer to parent mesh object and associated nodes when 
    *    this radial mesh is going to be dynamically constructed from 
    *    a groundwater multiscale simulation involving a dynamically 
    *    hidden mesh refinement to improve numerical accuracy locally. 
    *
    * 3) add methods to automatically synchronize a groundwater variable 
    *    between this mesh and its parent during an iterative simulation 
    *    process. 
    *
    *    First tests could just try to get a local flow simulation working 
    *    correctly when the coarse level variable is used as a boundary 
    *    condition to initialize the local flow simulation. 
    */
   class RadialMesh : public StrucMesh
   {
   private:
      /// outermost radius
      double rOuter;

      /// multiplier for geometric/log spacing 
      double multiplier; 

      /// calculate initial spacing dr of the radial mesh 
      double InitialSpacing();

      /// arc fraction along which to dispose theta-nodes 
      double arcFraction;

      /// is the radial mesh closed along i- or theta- direction 
      bool isClosed;

   public:
      /// domain along r, theta, et z directions 
      Tpoint3D length; 

      /// well or innermost radius 
      double rInner; 

      /// enumerated type for radial spacing
      enum spacing_
      {
         uniform, log, geometric 
      };

      /// type of radial spacing 
      spacing_ radialSpacing; 

      /// default constructor 
      RadialMesh(const unsigned int _nt, 
         const unsigned int _nr, 
         const unsigned int _nz,
         const Tpoint3D _length,
         const double _rInner,
         const spacing_ _radialSpacing)
      {
         ni = _nt;
         nj = _nr;
         nSlices = _nz;
         length = _length;
         rInner = _rInner;
         radialSpacing = _radialSpacing;
         rOuter = length.x - rInner;
         isRadial = true; 
         multiplier = 1.2;
         arcFraction = 1.0;
         isClosed = true;
      }

      /// destructor 
      virtual ~RadialMesh() {}

      inline void SetMultiplier(const double mult)
      {
         if (mult > 1.0) multiplier = mult;
         else if (mult < 1.0) multiplier = 1.0/mult; 
         else if (mult == 1.0) radialSpacing = spacing_::uniform;
      }

      inline double GetMultiplier() const {
         return multiplier;
      }

      inline void SetArcFraction(const double frac)
      {
         if (abs(frac) >= 1.0) {
            arcFraction = 1.0;
            isClosed = true;
         }
         else {
            arcFraction = abs(frac);
            isClosed = false;
         }
      }

      inline double GetArcFraction() const {
         return arcFraction; 
      }

      /// generate the mesh 
      virtual bool Generate(); 

      /// overloading of output operator
      friend ostream& operator <<(ostream& os, RadialMesh& mesh);

   };
}

#endif // HAPP_RADIAL_MESH
