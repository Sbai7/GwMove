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

#ifndef _PATH_LINE_H
#define _PATH_LINE_H

// Project headers
#include <point3d.h>

// C++ headers
#include <vector>
#include <list> 
using namespace std;
using namespace happ; 


/**
* Flow line class which characterize a three-dimensional particle
* path traced from a three-dimensional velocity field.
*/
class FlowLine
{
private:
   /// a specific label or ID
   unsigned int id;

   /// traced forwards or backwards?
   bool isForward;

protected:

public:

   /// Initial coordinates of the moving particle 
   Tpoint3D initPoint; 

   /// positions vector of points constructing the flowlines
   std::vector<Tpoint3D> point;

   /// particle travel times at 'points' positions 
   std::vector<double> time;

   /// materials vector traced by this flowline 
   std::vector<unsigned int> soilId;

   /// model layer where this function is situated 


   /// Default constructor 
   FlowLine(unsigned int ID, bool Forward, Tpoint3D & init) {
      id = ID;
      isForward = Forward;
      initPoint = init;
   };

   /// Destructor
   virtual ~FlowLine() {
      point.clear();
      time.clear();
      soilId.clear();
   };

   /// Gets flowline ID 
   unsigned int GetId() {
      return id;
   };

   /// Is tracing direction forwards?
   bool IsForwardDirection() {
      return isForward;
   };

   /// overloading of input operator
   friend istream & operator >>(istream & is, FlowLine & pathline)
   {
      char buff[120];
      double x, y, z;

      // line containing the initial position of the pathline 
      // and the particle tracking direction !
      is >> x >> y >> z >> pathline.isForward;
      is.getline(buff, 120, '\n');

      // update initial point variable 
      pathline.initPoint = Tpoint3D(x, y, z);

      return is;
   };

   /// overloading of output operator 
   friend ostream & operator <<(ostream & os, FlowLine & pathline)
   {
      os << std::setiosflags(ios::scientific) << std::setprecision(8);

      // stream all registered particle positions for this pathline  
      unsigned int i; 
      unsigned int np = pathline.point.size();
      os << np << std::endl; 
      for (i = 0; i < np; i++) {
         os << std::setw(15) << pathline.point[i].x 
            << std::setw(15) << pathline.point[i].y 
            << std::setw(15) << pathline.point[i].z 
            //<< std::setw(15) << time[i] 
            //<< std::setw(5)  << soilId[i]
            << std::endl;
      }

      return os;
   };

};

#endif // #ifndef _PATH_LINE_H
