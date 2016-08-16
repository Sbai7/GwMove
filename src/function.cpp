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
#include <function.h>

// C++ headers
#include <iomanip>


namespace happ
{
   TimeFunction::TimeFunction()
   {
      name.Empty();
      points.clear();
   }

   TimeFunction::~TimeFunction()
   {
      points.clear();
   }

   TimeFunction::TimeFunction(const TimeFunction &fct)
   {
      name = fct.name;
      for (unsigned int i = 0; i < fct.points.size(); i++)
         points.push_back(fct.points[i]);
   }

   TimeFunction& TimeFunction::operator=(const TimeFunction &fct)
   {
      name = fct.name;
      points.resize(fct.points.size());
      for (unsigned int i = 0; i < fct.points.size(); i++)
         points[i] = fct.points[i];

      return *this;
   }

   double TimeFunction::MinX()
   {
      double min = 1.0E+99;

      for (unsigned int i = 0; i < points.size(); i++)
         min = (points[i].x < min ? points[i].x : min);

      return min;
   }

   double TimeFunction::MaxX()
   {
      double max = -1.0E+99;

      for (unsigned int i = 0; i < points.size(); i++)
         max = (points[i].x > max ? points[i].x : max);

      return max;
   }

   double TimeFunction::MinY()
   {
      double min = 1.0E+99;

      for (unsigned int i = 0; i < points.size(); i++)
         min = (points[i].y < min ? points[i].y : min);

      return min;
   }

   double TimeFunction::MaxY()
   {
      double max = -1.0E+99;

      for (unsigned int i = 0; i < points.size(); i++)
         max = (points[i].y > max ? points[i].y : max);

      return max;
   }

   point TimeFunction::Max()
   {
      return point(MaxX(),MaxY());
   }

   void TimeFunction::Deriv(vector<double> & v)
   {
      unsigned int i; 
      v.clear();

      for (i = 1; i < points.size(); i++) {
         double dX = points[i].x - points[i - 1].x;
         double dY = points[i].y - points[i - 1].y;
         v.push_back(dY / dX);
      }
   }

   ostream& operator <<(ostream& os, TimeFunction& fct)
   {
      unsigned int i;

      os << setiosflags(ios::fixed) << setprecision(8) << setiosflags(ios::showpoint);
      os << "  -name          " << fct.name << endl;
      os << "  -points        " << fct.Size() << endl;
      for (i = 0; i < fct.Size(); i++) {
         os << setw(27) << fct.points[i].x << setw(20) << fct.points[i].y << endl;
      }

      return os;
   }

   istream& operator >>(istream& is, TimeFunction& fct)
   {
      unsigned int i;
      char ch[1024];
      point pt;

      for (i = 0; i < 3; i++) { // 3 is the number of allowed keywords ...
         is.getline(ch, 120, '\n');
         String strBuffer = String(ch);
         strBuffer.Trim();
         int pos = strBuffer.Find(' ');
         if (pos > 0) {
            String key /*= strBuffer.Left(pos)*/;
            String txt /*= strBuffer.Right(strBuffer.GetLength() - pos)*/;
            txt.Trim();
            if (key == "-name")
               fct.name = txt;
            else if (key == "-index")
               continue;
            else if (key == "-points")
               fct.Size(atoi(txt.GetBuffer()));
         }
      }

      for (i = 0; i < fct.Size(); i++) {
         is >> pt.x >> pt.y;
         fct.Point(i, pt);
         is.getline(ch, 120, '\n');
      }

      return is;
   }

}
