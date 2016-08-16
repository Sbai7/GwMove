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

#ifndef HAPP_TIME_FUNCTION
#define HAPP_TIME_FUNCTION

// Project headers
#include <AdString.h>

// C++ headers
#include <iostream>
#include <vector>
using namespace std;


namespace happ
{
   /**
    * Data structure holding a function couple t and f(t)
    */
   struct point
   {
      double x;
      double y;

      /// default constructor 
      point() { x = y = 0.; }
      /// 2nd constructor 
      point(double _x, double _y)
      { x = _x; y = _y; }
      /// += operators
      point & operator+=(const point &pt) {
         this->x += pt.x;
         this->y += pt.y;
         return (*this);
      }
      point & operator+=(double c) {
         this->x += c;
         this->y += c;
         return (*this);
      }
      /// -= operators 
      point & operator-=(const point &pt) {
         this->x -= pt.x;
         this->y -= pt.y;
         return (*this);
      }
      point & operator-=(double c) {
         this->x -= c;
         this->y -= c;
         return (*this);
      }
      /// *= operator 
      point & operator*=(double c) {
         this->x *= c;
         this->y *= c;
         return (*this);
      }
      /// /= operator 
      point & operator/=(double c) {
         this->x /= c;
         this->y /= c;
         return (*this);
      }
      /// = operators
      point & operator=(double value) {
         this->x = value;
         this->y = value;
         return (*this);
      }
      point & operator=(const point &pt) {
         this->x = pt.x;
         this->y = pt.y;
         return (*this);
      }
      /// () operator 
      inline double & operator() (unsigned int i) {
         if (i == 0) return x;
         if (i > 0)  return y;
      }
   };



   /**
    * Class for modeling a function of time object
    */
   class TimeFunction
   {
   private:
      /// name of the timing function
      String name;

      /// vector of data points
      vector<point> points;

   public:
      /// Default constructeur
      TimeFunction();

      /// Destructor
      virtual ~TimeFunction();

      /// Copy constructor
      TimeFunction(const TimeFunction &fct);

      /// Sets the function name
      inline void Name(String s) {
         name = s;
      }

      /// Returns the function name
      inline String Name() {
         return name;
      }

      /// Sets the function size
      inline void Size(unsigned int size) {
         points.resize(size);
      }

      /// Returns the function size
      inline unsigned int Size() {
         return points.size();
      }

      /// Sets a data point
      inline void Point(unsigned int i, point &pt) {
         points[i] = pt;
      }

      /// Sets a data point
      inline void Point(unsigned int i, double &x, double &y) {
         points[i].x = x;
         points[i].y = y;
      }

      /// Returns a data point
      inline point& Point(unsigned int i) {
         return points[i];
      }

      /// Add a data point
      inline void Add(point &pt) {
         points.push_back(pt);
      }

      /// Add a data point
      inline void Add(double &x, double &y) {
         point pt; pt.x = x; pt.y = y;
         points.push_back(pt);
      }

      /// Inserts a data point a specified location
      inline void Insert(unsigned int i, point &pt) {
         points.insert(points.begin() + i, pt);
      }

      /// Remove a data point
      inline void Remove(unsigned int i) {
         points.erase(points.begin() + i);
      }

      /// Returns minimun abscissa value
      double MinX();

      /// Returns maximum abscissa value
      double MaxX();

      /// Returns minimun ordonate value
      double MinY();

      /// Returns maximum ordonate value
      double MaxY();

      /// Returns max point 
      point Max();

      /// Returns min point 
      point Min();

      /// Returns derivatives vector (i.e. dY/dX) 
      void Deriv(vector<double>& d);

      /// Overloding the equality operator
      TimeFunction& operator=(const TimeFunction &fct);

      /// Overloading the output stream operator
      friend ostream& operator <<(ostream& os, TimeFunction& fct);

      /// Overloading the input stream operator
      friend istream& operator >>(istream& is, TimeFunction& fct);
   };
}

#endif // #ifndef HAPP_TIME_FUNCTION
