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

#ifndef HAPP_PERIODS
#define HAPP_PERIODS

// Project headers
#include<function.h>

// C++ headers
#include <iostream>
#include <valarray>
#include <vector>
using namespace std;


namespace happ
{
   /**
    * Data structure holding the number of time steps and their length
    */
   struct TimePeriod
   {
      /// total time simulated by this period
      double length;

      /// automatic time step control flag
      bool auto_step;

      /// output results -at the period end- flag
      bool output;

      /// number of time steps - when auto. control disabled-
      unsigned int steps;

      /// overloading equality operator
      TimePeriod&  operator=(const TimePeriod &p) {
         output = p.output;
         length = p.length;
         auto_step = p.auto_step;
         steps = p.steps;

         return *this;
      };

      // 1st constructor 
      TimePeriod() {}

      // 2nd constructor 
      TimePeriod(const double _length, 
                 const bool _auto_step, 
                 const bool _output, 
                 const unsigned int _steps=0) {
         length = _length;
         auto_step = _auto_step; 
         output = _output;
         steps = _steps;
      }

      // destructor 
      virtual ~TimePeriod() {} 

      /// Copy constructor 
      TimePeriod(const TimePeriod & tp) {
         length = tp.length; 
         auto_step = tp.auto_step;
         output = tp.output;
         steps = tp.steps;
      }

   };

   /// Period's output operator
   ostream& operator <<(ostream& os, TimePeriod& per);

   /// Period's input operator
   istream& operator>>(istream& is, TimePeriod& per);


   /**
    * Time stepping class for management of time steps during
    * each stress period.
    */
   class TimeStepper
   {
   private:
      /// initial simulation time
      double t0;

      /// weighting factor for time marching scheme
      double theta;

      /// Current simulation time
      double time;

      /// current simulation time step
      double timeStep;

      /// current step
      unsigned int step;

      /// current period
      unsigned int period;

      /// is this is the last step in this time period ?
      bool lastStep = false;

      /// is it the end of time stepping sequence 
      bool end = false;

      /// maximal allowed time step 
      double maxTimeStep;

      /// minimal allowed time step 
      double minTimeStep; 

      /// number of backstepping events 
      unsigned int nBackSteps;

      /// time stepping history array 
      TimeFunction timeStepHist;

      /// Calculate total time until the end of the current period
      double GetNextOutTime();

      /// Make time step in acceptable min- max- bounds
      void BoundTimeStep();

   public:
      /// stress periods vector
      valarray<TimePeriod> Periods;

      /// Default constructor
      TimeStepper();

      /// Copy constructor
      TimeStepper(const TimeStepper &sp);

      /// destructor
      virtual ~TimeStepper();

      /// Returns initial simulation time
      inline double GetT0() const {
         return t0;
      }

      /// Sets initial simulation time
      inline void SetT0(const double _t0) {
         t0 = _t0;
      }

      /// Returns weighting factor
      inline double GetTheta() const {
         return theta;
      }

      /// Sets weighting factor
      inline void SetTheta(const double _theta) {
         theta = _theta;
      }

      /// Returns current time
      inline double GetTime() const {
         return time;
      }

      /// Returns current time step 
      inline double GetTimeStep() const {
         return timeStep;
      }

      /// Sets current time step 
      inline void SetTimeStep(const double _timeStep) {
         timeStep = _timeStep;
      }

      /// Sets max time step 
      inline void SetMaxTimeStep(const double _maxDt) {
         maxTimeStep = _maxDt;
      }

      /// Sets min time step 
      inline void SetMinTimeStep(const double _minDt) {
         minTimeStep = _minDt;
      }

      /// Returns current step
      inline unsigned int GetStep() const {
         return step;
      }

      /// Returns current period
      inline unsigned int GetPeriod() const {
         return period;
      }

      /// Signals the end of the time stepping sequence 
      inline bool End() const {
         return end; 
      }

      /// Signals the end of the current time period 
      inline bool EndPeriod() {
         if (time == GetNextOutTime()) return true;
         else return false;
      }

      /// Signals if time period output is user-requested 
      inline bool IsOutput() {
         return (EndPeriod() && Periods[period - 1].output); 
      }

      /// Returns the number of backstepping events 
      inline unsigned int GetBackSteps() const {
         return nBackSteps;
      }

      /// Resets time stepping internal counters & pointers 
      void Reset(); 

      /// Returns the total length of transient simulation 
      double GetTotalTime() const;

      /// Returns cumulative number of time steps
      unsigned int GetNumberOfSteps();

      /// Adds a time period to time stepping scheme 
      void AddPeriod(TimePeriod & tp);

      /// inserts a stress period 
      void InsertPeriod(TimePeriod & tp, unsigned int position); 

      /// removes a stress period 
      void RemovePeriod(unsigned int position);

      /// splits stress period into two 
      void SplitPeriod(unsigned int position, double new_len);

      /// concatenate two successive stress periods starting from position  
      void MergePeriods(unsigned int position);

      /// Returns number of steps per period
      void GetStepsHistory(valarray<unsigned int> &hist);

      /// Returns time stepping history
      void GetTimeStepsHistory(TimeFunction &hist);

      /// overloading equality operator
      TimeStepper&  operator=(const TimeStepper &sp);

      /// overloading ++ operator such that we skip to the next time step
      void operator ++();

      /// overloading -- operator to do back-stepping
      void operator --();

      /// overloading *= operator to multiply time step
      void operator *=(double &multiplier);

      /// overloading /= operator to divide time step
      void operator /=(double &divider);

      /// overloading += to increase time step by a constant
      void operator +=(double &increment);

      /// overloading -= to decrease time step by a constant
      void operator -=(double &decrement);

      /// overloading [] operator to access a stress period
      TimePeriod& operator [](unsigned int &i);

      /// overloading of output stream operator
      friend ostream& operator <<(ostream& os, TimeStepper& temp);

      /// overloading of input stream operator
      friend istream& operator >>(istream& is, TimeStepper& temp);
   };
}

#endif // HAPP_PERIODS
