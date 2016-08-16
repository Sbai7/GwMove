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
#include <periods.h>

// C++ headers
#include <iomanip>
#include <algorithm>


namespace happ
{
   TimeStepper::TimeStepper()
   {
      t0 = 0.0;
      theta = 1.0;
      time = t0;
      step = 0;
      period = 1;
      timeStep = 0.;
      end = false;
      maxTimeStep = 1.E+99;
      minTimeStep = 1.E-06;
      nBackSteps = 0;
   }

   TimeStepper::~TimeStepper()
   {
      Periods.free();
   }

   TimeStepper::TimeStepper(const TimeStepper &ts)
   {
      unsigned int i;

      t0 = ts.t0;
      theta = ts.theta;
      Periods.resize(ts.Periods.size());
      time = ts.time;
      timeStep = ts.timeStep;
      step = ts.step;
      period = ts.period;
      end = ts.end;
      maxTimeStep = ts.maxTimeStep;
      minTimeStep = ts.minTimeStep;
      nBackSteps = ts.nBackSteps;

      unsigned int n = ts.Periods.size();
      for (i = 0; i < n; i++)
         Periods[i] = ts.Periods[i];
   }

   void TimeStepper::Reset()
   {
      time = t0; 
      step = 0; 
      period = 1;
      nBackSteps = 0;
      end = false;
   }

    double TimeStepper::GetTotalTime() const
   {
      unsigned int i, n = Periods.size();
      double sum = 0.;
      for (i = 0; i < n; i++) {
         TimePeriod p = Periods[i];
         sum += p.length;
      }
      return sum;
   }

   unsigned int TimeStepper::GetNumberOfSteps()
   {
      valarray<unsigned int> a;
      GetStepsHistory(a);
      return a.sum();
   }

  void TimeStepper::GetStepsHistory(valarray<unsigned int> &hist)
   {
      unsigned int i;
      unsigned int n = Periods.size();
      hist.resize(n);
      for (i = 0; i < n; i++) hist[i] = Periods[i].steps;
   }

   void TimeStepper::GetTimeStepsHistory(TimeFunction &hist)
   {
      hist = timeStepHist;
   }

   void TimeStepper::AddPeriod(TimePeriod & tp)
   {
      unsigned int i, n = Periods.size();
      // make a copy of the time periods 
      valarray<TimePeriod> pBack; 
      pBack.resize(n);
      for (i = 0; i < n; i++) pBack[i] = Periods[i]; 

      // free current time periods array 
      Periods.free(); 

      // make new array
      Periods.resize(n+1); 
      for (i = 0; i < n; i++) Periods[i] = pBack[i];
      Periods[n] = tp; 

      // free backup array 
      pBack.free();
   }

   void TimeStepper::InsertPeriod(TimePeriod & tp, unsigned int position)
   {
      unsigned int i, n = Periods.size();
      if (position >= n) return;

      // make a copy of the time periods 
      valarray<TimePeriod> pBack;
      pBack.resize(n);
      for (i = 0; i < n; i++) pBack[i] = Periods[i];

      // free current time periods array 
      Periods.free();

      // make new array
      Periods.resize(n + 1);
      if (position == 0) {
         Periods[0] = tp;
         for (i = 0; i < n; i++) Periods[i+1] = pBack[i];
         return;
      }
      else {
         for (i = 0; i < position; i++) Periods[i] = pBack[i];
         Periods[position] = tp;
         for (i = position + 1; i < n + 1; i++) Periods[i] = pBack[i - 1];
      }

      // free backup array 
      pBack.free();
   }

   void TimeStepper::RemovePeriod(unsigned int position)
   {
      unsigned int i, n = Periods.size();
      if (position >= n || n < 2) return;

      // make a copy of the time periods 
      valarray<TimePeriod> pBack;
      pBack.resize(n);
      for (i = 0; i < n; i++) pBack[i] = Periods[i];

      // free current time periods array 
      Periods.free();

      // make new array
      Periods.resize(n - 1);
      if (position == 0) {
         for (i = 0; i < n-1; i++) Periods[i] = pBack[i + 1];
         return;
      }
      else {
         for (i = 0; i < position; i++) Periods[i] = pBack[i];
         for (i = position; i < n-1; i++) Periods[i] = pBack[i + 1];
      }

      // free backup array 
      pBack.free();
   }

   void TimeStepper::SplitPeriod(unsigned int position, double new_len)
   {
      unsigned int i, n = Periods.size();
      if (position >= n || new_len >= Periods[position].length) return;

      // make a copy of the time periods 
      valarray<TimePeriod> pBack;
      pBack.resize(n);
      for (i = 0; i < n; i++) pBack[i] = Periods[i];

      // free current time periods array 
      Periods.free();

      // make new array
      Periods.resize(n + 1);
      TimePeriod tp = pBack[position]; 
      tp.length = new_len;
      pBack[position].length -= new_len;
      if (position == 0) {
         Periods[0] = tp;
         for (i = 0; i < n; i++) Periods[i + 1] = pBack[i];
         return;
      }
      else {
         for (i = 0; i < position; i++) Periods[i] = pBack[i];
         Periods[position] = tp;
         for (i = position + 1; i < n + 1; i++) Periods[i] = pBack[i - 1];
      }

      // free backup array 
      pBack.free();
   }

   void TimeStepper::MergePeriods(unsigned int position)
   {
      unsigned int i, n = Periods.size();
      if (position >= n - 1 || n < 2) return;

      // make a copy of the time periods 
      valarray<TimePeriod> pBack;
      pBack.resize(n);
      for (i = 0; i < n; i++) pBack[i] = Periods[i];

      // free current time periods array 
      Periods.free();

      // make new array
      Periods.resize(n - 1);

      // adjusts merged period length 
      pBack[position + 1].length += pBack[position].length;

      // decision on internal time stepping of the merged period 
      if (pBack[position + 1].auto_step || pBack[position].auto_step) 
         pBack[position + 1].auto_step = true; 
      if (!pBack[position].auto_step && !pBack[position + 1].auto_step) {
         pBack[position + 1].auto_step = false;
         pBack[position + 1].steps += pBack[position].steps;
      }

      if (position == 0) {
         for (i = 0; i < n - 1; i++) Periods[i] = pBack[i + 1];
         return;
      }
      else {
         for (i = 0; i < position; i++) Periods[i] = pBack[i];
         for (i = position; i < n - 1; i++) Periods[i] = pBack[i + 1];
      }

      // free backup array 
      pBack.free();
   }

   double TimeStepper::GetNextOutTime()
   {
      unsigned int i;
      double sum = 0.;
      for (i = 0; i < period; i++) {
         sum += Periods[i].length;
      }
      return sum;
   }

   void TimeStepper::BoundTimeStep()
   {
      timeStep = std::max(timeStep, minTimeStep);
      timeStep = std::min(timeStep, maxTimeStep);
   }

   TimeStepper& TimeStepper::operator=(const TimeStepper &sp)
   {
      t0 = sp.t0;
      theta = sp.theta;
      Periods.resize(sp.Periods.size());
      time = sp.time;
      timeStep = sp.timeStep;
      step = sp.step;
      period = sp.period;
      end = sp.end;
      maxTimeStep = sp.maxTimeStep;
      minTimeStep = sp.minTimeStep;
      nBackSteps = sp.nBackSteps;

      unsigned int i, n = sp.Periods.size();
      for (i = 0; i < n; i++) Periods[i] = sp.Periods[i];
      return *this;
   }

   void TimeStepper::operator ++()
   {
      bool inPeriod = (Periods[period - 1].auto_step && time < GetNextOutTime()) || 
         (!Periods[period - 1].auto_step && step < Periods[period - 1].steps);
      
      if (inPeriod) {
         // we're still in the current time period
         step++;
         if (Periods[period - 1].auto_step) {
            timeStep = std::min(timeStep, GetNextOutTime() - time);
            BoundTimeStep();
            Periods[period - 1].steps = step;
         }
         else {
            timeStep = Periods[period - 1].length / Periods[period - 1].steps;
            BoundTimeStep();
         }
         time += timeStep;
      }

      else {
         // check for the last stress period 
         if (period == Periods.size()) {
            end = true;
            return;
         }

         period++; // goto the next time period
         step = 1;
         if (Periods[period - 1].auto_step) {
            timeStep = std::min(timeStep, GetNextOutTime() - time);
            BoundTimeStep();
         } 
         else {
            timeStep = Periods[period - 1].length / Periods[period - 1].steps;
            BoundTimeStep();
         }
         time += timeStep;
      }
      
      timeStepHist.Add(time, timeStep);
   }

   void TimeStepper::operator --()
   {
      if (time == t0) return;

      time -= timeStep;
      //step++;
      nBackSteps++;
      //timeStepHist.Add(time, timeStep);
   }

   void TimeStepper::operator *=(double &multiplier)
   {
      if (time != GetNextOutTime()) {
         timeStep = std::min(fabs(timeStep*multiplier),
            GetNextOutTime() - time);
         BoundTimeStep();
      }
   }

   void TimeStepper::operator /=(double &divider)
   {
      if (time != GetNextOutTime()) {
         timeStep = std::min(fabs(timeStep / divider),
            GetNextOutTime() - time);
         BoundTimeStep();
      }
   }

   void TimeStepper::operator +=(double &increment)
   {
      if (time != GetNextOutTime())
         timeStep = std::min(fabs(timeStep + increment), 
                             GetNextOutTime() - time);
   }

   void TimeStepper::operator -=(double &decrement)
   {
      if (time != GetNextOutTime()) {
         timeStep = std::min(fabs(timeStep - decrement),
            GetNextOutTime() - time);
         BoundTimeStep();
      }
   }

   TimePeriod& TimeStepper::operator [](unsigned int &i)
   {
      return Periods[i];
   }

   ostream& operator <<(ostream& os, TimePeriod& per)
   {
      os << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(4);
      os << setw(15) << per.length 
         << setw(5) << per.auto_step
         << setw(5) << per.output
         << setw(8) << per.steps << endl;
      return os;
   }

   ostream& operator <<(ostream& os, TimeStepper& temp)
   {
      os << "t0    = " << temp.t0    << endl;
      os << "theta = " << temp.theta << endl;

      unsigned int i, n = temp.Periods.size();
      os << "time periods = " << n << endl;
      //     1234567890123456789012345678901234567890
      os << "# ------ length - auto - out - steps -" << endl;
      for (i = 0; i < n; i++) os << temp.Periods[i];
      os << endl;

      return os;
   }

   istream& operator >>(istream& is, TimePeriod& per)
   {
      int Output, AutoStep;
      char d[120];

      is >> per.length >> AutoStep >> Output >> per.steps;
      per.output = (Output == 1 ? true : false);
      per.auto_step = (AutoStep == 1 ? true : false);

      is.getline(d, 120, '\n');
      return is;
   }

   istream& operator >>(istream& is, TimeStepper& temp)
   {
      char d[120];
      TimePeriod per;
      unsigned int i, nper;

      is >> temp.t0 >> temp.theta;
      is.getline(d, 120, '\n');

      is.getline(d, 120, '\n');

      is >> nper;
      temp.Periods.resize(nper);
      is.getline(d, 120, '\n');

      for (i = 0; i < nper; i++) {
         is >> per;
         temp.Periods[i] = per;
      }

      return is;
   }

}
