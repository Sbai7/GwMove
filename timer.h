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

#ifndef HAPP_TIMER
#define HAPP_TIMER

// C headers
#include <time.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/timeb.h>

// C++ headers
#include <iostream>
using namespace std;


namespace happ
{
   /**
    * This class may be used to calculate execution times in seconds.
    */
   class Timer
   {
   private:
      time_t t0, t1, ct;
      bool st0, st1;

   public:

      /// Default constructor. Initialize time to zero
      inline Timer() {
         ct = time_t(0);
         st0 = st1 = false;
      }

      /// Start time computation
      inline void Start() {
         time(&t0);
         st0 = true;
      }

      /// Stop time computation
      inline void Stop() {
         if (!st0) {
            cerr << " >> Error: Method Start() must have been called before." << endl;
            exit(1);
         }
         time(&t1);
         ct += t1 - t0;
         st1 = true;
      }

      inline void Resume() {
         if (!st0) {
            cerr << " >> Error: Method Start() must have been called before." << endl;
            exit(1);
         }
         if (!st1) {
            cerr << " >> Error: Method Stop() must have been called before." << endl;
            cerr << "Call method 'Start()'\n";
            exit(1);
         }
         time(&t1);
      }

      /// Restart computation
      inline void Clear() {
         time(&t0);
         ct = time_t(0);
      }

      /// Return elapsed time in seconds
      inline time_t ClockTime() {
         if (!st0) {
            cerr << "  >> Error: Method Start() must have been called before." << endl;
            exit(1);
         }
         if (!st1) {
            cerr << " >> Error: Method Stop() must have been called before." << endl;
            exit(1);
         }
         return ct;
      }

   };
}

#endif // #ifndef HAPP_TIMER
