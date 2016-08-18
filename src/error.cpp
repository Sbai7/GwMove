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

#include <error.h>

#include <cstdlib>
#include <iostream>

#ifdef HAPP_USE_MPI
#include <mpi.h>
#endif

namespace happ
{

   void happ_error(const char *msg)
   {
      if (msg)
      {
         // NOTE: This endl also flushes the I/O stream, which can be a very bad
         // thing if all your processors try to do it at the same time.
         std::cerr << "\n\n" << msg << std::endl;
      }
#ifdef HAPP_USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#else
      std::abort(); // force crash by calling abort
#endif
   }

   void happ_warning(const char *msg)
   {
      if (msg)
      {
         std::cout << "\n\n" << msg << std::endl;
      }
   }

}
