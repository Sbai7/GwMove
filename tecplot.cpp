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
#include <tecplot.h>

// C++ headers
#include <iomanip>


namespace happ
{
   TecplotIO::TecplotIO()
   {
      m_datapacking = "POINT";         // POINT (column) nodal data by default 
      m_zonetype = ZONETYPE::FEBRICK;  // Hexahedron by default 
      m_binary   = false;              // ASCII file by default  
      m_zonetype = ZONETYPE::ORDERED;  // Ordered IJK zone by default 
      m_meshShared = false; 
      m_zCoordShared = false;
      m_fileOpened = false;
      m_zonetitle = "Zone I";          // Default title of first zone 
      m_filename = "Output.plt";       // Default output file 
   }

   TecplotIO::TecplotIO(StrucMesh *pmesh)
   {
      mesh = pmesh;
      m_nvariables = 3;
      m_svariables = "\"X\", \"Y\", \"Z\"";
      m_currentzone = 1;
      if (mesh->GetDim() == 2) m_zonetype = ZONETYPE::FEQUADRILATERAL;
      if (mesh->GetDim() == 3) m_zonetype = ZONETYPE::FEBRICK;
      m_datapacking = "POINT";
      m_binary = false;
      m_zonetype = ZONETYPE::ORDERED;
      m_meshShared = false;
      m_zCoordShared = false;
      m_fileOpened = false;
      m_zonetitle = "Zone I";
      m_filename = "Output.plt"; 
   }

   TecplotIO::~TecplotIO()
   {
      variables.clear();
   }

   bool TecplotIO::AddDataset(double *var, const String name)
   {
      variables.push_back(var);
      // add variable name to m_svariables string
      m_svariables += String(", \"");
      m_svariables += name;
      m_svariables += String("\"");
      m_nvariables++;
      return true;
   }

   void TecplotIO::UpdateDataset(double *var, const unsigned int j)
   {
      if (j >= variables.size()) return;
      double* thisVar = variables[j];
      unsigned int i, nn = mesh->GetNnodes(); 
      for (i = 0; i < nn; i++) thisVar[i] = var[i]; 
   }

   bool TecplotIO::RemoveDataset(const unsigned int i)
   {
      if (i >= variables.size())
         return false;
      variables.erase(variables.begin() + i);
      // remove variable title from m_svariables string
      int index = m_svariables.Find(',', -1);
      m_svariables = m_svariables.SubString(0, index);
      m_nvariables--;
      return true;
   }

   void TecplotIO::InitNextZone(String title)
   {
     m_currentzone++;
     m_zonetitle = title;       
   }

   bool TecplotIO::Write()
   {
      if (!m_fileOpened) outfile.open(m_filename.GetBuffer(), ios::out);
      ostream os(&outfile);
      os << *this; 
      
      return true;
   }

   bool TecplotIO::Write(String filename)
   {
      if (m_fileOpened) this->Dispose();
      m_filename = filename;
      return this->Write();
   }


   bool TecplotIO::Dispose()
   {
      if (m_fileOpened) {
         outfile.close();
         m_currentzone = 1;
         return true;
      }
      else
         return false;
   }

   bool TecplotIO::WriteWaterFile(String filename)
   {
      std::filebuf waterFile;
      waterFile.open(filename.GetBuffer(), ios::out);
      ostream os(&waterFile);

      if (variables.size() < 2) {
         waterFile.close();
         return false;
      }

      // number of nodes 
      os << setw(8) << mesh->GetNnodes() << "   # Number of nodes" << endl; 
      os << setw(8) << variables.size()  << "   # number of variables" << endl; 
      os << "# Variable 1 = groundwater head [L]" << endl; 
      os << "# Variable 2 = volumetric flow rate[L^3/T]" << endl;

      os << scientific;
      unsigned int i,j;
      for (i = 0; i < mesh->GetNnodes(); i++) {
         for (j = 0; j < variables.size(); j++) {
            double *var = variables[j];
            os << setw(16) << var[i];
         }
         os << endl;
      }

      waterFile.close();

      return true; 
   }


   ostream& operator << (std::ostream & os, TecplotIO& tec)
   {
      unsigned int i, j;
      if (!tec.m_binary) { // 1 -> ASCII TecFile

         // Write tecplot header file
         if (tec.m_currentzone == 1) {
            os << "TITLE = \"" << tec.m_title << "\"" << endl;
            os << "VARIABLES = " << tec.m_svariables << endl;
         }

         /* Write the unique zone to output stream */

         // 1- zone header 
         os << "ZONE T=\"" << tec.m_zonetitle << "\"";
         if (tec.m_zonetype != ZONETYPE::ORDERED) {
            os << ", NODES="   << tec.mesh->GetNnodes()
               << ", ELEMENTS=" << tec.mesh->GetNelements();
         }
         else {
            os << "I="   << tec.mesh->GetNi() 
               << ", J=" << tec.mesh->GetNj()
               << ", K=" << tec.mesh->GetNslices();
         }
         os << ", DATAPACKING=" << tec.m_datapacking
            << ", ZONETYPE=" << tec.GetZoneTypeStr();
         if (tec.m_meshShared) os << ", VARSHARELIST=([1-3]=1)"; 
         if (tec.m_zCoordShared) os << ", VARSHARELIST=([1-2]=1)";
         os << endl;

         // 2- zone core section: write mesh coordinates + variables
         if (tec.m_datapacking == "POINT") { // in POINT format ... 

            for (i = 0; i < tec.mesh->GetNnodes(); i++)  {
               double x[3];
               tec.mesh->GetCoord(i, x);
               os << scientific;

               if (!tec.m_meshShared) {
                  if (tec.m_zCoordShared)
                     os << setw(16) << x[0]
                        << setw(16) << x[1];
                  else 
                     os << setw(16) << x[0]
                        << setw(16) << x[1]
                        << setw(16) << x[2];
               }

               for (j = 0; j < tec.variables.size(); j++) {
                  double *var = tec.variables[j];
                  os << setw(16) << var[i];
               }

               os << endl;
            }
         }

         else { // in BLOCK format ... 
            
            double* X = new double[tec.mesh->GetNnodes()];
            double* Y = new double[tec.mesh->GetNnodes()];
            double* Z = new double[tec.mesh->GetNnodes()];
            
            tec.mesh->GetX(X); tec.mesh->GetY(Y); tec.mesh->GetZ(Z); 
            
            if (!tec.m_meshShared) {
               
               // X data block 
               for (i = 0; i < tec.mesh->GetNnodes(); i++) {
                  os << setw(16) << X[i];
                  if (i % 500 == 0) os << endl;
               }
               os << endl;
               
               // Y data block 
               for (i = 0; i < tec.mesh->GetNnodes(); i++) {
                  os << setw(16) << Y[i];
                  if (i % 500 == 0) os << endl;
               }
               os << endl;
                   
               // Z data block 
              if (!tec.m_zCoordShared) {
                  for (i = 0; i < tec.mesh->GetNnodes(); i++) {
                     os << setw(16) << Z[i];
                     if (i % 500 == 0) os << endl;
                  }
                  os << endl;
               }
            }
            
            delete[] X; 
            delete[] Y;
            delete[] Z;

            // then other variables ...
            for (j = 0; j < tec.variables.size(); j++) {
               double *var = tec.variables[j];
               for (i = 0; i < tec.mesh->GetNnodes(); i++)
               {
                  os << setw(16) << var[i];
                  if (i % 500 == 0) os << endl;
               }
               os << endl;
            }
         }

         // 3- write elements-nodes connectivity list
         if (tec.m_zonetype != ZONETYPE::ORDERED) {
            unsigned int np;
            if (tec.m_zonetype == ZONETYPE::FEQUADRILATERAL) np = 4;
            if (tec.m_zonetype == ZONETYPE::FEBRICK) np = 8;
            for (i = 0; i < tec.mesh->GetNelements(); i++) {
               unsigned int conn[8]; // WORKS ONLY IN 3D CASE NOW 
               tec.mesh->GetConnectivity(i, conn);
               for (j = 0; j < np; j++)
                  os << setw(8) << conn[j] + 1;
               os << endl;
            }
         }
         os << endl;
      }
      return os;
   }
}

