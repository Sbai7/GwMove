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
#include <bndcondition.h>

// C++ headers
#include <iomanip>


namespace happ
{
   BndCondition::BndCondition(StrucMesh *pMesh)
   {
      // no boundary condition is specified yet
      nBC = 0;

      // set mesh pointer
      mesh = pMesh;

      // set dimension
      unsigned int nn = mesh->GetNnodes();
      dimension = nn;

      // allocate memory for internal arrays
      bcType = new unsigned int[dimension];
      bcValue = new double[4 * dimension];
      fixedHead = new bool[dimension];
      active = new unsigned int[dimension];

      // initialize arrays
      unsigned int i;
      for (i = 0; i < dimension; i++) {
         bcType[i] = NO_BC;
         fixedHead[i] = false;
         active[i] = true;
      }
      for (i = 0; i < 4 * dimension; i++)
         bcValue[i] = 0.;

      // calculates nodal surfaces on top surface
      mesh->ComputeNodalSurfaces();
   }

   BndCondition::~BndCondition()
   {
      unsigned int i, n = bcZone.size();
      for (i = 0; i < n; i++) delete bcZone[i];
      bcZone.clear();
      if (bcType)    delete bcType;
      if (bcValue)   delete bcValue;
      if (fixedHead) delete fixedHead;
      if (active)    delete active;
   }

   void BndCondition::Add(const BcZone zone)
   {
      unsigned int i, j, k, n;
      unsigned int number;
      unsigned int ni = mesh->GetNi();
      unsigned int nj = mesh->GetNj();

      for (k = zone.K1; k <= zone.K2; k++) {
         for (j = zone.J1; j <= zone.J2; j++) {
            for (i = zone.I1; i <= zone.I2; i++) {
               number = i + (j - 1)*ni + (k - 1)*ni*nj;
               bcType[number - 1] = zone.bcType;

               // set fixed head/flow bc's
               for (n = 0; n < 4; n++) bcValue[4 * (number - 1) + n] = zone.dValues[n];

               // set flux for top surface bc's 
               // multiply by nodal surface area to obtain a flow rate 
               if (zone.bcType == FIXED_FLUX_BC || zone.bcType == INFILTRATION_BC) {
                  bcValue[4 * (number - 1)] *= mesh->surface[number - 1];
               }
               else if (zone.bcType == LEAKAGE_BC) {
                  bcValue[4 * (number - 1) + 1] *= mesh->surface[number - 1];
               }
               else if (zone.bcType == RECHARGE_BC || zone.bcType == EVAPOTRANSPIRATION_BC) {
                  bcValue[4 * (number - 1) + 1] *= mesh->surface[number - 1];
                  bcValue[4 * (number - 1) + 3] *= mesh->surface[number - 1];
               }
            }
         }
      }
   }

   void BndCondition::Add(const ijkZone zone, unsigned int type, double val[4])
   {
      unsigned int i, j, k, n;
      unsigned int number;
      unsigned int ni = mesh->GetNi();
      unsigned int nj = mesh->GetNj();

      for (k = zone.K1; k <= zone.K2; k++) {
         for (j = zone.J1; j <= zone.J2; j++) {
            for (i = zone.I1; i <= zone.I2; i++) {
               number = i + (j - 1)*ni + (k - 1)*ni*nj;
               bcType[number - 1] = type;

               // set fixed head/flow bc's
               for (n = 0; n < 4; n++) bcValue[4 * (number - 1) + n] = val[n];

               // set flux for top surface bc's 
               // multiply by nodal surface area to obtain a flow rate 
               if (type == FIXED_FLUX_BC || type == INFILTRATION_BC) {
                  bcValue[4 * (number - 1)] *= mesh->surface[number - 1];
               }
               else if (type == LEAKAGE_BC) {
                  bcValue[4 * (number - 1) + 1] *= mesh->surface[number - 1];
               }
               else if (type == RECHARGE_BC || type == EVAPOTRANSPIRATION_BC) {
                  bcValue[4 * (number - 1) + 1] *= mesh->surface[number - 1];
                  bcValue[4 * (number - 1) + 3] *= mesh->surface[number - 1];
               }
            }
         }
      }
   }

   ostream& operator <<(ostream& os, BndCondition& bc)
   {
      unsigned int i, j;

      os << bc.nBC << endl;

      for (i = 0; i < bc.dimension; i++) {
         if (bc.bcType[i] > NO_BC) {
            os << setw(8) << bc.bcType[i];
            for (j = 0; j < 4; j++)
               os << setw(16) << bc.bcValue[4 * i + j];
            os << endl;
         }
      }

      return os;
   }

   bool BndCondition::ExtractAtBC(Vector& v) const 
   {
      if (v.Size() != dimension) return false; 

      unsigned int m = this->GetnBC(); // number of BC nodes 
      double* newVec = new double[m];

      unsigned int i, j = 0;
      for (i = 0; i < dimension; i++) {
         if (bcType[i] > NO_BC && bcType[i] <= FLOOD_BC) {
            newVec[j] = v[i];
            j++;
         }
      }

      v.NewDataAndSize(newVec, m);
      return true;
   }

   bool BndCondition::Read(String projectFile)
   {
      filebuf file;
      String bcFile = projectFile + ".fbc";

      file.open(bcFile.GetBuffer(), ios::in);
      istream is(&file);
      FileReader fr(is);

      int j;
      String label;

      //int nTokens = fr.GetLine();
      //if (nTokens < 1) return false;
      //nBC = atoi(fr.GetToken(0));
      //bcZone.reserve(nBC);

      do {
         int nTokens = fr.GetLine();
         if (nTokens < 5) {
            // print error message indicating the line where it occurs
            return false;
         } 
         // end of file 
         else if (nTokens == 0 ) return true;

         nBC++;

         unsigned int type = atoi(fr.GetToken(0));
         unsigned int i1 = atoi(fr.GetToken(1));
         unsigned int j1 = atoi(fr.GetToken(2));
         unsigned int k1 = atoi(fr.GetToken(3));
         unsigned int i2, j2, k2, v_index;
         double value[4];
         for (j = 0; j < 4; j++) value[j] = 0.;

         i2 = atoi(fr.GetToken(4));
         if (i2 == 0) {
            j2 = j1;
            k2 = k1;
            v_index = 5;
         }
         else {
            j2 = atoi(fr.GetToken(5));
            k2 = atoi(fr.GetToken(6));
            v_index = 7;
         }

         if (type == FIXED_HEAD_BC || type == FIXED_FLOW_BC ||
            type == FIXED_FLUX_BC || type == DRAINAGE_BC) {
            value[0] = atof(fr.GetToken(v_index)); v_index++;
         }

         else if (type == RIVER_BC || type == LEAKAGE_BC ||
            type == INFILTRATION_BC || type == ABSTRACTION_BC) {
            value[0] = atof(fr.GetToken(v_index)); v_index++;
            value[1] = atof(fr.GetToken(v_index)); v_index++;
         }

         else if (type == RECHARGE_BC || type == EVAPOTRANSPIRATION_BC) {
            value[0] = atof(fr.GetToken(v_index)); v_index++;
            value[1] = atof(fr.GetToken(v_index)); v_index++; 
            value[2] = atof(fr.GetToken(v_index)); v_index++;
            value[3] = atof(fr.GetToken(v_index)); v_index++;
         }

         std::vector<double> val;
         for (j = 0; j < 4; j++) val.push_back(value[j]);

         // try to get the label of this BC if set by the user
         unsigned int n = fr.GetNtokens();
         label = "";
         if (n > v_index) {
            for (j = v_index; j < n; j++) {
               label += fr.GetToken(j);
               label += " ";
            }
         }
         else {
            label  = "BC zone # "; 
            label += (int) nBC; 
         }

         // add this oundary condition zone
         ijkZone ijkZ(i1, j1, k1, i2, j2, k2);
         ijkZ.Name(label);

         BcZone *zone = new BcZone(ijkZ, type, val);
         bcZone.push_back(zone);
         Add(*zone);
      } while (true);

      return true;
   }
}
