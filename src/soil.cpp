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
#include <soil.h>
#include <algorithm>


namespace happ
{
   Soil::Soil() {
      hydrParam = hydraulic_parameters();
      stParam = solute_transport_parameters();
      htParam = heat_transport_parameters();
   }

   Soil::Soil(unsigned int Id)
   {
      id = Id;
      hydrParam = hydraulic_parameters();
      stParam = solute_transport_parameters();
      htParam = heat_transport_parameters();
   }

   Soil::Soil(String Name, unsigned int Id,
      const double _Kx, const double _Ky, const double _Kz,
      const double _S0) {
      name = Name;
      id = Id;
      hydrParam = hydraulic_parameters(_Kx, _Ky, _Kz, _S0);
   }

   Soil::Soil(const Soil& another_soil) {
      name = another_soil.name;
      id = another_soil.id;
      density = another_soil.density;
      hydrParam = another_soil.hydrParam;
      stParam = another_soil.stParam;
      htParam = another_soil.htParam;
   }

   Soil::~Soil() {
   }

   void Soil::SetKx(const double _Kx) {
      hydrParam.Kx = _Kx;
   }

   double Soil::GetKx() {
      return hydrParam.Kx;
   }

   void Soil::SetKy(const double _Ky) {
      hydrParam.Ky = _Ky;
   }

   double Soil::GetKy() {
      return hydrParam.Ky;
   }

   void Soil::SetKz(const double _Kz) {
      hydrParam.Kz = _Kz;
   }

   double Soil::GetKz() {
      return hydrParam.Kz;
   }

   void Soil::SetK(const double _K[3]) {
      hydrParam.Kx = _K[0];
      hydrParam.Ky = _K[1];
      hydrParam.Kz = _K[2];
   }

   void Soil::GetK(double k[3]) {
      k[0] = hydrParam.Kx;
      k[1] = hydrParam.Ky;
      k[2] = hydrParam.Kz;
   }

   void Soil::SetS0(const double _S0) {
      hydrParam.S0 = _S0;
   }

   double Soil::GetS0() {
      return hydrParam.S0;
   }

   void Soil::SetSatWC(const double _wc)
   {
      satWC = std::min(_wc, 1.);
   }

   double Soil::GetSatWC()
   {
      return satWC;
   }

   void Soil::SetResidualWC(const double _ResWc)
   {
      residualWC = std::min(_ResWc, std::min(satWC, 1.));
   }

   double Soil::GetResidualWC()
   {
      return residualWC;
   }

   void Soil::SetId(const unsigned int _id) {
      id = _id;
   }

   unsigned int Soil::GetId() {
      return id;
   }

   void Soil::operator = (const Soil &another_soil) {
      name = another_soil.name;
      id = another_soil.id;
      density = another_soil.density;
      hydrParam = another_soil.hydrParam;
      stParam = another_soil.stParam;
      htParam = another_soil.htParam;
   }
}
