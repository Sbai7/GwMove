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

#ifndef HAPP_SOIL
#define HAPP_SOIL

// Project headers
#include <AdString.h>
//#include <UnsatSoilCurves.h>

// C++ headers 
#include <vector>
using namespace std;


namespace happ
{
   struct hydraulic_parameters {
      /// Entries of diagonal hydraulic conductivity tensor
      double Kx, Ky, Kz;
      /// Off-diagonals of symmetric hydraulic conductivity tensor
      double Kxy, Kxz, Kyz;
      /// Entries of diagonal absolute permeability tensor 
      double absKx, absKy, absKz;
      /// Off-diagonals of symmetric absolute permeability tensor
      double absKxy, absKxz, absKyz;
      /// Rock compressibility 
      double compressibility;
      /// Hydraulic conductivity variance for isotropic materials
      double variance;
      /// Saturated specific storage coefficient [1/m]
      double S0;
      /// Unconfined specific storage or specific yield 
      double Sy;

      /// Default constructor
      hydraulic_parameters() {
         // isotropic permeability of ~ 1 darcy 
         Kx = 1e-12; Ky = 1e-12; Kz = 1e-12;
         Kxy = 0.;    Kxz = 0.;    Kyz = 0.;
         // Absolute permeabilities calculated for freshwater 
         // @ standard pressure, temperature, salinity conditions 

         absKxy = 0.; absKxz = 0.; absKz = 0.;
         // 
         S0 = 1e-6; Sy = 0.1;
      }

      /// 2nd constructor 
      hydraulic_parameters(const double _Kx, const double _Ky, const double _Kz,
         const double _S0 = 0) {
         Kx = _Kx;  Ky = _Ky;   Kz = _Kz;
         S0 = _S0;
      }

      /// 3rd constructor 
      hydraulic_parameters(const double _Kx, const double _Ky, const double _Kz,
         const double _Kxy, const double _Kxz, const double _Kyz,
         const double _S0 = 0) {
         Kx = _Kx;  Ky = _Ky;   Kz = _Kz;
         Kxy = _Kxy; Kxz = _Kxz; Kyz = _Kyz;
         S0 = _S0;
      }

      /// 4th constructor 
      hydraulic_parameters(const std::vector<double> & _K, const double _S0 = 0) {
         unsigned int size = _K.size();
         if (size == 3) {
            Kx = _K[0]; Ky = _K[1]; Kz = _K[2];
         }
         else if (size == 6) {
            Kx = _K[0]; Ky = _K[1]; Kz = _K[2];
            Kxy = _K[3]; Kxz = _K[4]; Kyz = _K[5];
         }
         else {
            // assertion failure 
         }
         S0 = _S0;
      }

      /// Destructor 
      virtual ~hydraulic_parameters() {};

      /// Calculate absolute permeability components from 
      /// those of hydraulic conductivity 
      void Cond2Perm(const double _density, const double _viscosity) {
         const double g = 9.81; // en m.s-2 ... to refine later
         double f = _viscosity / (_density*g);
         absKx = f*Kx;  absKy = f*Ky;  absKz = f*Kz;
         absKxy = f*Kxy; absKxz = f*Kxz; absKyz = f*Kyz;
      }

      /// Calculate hydraulic conductivity components from 
      /// those of hydraulic conductivity 
      void Perm2Cond(const double _density, const double _viscosity) {
         const double g = 9.81; // en m.s-2 ... to refine later
         double f = (_density*g) / _viscosity;
         Kx = f*absKx;  Ky = f*absKy;  Kz = f*absKz;
         Kxy = f*absKxy; Kxz = f*absKxz; Kyz = f*absKyz;
      }

      // methods to identify isotropy of K ...  
   };

   struct solute_transport_parameters {
      /// Solute dispersivities in logitudinal & transverse directions
      double dispX, dispY, dispZ;
      /// Effective porosity 
      double porosity;
      /// Mobile fraction of pore water ( a value in ]0,1])
      double mobileFrac;
      /// Retardation factor for a decaying solute 
      double Rfactor;
      /// Molecular diffusion coefficient of a solute specie 
      double diffusion;

      /// Default constructor 
      solute_transport_parameters() {
         dispX = 10.;
         dispY = 1.;
         dispZ = 1.;
         porosity = 0.2;
         mobileFrac = 1.;
         Rfactor = 1.;
         diffusion = 1e-9;
      }

      /// 2nd constructor 
      solute_transport_parameters(const double _dispX, const double _dispY,
         const double _dispZ, const double _porosity, const double _diffusion = 0.,
         const double _mobileFrac = 0.) {
         dispX = _dispX; dispY = _dispY; dispZ = _dispZ;
         porosity = _porosity;
         diffusion = _diffusion;
         mobileFrac = _mobileFrac;
      }

      /// 3rd constructor 
      solute_transport_parameters(const double _disp[3], const double _porosity,
         const double _diffusion = 0., const double _mobileFrac = 0.) {
         dispX = _disp[0]; dispY = _disp[1]; dispZ = _disp[2];
         porosity = _porosity;
         diffusion = _diffusion;
         mobileFrac = _mobileFrac;
      }

      /// Destructor 
      virtual ~solute_transport_parameters() {};
   };

   struct heat_transport_parameters {
      /// Heat dispersivities in logitudinal & transverse directions
      double dispX, dispY, dispZ;
      /// Effective porosity 
      double porosity;
      /// Heat diffusion coefficient
      double diffusion;
      /// Specific heat capacity 
      double heatCapacity;

      /// Default constructor 
      heat_transport_parameters() {
         dispX = 10.;
         dispY = 1.;
         dispZ = 1.;
         porosity = 0.2;
         diffusion = 1e-6;
         heatCapacity = 2.5e+6; // ... verify 
      }

      /// 2nd constructor 
      heat_transport_parameters(const double _dispX, const double _dispY,
         const double _dispZ, const double _porosity, const double _diffusion = 0.,
         const double _heatCapacity = 0.) {
         dispX = _dispX; dispY = _dispY; dispZ = _dispZ;
         porosity = _porosity;
         diffusion = _diffusion;
         heatCapacity = _heatCapacity;
      }

      /// 3rd constructor 
      heat_transport_parameters(const double _disp[3], const double _porosity,
         const double _diffusion = 0., const double _heatCapacity = 0.) {
         dispX = _disp[0]; dispY = _disp[1]; dispZ = _disp[2];
         porosity = _porosity;
         diffusion = _diffusion;
         heatCapacity = _heatCapacity;
      }
   };

   /**
    * Soil type definition class.
    */
   class Soil
   {
   private:
      /// Soil type id used in FE computations
      unsigned int id;

      /// Hydraulic parameters
      hydraulic_parameters hydrParam;

      /// Solute transport parameters 
      solute_transport_parameters stParam;

      /// Heat transport parameters
      heat_transport_parameters htParam;

      /// Saturated water content
      double satWC;

      /// Residual water content
      double residualWC;

      /// pointer to water relative permeability and specific capacity
      /// relationship kr(P) & C(P) including possible hysteretic effects.
      // UnsatSoilCurves *curves;

      /// Specific soil density 
      double density;

   public:
      /// name of the soil type (i.e. Sands, Clay, etc.)
      String name;

      /// default constructor
      Soil();

      /// 2rd constructor from soil id only
      Soil(unsigned int Id);

      /// 3nd constructor
      Soil(String Name, unsigned int Id,
         const double _Kx, const double _Ky, const double _Kz,
         const double _S0 = 0);

      /// copy constructor
      Soil(const Soil& soil);

      /// destructor
      virtual ~Soil();

      /// Sets horizontal longitudinal hydraulic conductivity
      void SetKx(const double _Kx);

      /// Returns horizontal longitudinal hydraulic conductivity
      double GetKx();

      /// Sets horizontal transverse hydraulic conductivity
      void SetKy(const double _Ky);

      /// Gets horizontal transverse hydraulic conductivity
      double GetKy();

      /// Sets vertical hydraulic conductivity
      void SetKz(const double _Kz);

      /// Gets vertical hydraulic conductivity
      double GetKz();

      /// Sets diagonal hydraulic conductivity tensor at once
      void SetK(const double k[3]);

      /// Gets diagonal hydraulic conductivity tensor at once
      void GetK(double k[3]);

      /// Sets saturated specific storage coefficient
      void SetS0(const double _S0);

      /// Gets saturated specific storage coefficient
      double GetS0();

      /// Sets saturated water content or porosity
      void SetSatWC(const double _wc);

      /// Gets saturated water content
      double GetSatWC();

      /// Sets residual water content
      void SetResidualWC(const double _ResWc);

      /// Gets residual water content
      double GetResidualWC();

      /// Sets soil type id
      void SetId(const unsigned int _id);

      /// Gets soil type id
      unsigned int GetId();

      /// Equality operator
      void operator = (const Soil &another_soil);
   };

}
#endif // HAPP_SOIL

// Inline all basic Get & Set functions ... 
