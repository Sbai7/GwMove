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

#ifndef HAPP_BNDRY_CONDITION
#define HAPP_BNDRY_CONDITION

// Project headers
#include <FileReader.h>
#include <strucmesh.h>
#include <soil.h>
#include <Vector.h> 

// C++ headers
#include <iostream>
#include <vector>
using namespace std;

namespace happ
{/*
 * Types of boundary conditons
 * used in Flow/Transport solvers
 */

   /// No BC defined
#define   NO_BC                   0

   /* BC's for Flow processes */

   /// Fixed head boundary condition code
#define   FIXED_HEAD_BC           1

   /// Fixed flow boundary condition code
#define   FIXED_FLOW_BC           2

   /// Fixed f lux boundary condition code
#define   FIXED_FLUX_BC           3

   /// River boundary condition code
#define   RIVER_BC                4

   /// Leakage boundary condition code
#define   LEAKAGE_BC              5

   /// Infiltration boundary condition code
#define   INFILTRATION_BC         6

   /// Drainage boundary condition code
#define   DRAINAGE_BC             7

   /// Abstraction boundary condition code
#define   ABSTRACTION_BC          8

   /// Seepage boundary condition code
#define   SEEPAGE_BC              9

   /// Recharge boundary condition code
#define   RECHARGE_BC            10

   /// Evapotranspiration boundary condition code
#define   EVAPOTRANSPIRATION_BC  11

   /// Flooding boundary condition code
#define   FLOOD_BC               12

   /* BC's for mass Transport processes */

   /// Fixed concentration boundary condition code
#define   FIXED_CONC_BC          21

   /// Fixed input concentration boundary condition code
#define   FIXED_INPUT_CONC_BC    22 

   /// Fixed mass boundary condition code
#define   FIXED_MASS_BC          23 

   /// Fixed mass flux boundary condition code
#define   FIXED_MASS_FLUX_BC     24

   /// Mixed concentration boundary condition code
#define   MIXED_BC               25

   /// Mixed concentration flux boundary condition code
#define   MIXED_FLUX_BC          26

   /* BC's for heat Transport processes */

   /// Fixed temperature boundary condition code
#define   FIXED_TEMP_BC          31

   /// Fixed input temperature boundary condition code
#define   FIXED_INPUT_TEMP_BC    32

   /// Inactive node condition code 
#define   SINGLE_NODES           99


   /*
   * Types of IJK zones
   */

   /// Generic ijk zone 
#define	 GENERIC_ZONE            1

   /// Boundary condition zone 
#define   BC_ZONE                 2 

   /// Soil type zone 
#define   SOIL_ZONE               3 

   /// SubMesh zone 
#define   SUBMESH_ZONE            4 


   // Error codes related to incorrect boundary conditons


   /*
    * Defines an ordered IJK zone or a submesh of the finite
    * element mesh used for computations.
    * Move this class from here later
    */
   class ijkZone
   {
   protected:
      /// ...
      double zero = 0.;

      /// Zone name 
      String name;

      /// zone type
      unsigned int type;

      /// 
      std::vector<double> bcVar;

   public:
      unsigned int I1;
      unsigned int I2;
      unsigned int J1;
      unsigned int J2;
      unsigned int K1;
      unsigned int K2;
      bool oneNode = false;

      /// Default constructor
      ijkZone() { };

      /// Constructor from one node (0-dim zone)
      ijkZone(const unsigned int i1,
         const unsigned int j1,
         const unsigned int k1) {
         I1 = I2 = i1;
         J1 = J2 = j1;
         K1 = K2 = k1;
         oneNode = true;
      };

      /// Third constructor for a plain zone
      ijkZone(const unsigned int i1,
         const unsigned int j1,
         const unsigned int k1,
         const unsigned int i2,
         const unsigned int j2,
         const unsigned int k2)
      {
         I1 = i1; I2 = i2;
         J1 = j1; J2 = j2;
         K1 = k1; K2 = k2;
      };

      /// Destructor
      virtual ~ijkZone() { };

      /// Copy constructor 
      void operator = (const ijkZone &zone) {
         I1 = zone.I1; I2 = zone.I2;
         J1 = zone.J1; J2 = zone.J2;
         K1 = zone.K1; K2 = zone.K2;
         name = zone.name;
         type = zone.type;
         oneNode = zone.oneNode;
      }

      /// Sets the zone name 
      inline void Name(const String _name) { name = _name; }

      /// Gets the zone name 
      inline String Name() { return name; }

      /// Sets the zone type 
      inline void Type(const unsigned int _type) { type = _type; }

      /// Gets the zone type 
      inline unsigned int & Type() { return type; }

      /// Returns a reference to extracted BC's vector 
      inline std::vector<double>& GetBCVar() {
         return bcVar; 
      }

      /// Returns 0-starting node indices in IJK zone
      void GetNodeIndices(StrucMesh *mesh, std::vector<unsigned int> &index)
      {
         unsigned int i, j, k;
         unsigned int number;
         unsigned int ni = mesh->GetNi();
         unsigned int nj = mesh->GetNj();

         index.clear();
         for (k = K1; k <= K2; k++) {
            for (j = J1; j <= J2; j++) {
               for (i = I1; i <= I2; i++) {
                  number = i + (j - 1)*ni + (k - 1)*ni*nj;
                  index.push_back(number - 1);
               }
            }
         }
      }

      /// Extracts array for boundary condition nodes from global array  
      virtual void ExtractVar(const double *var, StrucMesh *mesh, bool all = true) {
         if (!mesh && !var) return;

         if (!bcVar.empty()) bcVar.clear();

         // ... 
         unsigned int dim;
         if (all) 
            dim = (K2 - K1 + 1)*(J2 - J1 + 1)*(I2 - I1 + 1);
         else 
            dim = (J2 - J1 + 1)*(I2 - I1 + 1);
         bcVar.reserve(dim);

         unsigned int i, j, k;
         unsigned int ni = mesh->GetNi();
         unsigned int nj = mesh->GetNj();
         unsigned int number;
         if (oneNode) {
            // calculate node number from i1,j1,k1 
            if (all)
               number = I1 + (J1 - 1)*ni + (K1 - 1)*ni*nj;
            else 
               number = I1 + (J1 - 1)*ni;
            bcVar.push_back(var[number - 1]);
         }
         else {
            unsigned int start, end;
            if (all) {
               start = K1; end = K2;
            }
            else {
               start = 1; end = 1; 
            }
            for (k = start; k <= end; k++) {
               for (j = J1; j <= J2; j++) {
                  for (i = I1; i <= I2; i++) {
                     number = i + (j - 1)*ni + (k - 1)*ni*nj;
                     bcVar.push_back(var[number - 1]);
                  }
               }
            }
         }
      }

      /// returns total (sum) of extracted variable in this zone 
      virtual double SumVar() const {
         if (bcVar.empty()) return zero;
         double sum = zero;
         unsigned int i, n = bcVar.size();
         for (i = 0; i < n; i++) sum += bcVar[i];
         return sum;
      }

      /// returns sum of positive entries from extracted variable 
      virtual double PositiveVar() const {
         if (bcVar.empty()) return zero;
         double sum = zero;
         unsigned int i, n = bcVar.size();
         for (i = 0; i < n; i++) {
            if (bcVar[i] > zero) sum += bcVar[i];
         }
         return sum;
      }

      /// returns sum of negative entries from extracted variable 
      virtual double NegativeVar() const {
         if (bcVar.empty()) return zero;
         double sum = zero;
         unsigned int i, n = bcVar.size();
         for (i = 0; i < n; i++) {
            if (bcVar[i] < zero) sum += bcVar[i];
         }
         return sum;
      }

      /// Virtual method for nodes dof activation 
      virtual void Activate(unsigned int) {}
      virtual void Activate(unsigned int*) {}
      virtual void Activate(std::vector<unsigned int>&) {}

      /// Virtual method for nodes dof deactivation 
      virtual void Deactivate(unsigned int) {}
      virtual void Deactivate(unsigned int*) {}
      virtual void Deactivate(std::vector<unsigned int>&) {}

      /// Sets soil type of this zone 
      virtual void SoilType(Soil & soil) const {}

      /// Virtual method for flow budget calculation
      virtual void ComputeFlowBudget() const {}

      /// builds a submesh zone with valid FE connectivity 
      virtual StrucMesh *SubMesh(StrucMesh *mesh) { return NULL; }
   };


   /**
    * Class for modelling a zone of a given boundary condition type.
    */
   class BcZone : public ijkZone
   {
   public:
      /// Boundary condition type 
      unsigned int bcType;

      /// Boundary condition 'double' values/parameters 
      std::vector<double> dValues;

      /// Boundary conditions 'int' parameters 
      std::vector<unsigned int> iValues;

      /// Default constructor 
      BcZone() { type = BC_ZONE; }

      /// 2nd constructor 
      BcZone(ijkZone &zone,
         const unsigned int _bcType,
         const std::vector<double> _dValues)
         : ijkZone(zone.I1, zone.J1, zone.K1,
         zone.I2, zone.J2, zone.K2) {
         bcType = _bcType;
         dValues = _dValues;
         type = BC_ZONE;
         name = zone.Name();
      }

      /// 3rd constructor 
      BcZone(ijkZone &zone,
         const unsigned int _bcType,
         const std::vector<double> _dValues,
         const std::vector<unsigned int> _iValues)
         : ijkZone(zone.I1, zone.J1, zone.K1,
         zone.I2, zone.J2, zone.K2) {
         bcType = _bcType;
         dValues = _dValues;
         iValues = _iValues;
         type = BC_ZONE;
         name = zone.Name();
      }

      /// class destructor 
      virtual ~BcZone() {
         dValues.clear();
         iValues.clear();
         bcVar.clear();
      }

      /// Copy constructor 
      void operator = (const BcZone &zone) {
         (ijkZone)*this = (ijkZone)zone;
         bcType = zone.bcType;
         dValues = zone.dValues;
         iValues = zone.iValues;
         name = zone.name;
      }
   };

   /* Other classes to implement:
    *
    * class SuperZone      : public std::vector<ijkZone>   // union of many ijk zones
    * class SoilZone       : public ijkZone
    * class RefinementZone : public SuperZone
    *
    * Other possible specializations:
    *
    * class DirichletZone : public
    * class WellZone      : public SuperZone
    * class RechargeZone  : public SuperZone
    * class FluxZone      : public SuperZone
    * class LeakageZone   : public SuperZone
    * class RiverZone     : public SuperZone
    * class DrainageZone  : public SuperZone
    * class SeepageZone   : public SuperZone
    * class SgdZone       : public SuperZone   // Sea groundwater discharge zone
    * class EvtZone       : public SuperZone   // Evapotranspiration zone
    */


   /**
    * Class representing flow boundary conditions.
    */
   class BndCondition
   {
   private:
      /// number of 'packed' IJK-type boundary conditions
      unsigned int nBC;

      /// dimension of boundary conditions tables
      unsigned int dimension;

      /// activation flags for BC's 
      unsigned int *active;


   public:
      /// pointer to finite element mesh
      StrucMesh *mesh;

      /// Array of ijk 'packed' boundary conditions objects
      std::vector<BcZone*> bcZone;

      /// type of nodal boundary conditions
      unsigned int *bcType;

      /// parameters of boundary conditions
      double *bcValue;

      /// logical flag indicating a fixed head boundary condition type
      bool *fixedHead;

      /// constructor
      BndCondition(StrucMesh *pMesh);

      /// destructor
      virtual ~BndCondition();

      /// Returns boundary condition type at node number n
      inline unsigned int GetType(const unsigned int n) {
         if (n < dimension) return bcType[n];
         else return 0;
      }

      /// Sets boundary condition type for node n
      inline void SetType(const unsigned int n,
         const unsigned int type) {
         if (n < dimension && type <= SINGLE_NODES) bcType[n] = type;
      }

      /// returns boundary condition value
      inline double GetValue(const unsigned int n,
         const unsigned int i) {
         return bcValue[4 * n + i - 1];
      }

      /// sets boundary condition value
      inline void SetValue(const unsigned int n,
         const unsigned int i,
         const double value) {
         bcValue[4 * n + i - 1] = value;
      }

      /// Is BC @ dof i active ? 
      inline bool isActive(unsigned int i) {
         return active[i];
      }

      /// deactivates a given boundary condition @ dof i 
      inline void Deactivate(unsigned int i) {
         if (active[i]) active[i] = false;
      }

      /// activates a given boundary condition @ dof i 
      inline void Activate(unsigned int i) {
         if (!active[i]) active[i] = true;
      }

      /// returns the number of node-free boundary condition 
      inline unsigned int GetNoBC() const {
         unsigned int i, nNoBC = 0;
         unsigned int nn = mesh->GetNnodes();
         for (i = 0; i < nn; i++) if (bcType[i] == NO_BC) nNoBC++;

         return nNoBC; 
      }

      /// returns the number of elemental boundary conditions 
      inline unsigned int GetnBC() const {
         unsigned int i, nNodeBC = 0;
         unsigned int nn = mesh->GetNnodes();
         for (i = 0; i < nn; i++)
            if (bcType[i] > NO_BC && bcType[i] <= FLOOD_BC) nNodeBC++;

         return nNodeBC;
      }

      /// returns the number of single/inactive nodes
      inline unsigned int GetnSingleNodes() const {
         unsigned int i, nDeactivated = 0;
         unsigned int nn = mesh->GetNnodes();
         for (i = 0; i < nn; i++) if (bcType[i] == SINGLE_NODES) nDeactivated++;

         return nDeactivated;
      }

      /// extracts the elements of a vector at imposed nodes 
      bool ExtractAtBC(Vector& v) const;

      /// reads a boundary conditions file
      bool Read(String projectFile);

      /// Adds a boundary condition from ijk zone
      void Add(const BcZone zone);
      void Add(const ijkZone zone, unsigned int type, double val[4]);

      /// overloading of output operator
      friend ostream& operator <<(ostream& os, BndCondition& bc);
   };
}

#endif // HAPP_BNDRY_CONDITION
