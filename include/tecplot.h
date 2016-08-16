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

#ifndef HAPP_TECPLOT
#define HAPP_TECPLOT

// Project headers
#include <strucmesh.h>
#include <AdString.h>

// C++ headers
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;


namespace happ
{
   // Types definitions

   enum ZONETYPE {
      ORDERED,
      FETRIANGLE, 
      FEQUADRILATERAL, 
      FETETRAHEDRON, 
      FEBRICK, 
      FEPOLYHEDRON
   };

   enum DATAPACKING {
      POINT,
      BLOCK
   };

   /**
    *  Writes model output results for post-processing with Tecplot.
    */
   class TecplotIO
   {
   private:
      /// pointer to mesh object 
      StrucMesh *mesh;

      /// model scalar or vector variables 
      std::vector<double*> variables;

      /// is the file to write to in binary format ?
      bool m_binary;

      /// index of current zone in this tecplot file instance
      unsigned int m_currentzone;

      /// number of active variables
      unsigned int m_nvariables;

      /// file name, including the path, of the tecplot file
      String m_filename;

      /// path of the tecplot file
      String m_path;

      /// tecplot file title
      String m_title;

      /// zone title 
      String m_zonetitle;

      /// variables string
      String m_svariables;

      /// data packing: either POINT or BLOCK
      String m_datapacking;

      /// zone type 
      unsigned int m_zonetype;

      /// is the last zone mesh (i.e. X,Y,Z coordinates) shared ?
      bool m_meshShared; 

      /// is the last zone Z-coordinate shared ? 
      bool m_zCoordShared;

      /// is the output file opened?
      bool m_fileOpened;

      /// output file buffer 
      std::filebuf outfile;


      /// Parse Tecplot input stream title line 
      bool ParseTitle();

      /// Parse Tecplot input stream variables line
      bool ParseVariables();

      /// Parse Tecplot input stream zone (first one for instance) header
      bool ParseZoneHeader();

      /// Parse Tecplot input stream zone data 
      bool ParseZoneData();

   public:
      /// Default constructor
      TecplotIO();

      /// Constructor from a mesh pointer
      TecplotIO(StrucMesh *pmesh);

      /// Destructor
      virtual ~TecplotIO();

      /// Sets the file name
      inline void SetFilename(const String filename) {
         m_filename = filename;
      }

      /// Sets the file title
      inline void SetTitle(const String title) {
         m_title = title;
      }

      /// Sets zone title 
      inline void SetZoneTitle(const String title) {
         m_zonetitle = title;
      }

      /// Sets the path
      inline void SetPath(const String path) {
         m_path = path;
      }

      /// Sets the variables string
      inline void SetVariables(const String var) {
         m_svariables = var;
      }

      /// Sets data packing format
      inline void SetDataPacking(const DATAPACKING datapacking) {
         m_datapacking = (datapacking == DATAPACKING::POINT ? "POINT" : "BLOCK");
      }

      /// Sets zone type 
      inline void SetZoneType(const ZONETYPE zonetype) {
         m_zonetype = zonetype;
      }

      /// Sets number of variables
      inline void SetnVariables(const unsigned int nvar) {
         m_nvariables = nvar;
      }

      /// Sets number of current zone 
      inline void SetZone(const unsigned int zone) {
         m_currentzone = zone;
      }

      /// Sets file format
      inline void SetFileFormat(const bool binary) {
         m_binary = binary;
      }

      /// Allow mesh sharing from last zone 
      inline void ShareMesh(bool meshshared) {
         m_meshShared = meshshared;
         if (m_meshShared) m_zCoordShared = false;
      }

      /// Allow Z-coordinate sharing from last mesh 
      inline void ShareZCoord(bool zshared) {
         m_zCoordShared = zshared;
         if (m_zCoordShared) m_meshShared = false;
      }

      /// Returns file name
      inline String GetFilename() const {
         return m_filename;
      }

      /// Returns the file path
      inline String GetPath() const {
         return m_path;
      }

      /// Returns the file header 
      inline String GetTitle() const {
         return m_title;
      }

      /// Returns the current zone title
      inline String GetZoneTitle() const {
         return m_zonetitle;
      }

      /// Returns the variables string
      inline String GetVariables() const {
         return m_svariables;
      }

      /// Returns data packing 
      inline String GetDataPacking() const {
         return m_datapacking;
      }

      /// Returns zone type 
      inline unsigned int GetZoneType() const {
         return m_zonetype;
      }

       inline String GetZoneTypeStr() const {
          if (m_zonetype == ZONETYPE::ORDERED) return "ORDERED";
          if (m_zonetype == ZONETYPE::FETRIANGLE) return "FETRIANGLE";
          if (m_zonetype == ZONETYPE::FEQUADRILATERAL) return "FEQUADRILATERAL";
          if (m_zonetype == ZONETYPE::FETETRAHEDRON) return "FETETRAHEDRON";
          if (m_zonetype == ZONETYPE::FEBRICK) return "FEBRICK";
          if (m_zonetype == ZONETYPE::FEPOLYHEDRON) return "FEPOLYHEDRON";
          else
            return "FEBRICK";
      }

     /// Returns number of variables
      inline unsigned int GetnVariables() const {
         return m_nvariables;
      }

      /// Returns number of current tecplot zone
      inline unsigned int GetZone() const {
         return m_currentzone;
      }

      /// Returns true if the mesh is being shared between zones 
      inline bool IsMeshShared() const {
         return m_meshShared; 
      }

      /// Returns true if Z-ccordinate of the mesh is being shared 
      inline bool IsZCoordShared() const {
         return m_zCoordShared;
      }
          
      /// Returns the file format
      inline bool GetFileFormat() const {
         return m_binary;
      }

      /// Adds a dataset to the file
      bool AddDataset(double *var, const String name);

      /// Updates dataset of index i 
      void UpdateDataset(double *var, const unsigned int i);

      /// Removes a dataset having index i from the file
      bool RemoveDataset(const unsigned int i);

      /// Initializes a new Tecplot zone with indices of shared variables
      void InitNextZone(String title);

      /// Reads variables from a tecplot file
      bool Read();

      /// Writes current zone to opened tecplot file 
      bool Write();
      bool Write(String filename);

      /// Writes current zone to old Water file format 
      bool WriteWaterFile(String filename);

      /// Closes the ouput file 
      bool Dispose();

      /// Overloading the output stream operator
      friend ostream& operator << (std::ostream& os, TecplotIO& tec);

      /// Overloading the input stream operator
      friend istream& operator >> (std::istream& is, TecplotIO& tec);
   };
}

#endif // HAPP_TECPLOT
