// ----------------------------------------------------------------------------
 // nexus | SiPMVUVnaked.cc
 //
 // 3x3 mm2 FBK SiPM geometry.
 //
 // The NEXT Collaboration
 // ----------------------------------------------------------------------------

 #ifndef SIPM_FBK_VUV_H
 #define SIPM_FBK_VUV_H

 #include "BaseGeometry.h"
 #include <G4ThreeVector.hh>

 class G4GenericMessenger;

 namespace nexus {

   class SiPMVUVnaked: public BaseGeometry
   {
   public:
     /// Constructor
     SiPMVUVnaked();
     /// Destructor
     ~SiPMVUVnaked();

     /// Invoke this method to build the volumes of the geometry
     void Construct();

     void SetVisibility(G4bool vis);
     void SetPDE(G4double eff);
     void SetTimeBinning (G4double time_binning);
     void SetSensorDepth (G4int sensor_depth);
     void SetMotherDepth (G4int mother_depth);
     void SetNamingOrder (G4int naming_order);


   private:
     G4bool visibility_;

     // PDE for the sensor
     G4double eff_;

     G4double time_binning_;
     G4int    sensor_depth_;
     G4int    mother_depth_;
     G4int    naming_order_;

   };

   inline void SiPMVUVnaked::SetVisibility(G4bool vis)
   { visibility_ = vis; }

   inline void SiPMVUVnaked::SetPDE(G4double eff)
   { eff_ = eff; }

   inline void SiPMVUVnaked::SetTimeBinning(G4double time_binning)
   { time_binning_ = time_binning; }

   inline void SiPMVUVnaked::SetSensorDepth(G4int sensor_depth)
   { sensor_depth_ = sensor_depth; }

   inline void SiPMVUVnaked::SetMotherDepth(G4int mother_depth)
   { mother_depth_ = mother_depth; }

   inline void SiPMVUVnaked::SetNamingOrder(G4int naming_order)
   { naming_order_ = naming_order; }


 } // end namespace nexus

 #endif
