// ----------------------------------------------------------------------------
//  Author:  <paola.ferrario@dipc.org>
//  Created: 5 March 2020
// ----------------------------------------------------------------------------

#ifndef SILICON_PM_pet_H
#define SILICON_PM_pet_H

#include "BaseGeometry.h"
#include <G4ThreeVector.hh>

class G4GenericMessenger;

namespace nexus {


  /// Geometry of a generic SiPM

  class SiPMpet: public BaseGeometry
  {
  public:
    /// Constructor
    SiPMpet();
    /// Destructor
    ~SiPMpet();

    /// Invoke this method to build the volumes of the geometry
    void Construct();

  private:

    // Visibility of the tracking plane
    G4bool visibility_;

    // Optical properties to be used for epoxy
    G4double refr_index_;

    // PDE of SiPMs
    G4double eff_;

    G4double time_binning_;
    G4double sipm_size_, sipm_z_, offset_;
    G4double active_depth_, active_size_;

    G4int mother_depth_, naming_order_;

    G4bool wls_coating_;
    G4double decay_time_; ///< decay time of WLS
    G4double qe_; ///< quantum efficiency of WLS

    // Messenger for the definition of control commands
    G4GenericMessenger* msg_;

  };


} // end namespace nexus

#endif
