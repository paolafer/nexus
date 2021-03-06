// ----------------------------------------------------------------------------
// nexus | IonizationElectron.cc
//
// Definition of the ionization electron.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------


#include "IonizationElectron.h"

#include <G4ParticleTable.hh>

#include "CLHEP/Units/PhysicalConstants.h"


namespace nexus {

  using namespace CLHEP;


  IonizationElectron* IonizationElectron::instance_ = 0;



  IonizationElectron::IonizationElectron(): G4ParticleDefinition()
  {
  }



  IonizationElectron::~IonizationElectron()
  {
  }



  IonizationElectron* IonizationElectron::Definition()
  {
    // Return the (only) instance if it has been defined already.
    if (instance_) return instance_;

    // Search for ionization electron definition in particle table.
    // If we get a null pointer, the particle has not been created yet.

    const G4String name = "ie-";

    G4ParticleDefinition* pdef =
      G4ParticleTable::GetParticleTable()->FindParticle(name);

    if (!pdef) {
      // Create the particle. Arguments for constructor are as follows:
      //    name         mass            width           charge
      //    2*spin       parity          C-conjugation
      //    2*Isospin    2*Isospin3      G-parity
      //    type         lepton number   baryon number   PDG encoding
      //    stable       lifetime        decay table
      //    shortlived   subType         anti_encoding
      pdef =
        new G4ParticleDefinition(name, electron_mass_c2, 0., -1.*eplus,
				  1, 0, 0,
				  0, 0, 0,
				  "lepton", 1, 0, 11,
				  true, -1.0, NULL,
				  false, "e", 0);
    }

    instance_ = (IonizationElectron*) pdef;
    return instance_;
  }


} // end namespace nexus
