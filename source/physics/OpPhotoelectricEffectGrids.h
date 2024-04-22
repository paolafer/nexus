// ----------------------------------------------------------------------------
// nexus | OpPhotoelectricEffectGrids.h
//
// Add photoelectric effect physics to optical photons hitting
// specific materials which are not dielectric.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef OPTICAL_PHOTOELECTRIC_GRIDS_H
#define OPTICAL_PHOTOELECTRIC_GRIDS_H


#include <G4VRestProcess.hh>

class G4Material;

namespace nexus {

  class OpPhotoelectricEffectGrids: public G4VRestProcess
  {
  public:

    /// Constructor
    OpPhotoelectricEffectGrids(const G4String& process_name="OpPhotoelectricEffectGrids",
			                   G4ProcessType type = fUserDefined);
    /// Destructor
    ~OpPhotoelectricEffectGrids();

    /// Whether a particle is sensitive to this physics.
    /// Only optical photons apply
    G4bool IsApplicable(const G4ParticleDefinition&) override;

    /// Calculates ionization electron emission probability
    /// and generates new particles if necessary.
    G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) override;

    //  protected:
    //    void InitialiseProcess(const G4ParticleDefinition*) override;

  private:
    G4double AtRestGetPhysicalInteractionLength(const G4Track&, G4ForceCondition*) override;

    /*
    /// Returns infinity; i. e. the process does not limit the step,
    /// but sets the 'StronglyForced' condition for the PostStepDoIt
    /// to be invoked at every step
    G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);
    */
    /// Returns infinity; i. e. the process does not limit the time,
    /// but sets the 'StronglyForced' condition for the AtRestDoIt
    /// to be invoked at every step
    G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*) override;

    G4bool ValidMaterial(const G4Material*);

  private:
    G4ParticleChange* particle_change_;

    G4bool is_initialised_;

  };

} // end namespace nexus

#endif
