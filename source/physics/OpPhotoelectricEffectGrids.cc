// ----------------------------------------------------------------------------
// nexus | OpPhotoelectricEffectGrids.h
//
// Add photoelectric effect physics to optical photons hitting
// specific materials which are not dielectric.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "OpPhotoelectricEffectGrids.h"

#include "IonizationElectron.h"
#include "OpticalMaterialProperties.h"

#include <G4Material.hh>
#include <G4ParticleDefinition.hh>
#include <G4OpticalPhoton.hh>
#include <G4RandomDirection.hh>
#include <G4Poisson.hh>

namespace nexus {


  using namespace CLHEP;

  OpPhotoelectricEffectGrids::OpPhotoelectricEffectGrids(const G4String& process_name,
                                             G4ProcessType type):
    G4VRestProcess(process_name, type), particle_change_(0)
  {
    // Create particle change object
    particle_change_ = new G4ParticleChange();
    pParticleChange = particle_change_;
  }



  OpPhotoelectricEffectGrids::~OpPhotoelectricEffectGrids()
  {
    delete particle_change_;
  }

  /*
  void OpPhotoelectricEffectGrids::InitialiseProcess(const G4ParticleDefinition*)
  {
    if(!is_initialised_) {
     is_initialised_ = true;
     if(nullptr == EmModel(0)) { SetEmModel(new G4EmStandardPhysics_option4()); }
     EmModel(0)->SetLowEnergyLimit(MinKinEnergy());
     EmModel(0)->SetHighEnergyLimit(MaxKinEnergy());
     AddEmModel(1, EmModel(0));
    }
  }

  */

  G4bool OpPhotoelectricEffectGrids::IsApplicable(const G4ParticleDefinition& pdef)
  {
    return pdef == *G4OpticalPhoton::Definition();
  }

  G4VParticleChange*
  OpPhotoelectricEffectGrids::AtRestDoIt(const G4Track& track, const G4Step& step)
  {
    G4cout << "here" << G4endl;
    // Initialize particle change with current track values
    particle_change_->Initialize(track);

    if (step.GetPostStepPoint()) {
      G4String volname = step.GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
      if ((volname == "EL_GRID_GATE" || volname == "EL_GRID_ANODE")) {

        
      /*
    const G4Material* material = track.GetMaterial();

    G4MaterialPropertiesTable* mpt = material->GetMaterialPropertiesTable();

    if (!mpt ||
        !mpt->ConstPropertyExists("WORK_FUNCTION") ||
        !mpt->ConstPropertyExists("OP_PHOTOELECTRIC_PROBABILITY"))
      particle_change_->ProposeTrackStatus(fStopAndKill);
      //return G4VRestProcess::PostStepDoIt(track, step);
      */   
     //G4cout << material->GetName() << G4endl;

        G4double photon_energy = track.GetDynamicParticle()->GetTotalEnergy();
        
        G4MaterialPropertiesTable* mpt = opticalprops::Steel(1.e-3);
        G4double work_function = mpt->GetConstProperty("WORK_FUNCTION");
        G4double probability   = mpt->GetConstProperty("OP_PHOTOELECTRIC_PROBABILITY");
        
        if (!work_function || !probability)
          particle_change_->ProposeTrackStatus(fStopAndKill);
        //return G4VRestProcess::PostStepDoIt(track, step);
        
        // We have to compare the energy with the work function here because
        // Geant4 doesn't deal with the vector of probabilities correctly.
        if ((photon_energy   <  work_function) ||
            (G4UniformRand() >= probability  ) ) {
          
          particle_change_->ProposeTrackStatus(fStopAndKill);
          //return G4VRestProcess::PostStepDoIt(track, step);
        }
        
        particle_change_->ProposeTrackStatus(fStopAndKill);
        
        // We force each photon to produce a maximum of one electron.
        particle_change_->SetNumberOfSecondaries(1);
        
        G4ThreeVector momentum_direction = G4RandomDirection();
        G4double      kinetic_energy     = photon_energy - work_function;
        
        G4DynamicParticle*
          ie = new G4DynamicParticle(IonizationElectron::Definition(),
                                     momentum_direction, kinetic_energy);
        
        // Setting the emission point to either the pre or the post
        // step points results in a warning about the new particle
        // having a negative local time. This seems to be caused by
        // the positions being in a boundary between two volumes.
        // Taking the midpoint avoids the warning without affecting
        // the physics.
        G4StepPoint*  pre  = step.GetPreStepPoint();
        G4StepPoint*  post = step.GetPostStepPoint();
        G4double      time = post->GetGlobalTime();
        G4ThreeVector pos  = (pre->GetPosition() + post->GetPosition()) / 2;
        
        G4Track* new_track = new G4Track(ie, time, pos);
        new_track->SetTouchableHandle(post->GetTouchableHandle());
        new_track->SetParentID(track.GetTrackID());
        
        particle_change_->AddSecondary(new_track);
        
      }
    }
    
    return particle_change_;
  }


  G4double OpPhotoelectricEffectGrids::AtRestGetPhysicalInteractionLength(const G4Track&, G4ForceCondition* condition)
  {
    *condition = StronglyForced;
    return DBL_MAX;
  }
  /*
  G4double OpPhotoelectricEffectGrids::GetMeanFreePath(const G4Track&, G4double,
                                                 G4ForceCondition* condition)
  {
    *condition = StronglyForced;
    return DBL_MAX;
  }
  */ 
  
  G4double OpPhotoelectricEffectGrids::GetMeanLifeTime(const G4Track&,
                                                 G4ForceCondition* condition)
  {
    *condition = Forced;
    return DBL_MAX;
  }
   
  
} // end namespace nexus
