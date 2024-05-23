// ----------------------------------------------------------------------------
// nexus | AnalysisSteppingAction.cc
//
// This class allows the user to print the total number of photons detected by
// all kinds of photosensors at the end of the run.
// It also shows examples of information that can be accessed at the stepping
// level, so it is useful for debugging.
//
// The  NEXT Collaboration
// ----------------------------------------------------------------------------

#include "AnalysisSteppingAction.h"
#include "FactoryBase.h"
#include "IonizationElectron.h"

#include <G4Step.hh>
#include <G4SteppingManager.hh>
#include <G4ProcessManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4VPhysicalVolume.hh>

using namespace nexus;

REGISTER_CLASS(AnalysisSteppingAction, G4UserSteppingAction)

AnalysisSteppingAction::AnalysisSteppingAction(): G4UserSteppingAction(),
  ph_on_grids_(0), photoel_ie_(0), refl_(0), abs_(0), refr_(0)
{
}



AnalysisSteppingAction::~AnalysisSteppingAction()
{
  /*
  G4double total_counts = 0;
  detectorCounts::iterator it = my_counts_.begin();
  while (it != my_counts_.end()) {
    G4cout << "Detector " << it->first << ": " << it->second << " counts" << G4endl;
    total_counts += it->second;
    it ++;
  }
  G4cout << "TOTAL COUNTS: " << total_counts << G4endl;
  */
  G4cout << "Photons on SS mesh: " << ph_on_grids_ << G4endl;
  G4cout << "Ionization e- created: " << photoel_ie_ << G4endl;
  G4cout << "Reflected photons: " << refl_ << G4endl;
  G4cout << "% reflected: " << G4double(refl_) / ph_on_grids_ << G4endl;
  G4cout << "% refracted: " << G4double(refr_) / ph_on_grids_ << G4endl;
  G4cout << "% absorbed: " << G4double(abs_) / ph_on_grids_ << G4endl;
  G4cout << "% ie-: " << G4double(photoel_ie_) / (ph_on_grids_ - refl_) << G4endl;
  
}



void AnalysisSteppingAction::UserSteppingAction(const G4Step* step)
{

  G4ParticleDefinition* pdef = step->GetTrack()->GetDefinition();

  //Check whether the track is an optical photon
  if (pdef != G4OpticalPhoton::Definition() && pdef != IonizationElectron::Definition()) return;

  /*
  // example of information one can access about optical photons

  G4Track* track = step->GetTrack();
  G4int pid = track->GetParentID();
  G4int tid = track->GetTrackID();
  G4StepPoint* point1 = step->GetPreStepPoint();
  G4StepPoint* point2 = step->GetPostStepPoint();
  G4TouchableHandle touch1 = point1->GetTouchableHandle();
  G4TouchableHandle touch2 = point2->GetTouchableHandle();
  G4String vol1name = touch1->GetVolume()->GetName();
  G4String vol2name = touch2->GetVolume()->GetName();

  G4int copy_no = step->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);
  */

  
  if (step->GetPostStepPoint()) {
    G4String volname = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
    G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    
    //    if ((volname == "EL_GRID_GATE" || volname == "EL_GRID_ANODE") && pdef == G4OpticalPhoton::Definition() &&  step->GetTrack()->GetKineticEnergy() > 4.3*CLHEP::eV) {
    if (volname == "STEEL" && pdef == G4OpticalPhoton::Definition() &&  step->GetTrack()->GetKineticEnergy() > 4.3*CLHEP::eV && (proc == "Transportation")) {
      ph_on_grids_++;
      G4cout << "Cuenta!" << G4endl;
    }

    if (step->GetTrack()->GetCreatorProcess() ) {
      G4String proc_name = step->GetTrack()->GetCreatorProcess()->GetProcessName();
      if (pdef == IonizationElectron::Definition() && proc_name == "OpPhotoelectricEffect") {
        photoel_ie_++;
      }
    }

    G4String proc_name = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    if ((proc_name == "OpAbsorption") && volname == "STEEL") {
      abs_ += 1;
    }
  }
  
  // Retrieve the pointer to the optical boundary process.
  // We do this only once per run defining our local pointer as static.
  static G4OpBoundaryProcess* boundary = 0;

  if (!boundary) { // the pointer is not defined yet
    // Get the list of processes defined for the optical photon
    // and loop through it to find the optical boundary process.
    G4ProcessVector* pv = pdef->GetProcessManager()->GetProcessList();
    for (size_t i=0; i<pv->size(); i++) {
      if ((*pv)[i]->GetProcessName() == "OpBoundary") {
	boundary = (G4OpBoundaryProcess*) (*pv)[i];
	break;
      }
    }
  }
  
  if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    if ((boundary->GetStatus() == FresnelReflection) ||
        (boundary->GetStatus() == SpikeReflection) ||
        (boundary->GetStatus() == LambertianReflection) ||
        (boundary->GetStatus() == TotalInternalReflection) ||
        (boundary->GetStatus() == LobeReflection) ||
        (boundary->GetStatus() == BackScattering)) {
      refl_ += 1;
    } else if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "STEEL") {
      G4cout << "Enters steel with " << boundary->GetStatus() << G4endl;
      // }
      refr_ += 1;
    } else {
      // G4cout << boundary->GetStatus() << G4endl;
      }

    /*
    if (boundary->GetStatus() == Detection ){
      G4String detector_name = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
      //G4cout << "##### Sensitive Volume: " << detector_name << G4endl;

      detectorCounts::iterator it = my_counts_.find(detector_name);
      if (it != my_counts_.end()) my_counts_[it->first] += 1;
      else my_counts_[detector_name] = 1;
    }
    */
  }
  
  return;
}
