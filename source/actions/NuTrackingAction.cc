#include "NuTrackingAction.h"

#include "Trajectory.h"
#include "TrajectoryMap.h"
#include "IonizationElectron.h"

#include <G4Track.hh>
#include <G4TrackingManager.hh>
#include <G4Trajectory.hh>
#include <G4ParticleDefinition.hh>
#include <G4NeutrinoE.hh>
#include <G4NeutrinoMu.hh>
#include <G4NeutrinoTau.hh>
#include <G4AntiNeutrinoE.hh>
#include <G4AntiNeutrinoMu.hh>
#include <G4AntiNeutrinoTau.hh>



using namespace nexus;



NuTrackingAction::NuTrackingAction(): G4UserTrackingAction()
{
}



NuTrackingAction::~NuTrackingAction()
{
}



void NuTrackingAction::PreUserTrackingAction(const G4Track* track)
{
  // Do nothing if the track is not a neutrino
  if (track->GetDefinition() != G4NeutrinoE::Definition()       & 
      track->GetDefinition() != G4NeutrinoMu::Definition()      &
      track->GetDefinition() != G4NeutrinoTau::Definition()     &
      track->GetDefinition() != G4AntiNeutrinoE::Definition()   & 
      track->GetDefinition() != G4AntiNeutrinoMu::Definition()  &
      track->GetDefinition() != G4AntiNeutrinoTau::Definition()) {
      fpTrackingManager->SetStoreTrajectory(false);    
      return;
  }

  // Create a new trajectory associated to the track.
  // N.B. If the processesing of a track is interrupted to be resumed
  // later on (to process, for instance, its secondaries) more than
  // one trajectory associated to the track will be created, but 
  // the event manager will merge them at some point.
  G4VTrajectory* trj = new Trajectory(track);
  
   // Set the trajectory in the tracking manager
  fpTrackingManager->SetStoreTrajectory(true);
  fpTrackingManager->SetTrajectory(trj);
 }



void NuTrackingAction::PostUserTrackingAction(const G4Track* track)
{
  // Do nothing if the track is an optical photon or an ionization electron
  if (track->GetDefinition() != G4NeutrinoE::Definition()       & 
      track->GetDefinition() != G4NeutrinoMu::Definition()      &
      track->GetDefinition() != G4NeutrinoTau::Definition()     &
      track->GetDefinition() != G4AntiNeutrinoE::Definition()   & 
      track->GetDefinition() != G4AntiNeutrinoMu::Definition()  &
      track->GetDefinition() != G4AntiNeutrinoTau::Definition()) return;
 

  Trajectory* trj = (Trajectory*) TrajectoryMap::Get(track->GetTrackID());

  // Do nothing if the track has no associated trajectory in the map
  if (!trj) return;

  // Record final time and position of the track
  trj->SetFinalPosition(track->GetPosition());
  trj->SetFinalTime(track->GetGlobalTime());
  trj->SetTrackLength(track->GetTrackLength());
  trj->SetDecayVolume(track->GetVolume()->GetName());
}
