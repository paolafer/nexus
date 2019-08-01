// ----------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------

#include "NeutrinoEventAction.h"
#include "Trajectory.h"
#include "PersistencyManager.h"
#include "IonizationHit.h"

#include <G4Event.hh>
#include <G4VVisManager.hh>
#include <G4Trajectory.hh>
#include <G4GenericMessenger.hh>
#include <G4HCofThisEvent.hh>
#include <G4SDManager.hh>
#include <G4HCtable.hh>
#include <globals.hh>


namespace nexus {


  NeutrinoEventAction::NeutrinoEventAction():
    G4UserEventAction(), _nevt(0), _nupdate(10)
  {
    // _msg = new G4GenericMessenger(this, "/Actions/NeutrinoEventAction/");
  }



  NeutrinoEventAction::~NeutrinoEventAction()
  {
  }



  void NeutrinoEventAction::BeginOfEventAction(const G4Event* /*event*/)
  {
    // Print out event number info
    if ((_nevt % _nupdate) == 0) {
      G4cout << " >> Event no. " << _nevt  << G4endl;
      if (_nevt  == (10 * _nupdate)) _nupdate *= 10;
    }
  }



  void NeutrinoEventAction::EndOfEventAction(const G4Event* event)
  {
    _nevt++;

    G4bool neutrino = false;

    // Get the trajectories stored for this event and loop through them
    // to find if there is any neutrino.

      G4TrajectoryContainer* tc = event->GetTrajectoryContainer();
      if (tc) {
        for (unsigned int i=0; i<tc->size(); ++i) {
          Trajectory* trj = dynamic_cast<Trajectory*>((*tc)[i]);
          if (trj->GetParticleName() == "nu_e" || trj->GetParticleName() == "nu_mu" ||
	      trj->GetParticleName() == "nu_tau" || trj->GetParticleName() == "anti_nu_e" ||
	      trj->GetParticleName() == "anti_nu_mu" || trj->GetParticleName() == "anti_nu_tau") {
	    neutrino = true;
	  }
          // Draw tracks in visual mode
          if (G4VVisManager::GetConcreteInstance()) trj->DrawTrajectory();
        }
      }

      PersistencyManager* pm = dynamic_cast<PersistencyManager*>
        (G4VPersistencyManager::GetPersistencyManager());

      if (!event->IsAborted() && neutrino) {
	pm->StoreCurrentEvent(true);
      } else {
	pm->InteractingEvent(false);
      }
    
  }


} // end namespace nexus
