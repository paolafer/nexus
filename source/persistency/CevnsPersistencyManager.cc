#include "CevnsPersistencyManager.h"

#include "Trajectory.h"
#include "TrajectoryMap.h"
#include "IonizationSD.h"
#include "PmtSD.h"
#include "NexusApp.h"
#include "DetectorConstruction.h"
#include "BaseGeometry.h"
#include "HDF5Writer.h"

#include <G4GenericMessenger.hh>
#include <G4Event.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4SDManager.hh>
#include <G4HCtable.hh>
#include <G4RunManager.hh>
#include <G4Run.hh>

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "CLHEP/Units/SystemOfUnits.h"

using namespace nexus;

CevnsPersistencyManager::CevnsPersistencyManager(G4String historyFile_init, G4String historyFile_conf):
  G4VPersistencyManager(), _msg(0),
  _ready(false), _store_evt(true),  event_type_("other"),
  _saved_evts(0),  _nevt(0), _start_id(0), _first_evt(true),
  _h5writer(0)
{

  _historyFile_init = historyFile_init;
  _historyFile_conf = historyFile_conf;

  _msg = new G4GenericMessenger(this, "/nexus/persistency/");
  _msg->DeclareMethod("outputFile", &CevnsPersistencyManager::OpenFile, "");
  _msg->DeclareProperty("eventType", event_type_, "Type of event: bb0nu, bb2nu or background.");
  _msg->DeclareProperty("start_id", _start_id, "Starting event ID for this job.");
}



CevnsPersistencyManager::~CevnsPersistencyManager()
{
  delete _msg;
  delete _h5writer;
}



void CevnsPersistencyManager::Initialize(G4String historyFile_init, G4String historyFile_conf)
{

  // Get a pointer to the current singleton instance of the persistency
  // manager using the method of the base class
  CevnsPersistencyManager* current = dynamic_cast<CevnsPersistencyManager*>
    (G4VPersistencyManager::GetPersistencyManager());

  // If no instance exists yet, create a new one.
  // (Notice that the above dynamic cast would also return 0 if an instance
  // of another G4VCevnsPersistencyManager-derived was previously set, resulting
  // in the leak of that object since the pointer will no longer be
  // accessible.)
  if (!current) current = new CevnsPersistencyManager(historyFile_init, historyFile_conf);
}


void CevnsPersistencyManager::OpenFile(G4String filename)
{
  _h5writer = new HDF5Writer();
  G4String hdf5file = filename + ".h5";
  _h5writer->Open(hdf5file);
  return;
}



void CevnsPersistencyManager::CloseFile()
{
  _h5writer->Close();
  G4cout << "Closing file" << G4endl;
  return;
}



G4bool CevnsPersistencyManager::Store(const G4Event* event)
{
  if (!_store_evt) {
    TrajectoryMap::Clear();
    return false;
  }

  _saved_evts++;


  if (_first_evt) {
    _first_evt = false;
    _nevt = _start_id;
  }


  // Store the trajectories of the event as Gate particles
  StoreTrajectories(event->GetTrajectoryContainer());

  _nevt++;

  TrajectoryMap::Clear();
  StoreCurrentEvent(true);

  return true;
}



void CevnsPersistencyManager::StoreTrajectories(G4TrajectoryContainer* tc)         
{
  // If the pointer is null, no trajectories were stored in this event
  if (!tc) return;

  // Loop through the trajectories stored in the container
  for (G4int i=0; i<tc->entries(); ++i) {
    Trajectory* trj = dynamic_cast<Trajectory*>((*tc)[i]);
    if (!trj) continue;

    G4int trackid = trj->GetTrackID();

    G4ThreeVector ini_xyz = trj->GetInitialPosition();
    G4double ini_t = trj->GetInitialTime();
    G4ThreeVector xyz = trj->GetFinalPosition();
    G4double t = trj->GetFinalTime();

    G4String ini_volume = trj->GetInitialVolume();
    G4String volume = trj->GetDecayVolume();

    G4double mass = trj->GetParticleDefinition()->GetPDGMass();
    G4ThreeVector mom = trj->GetInitialMomentum();
    G4double energy = sqrt(mom.mag2() + mass*mass);

    G4ThreeVector final_mom = trj->GetFinalMomentum();
    G4double final_energy = sqrt(final_mom.mag2() + mass*mass);
   
    float ini_pos[4] = {(float)ini_xyz.x(), (float)ini_xyz.y(), (float)ini_xyz.z(), (float)ini_t};
    float final_pos[4] = {(float)xyz.x(), (float)xyz.y(), (float)xyz.z(), (float)t};
    float momentum[3] = {(float)mom.x(), (float)mom.y(), (float)mom.z()};
    float final_momentum[3] = {(float)final_mom.x(), (float)final_mom.y(), (float)final_mom.z()};
    float kin_energy = energy - mass;
    
    char primary = 0;
    G4int mother_id = 0;
    if (!trj->GetParentID()) {
      primary = 1;
    } else {
      mother_id = trj->GetParentID();
    }
 
    _h5writer->WriteParticleInfo(_nevt, trackid, trj->GetParticleName().c_str(),
				 primary, mother_id,
				 ini_pos[0], ini_pos[1], ini_pos[2], ini_pos[3],
				 final_pos[0], final_pos[1], final_pos[2], final_pos[3],
				 ini_volume.c_str(), volume.c_str(),
				 momentum[0],  momentum[1], momentum[2],
                                 final_momentum[0], final_momentum[1], final_momentum[2],
				 kin_energy, trj->GetCreatorProcess().c_str());
  }
}


G4bool CevnsPersistencyManager::Store(const G4Run*)
{
  // Store the number of events to be processed
  NexusApp* app = (NexusApp*) G4RunManager::GetRunManager();
  G4int num_events = app->GetNumberOfEventsToBeProcessed();

  G4String key = "num_events";
  _h5writer->WriteRunInfo(key, std::to_string(num_events).c_str());
  key = "saved_events";
  _h5writer->WriteRunInfo(key, std::to_string(_saved_evts).c_str());
  SaveConfigurationInfo(_historyFile_init);
  SaveConfigurationInfo(_historyFile_conf);

  return true;
}

void CevnsPersistencyManager::SaveConfigurationInfo(G4String file_name)
{
  std::ifstream history(file_name, std::ifstream::in);
  while (history.good()) {

    std::string key, value;
    std::getline(history, key, ' ');
    std::getline(history, value);

    if (key != "") {
      _h5writer->WriteRunInfo(key.c_str(), value.c_str());
    }
    
  }
  
  history.close();
}
