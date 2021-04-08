// ----------------------------------------------------------------------------
// nexus | NexusApp.cc
//
// This class is the run manager of the nexus simulation. It takes care of
// setting up the simulation (geometry, physics lists, generators, actions),
// so that it is ready to be run.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "NexusApp.h"

#include "GeometryFactory.h"
#include "GeneratorFactory.h"
#include "ActionsFactory.h"
#include "PersistencyFactory.h"
#include "DetectorConstruction.h"
#include "PrimaryGeneration.h"
#include "PersistencyManager.h"
#include "BatchSession.h"

#include <G4GenericPhysicsList.hh>
#include <G4UImanager.hh>
#include <G4StateManager.hh>

using namespace nexus;



NexusApp::NexusApp(G4String init_macro): G4RunManager()
{
  // Create and configure a generic messenger for the app
  msg_ = new G4GenericMessenger(this, "/nexus/", "Nexus control commands.");

  // Define the command to register a configuration macro.
  // The user may invoke the command as many times as needed.
  msg_->DeclareMethod("RegisterMacro", &NexusApp::RegisterMacro, "");

  // Some commands, which we call 'delayed', only work if executed
  // after the initialization of the application. The user may include
  // them in configuration macros registered with the command defined below.
  msg_->DeclareMethod("RegisterDelayedMacro",
                      &NexusApp::RegisterDelayedMacro, "");

  // Define a command to set a seed for the random number generator.
  msg_->DeclareMethod("random_seed", &NexusApp::SetRandomSeed,
                      "Set a seed for the random number generator.");

  /////////////////////////////////////////////////////////

  // We will set now the user initialization class instances
  // in the run manager. In order to do so, we create first the factories
  // (the objects that construct the appropriate instances according
  // to user's input) so that the messenger commands are already defined
  // by the time we process the initialization macro.

  GeometryFactory  geomfctr;
  GeneratorFactory genfctr;
  ActionsFactory   actfctr;
  PersistencyFactory persfctr;

  // The physics lists are handled with Geant4's own 'factory'
  physicsList = new G4GenericPhysicsList();

  BatchSession* batch = new BatchSession(init_macro.c_str());
  batch->SessionStart();

  // Set the physics list in the run manager
  this->SetUserInitialization(physicsList);

  // Set the detector construction instance in the run manager
  DetectorConstruction* dc = new DetectorConstruction();
  dc->SetGeometry(geomfctr.CreateGeometry());
  this->SetUserInitialization(dc);

  // Set the primary generation instance in the run manager
  PrimaryGeneration* pg = new PrimaryGeneration();
  pg->SetGenerator(genfctr.CreateGenerator());
  this->SetUserAction(pg);

  // User interface
  G4UImanager* UI = G4UImanager::GetUIpointer();

  persfctr.CreatePersistencyManager(init_macro, macros_, delayed_);
  //PersistencyManager::Initialize(init_macro, macros_, delayed_);

  // Set the user action instances, if any, in the run manager

  if (UI->GetCurrentValues("/Actions/RegisterRunAction") != "")
    this->SetUserAction(actfctr.CreateRunAction());

  if (UI->GetCurrentValues("/Actions/RegisterEventAction") != "")
    this->SetUserAction(actfctr.CreateEventAction());

  if (UI->GetCurrentValues("/Actions/RegisterStackingAction") != "")
    this->SetUserAction(actfctr.CreateStackingAction());

  if (UI->GetCurrentValues("/Actions/RegisterTrackingAction") != "")
    this->SetUserAction(actfctr.CreateTrackingAction());

  if (UI->GetCurrentValues("/Actions/RegisterSteppingAction") != "")
    this->SetUserAction(actfctr.CreateSteppingAction());

  /////////////////////////////////////////////////////////

  // Set by default a random seed (system time) for the random
  // number generator
  SetRandomSeed(-1);

}



NexusApp::~NexusApp()
{
  // Close output file before finishing
  PersistencyManager* current = dynamic_cast<PersistencyManager*>
    (G4VPersistencyManager::GetPersistencyManager());
  current->CloseFile();

  delete msg_;
}



void NexusApp::RegisterMacro(G4String macro)
{
  // Store the name of the macro file
  macros_.push_back(macro);
}



void NexusApp::RegisterDelayedMacro(G4String macro)
{
  // Store the name of the macro file
  delayed_.push_back(macro);
}



void NexusApp::Initialize()
{
  // Execute all command macro files before initializing the app
  // so that all objects get configured
  // G4UImanager* UI = G4UImanager::GetUIpointer();

  for (unsigned int i=0; i<macros_.size(); i++) {
    ExecuteMacroFile(macros_[i].data());
  }

  G4RunManager::Initialize();

  for (unsigned int j=0; j<delayed_.size(); j++) {
    ExecuteMacroFile(delayed_[j].data());
  }
}



void NexusApp::ExecuteMacroFile(const char* filename)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4UIsession* batchSession = new BatchSession(filename, UI->GetSession());
  UI->SetSession(batchSession);
  G4UIsession* previousSession = UI->GetSession()->SessionStart();
  delete UI->GetSession();
  UI->SetSession(previousSession);
}



void NexusApp::SetRandomSeed(G4int seed)
{
  // Set the seed chosen by the user for the pseudo-random number
  // generator unless a negative number was provided, in which case
  // we will set as seed the system time.
  if (seed < 0) CLHEP::HepRandom::setTheSeed(time(0));
  else CLHEP::HepRandom::setTheSeed(seed);
}
