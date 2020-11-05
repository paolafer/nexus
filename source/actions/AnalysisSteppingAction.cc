// ----------------------------------------------------------------------------
//  $Id$
//
//  Author : <>    
//
//  Copyright (c) 2013 NEXT Collaboration. All rights reserved.
// ----------------------------------------------------------------------------

#include "AnalysisSteppingAction.h"

#include <G4Step.hh>
#include <G4SteppingManager.hh>
#include <G4ProcessManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4Electron.hh>
#include <G4GenericMessenger.hh>

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <iostream>
#include <sstream>

using namespace nexus;
using namespace CLHEP;


AnalysisSteppingAction::AnalysisSteppingAction(): G4UserSteppingAction(), file_no_(0) 
{
  _msg = new G4GenericMessenger(this, "/Actions/AnalysisSteppingAction/");
  _msg->DeclareProperty("file_number", file_no_, "");
  
  detected = 0;
  not_det = 0;
 
  times.clear();
  wavelengths.clear();
}



AnalysisSteppingAction::~AnalysisSteppingAction()
{


  G4cout << "Detected photons = " << detected << G4endl;
  G4cout << "Non detected photons = " << not_det << G4endl;
   
   
   double first = 1000.*second ;
   
   for (unsigned int i=0; i< times.size(); ++i) {
     if (times[i] < first) {
       first = times[i];
     }
   }
   //G4cout << first/picosecond << G4endl;
   //std::vector<double> times_sub;



   std::ostringstream file_number;
   file_number << file_no_;
   G4String filename = "PhotonVelocitiesCherLXe."+file_number.str()+".root";
  
  
 
  G4double total_counts = 0;
  detectorCounts::iterator it = my_counts.begin();
  while (it != my_counts.end()) {
    G4cout << "Detector " << it->first << ": " << it->second << " counts" << G4endl;
    total_counts += it->second;
    it ++;
  }
  G4cout << "TOTAL COUNTS: " << total_counts << G4endl;
}



void AnalysisSteppingAction::UserSteppingAction(const G4Step* step)
{ 
  G4ParticleDefinition* pdef = step->GetTrack()->GetDefinition();

  //  Check whether the track is an optical photon  
  if (pdef != G4OpticalPhoton::Definition()) return;


  G4StepPoint* point1 = step->GetPreStepPoint();
  G4StepPoint* point2 = step->GetPostStepPoint();
  G4TouchableHandle touch1 = point1->GetTouchableHandle();
  G4TouchableHandle touch2 = point2->GetTouchableHandle();
  G4String vol1name = touch1->GetVolume()->GetName();
  G4String vol2name = touch2->GetVolume()->GetName();
  //G4Track* track = step->GetTrack();

  //G4String proc_name = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  //G4int copy_no = step->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);

  // Retrieve the pointer to the optical boundary process.
  // We do this only once per run defining our local pointer as static.
  static G4OpBoundaryProcess* boundary = 0;
  
  if (!boundary) { // the pointer is not defined yet
    // Get the list of processes defined for the optical photon
    // and loop through it to find the optical boundary process.
    G4ProcessVector* pv = pdef->GetProcessManager()->GetProcessList();
    for (G4int i=0; i<pv->size(); i++) {
      if ((*pv)[i]->GetProcessName() == "OpBoundary") {
	boundary = (G4OpBoundaryProcess*) (*pv)[i];
	break;
      }
    }
  }

  if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {

    // if boundary->GetStatus() == 2 in SiPMpet refraction takes place
    if (boundary->GetStatus() == Detection) {
	detected = detected + 1;	
	double distance = 
	  std::pow(point2->GetPosition().getX() - point1->GetPosition().getX(), 2) + 
	  std::pow(point2->GetPosition().getY() - point1->GetPosition().getY(), 2)  + std::pow(point2->GetPosition().getZ() - point1->GetPosition().getZ(), 2) ;
	distance = std::sqrt(distance);
	G4double lambda = h_Planck*c_light/step->GetTrack()->GetKineticEnergy()/nanometer;
	times.push_back(step->GetDeltaTime());
	wavelengths.push_back(distance/ step->GetDeltaTime());

	G4String detector_name = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
	//G4cout << "##### Sensitive Volume: " << detector_name << G4endl;
	
	detectorCounts::iterator it = my_counts.find(detector_name);
	if (it != my_counts.end()) my_counts[it->first] += 1;
	else my_counts[detector_name] = 1;
      } else {
	not_det = not_det + 1;
      }
    
      //	G4cout << "check: " << velocity << ", " << track_velocity << G4endl;


  }

  return;
}
