// ----------------------------------------------------------------------------
//  $Id$
//
//  Author : <justo.martin-albo@ific.uv.es>
//  Created: 13 March 2013
//
//  Copyright (c) 2013 NEXT Collaboration. All rights reserved.
// ----------------------------------------------------------------------------

#include "DefaultRunAction.h"
#include "MyAnalysis.hh"

#include <G4Run.hh>
#include <G4GenericMessenger.hh>

using namespace nexus;



DefaultRunAction::DefaultRunAction(): G4UserRunAction(),
				      analysis_file_("analysis")
{
  msg_ = new G4GenericMessenger(this, "/Actions/DefaultRunAction/",
                                  "Control commands of default run action.");
  msg_->DeclareProperty("analysisFile", analysis_file_,
			"Name of analysis file");
}



DefaultRunAction::~DefaultRunAction()
{
}



void DefaultRunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  // Create/get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);

  // Open an output file
  analysisManager->OpenFile(analysis_file_);

  // Creation of ntuple
  analysisManager->CreateNtuple("MyNtuple", "Whatever");
  analysisManager->CreateNtupleDColumn("GammaE");
  analysisManager->FinishNtuple();
}


void DefaultRunAction::EndOfRunAction(const G4Run* run)
{
  G4cout << "### Run " << run->GetRunID() << " end." << G4endl;

  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Write and close the output file
  analysisManager->Write();
  analysisManager->CloseFile();
}
