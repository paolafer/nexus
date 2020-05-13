// ----------------------------------------------------------------------------
//  $Id$
//
//  Author : <justo.martin-albo@ific.uv.es>
//  Created: 13 March 2013
//
//  Copyright (c) 2013 NEXT Collaboration. All rights reserved.
// ----------------------------------------------------------------------------

#include "ActionsFactory.h"

#include <G4GenericMessenger.hh>
#include <G4UserRunAction.hh>
#include <G4UserEventAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserStackingAction.hh>

using namespace nexus;


ActionsFactory::ActionsFactory(): _msg(0)
{
  _msg = new G4GenericMessenger(this, "/Actions/");
  _msg->DeclareProperty("RegisterRunAction",      _runact_name, "");
  _msg->DeclareProperty("RegisterEventAction",    _evtact_name, "");
  _msg->DeclareProperty("RegisterTrackingAction", _trkact_name, "");
  _msg->DeclareProperty("RegisterSteppingAction", _stpact_name, "");
  _msg->DeclareProperty("RegisterStackingAction", _stkact_name, "");
}


ActionsFactory::~ActionsFactory()
{
  delete _msg;
}


//////////////////////////////////////////////////////////////////////
#include "DefaultRunAction.h"


G4UserRunAction* ActionsFactory::CreateRunAction() const
{
  G4UserRunAction* p = 0;

  if (_runact_name == "DEFAULT") p = new DefaultRunAction();

  return p;
}


//////////////////////////////////////////////////////////////////////
#include "DefaultEventAction.h"
#include "ELSimEventAction.h"
#include "AnalysisEventAction.h"

G4UserEventAction* ActionsFactory::CreateEventAction() const
{
  G4UserEventAction* p = 0;

  if      (_evtact_name == "DEFAULT") p = new DefaultEventAction();
  else if (_evtact_name == "EL_SIM") p = new ELSimEventAction();
  else if (_evtact_name == "ANALYSIS") p = new AnalysisEventAction();
  else {
    G4String err = "Unknown user event action: " + _evtact_name;
    G4Exception("CreateEventAction()", "[ActionsFactory]", JustWarning, err);
  }

  return p;
}


//////////////////////////////////////////////////////////////////////
#include "DefaultTrackingAction.h"
#include "AnalysisTrackingAction.h"
#include "OpticalTrackingAction.h"

G4UserTrackingAction* ActionsFactory::CreateTrackingAction() const
{
  G4UserTrackingAction* p = 0;

  if (_trkact_name == "DEFAULT") p = new DefaultTrackingAction();
  else if (_trkact_name == "ANALYSIS") p = new AnalysisTrackingAction();
  else if (_trkact_name == "OPTICAL") p = new OpticalTrackingAction();
  else {
    G4String err = "Unknown user tracking action: " + _trkact_name;
    G4Exception("CreateTrackingAction()", "[ActionsFactory]",
      JustWarning, err);
  }

  return p;
}


//////////////////////////////////////////////////////////////////////
#include "DefaultSteppingAction.h"
#include "AnalysisSteppingAction.h"


G4UserSteppingAction* ActionsFactory::CreateSteppingAction() const
{
  G4UserSteppingAction* p = 0;

  if (_stpact_name == "DEFAULT") p = new DefaultSteppingAction();
  else if (_stpact_name == "ANALYSIS") p = new AnalysisSteppingAction();

  else {
    G4String err = "Unknown user stepping action: " + _stpact_name;
    G4Exception("CreateSteppingAction()", "[ActionsFactory]",
      JustWarning, err);
  }

  return p;
}


//////////////////////////////////////////////////////////////////////
#include "DefaultStackingAction.h"
#include "NESTdevStackingAction.h"

G4UserStackingAction* ActionsFactory::CreateStackingAction() const
{
  G4UserStackingAction* p = 0;

  if (_stkact_name == "DEFAULT") p = new DefaultStackingAction();
  if (_stkact_name == "NEST") p = new NESTdevStackingAction();

  else {
    G4String err = "Unknown user stacking action: " + _stkact_name;
    G4Exception("CreateStackingAction()", "[ActionsFactory]",
      JustWarning, err);
  }

  return p;
}
