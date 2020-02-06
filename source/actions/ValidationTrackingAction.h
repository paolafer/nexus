// ----------------------------------------------------------------------------
///  \file   ValidationTrackingAction.h
///  \brief  General-purpose user tracking action.
///
///  \author   <justo.martin-albo@ific.uv.es>
///  \date     25 Mar 2013
///  \version  $Id$
///
///  Copyright (c) 2013 NEXT Collaboration. All rights reserved.
// ----------------------------------------------------------------------------

#ifndef VALIDATION_TRACKING_ACTION_H
#define VALIDATION_TRACKING_ACTION_H

#include <G4UserTrackingAction.hh>

class G4Track;

namespace nexus {

  // General-purpose user tracking action

  class ValidationTrackingAction: public G4UserTrackingAction
  {
  public:
    /// Constructor
    ValidationTrackingAction();
    /// Destructor
    virtual ~ValidationTrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);

  };

}

#endif
