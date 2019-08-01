#ifndef NU_TRACKING_ACTION_H
#define NU_TRACKING_ACTION_H

#include <G4UserTrackingAction.hh>

class G4Track;


namespace nexus {

  // General-purpose user tracking action

  class NuTrackingAction: public G4UserTrackingAction
  {
  public:
    /// Constructor
    NuTrackingAction();
    /// Destructor
    virtual ~NuTrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);
  };

}

#endif
