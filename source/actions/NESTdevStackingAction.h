#ifndef NEST_DEV_STACKING_ACTION_H
#define NEST_DEV_STACKING_ACTION_H

#include <G4UserStackingAction.hh>
#include <NESTStackingAction.hh>

namespace nexus {

  // General-purpose user stacking action

  class NESTdevStackingAction: public NESTStackingAction
  {
  public:
    /// Constructor
    NESTdevStackingAction();
    /// Destructor
    ~NESTdevStackingAction();
  };

} // end namespace nexus

#endif
