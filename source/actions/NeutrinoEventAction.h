#ifndef NEUTRINO_EVENT_ACTION_H
#define NEUTRINO_EVENT_ACTION_H

#include <G4UserEventAction.hh>
#include <globals.hh>

class G4Event;
class G4GenericMessenger;

namespace nexus {
    
  /// This class is a general-purpose event run action.
  
  class NeutrinoEventAction: public G4UserEventAction
  {
  public:
    /// Constructor
    NeutrinoEventAction();
    /// Destructor
    ~NeutrinoEventAction();
    
    /// Hook at the beginning of the event loop
    void BeginOfEventAction(const G4Event*);
    /// Hook at the end of the event loop
    void EndOfEventAction(const G4Event*);

  private:
    // G4GenericMessenger* _msg;
    G4int _nevt, _nupdate;
  };
  
} // namespace nexus

#endif
