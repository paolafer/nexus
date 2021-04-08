// ----------------------------------------------------------------------------
// nexus | PersistencyFactory.h
//
// This class instantiates the persistency manager
// that the user specifies via configuration parameters.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef PERSISTENCY_FACTORY_H
#define PERSISTENCY_FACTORY_H

#include <G4String.hh>

class G4VPersistencyManager;
class G4GenericMessenger;


//namespace nexus {

  class PersistencyFactory
  {
  public:
    PersistencyFactory();
    ~PersistencyFactory();

    G4VPersistencyManager* CreatePersistencyManager(G4String init,
                                                    std::vector<G4String>& macros,
                                                    std::vector<G4String>& delayed) const;

  private:
    G4GenericMessenger* msg_;

    G4String persman_name_; ///< Name of the user persistency manager
  };

//} // end namespace nexus

#endif
