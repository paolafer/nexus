// ----------------------------------------------------------------------------
// nexus | PersistencyFactory.cc
//
// This class instantiates the the persistency manager
// that the user specifies via configuration parameters.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "PersistencyFactory.h"

#include <G4GenericMessenger.hh>
#include <G4VPersistencyManager.hh>

using namespace nexus;


PersistencyFactory::PersistencyFactory(): msg_(0)
{
  msg_ = new G4GenericMessenger(this, "/Persistency/");
  msg_->DeclareProperty("RegisterPersistencyManager", persman_name_, "");
}


PersistencyFactory::~PersistencyFactory()
{
  delete msg_;
}


//////////////////////////////////////////////////////////////////////
#include "PersistencyManager.h"


G4VPersistencyManager* PersistencyFactory::CreatePersistencyManager(G4String init,
                                                                    std::vector<G4String>& macros,
                                                                    std::vector<G4String>& delayed) const
{
  G4VPersistencyManager* p = 0;

  if (persman_name_ == "NEXT") p = new PersistencyManager(init, macros, delayed);
  else {
    G4String err = "Unknown user persistency manager: " + persman_name_;
    G4Exception("[PersistencyFactory]", "CreatePersistencyManager()", JustWarning, err);
  }
  return p;
}


