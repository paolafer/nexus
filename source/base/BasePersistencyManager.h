#ifndef BASE_PERSISTENCY_MANAGER_H
#define BASE_PERSISTENCY_MANAGER_H

#include <G4VPersistencyManager.hh>
#include <map>
#include <vector>


class BasePersistencyManager: public G4VPersistencyManager
  {
  public:
    virtual void CloseFile() = 0;
    
    
  };

#endif
