// ----------------------------------------------------------------------------
///  \file   ActionsFactory.h
///  \brief  
/// 
///  \author   <justo.martin-albo@ific.uv.es>
///  \date     13 March 2013
///  \version  $Id$
///
///  Copyright (c) 2009-2013 NEXT Collaboration. All rights reserved.
// ----------------------------------------------------------------------------

#ifndef __DEFAULT_RUN_ACTION__
#define __DEFAULT_RUN_ACTION__ 

#include <G4UserRunAction.hh>
#include <G4String.hh>

class G4GenericMessenger;


namespace nexus {

  /// TODO. CLASS DESCRIPTION

  class DefaultRunAction: public G4UserRunAction
  {
  public:
    /// Constructor
    DefaultRunAction();
    /// Destructor
    ~DefaultRunAction();
  
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

  private:

    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;

    /// Name of files with ntuples for posterior plotting
    G4String analysis_file_;
  };

}

#endif

