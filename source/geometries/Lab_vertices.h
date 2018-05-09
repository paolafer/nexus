// ----------------------------------------------------------------------------
///  \file 
///  \brief
///
///  \author   <paolafer@ific.uv.es>
///  \date     2015
///  \version  $Id: Lab.h 9593 2014-02-13 16:38:56Z paola $
///
///  Copyright (c) 2015-2017 NEXT Collaboration. All rights reserved.
// ----------------------------------------------------------------------------

#ifndef LAB_V_H
#define LAB_V_H

#include "BaseGeometry.h"

#include <TFile.h>

class G4GenericMessenger;
namespace nexus {class Pet2boxes;}
namespace nexus {class PetLYSObox;}
namespace nexus {class PetLXeCellCheck;}
namespace nexus {class PetLYSOCell;}

namespace nexus {

  class Lab_vertices: public BaseGeometry
  {
  public:
    /// Constructor
    Lab_vertices();
    /// Destructor
    ~Lab_vertices();

    /// Return vertex within region <region> of the chamber
    virtual G4ThreeVector GenerateVertex(const G4String& region) const;
    std::pair<G4ThreeVector, G4ThreeVector> GenerateVertices(const G4String& /*region*/) const;

    virtual void Construct();

  private:
    /// Messenger for the definition of control commands
    G4GenericMessenger* _msg;

    PetLXeCellCheck* module_;

    TFile* file_;
    G4float px1_, py1_, pz1_, px2_, py2_, pz2_;
    mutable G4int index_;
    mutable std::vector<std::pair<G4ThreeVector, G4ThreeVector> > vertices_;

    G4int starting_point_;
    G4String filename_;
    G4String type_;
    
  };

} // end namespace nexus

#endif
