// ----------------------------------------------------------------------------
// nexus | NexusPhysics.cc
//
// This class registers any new physics process defined in nexus.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "NexusPhysics.h"
#include "WavelengthShifting.h"

#include "NESTProc.hh"
#include "PetaloDetector.hh"

#include <G4Scintillation.hh>
#include <G4GenericMessenger.hh>
#include <G4OpticalPhoton.hh>
#include <G4Gamma.hh>
#include <G4Electron.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessTable.hh>
#include <G4StepLimiter.hh>
#include <G4FastSimulationManagerProcess.hh>
#include <G4PhysicsConstructorFactory.hh>


namespace nexus {

  /// Macro that allows the use of this physics constructor
  /// with the generic physics list
  G4_DECLARE_PHYSCONSTR_FACTORY(NexusPhysics);


  NexusPhysics::NexusPhysics():
    G4VPhysicsConstructor("NexusPhysics"), risetime_(false), noCompt_(false),
    nest_(false), sc_yield_factor_(0.818)
  {
    msg_ = new G4GenericMessenger(this, "/PhysicsList/Nexus/",
      "Control commands of the nexus physics list.");

    msg_->DeclareProperty("scintRiseTime", risetime_,
      "True if LYSO is used");

    msg_->DeclareProperty("offCompt", noCompt_,
			  "Switch off Compton Scattering.");

    msg_->DeclareProperty("nest", nest_,
      "True if NEST is used for scintillation");
    msg_->DeclareProperty("sc_yield_factor", sc_yield_factor_,
      "Factor to reduce scintillation yield in NEST");
  }



  NexusPhysics::~NexusPhysics()
  {
    delete msg_;
    delete wls_;
  }



  void NexusPhysics::ConstructParticle()
  {
    G4OpticalPhoton::Definition();
    NEST::NESTThermalElectron::Definition();
    //G4OpticalPhoton::OpticalPhotonDefinition();
  }



  void NexusPhysics::ConstructProcess()
  {
    G4ProcessManager* pmanager = 0;

    // Add our own wavelength shifting process for the optical photon
    pmanager = G4OpticalPhoton::Definition()->GetProcessManager();
    if (!pmanager) {
      G4Exception("ConstructProcess()", "[NexusPhysics]", FatalException,
        "G4OpticalPhoton without a process manager.");
    }
    //WavelengthShifting* wls = new WavelengthShifting();
    wls_ = new WavelengthShifting();
    pmanager->AddDiscreteProcess(wls_);

    // Add rise time to scintillation
    if (risetime_) {
      pmanager  = G4Electron::Definition()->GetProcessManager();
      G4Scintillation*  theScintillationProcess =
	(G4Scintillation*)G4ProcessTable::GetProcessTable()->
	FindProcess("Scintillation", G4Electron::Definition());
      theScintillationProcess->SetFiniteRiseTime(true);
    }


    if (noCompt_) {
      pmanager  = G4Gamma::Definition()->GetProcessManager();
       G4VProcess* cs = G4ProcessTable::GetProcessTable()->
        FindProcess("compt", G4Gamma::Definition());
       pmanager->RemoveProcess(cs);
    }

    if (nest_) {
      PetaloDetector* petalo = new PetaloDetector();
      NEST::NESTcalc* petaloCalc = new NEST::NESTcalc(petalo);
      NEST::NESTProc* theNESTScintillationProcess =
	new NEST::NESTProc("S1", fElectromagnetic, petaloCalc, petalo);
      theNESTScintillationProcess->SetDetailedSecondaries(true); // this is to use the full scintillation spectrum of LXe.
      theNESTScintillationProcess->SetStackElectrons(false);
      theNESTScintillationProcess->SetScintillationYieldFactor(sc_yield_factor_);

      auto aParticleIterator = GetParticleIterator();
      aParticleIterator->reset();
      while ((*aParticleIterator)()) {
        G4ParticleDefinition* particle = aParticleIterator->value();
        if (theNESTScintillationProcess->IsApplicable(*particle)) {
          pmanager = particle->GetProcessManager();
          pmanager->AddProcess(theNESTScintillationProcess, ordDefault+1, ordInActive, ordDefault+1);
        }
      }
    }

  }



} // end namespace nexus
