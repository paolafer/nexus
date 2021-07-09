//
//
////////////////////////////////////////////////////////////////////////

#ifndef LXeScintillation_h
#define LXeScintillation_h

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"

#include "G4EmSaturation.hh"


/////////////////////
// Class Definition
/////////////////////

class LXeScintillation : public G4VRestDiscreteProcess
{

public:

	////////////////////////////////
	// Constructors and Destructor
	////////////////////////////////

	explicit LXeScintillation(const G4String& processName = "LXeScintillation",
                                 G4ProcessType type = fElectromagnetic);
	~LXeScintillation();

private:

        LXeScintillation(const LXeScintillation &right) = delete;

        //////////////
        // Operators
        //////////////

        LXeScintillation& operator=(const LXeScintillation &right) = delete;

public:

        ////////////
        // Methods
        ////////////

        // LXeScintillation Process has both PostStepDoIt (for energy
        // deposition of particles in flight) and AtRestDoIt (for energy
        // given to the medium by particles at rest)

        G4bool IsApplicable(
          const G4ParticleDefinition& aParticleType) override;
        // Returns true -> 'is applicable', for any particle type except
        // for an 'opticalphoton' and for short-lived particles

        void BuildPhysicsTable(
          const G4ParticleDefinition& aParticleType) override;
        // Build table at the right time

        G4double GetMeanFreePath(const G4Track& aTrack,
                                       G4double ,
                                       G4ForceCondition* ) override;
        // Returns infinity; i. e. the process does not limit the step,
        // but sets the 'StronglyForced' condition for the DoIt to be
        // invoked at every step.

        G4double GetMeanLifeTime(const G4Track& aTrack,
                                 G4ForceCondition* ) override;
        // Returns infinity; i. e. the process does not limit the time,
        // but sets the 'StronglyForced' condition for the DoIt to be
        // invoked at every step.

        G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                            const G4Step&  aStep) override;
        G4VParticleChange* AtRestDoIt (const G4Track& aTrack,
                                       const G4Step& aStep) override;

        G4double GetScintillationYieldByParticleType(const G4Track &aTrack,
						     const G4Step &aStep);
        // Returns the number of scintillation photons calculated when
        // scintillation depends on the particle type and energy
        // deposited (includes nonlinear dependendency)

        // These are the methods implementing the scintillation process.

        void SetTrackSecondariesFirst(const G4bool state);
        // If set, the primary particle tracking is interrupted and any
        // produced scintillation photons are tracked next. When all
        // have been tracked, the tracking of the primary resumes.

        G4bool GetTrackSecondariesFirst() const;
        // Returns the boolean flag for tracking secondaries first.

        void SetFiniteRiseTime(const G4bool state);
        // If set, the LXeScintillation process expects the user to have
        // set the constant material property FAST/SLOWSCINTILLATIONRISETIME.

        G4bool GetFiniteRiseTime() const;
        // Returns the boolean flag for a finite scintillation rise time.

        void SetScintillationYieldFactor(const G4double yieldfactor);
        // Called to set the scintillation photon yield factor, needed when
        // the yield is different for different types of particles. This
        // scales the yield obtained from the G4MaterialPropertiesTable.

        G4double GetScintillationYieldFactor() const;
        // Returns the photon yield factor.

        void SetScintillationExcitationRatio(const G4double ratio);
        // Called to set the scintillation exciation ratio, needed when
        // the scintillation level excitation is different for different
        // types of particles. This overwrites the YieldRatio obtained
        // from the G4MaterialPropertiesTable.

        G4double GetScintillationExcitationRatio() const;
        // Returns the scintillation level excitation ratio.

        G4PhysicsTable* GetFastIntegralTable() const;
        // Returns the address of the fast scintillation integral table.

        G4PhysicsTable* GetSlowIntegralTable() const;
        // Returns the address of the slow scintillation integral table.

        void AddSaturation(G4EmSaturation* sat);
        // Adds Birks Saturation to the process.

        void RemoveSaturation();
        // Removes the Birks Saturation from the process.

        G4EmSaturation* GetSaturation() const;
        // Returns the Birks Saturation.

        void SetScintillationByParticleType(const G4bool );
        // Called by the user to set the scintillation yield as a function
        // of energy deposited by particle type

        G4bool GetScintillationByParticleType() const;
        // Return the boolean that determines the method of scintillation
        // production

        void SetScintillationTrackInfo(const G4bool trackType);
        // Call by the user to set the G4ScintillationTrackInformation
        // to scintillation photon track

        G4bool GetScintillationTrackInfo() const;
        // Return the boolean for whether or not the
        // G4ScintillationTrackInformation is set to the scint. photon track

        void SetStackPhotons(const G4bool );
        // Call by the user to set the flag for stacking the scint. photons

        G4bool GetStackPhotons() const;
        // Return the boolean for whether or not the scint. photons are stacked

        G4int GetNumPhotons() const;
        // Returns the current number of scint. photons (after PostStepDoIt)

        void DumpPhysicsTable() const;
        // Prints the fast and slow scintillation integral tables.

protected:

        void BuildThePhysicsTable();
        // It builds either the fast or slow scintillation integral table;
        // or both.

        ///////////////////////
        // Class Data Members
        ///////////////////////

        G4PhysicsTable* fFastIntegralTable;
        G4PhysicsTable* fSlowIntegralTable;

private:

        G4bool fTrackSecondariesFirst;
        G4bool fFiniteRiseTime;

        G4double fYieldFactor;

        G4double fExcitationRatio;

        G4bool fScintillationByParticleType;

        G4bool fScintillationTrackInfo;

        G4bool fStackingFlag;

        G4int fNumPhotons;

        G4double single_exp(G4double t, G4double tau2);
        G4double bi_exp(G4double t, G4double tau1, G4double tau2);

        // emission time distribution when there is a finite rise time
        G4double sample_time(G4double tau1, G4double tau2);

        G4EmSaturation* fEmSaturation;

};

////////////////////
// Inline methods
////////////////////

inline
void LXeScintillation::SetTrackSecondariesFirst(const G4bool state)
{
        fTrackSecondariesFirst = state;
}

inline
G4bool LXeScintillation::GetTrackSecondariesFirst() const
{
        return fTrackSecondariesFirst;
}

inline
void LXeScintillation::SetFiniteRiseTime(const G4bool state)
{
        fFiniteRiseTime = state;
}

inline
G4bool LXeScintillation::GetFiniteRiseTime() const
{
        return fFiniteRiseTime;
}

inline
void LXeScintillation::SetScintillationYieldFactor(const G4double yieldfactor)
{
        fYieldFactor = yieldfactor;
}

inline
G4double LXeScintillation::GetScintillationYieldFactor() const
{
        return fYieldFactor;
}

inline
void LXeScintillation::SetScintillationExcitationRatio(const G4double ratio)
{
        fExcitationRatio = ratio;
}

inline
G4double LXeScintillation::GetScintillationExcitationRatio() const
{
        return fExcitationRatio;
}

inline
G4PhysicsTable* LXeScintillation::GetSlowIntegralTable() const
{
        return fSlowIntegralTable;
}

inline
G4PhysicsTable* LXeScintillation::GetFastIntegralTable() const
{
        return fFastIntegralTable;
}

inline
void LXeScintillation::AddSaturation(G4EmSaturation* sat)
{
        fEmSaturation = sat;
}

inline
void LXeScintillation::RemoveSaturation()
{
        fEmSaturation = nullptr;
}

inline
G4EmSaturation* LXeScintillation::GetSaturation() const
{
        return fEmSaturation;
}

inline
G4bool LXeScintillation::GetScintillationByParticleType() const
{
        return fScintillationByParticleType;
}

inline
void LXeScintillation::SetScintillationTrackInfo(const G4bool trackType)
{
        fScintillationTrackInfo = trackType;
}

inline
G4bool LXeScintillation::GetScintillationTrackInfo() const
{
        return fScintillationTrackInfo;
}

inline
void LXeScintillation::SetStackPhotons(const G4bool stackingFlag)
{
        fStackingFlag = stackingFlag;
}

inline
G4bool LXeScintillation::GetStackPhotons() const
{
        return fStackingFlag;
}

inline
G4int LXeScintillation::GetNumPhotons() const
{
        return fNumPhotons;
}


inline
G4double LXeScintillation::single_exp(G4double t, G4double tau2)
{
         return std::exp(-1.0*t/tau2)/tau2;
}

inline
G4double LXeScintillation::bi_exp(G4double t, G4double tau1, G4double tau2)
{
         return std::exp(-1.0*t/tau2)*(1-std::exp(-1.0*t/tau1))/tau2/tau2*(tau1+tau2);
}

#endif /* LXeScintillation_h */
