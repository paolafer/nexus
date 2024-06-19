// ----------------------------------------------------------------------------
// nexus | OpBoundaryProcess.h
//
// This class is a copy of G4OpBoundaryProcess, modified with the optical
// photoelectric effect.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------


#ifndef G4OpBoundaryProcess_h
#define G4OpBoundaryProcess_h 1

#include "G4OpticalPhoton.hh"
#include "G4OpticalSurface.hh"
#include "G4RandomTools.hh"
#include "G4VDiscreteProcess.hh"

namespace nexus {
enum OpBoundaryProcessStatus
{
  Undefined,
  Transmission,
  FresnelRefraction,
  FresnelReflection,
  TotalInternalReflection,
  LambertianReflection,
  LobeReflection,
  SpikeReflection,
  BackScattering,
  Absorption,
  Detection,
  NotAtBoundary,
  SameMaterial,
  StepTooSmall,
  NoRINDEX,
  PolishedLumirrorAirReflection,
  PolishedLumirrorGlueReflection,
  PolishedAirReflection,
  PolishedTeflonAirReflection,
  PolishedTiOAirReflection,
  PolishedTyvekAirReflection,
  PolishedVM2000AirReflection,
  PolishedVM2000GlueReflection,
  EtchedLumirrorAirReflection,
  EtchedLumirrorGlueReflection,
  EtchedAirReflection,
  EtchedTeflonAirReflection,
  EtchedTiOAirReflection,
  EtchedTyvekAirReflection,
  EtchedVM2000AirReflection,
  EtchedVM2000GlueReflection,
  GroundLumirrorAirReflection,
  GroundLumirrorGlueReflection,
  GroundAirReflection,
  GroundTeflonAirReflection,
  GroundTiOAirReflection,
  GroundTyvekAirReflection,
  GroundVM2000AirReflection,
  GroundVM2000GlueReflection,
  Dichroic,
  CoatedDielectricReflection,
  CoatedDielectricRefraction,
  CoatedDielectricFrustratedTransmission
};

class OpBoundaryProcess : public G4VDiscreteProcess
{
 public:
  explicit OpBoundaryProcess(const G4String& processName = "OpBoundary",
                             G4ProcessType type          = fOptical);
  virtual ~OpBoundaryProcess();

  virtual G4bool IsApplicable(
    const G4ParticleDefinition& aParticleType) override;
  // Returns true -> 'is applicable' only for an optical photon.

  virtual G4double GetMeanFreePath(const G4Track&, G4double,
                                   G4ForceCondition* condition) override;
  // Returns infinity; i. e. the process does not limit the step, but sets the
  // 'Forced' condition for the DoIt to be invoked at every step. However, only
  // at a boundary will any action be taken.

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                  const G4Step& aStep) override;
  // This is the method implementing boundary processes.

  virtual OpBoundaryProcessStatus GetStatus() const;
  // Returns the current status.

  virtual void SetInvokeSD(G4bool);
  // Set flag for call to InvokeSD method.

  virtual void PreparePhysicsTable(const G4ParticleDefinition&) override;

  virtual void Initialise();

  void SetVerboseLevel(G4int);

 private:
  OpBoundaryProcess(const OpBoundaryProcess& right) = delete;
  OpBoundaryProcess& operator=(const OpBoundaryProcess& right) = delete;

  G4bool G4BooleanRand(const G4double prob) const;

  G4ThreeVector GetFacetNormal(const G4ThreeVector& Momentum,
                               const G4ThreeVector& Normal) const;

  void DielectricMetal();
  void DielectricDielectric();

  void DielectricLUT();
  void DielectricLUTDAVIS();

  void DielectricDichroic();
  void CoatedDielectricDielectric();

  void ChooseReflection();
  void DoAbsorption();
  void DoOpPhotoelectric();
  void DoReflection();

  G4double GetIncidentAngle();
  // Returns the incident angle of optical photon

  G4double GetReflectivity(G4double E1_perp, G4double E1_parl,
                           G4double incidentangle, G4double RealRindex,
                           G4double ImaginaryRindex);
  // Returns the Reflectivity on a metallic surface

  G4double GetReflectivityThroughThinLayer(G4double sinTL, G4double E1_perp,
                                           G4double E1_parl, G4double wavelength,
                                           G4double cost1, G4double cost2);
  // Returns the Reflectivity on a coated surface

  void CalculateReflectivity();

  void BoundaryProcessVerbose() const;

  // Invoke SD for post step point if the photon is 'detected'
  G4bool InvokeSD(const G4Step* step);

  G4ThreeVector fOldMomentum;
  G4ThreeVector fOldPolarization;

  G4ThreeVector fNewMomentum;
  G4ThreeVector fNewPolarization;

  G4ThreeVector fGlobalNormal;
  G4ThreeVector fFacetNormal;

  const G4Material* fMaterial1;
  const G4Material* fMaterial2;

  G4OpticalSurface* fOpticalSurface;

  G4MaterialPropertyVector* fRealRIndexMPV;
  G4MaterialPropertyVector* fImagRIndexMPV;
  G4Physics2DVector* fDichroicVector;

  G4double fPhotonMomentum;
  G4double fRindex1;
  G4double fRindex2;

  G4double fSint1;

  G4double fReflectivity;
  G4double fEfficiency;
  G4double fTransmittance;
  G4double fSurfaceRoughness;

  G4double fProb_sl, fProb_ss, fProb_bs;
  G4double fCarTolerance;

  // Used by CoatedDielectricDielectric()
  G4double fCoatedRindex, fCoatedThickness;

  G4OpBoundaryProcessStatus fStatus;
  G4OpticalSurfaceModel fModel;
  G4OpticalSurfaceFinish fFinish;

  G4int f_iTE, f_iTM;

  G4int fNumSmallStepWarnings = 0; // number of times small step warning printed
  G4int fNumBdryTypeWarnings = 0;  // number of times boundary type warning printed

  size_t idx_dichroicX      = 0;
  size_t idx_dichroicY      = 0;
  size_t idx_rindex1        = 0;
  size_t idx_rindex_surface = 0;
  size_t idx_reflect        = 0;
  size_t idx_eff            = 0;
  size_t idx_trans          = 0;
  size_t idx_lobe           = 0;
  size_t idx_spike          = 0;
  size_t idx_back           = 0;
  size_t idx_rindex2        = 0;
  size_t idx_groupvel       = 0;
  size_t idx_rrindex        = 0;
  size_t idx_irindex        = 0;
  size_t idx_coatedrindex   = 0;

  // Used by CoatedDielectricDielectric()
  G4bool fCoatedFrustratedTransmission = true;

  G4bool fInvokeSD;
};

////////////////////
// Inline methods
////////////////////

inline G4bool OpBoundaryProcess::G4BooleanRand(const G4double prob) const
{
  // Returns a random boolean variable with the specified probability
  return (G4UniformRand() < prob);
}

inline G4bool OpBoundaryProcess::IsApplicable(
  const G4ParticleDefinition& aParticleType)
{
  return (&aParticleType == G4OpticalPhoton::OpticalPhoton());
}

inline OpBoundaryProcessStatus G4OpBoundaryProcess::GetStatus() const
{
  return fStatus;
}

inline void OpBoundaryProcess::ChooseReflection()
{
  G4double rand = G4UniformRand();
  if(rand < fProb_ss)
  {
    fStatus      = SpikeReflection;
    fFacetNormal = fGlobalNormal;
  }
  else if(rand < fProb_ss + fProb_sl)
  {
    fStatus = LobeReflection;
  }
  else if(rand < fProb_ss + fProb_sl + fProb_bs)
  {
    fStatus = BackScattering;
  }
  else
  {
    fStatus = LambertianReflection;
  }
}

inline void OpBoundaryProcess::DoAbsorption()
{
  fStatus = Absorption;

  if(G4BooleanRand(fEfficiency))
  {
    // EnergyDeposited =/= 0 means: photon has been detected
    fStatus = Detection;
    aParticleChange.ProposeLocalEnergyDeposit(fPhotonMomentum);
  }
  else
  {
    aParticleChange.ProposeLocalEnergyDeposit(0.0);
  }

  fNewMomentum     = fOldMomentum;
  fNewPolarization = fOldPolarization;

  aParticleChange.ProposeTrackStatus(fStopAndKill);
}


inline void OpBoundaryProcess::DoReflection()
{
  if(fStatus == LambertianReflection)
  {
    fNewMomentum = G4LambertianRand(fGlobalNormal);
    fFacetNormal = (fNewMomentum - fOldMomentum).unit();
  }
  else if(fFinish == ground)
  {
    fStatus = LobeReflection;
    if(!fRealRIndexMPV || !fImagRIndexMPV)
    {
      fFacetNormal = GetFacetNormal(fOldMomentum, fGlobalNormal);
    }
    // else
      // complex ref. index to be implemented
    fNewMomentum =
      fOldMomentum - (2. * fOldMomentum * fFacetNormal * fFacetNormal);
  }
  else
  {
    fStatus      = SpikeReflection;
    fFacetNormal = fGlobalNormal;
    fNewMomentum =
      fOldMomentum - (2. * fOldMomentum * fFacetNormal * fFacetNormal);
  }
  fNewPolarization =
    -fOldPolarization + (2. * fOldPolarization * fFacetNormal * fFacetNormal);
}

}
#endif /* G4OpBoundaryProcess_h */

