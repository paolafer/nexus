
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTypes.hh"
#include "G4EmProcessSubType.hh"
#include "G4ScintillationTrackInformation.hh"

#include "LXeScintillation.h"

LXeScintillation::LXeScintillation(const G4String& processName,
                                       G4ProcessType type)
                  : G4VRestDiscreteProcess(processName, type) ,
    fTrackSecondariesFirst(false),
    fFiniteRiseTime(false),
    fYieldFactor(1.0),
    fExcitationRatio(1.0),
    fScintillationByParticleType(false),
    fScintillationTrackInfo(false),
    fStackingFlag(true),
    fNumPhotons(0),
    fEmSaturation(nullptr)
{
        SetProcessSubType(fScintillation);

        fFastIntegralTable = nullptr;
        fSlowIntegralTable = nullptr;

        if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << G4endl;
        }
}

LXeScintillation::~LXeScintillation()
{
	if (fFastIntegralTable != nullptr) {
           fFastIntegralTable->clearAndDestroy();
           delete fFastIntegralTable;
        }
        if (fSlowIntegralTable != nullptr) {
           fSlowIntegralTable->clearAndDestroy();
           delete fSlowIntegralTable;
        }
}

G4bool LXeScintillation::IsApplicable(const G4ParticleDefinition& aParticleType)
{
       return (aParticleType.GetPDGCharge() == 0.0 ||
	       aParticleType.IsShortLived()) ? false : true;
}

void LXeScintillation::BuildPhysicsTable(const G4ParticleDefinition&)
{
    if (fFastIntegralTable != nullptr) {
       fFastIntegralTable->clearAndDestroy();
       delete fFastIntegralTable;
       fFastIntegralTable = nullptr;
    }
    if (fSlowIntegralTable != nullptr) {
       fSlowIntegralTable->clearAndDestroy();
       delete fSlowIntegralTable;
       fSlowIntegralTable = nullptr;
    }
    BuildThePhysicsTable();
}


// AtRestDoIt
// ----------
//
G4VParticleChange*
LXeScintillation::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine simply calls the equivalent PostStepDoIt since all the
// necessary information resides in aStep.GetTotalEnergyDeposit()

{
        return LXeScintillation::PostStepDoIt(aTrack, aStep);
}

// PostStepDoIt
// -------------
//
G4VParticleChange*
LXeScintillation::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine is called for each tracking step of a charged particle
// in a scintillator. A Poisson/Gauss-distributed number of photons is
// generated according to the scintillation yield formula, distributed
// evenly along the track segment and uniformly into 4pi.

{
        aParticleChange.Initialize(aTrack);
        fNumPhotons = 0;

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        const G4Material* aMaterial = aTrack.GetMaterial();

        G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
        G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

        G4ThreeVector x0 = pPreStepPoint->GetPosition();
        G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
        G4double      t0 = pPreStepPoint->GetGlobalTime();

        G4double TotalEnergyDeposit = aStep.GetTotalEnergyDeposit();
        G4double ResAt511keV        = 0.06;
        G4double ResAtStepEnergy    = ResAt511keV * sqrt(511 * keV / TotalEnergyDeposit);
        G4double FluctEnergyDeposit = G4RandGauss::shoot(TotalEnergyDeposit, TotalEnergyDeposit*ResAtStepEnergy);

        G4MaterialPropertiesTable* aMaterialPropertiesTable =
                               aMaterial->GetMaterialPropertiesTable();
        if (!aMaterialPropertiesTable)
             return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

        G4MaterialPropertyVector* Fast_Intensity =
                aMaterialPropertiesTable->GetProperty(kFASTCOMPONENT);
        G4MaterialPropertyVector* Slow_Intensity =
                aMaterialPropertiesTable->GetProperty(kSLOWCOMPONENT);

        if (!Fast_Intensity && !Slow_Intensity )
             return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

        G4int nscnt = 1;
        if (Fast_Intensity && Slow_Intensity) nscnt = 2;

        G4double ScintillationYield = 0.;

	// Scintillation depends on particle type, energy deposited
        if (fScintillationByParticleType) {

           ScintillationYield =
                          GetScintillationYieldByParticleType(aTrack, aStep);

           // The default linear scintillation process
        } else {

           ScintillationYield = aMaterialPropertiesTable->
                                      GetConstProperty(kSCINTILLATIONYIELD);

           // Units: [# scintillation photons / MeV]
           ScintillationYield *= fYieldFactor;
        }

        G4double ResolutionScale    = aMaterialPropertiesTable->
                                      GetConstProperty(kRESOLUTIONSCALE);

        // Birks law saturation:

        //G4double constBirks = 0.0;

        //constBirks = aMaterial->GetIonisation()->GetBirksConstant();

        G4double MeanNumberOfPhotons;

        // Birk's correction via fEmSaturation and specifying scintillation by
        // by particle type are physically mutually exclusive

        if (fScintillationByParticleType)
           MeanNumberOfPhotons = ScintillationYield;
        else if (fEmSaturation)
           MeanNumberOfPhotons = ScintillationYield*
             (fEmSaturation->VisibleEnergyDepositionAtAStep(&aStep));
        else
            MeanNumberOfPhotons = ScintillationYield*FluctEnergyDeposit;
           // MeanNumberOfPhotons = ScintillationYield*TotalEnergyDeposit;

        if (MeanNumberOfPhotons > 10.)
        {
          G4double sigma = ResolutionScale * std::sqrt(MeanNumberOfPhotons);
          fNumPhotons=G4int(G4RandGauss::shoot(MeanNumberOfPhotons,sigma)+0.5);
        }
        else
        {
          fNumPhotons = G4int(G4Poisson(MeanNumberOfPhotons));
        }

        if ( fNumPhotons <= 0 || !fStackingFlag )
        {
           // return unchanged particle and no secondaries

           aParticleChange.SetNumberOfSecondaries(0);

           return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
        }

        ////////////////////////////////////////////////////////////////

        aParticleChange.SetNumberOfSecondaries(fNumPhotons);

        if (fTrackSecondariesFirst) {
           if (aTrack.GetTrackStatus() == fAlive )
                  aParticleChange.ProposeTrackStatus(fSuspend);
        }

        ////////////////////////////////////////////////////////////////

        G4int materialIndex = aMaterial->GetIndex();

        // Retrieve the Scintillation Integral for this material
        // new G4PhysicsOrderedFreeVector allocated to hold CII's

        G4int Num = fNumPhotons;

        for (G4int scnt = 1; scnt <= nscnt; scnt++) {

            G4double ScintillationTime = 0.*ns;
            G4double ScintillationRiseTime = 0.*ns;
            G4PhysicsOrderedFreeVector* ScintillationIntegral = nullptr;
            G4ScintillationType ScintillationType = Slow;

            if (scnt == 1) {
               if (nscnt == 1) {
                 if (Fast_Intensity) {
                   ScintillationTime   = aMaterialPropertiesTable->
                                           GetConstProperty(kFASTTIMECONSTANT);
                   if (fFiniteRiseTime) {
                      ScintillationRiseTime = aMaterialPropertiesTable->
                                  GetConstProperty(kFASTSCINTILLATIONRISETIME);
                   }
                   ScintillationType = Fast;
                   ScintillationIntegral =
                   (G4PhysicsOrderedFreeVector*)
                                       ((*fFastIntegralTable)(materialIndex));
                 }
                 if (Slow_Intensity) {
                   ScintillationTime   = aMaterialPropertiesTable->
                                           GetConstProperty(kSLOWTIMECONSTANT);
                   if (fFiniteRiseTime) {
                      ScintillationRiseTime = aMaterialPropertiesTable->
                                  GetConstProperty(kSLOWSCINTILLATIONRISETIME);
                   }
                   ScintillationType = Slow;
                   ScintillationIntegral =
                   (G4PhysicsOrderedFreeVector*)
                                       ((*fSlowIntegralTable)(materialIndex));
                 }
               }
               else {
                 G4double yieldRatio = aMaterialPropertiesTable->
                                          GetConstProperty(kYIELDRATIO);
                 if ( fExcitationRatio == 1.0 || fExcitationRatio == 0.0) {
                    Num = G4int (std::min(yieldRatio,1.0) * fNumPhotons);
                 }
                 else {
                    Num = G4int (std::min(fExcitationRatio,1.0) * fNumPhotons);
                 }
                 ScintillationTime   = aMaterialPropertiesTable->
                                          GetConstProperty(kFASTTIMECONSTANT);
                 if (fFiniteRiseTime) {
                      ScintillationRiseTime = aMaterialPropertiesTable->
                                 GetConstProperty(kFASTSCINTILLATIONRISETIME);
                 }
                 ScintillationType = Fast;
                 ScintillationIntegral =
                 (G4PhysicsOrderedFreeVector*)
                                      ((*fFastIntegralTable)(materialIndex));
               }
            }
            else {
               Num = fNumPhotons - Num;
               ScintillationTime   =   aMaterialPropertiesTable->
                                          GetConstProperty(kSLOWTIMECONSTANT);
               if (fFiniteRiseTime) {
                    ScintillationRiseTime = aMaterialPropertiesTable->
                                 GetConstProperty(kSLOWSCINTILLATIONRISETIME);
               }
               ScintillationType = Slow;
               ScintillationIntegral =
               (G4PhysicsOrderedFreeVector*)
                                      ((*fSlowIntegralTable)(materialIndex));
            }

            if (!ScintillationIntegral) continue;

            // Max Scintillation Integral

            G4double CIImax = ScintillationIntegral->GetMaxValue();

            for (G4int i = 0; i < Num; i++) {

                // Determine photon energy

                G4double CIIvalue = G4UniformRand()*CIImax;
                G4double sampledEnergy =
                              ScintillationIntegral->GetEnergy(CIIvalue);

                if (verboseLevel>1) {
                   G4cout << "sampledEnergy = " << sampledEnergy << G4endl;
                   G4cout << "CIIvalue =        " << CIIvalue << G4endl;
                }

                // Generate random photon direction

                G4double cost = 1. - 2.*G4UniformRand();
                G4double sint = std::sqrt((1.-cost)*(1.+cost));

                G4double phi = twopi*G4UniformRand();
                G4double sinp = std::sin(phi);
                G4double cosp = std::cos(phi);

                G4double px = sint*cosp;
                G4double py = sint*sinp;
                G4double pz = cost;

                // Create photon momentum direction vector

                G4ParticleMomentum photonMomentum(px, py, pz);

                // Determine polarization of new photon

                G4double sx = cost*cosp;
                G4double sy = cost*sinp;
                G4double sz = -sint;

                G4ThreeVector photonPolarization(sx, sy, sz);

                G4ThreeVector perp = photonMomentum.cross(photonPolarization);

                phi = twopi*G4UniformRand();
                sinp = std::sin(phi);
                cosp = std::cos(phi);

                photonPolarization = cosp * photonPolarization + sinp * perp;

                photonPolarization = photonPolarization.unit();

                // Generate a new photon:

                G4DynamicParticle* aScintillationPhoton =
                  new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),
                                                         photonMomentum);
                aScintillationPhoton->SetPolarization
                                     (photonPolarization.x(),
                                      photonPolarization.y(),
                                      photonPolarization.z());

                aScintillationPhoton->SetKineticEnergy(sampledEnergy);

                // Generate new G4Track object:

                G4double rand;

                if (aParticle->GetDefinition()->GetPDGCharge() != 0) {
                   rand = G4UniformRand();
                } else {
                   rand = 1.0;
                }

                G4double delta = rand * aStep.GetStepLength();
                G4double deltaTime = delta / (pPreStepPoint->GetVelocity()+
                                      rand*(pPostStepPoint->GetVelocity()-
                                            pPreStepPoint->GetVelocity())/2.);

                // emission time distribution
                if (ScintillationRiseTime==0.0) {
                   deltaTime = deltaTime -
                          ScintillationTime * std::log( G4UniformRand() );
                } else {
                   deltaTime = deltaTime +
                          sample_time(ScintillationRiseTime, ScintillationTime);
                }

                G4double aSecondaryTime = t0 + deltaTime;

                G4ThreeVector aSecondaryPosition =
                                    x0 + rand * aStep.GetDeltaPosition();

                G4Track* aSecondaryTrack = new G4Track(aScintillationPhoton,
                                                       aSecondaryTime,
                                                       aSecondaryPosition);

                aSecondaryTrack->SetTouchableHandle(
                                 aStep.GetPreStepPoint()->GetTouchableHandle());
                // aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);

                aSecondaryTrack->SetParentID(aTrack.GetTrackID());

                if (fScintillationTrackInfo) aSecondaryTrack->
                   SetUserInformation(new G4ScintillationTrackInformation(ScintillationType));

                aParticleChange.AddSecondary(aSecondaryTrack);

            }
        }

        if (verboseLevel>0) {
        G4cout << "\n Exiting from LXeScintillation::DoIt -- NumberOfSecondaries = "
               << aParticleChange.GetNumberOfSecondaries() << G4endl;
        }

        return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// BuildThePhysicsTable for the scintillation process
// --------------------------------------------------
//

void LXeScintillation::BuildThePhysicsTable()
{
        if (fFastIntegralTable && fSlowIntegralTable) return;

        const G4MaterialTable* theMaterialTable =
                               G4Material::GetMaterialTable();
        G4int numOfMaterials = G4Material::GetNumberOfMaterials();

        // create new physics table

        if(!fFastIntegralTable)fFastIntegralTable =
                                            new G4PhysicsTable(numOfMaterials);
        if(!fSlowIntegralTable)fSlowIntegralTable =
                                            new G4PhysicsTable(numOfMaterials);

        // loop for materials

        for (G4int i=0 ; i < numOfMaterials; i++)
        {
                G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
					new G4PhysicsOrderedFreeVector();
                G4PhysicsOrderedFreeVector* bPhysicsOrderedFreeVector =
                                        new G4PhysicsOrderedFreeVector();

                // Retrieve vector of scintillation wavelength intensity for
                // the material from the material's optical properties table.

                G4Material* aMaterial = (*theMaterialTable)[i];

                G4MaterialPropertiesTable* aMaterialPropertiesTable =
                                aMaterial->GetMaterialPropertiesTable();

                if (aMaterialPropertiesTable) {

                   G4MaterialPropertyVector* theFastLightVector =
                   aMaterialPropertiesTable->GetProperty(kFASTCOMPONENT);

                   if (theFastLightVector) {

                      // Retrieve the first intensity point in vector
                      // of (photon energy, intensity) pairs

                      G4double currentIN = (*theFastLightVector)[0];

                      if (currentIN >= 0.0) {

                         // Create first (photon energy, Scintillation
                         // Integral pair

                         G4double currentPM = theFastLightVector->Energy(0);

                         G4double currentCII = 0.0;

                         aPhysicsOrderedFreeVector->
                                 InsertValues(currentPM , currentCII);

                         // Set previous values to current ones prior to loop

                         G4double prevPM  = currentPM;
                         G4double prevCII = currentCII;
                         G4double prevIN  = currentIN;

                         // loop over all (photon energy, intensity)
                         // pairs stored for this material

                         for (size_t ii = 1;
                              ii < theFastLightVector->GetVectorLength();
                              ++ii)
                         {
                                currentPM = theFastLightVector->Energy(ii);
                                currentIN = (*theFastLightVector)[ii];

                                currentCII = 0.5 * (prevIN + currentIN);

                                currentCII = prevCII +
                                             (currentPM - prevPM) * currentCII;

                                aPhysicsOrderedFreeVector->
                                    InsertValues(currentPM, currentCII);

                                prevPM  = currentPM;
                                prevCII = currentCII;
                                prevIN  = currentIN;
                         }

                      }
                   }

                   G4MaterialPropertyVector* theSlowLightVector =
                   aMaterialPropertiesTable->GetProperty(kSLOWCOMPONENT);

                   if (theSlowLightVector) {

                      // Retrieve the first intensity point in vector
                      // of (photon energy, intensity) pairs

                      G4double currentIN = (*theSlowLightVector)[0];

                      if (currentIN >= 0.0) {

                         // Create first (photon energy, Scintillation
                         // Integral pair

                         G4double currentPM = theSlowLightVector->Energy(0);

                         G4double currentCII = 0.0;

                         bPhysicsOrderedFreeVector->
                                 InsertValues(currentPM , currentCII);

                         // Set previous values to current ones prior to loop

                         G4double prevPM  = currentPM;
                         G4double prevCII = currentCII;
                         G4double prevIN  = currentIN;

                         // loop over all (photon energy, intensity)
                         // pairs stored for this material

                         for (size_t ii = 1;
                              ii < theSlowLightVector->GetVectorLength();
                              ++ii)
                         {
                                currentPM = theSlowLightVector->Energy(ii);
                                currentIN = (*theSlowLightVector)[ii];

                                currentCII = 0.5 * (prevIN + currentIN);

                                currentCII = prevCII +
                                             (currentPM - prevPM) * currentCII;

                                bPhysicsOrderedFreeVector->
                                    InsertValues(currentPM, currentCII);

                                prevPM  = currentPM;
                                prevCII = currentCII;
                                prevIN  = currentIN;
                         }

                      }
                   }
                }

        // The scintillation integral(s) for a given material
        // will be inserted in the table(s) according to the
        // position of the material in the material table.

        fFastIntegralTable->insertAt(i,aPhysicsOrderedFreeVector);
        fSlowIntegralTable->insertAt(i,bPhysicsOrderedFreeVector);

        }
}

void LXeScintillation::SetScintillationByParticleType(const G4bool scintType)
{
        if (fEmSaturation) {
           G4Exception("LXeScintillation::SetScintillationByParticleType", "Scint02",
                       JustWarning, "Redefinition: Birks Saturation is replaced by ScintillationByParticleType!");
           RemoveSaturation();
        }
        fScintillationByParticleType = scintType;
}

// GetMeanFreePath
// ---------------
//

G4double LXeScintillation::GetMeanFreePath(const G4Track&,
                                          G4double ,
                                          G4ForceCondition* condition)
{
        *condition = StronglyForced;

        return DBL_MAX;

}

// GetMeanLifeTime
// ---------------
//

G4double LXeScintillation::GetMeanLifeTime(const G4Track&,
                                          G4ForceCondition* condition)
{
        *condition = Forced;

        return DBL_MAX;

}

G4double LXeScintillation::sample_time(G4double tau1, G4double tau2)
{
// tau1: rise time and tau2: decay time

	// Loop checking, 07-Aug-2015, Vladimir Ivanchenko
        while(1) {
          // two random numbers
          G4double ran1 = G4UniformRand();
          G4double ran2 = G4UniformRand();
          //
          // exponential distribution as envelope function: very efficient
          //
          G4double d = (tau1+tau2)/tau2;
          // make sure the envelope function is
          // always larger than the bi-exponential
          G4double t = -1.0*tau2*std::log(1-ran1);
          G4double gg = d*single_exp(t,tau2);
          if (ran2 <= bi_exp(t,tau1,tau2)/gg) return t;
        }
        return -1.0;
}

G4double LXeScintillation::
GetScintillationYieldByParticleType(const G4Track &aTrack, const G4Step &aStep)
{
  ////////////////////////////////////////
  // Get the scintillation yield vector //
  ////////////////////////////////////////

  G4ParticleDefinition *pDef = aTrack.GetDynamicParticle()->GetDefinition();

  G4MaterialPropertyVector *Scint_Yield_Vector = nullptr;

  G4MaterialPropertiesTable *aMaterialPropertiesTable
    = aTrack.GetMaterial()->GetMaterialPropertiesTable();

  // Get the G4MaterialPropertyVector containing the scintillation
  // yield as a function of the energy deposited and particle type

  // Protons
  if(pDef==G4Proton::ProtonDefinition())
    Scint_Yield_Vector = aMaterialPropertiesTable->
      GetProperty(kPROTONSCINTILLATIONYIELD);

  // Deuterons
  else if(pDef==G4Deuteron::DeuteronDefinition())
    Scint_Yield_Vector = aMaterialPropertiesTable->
      GetProperty(kDEUTERONSCINTILLATIONYIELD);

  // Tritons
  else if(pDef==G4Triton::TritonDefinition())
    Scint_Yield_Vector = aMaterialPropertiesTable->
      GetProperty(kTRITONSCINTILLATIONYIELD);

  // Alphas
  else if(pDef==G4Alpha::AlphaDefinition())
    Scint_Yield_Vector = aMaterialPropertiesTable->
      GetProperty(kALPHASCINTILLATIONYIELD);

  // Ions (particles derived from G4VIon and G4Ions) and recoil ions
  // below the production cut from neutrons after hElastic
  else if(pDef->GetParticleType()== "nucleus" ||
	  pDef==G4Neutron::NeutronDefinition())
    Scint_Yield_Vector = aMaterialPropertiesTable->
      GetProperty(kIONSCINTILLATIONYIELD);

  // Electrons (must also account for shell-binding energy
  // attributed to gamma from standard photoelectric effect)
  else if(pDef==G4Electron::ElectronDefinition() ||
	  pDef==G4Gamma::GammaDefinition())
    Scint_Yield_Vector = aMaterialPropertiesTable->
      GetProperty(kELECTRONSCINTILLATIONYIELD);

  // Default for particles not enumerated/listed above
  else
    Scint_Yield_Vector = aMaterialPropertiesTable->
      GetProperty(kELECTRONSCINTILLATIONYIELD);

  // If the user has specified none of the above particles then the
  // default is the electron scintillation yield
  if(!Scint_Yield_Vector)
    Scint_Yield_Vector = aMaterialPropertiesTable->
      GetProperty(kELECTRONSCINTILLATIONYIELD);

  // Throw an exception if no scintillation yield vector is found
  if (!Scint_Yield_Vector) {
    G4ExceptionDescription ed;
    ed << "\nG4Scintillation::PostStepDoIt(): "
       << "Request for scintillation yield for energy deposit and particle\n"
       << "type without correct entry in MaterialPropertiesTable.\n"
       << "ScintillationByParticleType requires at minimum that \n"
       << "ELECTRONSCINTILLATIONYIELD is set by the user\n"
       << G4endl;
    G4String comments = "Missing MaterialPropertiesTable entry - No correct entry in MaterialPropertiesTable";
    G4Exception("LXeScintillation::PostStepDoIt","Scint01",
		FatalException,ed,comments);
  }

  ///////////////////////////////////////
  // Calculate the scintillation light //
  ///////////////////////////////////////
  // To account for potential nonlinearity and scintillation photon
  // density along the track, light (L) is produced according to:
  //
  // L_currentStep = L(PreStepKE) - L(PreStepKE - EDep)

  G4double ScintillationYield = 0.;

  G4double ResAt511keV        = 0.06;
  G4double StepEnergyDeposit  = aStep.GetTotalEnergyDeposit();
  G4double ResAtStepEnergy    = ResAt511keV * sqrt(511 * keV / StepEnergyDeposit);
  G4double FluctEnergyDeposit = G4RandGauss::shoot(StepEnergyDeposit, StepEnergyDeposit*ResAtStepEnergy);

  G4double PreStepKineticEnergy = aStep.GetPreStepPoint()->GetKineticEnergy();

  if(PreStepKineticEnergy <= Scint_Yield_Vector->GetMaxEnergy()){
    G4double Yield1 = Scint_Yield_Vector->Value(PreStepKineticEnergy);
    G4double Yield2 = Scint_Yield_Vector->
                               Value(PreStepKineticEnergy - FluctEnergyDeposit); // Value(PreStepKineticEnergy - StepEnergyDeposit);
    G4cout << Yield2 << ", " << Scint_Yield_Vector->Value(PreStepKineticEnergy - StepEnergyDeposit) << G4endl;
    ScintillationYield = Yield1 - Yield2;
  } else {
    G4ExceptionDescription ed;
    ed << "\nG4Scintillation::GetScintillationYieldByParticleType(): Request\n"
       <<   "for scintillation light yield above the available energy range\n"
       <<   "specifed in G4MaterialPropertiesTable. A linear interpolation\n"
       <<   "will be performed to compute the scintillation light yield using\n"
       <<   "(L_max / E_max) as the photon yield per unit energy."
       << G4endl;
    G4String cmt = "\nScintillation yield may be unphysical!\n";
    G4Exception("LXeScintillation::GetScintillationYieldByParticleType()",
		"Scint03", JustWarning, ed, cmt);

    G4double LinearYield = Scint_Yield_Vector->GetMaxValue()
                                         / Scint_Yield_Vector->GetMaxEnergy();

    // Units: [# scintillation photons]
    ScintillationYield = LinearYield * StepEnergyDeposit;
  }

  return ScintillationYield;
}

void LXeScintillation::DumpPhysicsTable() const
{
        if (fFastIntegralTable) {
           G4int PhysicsTableSize = fFastIntegralTable->entries();
           G4PhysicsOrderedFreeVector *v;

           for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
           {
        	v = (G4PhysicsOrderedFreeVector*)(*fFastIntegralTable)[i];
        	v->DumpValues();
           }
         }

        if (fSlowIntegralTable) {
           G4int PhysicsTableSize = fSlowIntegralTable->entries();
           G4PhysicsOrderedFreeVector *v;

           for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
           {
                v = (G4PhysicsOrderedFreeVector*)(*fSlowIntegralTable)[i];
                v->DumpValues();
           }
         }
}
