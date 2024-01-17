//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//

#ifndef G4Channeling_h
#define G4Channeling_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ChannelingMaterialData.hh"
#include "G4ExtendedMaterial.hh"
#include "G4LogicalCrystalVolume.hh"
#include "G4RunManager.hh"
#include "ExExChParticleUserInfo.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4hIonEffChargeSquare.hh"

class G4ChannelingTrackData;
class G4ChannelingMessenger;

class G4Channeling : public G4VDiscreteProcess
{
	public:
		G4Channeling();
		virtual ~G4Channeling();
		virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);
		virtual G4bool IsApplicable(const G4ParticleDefinition& aPD){
			return(aPD.GetPDGCharge() != 0.);
		};
		virtual void BuildPhysicsTable(const G4ParticleDefinition&){;};
		
	protected:
		virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

	private:
		G4hIonEffChargeSquare* eff_charge = new G4hIonEffChargeSquare("");
		G4ParticleDefinition* GetParticleDefinition(const G4Track& aTrack)
			{return const_cast<G4ParticleDefinition*>(aTrack.GetParticleDefinition());}
		ExExChParticleUserInfo* GetInfo(const G4Track&);	

	private:
		G4StepPoint* GetPre(const G4Track& aTrack){return aTrack.GetStep()->GetPreStepPoint();}
		G4StepPoint* GetPost(const G4Track& aTrack){return aTrack.GetStep()->GetPostStepPoint();}

	private:
		G4ChannelingMaterialData* GetMatData(const G4Track& aTrack) {
			G4LogicalVolume* aLV = aTrack.GetVolume()->GetLogicalVolume();
			if (aLV->IsExtended() == true) {
				G4ExtendedMaterial* aEM = (G4ExtendedMaterial*) aTrack.GetVolume()->GetLogicalVolume()->GetMaterial();
				return (G4ChannelingMaterialData*) aEM->RetrieveExtension("channeling");
			} else {
				return nullptr;
			}
		}

		G4double GetEffectiveCharge(const G4Track& aTrack) {
			const G4ParticleDefinition* particle = aTrack.GetParticleDefinition();
			G4Material* material = aTrack.GetMaterial();
			G4double kinetic_energy = aTrack.GetKineticEnergy();
			return sqrt(eff_charge->TheValue(particle,material,kinetic_energy));
		}
		//----------------------------------------
		// Functions for the calculations of
		// parameters related to channeling
		//----------------------------------------
	public:
		G4double GetCriticalAngle2(const G4Track& aTrack)
			{return (std::sqrt(GetMatData(aTrack)->GetPot()->GetMaxMin() / GetPre(aTrack)->GetKineticEnergy()));}

		G4double GetCriticalAngle4(const G4Track& aTrack) {
			G4double pot_diff = GetMatData(aTrack)->GetPot()->GetMax() - abs(GetMatData(aTrack)->GetPot()->GetMin());
			return (std::sqrt(pot_diff/GetPre(aTrack)->GetKineticEnergy()));
		}

		G4double GetCriticalAngle3(const G4Track& aTrack)
			{return std::sqrt((GetMatData(aTrack)->GetPot()->GetMax() * GetEffectiveCharge(aTrack)) / GetPre(aTrack)->GetKineticEnergy());}

		G4double GetCriticalAngle(const G4Track& aTrack) {
			return std::sqrt(2.0 * GetMatData(aTrack)->GetPot()->GetMaxMin() / GetPre(aTrack)->GetTotalEnergy());} 
		
		G4double GetOscillationPeriod(const G4Track& aTrack) {
			return (CLHEP::pi * GetMatData(aTrack)->GetPot()->GetIntSp(0) / GetCriticalAngle(aTrack));}

		G4double GetOscillationPeriod2(const G4Track& aTrack) {
			return (CLHEP::pi * GetMatData(aTrack)->GetPot()->GetIntSp(0) / GetCriticalAngle2(aTrack));}

		// axial critical energy
		// from Critical Parameters of Channeling, S. I. Matyukhin
		G4double GetCriticalEnergy(const G4Track& aTrack) {
			G4double Z1 = aTrack.GetDynamicParticle()->GetCharge();
			G4double Z2 = (aTrack.GetMaterial()->GetTotNbOfElectPerVolume() / aTrack.GetMaterial()->GetTotNbOfAtomsPerVolume());
			G4double d = GetMatData(aTrack)->GetPot()->GetIntSp(0);
			G4double e2 = 1.4399764 * (MeV *fermi); // MeV*fm
			G4double a_tf = GetScreeningDistance(Z1,Z2);
			return (2 * Z1 * Z2 * e2 * d / (std::pow(a_tf, 2.)));
		}
		G4double GetScreeningDistance(G4double Z1, G4double Z2){
			G4double a_bohr = 5.291772e-11*m;
			return (0.885 * a_bohr * std::pow(std::pow(Z1, (1. / 2.)) + std::pow(Z2, (1. / 2.)), (-2. / 3.)));
		}
		// Critical angle, when energy is higher than the critical
		G4double GetCritAngleLimitAxial(const G4Track& aTrack) {
			G4double Z1 = aTrack.GetDynamicParticle()->GetCharge();
			G4double Z2 = (aTrack.GetMaterial()->GetTotNbOfElectPerVolume() / aTrack.GetMaterial()->GetTotNbOfAtomsPerVolume());
			G4double d = GetMatData(aTrack)->GetPot()->GetIntSp(0);
			G4double e2 = 1.4399764 *(MeV*fermi); // MeV*fm
			G4double kin_en = GetPre(aTrack)->GetKineticEnergy() /MeV;
			return sqrt(2 * Z1 * Z2 * e2 / (d * kin_en));
		}
		// PLANAR case
		// Critical angle, when energy is higher than the critical
		G4double GetCritAngleLimitPlanar(const G4Track& aTrack) {
			G4double Z1 = aTrack.GetDynamicParticle()->GetCharge();
			G4double Z2 = (aTrack.GetMaterial()->GetTotNbOfElectPerVolume() / aTrack.GetMaterial()->GetTotNbOfAtomsPerVolume());
			G4double d = GetMatData(aTrack)->GetPot()->GetIntSp(0) /fermi;
			G4double a_tf = GetScreeningDistance(Z1, Z2) /fermi;
			G4double density = aTrack.GetMaterial()->GetTotNbOfAtomsPerVolume() / (1 / (fermi * fermi * fermi));
			G4double surf_density = density * d;
			G4double e2 = 1.4399764; // MeV*fm *(MeV*CLHEP::fermi)
			G4double kin_en = GetPre(aTrack)->GetKineticEnergy()/MeV;
			return sqrt((2 *CLHEP::pi * Z1 * Z2 * e2 * surf_density * a_tf) / kin_en);
		}

		void SetChannelingStep(G4double a)  {chan_step = a;}
		G4double GetChannelingStep()        {return chan_step;}
		//----------------------------------------
		// Channeling Auxiliary Track Information
		//----------------------------------------
	private:
		G4int fChannelingID;
		G4ChannelingTrackData* GetTrackData(const G4Track&);
		//----------------------------------------
		// Variables for the integration
		// of the particle trajectory
		//----------------------------------------
	private:
		G4bool UpdateIntegrationStep(G4ChannelingMaterialData* materialData, const G4Track&, G4ThreeVector&, G4double&);
		G4bool UpdateParameters(const G4Track&);

		void GetEF(const G4Track&,G4ThreeVector&,G4ThreeVector&);
		void GetEF(G4ChannelingMaterialData*,G4ThreeVector&,G4ThreeVector&);
	public:
		void PosToLattice(G4StepPoint* step,G4ThreeVector&);
		
	public:
		G4double GetTransverseVariationMax()             {return fTransverseVariationMax;};
		void SetTransverseVariationMax(G4double aDouble) {fTransverseVariationMax = aDouble;};
		
		G4double GetTimeStepMin()                        {return fTimeStepMin;};
		void SetTimeStepMin(G4double aDouble)            {fTimeStepMin = aDouble;};

		G4double GetMinimumEnergy()                      {return fMinimumEnergy;};
		void SetMinimumEnergy(G4double aDouble)          {fMinimumEnergy = aDouble;};
		
		G4double GetMaximumMomentumRatio()               {return fMaximumMomentumRatio;};
		void SetMaximumMomentumRatio(G4double aDouble)   {fMaximumMomentumRatio = aDouble;};

		void EnableRechanneling(bool a)                  {enable_rechanneling = a;}

		void SetMinimumDensityLimit(G4double aDouble)    {density_limit = aDouble;}
		void UseV11ChannelingAlgo(bool a)                {v11_algo = a;}
		void UseOrgStepSize(bool a)                      {org_step_size = a;}

		void SetStepSizeUnit(G4double aDouble)           {step_size_value = aDouble;}
		void SetMFPSizeUnit(G4double aDouble)            {mfp_value = aDouble;}
		void UseStepSize(bool a)                         {use_step_size = a;}
		
	private:
		G4double fMinimumEnergy, fMaximumMomentumRatio;
		G4double fTimeStepMin, fTimeStepMax;
		G4double fTransverseVariationMax;
		
		const G4ThreeVector k010;

		G4bool enable_rechanneling, v11_algo, org_step_size, use_step_size;
		G4double density_limit, step_size_value, mfp_value, chan_step;
		G4ChannelingMessenger* fChMessenger;
};

#endif