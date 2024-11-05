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

#include "G4Channeling.hh"
#include "G4ChannelingMessenger.hh"
#include "Randomize.hh"
 
#include "G4ChannelingTrackData.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"

#include "G4LambdacPlus.hh"
#include "G4PhysicsModelCatalog.hh"
#include "ExExChParticleUserInfo.hh"

G4Channeling::G4Channeling():
G4VDiscreteProcess("channeling"),
fChannelingID(G4PhysicsModelCatalog::GetModelID("model_channeling")),
fMinimumEnergy(1.*CLHEP::keV), fMaximumMomentumRatio(0.01), 
fTimeStepMin(0.), fTimeStepMax(0.), fTransverseVariationMax(2.E-2 * CLHEP::angstrom),
enable_rechanneling(0), v11_algo(0), org_step_size(0), use_step_size(0), print_debug(0),
dynamic_step_size(0), density_limit(0), step_size_value(1.*nm), mfp_value(1.*nm), chan_step(0.1),
fChMessenger(0)
{
	fChMessenger = new G4ChannelingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Channeling::~G4Channeling(){
	delete fChMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ChannelingTrackData* G4Channeling::GetTrackData(const G4Track& aTrack){
	G4ChannelingTrackData* trackdata = (G4ChannelingTrackData*)(aTrack.GetAuxiliaryTrackInformation(fChannelingID));
	if (trackdata == nullptr) {
		trackdata = new G4ChannelingTrackData();
		aTrack.SetAuxiliaryTrackInformation(fChannelingID,trackdata);
	}
	return trackdata;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Channeling::GetEF(const G4Track& aTrack, G4ThreeVector& pos, G4ThreeVector& out){
	out = G4ThreeVector((GetMatData(aTrack)->GetEFX()->GetEC(pos)),
						(GetMatData(aTrack)->GetEFY()->GetEC(pos)),
						0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Channeling::GetEF(G4ChannelingMaterialData* materialData, G4ThreeVector& pos, G4ThreeVector& out){
	out = G4ThreeVector((materialData->GetEFX()->GetEC(pos)),
						(materialData->GetEFY()->GetEC(pos)),
						0.);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Channeling::PosToLattice(G4StepPoint* step,G4ThreeVector& pos)
{
	G4TouchableHistory* theTouchable = (G4TouchableHistory*)(step->GetTouchable());
	pos -= theTouchable->GetTranslation();
	pos = ((*theTouchable->GetRotation()).inverse())(pos);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4Channeling::UpdateParameters(const G4Track& aTrack)
{
	G4ChannelingMaterialData* matData = GetMatData(aTrack);
	G4LogicalCrystalVolume* aLCV = (G4LogicalCrystalVolume*)(aTrack.GetVolume()->GetLogicalVolume());
	
	G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();
	G4StepPoint* preStepPoint = aTrack.GetStep()->GetPreStepPoint();

	G4ThreeVector posPre = preStepPoint->GetPosition();
	G4ThreeVector posPost = postStepPoint->GetPosition(); // position in world

	aLCV->RotateToLattice(posPre);
	aLCV->RotateToLattice(posPost);

	G4double integrationLimit = std::fabs(posPost.z() - posPre.z());

	if (integrationLimit > 0.) {
		if (print_debug)
			G4cout << " integration limit = " << integrationLimit/nm << "[nm]" << G4endl;
		//----------------------------------------
		// Initialize particle variables
		//---------
		G4double beta = aTrack.GetVelocity()/CLHEP::c_light;
		G4double mass = aTrack.GetDynamicParticle()->GetDefinition()->GetPDGMass();
		//----------------------------------------
		// Get the momentum in the world reference 
		// frame and rotate to the solid reference frame
		//----------------------------------------
		G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
		G4ThreeVector momWorld = aTrack.GetStep()->GetPreStepPoint()->GetMomentum();
		G4ThreeVector mom = (*theTouchable->GetRotation())(momWorld);
		// randomly generates a position of the particle to evaluate the transverse energy
		// if transverse energy (potential + kinetic energy) is lower than maximum potential
		G4ThreeVector preChanPos;

		// check if no channeling position is present
		if (GetTrackData(aTrack)->GetPosCh().x() == DBL_MAX){
		// check whether it's the first step in crystal, then randomly generate position for potential evaluation
			if (enable_rechanneling == 1) {
				if (GetInfo(aTrack)->GetLastChannelingPosition().x() == DBL_MAX) {
					G4double positX = G4UniformRand() * matData->GetPot()->GetIntSp(0);
					G4double positY = G4UniformRand() * matData->GetPot()->GetIntSp(1);
					preChanPos = G4ThreeVector(positX, positY, 0.);
				} else {
					// if the particle was in channeling, use the new position of step in World reference system minus the last World position channeled, 
					// and add the difference to the last know channeling position to set the particle in the crystal reference system
					G4ThreeVector chanPos = GetInfo(aTrack)->GetLastChannelingPosition();
					G4ThreeVector newWorldPos = posPre - GetInfo(aTrack)->GetLastChannelingWorldPositionPost();
					preChanPos = chanPos + newWorldPos;
				}
			} else {
				// rechanneling disabled, generate new position in channel
				G4double positX = G4UniformRand() * matData->GetPot()->GetIntSp(0);
				G4double positY = G4UniformRand() * matData->GetPot()->GetIntSp(1);
				preChanPos = G4ThreeVector(positX, positY, 0.);
			}
		} else
			preChanPos = GetTrackData(aTrack)->GetPosCh();

		G4double mom_xy = mom.x() * mom.x() + mom.y() * mom.y();
		G4double transverse_potential_en = GetMatData(aTrack)->GetPot()->GetEC(preChanPos) * GetEffectiveCharge(aTrack);
		G4double transverse_kinetic_en = mom_xy / (2 * mass);
		G4double total_transverse_energy = transverse_kinetic_en + transverse_potential_en;

		// check for channeling condition
		// maybe also check for distance to the string?
		G4double max_potential = (GetMatData(aTrack)->GetPot()->GetMax() * GetEffectiveCharge(aTrack));
		if (total_transverse_energy > max_potential) {
			if (print_debug) {
				G4cout << " *********** " << G4endl;
				G4cout << " Dechanneling occured! Information: " << G4endl;
				G4cout << " Depth = " << posPre.z()/nm << " [nm], kinetic energy = " << aTrack.GetKineticEnergy() /keV << " [keV] " << G4endl;
				G4cout << " Crystal potential = " << transverse_potential_en /eV << " [eV], transverse kinetic energy = " << transverse_kinetic_en /eV << " [eV] " << G4endl;
			}
			GetTrackData(aTrack)->Reset();
			GetInfo(aTrack)->Reset();
			return false;
		}

		//----------------------------------------
		// Get the momentum in the solid reference 
		// frame and rotate to the crystal reference frame
		//----------------------------------------
		aLCV->RotateToLattice(mom);
		//----------------------------------------
		// Take the position stored in the track data.
		// If the particle enters the crystal,
		// the position in the channel is randomly
		// generated using a uniform distribution
		//----------------------------------------

		G4ThreeVector pos = preChanPos;
		G4double step = 0., stepTot = 0.;
		G4double nud = 0., eld = 0.;
		G4double efx = 0., efy = 0.;
		G4double nud_temp = 0., eld_temp = 0.;
		G4double Z = GetParticleDefinition(aTrack)->GetPDGCharge();
		
		const G4double oneSixth = 1. / 6.;
		G4ThreeVector posk1,posk2,posk3,posk4,posk5,posk6;
		G4ThreeVector momk1,momk2,momk3,momk4,momk5,momk6;
		G4ThreeVector pos_temp, efxy;

		do {
			//----------------------------------------
			// Limit the variable step length for the
			// integration via the selected algorithm
			// and update variables for the integration
			//----------------------------------------
			UpdateIntegrationStep(matData,aTrack,mom,step);

			if (step + stepTot > integrationLimit) {
				step = integrationLimit - stepTot;
			}

			//----------------------------------------
			// Function integration algorithm
			// 4th Order Runge-Kutta
			//----------------------------------------
			if (print_debug) {
				G4cout << " step = " << step/nm << " nm " << G4endl;
			}
			GetEF(matData, pos, efxy);
			posk1 = step / mom.z() * mom;
			if(v11_algo)
				momk1 = step / beta * Z * efxy;
			else
			momk1 = step / beta * Z * efxy * 0.5;

			GetEF(matData,pos_temp = pos + posk1 * 0.5,efxy);
			posk2 = step / mom.z() * (mom + momk1 * 0.5);
			if(v11_algo)
				momk2 = step / beta * Z * efxy;
			else
				momk2 = step / beta * Z * efxy * 0.5;

			GetEF(matData,pos_temp = pos + posk2 * 0.5,efxy);
			posk3 = step / mom.z() * (mom + momk2 * 0.5);
			if(v11_algo)
				momk3 = step / beta * Z * efxy;
			else
				momk3 = step / beta * Z * efxy * 0.5;
			
			GetEF(matData,pos_temp = pos + posk3,efxy);
			posk4 = step / mom.z() * (mom + momk3);
			if(v11_algo)
				momk4 = step / beta * Z * efxy;
			else
				momk4 = step / beta * Z * efxy * 0.5;

			pos = pos + oneSixth * (posk1 + 2.*posk2 + 2.*posk3 + posk4);
			mom = mom + oneSixth * (momk1 + 2.*momk2 + 2.*momk3 + momk4);
			//----------------------------------------
			// Update the total step and the electron
			// and nuclei density experienced by
			// the particle during its motion
			//----------------------------------------

			stepTot += step;
			nud_temp = matData->GetNuD()->GetEC(pos);
			eld_temp = matData->GetElD()->GetEC(pos);

			if(nud_temp < 0.) {nud_temp = 0.;}
			if(eld_temp < 0.) {eld_temp = 0.;}

			nud += (step * nud_temp);
			eld += (step * eld_temp);

			efx += (step * matData->GetEFX()->GetEC(pos));
			efy += (step * matData->GetEFY()->GetEC(pos));
			//----------------------------------------
		} while (stepTot < integrationLimit);
 
		nud /= stepTot;
		eld /= stepTot;

		if(nud < density_limit)
			nud = density_limit;

		if(eld < density_limit)
			eld = density_limit;

		GetTrackData(aTrack)->SetNuD(nud);
		GetTrackData(aTrack)->SetElD(eld);

		GetInfo(aTrack)->SetNucleiDensity(nud);
		GetInfo(aTrack)->SetElectronDensity(eld);

		GetTrackData(aTrack)->SetEFX(efx);
		GetTrackData(aTrack)->SetEFY(efy);
		
		//----------------------------------------
		// Scattering on electrons
		//----------------------------------------
		
		G4double gamma = aTrack.GetTotalEnergy()/mass;
		G4double m_e = CLHEP::electron_mass_c2;

		G4double T = aTrack.GetStep()->GetTotalEnergyDeposit();
		T -= aTrack.GetStep()->GetNonIonizingEnergyDeposit();
		// there's no such thing in geant4 github, taken from one of Bagli's source codes
		if (T > 0.){
			G4double Tmax = 2. * m_e * gamma * gamma * beta * beta;
			Tmax /= (1. + 2. * gamma * m_e / mass + m_e * m_e / mass / mass);
			G4double mmod = aTrack.GetStep()->GetPreStepPoint()->GetMomentum().mag();
			G4double theta_sc_el = sqrt(2. * m_e * T * (1. -T / Tmax)) / mmod;
			if (theta_sc_el > 0.) {
				G4double rot_angle = G4UniformRand() * CLHEP::twopi;
				G4double rot_mod   = G4RandGauss::shoot(0., theta_sc_el * eld);
				mom.rotateX(cos(rot_angle) * rot_mod);
				mom.rotateY(sin(rot_angle) * rot_mod);
			}
		}

		GetTrackData(aTrack)->SetMomCh(mom);
		GetTrackData(aTrack)->SetPosCh(pos);

		GetInfo(aTrack)->SetMomentumChanneled(mom);
		GetInfo(aTrack)->SetPositionChanneled(pos);

		GetInfo(aTrack)->SetLastChannelingPosition(pos);
		GetInfo(aTrack)->SetLastChannelingMomentum(mom);

		GetInfo(aTrack)->SetLastChannelingWorldPositionPost(posPost);
		GetInfo(aTrack)->SetLastChannelingWorldPositionPre(posPre);
		GetInfo(aTrack)->SetLastChannelingWorldMomentum(momWorld);

		return true;
	}
	return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4Channeling::
UpdateIntegrationStep(G4ChannelingMaterialData* matData, const G4Track& aTrack, G4ThreeVector& mom, G4double& step)
{
	if (use_step_size) {
		step = step_size_value;
		if (step < fTimeStepMin)
			step = fTimeStepMin;
		return true;
	} else {
		if (mom.x() != 0.0 || mom.y() != 0.0) {
			G4double xy2 = mom.x() * mom.x() + mom.y() * mom.y();
			if (xy2 != 0.) {
				step = std::fabs(fTransverseVariationMax * GetPre(aTrack)->GetKineticEnergy() / std::pow(xy2, 0.5));
				if (step < fTimeStepMin)
					step = fTimeStepMin;
				else {
					fTimeStepMax = std::sqrt( fTransverseVariationMax * GetPre(aTrack)->GetKineticEnergy() / std::fabs(matData->GetEFX()->GetMax()));
					if (step > fTimeStepMax)
						step = fTimeStepMax;
				}
			} else {
				step = fTimeStepMin;
			}
			return true;
		} else {
			step = fTimeStepMin;
		}
		return false;
	}
	return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4Channeling::
GetMeanFreePath(const G4Track& aTrack, G4double /*previousStepSize */, G4ForceCondition* condition){
	//----------------------------------------
	// the condition is forced to check if
	// the volume has a lattice at each step.
	// if it hasn't, return DBL_MAX
	//----------------------------------------
	*condition = Forced;
	G4LogicalVolume* aLV = aTrack.GetVolume()->GetLogicalVolume();
	G4LogicalVolume* aNLV = aTrack.GetNextVolume()->GetLogicalVolume();

	if (G4LogicalCrystalVolume::IsLattice(aLV) == true && G4LogicalCrystalVolume::IsLattice(aNLV) == true && aTrack.GetKineticEnergy() > fMinimumEnergy) {
		GetInfo(aTrack)->SetInTheCrystal(true);
		
		if (use_step_size) {
			fTimeStepMin = step_size_value;
			if (print_debug)
				G4cout << " mfp = " << mfp_value/nm << " [nm], timestepmin = " << fTimeStepMin/nm << " [nm] " << G4endl;
			return mfp_value;
		}

		if (org_step_size) {
			G4double osc_per = GetOscillationPeriod(aTrack);
			fTimeStepMin = osc_per * 2.e-4;
			if (print_debug)
				G4cout << " mfp = " << osc_per * 0.01 << " [nm], timestepmin = " << fTimeStepMin/nm << " [nm] " << G4endl;
			return osc_per * 0.01;
		} else {
			G4double osc_per2 = GetOscillationPeriod2(aTrack);
			fTimeStepMin = osc_per2 * chan_step * 0.01;
			if (print_debug)
				G4cout << " mfp = " << osc_per2 * chan_step << " [nm], timestepmin = " << fTimeStepMin/nm << " [nm] " << G4endl;
			return osc_per2 * chan_step;
		}
	} else {
		GetTrackData(aTrack)->Reset();
		GetInfo(aTrack)->Reset();
		GetInfo(aTrack)->SetInTheCrystal(false);
		return DBL_MAX;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4Channeling::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
	//----------------------------------------
	// check if the volume has a lattice
	// and if the particle is in channeling.
	// If it is so, the particle is forced
	// to follow the channeling plane
	// direction. If the particle has
	// dechanneled or exited the crystal,
	// the outgoing angle is evaluated
	//----------------------------------------

	aParticleChange.Initialize(aTrack);
	GetInfo(aTrack)->StoreDensityPreviousStep();
	G4LogicalVolume* aLV = aTrack.GetVolume()->GetLogicalVolume();
	G4LogicalVolume* aNLV = aTrack.GetNextVolume()->GetLogicalVolume();

	if (G4LogicalCrystalVolume::IsLattice(aLV) == true && G4LogicalCrystalVolume::IsLattice(aNLV) == true) {
		G4bool bModifiedTraj = UpdateParameters(aTrack);
		if (bModifiedTraj == true) {
			//----------------------------------------
			// Get the momentum in the reference frame
			// solidal to the bent planes and rotate
			// to the reference frame
			//----------------------------------------
			G4LogicalCrystalVolume* aLCV = (G4LogicalCrystalVolume*)(aTrack.GetVolume()->GetLogicalVolume());
			G4ThreeVector momCh = GetTrackData(aTrack)->GetMomCh();
			G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();
			G4TouchableHistory* theTouchable = (G4TouchableHistory*)(postStepPoint->GetTouchable());
			//----------------------------------------
			// Get the momentum in the crystal reference
			// frame and rotate to the solid reference frame
			//----------------------------------------
			aLCV->RotateToSolid(momCh);
			//----------------------------------------
			// Get the momentum in the solid reference
			// frame and rotate to the world reference frame
			//----------------------------------------
			G4ThreeVector mom = ((*theTouchable->GetRotation()).inverse())(momCh);
			aParticleChange.ProposeMomentumDirection(mom.unit());
			aParticleChange.ProposePolarization(GetTrackData(aTrack)->GetSpinCh());
		}
	} else {
		// if the volume has no lattice it resets the density factors
		GetTrackData(aTrack)->Reset();
		GetInfo(aTrack)->Reset();
	}
	
	return &aParticleChange;
}
ExExChParticleUserInfo* G4Channeling::GetInfo(const G4Track& aTrack)
{
	return (ExExChParticleUserInfo*) aTrack.GetUserInformation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
