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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 76293 2013-11-08 13:11:23Z gcosmo $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "SensitiveDetectorHit.hh"
#include "G4Material.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4EmCalculator.hh"
#include "G4NistManager.hh"
#include "RBSHelpers.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* DE, PrimaryGeneratorAction* PRA )
:G4UserEventAction(),detector(DE),fPrimary(PRA),
 fTotalEnergyDeposit(0.), fTotalEnergyFlow(0.), fTotalNiel(0), theta(0),
 prim_step(0), trak_len_prim(0), trak_len_sec(0),
 TrueTrakLen(0.), ProjTrakLen(0.), projected_range(0),
 sdht_ID(-1), sd0_ID(-1),sd1_ID(-1),sd2_ID(-1),sd3_ID(-1),sd4_ID(-1)
{
	fVectorSi_total = nullptr;
	fVectorO_total = nullptr;
	fVectorNa_total = nullptr;
	fVectorN_total = nullptr;
	fVectorC_total = nullptr;
	fVectorF_total = nullptr;
	fVectorB_total = nullptr;
	fVectorNi_total = nullptr;
	fVectorCu_total = nullptr;
	//G4UAtomicDeexcitation* deExcitation = new G4UAtomicDeexcitation();
	//deExcitation->InitialiseAtomicDeexcitation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
	for (int i=0; i<5; i++) {
		delete[] Znumb[i];
		delete[] Mnumb[i];
		delete[] Adens[i];
		delete[] HistoBase[i];
		delete[] nud_el[i];
		delete[] mat_matrix[i];
	}
	delete[] NoOfElements;
	delete[] Znumb;
	delete[] Mnumb;
	delete[] Adens;
	delete[] thicknesses;
	delete[] positions;
	delete[] layer_start_pos;
	delete[] layer_end_pos;
	delete[] dist_between_l;
	delete[] dist_matrix;
	delete[] start_pos;
	delete[] HistoBase;
	delete[] nud_el;
	//
	delete fVectorB_total;
	delete fVectorC_total;
	delete fVectorCu_total;
	delete fVectorF_total;
	delete fVectorN_total;
	delete fVectorNa_total;
	delete fVectorNi_total;
	delete fVectorO_total;
	delete fVectorSi_total;
	//
	delete ElementVector;
	delete Element;
	delete Isotope;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
	fTotalEnergyDeposit = 0.;
	fTotalEnergyFlow = 0.;
	theta = 0.;
	TrueTrakLen = 0.;
	ProjTrakLen = 0.;
	for (G4int i=0; i<5; ++i) {
		StepsPrim[i] = 0;
		step_layer[i] = 0;
		ion_dep_layer[i] = 0;
		non_dep_layer[i] = 0;
	}
	fTotalNiel = 0;
	projected_range = 0;
	trak_len_prim = 0;
	prim_step = 0;
	trak_len_sec = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
	Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

	run->AddEdep(fTotalEnergyDeposit);
	run->AddEflow(fTotalEnergyFlow);
	run->AddTheta(theta);
	run->AddTrueRange(TrueTrakLen);
	run->AddProjRange(ProjTrakLen);
	run->AddProjectedRange(projected_range);
	run->AddNonIonEnergy(fTotalNiel);
	run->AddTrakLenPrim(trak_len_prim);
	run->AddTrakLenSec(trak_len_sec);
	run->AddNumberOfSteps(prim_step);
	
	for (G4int i=0; i<5; ++i) {
		run->absStepLayer(StepsPrim[i], i);
		run->absTrackLenLayer(step_layer[i], i);
		run->absIonLayer(ion_dep_layer[i], i);
		run->absNonLayer(non_dep_layer[i], i);
	}

	G4AnalysisManager::Instance()->FillH1(1, fTotalEnergyDeposit);
	G4AnalysisManager::Instance()->FillH1(3, fTotalEnergyFlow);
	
	// ******************************
	for (int i=0; i<NUMBER_OF_MAX_LAYERS; i++) {
		Znumb[i] = new G4double[4];
		Mnumb[i] = new G4double[4];
		Adens[i] = new G4double[4];
		HistoBase[i] = new G4int[4];
		nud_el[i] = new G4double[3];
	}

	for (int i=0; i<NUMBER_OF_MAX_LAYERS; ++i) {
		for (int j=0; j<3; ++j) {
			nud_el[i][j] = 0;
		}
	}

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	//detector energy resolution
	G4int rbs_eval = detector->GetRBSCalc();
	G4ParticleDefinition* fParticle = fPrimary->GetParticleGun()->GetParticleDefinition();
	G4double primary_energy = fPrimary->GetParticleGun()->GetParticleEnergy();
	run->PrimaryEnergy(primary_energy);
	// G4int track_histos = detector->GetTrackingHistos();
	//**************************************************************************
	// Incidence angle
	G4double angle_of_incidence = atan(fPrimary->GetParticleGun()->GetParticleMomentumDirection().x() / fPrimary->GetParticleGun()->GetParticleMomentumDirection().z());
	run->Inc_angle(angle_of_incidence);

	PrepareHistoBase();
	// prepare depth intervals
	detector->CopyArrayForDepth(depth_for_histos);
	
	//scattering angle
	Angle = detector->GetRBSAngle();
	angle_of_exit = pi - Angle;
	// G4EmCalculator emCalculator;
	G4String dead_material_name = detector->GetDeadLayer();
	// G4ThreeVector ch_pos;

	for (int i=0; i<NUMBER_OF_MAX_LAYERS; i++) {
		sample_material[i] = detector->GetMaterialM(i);
		NoOfElements[i]    = sample_material[i]->GetNumberOfElements();
	}
	G4int max_M = 0;

	for (int i=0; i<NUMBER_OF_MAX_LAYERS; i++) {
		for (int j = 0; j<NoOfElements[i]; j++) {
			Znumb[i][j] = sample_material[i]->GetElement(j)->GetZ();
			Mnumb[i][j] = sample_material[i]->GetElement(j)->GetA() /(g/mole);
			if (Mnumb[i][j] > max_M)
				max_M = Mnumb[i][j];
			const G4double *atomDensVector	= sample_material[i]->GetVecNbOfAtomsPerVolume();
			Adens[i][j] = atomDensVector[j] /(1/cm3);
		}
	}

	// RTR vector fillup
	if (detector->GetSigmaCalc() && rbs_eval) {
		FillSigmaCalcVectors();
	}

	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	
	G4ThreeVector ssd[4];
	ssd[0]= G4ThreeVector(0.,0.,0.);
	ssd[1]= G4ThreeVector(0.,0.,0.);
	ssd[2]= G4ThreeVector(0.,0.,0.);
	ssd[3]= G4ThreeVector(0.,0.,0.);

	G4double energy = 0., kinen = 0;

	if (sdht_ID == -1) {
		G4String sdName;
		if (SDman->FindSensitiveDetector(sdName="telescope",0)){
			sdht_ID = SDman->GetCollectionID(sdName="telescope/collection");
		}
	}

	SensitiveDetectorHitsCollection* sdht = 0;
	G4HCofThisEvent *hce = evt->GetHCofThisEvent();

	if (hce){
		if (sdht_ID != -1) {
			G4VHitsCollection* aHCSD = hce->GetHC(sdht_ID);
			sdht = (SensitiveDetectorHitsCollection*)(aHCSD);
		}
	}
	int bTotalHits = 0;
	if (sdht) {
		int n_hit_sd = sdht->entries();
		for (int i2=0; i2<4; i2++) {
			for (int i1=0; i1<n_hit_sd; i1++) {
				SensitiveDetectorHit* aHit = (*sdht)[i1];
				if (aHit->GetLayerID() == i2) {
					ssd[i2] = aHit->GetWorldPos();
					bTotalHits++; // checks, how many detectors got a hit
				}
				if (aHit->GetLayerID() == 2) { // info from the 3rd detector
					energy = aHit->GetKinE();
					analysisManager->FillH1(36, energy);
				}
				// particle hit in the last detector [no 4]
				if (aHit->GetLayerID() == 3) {
					run->AddDetHits(); // adds hits
					kinen = aHit->GetKinE(); // takes kinetinic energy
					run->AddKinEn(kinen); // add kinetic energy 
				}
			}
		}
	}

	// sensitive detectors
	G4String sdName;
	if (sd0_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector0",0))
			sd0_ID = SDman->GetCollectionID(sdName="crystaldetector0/collection");
	}
	if (sd1_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector1",0))
			sd1_ID = SDman->GetCollectionID(sdName="crystaldetector1/collection");
	}
	if (sd2_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector2",0))
			sd2_ID = SDman->GetCollectionID(sdName="crystaldetector2/collection");
	}
	if (sd3_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector3",0))
			sd3_ID = SDman->GetCollectionID(sdName="crystaldetector3/collection");
	}
	if (sd4_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector4",0))
			sd4_ID = SDman->GetCollectionID(sdName="crystaldetector4/collection");
	}

	CrystalDetectorHitsCollection* sd0 = 0;
	G4HCofThisEvent *hx0 = evt->GetHCofThisEvent();
	if (hx0) {
		if (sd0_ID != -1) {
			G4VHitsCollection* aHCSD0 = hx0->GetHC(sd0_ID);
			sd0 = (CrystalDetectorHitsCollection*)(aHCSD0);
		}
	}

	CrystalDetectorHitsCollection* sd1 = 0;
	G4HCofThisEvent *hx1 = evt->GetHCofThisEvent();
	if (hx1) {
		if(sd1_ID != -1) {
			G4VHitsCollection* aHCSD1 = hx1->GetHC(sd1_ID);
			sd1 = (CrystalDetectorHitsCollection*)(aHCSD1);
		}
	}

	CrystalDetectorHitsCollection* sd2 = 0;
	G4HCofThisEvent *hx2 = evt->GetHCofThisEvent();
	if (hx2) {
		if (sd2_ID != -1){
			G4VHitsCollection* aHCSD2 = hx2->GetHC(sd2_ID);
			sd2 = (CrystalDetectorHitsCollection*)(aHCSD2);
		}	
	}

	CrystalDetectorHitsCollection* sd3 = 0;
	G4HCofThisEvent *hx3 = evt->GetHCofThisEvent();
	if (hx3) {
		if (sd3_ID != -1) {
			G4VHitsCollection* aHCSD3 = hx3->GetHC(sd3_ID);
			sd3 = (CrystalDetectorHitsCollection*)(aHCSD3);
		}
	}

	CrystalDetectorHitsCollection* sd4 = 0;
	G4HCofThisEvent *hx4 = evt->GetHCofThisEvent();
	if (hx4) {
		if (sd4_ID != -1) {
			G4VHitsCollection* aHCSD4 = hx4->GetHC(sd4_ID);
			sd4 = (CrystalDetectorHitsCollection*)(aHCSD4);
		}
	}

	// prepare distances
	PrepareDistances();

	G4double steps = 0, sample_energy = 0, tot_step = 0;
	G4ThreeVector position, momDir;
	// chan params
	G4double efx = 0, efy = 0, eld = 0, nud = 0;
	G4double efxavg = 0, efyavg = 0, eldavg = 0, nudavg = 0;
	G4double nudavg_a = 0, nudavg_b = 0, nudavg_c = 0;

	// fill ntuples
	if (bTotalHits > 3 ) {
		G4double angXin  = ((ssd[1].x() - ssd[0].x()) / (ssd[1].z() - ssd[0].z()));
		G4double angYin  = ((ssd[1].y() - ssd[0].y()) / (ssd[1].z() - ssd[0].z()));
		analysisManager->FillNtupleDColumn(0, atan(angXin) /degree); 
		analysisManager->FillNtupleDColumn(1, atan(angYin) /degree); 
		G4double posXin = ssd[1].x() - atan(angXin) * ssd[1].z();
		G4double posYin = ssd[1].y() - atan(angYin) * ssd[1].z();
		analysisManager->FillNtupleDColumn(3, posXin / CLHEP::mm);
		analysisManager->FillNtupleDColumn(4, posYin / CLHEP::mm);
		// PARAMETERS FOR OUT
		G4double angXout = -9999.;
		G4double angYout = -9999.;
		G4double posXout = -9999.;
		G4double posYout = -9999.;
		
		if(bTotalHits == 4) {
			angXout = ((ssd[3].x() - ssd[2].x()) / (ssd[3].z() - ssd[2].z()));
			angYout = ((ssd[3].y() - ssd[2].y()) / (ssd[3].z() - ssd[2].z()));
			posXout = ssd[2].x() - atan(angXout) * ssd[2].z();
			posYout = ssd[2].y() - atan(angYout) * ssd[2].z();
		}
		analysisManager->FillNtupleDColumn(42, atan(angXout) /degree); 
		analysisManager->FillNtupleDColumn(43, atan(angYout) /degree); 
		analysisManager->FillNtupleDColumn(40, posXout / CLHEP::mm); 
		analysisManager->FillNtupleDColumn(41, posYout / CLHEP::mm);
	}

	ProcessDetHits(sd0, 0, NoOfElements[0]);
	ProcessDetHits(sd1, 1, NoOfElements[1]);
	ProcessDetHits(sd2, 2, NoOfElements[2]);
	ProcessDetHits(sd3, 3, NoOfElements[3]);
	ProcessDetHits(sd4, 4, NoOfElements[4]);

	analysisManager->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void EventAction::ProcessDetHits(CrystalDetectorHitsCollection* sd, G4int layer, G4int elements)
{
	if (!sd)
		return;

	Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	G4int track_histos = detector->GetTrackingHistos();
	G4double RBSROImin = detector->GetRBSROImin();
	G4ParticleDefinition* fParticle = fPrimary->GetParticleGun()->GetParticleDefinition();

	G4double efx = 0, efy = 0, eld = 0, nud = 0;
	G4double efxavg = 0, efyavg = 0, eldavg = 0, nudavg = 0;
	G4double nudavg_a = 0, nudavg_b = 0, nudavg_c = 0;
	G4double tot_step = 0;
	int n_hit_sd = sd->entries();
	run->add_entry_sd(n_hit_sd);

	for (int i = 0; i < n_hit_sd; ++i) {
		CrystalDetectorHit* aHit = (*sd)[i];
		G4double steps = aHit->GetStep();
		G4ThreeVector position = aHit->GetWorldPos();
		G4double sample_energy = aHit->GetKinECR();
		G4ThreeVector ch_pos = aHit->GetChPos();
		G4double nud = aHit->GetNud();
		G4double efx = aHit->GetEFX();
		G4double efy = aHit->GetEFY();
		G4double eld = aHit->GetEld();
		G4ThreeVector momDir = aHit->GetWorldMomentumDirection();
		G4double nud_el[3];
		if (nud == 1.)
			nud_el[0] = nud_el[1] = nud_el[2] = 1;
		else {
			nud_el[0] = aHit->GetNuD_a();
			nud_el[1] = aHit->GetNuD_b();
			nud_el[2] = aHit->GetNuD_c();
		}

		tot_step += steps;
		eldavg += eld * steps;
		nudavg += nud * steps;
		nudavg_a += nud_el[0] * steps;
		nudavg_b += nud_el[1] * steps;
		nudavg_c += nud_el[2] * steps;
		efxavg += efx * steps;
		efyavg += efy * steps;

		G4double z_pos = position.z() + detector->GetLength(0) / 2;
		if (track_histos) {
			FillH2ChannelingHistos(z_pos, ch_pos.x(), ch_pos.y(), position.x(), position.y());
			analysisManager->FillP1(2 * layer + 1, z_pos, nud);
			analysisManager->FillP1(2 * layer + 2, z_pos, eld);
		}

		if (sample_energy > RBSROImin) {
			CalculateRBS(layer, sample_energy, position, momDir, steps, fParticle, elements, nud_el);
		}
	}

	if (tot_step > 0) {
		nudavg /= tot_step;
		eldavg /= tot_step;
		efxavg /= tot_step;
		efyavg /= tot_step;
		nudavg_a /= tot_step;
		nudavg_b /= tot_step;
		nudavg_c /= tot_step;
	}
	run->add_total_step(tot_step);

	G4AnalysisManager::Instance()->FillNtupleDColumn(4 * layer + 1, efxavg / CLHEP::eV * CLHEP::angstrom);
	G4AnalysisManager::Instance()->FillNtupleDColumn(4 * layer + 2, efyavg / CLHEP::eV * CLHEP::angstrom);
	G4AnalysisManager::Instance()->FillNtupleDColumn(4 * layer + 3, nudavg);
	G4AnalysisManager::Instance()->FillNtupleDColumn(4 * layer + 4, eldavg);
}

G4double EventAction::GenerateGaussian(G4double x, G4double y, G4double sigma_sq)
{
	G4double inv_sqrt_2pi = 0.398942280401433;
	G4double a = (x - y);
	return (inv_sqrt_2pi / std::sqrt(sigma_sq)) * std::exp(-0.5 * a * a / sigma_sq);
}

G4Physics2DVector* EventAction::Get2DRTRVector(G4String element, G4int A)
{
	G4Physics2DVector* vector = new G4Physics2DVector;
	G4String start_file = "RTR_values/";
	G4String end_file = "_total.txt";
	G4String folder_def, part_def;

	if (A == 1) {
		folder_def = "Hydrogen/";
		part_def = "H_";
	} else if (A == 4) {
		folder_def = "Helium/";
		part_def = "He_";
	} else
		return 0;

	G4String tot_filename = start_file + folder_def + part_def + element + end_file;
	G4cout << " Reading RTR values from file: " << tot_filename << G4endl;
	std::ifstream vFile(tot_filename);
	if (!vFile.is_open()) {
		G4cout << " No RTR file " << tot_filename << G4endl;
		exit(1);
	}
	vector->Retrieve(vFile);
	return vector;
}

G4double EventAction::Get2DRTRValue(G4double energy, G4String elname, G4double angle)
{
	size_t fIndx(0);
	size_t fIndy(0);

	G4double value = 1.;
	if(elname == "Si")
		value = fVectorSi_total->Value(angle /degree, energy /keV, fIndx, fIndy);
	else if(elname == "O")
		value = fVectorO_total->Value(angle /degree, energy /keV, fIndx, fIndy);
	else if(elname == "Na")
		value = fVectorNa_total->Value(angle /degree, energy /keV, fIndx, fIndy);
	else if(elname == "N")
		value = fVectorN_total->Value(angle /degree, energy /keV, fIndx, fIndy);
	else if(elname == "C")
		value = fVectorC_total->Value(angle /degree, energy /keV, fIndx, fIndy);
	else if(elname == "F")
		value = fVectorF_total->Value(angle /degree, energy /keV, fIndx, fIndy);
	else if(elname == "B")
		value = fVectorB_total->Value(angle /degree, energy /keV, fIndx, fIndy);
	else if(elname == "Ni")
		value = fVectorNi_total->Value(angle /degree, energy /keV, fIndx, fIndy);
	else if(elname == "Cu")
		value = fVectorCu_total->Value(angle /degree, energy /keV, fIndx, fIndy);
	return value;
}
// calculates energy loss, energy fwhm
void EventAction::CalcEnergyLeft(G4double depth, G4double energy, G4double ptr_pars[2])
{
	G4int location = 0;
	G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();
	G4double start_distance = 0;
	G4double ROI_region = detector->GetRBSROImin();
	G4double step_number = detector->GetEnLossStep();
	ptr_pars[0] = ptr_pars[1] = 0.;

	if (depth < layer_start_pos[0])                                     location = 1;
	else if (depth > layer_start_pos[0] && depth < layer_end_pos[0])    location = 2;
	else if (depth > layer_end_pos[0] && depth < layer_start_pos[1])    location = 3;
	else if (depth > layer_start_pos[1] && depth < layer_end_pos[1])    location = 4;
	else if (depth > layer_end_pos[1] && depth < layer_start_pos[2])    location = 5;
	else if (depth > layer_start_pos[2] && depth < layer_end_pos[2])    location = 6;
	else if (depth > layer_end_pos[2] && depth < layer_start_pos[3])    location = 7;
	else if (depth > layer_start_pos[3] && depth < layer_end_pos[3])    location = 8;
	else if (depth > layer_end_pos[3])                                  location = 9;
	else                                                                exit(1);

	G4double final_fwhm = 0;
	G4double final_energy = energy;
	G4double distance = 0;
	G4double cosine = cos(angle_of_exit);
	start_distance = depth - start_pos[location-1];

	while (location > 0) {
		G4Material* mat = mat_matrix[location];
		if (start_distance == 0)
			distance = (dist_matrix[location] / cosine);
		else
			distance = (start_distance / cosine);
		if ((distance / step_number) /nm < 1.)
			step_number = ceil(distance / (1. *nm));
		final_fwhm += CalculateTotalBohrStraggling(final_energy, particle, mat, distance);
		final_energy = CalcTotEnLoss(final_energy, distance, step_number, particle, mat);
		location--;
		start_distance = 0;
		// if final energy less than ROI, stop evaluation
		if (final_energy < ROI_region) {
			final_energy = 0;
			final_fwhm = 0;
			break;
		}
	}
	
	ptr_pars[0] = final_energy;
	ptr_pars[1] = final_fwhm;
}

void EventAction::PrepareDistances()
{
	for (G4int i=0; i<4; ++i) {
		thicknesses[i]     = detector->GetLength(i + 1);
		positions[i]       = (detector->GetLength(0) / 2)+detector->GetPosition(i).z();
		layer_start_pos[i] = positions[i] - thicknesses[i] / 2;
		layer_end_pos[i]   = positions[i] + thicknesses[i] / 2;
	}
	for (G4int i=0; i<3; ++i) {
		dist_between_l[i]  = layer_start_pos[i + 1] - layer_end_pos[i];
	}
	// prepare distance matrices for energy loss calculations
	dist_matrix[0] = 0;
	dist_matrix[1] = layer_start_pos[0];
	start_pos[0] = 0;
	mat_matrix[0] = sample_material[0];

	uint8_t x = 0;
	uint8_t y = 0;

	for (uint8_t i=1; i<10; i++) {
		if (i % 2 == 0) {
			start_pos[i] = layer_end_pos[x];
			mat_matrix[i] = sample_material[x+1];
			x++;
		} else {
			start_pos[i] = layer_start_pos[y];
			mat_matrix[i] = sample_material[0];
			y++;
		}
	}

	x = 0;
	y = 0;

	for (uint8_t i=2; i<10; i++) {
		if (i % 2 == 0) {
			dist_matrix[i] = thicknesses[x];
			x++;
		} else {
			dist_matrix[i] = dist_between_l[y];
			y++;
		}
	}
	// for (int i=0; i<10; ++i) {
		//G4cout << " i = " << i << " dist matrix = " << dist_matrix[i]/nm << " material = " << mat_matrix[i]->GetName() << G4endl;
	// }
}

G4double EventAction::GetGaussSum(G4double energy, G4double sigma)
{
	G4int GK = detector->GetGaussCounter();
	G4double Gauss_sum = 0;
	G4double gaussian_denominator = GenerateGaussian(0, 0, sigma);
	for (int k=-GK; k<=GK; k++)	{
		G4double coef = k;
		G4double Gauss_en = energy * (1 - (coef / 1000));
		G4double Gauss_value = (GenerateGaussian(Gauss_en, energy, sigma)) / gaussian_denominator;
		Gauss_sum += Gauss_value;
	}
	return Gauss_sum;
}

void EventAction::FillGaussians(G4double energy, G4double sigma, G4double steps, G4double yield, G4int element, G4int layer)
{
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	G4double RBS_norm_dist = 0.1*nm;
	G4double ROI_region	= detector->GetRBSROImin();
	G4double Gauss_sum = GetGaussSum(energy, sigma);
	G4int GK = detector->GetGaussCounter();
	G4double normalized_gaussian = GenerateGaussian(0, 0, sigma);
	G4double steps_to_dist = (steps / RBS_norm_dist);
	G4double yield_x_steps = steps_to_dist * yield;

	for (int k=-GK; k<=GK; k++) {
		G4double coef = (G4double)k;
		G4double Gauss_en 	= energy * (1 - (coef / 1000));
		G4double Gauss_value = (GenerateGaussian(Gauss_en, energy, sigma)) / normalized_gaussian;
		G4double mod_G_value = Gauss_value / Gauss_sum;
		if (Gauss_en /MeV > ROI_region) {
			analysisManager->FillH1(HistoBase[layer][0], Gauss_en, (mod_G_value * yield_x_steps)); // total
			analysisManager->FillH1(HistoBase[layer][1], Gauss_en, (mod_G_value * yield_x_steps)); // substrate
			if (element == 0)
				analysisManager->FillH1(HistoBase[layer][2], Gauss_en, (mod_G_value * yield_x_steps)); // el1
			else if (element == 1)
				analysisManager->FillH1(HistoBase[layer][3], Gauss_en, (mod_G_value * yield_x_steps)); // el2
		}
	}
}

void EventAction::CalculateRBS(G4int layer, G4double energy, G4ThreeVector position, G4ThreeVector momDir, G4double steps, G4ParticleDefinition* fParticle, G4int elements, G4double *nud_el)
{
	G4double angle_of_incidence = atan(fPrimary->GetParticleGun()->GetParticleMomentumDirection().x() / fPrimary->GetParticleGun()->GetParticleMomentumDirection().z());
	G4double rbs_angle = detector->GetRBSAngle();
	G4double scattering_angle = (rbs_angle - angle_of_incidence);
	G4double xsecRTR = 0, RBS_yield = 0, tot_sigma = 0, primary_energy = energy, final_energy = 0, bohr_straggling = 0, diff_angle = 0;

	if (detector->GetConstAngle())
		diff_angle = scattering_angle;
	else {
		G4double curr_angle = atan(momDir.x() / momDir.z());
		diff_angle = rbs_angle - curr_angle;
	}
	//
	G4double ROI_region = detector->GetRBSROImin();
	const G4Material* dead_material = G4NistManager::Instance()->FindOrBuildMaterial(detector->GetDeadLayer());
	G4double det_length = detector->GetLength(0);
	if (detector->GetRBSCalc()) {
		if (layer < NUMB_MAX_LAYERS) {
			for (G4int i=0; i<elements; i++) {
				G4int A1 = fPrimary->GetParticleGun()->GetParticleDefinition()->GetAtomicMass();
				G4double RecEn = RecoilEnergy(primary_energy, diff_angle, A1, Mnumb[layer][i]);
				if (RecEn > ROI_region) {
					bool is_within_ROI = CalculateFinalEnergy(RecEn, position, final_energy, bohr_straggling, det_length, ROI_region);

					if (is_within_ROI) {
						ProcessRBSYield(primary_energy, final_energy, bohr_straggling, layer, i, diff_angle, A1, fParticle, dead_material, steps, angle_of_incidence, nud_el[i]);
					} else if (final_energy > 0.95 * ROI_region) {
						Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
						run->MaxRBSDepth(position.z() + det_length / 2);
						run->AddCount();
					}
				}
			}
		}
	}
}

bool EventAction::CalculateFinalEnergy(G4double recoil_energy, G4ThreeVector position, G4double& final_energy, G4double& bohr_straggling, G4double detector_length, G4double ROI_region)
{
	G4double ptr_pars[2];
	G4double newWorldPosition = (position.z() + detector->GetLength(0) / 2);
	CalcEnergyLeft(newWorldPosition, recoil_energy, ptr_pars);
	final_energy = ptr_pars[0];
	bohr_straggling = ptr_pars[1];

	return final_energy > ROI_region;
}

void EventAction::ProcessRBSYield(G4double primary_energy, G4double& final_energy, G4double& bohr_straggling, G4int layer, G4int element_idx, G4double diff_angle, G4int A1, G4ParticleDefinition* fParticle, const G4Material* dead_material, G4double steps, G4double angle_of_incidence, G4double nud_el)
{
	G4double xsecRTR = GetCrossSection(primary_energy, layer, element_idx, diff_angle, nud_el);
	G4double solidAngle = detector->GetSolidAngle();
	G4double dead_thickness = detector->GetDeadLayerThickness();
	G4int Z1 = fPrimary->GetParticleGun()->GetParticleDefinition()->GetAtomicNumber();
	G4double RBS_yield = CalculateTotalRBSYield(primary_energy, A1, Mnumb[layer][element_idx], Z1, Znumb[layer][element_idx], diff_angle, solidAngle, xsecRTR, Adens[layer][element_idx], angle_of_incidence);
	final_energy = CalculateDeadLayerEffect(final_energy, dead_material, dead_thickness, fParticle);
	bohr_straggling += CalculateTotalBohrStraggling(final_energy, fParticle, dead_material, dead_thickness);
	G4double tot_sigma = GetDetectorSigmaSq(Z1, final_energy) + bohr_straggling;

	FillGaussians(final_energy, tot_sigma, steps, RBS_yield, element_idx, layer);
}

G4double EventAction::GetCrossSection(G4double primary_energy, G4int layer, G4int element_idx, G4double diff_angle, G4double nud)
{
	if (detector->GetSigmaCalc()) {
		G4String el_name = sample_material[layer]->GetElement(element_idx)->GetName();
		G4double xsec_rtr = Get2DRTRValue(primary_energy, el_name, diff_angle);
		return xsec_rtr * nud;
	}
	return 1.0;
}

void EventAction::PrepareHistoBase(void)
{
	uint8_t x = 0;
	for (uint8_t i=0; i<NUMBER_OF_MAX_LAYERS; ++i) {
		HistoBase[i][0] = 20;
		for (uint8_t j=1; j<4; j++) {
			if(i == 0)
				x = 21;
			else if(i == 1)
				x = 24;
			else if(i == 2)
				x = 27;
			else if(i == 3)
				x = 30;
			else if(i == 4)
				x = 33;
			else
				exit(1);
			HistoBase[i][j] = x + (j - 1);
		}
	}
}

void EventAction::FillH2ChannelingHistos(G4double z_pos, G4double ch_x_pos, G4double ch_y_pos, G4double wo_x_pos, G4double wo_y_pos)
{
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->FillH2(1, wo_x_pos, z_pos);
	analysisManager->FillH2(2, wo_y_pos, z_pos);
	analysisManager->FillH2(3, ch_x_pos, z_pos);
	analysisManager->FillH2(4, ch_y_pos, z_pos);

	for (uint8_t i=0; i<11; i++) {
		if (z_pos > depth_for_histos[i] && z_pos < depth_for_histos[i + 1])
			analysisManager->FillH2(5 + i, ch_x_pos, ch_y_pos);
	}
}

void EventAction::FillSigmaCalcVectors(void)
{
	G4String mat_name;
	G4double A1 = fPrimary->GetParticleGun()->GetParticleDefinition()->GetAtomicMass();
	for (uint8_t i=0; i<NUMBER_OF_MAX_LAYERS; i++) {
		for (uint8_t j=0; j<NoOfElements[i]; j++) {
			mat_name = sample_material[i]->GetElement(j)->GetName();
			if (mat_name == "Si" && fVectorSi_total == nullptr) {
				fVectorSi_total = Get2DRTRVector("Si", A1);
				fVectorSi_total->SetBicubicInterpolation(true);
			} else if (mat_name == "O" && fVectorO_total == nullptr) {
				fVectorO_total = Get2DRTRVector("O",A1);
				fVectorO_total->SetBicubicInterpolation(true);
			} else if (mat_name == "B"  && fVectorB_total == nullptr) {
				fVectorB_total = Get2DRTRVector("B",A1);
				fVectorB_total->SetBicubicInterpolation(true);
			} else if (mat_name == "C"  && fVectorC_total == nullptr) {
				fVectorC_total = Get2DRTRVector("C",A1);
				fVectorC_total->SetBicubicInterpolation(true);
			} else if (mat_name == "F" && fVectorF_total == nullptr) {
				fVectorF_total = Get2DRTRVector("F",A1);
				fVectorF_total->SetBicubicInterpolation(true);
			} else if (mat_name == "N"  && fVectorN_total == nullptr) {
				fVectorN_total = Get2DRTRVector("N",A1);
				fVectorN_total->SetBicubicInterpolation(true);
			} else if (mat_name == "Na" && fVectorNa_total == nullptr) {
				fVectorNa_total = Get2DRTRVector("Na",A1);
				fVectorNa_total->SetBicubicInterpolation(true);
			} else if (mat_name == "Ni" && fVectorNi_total == nullptr) {
				fVectorNi_total = Get2DRTRVector("Ni",A1);
				fVectorNi_total->SetBicubicInterpolation(true);
			} else if (mat_name == "Cu" && fVectorCu_total == nullptr) {
				fVectorCu_total = Get2DRTRVector("Cu",A1);
				fVectorCu_total->SetBicubicInterpolation(true);
			}
		} // end of cycle throu elements
	} // end of cycle throu materials
}

G4double EventAction::GetDetectorSigmaSq(G4int Z1, G4double final_energy)
{
	G4double detector_fwhm = 0;
	G4double det_FWHM = detector->GetDetectorResolution();

	if (det_FWHM < 10. *keV)
		detector_fwhm = (10. *keV) / 2.355;
	else
		detector_fwhm = det_FWHM / 2.355;

	if (Z1 == 3 && detector->GetCalcFWHM())
		detector_fwhm = (CalcDetectorFWHM(final_energy, Z1) / 1000 *MeV) / 2.355;
	G4double sigma_det_sq = std::pow(detector_fwhm, 2);
	return sigma_det_sq;
}