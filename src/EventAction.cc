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
#include "Run.hh"
#include "HistoManager.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
//#include "G4TrajectoryContainer.hh"
//#include "G4Trajectory.hh"
//#include "G4VVisManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "CrystalDetectorHit.hh"
#include "SensitiveDetectorHit.hh"
#include "Analysis.hh"
// new addition
#include "G4Material.hh"
#include "G4Track.hh"
#include "DetectorConstruction.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4EmCalculator.hh"
//#include "G4LossTableManager.hh"
#include <cmath>
#include "G4NistManager.hh"
//test fluctuations
#include "G4IonYangFluctuationModel.hh"
#include "G4IonChuFluctuationModel.hh"
#include "G4hIonEffChargeSquare.hh"

#include "G4Pow.hh"
#include "RBSHelpers.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* DE, PrimaryGeneratorAction* PRA )
:G4UserEventAction(),detector(DE),fPrimary(PRA),
 fTotalEnergyDeposit(0.), fTotalEnergyFlow(0.), EnergyDeposit(0.), NonIonEnergyDeposit(0.),
  theta(0), TrakLenPrim(0.), nbStepsPrim(0.), TrakLenSec(0.), absStepPrim(0.), abs2StepPrim(0.), TrueTrakLen(0.), ProjTrakLen(0.), sdht_ID(-1), // paskutiniai is chan
  sd0_ID(-1),sd1_ID(-1),sd2_ID(-1),sd3_ID(-1),sd4_ID(-1)
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
	for (int i = 0; i<5; i++)
	{
		delete[] Znumb[i];
		delete[] Mnumb[i];
		delete[] Adens[i];
		delete[] HistoBase[i];
		delete[] nud_el[i];
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
	fTotalEnergyDeposit = 0.;
  	fTotalEnergyFlow = 0.;
 	// initialisation per event
 	EnergyDeposit  = 0.;
 	NonIonEnergyDeposit=0.;
 	theta=0.;
 	TrakLenPrim=0.;
 	TrakLenSec=0.;
 	nbStepsPrim=0.;
 	TrueTrakLen=0.;
	ProjTrakLen=0.;
 	absStepPrim = 0;
 	abs1StepPrim = 0;
 	abs2StepPrim = 0; 
 	abs3StepPrim = 0;
 	abs4StepPrim = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// end of functions relevant to RBS spectrum
void EventAction::AddEdep(G4double Edep)
{
  fTotalEnergyDeposit += Edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEflow(G4double Eflow)
{
  fTotalEnergyFlow += Eflow;
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
	Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());

	run->AddEdep (fTotalEnergyDeposit);
	run->AddEflow(fTotalEnergyFlow);
	run->AddTrakLenPrim(TrakLenPrim);
	run->AddNumberOfSteps(nbStepsPrim);
	run->AddTheta(theta);
	run->AddTrakLenSec(TrakLenSec);
	run->AddTrueRange(TrueTrakLen);
	run->AddProjRange(ProjTrakLen);
	
	run->absSTP(absStepPrim);
	run->abs1STP(abs1StepPrim);
	run->abs2STP(abs2StepPrim);
	run->abs3STP(abs3StepPrim);
	run->abs4STP(abs4StepPrim);    
	
	G4AnalysisManager::Instance()->FillH1(1,fTotalEnergyDeposit);
	G4AnalysisManager::Instance()->FillH1(3,fTotalEnergyFlow);  
	
	// ******************************
  	for (int i = 0; i<5; i++) {
		Znumb[i] = new G4double[4];
		Mnumb[i] = new G4double[4];
		Adens[i] = new G4double[4];
		HistoBase[i] = new G4int[4];
		nud_el[i] = new G4double[3];
	}

	for(int i=0; i<5; ++i) {
		for(int j=0; j<3; ++j) {
			nud_el[i][j] = 0;
		}
	}

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	//detector energy resolution
	G4int rbs_eval = detector->GetRBSCalc();
	G4ParticleDefinition*  fParticle = fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double A1 = fPrimary->GetParticleGun()->GetParticleDefinition()->GetAtomicMass();
   	G4double primary_energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    run->PrimaryEnergy(primary_energy);
	G4int track_histos = detector->GetTrackingHistos();
	//**************************************************************************
	// Incidence angle
	G4double angle_of_incidence = 0.;
	angle_of_incidence = atan(fPrimary->GetParticleGun()->GetParticleMomentumDirection().x()/fPrimary->GetParticleGun()->GetParticleMomentumDirection().z());
	run->Inc_angle(angle_of_incidence);

	PrepareHistoBase();
	// prepare depth intervals
	detector->CopyArrayForDepth(depth_for_histos);
	
	//scattering angle
	Angle = detector->GetRBSAngle();
	angle_of_exit = pi - Angle;
  	G4EmCalculator emCalculator;
    G4String dead_material_name = detector->GetDeadLayer();
    		
	G4ThreeVector ch_pos;

	for (int i=0;i<5;i++) {
		sample_material[i] 	= detector->GetMaterialM(i);
		NoOfElements[i]    	= sample_material[i]->GetNumberOfElements();
	}
	G4int max_M = 0;

	for (int i=0;i<5;i++) {
		for (int j = 0; j<NoOfElements[i]; j++) {
			Znumb[i][j] = sample_material[i]->GetElement(j)->GetZ();
			Mnumb[i][j] = sample_material[i]->GetElement(j)->GetA()/(g/mole);
			if(Mnumb[i][j] > max_M)
				max_M = Mnumb[i][j];
			const G4double *atomDensVector	= sample_material[i]->GetVecNbOfAtomsPerVolume();
			Adens[i][j] = atomDensVector[j]/(1/cm3);
		}
	}

	// RTR vector fillup
	if(detector->GetSigmaCalc() && rbs_eval) {
  		G4String mat_name;
  		for(int i=0;i<5;i++) {
  			for(int j=0;j<NoOfElements[i]; j++)	{
				mat_name = sample_material[i]->GetElement(j)->GetName();
			  	if 	(mat_name == "Si" && fVectorSi_total == nullptr) { 
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
  	}	// end of sigma calc

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    G4ThreeVector ssd[4];
    ssd[0]= G4ThreeVector(0.,0.,0.);
    ssd[1]= G4ThreeVector(0.,0.,0.);
    ssd[2]= G4ThreeVector(0.,0.,0.);
    ssd[3]= G4ThreeVector(0.,0.,0.);

    G4double energy = 0.;
    G4double kinen = 0.;

    if(sdht_ID == -1) {
        G4String sdName;
        if(SDman->FindSensitiveDetector(sdName="telescope",0)){
            sdht_ID = SDman->GetCollectionID(sdName="telescope/collection");
        }
    }
    
    SensitiveDetectorHitsCollection* sdht = 0;
    G4HCofThisEvent *hce = evt->GetHCofThisEvent();
    
    if(hce){
        if(sdht_ID != -1){
            G4VHitsCollection* aHCSD = hce->GetHC(sdht_ID);
            sdht = (SensitiveDetectorHitsCollection*)(aHCSD);
        }
    }
    int bTotalHits = 0;
    if(sdht){
        int n_hit_sd = sdht->entries();
		for(int i2=0;i2<4;i2++){
            for(int i1=0;i1<n_hit_sd;i1++) {
                SensitiveDetectorHit* aHit = (*sdht)[i1];
                if(aHit->GetLayerID()==i2) {
                    ssd[i2] = aHit->GetWorldPos();
                    bTotalHits++;	// checks, how many detectors got a hit
                }
                if(aHit->GetLayerID()==2) {	// info from the 3rd detector
                    energy  = aHit->GetKinE();
					analysisManager->FillH1(36, energy);
                }
				// particle hit in the last detector [no 4]
				if(aHit->GetLayerID()==3) {
					run->addh(); // adds hits
					kinen = aHit->GetKinE(); // takes kinetinic energy
					//analysisManager->FillH1(36, kinen);
					run->addkinen(kinen); // add kinetic energy 
				}
         	}
        }
    }

  	// sensitive detectors
  	G4String sdName;
	if(sd0_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector",0))
            sd0_ID = SDman->GetCollectionID(sdName="crystaldetector/collection");
   	 }
   	if(sd1_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector2",0))
            sd1_ID = SDman->GetCollectionID(sdName="crystaldetector2/collection");
   	 }
   	if(sd2_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector3",0))
            sd2_ID = SDman->GetCollectionID(sdName="crystaldetector3/collection");
   	 }
   	if(sd3_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector4",0))
            sd3_ID = SDman->GetCollectionID(sdName="crystaldetector4/collection");
   	 }
   	if(sd4_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector5",0))
            sd4_ID = SDman->GetCollectionID(sdName="crystaldetector5/collection");
   	 }

    CrystalDetectorHitsCollection* sd0 = 0;
    G4HCofThisEvent *hx0 = evt->GetHCofThisEvent();
    if(hx0){
       	if(sd0_ID != -1){
           		G4VHitsCollection* aHCSD0 = hx0->GetHC(sd0_ID);
           		sd0 = (CrystalDetectorHitsCollection*)(aHCSD0);
       		}
    	}

    CrystalDetectorHitsCollection* sd1 = 0;
    G4HCofThisEvent *hx1 = evt->GetHCofThisEvent();
    if(hx1){
       	if(sd1_ID != -1){
           		G4VHitsCollection* aHCSD1 = hx1->GetHC(sd1_ID);
           		sd1 = (CrystalDetectorHitsCollection*)(aHCSD1);
       		}
    	}
    
	CrystalDetectorHitsCollection* sd2 = 0;
    G4HCofThisEvent *hx2 = evt->GetHCofThisEvent();
    if(hx2){
       	if(sd2_ID != -1){
           		G4VHitsCollection* aHCSD2 = hx2->GetHC(sd2_ID);
           		sd2 = (CrystalDetectorHitsCollection*)(aHCSD2);
       		}
    	}
    
	CrystalDetectorHitsCollection* sd3 = 0;
    G4HCofThisEvent *hx3 = evt->GetHCofThisEvent();
    if(hx3){
       	if(sd3_ID != -1){
           		G4VHitsCollection* aHCSD3 = hx3->GetHC(sd3_ID);
           		sd3 = (CrystalDetectorHitsCollection*)(aHCSD3);
       		}
    	}
    
	CrystalDetectorHitsCollection* sd4 = 0;
    G4HCofThisEvent *hx4 = evt->GetHCofThisEvent();
    if(hx4){
       	if(sd4_ID != -1){
           		G4VHitsCollection* aHCSD4 = hx4->GetHC(sd4_ID);
           		sd4 = (CrystalDetectorHitsCollection*)(aHCSD4);
       		}
    	}

	// prepare distances
	PrepareDistances();

	G4double steps;
	G4ThreeVector position;
	G4double sample_energy =0.;
	G4double tot_step =0.;
	G4ThreeVector momDir;
	// chan params
	G4double efx = 0;
	G4double efy = 0;
	G4double eld = 0;
	G4double nud = 0;
	G4double efxavg = 0;
	G4double efyavg = 0;
	G4double eldavg = 0;
	G4double nudavg = 0;
	G4double nudavg_a = 0;
	G4double nudavg_b = 0;
	G4double nudavg_c = 0;

	// fill ntuples
	if(bTotalHits > 3 ) {
        G4double angXin  = ((ssd[1].x() - ssd[0].x()) / (ssd[1].z() - ssd[0].z()));
        G4double angYin  = ((ssd[1].y() - ssd[0].y()) / (ssd[1].z() - ssd[0].z()));
        analysisManager->FillNtupleDColumn(0, atan(angXin)/degree); 
        analysisManager->FillNtupleDColumn(1, atan(angYin)/degree); 
        double posXin = ssd[1].x() - atan(angXin) * ssd[1].z();
        double posYin = ssd[1].y() - atan(angYin) * ssd[1].z();
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
    	    posXout = ssd[2].x() - atan(angXout)*ssd[2].z();
    	    posYout = ssd[2].y() - atan(angYout)*ssd[2].z();
		}
        analysisManager->FillNtupleDColumn(42, atan(angXout)/degree); 
        analysisManager->FillNtupleDColumn(43, atan(angYout)/degree); 
        analysisManager->FillNtupleDColumn(40, posXout / CLHEP::mm); 
        analysisManager->FillNtupleDColumn(41, posYout / CLHEP::mm);
	}

/*
	for (int i = 10000; i >= 100; i--) {	
    	if(i%100 == 0) {
    		G4double en = i*keV;
    		//G4double stp_pwr = emCalculator.GetDEDX(en,fParticle,sample_material[0])/(eV/angstrom);
			G4double stp_pwr = emCalculator.GetDEDX(en,fParticle,sample_material[0])/(MeV/centimeter);
			G4double dens = sample_material[0]->GetDensity()/(milligram/cm3);
			//G4cout << dens << G4endl;
    		//G4double stp_pwr = emCalculator.ComputeTotalDEDX(en,fParticle,sample_material[0])/(MeV/(milligram/centimeter2));
    		G4cout << " " << en/keV << " " << stp_pwr/dens << G4endl;
    	}
	}*/


	// mother volume sensitive detector
    if(sd0)
    {
        int n_hit_sd = sd0->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++) {
        	CrystalDetectorHit* aHit = (*sd0)[i1];
           	steps = aHit->GetStep();
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();
			ch_pos = aHit->GetChPos();
			nud = aHit->GetNud();
			efx = aHit->GetEFX();
			efy = aHit->GetEFY();
			eld = aHit->GetEld();
			momDir = aHit->GetWorldMomentumDirection();
			if(nud == 1.)
				nud_el[0][0] = nud_el[0][1] = nud_el[0][2] = 1;
			else {
				nud_el[0][0] = aHit->GetNuD_a();
				nud_el[0][1] = aHit->GetNuD_b();
				nud_el[0][2] = aHit->GetNuD_c();
			}

			tot_step +=steps;
			eldavg += eld*steps;
			nudavg += nud*steps;
			nudavg_a += nud_el[0][0]*steps;
			nudavg_b += nud_el[0][1]*steps;
			nudavg_c += nud_el[0][2]*steps;
			efxavg += efx*steps;
			efyavg += efy*steps;

			G4double z_pos = position.z()+detector->GetLength(0)/2;
			if(track_histos) {
				FillH2ChannelingHistos(z_pos, ch_pos.x(), ch_pos.y(),position.x(),position.y());
				analysisManager->FillP1(1,z_pos,nud);
				analysisManager->FillP1(2,z_pos,eld);
			}
			CalculateRBS(0,sample_energy,position,momDir,steps,fParticle);
		}
		run->add_total_step(tot_step);
  	}
	if(tot_step > 0) {
		nudavg /= tot_step;
		eldavg /= tot_step;
		efxavg /= tot_step;
		efyavg /= tot_step;
		nudavg_a /= tot_step;
		nudavg_b /= tot_step;
		nudavg_c /= tot_step;
	} else {
		nudavg = 0;
		eldavg = 0;
		efxavg = 0;
		efyavg = 0;
		nudavg_a = 0;
		nudavg_b = 0;
		nudavg_c = 0;
	}
	analysisManager->FillNtupleDColumn(5, efxavg/ CLHEP::eV * CLHEP::angstrom);
	analysisManager->FillNtupleDColumn(6, efyavg/ CLHEP::eV * CLHEP::angstrom);
	analysisManager->FillNtupleDColumn(7, nudavg);
	analysisManager->FillNtupleDColumn(8, eldavg);
	// elements
	/*
	analysisManager->FillNtupleDColumn(28, nudavg_a);
	analysisManager->FillNtupleDColumn(29, nudavg_b);
	analysisManager->FillNtupleDColumn(30, nudavg_c);
	*/
	//analysisManager->AddNtupleRow();
	// end of mother sensitive detector

//======================================================
// LAYER 1
//======================================================
	efx = 0;
	efy = 0;
	eld = 0;
	nud = 0;
	efxavg = 0;
	efyavg = 0;
	eldavg = 0;
	nudavg = 0;
	nudavg_a = 0;
	nudavg_b = 0;
	nudavg_c = 0;

	// layer1 sensitive detector
    if(sd1){
    	tot_step =0.;
        int n_hit_sd = sd1->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            CrystalDetectorHit* aHit = (*sd1)[i1];
           	steps = aHit->GetStep();
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();
			ch_pos = aHit->GetChPos();
			nud = aHit->GetNud();
			efx = aHit->GetEFX();
			efy = aHit->GetEFY();
			eld = aHit->GetEld();
			momDir = aHit->GetWorldMomentumDirection();
			if(nud == 1.)
				nud_el[1][0] = nud_el[1][1] = nud_el[1][2] = 1;
			else {
				nud_el[1][0] = aHit->GetNuD_a();
				nud_el[1][1] = aHit->GetNuD_b();
				nud_el[1][2] = aHit->GetNuD_c();
			}

			tot_step +=steps;
			eldavg += eld*steps;
			nudavg += nud*steps;
			nudavg_a += nud_el[1][0]*steps;
			nudavg_b += nud_el[1][1]*steps;
			nudavg_c += nud_el[1][2]*steps;
			efxavg += efx*steps;
			efyavg += efy*steps;

			G4double z_pos = position.z()+detector->GetLength(0)/2;
			if(track_histos) {
				FillH2ChannelingHistos(z_pos, ch_pos.x(), ch_pos.y(),position.x(),position.y());
				analysisManager->FillP1(3,z_pos,nud);
				analysisManager->FillP1(4,z_pos,eld);
			}
			CalculateRBS(1,sample_energy,position,momDir,steps,fParticle);
		}
    	run->add_total_step(tot_step);
    }// end of layer1 sensitive detector
	if(tot_step > 0) {
		nudavg /= tot_step;
		eldavg /= tot_step;
		efxavg /= tot_step;
		efyavg /= tot_step;
		nudavg_a /= tot_step;
		nudavg_b /= tot_step;
		nudavg_c /= tot_step;
	} else {
		nudavg = 0;
		eldavg = 0;
		efxavg = 0;
		efyavg = 0;
		nudavg_a = 0;
		nudavg_b = 0;
		nudavg_c = 0;
	}
	analysisManager->FillNtupleDColumn(9, efxavg/ CLHEP::eV * CLHEP::angstrom);
	analysisManager->FillNtupleDColumn(10, efyavg/ CLHEP::eV * CLHEP::angstrom);
	analysisManager->FillNtupleDColumn(11, nudavg);
	analysisManager->FillNtupleDColumn(12, eldavg);
	// elements
	/*
	analysisManager->FillNtupleDColumn(31, nudavg_a);
	analysisManager->FillNtupleDColumn(32, nudavg_b);
	analysisManager->FillNtupleDColumn(33, nudavg_c);
	*/

//======================================================
// LAYER 2
//======================================================
	efx = 0;
	efy = 0;
	eld = 0;
	nud = 0;
	efxavg = 0;
	efyavg = 0;
	eldavg = 0;
	nudavg = 0;
	nudavg_a = 0;
	nudavg_b = 0;
	nudavg_c = 0;
	// layer2 sensitive detector
    if(sd2){
    	tot_step =0.;
        int n_hit_sd = sd2->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            CrystalDetectorHit* aHit = (*sd2)[i1];
           	steps = aHit->GetStep();
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();
			ch_pos = aHit->GetChPos();
			nud = aHit->GetNud();
			efx = aHit->GetEFX();
			efy = aHit->GetEFY();
			eld = aHit->GetEld();
			momDir = aHit->GetWorldMomentumDirection();
			if(nud == 1.)
				nud_el[2][0] = nud_el[2][1] = nud_el[2][2] = 1;
			else {
				nud_el[2][0] = aHit->GetNuD_a();
				nud_el[2][1] = aHit->GetNuD_b();
				nud_el[2][2] = aHit->GetNuD_c();
			}

			tot_step +=steps;
			eldavg += eld*steps;
			nudavg += nud*steps;
			nudavg_a += nud_el[2][0]*steps;
			nudavg_b += nud_el[2][1]*steps;
			nudavg_c += nud_el[2][2]*steps;
			efxavg += efx*steps;
			efyavg += efy*steps;

			G4double z_pos = position.z()+detector->GetLength(0)/2;
			if(track_histos) {
				FillH2ChannelingHistos(z_pos, ch_pos.x(), ch_pos.y(),position.x(),position.y());
				analysisManager->FillP1(5,z_pos,nud);
				analysisManager->FillP1(6,z_pos,eld);
			}
			CalculateRBS(2,sample_energy,position,momDir,steps,fParticle);
		}
		run->add_total_step(tot_step); //all hits		
    }// end of layer2 sensitive detector
	if(tot_step > 0) {
		nudavg /= tot_step;
		eldavg /= tot_step;
		efxavg /= tot_step;
		efyavg /= tot_step;
		nudavg_a /= tot_step;
		nudavg_b /= tot_step;
		nudavg_c /= tot_step;
	} else {
		nudavg = 0;
		eldavg = 0;
		efxavg = 0;
		efyavg = 0;
		nudavg_a = 0;
		nudavg_b = 0;
		nudavg_c = 0;
	}
	analysisManager->FillNtupleDColumn(13, efxavg/ CLHEP::eV * CLHEP::angstrom);
	analysisManager->FillNtupleDColumn(14, efyavg/ CLHEP::eV * CLHEP::angstrom);
	analysisManager->FillNtupleDColumn(15, nudavg);
	analysisManager->FillNtupleDColumn(16, eldavg);
	// elements
	/*
	analysisManager->FillNtupleDColumn(34, nudavg_a);
	analysisManager->FillNtupleDColumn(35, nudavg_b);
	analysisManager->FillNtupleDColumn(36, nudavg_c);
	*/
//======================================================
// LAYER 3
//======================================================
	efx = 0;
	efy = 0;
	eld = 0;
	nud = 0;
	efxavg = 0;
	efyavg = 0;
	eldavg = 0;
	nudavg = 0;
	nudavg_a = 0;
	nudavg_b = 0;
	nudavg_c = 0;
	// layer3 sensitive detector
    if(sd3){
    	tot_step =0.;
        int n_hit_sd = sd3->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            CrystalDetectorHit* aHit = (*sd3)[i1];
           	steps = aHit->GetStep();
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();
			ch_pos = aHit->GetChPos();
			nud = aHit->GetNud();
			efx = aHit->GetEFX();
			efy = aHit->GetEFY();
			eld = aHit->GetEld();
			momDir = aHit->GetWorldMomentumDirection();
			if(nud == 1.)
				nud_el[3][0] = nud_el[3][1] = nud_el[3][2] = 1;
			else {
				nud_el[3][0] = aHit->GetNuD_a();
				nud_el[3][1] = aHit->GetNuD_b();
				nud_el[3][2] = aHit->GetNuD_c();
			}

			tot_step +=steps;
			eldavg += eld*steps;
			nudavg += nud*steps;
			nudavg_a += nud_el[3][0]*steps;
			nudavg_b += nud_el[3][1]*steps;
			nudavg_c += nud_el[3][2]*steps;
			efxavg += efx*steps;
			efyavg += efy*steps;

			G4double z_pos = position.z()+detector->GetLength(0)/2;
			if(track_histos) {
				FillH2ChannelingHistos(z_pos, ch_pos.x(), ch_pos.y(),position.x(),position.y());
				analysisManager->FillP1(7,z_pos,nud);
				analysisManager->FillP1(8,z_pos,eld);
			}
			CalculateRBS(3,sample_energy,position,momDir,steps,fParticle);
		}
		run->add_total_step(tot_step);
    }
	if(tot_step > 0) {
		nudavg /= tot_step;
		eldavg /= tot_step;
		efxavg /= tot_step;
		efyavg /= tot_step;
		nudavg_a /= tot_step;
		nudavg_b /= tot_step;
		nudavg_c /= tot_step;
	} else {
		nudavg = 0;
		eldavg = 0;
		efxavg = 0;
		efyavg = 0;
		nudavg_a = 0;
		nudavg_b = 0;
		nudavg_c = 0;
	}
	analysisManager->FillNtupleDColumn(17, efxavg/ CLHEP::eV * CLHEP::angstrom);
	analysisManager->FillNtupleDColumn(18, efyavg/ CLHEP::eV * CLHEP::angstrom);
	analysisManager->FillNtupleDColumn(19, nudavg);
	analysisManager->FillNtupleDColumn(20, eldavg);
	// elements
	/*
	analysisManager->FillNtupleDColumn(37, nudavg_a);
	analysisManager->FillNtupleDColumn(38, nudavg_b);
	analysisManager->FillNtupleDColumn(39, nudavg_c);
	*/

//======================================================
// LAYER 4
//======================================================
	efx = 0;
	efy = 0;
	eld = 0;
	nud = 0;
	efxavg = 0;
	efyavg = 0;
	eldavg = 0;
	nudavg = 0;
	nudavg_a = 0;
	nudavg_b = 0;
	nudavg_c = 0;
	// layer4 sensitive detector
    if(sd4){
    	tot_step =0.;
        int n_hit_sd = sd4->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            CrystalDetectorHit* aHit = (*sd4)[i1];
           	steps = aHit->GetStep();
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();
			ch_pos = aHit->GetChPos();
			nud = aHit->GetNud();
			efx = aHit->GetEFX();
			efy = aHit->GetEFY();
			eld = aHit->GetEld();
			momDir = aHit->GetWorldMomentumDirection();
			if(nud == 1.)
				nud_el[4][0] = nud_el[4][1] = nud_el[4][2] = 1;
			else {
				nud_el[4][0] = aHit->GetNuD_a();
				nud_el[4][1] = aHit->GetNuD_b();
				nud_el[4][2] = aHit->GetNuD_c();
			}

			tot_step +=steps;
			eldavg += eld*steps;
			nudavg += nud*steps;
			nudavg_a += nud_el[4][0]*steps;
			nudavg_b += nud_el[4][1]*steps;
			nudavg_c += nud_el[4][2]*steps;
			efxavg += efx*steps;
			efyavg += efy*steps;

			G4double z_pos = position.z()+detector->GetLength(0)/2;
			if(track_histos) {
				FillH2ChannelingHistos(z_pos, ch_pos.x(), ch_pos.y(),position.x(),position.y());
				analysisManager->FillP1(9,z_pos,nud);
				analysisManager->FillP1(10,z_pos,eld);
			}
			CalculateRBS(4,sample_energy,position,momDir,steps,fParticle);
		}
		run->add_total_step(tot_step);
	}
	if(tot_step > 0) {
		nudavg /= tot_step;
		eldavg /= tot_step;
		efxavg /= tot_step;
		efyavg /= tot_step;
		nudavg_a /= tot_step;
		nudavg_b /= tot_step;
		nudavg_c /= tot_step;
	} else {
		nudavg = 0;
		eldavg = 0;
		efxavg = 0;
		efyavg = 0;
		nudavg_a = 0;
		nudavg_b = 0;
		nudavg_c = 0;
	}
	analysisManager->FillNtupleDColumn(21, efxavg/ CLHEP::eV * CLHEP::angstrom);
	analysisManager->FillNtupleDColumn(22, efyavg/ CLHEP::eV * CLHEP::angstrom);
	analysisManager->FillNtupleDColumn(23, nudavg);
	analysisManager->FillNtupleDColumn(24, eldavg);
	// elements
	/*
	analysisManager->FillNtupleDColumn(40, nudavg_a);
	analysisManager->FillNtupleDColumn(41, nudavg_b);
	analysisManager->FillNtupleDColumn(42, nudavg_c);
	*/
	analysisManager->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EventAction::GenerateGaussian(G4double x, G4double y, G4double sigma_sq)
{
	G4double inv_sqrt_2pi = 0.398942280401433;
	G4double a = (x-y);
	return (inv_sqrt_2pi/std::sqrt(sigma_sq)) * std::exp(-0.5*a*a/sigma_sq);
}

G4Physics2DVector* EventAction::Get2DRTRVector(G4String element, G4int A)
{
	G4Physics2DVector* vector = new G4Physics2DVector;
	G4String start_file = "RTR_values/";
	G4String end_file = "_total.txt";
	
	G4String folder_def;
	G4String part_def;

	if (A == 1) {
		folder_def = "Hydrogen/";
		part_def = "H_";
	} else if (A == 4) {
		folder_def = "Helium/";
		part_def = "He_";
	} else
		return 0;

	G4String tot_filename = start_file+folder_def+part_def+element+end_file;
	G4cout << " Reading RTR values from file: " << tot_filename << G4endl;
	std::ifstream vFile(tot_filename);
	if(!vFile.is_open()) {
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

	G4double value = 0.;
	if(elname == "Si")
		value = fVectorSi_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "O")
		value = fVectorO_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "Na")
		value = fVectorNa_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "N")
		value = fVectorN_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "C")
		value = fVectorC_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "F")
		value = fVectorF_total->Value(angle/degree,energy/keV,fIndx,fIndy);	
	else if(elname == "B")
		value = fVectorB_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "Ni")
		value = fVectorNi_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "Cu")
		value = fVectorCu_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else 
		value = 1.;
	return value;
}
// calculates energy loss, energy fwhm, multiple scattering
G4double* EventAction::CalcEnergyLeft(G4double depth, G4double energy)
{
	G4double* final_parameters 	= new G4double[2];
	G4int location = 0;
	G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();
	G4double start_distance = 0;

	if(depth < layer_start_pos[0])									location = 1;
	else if(depth > layer_start_pos[0] && depth < layer_end_pos[0])	location = 2;
	else if(depth > layer_end_pos[0] && depth < layer_start_pos[1])	location = 3;
	else if(depth > layer_start_pos[1] && depth < layer_end_pos[1])	location = 4;
	else if(depth > layer_end_pos[1] && depth < layer_start_pos[2])	location = 5;
	else if(depth > layer_start_pos[2] && depth < layer_end_pos[2])	location = 6;
	else if(depth > layer_end_pos[2] && depth < layer_start_pos[3])	location = 7;
	else if(depth > layer_start_pos[3] && depth < layer_end_pos[3])	location = 8;
	else if(depth > layer_end_pos[3])								location = 9;
	else
		exit(1);


	G4double test_fwhm = 0;
	G4double test_energy = energy;
	start_distance = depth - start_pos[location-1];

	while(location > 0) {
		test_fwhm = test_fwhm + CalcStragglingInLayer(location,start_distance,test_energy,particle);
		test_energy = CalcLossInLayer(location, start_distance,test_energy,particle);
		location--;
		start_distance = 0;
	}
	
	final_parameters[0] = test_energy;
	final_parameters[1] = test_fwhm;

	return final_parameters;
}

void EventAction::PrepareDistances()
{
	for(int i=0; i<4; ++i) {
		thicknesses[i]		= detector->GetLength(i+1);
		positions[i] 		= (detector->GetLength(0)/2)+detector->GetPosition(i).z();
		layer_start_pos[i] 	= positions[i] - thicknesses[i]/2;
		layer_end_pos[i] 	= positions[i] + thicknesses[i]/2;
	}
	for(int i = 0; i<3;++i){
		dist_between_l[i] 	= layer_start_pos[i+1] - layer_end_pos[i];
	}
	// prepare distance matrices for energy loss calculations
	dist_matrix[0] = 0;
	dist_matrix[1] = layer_start_pos[0];
	start_pos[0] = 0;
	mat_matrix[0] = sample_material[0];

	int x = 0;
	int y = 0;

	for(int i=1; i<10;i++) {
		if(i%2 == 0) {
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

	for(int i=2;i<10;i++) {
		if(i%2==0) {
			dist_matrix[i] = thicknesses[x];
			x++;
		} else {
			dist_matrix[i] = dist_between_l[y];
			y++;
		}
	}
	for(int i=0; i<10; ++i) {
		//G4cout << " i = " << i << " dist matrix = " << dist_matrix[i]/nm << " material = " << mat_matrix[i]->GetName() << G4endl;
	}
}

G4double EventAction::CalcLossInLayer(G4int i, G4double dist, G4double energy, G4ParticleDefinition* particle)
{
	G4double final_energy = 0;
	G4double ROI = detector->GetRBSROImin();
	G4double step_number = detector->GetEnLossStep();
	G4double distance = 0;
	 
	if(dist == 0)
		distance = (dist_matrix[i]/cos(angle_of_exit));
	else
		distance = dist/cos(angle_of_exit);
	if((distance/step_number)/nm < 1.)
		step_number = ceil(distance/(1.*nm));

	G4Material* mat = mat_matrix[i];
	if(energy > ROI)
		final_energy = CalcTotEnLoss(energy, distance, step_number, particle, mat);
	else
		final_energy = 0;
	return final_energy;
}

G4double EventAction::CalcStragglingInLayer(G4int i, G4double dist, G4double energy, G4ParticleDefinition* particle)
{
	G4double final_stragg = 0;
	G4double ROI = detector->GetRBSROImin();
	G4double step_number = detector->GetEnLossStep();
	G4double distance = 0;
	if(dist == 0)
		distance = (dist_matrix[i]/cos(angle_of_exit));
	else
		distance = dist/cos(angle_of_exit);
	if((distance/step_number)/nm < 1.)
		step_number = ceil(distance/(1.*nm));
	G4Material* mat = mat_matrix[i];
	if(energy > ROI)
		final_stragg = CalculateTotalBohrStraggling(energy,particle,mat,distance);
	else
		final_stragg = 0;
	return final_stragg;	
}

G4double EventAction::GetGaussSum(G4double energy, G4double sigma)
{
	G4int GK = detector->GetGaussCounter();
	G4double Gauss_sum = 0;
	for (int k=-GK; k<=GK; k++)	{
		G4double coef = k;
		G4double Gauss_en = energy*(1-(coef/1000));
		G4double Gauss_value = (GenerateGaussian(Gauss_en,energy,sigma))/GenerateGaussian(0,0,sigma);
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

	for(int k=-GK; k<=GK; k++) {
		G4double coef = (G4double)k;
		G4double Gauss_en 	= energy*(1-(coef/1000));
		G4double Gauss_value = (GenerateGaussian(Gauss_en,energy,sigma))/GenerateGaussian(0,0,sigma);
		G4double  mod_G_value = Gauss_value/Gauss_sum;
		if (Gauss_en/MeV > ROI_region) {
			//run->add_entry_reach(1);
			analysisManager->FillH1(HistoBase[layer][0], Gauss_en, (yield*mod_G_value*(steps/RBS_norm_dist)));	// total
			analysisManager->FillH1(HistoBase[layer][1], Gauss_en, (yield*mod_G_value*(steps/RBS_norm_dist)));	//substrate
			if (element == 0)
				analysisManager->FillH1(HistoBase[layer][2],Gauss_en, (yield*mod_G_value*(steps/RBS_norm_dist)));	//el1
			else if (element == 1)
				analysisManager->FillH1(HistoBase[layer][3],Gauss_en, (yield*mod_G_value*(steps/RBS_norm_dist)));
			}
	}
}

void EventAction::CalculateRBS(G4int layer, G4double energy, G4ThreeVector position, G4ThreeVector momDir, G4double steps, G4ParticleDefinition* fParticle)
{
	G4double A1 = fPrimary->GetParticleGun()->GetParticleDefinition()->GetAtomicMass();
	G4double Z1 = fPrimary->GetParticleGun()->GetParticleDefinition()->GetAtomicNumber();
	G4double angle_of_incidence = 0.;
	angle_of_incidence = atan(fPrimary->GetParticleGun()->GetParticleMomentumDirection().x()/fPrimary->GetParticleGun()->GetParticleMomentumDirection().z());
	//G4cout << " incidence angle = " << angle_of_incidence/degree << G4endl;
	//G4cout << " scattering angle = " << (Angle-angle_of_incidence)/degree << " exit angle = " << angle_of_exit/degree << G4endl;
	G4double scattering_angle = (Angle - angle_of_incidence);
    G4String dead_material_name = detector->GetDeadLayer();
	const G4Material* dead_material = G4NistManager::Instance()->FindOrBuildMaterial(dead_material_name);
    G4Material* d_mat = G4NistManager::Instance()->FindOrBuildMaterial(dead_material_name);

	G4double det_FWHM = detector->GetDetectorResolution();
	G4double xsecRTR = 0;
	G4double RBS_yield = 0;
	G4double* parameters;
	G4double tot_sigma;
	G4double primary_energy = energy;
	G4double final_energy = 0;
	G4double bohr_straggling = 0.;
	G4double detector_fwhm;
	if(det_FWHM < 10.*keV)
		detector_fwhm = (10.*keV)/2.355;
	else
		detector_fwhm = det_FWHM/2.355;

	G4int use_const_angle 		= detector->GetConstAngle();
	G4int rbs_eval 				= detector->GetRBSCalc();
	G4double ROI_region			= detector->GetRBSROImin();
	G4int fwhm_calc 			= detector->GetCalcFWHM();
    G4double RBS_norm_dist 		= 0.1*nm;
    G4double dead_thickness 	= detector->GetDeadLayerThickness();
	G4double solidAngle 		= detector->GetSolidAngle();

	G4double newWorldPosition = (position.z()+detector->GetLength(0)/2);
	G4double curr_angle = atan(momDir.x()/momDir.z());
	G4double diff_angle;
	
	if(use_const_angle == 1)
		//diff_angle = Angle;
		diff_angle = scattering_angle;
	else
		diff_angle = Angle-curr_angle;
	if(rbs_eval) {
		for(int i=0;i<NoOfElements[layer];i++) {
			G4double RecEn = RecoilEnergy(primary_energy,diff_angle,A1,Mnumb[layer][i]);
			G4String el_name = sample_material[layer]->GetElement(i)->GetName();
			if(detector->GetSigmaCalc())
				xsecRTR = Get2DRTRValue(primary_energy, el_name,diff_angle);
			else
				xsecRTR = 1.;
			if (RecEn/MeV > ROI_region)
			{
				xsecRTR = xsecRTR*nud_el[layer][i];
				RBS_yield = CalculateTotalRBSYield(primary_energy, A1, Mnumb[layer][i], Z1, Znumb[layer][i], diff_angle, RBS_norm_dist, solidAngle, xsecRTR, Adens[layer][i],angle_of_incidence);
				parameters = CalcEnergyLeft(newWorldPosition, RecEn);
				final_energy = parameters[0];
				bohr_straggling = parameters[1];
				parameters[0] = parameters[1] = 0;
				// filling of RBS histograms
				if (final_energy > ROI_region) {
					final_energy = CalculateDeadLayerEffect(final_energy, dead_material,dead_thickness,fParticle);
					G4double dead_layer_straggling = CalculateTotalBohrStraggling(final_energy, fParticle, d_mat, dead_thickness);
					bohr_straggling += dead_layer_straggling;
					if(Z1 == 3 && fwhm_calc)
						detector_fwhm = (CalcDetectorFWHM(final_energy,Z1)/1000*MeV)/2.355;
					G4double sigma_det_sq = std::pow(detector_fwhm,2);
					tot_sigma = sigma_det_sq + bohr_straggling;
					FillGaussians(final_energy,tot_sigma,steps,RBS_yield,i,layer);
				}
			}
		}// end of elements
	}	
}

void EventAction::PrepareHistoBase(void)
{
	G4int x = 0;
	for(int i=0; i<5; ++i) {
		HistoBase[i][0] = 20;
		for(int j=1; j<4; j++){
			if(i==0)
				x = 21;
			else if(i==1)
				x = 24;
			else if(i==2)
				x = 27;
			else if(i==3)
				x = 30;
			else if(i==4)
				x = 33;
			else
				exit(1);
			HistoBase[i][j] = x+(j-1);
			//G4cout << " [" << i << "][" << j << "] = " << HistoBase[i][j] << G4endl;
		}
	}
}

void EventAction::FillH2ChannelingHistos(G4double z_pos, G4double ch_x_pos, G4double ch_y_pos, G4double wo_x_pos, G4double wo_y_pos)
{
	/*
	for(int i=0; i<12; i++) {
		G4cout << " i = " << i << " histo = " << depth_for_histos[i]/nm << G4endl;
	}
	*/
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->FillH2(1,wo_x_pos,z_pos);
	analysisManager->FillH2(2,wo_y_pos,z_pos);
	analysisManager->FillH2(3,ch_x_pos,z_pos);
	analysisManager->FillH2(4,ch_y_pos,z_pos);

	if(z_pos > depth_for_histos[0] && z_pos < depth_for_histos[1])
		analysisManager->FillH2(5,ch_x_pos,ch_y_pos);
	else if(z_pos > depth_for_histos[1] && z_pos < depth_for_histos[2])
		analysisManager->FillH2(6,ch_x_pos,ch_y_pos);
	else if(z_pos > depth_for_histos[2] && z_pos < depth_for_histos[3])
		analysisManager->FillH2(7,ch_x_pos,ch_y_pos);
	else if(z_pos > depth_for_histos[3] && z_pos < depth_for_histos[4])
		analysisManager->FillH2(8,ch_x_pos,ch_y_pos);
	else if(z_pos > depth_for_histos[4] && z_pos < depth_for_histos[5])
		analysisManager->FillH2(9,ch_x_pos,ch_y_pos);
	else if(z_pos > depth_for_histos[5] && z_pos < depth_for_histos[6])
		analysisManager->FillH2(10,ch_x_pos,ch_y_pos);
	else if(z_pos > depth_for_histos[6] && z_pos < depth_for_histos[7])
		analysisManager->FillH2(11,ch_x_pos,ch_y_pos);
	else if(z_pos > depth_for_histos[7] && z_pos < depth_for_histos[8])
		analysisManager->FillH2(12,ch_x_pos,ch_y_pos);
	else if(z_pos > depth_for_histos[8] && z_pos < depth_for_histos[9])
		analysisManager->FillH2(13,ch_x_pos,ch_y_pos);
	else if(z_pos > depth_for_histos[9] && z_pos < depth_for_histos[10])
		analysisManager->FillH2(14,ch_x_pos,ch_y_pos);
	else if(z_pos > depth_for_histos[10] && z_pos < depth_for_histos[11])
		analysisManager->FillH2(15,ch_x_pos,ch_y_pos);
}