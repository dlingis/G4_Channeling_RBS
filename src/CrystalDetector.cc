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

#include "CrystalDetector.hh"
#include "CrystalDetectorHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"
#include "G4ChannelingTrackData.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"


CrystalDetector::CrystalDetector(G4String name,const G4String& fileName1,const G4String& fileName2,const G4String& fileName3,G4int amorphous):
//CrystalDetector::CrystalDetector(G4String name):
G4VSensitiveDetector(name){
    G4String HCname;
    collectionName.insert(HCname="collection");
    fHCID = -1;
    //fChannelingID = -1;

    if(!amorphous) {
        if (fileName1 != "") {
            G4String fileAtD_a = fileName1 + "_atd.txt";
            G4cout << "\n Added atomic density for sensitive detector " << name << " from file " << fileName1 << G4endl;
            fNucleiDensity_a =   new G4ChannelingECHARM(fileAtD_a,1.);
            set_name1 = 1; 
        } else {
            G4cout << "\n No atomic density for sensitive detector " << name << G4endl;
            set_name1 = 0;
        }

        if (fileName2 != "") {
            G4String fileAtD_b = fileName2 + "_atd.txt";
            G4cout << "\n Added atomic density for sensitive detector " << name << " from file " << fileName2 << G4endl;
            fNucleiDensity_b =   new G4ChannelingECHARM(fileAtD_b,1.);
            set_name2 = 1;
        } else {
            set_name2 = 0;
            G4cout << "\n No atomic density for sensitive detector " << name << G4endl;	
        }
        if (fileName3 != "") {
            G4String fileAtD_c = fileName3 + "_atd.txt";
            G4cout << "\n Added atomic density for sensitive detector " << name << " from file " << fileName3 << G4endl;
            fNucleiDensity_c =   new G4ChannelingECHARM(fileAtD_c,1.);
            set_name3 = 1;
        } else {
            set_name3 = 0;
            G4cout << "\n No atomic density for sensitive detector " << name << G4endl;	
        }
    } else {
        set_name1 = 0;
        set_name2 = 0;
        set_name3 = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

CrystalDetector::~CrystalDetector(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CrystalDetector::Initialize(G4HCofThisEvent*HCE){
    fHitsCollection =
        new CrystalDetectorHitsCollection(SensitiveDetectorName,
                                                        collectionName[0]);
    if(fHCID<0){
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
    }
    HCE->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool CrystalDetector::ProcessHits(G4Step*aStep,
                                                  G4TouchableHistory*){

    if(aStep->GetTrack()->GetTrackID()>1) return true;

    G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();


    if((postStepPoint->GetStepStatus() == fGeomBoundary)) return true;
    if((preStepPoint->GetStepStatus() == fGeomBoundary)) return true;
    
    fChannelingID = G4PhysicsModelCatalog::GetModelID("model_channeling");

    G4ChannelingTrackData* trackdata =
    (G4ChannelingTrackData*)(aStep->GetTrack()->GetAuxiliaryTrackInformation(fChannelingID));

    G4TouchableHistory* theTouchable =
        (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4VPhysicalVolume* thePhysical = theTouchable->GetVolume(0);
    G4int copyNo = thePhysical->GetCopyNo();

    G4ThreeVector momWorld = preStepPoint->GetMomentum();
    G4ThreeVector momDirect = preStepPoint->GetMomentumDirection();
    G4ThreeVector worldPos = preStepPoint->GetPosition();
    
    CrystalDetectorHit* aHit =
        new CrystalDetectorHit(copyNo);
    aHit->SetLayerID(copyNo);
    
    G4ThreeVector pos = trackdata->GetPosCh();
    
    if(aStep->GetTrack()->GetKineticEnergy()/MeV > 0.1) // new addition 2020-09-17
    {
        aHit->SetWorldMomentumDirection(momDirect);
        aHit->SetWorldPos(worldPos);
        aHit->SetStep(aStep->GetStepLength());
        aHit->SetKinECR(aStep->GetTrack()->GetKineticEnergy());
        aHit->SetWorldMomentum(momWorld);
        aHit->SetEFX(trackdata->GetEFX());
        aHit->SetEFY(trackdata->GetEFY());
        aHit->SetNud(trackdata->GetNuD());
        aHit->SetEld(trackdata->GetElD());
        aHit->SetChPos(pos);
        if (set_name1) 
    	    aHit->SetNuD_a(fNucleiDensity_a->GetEC(pos));	
	    else
		    aHit->SetNuD_a(1.);
	    if (set_name2) 
    	    aHit->SetNuD_b(fNucleiDensity_b->GetEC(pos));	
	    else
		    aHit->SetNuD_b(1.);
	    if (set_name3) 
    	    aHit->SetNuD_c(fNucleiDensity_c->GetEC(pos));	
	    else
		    aHit->SetNuD_c(1.);

        fHitsCollection->insert(aHit);
    }



    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void CrystalDetector::EndOfEvent(G4HCofThisEvent* /*HCE*/){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
