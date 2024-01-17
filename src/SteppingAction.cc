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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 71404 2013-06-14 16:56:38Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* event, DetectorConstruction* DET)
: G4UserSteppingAction(), fEventAction(event), detector(DET)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
	//safety
	G4double edep = aStep->GetTotalEnergyDeposit();
	if (edep <= 0.)
		return;
	G4VPhysicalVolume* preMAT = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	G4VPhysicalVolume* posMAT = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
	G4double charge = aStep->GetTrack()->GetDynamicParticle()->GetCharge();
	G4String parName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	G4ThreeVector end = aStep->GetPostStepPoint()->GetPosition();
	G4ThreeVector start = aStep->GetPreStepPoint()->GetPosition();
	G4double depth = (start + G4UniformRand() * (end - start)).z() + (detector->GetLength(0) / 2);
	// primary particles 
	G4int IDp =  aStep->GetTrack()->GetParentID();
	Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

	for (G4int i=0; i<NUMB_MAX_LAYERS; ++i) {
		if (preMAT == detector->GetIntAbsorber(i) && posMAT == detector->GetIntAbsorber(i)) {
			if (IDp == 0) {
				fEventAction->AddTrackLenLayer(aStep->GetStepLength(), i);
				fEventAction->AddIonDepLayer(aStep->GetTotalEnergyDeposit(), i);
				fEventAction->AddNonDepLayer(aStep->GetNonIonizingEnergyDeposit(), i);
				fEventAction->AddStepLayer(i);
				// totals
				fEventAction->AddTrakLenPrim(aStep->GetStepLength());
				fEventAction->AddPrimSteps();
				G4double edep = aStep->GetTotalEnergyDeposit();
				G4double niel = aStep->GetNonIonizingEnergyDeposit();
				fEventAction->AddEdep(edep);
				fEventAction->AddNiel(niel);

				if (niel > 0.) {
					// NIEL
					analysisManager->FillH1(17, depth, niel /keV);
					analysisManager->FillH1(18, niel);
					// IONiZING energy loss
					analysisManager->FillH1(2, depth, edep);
					analysisManager->FillH1(19, edep);
				}
			}
			// secondaries
			if (IDp > 0 && charge > 0) {
				fEventAction->AddTrakLenSec(aStep->GetStepLength());
			}
		}
	}

	const G4StepPoint* endPoint = aStep->GetPostStepPoint();
	const G4VProcess* process   = endPoint->GetProcessDefinedStep();
	run->CountProcesses(process);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


