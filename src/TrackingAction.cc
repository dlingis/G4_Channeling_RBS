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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
// $Id: TrackingAction.cc 69099 2013-04-18 12:25:19Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4StepStatus.hh"
#include "G4ParticleTypes.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "DetectorConstruction.hh"

#include "G4Channeling.hh"
#include "G4ChannelingTrackData.hh"
#include "G4ChannelingMaterialData.hh"


#include "ExExChParticleUserInfo.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(EventAction* event, DetectorConstruction* detector,PrimaryGeneratorAction* prim)
:G4UserTrackingAction(), fEventAction(event), fDetector(detector),fPrimary(prim)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{

	ExExChParticleUserInfo *trackInfo(static_cast< ExExChParticleUserInfo * > (track->GetUserInformation() ) );
	if (trackInfo){
		return;
	} else {
		G4Track *  theTrack( const_cast< G4Track * >( track ) );
		trackInfo = new ExExChParticleUserInfo();
		theTrack->SetUserInformation( trackInfo );
	}

	G4VPhysicalVolume* volume = track->GetTouchableHandle()->GetVolume();
	if (track->GetTrackID() == 1)
		return;
	
	for (uint8_t i=0; i<NUMB_MAX_LAYERS; ++i) {
		if (volume == fDetector->GetIntAbsorber(i)) {
			G4String name = track->GetDefinition()->GetParticleName();
			G4double energy = track->GetKineticEnergy();
			Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
			run->ParticleCount(name,energy);
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{

	Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
	G4VPhysicalVolume* volume = track->GetTouchableHandle()->GetVolume();

	G4String parName = track->GetDefinition()->GetParticleName();
	G4double Trleng = track->GetTrackLength() + fPrimary->GetParticleGun()->GetParticlePosition().z() + 0.5 * fDetector->GetLength(0);
	G4ThreeVector vertex = track->GetVertexPosition();
	G4ThreeVector position = track->GetPosition() + 0.5 * fDetector->GetDimensions(0);
	G4double z = track->GetPosition().z() + 0.5 * fDetector->GetLength(0);

	for (uint8_t i=0; i<5; ++i) {
		if (volume == fDetector->GetIntAbsorber(i)) {
			if (track->GetParentID() == 0) {
				fEventAction->AddTrueTrakLen(Trleng);
				fEventAction->AddProjTrakLen(position.z());
				//scattering angle of  primary particle
				// TODO check whether this actually works
				G4double theta = 0.0;
				G4double vz = (track->GetMomentumDirection()).z();
				//G4cout<<"Angle_Lab: "<<std::acos(vx)<<" Position_x: "<<px<<G4endl;
				if(vz <= -1.0)
					theta = -CLHEP::pi;
				else if(vz < 1.0)
					theta = std::acos(vz);
				fEventAction->AddTheta(theta);
			}
			if (track->GetTrackID() == 1) {
				fEventAction->AddProjectedRange(z);
				// run->AddProjectedRange(z);
				G4AnalysisManager* analysis = G4AnalysisManager::Instance();
				analysis->FillH1(16, z);
			}
		}
	}
	// keep only outgoing particle
	G4StepStatus status = track->GetStep()->GetPostStepPoint()->GetStepStatus();
	if (status != fWorldBoundary)
		return;

	if (track->GetParentID() == 0)
		run->AddEmerging();

	const G4ParticleDefinition* particle = track->GetParticleDefinition();
	G4String name = particle->GetParticleName();
	G4double energy = track->GetKineticEnergy();

	fEventAction->AddEflow(energy);
	run->ParticleFlux(name,energy);

	// histograms: enery flow
	G4AnalysisManager* analysis = G4AnalysisManager::Instance();

	G4int ih = 0;
	G4String type = particle->GetParticleType();
	G4double charge = particle->GetPDGCharge();
	if (charge > 3.)                             ih = 10;
	else if (particle == G4Gamma::Gamma())       ih = 4;
	else if (particle == G4Electron::Electron()) ih = 5;
	else if (particle == G4Positron::Positron()) ih = 5;
	else if (particle == G4Neutron::Neutron())   ih = 6;
	else if (particle == G4Proton::Proton())     ih = 7;
	else if (particle == G4Deuteron::Deuteron()) ih = 8;
	else if (particle == G4Alpha::Alpha())       ih = 9;
	else if (type == "nucleus")                  ih = 10;
	else if (type == "baryon")                   ih = 11;
	else if (type == "meson")                    ih = 12;
	else if (type == "lepton")                   ih = 13;
	if (ih > 0) analysis->FillH1(ih,energy);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......