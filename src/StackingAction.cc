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
// StackingAction.cc
// 
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include <iostream>
#include <iomanip>

#include "RunAction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "StackingMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Material.hh"
#include "G4RunManager.hh"

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(EventAction* EA, DetectorConstruction* DE )
:eventaction(EA), detector(DE)
{
	killSecondary = 0;
	stackMessenger = new StackingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
	delete stackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// old Lindhard - Robinson partition function for protons
G4double StackingAction::DamageEnergy(G4double T,G4double A, G4double Z)
{
	//.................. T in eV!!!!!!!!!!!!!
	G4double Z2 = Z;
	G4double M2 = A;
	G4double k_d, epsilon_d, g_epsilon_d, E_nu;

	// original, for protons
	k_d = 0.1334 * std::pow(Z2, (2. / 3.)) * std::pow(M2, (-1. / 2.));
	epsilon_d = 0.01014 * std::pow(Z2, (-7. / 3.)) * (T/eV);
	g_epsilon_d = epsilon_d + 0.40244 * std::pow(epsilon_d, (3. / 4.)) + 3.4008 * std::pow(epsilon_d, (1. / 6.));
	E_nu=1. / (1.+ k_d * g_epsilon_d);
	return E_nu;//partition fraction!!!
}

//...................................................................

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
	G4int IDp= aTrack->GetParentID();
	G4double charge = aTrack->GetDefinition()->GetPDGCharge();

	//keep primary particle
	if (IDp == 0)
		return fUrgent;
	Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

	G4double partition;
	if (IDp > 0 && charge > 0) {
		G4double energy = aTrack->GetKineticEnergy();
		G4double u = 931.49410242;
		G4double A1 = aTrack->GetParticleDefinition()->GetPDGMass()/u;
		G4double Z1 = aTrack->GetParticleDefinition()->GetPDGCharge();
		G4double Spoint  = (aTrack->GetPosition()).mag();
		G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		// MAIN LAYER
		for (uint8_t i=0; i<NUMB_MAX_LAYERS; ++i) {
			if (aTrack->GetTouchableHandle()->GetVolume() == detector->GetIntAbsorber(i)) {
				G4Material* material  = detector->GetMaterialM(i);
				G4double mA2 = material->GetDensity() / (material->GetTotNbOfAtomsPerVolume() *amu);
				G4double mZ2 = material->GetTotNbOfElectPerVolume() / material->GetTotNbOfAtomsPerVolume();
				partition = DamagePartitionFunction(energy, A1, mA2, Z1, mZ2);
				run->absSumTLLayer(partition * energy, i); // damage energy
				run->absRecLayer(1, i); // recoil
				run->absSumTLayer(energy, i); // kinetic energy
			}
		}
		analysisManager->FillH1(14, Spoint, partition * energy);
		run->NumberRec(1);
	}

	//stack or delete secondaries
	G4ClassificationOfNewTrack status = fUrgent;
	if (killSecondary == 1)
		status = fKill;

	return status;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double StackingAction::DamagePartitionFunction(G4double T, G4double M1, G4double M2, G4double Z1, G4double Z2)
{
	G4double E_nu, a;
	G4double e_ch_squared = 1.4399764e-6;//*(MeV*nm);
	//screening distance
	G4double a_0 = 0.0529*nm; // Bohr radius
	a = 0.8853 * a_0 / (std::pow(std::pow(Z1, (2. / 3.)) + std::pow(Z2, (2. / 3.)), (1. / 2.))) /nm; // Thomas-Fermi
	//G4double x = 0.23;
	//a = 0.8853 * a_0 / (std::pow(Z1, x) + std::pow(Z2, x)); // ZBL
	G4double E_l = (Z1 * Z2 * e_ch_squared / a) * (M1 + M2) / M2;
	G4double first_term = std::pow(M1 + M2, (3. / 2.)) * std::pow(Z1, (2. / 3.)) * std::pow(Z2, (1. / 2.));
	G4double second_term = std::pow(M1, (3. / 2.)) * std::pow(std::pow(Z1, (2. / 3.)) + std::pow(Z2, (2. / 3.)), (3. / 4.));
	G4double m_0 = 0.510998950; // electron rest mass
	G4double M2_x = 931.49410242 * M2; // target atom mass
	G4double k_l = (32 / (3 * pi)) * std::pow((m_0 / M2_x), (1. / 2.)) * first_term / second_term;
	G4double T_1 = T / eV; // kinetic energy of projectile in eV
	G4double T_El = T_1 / (E_l * 1000000);

	// Lindhard - Robinson partition function for any ion
	//G4double g_l = T_El + 0.40244 * std::pow(T_El, (3. / 4.)) + 3.4008 * std::pow(T_El, (1. / 6.));
	// Akerman partition function https://ieeexplore.ieee.org/document/4033183
	G4double g_l = 0.74422 * T_El + 1.6812 * std::pow(T_El, (0.75)) + 0.90565 * std::pow(T_El, (1./6.));
	E_nu = 1. / (1. + k_l * g_l);
	return E_nu;
}