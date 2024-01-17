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
/// \file Run.hh
/// \brief Definition of the Run class
//
// $Id: Run.hh 71375 2013-06-14 07:39:33Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include "G4THitsMap.hh"
#include <map>

#define NUMB_MAX_LAYERS          5

//test
class DetectorConstruction;
class G4ParticleDefinition;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
	public:
		Run(DetectorConstruction*);
		~Run();

	public:
	void SetPrimary(G4ParticleDefinition* particle, G4double energy);
	void CountProcesses(const G4VProcess* process);
	void ParticleCount(G4String, G4double); 
	void AddEdep (G4double edep) {fEnergyDeposit += edep; fEnergyDeposit2 += edep * edep;};
	void AddEflow (G4double eflow) {fEnergyFlow += eflow; fEnergyFlow2 += eflow * eflow;};
	void ParticleFlux(G4String, G4double);
	void AddNonIonEnergy (G4double enondep) {NonIonEnergyDeposit += enondep; NonIonEnergyDeposit2 += enondep * enondep;};
	void NumberRec (G4int i) {N_rec += i;};
	void AddNumberOfSteps (G4int i) {Nsteps += double(i); Nsteps2 += double(i) * double(i);};
	void AddTheta (G4double tet) {theta += tet; theta2 += tet * tet;};
	void AddTrakLenPrim (G4double length) {TrakLenPrim += length; TrakLenPrim2 += length * length;};
	void AddTrakLenSec (G4double length) {TrakLenSec += length; TrakLenSec2 += length * length;};
	void AddDetHits() {hit++;}; // sums detector hits
	void AddEmerging() {partEmerging++;}; // sums particles leaving material
	void AddKinEn(G4double energy) { detKinEn += energy; detKinEn2 += energy*energy;}; // kinetic energy on 4th detector

	void AddProjectedRange(G4double z) {projectedR += z; projectedR2 += z * z;};
	void PrimaryEnergy(G4double en) {primary_energy += en;};
	virtual void Merge(const G4Run*);
	void EndOfRun();

	void AddTrueRange (G4double l) { fTrueRange += l; fTrueRange2 += l * l;};
	void AddProjRange (G4double x) { fProjRange += x; fProjRange2 += x * x;};
	void absTrackLenLayer (G4double len, G4int l) {trackLenLayer[l][0] += len; trackLenLayer[l][1] += (len * len);};
	void absStepLayer (G4int step, G4int l)       {stepsLayer[l][0] += step; stepsLayer[l][1] += (step * step);};
	void absIonLayer (G4double edep, G4int l)     {edepLayer[l][0] += edep; edepLayer[l][1] += (edep * edep);};
	void absNonLayer (G4double niel, G4int l)     {nielLayer[l][0] += niel; nielLayer[l][1] += (niel * niel);};
	void absRecLayer (G4int rec, G4int l)         {numbRecLayer[l] += rec;};
	void absSumTLayer (G4double T, G4int l)       {secKinEnLayer[l][0] += T; secKinEnLayer[l][1] += (T * T);};
	void absSumTLLayer (G4double TL, G4int l)     {secDamEnLayer[l][0] += TL; secDamEnLayer[l][1] += (TL * TL);};

	G4double GetVariation(G4double val2, G4double val1)
	{
		return std::sqrt(std::fabs(val2 - (val1 * val1)));
	}
	// check for entries
	void add_entry_sd(G4int en) {entry_sd += en;};
	void add_total_step(G4double en) {total_step += en;};
	void MaxRBSDepth(G4double dist) {RBSDepth += dist; RBSDepth2 += dist*dist;};
	void AddCount() {counts++;};
	void Inc_angle(G4double a) {angle_of_incidence += a;};
	void add_entry_reach(G4double en) {entry_reach += en;};


	private:
		struct ParticleData {
			ParticleData()
			: fCount(0), fEmean(0.), fEmin(0.), fEmax(0.) {}
			ParticleData(G4int count, G4double ekin, G4double emin, G4double emax)
			: fCount(count), fEmean(ekin), fEmin(emin), fEmax(emax) {}
			G4int     fCount;
			G4double  fEmean;
			G4double  fEmin;
			G4double  fEmax;
		};

	private:
		DetectorConstruction*              fDetector;
		G4ParticleDefinition*              fParticle;
		G4double                           fEkin, primary_energy;
		
		G4double                           fEnergyDeposit, fEnergyDeposit2;
		G4double                           fEnergyFlow, fEnergyFlow2;
		std::map<G4String,G4int>           fProcCounter;
		std::map<G4String,ParticleData>    fParticleDataMap1;
		std::map<G4String,ParticleData>    fParticleDataMap2;

		G4double                           fTrueRange, fTrueRange2, fProjRange, fProjRange2;

		G4double trackLenLayer[NUMB_MAX_LAYERS][3]; // 0 - sum of vals, 1 - sum of squares, 2 - std dev
		G4double edepLayer[NUMB_MAX_LAYERS][3];
		G4double nielLayer[NUMB_MAX_LAYERS][3];
		G4int stepsLayer[NUMB_MAX_LAYERS][3];
		G4double secKinEnLayer[NUMB_MAX_LAYERS][3];
		G4double secDamEnLayer[NUMB_MAX_LAYERS][3];
		G4int numbRecLayer[NUMB_MAX_LAYERS];

		G4double RBSDepth, RBSDepth2;
		G4int counts, N_rec, entry_sd, hit, partEmerging;
		G4double NonIonEnergyDeposit,  NonIonEnergyDeposit2;
		G4double Th, Nsteps, Nsteps2, theta, theta2;

		G4double TrakLenPrim, TrakLenPrim2;
		G4double TrakLenSec, TrakLenSec2;
		G4double entry_reach, total_step;
		G4double detKinEn, detKinEn2, projectedR, projectedR2, angle_of_incidence;


		void Print(const std::vector<G4String>& title,
				const std::map< G4int, std::vector<G4double> >&out,
				G4String&) const;
	
		std::map<G4int, G4THitsMap<G4double>* > fMap;
		G4String fOutputFileSpec;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif