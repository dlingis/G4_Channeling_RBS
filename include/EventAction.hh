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
/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 76293 2013-11-08 13:11:23Z gcosmo $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4Physics2DVector.hh"
#include "G4ElementVector.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UAtomicDeexcitation.hh"

#define NUMBER_OF_LOCATIONS    10
#define NUMBER_OF_MAX_LAYERS    5


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class RunAction;
class DetectorConstruction; 
class PrimaryGeneratorAction;

class EventAction : public G4UserEventAction
{
	public:
		EventAction(DetectorConstruction*, PrimaryGeneratorAction*);
		~EventAction();

	public:
		virtual void BeginOfEventAction(const G4Event*);
		virtual void EndOfEventAction(const G4Event*);
		
		void AddEdep (G4double Edep)              {fTotalEnergyDeposit += Edep;};
		void AddEflow(G4double Eflow)             {fTotalEnergyFlow += Eflow;};
		void AddTheta(G4double tet)               {theta += tet;};
		void AddNiel (G4double niel)              {fTotalNiel += niel;};
		void AddStepLayer(G4int layer)            {StepsPrim[layer]++;};
		void AddTrueTrakLen(G4double trueLength)  {TrueTrakLen += trueLength;};
		void AddProjTrakLen(G4double projLength)  {ProjTrakLen += projLength;};
		void AddProjectedRange(G4double l)        {projected_range += l;};
		void AddTrackLenLayer(G4double s, G4int l){step_layer[l] += s;};
		void AddIonDepLayer(G4double e, G4int l)  {ion_dep_layer[l] += e;};
		void AddNonDepLayer(G4double n, G4int l)  {non_dep_layer[l] += n;};
		void AddTrakLenPrim(G4double l)           {trak_len_prim += l;};
		void AddPrimSteps(void)                   {prim_step++;};
		void AddTrakLenSec(G4double l)            {trak_len_sec += l;};
		// RBS
		G4Physics2DVector* Get2DRTRVector(G4String element, G4int Z);
		G4double Get2DRTRValue(G4double energy, G4String elname, G4double angle);
		// distribution functions
		G4double GenerateGaussian(G4double x, G4double y, G4double sigma_sq);
		G4double* CalcEnergyLeft(G4double depth, G4double energy);
		void PrepareDistances();
		G4double CalcLossInLayer(G4int i, G4double dist, G4double energy, G4ParticleDefinition* particle);
		G4double CalcStragglingInLayer(G4int i, G4double dist, G4double energy, G4ParticleDefinition* particle);
		G4double GetGaussSum(G4double energy, G4double sigma);
		void FillGaussians(G4double energy, G4double sigma, G4double steps, G4double yield, G4int element, G4int layer);
		void CalculateRBS(G4int layer, G4double energy, G4ThreeVector position, G4ThreeVector momDir, G4double steps, G4ParticleDefinition* fParticle);
		void PrepareHistoBase(void);
		void FillH2ChannelingHistos(G4double z_pos, G4double x_pos, G4double y_pos, G4double wo_x_pos, G4double wo_y_pos);
		void FillSigmaCalcVectors(void);

	private:
		DetectorConstruction* detector;
		PrimaryGeneratorAction* fPrimary;

		G4double fTotalEnergyDeposit, fTotalEnergyFlow, fTotalNiel, theta, prim_step, trak_len_prim, trak_len_sec;
		G4double step_layer[5], ion_dep_layer[5], non_dep_layer[5];


		G4int StepsPrim[5];
		G4double TrueTrakLen, ProjTrakLen, En, Angle, Ma1, Ma2, projected_range;
		G4int sdht_ID, sdct_ID, sdxt_ID, sd0_ID,sd1_ID,sd2_ID,sd3_ID,sd4_ID;
		//2D vectors of rtr values
		G4Physics2DVector* fVectorSi_total;
		G4Physics2DVector* fVectorO_total;
		G4Physics2DVector* fVectorNa_total;
		G4Physics2DVector* fVectorN_total;
		G4Physics2DVector* fVectorC_total;
		G4Physics2DVector* fVectorF_total;
		G4Physics2DVector* fVectorB_total;
		G4Physics2DVector* fVectorNi_total;
		G4Physics2DVector* fVectorCu_total;
		
		const G4ElementVector* ElementVector;
		const G4Element* Element;
		const G4Isotope* Isotope;

		G4UAtomicDeexcitation* deExcitation = new G4UAtomicDeexcitation();

		G4Material* sample_material[5];
		
		G4int* NoOfElements         = new G4int [5];
		G4double** Znumb            = new G4double* [5];
		G4double** Mnumb            = new G4double* [5];
		G4double** Adens            = new G4double* [5];
		G4int** HistoBase           = new G4int* [5];
		G4double** nud_el           = new G4double* [5];
	
		G4double* thicknesses       = new G4double[4];
		G4double* positions         = new G4double[4];
		G4double* layer_start_pos   = new G4double[4];
		G4double* layer_end_pos     = new G4double[4];
		G4double* dist_between_l    = new G4double[4];

		G4double* dist_matrix       = new G4double[NUMBER_OF_LOCATIONS];
		G4double angle_of_exit;
		G4double* start_pos         = new G4double[NUMBER_OF_LOCATIONS];
		G4Material* mat_matrix[NUMBER_OF_LOCATIONS];
		G4double  depth_for_histos[12];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif