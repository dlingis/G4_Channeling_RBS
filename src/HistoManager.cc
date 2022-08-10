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
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 67909 2013-03-12 18:51:09Z vnivanch $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>
#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("Hadr06")
{

G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	G4cout << "Using " << analysisManager->GetType() << G4endl;
	//** Create directories **//
    
    analysisManager->SetVerboseLevel(1);
  
    analysisManager->SetNtupleMerging(true);  //sumergina kiekvieno i bendra suma
    // Creating ntuple 
    analysisManager->CreateNtuple("sim", "Channeling");
    analysisManager->CreateNtupleDColumn("angXin");    //0
    analysisManager->CreateNtupleDColumn("angYin");    //1
    analysisManager->CreateNtupleDColumn("energy,MeV");//2
    analysisManager->CreateNtupleDColumn("posXin");    //3
    analysisManager->CreateNtupleDColumn("posYin");    //4
    analysisManager->CreateNtupleDColumn("l1_efx");    //5
    analysisManager->CreateNtupleDColumn("l1_efy");    //6
    analysisManager->CreateNtupleDColumn("l1_nud");    //7
    analysisManager->CreateNtupleDColumn("l1_eld");    //8
    analysisManager->CreateNtupleDColumn("l2_efx");    //9
    analysisManager->CreateNtupleDColumn("l2_efy");    //10
    analysisManager->CreateNtupleDColumn("l2_nud");    //11
    analysisManager->CreateNtupleDColumn("l2_eld");    //12
    analysisManager->CreateNtupleDColumn("l3_efx");    //13
    analysisManager->CreateNtupleDColumn("l3_efy");    //14
    analysisManager->CreateNtupleDColumn("l3_nud");    //15
    analysisManager->CreateNtupleDColumn("l3_eld");    //16
    analysisManager->CreateNtupleDColumn("l4_efx");    //17
    analysisManager->CreateNtupleDColumn("l4_efy");    //18
    analysisManager->CreateNtupleDColumn("l4_nud");    //19
    analysisManager->CreateNtupleDColumn("l4_eld");	   //20
    analysisManager->CreateNtupleDColumn("l5_efx");    //21
    analysisManager->CreateNtupleDColumn("l5_efy");    //22
    analysisManager->CreateNtupleDColumn("l5_nud");    //23
    analysisManager->CreateNtupleDColumn("l5_eld");    //24
    //analysisManager->FinishNtuple();
    analysisManager->CreateNtupleDColumn("l1_nud_a");  //25
    analysisManager->CreateNtupleDColumn("l1_nud_b");  //26
    analysisManager->CreateNtupleDColumn("l1_nud_c");  //27
    analysisManager->CreateNtupleDColumn("l2_nud_a");  //28
    analysisManager->CreateNtupleDColumn("l2_nud_b");  //29
    analysisManager->CreateNtupleDColumn("l2_nud_c");  //30
    analysisManager->CreateNtupleDColumn("l3_nud_a");  //31
    analysisManager->CreateNtupleDColumn("l3_nud_b");  //32
    analysisManager->CreateNtupleDColumn("l3_nud_c");  //33
    analysisManager->CreateNtupleDColumn("l4_nud_a");  //34
    analysisManager->CreateNtupleDColumn("l4_nud_b");  //35
    analysisManager->CreateNtupleDColumn("l4_nud_c");  //36
    analysisManager->CreateNtupleDColumn("l5_nud_a");  //37
    analysisManager->CreateNtupleDColumn("l5_nud_b");  //38
    analysisManager->CreateNtupleDColumn("l5_nud_c");	 //39
    //
    analysisManager->CreateNtupleDColumn("posXout");	 //40
    analysisManager->CreateNtupleDColumn("posYout");	 //41
    analysisManager->CreateNtupleDColumn("angXout");	 //42
    analysisManager->CreateNtupleDColumn("angYout");	 //43
    analysisManager->FinishNtuple();
	Book();
	BookH2();
	BookP1();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  //delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
	// Create or get analysis manager
	// The choice of analysis technology is done via selection of a namespace
	// in HistoManager.hh
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetDefaultFileType("root");
	analysisManager->SetFileName(fFileName);
	//analysisManager->SetFirstHistoId(1);
	analysisManager->SetVerboseLevel(1);
	analysisManager->SetActivation(true);     //enable inactivation of histograms

	// Define histograms start values
	const G4int kMaxHisto = 37;
	const G4String id[] = {"0",
  		"Total_energy_deposit",
  		"EDEP_(MeV/mm)",
  		"3",
  		"4",
  		"5",
  		"6",
  		"7",
  		"8",
  		"9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "Projected_range",
        "17",
        "18",
		"19",
		"RBS_total_energy",		//20
		"RBS_substrate_energy",	//21
		"RBS_substrate_el1_en",	//22
		"RBS_substrate_el2_en",	//23
		"RBS_la1_en",			//24
		"RBS_la1_el1_en",		//25
		"RBS_la1_el2_en",		//26
		"RBS_la2_en",			//27
		"RBS_la2_el1_en",		//28
		"RBS_la2_el2_en",		//29
		"RBS_la3_en",			//30
		"RBS_la3_el1_en",		//31
		"RBS_la3_el2_en",		//32
		"RBS_la4_en",			//33
		"RBS_la4_el1_en",		//34
		"RBS_la4_el2_en",		//35
		"energy_at_detector"	//36
};
	const G4String title[] = { 
    	"dummy",										//0
    	"total energy deposit",							//1
    	"Edep (MeV/mm) along absorber",					//2
    	"total kinetic energy flow",					//3
    	"gamma flux (dN/dE) at exit",					//4
    	"e+- flux (dN/dE) at exit",						//5
    	"neutrons flux (dN/dE) at exit",				//6
    	"protons flux (dN/dE) at exit",					//7
    	"deuterons flux (dN/dE) at exit",				//8
    	"alphas flux (dN/dE) at exit",					//9
    	"all others ions flux (dN/dE) at exit",			//10
    	"all others baryons flux (dN/dE) at exit",		//11
    	"all others mesons flux (dN/dE) at exit",		//12
    	"all others leptons flux (dN/dE) at exit",		//13 
 		"Frenkel pairs(#/cm3)",							//14
		" Sklaidos kampas ",							//15
		" Projected range distribution ",				//16	
		" NIEL by distance, keV	",						//17 
		" NIEL by energy",								//18
		" EDEP by energy ",								//19
		" Total RBS Spectrum ",							//20
		" RBS spectrum substrate en",					//21
		" RBS spectrum substrate el 1",					//22
		" RBS spectrum substrate el 2",					//23
		" RBS Layer 1 en",								//24
		" RBS Layer 1 el 1 ",							//25
		" RBS Layer 1 el 2 ",							//26
		" RBS Layer 2 en",								//27
		" RBS Layer 2 el 1 ",							//28
		" RBS Layer 2 el 2 ",							//29
		" RBS Layer 3 en",								//30
		" RBS Layer 3 el 1 ",							//31
		" RBS Layer 3 el 2 ",							//32	  
		" RBS Layer 4 en",								//33
		" RBS Layer 4 el 1 ",							//34
		" RBS Layer 4 el 2 ",							//35
    	"energy_at_detector"							//36
		};  

	// Default values (to be reset via /analysis/h1/set command) 
	G4int nbins = 100;
	G4double vmin = 0.;
	G4double vmax = 100.;

	// Create all histograms as inactivated 
	// as we have not yet set nbins, vmin, vmax
	for (G4int k=0; k<kMaxHisto; k++) {
		G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
		analysisManager->SetH1Activation(ih, false);
	}
}

void HistoManager::BookH2()
{
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	// Define histograms start values
	const G4int kMaxHisto = 16;
	const G4String id[] = {"H2_0",
		"H2_1",
		"H2_2",
		"H2_3",
		"H2_4",
		"H2_5",
		"H2_6",
		"H2_7",
		"H2_8",
		"H2_9",
		"H2_10",
		"H2_11",
		"H2_12",
		"H2_13",
		"H2_14",
		"H2_15"
	};

	const G4String title[] = { 
		"dummy",							//0
		"XZ position World system",			//1
		"YZ position World system",			//2
		"XZ position Crystal system",		//3
		"YZ position Crystal system",		//4
		"XY position Crystal interval 1",	//5
		"XY position Crystal interval 2",	//6
		"XY position Crystal interval 3",	//7
		"XY position Crystal interval 4",	//8
		"XY position Crystal interval 5",	//9
		"XY position Crystal interval 6",	//10
		"XY position Crystal interval 7",	//11
		"XY position Crystal interval 8",	//12
		"XY position Crystal interval 9",	//13
		"XY position Crystal interval 10",	//14
		"XY position Crystal interval 11"	//15
	};

	G4int nbinsx = 100;
	G4double xmin = 0;
	G4double xmax = 100;
	G4int nbinsy = 100;
	G4double ymin = 0;
	G4double ymax = 100;
	for(int k=0; k<kMaxHisto; ++k) {
		G4int ih = analysisManager->CreateH2(id[k],title[k],nbinsx,xmin,xmax,nbinsy,ymin,ymax);
		analysisManager->SetH2Activation(ih,false);
	}
}

void HistoManager::BookP1()
{
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	// Define histograms start values
	const G4int kMaxHisto = 11;
	const G4String id[] = {"P1_0",
		"P1_1",
		"P1_2",
		"P1_3",
		"P1_4",
		"P1_5",
		"P1_6",
		"P1_7",
		"P1_8",
		"P1_9",
		"P1_10"
	};

	const G4String title[] = { 
		"dummy",						//0
		"NUD sum along depth",			//1
		"ELD sum along depth",			//2
		"NUD layer1 along depth",		//3
		"ELD layer1 along depth",		//4
		"NUD layer2 along depth",		//5
		"ELD layer2 along depth",		//6
		"NUD layer3 along depth",		//7
		"ELD layer3 along depth",		//8
		"NUD layer4 along depth",		//9
		"ELD layer4 along depth"		//10
	};

	G4int nbins = 100;
	G4double min = 0;
	G4double max = 100;

	for(int k=0; k<kMaxHisto; ++k) {
		G4int ih = analysisManager->CreateP1(id[k],title[k],nbins,min,max);
		analysisManager->SetP1Activation(ih,false);
	}

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......