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
/// \file Hadr06.cc
/// \brief Main program of the hadronic/Hadr06 example
//
//
// $Id: TestEm1.cc,v 1.16 2010-04-06 11:11:24 maire Exp $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4RunManagerFactory.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "SteppingVerbose.hh"

#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BIC.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

// channe
#include "G4Types.hh"
#include "G4ScoringManager.hh"
#include "G4EmStandardPhysics_option4_channeling.hh"
#include "G4ChannelingPhysics.hh"


#include "G4GenericBiasingPhysics.hh"
//#include "G4HadronElasticPhysics.hh"
//#include "G4HadronPhysicsFTFP_BERT.hh"

#include "G4EmStandardPhysics_option4.hh"

// new 2022-03-02
#include "HadronElasticPhysicsHP.hh"
#include "G4IonElasticPhysics.hh"
#include "G4IonPhysicsXS.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
int main(int argc,char** argv) {
 
   //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);
 
  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
 
  //construct the run manager
  auto runManager = G4RunManagerFactory::CreateRunManager();  
  runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores() - 4);

    // Set mandatory initialization classes
    //G4VModularPhysicsList* physlist= new FTFP_BERT();
    G4VModularPhysicsList* physlist= new QGSP_BIC();
    physlist->RegisterPhysics(new G4ChannelingPhysics());
 // new addition 2022-03-02
    //physlist->RegisterPhysics(new HadronElasticPhysicsHP(1));
    //physlist->RegisterPhysics(new G4IonElasticPhysics(1));
    //physlist->RegisterPhysics(new G4IonPhysicsXS(1));
    //
    physlist->ReplacePhysics(new G4EmStandardPhysics_option4_channeling());
    runManager->SetUserInitialization(physlist);

/*
    // Set mandatory initialization classes
    G4VModularPhysicsList* physlist= new FTFP_BERT();
    G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
    physlist->RegisterPhysics(new G4ChannelingPhysics());
    physlist->ReplacePhysics(new G4EmStandardPhysics_option4());
    biasingPhysics->PhysicsBiasAllCharged();
    physlist->RegisterPhysics(biasingPhysics);
    runManager->SetUserInitialization(physlist);
*/
  // set mandatory initialization classes
  DetectorConstruction* det= new DetectorConstruction;
  runManager->SetUserInitialization(det);
    
  runManager->SetUserInitialization(new ActionInitialization(det));    
     
  //initialize visualization
  G4VisManager* visManager = nullptr;

  //get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui)  {
   //interactive mode
   visManager = new G4VisExecutive;
   visManager->Initialize();
   ui->SessionStart();
   delete ui;
  }
  else  {
   //batch mode
   G4String command = "/control/execute ";
   G4String fileName = argv[1];
   UImanager->ApplyCommand(command+fileName);
  }

  //job termination
  delete visManager;
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
