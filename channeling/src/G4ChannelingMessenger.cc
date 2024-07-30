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

#include "G4ChannelingMessenger.hh"
#include "G4Channeling.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"

#include "G4RunManager.hh"
#include "G4ios.hh"

G4ChannelingMessenger::
G4ChannelingMessenger(G4Channeling* mpgax)
:fTarget(mpgax)
{
	fChannelingDirectory = new G4UIdirectory("/chan/",true);
	fChannelingDirectory->SetGuidance("Channeling Process control commands.");
	
	fMinimumEnergyCmd = new G4UIcmdWithADoubleAndUnit("/chan/set_min_energy",this);
	fMinimumEnergyCmd->SetGuidance("Set minimum energy allowed.");
	fMinimumEnergyCmd->SetParameterName("min_energy", true);
	fMinimumEnergyCmd->SetDefaultValue(1.);
	fMinimumEnergyCmd->SetDefaultUnit("MeV");
	
	fMaximumMomentumRatioCmd = new G4UIcmdWithADouble("/chan/set_max_mom_ratio",this);
	fMaximumMomentumRatioCmd->SetGuidance("Set maximum momentum ratio allowed.");
	fMaximumMomentumRatioCmd->SetParameterName("max_mom_ratio", true);
	fMaximumMomentumRatioCmd->SetDefaultValue(0.01);

	fStepCmd = new G4UIcmdWithADouble("/chan/set_channeling_step",this);
	fStepCmd->SetGuidance("Set channeling simulation step size.");
	fStepCmd->SetParameterName("max_step",true);
	fStepCmd->SetDefaultValue(0.1);

	fRechannelingCmd = new G4UIcmdWithABool("/chan/set_allow_rechanneling",this);
	fRechannelingCmd->SetGuidance("Allow rechanneling to occur");
	fRechannelingCmd->SetParameterName("allow_rechannel",true);
	fRechannelingCmd->SetDefaultValue(false);

	fNudEldLimCmd = new G4UIcmdWithADouble("/chan/set_density_limit",this);
	fNudEldLimCmd->SetGuidance("Limit min limit of density");
	fNudEldLimCmd->SetParameterName("limit_density",true);
	fNudEldLimCmd->SetDefaultValue(0.01);

	fGeantV11AlgorithmCmd = new G4UIcmdWithABool("/chan/v11_algorithm",this);
	fGeantV11AlgorithmCmd->SetGuidance("Use GEANT4 V11 channeling algorithm");
	fGeantV11AlgorithmCmd->SetParameterName("v11_algo",true);
	fGeantV11AlgorithmCmd->SetDefaultValue(false);

	fOrgChanStepSizeCmd = new G4UIcmdWithABool("/chan/org_step_size",this);
	fOrgChanStepSizeCmd->SetGuidance("Use org chan step size");
	fOrgChanStepSizeCmd->SetParameterName("org_step_size",true);
	fOrgChanStepSizeCmd->SetDefaultValue(false);

	fManualStepSizeCmd = new G4UIcmdWithABool("/chan/use_step_size_unit",this);
	fManualStepSizeCmd->SetGuidance("Use defined step size");
	fManualStepSizeCmd->SetParameterName("use_step_size",true);
	fManualStepSizeCmd->SetDefaultValue(false);

	fManualStepSizeValue = new G4UIcmdWithADoubleAndUnit("/chan/set_channeling_step_with_unit",this);
	fManualStepSizeValue->SetGuidance("Set channeling step with unit ");
	fManualStepSizeValue->SetParameterName("step_size",true);
	fManualStepSizeValue->SetDefaultValue(1.);
	fManualStepSizeValue->SetDefaultUnit("nm");

	fMeanFPSizeValue = new G4UIcmdWithADoubleAndUnit("/chan/set_channeling_mean_free_path",this);
	fMeanFPSizeValue->SetGuidance("Set channeling mfp size with unit ");
	fMeanFPSizeValue->SetParameterName("mfp_size",true);
	fMeanFPSizeValue->SetDefaultValue(1.);
	fMeanFPSizeValue->SetDefaultUnit("nm");

	fPrintDebugInfo = new G4UIcmdWithABool("/chan/print_debug",this);
	fPrintDebugInfo->SetGuidance("Print debug info");
	fPrintDebugInfo->SetParameterName("debug_info",true);
	fPrintDebugInfo->SetDefaultValue(false);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ChannelingMessenger:: 
~G4ChannelingMessenger(){
	delete fMinimumEnergyCmd;
	delete fMaximumMomentumRatioCmd;
	delete fStepCmd;
	delete fRechannelingCmd;
	delete fNudEldLimCmd;
	delete fGeantV11AlgorithmCmd;
	delete fOrgChanStepSizeCmd;
	delete fManualStepSizeCmd;
	delete fManualStepSizeValue;
	delete fMeanFPSizeValue;
	delete fChannelingDirectory;
	delete fPrintDebugInfo;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ChannelingMessenger::SetNewValue(G4UIcommand *command, G4String newValue){
	if (command == fMinimumEnergyCmd)
		fTarget->SetMinimumEnergy(fMinimumEnergyCmd->GetNewDoubleValue(newValue));
	if (command == fMaximumMomentumRatioCmd)
		fTarget->SetMaximumMomentumRatio(fMaximumMomentumRatioCmd->GetNewDoubleValue(newValue));
	if (command == fStepCmd)
		fTarget->SetChannelingStep(fStepCmd->GetNewDoubleValue(newValue));
	if (command == fRechannelingCmd)
		fTarget->EnableRechanneling(true);
	if (command == fNudEldLimCmd)
		fTarget->SetMinimumDensityLimit(fNudEldLimCmd->GetNewDoubleValue(newValue));
	if (command == fGeantV11AlgorithmCmd) 
		fTarget->UseV11ChannelingAlgo(true);
	if (command == fOrgChanStepSizeCmd)
		fTarget->UseOrgStepSize(true);
	if (command == fManualStepSizeCmd)
		fTarget->UseStepSize(true);
	if (command == fManualStepSizeValue)
		fTarget->SetStepSizeUnit(fManualStepSizeValue->GetNewDoubleValue(newValue));
	if (command == fMeanFPSizeValue)
		fTarget->SetMFPSizeUnit(fMeanFPSizeValue->GetNewDoubleValue(newValue));
	if (command == fPrintDebugInfo)
		fTarget->SetPrintDebugInfo(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String G4ChannelingMessenger::GetCurrentValue(G4UIcommand* command) 
{
	G4String cv;

	if (command==fMinimumEnergyCmd) {
		cv = fMinimumEnergyCmd->ConvertToString(fTarget->GetMinimumEnergy(),"MeV");
	}
	if (command==fMaximumMomentumRatioCmd) {
		cv = fMaximumMomentumRatioCmd->ConvertToString(fTarget->GetMaximumMomentumRatio());
	}
	if (command==fStepCmd ) {
		cv = fStepCmd->ConvertToString(fTarget->GetChannelingStep());
	}
	return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....