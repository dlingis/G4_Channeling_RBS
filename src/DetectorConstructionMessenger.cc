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

#include "DetectorConstructionMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

DetectorConstructionMessenger::
DetectorConstructionMessenger(DetectorConstruction* mpga)
:fTarget(mpga) {
	fMyXtalDirectory = new G4UIdirectory("/xtal/");
	fMyXtalDirectory->SetGuidance("Crystal setup control commands.");
	
	fWorldMaterial = new G4UIcmdWithAString("/mydet/setWorldMaterial", this);
	fWorldMaterial->SetGuidance("Set world material.");
	fWorldMaterial->SetParameterName("worldMat",true);
	fWorldMaterial->SetDefaultValue("G4_Galactic");


	// material name matrix
	G4String defaultMaterials[5] = {"G4_Si","G4_Si","G4_Si","G4_Si","G4_Si"};
	for (uint8_t i=0; i<5; i++) {
		G4String command = "/xtal/setMaterial" + std::to_string(i + 1);
		fMaterialCmd[i] = new G4UIcmdWithAString(command, this);
		fMaterialCmd[i]->SetGuidance("Set materials.");
		fMaterialCmd[i]->SetParameterName("materials", true);
		fMaterialCmd[i]->SetDefaultValue(defaultMaterials[i]);
	}

	// material size matrix
	G4ThreeVector material_size = G4ThreeVector(1.,1.,1.);
	G4ThreeVector defaultSizes[5] = {material_size,material_size,material_size,material_size,material_size};
	for (uint8_t i=0; i<5; i++) {
		G4String command = "/xtal/setSize" + std::to_string(i + 1);
		fSizeCmd[i] = new G4UIcmdWith3VectorAndUnit(command, this);
		fSizeCmd[i]->SetGuidance("Set crystal size.");
		fSizeCmd[i]->SetParameterName("xtalSizeX", "xtalSizeY", "xtalSizeZ", +true);
		fSizeCmd[i]->SetDefaultValue(defaultSizes[i]);
		fSizeCmd[i]->SetDefaultUnit("mm");
	}
	
	// material position matrix
	G4ThreeVector material_pos = G4ThreeVector(1.,1.,1.);	
	G4ThreeVector defaultPos[4] = {material_pos,material_pos,material_pos,material_pos};
	for (uint8_t i=0; i<4; i++) {
		G4String command = "/xtal/setPos" + std::to_string(i + 1);
		fPosCmd[i] = new G4UIcmdWith3VectorAndUnit(command, this);
		fPosCmd[i]->SetGuidance("Set crystal position.");
		fPosCmd[i]->SetParameterName("crystalPosX", "crystalPosY", "crystalPosZ", true);
		fPosCmd[i]->SetDefaultValue(defaultPos[i]);
		fPosCmd[i]->SetDefaultUnit("mm");
	}

	//Crystal EC's
	for (uint8_t i=0; i<5; i++) {
		for (uint8_t j=0; j<3; j++) {
			G4String command = "/xtal/setEC"+std::to_string(i + 1)+std::to_string(j + 1);
			fXtalEC[i][j] = new G4UIcmdWithAString(command,this);
			fXtalEC[i][j]->SetGuidance("Set xtal ECs");
			fXtalEC[i][j]->SetParameterName("ECs", true);
			fXtalEC[i][j]->SetDefaultValue("data/Si001ax_100K");
		}
	}
	// amorphous crystals 
	for (uint8_t i=0; i<5; i++) {
		G4String command = "/xtal/setAmorphous"+std::to_string(i + 1);
		fXtalAmorphous[i] = new G4UIcmdWithABool(command, this);
		fXtalAmorphous[i]->SetGuidance("Set crystal amorphous");
		fXtalAmorphous[i]->SetParameterName("fCrystalAmorphous",true);
		fXtalAmorphous[i]->SetDefaultValue(false);
	}

	// Particle tracking into histograms
	fXtalTrackCmd = new G4UIcmdWithABool("/xtal/SetParticleTracking", this);
	fXtalTrackCmd->SetGuidance("Enable particle tracking for histograms");
	fXtalTrackCmd->SetParameterName("fTracking", true);
	fXtalTrackCmd->SetDefaultValue(false);
	// set material mixing
	fMatMixCmd = new G4UIcmdWithABool("/xtal/SetMaterialMixing", this);
	fMatMixCmd->SetGuidance("Enable mixing of material");
	fMatMixCmd->SetParameterName("fMatMix", true);
	fMatMixCmd->SetDefaultValue(false);
	// mixing ratio
	fMatMixRatioCmd = new G4UIcmdWithADouble("/xtal/SetMaterialMixingRatio", this);
	fMatMixRatioCmd->SetGuidance("Set material mixing ratio");
	fMatMixRatioCmd->SetParameterName("MIXING_RATIO", true);
	fMatMixRatioCmd->SetDefaultValue(0.);
	// detector energy resolution
	fDetResCmd = new G4UIcmdWithADoubleAndUnit("/xtal/SetDetectorResolution", this);
	fDetResCmd->SetGuidance("Set detector energy resolution.");
	fDetResCmd->SetParameterName("RES", true);
	fDetResCmd->SetDefaultValue(10.);
	fDetResCmd->SetDefaultUnit("keV");
	// Number of steps for Gauss calculation
	fGaussCountCmd = new G4UIcmdWithAnInteger("/xtal/SetGaussSteps", this);
	fGaussCountCmd->SetGuidance("Set # of steps for Gauss Calculation.");
	fGaussCountCmd->SetParameterName("GAUSS", true);
	fGaussCountCmd->SetDefaultValue(5.);
	// crystal rotation angle
	fXtalAngleCmd = new G4UIcmdWith3VectorAndUnit("/xtal/setAngle", this);
	fXtalAngleCmd->SetGuidance("Set crystal angles.");
	fXtalAngleCmd->SetParameterName("xtalAngleX", "xtalAngleY", "xtalAngleZ", true);
	fXtalAngleCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
	fXtalAngleCmd->SetDefaultUnit("rad");
	// max step size
	fMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/xtal/setMaxStep", this);
	fMaxStepCmd->SetGuidance("Set max step size for crystal and amorphous.");
	fMaxStepCmd->SetParameterName("maxStep", true);
	fMaxStepCmd->SetDefaultValue(10.);
	fMaxStepCmd->SetDefaultUnit("nm");
	// en loss steps for rbs calc
	fEnLossStepCmd = new G4UIcmdWithAnInteger("/xtal/SetEnLossStep", this);
	fEnLossStepCmd->SetGuidance("Set # of steps for RBS evaluation.");
	fEnLossStepCmd->SetParameterName("ENLOSS", true);
	fEnLossStepCmd->SetDefaultValue(10.);
	// set RBS angle
	fRBSAngleCmd = new G4UIcmdWithADoubleAndUnit("/xtal/setRBSAngle", this);
	fRBSAngleCmd->SetGuidance("Set RBS angle.");
	fRBSAngleCmd->SetParameterName("ang", true);
	fRBSAngleCmd->SetDefaultValue(160.);
	fRBSAngleCmd->SetDefaultUnit("degree");
	// use of sigmacalc xsec modificators
	fSigmaCalcCmd = new G4UIcmdWithABool("/xtal/SetSigmaCalc", this);
	fSigmaCalcCmd->SetGuidance("Use sigma calc xsec modificators");
	fSigmaCalcCmd->SetParameterName("fSigmaCalc", true);
	fSigmaCalcCmd->SetDefaultValue(false);
	// evaluate RBS spectra
	fRBSCmd = new G4UIcmdWithABool("/xtal/SetRBSEvaluation", this);
	fRBSCmd->SetGuidance("Calculate RBS spectrum");
	fRBSCmd->SetParameterName("fRBSEval", true);
	fRBSCmd->SetDefaultValue(false);
	// set minimum RBS energy ROI
	fRBSROICmd = new G4UIcmdWithADoubleAndUnit("/xtal/setRBSROI", this);
	fRBSROICmd->SetGuidance("Set minimum ROI.");
	fRBSROICmd->SetParameterName("minRoi", true);
	fRBSROICmd->SetDefaultValue(100.);
	fRBSROICmd->SetDefaultUnit("keV");
	// custom material definition
	fCustomMatCmd = new G4UIcmdWithABool("/xtal/EnableCustomMaterial", this);
	fCustomMatCmd->SetGuidance("Enable custom material");
	fCustomMatCmd->SetParameterName("fCustomMat", true);
	fCustomMatCmd->SetDefaultValue(false);
	// custom density
	fCustomMatDensityCmd = new G4UIcmdWithADoubleAndUnit("/xtal/SetCustomMaterialDensity", this);
	fCustomMatDensityCmd->SetGuidance("Set density for custom material.");
	fCustomMatDensityCmd->SetParameterName("DENSITY", true);
	fCustomMatDensityCmd->SetDefaultValue(1.);
	fCustomMatDensityCmd->SetDefaultUnit("g/cm3");
	//element1 part
	fCustomMatElement1PartCmd = new G4UIcmdWithADouble("/xtal/SetCustomMaterialElement1Part", this);
	fCustomMatElement1PartCmd->SetGuidance("Set part of element1 in custom material ");
	fCustomMatElement1PartCmd->SetParameterName("Part_of_el1", true);
	fCustomMatElement1PartCmd->SetDefaultValue(0.);
	//element2 part
	fCustomMatElement2PartCmd = new G4UIcmdWithADouble("/xtal/SetCustomMaterialElement2Part", this);
	fCustomMatElement2PartCmd->SetGuidance("Set part of element2 in custom material ");
	fCustomMatElement2PartCmd->SetParameterName("Part_of_el2", true);
	fCustomMatElement2PartCmd->SetDefaultValue(0.);
	//element3 part
	fCustomMatElement3PartCmd = new G4UIcmdWithADouble("/xtal/SetCustomMaterialElement3Part", this);
	fCustomMatElement3PartCmd->SetGuidance("Set part of element3 in custom material ");
	fCustomMatElement3PartCmd->SetParameterName("Part_of_el3", true);
	fCustomMatElement3PartCmd->SetDefaultValue(0.);
	// element1 name
	fCustomMatElement1Cmd = new G4UIcmdWithAString("/xtal/SetCustomMaterialElement1", this);
	fCustomMatElement1Cmd->SetGuidance("Set element1 of custom material.");
	fCustomMatElement1Cmd->SetParameterName("cMat_el1", true);
	fCustomMatElement1Cmd->SetDefaultValue("G4_Si");
	// element2 name
	fCustomMatElement2Cmd = new G4UIcmdWithAString("/xtal/SetCustomMaterialElement2", this);
	fCustomMatElement2Cmd->SetGuidance("Set element2 of custom material.");
	fCustomMatElement2Cmd->SetParameterName("cMat_el2", true);
	fCustomMatElement2Cmd->SetDefaultValue("G4_Si");
	// element3 name
	fCustomMatElement3Cmd = new G4UIcmdWithAString("/xtal/SetCustomMaterialElement3", this);
	fCustomMatElement3Cmd->SetGuidance("Set element3 of custom material.");
	fCustomMatElement3Cmd->SetParameterName("cMat_el3", true);
	fCustomMatElement3Cmd->SetDefaultValue("G4_Si");
	// dead material definition
	fDeadMaterialCmd = new G4UIcmdWithAString("/xtal/SetDeadMaterialName", this);
	fDeadMaterialCmd->SetGuidance("Set dead material.");
	fDeadMaterialCmd->SetParameterName("dead_mat", true);
	fDeadMaterialCmd->SetDefaultValue("G4_Si");
	// thickness of dead layer
	fDeadThickCmd = new G4UIcmdWithADoubleAndUnit("/xtal/SetDeadLayerThick", this);
	fDeadThickCmd->SetGuidance("Set thickness of dead layer.");
	fDeadThickCmd->SetParameterName("dead_thick", true);
	fDeadThickCmd->SetDefaultValue(1.);
	fDeadThickCmd->SetDefaultUnit("nm");
	// detector solid angle
	fSolidAngleCmd = new G4UIcmdWithADouble("/xtal/SetSolidAngle", this);
	fSolidAngleCmd->SetGuidance("Set detector solid angle");
	fSolidAngleCmd->SetParameterName("Solid_angle", true);
	fSolidAngleCmd->SetDefaultValue(1.);
	// number of material for mixing
	fNoMatMixCmd = new G4UIcmdWithAnInteger("/xtal/SetNoOfMaterialMix", this);
	fNoMatMixCmd->SetGuidance("Set no of material for mix.");
	fNoMatMixCmd->SetParameterName("mix", true);
	fNoMatMixCmd->SetDefaultValue(1);
	// set mixing material name
	fMixMatNameCmd = new G4UIcmdWithAString("/xtal/SetMixingMaterialName", this);
	fMixMatNameCmd->SetGuidance("Set material name for mixing.");
	fMixMatNameCmd->SetParameterName("mat_mix", true);
	fMixMatNameCmd->SetDefaultValue("G4_Si");
	// set detector fwhm calculation
	fFwhmCalcCmd = new G4UIcmdWithABool("/xtal/SetFWHMCalc", this);
	fFwhmCalcCmd->SetGuidance("Set for calculation of fwhm");
	fFwhmCalcCmd->SetParameterName("fFwhmCalc", true);
	fFwhmCalcCmd->SetDefaultValue(false);
	// set to use constant scattering angle (set by RBS angle command)	
	fConstScatCmd = new G4UIcmdWithABool("/xtal/UseConstScatAngle", this);
	fConstScatCmd->SetGuidance("Use constant scattering angle");
	fConstScatCmd->SetParameterName("fConstScatAng", true);
	fConstScatCmd->SetDefaultValue(false);
	// use xsec transformation
	fUseXStranformation = new G4UIcmdWithABool("/xtal/use_xs_transformation", this);
	fUseXStranformation->SetGuidance("Use XS tranformation for channeling");
	fUseXStranformation->SetParameterName("fUseXST", true);
	fUseXStranformation->SetDefaultValue(false);
	//********************************************
	fMyDetDirectory = new G4UIdirectory("/mydet/");
	fMyDetDirectory->SetGuidance("Detector setup control commands.");
	fDetectorsMaterialCmd = new G4UIcmdWithAString("/mydet/setDetMaterial", this);
	fDetectorsMaterialCmd->SetGuidance("Set detector material.");
	fDetectorsMaterialCmd->SetParameterName("xDetMat", true);
	fDetectorsMaterialCmd->SetDefaultValue("");
	// set mat size
	fDetSizesCmd = new G4UIcmdWith3VectorAndUnit("/mydet/setSize", this);
	fDetSizesCmd->SetGuidance("Set detector size.");
	fDetSizesCmd->SetParameterName("detSizeX", "detSizeY", "detSizeZ", true);
	fDetSizesCmd->SetDefaultValue(G4ThreeVector(50.,50.,1.));
	fDetSizesCmd->SetDefaultUnit("mm");

	G4double defaultDistances[4] = {-20.,-19.,+19,+20.};
	for (uint8_t i=0; i<4; i++) {
		G4String command = "/mydet/setDistance" + std::to_string(i + 1);
		fDetDistanceCmd[i] = new G4UIcmdWithADoubleAndUnit(command, this);
		fDetDistanceCmd[i]->SetGuidance("Set detector size.");
		fDetDistanceCmd[i]->SetParameterName("detSize", true);
		fDetDistanceCmd[i]->SetDefaultValue(defaultDistances[i]);
		fDetDistanceCmd[i]->SetDefaultUnit("mm");
	}

	G4double defaultDepths[12] = {0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,300.};
	for (uint8_t i=0; i<12; i++) {
		G4String command = "/mydet/setHistoDepthIntervals" + std::to_string(i + 1);
		fDepthForHistos[i] = new G4UIcmdWithADoubleAndUnit(command, this);
		fDepthForHistos[i]->SetGuidance("Set depth intervals for channeling histos");
		fDepthForHistos[i]->SetParameterName("depthInt", true);
		fDepthForHistos[i]->SetDefaultValue(defaultDepths[i]);
		fDepthForHistos[i]->SetDefaultUnit("nm");
	}
	fHistoTrackingCmd = new G4UIcmdWithABool("/mydet/trackChannelingHistos", this);
	fHistoTrackingCmd->SetGuidance(" Enable tracking of particles for histograms ");
	fHistoTrackingCmd->SetParameterName("fTrackHisto", true);
	fHistoTrackingCmd->SetDefaultValue(false);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstructionMessenger::
~DetectorConstructionMessenger() {
	delete fWorldMaterial;
	
	delete fMatMixCmd;
	delete fMatMixRatioCmd;
	delete fDetResCmd;
	delete fGaussCountCmd;

	delete fXtalAngleCmd;

	delete fMaxStepCmd;
	delete fRBSAngleCmd;
	delete fSigmaCalcCmd;

	delete fRBSCmd;
	delete fRBSROICmd;
	delete fEnLossStepCmd;

	// size matrix
	for (uint8_t i=0; i<5; i++) {
		delete fSizeCmd[i];
	}
	// material name matrix
	for (uint8_t i=0; i<5; i++) {
		delete fMaterialCmd[i];
	}
	// material position matrix
	for (uint8_t i=0; i<4; i++) {
		delete fPosCmd[i];
	}
	// amo and crystal characs
	for (uint8_t i=0; i<5; i++) {
		delete fXtalAmorphous[i];
		for (uint8_t j=0; j<3; j++) {
			delete fXtalEC[i][j];
		}
	}

	delete fXtalTrackCmd;
	delete fCustomMatCmd;
	delete fCustomMatDensityCmd;
	delete fCustomMatElement1PartCmd;
	delete fCustomMatElement2PartCmd;
	delete fCustomMatElement3PartCmd;
	delete fCustomMatElement1Cmd;
	delete fCustomMatElement2Cmd;
	delete fCustomMatElement3Cmd;
	
	delete fDeadMaterialCmd;
	delete fDeadThickCmd;
	delete fSolidAngleCmd;
	
	delete fNoMatMixCmd;
	delete fMixMatNameCmd;
	delete fFwhmCalcCmd;
	delete fConstScatCmd;
	delete fUseXStranformation;

	delete fDetectorsMaterialCmd;
	delete fDetSizesCmd;
	for (uint8_t i=0; i<5; i++) {
		delete fDetDistanceCmd[i];
	}

	for (uint8_t i=0; i<12; i++) {
		delete fDepthForHistos[i];
	}

	delete fMyXtalDirectory;
	delete fMyDetDirectory;
	delete fHistoTrackingCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstructionMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
	if (command == fWorldMaterial)
		fTarget->SetWorldMaterial(newValue);
	if (command == fXtalTrackCmd)
		fTarget->SetTracking(true);
	if (command == fMatMixCmd)
		fTarget->SetMixing(true);
	if (command == fMatMixRatioCmd)
		fTarget->SetMixingRatio(fMatMixRatioCmd->GetNewDoubleValue(newValue));
	if (command == fDetResCmd)
		fTarget->SetDetectorResolution(fDetResCmd->GetNewDoubleValue(newValue));
	if (command == fGaussCountCmd)
		fTarget->SetGaussCounter(fGaussCountCmd->GetNewIntValue(newValue));
	if (command == fXtalAngleCmd )
		fTarget->SetAngles(fXtalAngleCmd->GetNew3VectorValue(newValue));
	if (command == fEnLossStepCmd)
		fTarget->SetEnLossStep(fEnLossStepCmd->GetNewIntValue(newValue));
	if (command == fSigmaCalcCmd)
		fTarget->SetSigmaCal(true);
	if (command == fRBSCmd)
		fTarget->SetRBSCalc(true);
	if (command == fRBSROICmd)
		fTarget->SetRBSROImin(fRBSROICmd->GetNewDoubleValue(newValue));
	if (command == fMaxStepCmd)
		fTarget->SetMaxStep(fMaxStepCmd->GetNewDoubleValue(newValue));
	if (command == fRBSAngleCmd)
		fTarget->SetRBSAngle(fRBSAngleCmd->GetNewDoubleValue(newValue));
	if (command == fCustomMatCmd)
		fTarget->SetCustomMaterial(true);
	if (command == fCustomMatDensityCmd)
		fTarget->SetCustomMaterialDensity(fCustomMatDensityCmd->GetNewDoubleValue(newValue));
	if (command == fCustomMatElement1PartCmd)
		fTarget->SetCustomElement1Part(fCustomMatElement1PartCmd->GetNewDoubleValue(newValue));
	if (command == fCustomMatElement2PartCmd)
		fTarget->SetCustomElement2Part(fCustomMatElement2PartCmd->GetNewDoubleValue(newValue));
	if (command == fCustomMatElement3PartCmd)
		fTarget->SetCustomElement3Part(fCustomMatElement3PartCmd->GetNewDoubleValue(newValue));
	if (command == fCustomMatElement1Cmd)
		fTarget->SetCustomElement1(newValue);
	if (command == fCustomMatElement2Cmd)
		fTarget->SetCustomElement2(newValue);
	if (command == fCustomMatElement3Cmd)
		fTarget->SetCustomElement3(newValue);
	if (command == fDeadThickCmd)
		fTarget->SetDeadLayerThickness(fDeadThickCmd->GetNewDoubleValue(newValue));
	if (command == fDeadMaterialCmd)
		fTarget->SetDeadLayer(newValue);
	if (command == fSolidAngleCmd)
		fTarget->SetSolidAngle(fSolidAngleCmd->GetNewDoubleValue(newValue));
	if (command == fNoMatMixCmd)
		fTarget->SetMaterialForMix(fNoMatMixCmd->GetNewIntValue(newValue));
	if (command == fMixMatNameCmd)
		fTarget->SetMixingMaterial(newValue);
	if (command == fFwhmCalcCmd)
		fTarget->SetCalcFWHM(newValue);
	if (command == fConstScatCmd)
		fTarget->UseConstAngle(newValue);
	if (command == fDetectorsMaterialCmd )
		fTarget->SetDetectorMaterial(newValue);
	if (command == fDetSizesCmd )
		fTarget->SetDetectorSizes(fDetSizesCmd->GetNew3VectorValue(newValue));
	if (command == fUseXStranformation)
		fTarget->UseXStranformation(newValue);
	if (command == fHistoTrackingCmd)
		fTarget->SetTrackingHistos(true);

	// material setting
	for (uint8_t i=0; i<5; i++) {
		if (command == fMaterialCmd[i]) {
			fTarget->SetMaterial(i, newValue);
		}
	}
	// material size setting
	for (uint8_t i=0; i<5; i++) {
		if (command == fSizeCmd[i]) {
			fTarget->SetSize(i,fSizeCmd[i]->GetNew3VectorValue(newValue));
		}
	}
	// material position setting
	for (uint8_t i=0; i<4; i++) {
		if (command == fPosCmd[i]) {
			fTarget->SetPosition(i, fPosCmd[i]->GetNew3VectorValue(newValue));
		}
	}
	// amorphous and xtal ec setting
	for (uint8_t i=0;i<5;i++) {
		if (command == fXtalAmorphous[i]) {
			fTarget->SetAmorphous(i, true);
		} else {
			for (uint8_t j=0; j<3; j++) {
				if (command == fXtalEC[i][j])
					fTarget->SetEC(i, j, newValue);
			}
		}
	}
	// detector distance
	for (uint8_t i=0; i<5; i++) {
		if (command == fDetDistanceCmd[i]) {
			fTarget->SetDetectorDistance(i, fDetDistanceCmd[i]->GetNewDoubleValue(newValue));
		}
	}
	// depth limits for histograms
	for (uint8_t i=0; i<12; i++) {
		if (command == fDepthForHistos[i]) {
			fTarget->SetArrayForHistos(i, fDepthForHistos[i]->GetNewDoubleValue(newValue));
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String DetectorConstructionMessenger::GetCurrentValue(G4UIcommand * command)
{
	G4String cv;

	if (command == fWorldMaterial) {
		cv = fTarget->GetWorldMaterial();
	}
	if (command == fDetResCmd) {
		cv = fTarget->GetDetectorResolution();
	}
	if (command == fXtalAngleCmd) {
		cv = fXtalAngleCmd->ConvertToString(fTarget->GetAngles(), "rad");
	}
	if (command ==  fEnLossStepCmd) {
		cv = fTarget->GetEnLossStep();
	}
	if (command == fMaxStepCmd) {
		cv = fTarget->GetMaxStep();
	}
	if (command == fRBSAngleCmd) {
		cv = fTarget->GetRBSAngle();
	}
	// material set
	for (uint8_t i=0; i<5; i++) {
		if (command == fMaterialCmd[i]) {
			cv = fTarget->GetMaterial(i);
		}
	}
	// size set
	for (uint8_t i=0; i<5; i++) {
		if (command == fSizeCmd[i]) {
			cv = fSizeCmd[i]->ConvertToString(fTarget->GetSize(i), "mm");
		}
	}
	// layer positions
	for (uint8_t i=0; i<4; i++) {
		if (command == fPosCmd[i]) {
			cv = fPosCmd[i]->ConvertToString(fTarget->GetPosition(i), "mm");
		}
	}
	// detector distance
	for (uint8_t i=0; i<5; i++) {
		if (command == fDetDistanceCmd[i]) {
			cv = fDetDistanceCmd[i]->ConvertToString(fTarget->GetDetectorDistance(i), "mm");
		}
	}
	return cv;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....