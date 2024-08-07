//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
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
// --------------------------------------------------------------
//

#ifndef DetectorConstructionMessenger_h
#define DetectorConstructionMessenger_h 1

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAIntAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstructionMessenger: public G4UImessenger
{
	public:
		DetectorConstructionMessenger(DetectorConstruction* mpga);
		~DetectorConstructionMessenger();

		virtual void SetNewValue(G4UIcommand * command,G4String newValues);
		virtual G4String GetCurrentValue(G4UIcommand * command);

	private:
		DetectorConstruction *          fTarget;

		G4UIdirectory*                  fMyXtalDirectory;
		G4UIdirectory*                  fMyDetDirectory;

		G4UIcmdWith3VectorAndUnit*      fSizeCmd[5];
		G4UIcmdWith3VectorAndUnit*      fPosCmd[4];
		G4UIcmdWith3VectorAndUnit*      fXtalAngleCmd;
		G4UIcmdWith3VectorAndUnit*      fDetSizesCmd;

		G4UIcmdWithAnInteger*           fGaussCountCmd;
		G4UIcmdWithAnInteger*           fEnLossStepCmd;
		G4UIcmdWithAnInteger*           fNoMatMixCmd;

		G4UIcmdWithAString*             fCustomMatElement1Cmd;
		G4UIcmdWithAString*             fCustomMatElement2Cmd;
		G4UIcmdWithAString*             fCustomMatElement3Cmd;
		G4UIcmdWithAString*             fDeadMaterialCmd;
		G4UIcmdWithAString*             fMixMatNameCmd;
		G4UIcmdWithAString*             fDetectorsMaterialCmd;
		G4UIcmdWithAString*             fWorldMaterial;
		G4UIcmdWithAString*             fDetMaterialCmd;
		G4UIcmdWithAString*             fXtalEC[5][3];
		G4UIcmdWithAString*             fMaterialCmd[5];

		G4UIcmdWithADoubleAndUnit*      fDetResCmd;
		G4UIcmdWithADoubleAndUnit*      fMaxStepCmd;
		G4UIcmdWithADoubleAndUnit*      fRBSAngleCmd;
		G4UIcmdWithADoubleAndUnit*      fRBSROICmd;
		G4UIcmdWithADoubleAndUnit*      fCustomMatDensityCmd;
		G4UIcmdWithADoubleAndUnit*      fDeadThickCmd;
		G4UIcmdWithADoubleAndUnit*      fDetDistanceCmd[5];
		G4UIcmdWithADoubleAndUnit*      fDepthForHistos[12];

		G4UIcmdWithADouble*             fSolidAngleCmd;
		G4UIcmdWithADouble*             fCustomMatElement1PartCmd;
		G4UIcmdWithADouble*             fCustomMatElement2PartCmd;
		G4UIcmdWithADouble*             fCustomMatElement3PartCmd;
		G4UIcmdWithADouble*             fMatMixRatioCmd;

		G4UIcmdWithABool*               fCustomMatCmd;
		G4UIcmdWithABool*               fSigmaCalcCmd;
		G4UIcmdWithABool*               fXtalAmorphous[5];
		G4UIcmdWithABool*               fXtalTrackCmd;
		G4UIcmdWithABool*               fMatMixCmd;
		G4UIcmdWithABool*               fRBSCmd;
		G4UIcmdWithABool*               fFwhmCalcCmd;
		G4UIcmdWithABool*               fConstScatCmd;
		G4UIcmdWithABool*               fUseXStranformation;
		G4UIcmdWithABool*               fHistoTrackingCmd;
};

#endif


