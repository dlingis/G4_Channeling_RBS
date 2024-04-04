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

#ifndef G4ChannelingMessenger_h
#define G4ChannelingMessenger_h 1

class G4Channeling;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ChannelingMessenger: public G4UImessenger
{
	public:
		G4ChannelingMessenger(G4Channeling* mpgax);
		~G4ChannelingMessenger();
		
		virtual void SetNewValue(G4UIcommand* command,G4String newValues);
		virtual G4String GetCurrentValue(G4UIcommand* command);
		
	private:
		G4UIdirectory* fChannelingDirectory;
		G4Channeling * fTarget;

		G4UIcmdWithADouble*  fMaximumMomentumRatioCmd;
		G4UIcmdWithADouble* fStepCmd;
		G4UIcmdWithADouble* fNudEldLimCmd;

		G4UIcmdWithABool* fRechannelingCmd;
		G4UIcmdWithABool* fGeantV11AlgorithmCmd;
		G4UIcmdWithABool* fOrgChanStepSizeCmd;
		G4UIcmdWithABool* fManualStepSizeCmd;
		G4UIcmdWithABool* fPrintDebugInfo;

		G4UIcmdWithADoubleAndUnit* fMinimumEnergyCmd;
		G4UIcmdWithADoubleAndUnit* fManualStepSizeValue;
		G4UIcmdWithADoubleAndUnit* fMeanFPSizeValue;
};

#endif
