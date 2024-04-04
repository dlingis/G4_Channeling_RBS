#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"


PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* pga)
:fTarget(pga)
{
	fPgaDir = new G4UIdirectory("/pga/", true);
	fPgaDir->SetGuidance("Primary generator action control commands.");

	useDirAngle = new G4UIcmdWithABool("/pga/use_provided_angle", this);
	useDirAngle->SetGuidance("Set if direction vector is to be used.");
	useDirAngle->SetParameterName("use_dir", true);
	useDirAngle->SetDefaultValue(false);

	fDirAngle = new G4UIcmdWithADouble("/pga/direction_angle", this);
	fDirAngle->SetGuidance("Set direction angle.");
	fDirAngle->SetParameterName("angle", true);
	fDirAngle->SetDefaultValue(0.);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
	delete useDirAngle;
	delete fDirAngle;
	delete fPgaDir;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
	if (command == useDirAngle)
		fTarget->SetUseDirAngle(true);
	if (command == fDirAngle)
		fTarget->SetParticleDirectionAngle(fDirAngle->GetNewDoubleValue(newValue));
}