#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;

class PrimaryGeneratorMessenger: public G4UImessenger
{
	public:
		PrimaryGeneratorMessenger(PrimaryGeneratorAction *pga);
		~PrimaryGeneratorMessenger();

		void SetNewValue(G4UIcommand * command, G4String newValues);

	private:
		PrimaryGeneratorAction*     fTarget;
		G4UIdirectory*              fPgaDir;

		G4UIcmdWithABool*           useDirAngle;
		G4UIcmdWithADouble*         fDirAngle;
};

#endif