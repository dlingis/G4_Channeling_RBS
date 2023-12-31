#include "globals.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// functions for RBS evaluations
G4double RecoilEnergy(G4double E, G4double angle, G4double M1, G4double M2); //evaluates recoiled energy
G4double KinematicFactor(G4double angle, G4double M1, G4double M2);
G4double RandomEnLoss(G4double E, G4double dedx, G4double position, G4double angle); //calculates the energy lost when particle reaches outside volume after the scattering event
G4double CalcDiffRuthXsec(G4double E, G4double M1, G4double M2, G4double angle,G4double Z1, G4double Z2); // calculates the differential Rutherford xsec in laboratory system
//G4double CalcRBSYield(G4double xsec, G4double dist, G4double solidAngle, G4double atomDens);
G4double CalcRBSYield(G4double xsec, G4double dist, G4double solidAngle, G4double atomDens, G4double inc_angle);
//function to calculate Bohr energy straggling
G4double CalcBohrStrag(G4double Z1, G4double Z2, G4double atomDens, G4double dist);
G4double CalcRuthXsecMod(G4double xsec, G4double Z1, G4double Z2, G4double energy);
//function to calculate scattering angle in the CM reference frame
G4double CalcAngleCMFrame(G4double angle, G4double M1, G4double M2);
//function to calculate energy in the CM reference frame
G4double CalcEnergyCMFrame(G4double energy, G4double M1, G4double M2);
// RBS xsec in the CM reference fram
G4double CalcDiffRuthXsecCM(G4double E, G4double angleCM ,G4double Z1, G4double Z2); 
// RBS xsec in the Lab frame from the CM reference fram
G4double CalcDiffRuthXsecLAB(G4double M1, G4double M2, G4double angle, G4double xsection);
// integrates energy loss on the "particle way out" 
G4double CalcTotEnLoss(G4double E, G4double distance, G4int steps, G4ParticleDefinition* fParticle, G4Material* mat);
// function for Total RBS yield, combining other functions into single one
G4double CalculateTotalRBSYield(G4double energy, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double angle, G4double dist,G4double solidAngle, G4double xsecmod, G4double atomDensity, G4double inc_angle);
// Total energy straggling, both electronic and nuclear included
G4double CalculateTotalBohrStraggling(G4double energy, G4ParticleDefinition* particle, G4Material* mat, G4double distance);
// dead layer energy loss
G4double CalculateDeadLayerEffect(G4double energy, const G4Material* dead_mat, G4double thickness,G4ParticleDefinition* particle);
//https://www.sciencedirect.com/science/article/pii/S0168583X97006642
G4double CalcDetectorFWHM(G4double energy, G4double Z1);
G4double CalcAndersenScreening(G4double energy_cm, G4double angle_cm, G4double Z1, G4double Z2);
// from SIMNRA user's guide
G4double CalcNuclEnStraggling(G4double Z1, G4double Z2, G4double M1, G4double M2, G4double atdens,G4double distance);
G4double CalcScreening_TF(G4double Z1, G4double Z2);
G4double CalcScreening_ZBL(G4double Z1, G4double Z2);