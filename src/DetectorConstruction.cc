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
//

#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4CrystalExtension.hh"
#include "G4ExtendedMaterial.hh"
#include "G4LogicalCrystalVolume.hh"

#include "G4ChannelingMaterialData.hh"
#include "G4ChannelingOptrMultiParticleChangeCrossSection.hh"

#include "CrystalDetector.hh"
#include "SensitiveDetector.hh"

#include "G4SDManager.hh"

#define EC_SAMPLE     "data/Si001ax_100K"
#define SZ_SAMPLE     G4ThreeVector(1.,1.,1.)
#define POS_SAMPLE    G4ThreeVector(0.,0.,0.)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():
fECfileName{EC_SAMPLE, EC_SAMPLE, EC_SAMPLE,
	EC_SAMPLE, EC_SAMPLE, EC_SAMPLE,
	EC_SAMPLE, EC_SAMPLE, EC_SAMPLE,
	EC_SAMPLE, EC_SAMPLE, EC_SAMPLE,
	EC_SAMPLE, EC_SAMPLE, EC_SAMPLE},
stepLimit(0),
fMaterialName{"G4_Si","G4_Si","G4_Si","G4_Si","G4_Si"},
fSizes{SZ_SAMPLE, SZ_SAMPLE, SZ_SAMPLE, SZ_SAMPLE, SZ_SAMPLE},
fAngles(POS_SAMPLE),
fWorldMaterial("G4_Galactic"),
fDetectorMaterialName(""),
fDetectorSizes(G4ThreeVector(50. * CLHEP::mm,50. * CLHEP::mm,1 * CLHEP::mm)),
fDetectorDistance{-20. * CLHEP::cm,-19. * CLHEP::cm,+19. * CLHEP::cm,+20. * CLHEP::cm,+40. * CLHEP::cm}, 
fCrystalAmorphous{0,0,0,0,0},
position{POS_SAMPLE, POS_SAMPLE, POS_SAMPLE, POS_SAMPLE},
rbs_angle(160.*degree),sigma_calc(0),rbs_calc(0),rbs_step(0),/*add element or material for mix*/material_mixing(0.), sec_material_ratio(0.),
detector_resolution(10.*keV),gauss_counter(5),
//custom material definition
enable_custom_material(0),element1(""),element2(""),element3(""),custom_density(0.),part1(0),part2(0),part3(0),
// other parameters
dead_thickness(100.),solidAngle(1.),
material_for_mix(0.),mix_material_name("G4_Si"),enable_fwhm_calc(0.),rbs_roi_min(50.*keV),use_const_angle(0.),
use_xs_transformation(0),
histo_tracking(false)
{
	fMessenger = new DetectorConstructionMessenger(this); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::DefineMaterials(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){

	//** World **//
	G4Material* worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial(fWorldMaterial);

	G4double worldSizeXY = 2. * CLHEP::meter;
	G4double worldSizeZ = 30. * CLHEP::meter; // buvo 30 metru
	G4Box* worldSolid = new G4Box("world.solid",
	                              worldSizeXY/2.,
	                              worldSizeXY/2.,
	                              worldSizeZ/2.);

	G4LogicalVolume* worldLogic = new G4LogicalVolume(worldSolid,
	                                                  worldMaterial,
	                                                  "world.logic");

	worldPhysical = new G4PVPlacement(0,
                                    G4ThreeVector(),
                                    worldLogic,
                                    "world.physic",
                                    0,
                                    false,
                                    0);

	if (position[0].z()>position[1].z() || position[1].z() > position[2].z() || position[2].z() > position[3].z()) {
		G4cout << " ************************ " << G4endl;
		G4cout << " ERROR in layer placement " << G4endl;
		G4cout << " ************************ " << G4endl;
		exit(1);
	}

	for (G4int i=1; i<12; i++) {
		if (array_for_histos[i-1] > array_for_histos[i]) {
			G4cout << " ************************************** " << G4endl;
			G4cout << " Channeling profile intervals incorrect " << G4endl;
			G4cout << " ************************************** " << G4endl;
			exit(1);
		}
	}

	//** Detectors instantiation **//
	G4Box* ssdSolid = new G4Box("ssd.solid",
	                            fDetectorSizes.x()/2.,
	                            fDetectorSizes.y()/2.,
	                            fDetectorSizes.z()/2.);

	G4Material* detectorMaterial;
	if (fDetectorMaterialName == "")
		detectorMaterial = worldMaterial;
	else
		detectorMaterial = G4NistManager::Instance()->FindOrBuildMaterial(fDetectorMaterialName);

	G4LogicalVolume* ssdLogic = new G4LogicalVolume(ssdSolid, detectorMaterial, "ssd.logic");

	for (size_t i1=0; i1<4; i1++) {
		new G4PVPlacement(0,
	                      G4ThreeVector(0.,0.,fDetectorDistance[i1]),
	                      ssdLogic,
	                      "ssd.physic",
	                      worldLogic,
	                      false,
	                      i1);
	}

	G4Box* crystalSolid = new G4Box("crystal.Solid",
	                                fSizes[0].x()/2.,
	                                fSizes[0].y()/2.,
	                                fSizes[0].z()/2.);

	G4Box* crystalSolid2 = new G4Box("crystal.solid2",
	                                fSizes[1].x()/2.,
	                                fSizes[1].y()/2.,
	                                fSizes[1].z()/2.);

	G4Box* crystalSolid3 = new G4Box("crystal.solid3",
	                                fSizes[2].x()/2.,
	                                fSizes[2].y()/2.,
	                                fSizes[2].z()/2.);

	G4Box* crystalSolid4 = new G4Box("crystal.solid4",
	                                fSizes[3].x()/2.,
	                                fSizes[3].y()/2.,
	                                fSizes[3].z()/2.);

	G4Box* crystalSolid5 = new G4Box("crystal.solid5",
	                                fSizes[4].x()/2.,
	                                fSizes[4].y()/2.,
	                                fSizes[4].z()/2.);

	G4String name, symbol;

	// Elements    
	G4Element* elHfn = G4NistManager::Instance()->FindOrBuildElement(72);
	G4Element* elOx = G4NistManager::Instance()->FindOrBuildElement(8);
	G4Element* elNb = G4NistManager::Instance()->FindOrBuildElement(41);
	G4Element* elIn = G4NistManager::Instance()->FindOrBuildElement(49);
	G4Element* elP = G4NistManager::Instance()->FindOrBuildElement(15);
	G4Element* elSi = G4NistManager::Instance()->FindOrBuildElement(14);
	G4Element* elC = G4NistManager::Instance()->FindOrBuildElement(6);
	//G4Element* elAr = G4NistManager::Instance()->FindOrBuildElement(18);

	G4int ncomponents, natoms;

	// Nb2O5
	G4Material* Nb2O5 = new G4Material("Nb2O5", 4.6 *g/cm3, 2);
	Nb2O5->AddElement(elNb, 2);
	Nb2O5->AddElement(elOx, 5);

	G4Material* InP = new G4Material("InP", 4.81 *g/cm3, 2);
	InP->AddElement(elIn, 1);
	InP->AddElement(elP, 1);

	G4Material* HfO2 = new G4Material("HfO2",9.68 *g/cm3, 2);
	HfO2->AddElement(elHfn, 1);
	HfO2->AddElement(elOx, 2); 

	//https://srdata.nist.gov/CeramicDataPortal/Pds/Scdscs
	G4Material* SiC = new G4Material("SiC", 3.16 *g/cm3, 2);
	SiC->AddElement(elSi, 1);
	SiC->AddElement(elC, 1);

	G4Material* diamond = new G4Material("Diamond", 3.5 *g/cm3, 1);
	diamond->AddElement(elC,1);

	G4Material* graphite = new G4Material("Graphite", 2.25 *g/cm3, 1);
	graphite->AddElement(elC,1);


	G4String plus = " + ";
	G4String gap  = " ";


	for (G4int i=0; i<5; ++i) {
		G4cout << " Material " << i << " name: " << fMaterialName[i] << G4endl;
		if (fMaterialName[i] == "")
			mat[i] = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
		else if (fMaterialName[i] == "Nb2O5")
			mat[i] = Nb2O5;
		else if (fMaterialName[i] == "InP")
			mat[i] = InP;
		else if (fMaterialName[i] == "HfO2")
			mat[i] = HfO2;
		else if (fMaterialName[i] == "SiC")
			mat[i] = SiC;
		else if (fMaterialName[i] == "Diamond")
			mat[i] = diamond;
		else if (fMaterialName[i] == "Graphite")
			mat[i] = graphite;
		else
			mat[i] = G4NistManager::Instance()->FindOrBuildMaterial(fMaterialName[i]);
	}


	// detector dead layer
	if (dead_material_name == "")
		dead_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
	else
	dead_material = G4NistManager::Instance()->FindOrBuildMaterial(dead_material_name);

	G4cout << "Detector dead layer: " << dead_thickness/nm << " nm of " << dead_material_name << " Z " << dead_material->GetZ() << G4endl;

	G4double a = 1.*um, b = 20.*cm;;
	G4Box* dead_layer_solid = new G4Box("dead_layer_solid",a,a,dead_thickness);
	G4LogicalVolume* dead_layer_log = new G4LogicalVolume (dead_layer_solid, dead_material, "dead_layer_logic");
	G4VPhysicalVolume* dead_layer_phys = new G4PVPlacement(0, G4ThreeVector(b, b, b), dead_layer_log, "dlayer_physic", worldLogic, false, 0,true);
	// end

	if (material_mixing != 0) {
		G4cout << " Mixing of materials enabled " << G4endl;
		G4int no_of_material_mix = GetMaterialForMix();	
		G4cout << " Material " << no_of_material_mix << " name " << fMaterialName[no_of_material_mix] << G4endl;
		G4cout << " To be mixed with " << mix_material_name << G4endl;
		mixing_material = G4NistManager::Instance()->FindOrBuildMaterial(mix_material_name);
		G4double den1 = mixing_material->GetDensity() /(g/cm3);
		G4double den2 = mat[no_of_material_mix]->GetDensity() /(g/cm3);
		G4String ratio = std::to_string(sec_material_ratio);
		G4String ratio2 = std::to_string(1-sec_material_ratio);
		G4String mix_name = ratio + gap + fMaterialName[no_of_material_mix] + plus + ratio2 + gap + mix_material_name;
		G4Material* material_mix = new G4Material(mix_name, (den1 * (1 - sec_material_ratio) + den2 * sec_material_ratio) *g/cm3, 2);

		material_mix->AddMaterial(mat[no_of_material_mix], sec_material_ratio);
		material_mix->AddMaterial(mixing_material, 1 - sec_material_ratio);

		mat[no_of_material_mix] = material_mix;
		G4cout << " Mixed material density = " << material_mix->GetDensity() /(g/cm3) << G4endl;
	}

	if (enable_custom_material) {
		G4cout << " Enabled custom material " << G4endl;
		//part2 = 1-part1;
		G4int no_of_material_mix = GetMaterialForMix();
		G4String ratio1 = std::to_string(part1);
		G4String ratio2 = std::to_string(part2);
		G4String ratio3 = std::to_string(part3);
		G4String custom_material_name = ratio1 + element1 + plus + ratio2 + element2 + plus + ratio3 + element3;
		// addition of custom material
		G4Element* el1 = G4NistManager::Instance()->FindOrBuildElement(element1);
		G4cout << " el 1: " << el1->GetName() << " Z " << el1->GetZ() << " S " << el1->GetSymbol() << " Mass % " << part1 << G4endl;
		G4Element* el2 = G4NistManager::Instance()->FindOrBuildElement(element2);
		G4cout << " el 2:  " << el2->GetName() << " Z " << el2->GetZ() << " S " << el2->GetSymbol() << " Mass % " << part2 <<  G4endl;
		G4Element* el3 = G4NistManager::Instance()->FindOrBuildElement(element3);
		G4cout << " el 3:  " << el3->GetName() << " Z " << el3->GetZ() << " S " << el3->GetSymbol() << " Mass % " << part3 <<  G4endl;
		G4Material* custom_material = new G4Material(custom_material_name,custom_density,3);
		custom_material->AddElement(el1,part1);
		custom_material->AddElement(el2,part2);
		custom_material->AddElement(el3,part3);
		mat[no_of_material_mix] = custom_material;
	}

	// separation of materials into components
	G4String mat_name = "material";
	G4String el_name = "element";
	G4double avogadro = 6.022e+23;

	for (G4int i=0; i<5; i++) {
		G4int no_of_elements = mat[i]->GetNumberOfElements();
		for(int j=0; j<no_of_elements; j++) {
			G4String total_name = mat_name + std::to_string(i) + el_name + std::to_string(j);
			G4double Z = mat[i]->GetElement(j)->GetZ();
			G4double M = mat[i]->GetElement(j)->GetA() /(g/mole);
			const G4double *atomDensVector = mat[i]->GetVecNbOfAtomsPerVolume();
			G4double aDensity = atomDensVector[j] /(1/cm3);
			G4double mDensity = aDensity * M / avogadro;
			mat_components[i][j] = new G4Material(total_name, Z, M, mDensity *(g/cm3));
		}
	}
	/*
	// creation of phantom volumes for G4EmCalculator error
	// regarding couples
	G4double test_distance = -2*m;
	G4String name2 = "phantom2";
	G4String type = "logic";
	G4String type2 = "physical";
	
  G4Box* phantom2 = new G4Box("phantom2",
                                1*nm,
                                1*nm,
                                1*nm);
    
  // -20 m
  for(int i = 0; i<5; i++) {
    G4int no_of_elements = mat[i]->GetNumberOfElements();
    	for(int k = 0; k< no_of_elements; k++) {
    	G4LogicalVolume* phantom2_logic = new G4LogicalVolume(phantom2,
                                                    mat_components[i][k],
                                                    name2+type+std::to_string(i+k));
    
        G4VPhysicalVolume* phantom2_physic = new G4PVPlacement(0,
                          G4ThreeVector(0.,0.,test_distance+(i*nm)),
                          phantom2_logic,
                          name2+type2+std::to_string(i+k),
                          worldLogic,
                          false,
                          0);
        }
    }
*/
	stepLimit = new G4UserLimits(maxStep);

	G4RotationMatrix* rot = new G4RotationMatrix;
	if(fAngles.x()!=0.)
		rot->rotateX(fAngles.x());
	if(fAngles.y()!=0.)
		rot->rotateY(fAngles.y());
	if(fAngles.z()!=0.)
		rot->rotateZ(fAngles.z());


	// main layer crystal
	Crystal = new G4ExtendedMaterial("crystal.material",mat[0]);
	Crystal->RegisterExtension(std::unique_ptr<G4CrystalExtension>(new G4CrystalExtension(Crystal)));
	G4CrystalExtension* crystalExtension = (G4CrystalExtension*)Crystal->RetrieveExtension("crystal");


	crystalExtension->SetUnitCell(new G4CrystalUnitCell(5.43 * CLHEP::angstrom,
	                                                    5.43 * CLHEP::angstrom,
	                                                    5.43 * CLHEP::angstrom,
	                                                    CLHEP::halfpi,
	                                                    CLHEP::halfpi,
	                                                    CLHEP::halfpi,
	                                                    227));

	Crystal->RegisterExtension(std::unique_ptr<G4ChannelingMaterialData>(new G4ChannelingMaterialData("channeling")));
	G4ChannelingMaterialData* crystalChannelingData = (G4ChannelingMaterialData*)Crystal->RetrieveExtension("channeling");



	crystalLogic_amo = new G4LogicalVolume (crystalSolid, mat[0], "crystal.logic_amo");
	crystalLogic_cry = new G4LogicalCrystalVolume(crystalSolid, Crystal,"crystal.logic_cry"); 

	if (fCrystalAmorphous[0]) {
		//crystalLogic_amo = new G4LogicalVolume (crystalSolid, mat[0], "crystal.logic_amo"); // new
		physAbsor = new G4PVPlacement(rot,
		                              G4ThreeVector(),
		                              crystalLogic_amo,
		                              "crystal.physic_amo",
		                              worldLogic,
		                              false,
		                              0,
		                              true);

		crystalLogic_amo->SetUserLimits(stepLimit);
	} else {
		crystalChannelingData->SetFilename(fECfileName[0][0]);
		crystalLogic_cry->SetVerbose(1);

		//crystalLogic_cry = new G4LogicalCrystalVolume(crystalSolid, Crystal,"crystal.logic_cry");  // new
		physAbsor = new G4PVPlacement(rot,
		                              G4ThreeVector(),
		                              crystalLogic_cry,
		                              "crystal.physic_cry",
		                              worldLogic,
		                              false,
		                              0,
		                              true);
	}
	// layer1 crystal
	Crystal2 = new G4ExtendedMaterial("crystal2.material",mat[1]);
	Crystal2->RegisterExtension(std::unique_ptr<G4CrystalExtension>(new G4CrystalExtension(Crystal2)));
	G4CrystalExtension* crystalExtension2 = (G4CrystalExtension*)Crystal2->RetrieveExtension("crystal");


	crystalExtension2->SetUnitCell(new G4CrystalUnitCell(5.43 * CLHEP::angstrom,
	                                                    5.43 * CLHEP::angstrom,
	                                                    5.43 * CLHEP::angstrom,
	                                                    CLHEP::halfpi,
	                                                    CLHEP::halfpi,
	                                                    CLHEP::halfpi,
	                                                    227));

	Crystal2->RegisterExtension(std::unique_ptr<G4ChannelingMaterialData>(new G4ChannelingMaterialData("channeling")));
	G4ChannelingMaterialData* crystalChannelingData2 = (G4ChannelingMaterialData*)Crystal2->RetrieveExtension("channeling");

	if (fCrystalAmorphous[1]) {
		crystalLogic2_amo = new G4LogicalVolume (crystalSolid2, mat[1], "crystal2.logic_amo");
		if(fCrystalAmorphous[0]) {
			physAbsor2 = new G4PVPlacement(0,
			                              position[0],
			                              crystalLogic2_amo,
			                              "crystal2.physic_amo",
			                              crystalLogic_amo,
			                              false,
			                              0,
			                              true);
			crystalLogic2_amo->SetUserLimits(stepLimit);
		} else {
			physAbsor2 = new G4PVPlacement(0,
			                              position[0],
			                              crystalLogic2_amo,
			                              "crystal2.physic_amo",
			                              crystalLogic_cry,
			                              false,
			                              0,
			                              true);
		}
	} else {
		crystalLogic2_cry = new G4LogicalCrystalVolume(crystalSolid2, Crystal2,"crystal2.logic_cry"); 
		crystalChannelingData2->SetFilename(fECfileName[1][0]);
		crystalLogic2_cry->SetVerbose(1);

		if (fCrystalAmorphous[0]) {
			physAbsor2 = new G4PVPlacement(0,
			                               position[0],
			                               crystalLogic2_cry,
			                               "crystal2.physic_cry",
			                               crystalLogic_amo,
			                               false,
			                               0,
			                               true);
		} else {
			physAbsor2 = new G4PVPlacement(0,
			                               position[0],
			                               crystalLogic2_cry,
			                               "crystal2.physic_cry",
			                               crystalLogic_cry,
			                               false,
			                               0,
			                               true);
		}
	}

	// layer2 crystal
	Crystal3 = new G4ExtendedMaterial("crystal3.material",mat[2]);
	Crystal3->RegisterExtension(std::unique_ptr<G4CrystalExtension>(new G4CrystalExtension(Crystal3)));
	G4CrystalExtension* crystalExtension3 = (G4CrystalExtension*)Crystal3->RetrieveExtension("crystal");


	crystalExtension3->SetUnitCell(new G4CrystalUnitCell(5.43 * CLHEP::angstrom,
	                                                    5.43 * CLHEP::angstrom,
	                                                    5.43 * CLHEP::angstrom,
	                                                    CLHEP::halfpi,
	                                                    CLHEP::halfpi,
	                                                    CLHEP::halfpi,
	                                                    227));

	Crystal3->RegisterExtension(std::unique_ptr<G4ChannelingMaterialData>(new G4ChannelingMaterialData("channeling")));
	G4ChannelingMaterialData* crystalChannelingData3 = (G4ChannelingMaterialData*)Crystal3->RetrieveExtension("channeling");


	crystalLogic3_amo = new G4LogicalVolume (crystalSolid3, mat[2], "crystal3.logic_amo");
	crystalLogic3_cry = new G4LogicalCrystalVolume(crystalSolid3, Crystal3,"crystal3.logic_cry"); 

	if (fCrystalAmorphous[2]) {
		if (fCrystalAmorphous[0]) {
			physAbsor3 = new G4PVPlacement(0,
			                               position[1],
			                               crystalLogic3_amo,
			                               "crystal3.physic_amo",
			                               crystalLogic_amo,
			                               false,
			                               0,
			                               true);
		} else {
			physAbsor3 = new G4PVPlacement(0,
			                               position[1],
			                               crystalLogic3_amo,
			                               "crystal3.physic_amo",
			                               crystalLogic_cry,
			                               false,
			                               0,
			                               true);
		}
		crystalLogic3_amo->SetUserLimits(stepLimit);
	} else {
		crystalChannelingData3->SetFilename(fECfileName[2][0]);
		crystalLogic3_cry->SetVerbose(1);

		if (fCrystalAmorphous[0]) {
			physAbsor3 = new G4PVPlacement(0,
			                               position[1],
			                               crystalLogic3_cry,
			                               "crystal3.physic_cry",
			                               crystalLogic_amo,
			                               false,
			                               0,
			                               true);
		} else {
			physAbsor3 = new G4PVPlacement(0,
			                               position[1],
			                               crystalLogic3_cry,
			                               "crystal3.physic_cry",
			                               crystalLogic_cry,
			                               false,
			                               0,
			                               true);
		}
	}

	// layer3 crystal
	Crystal4 = new G4ExtendedMaterial("crystal4.material",mat[3]);
	Crystal4->RegisterExtension(std::unique_ptr<G4CrystalExtension>(new G4CrystalExtension(Crystal4)));
	G4CrystalExtension* crystalExtension4 = (G4CrystalExtension*)Crystal4->RetrieveExtension("crystal");


	crystalExtension4->SetUnitCell(new G4CrystalUnitCell(5.43 * CLHEP::angstrom,
	                                                    5.43 * CLHEP::angstrom,
	                                                    5.43 * CLHEP::angstrom,
	                                                    CLHEP::halfpi,
	                                                    CLHEP::halfpi,
	                                                    CLHEP::halfpi,
	                                                    227));

	Crystal4->RegisterExtension(std::unique_ptr<G4ChannelingMaterialData>(new G4ChannelingMaterialData("channeling")));
	G4ChannelingMaterialData* crystalChannelingData4 = (G4ChannelingMaterialData*)Crystal4->RetrieveExtension("channeling");

	crystalLogic4_amo = new G4LogicalVolume (crystalSolid4, mat[3], "crystal4.logic_amo");
	crystalLogic4_cry = new G4LogicalCrystalVolume(crystalSolid4, Crystal4,"crystal4.logic_cry"); 

	if(fCrystalAmorphous[3]) {
		if(fCrystalAmorphous[0]) {
			physAbsor4 = new G4PVPlacement(0,
			                               position[2],
			                               crystalLogic4_amo,
			                               "crystal4.physic_amo",
			                               crystalLogic_amo,
			                               false,
			                               0,
			                               true);
		} else {
			physAbsor4 = new G4PVPlacement(0,
			                               position[2],
			                               crystalLogic4_amo,
			                               "crystal4.physic_amo",
			                               crystalLogic_cry,
			                               false,
			                               0,
			                               true);
		}
		crystalLogic4_amo->SetUserLimits(stepLimit);
	} else {
		crystalChannelingData4->SetFilename(fECfileName[3][0]);
		crystalLogic4_cry->SetVerbose(1);
		if (fCrystalAmorphous[0]) {
			physAbsor4 = new G4PVPlacement(0,
			                               position[2],
			                               crystalLogic4_cry,
			                               "crystal4.physic_cry",
			                               crystalLogic_amo,
			                               false,
			                               0,
			                               true);
		} else {
			physAbsor4 = new G4PVPlacement(0,
			                               position[2],
			                               crystalLogic4_cry,
			                               "crystal4.physic_cry",
			                               crystalLogic_cry,
			                               false,
			                               0,
			                               true);
		}
	}
	// layer4 crystal
	Crystal5 = new G4ExtendedMaterial("crystal5.material",mat[4]);
	Crystal5->RegisterExtension(std::unique_ptr<G4CrystalExtension>(new G4CrystalExtension(Crystal5)));
	G4CrystalExtension* crystalExtension5 = (G4CrystalExtension*)Crystal5->RetrieveExtension("crystal");


	crystalExtension5->SetUnitCell(new G4CrystalUnitCell(5.43 * CLHEP::angstrom,
	                                                    5.43 * CLHEP::angstrom,
	                                                    5.43 * CLHEP::angstrom,
	                                                    CLHEP::halfpi,
	                                                    CLHEP::halfpi,
	                                                    CLHEP::halfpi,
	                                                    227));

	Crystal5->RegisterExtension(std::unique_ptr<G4ChannelingMaterialData>(new G4ChannelingMaterialData("channeling")));
	G4ChannelingMaterialData* crystalChannelingData5 = (G4ChannelingMaterialData*)Crystal5->RetrieveExtension("channeling");

	crystalLogic5_amo = new G4LogicalVolume (crystalSolid5, mat[4], "crystal5.logic_amo");
	crystalLogic5_cry = new G4LogicalCrystalVolume(crystalSolid4, Crystal5,"crystal5.logic_cry"); 

	if (fCrystalAmorphous[4]) {
		if (fCrystalAmorphous[0]) {
			physAbsor5 = new G4PVPlacement(0,
			                               position[3],
			                               crystalLogic5_amo,
			                               "crystal5.physic_amo",
			                               crystalLogic_amo,
			                               false,
			                               0,
			                               true);
		} else {
			physAbsor5 = new G4PVPlacement(0,
			                               position[3],
			                               crystalLogic5_amo,
			                               "crystal5.physic_amo",
			                               crystalLogic_cry,
			                               false,
			                               0,
			                               true);
		}
		crystalLogic5_amo->SetUserLimits(stepLimit);
	} else {
		crystalChannelingData5->SetFilename(fECfileName[4][0]);
		crystalLogic5_cry->SetVerbose(1);
		if (fCrystalAmorphous[0]) {
			physAbsor5 = new G4PVPlacement(0, 
			                               position[3],
			                               crystalLogic5_cry,
			                               "crystal5.physic_amo",
			                               crystalLogic_amo,
			                               false,
			                               0,
			                               true);
		} else {
			physAbsor5 = new G4PVPlacement(0,
			                               position[3],
			                               crystalLogic5_cry,
			                               "crystal5.physic_amo",
			                               crystalLogic_cry,
			                               false,
			                               0,
			                               true);
		}
	}

	G4ChannelingMaterialData* matData;
	G4double LatConst1, LatConst2, LatConst3;
	if (!fCrystalAmorphous[0]) {
		matData = GetMatData(crystalLogic_cry);
		LatConst1 = matData->GetPot()->GetIntSp(0);
		LatConst2 = matData->GetPot()->GetIntSp(1);
		LatConst3 = matData->GetPot()->GetIntSp(2);
		SetLatticeConstants(0,G4ThreeVector(LatConst1,LatConst2,LatConst3));
	}
	if (!fCrystalAmorphous[1]) {
		matData = GetMatData(crystalLogic2_cry);
		LatConst1 = matData->GetPot()->GetIntSp(0);
		LatConst2 = matData->GetPot()->GetIntSp(1);
		LatConst3 = matData->GetPot()->GetIntSp(2);
		SetLatticeConstants(1,G4ThreeVector(LatConst1,LatConst2,LatConst3));
	}
	if (!fCrystalAmorphous[2]) {
		matData = GetMatData(crystalLogic3_cry);
		LatConst1 = matData->GetPot()->GetIntSp(0);
		LatConst2 = matData->GetPot()->GetIntSp(1);
		LatConst3 = matData->GetPot()->GetIntSp(2);
		SetLatticeConstants(2,G4ThreeVector(LatConst1,LatConst2,LatConst3));
	}
	if (!fCrystalAmorphous[3]) {
		matData = GetMatData(crystalLogic4_cry);
		LatConst1 = matData->GetPot()->GetIntSp(0);
		LatConst2 = matData->GetPot()->GetIntSp(1);
		LatConst3 = matData->GetPot()->GetIntSp(2);
		SetLatticeConstants(3,G4ThreeVector(LatConst1,LatConst2,LatConst3));
	}
	if (!fCrystalAmorphous[4]) {
		matData = GetMatData(crystalLogic5_cry);
		LatConst1 = matData->GetPot()->GetIntSp(0);
		LatConst2 = matData->GetPot()->GetIntSp(1);
		LatConst3 = matData->GetPot()->GetIntSp(2);
		SetLatticeConstants(4,G4ThreeVector(LatConst1,LatConst2,LatConst3));
	}
	return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{

	if (fCrystalAmorphous[0]) 
		crystalLogic = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal.logic_amo");
	else
		crystalLogic = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal.logic_cry");

	G4VSensitiveDetector* crystaldetector = new CrystalDetector("/crystaldetector",fECfileName[0][0],fECfileName[0][1],fECfileName[0][2],fCrystalAmorphous[0]);
	G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector);
	crystalLogic->SetSensitiveDetector(crystaldetector);


	if (fCrystalAmorphous[1]) 
		crystalLogic2 = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal2.logic_amo");
	else
		crystalLogic2 = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal2.logic_cry");

	G4VSensitiveDetector* crystaldetector2 = new CrystalDetector("/crystaldetector2",fECfileName[1][0],fECfileName[1][1],fECfileName[1][2],fCrystalAmorphous[1]);
	G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector2);
	crystalLogic2->SetSensitiveDetector(crystaldetector2);

	if (fCrystalAmorphous[2]) 
		crystalLogic3 = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal3.logic_amo");
	else
		crystalLogic3 = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal3.logic_cry");

	G4VSensitiveDetector* crystaldetector3 = new CrystalDetector("/crystaldetector3",fECfileName[2][0],fECfileName[2][1],fECfileName[2][2],fCrystalAmorphous[2]);
	G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector3);
	crystalLogic3->SetSensitiveDetector(crystaldetector3);

	if (fCrystalAmorphous[3]) 
		crystalLogic4 = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal4.logic_amo");
	else
		crystalLogic4 = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal4.logic_cry");

	G4VSensitiveDetector* crystaldetector4 = new CrystalDetector("/crystaldetector4",fECfileName[3][0],fECfileName[3][1],fECfileName[3][2],fCrystalAmorphous[3]);
	G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector4);
	crystalLogic4->SetSensitiveDetector(crystaldetector4);

	if (fCrystalAmorphous[4]) 
		crystalLogic5 = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal5.logic_amo");
	else
		crystalLogic5 = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal5.logic_cry");

	G4VSensitiveDetector* crystaldetector5 = new CrystalDetector("/crystaldetector5",fECfileName[4][0],fECfileName[4][1],fECfileName[4][2],fCrystalAmorphous[4]);
	G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector5);
	crystalLogic5->SetSensitiveDetector(crystaldetector5);

	G4LogicalVolume* ssdLogic = G4LogicalVolumeStore::GetInstance()->GetVolume("ssd.logic");
	G4VSensitiveDetector* telescope = new SensitiveDetector("/telescope");
	G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
	for (G4int i1=0; i1<4; i1++) {
		ssdLogic->SetSensitiveDetector(telescope);
	}

	if (use_xs_transformation) {
		G4cout << " Using XS TRANFORMATION " << G4endl;
		G4ChannelingOptrMultiParticleChangeCrossSection* xsecopr = new G4ChannelingOptrMultiParticleChangeCrossSection();
		xsecopr->AttachTo(crystalLogic);
		xsecopr->AttachTo(crystalLogic2);
		xsecopr->AttachTo(crystalLogic3);
		xsecopr->AttachTo(crystalLogic4);
		xsecopr->AttachTo(crystalLogic5);
	}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......