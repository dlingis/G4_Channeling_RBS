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
#define MAT_SAMPLE    "G4_Si"
#define DETSZ_SAMPLE  G4ThreeVector(50.*mm, 50.*mm, 1.*mm)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():
fECfileName{EC_SAMPLE, EC_SAMPLE, EC_SAMPLE, EC_SAMPLE, EC_SAMPLE, EC_SAMPLE,
			EC_SAMPLE, EC_SAMPLE, EC_SAMPLE, EC_SAMPLE, EC_SAMPLE, EC_SAMPLE,
			EC_SAMPLE, EC_SAMPLE, EC_SAMPLE},
stepLimit(0), fMaterialName{MAT_SAMPLE, MAT_SAMPLE, MAT_SAMPLE, MAT_SAMPLE, MAT_SAMPLE},
fSizes{SZ_SAMPLE, SZ_SAMPLE, SZ_SAMPLE, SZ_SAMPLE, SZ_SAMPLE},
fAngles(POS_SAMPLE), fWorldMaterial("G4_Galactic"),
fDetectorMaterialName(""), fDetectorSizes(DETSZ_SAMPLE), fDetectorDistance{-20. *cm, -19. *cm, +19. *cm, +20. *cm, +40. *cm}, 
fCrystalAmorphous{0,0,0,0,0}, maxStep(0.), rbs_angle(160.*degree), rbs_step(0), sec_material_ratio(0.),
detector_resolution(10.*keV), position{POS_SAMPLE, POS_SAMPLE, POS_SAMPLE, POS_SAMPLE}, 
sigma_calc(0), rbs_calc(0), material_mixing(0), enable_custom_material(0), use_const_angle(false), enable_fwhm_calc(0.),
use_xs_transformation(0), histo_tracking(false), gauss_counter(5), material_for_mix(0),
element1(""), element2(""), element3(""), dead_material_name(MAT_SAMPLE), mix_material_name(MAT_SAMPLE),
custom_density(0.), part1(0), part2(0), part3(0), dead_thickness(100.), solidAngle(1.), rbs_roi_min(50. *keV),
array_for_histos{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
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

	G4Box *crystalSolid[LAYERS];
	G4String crname;
	for (uint8_t i=0; i<LAYERS; ++i) {
		crname = "crystal.solid" + std::to_string(i);
		crystalSolid[i] = new G4Box(crname,
	                                fSizes[i].x()/2.,
	                                fSizes[i].y()/2.,
	                                fSizes[i].z()/2.);
	}

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

	for (G4int i=0; i<LAYERS; ++i) {
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

	for (G4int i=0; i<LAYERS; i++) {
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
  for(int i = 0; i<LAYERS; i++) {
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

	// setup crystalline extension
	G4ExtendedMaterial* crystals[LAYERS];
	G4ChannelingMaterialData* crystalChannelingDatas[LAYERS];
	G4CrystalExtension* cry_ext[LAYERS];
	G4double sz = 5.43 * CLHEP::angstrom;
	G4double half_pi = CLHEP::halfpi;
	G4String cryname;
	for (uint8_t i=0; i<LAYERS; ++i) {
		// setup only if material is crystalline
		if (!fCrystalAmorphous[i]) {
			cryname = "crystal.material" + std::to_string(i);
			crystals[i] = new G4ExtendedMaterial(cryname, mat[i]);
			crystals[i]->RegisterExtension(std::unique_ptr<G4CrystalExtension>(new G4CrystalExtension(crystals[i])));
			cry_ext[i] = (G4CrystalExtension*)crystals[i]->RetrieveExtension("crystal");
			cry_ext[i]->SetUnitCell(new G4CrystalUnitCell(sz, sz, sz, half_pi, half_pi, half_pi, 227));
			crystals[i]->RegisterExtension(std::unique_ptr<G4ChannelingMaterialData>(new G4ChannelingMaterialData("channeling")));
			crystalChannelingDatas[i] = (G4ChannelingMaterialData*)crystals[i]->RetrieveExtension("channeling");
		}
	}

	// building logical and physical volumes
	G4String log_v_name, phys_v_name;
	for (uint8_t i=0; i<LAYERS; ++i) {
		log_v_name = "crystal.logic" + std::to_string(i);
		phys_v_name = "crystal.physic" + std::to_string(i);
		// construct logicals and link to crystal or amorphous
		if (fCrystalAmorphous[i]) {
			log_v_name += "_amo";
			crystal_logic_v_amo[i] = new G4LogicalVolume(crystalSolid[i], mat[i], log_v_name);
			crystal_logic_v[i] = crystal_logic_v_amo[i];
		} else {
			log_v_name += "_cry";
			crystal_logic_v_cry[i] = new G4LogicalCrystalVolume(crystalSolid[i], crystals[i], log_v_name);
			crystal_logic_v[i] = crystal_logic_v_cry[i];
		}
		// construct physicals
		if (i == 0)
			crystal_physical_v[i] = new G4PVPlacement(rot, G4ThreeVector(), crystal_logic_v[i], phys_v_name, worldLogic, false, 0, true);
		else
			crystal_physical_v[i] = new G4PVPlacement(rot, position[i - 1], crystal_logic_v[i], phys_v_name, crystal_logic_v[0], false, 0, true);
		// set limitations
		if (fCrystalAmorphous[i])
			crystal_logic_v[i]->SetUserLimits(stepLimit);
		else
			crystalChannelingDatas[i]->SetFilename(fECfileName[i][0]);
		G4cout << " Logical name = " << crystal_logic_v[i]->GetName() << G4endl;
		G4cout << " Physical name = " << crystal_physical_v[i]->GetName() << G4endl;
		// G4cout << " Pot = " << crystalChannelingDatas[i]->GetPot()->GetMax() << G4endl;
	}

	return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void PrintMap(const std::map<G4String, std::vector<G4LogicalVolume*> >& bmap) {
//     // Iterate through the map
//     for (const auto& pair : bmap) {
//         const G4String& key = pair.first;
//         const std::vector<G4LogicalVolume*>& volumes = pair.second;

//         // Print the key (G4String)
//         std::cout << "Key: " << key << std::endl;

//         // Print the logical volumes associated with the key
//         for (size_t i = 0; i < volumes.size(); ++i) {
//             if (volumes[i] != nullptr) {
//                 std::cout << "  Volume " << i << ": " << volumes[i]->GetName() << " (Address: " << volumes[i] << ")" << std::endl;
//             } else {
//                 std::cout << "  Volume " << i << ": nullptr" << std::endl;
//             }
//         }
//     }
// }

void DetectorConstruction::ConstructSDandField()
{
	G4String name, det_name;
	G4VSensitiveDetector* crystal_det[LAYERS];
	// PrintMap(G4LogicalVolumeStore::GetInstance()->GetMap());
	for (uint8_t i=0; i<LAYERS; ++i) {
		name = "crystal.logic" + std::to_string(i);
		det_name = "/crystaldetector" + std::to_string(i);
		if (fCrystalAmorphous[i])
			name += "_amo";
		else
			name += "_cry";
		crystal_logic_v[i] = G4LogicalVolumeStore::GetInstance()->GetVolume(name);
		crystal_det[i] = new CrystalDetector(det_name, fECfileName[i][0],fECfileName[i][1],fECfileName[i][2],fCrystalAmorphous[i]);
		G4SDManager::GetSDMpointer()->AddNewDetector(crystal_det[i]);
		crystal_logic_v[i]->SetSensitiveDetector(crystal_det[i]);
	}

	G4LogicalVolume* ssdLogic = G4LogicalVolumeStore::GetInstance()->GetVolume("ssd.logic");
	G4VSensitiveDetector* telescope = new SensitiveDetector("/telescope");
	G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
	for (G4int i1=0; i1<4; i1++) {
		ssdLogic->SetSensitiveDetector(telescope);
	}

	if (use_xs_transformation) {
		G4cout << " Using XS TRANFORMATION " << G4endl;
		G4ChannelingOptrMultiParticleChangeCrossSection* xsecopr = new G4ChannelingOptrMultiParticleChangeCrossSection();
		for (uint8_t i=0; i<LAYERS; ++i)
			xsecopr->AttachTo(crystal_logic_v[i]);
	}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......