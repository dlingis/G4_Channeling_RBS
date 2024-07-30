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

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1
#endif

#define LAYERS 5

#include "G4VUserDetectorConstruction.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "DetectorConstructionMessenger.hh"

#include "G4LogicalCrystalVolume.hh"

#include "globals.hh"

#include "G4UserLimits.hh"
#include "G4ExtendedMaterial.hh"
#include "G4ChannelingMaterialData.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class DetectorConstruction : public G4VUserDetectorConstruction
{
	public:
		DetectorConstruction();
		~DetectorConstruction();
		
		void DefineMaterials();
		G4VPhysicalVolume* Construct();

	private:
		void ConstructSDandField();

	private:
		DetectorConstructionMessenger* fMessenger;
		G4VPhysicalVolume* physAbsor;
		G4VPhysicalVolume* physAbsor2;
		G4VPhysicalVolume* physAbsor3;
		G4VPhysicalVolume* physAbsor4;
		G4VPhysicalVolume* physAbsor5;
		G4VPhysicalVolume* worldPhysical; 

		G4String fECfileName[LAYERS][3];
		G4UserLimits* stepLimit; // pointer to user step limits
		G4String fMaterialName[LAYERS];
		G4ThreeVector fSizes[LAYERS];
	
		G4ThreeVector fAngles;
		G4String fWorldMaterial;

		G4LogicalVolume*        crystalLogic_amo;
		G4LogicalVolume*        crystalLogic2_amo;
		G4LogicalVolume*        crystalLogic3_amo;
		G4LogicalVolume*        crystalLogic4_amo;
		G4LogicalVolume*        crystalLogic5_amo;

		G4LogicalCrystalVolume* crystalLogic_cry;
		G4LogicalCrystalVolume* crystalLogic2_cry;
		G4LogicalCrystalVolume* crystalLogic3_cry;
		G4LogicalCrystalVolume* crystalLogic4_cry;
		G4LogicalCrystalVolume* crystalLogic5_cry;

		G4ExtendedMaterial*     Crystal;
		G4ExtendedMaterial*     Crystal2;
		G4ExtendedMaterial*     Crystal3;
		G4ExtendedMaterial*     Crystal4;
		G4ExtendedMaterial*     Crystal5;

		G4LogicalVolume*        crystalLogic;
		G4LogicalVolume*        crystalLogic2;
		G4LogicalVolume*        crystalLogic3;
		G4LogicalVolume*        crystalLogic4;
		G4LogicalVolume*        crystalLogic5;
	
	public:
		G4String GetEC(G4int la, G4int el) {return fECfileName[la][el];}
		void SetEC(G4int la, G4int el, G4String name) {fECfileName[la][el] = name;}

		G4String GetMaterial(G4int i) {return fMaterialName[i];}
		void SetMaterial(G4int i,G4String aString) {fMaterialName[i] = aString;}

		G4ThreeVector GetSize(G4int i) {return fSizes[i];}
		void SetSize(G4int i, G4ThreeVector a3vec) {fSizes[i] = a3vec;}

		G4ThreeVector GetAngles() {return fAngles;}
		void SetAngles(G4ThreeVector a3vec) {fAngles = a3vec;}

		G4VPhysicalVolume* GetIntAbsorber(G4int i) {
								if(i == 0) return physAbsor;
								else if(i == 1) return physAbsor2;
								else if(i == 2) return physAbsor3;
								else if(i == 3) return physAbsor4;
								else if(i == 4) return physAbsor5;
								else return physAbsor;
								};

		const G4VPhysicalVolume* GetWorld()                 {return worldPhysical;};
		G4Material* GetMaterialM(G4int i)                   {return mat[i];};
		G4Material* GetMaterialComponents(G4int i, G4int j) {return mat_components[i][j];};
		G4double GetLength(G4int i)                         {return fSizes[i].z();};
		G4ThreeVector GetDimensions(G4int i )               {return fSizes[i];};

		void SetPosition(G4int i, G4ThreeVector a3vec)      {position[i] = a3vec;};
		G4ThreeVector GetPosition(G4int i)                  {return position[i];};

		void SetWorldMaterial(G4String aString)             {fWorldMaterial = aString;}
		G4String GetWorldMaterial()                         {return fWorldMaterial;}

		void SetDetectorMaterial(G4String aString)          {fDetectorMaterialName = aString;}
		G4String GetDetectorMaterial()                      {return fDetectorMaterialName;}

		void SetDetectorSizes(G4ThreeVector a3vec)          {fDetectorSizes = a3vec;}
		G4ThreeVector GetDetectorSizes()                    {return fDetectorSizes;}

		void SetDetectorDistance(G4int aInt,G4double aDouble){fDetectorDistance[aInt] = aDouble;}
		G4double GetDetectorDistance(G4int aInt)             {return fDetectorDistance[aInt];}

		void SetAmorphous(G4int la, bool a)                 {fCrystalAmorphous[la] = a;}
		G4int GetAmorphous(G4int la)                        {return fCrystalAmorphous[la];}
	
		G4double GetMaxStep()                               {return maxStep;}
		void SetMaxStep(G4double stp)                       {maxStep = stp;}

		G4double GetRBSAngle()                              {return rbs_angle;}
		void SetRBSAngle(G4double ang)                      {rbs_angle = ang;}

		G4bool GetSigmaCalc()                               {return sigma_calc;}
		void SetSigmaCal(bool a)                            {sigma_calc = a;}

		G4bool GetRBSCalc()                                 {return rbs_calc;}
		void SetRBSCalc(bool a)                             {rbs_calc = a;}

		void SetEnLossStep(G4double a)                      {rbs_step = a;}
		G4double GetEnLossStep()                            {return rbs_step;}
		
		void SetMixing(bool a)                              {material_mixing = a;}
		G4bool GetMaterialMixing()                          {return material_mixing;}
		
		void SetMixingRatio(G4double a)                     {sec_material_ratio = a;}
		G4double GetMixingRatio()                           {return sec_material_ratio;}
		
		void SetDetectorResolution(G4double a)              {detector_resolution = a;}
		G4double GetDetectorResolution()                    {return detector_resolution;}
		
		void SetGaussCounter(G4int a)                       {gauss_counter = a;}
		G4int GetGaussCounter()                             {return gauss_counter;}
		// custom material
		void SetCustomMaterial(bool a)                      {enable_custom_material = a;}
		G4bool GetCustomMaterial()                          {return enable_custom_material;}
		
		void SetCustomElement1(G4String a)                  {element1 = a;}
		G4String GetCustomElement1()                        {return element1;}
		
		void SetCustomElement2(G4String a)                  {element2 = a;}
		G4String GetCustomElement2()                        {return element2;}
		
		void SetCustomElement3(G4String a)                  {element3 = a;}
		G4String GetCustomElement3()                        {return element3;}
		
		void SetCustomMaterialDensity(G4double a)           {custom_density = a;}
		G4double GetCustomMaterialDensity()                 {return custom_density;}
		
		void SetCustomElement1Part(G4double a)              {part1 = a;}
		G4double GetCustomElement1Part()                    {return part1;}
		
		void SetCustomElement2Part(G4double a)              {part2 = a;}
		G4double GetCustomElement2Part()                    {return part2;}
		
		void SetCustomElement3Part(G4double a)              {part3 = a;}
		G4double GetCustomElement3Part()                    {return part3;}
		
		void SetDeadLayer(G4String a)                       {dead_material_name = a;};
		G4String GetDeadLayer()                             {return dead_material_name;};

		void SetDeadLayerThickness(G4double a)              {dead_thickness = a;};
		G4double GetDeadLayerThickness()                    {return dead_thickness;};
		G4Material* GetDeadLayerMaterial()                  {return dead_material;};
		
		void SetSolidAngle(G4double a)                      {solidAngle = a;};
		G4double GetSolidAngle()                            {return solidAngle;};
		
		// set the number of material which should be mixed
		void SetMaterialForMix(G4int a)                     {material_for_mix = a;};
		G4int GetMaterialForMix()                           {return material_for_mix;};
		
		// set the material 2 for mix
		G4String GetMixingMaterial()                        {return mix_material_name;}
		void SetMixingMaterial(G4String a)                  {mix_material_name = a;}
		
		// calculate detector fwhm
		G4bool GetCalcFWHM()                                {return enable_fwhm_calc;}
		void SetCalcFWHM(bool a)                            {enable_fwhm_calc = a;}
		
		G4double GetRBSROImin()                             {return rbs_roi_min;}
		void SetRBSROImin(G4double a)                       {rbs_roi_min = a;}
		
		G4bool GetConstAngle()                              {return use_const_angle;}
		void UseConstAngle(bool a)                          {use_const_angle = a;}

		void UseXStranformation(bool a)                     {use_xs_transformation = a;}

		void SetArrayForHistos(G4int aInt,G4double aDouble) {array_for_histos[aInt] = aDouble;}
		void CopyArrayForDepth(G4double array[])            {memcpy(array, array_for_histos, sizeof(array_for_histos));}

		void SetTrackingHistos(bool a)                      {histo_tracking = a;}
		G4bool GetTrackingHistos()                          {return histo_tracking;}
	private:
		G4String                fDetectorMaterialName;
		G4ThreeVector           fDetectorSizes;
		G4double                fDetectorDistance[5];

		G4bool                  fCrystalAmorphous[5];

		G4double                maxStep, rbs_angle, rbs_step, sec_material_ratio, detector_resolution;
		G4ThreeVector           position[4];

		G4bool                  sigma_calc, rbs_calc, material_mixing, enable_custom_material, use_const_angle, enable_fwhm_calc, use_xs_transformation, histo_tracking;
		G4int                   gauss_counter, material_for_mix;

		G4String                element1, element2, element3, dead_material_name, mix_material_name;
		G4double                custom_density, part1, part2, part3, dead_thickness, solidAngle, rbs_roi_min;
		
		G4Material*             mat[5];
		G4Material*             dead_material;
		G4Material*             mixing_material;
		G4Material*             mat_components[5][5];
		G4double                array_for_histos[12];

	G4ChannelingMaterialData* GetMatData(G4LogicalVolume* aLV) {
		if(aLV->IsExtended() == true){
			G4ExtendedMaterial* aEM = (G4ExtendedMaterial*) aLV->GetMaterial();
			return (G4ChannelingMaterialData*) aEM->RetrieveExtension("channeling");
		}
		else{
			return nullptr;
		}
	}
};