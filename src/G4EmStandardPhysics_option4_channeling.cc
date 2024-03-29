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
//---------------------------------------------------------------------------
//
// ClassName:   G4EmStandardPhysics_option4_channeling
//
// Author:      V.Ivanchenko 28.09.2012 from Option3 physics constructor
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4EmStandardPhysics_option4_channeling.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"
#include "G4UrbanMscModel.hh"
#include "G4DummyModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4Generator2BS.hh"
#include "G4Generator2BN.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eplusAnnihilation.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include "G4ePairProduction.hh"

#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4hBremsstrahlungModel.hh"
#include "G4hPairProductionModel.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmModelActivator.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"



// jomajo

#include "G4StepLimiter.hh"
#include "G4ProcessManager.hh"

#include "G4IonCoulombScatteringModel.hh"


// Wrapper
#include "XWrapperDiscreteProcess.hh"
#include "XWrapperContinuousDiscreteProcess.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronElastic.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4ChipsProtonElasticXS.hh"
#include "G4ChipsElasticModel.hh"
#include "G4IonPhysics.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4FTFBuilder.hh"
#include "G4HadronicInteraction.hh"
#include "G4BuilderType.hh"
//#include "G4ProtonInelasticProcess.hh"

// ATIMA
#include "G4AtimaEnergyLossModel.hh"
#include "G4AtimaFluctuations.hh"
#include "G4BraggIonModel.hh"
#include "G4LindhardSorensenIonModel.hh"

#include "G4BetheBlochModel.hh"
#include "G4BraggModel.hh"


// modified atima
#include "G4ionIonisation_mod.hh"


#include "G4PhysicsFreeVector.hh"

#include <fstream>
#include <sstream>
#include <iomanip>

#include "G4eSingleCoulombScatteringModel.hh"
#include "G4IonParametrisedLossModel_mod.hh"


//#include "G4IonParametrisedLossModel_mod.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysics_option4_channeling::G4EmStandardPhysics_option4_channeling(G4int ver, const G4String&)
	: G4VPhysicsConstructor("G4EmStandard_opt4_channeling"), verbose(ver)
{
	G4EmParameters* param = G4EmParameters::Instance();
	param->SetDefaults();
	param->SetVerbose(verbose);
	param->SetMinEnergy(100*eV);// orgas
	param->SetMaxEnergy(100*MeV);
	param->SetLowestElectronEnergy(10*eV);
	//param->NumberOfBins(40000);//pridetas -> nuimti?
	param->SetNumberOfBinsPerDecade(200);
	param->ActivateAngularGeneratorForIonisation(true);
	param->SetMscThetaLimit(0.0);
	param->SetFluo(true);
	param->SetAuger(true);
	param->SetPixe(true);
	param->SetBuildCSDARange(true); // nebuvo, mano pridetas 
	param->SetLossFluctuations(true);
	SetPhysicsType(bElectromagnetic);
	/*
	param->SetDefaults();
	param->SetMinEnergy(1000*eV);
	param->SetMaxEnergy(100*MeV);
	param->SetNumberOfBinsPerDecade(10);
	param->SetBuildCSDARange(true);
	param->SetLossFluctuations(true);  
	param->SetMaxEnergyForCSDARange(20*MeV);
	param->SetStepFunctionLightIons(0.1, 20*CLHEP::um);
	param->SetStepFunctionIons(0.1, 1*CLHEP::um);
	param->SetUseMottCorrection(true);  
	SetPhysicsType(bElectromagnetic);  
	*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysics_option4_channeling::~G4EmStandardPhysics_option4_channeling()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysics_option4_channeling::ConstructParticle()
{
	// gamma
	G4Gamma::Gamma();
	// leptons
	G4Electron::Electron();
	G4Positron::Positron();
	G4MuonPlus::MuonPlus();
	G4MuonMinus::MuonMinus();
	// mesons
	G4PionPlus::PionPlusDefinition();
	G4PionMinus::PionMinusDefinition();
	G4KaonPlus::KaonPlusDefinition();
	G4KaonMinus::KaonMinusDefinition();
	// barions
	G4Proton::Proton();
	G4AntiProton::AntiProton();
	// ions
	G4Deuteron::Deuteron();
	G4Triton::Triton();
	G4He3::He3();
	G4Alpha::Alpha();
	G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysics_option4_channeling::ConstructProcess()
{
	if(verbose > 1) {
		G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
	}
	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
	// muon & hadron bremsstrahlung and pair production
	G4MuBremsstrahlung* mub = new G4MuBremsstrahlung();
	G4MuPairProduction* mup = new G4MuPairProduction();
	G4hBremsstrahlung* pib = new G4hBremsstrahlung();
	G4hPairProduction* pip = new G4hPairProduction();
	G4hBremsstrahlung* kb = new G4hBremsstrahlung();
	G4hPairProduction* kp = new G4hPairProduction();
	G4ePairProduction* ee = new G4ePairProduction();
	// muon & hadron multiple scattering
	G4CoulombScattering* muss = new G4CoulombScattering();
	G4CoulombScattering* piss = new G4CoulombScattering();
	G4CoulombScattering* kss = new G4CoulombScattering();
	// energy limits for e+- scattering models
	//G4double highEnergyLimit = 100*MeV;
	// energy limits for e+- ionisation models
	G4double penEnergyLimit = 1*MeV;
	// nuclear stopping
	G4NuclearStopping* pnuc = new G4NuclearStopping();

	// Add standard EM Processes
	auto myParticleIterator=GetParticleIterator();
	myParticleIterator->reset();
	while( (*myParticleIterator)() ){
		G4ParticleDefinition* particle = myParticleIterator->value();
		G4String particleName = particle->GetParticleName();
		G4ProcessManager* pmanager = particle->GetProcessManager();

		if (particleName == "gamma") {
			// Photoelectric
			G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
			G4VEmModel* theLivermorePEModel = new G4LivermorePhotoElectricModel();
			pe->SetEmModel(theLivermorePEModel,1);
			ph->RegisterProcess(pe, particle);
			// Compton scattering
			G4ComptonScattering* cs = new G4ComptonScattering;
			cs->SetEmModel(new G4KleinNishinaModel(),1);
			G4VEmModel* theLowEPComptonModel = new G4LowEPComptonModel();
			theLowEPComptonModel->SetHighEnergyLimit(20*MeV);
			cs->AddEmModel(0, theLowEPComptonModel);
			ph->RegisterProcess(cs, particle);
			// Gamma conversion
			G4GammaConversion* gc = new G4GammaConversion();
			G4VEmModel* thePenelopeGCModel = new G4PenelopeGammaConversionModel();
			thePenelopeGCModel->SetHighEnergyLimit(1*GeV);
			gc->SetEmModel(thePenelopeGCModel,1);
			ph->RegisterProcess(gc, particle);
			// Rayleigh scattering
			ph->RegisterProcess(new G4RayleighScattering(), particle);
		} else if (particleName == "e-") {
			// multiple scattering
			G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
			G4CoulombScattering* ss = new G4CoulombScattering();
			ss->SetEmModel(ssm, 1); 
			// ionisation
			G4eIonisation* eIoni = new G4eIonisation();
			eIoni->SetStepFunction(0.2, 100*um);
			G4PenelopeIonisationModel* pen = new G4PenelopeIonisationModel();
			pen->SetHighEnergyLimit(penEnergyLimit);
			eIoni->AddEmModel(0, pen, new G4UniversalFluctuation());
			// bremsstrahlung
			G4eBremsstrahlung* brem = new G4eBremsstrahlung();
			G4SeltzerBergerModel* br1 = new G4SeltzerBergerModel();
			G4eBremsstrahlungRelModel* br2 = new G4eBremsstrahlungRelModel();
			br1->SetAngularDistribution(new G4Generator2BS());
			br2->SetAngularDistribution(new G4Generator2BS());
			brem->SetEmModel(br1,1);
			brem->SetEmModel(br2,2);
			br2->SetLowEnergyLimit(GeV);
			// register processes
			ph->RegisterProcess(eIoni, particle);
			ph->RegisterProcess(brem, particle);
			ph->RegisterProcess(ee, particle);
			ph->RegisterProcess(ss, particle);
		} else if (particleName == "e+") {
			ph->RegisterProcess(new G4eMultipleScattering(), particle);
			G4eIonisation* eIoni = new G4eIonisation();
			eIoni->SetStepFunction(0.1, 100*um);
			ph->RegisterProcess(eIoni, particle);
			ph->RegisterProcess(new G4eBremsstrahlung(), particle);
			ph->RegisterProcess(new G4eplusAnnihilation(), particle);
		} else if (particleName == "mu+" || particleName == "mu-") {
			G4MuIonisation* muIoni = new G4MuIonisation();
			muIoni->SetStepFunction(0.2, 50*um);
			ph->RegisterProcess(muIoni, particle);
			ph->RegisterProcess(mub, particle);
			ph->RegisterProcess(mup, particle);
			ph->RegisterProcess(muss, particle);
		} else if (particleName == "alpha" || particleName == "He3") {
			// ionisation
			G4ionIonisation* ionIoni = new G4ionIonisation();
			ionIoni->SetEmModel(new G4IonParametrisedLossModel());
			ionIoni->SetStepFunction(0.1, 0.02*mm); // org 10 um
			XWrapperContinuousDiscreteProcess *ionIoni_wrapper = new XWrapperContinuousDiscreteProcess("Alpha_Ionization");
			ionIoni_wrapper->RegisterProcess(ionIoni, -1);
			//ionIoni_wrapper->RegisterProcess(ionIoni,-1,1); // same
			//ionIoni_wrapper->RegisterProcess(ionIoni,0,2); // new addition, kad pridetu ir amorfiniame
			ph->RegisterProcess(ionIoni_wrapper, particle);
			// coulomb scattering
			G4CoulombScattering* ecs = new G4CoulombScattering();
			G4IonCoulombScatteringModel* mod = new G4IonCoulombScatteringModel();
			ecs->AddEmModel(0,mod);
			ecs->SetBuildTableFlag(false);
			XWrapperDiscreteProcess *ecs_wrapper = new XWrapperDiscreteProcess("Alpha_CoulombScattering");
			ecs_wrapper->RegisterProcess(ecs, 1, 1);
			ecs_wrapper->RegisterProcess(ecs, 0, 2);
			ph->RegisterProcess(ecs_wrapper, particle);
			ph->RegisterProcess(pnuc, particle);
			/*
			G4ionIonisation* ionIoni = new G4ionIonisation();
			ionIoni->SetEmModel(new G4IonParametrisedLossModel());
			ionIoni->SetStepFunction(0.1, 0.02*mm); // org 10 um
			ph->RegisterProcess(ionIoni, particle);
			// coulomb scattering
			G4CoulombScattering* ecs = new G4CoulombScattering();
			G4IonCoulombScatteringModel* mod = new G4IonCoulombScatteringModel();
			ecs->AddEmModel(0,mod);
			ecs->SetBuildTableFlag(false);
			ph->RegisterProcess(ecs, particle);
			ph->RegisterProcess(pnuc, particle);
			*/
			pmanager->AddProcess(new G4StepLimiter, -1, -1, 3);
		} else if (particleName == "GenericIon") {
			// Original ICRU
			/*
			G4ionIonisation_mod* ionIonis = new G4ionIonisation_mod();
			ionIonis->SetEmModel(new G4IonParametrisedLossModel());
			ionIonis->SetStepFunction(0.1, 0.02*mm);
			ionIonis->SetMinKinEnergy(100.*keV);
			ph->RegisterProcess(ionIonis,particle);
			*/
			G4IonParametrisedLossModel_mod* ionpar = new G4IonParametrisedLossModel_mod();
			ionpar->UseSRIMData(1);	// set to 0 to use ICRU data
			G4ionIonisation* ionIonis = new G4ionIonisation();
			ionIonis->SetEmModel(ionpar,0);
			ionIonis->SetStepFunction(0.1, 0.02*mm); // org 10 um
			ionIonis->SetMinKinEnergy(100.*keV);
			//ph->RegisterProcess(ionIonis, particle); 
			// chann
			XWrapperContinuousDiscreteProcess *ionIoni_wrapper = new XWrapperContinuousDiscreteProcess("Ion_Ionization");
			ionIoni_wrapper->RegisterProcess(ionIonis,-1); // test
			ph->RegisterProcess(ionIoni_wrapper, particle);
			//*****************************************
			// SCATTERING 
			//*****************************************	
			// coulomb scattering
			G4CoulombScattering* ecs = new G4CoulombScattering();
			G4IonCoulombScatteringModel* mod = new G4IonCoulombScatteringModel();
			mod->SetLowEnergyLimit(100.*keV);
			mod->SetRecoilThreshold(28.*eV);
			ecs->AddEmModel(0,mod);
			//ecs->SetBuildTableFlag(false);
			//ph->RegisterProcess(ecs, particle);
			XWrapperDiscreteProcess *ecs_wrapper = new XWrapperDiscreteProcess("Ion_CoulombScattering");
			ecs_wrapper->RegisterProcess(ecs, 1, 1);
			ecs_wrapper->RegisterProcess(ecs, 0, 2);
			ph->RegisterProcess(ecs_wrapper, particle);
			//ph->RegisterProcess(pnuc,particle);
			//ph->RegisterProcess(new G4hMultipleScattering(), particle);
			pmanager->AddProcess(new G4StepLimiter, 1, -1, 3);
		} else if (particleName == "pi+" || particleName == "pi-" ) {
			//G4hMultipleScattering* pimsc = new G4hMultipleScattering();
			G4hIonisation* hIoni = new G4hIonisation();
			hIoni->SetStepFunction(0.2, 50*um);
			ph->RegisterProcess(hIoni, particle);
			ph->RegisterProcess(pib, particle);
			ph->RegisterProcess(pip, particle);
			ph->RegisterProcess(piss, particle);
		} else if (particleName == "kaon+" || particleName == "kaon-" ) {
			//G4hMultipleScattering* kmsc = new G4hMultipleScattering();
			G4hIonisation* hIoni = new G4hIonisation();
			hIoni->SetStepFunction(0.2, 50*um);
			ph->RegisterProcess(hIoni, particle);
			ph->RegisterProcess(kb, particle);
			ph->RegisterProcess(kp, particle);
			ph->RegisterProcess(kss, particle);

		} else if (particleName == "proton" || particleName == "anti_proton") {

			// ionisation
			G4hIonisation* hIoni = new G4hIonisation();
			hIoni->SetStepFunction(0.1, 2*um); //org 2 um
			XWrapperContinuousDiscreteProcess *hIoni_wrapper = new XWrapperContinuousDiscreteProcess("prot_Ionization");
			hIoni_wrapper->RegisterProcess(hIoni,-1);
			ph->RegisterProcess(hIoni_wrapper, particle);
			// bremsstrahlung
			G4hBremsstrahlung* hBrem = new G4hBremsstrahlung();
			XWrapperContinuousDiscreteProcess *hBrem_wrapper = new XWrapperContinuousDiscreteProcess("prot_Bremsstrahlung");
			hBrem_wrapper->RegisterProcess(hBrem,-1);
			ph->RegisterProcess(hBrem_wrapper, particle);
			// pair production
			G4hPairProduction* hPair = new G4hPairProduction();
			XWrapperContinuousDiscreteProcess* hPair_wrapper = new XWrapperContinuousDiscreteProcess("prot_PairProduction");
			hPair_wrapper->RegisterProcess(hPair,-1);
			ph->RegisterProcess(hPair_wrapper, particle);
			// coulomb scattering
			G4CoulombScattering* ecs = new G4CoulombScattering();
			ecs->SetBuildTableFlag(false);
			//ph->RegisterProcess(ecs,particle);
			XWrapperDiscreteProcess *ecs_wrapper = new XWrapperDiscreteProcess("prot_CoulombScattering");
			ecs_wrapper->RegisterProcess(ecs, 1, 1);
			ecs_wrapper->RegisterProcess(ecs, 0, 2);
			ph->RegisterProcess(ecs_wrapper, particle);
			//ph->RegisterProcess(pnuc, particle); // pridetas 2022-02-18
			pmanager->AddProcess(new G4StepLimiter,-1,-1,3);
		} else if (particleName == "deuteron" || particleName == "anti_deuteron"){
			// ionisation
			G4hIonisation* hIoni = new G4hIonisation();
			//hIoni->SetStepFunction(0.1, 2*um);
			XWrapperContinuousDiscreteProcess *hIoni_wrapper = new XWrapperContinuousDiscreteProcess("deuteron_Ionization");
			hIoni_wrapper->RegisterProcess(hIoni,-1);
			ph->RegisterProcess(hIoni_wrapper, particle);
			// coulomb scattering
			G4CoulombScattering* ecs = new G4CoulombScattering();
			ecs->SetBuildTableFlag(false);
			//ph->RegisterProcess(ecs,particle);
			XWrapperDiscreteProcess *ecs_wrapper = new XWrapperDiscreteProcess("deuteron_CoulombScattering");
			ecs_wrapper->RegisterProcess(ecs, 1, 1);
			ecs_wrapper->RegisterProcess(ecs, 0, 2);
			ph->RegisterProcess(ecs_wrapper, particle);
			/*
			// bremsstrahlung
			G4hBremsstrahlung* hBrem = new G4hBremsstrahlung();
			XWrapperContinuousDiscreteProcess *hBrem_wrapper =
			new XWrapperContinuousDiscreteProcess("deuteron_Bremsstrahlung");
			hBrem_wrapper->RegisterProcess(hBrem,-1);
			ph->RegisterProcess(hBrem_wrapper, particle);
			ph->RegisterProcess(hBrem, particle);

			// pair production
			G4hPairProduction* hPair = new G4hPairProduction();
			XWrapperContinuousDiscreteProcess* hPair_wrapper =
			new XWrapperContinuousDiscreteProcess("deuteron_PairProduction");
			hPair_wrapper->RegisterProcess(hPair,-1);
			ph->RegisterProcess(hPair_wrapper, particle);
			*/
			//ph->RegisterProcess(pnuc, particle);
			pmanager->AddProcess(new G4StepLimiter, -1, -1, 3);
		} else if (particleName == "B+" ||
				particleName == "B-" ||
				particleName == "D+" ||
				particleName == "D-" ||
				particleName == "Ds+" ||
				particleName == "Ds-" ||
				particleName == "anti_He3" ||
				particleName == "anti_alpha" ||
				//particleName == "anti_deuteron" ||
				particleName == "anti_lambda_c+" ||
				particleName == "anti_omega-" ||
				particleName == "anti_sigma_c+" ||
				particleName == "anti_sigma_c++" ||
				particleName == "anti_sigma+" ||
				particleName == "anti_sigma-" ||
				particleName == "anti_triton" ||
				particleName == "anti_xi_c+" ||
				particleName == "anti_xi-" ||
				//particleName == "deuteron" ||
				particleName == "lambda_c+" ||
				particleName == "omega-" ||
				particleName == "sigma_c+" ||
				particleName == "sigma_c++" ||
				particleName == "sigma+" ||
				particleName == "sigma-" ||
				particleName == "tau+" ||
				particleName == "tau-" ||
				particleName == "triton" ||
				particleName == "xi_c+" ||
				particleName == "xi-" ) {
			ph->RegisterProcess(new G4CoulombScattering(), particle);
			ph->RegisterProcess(new G4hIonisation(), particle);
			//ph->RegisterProcess(pnuc, particle);
		}
	}
	// Nuclear stopping
	pnuc->SetMaxKinEnergy(MeV);
/*
	G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
	de->SetFluo(true);
	de->SetAuger(true);
	de->SetPIXE(true);
	G4LossTableManager::Instance()->SetAtomDeexcitation(de);
*/
	G4EmModelActivator mact(GetPhysicsName());

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
