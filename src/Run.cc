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
/// \file Run.cc
/// \brief Implementation of the Run class
//
// $Id: Run.cc 71376 2013-06-14 07:44:50Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(),
	fDetector(det), fParticle(0), fEkin(0.),  fTrueRange(0.), fTrueRange2(0.),
	fProjRange(0.), fProjRange2(0.), fMap()
{
	fEnergyDeposit = fEnergyDeposit2 = 0.;
	fEnergyFlow = fEnergyFlow2 = 0.;
	NonIonEnergyDeposit = NonIonEnergyDeposit2 = 0.;
	Nsteps = Nsteps2 = 0.;
	theta = theta2 = 0.;
	TrakLenPrim = TrakLenPrim2 = 0.;
	TrakLenSec  = TrakLenSec2  = 0.;
	N_rec  = 0.;
	Th = 28.*eV;
	fTrueRange = fTrueRange2 = 0.;
	fProjRange = fProjRange2 = 0.;
	hit = 0.;
	partEmerging = 0.;
	detKinEn = detKinEn2 = 0.;

	for (G4int i=0; i<NUMB_MAX_LAYERS; ++i) {
		for (G4int j=0; j<3; ++j) {
			trackLenLayer[i][j] = 0;
			edepLayer[i][j] = 0;
			nielLayer[i][j] = 0;
			stepsLayer[i][j] = 0;
			secKinEnLayer[i][j] = 0;
			secDamEnLayer[i][j] = 0;
		}
		numbRecLayer[i] = 0;
	}

	projectedR = projectedR2 = 0.;
	entry_sd = 0.;
	total_step = 0.;
	entry_reach = 0.;
	RBSDepth = RBSDepth2 = 0.;
	counts = 0;
	angle_of_incidence = 0;
	primary_energy = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{
	// Important to clean up the map
	std::map<G4int, G4THitsMap<G4double>* >::iterator iter = fMap.begin();
	
	while (iter != fMap.end()) {
		delete iter->second;
		iter++;}
}
//========================================================
void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{
	fParticle = particle;
	fEkin = energy;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(const G4VProcess* process) 
{
	G4String procName = process->GetProcessName();
	std::map<G4String,G4int>::iterator it = fProcCounter.find(procName);
	if ( it == fProcCounter.end()) {
		fProcCounter[procName] = 1;
	} else {
		fProcCounter[procName]++; 
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleCount(G4String name, G4double Ekin)
{
	std::map<G4String, ParticleData>::iterator it = fParticleDataMap1.find(name);
	if (it == fParticleDataMap1.end()) {
		fParticleDataMap1[name] = ParticleData(1, Ekin, Ekin, Ekin);
	} else {
		ParticleData& data = it->second;
		data.fCount++;
		data.fEmean += Ekin;
		//update min max
		G4double emin = data.fEmin;
		if (Ekin < emin)
			data.fEmin = Ekin;
		G4double emax = data.fEmax;
		if (Ekin > emax)
			data.fEmax = Ekin;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleFlux(G4String name, G4double Ekin)
{
	std::map<G4String, ParticleData>::iterator it = fParticleDataMap2.find(name);
	if (it == fParticleDataMap2.end()) {
		fParticleDataMap2[name] = ParticleData(1, Ekin, Ekin, Ekin);
	} else {
		ParticleData& data = it->second;
		data.fCount++;
		data.fEmean += Ekin;
		//update min max
		G4double emin = data.fEmin;
		if (Ekin < emin)
			data.fEmin = Ekin;
		G4double emax = data.fEmax;
		if (Ekin > emax)
			data.fEmax = Ekin;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
	const Run* localRun = static_cast<const Run*>(run);
	//primary particle info
	fParticle = localRun->fParticle;
	fEkin     = localRun->fEkin;
	// accumulate sums
	fEnergyDeposit       += localRun->fEnergyDeposit;
	fEnergyDeposit2      += localRun->fEnergyDeposit2;
	fEnergyFlow          += localRun->fEnergyFlow;
	fEnergyFlow2         += localRun->fEnergyFlow2;
	NonIonEnergyDeposit  += localRun->NonIonEnergyDeposit;
	NonIonEnergyDeposit2 += localRun->NonIonEnergyDeposit2;
	N_rec                += localRun->N_rec; 
	Nsteps               += localRun->Nsteps;
	Nsteps2              += localRun->Nsteps2;
	theta                += localRun->theta;
	theta2               += localRun->theta2;
	TrakLenPrim          += localRun->TrakLenPrim;
	TrakLenPrim2         += localRun->TrakLenPrim2;
	TrakLenSec           += localRun->TrakLenSec;
	TrakLenSec2          += localRun->TrakLenSec2;

	fTrueRange           += localRun->fTrueRange;
	fTrueRange2          += localRun->fTrueRange2;
	fProjRange           += localRun->fProjRange;
	fProjRange2          += localRun->fProjRange2;

	for (G4int i=0; i<NUMB_MAX_LAYERS; ++i) {
		for (G4int j=0; j<2; ++j) {
			trackLenLayer[i][j] += localRun->trackLenLayer[i][j];
			edepLayer[i][j]     += localRun->edepLayer[i][j];
			nielLayer[i][j]     += localRun->nielLayer[i][j];
			stepsLayer[i][j]    += localRun->stepsLayer[i][j];
			secKinEnLayer[i][j] += localRun->secKinEnLayer[i][j];
			secDamEnLayer[i][j] += localRun->secDamEnLayer[i][j];
		}
		numbRecLayer[i]   += localRun->numbRecLayer[i];
	}

	primary_energy        += localRun->primary_energy;
	angle_of_incidence    += localRun->angle_of_incidence;
	hit                   += localRun->hit;
	partEmerging          += localRun->partEmerging;
	detKinEn              += localRun->detKinEn;
	detKinEn2             += localRun->detKinEn2;
	projectedR            += localRun->projectedR;
	projectedR2           += localRun->projectedR2;
	entry_sd              += localRun->entry_sd;
	entry_reach           += localRun->entry_reach;
	total_step            += localRun->total_step;
	RBSDepth              += localRun->RBSDepth;
	counts                += localRun->counts;
	// END

	//map: processes count
	std::map<G4String,G4int>::const_iterator itp;
	for (itp = localRun->fProcCounter.begin(); itp != localRun->fProcCounter.end(); ++itp ) {
		G4String procName = itp->first;
		G4int localCount = itp->second;
		if (fProcCounter.find(procName) == fProcCounter.end()) {
			fProcCounter[procName] = localCount;
		} else {
			fProcCounter[procName] += localCount;
		}
	}
	
	//map: created particles count
	std::map<G4String,ParticleData>::const_iterator itc;
	for (itc = localRun->fParticleDataMap1.begin(); itc != localRun->fParticleDataMap1.end(); ++itc) {
		G4String name = itc->first;
		const ParticleData& localData = itc->second;
		if (fParticleDataMap1.find(name) == fParticleDataMap1.end()) {
			fParticleDataMap1[name] = ParticleData(localData.fCount, localData.fEmean, localData.fEmin, localData.fEmax);
		} else {
			ParticleData& data = fParticleDataMap1[name];
			data.fCount += localData.fCount;
			data.fEmean += localData.fEmean;
			G4double emin = localData.fEmin;
			if (emin < data.fEmin)
				data.fEmin = emin;
			G4double emax = localData.fEmax;
			if (emax > data.fEmax)
				data.fEmax = emax; 
		}
	}
	
	//map: particles flux count
	std::map<G4String,ParticleData>::const_iterator itn;
	for (itn = localRun->fParticleDataMap2.begin(); itn != localRun->fParticleDataMap2.end(); ++itn) {
		G4String name = itn->first;
		const ParticleData& localData = itn->second;
		if (fParticleDataMap2.find(name) == fParticleDataMap2.end()) {
			fParticleDataMap2[name] = ParticleData(localData.fCount, localData.fEmean, localData.fEmin, localData.fEmax);
		} else {
			ParticleData& data = fParticleDataMap2[name];
			data.fCount += localData.fCount;
			data.fEmean += localData.fEmean;
			G4double emin = localData.fEmin;
			if (emin < data.fEmin)
				data.fEmin = emin;
			G4double emax = localData.fEmax;
			if (emax > data.fEmax)
				data.fEmax = emax; 
		}
	}
	const std::map< G4int, G4THitsMap<G4double>* >& localMap = localRun->fMap;
	std::map< G4int, G4THitsMap<G4double>* >::const_iterator iter = localMap.begin();
	for (; iter != localMap.end() ; ++iter)
		(*(fMap[iter->first])) += (*(iter->second));

	G4Run::Merge(run);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun() 
{
	G4int prec = 5, wid = prec + 2;  
	G4int dfprec = G4cout.precision(prec);

	G4Material* mat_matrix[5];
	G4double dens_matrix[5];

	for (G4int i=0; i<5; ++i) {
		mat_matrix[i] = fDetector->GetMaterialM(i);
		dens_matrix[i] = mat_matrix[i]->GetDensity();
	}

	G4cout << " ********************************************************** " << G4endl;
	G4cout << " ******************** FINISH OF RUN *********************** " << G4endl;

	if (numberOfEvent == 0) {
		G4cout.precision(dfprec);
		return;
	}

	G4String Particle = fParticle->GetParticleName();
	G4cout << "\n The run is " << numberOfEvent << " "<< Particle << " of "	<< G4BestUnit(primary_energy / numberOfEvent,"Energy") << " through "  
		<< G4BestUnit(fDetector->GetLength(0),"Length") << " length of " << mat_matrix[0]->GetName() << " (density: " << G4BestUnit(dens_matrix[0], "Volumic Mass") << ")" << G4endl;

	//frequency of processes
	G4cout << "\n Process calls frequency :" << G4endl;
	G4int index = 0;
	std::map<G4String,G4int>::iterator it;    
	for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
		G4String procName = it->first;
		G4int    count    = it->second;
		G4String space = " "; if (++index%3 == 0) space = "\n";
		G4cout << " " << std::setw(20) << procName << "="<< std::setw(7) << count << space;
	}
	G4cout << G4endl;
  
	//particles count
	G4cout << "\n List of generated particles:" << G4endl;
	std::map<G4String,ParticleData>::iterator itc;
	for (itc = fParticleDataMap1.begin(); itc != fParticleDataMap1.end(); itc++) { 
		G4String name = itc->first;
		ParticleData data = itc->second;
		G4int count = data.fCount;
		G4double eMean = data.fEmean/count;
		G4double eMin = data.fEmin;
		G4double eMax = data.fEmax;
			
		G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
			<< "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
			<< "\t( "  << G4BestUnit(eMin, "Energy")
			<< " --> " << G4BestUnit(eMax, "Energy") 
			<< ")" << G4endl;
	}

	// compute mean Energy deposited and rms
	G4int TotNbofEvents = numberOfEvent;
	fEnergyDeposit /= TotNbofEvents; fEnergyDeposit2 /= TotNbofEvents;
	G4double rmsEdep = GetVariation(fEnergyDeposit2, fEnergyDeposit);
	G4cout << "\n Mean energy deposit per event = "
			<< G4BestUnit(fEnergyDeposit,"Energy") << ";  rms = "
			<< G4BestUnit(rmsEdep,      "Energy") 
			<< G4endl;
	// compute mean Energy flow and rms
	fEnergyFlow /= TotNbofEvents; fEnergyFlow2 /= TotNbofEvents;
	G4double rmsEflow = GetVariation(fEnergyFlow2, fEnergyFlow);
	G4cout << " Mean energy flow per event    = "
			<< G4BestUnit(fEnergyFlow,"Energy") << ";  rms = "
			<< G4BestUnit(rmsEflow,   "Energy") 
			<< G4endl;

	//particles flux
	G4cout << "\n List of particles emerging from the absorber :" << G4endl;
	std::map<G4String,ParticleData>::iterator itn;
	for (itn = fParticleDataMap2.begin(); itn != fParticleDataMap2.end(); itn++) { 
		G4String name = itn->first;
		ParticleData data = itn->second;
		G4int count = data.fCount;
		G4double eMean = data.fEmean/count;
		G4double eMin = data.fEmin;
		G4double eMax = data.fEmax;
		G4double Eflow = data.fEmean/TotNbofEvents;
		G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
			<< "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
			<< "\t( "  << G4BestUnit(eMin, "Energy")
			<< " --> " << G4BestUnit(eMax, "Energy") 
			<< ") \tEflow/event = " << G4BestUnit(Eflow, "Energy") << G4endl;
	}

	//mean number of steps:
	Nsteps/=TotNbofEvents;  Nsteps2/=TotNbofEvents;
	G4double rmsSteps = GetVariation(Nsteps2, Nsteps);
	//scattering angle
	theta/=TotNbofEvents ; theta2/=TotNbofEvents;
	G4double rmsTheta = GetVariation(theta2, theta);
	//track length
	TrakLenPrim /= TotNbofEvents; TrakLenPrim2 /= TotNbofEvents;
	G4double rmsTLPrim = GetVariation(TrakLenPrim2, TrakLenPrim);
	G4double rmsTLSec = 0;
	// secondaries track length
	if (N_rec > 0) {
		TrakLenSec /= N_rec; TrakLenSec2 /= N_rec;
		rmsTLSec = GetVariation(TrakLenSec2, TrakLenSec);
	}
	// non ionion edep
	NonIonEnergyDeposit /= TotNbofEvents; NonIonEnergyDeposit2 /= TotNbofEvents;
	G4double rmsNonIon = GetVariation(NonIonEnergyDeposit2, NonIonEnergyDeposit);
	// AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
	//************************************************
	// --------------------------------------------
	G4double sum_Tl, sum_Tl2;
	for (G4int i=0; i<NUMB_MAX_LAYERS; ++i) {
		for (G4int j=0; j<2; ++j) {
			edepLayer[i][j] /= TotNbofEvents;
			trackLenLayer[i][j] /= TotNbofEvents;
			nielLayer[i][j] /= TotNbofEvents;
			stepsLayer[i][j] /= TotNbofEvents;
			secDamEnLayer[i][j] /= TotNbofEvents;
		}
		edepLayer[i][2] = GetVariation(edepLayer[i][1], edepLayer[i][0]);
		trackLenLayer[i][2] = GetVariation(trackLenLayer[i][1], trackLenLayer[i][0]);
		nielLayer[i][2] = GetVariation(nielLayer[i][1], nielLayer[i][0]);
		stepsLayer[i][2] = GetVariation(stepsLayer[i][1], stepsLayer[i][0]);
		secDamEnLayer[i][2] = GetVariation(secDamEnLayer[i][1], secDamEnLayer[i][0]);
		if (numbRecLayer[i] > 0) {
			secKinEnLayer[i][0] /= numbRecLayer[i];
			secKinEnLayer[i][1] /= numbRecLayer[i];
			secKinEnLayer[i][2] = GetVariation(secKinEnLayer[i][1], secKinEnLayer[i][0]);
		}
		sum_Tl += secDamEnLayer[i][0];
		sum_Tl2 += secDamEnLayer[i][1];
	}
	//*************************************************
	RBSDepth /= TotNbofEvents; RBSDepth2 /= TotNbofEvents;
	G4double rmsRBSDepth = GetVariation(RBSDepth2, RBSDepth);
	G4double rmsSum_Tl = GetVariation(sum_Tl2, sum_Tl);
	//..............................................................
	// effective length
	G4double length = TrakLenPrim;
	// total energy loss  
	G4double meandEdx  = fEnergyDeposit / length;
	// nuclear energy loss
	G4double meandEdx_nucl  = NonIonEnergyDeposit / length;
	// NIEL 
	G4double meandEdx_sumTL = sum_Tl / length;
	//[MeVcm2/g]
	G4double stopPower = meandEdx / dens_matrix[0];
	G4double stopPower_nucl = meandEdx_nucl / dens_matrix[0];
	G4double stopPower_sumTL = meandEdx_sumTL / dens_matrix[0];
	//mean free path & corss section 
	G4double freepath = TrakLenPrim / Nsteps;
	G4double er1 = rmsTLPrim / Nsteps;
	G4double er2 = freepath * rmsSteps / Nsteps;
	G4double rmsFreepath = std::sqrt(er1 * er1 + er2 * er2);

	G4double NA = mat_matrix[0]->GetTotNbOfAtomsPerVolume();
	G4double CrossSection = 1. / (NA * freepath); 
	G4double rmsCrossSection = rmsFreepath * CrossSection / freepath;
	// true and projected range from TestEm1
	fTrueRange /= TotNbofEvents; fTrueRange2 /= TotNbofEvents;
	G4double trueRms = GetVariation(fTrueRange2, fTrueRange);
	fProjRange /= TotNbofEvents; fProjRange2 /= TotNbofEvents;
	G4double projRms = GetVariation(fProjRange2, fProjRange);
	// projected range from testem11
	projectedR /= TotNbofEvents; projectedR2 /= TotNbofEvents;
	G4double rmsPR = GetVariation(projectedR2, projectedR);
	G4cout << "\n ============= Simulation statistics ==============\n";
	G4cout << "\n Primary Total track length in absorber:\n "
			<< G4BestUnit(TrakLenPrim,"Length") << " +- "
			<< G4BestUnit(rmsTLPrim,       "Length") << G4endl;

	G4cout << "\n Secondaries total track length in absorber:\n "
			<< G4BestUnit(TrakLenSec,"Length") << " +- "
			<< G4BestUnit(rmsTLSec,       "Length") << G4endl;

	G4cout << "\n ==================================================\n ";
	G4cout << "\n Primary particle statistics\n ";
	G4cout << "\n Mean Number of Steps:\n "<<Nsteps<< " +/- "<< rmsSteps<<G4endl;
									
	G4cout << "\n Mean Free Path :\n "<<G4BestUnit(freepath,"Length")<<  
					" +/- "<< G4BestUnit(rmsFreepath,"Length")<<G4endl;

	G4cout << "\n Mean Cross Section :\n "<<G4BestUnit(CrossSection,"Surface")<<
				" +/- "<<G4BestUnit(rmsCrossSection,"Surface")<<G4endl;

	G4cout << "\n Mean scattering angle :\n "<<G4BestUnit(theta,"Angle")<< " +/- "
				<< G4BestUnit(rmsTheta,"Angle")<<G4endl;
	
	G4cout << "\n Total energy deposit in absorber:\n "
			<< G4BestUnit(fEnergyDeposit,"Energy") << " +/- "
			<< G4BestUnit(rmsEdep,      "Energy") 
			<< G4endl;
	G4cout << "-----> dE/dx total= " << meandEdx/(MeV/cm) << " MeV/cm"
			<< "\t(" << stopPower/(MeV*cm2/g) << " MeV*cm2/g)"
			<< G4endl;


	G4cout << "\n Nuclear energy deposit in absorber:\n "
			<< G4BestUnit(NonIonEnergyDeposit,"Energy") << " +/- "
			<< G4BestUnit(rmsNonIon,      "Energy")
			<< G4endl;
	G4cout << "-----> dE/dx  nucl = " << meandEdx_nucl/(MeV/cm) << " MeV/cm"
			<< "\t(" << stopPower_nucl/(MeV*cm2/g) << " MeV*cm2/g)"
			<< G4endl;

	G4cout <<"\n NIEL in absorber (Th>"<<Th/eV <<" eV):\n "
			<< G4BestUnit(sum_Tl,"Energy") << " +/- "
			<< G4BestUnit(rmsSum_Tl,      "Energy")
			<< G4endl;
	G4cout << "-----> NIEL = " << meandEdx_sumTL/(MeV/cm) << " MeV/cm"
			<< "\t(" << stopPower_sumTL/(MeV*cm2/g) << " MeV*cm2/g)"
			<< G4endl;

	G4cout << "\n ===========================================================\n";
	G4cout << "\n true Range = " << G4BestUnit(fTrueRange,"Length")
			<< "   rms = "        << G4BestUnit(trueRms,  "Length");

	G4cout << "\n proj Range = " << G4BestUnit(fProjRange,"Length")
			<< "   rms = "        << G4BestUnit(projRms,  "Length");
	G4cout << "\n ===========================================================\n";
	G4cout << " ====================== Layers =============================";
	G4cout << "\n ===========================================================\n";
	G4cout << " \n"; 

	G4double widths[NUMB_MAX_LAYERS - 1], lengths[NUMB_MAX_LAYERS - 1], neigh_dist[NUMB_MAX_LAYERS - 2];
	for (G4int i=0; i<NUMB_MAX_LAYERS - 1; ++i) {
		widths[i] = fDetector->GetLength(i + 1);
		lengths[i] = fDetector->GetPosition(i).z();
	}
	G4double half_tot_len = fDetector->GetLength(0) / 2;

	for (G4int i=0; i<NUMB_MAX_LAYERS - 2; ++i) {
		neigh_dist[i] = (lengths[i + 1] - widths[i + 1] / 2) - (lengths[i] + widths[i]/2);
	}

	G4double distance_to_surface[NUMB_MAX_LAYERS - 1];
	for (G4int i=0; i<NUMB_MAX_LAYERS - 1; ++i) {
		distance_to_surface[i] = half_tot_len + lengths[i] - (widths[i] / 2);
	}
	if (distance_to_surface[0] < 1e-5 *nm)
		distance_to_surface[0] = 0;

	G4double final_dist = half_tot_len - lengths[3] - widths[3] / 2;
	G4double TFU_main[NUMB_MAX_LAYERS];
	G4double tfu_norm = 1e+15;
	for (G4int i=0; i<NUMB_MAX_LAYERS - 1; ++i) {
		TFU_main[i] = (mat_matrix[i + 1]->GetTotNbOfAtomsPerVolume() /(1/cm3) * fDetector->GetLength(i + 1) /cm) / tfu_norm;
	}
	TFU_main[4] = (mat_matrix[0]->GetTotNbOfAtomsPerVolume() /(1/cm3) * final_dist /cm) / tfu_norm;

	G4double TFU_mid[NUMB_MAX_LAYERS - 1];
	TFU_mid[0] = (mat_matrix[0]->GetTotNbOfAtomsPerVolume() /(1/cm3) * distance_to_surface[0]) / tfu_norm;
	for (G4int i=1; i<NUMB_MAX_LAYERS - 1; ++i) {
		TFU_mid[i] = (mat_matrix[0]->GetTotNbOfAtomsPerVolume() /(1/cm3) * neigh_dist[i - 1] /cm) / tfu_norm;
	}

	G4cout << " >>>>>>>>>>>>>>>>> MAIN LAYER <<<<<<<<<<<<<<<<<<<<<<<<\n";
	G4cout << " Total thickness = " << G4BestUnit(half_tot_len * 2, "Length") << " and surface layer thickness " << TFU_mid[0] << " [TFU] " << G4endl;
	for (G4int i=0; i<NUMB_MAX_LAYERS; ++i) {
		G4cout << " ************************************* " << G4endl;
		G4cout << " Layer # " << i << G4endl;
		G4cout << " Thickness " << G4BestUnit(fDetector->GetLength(i), "Length") << " or " << TFU_main[i] << " TFU " << G4endl; // pakeisti del pirmo layeriaus i mid, nuo pirmo palikti main???
		G4cout << " Layer dimensions : " << G4BestUnit(fDetector->GetSize(i), "Length") << G4endl;
		if (i != 0) {
			G4cout << " Position of the layer: " << G4BestUnit(fDetector->GetPosition(i - 1), "Length") << G4endl;
			G4cout << " Distance from the surface: " << G4BestUnit(distance_to_surface[i - 1], "Length") << G4endl;
		}
		G4cout << " Material: " << fDetector->GetMaterial(i) << G4endl;
		G4cout << " Density = " << dens_matrix[i] /(g/cm3) << " g/cm3 "  << ", or " << mat_matrix[i]->GetTotNbOfAtomsPerVolume() /(1/cm3) << " [at3] "  << G4endl;
		G4cout << " Composition: " << G4endl;
		const G4double *atomDensVector = mat_matrix[i]->GetVecNbOfAtomsPerVolume();
		G4double NoOfElements_material = mat_matrix[i]->GetNumberOfElements();
		for (int j=0; j<NoOfElements_material; j++) {
			G4cout << " " << mat_matrix[i]->GetElement(j)->GetName() << " = " << atomDensVector[j] /(1/cm3) << " [cm-3] " << G4endl;
		}
		G4cout << " -------------------- Particle parameters -------------------- " << G4endl;
		G4cout << " Ionizing energy loss : \t\t" << G4BestUnit(edepLayer[i][0], "Energy") << " +/- " << G4BestUnit(edepLayer[i][2], "Energy") << G4endl;
		G4cout << " Non ionizing energy loss [NIEL] : \t" << G4BestUnit(nielLayer[i][0], "Energy") << " +/- " << G4BestUnit(nielLayer[i][2], "Energy")  << G4endl;
		G4cout << " Mean number of steps : \t\t" << stepsLayer[i][0] << " +/- " << stepsLayer[i][2] << G4endl;
		G4cout << " Primary track length : \t\t" << G4BestUnit(trackLenLayer[i][0], "Length") << " +/- " << G4BestUnit(trackLenLayer[i][2],"Length") << G4endl;
		G4cout << " Mean dE/dx : \t\t\t\t" << edepLayer[i][0] / trackLenLayer[i][0] /(MeV/cm) << " MeV/cm " << G4endl;
		G4cout << " Mean NIEL dE/dx: \t\t\t" << nielLayer[i][0] / trackLenLayer[i][0] /(MeV/cm) << " MeV/cm " << G4endl;
		G4cout << " Mean free path: \t\t\t" <<G4BestUnit(trackLenLayer[i][0] / stepsLayer[i][0],"Length")<<  G4endl;
		G4cout << " Steps per length [nm-1]: \t\t" << stepsLayer[i][0] / ((half_tot_len * 2) /nm) << G4endl;
		G4cout << " # of secondaries: \t\t\t" << numbRecLayer[i] << " \t per primary : " << (double)numbRecLayer[i] / TotNbofEvents << G4endl;
		G4cout << " Mean E of secondaries: \t\t" << G4BestUnit(secKinEnLayer[i][0], "Energy") << " +/- " << G4BestUnit(secKinEnLayer[i][2], "Energy") << G4endl;
		G4cout << " Mean damage en of secondaries: \t" << G4BestUnit(secDamEnLayer[i][0], "Energy") << " +/- " << G4BestUnit(secDamEnLayer[i][2], "Energy") << G4endl;
		G4cout << " \n"; 
	}

	G4cout << " ==================================================\n";
	G4cout << " Max Step size : \t\t\t" <<G4BestUnit(fDetector->GetMaxStep(),"Length")<<  G4endl;
	G4double rbs_angle = fDetector->GetRBSAngle() /degree;
	G4cout << " RBS detector angle: \t\t\t" << rbs_angle << " degrees " << G4endl;
	G4cout << " Average Max RBS depth: \t\t" << G4BestUnit(RBSDepth / counts, "Length") << " +/- " << G4BestUnit(rmsRBSDepth / counts, "Length") << G4endl;
	G4cout << " Use of sigmacalc: \t\t\t" ;
	if (fDetector->GetSigmaCalc()) 
		G4cout << "ENABLED " << G4endl;
	else
		G4cout << "DISABLED " << G4endl;
	G4cout << " RBS evaluation was: \t\t\t" ;
	if (fDetector->GetRBSCalc()) 
		G4cout << "ENABLED " << G4endl;
	else
		G4cout << "DISABLED " << G4endl;
	G4cout << " Constant scattering angle: \t\t" ;
	if (fDetector->GetConstAngle()) 
		G4cout << "ENABLED " << G4endl;
	else
		G4cout << "DISABLED " << G4endl;
	G4cout << " RBS minimum ROI : \t\t\t" << G4BestUnit(fDetector->GetRBSROImin(), "Energy") << G4endl;
	G4cout << " \n";

	G4double chan = hit / 4;
	G4double totchan = chan / TotNbofEvents;
	G4double totem = chan / partEmerging;
	// kinetic energy of particles reaching the 4th detector
	// TODO check whether this is correct!!!!
	G4double rmssum_kinen = 0.0;
	if (chan > 0) {
		detKinEn /= hit;
		detKinEn2 /= hit;
		rmssum_kinen = GetVariation(detKinEn2, detKinEn);
	}
	G4cout << " ==================================================\n";
	G4cout << " Projected range = " << G4BestUnit(projectedR, "Length") << " +- " << G4BestUnit(rmsPR, "Length")<< G4endl;
	G4cout << " Rotation of the sample: " << fDetector->GetAngles() /degree << " degrees " << G4endl;
	G4cout << " # of particles that reaches last detector " << hit / 4 << G4endl;
	G4cout << " part of total particles " << std::fixed << std::setprecision(2) << totchan * 100 << " % " <<G4endl;
	G4cout << " part of emitted particles " << std::setprecision(2) << totem * 100 << " % "<< G4endl;
	G4cout << " Mean KinEn of particles reaching 4th detector = "<<G4BestUnit(detKinEn, "Energy") <<" +/- "  <<G4BestUnit(rmssum_kinen, "Energy") << G4endl; 
	G4double no_of_steps_per_particle = entry_sd / numberOfEvent;
	G4double no_of_reach_per_particle = entry_reach / numberOfEvent;
	G4cout << " Number of steps : " << no_of_steps_per_particle << G4endl; 
	G4cout << " Number of entries reached histo : " << no_of_reach_per_particle << G4endl;
	G4cout << " Average number of entries per step : " << no_of_reach_per_particle / no_of_steps_per_particle << G4endl;
	G4double total_step_length = total_step / numberOfEvent;
	G4cout << " Total step : " << G4BestUnit(total_step_length, "Length") << G4endl;
	G4cout << " Average step length: " << G4BestUnit(total_step_length / no_of_steps_per_particle, "Length") << G4endl;
	G4cout << " ****************************************************** " << G4endl;
	G4cout << " ******************** GEOMETRY ************************ " << G4endl;
	G4cout << " ****************************************************** " << G4endl;

	G4double inc_angle = (angle_of_incidence / numberOfEvent);
	G4cout << " Angles in degrees: incidence = " << inc_angle / degree << "; detection = " << fDetector->GetRBSAngle() / degree 
		<< "; scattering angle = " << (fDetector->GetRBSAngle() - inc_angle) /degree << "; exit angle = " << (M_PI - fDetector->GetRBSAngle() - inc_angle) /degree << G4endl;

	G4cout << " Materials: \n " << G4endl;
	G4cout << " X = " << mat_matrix[0]->GetName() << G4endl;
	G4cout << " A = " << mat_matrix[1]->GetName() << G4endl;
	G4cout << " B = " << mat_matrix[2]->GetName() << G4endl;
	G4cout << " C = " << mat_matrix[3]->GetName() << G4endl;
	G4cout << " D = " << mat_matrix[4]->GetName() << G4endl;
	G4cout << " \n" << G4endl;

	G4cout << " TFU UNITS " << G4endl;
	G4cout << " " << TFU_mid[0] << " | " << TFU_main[0] << " | " << TFU_mid[1] << " | " << TFU_main[1] << " | " << TFU_mid[2] << " | " 
	       << TFU_main[2] << " | " << TFU_mid[3] << " | " << TFU_main[3] << " | " << TFU_main[4] << G4endl;
	
	G4cout << " " << G4BestUnit(distance_to_surface[0], "Length") <<  " | " << G4BestUnit(widths[0], "Length") 
		   << " | " << G4BestUnit(neigh_dist[0], "Length") <<  " | " << G4BestUnit(widths[1], "Length") 
		   << " | " << G4BestUnit(neigh_dist[1], "Length") <<  " | " << G4BestUnit(widths[2], "Length") 
		   << " | " << G4BestUnit(neigh_dist[2], "Length") <<  " | " << G4BestUnit(widths[3], "Length") 
		   << " | " << G4BestUnit(final_dist, "Length") << G4endl;
	G4cout << " XXXXXX->|<-AAAAAA->|<-XXXXXX->|<-BBBBBB->|<-XXXXXX->|<-CCCCCC->|<-XXXXXX->|<-DDDDDD->|<-XXXXXX" << G4endl;
	G4cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " << G4endl;
	G4cout << " |<-------------- Substrate " << G4BestUnit(half_tot_len * 2, "Length") << " --------------------------------------------------->| " << G4endl;
	G4cout << "                                                        " << G4endl;

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
	if (analysisManager->IsActive()) {
		G4double fac; 
		G4double unit = analysisManager->GetH1Unit(2); 
		G4double binWidth = analysisManager->GetH1Width(2) * unit;
		fac = (1. / (numberOfEvent * binWidth)) *(mm/MeV);
		analysisManager->ScaleH1(2, fac);

		for (G4int ih=3; ih<15; ih++) {
			unit = analysisManager->GetH1Unit(ih);
			binWidth = analysisManager->GetH1Width(ih) * unit;
		}

		fac = (1. / TrakLenPrim);
		analysisManager->ScaleH1(15, fac);
		fac = (1. / numberOfEvent);
		analysisManager->ScaleH1(18, fac);
	
		binWidth = analysisManager->GetH1Width(17) * analysisManager->GetH1Unit(17);
		fac = (1. / (numberOfEvent * binWidth)) *(mm/keV);
		analysisManager->ScaleH1(17, fac);
		binWidth = analysisManager->GetH1Width(19) * analysisManager->GetH1Unit(19);
		fac = (1. / (TrakLenPrim * binWidth));
		analysisManager->ScaleH1(19, fac);

		for (G4int ih=20; ih<36; ih++) {
			G4double binW = analysisManager->GetH1Width(ih);
			G4double ave_step = total_step_length / no_of_steps_per_particle;
			G4double RBS_norm_dist = 0.1 *nm;
			G4double norm = ave_step / RBS_norm_dist;
			G4double exponent = exp(1);
			fac = (1 / (exponent * norm * no_of_steps_per_particle * numberOfEvent * binW)); // latest, 2021-02-23
			analysisManager->ScaleH1(ih, fac);
		}
		binWidth = analysisManager->GetH1Width(16) * analysisManager->GetH1Unit(16);
		G4double fcc = (1. / (binWidth * numberOfEvent)) *(cm);
		analysisManager->ScaleH1(16,fcc);

		binWidth = analysisManager->GetH1Width(17) * analysisManager->GetH1Unit(17);
			fcc = (1. / (binWidth * numberOfEvent)) *(cm/eV);	
		analysisManager->ScaleH1(17,fcc);
	}
	//remove all contents in fProcCounter, fCount 
	fProcCounter.clear();
	fParticleDataMap2.clear();

	//restore default format
	G4cout.precision(dfprec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
