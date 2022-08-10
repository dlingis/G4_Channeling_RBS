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

#include "G4Channeling.hh"
#include "G4ChannelingMessenger.hh"
#include "Randomize.hh"
 
#include "G4ChannelingTrackData.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"

#include "G4LambdacPlus.hh"
#include "G4PhysicsModelCatalog.hh"
#include "ExExChParticleUserInfo.hh"

G4Channeling::G4Channeling():
G4VDiscreteProcess("channeling"),
fChannelingID(G4PhysicsModelCatalog::GetModelID("model_channeling")),
fMinimumEnergy(1.*CLHEP::keV),  
fMaximumMomentumRatio(0.01), 
fTimeStepMin(0.),
fTimeStepMax(0.),
fTransverseVariationMax(2.E-2 * CLHEP::angstrom),
chan_step(0),//,fCritAngle(0*degree)
density_limit(0),
v11_algo(0),
org_step_size(0),
use_step_size(0),
step_size_value(1.*nm),
mfp_value(1.*nm),
fChMessenger(0)
{
    fChMessenger = new G4ChannelingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Channeling::~G4Channeling(){
    //delete fChMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ChannelingTrackData* G4Channeling::GetTrackData(const G4Track& aTrack){
    G4ChannelingTrackData* trackdata =
        (G4ChannelingTrackData*)(aTrack.GetAuxiliaryTrackInformation(fChannelingID));
    if(trackdata == nullptr){
        trackdata = new G4ChannelingTrackData();
        aTrack.SetAuxiliaryTrackInformation(fChannelingID,trackdata);
    }
    return trackdata;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Channeling::GetEF(const G4Track& aTrack,
                         G4ThreeVector& pos,
                         G4ThreeVector& out){
    out = G4ThreeVector((GetMatData(aTrack)->GetEFX()->GetEC(pos)),
                        (GetMatData(aTrack)->GetEFY()->GetEC(pos)),
                        0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Channeling::GetEF(G4ChannelingMaterialData* materialData,
                         G4ThreeVector& pos,
                         G4ThreeVector& out){
    out = G4ThreeVector((materialData->GetEFX()->GetEC(pos)),
                        (materialData->GetEFY()->GetEC(pos)),
                        0.);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Channeling::PosToLattice(G4StepPoint* step,G4ThreeVector& pos){
    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(step->GetTouchable());

    pos -= theTouchable->GetTranslation();
    pos = ((*theTouchable->GetRotation()).inverse())(pos);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



G4bool G4Channeling::UpdateParameters(const G4Track& aTrack){

    G4ChannelingMaterialData* matData = GetMatData(aTrack);
    G4LogicalCrystalVolume* aLCV = (G4LogicalCrystalVolume*)(aTrack.GetVolume()->GetLogicalVolume());
    
    G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();
    G4StepPoint* preStepPoint = aTrack.GetStep()->GetPreStepPoint();

    G4ThreeVector posPre = preStepPoint->GetPosition();
    G4ThreeVector posPost = postStepPoint->GetPosition();	//pozicija world atskaitos sistemoje

    aLCV->RotateToLattice(posPre);
    aLCV->RotateToLattice(posPost);
    /*
    G4cout << " CRIT ANGLE = " << GetCriticalAngle2(aTrack)/degree << " crit angle 2 = " << GetCritAngleLimitAxial(aTrack)/degree<< G4endl;
    G4cout << " Crit angle v3 = " << GetCriticalAngle3(aTrack)/degree << G4endl;
    G4cout << " CRIT ENERGY = " << GetCriticalEnergy(aTrack)/MeV << G4endl;
    */
    G4double integrationLimit = std::fabs(posPost.z() - posPre.z());

    if(integrationLimit > 0.){

        //----------------------------------------
        // Initialize particle variables
        //---------
        G4double beta = aTrack.GetVelocity()/CLHEP::c_light;
        G4double mass = aTrack.GetDynamicParticle()->GetDefinition()->GetPDGMass();
        //----------------------------------------
        // Get the momentum in the world reference 
        // frame and rotate to the solid reference frame
        //----------------------------------------
        G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
        G4ThreeVector momWorld = aTrack.GetStep()->GetPreStepPoint()->GetMomentum();
        G4ThreeVector mom = (*theTouchable->GetRotation())(momWorld);

        // randomly generates a position of the particle to evatuate the transverse energy
        // if transverse energy (potential + kinetic energy) is lower than maximum potential
        G4ThreeVector preChanPos;
        /*
        G4double kinetic_energy = (mom.x()*mom.x()+mom.y()*mom.y()+mom.z()*mom.z())/(2*aTrack.GetParticleDefinition()->GetPDGMass());
        G4cout << " Kinetic energy = " << kinetic_energy/MeV << G4endl;
        */

        if(GetTrackData(aTrack)->GetPosCh().x() == DBL_MAX){
		// check whether it's the first step in crystal, then randomly generate position for potential evaluation
            if(enable_rechanneling == 1) {
                if(GetInfo(aTrack)->GetLastChannelingPosition().x() == DBL_MAX) {
                        G4double positX = G4UniformRand() * matData->GetPot()->GetIntSp(0);
                        G4double positY = G4UniformRand() * matData->GetPot()->GetIntSp(1);
                        preChanPos = G4ThreeVector(positX,positY,0.);
                }
                // if the particle was in channeling, use the new position of step in World reference system minus the last World position channeled, and add the difference to the last know channeling position to set the particle in the crystal reference system
                else {
                    G4ThreeVector chanPos = GetInfo(aTrack)->GetLastChannelingPosition();
                    G4ThreeVector newWorldPos = posPre - GetInfo(aTrack)->GetLastChannelingWorldPositionPost();
                    preChanPos = chanPos + newWorldPos;
                }
            } else {
                G4double positX = G4UniformRand() * matData->GetPot()->GetIntSp(0);
                G4double positY = G4UniformRand() * matData->GetPot()->GetIntSp(1);
                preChanPos = G4ThreeVector(positX,positY,0.);
            }
        } else
            preChanPos = GetTrackData(aTrack)->GetPosCh();

        G4double momentas = mom.x()*mom.x()+mom.y()*mom.y();
        G4double momentas_z = mom.z()*mom.z();

    	//G4double vTotalEnergy = aTrack.GetStep()->GetPreStepPoint()->GetTotalEnergy();
	    //G4double kinetic_energy = ((mom.x()*mom.x()+mom.y()*mom.y())/(2*mass));

	    //G4double transverse_energy = potential + kinetic_energy; 

        //G4double potential_maximum = GetMatData(aTrack)->GetPot()->GetMax();
        //G4double potential_at_pos = GetMatData(aTrack)->GetPot()->GetEC(preChanPos);

        G4ThreeVector mom_dir = aTrack.GetStep()->GetPreStepPoint()->GetMomentumDirection();
        G4ThreeVector mom_rotated = (*theTouchable->GetRotation())(mom_dir);

        //G4double momentas_dir = mom_dir.x()*mom_dir.x() + mom_dir.y()*mom_dir.y();
        //G4double momentas_dir_z = mom_dir.z()*mom_dir.z();

        //G4double angle_non_rot = atan(sqrt(momentas_dir/momentas_dir_z));

        //G4int potential_limit = (potential_at_pos > 0.90*potential_maximum);
        G4int angle_lim = (atan(sqrt(momentas/momentas_z)) > 1.1*GetCriticalAngle2(aTrack));

        //G4cout << " Critical angle = " << GetCriticalAngle2(aTrack)/degree << G4endl;

        if(angle_lim ){//&& potential_limit) { //potential_limit &&
            //G4cout << " angle " << angle_non_rot/degree << G4endl;
            GetTrackData(aTrack)->Reset();
	        GetInfo(aTrack)->Reset();
            return false;
        }
        
        
        /*
	    if (atan(sqrt(momentas/momentas_z))/degree > GetCriticalAngle2(aTrack)/degree ) {
            GetTrackData(aTrack)->Reset();
	        GetInfo(aTrack)->Reset();
            return false;
        } 
        */
        //----------------------------------------
        // Get the m omentum in the solid reference 
        // frame and rotate to the crystal reference frame
        //----------------------------------------
        aLCV->RotateToLattice(mom);
        //----------------------------------------
        // Take the position stored in the track data.
        // If the particle enters the crystal,
        // the position in the channel is randomly
        // generated using a uniform distribution
        //----------------------------------------

        G4ThreeVector pos = preChanPos;

        G4double step=0., stepTot=0.;
        G4double nud =0., eld    =0.;
        G4double efx =0., efy    =0.;
        G4double nud_temp =0., eld_temp    =0.;

        G4double Z = GetParticleDefinition(aTrack)->GetPDGCharge();
        
        const G4double oneSixth = 1./6.;
        G4ThreeVector posk1,posk2,posk3,posk4,posk5,posk6;
        G4ThreeVector momk1,momk2,momk3,momk4,momk5,momk6;
        G4ThreeVector pos_temp, efxy;

        do{
            //----------------------------------------
            // Limit the variable step length for the
            // integration via the selected algorithm
            // and update variables for the integration
            //----------------------------------------
            UpdateIntegrationStep(matData,aTrack,mom,step);

            if(step + stepTot > integrationLimit){
                step = integrationLimit - stepTot;
            }

            //----------------------------------------
            // Function integration algorithm
            // 4th Order Runge-Kutta
            //----------------------------------------
            //G4cout << " step = " << step/nm << G4endl;
            GetEF(matData,pos,efxy);
            posk1 = step / mom.z() * mom;
            if(v11_algo)
                momk1 = step / beta * Z * efxy;
            else
                momk1 = step / beta * Z * efxy * 0.5;
            //if(isBent) momk1.setX(momk1.x() - step * mom.z() * beta / (matData->GetBR(pos)).x());
            
            GetEF(matData,pos_temp = pos + posk1 * 0.5,efxy);
            posk2 = step / mom.z() * (mom + momk1 * 0.5);
            if(v11_algo)
                momk2 = step / beta * Z * efxy;
            else
                momk2 = step / beta * Z * efxy * 0.5;
            //if(isBent) momk2.setX(momk2.x() - step * mom.z() * beta / (matData->GetBR(pos_temp)).x());

            GetEF(matData,pos_temp = pos + posk2 * 0.5,efxy);
            posk3 = step / mom.z() * (mom + momk2 * 0.5);
            if(v11_algo)
                momk3 = step / beta * Z * efxy;
            else
                momk3 = step / beta * Z * efxy * 0.5;
            //if(isBent) momk3.setX(momk3.x() - step * mom.z() * beta / (matData->GetBR(pos_temp)).x());
            
            GetEF(matData,pos_temp = pos + posk3,efxy);
            posk4 = step / mom.z() * (mom + momk3);
            if(v11_algo)
                momk4 = step / beta * Z * efxy;
            else
                momk4 = step / beta * Z * efxy * 0.5;
            //if(isBent) momk4.setX(momk4.x() - step * mom.z() * beta / (matData->GetBR(pos_temp)).x());

            pos = pos + oneSixth * (posk1 + 2.*posk2 + 2.*posk3 + posk4);
            mom = mom + oneSixth * (momk1 + 2.*momk2 + 2.*momk3 + momk4);
       
            //----------------------------------------
            // Update the total step and the electron
            // and nuclei density experienced by
            // the particle during its motion
            //----------------------------------------

            stepTot += step;

            nud_temp = matData->GetNuD()->GetEC(pos);
            eld_temp = matData->GetElD()->GetEC(pos);
				
            if(nud_temp < 0.) {nud_temp = 0.;}
            if(eld_temp < 0.) {eld_temp = 0.;}

            nud += (step * nud_temp);
            eld += (step * eld_temp);

            efx += (step * matData->GetEFX()->GetEC(pos));
            efy += (step * matData->GetEFY()->GetEC(pos));
            //----------------------------------------
        } while(stepTot<integrationLimit);
 
        nud /= stepTot;
        eld /= stepTot;
/*
        if(nud < 1.e-2) nud = 1.e-2;
        if(eld < 1.e-2) eld = 1.e-2;
*/
        if(nud < density_limit) nud = density_limit;
        if(eld < density_limit) eld = density_limit;

        GetTrackData(aTrack)->SetNuD(nud);
        GetTrackData(aTrack)->SetElD(eld);

        GetInfo(aTrack)->SetNucleiDensity(nud);
        GetInfo(aTrack)->SetElectronDensity(eld);

        GetTrackData(aTrack)->SetEFX(efx);
        GetTrackData(aTrack)->SetEFY(efy);
        
        //----------------------------------------
        // Scattering on electrons
        //----------------------------------------
        
        G4double gamma = aTrack.GetTotalEnergy()/mass;
        G4double m_e = CLHEP::electron_mass_c2;

        G4double T = aTrack.GetStep()->GetTotalEnergyDeposit();
        T -= aTrack.GetStep()->GetNonIonizingEnergyDeposit();
	    // default githube nera . 
        if(T>0.){
            G4double Tmax = 2. * m_e * gamma * gamma * beta * beta;
            Tmax /= (1. + 2.*gamma*m_e/mass + m_e*m_e/mass/mass);
            G4double mmod = aTrack.GetStep()->GetPreStepPoint()->GetMomentum().mag();
            G4double theta_sc_el = sqrt(2.*m_e*T*(1.-T/Tmax))/mmod;
            if(theta_sc_el>0.){
                G4double rot_angle = G4UniformRand()*CLHEP::twopi;
                G4double rot_mod   = G4RandGauss::shoot(0.,theta_sc_el*eld);
                mom.rotateX(cos(rot_angle)*rot_mod);
                mom.rotateY(sin(rot_angle)*rot_mod);
            }
        }

        GetTrackData(aTrack)->SetMomCh(mom);
        GetTrackData(aTrack)->SetPosCh(pos);

    	GetInfo(aTrack)->SetMomentumChanneled(mom);	
    	GetInfo(aTrack)->SetPositionChanneled(pos);

        GetInfo(aTrack)->SetLastChannelingPosition(pos);
        GetInfo(aTrack)->SetLastChannelingMomentum(mom);

        GetInfo(aTrack)->SetLastChannelingWorldPositionPost(posPost);
        GetInfo(aTrack)->SetLastChannelingWorldPositionPre(posPre);
        GetInfo(aTrack)->SetLastChannelingWorldMomentum(momWorld);

        return true;
    } else
        return false;
    
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4Channeling::
UpdateIntegrationStep(G4ChannelingMaterialData* matData,
                      const G4Track& aTrack,
                      G4ThreeVector& mom,
                      G4double& step){
    

    if(mom.x() != 0.0 || mom.y() != 0.0){
        double xy2 = mom.x() * mom.x() + mom.y()*mom.y();
        //G4cout << " use step size = " << use_step_size << " step size = " << step_size_value/nm << G4endl;
        if(use_step_size)
            step = step_size_value;
        else {
            if(xy2!=0.){
                step = std::fabs(fTransverseVariationMax * GetPre(aTrack)->GetKineticEnergy() / std::pow(xy2,0.5));
                if(step < fTimeStepMin) step = fTimeStepMin;
                else{
                    fTimeStepMax = std::sqrt( fTransverseVariationMax * GetPre(aTrack)->GetKineticEnergy()
                                        / std::fabs(matData->GetEFX()->GetMax()));
                    
                    if(step > fTimeStepMax) step = fTimeStepMax;
                }
            }
            else{
                step = fTimeStepMin;
            }
        }
        
        return true;
    }
    else{
        if(use_step_size)
            step = step_size_value;
        else
            step = fTimeStepMin;
    }
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4Channeling::
GetMeanFreePath(const G4Track& aTrack,
                G4double, // previousStepSize
                G4ForceCondition* condition){
    
    //----------------------------------------
    // the condition is forced to check if
    // the volume has a lattice at each step.
    // if it hasn't, return DBL_MAX
    //----------------------------------------
    
    *condition = Forced;

    G4LogicalVolume* aLV = aTrack.GetVolume()->GetLogicalVolume(); 
    G4LogicalVolume* aNLV = aTrack.GetNextVolume()->GetLogicalVolume();  
      
    if(G4LogicalCrystalVolume::IsLattice(aLV) == true &&
       G4LogicalCrystalVolume::IsLattice(aNLV) == true &&
       aTrack.GetKineticEnergy() > fMinimumEnergy){
        GetInfo(aTrack)->SetInTheCrystal(true);
        G4double osc_per = GetOscillationPeriod(aTrack);
        G4double osc_per2 = GetOscillationPeriod2(aTrack);

        if(use_step_size) {
            fTimeStepMin = step_size_value;
            return mfp_value;
        }

        if(org_step_size) {
            fTimeStepMin = osc_per * 2.e-4;
            return osc_per * 0.01;
        } else {
            fTimeStepMin = osc_per2*chan_step*1e-4;
            return osc_per2*chan_step;
        }
    }
    else{
        GetTrackData(aTrack)->Reset();
	    GetInfo(aTrack)->Reset();
        GetInfo(aTrack)->SetInTheCrystal(false);
        return DBL_MAX;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4Channeling::
PostStepDoIt(const G4Track& aTrack,
             const G4Step&){
    
    //----------------------------------------
    // check if the volume has a lattice
    // and if the particle is in channeling.
    // If it is so, the particle is forced
    // to follow the channeling plane
    // direction. If the particle has
    // dechanneled or exited the crystal,
    // the outgoing angle is evaluated
    //----------------------------------------
    
    aParticleChange.Initialize(aTrack);

    GetInfo(aTrack)->StoreDensityPreviousStep();

    G4LogicalVolume* aLV = aTrack.GetVolume()->GetLogicalVolume();
    G4LogicalVolume* aNLV = aTrack.GetNextVolume()->GetLogicalVolume();

    
    if(G4LogicalCrystalVolume::IsLattice(aLV) == true &&
       G4LogicalCrystalVolume::IsLattice(aNLV) == true){

        G4bool bModifiedTraj = UpdateParameters(aTrack);

        if(bModifiedTraj==true){
            //----------------------------------------
            // Get the momentum in the reference frame
            // solidal to the bent planes and rotate
            // to the reference frame
            //----------------------------------------
            G4LogicalCrystalVolume* aLCV = (G4LogicalCrystalVolume*)(aTrack.GetVolume()->GetLogicalVolume());
            G4ThreeVector momCh = GetTrackData(aTrack)->GetMomCh();
            //G4ThreeVector momCh = GetInfo(aTrack)->GetMomentumChanneled();
	    

            G4StepPoint* postStepPoint = aTrack.GetStep()->GetPostStepPoint();
            G4TouchableHistory* theTouchable = (G4TouchableHistory*)(postStepPoint->GetTouchable());
            

            //----------------------------------------
            // Get the momentum in the crystal reference
            // frame and rotate to the solid reference frame
            //----------------------------------------
            aLCV->RotateToSolid(momCh);
            //----------------------------------------
            // Get the momentum in the solid reference
            // frame and rotate to the world reference frame
            //----------------------------------------
            G4ThreeVector mom = ((*theTouchable->GetRotation()).inverse())(momCh);

            aParticleChange.ProposeMomentumDirection(mom.unit());
            aParticleChange.ProposePolarization(GetTrackData(aTrack)->GetSpinCh());
        }
    }
    else{
        // if the volume has no lattice it resets the density factors
        GetTrackData(aTrack)->Reset();
        GetInfo(aTrack)->Reset();
    }
    
    return &aParticleChange;
}
ExExChParticleUserInfo* G4Channeling::
GetInfo(const G4Track& aTrack){
    return (ExExChParticleUserInfo*) aTrack.GetUserInformation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
