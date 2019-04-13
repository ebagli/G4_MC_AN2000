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

#include "EventAction.hh"

#include "G4RunManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "CrystalDetectorHit.hh"
#include "CrystalDetectorHit.hh"
#include "Analysis.hh"

#include "G4GenericMessenger.hh"

EventAction::EventAction():
sdht_ID(-1),
sdct_ID(-1),
selectedReaction(0){
    fVectorEC = new G4PhysicsFreeVector(23);
    fVectorEC->PutValue(0,0.*CLHEP::keV,0.);
    fVectorEC->PutValue(1,258.75*CLHEP::keV,13.85);
    fVectorEC->PutValue(2,517.5*CLHEP::keV,27.7);
    fVectorEC->PutValue(3,542.5*CLHEP::keV,40.46);
    fVectorEC->PutValue(4,572.5*CLHEP::keV,66.98);
    fVectorEC->PutValue(5,592.5*CLHEP::keV,92.22);
    fVectorEC->PutValue(6,612.5*CLHEP::keV,136.58);
    fVectorEC->PutValue(7,622.5*CLHEP::keV,209.26);
    fVectorEC->PutValue(8,623.5*CLHEP::keV,239.1);
    fVectorEC->PutValue(9,624.5*CLHEP::keV,248.72);
    fVectorEC->PutValue(10,625.5*CLHEP::keV,300.18);
    fVectorEC->PutValue(11,626.5*CLHEP::keV,374.9);
    fVectorEC->PutValue(12,627.5*CLHEP::keV,920.2);
    fVectorEC->PutValue(13,628.5*CLHEP::keV,996.64);
    fVectorEC->PutValue(14,629.5*CLHEP::keV,608.64);
    fVectorEC->PutValue(15,630.5*CLHEP::keV,195.4);
    fVectorEC->PutValue(16,631.5*CLHEP::keV,109.94);
    fVectorEC->PutValue(17,632.5*CLHEP::keV,104.9);
    fVectorEC->PutValue(18,637.5*CLHEP::keV,117.74);
    fVectorEC->PutValue(19,642.5*CLHEP::keV,126.18);
    fVectorEC->PutValue(20,647.5*CLHEP::keV,142.34);
    fVectorEC->PutValue(21,652.5*CLHEP::keV,155.330324386);
    fVectorEC->PutValue(22,652.5*CLHEP::keV,150.94);
    
    
    fVectorEC_Y89 = new G4PhysicsFreeVector(45);
    fVectorEC_Y89->PutValue(0,1.00000E+00 *CLHEP::MeV,0.00000E+00);
    fVectorEC_Y89->PutValue(1,2.00000E+00 *CLHEP::MeV,0.00000E+00);
    fVectorEC_Y89->PutValue(2,3.00000E+00 *CLHEP::MeV,0.00000E+00);
    fVectorEC_Y89->PutValue(3,4.00000E+00 *CLHEP::MeV,2.27052E+00);
    fVectorEC_Y89->PutValue(4,5.00000E+00 *CLHEP::MeV,5.89439E+01);
    fVectorEC_Y89->PutValue(5,6.00000E+00 *CLHEP::MeV,1.79371E+02);
    fVectorEC_Y89->PutValue(6,7.00000E+00 *CLHEP::MeV,3.43627E+02);
    fVectorEC_Y89->PutValue(7,8.00000E+00 *CLHEP::MeV,5.02739E+02);
    fVectorEC_Y89->PutValue(8,9.00000E+00 *CLHEP::MeV,6.21005E+02);
    fVectorEC_Y89->PutValue(9,1.00000E+01 *CLHEP::MeV,7.02344E+02);
    fVectorEC_Y89->PutValue(10,1.10000E+01 *CLHEP::MeV,7.62526E+02);
    fVectorEC_Y89->PutValue(11,1.20000E+01 *CLHEP::MeV,8.14053E+02);
    fVectorEC_Y89->PutValue(12,1.30000E+01 *CLHEP::MeV,8.59700E+02);
    fVectorEC_Y89->PutValue(13,1.40000E+01 *CLHEP::MeV,7.95993E+02);
    fVectorEC_Y89->PutValue(14,1.50000E+01 *CLHEP::MeV,6.46547E+02);
    fVectorEC_Y89->PutValue(15,1.60000E+01 *CLHEP::MeV,4.70527E+02);
    fVectorEC_Y89->PutValue(16,1.70000E+01 *CLHEP::MeV,3.12452E+02);
    fVectorEC_Y89->PutValue(17,1.80000E+01 *CLHEP::MeV,2.09222E+02);
    fVectorEC_Y89->PutValue(18,1.90000E+01 *CLHEP::MeV,1.47963E+02);
    fVectorEC_Y89->PutValue(19,2.00000E+01 *CLHEP::MeV,1.14786E+02);
    fVectorEC_Y89->PutValue(20,2.20000E+01 *CLHEP::MeV,8.47757E+01);
    fVectorEC_Y89->PutValue(21,2.40000E+01 *CLHEP::MeV,7.22970E+01);
    fVectorEC_Y89->PutValue(22,2.60000E+01 *CLHEP::MeV,6.52142E+01);
    fVectorEC_Y89->PutValue(23,2.80000E+01 *CLHEP::MeV,5.96815E+01);
    fVectorEC_Y89->PutValue(24,3.00000E+01 *CLHEP::MeV,5.50719E+01);
    fVectorEC_Y89->PutValue(25,3.50000E+01 *CLHEP::MeV,4.49117E+01);
    fVectorEC_Y89->PutValue(26,4.00000E+01 *CLHEP::MeV,3.68834E+01);
    fVectorEC_Y89->PutValue(27,4.50000E+01 *CLHEP::MeV,3.05675E+01);
    fVectorEC_Y89->PutValue(28,5.00000E+01 *CLHEP::MeV,2.59446E+01);
    fVectorEC_Y89->PutValue(29,5.50000E+01 *CLHEP::MeV,2.18940E+01);
    fVectorEC_Y89->PutValue(30,6.00000E+01 *CLHEP::MeV,1.89228E+01);
    fVectorEC_Y89->PutValue(31,6.50000E+01 *CLHEP::MeV,1.60827E+01);
    fVectorEC_Y89->PutValue(32,7.00000E+01 *CLHEP::MeV,1.42660E+01);
    fVectorEC_Y89->PutValue(33,7.50000E+01 *CLHEP::MeV,1.24970E+01);
    fVectorEC_Y89->PutValue(34,8.00000E+01 *CLHEP::MeV,1.12482E+01);
    fVectorEC_Y89->PutValue(35,9.00000E+01 *CLHEP::MeV,8.83529E+00);
    fVectorEC_Y89->PutValue(36,1.00000E+02 *CLHEP::MeV,7.10967E+00);
    fVectorEC_Y89->PutValue(37,1.10000E+02 *CLHEP::MeV,5.82436E+00);
    fVectorEC_Y89->PutValue(38,1.20000E+02 *CLHEP::MeV,4.78898E+00);
    fVectorEC_Y89->PutValue(39,1.30000E+02 *CLHEP::MeV,4.07154E+00);
    fVectorEC_Y89->PutValue(40,1.40000E+02 *CLHEP::MeV,3.47351E+00);
    fVectorEC_Y89->PutValue(41,1.50000E+02 *CLHEP::MeV,3.05016E+00);
    fVectorEC_Y89->PutValue(42,1.60000E+02 *CLHEP::MeV,2.68101E+00);
    fVectorEC_Y89->PutValue(43,1.80000E+02 *CLHEP::MeV,2.20578E+00);
    fVectorEC_Y89->PutValue(44,2.00000E+02 *CLHEP::MeV,1.84821E+00);
    // -- Define messengers:
    fChangeReactionType =
    new G4GenericMessenger(this, "/reaction/","Change Reaction" );
    
    fChangeReactionType->DeclareProperty("set",
                                         selectedReaction,
                                         "Set reaction." );

    // -- Define messengers:
    fChangeStepSave =
    new G4GenericMessenger(this, "/reaction/","Change Reaction" );
    
    fChangeStepSave->DeclareProperty("step",
                                     selectedStep,
                                     "Set step." );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::BeginOfEventAction(const G4Event*){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event* evt){
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    G4double nudavg = 0.;
    G4double eldavg = 0.;
    G4double nudOavg = 0.;
    G4double steptot = 0.;
    G4double kinE = 0.;
    G4double kinEtotOxygen = 0.;
    G4double kinEtot = 0.;
    G4double xs = 0.;

    G4double nudAlavg = 0.;

    G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    //G4int eventID = G4EventManager::GetEventManager()->GetTrackingManager()->GetTrack()->GetTrackID();
    
    G4double nud;
    G4double eld;
    G4double nudO;
    G4double nudAl;

    G4double step;

    if(sdct_ID == -1) {
        G4String sdName;
        if(SDman->FindSensitiveDetector(sdName="crystaldetector",0)){
            sdct_ID = SDman->GetCollectionID(sdName="crystaldetector/collection");
        }
    }
    
    CrystalDetectorHitsCollection* sdct = 0;
    G4HCofThisEvent *hce1 = evt->GetHCofThisEvent();
    
    if(hce1){
        if(sdct_ID != -1){
            G4VHitsCollection* aHCSD1 = hce1->GetHC(sdct_ID);
            sdct = (CrystalDetectorHitsCollection*)(aHCSD1);
        }
    }
    
    // Select the reaction
    //0 -> O18(p,a)
    //1 -> P31(a,p)
    //2 -> Al(p,gamma)

    G4double stepSave    = selectedStep * CLHEP::micrometer;
    G4double steptotPre    = 0.;
    G4double stepTemp    = 0.;
    G4int    stepSaveBin = 0;
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4ThreeVector pos;
    if(sdct){
        int n_hit_sd = sdct->entries();
        G4double kinEprev = 0.;
        for(int i1=0;i1<n_hit_sd;i1++){
            CrystalDetectorHit* aHit = (*sdct)[i1];
            step = aHit->GetStep();
            nud = aHit->GetNuD();
            eld = aHit->GetElD();
            nudO = aHit->GetNuD2();
            kinE = aHit->GetKinE();
            if(kinEprev == 0.) kinEprev = kinE;
            nudAl = aHit->GetNuD3();
            pos = aHit->GetChPos();

            steptot += step;
            nudavg  += nud * step;
            eldavg  += eld * step;
            nudOavg  += nudO * step;
            kinEtot  += kinE * step;
            nudAlavg  += nudAl * step;

            xs = 0;
            
            if(selectedReaction==0){
                xs = fVectorEC->Value(kinE);
            }
            else if(selectedReaction==3){
                xs = fVectorEC_Y89->Value(kinE);
            }
            else if(selectedReaction==1){
                G4double gamma = 7. * CLHEP:: keV;
                G4double x0 = 5. * CLHEP:: MeV;
                G4double x = kinE;
                G4double cauchy = gamma*gamma;
                cauchy /= ( (x-x0)*(x-x0) + gamma*gamma );
                cauchy /= (CLHEP::pi * gamma);
                xs = cauchy;
            }
            /*
            else if(selectedReaction==2){
                G4double gamma = 0.07 * CLHEP:: keV;
                G4double x0 = 991.8 * CLHEP:: keV;
                G4double x = kinE;
                G4double cauchy = gamma*gamma;
                cauchy /= ( (x-x0)*(x-x0) + gamma*gamma );
                cauchy /= (CLHEP::pi * gamma);
                xs = cauchy;
            }
            else if(selectedReaction==2){
                G4double gamma = 0.07 * CLHEP:: keV;
                G4double x0 = 991.8 * CLHEP:: keV;
                G4double cauchy_1 = std::atan2((kinEprev-x0),gamma)/CLHEP::pi + 0.5;
                G4double cauchy_2 = std::atan2((kinE-x0),gamma)/CLHEP::pi + 0.5;
                xs = fabs(cauchy_2 - cauchy_1);
                
                //G4cout << kinEprev/CLHEP:: keV << " " << kinE/CLHEP:: keV << G4endl;
                //G4cout << cauchy_1 << " " << cauchy_2 << G4endl;
                //G4cout << xs << G4endl;
                
            }
            */
            else if(selectedReaction==2){
                G4double gamma = 0.07 * CLHEP:: keV;
                G4double x0 = 991.8 * CLHEP:: keV;
                G4double x1 = x0 - gamma * 0.5;
                G4double x2 = x0 + gamma * 0.5;
                G4double Emin = kinE;
                G4double Emax = kinEprev;
                if(Emin < x1 && Emax > x2){
                    //G4cout << "case 1" << G4endl;
                    xs = (Emax - Emin)/gamma;
                }
                else if(Emin < x2 && Emax > x2){
                    //G4cout << "case 2" << G4endl;
                    xs = (x2 - Emin)/(Emax - Emin);
                }
                else if(Emin < x1 && Emax > x1){
                    //G4cout << "case 3" << G4endl;
                    xs = (Emax - x1)/(Emax - Emin);
                }
                else{
                    //G4cout << "None" << G4endl;
                }
                //G4cout << kinEprev/CLHEP:: keV << " " << kinE/CLHEP:: keV << G4endl;
                //G4cout << x1/CLHEP:: keV << " " << x2/CLHEP:: keV << G4endl;

            }
            kinEtotOxygen  += nudO * step * xs;
            //G4cout << nudO << " " << step << " " << xs << " " << kinEtotOxygen << G4endl;
            //while(!getchar());
			/*
            if((steptot > stepSave * stepSaveBin) || (i1 == n_hit_sd - 1)){
                
                stepTemp = (steptot-steptotPre);
                steptotPre = steptot;
                stepSaveBin += 1;
                if( stepTemp != 0.0 ){
                    
                    analysisManager->FillNtupleDColumn(0, nudavg / stepTemp);
                    analysisManager->FillNtupleDColumn(1, eldavg / stepTemp);
                    analysisManager->FillNtupleDColumn(2, nudOavg / stepTemp);
                    analysisManager->FillNtupleDColumn(3, kinEtotOxygen / stepTemp);
                    analysisManager->FillNtupleDColumn(4, steptot / CLHEP::micrometer);
                    analysisManager->FillNtupleDColumn(5, stepTemp / CLHEP::micrometer);
                    analysisManager->FillNtupleDColumn(6, xs);
                    analysisManager->FillNtupleDColumn(7, kinEtot / stepTemp);
                    analysisManager->FillNtupleIColumn(8, eventID);
                    analysisManager->FillNtupleDColumn(9, nudAlavg / stepTemp);
                    analysisManager->FillNtupleDColumn(10, pos.x()/CLHEP::angstrom);
                    analysisManager->FillNtupleDColumn(11, pos.y()/CLHEP::angstrom);
                    analysisManager->AddNtupleRow();
                   
                    
                    nudavg = 0.;
                    eldavg = 0.;
                    nudOavg = 0.;
                    kinEtot = 0.;
                    kinEtotOxygen = 0.;
                    nudAlavg = 0.;

                }
				*/
             
            }
            
            kinEprev = kinE;
        }
        
        analysisManager->FillNtupleDColumn(0, nudavg / stepTemp);
        analysisManager->FillNtupleDColumn(1, eldavg / stepTemp);
        analysisManager->FillNtupleDColumn(2, nudOavg / stepTemp);
        analysisManager->FillNtupleDColumn(3, kinEtotOxygen / stepTemp);
        analysisManager->FillNtupleDColumn(4, steptot / CLHEP::micrometer);
        analysisManager->FillNtupleDColumn(5, (steptotPre + stepTemp * 0.5) / CLHEP::micrometer);
        analysisManager->FillNtupleDColumn(6, xs);
        analysisManager->FillNtupleDColumn(7, kinEtot / stepTemp);
        analysisManager->FillNtupleIColumn(8, eventID);
        analysisManager->AddNtupleRow();
        
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
