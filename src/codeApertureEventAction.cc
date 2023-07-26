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
/// \file codeApertureEventAction.cc
/// \brief Implementation of the codeApertureEventAction class

#include "codeApertureEventAction.hh"
#include "codeApertureRunAction.hh"
#include "codeApertureSteppingAction.hh"
#include "codeApertureAnalysis.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4THitsMap.hh"

#include "codeApertureLayer1Hit.hh"
#include "codeApertureLayer1SD.hh"

#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"

#include "random"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {

// Utility function which finds a hit collection with the given Id
// and print warnings if not found
G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
  auto hce = event->GetHCofThisEvent();
  if (!hce) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl;
      G4Exception("codeApertureEventAction::EndOfEventAction()",
                  "codeApertureCode001", JustWarning, msg);
      return nullptr;
  }

  auto hc = hce->GetHC(collId);
  if ( ! hc) {
    G4ExceptionDescription msg;
    msg << "Hits collection " << collId << " of this event not found." << G4endl;
    G4Exception("codeApertureEventAction::EndOfEventAction()",
                "codeApertureCode001", JustWarning, msg);
  }
  return hc;
}
}

codeApertureEventAction::codeApertureEventAction(codeApertureRunAction* runAction)
: G4UserEventAction(),
  fMessenger(nullptr),
  fRunAction(runAction),
  fEdep(0.),
  f3LayerCoincidence(0),
  fTotalNumber(0),
  fOrdered3(0),
  fOrdered2(0),
  f2LayerCoincidence(0),
  fLayer1ID(-1),
  fLayer2ID(-1),
  fLayer3ID(-1),
  f2TotalE_low(0.6330),
  f2TotalE_high(0.693),
  f2ScatterE_low(0.020),
  f2ScatterE_high(0.16)
{
    DefineCommands();
     G4RunManager::GetRunManager()->SetPrintProgress(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

codeApertureEventAction::~codeApertureEventAction()
{ delete fMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void codeApertureEventAction::BeginOfEventAction(const G4Event*)
{    
//   fEdep = 0.;
//   f3LayerCoincidence = 0;
//   fTotalNumber = 0;
//   fOrdered3 = 0;
//   fOrdered2 = 0;
//   f2LayerCoincidence = 0;

//   if(fLayer1ID == -1 ){
//       auto sdManager = G4SDManager::GetSDMpointer();

//       fLayer1ID = sdManager->GetCollectionID("Layer1/LayerColl");
//       fLayer2ID = sdManager->GetCollectionID("Layer2/LayerColl");
//       fLayer3ID = sdManager->GetCollectionID("Layer3/LayerColl");
//   }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void codeApertureEventAction::EndOfEventAction(const G4Event* event)
{   
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();


    // auto hcLayer1 = GetHC(event, fLayer1ID);
    // if(!hcLayer1) return;
    // auto hcLayer2 = GetHC(event, fLayer2ID);
    // if(!hcLayer2) return;
    // auto hcLayer3 = GetHC(event, fLayer3ID);
    // if(!hcLayer3) return;



  fRunAction->AddTotalNumber(fTotalNumber);
  //G4cout << "the total number is " << totalNumber << G4endl;

    auto analysisManager = G4AnalysisManager::Instance();

    auto HCE = event->GetHCofThisEvent();
    if(!HCE) return;

    if(fLayer1ID<0){
        auto SDMan = G4SDManager::GetSDMpointer();
        fLayer1ID = SDMan->GetCollectionID("detectorEnergy/energyDeposition");
    }

    G4double eThreshold = 300*eV;

    G4THitsMap<G4double>* evtMap =
                       static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fLayer1ID));

    std::map<G4int,G4double*>::iterator itr;

    for(itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); ++itr){

        G4double CopyNo = static_cast<G4double>(itr->first);
        G4double eDep = *(itr->second);

        fEdep += eDep;

        if(eDep>100*keV && eDep<500*keV){
            analysisManager->FillH1(1,CopyNo);
           // G4cout << "CopyNo is " << CopyNo << G4endl;
        }
    }


  if(fEdep > 0.02*MeV){

      // analysisManager->FillH1(0,AddFluction(fEdep));
      analysisManager->FillH1(0,fEdep);
  }


}

G4double codeApertureEventAction::AddFluction(G4double val){

    //calculate the coefficient of the resolution, based on 662keV
    G4double Resolution662 =0.08,
             coefficient = Resolution662*std::sqrt(0.662);

    G4double valVar = coefficient*(std::sqrt(val))/(2*std::sqrt(2*std::log(2)));
    //G4double valVar = 0.08*val/(2*std::sqrt(2*std::log(2)));


    std::random_device rd;
    std::mt19937 generator(rd());

    std::normal_distribution<G4double> distribution(val,valVar);
    G4double var_new = distribution(generator);

    return var_new;
}


void codeApertureEventAction::Set2TotalELow(G4double val){f2TotalE_low = val;}
void codeApertureEventAction::Set2TotalEHigh(G4double val){f2TotalE_high = val;}
void codeApertureEventAction::Set2ScatterELow(G4double val){f2ScatterE_low = val;}
void codeApertureEventAction::Set2ScatterEHigh(G4double val){f2ScatterE_high = val;}

void codeApertureEventAction::DefineCommands(){
    fMessenger = new G4GenericMessenger(this,
                                        "/codeApertureCamera/event/",
                                        "Event Control");


    auto& ThresTotalEHighCmd
      = fMessenger->DeclareMethodWithUnit("ThresTotalEHigh","MeV",
                                  &codeApertureEventAction::Set2TotalEHigh,
                                  "Set high threshold of total edep in a 2coincident event.");
    ThresTotalEHighCmd.SetParameterName("Energy", true);
    ThresTotalEHighCmd.SetRange("Energy>=0. && Energy<100.");
    ThresTotalEHighCmd.SetDefaultValue("0.710");

    auto& ThresTotalELowCmd
      = fMessenger->DeclareMethodWithUnit("ThresTotalELow","MeV",
                                  &codeApertureEventAction::Set2TotalELow,
                                  "Set low threshold of total edep in a 2coincident event.");
    ThresTotalELowCmd.SetParameterName("Energy", true);
    ThresTotalELowCmd.SetRange("Energy>=0. && Energy<100.");
    ThresTotalELowCmd.SetDefaultValue("0.615");

    auto& ThresScatterEHighCmd
      = fMessenger->DeclareMethodWithUnit("ThresScatterEHigh","MeV",
                                  &codeApertureEventAction::Set2ScatterEHigh,
                                  "Set high threshold of Scatter layer edep in a 2coincident event.");
    ThresScatterEHighCmd.SetParameterName("Energy", true);
    ThresScatterEHighCmd.SetRange("Energy>=0. && Energy<100.");
    ThresScatterEHighCmd.SetDefaultValue("0.165");

    auto& ThresScatterELowCmd
      = fMessenger->DeclareMethodWithUnit("ThresScatterELow","MeV",
                                  &codeApertureEventAction::Set2ScatterELow,
                                  "Set Low threshold of Scatter layer edep in a 2coincident event.");
    ThresScatterELowCmd.SetParameterName("Energy", true);
    ThresScatterELowCmd.SetRange("Energy>=0. && Energy<100.");
    ThresScatterELowCmd.SetDefaultValue("0.01");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
