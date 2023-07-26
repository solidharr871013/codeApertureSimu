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
/// \file B3PhysicsList.cc
/// \brief Implementation of the B3PhysicsList class


/*********************
#include "codeAperturePhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//#include "G4EmLivermorePhysics.hh"
//#include "G4LivermorecodeApertureModel.hh"
//#include "G4LivermorePhotoElectricModel.hh"
#include "G4PenelopecodeApertureModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4PenelopeRayleighModel.hh"
//#include "G4KleinNishinaModel.hh"
#include "G4codeApertureScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4StepLimiter.hh"

#include "G4SystemOfUnits.hh"

//#include "G4DecayPhysics.hh"
//#include "G4EmLivermorePhysics.hh"
//#include "G4EmStandardPhysics_option4.hh"

//#include "G4EmPenelopePhysics.hh"
//#include "G4RadioactiveDecayPhysics.hh"
//#include "G4OpticalPhysics.hh"
***************************/


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*************************
codeAperturePhysicsList::codeAperturePhysicsList()
: G4VUserPhysicsList(){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

codeAperturePhysicsList::~codeAperturePhysicsList()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void codeAperturePhysicsList::ConstructParticle(){
    ConstructBosons();
    ConstructLeptons();

    G4GenericIon::GenericIonDefinition();
}

void codeAperturePhysicsList::ConstructBosons(){
    G4Geantino::GeantinoDefinition();
    G4ChargedGeantino::ChargedGeantinoDefinition();

    G4Gamma::GammaDefinition();
}

void codeAperturePhysicsList::ConstructLeptons(){
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();

    G4NeutrinoE::NeutrinoEDefinition();
    G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void codeAperturePhysicsList::ConstructProcess(){
    AddTransportation();
    ConstructEM();
    ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void codeAperturePhysicsList::ConstructEM(){
    auto particleIterator = GetParticleIterator();
    particleIterator->reset();

    while ((*particleIterator)()) {
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        if(particleName == "gamma"){

            G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
            //theComptonScattering->SetEmModel(new G4LivermoreComptonModel(),1);
            theComptonScattering->SetEmModel(new G4PenelopeComptonModel(),1);

            G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
            thePhotoElectricEffect->SetEmModel(new G4PenelopePhotoElectricModel(),1);

            G4GammaConversion* theGammaConversion = new G4GammaConversion();
            theGammaConversion->SetEmModel(new G4PenelopeGammaConversionModel(),1);

            G4RayleighScattering* theRayleighScattering = new G4RayleighScattering();
            theRayleighScattering->SetEmModel(new G4PenelopeRayleighModel(),1);
            //G4PhotoElectricEffect* thePhotoEffect = new G4PhotoElectricEffect();
            //thePhotoEffect->SetEmModel(new G4LivermoreComptonModel(),1);

            pmanager->AddDiscreteProcess(thePhotoElectricEffect);
            //pmanager->AddDiscreteProcess(thePhotoEffect);
            pmanager->AddDiscreteProcess(theComptonScattering);
            pmanager->AddDiscreteProcess(theGammaConversion);
            pmanager->AddDiscreteProcess(theRayleighScattering);

        }
        else if(particleName == "e-"){
            pmanager->AddProcess(new G4eMultipleScattering, -1,1,1);
            pmanager->AddProcess(new G4eIonisation, -1,2,2);
            pmanager->AddProcess(new G4eBremsstrahlung, -1,3,3);
            pmanager->AddProcess(new G4StepLimiter, -1,-1,4);
        }
        else if (particleName == "e+") {
            pmanager->AddProcess(new G4eMultipleScattering, -1,1,1);
            pmanager->AddProcess(new G4eIonisation, -1,2,2);
            pmanager->AddProcess(new G4eBremsstrahlung,-1,3,3);
            pmanager->AddProcess(new G4eplusAnnihilation, 0,-1,4);
            pmanager->AddProcess(new G4StepLimiter, -1,-1,5);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4Decay.hh"
void codeAperturePhysicsList::ConstructGeneral(){
    G4Decay* theDecayProcess = new G4Decay();
    auto particleIterator=GetParticleIterator();
    particleIterator->reset();
    while ((*particleIterator)()){
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        if (theDecayProcess->IsApplicable(*particle)) {
          pmanager ->AddProcess(theDecayProcess);
          // set ordering for PostStepDoIt and AtRestDoIt
          pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
          pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
        }
    }
}
*****************************/
#include "codeAperturePhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

codeAperturePhysicsList::codeAperturePhysicsList()
: G4VModularPhysicsList(){
  SetVerboseLevel(1);

  // Default physics
  RegisterPhysics(new G4DecayPhysics());

  // EM physics
  RegisterPhysics(new G4EmPenelopePhysics());
  //RegisterPhysics(new G4EmLivermorePhysics());

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

codeAperturePhysicsList::~codeAperturePhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void codeAperturePhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCuts();
}
