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
/// \file codeAperturePrimaryGeneratorAction.cc
/// \brief Implementation of the codeAperturePrimaryGeneratorAction class

#include "codeAperturePrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "random"
#include "G4GenericMessenger.hh"
#include "cmath"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

codeAperturePrimaryGeneratorAction::codeAperturePrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fMessenger(nullptr),
  fParticleGun(nullptr),
  xpos(0.),
  ypos(0.),
  zpos(-40.)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  DefineCommands();

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(500*keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

codeAperturePrimaryGeneratorAction::~codeAperturePrimaryGeneratorAction()
{
    delete  fMessenger;
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void codeAperturePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.






    //the following code is for generating two seperated source, the source is aparted in
    // the distance of source_divergence in x-axis.


    //G4int randNumber = rand()%2;
    G4double randNumber = G4UniformRand();


    // G4double cellNo = 16, cellSide = 2.*mm;

    // G4double detectDim = 0.5*cellNo*cellSide;

    // G4double targetX = detectDim*(4*G4UniformRand()-2),
    //                     targetY = detectDim*(4*G4UniformRand()-2),
    //                     targetZ = -15*mm;
    // if(randNumber<0.5){

    //     G4double source_x= 10*cm;
    //     G4double source_y= 0*cm;
    //     G4double source_z= -200*cm;

    //     G4double x_direction = targetX-source_x,
    //              y_direction = targetY-source_y,
    //              z_direction = targetZ-source_z;

    //     fParticleGun->SetParticleEnergy(500*keV);

    //     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x_direction,y_direction,z_direction));

    //     fParticleGun->SetParticlePosition(G4ThreeVector(source_x,source_y,source_z));

    //     fParticleGun->GeneratePrimaryVertex(anEvent);

    // }

    // else if(randNumber>=0.5){

    //     G4double source_x=-10*cm;
    //     G4double source_y=0*cm;
    //     G4double source_z=-200*cm;

    //     G4double x_direction = targetX-source_x,
    //              y_direction = targetY-source_y,
    //              z_direction = targetZ-source_z;

    //     fParticleGun->SetParticleEnergy(500*keV);

    //     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x_direction,y_direction,z_direction));

    //     fParticleGun->SetParticlePosition(G4ThreeVector(source_x,source_y,source_z));

    //     fParticleGun->GeneratePrimaryVertex(anEvent);
    // }


/*************multi-source end**********************/

/********* single source ************/


   G4double source_x = 0*cm, 
            source_y = 0*cm, 
            source_z = -200*cm;
  
   G4double cos_phi = 1-0.03*randNumber, sin_phi = std::sqrt(1-cos_phi*cos_phi);

   G4double theta = G4UniformRand()*2*pi,
            cos_theta = std::cos(theta), sin_theta = std::sin(theta);

    G4double x_direction = sin_phi*cos_theta,
             y_direction = sin_phi*sin_theta,
             z_direction = cos_phi;
//    G4double x_direction = targetX-source_x,
//             y_direction = targetY-source_y,
//             z_direction = targetZ-source_z;

    fParticleGun->SetParticleEnergy(511*keV);

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x_direction,y_direction,z_direction));

    fParticleGun->SetParticlePosition(G4ThreeVector(source_x,source_y,source_z));

    fParticleGun->GeneratePrimaryVertex(anEvent);


/**********single source end******************/

}

void codeAperturePrimaryGeneratorAction::DefineCommands(){
    fMessenger = new G4GenericMessenger(this,
                                        "/codeAperture/gunpos/",
                                        "gun Control");


    auto& xposCmd
      = fMessenger->DeclareMethodWithUnit("xpos","cm",
                                  &codeAperturePrimaryGeneratorAction::setPosX,
                                  "Set x position for source.");
    xposCmd.SetParameterName("xPosition", true);
    xposCmd.SetRange("xPosition>=-10000. && xPosition<10000.");
    xposCmd.SetDefaultValue("0");

    auto& yposCmd
      = fMessenger->DeclareMethodWithUnit("ypos","cm",
                                  &codeAperturePrimaryGeneratorAction::setPosY,
                                  "Set y position for source.");
    yposCmd.SetParameterName("yPosition", true);
    yposCmd.SetRange("yPosition>=-10000. && yPosition<10000.");
    yposCmd.SetDefaultValue("40");

    auto& zposCmd
      = fMessenger->DeclareMethodWithUnit("zpos","cm",
                                  &codeAperturePrimaryGeneratorAction::setPosZ,
                                  "Set z position for source.");
    zposCmd.SetParameterName("zPosition", true);
    zposCmd.SetRange("zPosition>=-10000. && zPosition<10000.");
    zposCmd.SetDefaultValue("-40.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

