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
/// \file codeApertureDetectorConstruction.cc
/// \brief Implementation of the codeApertureDetectorConstruction class

#include "codeApertureDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"

#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4GenericMessenger.hh"

#include "codeApertureLayer1SD.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

codeApertureDetectorConstruction::codeApertureDetectorConstruction()
: G4VUserDetectorConstruction(),
  fMessenger(nullptr),
  fLogMask(nullptr),
  fLogDet(nullptr),
  fLogCell(nullptr)
{
    //DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

codeApertureDetectorConstruction::~codeApertureDetectorConstruction()
{
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* codeApertureDetectorConstruction::Construct()
{
   //
   // Construct the material to be used
   //
   ConstructMaterial();
   G4Material* world_mat = G4Material::GetMaterial("G4_AIR");
   G4Material* mask_mat = G4Material::GetMaterial("tungsten");
   G4Material* det_mat = G4Material::GetMaterial("G4_SODIUM_IODIDE");
   G4Material* detLayer_mat = G4Material::GetMaterial("BaSO4");

  //
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 10*m;
  G4double world_sizeZ  = 10*m;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  fLogicalWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(nullptr,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      fLogicalWorld,            //its logical volume
                      "World",               //its name
                      nullptr,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  G4RotationMatrix* rotation = new G4RotationMatrix();
  rotation->rotateY(0.*deg);

  // Detector layer define
  G4int cellNO = 19;
  G4double cellSide = 2*mm, thicknessDet = 20*mm, thicknessCell = thicknessDet, maskCellThick = 2*cm;
  G4double gap = 0.02*mm, det_xy = cellSide*cellNO + gap*(cellNO-1);

  G4ThreeVector detectorPos = G4ThreeVector(0,0,det_xy+0.5*(thicknessCell+maskCellThick));

  fSolidDet = new G4Box("soliddet", 0.5*det_xy, 0.5*det_xy, 0.5*thicknessDet);

  fLogDet = new G4LogicalVolume(fSolidDet, detLayer_mat, "detector");

  fPhysDet = new G4PVPlacement(rotation, detectorPos, fLogDet, "detector", fLogicalWorld, false, 0, checkOverlaps);

  // cell define

  fSolidCell = new G4Box("solidcell", 0.5*cellSide, 0.5*cellSide, 0.5*thicknessCell);

  fLogCell = new G4LogicalVolume(fSolidCell, det_mat, "cell");

  G4double adjust_displacement_x = (cellNO*cellSide+(cellNO+1)*gap)*0.5,
           adjust_displacement_y = (cellNO*cellSide+(cellNO+1)*gap)*0.5;
  G4int layer_NOCopy = 0;

  for (G4int i=0; i!=cellNO; ++i) {
    for (G4int j=0; j!=cellNO; ++j) {


        G4VPhysicalVolume* cellphys = new G4PVPlacement(nullptr,
                            G4ThreeVector((i+0.5)*cellSide+(i+1)*gap
                                           -adjust_displacement_x,
                                            ((j+0.5)*cellSide+(j+1)*gap
                                           -adjust_displacement_y),
                                            0),
                            fLogCell,
                            "cells",
                            fLogDet,
                            false,
                            layer_NOCopy,
                            true);
        fPhysCell.push_back(cellphys);
            ++layer_NOCopy;
        }

    }

  //Mask define

  G4int maskCellNo = 19;
  G4double maskCellSide = cellSide, maskLayerSide = maskCellSide*maskCellNo, maskLayerThick = maskCellThick;

  G4double maskMotherLayerSide = 2*maskLayerSide, maskMotherLayerThick = maskLayerThick;

  G4Box* maskMotherLayerSolid = new G4Box("maskMotherLayerSolid",0.5*maskMotherLayerSide,0.5*maskMotherLayerSide,0.5*maskMotherLayerThick);
  G4LogicalVolume* maskMotherLayerLog = new G4LogicalVolume(maskMotherLayerSolid,world_mat,"maskMotherLayerLog");

  G4Box* maskLayerSolid = new G4Box("maskLayerSolid",0.5*maskLayerSide,0.5*maskLayerSide,0.5*maskLayerThick);
  G4LogicalVolume* maskLayerLog = new G4LogicalVolume(maskLayerSolid,world_mat,"maskLayerLog");

  G4RotationMatrix* rotation_mask = new G4RotationMatrix();
  rotation_mask->rotateZ(90.*deg);

  for(G4int xitr=0; xitr!=2; ++xitr){
    for(G4int yitr=0; yitr!=2; ++yitr){

      G4ThreeVector maskLayerPos = G4ThreeVector((xitr-0.5)*maskLayerSide,(yitr-0.5)*maskLayerSide,0);

      new G4PVPlacement(
        rotation_mask,
        maskLayerPos,
        maskLayerLog,
        "maskLayerPhy",
        maskMotherLayerLog,
        false,
        0,
        checkOverlaps);

    }
  }

  new G4PVPlacement(
    nullptr,
    G4ThreeVector(0,0,0),
    maskMotherLayerLog,
    "maskMotherLayerPhy",
    fLogicalWorld,
    false,
    0,
    checkOverlaps
  );

  //place the coded unit to the mask

  G4Box* maskCellSolid = new G4Box("maskCellSolid",maskCellSide*0.5,maskCellSide*0.5,maskCellThick*0.5);
  G4LogicalVolume* maskCellLog = new G4LogicalVolume(maskCellSolid,mask_mat,"maskCellLog");

  int endCount = round((maskCellNo-1)/2);

  for(G4int xitr=0; xitr!=maskCellNo; ++xitr){
    for(G4int yitr=0; yitr!=maskCellNo; ++yitr){

      if(xitr!=0&&yitr==0){
        new G4PVPlacement(
          nullptr,
          G4ThreeVector((xitr-endCount)*maskCellSide,(yitr-endCount)*maskCellSide,0),
          maskCellLog,
          "maskCellPhy",
          maskLayerLog,
          false,
          0,
          checkOverlaps
        );
      }
    }
  }



 std::vector<int> QRnum, NonQRnum;

 for(int itr=1; itr!=endCount+1; ++itr){
     int label = (itr*itr)%maskCellNo;
     QRnum.push_back(label);
 }

 for(int itr=1; itr!=maskCellNo; ++itr){
     bool exists = std::find(std::begin(QRnum), std::end(QRnum), itr) != std::end(QRnum);
     if(exists==false){
         NonQRnum.push_back(itr);
     }
 }


  for(unsigned long xitr=0; xitr!=QRnum.size(); ++xitr){
    for(unsigned long yitr=0; yitr!=QRnum.size(); ++yitr){

        G4int xlabel = QRnum[xitr], ylabel = QRnum[yitr];

        new G4PVPlacement(
          nullptr,
          G4ThreeVector((xlabel-endCount)*maskCellSide,(ylabel-endCount)*maskCellSide,0),
          maskCellLog,
          "maskCellPhy",
          maskLayerLog,
          false,
          0,
          checkOverlaps
        );
           
    }
  }

  for(unsigned long xitr=0; xitr!=NonQRnum.size(); ++xitr){
    for(unsigned long yitr=0; yitr!=NonQRnum.size(); ++yitr){

        G4int xlabel = NonQRnum.at(xitr), ylabel = NonQRnum.at(yitr);

        new G4PVPlacement(
          nullptr,
          G4ThreeVector((xlabel-endCount)*maskCellSide,(ylabel-endCount)*maskCellSide,0),
          maskCellLog,
          "maskCellPhy",
          maskLayerLog,
          false,
          0,
          checkOverlaps
        );

    }
  }


  fLogDet->SetVisAttributes(G4VisAttributes::GetInvisible());


  //
  //always return the physical World
  //
  return physWorld;
}

void codeApertureDetectorConstruction::ConstructMaterial(){
    // Get NistManager
    G4NistManager* man = G4NistManager::Instance();

    man->FindOrBuildMaterial("G4_AIR");

    // GAGG material
    G4Element*  O = man->FindOrBuildElement( "O");
    G4Element*  S = man->FindOrBuildElement( "S");
    G4Element* Si = man->FindOrBuildElement("Si");
    G4Element* Lu = man->FindOrBuildElement("Lu");
    G4Element* Gd = man->FindOrBuildElement("Gd");
    G4Element* Ca = man->FindOrBuildElement("Ca");
    G4Element* Ba = man->FindOrBuildElement("Ba");
    G4Element* Al = man->FindOrBuildElement("Al");
    G4Element* Br = man->FindOrBuildElement("Br");
    G4Element* La = man->FindOrBuildElement("La");
    G4Element* Y = man->FindOrBuildElement("Y");
    G4Element* W = man->FindOrBuildElement("W");

    G4Material* LaBr = new G4Material("LaBr",5.06*g/cm3,2);
    LaBr->AddElement(La,1);
    LaBr->AddElement(Br,3);

    G4Material* LSO = new G4Material("Lu2SiO5",7.1*g/cm3,3);
    LSO->AddElement(Lu,2);
    LSO->AddElement(Si,1);
    LSO->AddElement(O,5);

    G4Material* YSO = new G4Material("Y2SiO5",4.5*g/cm3,3);
    YSO->AddElement(Y,2);
    YSO->AddElement(Si,1);
    YSO->AddElement(O,5);

    G4Material* GAGG = new G4Material("Gd3CaAl4O12", 6.63*g/cm3,4);
    GAGG->AddElement(Gd,3);
    GAGG->AddElement(Ca,1);
    GAGG->AddElement(Al,4);
    GAGG->AddElement(O,12);

    G4Material* BaSO4 = new G4Material("BaSO4", 4.49*g/cm3,3);
    BaSO4->AddElement(Ba,1);
    BaSO4->AddElement(S,1);
    BaSO4->AddElement(O,4);

    G4Material* Siboard = new G4Material("Siboard", 2.33*g/cm3,1);
    Siboard->AddElement(Si,1);

    G4Material* tungsten = new G4Material("tungsten", 19.35*g/cm3,1);
    tungsten->AddElement(W,1);

    man->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    man->FindOrBuildMaterial("G4_Pb");

}

void codeApertureDetectorConstruction::ConstructSDandField(){

    // auto sdManger = G4SDManager::GetSDMpointer();
    // G4String SDname;

    // auto Layer1 = new codeApertureLayer1SD(SDname = "/Layer1");
    // sdManger->AddNewDetector(Layer1);
    // fLogCell->SetSensitiveDetector(Layer1);

    auto sdManger = G4SDManager::GetSDMpointer();
    sdManger->SetVerboseLevel(1);
    G4String SDname;

    auto Layer1 = new codeApertureLayer1SD(SDname = "/Layer1");
    sdManger->AddNewDetector(Layer1);
    fLogCell->SetSensitiveDetector(Layer1);

    auto detectorEnergy = new G4MultiFunctionalDetector("detectorEnergy");
    sdManger->AddNewDetector(detectorEnergy);
    G4VPrimitiveScorer* primitive = new G4PSEnergyDeposit("energyDeposition");
    detectorEnergy->RegisterPrimitive(primitive);
    SetSensitiveDetector(fLogCell,detectorEnergy);



}

// void codeApertureDetectorConstruction::SetMaskThickness(G4double val){
//     fThicknessMask = val;
//     fSolidDet->SetZHalfLength(0.5*fThicknessDet);
//     fSolidCell->SetZHalfLength(0.5*fThicknessCell);

//     G4RunManager::GetRunManager()->ReinitializeGeometry();
//     // tell G4RunManager that we change the geometry
//     G4RunManager::GetRunManager()->GeometryHasBeenModified();
// // }



// /************** Setting the Layer position ***************/
// void codeApertureDetectorConstruction::SetDetZPosition(G4double val){

//     fPhysDet->SetTranslation(G4ThreeVector(0,0,val));

//     G4RunManager::GetRunManager()->GeometryHasBeenModified();
// }


// /****************** Setting the Cell size **********/
// void codeApertureDetectorConstruction::SetCellSize(G4double val){
//     fCellSide = val;
//     G4int cellNO = 16;
//     G4double gap = 0.*mm;

//     //change the cell size
//     fSolidCell->SetXHalfLength(0.5*fCellSide);
//     fSolidCell->SetYHalfLength(0.5*fCellSide);

//     //new layer size
//     fSolidDet->SetXHalfLength(0.5*((cellNO+1)*gap+cellNO*fCellSide));
//     fSolidDet->SetYHalfLength(0.5*((cellNO+1)*gap+cellNO*fCellSide));

//     //replace the cells, fisrt place back, then place them to the new places
//     G4double adjust_displacement_x = (cellNO*fCellSide+(cellNO+1)*gap)*0.5,
//              adjust_displacement_y = (cellNO*fCellSide+(cellNO+1)*gap)*0.5;
//     G4int cellItr = 255, rowNo = 0, columnNo = 0;

//     for(G4int i=0; i!=cellItr+1; ++i){

//         rowNo = i%cellNO;
//         columnNo = i/cellNO;

//         G4ThreeVector newTrans = G4ThreeVector((rowNo+0.5)*fCellSide+(rowNo+1)*gap
//                                                -adjust_displacement_x,
//                                                 ((columnNo+0.5)*fCellSide+(columnNo+1)*gap
//                                                -adjust_displacement_y),
//                                                 0);

//         fPhysCell[i]->SetTranslation(newTrans);

//     }
//     for(G4int i=0; i!=cellItr+1; ++i){
//         fPhysCell[i]->CheckOverlaps();
//     }

//     G4RunManager::GetRunManager()->GeometryHasBeenModified();
// }



// /******************** setting the physical volume of the cells ************/

// void codeApertureDetectorConstruction::DefineCommands(){

//     fMessenger = new G4GenericMessenger(this,
//                                         "/codeApertureCamera/detector/",
//                                         "Detector control");

//     //  command setting the thickness of the Mask
//     auto& Layer1ThicknessCmd
//       = fMessenger->DeclareMethodWithUnit("MaskThickness","mm",
//                                   &codeApertureDetectorConstruction::SetMaskThickness,
//                                   "Set thickness of the Mask.");
//     Layer1ThicknessCmd.SetParameterName("Thickness", true);
//     Layer1ThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
//     Layer1ThicknessCmd.SetDefaultValue("8.");

//     // command setting the thickness of the Detector
//     auto& Layer2ThicknessCmd
//       = fMessenger->DeclareMethodWithUnit("DetThickness","mm",
//                                   &codeApertureDetectorConstruction::SetDetThickness,
//                                   "Set thickness of the detector.");
//     Layer2ThicknessCmd.SetParameterName("Thickness", true);
//     Layer2ThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
//     Layer2ThicknessCmd.SetDefaultValue("8.");

//     // command setting the side width of the cell of the detector

//     auto& Cell1SizeCmd
//       = fMessenger->DeclareMethodWithUnit("CellSize","mm",
//                                   &codeApertureDetectorConstruction::SetCellSize,
//                                   "Set size of the Cell.");
//     Cell1SizeCmd.SetParameterName("cellsize", true);
//     Cell1SizeCmd.SetRange("cellsize>=0. && cellsize<100.");
//     Cell1SizeCmd.SetDefaultValue("2.");


//     auto& Layer1ZPositionCmd
//       = fMessenger->DeclareMethodWithUnit("Layer1Z","mm",
//                                   &codeApertureDetectorConstruction::SetDetZPosition,
//                                   "Set z position of the Layer1.");
//     Layer1ZPositionCmd.SetParameterName("ZPosition", true);
//     Layer1ZPositionCmd.SetRange("ZPosition>=-100. && ZPosition<100.");
//     Layer1ZPositionCmd.SetDefaultValue("-18");

// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
