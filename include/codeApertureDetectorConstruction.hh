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
/// \file codeApertureDetectorConstruction.hh
/// \brief Definition of the codeApertureDetectorConstruction class

#ifndef codeApertureDetectorConstruction_h
#define codeApertureDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;

class G4GenericMessenger;
/// Detector construction class to define materials and geometry.

class codeApertureDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    codeApertureDetectorConstruction();
    virtual ~codeApertureDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // void SetMaskThickness(G4double val);
    // G4double GetMaskThickness() const { return fThicknessMask; }

    // void SetDetThickness(G4double val);
    // G4double GetDetThickness() const { return fThicknessDet; }

    // void SetCellThickness(G4double val);
    // G4double GetCellThickness() const { return fThicknessCell; }

    // void SetDetZPosition(G4double val);


    // void SetCellSize(G4double val);


    // G4LogicalVolume* GetMaskLogicalVolume() const { return fLogMask; }
    // G4LogicalVolume* GetDetLogicalVolume() const { return fLogDet; }

    void ConstructMaterial();

  private:
    //void DefineCommands();
    G4GenericMessenger* fMessenger;

    G4Box* fSolidMask, *fSolidDet;

    G4Box* fSolidCell;

    G4VPhysicalVolume* fPhysMask, *fPhysDet;

    std::vector<G4VPhysicalVolume*> fPhysCell;

    G4LogicalVolume*  fLogMask, *fLogDet;

    G4LogicalVolume* fLogCell;

    G4LogicalVolume* fLogicalWorld;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

