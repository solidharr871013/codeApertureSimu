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
/// \file codeApertureLayer1Hit.hh
/// \brief Definition of the codeApertureLayer1Hit class

#ifndef codeApertureLayer1Hit_h
#define codeApertureLayer1Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Drift chamber hit
///
/// It records:
/// - the layer ID
/// - the particle time
/// - the particle local and global positions

class codeApertureLayer1Hit : public G4VHit
{
  public:
    codeApertureLayer1Hit();
    //codeApertureLayer1Hit(G4int layerID);
    codeApertureLayer1Hit(const codeApertureLayer1Hit &right);
    virtual ~codeApertureLayer1Hit();

    const codeApertureLayer1Hit& operator=(const codeApertureLayer1Hit &right);
    G4bool operator==(const codeApertureLayer1Hit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();

    void SetCellID(G4int z) { fCellID = z; }
    G4int GetCellID() const { return fCellID; }

    void SetEdep(G4double edep) { fEdep = edep; }
    void AddEdep(G4double edep){fEdep += edep;}
    G4double GetEdep() const { return fEdep; }

    void SetScatteringPos(G4ThreeVector xyz) { fScatteringPos = xyz; }
    G4ThreeVector GetScatteringPos() const { return fScatteringPos; }

    void SetProcessName(G4String procName) { fProcessName = procName; }
    G4String GetProcessName() const { return fProcessName; }

    void SetParticleName(G4String pName) { fParticleName = pName; }
    G4String GetParticleName() const { return fParticleName; }

    void SetVolumeName(G4String volName) { fVolumeName = volName; }
    G4String GetVolumeName() const { return fVolumeName; }

    void SetTime(G4double time){fTime = time;}
    G4double GetTime() const {return fTime;}

    void SetParentID(G4int PID){fParentID = PID;}
    G4int GetParentID() const {return  fParentID;}
    
  private:
    G4double fEdep;
    G4double fTime;
    G4ThreeVector fScatteringPos;
    G4String fProcessName;
    G4String fVolumeName;
    G4String fParticleName;
    G4int fCellID;
    G4int fParentID;

};

using codeApertureLayer1HitsCollection = G4THitsCollection<codeApertureLayer1Hit>;

extern G4ThreadLocal G4Allocator<codeApertureLayer1Hit>* codeApertureLayer1HitAllocator;

inline void* codeApertureLayer1Hit::operator new(size_t)
{
  if (!codeApertureLayer1HitAllocator) {
       codeApertureLayer1HitAllocator = new G4Allocator<codeApertureLayer1Hit>;
  }
  return (void*)codeApertureLayer1HitAllocator->MallocSingle();
}

inline void codeApertureLayer1Hit::operator delete(void* aHit)
{
  codeApertureLayer1HitAllocator->FreeSingle((codeApertureLayer1Hit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
