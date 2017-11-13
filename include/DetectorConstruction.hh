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
/// \file OpNovice/include/OpNoviceDetectorConstruction.hh
/// \brief Definition of the OpNoviceDetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialPropertiesTable.hh"
class DetectorMessenger;
class G4LogicalVolume;
class G4Material;
class G4Box;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    void SetSize  (G4double);
    void SetTargetMaterial(G4String);
    void SetWorldMaterial(G4String);
    void SetPMTPlacement(G4ThreeVector);
    void UpdateGeometry(void);
    void SetPropertyTable(G4Material* mat, G4MaterialPropertiesTable* tab);
    void SetDetectorType(G4int);
    G4int GetDetectorType(){return fDetectorType;};


  public:
    const G4VPhysicalVolume* GetWorld() {return fPBox;};

    G4double GetSize()  {return fThickness;};
    G4double GetTargetSize()  {return fTargetThickness;};
    G4String GetTargetMaterialName(){return fTargetName;};
    G4Material*        GetWorldMaterial()   {return fWorldMaterial;};
    G4Material*        GetTargetMaterial()   {return fTargetMaterial;};

    void               DefineMaterials();

  private:
    G4double fExpHall_x;
    G4double fExpHall_y;
    G4double fExpHall_z;

    G4double fThickness;
    G4double fTargetThickness;
    G4MaterialPropertiesTable* fTargetMPT;
    G4MaterialPropertiesTable* fWorldMPT;

    G4Material* fWorldMaterial;
    G4Material* fTargetMaterial;
    G4String fTargetName;
    G4VPhysicalVolume* fPBox;
    G4LogicalVolume* fLBox;
    G4Box* fBox;

    G4VPhysicalVolume* fWPBox;
    G4LogicalVolume* fWLBox;
    G4Box* fWorldBox;
    G4ThreeVector fPMTPlacement = G4ThreeVector(0,0,63.5*mm+fTargetThickness);

    G4int fDetectorType;

    DetectorMessenger* fDetectorMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*DetectorConstruction_h*/
