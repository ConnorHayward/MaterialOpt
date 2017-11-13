#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Cons.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include <iostream>
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),fPBox(nullptr), fLBox(nullptr),
  fBox(nullptr)
{
  fDetectorMessenger = new DetectorMessenger(this);
  fExpHall_x = fExpHall_y = fExpHall_z = 0.5*m;
  fTargetName = "holder";
  fThickness = 1*mm;
  fTargetThickness = 10*mm;
  fDetectorType = 0;
  DefineMaterials();
  SetTargetMaterial("PEN");
  SetWorldMaterial("air");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value){
  fTargetThickness=value;
  if(fBox){
    fBox->SetZHalfLength(fTargetThickness/2);
  }
  UpdateGeometry();

  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetDetectorType(G4int value){
  fDetectorType=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fTargetMaterial = pttoMaterial;
    fTargetName = fTargetMaterial->GetName();
    if ( fLBox ) { fLBox->SetMaterial(fTargetMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if ( fWLBox ) { fWLBox->SetMaterial(fWorldMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::DefineMaterials(){// ------------- Materials -------------
  G4double a, z, density;
  G4int nelements;

  // air
  //
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  G4Material* air = new G4Material("air", density=1.29*mg/cm3, nelements=2);
  air->AddElement(N, 70.*perCent);
  air->AddElement(O, 30.*perCent);

  G4NistManager* man = G4NistManager::Instance();

  G4Material* glass = man->FindOrBuildMaterial("G4_Pyrex_Glass");
  G4Material* teflon = man->FindOrBuildMaterial("G4_TEFLON");

  G4MaterialPropertiesTable* glassMPT = new G4MaterialPropertiesTable();
  G4double photonEnergy2[] = {2.479684*eV, 2.610194*eV, 2.755204*eV, 2.883353*eV, 2.917275*eV, 3.099605*eV};
  G4double refractiveIndex4[] = {1.52, 1.52, 1.52, 1.52, 1.52, 1.52};
  G4double absorption[] = {10*m,10*m,10*m,10*m,10*m,10*m};
  assert(sizeof(refractiveIndex4) == sizeof(photonEnergy2));
  glassMPT->AddProperty("RINDEX",photonEnergy2,refractiveIndex4,6)->SetSpline(true);
  glassMPT->AddProperty("ABSLENGTH",photonEnergy2,absorption,6);
  glass->SetMaterialPropertiesTable(glassMPT);

  G4Material* aluminium = man->FindOrBuildMaterial("G4_Al");
  G4Material* Si = man->FindOrBuildMaterial("G4_Si");
  G4Material* Pb = man->FindOrBuildMaterial("G4_Pb");
  G4Material* Cu = man->FindOrBuildMaterial("G4_Cu");

  // Water
  //
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);

  G4Element* C = new G4Element("Carbon", "C", z=12, a=12*g/mole);

  G4Material* PEN = new G4Material("PEN", density= 1.3*g/cm3, nelements=3);
  G4int number_of_atoms;
  PEN->AddElement(O, number_of_atoms=4);
  PEN->AddElement(H, number_of_atoms=10);
  PEN->AddElement(C, number_of_atoms=14);

  G4double wavelength;
  char filler;
  G4double varabsorlength;
  G4double ems;
  G4double rindex;

  G4double absEnergy[101]  = {0};
  G4double abs[101]={0};
  G4double emission[101]={0};
  G4double rIndex[101]={0};
  G4double rIndex_air[101]={0};

  G4int absEntries = 0;
  ifstream ReadAbs;

  G4String abs_file = "../input_files/Exp4.csv";
  ReadAbs.open(abs_file);

  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>varabsorlength>>filler>>ems>>filler>>rindex;
      if(ReadAbs.eof()){
        break;
      }
      absEnergy[absEntries] = (1240/wavelength)*eV;
      abs[absEntries] = varabsorlength*mm;
      emission[absEntries] = ems;
      rIndex[absEntries] = rindex;
      rIndex_air[absEntries]=1.0;
      absEntries++;
    }
  }

  else G4cout<<"Error opening file: " <<abs_file<<G4endl;
  ReadAbs.close();
  absEntries--;

  const G4int nEntries1 = sizeof(absEnergy)/sizeof(G4double);
  assert(sizeof(rIndex) == sizeof(absEnergy));
  assert(sizeof(abs) == sizeof(absEnergy));
  assert(sizeof(emission) == sizeof(absEnergy));
  assert(sizeof(rIndex_air == sizeof(absEnergy)));

  G4MaterialPropertiesTable* penMPT = new G4MaterialPropertiesTable();

  penMPT->AddProperty("RINDEX",       absEnergy, rIndex, nEntries1)->SetSpline(true);
  penMPT->AddProperty("ABSLENGTH",    absEnergy, abs, nEntries1)->SetSpline(true);
  penMPT->AddProperty("FASTCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);
  penMPT->AddProperty("SLOWCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);

  penMPT->AddConstProperty("SCINTILLATIONYIELD",10500./MeV);
  penMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
  penMPT->AddConstProperty("FASTTIMECONSTANT", 5.198*ns);
  penMPT->AddConstProperty("SLOWTIMECONSTANT",24.336*ns);
  penMPT->AddConstProperty("YIELDRATIO",0.05);

  PEN->SetMaterialPropertiesTable(penMPT);
  density = universe_mean_density;    //from PhysicalConstants.h
  G4Material* vacuum = new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                           kStateGas,2.73*kelvin,3.e-18*pascal);
  //
  // air
  //

  G4MaterialPropertiesTable* worldMPT = new G4MaterialPropertiesTable();
  worldMPT->AddProperty("RINDEX", absEnergy, rIndex_air, nEntries1)->SetSpline(true);

  air->SetMaterialPropertiesTable(worldMPT);
  vacuum->SetMaterialPropertiesTable(worldMPT);
}

void DetectorConstruction::SetPropertyTable(G4Material* material, G4MaterialPropertiesTable* table){
  material->SetMaterialPropertiesTable(table);
}
void DetectorConstruction::UpdateGeometry(){
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    DefineMaterials();
//
// ------------- Volumes --------------

// The experimental Hall
//
  fWorldBox = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);

  fWLBox = new G4LogicalVolume(fWorldBox,fWorldMaterial,"World",0,0,0);

  fWPBox = new G4PVPlacement(0,G4ThreeVector(),fWLBox,"World",0,false,0);

//  fThickness=1*mm;
  double thickness_foil = 0.5*mm;
  double holeHeight = fThickness;
  G4Box* absorberPb_box = new G4Box("absorber",50*mm,50*mm, fThickness);
  G4Tubs* hole = new G4Tubs("hole",0*mm,1.5*mm,holeHeight, 0.*deg, 360.*deg);
  G4SubtractionSolid* collimator = new G4SubtractionSolid("collimator",absorberPb_box,hole);
  // G4LogicalVolume* absorberPb_log = new G4LogicalVolume(collimator,teflon,"absorber",0,0,0);
//  G4VPhysicalVolume* absorberPb_phys = new G4PVPlacement(0,G4ThreeVector(),absorberPb_log,"absorber",fWLBox,false,0);

  G4Box* absorber_foil = new G4Box("foil",50*mm,50*mm,thickness_foil);
  G4Tubs* hole_foil = new G4Tubs("hole_foil",0*mm,1.5*mm,thickness_foil, 0.*deg, 360.*deg);
  G4SubtractionSolid* collimator_foil = new G4SubtractionSolid("collimator_foil",absorber_foil,hole_foil);
  // G4LogicalVolume* absorberCu_log = new G4LogicalVolume(collimator_foil,aluminium,"foil",0,0,0);
  //G4VPhysicalVolume* absorberCu_phys = new G4PVPlacement(0,G4ThreeVector(0,0,fThickness+(thickness_foil)),absorberCu_log,"foil",absorberPb_log,false,0);

//  fTargetThickness = 2*cm;
  double target_width = 15*mm;
  fBox = new G4Box("target", target_width, target_width, fTargetThickness);
  fLBox = new G4LogicalVolume(fBox,fTargetMaterial, "target",0,0,0);
  double position = 0;
  fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,position),fLBox,"target",fWLBox,false,0);

  // --------------Detectors--------------

  G4NistManager* man = G4NistManager::Instance();

  G4Material* glass = man->FindOrBuildMaterial("G4_Pyrex_Glass");
  G4Material* aluminium = man->FindOrBuildMaterial("G4_Al");
  G4Material* Si = man->FindOrBuildMaterial("G4_Si");
  G4Material* Pb = man->FindOrBuildMaterial("G4_Pb");
  G4Material* Cu = man->FindOrBuildMaterial("G4_Cu");

  G4MaterialPropertiesTable* glassMPT = new G4MaterialPropertiesTable();
  G4double photonEnergy2[] = {2.479684*eV, 2.610194*eV, 2.755204*eV, 2.883353*eV, 2.917275*eV, 3.099605*eV};
  G4double refractiveIndex4[] = {1.52, 1.52, 1.52, 1.52, 1.52, 1.52};
  G4double absorption[] = {10*m,10*m,10*m,10*m,10*m,10*m};
  assert(sizeof(refractiveIndex4) == sizeof(photonEnergy2));
  glassMPT->AddProperty("RINDEX",photonEnergy2,refractiveIndex4,6)->SetSpline(true);
  glassMPT->AddProperty("ABSLENGTH",photonEnergy2,absorption,6);
  glass->SetMaterialPropertiesTable(glassMPT);

  // PMT
  char filler;
  G4double wavelength;
  G4double innerRadius_pmt = 0.*cm;
  G4double outerRadius_pmt = 53*mm;
  G4double outerRadius_cath = 46*mm;
  G4double height_pmt = 63.5*mm;
  G4double height_cath = 62.*mm;
  G4double startAngle_pmt = 0.*deg;
  G4double spanningAngle_pmt = 360.*deg;

  G4Tubs* pmt = new G4Tubs("pmt_tube",innerRadius_pmt,outerRadius_pmt, height_pmt,startAngle_pmt,spanningAngle_pmt);
  G4LogicalVolume* pmt_log = new G4LogicalVolume(pmt,glass, "pmt_log");
  G4Tubs* Photocath = new G4Tubs("photocath_tube",innerRadius_pmt,outerRadius_cath,
                          height_cath,startAngle_pmt,spanningAngle_pmt);
  G4LogicalVolume* photocath_log = new G4LogicalVolume(Photocath, aluminium,"photocath_log");

  G4double cath_eff;
  G4double photocath_energy[57];
  G4double photocath_EFF[57];
  G4double photocath_REFL[57]={0};
  G4String pmt_file = "../input_files/pmtQE.csv";

  ifstream ReadEff;
  G4int effCounter = 0;
  ReadEff.open(pmt_file);

  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>cath_eff;
      if(ReadEff.eof()){
        break;
      }
      photocath_energy[57-effCounter] = (1240/wavelength)*eV;
      photocath_EFF[57-effCounter] = cath_eff;
      effCounter++;
    }
  }
  else G4cout<<"Error opening file: " <<pmt_file<<G4endl;
  ReadEff.close();
  effCounter--;

  const G4int nPMT_EFF = sizeof(photocath_energy)/sizeof(G4double);

  G4OpticalSurface* photocath_optsurf = new G4OpticalSurface("photocath_opsurf",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* photocath_MT = new G4MaterialPropertiesTable();
  photocath_MT->AddProperty("EFFICIENCY", photocath_energy, photocath_EFF,nPMT_EFF);
  photocath_MT->AddProperty("REFLECTIVITY", photocath_energy, photocath_REFL,nPMT_EFF);
  photocath_optsurf->SetMaterialPropertiesTable(photocath_MT);
  new G4LogicalSkinSurface("photocath_surf",photocath_log,photocath_optsurf);

  // SiPM
  double siPM_thickness = 0.25*mm;
  G4Box* siPM_box = new G4Box("siPM_box", 3*mm, 3*mm, siPM_thickness);
  G4LogicalVolume* siPM_log = new G4LogicalVolume(siPM_box, Si, "siPM_log");

  G4double sipm_eff;
  G4double siPM_energy[56];
  G4double siPM_EFF[56];
  G4double siPM_REFL[56]={0};
  G4String siPM_file = "../input_files/sipmQE.csv";

  effCounter = 0;

  ReadEff.open(siPM_file);

  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>cath_eff;
      if(ReadEff.eof()){
        break;
      }
      siPM_energy[56-effCounter] = (1240/wavelength)*eV;
      siPM_EFF[56-effCounter] = cath_eff;
      effCounter++;
    }
  }
  else G4cout<<"Error opening file: " <<siPM_file<<G4endl;
  ReadEff.close();
  effCounter--;

  const G4int nSIPM_EFF = sizeof(siPM_energy)/sizeof(G4double);

  G4OpticalSurface* siPM_optsurf = new G4OpticalSurface("siPM_optsurf",glisur, polished, dielectric_metal);
  G4MaterialPropertiesTable* siPM_MT = new G4MaterialPropertiesTable();
  siPM_MT->AddProperty("EFFICIENCY", siPM_energy, siPM_EFF,nSIPM_EFF);
  siPM_MT->AddProperty("REFLECTIVITY", siPM_energy, siPM_REFL,nSIPM_EFF);
  siPM_optsurf->SetMaterialPropertiesTable(siPM_MT);
  new G4LogicalSkinSurface("siPM_surf", siPM_log, siPM_optsurf);

  G4VPhysicalVolume* pmtPlacement;
  G4VPhysicalVolume* cathPlacement;
  G4VPhysicalVolume* siPMPlacement1;
  G4VPhysicalVolume* siPMPlacement2;
  G4RotationMatrix* rotm = new G4RotationMatrix();
  G4ThreeVector* pmtVector = new G4ThreeVector(0,0,height_pmt+fTargetThickness);
  switch (fDetectorType){
    case 0:
      pmtPlacement = new G4PVPlacement(0,*pmtVector,pmt_log,"pmt",fWLBox,false,0);
      cathPlacement = new G4PVPlacement(0,G4ThreeVector(),photocath_log,"detector",pmt_log,false,0);
      break;
    case 1:
      rotm->rotateX(90*deg);
      siPMPlacement1 = new G4PVPlacement(rotm,G4ThreeVector(0,siPM_thickness+target_width,0),siPM_log,"detector",fWLBox,false,0);
      siPMPlacement2 = new G4PVPlacement(rotm,G4ThreeVector(0,-1*(siPM_thickness+target_width),0),siPM_log,"detector",fWLBox,false,0);
      break;
  }
  return fWPBox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
