#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Transform3D.hh"
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

/*
Constructs DetectorConstruction, defines default values.
*/
DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),fPBox(nullptr), fLBox(nullptr),
  fBox(nullptr)
{
  fDetectorMessenger = new DetectorMessenger(this);
  fTargetMPT = new G4MaterialPropertiesTable();
  fExpHall_x = fExpHall_y = fExpHall_z = 1*m;
  fTargetName = "holder";
  fThickness = 1*mm;
  fTargetThickness = 10*mm;
  fDetectorType = 0;
  DefineMaterials();
  SetTargetMaterial("PEN");
  SetWorldMaterial("Air");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
Sets thickness of target.
*/
void DetectorConstruction::SetSize(G4double value){
  fTargetThickness=value;
  if(fBox){
    fBox->SetZHalfLength(fTargetThickness/2);
  }
  UpdateGeometry();

  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets which detector geometry is used.
*/
void DetectorConstruction::SetDetectorType(G4int value){
  fDetectorType=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets material of target.
*/
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

/*
Sets material of world volume.
*/
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

/*
Defines materials used in simulation. Sets material properties for PEN and other optical components.
*/
void DetectorConstruction::DefineMaterials(){// ------------- Materials -------------
  G4double a, z, density;
  G4int nelements;

  // fAir
  //
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  fAir = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  fAir->AddElement(N, 70.*perCent);
  fAir->AddElement(O, 30.*perCent);

  G4NistManager* man = G4NistManager::Instance();

  fGlass = man->FindOrBuildMaterial("G4_Pyrex_Glass");
  fTeflon = man->FindOrBuildMaterial("G4_TEFLON");
  fLAr = man->FindOrBuildMaterial("G4_lAr");

  G4MaterialPropertiesTable* glassMPT = new G4MaterialPropertiesTable();
  G4double photonEnergy2[] = {2.479684*eV, 2.610194*eV, 2.755204*eV, 2.883353*eV, 2.917275*eV, 3.099605*eV};
  G4double refractiveIndex4[] = {1.52, 1.52, 1.52, 1.52, 1.52, 1.52};
  G4double absorption[] = {10*m,10*m,10*m,10*m,10*m,10*m};
  assert(sizeof(refractiveIndex4) == sizeof(photonEnergy2));
  glassMPT->AddProperty("RINDEX",photonEnergy2,refractiveIndex4,6)->SetSpline(true);
  glassMPT->AddProperty("ABSLENGTH",photonEnergy2,absorption,6);
  fGlass->SetMaterialPropertiesTable(glassMPT);

  fAluminium = man->FindOrBuildMaterial("G4_Al");
  fSi = man->FindOrBuildMaterial("G4_Si");
  fPb = man->FindOrBuildMaterial("G4_Pb");
  fCu = man->FindOrBuildMaterial("G4_Cu");
  fGe = man->FindOrBuildMaterial("G4_Ge");

  // Water
  //
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);

  G4Element* C = new G4Element("Carbon", "C", z=12, a=12*g/mole);

  fPEN = new G4Material("PEN", density= 1.3*g/cm3, nelements=3);
  G4int number_of_atoms;
  fPEN->AddElement(O, number_of_atoms=4);
  fPEN->AddElement(H, number_of_atoms=10);
  fPEN->AddElement(C, number_of_atoms=14);

  G4double wavelength;
  char filler;
  G4double varabsorlength;
  G4double ems;
  G4double rindex;

  G4double absEnergy[102]  = {0};
  G4double abs[102]={0};
  G4double emission[102]={0};
  G4double rIndex[102]={0};
  G4double rIndex_fAir[102]={0};

  G4int absEntries = 0;
  ifstream ReadAbs;

//  G4String abs_file = "../input_files/lancSpec.csv";
//  G4String abs_file = "../input_files/BC408.csv";
//  G4String abs_file = "../input_files/OldPen.csv";
  G4String abs_file = "../input_files/Exp4.csv";
//  G4String abs_file = "../input_files/flatAbs.csv";
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
      rIndex_fAir[absEntries]=1.0;
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
  assert(sizeof(rIndex_fAir == sizeof(absEnergy)));

  fTargetMPT->AddProperty("RINDEX",       absEnergy, rIndex, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("ABSLENGTH",    absEnergy, abs, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("FASTCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("SLOWCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);

  fTargetMPT->AddConstProperty("SCINTILLATIONYIELD",10500./MeV);
  fTargetMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
  fTargetMPT->AddConstProperty("FASTTIMECONSTANT", 5.198*ns);
  fTargetMPT->AddConstProperty("SLOWTIMECONSTANT",24.336*ns);
  fTargetMPT->AddConstProperty("YIELDRATIO",0.05);

  fPEN->SetMaterialPropertiesTable(fTargetMPT);
  density = universe_mean_density;    //from PhysicalConstants.h
  G4Material* fVacuum = new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                           kStateGas,2.73*kelvin,3.e-18*pascal);
  //
  // fAir
  //

  G4MaterialPropertiesTable* worldMPT = new G4MaterialPropertiesTable();
  worldMPT->AddProperty("RINDEX", absEnergy, rIndex_fAir, nEntries1)->SetSpline(true);

  fAir->SetMaterialPropertiesTable(worldMPT);
  fVacuum->SetMaterialPropertiesTable(worldMPT);

  G4MaterialPropertiesTable* lARMPT = new G4MaterialPropertiesTable();
  abs_file = "../input_files/lArScint.csv";

  ReadAbs.open(abs_file);
  absEntries = 0;
  G4double absEnergylArScint[43] = {};
  G4double emissionlArScint[43] = {};

  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>ems;
      if(ReadAbs.eof()){
        break;
      }
      absEnergylArScint[absEntries] = (1240/wavelength)*eV;
      emissionlArScint[absEntries] = ems;
      absEntries++;
    }
  }
  else G4cout<<"Error opening file: " <<abs_file<<G4endl;
  ReadAbs.close();
  absEntries--;

  G4double absEnergyLAr[]  = {1.9255*eV, 2.145*eV, 2.2704*eV, 2.4378*eV, 2.6085*eV,2.845*eV,3.0515*eV,3.397*eV};
  G4double absLAr[]={10*m,10*m,10*m,10*m,10*m,10*m,10*m,10*m};
  G4double rIndexLAr[]={1.2295, 1.2303, 1.2308,1.2316,1.2324,1.2336,1.2347,1.2367};


  lARMPT->AddProperty("RINDEX",       absEnergyLAr, rIndexLAr, 8)->SetSpline(true);
  lARMPT->AddProperty("ABSLENGTH",    absEnergyLAr, absLAr, 8)->SetSpline(true);
  lARMPT->AddProperty("FASTCOMPONENT",absEnergylArScint, emissionlArScint, 43)->SetSpline(true);
  lARMPT->AddProperty("SLOWCOMPONENT",absEnergylArScint, emissionlArScint, nEntries1)->SetSpline(true);
  lARMPT->AddConstProperty("SCINTILLATIONYIELD", 51.*MeV);
  lARMPT->AddConstProperty("FASTTIMECONSTANT", 6.2*us);
  lARMPT->AddConstProperty("SLOWTIMECONSTANT",1300*us);
  lARMPT->AddConstProperty("YIELDRATIO",0.05);
  lARMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
  fLAr->SetMaterialPropertiesTable(lARMPT);
}

void DetectorConstruction::SetPropertyTable(G4Material* material, G4MaterialPropertiesTable* table){
  material->SetMaterialPropertiesTable(table);
}

void DetectorConstruction::UpdateGeometry(){
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

/*
Clears stored geometry, then constructs all volumes that can be used in the simulation.

Builds and places volumes in world.

Defines detector sensitivities and properties.
*/
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
//    DefineMaterials();

//
// ------------- Volumes --------------

// The experimental Hall
//
  fWorldBox = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);

  fWLBox = new G4LogicalVolume(fWorldBox,fWorldMaterial,"World",0,0,0);

  fWPBox = new G4PVPlacement(0,G4ThreeVector(),fWLBox,"World",0,false,0);

//  fThickness=1*mm;
  double thickness_foil = 2.5*mm;
  double holeHeight = fThickness;
  G4Box* absorberPb_box = new G4Box("absorber",15*mm,15*mm, fThickness);
  G4Tubs* hole = new G4Tubs("hole",0*mm,1.5*mm,holeHeight, 0.*deg, 360.*deg);
  G4SubtractionSolid* collimator = new G4SubtractionSolid("collimator",absorberPb_box,hole);

  G4Box* absorber_foil = new G4Box("foil",15*mm,15*mm,thickness_foil);
  G4Tubs* hole_foil = new G4Tubs("hole_foil",0*mm,1.5*mm,thickness_foil, 0.*deg, 360.*deg);
  G4SubtractionSolid* collimator_foil = new G4SubtractionSolid("collimator_foil",absorber_foil,hole_foil);


//  fTargetThickness = 2*cm;
  double target_width = 15*mm;
  fBox = new G4Box("target", target_width, target_width, fTargetThickness);
  fLBox = new G4LogicalVolume(fBox,fTargetMaterial, "target",0,0,0);
  double position = 0;
//  fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,position),fLBox,"target",fWLBox,false,0);


G4Sphere* penShell = new G4Sphere( "penShell",  320*mm, 350*mm,  0*deg,  360*deg,  0*deg,  360*deg );
G4Sphere* dimpleSphere = new G4Sphere( "dimple",  0*mm, 30*mm,  0*deg,  360*deg,  0*deg,  360*deg );

double shift = 370*mm;
G4RotationMatrix* dimpleRot = new G4RotationMatrix();
G4ThreeVector* translation1 = new G4ThreeVector(0,0,shift);
G4ThreeVector* translation2 = new G4ThreeVector(0,shift,0);
G4ThreeVector* translation3 = new G4ThreeVector(shift,0,0);
G4SubtractionSolid* golfBall = new G4SubtractionSolid("golfBall",penShell,dimpleSphere, dimpleRot, *translation1);
double i=0;
double angle = 11.25*deg;
  // --------------Detectors--------------

  G4NistManager* man = G4NistManager::Instance();

  // Define Bottle

  double cupHeight = 40*mm;
  double coneHeight = 25*mm;
  shift = 21*mm;
  double offset = -35*mm;
  G4Tubs* cup = new G4Tubs("cup", 16*mm, 17.5*mm, cupHeight, 0*deg, 360*deg);
  G4Sphere* smallDimpleSphere = new G4Sphere( "dimple",  0*mm, 4*mm,  0*deg,  360*deg,  0*deg,  360*deg );
  G4ThreeVector* translation = new G4ThreeVector(0,shift,0);
  G4SubtractionSolid* dimpleTube = new G4SubtractionSolid("dimpleTube",cup,smallDimpleSphere, dimpleRot, G4ThreeVector());
  double j = 0;


  G4Tubs* gePuck = new G4Tubs("puck", 0*mm, 15*mm, 20*mm, 0*deg, 360*deg);
  G4Cons* guide = new G4Cons("guide", 0*mm, 3*mm, 16*mm, 17.5*mm, coneHeight, 0*deg, 360*deg);
  //G4Tubs* spacer = new G4Tubs("space",0*mm, 3*mm, cupHeight, 0*deg, 360*deg);

  //G4SubtractionSolid* bottle = new G4SubtractionSolid("bottle", cup, spacer);
  G4LogicalVolume* bottle_log = new G4LogicalVolume(cup, fTargetMaterial, "bottle_log");
  G4LogicalVolume* guide_log = new G4LogicalVolume(guide, fTargetMaterial,"guide_log");
  G4LogicalVolume* gePuck_log = new G4LogicalVolume(gePuck, fGe,"puck_log");
  G4LogicalVolume* dimpleTube_log = new G4LogicalVolume(dimpleTube,fTargetMaterial,"tube_log");

  G4LogicalVolume* collimator_logic = new G4LogicalVolume(collimator,fAluminium,"collimator",0,0,0);

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
  G4LogicalVolume* pmt_log = new G4LogicalVolume(pmt,fGlass, "pmt_log");
  G4Tubs* Photocath = new G4Tubs("photocath_tube",innerRadius_pmt,outerRadius_cath,
                          height_cath,startAngle_pmt,spanningAngle_pmt);
  G4LogicalVolume* photocath_log = new G4LogicalVolume(Photocath, fAluminium,"photocath_log");

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

   G4OpticalSurface* bottleSurf = new G4OpticalSurface("bottle_opsurf",glisur,polished,dielectric_dielectric);
  // bottleSurf->SetSigmaAlpha(0.1);
   new G4LogicalSkinSurface("bottle_surf",fLBox,bottleSurf);
  // new G4LogicalSkinSurface("guide_surf", guide_log, bottleSurf);

  G4OpticalSurface* photocath_optsurf = new G4OpticalSurface("photocath_opsurf",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* photocath_MT = new G4MaterialPropertiesTable();
  photocath_MT->AddProperty("EFFICIENCY", photocath_energy, photocath_EFF,nPMT_EFF);
  photocath_MT->AddProperty("REFLECTIVITY", photocath_energy, photocath_REFL,nPMT_EFF);
  photocath_optsurf->SetMaterialPropertiesTable(photocath_MT);
  new G4LogicalSkinSurface("photocath_surf",photocath_log,photocath_optsurf);

  G4Sphere* shadowShell = new G4Sphere( "penShell", 0*mm, 200*mm,  0*deg,  360*deg,  0*deg,  360*deg );
  G4LogicalVolume* shadowShell_log = new G4LogicalVolume(shadowShell, fGe,"shadow",0,0,0);

  G4MaterialPropertiesTable* shadowMPT = new G4MaterialPropertiesTable();
  shadowMPT->AddProperty("REFLECTIVITY", photocath_energy, photocath_REFL,nPMT_EFF);
  fGe->SetMaterialPropertiesTable(shadowMPT);

  // Spectrometer

  G4Box* spec_box = new G4Box("spec_box",12.8*mm,3.3*mm,0.25*mm);
  G4LogicalVolume* spec_log = new G4LogicalVolume(spec_box,fSi,"spec_log");

  G4double spec_eff;
  G4double spec_energy[65];
  G4double spec_EFF[65];
  G4double spec_REFL[65]={0};
  G4String spec_file = "../input_files/specQE.csv";

  effCounter = 0;
  ReadEff.open(spec_file);

  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>spec_eff;
      if(ReadEff.eof()){
        break;
      }
      spec_energy[65-effCounter]=(1240/wavelength)*eV;
      spec_EFF[65-effCounter]=spec_eff;
      effCounter++;
    }
  }
  else G4cout<<"Error opening file: " <<spec_file<<G4endl;
  ReadEff.close();
  effCounter--;

  const G4int nSPEC_EFF = sizeof(spec_energy)/sizeof(G4double);

  G4OpticalSurface* spec_optsurf = new G4OpticalSurface("spectrometer_opsurf",glisur, polished, dielectric_metal);
  G4MaterialPropertiesTable* spec_MT = new G4MaterialPropertiesTable();
  spec_MT->AddProperty("EFFICIENCY", spec_energy,spec_EFF,nSPEC_EFF);
  spec_MT->AddProperty("REFLECTIVITY",spec_energy,spec_REFL,nSPEC_EFF);
  spec_optsurf->SetMaterialPropertiesTable(spec_MT);
  new G4LogicalSkinSurface("spec_surf", spec_log,spec_optsurf);

  G4Box* slit_box = new G4Box("slit_box",12.8*mm,3.3*mm,1*mm);
  G4Box* slit_gap = new G4Box("slit_gap", 0.9*mm, 3.3*mm, 1*mm);
  G4SubtractionSolid* slits = new G4SubtractionSolid("slit",slit_box,slit_gap);
  G4LogicalVolume* slit_log = new G4LogicalVolume(slits,fTeflon,"slit_log");

  // SiPM
  double siPM_thickness = 0.25*mm;
  G4Box* siPM_box = new G4Box("siPM_box", 1.5*mm, 1.5*mm, siPM_thickness);
  G4LogicalVolume* siPM_log = new G4LogicalVolume(siPM_box, fSi, "siPM_log");

  G4Sphere* detShell = new G4Sphere( "detShell",  450*mm, 500*mm,  0*deg,  360*deg,  0*deg,  360*deg );
  G4LogicalVolume* detShell_log =new G4LogicalVolume(detShell, fSi,"target",0,0,0);

  G4double sipm_eff;
  G4double siPM_energy[56];
  G4double siPM_EFF[56];
  G4double siPM_REFL[56]={0};
  G4double siMax[56];
  G4double tefRI[56];
  G4double tefRef[56];
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
      siMax[effCounter]=1;
      tefRI[effCounter]=1.3;
      tefRef[effCounter]=1;
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

  G4OpticalSurface* max_optsurf = new G4OpticalSurface("max_optsurf",glisur, polished, dielectric_metal);
  G4MaterialPropertiesTable* max_MT = new G4MaterialPropertiesTable();
  max_MT->AddProperty("EFFICIENCY", siPM_energy, siMax,nSIPM_EFF);
  max_MT->AddProperty("REFLECTIVITY", siPM_energy, siPM_REFL,nSIPM_EFF);
  max_optsurf->SetMaterialPropertiesTable(max_MT);
  new G4LogicalSkinSurface("detShell_surf", detShell_log, max_optsurf);

  G4RotationMatrix* rotm = new G4RotationMatrix();

  // Reflectors

  G4Box* detectorBox = new G4Box("detector",target_width,target_width,siPM_thickness);
  G4LogicalVolume* detectorLogical = new G4LogicalVolume(detectorBox,fSi, "detector",0,0,0);
  G4Box* reflectorBox = new G4Box("reflector",target_width+3*mm, target_width+2*mm, fTargetThickness+1*mm);
  G4SubtractionSolid* foil = new G4SubtractionSolid("foil", reflectorBox, fBox);
  G4SubtractionSolid* reflector = new G4SubtractionSolid("foil",foil, siPM_box,rotm,G4ThreeVector(0,0,siPM_thickness+fTargetThickness));
  G4LogicalVolume* ref_log = new G4LogicalVolume(reflector, fTeflon, "reflector",0,0,0);
  G4SubtractionSolid* reflector2 = new G4SubtractionSolid("foil",reflectorBox,fBox, rotm, G4ThreeVector(0,-2*mm, 0));
  G4LogicalVolume* ref_log2 = new G4LogicalVolume(reflector2, fTeflon, "reflector",0,0,0);


  G4Box* holder = new G4Box("siPM_box", 1.5*mm, 1.5*mm, 1*mm);
  G4SubtractionSolid* siPM_reflector = new G4SubtractionSolid("foil",foil, holder,rotm,G4ThreeVector(0,0,2*mm));
  G4LogicalVolume* siPM_ref = new G4LogicalVolume(siPM_reflector, fTeflon, "reflector",0,0,0);

  G4OpticalSurface* tef_optsurf = new G4OpticalSurface("tef_optsurf", glisur, polished, dielectric_metal);
  G4MaterialPropertiesTable* tefMPT = new G4MaterialPropertiesTable();
  tefMPT->AddProperty("REFLECTIVITY", siPM_energy,tefRef, nSIPM_EFF);
  tefMPT->AddProperty("RINDEX",siPM_energy, tefRI, nSIPM_EFF);
  tef_optsurf->SetMaterialPropertiesTable(tefMPT);
  new G4LogicalSkinSurface("tefSurf", ref_log, tef_optsurf);
  new G4LogicalSkinSurface("tefSurf", siPM_ref,tef_optsurf);


  G4OpticalSurface* detector_optsurf = new G4OpticalSurface("detector_optsurf",glisur, polished, dielectric_metal);
  G4MaterialPropertiesTable* detector_MT = new G4MaterialPropertiesTable();
  detector_MT->AddProperty("EFFICIENCY", siPM_energy, siMax,nSIPM_EFF);
  detector_MT->AddProperty("REFLECTIVITY", siPM_energy, siPM_REFL,nSIPM_EFF);
  detector_optsurf->SetMaterialPropertiesTable(detector_MT);
  new G4LogicalSkinSurface("siPM_surf", detectorLogical, detector_optsurf);

  G4VPhysicalVolume* pmtPlacement;
  G4VPhysicalVolume* cathPlacement;
  G4VPhysicalVolume* siPMPlacement1;
  G4VPhysicalVolume* siPMPlacement2;
  G4VPhysicalVolume* reflectorPlacement;
  G4VPhysicalVolume* bottlePlacement;
  G4VPhysicalVolume* slitPlacement;
  G4VPhysicalVolume* collimator_phys;
  G4VPhysicalVolume* coll2;
  G4VPhysicalVolume* spec_phys;
  G4VPhysicalVolume* slit_phys;

  G4ThreeVector* pmtVector = new G4ThreeVector(0,0,height_pmt+fTargetThickness);
  G4ThreeVector* siPMVector = new G4ThreeVector(0,0,siPM_thickness+fTargetThickness);

  /*
  0 - PMT on base of tile, collimator included.
  1 - 2 SiPMs attached to the sides of a tile.
  2 - Prague PMT set up. Wrapped tile, coupled to PMT on the side. Collimator included above and below.
  3 - Basic bottle shape. No detector, for visual tests.
  4 - PMT coupled to tile on side. Tile is wrapped.
  5 - SiPM coupled to tile on side.
  6 - MPI Spectrometer set up. SiPM with spectrometer settings placed a short distance from edge of tile. Slits are included.
  7 - Wrapped tile with SiPM on base.
  8 - Bottle shape with SiPM at one end. Test case for shape.
  9 - 'Golfball' shape, surrounded by ideal spherical detector. A spherical shadow is placed in the centre of the PEN shell.
  10 - Cylinder with spherical sections removed.
  */
  switch (fDetectorType){
    case 0:
      fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,position),fLBox,"target",fWLBox,false,0);
      pmtPlacement = new G4PVPlacement(0,*pmtVector,pmt_log,"pmt",fWLBox,false,0);
      cathPlacement = new G4PVPlacement(0,G4ThreeVector(),photocath_log,"detector",pmt_log,false,0);
      collimator_phys = new G4PVPlacement(0,G4ThreeVector(0,0,-0.75*cm),collimator_logic,"collimator",fWLBox, false, 0);
      break;
    case 1:
      rotm->rotateX(90*deg);
      fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,position),fLBox,"target",fWLBox,false,0);
      siPMPlacement1 = new G4PVPlacement(rotm,G4ThreeVector(0,siPM_thickness+target_width,0),siPM_log,"detector",fWLBox,false,0);
      siPMPlacement2 = new G4PVPlacement(rotm,G4ThreeVector(0,-1*(siPM_thickness+target_width),0),siPM_log,"detector",fWLBox,false,0);
      break;
    case 2:
      rotm->rotateX(90*deg);
      fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,position),fLBox,"target",fWLBox,false,0);
      reflectorPlacement = new G4PVPlacement(0,G4ThreeVector(0,-3*mm,0), ref_log,"reflector",fWLBox, false, 0);
      pmtPlacement = new G4PVPlacement(rotm,G4ThreeVector(0,height_pmt+target_width,0),pmt_log,"pmt",fWLBox,false,0);
      cathPlacement = new G4PVPlacement(0,G4ThreeVector(),photocath_log,"detector",pmt_log,false,0);
      collimator_phys = new G4PVPlacement(0,G4ThreeVector(0,0,-0.75*cm),collimator_logic,"collimator",fWLBox, false, 0);
      coll2 = new G4PVPlacement(0,G4ThreeVector(0,0,0.75*cm),collimator_logic,"collimator",fWLBox, false, 0);
      break;
    case 3:
      bottlePlacement = new G4PVPlacement(0,G4ThreeVector(0,0,0),bottle_log,"target",fWLBox,false,0);
      siPMPlacement1 = new G4PVPlacement(0,G4ThreeVector(), detShell_log,"detector",fWLBox,false,0);
      break;
    case 4:
      rotm->rotateX(90*deg);
      fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,position),fLBox,"target",fWLBox,false,0);
      pmtPlacement = new G4PVPlacement(0,G4ThreeVector(0,0,height_pmt+fTargetThickness),pmt_log,"pmt",fWLBox,false,0);
      cathPlacement = new G4PVPlacement(0,G4ThreeVector(),photocath_log,"detector",pmt_log,false,0);
      reflectorPlacement = new G4PVPlacement(0,G4ThreeVector(0,0,0), ref_log,"reflector",fWLBox, false, 0);
      break;
    case 5:
      rotm->rotateX(90*deg);
      fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,position),fLBox,"target",fWLBox,false,0);
      siPMPlacement1 = new G4PVPlacement(rotm, G4ThreeVector(0,40*mm+target_width+siPM_thickness,0),spec_log,"detector",fWLBox,false,0);
      break;
    case 6:
      rotm->rotateX(90*deg);
      spec_phys = new G4PVPlacement(rotm,G4ThreeVector(0,siPM_thickness+target_width+1*cm,0),spec_log,"detector",fWLBox,false,0);
      break;
    case 7:
      reflectorPlacement = new G4PVPlacement(0,G4ThreeVector(0,0,0), ref_log,"reflector",fWLBox, false, 0);
      siPMPlacement1 = new G4PVPlacement(rotm, G4ThreeVector(0,0,siPM_thickness+fTargetThickness),siPM_log,"detector",fWLBox,false,0);
      break;
    case 8:
      fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,position),bottle_log,"target",fWLBox,false,0);
      bottlePlacement = new G4PVPlacement(0,G4ThreeVector(),gePuck_log,"blocker", fWLBox,false, 0);
      coll2 = new G4PVPlacement(0,G4ThreeVector(0,0,-1*(cupHeight+coneHeight)),guide_log,"guide",fWLBox,false,0);
      siPMPlacement1 = new G4PVPlacement(rotm, G4ThreeVector(0,0,-1*(coneHeight+siPM_thickness)),siPM_log,"detector",guide_log,false,0);
      break;
    case 9:
      angle = 11.25*deg;
      shift = 370*mm;
      G4cout<<"Test"<<G4endl;
      i=0;
      golfBall = new G4SubtractionSolid("golfBall",penShell,dimpleSphere, dimpleRot, *translation1);
      while(i<15){
        translation1->rotateX(22.5*deg);
        golfBall = new G4SubtractionSolid("golfBall",golfBall,dimpleSphere, dimpleRot, *translation1);
        i++;
      }
      translation1->rotateX(22.5*deg);
      translation1->rotateY(90*deg);
      golfBall = new G4SubtractionSolid("golfBall",golfBall,dimpleSphere, dimpleRot, *translation1);
      golfBall = new G4SubtractionSolid("golfBall",golfBall,dimpleSphere, dimpleRot, *translation1*-1);

      i=0;

      while(i<16){
        *translation1 = G4ThreeVector(0,0,shift);

        if((i*22.5) == 0 ||(i*22.5) == 90 || (i*22.5) == 180 || (i*22.5) == 270){
          i++;
        }

        translation1->rotateY(i*22.5*deg);
        golfBall = new G4SubtractionSolid("golfBall",golfBall,dimpleSphere, dimpleRot, *translation1);
        i++;
      }
      *translation1 = G4ThreeVector(0,0,shift);

      i=0;

      while(i<16){
        *translation2 =  G4ThreeVector(0,shift,0);

        if((i*22.5) == 0 ||(i*22.5) == 90 || (i*22.5) == 180 || (i*22.5) == 270){
          i++;
        }

        translation2->rotateZ(i*22.5*deg);
        golfBall = new G4SubtractionSolid("golfBall",golfBall,dimpleSphere, dimpleRot, *translation2);
        i++;
      }

      i=0;

      while(i<16){
        translation3 = new G4ThreeVector(shift,0,0);
        translation3->rotateZ(22.5*deg);
        if((i*22.5) == 0 ||(i*22.5) == 90 || (i*22.5) == 180 || (i*22.5) == 270){
          i++;
        }
        translation3->rotateX(i*22.5*deg);
        golfBall = new G4SubtractionSolid("golfBall",golfBall,dimpleSphere, dimpleRot, *translation3);
        golfBall = new G4SubtractionSolid("golfBall",golfBall,dimpleSphere, dimpleRot, *translation3*-1);
        i++;
      }

      i=0;

      while(i<16){
        translation3 = new G4ThreeVector(shift,0,0);
        translation3->rotateZ(45*deg);
        if((i*22.5) == 0 ||(i*22.5) == 90 || (i*22.5) == 180 || (i*22.5) == 270){
          i++;
        }

        translation3->rotateX(i*22.5*deg);
        golfBall = new G4SubtractionSolid("golfBall",golfBall,dimpleSphere, dimpleRot, *translation3);
        golfBall = new G4SubtractionSolid("golfBall",golfBall,dimpleSphere, dimpleRot, *translation3*-1);
        i++;
      }

      i=0;

      while(i<16){
        translation3 = new G4ThreeVector(shift,0,0);
        translation3->rotateZ(67.5*deg);
        if((i*22.5) == 0 ||(i*22.5) == 90 || (i*22.5) == 180 || (i*22.5) == 270){
          i++;
        }

        translation3->rotateX(i*22.5*deg);
        golfBall = new G4SubtractionSolid("golfBall",golfBall,dimpleSphere, dimpleRot, *translation3);
        golfBall = new G4SubtractionSolid("golfBall",golfBall,dimpleSphere, dimpleRot, *translation3*-1);
        i++;
      }

      fLBox = new G4LogicalVolume(golfBall, fTargetMaterial,"target",0,0,0);
      fPBox = new G4PVPlacement(0,G4ThreeVector(0,0,0), fLBox,"target", fWLBox,false, 0);
      siPMPlacement1 = new G4PVPlacement(0,G4ThreeVector(), detShell_log,"detector",fWLBox,false,0);
      siPMPlacement2 = new G4PVPlacement(0,G4ThreeVector(), shadowShell_log,"shadow",fWLBox,false,0);
      break;

    case 10:
      i = 0;
      j = 0;
      while(j<7){
        i=0;
        while(i<8){
          *translation = G4ThreeVector(0,21*mm,offset+j*10*mm);
          translation->rotateZ(i*45*deg);
          dimpleTube = new G4SubtractionSolid("dimpleTube",dimpleTube,smallDimpleSphere, dimpleRot, *translation);
          i++;
        }
      j++;
      }

      fLBox = new G4LogicalVolume(dimpleTube,fTargetMaterial,"tube_log");
      fPBox = new G4PVPlacement(0,G4ThreeVector(0,0,0), fLBox,"target", fWLBox,false, 0);
      siPMPlacement1 = new G4PVPlacement(0,G4ThreeVector(), detShell_log,"detector",fWLBox,false,0);
      break;
  }
  return fWPBox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
