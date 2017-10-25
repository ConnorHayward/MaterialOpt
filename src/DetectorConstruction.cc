#include "DetectorConstruction.hh"

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

#include <iostream>
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
  fExpHall_x = fExpHall_y = fExpHall_z = 0.5*m;
  fTank_x    = fTank_y    = fTank_z    =  0.5*m;
  fBubble_x  = fBubble_y  = fBubble_z  =  0.5*m;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

// ------------- Materials -------------

  G4double a, z, density;
  G4int nelements;

// Air
//
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  G4Material* air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
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

  //
  // PEN
  //

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

  G4int absEntries = 0;
  ifstream ReadAbs;

  G4int penType = 1;
  G4String abs_file;

  /*
  Switch for type of PEN
  0 - PEN8050 from TU Dortmund (2016)
  1 - Experiment 4 from TU Dortmund
  2 - Constant ABS with absorption length of 2m
  3 - Constant ABS with absorption of 1.4cm
  */
  switch (penType){

    case 0:
    abs_file = "../input_files/OldPen.csv";
    break;

    case 1:
    abs_file = "../input_files/Exp4.csv";
    break;

    case 2:
    abs_file = "../input_files/highAbs.csv";
    break;

    case 3:
    abs_file = "../input_files/flatAbs.csv";
    break;
  }

  ReadAbs.open(abs_file);

  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>varabsorlength>>filler>>ems>>filler>>rindex;
      if(ReadAbs.eof()){
        break;
      }
      //cout<<absEntries << " " <<wavelength << " " << varabsorlength << " " << ems <<endl;
      absEnergy[absEntries] = (1240/wavelength)*eV;
      abs[absEntries] = varabsorlength*mm;
      emission[absEntries] = ems;
      rIndex[absEntries] = rindex;
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

  //
  // Air
  //
    G4double refractiveIndex2[] =
            { 1.00,1.00,1.00,1.00,1.00,1.00 };
        //  { 1.50,1.50,1.50,1.50,1.50,1.50 };
    G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
    myMPT2->AddProperty("RINDEX", absEnergy, refractiveIndex2, 6);

    air->SetMaterialPropertiesTable(myMPT2);

//
// ------------- Volumes --------------

// The experimental Hall
//
  G4Box* expHall_box = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,air,"World",0,0,0);

  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

// ------------- Scintillator Volumes --------------

G4double startAngle_pmt = 0.*deg;
G4double spanningAngle_pmt = 360.*deg;
// PEN Tile

  G4Box* penTile_box = new G4Box("scintillator", 17.5*mm,17.5*mm, 1.5*mm);

  G4LogicalVolume* penTile_log = new G4LogicalVolume(penTile_box,PEN, "scintillator",0,0,0);

// Reflective wrap for Tile

  G4Box* refBox = new G4Box("reflector",18*mm,18.5*mm, 3*mm);
//  G4SubtractionSolid* wrapBox = new G4SubtractionSolid("teflon_tape", refBox, penTile_box);
  G4LogicalVolume* refBox_log = new G4LogicalVolume(refBox, teflon,"reflector",0,0,0);

// SiPM
  G4Box* siPM_box = new G4Box("siPM_box", 0.3*cm, 0.3*cm, 0.25*mm);
  G4LogicalVolume* siPM_log = new G4LogicalVolume(siPM_box, Si, "siPM_log");

  G4double siPM_energy[] = {2.479684*eV, 2.610194*eV, 2.755204*eV, 2.883353*eV , 2.917275*eV, 3.099605*eV};
  G4double siPM_EFF[]={0.47, 0.48, 0.5, 0.48,0.47,0,45};
  assert(sizeof(siPM_EFF)==sizeof(siPM_energy));
  G4double siPM_REFL[] = {0.,0.,0.,0.,0.,0.};

  G4OpticalSurface* siPM_optsurf = new G4OpticalSurface("siPM_optsurf",glisur, polished, dielectric_metal);
  G4MaterialPropertiesTable* siPM_MT = new G4MaterialPropertiesTable();
  siPM_MT->AddProperty("EFFICIENCY", siPM_energy, siPM_EFF,6);
  siPM_MT->AddProperty("REFLECTIVITY", siPM_energy, siPM_REFL,6);
  siPM_optsurf->SetMaterialPropertiesTable(siPM_MT);
  new G4LogicalSkinSurface("siPM_surf", siPM_log, siPM_optsurf);

  G4Box* absorberPb_box = new G4Box("absorber",50*mm,50*mm, 5*mm);
  G4LogicalVolume* absorberPb_log = new G4LogicalVolume(absorberPb_box,Pb,"absorber",0,0,0);
  G4VPhysicalVolume* absorberPb_phys = new G4PVPlacement(0,G4ThreeVector(),absorberPb_log,"absorber",expHall_log,false,0);

  G4Box* absorberCu_box = new G4Box("foil",50*mm,50*mm,1.5*mm);
  G4LogicalVolume* absorberCu_log = new G4LogicalVolume(absorberCu_box,Cu,"foil",0,0,0);
  G4VPhysicalVolume* absorberCu_phys = new G4PVPlacement(0,G4ThreeVector(0,0,6.5*mm),absorberCu_log,"absorber",absorberPb_log,false,0);



  // Teflon Wrap Surface
  G4MaterialPropertiesTable* wrapSurface_mpt = new G4MaterialPropertiesTable();
  G4double wrapReflectivity[]={0.99,0.99,0.99,0.99,0.99,0.99};
  G4double wrapEff[]={0.1,0.1,0.1,0.1,0.1,0.1};
  G4double specLobe[]={0.,0.,0.,0.,0.,0.};
  G4double backScatter[]={0.,0.,0.,0.,0.,0.};
  G4double specSpike[]={1.0,1.0,1.0,1.0,1.0,1.0};
  G4double RItef[]={1.35,1.35,1.35,1.35,1.35,1.35};
  G4double sigAlpha[]={2*degree,2*degree,2*degree,2*degree,2*degree,2*degree};

  G4OpticalSurface* wrapSurface = new G4OpticalSurface("wrapSurface");
  wrapSurface->SetType(dielectric_metal);
  wrapSurface->SetModel(unified);
  wrapSurface->SetFinish(polished);

  wrapSurface_mpt->AddProperty("RINDEX",siPM_energy,RItef,6);
  wrapSurface_mpt->AddProperty("REFLECTIVITY", siPM_energy, wrapReflectivity, 6);
  wrapSurface_mpt->AddProperty("EFFICIENCY", siPM_energy, wrapEff,6);
  wrapSurface_mpt->AddProperty("SPECULARLOBECONSTANT", siPM_energy, specLobe,6);
  wrapSurface_mpt->AddProperty("SPECULARSPIKECONSTANT",siPM_energy,specSpike,6);
  wrapSurface_mpt->AddProperty("BACKSCATTERCONSTANT",siPM_energy,backScatter,6);
  wrapSurface_mpt->AddProperty("SIGMA_ALPHA", siPM_energy, sigAlpha, 6);

  wrapSurface->SetMaterialPropertiesTable(wrapSurface_mpt);

//-------------- Define Detector Construction --------------
G4int geom = 0;

G4VPhysicalVolume* penTile_phys;
G4RotationMatrix* rot = new G4RotationMatrix();
G4RotationMatrix* rotm;
G4VPhysicalVolume* siPM_phys1;
G4VPhysicalVolume* siPM_phys2;
G4VPhysicalVolume* wrap_phys;
G4SubtractionSolid* refBox1;

/*

  0 - Standard tile coupled to SiPM on base.
  1 - Standard tile coupled to SiPM on side.

*/
/*  switch(geom){

//  Tile coupled to SiPM on base

  case 0:
    penTile_phys = new G4PVPlacement(0,G4ThreeVector(),penTile_log,"scintillator",expHall_log,false,0);
    siPM_phys1 = new G4PVPlacement(0,G4ThreeVector(0,0,2.75*mm),siPM_log,"detector",expHall_log,false,0);
    break;


//  Tile coupled to SiPM on side

  case 1:
    rotm = new G4RotationMatrix();
    rotm->rotateY(90*deg);
    penTile_phys = new G4PVPlacement(rotm,G4ThreeVector(),penTile_log,"scintillator",expHall_log,false,0);
    siPM_phys1 = new G4PVPlacement(0,G4ThreeVector(0,0,17.75*mm),siPM_log,"detector",expHall_log,false,0);
    break;

}*/

  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
