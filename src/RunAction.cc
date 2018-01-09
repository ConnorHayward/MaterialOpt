// Make this appear first!
#include "G4Timer.hh"

#include "RunAction.hh"

#include "G4Run.hh"

#include "Analysis.hh"
//#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//#include "Run.hh"

#include "EventAction.hh"

#include <string>

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
 : G4UserRunAction(),fDetector(det),fPrimary(kin),
   fTimer(0)//fRun(0),//fHistoManager(0)
{
  fTimer = new G4Timer;
  fMan = G4AnalysisManager::Instance();
  //fHistoManager = new HistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
//  delete fHistoManager
  delete fTimer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  fMan = G4AnalysisManager::Instance();

  G4String targetName = fDetector->GetTargetMaterialName();
  G4String targetThickness = G4BestUnit(fDetector->GetTargetSize(),"Length");
  G4int sourceName = fPrimary->GetSourceType();
  G4String sourceString;

  if(sourceName == 0 || sourceName==7){
    G4String sourceEnergy = G4BestUnit(fPrimary->GetSourceEnergy(),"Energy");
    sourceString = "gamma-"+ sourceEnergy;
  }
  else{
    sourceString = fPrimary->GetParticleName();
  }

  G4int directory = fDetector->GetDetectorType();
  G4String directorName;

  G4String position = G4BestUnit((fDetector->GetTargetSize()-fPrimary->GetPosition()),"Length");

  //fFileName = "Data/"+targetName + "-" + targetThickness+"-"+sourceString+"-"+position+"-1mm";
  fFileName = "dimpleCylinder";
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fTimer->Start();

  fMan->SetVerboseLevel(0);
  fMan->OpenFile(fFileName);

  fMan->CreateH1("Light Output","N Photons Detected",250,0,100000);
  fMan->CreateH1("Energy in Target","Deposited energy target [MeV]",500,0,2);
  fMan->CreateH1("Light Yield","N Photons Produced",5000,0,10000);
  fMan->CreateH1("Yield/Energy Dep [N/MeV]", "N Photons / E Dep[MeV]",1000,0,15000);
  fMan->CreateH1("Detected/Produced","Detected/Produced", 100,0,1);
  fMan->CreateH1("N Escaped Photons", "N Escaped Photons",1000,0,10000);
  fMan->CreateH1("Detected Wavelength [nm]", "Wavelgnth of detected photon [nm]",50,400,500);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  fTimer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent()
         << " " << *fTimer << G4endl;
    G4cout << "End Run" << G4endl;

    fMan->Write();
    fMan->CloseFile();
    delete fMan;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
