#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det),
 fPENDir(0),
 fDetDir(0),
 fTMaterCMD(0),
 fWMaterCMD(0),
 fSizeCMD(0),
 fTypeCMD(0)
{
  fDetDir = new G4UIdirectory("/PEN/det/");
  fDetDir->SetGuidance("detector construction commands");

  fTMaterCMD = new G4UIcmdWithAString("/PEN/det/setTargetMat",this);
  fTMaterCMD->SetGuidance("Select material of the target.");
  fTMaterCMD->SetParameterName("choice",false);
  fTMaterCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTMaterCMD->SetToBeBroadcasted(false);

  fWMaterCMD = new G4UIcmdWithAString("/PEN/det/setWorldMat",this);
  fWMaterCMD->SetGuidance("Select material of the world.");
  fWMaterCMD->SetParameterName("choice",false);
  fWMaterCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWMaterCMD->SetToBeBroadcasted(false);

  fSizeCMD = new G4UIcmdWithADoubleAndUnit("/PEN/det/setSize",this);
  fSizeCMD->SetGuidance("Set size of the box");
  fSizeCMD->SetParameterName("Size",false);
  fSizeCMD->SetRange("Size>0.");
  fSizeCMD->SetUnitCategory("Length");
  fSizeCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSizeCMD->SetToBeBroadcasted(false);

  fTypeCMD = new G4UIcmdWithAnInteger("/PEN/det/setDetectorType",this);
  fTypeCMD->SetGuidance("Set detector type");
  fTypeCMD->SetRange("Type 0-1");
  fTypeCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTypeCMD->SetToBeBroadcasted(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fTMaterCMD;
  delete fWMaterCMD;
  delete fSizeCMD;
  delete fTypeCMD;
  delete fDetDir;
  delete fPENDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
   if( command == fTMaterCMD )
    { fDetector->SetTargetMaterial(newValue);}

   if( command == fWMaterCMD )
    { fDetector->SetWorldMaterial(newValue);}

   if( command == fSizeCMD )
    { fDetector->SetSize(fSizeCMD->GetNewDoubleValue(newValue));}

  if( command == fTypeCMD )
    { fDetector->SetDetectorType(fTypeCMD->GetNewIntValue(newValue));}

}
