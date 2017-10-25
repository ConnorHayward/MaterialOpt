#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction* );
    virtual ~PrimaryGeneratorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    PrimaryGeneratorAction* 		fPENAction;
    G4UIdirectory*                  fGunDir;
    //G4UIcmdWithADoubleAndUnit*      fPolarCmd;
    G4UIcmdWithAnInteger*			fSourceType;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
