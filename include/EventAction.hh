#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class RunAction;

/// Event action class
///

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    void    AddDetectedPhoton(void){fDetectedPhotons++;};
    G4double   GetNumberDetectedPhotons(void){return fDetectedPhotons;};
    void    AddDepositedEnergy(G4double newEnergy){fDepositedEnergy += newEnergy;};
    void    AddDepositedEnergyA(G4double newEnergy){fDepositedEnergyA += newEnergy;};
    void    AddWavelength(G4double newWavelength);
    void    AddIWavelength(G4double startWavelength);
    void    AddBAWavelength(G4double boundaryAbsorbedWavelength);
    void    AddProducedPhoton(void) {fProducedPhotons++;};
    G4double   GetNumerProducedPhotons(void){return fProducedPhotons;};
    void AddEscapedPhoton(void){fEscapedPhoton++;};
    G4double GetNumberEscapedPhoton(void){return fEscapedPhoton;};

  private:
    RunAction* 	fRunAction;
    G4double fEscapedPhoton;
    G4double fProducedPhotons;
    G4double fDetectedPhotons;
    G4double        fDepositedEnergy;
    G4double        fDepositedEnergyA;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
