#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* myEvent)
{
	fDetectedPhotons = 0;
	fDepositedEnergy = 0;
	fProducedPhotons = 0;
	fEscapedPhoton = 0;
	if (myEvent->GetEventID() % 1000 == 0)
		G4cout << "event no.: " << myEvent->GetEventID() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddWavelength(G4double wavelength){
	auto analysisManager = G4AnalysisManager::Instance();

	analysisManager->FillH1(6, wavelength);
}

void EventAction::EndOfEventAction(const G4Event* myEvent)
{

	auto analysisManager = G4AnalysisManager::Instance();
	if (fDetectedPhotons > 0)
		analysisManager->FillH1(0,fDetectedPhotons);
	if (fDepositedEnergy > 0)
		analysisManager->FillH1(1,fDepositedEnergy);
	if (fProducedPhotons > 0)
		analysisManager->FillH1(2,fProducedPhotons);
	if (fDepositedEnergy>0 && fDetectedPhotons>0) {
		analysisManager->FillH1(3,fDetectedPhotons/fDepositedEnergy);
	}
	if (fDetectedPhotons > 0 && fProducedPhotons >0){
		analysisManager->FillH1(4,fDetectedPhotons/fProducedPhotons);
	}
	if (fEscapedPhoton >0 ){
		analysisManager->FillH1(5,fEscapedPhoton);
	}
}
