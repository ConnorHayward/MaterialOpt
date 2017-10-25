#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	PrimaryGeneratorAction();
	virtual ~PrimaryGeneratorAction();

public:
	virtual void GeneratePrimaries(G4Event*);

	void  SetSourceType(G4int newType);
	G4int GetSourceType(void){return fSourceType;};
	void SetPhotonWavelength(G4double newValue);

private:
	G4ParticleGun*  fParticleGun;
	G4int           fSourceType;
	G4double		fPhotonWavelength;
	PrimaryGeneratorMessenger* fGunMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*PrimaryGeneratorAction_h*/
