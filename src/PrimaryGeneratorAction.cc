#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4DecayTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0)
{
	G4int n_particle = 1;
	fParticleGun = new G4ParticleGun(n_particle);

	//create a messenger for this class
  	fGunMessenger = new PrimaryGeneratorMessenger(this);

 	 //default kinematic
 	 //
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("e-");

	fSourceType = 3;
	fPhotonWavelength = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	// Particle Types
	// 	0 - perpendicular mono-energetic electrons with fixed position
	//  1 - 207Bi source
	// 	2 - 90Sr source with isotropic continuous electron spectrum
	//  3 - cosmic muons, perpendicular, homogeneously distributed
	//  4 - monochromatic light
	//	5 - UV LED source
	// 	6 - Xe lamp

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.*keV;


	switch (fSourceType) {
		case 0:
			fParticleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
			fParticleGun->SetParticleEnergy(3*MeV);
			fParticleGun->SetParticlePosition(G4ThreeVector(1*cm,0*cm,0*cm));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-1.,0.,0));
			break;
		case 1:
			fParticleGun->SetParticleEnergy(0*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-1.0*cm));
		//	fParticleGun->SetParticlePosition(G4ThreeVector(30*cm,0.,0));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
			fParticleGun->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(83,207,excitEnergy));
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 2:
			fParticleGun->SetParticleEnergy(0*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-1.0*cm));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
			fParticleGun->SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(38,90,excitEnergy));
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 3:
			fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu-"));
			break;
		case 4:
			if (fPhotonWavelength == 0)
			{
				G4cerr << "Default value (420 nm) of primary photons was set!" << G4endl;
				fPhotonWavelength = 420;
			}
			fParticleGun->SetParticleDefinition(particleTable->FindParticle("opticalphoton"));
			fParticleGun->SetParticleEnergy(1.24/(fPhotonWavelength/1000));
			break;
		case 5:
			fParticleGun->SetParticleDefinition(particleTable->FindParticle("opticalphoton"));
			break;
		case 6:
			fParticleGun->SetParticleDefinition(particleTable->FindParticle("opticalphoton"));
			break;
		default:
			break;
	}


	// fParticleGun->SetParticlePosition(G4ThreeVector(0.0*cm,0.0*cm,0.0*cm));
	// fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
	// fParticleGun->SetParticleEnergy(500.0*keV);

	fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetSourceType(G4int newType)
{
	if (newType < 6 && newType > 0)
	{
		fSourceType = newType;
	}
	else
	{
		G4cerr << "The option is out of the possible values (0-5)!" << G4endl;
		G4cerr << "The default option (0) is set!" << G4endl;
		fSourceType = 0;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetPhotonWavelength(G4double newValue)
{
	if (newValue > 200 && newValue < 700)
	{
		fPhotonWavelength = newValue;
	}
	else
	{
		G4cerr << "The new desired wavelength is out of range (200-700 nm)!" << G4endl;
		G4cerr << "The photon wavelength is set to default value (420 nm)!" << G4endl;
		fPhotonWavelength = 420;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
