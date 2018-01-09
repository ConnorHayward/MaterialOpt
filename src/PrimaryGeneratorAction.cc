#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"
#include "G4DecayTable.hh"
#include "G4OpticalPhoton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0),
fDetector(det)
{
	G4int n_particle = 1;
	fPrimaryMessenger = new PrimaryGeneratorMessenger(this);
	fParticleGun = new G4ParticleGun(n_particle);

 	//default kinematic
 	//
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("e-");
	fSourceType = 0;
	fSourceEnergy = 1500*keV;
	fPhotonWavelength = 0;
	fParticleName = "void";
	DefineParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
}

void PrimaryGeneratorAction::DefineParticle(){

	// Particle Types
	// 	0 - perpendicular mono-energetic electrons with fixed position
	//  1 - 60Co source
	// 	2 - 137Cs source
	//	3 - 90Sr source
	//	4 - 241Am source
	//	5 - 106Ru source
	//  6 - 42Ar
	//	7 - 100 keV point gamma, for position scans.
	//	8 - 430 nm optical photon source, direction randomized.
	//	9 - 207Bi source
	//	10 - Very low energy gamma, Not used.

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.*keV;
		G4int Z=0, A=0;
		G4double x,y,z, sum;
		G4ParticleDefinition* ion;
		double position = 1*cm+(fDetector->GetTargetSize())/2;
	switch (fSourceType) {
		case 0:
			fParticleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
			fParticleGun->SetParticleEnergy(fSourceEnergy);
			fParticleGun->SetParticlePosition(G4ThreeVector(0*cm,0*mm,-15*mm));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2*(rand()-1),2*(rand()-1),2*(rand()-1)));
			break;
		case 1:
			Z = 27;
			A = 60;
			ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
			fParticleGun->SetParticleEnergy(0*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-position));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1));
			fParticleGun->SetParticleDefinition(ion);
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 2:
			Z = 55;
			A = 137;
			ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
			fParticleGun->SetParticleEnergy(0*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-position));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1));
			fParticleGun->SetParticleDefinition(ion);
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 3:
			Z=38;
			A=90;
			ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
			fParticleGun->SetParticleEnergy(0*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-15*mm));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1));
			fParticleGun->SetParticleDefinition(ion);
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 4:
			Z=95;
			A=241;
			ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
			fParticleGun->SetParticleEnergy(0*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-position));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1));
			fParticleGun->SetParticleDefinition(ion);
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 5:
			Z=44;
			A=106;
			ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
			fParticleGun->SetParticleEnergy(0*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-position));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1));
			fParticleGun->SetParticleDefinition(ion);
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 6:
			Z=8;
			A=42;
			ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
			fParticleGun->SetParticleEnergy(0*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-position));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1));
			fParticleGun->SetParticleDefinition(ion);
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 7:
			fParticleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
			fParticleGun->SetParticleEnergy(100*keV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0*cm,-1.25*cm,0));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,1,0));
			break;
		case 8:
			fParticleGun->SetParticleDefinition(G4OpticalPhoton::Definition());
			fParticleGun->SetParticlePolarization(G4ThreeVector(2*(rand()-1),2*(rand()-1),2*(rand()-1)));
			x = 2*(rand()-1);
			y= 2*(rand()-1);
			z=2*(rand()-1);
			sum = 1 / sqrt(x*x + y*y +z*z);
			fParticleGun->SetParticleEnergy(2.8*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(325*mm*x*sum,325*mm*y*sum,325*mm*z*sum));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2*(rand()-1),2*(rand()-1),2*(rand()-1)));
			break;
		case 9:
			Z=83;
			A=207;
			ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
			fParticleGun->SetParticleEnergy(0*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-15*mm));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1));
			fParticleGun->SetParticleDefinition(ion);
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 10:
			fParticleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
			fParticleGun->SetParticleEnergy(4.95*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0*cm,32.5*mm,0));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2*(rand()-1),2*(rand()-1),2*(rand()-1)));
			break;
		case 11:
			Z=19;
			A=40;
			ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
			fParticleGun->SetParticleEnergy(0*eV);
			fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1));
			fParticleGun->SetParticleDefinition(ion);
			fParticleGun->SetParticleCharge(ionCharge);
			break;
		case 12:
			fParticleGun->SetParticleDefinition(G4OpticalPhoton::Definition());
			fParticleGun->SetParticlePolarization(G4ThreeVector(2*(rand()-1),2*(rand()-1),2*(rand()-1)));
			x = 2*(rand()-1);
			y= 2*(rand()-1);
			z=2*(rand()-1);
			sum = 1 / sqrt(x*x + y*y);
			fParticleGun->SetParticleEnergy(2.8*eV);
			//fParticleGun->SetParticlePosition(G4ThreeVector(16.5*mm*x*sum,16.5*mm*y*sum,40*mm*z));
			fParticleGun->SetParticlePosition(G4ThreeVector(x*sum*16.25*mm,y*sum*16.25*mm,0));
			fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2*(rand()-1),2*(rand()-1),2*(rand()-1)));
			break;
		}
	SetParticleName(Z,A,excitEnergy);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Randomises placement and momentum vectors for required sources.

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	double x,y,z,sum;
	if(fSourceType==0){
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2*(rand()-1),2*(rand()-1),2*(rand()-1)));}
	else if(fSourceType==8){
		x = 2*(rand()-1);
		y= 2*(rand()-1);
		z=2*(rand()-1);
		sum = 1 / sqrt(x*x + y*y +z*z);
		fParticleGun->SetParticleEnergy(2.8*eV);
		fParticleGun->SetParticlePosition(G4ThreeVector(330*mm*x*sum,330*mm*y*sum,330*mm*z*sum));
		fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2*(rand()-1),2*(rand()-1),2*(rand()-1)));
	}
	else if (fSourceType==12){
		x = 2*(rand()-1);
		y= 2*(rand()-1);
		sum = 1 / sqrt(x*x + y*y);
		fParticleGun->SetParticleEnergy(2.8*eV);
		//fParticleGun->SetParticlePosition(G4ThreeVector(16.5*mm*x*sum,16.5*mm*y*sum,40*mm*z));
		fParticleGun->SetParticlePosition(G4ThreeVector(x*sum*16.25*mm,y*sum*16.25*mm,0));
		fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2*(rand()-1),2*(rand()-1),2*(rand()-1)));
	}

	fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticleName(G4int Z, G4int A, G4double excitEnergy)
{
	fParticleName = G4IonTable::GetIonTable()->GetIonName(Z,A,excitEnergy);
}

void PrimaryGeneratorAction::SetPosition(G4double newPosition){
	fPosition = newPosition;
	DefineParticle();
}

void PrimaryGeneratorAction::SetSourceType(G4int newType)
{
	if (newType <= 12 && newType >= 0)
	{
		fSourceType = newType;
	}
	else
	{
		G4cerr << "The option is out of the possible values (0-5)!" << G4endl;
		G4cerr << "The default option (0) is set!" << G4endl;
		fSourceType = 0;
	}
	DefineParticle();
}

void PrimaryGeneratorAction::SetSourceEnergy(G4double newEnergy)
{
	if (newEnergy>0)
	{
		fSourceEnergy = newEnergy;
	}
	else{
		G4cerr << "New energy is < 0." << G4endl;
		G4cerr << "The default option 60 keV is set!" << G4endl;
		fSourceEnergy = 60*keV;
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
