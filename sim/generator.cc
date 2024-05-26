#include "generator.hh"

MyPrimaryGenerator::MyPrimaryGenerator() {

    fParticleGun = new G4ParticleGun(1);

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

    G4String particleName = "mu-";

    G4ParticleDefinition *particle = particleTable->FindParticle(particleName);

    xWorld = 10;
    yWorld = 10;
    zWorld = 10;

    phi = G4UniformRand()*2*M_PI;
    theta = d(gen);

    randX = xWorld*(G4UniformRand() - G4UniformRand())*m;
    randZ = zWorld*(G4UniformRand() - G4UniformRand())*m;
    // randX = 0*m;
    // randZ = 0*m;

    dist = pow(randX, 2) + pow(randZ, 2);
    randPos = G4ThreeVector(randX, yWorld*m, randZ);
    tg = tan(theta);

    v = G4ThreeVector(tg*sin(phi), -1, tg*cos(phi));


    fParticleGun->SetParticleMomentumDirection(v);
    fParticleGun->SetParticlePosition(randPos);
    fParticleGun->SetParticleMomentum(100.*GeV);
    fParticleGun->SetParticleDefinition(particle);
}

MyPrimaryGenerator::~MyPrimaryGenerator() {
    delete fParticleGun;

}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent) {

    if(counter % (1000000000/100) == 0 ) {
        G4cout << counter*100/1000000000 << "%" << G4endl;
    }
    counter++;

    
    phi = G4UniformRand()*2*M_PI;
    theta = d(gen);

    randX = xWorld*(G4UniformRand() - G4UniformRand())*m;
    randZ = zWorld*(G4UniformRand() - G4UniformRand())*m;
    // randX = 0*m;
    // randZ = 0*m;


    dist = pow(randX, 2) + pow(randZ, 2);
    randPos = G4ThreeVector(randX, yWorld*m, randZ);
    tg = tan(theta);

    v = G4ThreeVector(tg*sin(phi), -1, tg*cos(phi));


    fParticleGun->SetParticleMomentumDirection(v);
    fParticleGun->SetParticlePosition(randPos);
    fParticleGun->GeneratePrimaryVertex(anEvent);


    
}