#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4RandomDirection.hh"
#include <random>


class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
    MyPrimaryGenerator();
    ~MyPrimaryGenerator();

    virtual void GeneratePrimaries(G4Event*);

private:
    G4ParticleGun *fParticleGun;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> d{0., M_PI/5};
    G4double phi, theta, randX, randZ, dist, tg;
    G4double xWorld, yWorld, zWorld;
    G4ThreeVector randPos, v;
    G4int counter;

};



#endif