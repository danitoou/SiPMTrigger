#ifndef CONSTRCUTION_HH
#define CONSTRUCTION_HH

#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4LogicalVolume.hh"
#include "G4GenericMessenger.hh"

#include "detector.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction {

public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();

    virtual G4VPhysicalVolume *Construct();
        

private:
    G4LogicalVolume *logicDetector;
    virtual void ConstructSDandField();

    G4Box *solidWorld, *solidDetector;
    G4LogicalVolume *logicWorld;
    G4VPhysicalVolume *physWorld, *physDetector;

    G4GenericMessenger *fMessenger;

    G4Material *Vacuum, *scintMat;
    G4Element *H;

    void DefineMaterials();
};


#endif
