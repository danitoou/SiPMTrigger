#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction() {}

MyDetectorConstruction::~MyDetectorConstruction() {}


G4VPhysicalVolume *MyDetectorConstruction::Construct() {
    
    G4NistManager *nist = G4NistManager::Instance();


    H = nist->FindOrBuildElement("H");

    Vacuum = new G4Material("Vacuum", pow(10, -25)*g/cm3, 1);
    Vacuum->AddElement(H, 100*perCent);

    G4double xWorld = 10; // meters
    G4double yWorld = 10;
    G4double zWorld = 10;


    solidWorld = new G4Box("solidWorld", xWorld*m, yWorld*m, zWorld*m);

    logicWorld = new G4LogicalVolume(solidWorld, Vacuum, "logicWorld");

    physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);


    solidDetector = new G4Box("solidDetector", 2*cm, 0.5*cm, 0.5*cm);

    logicDetector = new G4LogicalVolume(solidDetector, Vacuum, "logicDetector");


    G4RotationMatrix *rot = new G4RotationMatrix();

    rot->rotateZ(90*deg);

    physDetector = new G4PVPlacement(0, G4ThreeVector(0*cm, -(yWorld-0.005)*m, 0.5*cm), logicDetector, "physDetector", logicWorld, false, 0, true);
    physDetector = new G4PVPlacement(0, G4ThreeVector(0*cm, -(yWorld-0.015)*m, 0.5*cm), logicDetector, "physDetector", logicWorld, false, 1, true);
    physDetector = new G4PVPlacement(0, G4ThreeVector(0*cm, -(yWorld-0.025)*m, 0.5*cm), logicDetector, "physDetector", logicWorld, false, 2, true);
    physDetector = new G4PVPlacement(0, G4ThreeVector(0*cm, -(yWorld-0.035)*m, 0.5*cm), logicDetector, "physDetector", logicWorld, false, 3, true);
    physDetector = new G4PVPlacement(rot, G4ThreeVector(1.5*cm, -(yWorld-0.02)*m, -0.5*cm), logicDetector, "physDetector", logicWorld, false, 4, true);
    physDetector = new G4PVPlacement(rot, G4ThreeVector(0.5*cm, -(yWorld-0.02)*m, -0.5*cm), logicDetector, "physDetector", logicWorld, false, 5, true);
    physDetector = new G4PVPlacement(rot, G4ThreeVector(-0.5*cm, -(yWorld-0.02)*m, -0.5*cm), logicDetector, "physDetector", logicWorld, false, 6, true);
    physDetector = new G4PVPlacement(rot, G4ThreeVector(-1.5*cm, -(yWorld-0.02)*m, -0.5*cm), logicDetector, "physDetector", logicWorld, false, 7, true);

    


    return physWorld;
    
}

void MyDetectorConstruction::ConstructSDandField() {

    MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

    logicDetector->SetSensitiveDetector(sensDet);
}



