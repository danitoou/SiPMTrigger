#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction() {}

MyDetectorConstruction::~MyDetectorConstruction() {}


G4VPhysicalVolume *MyDetectorConstruction::Construct() {
    
    G4NistManager *nist = G4NistManager::Instance();


    H = nist->FindOrBuildElement("H");

    scintMat = nist->FindOrBuildMaterial("G4_POLYSTYRENE");

    Vacuum = new G4Material("Vacuum", pow(10, -25)*g/cm3, 1);
    Vacuum->AddElement(H, 100*perCent);

    G4double xWorld = 10; // meters
    G4double yWorld = 10;
    G4double zWorld = 10;


    solidWorld = new G4Box("solidWorld", xWorld*m, yWorld*m, zWorld*m);

    logicWorld = new G4LogicalVolume(solidWorld, Vacuum, "logicWorld");

    physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);


    solidDetector = new G4Box("solidDetector", 2*cm, 0.5*cm, 0.5*cm);

    logicDetector = new G4LogicalVolume(solidDetector, scintMat, "logicDetector");


    G4RotationMatrix *rot = new G4RotationMatrix();

    rot->rotateZ(90*deg);

    for(int i = 0; i < 4; i++) {
        physDetector = new G4PVPlacement(0, G4ThreeVector(0*cm, -(yWorld-0.01*(i+0.5))*m + 0.01*cm, 0.5*cm), logicDetector, "physDetector", logicWorld, false, i, true);
        physDetector = new G4PVPlacement(rot, G4ThreeVector((1.5-i)*cm, -(yWorld-0.02)*m + 0.01*cm, -0.5*cm), logicDetector, "physDetector", logicWorld, false, i+4, true);
    }

    
    return physWorld;
    
}

void MyDetectorConstruction::ConstructSDandField() {

    MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

    logicDetector->SetSensitiveDetector(sensDet);
}



