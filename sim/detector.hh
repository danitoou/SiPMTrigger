#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4String.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include <cmath>

class MySensitiveDetector : public G4VSensitiveDetector {

public:
    MySensitiveDetector(G4String);
    ~MySensitiveDetector();

private:
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
    G4double entryEnergy[4];
    G4int preStepCopyNo[4];
    bool prevEvtNum;

};


#endif