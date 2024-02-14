#include "detector.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name){}

MySensitiveDetector::~MySensitiveDetector(){}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist) {

    G4Track *track = aStep->GetTrack();


    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();

    G4int copyNo = touchable->GetCopyNumber();


    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->FillNtupleIColumn(0, copyNo, evt);
    // man->FillNtupleDColumn(copyNo, 1, posDetector[0]);
    // man->FillNtupleDColumn(copyNo, 2, posDetector[1]);
    // man->FillNtupleDColumn(copyNo, 3, posDetector[2]);
    man->AddNtupleRow(0);

    return true;

}