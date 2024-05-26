#include "detector.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name){}

MySensitiveDetector::~MySensitiveDetector(){}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist) {

    G4AnalysisManager *man = G4AnalysisManager::Instance();

    if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {
        const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
        G4int copyNo = touchable->GetCopyNumber();

        G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

        if(prevEvtNum != evt*10 + copyNo) {

        
            man->FillNtupleIColumn(0, copyNo, evt);
            man->AddNtupleRow(0);
            prevEvtNum = evt*10 + copyNo;

            G4ThreeVector preStepPointPos = aStep->GetPreStepPoint()->GetPosition();

            // G4cout << "---------------PreStepPoint---------------" << G4endl;
            // G4cout << "Event ID: " << evt << G4endl;
            // G4cout << "Detector ID: " << copyNo << G4endl;
            // G4cout << "X: " << preStepPointPos.getX() << G4endl;
            // G4cout << "Y: " << preStepPointPos.getY() << G4endl;
            // G4cout << "Z: " << preStepPointPos.getZ() << G4endl;
            

            G4int thr = G4Threading::G4GetThreadId();

            entryEnergy[thr] = aStep->GetPreStepPoint()->GetKineticEnergy();
            preStepCopyNo[thr] = copyNo;
        } else {
            G4cout << "Skipped: " << prevEvtNum << ", " << evt << ", " << copyNo << G4endl;
        }
    }

    if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {

        G4int thr = G4Threading::G4GetThreadId();
        const G4VTouchable *touchable = aStep->GetPostStepPoint()->GetTouchable();
        G4int copyNo = touchable->GetCopyNumber();

        G4double energyDep = entryEnergy[thr] - aStep->GetPostStepPoint()->GetKineticEnergy();
        G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

        G4ThreeVector postStepPointPos = aStep->GetPostStepPoint()->GetPosition();

        // G4cout << "---------------PostStepPoint---------------" << G4endl;
        // G4cout << "Event ID: " << evt << G4endl;
        // G4cout << "Detector ID: " << preStepCopyNo[thr] << G4endl;
        // G4cout << "Energy Deposited: " << energyDep << G4endl;
        // G4cout << "X: " << postStepPointPos.getX() << G4endl;
        // G4cout << "Y: " << postStepPointPos.getY() << G4endl;
        // G4cout << "Z: " << postStepPointPos.getZ() << G4endl;

        bool breakCheck = false;

        if(energyDep > 10000) {
            G4cout << "Not counted" << G4endl;
            breakCheck = true;
        }

        if(!breakCheck) {
            man->FillNtupleIColumn(1, preStepCopyNo[thr], energyDep);
            man->AddNtupleRow(1);
            man->FillNtupleDColumn(2, preStepCopyNo[thr], energyDep);
            man->AddNtupleRow(2);
            man->FillNtupleFColumn(3, preStepCopyNo[thr], energyDep);
            man->AddNtupleRow(3);
        }

    }

    return true;

    // G4Track *track = aStep->GetTrack();



    // G4double energyDep = aStep->GetTotalEnergyDeposit()*100000000;

    // G4ThreeVector preStepPointPos = aStep->GetPreStepPoint()->GetPosition();
    // G4ThreeVector postStepPointPos = aStep->GetPostStepPoint()->GetPosition();

    // // G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    // G4double dist = 1000*sqrt(pow(preStepPointPos.getX() - postStepPointPos.getX(), 2) + pow(preStepPointPos.getY() - postStepPointPos.getY(), 2) + pow(preStepPointPos.getZ() - postStepPointPos.getZ(), 2));

    // // G4cout << copyNo << " " << evt << G4endl;
    // // G4cout << "Pre Step Point:" << preStepPointPos.getX() << ", " << preStepPointPos.getY() << ", " << preStepPointPos.getZ() << ", " << G4endl;
    // // G4cout << "Post Step Point:" << postStepPointPos.getX() << ", " << postStepPointPos.getY() << ", " << postStepPointPos.getZ() << ", " << G4endl;
    // // G4cout << "Difference:" << preStepPointPos.getX() - postStepPointPos.getX() << ", " << preStepPointPos.getY() - postStepPointPos.getY() << ", " <<  preStepPointPos.getZ() - postStepPointPos.getZ() << ", " << G4endl;
    


    
    // // man->FillNtupleDColumn(copyNo, 1, posDetector[0]);
    // // man->FillNtupleDColumn(copyNo, 2, posDetector[1]);
    // // man->FillNtupleDColumn(copyNo, 3, posDetector[2]);
    
    // man->FillNtupleIColumn(1, copyNo, energyDep);
    // man->AddNtupleRow(1);
    // man->FillNtupleIColumn(2, copyNo, dist);
    // man->AddNtupleRow(2);

    // return true;

}