#include "run.hh"

MyRunAction::MyRunAction() {
    G4AnalysisManager *man = G4AnalysisManager::Instance();

    man->CreateNtuple("Hits", "Hits");
    man->CreateNtupleIColumn("Ch0");
    man->CreateNtupleIColumn("Ch1");
    man->CreateNtupleIColumn("Ch2");
    man->CreateNtupleIColumn("Ch3");
    man->CreateNtupleIColumn("Ch4");
    man->CreateNtupleIColumn("Ch5");
    man->CreateNtupleIColumn("Ch6");
    man->CreateNtupleIColumn("Ch7");
    man->FinishNtuple(0);
}

MyRunAction::~MyRunAction(){}

void MyRunAction::BeginOfRunAction(const G4Run* run) {

    G4AnalysisManager *man = G4AnalysisManager::Instance();

    // G4String file = "output.root";

    G4int runNumber = run->GetRunID();

    std::stringstream strRunID;
    strRunID << runNumber;

    man->OpenFile("output" + strRunID.str() + ".root");


}


void MyRunAction::EndOfRunAction(const G4Run*) {

    G4AnalysisManager *man = G4AnalysisManager::Instance();

    // G4String file = "output.root";

    man->Write();
    man->CloseFile();

}