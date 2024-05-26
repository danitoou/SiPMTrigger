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

    man->CreateNtuple("EnergyI", "EnergyI");
    man->CreateNtupleIColumn("Ch0");
    man->CreateNtupleIColumn("Ch1");
    man->CreateNtupleIColumn("Ch2");
    man->CreateNtupleIColumn("Ch3");
    man->CreateNtupleIColumn("Ch4");
    man->CreateNtupleIColumn("Ch5");
    man->CreateNtupleIColumn("Ch6");
    man->CreateNtupleIColumn("Ch7");
    man->FinishNtuple(1);

    man->CreateNtuple("EnergyD", "EnergyD");
    man->CreateNtupleDColumn("Ch0");
    man->CreateNtupleDColumn("Ch1");
    man->CreateNtupleDColumn("Ch2");
    man->CreateNtupleDColumn("Ch3");
    man->CreateNtupleDColumn("Ch4");
    man->CreateNtupleDColumn("Ch5");
    man->CreateNtupleDColumn("Ch6");
    man->CreateNtupleDColumn("Ch7");
    man->FinishNtuple(2);

    man->CreateNtuple("EnergyF", "EnergyF");
    man->CreateNtupleFColumn("Ch0");
    man->CreateNtupleFColumn("Ch1");
    man->CreateNtupleFColumn("Ch2");
    man->CreateNtupleFColumn("Ch3");
    man->CreateNtupleFColumn("Ch4");
    man->CreateNtupleFColumn("Ch5");
    man->CreateNtupleFColumn("Ch6");
    man->CreateNtupleFColumn("Ch7");
    man->FinishNtuple(3);

    // man->CreateNtuple("Distance", "Distance");
    // man->CreateNtupleIColumn("Ch0");
    // man->CreateNtupleIColumn("Ch1");
    // man->CreateNtupleIColumn("Ch2");
    // man->CreateNtupleIColumn("Ch3");
    // man->CreateNtupleIColumn("Ch4");
    // man->CreateNtupleIColumn("Ch5");
    // man->CreateNtupleIColumn("Ch6");
    // man->CreateNtupleIColumn("Ch7");
    // man->FinishNtuple(2);
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