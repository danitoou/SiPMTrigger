#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>

void readTree() {

    TFile *input = new TFile("output.root", "READ");

    TFile *output = new TFile("histos.root", "RECREATE");
    TCanvas *c1 = new TCanvas();

    TTree *tree = (TTree*)input->Get("Hits");

    int ch0, ch1, ch2, ch3, ch4, ch5, ch6, ch7 = 1;
    int arr[8] = {ch0, ch1, ch2, ch3, ch4, ch5, ch6, ch7};

    tree->SetBranchAddress("Ch0", &ch0);
    tree->SetBranchAddress("Ch1", &ch1);
    tree->SetBranchAddress("Ch2", &ch2);
    tree->SetBranchAddress("Ch3", &ch3);
    tree->SetBranchAddress("Ch4", &ch4);
    tree->SetBranchAddress("Ch5", &ch5);
    tree->SetBranchAddress("Ch6", &ch6);
    tree->SetBranchAddress("Ch7", &ch7);


    int entries = tree->GetEntries();

    TH1F* histos[8];

    char str[20];

    for(int i = 0; i < 8; i++) {
        sprintf(str, "Ch%d_Hits", i);
        histos[i] = new TH1F(str, str, 100, 0, 1000000000);
    }

    FILE * file;

    file = fopen("data.txt", "w+");


    for(int i = 0; i < entries; i++) {

        tree->GetEntry(i);

        // for(int j = 0; j < 8; j++) {
        //     histos[j]->Fill(arr[j]);
        // }

        arr[0] = ch0;
        arr[1] = ch1;
        arr[2] = ch2;
        arr[3] = ch3;
        arr[4] = ch4;
        arr[5] = ch5;
        arr[6] = ch6;
        arr[7] = ch7;

        for(int j = 0; j < 8; j++) {
            if(arr[j] == 0) continue;
            for(int k = j+1; k < 8; k++) {
                // cout << ch0 << " " << ch1 << " " << ch2 << " " << ch3 << " " << ch4 << " " << ch5 << " " << ch6 << " " << ch7 << endl;
                if(arr[k] != 0) {
                    cout << j << "x" << k << " " << arr[j] << " " << arr[k] << endl;
                }
            }
        }

        histos[0]->Fill(ch0);
        histos[1]->Fill(ch1);
        histos[2]->Fill(ch2);
        histos[3]->Fill(ch3);
        histos[4]->Fill(ch4);
        histos[5]->Fill(ch5);
        histos[6]->Fill(ch6);
        histos[7]->Fill(ch7);

        fprintf(file, "%d,%d,%d,%d,%d,%d,%d,%d\n", ch0, ch1, ch2, ch3, ch4, ch5, ch6, ch7);

    }

    fclose(file);

    tree = (TTree*)input->Get("EnergyD");

    char name[20];

    for(int i = 0; i < 8; i++) {
        sprintf(name, "Ch%d_Energy", i);
        TH1F *histo = new TH1F(name, name, 100, -5, 30);
        sprintf(name, "Ch%d>>Ch%d_Energy", i, i);
        tree->Draw(name);
    }

    for(int i = 0; i < 8; i++) {
        histos[0]->Draw();
    }



    output->Write();
    output->Close();

    input->Close();




}
