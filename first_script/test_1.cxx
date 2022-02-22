#include <iostream>

#include "useful_functions.h"

#define PI 3.14159265

using namespace std;

const int number_of_daughters = 2;
const int D0_PDG = 421;

const TString script_name = "test_1";

void test_1(){

    vector<double> v1 = {1., 0, 0};
    vector<double> v2 = {-1, 0, 0.0};

    auto angle = GetAngleBetweenVectors(v1, v2);

    cout << "Angle: " << angle <<  " rad (" << angle * 180 / PI << " grad)" << endl;


      // check the field
    KFParticle::SetField(0.);
    {
        float xyz[3] = {0,0,0};
        float B[3];
        KFParticle p;
        p.GetFieldValue(xyz,B);
        std::cout<<"Field is set to " <<B[0]<<" "<<B[1]<<" "<<B[2]<<std::endl;
    }
    
    //TestFilterBit(16) was used (track was used to fit the primary vertex)
    double cov[21]; //Kovarianz Matrix
    cov[ 0 ] = 0.004;
    cov[ 1 ] = 5.e-5;   cov[ 2 ] = 0.003;
    cov[ 3 ] = 1.7e-6; cov[ 4 ] = -5.5e-6; cov[ 5 ] = 0.014;
    cov[ 6 ] = -3.e-4;  cov[ 7 ] = 4.5e-5;  cov[ 8 ] = -5e-6;  cov[ 9 ] = 0.05;
    cov[10] = -1.e-5;  cov[11] = -2.e-4;   cov[12] = -3.e-6; cov[13] = 0.0013; cov[14] = 0.04;
    cov[15] = -3.e-5;  cov[16] = 5.e-5;    cov[17] = -6.e-4; cov[18] = 0.004;   cov[19] = 0.01; cov[20] = 0.1;

    double errors[6]; //x, y, z in mm und p_X, p_y, p_z in GeV
    errors[0] = sqrt(cov[0]); errors[1] = sqrt(cov[2]); errors[2] = sqrt(cov[5]);
    errors[3] = sqrt(cov[9]); errors[4] = sqrt(cov[14]); errors[5] = sqrt(cov[20]);


    //ROOT file reading...
    TChain chain("Events");
    TString logic_folder = "D0_decays";
    //TString logic_folder = "D0_decays_randomized";
    //TString path_to_file = "/Users/vlad/HEP/GEANT4/" + logic_folder + "/build/D0_decay_output.root";
     TString path_to_file = "D0_decay_output.root";
    cout << path_to_file << endl;
    chain.Add(path_to_file);

    TString output_name = "output_" + logic_folder + "_" + script_name + ".root";
    TFile output_file(output_name,"RECREATE");

    TH1F *x_res = new TH1F("x_res", "X residual", 1000, -10, 10);   // x, y, z , mm
    TH1F *x_pull = new TH1F("x_pull", "X pull", 1000, -10, 10);
    TH1F *y_res = new TH1F("y_res", "Y residual", 1000, -10, 10);
    TH1F *y_pull = new TH1F("y_pull", "Y pull", 1000, -10, 10);
    TH1F *z_res = new TH1F("z_res", "Z residual", 1000, -50, 50);
    TH1F *z_pull = new TH1F("z_pull", "Z pull", 1000, -10, 10);
    
   // TH1F *length_res = new TH1F("length_res", "Track length resolution", 1000, -1, 1);
    //TH1F *length_pull = new TH1F("length_pull", "Track length pull", 1000, -0.001, 0.001);
    TH1F *delta_phiHist = new TH1F( "angle", "Between daughter angle", 1000, 0, TMath::TwoPi()  );
    TH2F *x_res_angle = new TH2F( "x_res_angle", "X Res and angle", 1000, -10, 10, 1000, 0, TMath::TwoPi()  );

    TH1F *Px_res = new TH1F("Px_res", "PX residual", 1000, -10, 10);  // Px, Py, Pz, GeV
    TH1F *Px_pull = new TH1F("Px_pull", "PX pull", 1000, -10, 10);
    TH1F *Py_res = new TH1F("Py_res", "PY residual", 1000, -10, 10);
    TH1F *Py_pull = new TH1F("Py_pull", "PY pull", 1000, -10, 10);
    TH1F *Pz_res = new TH1F("Pz_res", "PZ residual", 1000, -10, 10);
    TH1F *Pz_pull = new TH1F("Pz_pull", "PZ pull", 1000, -10, 10);

    TH1F *E_res = new TH1F("e_res", "E residual", 1000, -10, 10);  // E
    TH1F *E_pull = new TH1F("e_pull", "E pull", 1000, -10, 10);

    TH1F *Chi2 = new TH1F("chi2/ndf", "Chi2/NDF", 1000, -8, 8);
    
    TH1F *VertDevHist = new TH1F("VertDevHist", "Distance between vertex and daughter tracks coord.", 1000, -100, 20000);

    const int NMaxTrack = 100;

    Int_t nTracks;
    Int_t trackID[NMaxTrack];
    Int_t parentID[NMaxTrack];
    Int_t pdg[NMaxTrack];
    Int_t charge[NMaxTrack];
    Double_t trackX[NMaxTrack];
    Double_t trackY[NMaxTrack];
    Double_t trackZ[NMaxTrack];
    Double_t vertexX[NMaxTrack];
    Double_t vertexY[NMaxTrack];
    Double_t vertexZ[NMaxTrack];
    Double_t trackLength[NMaxTrack];
    Double_t initialPX[NMaxTrack];
    Double_t initialPY[NMaxTrack];
    Double_t initialPZ[NMaxTrack];
    Double_t finalPX[NMaxTrack];
    Double_t finalPY[NMaxTrack];
    Double_t finalPZ[NMaxTrack];
    Double_t mass[NMaxTrack];

    chain.SetBranchAddress( "nTracks", &nTracks);
    chain.SetBranchAddress( "trackID", trackID);
    chain.SetBranchAddress( "parentID", parentID);
    chain.SetBranchAddress( "pdg", pdg);
    chain.SetBranchAddress( "charge", charge);
    chain.SetBranchAddress( "trackX", trackX);
    chain.SetBranchAddress( "trackY", trackY);
    chain.SetBranchAddress( "trackZ", trackZ);
    chain.SetBranchAddress( "vertexX", vertexX);
    chain.SetBranchAddress( "vertexY", vertexY);
    chain.SetBranchAddress( "vertexZ", vertexZ);
    chain.SetBranchAddress( "trackLength", trackLength);
    chain.SetBranchAddress( "initialPX", initialPX);
    chain.SetBranchAddress( "initialPY", initialPY);
    chain.SetBranchAddress( "initialPZ", initialPZ);
    chain.SetBranchAddress( "finalPX", finalPX);
    chain.SetBranchAddress( "finalPY", finalPY);
    chain.SetBranchAddress( "finalPZ", finalPZ);
    chain.SetBranchAddress( "mass", mass);

    int nEvents = chain.GetEntries();

    size_t daughterCounter = 0;
    size_t motherCounter = 0;



    for (int iEvent = 0; iEvent < nEvents; iEvent++){
        chain.GetEntry(iEvent);

        if(iEvent % 1 == 0) {
            cout << "processing " << iEvent << " event\r";
            cout << flush;
        }
        
        if (nTracks < number_of_daughters+1)
            continue;
                
        motherCounter = 0;
        daughterCounter = 0;

        KFParticle motherParticle;
        vector<KFParticle> daughterParticles(number_of_daughters);
        vector<vector<double>> daughterVectors(number_of_daughters);

        for (int iTrack = 0; iTrack < nTracks; iTrack++){
            if (daughterCounter >= number_of_daughters && motherCounter > 0)
                break;
            if (parentID[iTrack]==0){
                motherCounter++;
                if (pdg[iTrack] != D0_PDG){
                    cout << "1st particle is not D0... error..." << endl;
                    break;
                } else {
                    double motherParam[6] = {trackX[iTrack], trackY[iTrack], trackZ[iTrack], finalPX[iTrack]/1000, finalPY[iTrack]/1000, finalPZ[iTrack]/1000};
                    motherParticle.Create(motherParam, cov, charge[iTrack], mass[iTrack]/1000);
                }
            }
            if (parentID[iTrack] == 1){
                daughterVectors[daughterCounter] = {initialPX[iTrack]/1000,
                                                                            initialPY[iTrack]/1000,
                                                                            initialPZ[iTrack]/1000};
                double daughterParam[6] = {vertexX[iTrack], vertexY[iTrack], vertexZ[iTrack], initialPX[iTrack]/1000, initialPY[iTrack]/1000, initialPZ[iTrack]/1000};
                for( int i=0; i<6; i++){
                    daughterParam[i] = daughterParam[i] + gRandom->Gaus(0,errors[i]);
                }
                daughterParticles[daughterCounter].Create(daughterParam, cov, charge[iTrack], mass[iTrack]/1000);
                
                double devDist = sqrt (pow(vertexX[iTrack] - daughterParam[0], 2) +
                                                     pow(vertexY[iTrack] - daughterParam[1], 2) +
                                                     pow(vertexZ[iTrack] - daughterParam[2], 2) );

                VertDevHist -> Fill(devDist);

                daughterCounter++;
            }
        }
        //processing...
        
        KFParticle motherRecovered( daughterParticles[0], daughterParticles[1]);

        x_res->Fill(motherRecovered.X() - motherParticle.X());
        x_pull->Fill( (motherRecovered.X() - motherParticle.X())/motherRecovered.GetErrX() );
        y_res->Fill(motherRecovered.Y() - motherParticle.Y());
        y_pull->Fill( (motherRecovered.Y() - motherParticle.Y())/motherRecovered.GetErrY() );
        z_res->Fill(motherRecovered.Z() - motherParticle.Z());
        z_pull->Fill( (motherRecovered.Z() - motherParticle.Z())/motherRecovered.GetErrZ() );
        
      /*  length_res->Fill(motherRecovered.GetDecayLength() - trackLength[0]);
        cout << motherRecovered.GetDecayLength() << endl;
        cout << trackLength[0] << endl;
        cout << endl;
        length_pull->Fill( (motherRecovered.GetDecayLength() - trackLength[0] )/motherRecovered.GetErrDecayLength() );*/

        //vector<double> motherOrigCoord = {motherParticle.X(),motherParticle.Y(), motherParticle.Z()};
        //vector<double> daughterOrigCoord_1 = {daughterParticles[0].X(), daughterParticles[0].Y(), daughterParticles[0].Z()};
       // vector<double> daughterOrigCoord_2 = {daughterParticles[1].X(), daughterParticles[1].Y(), daughterParticles[1].Z()};
        //double delta_phi = GetAngleBetweenVectors(daughterOrigCoord_1,daughterOrigCoord_2);
        double delta_phi = GetAngleBetweenVectors(daughterVectors[0],daughterVectors[1]);
        delta_phiHist -> Fill(delta_phi);
        x_res_angle->Fill( (motherRecovered.X() - motherParticle.X()), delta_phi) ; 

        Px_res->Fill(motherRecovered.Px() - motherParticle.Px());
        Px_pull->Fill( (motherRecovered.Px() - motherParticle.Px())/motherRecovered.GetErrPx() );
        Py_res->Fill(motherRecovered.Py() - motherParticle.Py());
        Py_pull->Fill( (motherRecovered.Py() - motherParticle.Py())/motherRecovered.GetErrPy() );
        Pz_res->Fill(motherRecovered.Pz() - motherParticle.Pz());
        Pz_pull->Fill( (motherRecovered.Pz() - motherParticle.Pz())/motherRecovered.GetErrPz() );

        E_res->Fill(motherRecovered.E() - motherParticle.E());
        E_pull->Fill( (motherRecovered.E() - motherParticle.E())/motherRecovered.GetErrE() );

        Chi2->Fill( motherRecovered.GetChi2() / motherRecovered.GetNDF() );

    }
    cout << endl;


    output_file.cd();
    x_res->Write();
    x_pull->Write();
    y_res->Write();
    y_pull->Write();
    z_res->Write();
    z_pull->Write();
    
    //length_res->Write();
   // length_pull->Write();
    delta_phiHist->Write();
    x_res_angle->Write();

    Px_res->Write();
    Px_pull->Write();
    Py_res->Write();
    Py_pull->Write();
    Pz_res->Write();
    Pz_pull->Write();

    E_res->Write();
    E_pull->Write();

    Chi2->Write();
    
    VertDevHist -> Write();

    output_file.Close();


    cout << "OK " << endl;
    gROOT->ProcessLine("new TBrowser");

}




//############################################

