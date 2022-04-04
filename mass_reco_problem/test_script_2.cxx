#include <iostream>
#include <string>
#include <sstream>

//#include "useful_functions.h"

#include "CovMatDatabase/CovMatDatabase.h"
#include "CovMatDatabase/CovMatDatabase.cc"

#define PI 3.14159265

const int number_of_daughters = 2;
const int D0_PDG = 421;
const int KMINUS_PDG = -321;
const int PIPLUS_PDG = 211;

const float ZERO_COVMAT[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
const float ALMOST_ZERO_COVMAT[21] = {0,0,0,0,0,0,0,0,0,0.,0,0,0,0.,0.,0,0,0,0.,0.,0.1};
const float TEST_COVMAT[21] = {22745.8, -2925.62, 25838.8, -197.81, 1554.37, 101387, -0.483599, -0.102425, -0.0317067, 3.73148e-05, 0.208339, -0.24716, -0.0746562, 1.61557e-06, 3.11066e-05, 0.061912, 0.00167607, -0.696451, -6.09529e-07, -3.39482e-06, 2.28536e-05};

// structure to maintain MC particle quantities
// implemented to laconize the code
struct MCParticle{
    float X;
    float Y;
    float Z;
    float Px;
    float Py;
    float Pz;
    float M;
    float Charge;
    float PDG;
    float NDF;
    float Chi2;

    float Pt(){return Px*Px+Py*Py;}
    float P(){return Px*Px+Py*Py+Pz*Pz;}
    float E(){return Px*Px+Py*Py+Pz*Pz+M*M;}
    float Eta(){return 0.5*log((P()+Pz)/(P()-Pz));}
    float Phi(){return atan2(Py, Pz);}

};

// 2 examples of functions which will extract covariance matrix correspondong to track Pt
// (defined below the main script)
// this functions are written in accordance to LHC17m_278915_database.tsv database
// in case you use a different database, you have to rewrite this functions
// or better write your own function(s)
vector<float> GetCovarianceMatrixByPt1(CovMatDatabase database, double track_pt);
vector<float> GetCovarianceMatrixByPt2(CovMatDatabase database, double track_pt);

// ... ADD DESCRIPTION !!!
template<typename T>
TGraphErrors *MakeMeanGraphFromTH2(T fHisto2D);
template<typename T>
TGraphErrors *MakeSigmaGraphFromTH2(T fHisto2D);

const TString script_name = "test_script_2";

void test_script_2(){
    //////////////////////////////////////////////////////
    //First of all, lets read covariance matrix database
    // and check all the matrix that it is really
    // covariance matrixes
    //////////////////////////////////////////////////////
    CovMatDatabase covmat_data("CovMatDatabase/LHC17m_278915_database.tsv");
    //covmat_data.PrintListOfMatNames();
    bool is_database_ok = covmat_data.CheckDatabase();
    if (!is_database_ok){
        cerr << "\nBad database. Exit.\n";
        return;
    }

    cout << "\n";
    //////////////////////////////////////////////////////
    // Lets setup magnetic field, used in MC simulations.
    // According to ALICE experiment, homogenious mag. field
    // of 0.5 Tl was used.
    // You have to setup your own parameters
    ///////////////////////////////////////////////////////
    #define HomogeneousField

    #ifdef HomogeneousField //ALICE
    std::cout << "HomogeneousField option is set" << std::endl;
    #endif
    #ifdef NonhomogeneousField //cbm
    std::cout << "NonhomogeneousField option is set" << std::endl;
    #endif

    // check the field
    KFParticle::SetField(0.5);
    {
        float xyz[3] = {0,0,0};
        float B[3];
        KFParticle p;
        p.GetFieldValue(xyz,B);
        std::cout<<"Field is set to " <<B[0]<<" "<<B[1]<<" "<<B[2]<<std::endl;
    }

    ///////////////////////////////////////////////////////
    // Lets read data from the input file and create output file
    //ROOT file reading...
    ///////////////////////////////////////////////////////
    TChain chain("Events");
    TString logic_folder = "D0_decays";
    TString path_to_file = "./data/D0_decay_output.root";
    cout << path_to_file << endl;
    chain.Add(path_to_file);
    // output creating
    TString output_name = "output_" + logic_folder + "_" + script_name + ".root";
    TFile output_file(output_name,"RECREATE");



    ///////////////////////////////////////////////////////////
    // Next step is to create a directories structure
    // and make a histograms corresponding to this directories
    //////////////////////////////////////////////////////////

    TDirectory* QAHistos = output_file.mkdir("QAHistos");
        TDirectory* GeantMC = QAHistos->mkdir("GeantMC");
            TDirectory* MotherMC_D0 = GeantMC->mkdir("Mother_D0");
                TH1F *fHistXMotherMC_D0 = new TH1F("fHistXMotherMC_D0", "X MC", 1000, -1e5, 1e5);
                TH1F *fHistYMotherMC_D0 = new TH1F("fHistYMotherMC_D0", "Y MC",1000, -1e5, 1e5);
                TH1F *fHistZMotherMC_D0 = new TH1F("fHistZMotherMC_D0", "Z MC",1000, -3e5, 3e5);
                TH1F *fHistPxMotherMC_D0 = new TH1F("fHistPxMotherMC_D0", "Px MC", 1000, -10, 10);
                TH1F *fHistPyMotherMC_D0 = new TH1F("fHistPyMotherMC_D0", "Py MC",1000, -10, 10);
                TH1F *fHistPzMotherMC_D0 = new TH1F("fHistPzMotherMC_D0", "Pz MC",1000, -10, 10);
                TH1I *fHistChargeMotherMC_D0 = new TH1I("fHistChargeMotherMC_D0", "Charge MC",50, -25, 25);
                TH1I *fHistPDGMotherMC_D0 = new TH1I("fHistPDGMotherMC_D0", "PDG MC",2000, -1000, 1000);
                TH1F *fHistPtMotherMC_D0 = new TH1F("fHistPtMotherMC_D0", "Pt MC",1000,0.,10.);
                TH1F *fHistPMotherMC_D0 = new TH1F("fHistPMotherMC_D0", "P MC",1000,0.,30.);
                TH1F *fHistOneOverPtMotherMC_D0 = new TH1F("fHistOneOverPtMotherMC_D0", "1/Pt MC",1000,0.,20.);
                TH1F *fHistEMotherMC_D0 = new TH1F("fHistEMotherMC_D0", "E MC",1000,0.,30.);
                TH1F *fHistMassMotherMC_D0 = new TH1F("fHistMassMotherMC_D0", "Mass MC",1000,0.,10.);
                TH1F *fHistChi2MotherMC_D0 = new TH1F("fHistChi2MotherMC_D0", "Chi2 MC",1000, -10, 10);

    TH1F *fHistMassDiffZeroCovMat = new TH1F("fHistMassDiffZeroCovMat", "MC - KF",100, -10, 10);
    TH1F *fHistMassDiffNonZeroCovMat = new TH1F("fHistMassDiffNonZeroCovMat", "MC - KF",100, -10, 10);



    //////////////////////////////////////////////////////////////////////
    //  Now we are reading event and track info (tree) from input file
    //////////////////////////////////////////////////////////////////////
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


    //////////////////////////////////////////////////////////
    //      Start of the analysis code
    //////////////////////////////////////////////////////////

    int nEvents = chain.GetEntries(); //read number of MC events to analyse


    ////////////////////////////////////////////////////////
    //          START OF THE EVENT LOOP
    ////////////////////////////////////////////////////////
    for (int iEvent = 0; iEvent < nEvents; iEvent++){
        chain.GetEntry(iEvent);

        // print a processing progress
        if(iEvent % 10000 == 0) {
            //cout << "processing " << iEvent << " event\r";
            //cout << flush;
        }

        // create MC particles structure
        MCParticle motherMC_D0;
        //MCParticle daughterMC_Kminus;
        //MCParticle daughterMC_PIplus;

        // create KFParticle structure for particles
        KFParticle motherKF_D0_ZEROCOVMAT;
        KFParticle motherKF_D0_NONZEROCOVMAT;

        // loop over the tracks
        for (int iTrack = 0; iTrack < nTracks; iTrack++){
            if (parentID[iTrack]!=0) continue;

            motherMC_D0.X = (float)trackX[iTrack];
            motherMC_D0.Y = (float)trackY[iTrack];
            motherMC_D0.Z = (float)trackZ[iTrack];
            motherMC_D0.Px = (float)finalPX[iTrack];
            motherMC_D0.Py = (float)finalPY[iTrack];
            motherMC_D0.Pz = (float)finalPZ[iTrack];
            motherMC_D0.M = (float)mass[iTrack];
            motherMC_D0.Charge = (float)charge[iTrack];
            motherMC_D0.PDG = (float)pdg[iTrack];
            motherMC_D0.NDF = 6;
            motherMC_D0.Chi2  = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/6.));

            float motherParam_D0[6] = {motherMC_D0.X,motherMC_D0.Y,motherMC_D0.Z,motherMC_D0.Px,motherMC_D0.Py,motherMC_D0.Pz};
            motherKF_D0_ZEROCOVMAT.Create(motherParam_D0, ZERO_COVMAT, motherMC_D0.Charge, motherMC_D0.Chi2, motherMC_D0.NDF, motherMC_D0.M);
            motherKF_D0_NONZEROCOVMAT.Create(motherParam_D0, ALMOST_ZERO_COVMAT, motherMC_D0.Charge, motherMC_D0.Chi2, motherMC_D0.NDF, motherMC_D0.M);

            cout << "Mass Diff nonzero cov mat: " << motherMC_D0.M - motherKF_D0_NONZEROCOVMAT.GetMass() << "\n";
            fHistMassDiffNonZeroCovMat -> Fill( motherMC_D0.M - motherKF_D0_NONZEROCOVMAT.GetMass());
            fHistMassDiffZeroCovMat -> Fill( motherMC_D0.M - motherKF_D0_ZEROCOVMAT.GetMass());

        }
        // end of track loop

        /////////////////////////////////////////
        //          Fill MC QA histos
        /////////////////////////////////////////

        fHistXMotherMC_D0 -> Fill(motherMC_D0.X);
        fHistYMotherMC_D0 -> Fill(motherMC_D0.Y);
        fHistZMotherMC_D0 -> Fill(motherMC_D0.Z);
        fHistPxMotherMC_D0 -> Fill(motherMC_D0.Px);
        fHistPyMotherMC_D0 -> Fill(motherMC_D0.Py);
        fHistPzMotherMC_D0 -> Fill(motherMC_D0.Pz);
        fHistChargeMotherMC_D0 -> Fill(motherMC_D0.Charge);
        fHistPDGMotherMC_D0 -> Fill(motherMC_D0.PDG);
        fHistPtMotherMC_D0 -> Fill(motherMC_D0.Pt());
        fHistPMotherMC_D0 -> Fill(motherMC_D0.P());
        fHistOneOverPtMotherMC_D0 -> Fill(1./motherMC_D0.Pt());
        fHistEMotherMC_D0 -> Fill(motherMC_D0.E());
        fHistMassMotherMC_D0 -> Fill(motherMC_D0.M);
        fHistChi2MotherMC_D0 -> Fill(motherMC_D0.Chi2);

        ///////////////////////////////////////////
        //      Create KFParticles...
        ///////////////////////////////////////////

        /*float motherParam_D0[6] = {motherMC_D0.X,motherMC_D0.Y,motherMC_D0.Z, motherMC_D0.Px,motherMC_D0.Py,motherMC_D0.Pz};
        vector<float> covVec_D0 = GetCovarianceMatrixByPt2(covmat_data, motherMC_D0.Pt());
        float *cov_D0 = &covVec_D0[0];
        motherKF_D0.Create(motherParam_D0, ALMOST_ZERO_COVMAT, motherMC_D0.Charge, motherMC_D0.Chi2, motherMC_D0.NDF, motherMC_D0.M);

        cout << motherMC_D0.M - motherKF_D0.GetMass() << "\n";

        float daughterParam_Kminus[6] = {daughterMC_Kminus.X,daughterMC_Kminus.Y,daughterMC_Kminus.Z, daughterMC_Kminus.Px,daughterMC_Kminus.Py,daughterMC_Kminus.Pz};
        vector<float> covVec_Kminus = GetCovarianceMatrixByPt2(covmat_data, daughterMC_Kminus.Pt());
        float *cov_Kminus = &covVec_Kminus[0];
        float errors_Kminus[6] = {sqrt(cov_Kminus[0]), sqrt(cov_Kminus[2]), sqrt(cov_Kminus[5]), sqrt(cov_Kminus[9]), sqrt(cov_Kminus[14]), sqrt(cov_Kminus[20]) };
        for( int i=0; i<6; i++)
            daughterParam_Kminus[i] += (float)gRandom->Gaus(0,errors_Kminus[i]);
        daughterKF_Kminus.Create(daughterParam_Kminus, cov_Kminus, daughterMC_Kminus.Charge, daughterMC_Kminus.Chi2, daughterMC_Kminus.NDF, daughterMC_Kminus.M);

        float daughterParam_PIplus[6] = {daughterMC_PIplus.X,daughterMC_PIplus.Y,daughterMC_PIplus.Z, daughterMC_PIplus.Px,daughterMC_PIplus.Py,daughterMC_PIplus.Pz};
        vector<float> covVec_PIplus = GetCovarianceMatrixByPt2(covmat_data, daughterMC_PIplus.Pt());
        float *cov_PIplus = &covVec_PIplus[0];
        float errors_PIplus[6] = {sqrt(cov_PIplus[0]), sqrt(cov_PIplus[2]), sqrt(cov_PIplus[5]), sqrt(cov_PIplus[9]), sqrt(cov_PIplus[14]), sqrt(cov_PIplus[20]) };
        for( int i=0; i<6; i++)
            daughterParam_PIplus[i] += (float)gRandom->Gaus(0,errors_PIplus[i]);
        daughterKF_PIplus.Create(daughterParam_PIplus, cov_PIplus, daughterMC_PIplus.Charge, daughterMC_PIplus.Chi2, daughterMC_PIplus.NDF, daughterMC_PIplus.M);*/

    }
    //// END OF EVENT LOOP
    cout << endl;


    //////////////////////////////////////////////////////////
    //     Write the histograms into the output file
    //////////////////////////////////////////////////////////
    fHistMassDiffZeroCovMat -> Write();
    fHistMassDiffNonZeroCovMat -> Write();

    MotherMC_D0 -> cd();
        fHistXMotherMC_D0 -> Write();
        fHistYMotherMC_D0 -> Write();
        fHistZMotherMC_D0 -> Write();
        fHistPxMotherMC_D0 -> Write();
        fHistPyMotherMC_D0 -> Write();
        fHistPzMotherMC_D0 -> Write();
        fHistChargeMotherMC_D0 -> Write();
        fHistPDGMotherMC_D0 -> Write();
        fHistPtMotherMC_D0 -> Write();
        fHistPMotherMC_D0 -> Write();
        fHistOneOverPtMotherMC_D0 -> Write();
        fHistEMotherMC_D0 -> Write();
        fHistMassMotherMC_D0 -> Write();
        fHistChi2MotherMC_D0 -> Write();

    output_file.Close();

    //cout << "OK " << endl;

    // open TBrowser
    //gROOT->ProcessLine("new TBrowser");
    gROOT->ProcessLine(".q");

}


//#######################################################
//############# FUNCTIONS DETERMINATION #################
//#######################################################

// instance 1 - parsing of matrixes names
vector<float> GetCovarianceMatrixByPt1(CovMatDatabase database, double track_pt){
    //vector<float> covmat(21,0.);
    vector<string> matnames = database.GetListOfMatNames();
    for (auto& name:matnames){
        name.erase(0,3);
        int size=name.size();
        for(int i=0;i<size;i++){
            if( name[i]=='_')
                name[i]=' ';
        }
    }

    size_t num_of_matrix = 0;
    double pt_min = 0., pt_max = 1000.;
    for (int i = 0; i < matnames.size(); i++){
        string s = matnames.at(i);
        istringstream ss(s);
        ss >> pt_min >> pt_max;
        if (track_pt<pt_min || track_pt>pt_max) continue;
        num_of_matrix = i;
        break;
    }

    string res_matname = database.GetListOfMatNames().at(num_of_matrix);
    return database.GetCovMat(res_matname);
}

// instance 2 - if..else.. conditions loop over pt-intervals
vector<float> GetCovarianceMatrixByPt2(CovMatDatabase database, double track_pt){
    vector<float> zero_covmat(21,0.);
    vector<float> res(21);
    const vector<float> borders = {0., 0.3, 0.5, 1., 2., 5., 10., 15., 20., 1000.};

    if (track_pt<borders.at(0) || track_pt>borders.at(borders.size()-1))
        return zero_covmat;

    for (size_t i = 1; i < borders.size(); i++){
        if (track_pt>=borders.at(i-1) && track_pt<=borders.at(i)){
            string res_matname = database.GetListOfMatNames().at(i-1);
            res = database.GetCovMat(res_matname);
        } else continue;
    }

    return res;

}

template<typename T>
TGraphErrors *MakeMeanGraphFromTH2(T fHisto2D){
    int num_of_bins = fHisto2D->GetXaxis()->GetLast();
    float x[num_of_bins], y[num_of_bins], xerr[num_of_bins], yerr[num_of_bins];
    for (int ibin = 1; ibin <= num_of_bins; ibin++){
        x[ibin-1] = fHisto2D->GetXaxis()->GetBinCenter(ibin);
        xerr[ibin-1] = fHisto2D->GetXaxis()->GetBinWidth(ibin)/2;
        auto fProjHisto2D = fHisto2D->ProjectionY("", ibin, ibin) ;
        y[ibin-1] = fProjHisto2D->GetMean();
        yerr[ibin-1] = fProjHisto2D->GetMeanError();
    }
    TGraphErrors *fGraphOut = new TGraphErrors(num_of_bins,x,y,xerr,yerr);
    return fGraphOut;
}

template<typename T>
TGraphErrors *MakeSigmaGraphFromTH2(T fHisto2D){
    int num_of_bins = fHisto2D->GetXaxis()->GetLast();
    float x[num_of_bins], y[num_of_bins], xerr[num_of_bins], yerr[num_of_bins];
    for (int ibin = 1; ibin <= num_of_bins; ibin++){
        x[ibin-1] = fHisto2D->GetXaxis()->GetBinCenter(ibin);
        xerr[ibin-1] = fHisto2D->GetXaxis()->GetBinWidth(ibin)/2;
        auto fProjHisto2D = fHisto2D->ProjectionY("", ibin, ibin) ;
        y[ibin-1] = fProjHisto2D->GetStdDev();
        yerr[ibin-1] = fProjHisto2D->GetStdDevError();
    }
    TGraphErrors *fGraphOut = new TGraphErrors(num_of_bins,x,y,xerr,yerr);
    return fGraphOut;
}
