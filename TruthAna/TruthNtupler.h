#ifndef TruthAna_TruthNtupler_h
#define TruthAna_TruthNtupler_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include <iostream>
#include "TLorentzVector.h"


// SusyNtuple
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNt.h"
#include "SusyNtuple/SusyNtObject.h"
#include "SusyNtuple/SusyNtTools.h"

namespace Susy {
class TruthNtupler : public SusyNtTools
{
    public:
        TruthNtupler(TString MCID, TString suffix); 
        ~TruthNtupler(){};

        void BookTree();        ///< link objects & vars to the branches
        void setSumW(double sumw);
        void SaveTree();        ///< save ntuple to output file
        void WriteTree();       ///< write output tree and close output file

        // group fillings of branches
        void fillEventVars(int run, int event, int mcid,
                    double w);
        void fillLeptonVars(const TruthParticleVector& electrons, const TruthParticleVector& muons,
                            const TruthMet& met);
        void fillMetVars(const TruthMet& met);
        void fillJetVars(const TruthJetVector& sigJets, const TruthJetVector& centralJets, const TruthJetVector& forwardJets,
                        const TruthJetVector& bJets, const std::vector<int> clIdx, const std::vector<int> bIdx, const std::vector<int> fIdx);
        void fillRazorVars(double MDR, double DPB, double R2);



    ///////////////////////////////////////
    // BRANCHES
    ///////////////////////////////////////
        
        // event variables
        int b_runNumber;
        int b_eventNumber;
        int b_mcid;
        double b_eventweight;
        
        // lepton variables
        static const unsigned int nLepMax=4;
        int b_nlep;
        float b_lept1Pt;
        float b_lept1Eta;
        float b_lept1Phi;
        float b_lept1E;
        float b_lept1q;
        int b_lept1Flav;

        float b_lept2Pt;
        float b_lept2Eta;
        float b_lept2Phi;
        float b_lept2E;
        float b_lept2q;
        int b_lept2Flav;

        // dilepton variables
        bool b_isOS;
        bool b_isElEl;
        bool b_isMuMu;
        bool b_isElMu;
        float b_dphi_ll;
        float b_deta_ll;
        float b_dR_ll;
        float b_pTll;
        float b_mll;
        
        // MET vars
        float b_met;
        float b_metPhi;
        float b_metrel;

        // jet variables
        int b_njets;
        int b_nCL20jets;
        int b_nB20jets;
        int b_nF30jets;
        bool b_j_isCL20[25];
        bool b_j_isB20[25];
        bool b_j_isF30[25];
        float b_j_pt[25];
        float b_j_eta[25];
        float b_j_phi[25];
        float b_j_e[25];

        // these jet variables refer to central-light jets!
        float b_jet1Pt;
        float b_jet1Eta;
        float b_jet1Phi;
        float b_jet2Pt;
        float b_jet2Eta;
        float b_jet2Phi;

        // Razor ana
        float   b_dphi_ll_vBetaT;
        float b_mDeltaR;
        float b_R2;


        ClassDef(TruthNtupler,1);
        


    private:
        bool m_initialize;
        TString m_filename;
        TFile*  m_file;
        TTree*  m_tree;
        void clearOutputBranches();
        TString m_mcid;
}; // class TruthNtupler

}; // namespace Susy


#endif
