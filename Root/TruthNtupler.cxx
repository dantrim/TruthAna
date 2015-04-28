#include "TruthAna/TruthNtupler.h"

using namespace std;
using namespace Susy;


struct SortByPt {
    bool operator()(TLorentzVector* a, TLorentzVector* b) const {
        return a->Pt() > b->Pt();
    }
};

/* --------------------------------------------------- */
// Constructor
/* --------------------------------------------------- */
TruthNtupler::TruthNtupler(TString MCID, TString suffix) : 
    SusyNtTools()
{
    m_filename = MCID+"_"+suffix+".root";
    m_file     = TFile::Open(m_filename,"recreate");
    m_tree     = new TTree("truthNt","truthNt");
    (m_tree)->SetAutoSave(1000000);
    (m_tree)->SetDirectory(m_file);
}
/* --------------------------------------------------- */
// SetBranches
/* --------------------------------------------------- */
void TruthNtupler::BookTree()
{
    // event variables
    m_tree->Branch("runNumber",   &b_runNumber, "runNumber/I");
    m_tree->Branch("eventNumber", &b_eventNumber, "eventNumber/I");
    m_tree->Branch("mcid",        &b_mcid, "mcid/I");
    m_tree->Branch("eventweight", &b_eventweight, "eventweight/D");
    
    // lepton variables
    m_tree->Branch("nlep",      &b_nlep, "nlep/I");
    m_tree->Branch("lept1Pt",   &b_lept1Pt,     "lept1Pt/F");
    m_tree->Branch("lept1Eta",  &b_lept1Eta,    "lept1Eta/F");
    m_tree->Branch("lept1Phi",  &b_lept1Phi,    "lept1Phi/F");
    m_tree->Branch("lept1E",    &b_lept1E,      "lept1E/F");
    m_tree->Branch("lept2Pt",   &b_lept2Pt,     "lept2Pt/F");
    m_tree->Branch("lept2Eta",  &b_lept2Eta,    "lept2Eta/F");
    m_tree->Branch("lept2Phi",  &b_lept2Phi,    "lept2Phi/F");
    m_tree->Branch("lept2E",    &b_lept2E,      "lept2E/F");

    
    // dilepton variables
    m_tree->Branch("isOS",    &b_isOS,      "isOS/O");
    m_tree->Branch("isElEl",  &b_isElEl,    "isElEl/O");
    m_tree->Branch("isMuMu",  &b_isMuMu,    "isMuMu/O");
    m_tree->Branch("isElMu",  &b_isElMu,    "isElMu/O");
    m_tree->Branch("dphi_ll", &b_dphi_ll,   "dphi_ll/F");
    m_tree->Branch("deta_ll", &b_deta_ll,   "deta_ll/F");
    m_tree->Branch("dR_ll",   &b_dR_ll,     "dR_ll/F");
    m_tree->Branch("pTll",    &b_pTll,     "pTll/F");
    m_tree->Branch("dphi_ll", &b_dphi_ll,   "dphi_ll/F");
    m_tree->Branch("mll",     &b_mll,      "mll/F");
    
    // MET vars
    m_tree->Branch("met",     &b_met,       "met/F");
    m_tree->Branch("metphi",  &b_metPhi,   "metphi/F");
    m_tree->Branch("metrel",  &b_metrel,    "metrel/F");
    
    // jet vars
    m_tree->Branch("njets",       &b_njets, "njets/I");
    m_tree->Branch("nCentralLightJets",   &b_nCL20jets, "nCL20jets/I");
    m_tree->Branch("nCentralBJets",    &b_nB20jets, "nB20jets/I");
    m_tree->Branch("nForwardJets",    &b_nF30jets, "nF30jets/I");
    m_tree->Branch("j_isCL20",    b_j_isCL20, "j_isCL20[njets]/O");
    m_tree->Branch("j_isB20",     b_j_isB20,  "j_isB20[njets]/O");
    m_tree->Branch("j_isF30",     b_j_isF30,  "j_isF30[njets]/O");
    m_tree->Branch("j_pt",        b_j_pt, "j_pt[njets]/F");
    m_tree->Branch("j_eta",       b_j_eta, "j_eta[njets]/F");
    m_tree->Branch("j_phi",       b_j_phi, "j_phi[njets]/F");
    m_tree->Branch("j_e",         b_j_e,  "j_e[njets]/F");

    m_tree->Branch("jet1Pt",        &b_jet1Pt,  "jet1Pt/F");
    m_tree->Branch("jet1Eta",       &b_jet1Eta, "jet1Eta/F");
    m_tree->Branch("jet1Phi",       &b_jet1Phi, "jet1Phi/F");
    m_tree->Branch("jet2Pt",        &b_jet2Pt,  "jet2Pt/F");
    m_tree->Branch("jet2Eta",       &b_jet2Eta, "jet2Eta/F");
    m_tree->Branch("jet2Phi",       &b_jet2Phi, "jet2Phi/F");

    // Razor
    m_tree->Branch("dphi_ll_vBetaT",     &b_dphi_ll_vBetaT, "dphi_ll_vBetaT/F");
    m_tree->Branch("mDeltaR",            &b_mDeltaR,        "mDeltaR/F");
    m_tree->Branch("R2",                 &b_R2,             "R2/F");


    clearOutputBranches();
    return;
}
    
/* --------------------------------------------------- */
// Clear branches
/* --------------------------------------------------- */
void TruthNtupler::clearOutputBranches()
{
    // event var
    b_runNumber=-1;
    b_eventNumber=-1;
    b_mcid=-1;
    b_eventweight=1;

    // lepton var
    b_nlep=0;
    b_lept1Pt=-999;
    b_lept1Eta=-999;
    b_lept1Phi=-999;
    b_lept1E=-999;
    b_lept1q=0;
    b_lept1Flav=-999;
    b_lept2Pt=-999;
    b_lept2Eta=-999;
    b_lept2Phi=-999;
    b_lept2E=-999;
    b_lept2q=0;
    b_lept2Flav=-999;
    
    // dilepton var
    b_isOS=false;
    b_isElEl=false;
    b_isMuMu=false;
    b_isElMu=false;
    b_dphi_ll=-999;
    b_deta_ll=-999;
    b_dR_ll=-999;
    b_pTll=-999;
    b_dphi_ll=-999;
    b_mll=-999;

    // met var
    b_met=-999;
    b_metPhi=-999;
    b_metrel=-999;
    
    // jet var
    b_njets=0;
    b_nCL20jets=0;
    b_nB20jets=0;
    b_nF30jets=0;
    for(uint i=0;i<25; i++){
        b_j_isCL20[i]=false;
        b_j_isB20[i]=false;
        b_j_isF30[i]=false;
        b_j_pt[i]=-999;
        b_j_eta[i]=-999;
        b_j_phi[i]=-999;
        b_j_e[i]=-999;
    }
    b_jet1Pt=-999;
    b_jet1Eta=-999;
    b_jet1Phi=-999;
    b_jet2Pt=-999;
    b_jet2Eta=-999;
    b_jet2Phi=-999;

    // razor
    b_dphi_ll_vBetaT=-999;
    b_mDeltaR=-999;
    b_R2=-999;

    return;
}
/* --------------------------------------------------- */
// Save tree
/* --------------------------------------------------- */
void TruthNtupler::SaveTree()
{
    m_file->cd();
    m_tree->SetDirectory(m_file);
    m_tree->Write();
    m_tree->SetDirectory(0);
    delete m_tree;
    m_file->Close();
    delete m_file;
}
void TruthNtupler::WriteTree()
{
    m_tree->Fill();
    clearOutputBranches();
}

void TruthNtupler::setSumW(double sumw)
{
    TH1D *h_sumw = new TH1D("sumOfMCWeights_"+m_mcid, "sumOfMCWeights_"+m_mcid,1,0.,1.);
    h_sumw->Fill(0.,sumw);
    m_file->cd();
    h_sumw->SetDirectory(m_file);
    h_sumw->Write();
    h_sumw->SetDirectory(0);
    delete h_sumw;
    return;
}

/* --------------------------------------------------- */
//  Methods for filling the output ntuple
/* --------------------------------------------------- */

// event vars
void TruthNtupler::fillEventVars(int run, int event, int mcid, double weight)
{
    b_runNumber     = run;
    b_eventNumber   = event;
    b_mcid          = mcid;
    b_eventweight   = weight;
}

// lepton vars
void TruthNtupler::fillLeptonVars(const TruthParticleVector& electrons, const TruthParticleVector& muons,
                            const TruthMet& met)
{
    if( (electrons.size()+muons.size()) < nLepMax) b_nlep=(electrons.size()+muons.size());
    else{
        cerr<<"fillLeptonVars() error: too many input leptons (" << electrons.size()+muons.size() << endl;
        abort();
    }
    TruthParticleVector leptons;
    for(uint ie=0; ie<electrons.size(); ++ie){
        leptons.push_back(electrons.at(ie));
    }
    for(uint im=0; im<muons.size(); ++im){
        leptons.push_back(muons.at(im));
    }
    int nlep = leptons.size();
    int nmu=muons.size();
    int nel=electrons.size();
    if(nlep==2 && nel==2){
        b_isElEl=1;
        b_isMuMu=0;
        b_isElMu=0;
    }
    if(nlep==2 && nmu==2){
        b_isElEl=0;
        b_isMuMu=1;
        b_isElMu=0;
    }
    if(nlep==2 && nmu==1 && nel==1){
        b_isElEl=0;
        b_isMuMu=0;
        b_isElMu=1;
    }
    // sort leptons
    sort(leptons.begin(), leptons.end(),SortByPt());
    
    b_nlep = nlep;
    b_isOS      = (leptons.at(0)->charge * leptons.at(1)->charge)<0 ? true : false;
    b_lept1Pt   = leptons.at(0)->Pt()*1000.;
    b_lept1Eta  = leptons.at(0)->Eta();
    b_lept1Phi  = leptons.at(0)->Phi();
    b_lept1E    = leptons.at(0)->E()*1000.;
    b_lept1q    = leptons.at(0)->charge;
    b_lept1Flav =abs(leptons.at(0)->pdgId)==11 ? 1 : 0;
    b_lept2Pt   = leptons.at(1)->Pt()*1000.;
    b_lept2Eta  = leptons.at(1)->Eta();
    b_lept2Phi  = leptons.at(1)->Phi();
    b_lept2E    = leptons.at(1)->E()*1000.;
    b_lept2q    = leptons.at(1)->charge;
    b_lept2Flav =abs(leptons.at(1)->pdgId)==11 ? 1 : 0;

    TLorentzVector ll_tlv = (*leptons.at(0) + *leptons.at(1));
    b_mll = ll_tlv.M()*1000.;
    b_pTll = ll_tlv.Pt()*1000.;
    b_deta_ll = abs(leptons.at(0)->Eta()-leptons.at(1)->Eta());
    b_dphi_ll = abs(leptons.at(0)->Phi()-leptons.at(1)->Phi());
    
    

 //   TLorentzVector ll_tlv;
 //   int qqType=1;
 //   for(uint ilep=0; ilep<leptons.size();ilep++){
 //       const Susy::TruthParticle* lep = leptons.at(ilep);
 //       if(ilep<2){
 //           ll_tlv += (*lep);
 //           qqType *= lep->charge;
 //       }
 //   }
 //   if(leptons.size()>1){
 //       b_isOS = (qqType<0) ? true : false;
 //       b_dphi_ll = fabs(leptons.at(0)->DeltaPhi(*leptons.at(1)));
 //       b_deta_ll = leptons.at(0)->Eta()-leptons.at(1)->Eta();
 //       b_dR_ll = leptons.at(0)->DeltaR(*leptons.at(1));
 //       b_pT_ll = ll_tlv.Pt();
 //       b_m_ll = ll_tlv.M();
 //   }
}
void TruthNtupler::fillMetVars(const TruthMet& met)
{
    b_met = met.lv().Pt()*1000.;
    b_metPhi = met.lv().Phi();
}        
void TruthNtupler::fillJetVars(const TruthJetVector& signalJets, const TruthJetVector& centralJets,
                    const TruthJetVector& forwardJets, const TruthJetVector& centralBJets,
                    const std::vector<int> clIdx, const std::vector<int> b20Idx, const std::vector<int> fIdx)
{
    b_njets = signalJets.size();
    b_nCL20jets = centralJets.size();
    b_nB20jets = centralBJets.size();
    b_nF30jets = forwardJets.size();

    for(uint ij=0;ij<signalJets.size();++ij){
        if(ij>25) continue;
        b_j_pt[ij] = signalJets.at(ij)->Pt()*1000.;
        b_j_eta[ij] = signalJets.at(ij)->Eta();
        b_j_phi[ij] = signalJets.at(ij)->Phi();

        bool isCentral=false;
        bool isBjet=false;
        bool isForward=false;
        for(uint ic=0; ic<clIdx.size();ic++){
            if(ij==clIdx[ic]) { b_j_isCL20[ij] = true; isCentral = true; }
        }
        for(uint ifj=0; ifj<fIdx.size();ifj++){
            if(ij==fIdx[ifj]) { b_j_isF30[ij] = true; isForward = true; }
        }
        for(uint ibj=0; ibj<b20Idx.size();ibj++){
            if(ij==b20Idx[ibj]) { b_j_isB20[ij] = true; isBjet = true; }
        }
        if(isCentral && isBjet && isForward) {
            cout << "fillJetVars error: Jet is central && bjet && forward ! " << endl;
        }
    }

    if(centralJets.size()>0){
        b_jet1Pt  = centralJets.at(0)->Pt()*1000.;
        b_jet1Eta = centralJets.at(0)->Eta();
        b_jet1Phi = centralJets.at(0)->Phi();
    }
    if(centralJets.size()>1){
        b_jet2Pt  = centralJets.at(1)->Pt()*1000.;
        b_jet2Eta = centralJets.at(1)->Eta();
        b_jet2Phi = centralJets.at(1)->Phi(); 
    }
}
void TruthNtupler::fillRazorVars(double MDR, double DPB, double R2)
{
    b_dphi_ll_vBetaT = DPB;
    b_mDeltaR = MDR*1000.;
    b_R2 = R2;
}
