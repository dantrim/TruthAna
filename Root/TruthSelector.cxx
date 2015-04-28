#include "Mt2/mt2_bisect.h"
#include "TruthAna/TruthSelector.h"

#include <iomanip>
#include <string>

using namespace std;
using namespace Susy;

/*--------------------------------------------------------------------------------*/
// TruthSelector Constructor
/*--------------------------------------------------------------------------------*/
TruthSelector::TruthSelector() :
    m_truthNtupler(NULL),
    m_initialize(true) 
{
    cout<<"TruthSelector"<<endl;
}

/*--------------------------------------------------------------------------------*/
// Begin() - called at the start of the query
/*--------------------------------------------------------------------------------*/
void TruthSelector::Begin(TTree* /*tree*/)
{
}
void TruthSelector::SlaveBegin(TTree* /*tree*/)
{
    SusyNtTruthAna::Begin(0);
    
}
void TruthSelector::initializeNtupler()
{
    if(m_initialize==false) return; // check whether ntupler is already initialized  
    m_initialize=false;
    cout << "Initializing TruthNtupler" << endl;
    if(nt.evt()->isMC){
        TString mcid="";
        mcid.Form("%i",nt.evt()->mcChannel);
        m_truthNtupler = new TruthNtupler(mcid, TString("truthRazor"));
    }
    else {
        cerr << "TruthNtupler\t Need to intialize on MC! Exitting." << endl;
        abort();
    }
    m_truthNtupler->BookTree();
}  

/*--------------------------------------------------------------------------------*/
// Process - called at each event in the chain
/*--------------------------------------------------------------------------------*/
Bool_t TruthSelector::Process(Long64_t entry)
{
    GetEntry(entry); // c.f. SusyNtAna.h
    clearTruthObjects(); // c.f. SusyNtTruthAna.cxx
    m_chainEntry++;

    TruthSelector::initializeNtupler();

    if(m_dbg || m_chainEntry%5000==0){
        cout << "**** Processing entry " << setw(6) << m_chainEntry
             << " run " << setw(6) << nt.evt()->run
             << " event " << setw(6) << nt.evt()->event
             << " sumw= " << nt.evt()->sumw << " ****" << endl;
    }

    const Susy::TruthMet* truthMet = getTruthMet(&nt);
    TLorentzVector truthMet_TLV = truthMet->lv();
    TruthParticleVector truthElectrons = getPreTruthLeptons(&nt,11); // 10 GeV, |eta|<2.47
    TruthParticleVector truthMuons     = getPreTruthLeptons(&nt,13); // 10 GeV, |eta|<2.4
    TruthParticleVector truthTaus      = getPreTruthLeptons(&nt,15); // 20 GeV, |eta<2.5
    TruthJetVector      truthJets      = getPreTruthJets(&nt); // 20 GeV, |eta|<4.5
    
    // sort
    sort(truthElectrons.begin(), truthElectrons.end(), SortByPt());
    sort(truthMuons.begin(), truthMuons.end(), SortByPt());
    sort(truthTaus.begin(), truthTaus.end(), SortByPt());
    sort(truthJets.begin(), truthJets.end(), SortByPt());

    performOverlap(truthElectrons, truthMuons, truthTaus, truthJets);
    removeSFOSPair(truthElectrons);
    removeSFOSPair(truthMuons);

    m_cl20Idx.clear();
    m_f30Idx.clear();
    m_b20Idx.clear();
    TruthParticleVector signal_truthElectrons = truthElectrons;
    TruthParticleVector signal_truthMuons = truthMuons;
    TruthParticleVector signal_truthTaus = truthTaus;
    TruthJetVector signal_CL20truthJets = getTruthJetsCL20(truthJets); //jets: pt>20 GeV, |eta|<=2.4, not matched to b-quark
    TruthJetVector signal_F30truthJets = getTruthJetsF30(truthJets); //jets: pt>30 GeV, |eta|>2.4
    TruthJetVector signal_B20truthJets = getTruthJetsB20(truthJets); //jets: pt>20 GeV, |eta|<=2.4, matched to b-quark
    TruthJetVector signal_truthJets = getTruthSignalJets(truthJets); //jets: pt>20 GeV, |eta|<=2.4 OR pt>30 GeV, eta>2.5

    int cutflags=nt.evt()->cutFlags[NtSys_NOM];
    if(selectEvent(truthElectrons, truthMuons, truthTaus, truthJets, truthMet, cutflags)){
  /*      cout<<"Pass event selection"<<endl;
        cout<<"\t no. leptons: " << truthElectrons.size()+truthMuons.size()<<endl;
        cout<<"\t no. taus:    " << truthTaus.size() << endl;
        if(truthMuons.size()==2){
        cout<<"\t l0 pt:       " << truthMuons.at(0)->Pt() << " l1 pt:      " << truthMuons.at(1)->Pt() << endl;
        cout<<"\t sign:        " << truthMuons.at(0)->charge * truthMuons.at(1)->charge << endl;
        }  
        else if(truthElectrons.size()==2){
        cout<<"\t l0 pt:       " << truthElectrons.at(0)->Pt() << " l1 pt:     " << truthElectrons.at(1)->Pt()<<endl;
        cout<<"\t sign:        " << truthElectrons.at(0)->charge * truthElectrons.at(1)->charge <<endl;
        }
        else {
        cout<<"\t l0 pt:       " << truthElectrons.at(0)->Pt() << " l1 pt:     " << truthMuons.at(0)->Pt() << endl;
        cout<<"\t sign:        " << truthElectrons.at(0)->charge * truthMuons.at(0)->charge << endl;
        }
        for(uint ic=0; ic<signal_CL20truthJets.size(); ic++){
        cout<<"\t jCL20["<<ic<<"].pt: " << signal_CL20truthJets.at(ic)->Pt() << endl;
        }
        for(uint ib=0; ib<signal_B20truthJets.size();ib++){
        cout<<"\t jB20["<<ib<<"].pt: " << signal_B20truthJets.at(ib)->Pt() << endl;
        }
        for(uint i=0; i<signal_F30truthJets.size();i++){
        cout<<"\t jF30["<<i<<"].pt: " << signal_F30truthJets.at(i)->Pt() << endl;
        }
        for(uint is=0; is<signal_truthJets.size(); is++){
        cout<<"\t sigJet["<<is<<"].pt: " << signal_truthJets.at(is)->Pt() << endl;
        for(uint x=0; x<m_cl20Idx.size(); x++){
            if(m_cl20Idx[x]==is) cout<<"\t\t is CL20"<<endl;
        }
        for(uint x=0; x<m_f30Idx.size(); x++){
            if(m_f30Idx[x]==is) cout<<"\t\t is F30"<<endl;
        }
        for(uint x=0; x<m_b20Idx.size(); x++){
            if(m_b20Idx[x]==is) cout<<"\t\t is B20"<<endl;
        }
        }
  */      
        //////////////////////////////////
        // Fill output ntuple
        //////////////////////////////////

        // event var
        int run = nt.evt()->run;
        int eventno = nt.evt()->event;
        double weight = nt.evt()->w;
        int mcid = nt.evt()->mcChannel;
        m_truthNtupler->fillEventVars(run, eventno, mcid, weight);

        // lepton var
        m_truthNtupler->fillLeptonVars(truthElectrons, truthMuons, *truthMet);
        // met vars
        m_truthNtupler->fillMetVars(*truthMet);
        // jet vars
        m_truthNtupler->fillJetVars(signal_truthJets,signal_CL20truthJets, signal_F30truthJets, signal_B20truthJets,
                                    m_cl20Idx, m_b20Idx, m_f30Idx);

        // razor vars
        double mDeltaR=0.0;
        double shatr=0.0;
        double cosThetaRp1=0.0;
        double dphi_ll_vbeta_t=0.0;
        double dphi_l1_l2=0.0;
        double gamma_r=0.0;
        double dphi_vBeta_R_vBeta_T=0.0;
        double cosTheta_b=0.0;
        TVector3 vBeta_z;
        TVector3 pt_CM;
        TVector3 vBeta_T_CMtoR;
        TVector3 vBeta_R;
        TruthParticleVector leptons;
        for(uint ie=0;ie<truthElectrons.size();ie++){
            leptons.push_back(truthElectrons.at(ie));
        }
        for(uint im=0;im<truthMuons.size();im++){
            leptons.push_back(truthMuons.at(im));
        }
        sort(leptons.begin(), leptons.end(), SortByPt());
        TruthSelector::getSuperRazor(leptons, truthMet, vBeta_z, pt_CM,
                    vBeta_T_CMtoR, vBeta_R, shatr, dphi_ll_vbeta_t,
                    dphi_l1_l2, gamma_r, dphi_vBeta_R_vBeta_T,
                    mDeltaR, cosThetaRp1, cosTheta_b);
        
        double r2 = (truthMet->lv().Pt()) / (truthMet->lv().Pt() + leptons.at(0)->Pt() + leptons.at(1)->Pt());
        m_truthNtupler->fillRazorVars(mDeltaR, dphi_ll_vbeta_t, r2);
    
        

    
        m_truthNtupler->WriteTree();
        

    }
    return kTRUE;
    
}

/*--------------------------------------------------------------------------------*/
// Terminate methods
/*--------------------------------------------------------------------------------*/
void TruthSelector::Terminate()
{
}
void TruthSelector::SlaveTerminate()
{
    m_truthNtupler->setSumW(nt.evt()->sumw);
    m_truthNtupler->SaveTree();
}
/*--------------------------------------------------------------------------------*/
// Cleaning Methods
/*--------------------------------------------------------------------------------*/
void TruthSelector::performOverlap(TruthParticleVector& electrons, TruthParticleVector& muons,
                            TruthParticleVector& taus, TruthJetVector& jets)
{
    // remove electrons from electrons
    e_e_overlap(electrons, E_E_DR);
    // remove jets from electrons
    e_j_overlap(electrons, jets, J_E_DR,true);
    // remove taus from electrons
    t_e_overlap(taus, electrons, T_E_DR);
    // remove taus from muons
    t_m_overlap(taus, muons, T_M_DR);
    // remove electrons from jets
    e_j_overlap(electrons, jets, E_J_DR, false);
    // remove muons from jets
    m_j_overlap(muons, jets, M_J_DR);
    // remove electrons and muons that overlap
    e_m_overlap(electrons,muons, E_M_DR);
    // remove muons from muons
    m_m_overlap(muons, M_M_DR);
    // remove jets from taus
    t_j_overlap(taus, jets, J_T_DR, true);
}
void TruthSelector::e_e_overlap(TruthParticleVector& elecs, float minDr)
{
    uint nEl = elecs.size();
    if(nEl < 2) return;
    // Find all possible pairings
    static std::set<const TruthParticle*> removeElecs;
    removeElecs.clear();
    for(uint iEl=0; iEl<nEl; iEl++){
        const TruthParticle* e1 = elecs[iEl];
        for(uint jEl=iEl+1; jEl<nEl; jEl++){
            const TruthParticle* e2 = elecs[jEl];
            if(e1->DeltaR(*e2) < minDr){
                if(e1->Pt() < e2->Pt()){
                    removeElecs.insert(e1);
                    break;
                }else{
                    removeElecs.insert(e2);
                }
            } // dR
        } // e2 loop
    } // e1 loop

    // Remove electrons that overlap
    for(int iEl=nEl-1; iEl>=0; iEl--){
        if(removeElecs.find(elecs[iEl]) != removeElecs.end()){
            elecs.erase( elecs.begin() + iEl );
        }
    }
}
void TruthSelector::e_j_overlap(TruthParticleVector& elecs, TruthJetVector& jets, float minDr, bool removeJets)
{
//   cout << nt.evt()->event << "  removeJets= " << removeJets << " minDr= " << minDr << endl;
  if(elecs.size()==0 || jets.size()==0) return;

  for(int ie=elecs.size()-1; ie>=0; ie--){
    const TruthParticle* e = elecs.at(ie);
    for(int ij=jets.size()-1; ij>=0; ij--){
        const TruthJet* j = jets.at(ij);
        if(e->DeltaR(*j) > minDr) continue;       
        if(removeJets){
            jets.erase( jets.begin() + ij );
        }else{
            elecs.erase( elecs.begin() + ie );
            break;
        }
    }// end loop over jets
  }// end loop over electrons
}
void TruthSelector::m_j_overlap(TruthParticleVector& muons, TruthJetVector& jets, float minDr)
{
  if(muons.size()==0 || jets.size()==0) return;
  for(int im=muons.size()-1; im>=0; im--){
      const TruthParticle* mu = muons.at(im);
      for(int ij=jets.size()-1; ij>=0; ij--){
          const TruthJet* j = jets.at(ij);
          if(mu->DeltaR(*j) > minDr) continue;
          muons.erase( muons.begin() + im );
          break;
      }// end loop over jets
  }// end loop over muons
}
void TruthSelector::e_m_overlap(TruthParticleVector& elecs,
                                             TruthParticleVector& muons,
                                             float minDr)
{
    uint nEl = elecs.size();
    uint nMu = muons.size();
    if(nEl==0 || nMu==0) return;
    // Electron muon overlap should be pretty rare,
    // so we can take advantage of that and optimize
    static std::set<const TruthParticle*> removeElecs;
    static std::set<const TruthParticle*> removeMuons;
    removeElecs.clear();
    removeMuons.clear();
    // In this case we will want to remove both the electron and the muon
    for(uint iEl=0; iEl<nEl; iEl++){
        const TruthParticle* e = elecs[iEl];
        for(uint iMu=0; iMu<nMu; iMu++){
            const TruthParticle* mu = muons[iMu];
            if(e->DeltaR(*mu) < minDr){
                removeElecs.insert(e);
                removeMuons.insert(mu);
            }
        }
    }
    // Remove those electrons flagged for removal
    if(removeElecs.size()){
        for(int iEl=nEl-1; iEl>=0; iEl--){
            if(removeElecs.find(elecs[iEl])!=removeElecs.end()){
                elecs.erase( elecs.begin() + iEl );
            }
        }
    }
    // Remove those muons flagged for removal
    if(removeMuons.size()){
        for(int iMu=nMu-1; iMu>=0; iMu--){
            if(removeMuons.find(muons[iMu])!=removeMuons.end()){
                muons.erase( muons.begin() + iMu );
            }
        }
    }
}
void TruthSelector::m_m_overlap(TruthParticleVector& muons, float minDr)
{
    uint nMu = muons.size();
    if(nMu < 2) return;
    // If 2 muons overlap, toss them both!
    static std::set<const TruthParticle*> removeMuons;
    removeMuons.clear();
    for(uint iMu=0; iMu<nMu; iMu++){
        const TruthParticle* mu1 = muons[iMu];
        for(uint jMu=iMu+1; jMu<nMu; jMu++){
            const TruthParticle* mu2 = muons[jMu];
            if(mu1->DeltaR(*mu2) < minDr){
                removeMuons.insert(mu1);
                removeMuons.insert(mu2);
            }
        }
    }
    for(int iMu=nMu-1; iMu>=0; iMu--){
        const TruthParticle* mu = muons[iMu];
        if(removeMuons.find(mu) != removeMuons.end()){
            muons.erase( muons.begin() + iMu );
        }
    }
}
void TruthSelector::t_e_overlap(TruthParticleVector& taus,
                                             TruthParticleVector& elecs,
                                             float minDr)
{
    uint nTau = taus.size();
    uint nEle = elecs.size();
    if(nTau==0 || nEle==0) return;    
    for(int iTau=nTau-1; iTau>=0; iTau--){
        const TruthParticle* tau = taus[iTau];
        for(int iEl=nEle-1; iEl>=0; iEl--){
            const TruthParticle* e = elecs[iEl];
            if(tau->DeltaR(*e) < minDr){
                taus.erase( taus.begin() + iTau );
                break;
            }            
        }
    }
}
void TruthSelector::t_m_overlap(TruthParticleVector& taus,
                                             TruthParticleVector& muons,
                                             float minDr)
{
    uint nTau = taus.size();
    uint nMuo = muons.size();
    if(nTau==0 || nMuo==0) return;
    for(int iTau=nTau-1; iTau>=0; iTau--){
        const TruthParticle* tau = taus[iTau];
        for(int iMu=nMuo-1; iMu>=0; iMu--){
            const TruthParticle* mu = muons[iMu];            
            if(tau->DeltaR(*mu) < minDr){
                taus.erase( taus.begin() + iTau );
                break;
            }            
        }
    }
}
void TruthSelector::t_j_overlap(TruthParticleVector& taus,
                                             TruthJetVector& jets,
                                             float minDr,
                                             bool removeJets)
{
    uint nTau = taus.size();
    uint nJet = jets.size();
    if(nTau==0 || nJet==0) return;    
    for(int iTau=nTau-1; iTau>=0; iTau--){
        const TruthParticle* tau = taus.at(iTau);
        for(int iJet=jets.size()-1; iJet>=0; iJet--){
            const TruthJet* jet = jets.at(iJet);            
            if(tau->DeltaR(*jet) < minDr){
                if(removeJets)
                    jets.erase( jets.begin() + iJet );
                else{
                    taus.erase( taus.begin() + iTau );
                    break;
                }
            }
        } // end loop over jets
    } // end loop over electrons
}
/*--------------------------------------------------------------------------------*/
// Building methods
/*--------------------------------------------------------------------------------*/
TruthJetVector TruthSelector::getTruthJetsCL20(TruthJetVector& truthBaseJets)
{
    TruthJetVector truthSignalJets;
    for(uint i=0;i<truthBaseJets.size(); ++i){
        TruthJet* jet = truthBaseJets.at(i);
        bool bmatch = false;
        TruthJet* truthjet = truthBaseJets.at(i);
        for(uint ij=0;ij<nt.jet()->size();++ij){
            Jet* recoj = & nt.jet()->at(ij);
            if(recoj->truthLabel==5){
                if(recoj->DeltaR(*truthjet)<0.4) bmatch=true;
            }
        }
        if(!bmatch && jet->Pt()>=20. && fabs(jet->Eta())<=2.4) {
            truthSignalJets.push_back(jet);
            m_cl20Idx.push_back(i);
        }
    }
    return truthSignalJets;
}
TruthJetVector TruthSelector::getTruthJetsF30(TruthJetVector& truthBaseJets)
{
    TruthJetVector truthSignalJets;
    for(uint i=0; i<truthBaseJets.size();++i){
        TruthJet* jet=truthBaseJets.at(i);
        if(jet->Pt()>=30. && fabs(jet->Eta())>2.4){
            truthSignalJets.push_back(jet);
            m_f30Idx.push_back(i);
        }
    }
    return truthSignalJets;
}
TruthJetVector TruthSelector::getTruthJetsB20(TruthJetVector& truthBaseJets)
{
    TruthJetVector truthSignalJets;
    for(uint i=0; i<truthBaseJets.size(); ++i){
        bool bmatch=false;
        TruthJet* truthjet = truthBaseJets.at(i);
        for(uint ij=0; ij<nt.jet()->size();++ij){
            Jet* recoj = &nt.jet()->at(ij);
            if(recoj->truthLabel==5){
                if(recoj->DeltaR(*truthjet)<0.4) bmatch=true;
            }
        }
        if(bmatch && truthjet->Pt()>20. && fabs(truthjet->Eta())<=2.4){ 
            truthSignalJets.push_back(truthjet);
            m_b20Idx.push_back(i);
        }
    }
    return truthSignalJets;
}
TruthJetVector TruthSelector::getTruthSignalJets(TruthJetVector& truthBaseJets)
{
    TruthJetVector truthSignalJets;
    for(uint i=0; i<truthBaseJets.size(); ++i){
        TruthJet* jet = truthBaseJets.at(i);
        if(jet->Pt()>=20. && fabs(jet->Eta())<=2.4) truthSignalJets.push_back(jet);
        if(jet->Pt()>=30. && fabs(jet->Eta())>2.4)  truthSignalJets.push_back(jet);
    }
    return truthSignalJets;
}

/*--------------------------------------------------------------------------------*/
bool TruthSelector::selectEvent(const TruthParticleVector& electrons, const TruthParticleVector& muons,
                        const TruthParticleVector& taus, const TruthJetVector& jets, const TruthMet* met, int cutflags)
{
    if( !(cutflags & ECut_GRL) ) return false;
    if( !(cutflags & ECut_TileTrip) ) return false;
    if( !(cutflags & ECut_LarErr) ) return false;

    if(nt.evt()->hfor==4) return false;
    
    if(!(muons.size()==2 || electrons.size()==2 ||
            (electrons.size()+muons.size())==2)) return false;
    
    if(muons.size()==2 && SusyNtAna::Mll(muons[0], muons[1])<20) return false;
    if(electrons.size()==2 && SusyNtAna::Mll(electrons[0], electrons[1])<20) return false;
    if((electrons.size()==1 && muons.size()==1) &&
            SusyNtAna::Mll(electrons[0], muons[0])<20) return false;
    if(!taus.size()==0) return false;

    if(muons.size()==2 && ( (muons[0]->charge * muons[1]->charge)>0)) return false;
    if(electrons.size()==2 && ( (electrons[0]->charge * electrons[1]->charge)>0)) return false;
    if(electrons.size()==1 && muons.size()==1 && 
            (electrons.size()+muons.size())==2 && ((electrons[0]->charge * muons[0]->charge)>0)) return false;
    
    return true;
}
void TruthSelector::getSuperRazor(const TruthParticleVector& leptons, const TruthMet* met,
			     TVector3& vBETA_z, TVector3& pT_CM,
			     TVector3& vBETA_T_CMtoR, TVector3& vBETA_R,
			     double& SHATR, double& dphi_LL_vBETA_T, double& dphi_L1_L2,
			     double& gamma_R, double&  dphi_vBETA_R_vBETA_T,
			     double& MDELTAR, double& costhetaRp1,
                             double& cosTheta_b)
{
// MDR CALCULATION 
//
// Code written by Christopher Rogan <crogan@cern.ch>, 04-23-13
// Details given in paper (http://arxiv.org/abs/1310.4827) written by 
// Matthew R. Buckley, Joseph D. Lykken, Christopher Rogan, Maria Spiropulu
//
  if( leptons.size() < 2 ) return;
  
  // necessary variables
  TLorentzVector metlv = met->lv();
  TLorentzVector l0    = *leptons.at(0);
  TLorentzVector l1    = *leptons.at(1);

  //
  // Lab frame
  //
  //Longitudinal boost
  vBETA_z = (1./(l0.E()+l1.E()))*(l0+l1).Vect(); 
  vBETA_z.SetX(0.0);         
  vBETA_z.SetY(0.0);
  
  l0.Boost(-vBETA_z);
  l1.Boost(-vBETA_z);

  //pT of CM frame
  pT_CM = (l0+l1).Vect() + metlv.Vect();
  pT_CM.SetZ(0.0);     
  
  TLorentzVector ll = l0+l1;
  //invariant mass of the total event
  SHATR = sqrt( 2.*(ll.E()*ll.E() - ll.Vect().Dot(pT_CM) 
		   + ll.E()*sqrt( ll.E()*ll.E() + pT_CM.Mag2() - 2.*ll.Vect().Dot(pT_CM) )));
  
  vBETA_T_CMtoR = (1./sqrt(pT_CM.Mag2() + SHATR*SHATR))*pT_CM;

  // FAR_WW variable cosTheta_b
  cosTheta_b = l0.CosTheta();
  
  l0.Boost(-vBETA_T_CMtoR);
  l1.Boost(-vBETA_T_CMtoR);
  ll.Boost(-vBETA_T_CMtoR);  

  //
  //R-frame
  //
  dphi_LL_vBETA_T = fabs((ll.Vect()).DeltaPhi(vBETA_T_CMtoR));
  
  dphi_L1_L2 = fabs(l0.Vect().DeltaPhi(l1.Vect()));
  
  vBETA_R = (1./(l0.E()+l1.E()))*(l0.Vect() - l1.Vect());
  
  gamma_R = 1./sqrt(1.-vBETA_R.Mag2());
  
  dphi_vBETA_R_vBETA_T = fabs(vBETA_R.DeltaPhi(vBETA_T_CMtoR));
  
  l0.Boost(-vBETA_R);
  l1.Boost(vBETA_R);
 
  //
  //R+1 frame
  //
  MDELTAR = 2.*l0.E();
  costhetaRp1 = l0.Vect().Dot(vBETA_R)/(l0.Vect().Mag()*vBETA_R.Mag());

  return;
}

