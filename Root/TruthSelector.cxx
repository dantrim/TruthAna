#include "Mt2/mt2_bisect.h"
#include "TruthAna/TruthSelector.h"

#include <iomanip>
#include <string>

using namespace std;
using namespace Susy;

/*--------------------------------------------------------------------------------*/
// TruthSelector Constructor
/*--------------------------------------------------------------------------------*/
TruthSelector::TruthSelector()
{
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
/*--------------------------------------------------------------------------------*/
// Process - called at each event in the chain
/*--------------------------------------------------------------------------------*/
Bool_t TruthSelector::Process(Long64_t entry)
{
    GetEntry(entry); // c.f. SusyNtAna.h
    clearTruthObjects(); // c.f. SusyNtTruthAna.cxx
    m_chainEntry++;

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
    performOverlap(truthElectrons, truthMuons, truthTaus, truthJets);
    removeSFOSPair(truthElectrons);
    removeSFOSPair(truthMuons);

    TruthParticleVector signal_truthElectrons = truthElectrons;
    TruthParticleVector signal_truthMuons = truthMuons;
    TruthParticleVector signal_truthTaus = truthTaus;
    TruthJetVector signal_CL20truthJets = getTruthJetsCL20(truthJets); //jets: pt>20 GeV, |eta|<=2.4, not matched to b-quark
    TruthJetVector signal_F30truthJets = getTruthJetsF30(truthJets); //jets: pt>30 GeV, |eta|>2.4
    TruthJetVector signal_B20truthJets = getTruthJetsB20(truthJets); //jets: pt>20 GeV, |eta|<=2.4, matched to b-quark

    int cutflags=nt.evt()->cutFlags[NtSys_NOM];
    if(selectEvent(truthElectrons, truthMuons, truthTaus, truthJets, truthMet, cutflags)) cout<<"Event clean"<<endl;

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
        if(!bmatch && jet->Pt()>=20. && fabs(jet->Eta())<=2.4) truthSignalJets.push_back(jet);
    }
    return truthSignalJets;
}
TruthJetVector TruthSelector::getTruthJetsF30(TruthJetVector& truthBaseJets)
{
    TruthJetVector truthSignalJets;
    for(uint i=0; i<truthBaseJets.size();++i){
        TruthJet* jet=truthBaseJets.at(i);
        if(jet->Pt()>=30. && fabs(jet->Eta())>2.4) truthSignalJets.push_back(jet);
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
        if(bmatch && truthjet->Pt()>20. && fabs(truthjet->Eta())<=2.4) 
                truthSignalJets.push_back(truthjet);
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
    if(electrons.size()==1 && muons.size()==1 && ((electrons[0]->charge * muons[0]->charge)>0)) return false;
    
    return true;
}
    


