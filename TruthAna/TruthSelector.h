#ifndef truthRazor_TruthSelector_h
#define truthRazor_TruthSelector_h

//////////////////////////////////////////////////////
// I have realized that I do not have a truth       //
// analysis script.                                 //
//////////////////////////////////////////////////////

#include "SusyNtuple/SusyNtTruthAna.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/DilTrigLogic.h"
#include "SusyNtuple/SusyDefs.h"

#include <fstream>

// ROOT
#include "TTree.h"


using Susy::Lepton;
using Susy::Jet;
using Susy::Met;

class TruthSelector : public SusyNtTruthAna
{
    public :
        TruthSelector();
        virtual ~TruthSelector(){};
        
        // TSelector methods (inheriting TSelector through SusyNtTruthAna)
        virtual void Begin(TTree *tree);
        virtual void SlaveBegin(TTree *tree);
        virtual Bool_t Process(Long64_t entry);
        virtual void Terminate();
        virtual void SlaveTerminate(); 

        virtual void performOverlap(TruthParticleVector& ele, TruthParticleVector& mou,
                            TruthParticleVector& tau, TruthJetVector& jets);
        virtual void e_e_overlap(TruthParticleVector& electrons, float minDr);
        virtual void e_j_overlap(TruthParticleVector& electrons, TruthJetVector& jets, float minDr, bool removeJets);
        virtual void t_e_overlap(TruthParticleVector& taus, TruthParticleVector& electrons, float minDr);
        virtual void t_m_overlap(TruthParticleVector& taus, TruthParticleVector& muons, float minDr);
        virtual void m_j_overlap(TruthParticleVector& muons, TruthJetVector& jets, float minDr);
        virtual void e_m_overlap(TruthParticleVector& electrons, TruthParticleVector& muons, float minDr);
        virtual void m_m_overlap(TruthParticleVector& muons, float minDr);
        virtual void t_j_overlap(TruthParticleVector& taus, TruthJetVector& jets, float minDr, bool removeJets);

        TruthJetVector getTruthJetsCL20(TruthJetVector& truthBaseJets);
        TruthJetVector getTruthJetsF30(TruthJetVector& truthBaseJets);
        TruthJetVector getTruthJetsB20(TruthJetVector& truthBaseJets);

        bool selectEvent(const TruthParticleVector& baseElectrons, const TruthParticleVector& baseMuons,
                    const TruthParticleVector& baseTaus, const TruthJetVector& jets, const TruthMet* met, int cutflags);
                        

        ClassDef(TruthSelector,1);


}; // class TruthSelector



#endif
