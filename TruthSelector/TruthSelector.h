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
        virtual void Process(Long64_t entry);
        virtual void Terminate();
        virtual void SlaveTerminate(); 


}; // class TruthSelector



#endif
