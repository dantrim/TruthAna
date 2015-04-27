#include "Mt2/mt2_bisect.h"
#include "TruthSelector/TruthSelector.h"

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
