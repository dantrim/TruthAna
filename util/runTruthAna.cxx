#include "TruthAna/TruthSelector.h"
#include "SusyNtuple/ChainHelper.h"
#include "TChain.h"
#include "Cintex/Cintex.h"
#include <iostream>
#include <string>

using namespace std;


//----------------------------------------------------------
// run TruthSelector
//----------------------------------------------------------
void printHelp(const char *exeName)
{
    cout<<"Usage :"<<endl
        <<exeName<<endl
        <<"\t -i input (file,list, or dir)"<<endl
        <<"\t -n number of events to process"<<endl
        <<"\t -s samplename"<<endl
        <<"\t -v set verbose"<<endl
        <<"\t -h print this help"<<endl
        <<endl;
}

int main(int argc, char **argv)
{
    ROOT::Cintex::Cintex::Enable();
    
    int dbg=0;
    int nEvt=-1;
    string sample;
    string input;
    

    for(int i = 1; i < argc; i++) {
        if      (strcmp(argv[i], "-d")==0) dbg=atoi(argv[++i]);
        else if (strcmp(argv[i], "-i")==0) input=argv[++i];
        else if (strcmp(argv[i], "-s")==0) sample=argv[++i];
        else if (strcmp(argv[i], "-n")==0) nEvt=atoi(argv[++i]);
        else {
            printHelp(argv[0]);
            return 0;
        }
    }
    if(input.empty()){
        cout << "No input specified. Seek help." << endl;
        printHelp(argv[0]);
        return 1;
    }
    
    TChain* chain = new TChain("susyNt");   
    ChainHelper::addInput(chain, input, dbg>0);
    Long64_t nEntries = chain->GetEntries();
    chain->ls();
    
    TruthSelector* susyAna = new TruthSelector();
    susyAna->setDebug(dbg);
    if(nEvt<0) nEvt=nEntries;
    cout<<endl;
    cout<<"Total entries:\t"<<nEntries<<endl;
    cout<<"Processing:\t"<<nEvt<<endl;
    if(nEvt>0) chain->Process(susyAna,sample.c_str(),nEvt);
    cout<<endl;
    cout<<" ------------------------- "<<endl;
    cout<<" TruthSelector job complete"<<endl;
    cout<<" ------------------------- "<<endl;

    delete chain;
    return 0;
}
    
    
    
