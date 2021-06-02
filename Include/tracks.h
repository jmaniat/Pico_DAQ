//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 13 18:55:00 2017 by ROOT version 6.08/04
// from TTree tracks/contains info on tracks in 3d space
// found on file: Data/Tracks/anaRun0209.root
//////////////////////////////////////////////////////////

#ifndef tracks_h
#define tracks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class tracks {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           srstriggerctr;
   Int_t           srstimestamp;
   Int_t           ntracks;
   Int_t           tracknumber;
   Double_t        trackchi2;
   std::vector<std::vector<double> > *hits;
   std::vector<double>  *distnextcluster;
   std::vector<double>  *totchanextcluster;

   // List of branches
   TBranch        *b_srstriggerctr;   //!
   TBranch        *b_srstimestamp;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_tracknumber;   //!
   TBranch        *b_trackchi2;   //!
   TBranch        *b_hits;   //!
   TBranch        *b_distnextcluster;   //!
   TBranch        *b_totchanextcluster;   //!

   tracks(TTree *tree=0);
   virtual ~tracks();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tracks_cxx
tracks::tracks(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("tracks.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("tracks.root");
      }
      f->GetObject("tracks",tree);

   }
   Init(tree);
}

tracks::~tracks()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tracks::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tracks::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tracks::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   hits = 0;
   distnextcluster = 0;
   totchanextcluster = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("srstriggerctr", &srstriggerctr, &b_srstriggerctr);
   fChain->SetBranchAddress("srstimestamp", &srstimestamp, &b_srstimestamp);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("tracknumber", &tracknumber, &b_tracknumber);
   fChain->SetBranchAddress("trackchi2", &trackchi2, &b_trackchi2);
   fChain->SetBranchAddress("hits", &hits, &b_hits);
   fChain->SetBranchAddress("distnextcluster", &distnextcluster, &b_distnextcluster);
   fChain->SetBranchAddress("totchanextcluster", &totchanextcluster, &b_totchanextcluster);
   Notify();
}

Bool_t tracks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tracks::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tracks::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tracks_cxx
