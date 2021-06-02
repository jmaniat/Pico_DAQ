#ifndef __trackinfo_h__
#define __trackinfo_h__

#include <vector>

#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>
#include <TROOT.h>
#include <TDirectory.h>

#include "tracks.h"

struct Track
{
  double Chi2=-1.;
   std::vector<double> Hit1;
double Hit2;
  //std::vector<double> Hit2Y;
  //std::vector<double> Hit3X;
  //std::vector<double> Hit3Y;
  double SlopeXZ = 0.;
  double SlopeYZ = 0.;
    double SRSnumber = 0.;
};

class TrackInfo
{
    public:
        TrackInfo(const char* rootfilename);        
        ~TrackInfo();
        std::vector< Track > *TrackBank;
        int *NumberofTracks;
    private:
        void FillTracks();
        int max_trigger_number;
        int NofEntries;
        tracks *tracks_class;
        TFile* root_tfile;
};

#endif
