#ifndef __detectorsetup_h__
#define __detectorsetup_h__

#include <string>
#include <vector>

#include <TFile.h>
#include <TVectorD.h>
#include <TObject.h>
#include <TMath.h>
#include <TTree.h>

#include "AverageTool.h"
#include "Detector.h"
#include "TrackInfo.h"
class Detector;

class DetectorSetup
{
    public:
        ~DetectorSetup();
        void CreateMCP();
        void CreateMM();
        void SetDetectorWaveform(int i, std::vector<double> x, std::vector<double> y);
        void SetTriggerNumber(int TrigNum);
        //void DumpOld();
        void Fill_Tree();
        void Finalize_Tree(const char* outrootfile);
        void ScaleAndShiftTimes();
        void SetWaveformToAverage(AverageTool &aver);
        int TriggerNumber;
    protected:
        std::vector<Detector*> Detectors;
        int baseline_region_end;
        int max_region_end;
        int NofDetectors;

        TFile* OutFile;//("outrootfile.root","recreate");
        TTree* OutTree;//("pico","Analysis Output");
};
#endif
