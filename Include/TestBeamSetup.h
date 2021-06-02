#ifndef __lasersetup_h__
#define __lasersetup_h__

#include "DetectorSetup.h"
#include "TrackInfo.h"


class TestBeamSetup : public DetectorSetup
{
    public:
        void WriteFile();
        void TestBeamAnalysis();
        void SetWaveformToAverage(AverageTool &aver);
        void Dump();
        void init();
        void init_tree();
        void tree_add_tracks(int& NumberofTracks, Track& OneTrack);
        void init(int& NumberofTracks, Track& OneTrack);

};
#endif
