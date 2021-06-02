#include "TrackInfo.h"
#include <iostream>

TrackInfo::TrackInfo(const char* rootfile_name)
{

  //root_tfile= new TFile(rootfile_name,"read");
    root_tfile= new TFile("anaRun0659.root","read");
   
    tracks_class = new tracks((TTree*)root_tfile->Get("tracks"));

    std::cout << "Opening Tracks rootfile: " << root_tfile << std::endl;

    NofEntries = (tracks_class->fChain)->GetEntries();
    tracks_class->GetEntry(NofEntries-1);
    max_trigger_number = tracks_class->srstriggerctr +1;
    std::cout << "step2  " << max_trigger_number <<std::endl;

    TrackBank = new std::vector< Track >[max_trigger_number];
    NumberofTracks = new int[max_trigger_number];
    std::cout << "Number of tracks in rootfile: " << NofEntries << std::endl << std::endl;

    FillTracks();
    root_tfile->Close();
    delete tracks_class;
    delete root_tfile;
}

TrackInfo::~TrackInfo()
{
    delete[] TrackBank;
    delete[] NumberofTracks;
}

void TrackInfo::FillTracks()
{
    int previous_trig = 0;
    for(int i = 0; i < NofEntries; ++i)
    {
        tracks_class->GetEntry(i);
        int trignum = tracks_class->srstriggerctr;
        NumberofTracks[trignum] = tracks_class->ntracks;
        if(previous_trig > trignum)
        {
            std::cout << "****WARNING**** Tracks in rootfile not in an increasing trigger numer order" << std::endl;
            std::cout << "Some tracks are bound to be lost" << std::endl << std::endl;
        }

        for(int j = previous_trig + 1; j < trignum; ++j)
            TrackBank[j].clear();

        Track OneTrack;
        OneTrack.Chi2 = tracks_class->trackchi2;
        //old: 6-MCP,7-MCP,8-MM ; 6-Multipad, 7-Resistive
        OneTrack.Hit1 = tracks_class->hits->at(6);
	//	std::cout << "tracknr " << previous_trig << " Hit X " << OneTrack.Hit1X <<" Hit y " << OneTrack.Hit1Y <<std::endl;

		OneTrack.Hit2 = tracks_class->srstriggerctr;
		// OneTrack.Hit3 = tracks_class->hits->at(8);
        OneTrack.SlopeXZ = (tracks_class->hits->at(0).at(0)-tracks_class->hits->at(4).at(0))/(tracks_class->hits->at(0).at(2)-tracks_class->hits->at(4).at(2));
        
        OneTrack.SlopeYZ = (tracks_class->hits->at(1).at(1)-tracks_class->hits->at(5).at(1))/(tracks_class->hits->at(1).at(2)-tracks_class->hits->at(5).at(2));
        TrackBank[trignum].push_back(OneTrack);
//        std::cout<<"OneTrack.Hit="<<tracks_class->hits->at(8).at(1)<<std::endl;
//         std::cout<<"OneTrack.SRSnumber="<<tracks_class->srstriggerctr<<std::endl;
//        OneTrack.SRSnumber = tracks_class->srstriggerctr;
//                 std::cout<<"OneTrack.SRSnumber="<<OneTrack.SRSnumber<<std::endl;

        
        previous_trig = trignum;

    }

    int counter = 0;
    std::cout << "=== List of Triggers with irregular number of tracks (not equal to 1)=======" << std::endl;
    for(int i= 0; i < max_trigger_number; ++i) 
    {
        NumberofTracks[i] = TrackBank[i].size();
        if(NumberofTracks[i] != 1)
        {
            std::cout << i << " \t";
            counter++;
            if(counter%10 == 0) std::cout << std::endl;
            if(NumberofTracks[i] == 0)
            {
                Track OneTrack;
                std::vector<double> dummy_hit;
                dummy_hit.push_back(1000);
                dummy_hit.push_back(1000);
                dummy_hit.push_back(1000);


	    
		OneTrack.Hit1 = dummy_hit;
                OneTrack.Hit2 = 1000000;
		//	OneTrack.Hit3 = dummy_hit;

                TrackBank[i].push_back(OneTrack);
            }
        }
    }
    std::cout << std::endl << "============================================================================" << std::endl << std::endl;
}

