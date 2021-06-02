#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <chrono>

#include <TBenchmark.h>

#include "TRC_FileReader.h"
#include "Detector.h"
#include "TestBeamSetup.h"
#include "TrackInfo.h"

using namespace std::chrono;
using namespace std;

int main(int argc, char* argv[])
{
  auto start = high_resolution_clock::now();

    //1st input argument is is directory path with TRC files
    //2nd input argument is number of files for each channel in the directory
    //3rd input argument is name the output root file will have
    //4th input argument (optional), if it is 1 then dumping will be enabled

    TBenchmark bench;
    bench.Start("full");

    TrackInfo *Tracker;

    int enable_dumping = 0;
    if(argc > 4)
    {
      // std::cout << "step1: dump" << std::endl;
        if(atoi(argv[4]) == 1)
            enable_dumping = 1;
    }

    //Specify the channels written in the files
    std::vector<int> channel_IDs;
//    channel_IDs.push_back(1); //<- Channel 1 (mcp)
    channel_IDs.push_back(2); //<- Channel 2 (mm)
     channel_IDs.push_back(3);
//    channel_IDs.push_back(4);
    //channel_IDs.push_back(5);
    // channel_IDs.push_back(6);
    // channel_IDs.push_back(7);


//    TRC_FileReader myfile(channel_IDs,argv[1],atoi(argv[2]), "Trace");
     TRC_FileReader myfile(channel_IDs,argv[1],atoi(argv[2]), "--Trace--",2);
//      myfile.OpenTriggerChannel(argv[1],atoi(argv[2]), "--Trace--", 3);//channel 3


    
    TestBeamSetup mysetup;

//    mysetup.CreateMCP();
    // mysetup.CreateMCP();
//     mysetup.CreateMM();
     //	 mysetup.CreateMCP();
     //	 mysetup.CreateMCP();
    	  mysetup.CreateMM();
          mysetup.CreateMM();


    
    int NumberofTracks;
    Track OneTrack;
    if( atoi(argv[5]) != 0)
    {
        Tracker = new TrackInfo(argv[6]);
        mysetup.init(NumberofTracks, OneTrack);
    }
    else
   
      mysetup.init();

    myfile.SetDetectorSetup(mysetup);

    const char* output_rootfile_name = argv[3];
    
    int dumping_id = 0;

    auto stopread = high_resolution_clock::now(); 
    auto durationread = duration_cast<milliseconds>(stopread - start);
    double durread=durationread.count();
//    cout<<setprecision(4) << "It took "<<durread/1000 << " s read the files"<<endl;
    int i = 0;
    while( myfile.GetNextEvent() )
    {
        i++;
        printf("\rEvent is: %.5d ===============================================",i-1);//magic
	

	//	std::cout << "trignum = " << mysetup.TriggerNumber << std::endl;
        if( atoi(argv[5]) != 0 )
        {
	
	NumberofTracks = Tracker->NumberofTracks[mysetup.TriggerNumber];
	OneTrack = Tracker->TrackBank[mysetup.TriggerNumber].at(0);
	

        }

        mysetup.TestBeamAnalysis();
        

   if(i==1){
	auto stopcalc = high_resolution_clock::now(); 
	auto durationcalc = duration_cast<milliseconds>(stopcalc - stopread);
	double durcalc=durationcalc.count();
	      }
   
        if(i-1==dumping_id) // Specify Event to be dumped to dumpfile.root (read by executing ` root -l draw_dump.C `)
            mysetup.Dump();

        if(dumping_id >=0 && enable_dumping) 
            if(dumping_id<=i-1)
            {
                printf("\nNext event to dump: (give negative to stop dumping)");
                int ok = scanf("%d",&dumping_id);
                if(dumping_id < 0) break;
            }

        mysetup.Fill_Tree();
//        if(i==1) break;
    }
    mysetup.Finalize_Tree(output_rootfile_name);
    mysetup.WriteFile();
    bench.Show("full");
    
}

