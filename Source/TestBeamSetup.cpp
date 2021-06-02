#include "TestBeamSetup.h"
#include "TRC_FileReader.h"
#include <iostream>


void TestBeamSetup::WriteFile(){
std::ofstream tot_val;
tot_val.open("./s.txt");
    
    std::map<std::string, std::vector<float>> output_map;
    output_map = Detectors.at(0)->GetToT();
    
    std::cout << output_map["first_0"].size()  << std::endl;
    
    for(int i = 0; i < output_map["first_0"].size(); i++){
        tot_val << output_map["first_0"].at(i) << "\t" << output_map["second_0"].at(i) << "\t" << output_map["first_1"].at(i) << "\t" << output_map["second_1"].at(i) << "\t" << output_map["first_2"].at(i) << "\t" << output_map["second_2"].at(i) << "\t" << output_map["first_3"].at(i) << "\t" << output_map["second_3"].at(i) << "\t" << output_map["first_4"].at(i) << "\t" << output_map["second_4"].at(i) << "\t" << output_map["first_5"].at(i) << "\t" << output_map["second_5"].at(i) << std::endl;
    }

    
}

void TestBeamSetup::TestBeamAnalysis()
{

    ScaleAndShiftTimes();

    for(int i = 0; i < NofDetectors; ++i)
    {
           if(i == 0) Detectors.at(i)->InvertY();
    }
    //Dynamically find the end for the baseline calculation region
    baseline_region_end=2000000;
    
    for(int i = 0; i < NofDetectors; ++i)
    {
        Detectors.at(i)->FindGlobalMaximum(0, Detectors.at(i)->waveform_y.size()-1);
        if(Detectors.at(i)->global_maximum.position < baseline_region_end)
            baseline_region_end = Detectors.at(i)->global_maximum.position;
    }
    baseline_region_end-=200;
    if(baseline_region_end < 0 )
        baseline_region_end=100.;

    max_region_end = 2000+baseline_region_end;
    //===

    for(int i = 0; i < NofDetectors; ++i)
    {
        Detectors.at(i)->SubstractBaseline(baseline_region_end);

	// Detectors.at(i)->FindGlobalMaximum(baseline_region_end, max_region_end);
	//Detectors.at(i)->TimeSigmoid();
	Detectors.at(i)->FindGlobalMaximum(0,20000);
	Detectors.at(i)->FindStartPoint(baseline_region_end);
	Detectors.at(i)->FindEndPoint(max_region_end);
        Detectors.at(i)->FindElectronPeakEndPoint();
	 Detectors.at(i)->CalculateCharges();
	Detectors.at(i)->FindNaiveTiming();
	// Detectors.at(i)->FindRiseTime();
	// Detectors.at(i)->FindFirstPeak();
        //Detectors.at(i)->FindMaxDerivative();
	 //	Detectors.at(i)->FindCFDTiming();
	 //	Detectors.at(i)->Findmultiplepeak();
	 //Detectors.at(i)->Findmultipleinverse();
        Detectors.at(i)->TimeSigmoid();

        if(Detectors.at(i)->type > 0){
            if(i == 0){
                Detectors.at(i)->ToT();
                Detectors.at(i)->ExpoFit();
//                std::cout << "detector: " << i <<  " max = " << Detectors.at(i)->global_maximum.position << " peak point = " << std::endl;
            }
    Detectors.at(i)->FullSigmoidFit();


//            std::cin.ignore();
        }
            //	    std::cout<<"Mark 1"<<std::endl;


        if(Detectors.at(i)->type > 0)
        {
            
	  // Detectors.at(i)->TimeSigmoid();
	  // Detectors.at(i)->TimeSigmoid();
          //  Detectors.at(i)->FilterWaveformFFT(baseline_region_end, 8192, 2.5);
          //  Detectors.at(i)->FindGlobalMaximum(0, 8191);
	  // Detectors.at(i)->FindGlobalMaximum(0, 20000);

	  // Detectors.at(i)->TimeSigmoid();
	  // Detectors.at(i)->FindStartPoint(0);
	  // Detectors.at(i)->FindEndPoint(max_region_end-baseline_region_end);
	  // Detectors.at(i)->FindElectronPeakEndPoint();
	  // Detectors.at(i)->CalculateCharges();
	  // Detectors.at(i)->FindNaiveTiming();
	    //Detectors.at(i)->TimeSigmoid();
            //Detectors.at(i)->FindRiseTime();
            //Detectors.at(i)->FindFirstPeak();
            //Detectors.at(i)->FindMaxDerivative();
	  //  Detectors.at(i)->TimeSigmoid();
	  // Detectors.at(i)->FullSigmoidFit();
	    //Detectors.at(i)->FFTWaveform();

	    //Detectors.at(i)->MicroSigmoid();
	    //Detectors.at(i)->Findmultipleinverse();
	      //  Detectors.at(i)->Findmultiplepeak();
	    //  std::cout<<"Mark 2 Sigmoid: "<<Detectors.at(i)->Sigmoid.timepoint<<" cfd: "<<Detectors.at(i)->cfd_time<<std::endl;
        }
	// Detectors.at(i)->TimeInflection();
        //Detectors.at(i)->TimeTwentyPercent();
    }
}


  std::fstream dumptxt;


void TestBeamSetup::Dump()
{

  dumptxt.open ("waveform_dump.txt");
    TFile dumpfile("testbeam_dumpfile.root","recreate");
    //TTree *newtree = OutTree->CloneTree(0);
    //newtree->Fill();

    TObject integer;

    integer.SetUniqueID(Detectors.size());
    integer.Write("n");

    TGraph gr;
    for(int i = 0; i < Detectors.size() ; ++i)
    {
      for(int j=0; j<Detectors.at(i)->waveform_x.size();j++) {
	dumptxt<<&Detectors.at(i)->waveform_x[j]<<"   "<<&Detectors.at(i)->waveform_y[j]<<std::endl;
	}
      
        gr = TGraph(Detectors.at(i)->waveform_x.size(), &Detectors.at(i)->waveform_x[0], &Detectors.at(i)->waveform_y[0]);
        char str1[20];
        sprintf(str1,"graph_%d",i);
        gr.Write(str1);
        if( Detectors.at(i)->pre_filter_backup ) 
        {
            gr = TGraph(Detectors.at(i)->pre_filter_backup->waveform_x.size(), &Detectors.at(i)->pre_filter_backup->waveform_x[0], &Detectors.at(i)->pre_filter_backup->waveform_y[0]);
            sprintf(str1,"unfiltered_graph_%d",i);
            gr.Write(str1);
        }
      
    }

    
    
    //Detectors.at(0)->pre_filter_backup->Sigmoid.fit_func.Write("rawsigmoid");
    //newtree->Write();
    dumpfile.Close();

}

//void TestBeamSetup::SetWaveformToAverage(AverageTool &aver)
//{
//    ScaleAndShiftTimes();
//    for(int i = 0; i < NofDetectors; ++i)
//    {
//        if( Detectors.at(i)->type > 0 )
//            Detectors.at(i)->InvertY();
//    }
//    
//    baseline_region_end=2000000;
//    for(int i = 0; i < NofDetectors; ++i)
//    {
//        Detectors.at(i)->FindGlobalMaximum(0, Detectors.at(i)->waveform_y.size()-1);
//        if(Detectors.at(i)->global_maximum.position < baseline_region_end)
//            baseline_region_end = Detectors.at(i)->global_maximum.position;
//    }
//    baseline_region_end-=200;
//    if(baseline_region_end < 2) baseline_region_end = 2;
//    //===
//    //
//    int j = 0;
//    double ref_time=0;
//    double n = 0;
//    for(int i = 0; i < NofDetectors; ++i)
//    {
//        Detectors.at(i)->SubstractBaseline(baseline_region_end);
//        if( Detectors.at(i)->type > 0 )
//        {
//            j = i;
//            Detectors.at(i)->FilterWaveformFFT(0, 4096, 1.0);
//            Detectors.at(i)->FindGlobalMaximum(0, 4095); 
//            Detectors.at(i)->FindStartPoint(baseline_region_end);
//            Detectors.at(i)->FindEndPoint(max_region_end-baseline_region_end);
//            Detectors.at(i)->FindElectronPeakEndPoint();
//            Detectors.at(i)->CalculateCharges();
//            Detectors.at(i)->FindNaiveTiming();
//        }
//        else
//        {
//            Detectors.at(i)->FindGlobalMaximum(0, Detectors.at(i)->waveform_y.size()-1); 
//            Detectors.at(i)->FindStartPoint(baseline_region_end);
//            Detectors.at(i)->FindEndPoint(max_region_end-baseline_region_end);
//            Detectors.at(i)->FindNaiveTiming();
//            Detectors.at(i)->TimeTwentyPercent();
//            Detectors.at(i)->FindFirstPeak();
//            Detectors.at(i)->FindMaxDerivative();
//            Detectors.at(i)->TimeInflection();
//            //ref_time += Detectors.at(i)->Inflection.timing;
//            ref_time += Detectors.at(i)->TwentyPercent.x;
//            n++;
//        }
//    }
//
//    ref_time/= n; 
//
//    double normalization = 1./Detectors.at(j)->global_maximum.y;
//    double normalization2 = 1./Detectors.at(j)->charge_e_peak;
//    //double normalization = 0; 
//    //for(int i = 0; i < Detectors.at(j)->waveform_y.size(); ++i)
//    //{
//    //    normalization += Detectors.at(j)->waveform_y.at(i);
//    //}
//    //normalization /= Detectors.at(j)->waveform_y.size();
//
//    //double ph = Detectors.at(j)->charge_e_peak;
//    double tt = Detectors.at(0)->naive_time - Detectors.at(1)->TwentyPercent.x; 
//    double ph = Detectors.at(j)->global_maximum.y;
//    double ch = Detectors.at(j)->charge_e_peak;
//    //if (ph <0.3|| ph>=0.42) baseline_region_end = -1;
//    if (ch < 1.9|| ch>= 2.2) baseline_region_end = -1;
//    //if (tt<7.5) baseline_region_end = -1;
//    aver.SetWaveform(Detectors.at(j)->waveform_x, Detectors.at(j)->waveform_y, ref_time, normalization, baseline_region_end, normalization2);
//    //aver.SetWaveform(Detectors.at(j)->waveform_x, Detectors.at(j)->waveform_y, ref_time, normalization, baseline_region_end);
//}

void TestBeamSetup::SetWaveformToAverage(AverageTool &aver)
{
    ScaleAndShiftTimes();
    for(int i = 0; i < NofDetectors; ++i)
    {
            if(i == 0) Detectors.at(i)->InvertY();
    }
    baseline_region_end=2000000;
    for(int i = 0; i < NofDetectors; ++i)
    {
        Detectors.at(i)->FindGlobalMaximum(0, Detectors.at(i)->waveform_y.size()-1);
        if(Detectors.at(i)->global_maximum.position < baseline_region_end)
            baseline_region_end = Detectors.at(i)->global_maximum.position;
    }
    baseline_region_end-=200;
    if(baseline_region_end < 2) baseline_region_end = 2;
    max_region_end = 2000+baseline_region_end;
    //===
    //
    int j = 0;
    double ref_time=0;
    double n = 0;
    for(int i = 0; i < NofDetectors; ++i)
    {
        Detectors.at(i)->SubstractBaseline(baseline_region_end);
        Detectors.at(i)->FindGlobalMaximum(0, Detectors.at(i)->waveform_y.size()-1); 
        if( Detectors.at(i)->type > 0 )
        {
            //j = i;
        }
        else
        {

            Detectors.at(i)->FindStartPoint(baseline_region_end);
            Detectors.at(i)->FindEndPoint(max_region_end);
            Detectors.at(i)->FindElectronPeakEndPoint();
            Detectors.at(i)->CalculateCharges();
            Detectors.at(i)->FindNaiveTiming();
            Detectors.at(i)->FindRiseTime();
            Detectors.at(i)->FindFirstPeak();
            Detectors.at(i)->FindMaxDerivative();

            Detectors.at(i)->TimeTwentyPercent();
            Detectors.at(i)->TimeInflection();
            //ref_time += Detectors.at(i)->Inflection.timing;
            ref_time += Detectors.at(i)->TwentyPercent.x;
            n++;
        }
    }

    ref_time/= n; 
    if(ref_time > 1.e+40)
        ref_time = 15.;

    double normalization = 1./Detectors.at(j)->global_maximum.y;
    //double normalization = 0; 
    //for(int i = 0; i < Detectors.at(j)->waveform_y.size(); ++i)
    //{
    //    normalization += Detectors.at(j)->waveform_y.at(i);
    //}
    //normalization /= Detectors.at(j)->waveform_y.size();

    aver.SetWaveform(Detectors.at(j)->waveform_x, Detectors.at(j)->waveform_y, ref_time, normalization, baseline_region_end);

}

void TestBeamSetup::init()
{
    NofDetectors=Detectors.size();
    //max_region_end = 20001;
    init_tree();
}

void TestBeamSetup::tree_add_tracks(int& NumberofTracks, Track& OneTrack)
{
  

    OutTree->Branch("NumberofTracks",&NumberofTracks,"NumberofTracks/I");
  OutTree->Branch("POS0_track_hit",&OneTrack.Hit1);
   OutTree->Branch("SRSnumber",&OneTrack.Hit2);
    //OutTree->Branch("POS2_track_hit",&OneTrack.Hit3);
    OutTree->Branch("TrackChi2",&OneTrack.Chi2,"TrackChi2/D");
    OutTree->Branch("SlopeXZ",&OneTrack.SlopeXZ,"SlopeXZ/D");
    OutTree->Branch("SlopeYZ",&OneTrack.SlopeYZ,"SlopeYZ/D");
//    OutTree->Branch("SRSnumber",&OneTrack.SRSnumber,"SRSNumber/I");
    
}

void TestBeamSetup::init(int& NumberofTracks, Track& OneTrack)
{
    init();
    tree_add_tracks(NumberofTracks, OneTrack);
}

void TestBeamSetup::init_tree()
{
    //OutFile = new TFile("outrootfile.root","recreate");
    OutTree = new TTree("Pico","Analysis Output");
    int mcp_no = 0;
    int mm_no = 0;
    int filter = 0;
    std::cout << "===List of different types of characteristics===" << std::endl;
    OutTree->Branch("Baseline_Window",&baseline_region_end,"Baseline_Window/I");
    //OutTree->Branch("EvntCtr",&EvntCtr,"EventCounter/I");
    for(int i = 0; i < Detectors.size(); ++i)
    {
        Detector* det;
        int type = Detectors.at(i)->type;
        std::string typestr;
        if(type > 0)
        {
            mm_no++;
            char str[20];
            sprintf(str,"MM%d_",mm_no);
            typestr = str;
        }
        else
        {
            mcp_no++;
            char str[20];
            sprintf(str,"MCP%d_",mcp_no);
            typestr = str;
        }
        bool is_filter = 0; 
        if(type > 2)
        {
            if(filter == mm_no - 1)
            {
                det = Detectors.at(i);
                //std::cout << "fmm" << i << std::endl;
                filter++;
                i--;
                mm_no--;
                is_filter = 1;
                typestr = "Filtered_" + typestr;
            }
            else
            {
                det = Detectors.at(i)->pre_filter_backup;
                //std::cout << "mm" << i << std::endl;
            }
       
        }
        else
        {
            det = Detectors.at(i);
            //std::cout << "mcp" << i << std::endl;
        }
        
        std::string varname;
        std::string leafname;

        varname = typestr + "baseline_level";
        leafname = varname + "/D";
        if(!is_filter)
            OutTree->Branch(varname.c_str(),&det->baseline_level,leafname.c_str());

        varname = typestr + "baseline_rms";
        leafname = varname + "/D";
        if(!is_filter)
            OutTree->Branch(varname.c_str(),&det->baseline_rms,leafname.c_str());

        varname = typestr + "global_maximum_y";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(),&det->global_maximum.y,leafname.c_str());

        varname = typestr + "global_maximum_x";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(),&det->global_maximum.x,leafname.c_str());

        varname = typestr + "start_x";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(),&det->start_point.x,leafname.c_str());

        
	varname = typestr + "e_peak_end_x";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(),&det->e_peak_end.x,leafname.c_str());

	/*
        varname = typestr + "inflection_timing";
        leafname = varname + "/D";
        if(is_filter || type < 0)
            OutTree->Branch(varname.c_str(),&det->Inflection.timing,leafname.c_str());

        varname = typestr + "inflection_failed";
        leafname = varname + "/O";
        if(is_filter || type < 0)
            OutTree->Branch(varname.c_str(),&det->Inflection.failed,leafname.c_str());

        varname = typestr + "half_charge";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(),&det->charge_leading_edge,leafname.c_str());
	*/

        varname = typestr + "e_charge";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(),&det->charge_e_peak,leafname.c_str());

        varname = typestr + "all_charge";
        leafname = varname + "/D";
        OutTree->Branch(varname.c_str(),&det->charge_all,leafname.c_str());

        // varname = typestr + "rise_time";
        // leafname = varname + "/D";
        // OutTree->Branch(varname.c_str(),&det->rise_time,leafname.c_str());

        // varname = typestr + "twentypercent_x";
        // leafname = varname + "/D";
        // if(is_filter || type < 0)
        //     OutTree->Branch(varname.c_str(),&det->TwentyPercent.x,leafname.c_str());

	//   varname = typestr + "twentypercent_time";
        // leafname = varname + "/D";
        // if(is_filter || type < 0)
	//    OutTree->Branch(varname.c_str(),&det->TwentyPercent.timing,leafname.c_str());
	

        // varname = typestr + "twentypercent_failed";
        // leafname = varname + "/O";
        // if(is_filter || type < 0)
        //     OutTree->Branch(varname.c_str(),&det->TwentyPercent.failed,leafname.c_str());

	 	 varname = typestr + "naive_time";
	  leafname = varname + "/D";
	 OutTree->Branch(varname.c_str(),&det->naive_time,leafname.c_str());
	
	 varname = typestr + "naive_x";
	 leafname = varname + "/D";
	 OutTree->Branch(varname.c_str(),&det->naive_point.x,leafname.c_str());

	/*
	varname = typestr + "cfd_time";
	leafname = varname + "/D";
	OutTree->Branch(varname.c_str(),&det->cfd_time,leafname.c_str());
	

			varname = typestr + "inv_time1";
	leafname = varname + "/D";
	OutTree->Branch(varname.c_str(),&det->inversetime1,leafname.c_str());

	varname = typestr + "inv_time2";
	leafname = varname + "/D";
	OutTree->Branch(varname.c_str(),&det->invtime2,leafname.c_str());
	*/
	/*		varname = typestr + "inv_amp1";
	leafname = varname + "/D";
	OutTree->Branch(varname.c_str(),&det->ampfirst,leafname.c_str());
	
	varname = typestr + "inv_amp2";
	leafname = varname + "/D";
	OutTree->Branch(varname.c_str(),&det->Inverse.amp2,leafname.c_str());
	*/
	// 	varname = typestr + "cfd_maxval";
	// leafname = varname + "/D";
	// OutTree->Branch(varname.c_str(),&det->globmaxval,leafname.c_str());

	// 	varname = typestr + "cfd_maxpos";
	// leafname = varname + "/D";
	// OutTree->Branch(varname.c_str(),&det->globmaxpos,leafname.c_str());

	// 	varname = typestr + "cfd_minval";
	// leafname = varname + "/D";
	// OutTree->Branch(varname.c_str(),&det->globminval,leafname.c_str());

	// 	varname = typestr + "cfd_minpos";
	// leafname = varname + "/D";
	// OutTree->Branch(varname.c_str(),&det->globminpos,leafname.c_str());

	// 	varname = typestr + "cfd_zeropos";
	// leafname = varname + "/D";
	// OutTree->Branch(varname.c_str(),&det->zeropos,leafname.c_str());


        if(type > 0)
        {

	  
            varname = typestr + "sigmoid_parameters";
            leafname = varname + "[4]/D";
            OutTree->Branch(varname.c_str(),det->Sigmoid.parameters,leafname.c_str());

            varname = typestr + "sigmoid_chi_square";
            leafname = varname + "/D";
            OutTree->Branch(varname.c_str(),&det->Sigmoid.chisquare,leafname.c_str());

            varname = typestr + "sigmoid_failed";
            leafname = varname + "/O";
            OutTree->Branch(varname.c_str(),&det->Sigmoid.failed,leafname.c_str());
	    
	    varname = typestr + "sigmoid_timepoint";
            leafname = varname + "/D";
            OutTree->Branch(varname.c_str(),&det->Sigmoid.timepoint,leafname.c_str());
	  

	    varname = typestr + "Fullsigmoid_charge";
            leafname = varname + "/D";
            OutTree->Branch(varname.c_str(),&det->FullSigmoid.charge,leafname.c_str());
	    
	   varname = typestr + "Fullsigmoid_chi2";
            leafname = varname + "/D";
            OutTree->Branch(varname.c_str(),&det->FullSigmoid.chisquare,leafname.c_str());

	    //    varname = typestr + "Fullsigmoid_fail";
            //leafname = varname + "/D";
	    // OutTree->Branch(varname.c_str(),&det->FullSigmoid.failure,leafname.c_str());
	    
	    varname = typestr + "Fullsigmoid_parameters";
            leafname = varname + "[7]/D";
            OutTree->Branch(varname.c_str(),det->FullSigmoid.parameters,leafname.c_str());

            if(i == 0){
                varname = typestr + "Expofit";
                leafname = varname + "/D";
                OutTree->Branch(varname.c_str(),&det->charge_expo_full,leafname.c_str());
            }
	    /*
	       varname = typestr + "Fullsigmoid_peakx";
            leafname = varname + "/D";
            OutTree->Branch(varname.c_str(),&det->FullSigmoid.PeakX,leafname.c_str());
	       varname = typestr + "Fullsigmoid_peaky";
            leafname = varname + "/D";
            OutTree->Branch(varname.c_str(),&det->FullSigmoid.PeakY,leafname.c_str());
	    
	    varname = typestr + "micro_timepoint";
            leafname = varname + "/D";
            OutTree->Branch(varname.c_str(),&det->Sigmoid.microbulk,leafname.c_str());

	    varname = typestr + "micro_peak";
            leafname = varname + "/D";
            OutTree->Branch(varname.c_str(),&det->Sigmoid.micropeak,leafname.c_str());

	    varname = typestr + "micro_peakpoint";
            leafname = varname + "/D";
            OutTree->Branch(varname.c_str(),&det->Sigmoid.micropeakpos,leafname.c_str());
	    
            varname = typestr + "first_peak_y";
            leafname = varname + "/D";
            OutTree->Branch(varname.c_str(),&det->first_peak.y,leafname.c_str());

            varname = typestr + "first_peak_x";
            leafname = varname + "/D";
            OutTree->Branch(varname.c_str(),&det->first_peak.x,leafname.c_str());
	    */

        }

        std::cout << typestr << std::endl;
    }
}

