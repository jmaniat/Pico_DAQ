#include "DetectorSetup.h"

DetectorSetup::~DetectorSetup()
{
    for(int i = 0 ; i < Detectors.size(); ++i)
        delete Detectors.at(i);
}

void DetectorSetup::CreateMCP()
{
    Detector* MCP = new Detector();
    MCP->SetMCP();
    Detectors.push_back(MCP);
}

void DetectorSetup::CreateMM()
{
    Detector* MM = new Detector();
    MM->SetMM();
    Detectors.push_back(MM);
}

void DetectorSetup::Fill_Tree()
{
    OutTree->Fill();
}

void DetectorSetup::Finalize_Tree(const char* outrootfile)
{
    OutFile = new TFile(outrootfile,"recreate");
    OutTree->Write();
    OutFile->Write();
    OutFile->Close();
    delete OutFile;
    delete OutTree;
}

void DetectorSetup::SetDetectorWaveform(int i, std::vector<double> x, std::vector<double> y)
{
    Detectors.at(i)->SetWaveX(x);

    Detectors.at(i)->SetWaveY(y);
}

void DetectorSetup::SetTriggerNumber(int TrigNum)
{
    TriggerNumber = TrigNum;
}


/*
void DetectorSetup::DumpOld()
{
    TFile dumpfile("dumpfile.root","recreate");
    int N = Detectors.at(0)->waveform_y.size();
    TObject integer;
    integer.SetUniqueID(N);
    integer.Write("N");

    integer.SetUniqueID(Detectors.size());
    integer.Write("n");

    TVectorD dump_x(N);
    TVectorD dump_y(N);
    for(int i = 0; i < Detectors.size() ; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            dump_x[j] = Detectors.at(i)->waveform_x.at(j);
            dump_y[j] = Detectors.at(i)->waveform_y.at(j);
        }
        char str1[20];
        sprintf(str1,"x_wave_%d",i);
        dump_x.Write(str1);
        char str2[20];
        sprintf(str2,"y_wave_%d",i);
        dump_y.Write(str2);
        if( Detectors.at(i)->pre_filter_backup ) 
        {
        for(int j = 0; j < N; ++j)
        {
            dump_x[j] = Detectors.at(i)->pre_filter_backup->waveform_x.at(j);
            dump_y[j] = Detectors.at(i)->pre_filter_backup->waveform_y.at(j);
        }
        char str1[20];
        sprintf(str1,"x_wave_%d",10);
        dump_x.Write(str1);
        char str2[20];
        sprintf(str2,"y_wave_%d",10);
        dump_y.Write(str2);


        }
    }
    Detectors.at(0)->Sigmoid.fit_func.Write("rawsigmoid");
    dumpfile.Close();
}
*/

void DetectorSetup::ScaleAndShiftTimes()
{
  int t=0;
    double earliest_time = 10000000000;
    for(int i = 0; i < NofDetectors; ++i)
       if( Detectors.at(i)->waveform_x.at(0) < earliest_time )  {
           earliest_time = Detectors.at(i)->waveform_x.at(0);
	   t=i;

	     }
  //  std::cout<<"time "<<earliest_time<<" at  "<<t<<std::endl;

    for(int i = 0; i < NofDetectors; ++i)
    {

     
      //std::cout<<"before: "<<Detectors.at(i)->waveform_x.at(0)<<"  "<<Detectors.at(i)->waveform_x.at(1)<<"  "<<Detectors.at(i)->waveform_x.at(2)<<std::endl;
	
     
      
        std::transform(Detectors.at(i)->waveform_x.begin(),
                       Detectors.at(i)->waveform_x.end(),
                       Detectors.at(i)->waveform_x.begin(),
                       bind2nd(std::plus<double>(), -earliest_time));




	
      
        std::transform(Detectors.at(i)->waveform_x.begin(),
                       Detectors.at(i)->waveform_x.end(),
                       Detectors.at(i)->waveform_x.begin(),
                       [](double x){ return x*1.e+9; }); //Convert s to ns

	// std::cout<<"after: "<<Detectors.at(i)->waveform_x.at(0)<<"  "<<Detectors.at(i)->waveform_x.at(1)<<"  "<<Detectors.at(i)->waveform_x.at(2)<<std::endl;

	
    }
}

void DetectorSetup::SetWaveformToAverage(AverageTool &aver)
{
    /*
    ScaleAndShiftTimes();
    for(int i = 0; i < NofDetectors; ++i)
    {
        if( Detectors.at(i)->type > 0 )
            Detectors.at(i)->InvertY();
    }
    */
    baseline_region_end=2000000;
    for(int i = 0; i < NofDetectors; ++i)
    {
        Detectors.at(i)->FindGlobalMaximum(0, Detectors.at(i)->waveform_y.size()-1);
        if(Detectors.at(i)->global_maximum.position < baseline_region_end)
            baseline_region_end = Detectors.at(i)->global_maximum.position;
    }
    baseline_region_end-=200;
    if(baseline_region_end < 2) baseline_region_end = 2;
    //===
    //
    int j = 0;
    double ref_time=0;
    double n = 0;
    for(int i = 0; i < NofDetectors; ++i)
    {
        if( Detectors.at(i)->type > 0 )
        {
            j = i;
        }
        else
        {
            Detectors.at(i)->SubstractBaseline(baseline_region_end);
            Detectors.at(i)->FindGlobalMaximum(0, 4000); 
            Detectors.at(i)->FindStartPoint(baseline_region_end);
            Detectors.at(i)->FindEndPoint(max_region_end-baseline_region_end);
            Detectors.at(i)->FindNaiveTiming();
            Detectors.at(i)->TimeTwentyPercent();
            Detectors.at(i)->FindFirstPeak();
            Detectors.at(i)->FindMaxDerivative();
            Detectors.at(i)->TimeInflection();
	    Detectors.at(i)->TimeSigmoid();
//        Detectors.at(i)->FullSigmoidFit();
//            std::cout << "detector: " << i << std::endl;
//            std::cin.ignore();
            //ref_time += Detectors.at(i)->Inflection.timing;
            ref_time += Detectors.at(i)->TwentyPercent.x;
            n++;
        }
    }

    ref_time/= n; 

    double normalization = 1./Detectors.at(j)->global_maximum.y;
    //double normalization = 0; 
    //for(int i = 0; i < Detectors.at(j)->waveform_y.size(); ++i)
    //{
    //    normalization += Detectors.at(j)->waveform_y.at(i);
    //}
    //normalization /= Detectors.at(j)->waveform_y.size();

    aver.SetWaveform(Detectors.at(j)->waveform_x, Detectors.at(j)->waveform_y, ref_time, normalization, baseline_region_end);
}
