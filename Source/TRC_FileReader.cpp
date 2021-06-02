#include "TRC_FileReader.h"
#include "TrackInfo.h"


std::vector<int> trignumb(1,0);
int overflow=0;

TRC_FileReader::TRC_FileReader(std::vector<int> channel_ids, const char* dir_path, int max_n_files, const char* file_midname, int first_file_counter)
{

    First_File = first_file_counter;
    for(int i = 0; i < channel_ids.size(); ++i)
    {
        TRC_channel C(dir_path, max_n_files, channel_ids.at(i), file_midname, First_File);
        Channels.push_back(C);
    }
}

void TRC_FileReader::OpenTriggerChannel(const char* dir_path, int max_n_files, const char* file_midname, int ch_id)
{
    trigger_channel_is_open = 1;
    TriggerChannel = new TRC_channel(dir_path, max_n_files, ch_id, file_midname, First_File);
}


bool TRC_FileReader::GetNextEvent()
{
  //   std::cout << "--Debug1--" << std::endl;

    for(int i = 0; i < Channels.size(); ++i)
    {
      //  std::cout << "--Debug2--" << Channels.size() << std::endl;

      //  std::cout << "--Debug3--" << Channels.at(i).GetWaveformX().size()<< std::endl;

								  
        if(Channels.at(i).GetNextEvent())
        {
	  //                   std::cout << "--Debug3--" << std::endl;

          //  std::cout << "End of events" << std::endl;
            return 0;
        }
	//      std::cout << "--Debug3.2--" << std::endl;

								 

        Detectors->SetDetectorWaveform(i,Channels.at(i).GetWaveformX(),Channels.at(i).GetWaveformY());
	//std::cout<<i<<std::endl;

	
    }


         if( trigger_channel_is_open)
        {
	  //  std::cout << "End of trigger" << std::endl;
            TriggerChannel->GetNextEvent();
            Detectors->SetTriggerNumber( CalculateTriggerNumber() );

        }
	 //     std::cout << "--Debug4--" << std::endl;



    
    return 1;
}

void TRC_FileReader::SetDetectorSetup(DetectorSetup &det)
{
    Detectors = &det;
}

int TRC_FileReader::CalculateTriggerNumber()
{

  
    std::fstream myfile;
   double sample;
   int number, bit, kk, inipulse, endpulse, nsamples, step, lastnumb, numbertest;
    std::vector<double> trigger_wave_x = TriggerChannel->GetWaveformX();
    std::vector<double> trigger_wave_y = TriggerChannel->GetWaveformY();
    double fSampleTime = (trigger_wave_x.at(1) - trigger_wave_x.at(0))*1.e+9;
    number = 0; bit = 1;

    
    nsamples = trigger_wave_x.size();

    
    //   endpulse = nsamples/2;
    // std::cout<<nsamples<<" NSAMPLES"<<std::endl;
    //endpulse = 600;//20 ks endpulse=10250; 10 kS 5125: 40 ks 20500 desquewPool 550


    if(nsamples==40002){
      inipulse = (int)(29250) ;
      endpulse = (int) (12250);
      step = (int)(500);
    }
    else if(nsamples==20002){
      endpulse = (int)(10250) ;
      inipulse = (int) (18250);
      step = (int)(500);
    }
    else if(nsamples==10002){
      endpulse = (int)(5125) ;//9250
      inipulse = (int) (9125);//1250
      step = (int)(250);//500

      //endpulse = (int)(1250) ;
      //inipulse = (int) (9250);
      // step = (int)(500);
      
    }
    else if(nsamples==5002){
      endpulse = (int)(600) ;
      inipulse = (int) (4600);
      step = (int)(250);
    }
    else std::cout<<"ERROR in calculating Eventnumber"<<std::endl;

   //std::cout << "  " << nsamples<< "  "  << std::endl;

    lastnumb=trignumb.back();

    //std::cout << "  " << lastnumb<< "  "  << std::endl;
    
    //     inipulse = (int)(1480.0 / fSampleTime); //1480=> long bit;
    //inipulse = (int)(4600); //20 kS => 29600=> long bit; 18250 short bit; 10 kS => 9125: 40 ks 36500  desquewPool 4550
    //step = (int)(250); //20 kS 500; 10 kS 250: 40 ks 1000
    //     step=int)(25.0/fSampleTime);

    kk = inipulse;

    //   std::cout<<kk<< " > " << endpulse << " end "  << std::endl;

    while(kk>endpulse)
      {


      sample = trigger_wave_y.at(kk);
      // std::cout<<kk<< " kk " << sample << " sample "  << std::endl;

      if(sample < -0.4) // Last year 1.2, now is 0.5
        number += bit;
      
      bit = bit*2;
      kk = kk - step;

    }

    numbertest = number +  (overflow) * 65536;
    //std::cout << "  " << numbertest<< "  "  << std::endl;

       if(numbertest<lastnumb)overflow++;

       number +=  (overflow) * 65536;
      
      trignumb.push_back(number);
    
        std::cout << "SRS " << number<<" Overflow "<<overflow<<std::endl;
     //The SRS event number is returned.
      
       return ((number));
    
  //   myfile.open ("waveform.txt",std::fstream::app);
  //  myfile<<"TRACKNBR."<<std::endl;
  // myfile<<((number-1)/2)<<std::endl;
  //myfile<<"TRACK"<<std::endl;
    
  //  myfile.close();
    
       //      return ((number-1)/2);
}


///////////
//////////           OLD Trigger Signal 16 Bit
//////////
// int TRC_FileReader::CalculateTriggerNumber()
// {

//   std::cout << "SRS "  << std::endl;
//    double sample;
//     int number, bit, kk, inipulse, endpulse, nsamples, step;
//     std::vector<double> trigger_wave_x = TriggerChannel->GetWaveformX();
//     std::vector<double> trigger_wave_y = TriggerChannel->GetWaveformY();
//     double fSampleTime = (trigger_wave_x.at(1) - trigger_wave_x.at(0))*1.e+9;
//     number = 0; bit = 1;

//     nsamples = trigger_wave_x.size();
//     endpulse = nsamples / 2;

//     inipulse = (int)(740.0 / fSampleTime); //1480=> long bit;
//     step = (int)(25.0/fSampleTime);

//     kk = inipulse;
//     while(kk > endpulse)
//     {
//       //	   std::cout << "step1" << std::endl;

//       sample = trigger_wave_y.at(kk);
//       //  	   std::cout << "step2" << std::endl;
//       if(sample < -0.5) // Last year 1.2, now is 0.5
//         number += bit;

//       bit = bit*2;
//       kk = kk - step;

//     }
//     // std::cout << "SRS " << (number-1)/2 << std::endl;
//     // The SRS event number is returned.
//     return ((number-1)/2); 

// }
