#include "FileReader.h"

FileReader::FileReader(const char* a)
{
    infile = fopen(a,"r");
    if(!infile) 
        std::cout << "Invalid file" << std::endl;

    std::cout << "Opening Waveforms file: " << a << std::endl;
}

FileReader::~FileReader()
{
    fclose(infile);
}

bool FileReader::GetNextEvent(int n, int &TrigNum)
{
    double EventNumber, TriggerNumber, NumberofSamples;
    if( GetNextLine(EventNumber, TriggerNumber, NumberofSamples) )
        return 1;

    TrigNum = (int)TriggerNumber; //extract trigger number

    //read trigger times
    double time_zero[N_detectors];
    double min_time= 11111111111111111111.;
    if( GetNextLine(time_zero[0], time_zero[1], time_zero[2]) )
        return 1;
    for(int i = 0; i < N_detectors; ++i)
    {
        time_zero[i]*=1.e+9; // convert s to ns
        if(min_time > time_zero[i])
            min_time = time_zero[i];

	//	std::cout<<time_zero[i]<<std::endl;
    }

    // std::cout<<"BREAK"<<std::endl;

    //readjust offset
    for(int i = 0; i < N_detectors; ++i)
    {
        time_zero[i]-=min_time;
    }

    //read waveform
    std::vector<double> waveform_y[N_detectors];
    std::vector<double> waveform_x[N_detectors];
    for(int i = 0 ; i < (int)NumberofSamples; ++i)
    {
        double temp[N_detectors];
        if( GetNextLine(temp[0], temp[1], temp[2]) ) 
            return 1;
        for(int j = 0; j < N_detectors; ++j)
        {
            waveform_y[j].push_back(-temp[j]); //invert waveform
            double time = time_zero[j] + time_step*i;
            waveform_x[j].push_back(time);
        }
    }

    for(int i = 0; i < N_detectors; ++i)
    {
        Detectors->SetDetectorWaveform(i,waveform_x[i],waveform_y[i]);      
    }
    return 0;
}

void FileReader::SetDetectorSetup(DetectorSetup &det)
{
    Detectors = &det;
}

bool FileReader::GetNextLine(double &p1, double &p2, double &p3)
{
    char text[80];
    if(fgets(text, sizeof(text), infile))
    {
        //if(text[strlen(text)-1] == "\n")
        //    text[strlen(text)-1] = 0;
        p1 = atof(strtok(text," "));
        p2 = atof(strtok(NULL," "));
        p3 = atof(strtok(NULL," "));
        return 0;
    }
    else
    {
        //std::cout << "EOF" << std::endl;
        return 1;
    }
}
