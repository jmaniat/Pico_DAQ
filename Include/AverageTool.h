#ifndef __averagetool_h__
#define __averagetool_h__

#include <iostream>
#include <vector>

#include <TFile.h>
#include <TGraph.h>
#include <TMath.h>
#include <TVirtualFFT.h>

class AverageTool
{
    public:
        AverageTool(int n);
        ~AverageTool();
        
        void SetWaveform(std::vector<double> x, std::vector<double> y, double ref_time, double norm, int baseline_end);
        void SetWaveform(std::vector<double> x, std::vector<double> y, double ref_time, double norm, int baseline_end, double charge_norm);
        void StandardAverage();
        void GetNoiseSpectrum();
        void GetSpectrum();
        void GetSpectrumOfAverage();
        void AddSpectrumToAverage();
        void AddNoiseSpectrumToAverage();
        void AddWaveformToAverage();
        void Finalize();
        void Write(const char* outfile_name);

    private:
        std::vector<double> waveform_y;
        std::vector<double> waveform_x;

        double* spectrum_real;
        double* spectrum_imag;

        double* noise_spectrum_real;
        double* noise_spectrum_imag;

        double* average_of_waveform_x;
        double* average_of_waveform_y;

        double* average_of_spectrum_freq;
        double* average_of_spectrum_power;
        double* average_of_spectrum_real;
        double* average_of_spectrum_imag;

        double* average_of_spectrum_of_noise_freq;
        double* average_of_spectrum_of_noise_power;
        double* average_of_spectrum_of_noise_real;
        double* average_of_spectrum_of_noise_imag;

        double* spectrum_of_average_freq;
        double* spectrum_of_average_power;
        double* spectrum_of_average_real;
        double* spectrum_of_average_imag;

        int N;
        int N2;
        int m;
        int number_of_events;
        double Reference_time;
        double Normalization_factor;
        int baseline_region_end;
        double Pulse_height;
};

#endif
