#include "AverageTool.h"


AverageTool::AverageTool(int n)
{
    number_of_events = 0;
    N = n;
    N2 = n;
    //waveform_y = new double[N];
    //waveform_x = new double[N];

    spectrum_real = new double[N];
    spectrum_imag = new double[N];

    noise_spectrum_real = new double[N];
    noise_spectrum_imag = new double[N];

    average_of_waveform_x = new double[N];
    average_of_waveform_y = new double[N];

    average_of_spectrum_freq = new double[N];
    average_of_spectrum_power = new double[N];
    average_of_spectrum_real = new double[N];
    average_of_spectrum_imag = new double[N];

    average_of_spectrum_of_noise_freq = new double[N2];
    average_of_spectrum_of_noise_power = new double[N2];
    average_of_spectrum_of_noise_real = new double[N2];
    average_of_spectrum_of_noise_imag = new double[N2];

    spectrum_of_average_freq = new double[N];
    spectrum_of_average_power = new double[N];
    spectrum_of_average_real = new double[N];
    spectrum_of_average_imag = new double[N];
    for(int i = 0; i < N; ++i)
    {
        average_of_spectrum_freq[i] = 0;
        average_of_spectrum_real[i] = 0;
        average_of_spectrum_power[i] = 0;
        average_of_spectrum_imag[i] = 0;

        //average_of_waveform_x[i] = 0;
        average_of_waveform_y[i] = 0;
    }
    for(int i = 0; i < N2; ++i)
    {
        average_of_spectrum_of_noise_freq[i] = 0;
        average_of_spectrum_of_noise_real[i] = 0;
        average_of_spectrum_of_noise_imag[i] = 0;
        average_of_spectrum_of_noise_power[i] = 0;
    }
}

AverageTool::~AverageTool()
{
    //delete[] waveform_y;
    //delete[] waveform_x;

    delete[] spectrum_real;
    delete[] spectrum_imag;

    delete[] noise_spectrum_real;
    delete[] noise_spectrum_imag;

    delete[] average_of_waveform_x;
    delete[] average_of_waveform_y;

    delete[] average_of_spectrum_freq;
    delete[] average_of_spectrum_power;
    delete[] average_of_spectrum_real;
    delete[] average_of_spectrum_imag;

    delete[] average_of_spectrum_of_noise_freq;
    delete[] average_of_spectrum_of_noise_power;
    delete[] average_of_spectrum_of_noise_real;
    delete[] average_of_spectrum_of_noise_imag;

    delete[] spectrum_of_average_freq;
    delete[] spectrum_of_average_power;
    delete[] spectrum_of_average_real;
    delete[] spectrum_of_average_imag;
}


void AverageTool::SetWaveform(std::vector<double> x, std::vector<double> y, double ref_time, double norm, int baseline_end)
{
    waveform_x = x;
    waveform_y = y;
    m = x.size();
    for(int i = 0; i < m; ++i)
        waveform_y[i]*=norm;
    //for(int i = 0; i < N; ++i)
    //{
    //    waveform_x[i] = x.at((i+start)%x.size());
    //    waveform_y[i] = y.at((i+start)%y.size());
    //}
    Reference_time = ref_time;
    //Reference_time = 0;//cividec
    baseline_region_end = baseline_end;
    Normalization_factor = norm;
    Pulse_height = norm;
}

void AverageTool::SetWaveform(std::vector<double> x, std::vector<double> y, double ref_time, double norm, int baseline_end, double charge_norm)
{
    waveform_x = x;
    waveform_y = y;
    m = x.size();
    for(int i = 0; i < m; ++i)
        waveform_y[i]*=charge_norm;
    //for(int i = 0; i < N; ++i)
    //{
    //    waveform_x[i] = x.at((i+start)%x.size());
    //    waveform_y[i] = y.at((i+start)%y.size());
    //}
    Reference_time = ref_time;
    baseline_region_end = baseline_end;
    Pulse_height = norm;
    Normalization_factor = charge_norm;
}

void AverageTool::StandardAverage()
{
    double pulse_height = 1./Normalization_factor;
    if(Reference_time != Reference_time) return; 
    if(baseline_region_end < 0) return;
    //if(pulse_height < 0.03) return;
    //if(pulse_height > 0.4) return;
    GetNoiseSpectrum();
    AddNoiseSpectrumToAverage();
    GetSpectrum();
    AddSpectrumToAverage();
    AddWaveformToAverage();
    number_of_events++;
}

void AverageTool::GetNoiseSpectrum()
{
    double temp_noise[N2];
    for(int i = N2 - 1; i >= 0; i--)
    {
        int pos = i - N2 + 1 + baseline_region_end;
        if( pos < 0 )
            pos += m;
        temp_noise[i] = waveform_y.at(pos%m);
    }
    
    TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N2, "R2C M");
    fft_forward->SetPoints(temp_noise);
    fft_forward->Transform();

    fft_forward->GetPointsComplex(noise_spectrum_real, noise_spectrum_imag);
}

void AverageTool::GetSpectrum()
{
    double temp_wave[N];
    for(int i = 0; i < N; ++i)
        temp_wave[i]= waveform_y.at((i+baseline_region_end)%m);

    TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N, "R2C M");
    fft_forward->SetPoints(temp_wave);
    fft_forward->Transform();

    fft_forward->GetPointsComplex(spectrum_real, spectrum_imag);
}

void AverageTool::GetSpectrumOfAverage()
{
    double temp_wave[N];
    for(int i = 0; i < N; ++i)
        temp_wave[i]= average_of_waveform_y[i%m];

    TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N, "R2C M");
    fft_forward->SetPoints(temp_wave);
    fft_forward->Transform();

    fft_forward->GetPointsComplex(spectrum_of_average_real, spectrum_of_average_imag);
}

void AverageTool::AddSpectrumToAverage()
{
    for(int i = 0; i < N; ++i)
    {
        average_of_spectrum_power[i] += spectrum_real[i]*spectrum_real[i] + 
                                       spectrum_imag[i]*spectrum_imag[i] ;
        average_of_spectrum_real[i] += spectrum_real[i];
        average_of_spectrum_imag[i] += spectrum_imag[i];
    }
}

void AverageTool::AddNoiseSpectrumToAverage()
{
    for(int i = 0; i < N2; ++i)
    {
        average_of_spectrum_of_noise_power[i] += noise_spectrum_real[i]*noise_spectrum_real[i] + noise_spectrum_imag[i]*noise_spectrum_imag[i];
        average_of_spectrum_of_noise_real[i] += noise_spectrum_real[i];
        average_of_spectrum_of_noise_imag[i] += noise_spectrum_imag[i];
    }
}

void AverageTool::AddWaveformToAverage()
{
    double step = waveform_x[1]-waveform_x[0];
    int Reference_pos = (int)(Reference_time/step);
    //int Move_pos = 200 - Reference_pos;//NotTestbeam
    int Move_pos = 200 - Reference_pos;//Testbeam
    //Move_pos=0; //cividec
    double x_after_move = (step - (Reference_time - Reference_pos*step) )/step;

    // slewing correction
    //Move_pos-= int(0.02582*TMath::Sqrt(Pulse_height)/step);
    //std::cout <<0.02582*TMath::Sqrt(Normalization_factor) << std::endl; 
    //

    //careful
    for(int i = 0; i < N; ++i)
    {
        int pos = i - Move_pos;
        if( pos < 0 ) pos += m;
        average_of_waveform_y[i%N] +=  ( waveform_y.at((pos+1)%m)-waveform_y.at(pos%m) )*x_after_move + waveform_y.at(pos%m);
    }
}

void AverageTool::Finalize()
{
    for(int i = 0 ; i < N; ++i)
    {
        average_of_spectrum_power[i] /= number_of_events;
        average_of_spectrum_real[i] /= number_of_events;
        average_of_spectrum_imag[i] /= number_of_events;
        average_of_spectrum_freq[i] = i/(waveform_x[1]-waveform_x[0])/N;
        average_of_waveform_y[i] /= number_of_events;
        average_of_waveform_x[i] = (waveform_x[1]-waveform_x[0])*i;
    }
    for(int i = 0 ; i < N2; ++i)
    {
        average_of_spectrum_of_noise_power[i] /= number_of_events;
        average_of_spectrum_of_noise_real[i] /= number_of_events;
        average_of_spectrum_of_noise_imag[i] /= number_of_events;
        //average_of_spectrum_of_noise_power[i] *= 2.; //Multiply by two to correct for the different
        //average_of_spectrum_of_noise_real[i] *= 2.; // number of points in the transform
        average_of_spectrum_of_noise_freq[i] = i/(waveform_x[1]-waveform_x[0])/N2;
    }
    //for(int i = 0 ; i < N; ++i)
    //{
    //    average_of_spectrum_power[i] = average_of_spectrum_real[i]*average_of_spectrum_real[i]
    //                                 + average_of_spectrum_imag[i]*average_of_spectrum_imag[i];
    //    average_of_spectrum_of_noise_power[i] = average_of_spectrum_of_noise_real[i]*average_of_spectrum_of_noise_real[i] + average_of_spectrum_of_noise_imag[i]*average_of_spectrum_of_noise_imag[i];
    //}

    //disable for cividec
    //renormalize average waveform to unity pulse height
    double average_of_waveform_height = -10;
    for(int i = 0; i < N; ++i)
        if(average_of_waveform_y[i] > average_of_waveform_height)
            average_of_waveform_height = average_of_waveform_y[i];
    double norm_factor = 1./average_of_waveform_height;
    for(int i = 0; i < N; ++i)
        average_of_waveform_y[i]*=norm_factor;
    //==================================================

    GetSpectrumOfAverage();
    for(int i = 0 ; i < N; ++i)
        spectrum_of_average_power[i] = spectrum_of_average_real[i]*spectrum_of_average_real[i]+spectrum_of_average_imag[i]*spectrum_of_average_imag[i];
}

void AverageTool::Write(const char* outfile_name)
{
    TFile *f = new TFile(outfile_name,"recreate");
    TGraph *gr = new TGraph(N/2, average_of_spectrum_freq, average_of_spectrum_power); 
    gr->SetMarkerStyle(20);
    gr->SetLineColor(kRed);
    gr->SetLineWidth(3);
    gr->Write("average_spectrum");

    TGraph *gr1 = new TGraph(N, average_of_waveform_x, average_of_waveform_y); 
    gr1->SetMarkerStyle(20);
    gr1->Write("average_waveform");

    TGraph *gr2 = new TGraph(N/2, average_of_spectrum_freq, spectrum_of_average_power); 
    gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kBlack);
    gr2->SetLineWidth(3);
    gr2->Write("spectrum_of_average");

    TGraph *gr3 = new TGraph(N2/2, average_of_spectrum_of_noise_freq, average_of_spectrum_of_noise_power); 
    gr3->SetMarkerStyle(20);
    gr3->SetLineColor(kBlue);
    gr3->SetLineWidth(3);
    gr3->Write("average_spectrum_noise");

    f->Close();
    delete gr;
    delete gr1;
    delete gr2;
    delete gr3;
    delete f;
}

