#ifndef __detector_h__
#define __detector_h__

#include <vector>

#include <TMath.h>
#include <TVirtualFFT.h>
#include <TMatrixD.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TGraph.h>


#include <DetectorSetup.h>

struct WaveformPoint
{
    double y;
    double x;
    int position;
};

struct TimingInfo
{
    double y;
    double x;
    double slope;
    double intersect;
    double timing;
    double parameters[4];//pol3 parameters
    bool failed = 0;
};

struct InverseInfo
{
    double time1;
    double time2;
    double amp1;
    double amp2;
    
};

struct SigmoidInfo
{
    TF1 fit_func;
    double parameters[4];
    bool failed;
    double chisquare;
    int degrees_freedom;
    double timepoint;
    double microbulk;
    double micropeakpos;
    double micropeak;
    
};

struct FullSigmoidInfo
{
    TF1 fit_func;
    double parameters[7];
    bool failed;
    double chisquare;
    int degrees_freedom;
    double timepoint;
    double microbulk;
    double micropeakpos;
    double micropeak;
    double charge;
    double PeakX;
    double PeakY;
    int failure;
    
};

class Detector
{
    friend class DetectorSetup;
    friend class MuonSetup;
    friend class LaserSetup;
    friend class CalibrationSetup;
    friend class SimulationSetup;
    friend class TestBeamSetup;
    
public:
    Detector();
    ~Detector();
    void SetWaveY(std::vector<double> wave);
    void SetWaveX(std::vector<double> wave);
    void SetMCP();
    void SetMM();
    void InvertY();
    void SubstractBaseline(int base_region_end);
    void FindGlobalMaximum(int start, int end);
    void FindStartPoint(int start);
    void FindEndPoint(int start);
    void FindElectronPeakEndPoint();
    void FindElectronPeakEndPointSimulation();
    void CalculateCharges();
    void FindNaiveTiming();
    void FindCFDTiming();
    double linear_interpolation(double x1, double x2, double y1, double y2, double y);
    void FindFirstPeak();
    void FindMaxDerivative();
    void TimeInflection();
    bool FitPol3(double* x, double* y, double* fit_parameters);
    void TimeTwentyPercent();
    void TimeSigmoid();
    void FullSigmoidFit();
    void ExpoFit();
    
    //// ToT  ////
    void Threshold_Values();
    void ToT(int detector_type);
    std::map<std::string, std::vector<double>> GetToT();
    std::map<std::string, std::vector<double>> AddRef();
    void MicroSigmoid();
    void FindRiseTime();
    void Findmultiplepeak();
    void Findmultipleinverse();
    void FindWidth();
    //void FilterWaveformFFT_test();
    void FFTWaveform();
    void FilterWaveformFFT(int start, int N, double biggest_frequency_not_to_cut_GHz);
private:
    int type;
    
    std::vector<double> waveform_y;
    std::vector<double> waveform_x;
    
    std::vector<double> waveformFFT_y;
    std::vector<double> waveformFFT_x;
    
    Detector* pre_filter_backup;
    
    double baseline_level;
    double baseline_rms;
    
    WaveformPoint global_maximum;
    WaveformPoint start_point;
    WaveformPoint end_point;
    WaveformPoint e_peak_end;
    WaveformPoint naive_point;
    WaveformPoint first_peak;
    WaveformPoint max_derivative;
    
    WaveformPoint start_exp;
    WaveformPoint end_exp;
    
    double charge_leading_edge;
    double charge_e_peak;
    double charge_all;
    
    double ref_naive_time; //ns
    double naive_time;//ns
    double rise_time;// defined from 20% to 80% height of the pulse
    double width;
    double charge_expo_full;
    double cfd_time;
    double globmaxval;
    double globminval;
    int globmaxpos;
    int globminpos;
    int zeropos;
    
    double     inversetime1;
    double     invtime2;
    //    double somevariable;
    //    double     inverseamp2;
    
    TimingInfo Inflection;
    TimingInfo TwentyPercent;
    SigmoidInfo Sigmoid;
    FullSigmoidInfo FullSigmoid;
    
    InverseInfo Inverse;
    
    
    //// ToT ////
    std::vector<float> threshold;
    std::map<std::string, std::vector<double>> ToT_values;
    std::vector<double> ToT_fit;
    int numb_thres;
    int input_thres = 0;
};

double fermi_dirac(double *x, double *par);
double fermi_dirac_general(double *x, double *par);
double fermi_dirac_generalsub(double *x, double *par);
double fermi_dirac_simple(double x);

//  double pars[7];

#endif
