#include "Detector.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include "TCanvas.h"
#include "TAxis.h"
#include "TApplication.h"
#include "TSystem.h"


double pars[7];
TApplication *myapp=new TApplication("myapp",0,0);
auto wave_can = new TCanvas("wave_can", "wave_can", 0 , 0, 1600, 1600);
//auto wave_can2 = new TCanvas("wave_can2", "wave_can2", 0 , 0, 1600, 1600);


Detector::Detector()
{
    pre_filter_backup = 0;
}

Detector::~Detector()
{
    if(pre_filter_backup)
        delete pre_filter_backup;
}

void Detector::SetWaveY(std::vector<double> wave)
{
    waveform_y = wave;
}

void Detector::SetWaveX(std::vector<double> wave)
{
    waveform_x = wave;
}

void Detector::SetMCP()
{
    type = -1;
}

void Detector::SetMM()
{
    type = 1;
    pre_filter_backup = new Detector;
}

void Detector::InvertY()
{
    //Inverts waveform_y
    std::transform(waveform_y.begin(), waveform_y.end(), waveform_y.begin(),[](double y){ return -y; });
}

void Detector::SubstractBaseline(int base_region_end)
{
    //Find baseline
    double baseline_sum=0;
    double baseline_square_sum=0;
    for(int i = 0; i < base_region_end; ++i)
    {
        double y = waveform_y.at(i);
        baseline_sum+=y;
        baseline_square_sum+=y*y;
    }
    baseline_level = baseline_sum/base_region_end;
    double baseline_rms_squared = baseline_square_sum/base_region_end - baseline_level*baseline_level ;
    baseline_rms = baseline_rms_squared > 0 ? TMath::Sqrt(baseline_rms_squared) : 0;

    //Subtract baseline
    std::transform(waveform_y.begin(), waveform_y.end(), waveform_y.begin(),
                      bind2nd(std::plus<double>(), -baseline_level));
}


void Detector::FindGlobalMaximum(int start, int end)
{
    global_maximum.y=-111111;
    global_maximum.position=5;

    for(int i = start ; i <= end && i < waveform_y.size() ; i++)
    { 
        if( waveform_y.at(i) > global_maximum.y )
        {
            global_maximum.y = waveform_y.at(i);
            global_maximum.position = i;
        }
    }
    global_maximum.x = waveform_x.at(global_maximum.position);
}

void Detector::FindStartPoint(int start)
{
    start_point.y = 0;
    start_point.position = 1;
    for(int i = global_maximum.position; i >= start; --i)
    {
        if(waveform_y.at(i) - baseline_rms < 0 )
        {
            start_point.y = waveform_y.at(i);
            start_point.position = i;
            start_point.x = waveform_x.at(i);
            break;
        }
    }
}

void Detector::FindEndPoint(int start)
{
    int N = waveform_y.size();
    end_point.y = 0;
    end_point.position = N-1;

    int j = type > 0 ? start : global_maximum.position;
    for(int i = j ; i < N; ++i)
    {
        if( waveform_y.at(i) - baseline_rms < 0 )
        {
            end_point.y = waveform_y.at(i);
            end_point.position = i;
            end_point.x = waveform_x.at(i);
            break;
        }
    }
}

void Detector::FindElectronPeakEndPoint()
{
    if(type < 0) // is MCP
    {
        e_peak_end = end_point;
    }
    else  //is MM
    {
        e_peak_end.y = 11111111111.;
        e_peak_end.position = global_maximum.position;

//        int j = global_maximum.position + 60;
        int j = global_maximum.position + 120;

        if( j > waveform_y.size()-1 ) j = waveform_y.size()-2;

        for(int i = global_maximum.position + 3; i < j && i+1 < waveform_y.size() ; ++i)
        {
            if( waveform_y.at(i) < e_peak_end.y )
            {
                e_peak_end.y = waveform_y.at(i);
                e_peak_end.position = i;
                if(waveform_y.at(i+1) < 0) 
                    break;
//                                if(waveform_y.at(i+15) > e_peak_end.y ) break;

            }
        }
        e_peak_end.x = waveform_x.at(e_peak_end.position);
    }
}

void Detector::FindElectronPeakEndPointSimulation()
{
    if(type < 0) // is MCP
    {
        e_peak_end = end_point;
    }
    else  //is MM
    {
        int n = 15;
        int j = global_maximum.position + 60;
        if( j >= waveform_y.size()-n ) j = waveform_y.size()-1-n;

        e_peak_end.position = j;
        e_peak_end.y = waveform_y.at(j);
        

        for(int i = global_maximum.position + 3; i < j; ++i)
        {
            double mean_x = 0, mean_y = 0;
            for(int k = 0; k < n; k++)
            {
                mean_y += waveform_y.at(i+k);
                mean_x += waveform_x.at(i+k);
            }
            mean_y /= n;
            mean_x /= n;

            double num = 0, den = 0;
            for(int k = 0; k < n; k++)
            {
                num += (waveform_x.at(i+k) - mean_x)*(waveform_y.at(i+k) - mean_y);
                den += (waveform_x.at(i+k) - mean_x)*(waveform_x.at(i+k) - mean_x);
            }

            double slope = num/den;
            if ( abs(slope) < 0.001 ) 
            {
                e_peak_end.position = i;
                e_peak_end.y = waveform_y.at(i);
                break;
            }
        }
        e_peak_end.x = waveform_x.at(e_peak_end.position);
    }
}

void Detector::CalculateCharges()
{
    double Ohms = 50;
    double step = waveform_x.at(1)- waveform_x.at(0);
    double conversion = step/Ohms*1000;

    charge_leading_edge=0;
    for(int i = start_point.position; i <= global_maximum.position; ++i)
        charge_leading_edge += waveform_y.at(i);
    charge_leading_edge *= conversion;

    charge_e_peak=0;
    for(int i = start_point.position; i <= e_peak_end.position; ++i)
        charge_e_peak += waveform_y.at(i);
    charge_e_peak *= conversion;

    charge_all=0;
    for(int i = start_point.position; i <= end_point.position; ++i)
        charge_all += waveform_y.at(i);
    charge_all *= conversion;
}

void Detector::FindNaiveTiming()
{
    double cf = 0.2*global_maximum.y;
    for(int i = global_maximum.position ; i > start_point.position && i>1; --i)
    {
        if( waveform_y.at(i) - cf > 0 && waveform_y.at(i-1) - cf < 0 )
        {
            double x1 = waveform_x.at(i-1);
            double x2 = waveform_x.at(i);
            double y1 = waveform_y.at(i-1);
            double y2 = waveform_y.at(i);
            naive_time = linear_interpolation(x1,x2,y1,y2,cf);
            naive_point.y = waveform_y.at(i);
            naive_point.position = i;
            naive_point.x = waveform_x.at(i);

	    // std::cout<<global_maximum.y<<"  "<<global_maximum.x<<"  "<<start_point.position<<"  "<<start_point.x<<"  "<<naive_point.x <<std::endl;
	  //  std::cout<<"naive time"<<naive_time<<std::endl;
            break;
        }
    }
}

void Detector::FindCFDTiming()
{
  std::vector<double> CFD(waveform_x.size(),0);
  //std::fstream cfddump;
  //cfddump.open ("cfddump.txt",std::fstream::app);
   //Superposition of the waveform
   for(int i=0; i<(waveform_x.size()-3);i++){
     CFD.at(i)=waveform_y.at(i)-0.8*waveform_y.at(i+3);
    //  cfddump<<waveform_y.at(i)<<"  "<<CFD.at(i)<<"  "<<i<<std::endl;
   }

   //Find global max
   globmaxval=-111111;
   globmaxpos=0;

   for(int i = 0 ; i < CFD.size() ; i++){ 
     if( CFD.at(i) > globmaxval){
       globmaxval = CFD.at(i);
       globmaxpos = i;
     }
   }
   //std::cout<<"globmaxval "<<globmaxval<<" globmaxpos "<<globmaxpos<<std::endl;
   
   //Find global min
   globminval=111111;
   globminpos=0;

   for(int i = 0 ; i < CFD.size() ; i++){ 
     if( CFD.at(i) < globminval){
       globminval = CFD.at(i);
       globminpos = i;
     }
   }
   //std::cout<<"globminval "<<globminval<<" globminpos "<<globminpos<<std::endl;

  
   //Find zero crossing
   zeropos=0;
   for( int i=globminpos; i<=globmaxpos; i++){
     //  std::cout<<i<<"  "<<CFD.at(i)<<"  "<<waveform_x.at(i)<<std::endl;
     if( CFD.at(i) < 0 && CFD.at(i+1)>0){
       zeropos = i;
     }
     //std::cout<<zeropos<<" Zeropos "<<std::endl;
   }

    double x1 = waveform_x.at(zeropos);
    double x2 = waveform_x.at(zeropos+1);
    double y1 = CFD.at(zeropos);
    double y2 = CFD.at(zeropos+1);
    cfd_time = linear_interpolation(x1,x2,y1,y2,0);
   
  //  std::cout<<"CFD time"<<cfd_time<<std::endl;
  return;
}

double Detector::linear_interpolation(double x1, double x2, double y1, double y2, double y)
{
    double b = (y2-y1)/(x2-x1);
    double a = y2 - b*x2;
    double x = (y-a)/b;
    return x;
}

void Detector::FindFirstPeak()
{
    first_peak = global_maximum;
    int j = global_maximum.position;
    while(j+3 > waveform_y.size()) j--;
    if(type > 0) //is MM
    {
        for(int i = naive_point.position; i <= j && i+3 < waveform_y.size(); ++i)
        {
            while( i < 3 ) i++;

            double dm2 = waveform_y.at(i-2) - waveform_y.at(i-3); 
            double dm1 = waveform_y.at(i-1) - waveform_y.at(i-2); 
            double dm = waveform_y.at(i) - waveform_y.at(i-1); 
            double dp = waveform_y.at(i+1) - waveform_y.at(i); 
            double dp1 = waveform_y.at(i+2) - waveform_y.at(i+1); 
            double dp2 = waveform_y.at(i+3) - waveform_y.at(i+2);
            if(dm2 > 0 && dm1 >= 0 && dm >= 0 && dp <= 0 && dp1 <= 0 && dp2 < 0 && waveform_y.at(i) >= 0.4*global_maximum.y)
            {
                first_peak.y = waveform_y.at(i);
                first_peak.position = i;
                first_peak.x = waveform_x.at(i);
                break;
            }
        }
    }
}

void Detector::FindMaxDerivative()
{
    max_derivative.y = -1111111111;
    max_derivative.position = first_peak.position;

    for(int i = naive_point.position - 2; i < first_peak.position - 2 && i+1 < waveform_y.size(); ++i)
    {
        while (i < 0 ) i++;

        double dy = waveform_y.at(i+1) - waveform_y.at(i); 
        if( dy > max_derivative.y )
        {
            max_derivative.y = dy;
            max_derivative.position = i;
        }
    }
    max_derivative.x = waveform_x.at(max_derivative.position);
}

void Detector::TimeInflection()
{
    bool use_filtered_waveform = 0;
    int k = max_derivative.position;
    while(k-2 < 0) k++;
    while( k+2 >= waveform_y.size() ) k--;
    
    double y[4],x[4];
    double x_offset = waveform_x.at(k-2); //scale for numerical purposes
    for(int i = 0 ; i < 4; ++i)
    {
        x[i] = waveform_x.at(k-1+i) - x_offset;
        y[i] = waveform_y.at(k-1+i);
    }
    double par[4];
    bool success = FitPol3(x,y,par);
    if(!success)
    {
        Inflection.failed = 1;
        return;
    }

    Inflection.y = par[0] - par[1]*par[2]/3./par[3] + 2*par[2]*par[2]*par[2]/27./par[3]/par[3];
    Inflection.slope = par[1] - par[2]*par[2]/3./par[3];
    Inflection.x = -par[2]/3./par[3] + x_offset;
    Inflection.intersect = Inflection.y - Inflection.slope*Inflection.x;
    Inflection.timing = -Inflection.intersect/Inflection.slope;
    for(int i = 0; i < 4; ++i) Inflection.parameters[i] = par[i];
}

bool Detector::FitPol3(double* x, double* y, double* fit_parameters)
{
    TMatrixD X(4,4);
    for(int i = 0 ; i < 4; ++i)
    {
        X[i][0] = 1.;        
        X[i][1] = x[i];        
        X[i][2] = x[i]*x[i];        
        X[i][3] = x[i]*x[i]*x[i];        
    }
    TVectorD Y(4,y);
    double det;
    TVectorD pars = X.Invert(&det)*Y;

    if(det == 0) return 0;

    for(int i = 0; i < 4; ++i)
        fit_parameters[i] = pars[i];
    
    return 1;
}

void Detector::TimeTwentyPercent()
{
    //if(type > 0) return;
    bool use_filtered_waveform = 0;
    int k = naive_point.position;
    while(k-3 < 0) k++;
    while( k+3 >= waveform_y.size() ) k--;
    
    double y[4],x[4];
    double x_offset = waveform_x.at(k-2); //scale for numerical purposes
    for(int i = 0 ; i < 4; ++i)
    {
        x[i] = waveform_x.at(k-1+i) - x_offset;

	y[i] = waveform_y.at(k-1+i);
    }
    double par[4];
    bool success = FitPol3(x,y,par);
    if(!success)
    {
        TwentyPercent.failed = 1;
        return;
    }
    
    double par_inverse[4];
    success = FitPol3(y,x,par_inverse);
    if(!success)
    {
        TwentyPercent.failed = 1;
        return;
    }

    double x20 = par_inverse[0];
    double y20 = 0.2*global_maximum.y;
    double y20_pow = 1;
    for(int i = 1 ; i < 4; ++i)
    {
        y20_pow*=y20; 
        x20+=par_inverse[i]*y20_pow;
    }
    double y20e = par[0];
    double x20_pow = 1;
    for(int i = 1; i < 4; ++i)
    {
        x20_pow *= x20;
        y20e += par[i]*x20_pow; 
    }

    TwentyPercent.x = x20 + x_offset;
    TwentyPercent.y = y20e;
    TwentyPercent.slope = par[1] + 2*par[2]*x20 + 3*par[3]*x20*x20;
    TwentyPercent.intersect = TwentyPercent.y - TwentyPercent.slope*TwentyPercent.x;
    TwentyPercent.timing = -TwentyPercent.intersect/TwentyPercent.slope;
    for(int i = 0; i < 4; ++i) TwentyPercent.parameters[i] = par[i];

}

void Detector::Threshold_Values()
{   float thresholds_array[6] = {0.2, 0.4, 0.6, 0.8, 1., 1.2};
    numb_thres = sizeof(thresholds_array)/sizeof(*thresholds_array);
    
    for(int i = 0; i < numb_thres; i++){
        threshold.push_back(thresholds_array[i]);
    }
}


void Detector::ExpoFit()
{
    start_exp.position = global_maximum.position + 20;
    start_exp.y = waveform_y.at(start_exp.position);
    start_exp.x = waveform_x.at(start_exp.position);
    
    end_exp.position = start_exp.position + 30;
    end_exp.y = waveform_y.at(end_exp.position);
    end_exp.x = waveform_x.at(end_exp.position);
    
    TF1 exp_fit("exp_fit", "expo", start_exp.x, end_exp.x);
    int Npoints = end_exp.position - start_exp.position;
    double x[Npoints], y[Npoints], x_err[Npoints], y_err[Npoints];
    int counter = 0;
    for(int i = start_exp.position; i < end_exp.position; i++){
        x[counter] = waveform_x.at(i);
        y[counter] = waveform_y.at(i);
        counter++;
    }
    TGraph  *waveform = new TGraph(waveform_x.size(), &(waveform_x[0]), &(waveform_y[0]));
    waveform->GetXaxis()->SetRangeUser(waveform_x.at(start_exp.position - 50), waveform_x.at(end_exp.position + 150));
    waveform->SetMarkerStyle(20);
    waveform->SetMarkerSize(1.2);
    double pars[2];

    wave_can->SetLogy();
//    wave_can->cd();
//    waveform->Draw("AP");
    waveform->Fit("exp_fit","qR");
    pars[0] = exp_fit.GetParameter(0);
    pars[1] = exp_fit.GetParameter(1);
//    wave_can->Update();
//    wave_can->SaveAs("laser.pdf");

    TF1 inter_exp("inter_exp", "expo",waveform_x.at(start_exp.position - 50), waveform_x.at(end_exp.position + 350) );
    TF1 up_line("upline", "[0]",waveform_x.at(start_point.position), waveform_x.at(end_exp.position + 350));
    up_line.SetParameter(0,start_exp.y);
    TF1 down_line("downline", "[0]", waveform_x.at(start_point.position),waveform_x.at(end_exp.position+350));
    down_line.SetParameter(0,end_exp.y);
    inter_exp.SetParameters(pars[0],pars[1]);
//    wave_can2->cd();
//    waveform->Draw("AP");
//    inter_exp.Draw("same");
//    up_line.Draw("same");
//    down_line.Draw("same");
//    wave_can2->Update();
    
//    wave_can2->SaveAs("laser_fit.pdf");
//    std::cout << "\n par0 = " << pars[0]  << ", par1 = " << pars[1] << std::endl;
    
    
//    gSystem->ProcessEvents() ;
//    std::cin.ignore();

    float Ohms = 50.;
    float step = waveform_x.at(1) -waveform_x.at(0);
    float conversion = step/Ohms*1000;
    
    float charge_before_expo = 0;
    for(int i = start_point.position; i <= start_exp.position; i++){
        charge_before_expo = charge_before_expo + waveform_y.at(i);
    }
    charge_before_expo = charge_before_expo*conversion;
    
    float integration_limit  = (std::log(0.001) - pars[0])/pars[1];
    TF1 exp_integral("exp_integral", "expo",waveform_x.at(start_exp.position), integration_limit);
    auto charge_expo = (inter_exp.Integral(start_exp.x, integration_limit));
    charge_expo = (charge_expo/Ohms)*1000;
    
    charge_expo_full = charge_before_expo + charge_expo;
    
//    std::cout << "Fit ends on " << integration_limit << ", integral result = " << i << std::endl;
//    myapp->SetReturnFromRun(true);
}


void Detector::TimeSigmoid()
{

  // std::cout<<"SIGMOID"<<std::endl;
    int start = start_point.position - 20; 
    while(start < 0) start++;
    int end = global_maximum.position + 20;
    while(end >= waveform_y.size()) end--;
    int Npoints = end-start+1;
    if( Npoints > 100 || Npoints <=0 )
    {
        Sigmoid.failed = 1;
        //return;
    }
    

    double x[Npoints], y[Npoints], erx[Npoints], ery[Npoints];
    for(int i = start; i <= end; ++i)
    {
        x[i-start] = waveform_x.at(i);
        y[i-start] = waveform_y.at(i);
        erx[i-start] = 0;
        ery[i-start] = baseline_rms;
        
	//    std::cout<<waveform_x.at(i)<<" "<<waveform_y.at(i)<<"  "<<baseline_rms<<std::endl;
    }

    TGraphErrors waveform_graph(Npoints, x, y, erx, ery);

    double pars[4];
    pars[0]=global_maximum.y;
    pars[1]=(x[end-start] + x[1])/2;
    pars[2]=5./(x[end-start] - x[1]);
    pars[3]=0.;
    
    TF1 fd_fit("fd_fit", fermi_dirac, x[0],x[global_maximum.position-start], 4);
    fd_fit.SetParameters(pars[0],pars[1],pars[2],pars[3]);
    waveform_graph.Fit("fd_fit","qR");
    waveform_graph.Fit("fd_fit","qR");
    waveform_graph.Fit("fd_fit","qR");
    waveform_graph.Fit("fd_fit","qR");
    
    // TCanvas *c1 = new TCanvas("c1","c1",1600,1000);
    // waveform_graph.Draw();
    // c1->Update();
    // c1 -> SaveAs("sigmoid.pdf");
      
    
    Sigmoid.chisquare = fd_fit.GetChisquare();
    Sigmoid.parameters[0] = fd_fit.GetParameter(0);
    Sigmoid.parameters[1] = fd_fit.GetParameter(1);
    Sigmoid.parameters[2] = fd_fit.GetParameter(2);
    Sigmoid.parameters[3] = fd_fit.GetParameter(3);
    Sigmoid.degrees_freedom = fd_fit.GetNDF();
    Sigmoid.chisquare = fd_fit.GetChisquare()/fd_fit.GetNDF();
    Sigmoid.fit_func = fd_fit;
    Sigmoid.failed = Sigmoid.chisquare < 1000 ? 0 : 1;

    Sigmoid.timepoint=Sigmoid.parameters[1]-(TMath::Log((Sigmoid.parameters[0]/(0.2*global_maximum.y-Sigmoid.parameters[3]))-1)/Sigmoid.parameters[2]);
    
    if(input_thres == 0) Threshold_Values();
    input_thres++;
    ToT_fit.clear();
    for(int i = 0; i < threshold.size(); i++){
        ToT_fit.push_back(Sigmoid.parameters[1]-(TMath::Log((Sigmoid.parameters[0]/(threshold.at(i)*global_maximum.y-Sigmoid.parameters[3]))-1)/Sigmoid.parameters[2]));
//        std::cout << threshold.at(i) << ", " << Sigmoid.parameters[1]-(TMath::Log((Sigmoid.parameters[0]/(threshold.at(i)*global_maximum.y-Sigmoid.parameters[3]))-1)/Sigmoid.parameters[2]) << std::endl;
    }
    
}

void Detector::ToT(){
    float ampl_val;
    for(int i = 0; i < threshold.size(); i++){
        ToT_values[Form("first_%i",i)].push_back(ToT_fit.at(i)== ToT_fit.at(i) ? ToT_fit.at(i) : -1);
        for(int j = global_maximum.position; j < waveform_y.size(); j++){
            if(waveform_y.at(j) > threshold.at(i)) continue;
            ToT_values[Form("second_%i",i)].push_back(ToT_values[Form("first_%i",i)].back() != -1 ? waveform_x.at(j) : -1);
//            std::cout << threshold.at(i) << ", " << waveform_x.at(j) << std::endl;
            break;
        }
//        std::cout << i << ",  " << threshold.at(i) <<  ", first = " << ToT_values[Form("first_%i",i)].back() << ", second = " << ToT_values[Form("second_%i",i)].back() << std::endl;
    }
}

    std::map<std::string, std::vector<float>> Detector::GetToT() {

    std::map<std::string, std::vector<float>> output_map;
    output_map = ToT_values;
    std::cout << "inside GET " << ToT_values["first_0"].size() << std::endl;
    return output_map;
}

void Detector::FullSigmoidFit()
{
  
  int start = start_point.position - 125;
  int sigstart= start_point.position-start;
  int sigpeak=global_maximum.position-start;
  int sigpeakVal=global_maximum.y;
  
  while(start < 0) start++;

  int end = global_maximum.position + 20* ( global_maximum.position - start_point.position);

  if(end<= start)end= start +200;

  if(end>= waveform_y.size()){
    end= waveform_y.size()-1;
    // return;
  }
  //while(end >= waveform_y.size()) end--;

  int Npoints = end-start+1;

  // if( Npoints > 500 || Npoints <=0 )
  // {
      //  FullSigmoid.failed = 1;
      //     return;
  // }

  int sigend=0;
  
  for(int i=global_maximum.position; i<end; i++){
    if( waveform_y.at(i)<0.4*global_maximum.y){
      sigend=i;
      break;
    }
  }

  sigend-=start;
  
  //std::cout<<"ONE"<<std::endl;
  //std::cout<<start<<" "<<end<<std::endl;
 
    double x[Npoints], y[Npoints], erx[Npoints], ery[Npoints];
    for(int i = start; i <= end; ++i)
    {
        x[i-start] = waveform_x.at(i);
        y[i-start] = waveform_y.at(i);
        erx[i-start] = 0;
        ery[i-start] = baseline_rms;
	
       // std::cout<<waveform_x.at(i)<<" "<<waveform_y.at(i)<<"  "<<baseline_rms<<std::endl;
    }

     TGraphErrors waveform_graph(Npoints, x, y, erx, ery);
    

      pars[0]=sigpeakVal;
      pars[1]=x[sigstart+((sigstart-sigpeak)/2)];
      pars[2]=1;
      pars[3]=1.;
    
   
      TF1 fd_fit("fd_fit", fermi_dirac_general, x[0], x[sigpeak], 4);
    fd_fit.SetParameters(pars[0],pars[1],pars[2],pars[3]);

    waveform_graph.Fit("fd_fit","qR");
      waveform_graph.Fit("fd_fit","qR");
     waveform_graph.Fit("fd_fit","qR");
     waveform_graph.Fit("fd_fit","qR");
    // std::cout<<"TWO"<<std::endl;

   pars[0]=fd_fit.GetParameter(0);
   pars[1]=fd_fit.GetParameter(1);
   pars[2]=fd_fit.GetParameter(2);
   pars[3]=fd_fit.GetParameter(3);

    TF1 fd_fitbis("fd_fit", fermi_dirac_general, x[0], x[sigpeak+100], 4);
   fd_fitbis.SetParameters(pars[0],pars[1],pars[2],pars[3]);
  
   // TGraphErrors waveform_graph2(Npoints, x, y, erx, ery);

  
   pars[4]=fd_fit.GetParameter(1)+1;
   pars[5]=fd_fit.GetParameter(2);
   pars[6]=fd_fit.GetParameter(3);

   TF1 fd_fit2("fd_fit2", fermi_dirac_generalsub, x[0], x[sigend],7);
       fd_fit2.SetParameters(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6]);
    waveform_graph.Fit("fd_fit2","qR");
    waveform_graph.Fit("fd_fit2","qR");
     waveform_graph.Fit("fd_fit2","qR");
     waveform_graph.Fit("fd_fit2","qR");


   pars[0]=fd_fit2.GetParameter(0);
   pars[1]=fd_fit2.GetParameter(1);
   pars[2]=fd_fit2.GetParameter(2);
   pars[3]=fd_fit2.GetParameter(3);
   pars[4]=fd_fit2.GetParameter(4);
   pars[5]=fd_fit2.GetParameter(5);
   pars[6]=fd_fit2.GetParameter(6);   
   
   
   
   // /*
    TF1 fd_fit3("fd_fit3", fermi_dirac_generalsub, x[0], x[end-start],7);
   fd_fit3.SetParameters(pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6]);

   TF1 fd_fit3L("fd_fit", fermi_dirac_general, x[0], x[end-start], 4);
    fd_fit3L.SetParameters(pars[0],pars[1],pars[2],pars[3]);
    fd_fit3L.SetLineStyle(2);
    fd_fit3L.SetLineColor(kBlue);
    
      TF1 fd_fit3R("fd_fit", fermi_dirac_general, x[0], x[end-start], 4);
    fd_fit3R.SetParameters(pars[0],pars[4],pars[5],pars[6]);
      fd_fit3R.SetLineStyle(2);
    fd_fit3R.SetLineColor(kGreen);
  // */
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

   // /*
//     TCanvas *c1 = new TCanvas("c1","c1",1600,1000);
////    c1->Divide(1,2);
////    c1->cd(1);
//    waveform_graph.Draw();
//     fd_fitbis.Draw("same");
//    c1->Update();

//    c1->cd(2);
   
//   waveform_graph.Draw();
//
//    fd_fit3.Draw("same");
//       fd_fit3L.Draw("same");
//    fd_fit3R.Draw("same");

//   c1->Update();
   
    // */
   
   //
//    c1 -> SaveAs("Fullsigmoid2.pdf");

   
   

   // std::cout<<"FULL"<<std::endl;

   //std::cout<<  start<<"   "<<end<<std::endl;
   // std::cout<<  sigpeak<<"   "<<sigpeakVal<<std::endl;

   //  std::cout<<  <<"   "<<sigpeakVal<<std::endl;


   //std::cout<<  fd_fit2.GetChisquare()<<"   "<<std::endl;
   //std::cout<<  fd_fit3.GetParameter(0)<<"   "<<std::endl;
   // std::cout<<  fd_fit3.GetParameter(1)<<"   "<< std::endl;
   // std::cout<< fd_fit3.GetParameter(2)<<"   "<<std::endl;
   // std::cout<< fd_fit3.GetParameter(3)<<"   "<<std::endl;
   // std::cout<< fd_fit3.GetParameter(4)<<"   "<<std::endl;
    // std::cout<< fd_fit3.GetParameter(5)<<"   "<<std::endl;
    // std::cout<< fd_fit3.GetParameter(6)<<"   "<<std::endl;
    // std::cout<< fd_fit2.GetNDF()<<"   "<< std::endl;
   // std::cout<<  fd_fit2.GetChisquare()/fd_fit2.GetNDF()<<"   "<< std::endl;
    // std::cout<<  fd_fit3<<std::endl;
    
    
     //INTEGRATION

     double IntegrateSigmoid=0;
     double binNo;
     double Peakx;
     double Peaky;
     int intstart=sigstart+start;

     while(intstart>start){
       if(fermi_dirac_simple(waveform_x.at(intstart))<0||fermi_dirac_simple(waveform_x.at(intstart))==0.)break;
       intstart--;
     }

     //   std::cout<<start<<"  "<<sigstart+start<<"  "<<intstart<<"   "<<waveform_x.at(intstart)<<std::endl;
     

    double Ohms = 50;
    double step = waveform_x.at(intstart+1)- waveform_x.at(intstart);
    double conversion = step/Ohms*1000;
    

     // for(int i=x[0]-100;  i<x[end-start]+100; i++){
     for(int i=0;  i<(Npoints+start-intstart); i++){
       
       //   binNo=(double) i;
       binNo=waveform_x.at(intstart+i);
       
       
       if(fermi_dirac_simple(binNo)>0){
	 IntegrateSigmoid+=fermi_dirac_simple(binNo)*conversion;
	 //     std::cout<<binNo<<"  "<<fermi_dirac_simple(binNo)<<"  "<<fermi_dirac_simple(binNo)*conversion<<"   "<<IntegrateSigmoid<<std::endl;

	 
       }

     }
     //std::cout<<IntegrateSigmoid<<"   "<<PeakY<<"  "<<PeakX<< std::endl;
     //        std::cout<<IntegrateSigmoid<<"   "<<global_maximum.y<<"  "<<global_maximum.x<< std::endl;



     FullSigmoid.charge=IntegrateSigmoid;

     FullSigmoid.failure=0;


     FullSigmoid.parameters[0] = fd_fit2.GetParameter(0);
     FullSigmoid.parameters[1] = fd_fit2.GetParameter(1);
     FullSigmoid.parameters[2] = fd_fit2.GetParameter(2);
     FullSigmoid.parameters[3] = fd_fit2.GetParameter(3);
     FullSigmoid.parameters[4] = fd_fit2.GetParameter(4);
     FullSigmoid.parameters[5] = fd_fit2.GetParameter(5);
     FullSigmoid.parameters[6] = fd_fit2.GetParameter(6);

     /*
     if(fd_fit2.GetParameter(0)/global_maximum.y>1000 || (fd_fit2.GetParameter(0)||fd_fit2.GetParameter(1)||fd_fit2.GetParameter(2)||fd_fit2.GetParameter(3)||fd_fit2.GetParameter(4)||fd_fit2.GetParameter(5)||fd_fit2.GetParameter(6))<0){
       FullSigmoid.failure=1;
       FullSigmoid.chisquare =0;
       
        FullSigmoid.charge=0;
       
     }
     */



     
     FullSigmoid.chisquare = fd_fit2.GetChisquare()/fd_fit2.GetNDF();


     // FullSigmoid.PeakX=Peakx;
     //FullSigmoid.PeakY=Peaky;
     
     //   std::cout<<FullSigmoid.failure<<" "<< FullSigmoid.charge<<"  <-- "<< std::endl;
     
    
     for(int i=0; i<7;i++){
       pars[i]=0;
     }
     
     
     
     /*   FullSigmoid.chisquare = fd_fit2.GetChisquare();
     FullSigmoid.parameters[0] = fd_fit2.GetParameter(0);
     FullSigmoid.parameters[1] = fd_fit2.GetParameter(1);
     FullSigmoid.parameters[2] = fd_fit2.GetParameter(2);
     FullSigmoid.parameters[3] = fd_fit2.GetParameter(3);
     FullSigmoid.parameters[3] = fd_fit2.GetParameter(4);
     FullSigmoid.parameters[3] = fd_fit2.GetParameter(5);
     FullSigmoid.parameters[3] = fd_fit2.GetParameter(6);
     
     
     FullSigmoid.degrees_freedom = fd_fit2.GetNDF();
     */     
     // FullSigmoid.fit_func = fd_fit2;
     //FullSigmoid.failed = FullSigmoid.chisquare < 1000 ? 0 : 1;

     /*

     if(FullSigmoid.charge>0&&global_maximum.y>0.){
 std::cout<<"FULL"<<std::endl;

 TCanvas *c1 = new TCanvas("c1","c1",1600,1000);
    c1->Divide(1,2);
    c1->cd(1);
    waveform_graph.Draw();
     fd_fitbis.Draw("same");
    c1->Update();

    c1->cd(2);
   
   waveform_graph2.Draw();
   
    fd_fit3.Draw("same");
       fd_fit3L.Draw("same");
    fd_fit3R.Draw("same");

   c1->Update();
   
   c1 -> SaveAs("Fullsigmoid2.pdf");

   std::cout<< "--------------"<<std::endl;
   std::cout<<global_maximum.y  <<"  "<<fd_fit3.GetParameter(0)<<"  "<<fd_fit3.GetParameter(0)/global_maximum.y<<std::endl;
     std::cout<< FullSigmoid.charge<<"   "<< std::endl;
   std::cout<< "--------------"<<std::endl;

    std::cout<< fd_fit3.GetParameter(4)<<"   "<<std::endl;
     std::cout<< fd_fit3.GetParameter(5)<<"   "<<std::endl;
     std::cout<< fd_fit3.GetParameter(6)<<"   "<<std::endl;
   std::cout<< "--------------"<<std::endl;

   std::cout<<  start<<"   "<<end<<std::endl;
   std::cout<<  sigpeak<<"   "<<sigpeakVal<<std::endl;

      std::cout<<global_maximum.y  <<std::endl;
      
     std::cout<< FullSigmoid.charge<<"   "<< std::endl;

     std::cout<<  fd_fit2.GetChisquare()/fd_fit2.GetNDF()<<"   "<< std::endl;

     std::cout<<  fd_fit3.GetParameter(0)<<"   "<<std::endl;
   std::cout<<  fd_fit3.GetParameter(1)<<"   "<< std::endl;
    std::cout<< fd_fit3.GetParameter(2)<<"   "<<std::endl;
    std::cout<< fd_fit3.GetParameter(3)<<"   "<<std::endl;
    std::cout<< fd_fit3.GetParameter(4)<<"   "<<std::endl;
     std::cout<< fd_fit3.GetParameter(5)<<"   "<<std::endl;
     std::cout<< fd_fit3.GetParameter(6)<<"   "<<std::endl;
     }

     */
     
}


double fermi_dirac(double *x, double *par)
{
    double fdreturn = par[0]/(1+TMath::Exp(-(x[0]-par[1])*par[2]))+par[3];
    return fdreturn;
}

double fermi_dirac_general(double *x, double *par)
{
   double fdreturn = par[0]/TMath::Power((1+TMath::Exp(-(x[0]-par[1])*par[2])),par[3]);
   //double fdreturn = par[0]/(1+par[3]*TMath::Exp(-(x[0]-par[1])*par[2]));
  
    return fdreturn;
}
double fermi_dirac_generalsub(double *x, double *par)
{
    double fdreturn = (par[0]/TMath::Power((1+TMath::Exp(-(x[0]-par[1])*par[2])),par[3])) - (par[0]/TMath::Power((1+TMath::Exp(-(x[0]-par[4])*par[5])),par[6]));
    //double fdreturn = (par[0]/(1+par[3]*TMath::Exp(-(x[0]-par[1])*par[2]))) - (par[0]/(1+par[6]*TMath::Exp(-(x[0]-par[4])*par[5])));
    return fdreturn;
}

double fermi_dirac_simple(double x)
{
    double fdreturn = (pars[0]/TMath::Power((1+TMath::Exp(-(x-pars[1])*pars[2])),pars[3])) - (pars[0]/TMath::Power((1+TMath::Exp(-(x-pars[4])*pars[5])),pars[6]));
    //double fdreturn = (par[0]/(1+par[3]*TMath::Exp(-(x[0]-par[1])*par[2]))) - (par[0]/(1+par[6]*TMath::Exp(-(x[0]-par[4])*par[5])));
    return fdreturn;
}




void Detector::FindRiseTime()
{
    double eightypercent_naive_time;
    double cf = 0.8*global_maximum.y;
    for(int i = global_maximum.position ; i >= start_point.position && i > 1; --i)
    {
        if( waveform_y.at(i) - cf > 0 && waveform_y.at(i-1) - cf < 0 )
        {
            double x1 = waveform_x.at(i-1);
            double x2 = waveform_x.at(i);
            double y1 = waveform_y.at(i-1);
            double y2 = waveform_y.at(i);
            eightypercent_naive_time = linear_interpolation(x1,x2,y1,y2,cf);
            break;
        }
    }
    rise_time = eightypercent_naive_time - naive_time;
}


void Detector::MicroSigmoid()
{
  /*
   int peakpos,startpos,difpeakpos,dif2peakpos;
   double peakval,valatdifpeak,valatdif2peak,peak,start,difpeak,dif2peak,sigtime;
std::vector<Double_t> wavedif,wave2dif;

 
       //Find PEAK and Startpoint
       int k=0;
       
       peakval=0;
       for(int i=0;i<waveform_y.size();i++){
	 if(peakval<waveform_y.at(i)){
	   peakval=waveform_y.at(i);
	   peakpos=i;
	 }
       }

       k=peakpos;


       while(waveform_y.at(k)>0){
         startpos=k;
         k--;

       }
               

       if(startpos<1000){
	 startpos=1210;
	 peakpos=3210;
       }
        
       peak=waveform_x.at(peakpos);
       start=waveform_x.at(startpos);

       
        //Differentiate Waveform (PEAK is the center of the slope)
       int cut=50;
       double yl[cut],xl[cut],par[2];
	 double peakdif=0;

       for(int i=startpos-100; i< peakpos+100; i++){
	 for(int k=0; k<cut; k++){
	   yl[k]=waveform_y.at(i-(cut/2)+k);
	   xl[k]= i-(cut/2)+k ;
	 }

	 TGraph gr(cut, xl, yl);
	 TF1 f1("f1", "[0]+[1]*x",i-(cut/2) ,i+(cut/2));
	 f1.SetParameters(waveform_y.at(i),waveform_y.at(i)-waveform_y.at(i+1));
	 //gr.Fit(f1,"R+");
	 gr.Fit("f1", "S", "",i-(cut/2) ,i+(cut/2));

	 f1.GetParameters(&par[0]);
	 wavedif.push_back(par[1]);
	 	      

	 
	 if(peakdif<par[1]){
	   peakdif=par[1];
	   difpeakpos=i; 
 	 }
       }

            
	      difpeakpos>peakpos? difpeakpos=peakpos: difpeakpos=difpeakpos;
	      
       valatdifpeak=waveform_y.at(difpeakpos);
       difpeak=waveform_x.at(difpeakpos);
   
       //2nd differentiation to find starting and endpoint

           int cut2=2;
       
       double yl2[cut2],xl2[cut2],par2[2];
       double peakdif2=0;
       
       for(int i=cut2/2; i< wavedif.size()-(cut2/2); i++){
	 for(int k=0; k<cut2; k++){
	   yl2[k]=wavedif.at(i-(cut2/2)+k);
	   xl2[k]= i-(cut2/2)+k ;
	 }

	

	 TGraph gr2(cut2, xl2, yl2);

	 TF1 f12("f12", "[0]+[1]*x",i-(cut2/2) ,i+(cut2/2));
	 f12.SetParameters(wavedif.at(i),wavedif.at(i)-wavedif.at(i+1));
	 // gr2.Fit(f12,"R+");
	 gr2.Fit("f12", "S", "",i-(cut2/2) ,i+(cut2/2));
	 
	 f12.GetParameters(&par2[0]);
	 wave2dif.push_back(par2[1]);
	 
       }
       
       double peak2dif=0;
       for(int i=difpeakpos-startpos+100; i<wave2dif.size();i++){
	 if(peak2dif<wave2dif.at(i)){
	   peak2dif=wave2dif.at(i);
	   dif2peakpos=i; 
	 }
       }
       
       
       
       dif2peakpos=dif2peakpos+startpos-100;
       valatdif2peak=waveform_y.at(dif2peakpos);
       dif2peak=waveform_x.at(dif2peakpos);
       
       



	 //create tgrapherror
	 
       int wavesize=waveform_y.size();
       double wx[wavesize], wy[wavesize],wxerr[wavesize],wyerr[wavesize];
       
       for(int i=0; i<wavesize;i++){    
	 wx[i]=i;
	 wy[i]=waveform_y.at(i);
       	 wxerr[i]=0;
	 wyerr[i]=baseline_rms;
	 }
       
       // TGraphErrors *waveerr= new TGraph(wavesize,wx,wy);

       TGraphErrors waveerr(wavesize,wx,wy,wxerr,wyerr);
   
       //Sigmoid fit
       double ps[4];

       startpos=startpos-50;
       dif2peakpos=dif2peakpos+0;
       
       TF1 f2("f2","([0]/(1+exp(-[1]*(x-[2]))))-[3]",startpos,dif2peakpos);
       f2.SetParameters(waveform_y.at(difpeakpos), 1,difpeakpos ,0);
      waveerr.Fit("f2", "S", "",startpos,dif2peakpos);
       
       
       ps[0] = f2.GetParameter(0);
       ps[1] = f2.GetParameter(1);
       ps[2] = f2.GetParameter(2);
       ps[3] = f2.GetParameter(3);


       //scale time

       double timepoint;
       
       timepoint=waveform_x.at(0)+ps[2]*(sqrt((waveform_x.at(2)-waveform_x.at(1))*(waveform_x.at(2)-waveform_x.at(1))));

       Sigmoid.microbulk= timepoint;
  Sigmoid.micropeakpos= dif2peakpos;
  Sigmoid.micropeak=waveform_y.at(dif2peakpos);
*/
}

std::fstream myfile2;
void Detector::Findmultiplepeak()
{
  /*
    myfile2.open ("multipeak.txt",std::fstream::app);
    std::cout<<"One"<<std::endl;
  int count1=0;
  int i=0;
   do{
  
     if(waveform_y.at(i)>0){
       for(int j=i;j < waveform_x.size(); j++){
	 if(waveform_y.at(j)<0){
	   break;
	   //  std::cout<<"One"<<i<<" "<<count1<<std::endl;
	 } 
	 count1=j;
       }
       
        if((count1-i)>8){
	 double maxloc=0;
	 //	 int maxlocpos=0;
	 //	 double maxloctime=0;
	 for(int k=i;k<count1;k++){
	   if(waveform_y.at(k)>maxloc){
	     maxloc=waveform_y.at(k);
	     //     maxlocpos=k;
	     //	     maxloctime=waveform_x.at(k);
	   }
	    }
	   myfile2<<maxloc<<std::endl;
	  std::cout<<maxloc<<" MAXLOC"<<std::endl;
	  	   }
      i=count1+8;
     //  std::cout<<i<<std::endl;
     }
     i++;
     if(i>=waveform_x.size())break;
   }while(i< waveform_x.size());
   myfile2.close();
  */
  
}


//For multiple inverse peaks like laser trigger 
void Detector::Findmultipleinverse()
{
  /* double test1=0;
  double timefirst=0;
    double timesec=0;
    int end=0;
    int Npoints=15;

   for(int i=0;i < waveform_x.size()/2; i++){
     if(waveform_y.at(i)<test1){
       test1=waveform_y.at(i);
       timefirst=waveform_x.at(i);
       end=i;
     }
   }

   int start=end-15;

 double x[Npoints], y[Npoints], erx[Npoints], ery[Npoints];
    for(int i = start; i <= end; ++i)
    {
        x[i-start] = waveform_x.at(i);
        y[i-start] = waveform_y.at(i);
        erx[i-start] = 0;
        ery[i-start] = baseline_rms;
        
       // std::cout<<waveform_x.at(i)<<" "<<waveform_y.at(i)<<"  "<<baseline_rms<<std::endl;
    }




    TGraphErrors waveform_graph(Npoints, x, y, erx, ery);

    double pars[4];
    pars[0]=test1;
    pars[1]=(x[end-start] + x[1])/2;
    pars[2]=5./(x[end-start] - x[1]);
    pars[3]=0.;
    
    TF1 fd_fit("fd_fit", fermi_dirac, x[0], x[Npoints-1], 4);
    fd_fit.SetParameters(pars[0],pars[1],pars[2],pars[3]);
    waveform_graph.Fit("fd_fit","qR");
    waveform_graph.Fit("fd_fit","qR");
    waveform_graph.Fit("fd_fit","qR");
    waveform_graph.Fit("fd_fit","qR");

    
    inversetime1=fd_fit.GetParameter(1)-(TMath::Log((fd_fit.GetParameter(0)/(0.2*test1-fd_fit.GetParameter(3)))-1)/fd_fit.GetParameter(2));
    // std::cout<<fd_fit.GetParameter(1)-(TMath::Log((fd_fit.GetParameter(0)/(0.2*test1-fd_fit.GetParameter(3)))-1)/fd_fit.GetParameter(2))<<std::endl;


    double test2=0;
    test1=0;
    end=0;
   for(int i=waveform_x.size()/2;i < waveform_x.size(); i++){
     if(waveform_y.at(i)<test1){
       test2=waveform_y.at(i);
       timesec=waveform_x.at(i);
       end=i;

     }
   }

   start=end-15;
 for(int i = start; i <= end; ++i)
    {
        x[i-start] = waveform_x.at(i);
        y[i-start] = waveform_y.at(i);
        erx[i-start] = 0;
        ery[i-start] = baseline_rms;
        
       // std::cout<<waveform_x.at(i)<<" "<<waveform_y.at(i)<<"  "<<baseline_rms<<std::endl;
    }
 
  TGraphErrors waveform_graph2(Npoints, x, y, erx, ery);

  double pars2[4];
   pars2[0]=test2;
   pars2[1]=(x[end-start] + x[1])/2;
   pars2[2]=5./(x[end-start] - x[1]);
   pars2[3]=0.;
    
    TF1 fd_fit2("fd_fit2", fermi_dirac, x[0], x[Npoints-1], 4);
   fd_fit2.SetParameters(pars2[0],pars2[1],pars2[2],pars2[3]);
   waveform_graph2.Fit("fd_fit2","qR");
   waveform_graph2.Fit("fd_fit2","qR");
   waveform_graph2.Fit("fd_fit2","qR");
   waveform_graph2.Fit("fd_fit2","qR");

    
    invtime2=fd_fit2.GetParameter(1)-(TMath::Log((fd_fit2.GetParameter(0)/(0.2*test1-fd_fit2.GetParameter(3)))-1)/fd_fit2.GetParameter(2));

    // std::cout<<fd_fit2.GetParameter(1)-(TMath::Log((fd_fit2.GetParameter(0)/(0.2*test2-fd_fit2.GetParameter(3)))-1)/fd_fit2.GetParameter(2))<<std::endl;
   //std:: cout<<"X"<<std::endl;
    //invtime2=timefirst;
    */
}
  
void Detector::FindWidth()
{
    double left_time;
    double right_time;
    double cf = 0.5*global_maximum.y;
    for(int i = global_maximum.position ; i >= start_point.position; --i)
    {
        if( waveform_y.at(i) - cf > 0 && waveform_y.at(i-1) - cf < 0 )
        {
            double x1 = waveform_x.at(i-1);
            double x2 = waveform_x.at(i);
            double y1 = waveform_y.at(i-1);
            double y2 = waveform_y.at(i);
            left_time = linear_interpolation(x1,x2,y1,y2,cf);
            break;
        }
    }

    
    bool ion_tail_is_too_high = 1;
    for(int i = global_maximum.position ; i <= e_peak_end.position; ++i)
    {
        if( waveform_y.at(i) - cf > 0 && waveform_y.at(i+1) - cf < 0 )
        {
            double x1 = waveform_x.at(i);
            double x2 = waveform_x.at(i+1);
            double y1 = waveform_y.at(i);
            double y2 = waveform_y.at(i+1);
            right_time = linear_interpolation(x1,x2,y1,y2,cf);
            ion_tail_is_too_high = 0;
            break;
        }
    }
    if(!ion_tail_is_too_high)
    {
        width = right_time - left_time;
    }
    else
    {
        std::cout << "ERROR*** Ion tail is too high and width cannot be properly calculated" << std::endl;
        width  = e_peak_end.x - left_time;
    }
}


void Detector::FilterWaveformFFT(int start, int N, double biggest_frequency_not_to_cut_GHz)
{
    *pre_filter_backup = *this;
    pre_filter_backup->pre_filter_backup = 0;
    
    double temp_wave[N];
    double temp_x[N];
    int end = start + N;
    double delta = waveform_x.at(1) -  waveform_x.at(0);
    int waveform_size = waveform_y.size();
    for(int i = start ; i < end; ++i)
        temp_wave[i-start] = waveform_y.at(i%waveform_size);


    TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N, "R2C M");
    fft_forward->SetPoints(temp_wave);
    fft_forward->Transform();

    
    double fourier_real[N];
    double fourier_imag[N];
    fft_forward->GetPointsComplex(fourier_real, fourier_imag);

    for(int i = 0; i < N; ++i)
    {   
        double frequency = i/delta/N;

        if(frequency > biggest_frequency_not_to_cut_GHz)// Cut at this frequency
        {
            fourier_real[i] = 0;
            fourier_imag[i] = 0;
        }
        //double filter = 1.;
        //if(frequency > 0.8)
        //    filter = TMath::Exp(2.08-2.6*frequency);
        //fourier_real[i] *= filter;
        //fourier_imag[i] *= filter;
    }

    TVirtualFFT *fft_backward = TVirtualFFT::FFT(1, &N, "C2R M");
    fft_backward->SetPointsComplex(fourier_real, fourier_imag);
    fft_backward->Transform();

    std::vector<double> temp_waveform_y;
    std::vector<double> temp_waveform_x;
    for(int i = 0; i < N; ++i)
    {
        temp_waveform_y.push_back(fft_backward->GetPointReal(i)/N);
        temp_waveform_x.push_back(waveform_x.at(start)+i*delta);
    }
    waveformFFT_y = temp_waveform_y;
    waveformFFT_x = temp_waveform_x;
}

void Detector:: FFTWaveform(){
  //std::cout<<"FFTWAVEFORM"<<std::endl;
  
  /*
  int start = start_point.position - 100; 

  int end =start_point.position +500;

  int Npoints = end-start+1;
  
  double x[Npoints], y[Npoints], erx[Npoints], ery[Npoints];
  for(int i = start; i <= end; ++i)
    {
      x[i-start] = waveform_x.at(i);
      y[i-start] = waveform_y.at(i);
      erx[i-start] = 0;
      ery[i-start] = baseline_rms;
      
      // std::cout<<waveform_x.at(i)<<" "<<waveform_y.at(i)<<"  "<<baseline_rms<<std::endl;
    }
  
  TGraphErrors waveform_graph(Npoints, x, y, 0, 0);


  TCanvas *c1 = new TCanvas("c1","c1",1600,1000);
  c1->Divide(2,1);
  c1->cd(1);
  waveform_graph.Draw();
  c1->Update();

  
   double delta = waveform_x.at(1) -  waveform_x.at(0);

   TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &Npoints, "R2C");
   fftr2c->SetPoints(y);
   fftr2c->Transform();
   
     Double_t re[Npoints], im[Npoints];
   fftr2c->GetPointsComplex(re, im);
  
   for(int i = 0; i < Npoints; ++i)
    {   
        double frequency = i/delta/Npoints;

        if(frequency > 1.8)// Cut at this frequency in GHz
        {
            re[i] = 0;
            im[i] = 0;
        }
 
    }


   
 TVirtualFFT *fft_backward = TVirtualFFT::FFT(1, &Npoints, "C2R M");
    fft_backward->SetPointsComplex(re, im);
    fft_backward->Transform();

     double  yFFT[Npoints];

     for(int i = 0; i < Npoints; ++i)
    {
      yFFT[i]=fft_backward->GetPointReal(i)/Npoints;
       
    }
    
  
  TGraphErrors waveform_graph2(Npoints, x, yFFT, 0, 0);

   
   c1->cd(2);
   waveform_graph2.Draw();
   c1->Update();
   c1 -> SaveAs("FFT.pdf");

   

*/
   
}

/* unsafe
void Detector::FilterWaveformFFT_test()
{
    *pre_filter_backup = *this;
    pre_filter_backup->pre_filter_backup = 0;
    
    int N = waveform_y.size();
    //round to next power of two
    N--;
    N |=  N >> 1;
    N |=  N >> 2;
    N |=  N >> 4;
    N |=  N >> 1;
    N |=  N >> 16;
    N++;
    //================
    double temp_wave[N];
    double temp_x[N];
    double delta = waveform_x.at(1) -  waveform_x.at(0);

    for(int i = 0 ; i < waveform_y.size(); ++i)
        temp_wave[i] = waveform_y.at(i);
    for(int i = waveform_y.size() ; i < N; ++i)
        temp_wave[i] = 0.;

    TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N, "R2C M");
    fft_forward->SetPoints(temp_wave);
    fft_forward->Transform();
    
    double fourier_real[N];
    double fourier_imag[N];
    fft_forward->GetPointsComplex(fourier_real, fourier_imag);

    for(int i = 0; i < N; ++i)
    {   
        double frequency = i/delta/N;
        if(frequency > 1)//in GHz
        {
            fourier_real[i] = 0;
            fourier_imag[i] = 0;
        }
    }

    TVirtualFFT *fft_backward = TVirtualFFT::FFT(1, &N, "C2R M");
    fft_backward->SetPointsComplex(fourier_real, fourier_imag);
    fft_backward->Transform();

    waveform_y.clear();
    for(int i = 0; i < waveform_x.size(); ++i)
    {
        waveform_y.push_back(fft_backward->GetPointReal(i)/N);
    }
}
*/
