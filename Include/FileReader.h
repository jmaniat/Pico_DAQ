#ifndef __filereader_h__
#define __filereader_h__

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>

#include "Detector.h"
#include "DetectorSetup.h"

class FileReader
{
    public:
        FileReader(const char* a);
        ~FileReader();
        bool GetNextEvent(int n, int &TrigNum);
        void SetDetectorSetup(DetectorSetup &det);
    private:
        DetectorSetup* Detectors;
        FILE* infile;
        int N_detectors = 3;
        double time_step = 0.05; //ns
        bool GetNextLine(double& p1, double& p2, double& p3);
};
#endif
