#ifndef __trc_filereader_h__
#define __trc_filereader_h__

#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <string.h>

#include "Detector.h"
#include "DetectorSetup.h"
#include "TRC_channel.h"

class TRC_FileReader
{
    public:
        TRC_FileReader(std::vector<int> channel_ids, const char* dir_path, int max_n_files, const char* file_midname, int first_file_counter = 0);
        bool GetNextEvent();
        void SetDetectorSetup(DetectorSetup &det);
        void OpenTriggerChannel(const char* dir_path, int max_n_files, const char* file_midname, int ch_id);
        int CalculateTriggerNumber();
    private:
        int First_File;
        std::vector<TRC_channel> Channels;
        bool trigger_channel_is_open = 0;
        TRC_channel *TriggerChannel;
        DetectorSetup* Detectors;
};
#endif




