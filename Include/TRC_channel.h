#ifndef __trc_channel_h__
#define __trc_channel_h__

#include <stdio.h>
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

#pragma pack(push, 1)

typedef struct{
    uint8_t descriptor_name[16];
    uint8_t template_name[16];
    int8_t comm_type1;
    int8_t comm_type2;
    int8_t comm_order1;
    int8_t comm_order2;
    int32_t wave_descriptor;
    int32_t user_text;
    int32_t res_desc1;
    int32_t trigtime_array;
    int32_t ris_time_array;
    int32_t res_array1;
    int32_t wave_array_1;
    int32_t wave_array_2;
    int32_t res_array2;
    int32_t res_array3;
    uint8_t instrument_name[16];
    int32_t instrument_number;
    uint8_t trace_label[16];
    int16_t reserved1;
    int16_t reserved2;
    int32_t wave_array_count;
    int32_t pnts_per_screen;
    int32_t first_valid_pnt;
    int32_t last_valid_pnt;
    int32_t first_point;
    int32_t sparsing_factor;
    int32_t segment_index;
    int32_t subarray_count;
    int32_t sweeps_per_acq;
    int16_t points_per_pair;
    int16_t pair_offset;
    float vertical_gain;
    float vertical_offset;
    float max_value;
    float min_value;
    int16_t nominal_bits;
    int16_t nom_subarray_count;
    float horizontal_interval;
    double horizontal_offset;
    double pixel_offset;
    uint8_t vertunit[48];
    uint8_t horunit[48];
    float horiz_uncertainty;
    double trigger_time_s;
    int8_t trigger_time_m;
    int8_t trigger_time_h;
    int8_t trigger_time_D;
    int8_t trigger_time_M;
    int16_t trigger_time_Y;
    int16_t trigger_time_r;
    float acq_duration;
    int16_t record_type;
    int16_t processing_done;
    int16_t reserved5;
    int16_t ris_sweeps;
    int16_t timebase;
    int16_t vert_coupling;
    float probe_att;
    int16_t fixed_vert_gain;
    int16_t bandwidth_limit;
    float vertical_vernier;
    float acq_vert_offset;
    int16_t wave_source;
} wavedesc;

#pragma pack(pop)

class TRC_channel
{
    public:
        TRC_channel(const char* dir, int n, int ch, const char* file_midname, int first_file_counter);

        std::vector<double> GetWaveformX(); 
        std::vector<double> GetWaveformY(); 

        bool GetNextEvent();
        int OpenFile(int i);
    private:
        wavedesc desc;
        std::vector<double> triggerTime;
        std::vector<double> triggerOffset;
        std::vector<int> wave;

        std::vector<double> waveform_x;
        std::vector<double> waveform_y;

        std::string directory_path;
        std::string filename_header; // i.e. 'C2Seq100traces' of 'C2Seq100traces00000.trc'
        int event_number;
        int file_counter;
        int max_number_of_files;
        int Number_of_events_in_file;

        int readLeCroyBinary(FILE *, wavedesc *, std::vector<double> * , std::vector<double> * , std::vector<int> *);
        int ffread(unsigned short int byteOrder, void * val, unsigned int nByte, FILE * fileIn);
        void printWavedesc(wavedesc);
};



#endif
