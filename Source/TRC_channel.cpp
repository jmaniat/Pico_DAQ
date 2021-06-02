#include "TRC_channel.h"
#include <iostream>
#include <fstream>
#include <sstream>

TRC_channel::TRC_channel(const char* dir, int n, int ch, const char* file_midname, int first_file_counter )
{
    max_number_of_files = n;
    directory_path = std::string(dir);
    std::string middle_name = std::string(file_midname);
    filename_header = std::string("C") + std::to_string(ch) + middle_name;
    if( directory_path.back() != '/' )
        directory_path.push_back('/');


    //reserve vector sizes to prevent reallocations
    //triggerTime.reserve(100);
    // triggerOffset.reserve(100);
    // wave.reserve(400200);
    // waveform_x.reserve(4002);
    // waveform_y.reserve(4002);

    triggerTime.reserve(200);
    triggerOffset.reserve(200);
    wave.reserve(800400);
    waveform_x.reserve(8004);
    waveform_y.reserve(8004);

    file_counter = first_file_counter+0;
    event_number = 0;
    std::cout << "First file path of channel: " << directory_path << filename_header;
    std::cout  << "00000.trc = " << file_counter << std::endl;
    OpenFile(file_counter);
}


std::vector<double> TRC_channel::GetWaveformX() { return waveform_x; } 
std::vector<double> TRC_channel::GetWaveformY() { return waveform_y; } 
std::fstream myfile;
std::fstream myfilex;
bool TRC_channel::GetNextEvent()
{
  //std::cout<<event_number<<"  event nbr     nbr of events  "<<Number_of_events_in_file<<std::endl;
    if(event_number >= Number_of_events_in_file) 
    {

        event_number = 0;
        file_counter++;
        if(file_counter >= max_number_of_files)
            return 1;
        OpenFile(file_counter);
    }

    waveform_x.clear();
    waveform_y.clear();

    //  myfilex.open ("waveformx.txt",std::fstream::app);
    //  myfile.open ("waveform.txt",std::fstream::app);

    int waveform_size = desc.wave_array_count/desc.subarray_count;
    // std::cout<<triggerOffset.at(event_number)<<std::endl;
    //    std::cout<<(0*desc.horizontal_interval + triggerOffset.at(event_number))<<"  "<<(wave.at(0+event_number*waveform_size)*desc.vertical_gain - desc.vertical_offset)<<std::endl;

    for(int i = 0; i < waveform_size; ++i)
      {
	//waveform_x.push_back(i*desc.horizontal_interval);
	

	 
	//   myfile<<(i*desc.horizontal_interval + triggerOffset.at(event_number))<<"  "<<(wave.at(i+event_number*waveform_size)*desc.vertical_gain - desc.vertical_offset)<<std::endl;
	//	std::cout<<(i*desc.horizontal_interval + triggerOffset.at(event_number))<<"  "<<(wave.at(i+event_number*waveform_size)*desc.vertical_gain - desc.vertical_offset)<<std::endl;


	waveform_x.push_back(i*desc.horizontal_interval + triggerOffset.at(event_number));
	waveform_y.push_back(wave.at(i+event_number*waveform_size)*desc.vertical_gain - desc.vertical_offset);

	

//        myfile<<wave.at(i+event_number*waveform_size)<<" time "<<desc.vertical_offset<<"  minus  "<<desc.vertical_gain<<std::endl;
// myfile<<wave.at(i+event_number*waveform_size)*desc.vertical_gain - desc.vertical_offset<<std::endl;
// myfilex<<i*desc.horizontal_interval + triggerOffset.at(event_number)<<std::endl;
    }
    //       myfile<<"FAFA"<<std::endl;
	 //   myfilex<<"FAFA"<<std::endl;
    //        myfile.close();
    //    myfilex.close();
     
    event_number++;
    return 0;
}
int TRC_channel::OpenFile(int i)
{

    std::string file_counter_string = std::to_string(i);
    file_counter_string = std::string(5 - file_counter_string.length(), '0') + std::to_string(file_counter);
    std::string full_path = directory_path + filename_header + file_counter_string + std::string(".trc");
    triggerTime.clear();
    triggerOffset.clear();
    wave.clear();

    FILE* file;
    file = fopen(full_path.c_str(),"r");
    if(int tmp = readLeCroyBinary(file, &desc, &triggerTime, &triggerOffset, &wave))
    {
        std::cerr << "Unknown format: "<< tmp << std::endl << std::endl;
        return 1;
    }

    Number_of_events_in_file = desc.subarray_count;

    //printWavedesc(desc);
    fclose(file);

    return 0;

}

std::fstream lecroyfile;

int TRC_channel::readLeCroyBinary(FILE * file, wavedesc * desc, std::vector<double> * triggerTime, std::vector<double> * triggerOffset, std::vector<int> * wave)
{

  //lecroyfile.open ("binary_to.txt",std::fstream::app | std::fstream::out);
    wavedesc header;
    uint8_t tmp[9];
    do{

        if(fread(tmp, 8, 1, file) != 1) return 1;
        tmp[8] = 0;
        fseek(file, -7, SEEK_CUR);
    } while(strcmp((char*)tmp, "WAVEDESC") != 0);
    fseek(file, -1, SEEK_CUR);
    
    //BEGIN OF WAVEDESC BLOCK

    if(fread(&header, sizeof(wavedesc), 1, file) != 1) return 2;
    *desc = (wavedesc)header;
    //END OF WAVEDESC BLOCK
    
    //BEGIN OF USERTEXT BLOCK

    fseek(file, header.user_text, SEEK_CUR);
    //END OF USERTEXT BLOCK
    
    //BEGIN OF TRIGTIME ARRAY BLOCK

    for(int ii = 0; ii < header.subarray_count && header.subarray_count > 1; ii++){
        double val0, val1;
        if(ffread(header.comm_order1, &val0, 8, file) != 1) return 3;
        if(ffread(header.comm_order1, &val1, 8, file) != 1) return 4;
        triggerTime -> push_back(val0);
        triggerOffset -> push_back(val1);

	//	std::cout<<"triggerTime  "<<val0<<"   triggerOffset   "<<val1<<std::endl;
	//	lecroyfile<<val0<< "   "<<val1<<std::endl;
	
    }
    //END OF TRIGTIME ARRAY BLOCK
    
    //BEGIN OF RISTIME ARRAY BLOCK

    fseek(file, header.ris_time_array, SEEK_CUR);
    //END OF RISTIME ARRAY BLOCK
    
    //BEGIN OF DATA_ARRAY_1 ARRAY BLOCK

    
    for(int ii = 0; ii < header.wave_array_count; ii++){
        int comm_type;
        header.comm_order1 == 1 ? comm_type = header.comm_type1 : comm_type = header.comm_type2;
        switch(comm_type){
            case 0:
                int8_t val0;
                if(ffread(header.comm_order1, &val0, comm_type + 1, file) != 1) return 5;
                //if(fread(&val0, comm_type + 1, 1, file) != 1) return 1;
                wave -> push_back((int)val0);
		//	std::cout<<val0<<std::endl;
		break;
            case 1:
                int16_t val1;
		if(ffread(header.comm_order1, &val1, comm_type + 1, file) != 1) return 6;
                //if(fread(&val1, comm_type + 1, 1, file) != 1) return 1;
                wave -> push_back((int)val1);
		//		std::cout<<val1<<std::endl;
		//			lecroyfile<<val1<<std::endl;
		break;
        }
    }

    //END OF DATA_ARRAY_1 ARRAY BLOCK
    
    //BEGIN OF DATA_ARRAY_2 ARRAY BLOCK
    fseek(file, header.wave_array_2, SEEK_CUR);

    //END OF DATA_ARRAY_2 ARRAY BLOCK
    
    return 0;
}

int TRC_channel::ffread(unsigned short int byteOrder, void * val, unsigned int nByte, FILE * fileIn)
{

    const uint16_t one = 1;
    uint8_t ii;
    uint8_t * cc;
    uint8_t * tmp;
    if(fread(val, nByte, 1, fileIn) != 1) return -1;
    if(* ((uint8_t *)(&one)) == byteOrder) return 1;
    else{

        cc = (uint8_t *)val;
        tmp = (uint8_t *)malloc(nByte);
        for(ii = 0; ii < nByte; ii++){
            tmp[ii] = cc[ii];
        }
        for(ii = 0; ii < nByte; ii++){
            cc[ii] = tmp[nByte - ii - 1];
        }
        free(tmp);
        return 0;
    }
}

void TRC_channel::printWavedesc(wavedesc desc)
{
    std::cerr << "descriptor_name\t" << desc.descriptor_name << "\n";
    std::cerr << "template_name\t" << desc.template_name << "\n";
    std::cerr << "comm_type1\t" << (int)desc.comm_type1 << "\n";
    std::cerr << "comm_type2\t" << (int)desc.comm_type2 << "\n";
    std::cerr << "comm_order1\t" << (int)desc.comm_order1 << "\n";
    std::cerr << "comm_order2\t" << (int)desc.comm_order2 << "\n";
    std::cerr << "wave_descriptor\t" << desc.wave_descriptor << "\n";
    std::cerr << "user_text\t" << desc.user_text << "\n";
    std::cerr << "res_desc1\t" << desc.res_desc1 << "\n";
    std::cerr << "trigtime_array\t" << desc.trigtime_array << "\n";
    std::cerr << "ris_time_array\t" << desc.ris_time_array << "\n";
    std::cerr << "res_array1\t" << desc.wave_array_1 << "\n";
    std::cerr << "wave_array_1\t" << desc.wave_array_1 << "\n";
    std::cerr << "wave_array_2\t" << desc.wave_array_2 << "\n";
    std::cerr << "res_array2\t" << desc.res_array2 << "\n";
    std::cerr << "res_array3\t" << desc.res_array3 << "\n";
    std::cerr << "instrument_name\t" << desc.instrument_name << "\n";
    std::cerr << "instrument_number\t" << desc.instrument_number << "\n";
    std::cerr << "trace_label\t" << desc.trace_label << "\n";
    std::cerr << "reserved1\t" << desc.reserved1 << "\n";
    std::cerr << "reserved2\t" << desc.reserved2 << "\n";
    std::cerr << "wave_array_count\t" << desc.wave_array_count << "\n";
    std::cerr << "pnts_per_screen\t" << desc.pnts_per_screen << "\n";
    std::cerr << "first_valid_pnt\t" << desc.first_valid_pnt << "\n";
    std::cerr << "last_valid_pnt\t" << desc.last_valid_pnt << "\n";
    std::cerr << "first_point\t" << desc.first_point << "\n";
    std::cerr << "sparsing_factor\t" << desc.sparsing_factor << "\n";
    std::cerr << "segment_index\t" << desc.segment_index << "\n";
    std::cerr << "subarray_count\t" << desc.subarray_count << "\n";
    std::cerr << "sweeps_per_acq\t" << desc.sweeps_per_acq << "\n";
    std::cerr << "points_per_pair\t" << desc.points_per_pair << "\n";
    std::cerr << "pair_offset\t" << desc.pair_offset << "\n";
    std::cerr << "vertical_gain\t" << desc.vertical_gain << "\n";
    std::cerr << "vertical_offset\t" << desc.vertical_offset << "\n";
    std::cerr << "max_value\t" << desc.max_value << "\n";
    std::cerr << "min_value\t" << desc.min_value << "\n";
    std::cerr << "nominal_bits\t" << desc.nominal_bits << "\n";
    std::cerr << "nom_subarray_count\t" << desc.nom_subarray_count << "\n";
    std::cerr << "horizontal_interval\t" << desc.horizontal_interval << "\n";
    std::cerr << "horizontal_offset\t" << desc.horizontal_offset << "\n";
    std::cerr << "pixel_offset\t" << desc.pixel_offset << "\n";
    std::cerr << "vertunit\t" << desc.vertunit << "\n";
    std::cerr << "horunit\t" << desc.horunit << "\n";
    std::cerr << "horiz_uncertainty\t" << desc.horiz_uncertainty << "\n";
    std::cerr << "trigger_time\t" << desc.trigger_time_Y << "/" << (int)desc.trigger_time_M << "/" << (int)desc.trigger_time_D << " " << (int)desc.trigger_time_h << ":" << (int)desc.trigger_time_m << ":";
    fprintf(stderr, "%.13lf\n", desc.trigger_time_s);
    std::cerr << "trigger_time_r\t" << desc.trigger_time_r << "\n";
    std::cerr << "acq_duration\t" << desc.acq_duration << "\n";
    std::cerr << "record_type\t" << desc.record_type << "\n";
    std::cerr << "processing_done\t" << desc.processing_done << "\n";
    std::cerr << "reserved5\t" << desc.reserved5 << "\n";
    std::cerr << "ris_sweeps\t" << desc.ris_sweeps << "\n";
    std::cerr << "timebase\t" << desc.timebase << "\n";
    std::cerr << "vert_coupling\t" << desc.vert_coupling << "\n";
    std::cerr << "probe_att\t" << desc.probe_att << "\n";
    std::cerr << "fixed_vert_gain\t" << desc.fixed_vert_gain << "\n";
    std::cerr << "bandwidth_limit\t" << desc.bandwidth_limit << "\n";
    std::cerr << "vertical_vernier\t" << desc.vertical_vernier << "\n";
    std::cerr << "acq_vert_offset\t" << desc.acq_vert_offset << "\n";
    std::cerr << "wave_source\t" << desc.wave_source << "\n";
    std::cerr << "\n\n";
}
