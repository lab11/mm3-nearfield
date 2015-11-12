/* -*- c++ -*- */
/* 
 * Copyright 2014 <+YOU OR YOUR COMPANY+>.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define NTAPS 8

#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/io_signature.h>
#include <volk/volk.h>
#include "nearfield_demod_impl.h"
#include <numeric>
#include <iterator>
#include <queue>
#include <ctime>
#include <iostream>
#include <pthread.h>
namespace gr {
namespace nearfield {

void* rake_filter_process_helper(void *obj) {
    object *mm3 = (object *)obj;
    mm3 -> C -> rake_filter_process(mm3->start_num, mm3->end_num, mm3->thread_num);
    return NULL;
}
float mean(std::vector<float> in_vec){
	float sum = std::accumulate(in_vec.begin(), in_vec.end(), 0.0);
	float m = sum / in_vec.size();
	return m;
}

nearfield_demod::sptr nearfield_demod::make(float sample_rate, float bitrate, float bitrate_accuracy, float post_bitrate_accuracy, float pulse_len, float pulse_len_accuracy, float post_pulse_len_accuracy, int packet_len, int header_len, float bb_freq, const std::string gatd_id) {
		return gnuradio::get_initial_sptr
			(new nearfield_demod_impl(sample_rate, bitrate, bitrate_accuracy, post_bitrate_accuracy, pulse_len, pulse_len_accuracy, post_pulse_len_accuracy, packet_len, header_len, bb_freq, gatd_id));
	}

/*
 * The private constructor
 */
nearfield_demod_impl::nearfield_demod_impl(float sample_rate, float bitrate, float bitrate_accuracy, float post_bitrate_accuracy, float pulse_len, float pulse_len_accuracy, float post_pulse_len_accuracy, int packet_len, int header_len, float bb_freq, const std::string gatd_id)
	: gr::sync_block("nearfield_demod",
			gr::io_signature::make(1, 1, sizeof(float)),
			gr::io_signature::make(0, 1, sizeof(float))),
	  d_log_file("nearfield_log.txt"), d_gatd_id(gatd_id) {

	//Prepare FIR filter for internal use
	float taps[NTAPS] = {0.8, 1.34, 1.35, 1.38, 1.35, 0.91, 0.8, 0.91};
	d_fir = new gr::filter::kernel::fir_filter_fff(1, std::vector<float>(taps,taps+NTAPS));
	set_history(d_fir->ntaps());

	const int alignment_multiple = volk_get_alignment() / sizeof(float);
	set_alignment(std::max(1, alignment_multiple));

	// variables
	threshold = 8;            // threshold set after observing data
	setSampleRate(sample_rate);
	setPulseLen(pulse_len);
	setPulseLenAccuracy(pulse_len_accuracy);
	setPostPulseLenAccuracy(post_pulse_len_accuracy);
	setBitrate(bitrate);
	setBitrateAccuracy(bitrate_accuracy);
	setPostBitrateAccuracy(post_bitrate_accuracy);
	setPacketLen(packet_len);
	setHeaderLen(header_len);


	sample_ctr = 0;
	last_prf = 0.0;
	last_pulse = 0.0;
    sample_counter = 0;	
	
	// Treating this like it's streaming data so I won't use some of the 
	// efficient Matlab functions that would speed this up.
	last_pulse_length = 0;
	last_pulse_distance = 0;
	last_data = 0;      // previous received data point
	pulse_count = 0;    // length of current pulse (1s)
	prf_count = 0;      // length of current prf (0s)
	sync = 0;           // current consecutive sync pulses
	error = 0;
	fuzz = 1;
	last_prf_count = 0;
	valid_count = 0;
	sync_prf = 0;       // prf counter in sync part of loop
	sync_prf2 = 0;      // prf counter in sync part of loop (prf_window)
	sync_pulse = 0;     // pulse counter in sync part of loop
	prf_win_cnt = 0;    // counter for window of pulses between prf_min and prf_max in sync routine
	valid_pulse = 0;    // flag for pulse detection
	n = 0;              // counter for N
	max_sample = 0;
    subsample_rate = 250.0;
    data_window = 7.0;
	
	threshold_sync = 0.7;
	max_header_response = 0;
	last_peak_response = 0;
	process_counter = 0;
    time_offset = 0;
    data_acquired = 0;
	pos = 0;
    target_one_pos = 0;
    target_zero_pos = 0;


    for(int i = 0; i < 100; i++){
	    data_energy_0[i] = 0; 
		data_energy_1[i] = 0; 
	    data_pos_0[i] = 0; 
		data_pos_1[i] = 0; 
    }
    //unit_time = 103000.0/(16 * subsample_rate);
    //tuning params
    //v7
    //data_distance = 12033.0/(subsample_rate * 1.003);
    //header_data_distance = 11488.0/(subsample_rate * 1.003);
    //one_zero_rate = 0.185;
    //A = 546.6/float(subsample_rate * 1.003);
    //B = 7669.5/float(subsample_rate * 1.003);
    //dis_one_to_zero = 9300.0/(subsample_rate * 1.003);
    //dis_one_to_one = 11488.0/(subsample_rate * 1.003);
    //dis_zero_to_zero = 12033.0/(subsample_rate * 1.003);
    //dis_zero_to_one = 14200.0/(subsample_rate * 1.003);



    ////data_distance = 24300/subsample_rate;
    ////header_data_distance = 20500/subsample_rate;
    ////one_zero_rate = 0.185;
    ////A = 1106.0/float(subsample_rate);
    ////B = 15418/float(subsample_rate);
    ////dis_one_to_zero = 18500.0/subsample_rate;
    ////dis_one_to_one = 23166.0/subsample_rate;
    ////dis_zero_to_zero = 24300.0/subsample_rate;
    ////dis_zero_to_one = 29000.0/subsample_rate;
    //one_zero_rate = 0.185;
    //A = 1381.3/float(subsample_rate);
    //B = 19311/float(subsample_rate);
    //data_distance = 30360/subsample_rate;
    //header_data_distance = 26000/subsample_rate;
    //dis_one_to_zero = (29500 * 0.797)/subsample_rate;
    //dis_one_to_one = 29500/subsample_rate;
    //dis_zero_to_zero = 24250/subsample_rate;
    //dis_zero_to_one = (24250 * (1 + one_zero_rate))/subsample_rate;
    reset_data = 0;
    peak_distance = 0;
    //jitter = data_distance * 0.01;
    jitter = 1.0;
	delay = 0;
    //jitter = 1.1;
    //std::cout << jitter << std::endl;
	//jitter = 1;
    num_rake_filter = 30;

	sub_sample_counter = 0;	
	max_current = 0;
	avg_current = 0;
	last_time = time(0);
	pulse_vec.clear();
	prf_vec.clear();
	demod_data.clear();
	noise_power = 0;
    for (int k = 0; k < 100; k++ ) {
        lastpulses.push(0);
	}
    //init tables
    seed_table[0] = 0;
    seed_table[1] = 16;
    seed_table[2] = 32;
    seed_table[3] = 3;
    seed_table[4] = 6;
    seed_table[5] = 12;
    seed_table[6] = 24;
    seed_table[7] = 48;
    seed_table[8] = 35;
    seed_table[9] = 5;
    seed_table[10] = 10;
    seed_table[11] = 20;
    seed_table[12] = 40;
    seed_table[13] = 19;
    seed_table[14] = 38;
    seed_table[15] = 15;
 

	setbbfreq(bb_freq);

    /*
    distance_table[0] = 0;

    std::cout << "distance table" << std::endl;
    for(int i = 1; i < 16; i++){
        distance_table[i] = seed_table[i] * A + B;
        std::cout << distance_table[i] << std::endl;
    }

    sum_table[0] = 0;
    std::cout << "sum table" << std::endl;
    for(int i = 1; i < 16; i++){
        sum_table[i] = sum_table[i-1] + distance_table[i];
        std::cout << sum_table[i] << std::endl;
    }


	for (int k = 0; k < num_rake_filter; k++){

		//std::cout << "size of buffer[" << k << "] = " << (345 * unit_time * (1+unit_offset*k) + 2 * jitter + 345 * unit_offset * unit_time) << std::endl; 
 		if(k == num_rake_filter - 1){
            for(int j = 0; j < (sum_table[15] * (1+unit_offset*k) + 2 * jitter + sum_table[15] * unit_offset); j++){
			    matched_pulses.push_front(0);
            }
        }
        rake_offset[k] = (sum_table[15] * (1+unit_offset*(num_rake_filter-1)) + 2 * jitter + sum_table[15] * unit_offset) -
                            (sum_table[15] * (1+unit_offset*k) + 2 * jitter + sum_table[15] * unit_offset);


		all_pulse_energy[k] = 0;

		//std::cout << "size of buffer[" << k << "] = " << matched_pulses[k].size() << std::endl; 
	}

    */
    //std::cout << N << std::endl;
    /*
	for(int i = 0; i< num_rake_filter; i++){
		for(int j = 0; j < 16; j++) {
			if(j == 0) {
				window_length[i] = 2*jitter + 1;
			} else {
				window_length[i] += int((1+unit_offset*i) * sum_table[j] + (2 * jitter - 1) + 
						unit_offset * sum_table[j]) - int((1+unit_offset*i) * sum_table[j]);
			}
		}
	}
    */
    //std::cout << "jitter: " << jitter << std::endl;



    /*
    Objs_0.C = this;
    Objs_0.start_num =0;
    Objs_0.end_num =5;
    Objs_0.thread_num =0;
    Objs_1.C = this;
    Objs_1.start_num =5;
    Objs_1.end_num =10;
    Objs_1.thread_num =1;
    Objs_2.C = this;
    Objs_2.start_num =10;
    Objs_2.end_num =15;
    Objs_2.thread_num =2;
    Objs_3.C = this;
    Objs_3.start_num =15;
    Objs_3.end_num=20;
    Objs_3.thread_num =3;
    Objs_4.C = this;
    Objs_4.start_num =20;
    Objs_4.end_num =25;
    Objs_4.thread_num =4;
    Objs_5.C = this;
    Objs_5.start_num =25;
    Objs_5.end_num =30;
    Objs_5.thread_num =5;
    Objs_6.C = this;
    Objs_6.start_num =30;
    Objs_6.end_num =35;
    Objs_6.thread_num =6;
    Objs_7.C = this;
    Objs_7.start_num =35;
    Objs_7.end_num=40;
    Objs_7.thread_num =7;



    count = 0;
    //std::cerr << "start" << std::endl;
    shared_lock = PTHREAD_MUTEX_INITIALIZER;
    //pthread_mutex_lock(&shared_lock);
    //std::cerr << "1" << std::endl;
    shared_cond = PTHREAD_COND_INITIALIZER;
    //std::cerr << "2" << std::endl;



    locks_0 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_0);
    locks_1 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_1);
    locks_2 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_2);
    locks_3 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_3);
    locks_4 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_4);
    locks_5 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_5);
    locks_6 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_6);
    locks_7 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_7);



    //std::cerr << "3" << std::endl;
    
    irets[0] = pthread_create(&threads_0, NULL, &rake_filter_process_helper, &Objs_0);
    irets[1] = pthread_create(&threads_1, NULL, &rake_filter_process_helper, &Objs_1);
    irets[2] = pthread_create(&threads_2, NULL, &rake_filter_process_helper, &Objs_2);
    irets[3] = pthread_create(&threads_3, NULL, &rake_filter_process_helper, &Objs_3);
    
    irets[4] = pthread_create(&threads_4, NULL, &rake_filter_process_helper, &Objs_4);
    irets[5] = pthread_create(&threads_5, NULL, &rake_filter_process_helper, &Objs_5);
    irets[6] = pthread_create(&threads_6, NULL, &rake_filter_process_helper, &Objs_6);
    irets[7] = pthread_create(&threads_7, NULL, &rake_filter_process_helper, &Objs_7);

    //std::cerr << "4" << std::endl;



    */
	message_port_register_out(pmt::mp("frame_out"));

}

void nearfield_demod_impl::rake_filter_process(int start_num, int end_num, int thread_num) {
    int start_counter = start_num;
    int thread_counter = thread_num;
    //int *end_counter = (int *)end_num;

    //std::cout << "start thread: " << thread_counter << std::endl;


    while(1){
        //wait for mutex()

        if(thread_num == 0){
            pthread_mutex_lock(&locks_0);
        }else if(thread_num == 1){
            pthread_mutex_lock(&locks_1);
        }else if(thread_num == 2){
            pthread_mutex_lock(&locks_2);
        }else if(thread_num == 3){
            pthread_mutex_lock(&locks_3);
        }else if(thread_num == 4){
            pthread_mutex_lock(&locks_4);
        }else if(thread_num == 5){
            pthread_mutex_lock(&locks_5);
        }else if(thread_num == 6){
            //pthread_mutex_lock(&locks_6);
        }else if(thread_num == 7){
            //pthread_mutex_lock(&locks_7);
        }
            
        for(int k = 0; k < 16; k++) {
            //std::cout << start_counter << std::endl;
            for(int i = start_counter; i < end_num; i++){
            //for(int j = 0; j < matched_pulses[i].size(); j++){
                if(k == 0) {
                    float max_pulse = 0;     
                    for(int loc = rake_offset[i]; loc <= rake_offset[i] + 2 * jitter - 1; loc++){
                        
                       /* 
                        if(thread_num == 0 && i == 0) {
                            std::cout << "thread: " << thread_num << ", loc_0: " << loc << ", pulse: " << matched_pulses[loc] << std::endl;
                        }
                        */
                        
                        
                        if(matched_pulses[loc] > max_pulse) {
                            max_pulse = matched_pulses[loc];
                        }
                    }
                    long_matched_out[i] += max_pulse;


		    /*
		    if(matched_pulses[rake_offset[i] + 2 * jitter] > 1) {
			scores[i]++;
			}
		    */
                }
                else if(k > 0 && k <= 15) {
            
                    /*
                    long_matched_out[i] = long_matched_out[i] - 
                            matched_pulses[rake_offset[i]+int((1+unit_offset*i) * unit_time * sum_table[k] - 1)] * 
                            matched_pulses[rake_offset[i]+int((1+unit_offset*i) * unit_time * sum_table[k] - 1)] +
                            (matched_pulses[rake_offset[i]+int((1+unit_offset*i) * unit_time * sum_table[k] + (2 * jitter - 1) + 
                            unit_offset * unit_time * sum_table[k])]) * 
                            (matched_pulses[rake_offset[i]+int((1+unit_offset*i) * unit_time * sum_table[k] + (2 * jitter - 1) + 
                            unit_offset * unit_time * sum_table[k])]);

                    */

                    float max_pulse = 0;     
                    for(int loc = rake_offset[i]+int((1+unit_offset*i) * sum_table[k]); 
                            loc <= rake_offset[i]+int((1+unit_offset*i) * sum_table[k] + (2 * jitter - 1) + unit_offset * sum_table[k]);
                            loc++){
                       
                        
                        if(matched_pulses[loc] > max_pulse) {
                            max_pulse = matched_pulses[loc];

                            //pulse_loc_queue[i][k] = loc;
                            //pulse_queue[i][k] = max_pulse;
                        }
                        /*
                        if((i == 15 || i == 16)) {
                            std::cout << "filter: " << i << ", loc: " << loc << ", pulse: " << matched_pulses[loc] << std::endl;
                        }
                        */


                        
                        
                    }
                    long_matched_out[i] += max_pulse;
		    /*
		    if((matched_pulses[rake_offset[i]+int((1+unit_offset*i) * unit_time * sum_table[k] + (2 * jitter - 1) + 
                            unit_offset * unit_time * sum_table[k])]) > 1) {
			scores[i]++;
			}
		    */
                }
            }
	                
        }
        //pthread_cond_signal(&conds[thread_num]);
        //pthread_mutex_unlock(&locks[thread_num]);
        //Lock mutex()
        pthread_mutex_lock(&shared_lock);
        count = count + 1;
        if(count == (num_rake_filter/5)){
            pthread_cond_signal(&shared_cond);
        }
        pthread_mutex_unlock(&shared_lock);
    }
}


/*
 * Our virtual destructor.
 */
nearfield_demod_impl::~nearfield_demod_impl() {
	delete d_fir;

	d_log_file.close();
    pthread_mutex_destroy(&locks_0);
    pthread_mutex_destroy(&locks_1);
    pthread_mutex_destroy(&locks_2);
    pthread_mutex_destroy(&locks_3);
    pthread_mutex_destroy(&locks_4);
    pthread_mutex_destroy(&locks_5);
    pthread_mutex_destroy(&locks_6);
    pthread_mutex_destroy(&locks_7);
       //pthread_join(thread[i], NULL);

    pthread_mutex_destroy(&shared_lock);
    pthread_cond_destroy(&shared_cond);
           
}
float nearfield_demod_impl::getLastObservedBitrate(){
	return 1.0/last_prf;
}

float nearfield_demod_impl::getLastObservedPulseLen(){
	return last_pulse;
}

float nearfield_demod_impl::getThreshold(){
	return threshold;
}

void nearfield_demod_impl::setHeaderLen(int header_len_in){
	header = header_len_in;
}

void nearfield_demod_impl::setPacketLen(int packet_len_in){
	N = packet_len_in;
}

void nearfield_demod_impl::setbbfreq(float bb_freq_in){
	sample_ctr = 0;
	last_prf = 0.0;
	last_pulse = 0.0;
    sample_counter = 0;	
	
	// Treating this like it's streaming data so I won't use some of the 
	// efficient Matlab functions that would speed this up.
	last_pulse_length = 0;
	last_pulse_distance = 0;
	last_data = 0;      // previous received data point
	pulse_count = 0;    // length of current pulse (1s)
	prf_count = 0;      // length of current prf (0s)
	sync = 0;           // current consecutive sync pulses
	error = 0;
	fuzz = 1;
	last_prf_count = 0;
	valid_count = 0;
	sync_prf = 0;       // prf counter in sync part of loop
	sync_prf2 = 0;      // prf counter in sync part of loop (prf_window)
	sync_pulse = 0;     // pulse counter in sync part of loop
	prf_win_cnt = 0;    // counter for window of pulses between prf_min and prf_max in sync routine
	valid_pulse = 0;    // flag for pulse detection
	n = 0;              // counter for N
	max_sample = 0;
    subsample_rate = 250.0;
	
	threshold_sync = 0.7;
	max_header_response = 0;
	last_peak_response = 0;
	process_counter = 0;
    time_offset = 0;
    data_acquired = 0;
	pos = 0;
    target_one_pos = 0;
    target_zero_pos = 0;


    //unit_time = 103000.0/(16 * subsample_rate);
    //tuning params
    //v7
    //data_distance = 12033.0/(subsample_rate * 1.003);
    //header_data_distance = 11488.0/(subsample_rate * 1.003);
    //one_zero_rate = 0.185;
    //A = 546.6/float(subsample_rate * 1.003);
    //B = 7669.5/float(subsample_rate * 1.003);
    //dis_one_to_zero = 9300.0/(subsample_rate * 1.003);
    //dis_one_to_one = 11488.0/(subsample_rate * 1.003);
    //dis_zero_to_zero = 12033.0/(subsample_rate * 1.003);
    //dis_zero_to_one = 14200.0/(subsample_rate * 1.003);



    ////data_distance = 24300/subsample_rate;
    ////header_data_distance = 20500/subsample_rate;
    ////one_zero_rate = 0.185;
    ////A = 1106.0/float(subsample_rate);
    ////B = 15418/float(subsample_rate);
    ////dis_one_to_zero = 18500.0/subsample_rate;
    ////dis_one_to_one = 23166.0/subsample_rate;
    ////dis_zero_to_zero = 24300.0/subsample_rate;
    ////dis_zero_to_one = 29000.0/subsample_rate;
    //one_zero_rate = 0.185;
    //A = 1381.3/float(subsample_rate);
    //B = 19311/float(subsample_rate);
    //data_distance = 30360/subsample_rate;
    //header_data_distance = 26000/subsample_rate;
    //dis_one_to_zero = (29500 * 0.797)/subsample_rate;
    //dis_one_to_one = 29500/subsample_rate;
    //dis_zero_to_zero = 24250/subsample_rate;
    //dis_zero_to_one = (24250 * (1 + one_zero_rate))/subsample_rate;
    reset_data = 0;
    peak_distance = 0;
    //jitter = data_distance * 0.01;
    jitter = 1.0;
	delay = 0;
    //jitter = 1.1;
    //std::cout << jitter << std::endl;
	//jitter = 1;
    num_rake_filter = 30;

	sub_sample_counter = 0;	
	max_current = 0;
	avg_current = 0;
	last_time = time(0);
	pulse_vec.clear();
	prf_vec.clear();
	demod_data.clear();
	noise_power = 0;
    for (int k = 0; k < 100; k++ ) {
        lastpulses.push(0);
	}
    //init tables
    seed_table[0] = 0;
    seed_table[1] = 16;
    seed_table[2] = 32;
    seed_table[3] = 3;
    seed_table[4] = 6;
    seed_table[5] = 12;
    seed_table[6] = 24;
    seed_table[7] = 48;
    seed_table[8] = 35;
    seed_table[9] = 5;
    seed_table[10] = 10;
    seed_table[11] = 20;
    seed_table[12] = 40;
    seed_table[13] = 19;
    seed_table[14] = 38;
    seed_table[15] = 15;
 

	bb_offset = bb_freq_in;

    unit_offset = 0.0015;
    factor = 1 + unit_offset * bb_offset;//26.5
    //v8 fast
    /*
    data_distance = 11925.75/(subsample_rate * 1.003);
    header_data_distance = 11456.0/(subsample_rate * 1.003);
    one_zero_rate = 0.185;
    A = 542.4/float(subsample_rate * 1.003);
    B = 7589.7/float(subsample_rate * 1.003);
    dis_one_to_zero = 9223.75/(subsample_rate * 1.003);
    dis_one_to_one = 11386.75/(subsample_rate * 1.003);
    dis_zero_to_zero = 11925.875/(subsample_rate * 1.003);
    dis_zero_to_one = 14087.71/(subsample_rate * 1.003);
    */
    
    /*
    //v8 fast
    data_distance = 11925.75 /(subsample_rate * factor);
    header_data_distance = 11456.0 /(subsample_rate * factor);
    one_zero_rate = 0.185;
    A = 542.4 /float(subsample_rate * factor);
    B = 7589.7 /float(subsample_rate * factor);
    dis_one_to_zero = 9223.75 /(subsample_rate * factor);
    dis_one_to_one = 11386.75 /(subsample_rate * factor);
    dis_zero_to_zero = 11925.875 /(subsample_rate * factor);
    dis_zero_to_one = 14087.71 /(subsample_rate * factor);
    */


    //v8 slow
    data_distance = 11925.75 * 2/(subsample_rate * factor);
    header_data_distance = 11456.0 * 2/(subsample_rate * factor);
    one_zero_rate = 0.185;
    A = 542.4 * 2/float(subsample_rate * factor);
    B = 7589.7 * 2/float(subsample_rate * factor);
    dis_one_to_zero = 9223.75 * 2/(subsample_rate * factor);
    dis_one_to_one = 11386.75 * 2/(subsample_rate * factor);
    dis_zero_to_zero = 11925.875 * 2/(subsample_rate * factor);
    dis_zero_to_one = 14087.71 * 2/(subsample_rate * factor);
    

    distance_table[0] = 0;

    //std::cout << "distance table" << std::endl;
    for(int i = 1; i < 16; i++){
        distance_table[i] = seed_table[i] * A + B;
        //std::cout << distance_table[i] << std::endl;
    }

    sum_table[0] = 0;
    //std::cout << "sum table" << std::endl;
    for(int i = 1; i < 16; i++){
        sum_table[i] = sum_table[i-1] + distance_table[i];
        //std::cout << sum_table[i] << std::endl;
    }

	for (int k = 0; k < num_rake_filter; k++){

		//std::cout << "size of buffer[" << k << "] = " << (345 * unit_time * (1+unit_offset*k) + 2 * jitter + 345 * unit_offset * unit_time) << std::endl; 
 		if(k == num_rake_filter - 1){
            for(int j = 0; j < (sum_table[15] * (1+unit_offset*k) + 2 * jitter + sum_table[15] * unit_offset); j++){
			    matched_pulses.push_front(0);
            }
        }
        rake_offset[k] = (sum_table[15] * (1+unit_offset*(num_rake_filter-1)) + 2 * jitter + sum_table[15] * unit_offset) -
                            (sum_table[15] * (1+unit_offset*k) + 2 * jitter + sum_table[15] * unit_offset);


		all_pulse_energy[k] = 0;

		//std::cout << "size of buffer[" << k << "] = " << matched_pulses[k].size() << std::endl; 
	}
	energy = 0;
    for(int i = 0; i < num_rake_filter; i++){
		long_matched_out[i] = 0;     
        aggregated_header[i] = 0;
	//scores[i] = 0;
        last[i] = 0;
	}

        /*
        for(int j = 0; j < unit_time * 16 * 0.1; j++){
            data_queue[i].push_front(0);
        }
        */
        
    for(int j = 0; j < num_rake_filter; j++){
        for( int i = 0; i < 16; i++)
            pulse_queue[j].push_front(0);
            pulse_loc_queue[j].push_front(0);
    }
    

    data_energy_out = 0;
  
    last_offset = 0;
	correct_offset = 0;
	start = 0;
    last_max_response = 0;

    pthread_mutex_destroy(&locks_0);
    pthread_mutex_destroy(&locks_1);
    pthread_mutex_destroy(&locks_2);
    pthread_mutex_destroy(&locks_3);
    pthread_mutex_destroy(&locks_4);
    pthread_mutex_destroy(&locks_5);
    //pthread_mutex_destroy(&locks_6);
    //pthread_mutex_destroy(&locks_7);
    //pthread_join(thread[i], NULL);

    pthread_mutex_destroy(&shared_lock);
    pthread_cond_destroy(&shared_cond);

    Objs_0.C = this;
    Objs_0.start_num =0;
    Objs_0.end_num =5;
    Objs_0.thread_num =0;
    Objs_1.C = this;
    Objs_1.start_num =5;
    Objs_1.end_num =10;
    Objs_1.thread_num =1;
    Objs_2.C = this;
    Objs_2.start_num =10;
    Objs_2.end_num =15;
    Objs_2.thread_num =2;
    Objs_3.C = this;
    Objs_3.start_num =15;
    Objs_3.end_num=20;
    Objs_3.thread_num =3;
    Objs_4.C = this;
    Objs_4.start_num =20;
    Objs_4.end_num =25;
    Objs_4.thread_num =4;
    Objs_5.C = this;
    Objs_5.start_num =25;
    Objs_5.end_num =30;
    Objs_5.thread_num =5;
    Objs_6.C = this;
    Objs_6.start_num =30;
    Objs_6.end_num =35;
    Objs_6.thread_num =6;
    Objs_7.C = this;
    Objs_7.start_num =35;
    Objs_7.end_num=40;
    Objs_7.thread_num =7;



    count = 0;
    //std::cerr << "start" << std::endl;
    shared_lock = PTHREAD_MUTEX_INITIALIZER;
    //pthread_mutex_lock(&shared_lock);
    //std::cerr << "1" << std::endl;
    shared_cond = PTHREAD_COND_INITIALIZER;
    //std::cerr << "2" << std::endl;



    locks_0 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_0);
    locks_1 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_1);
    locks_2 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_2);
    locks_3 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_3);
    locks_4 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_4);
    locks_5 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_5);
    locks_6 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_6);
    locks_7 = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&locks_7);



    //std::cerr << "3" << std::endl;
    
    irets[0] = pthread_create(&threads_0, NULL, &rake_filter_process_helper, &Objs_0);
    irets[1] = pthread_create(&threads_1, NULL, &rake_filter_process_helper, &Objs_1);
    irets[2] = pthread_create(&threads_2, NULL, &rake_filter_process_helper, &Objs_2);
    irets[3] = pthread_create(&threads_3, NULL, &rake_filter_process_helper, &Objs_3);
    
    irets[4] = pthread_create(&threads_4, NULL, &rake_filter_process_helper, &Objs_4);
    irets[5] = pthread_create(&threads_5, NULL, &rake_filter_process_helper, &Objs_5);
    //irets[6] = pthread_create(&threads_6, NULL, &rake_filter_process_helper, &Objs_6);
    //irets[7] = pthread_create(&threads_7, NULL, &rake_filter_process_helper, &Objs_7);

    //std::cerr << "4" << std::endl;
}

void nearfield_demod_impl::setPulseLen(float pulse_len_in){
	d_pulse_len = pulse_len_in;
	pulse_max = d_pulse_len+d_pulse_len*d_pulse_len_accuracy/1e2;
	pulse_min = d_pulse_len-d_pulse_len*d_pulse_len_accuracy/1e2;
}

void nearfield_demod_impl::setPulseLenAccuracy(float pulse_len_accuracy_in){
	d_pulse_len_accuracy = pulse_len_accuracy_in;
	setPulseLen(d_pulse_len);
}

void nearfield_demod_impl::setBitrate(float bitrate_in){
	d_bitrate = bitrate_in;
	prf_max = 1/d_bitrate+1/d_bitrate*d_bitrate_accuracy/1e2;            // maximum pulse repetition frequency
	prf_min = 1/d_bitrate-1/d_bitrate*d_bitrate_accuracy/1e2;            // minimum pulse repetition freuqnecy
}

void nearfield_demod_impl::setBitrateAccuracy(float bitrate_accuracy_in){
	d_bitrate_accuracy = bitrate_accuracy_in;
	setBitrate(d_bitrate);
}

void nearfield_demod_impl::setPostPulseLenAccuracy(float post_pulse_len_accuracy_in){
	d_post_pulse_len_accuracy = post_pulse_len_accuracy_in;
}

void nearfield_demod_impl::setPostBitrateAccuracy(float post_bitrate_accuracy_in){
	d_post_bitrate_accuracy = post_bitrate_accuracy_in;
}

void nearfield_demod_impl::setSampleRate(float sample_rate){
	SDR_sample_rate = sample_rate;   // sample rate of software defined radio receiver
	sample_period = 1/SDR_sample_rate;
}

int nearfield_demod_impl::work(int noutput_items,
		gr_vector_const_void_star &input_items,
		gr_vector_void_star &output_items) {

	const float *in = (const float *) input_items[0];
	float *out = (float *) output_items[0];

	int nn = 0;
	float current = 0;

	//int valid_count = 0;
	//std::cout << "reset at beginning" << std::endl;
	while(nn < noutput_items){
		//Local non-persistent variables
		float rx_data;
		float transition;
        sample_counter++;

        	//do matched filter first
        /*
                float past = lastpulses.front();
                lastpulses.pop();
		float in_sqr = in[nn+NTAPS-1] * in[nn+NTAPS-1];
                energy = energy + in_sqr - past*past;
                lastpulses.push(in[nn]);
		float matched;
		if(sample_counter < 101) {
			matched = 0.01;
		} else { 
			matched = d_fir->filter(&in[nn])/sqrt(energy);
		}
        */

	    //std::cout << "get: " << in[nn] << std::endl;
    	if(in[nn] > max_current){
			//max_current = matched;
			max_current = in[nn];
		    //std::cout << "change into: " << max_current << std::endl;
	    	//avg_current += in[nn];
        }

		    //current = average_current;
	    sub_sample_counter++;
            //insert into the deque
		    //std::cout << "pushing: " << current << ", poping: " << matched_pulses[0].back() << std::endl;

    	if(sub_sample_counter == subsample_rate){
            //std::cout << "sample counter: " << sample_counter << std::endl;
            
            
           /* 
            if(max_current > 0.5) {
                std::cout << "current pulse: " << max_current << ", position: " << sample_counter << std::endl;
            }
            */
            
            
            
		/*
		for(int i=0; i<40; i++){
			scores[i] = 0;
		}
		*/	

    		sub_sample_counter = 0;
	    	current = max_current;
	    	//current = avg_current/subsample_rate;
		    max_current = 0;
		//std::cout << current << std::endl;
		    //std::cout << "finding header" << std::endl;
            for(int i = 0; i< num_rake_filter; i++){
	    		last[i] = matched_pulses[rake_offset[i]];
		    	all_pulse_energy[i] = all_pulse_energy[i] + current * current - last[i] * last[i];
                long_matched_out[i] = 0;

    		}
		    matched_pulses.pop_front();
	        matched_pulses.push_back(current);
            
            
            //std::cout << "enable 8 threads" << std::endl;
            pthread_mutex_unlock(&locks_0);
            pthread_mutex_unlock(&locks_1);
            pthread_mutex_unlock(&locks_2);
            pthread_mutex_unlock(&locks_3);
            pthread_mutex_unlock(&locks_4);
            pthread_mutex_unlock(&locks_5);
            pthread_mutex_unlock(&locks_6);
            pthread_mutex_unlock(&locks_7);
                    
            
            pthread_mutex_lock(&shared_lock);
            while(count < (num_rake_filter/5)){
                pthread_cond_wait(&shared_cond, &shared_lock);
            }
            count = 0;
            pthread_mutex_unlock(&shared_lock);
            

            //pthread_mutex_lock(&locks_0);
            //pthread_mutex_lock(&locks_1);
            //pthread_mutex_lock(&locks_2);
            //pthread_mutex_lock(&locks_3);
            //std::cout << "4 threads done" << std::endl;
/*
        for(int k = 0; k < 16; k++) {
            //std::cout << start_counter << std::endl;
            for(int i = 0; i < 40; i++){
            //for(int j = 0; j < matched_pulses[i].size(); j++){
                if(k == 0) {
                    
                    long_matched_out[i] = long_matched_out[i] - 
                                last[i] * last[i] + 
                                matched_pulses[rake_offset[i] + 2 * jitter] * matched_pulses[rake_offset[i] + 2 * jitter];
                }
                else if(k > 0 && k <= 15) {
            
                //std::cout << "k = " << k << ", insert: " << (1+unit_offset*i) * unit_time * sum_table[k] 
                ////		<< ", delete: " << (1+unit_offset*i) * unit_time * sum_table[k] + (2 * jitter - 1) + unit_offset * unit_time * sum_table[k] << std::endl;
                    long_matched_out[i] = long_matched_out[i] - 
                            matched_pulses[rake_offset[i]+int((1+unit_offset*i) * unit_time * sum_table[k] - 1)] * 
                            matched_pulses[rake_offset[i]+int((1+unit_offset*i) * unit_time * sum_table[k] - 1)] +
                            (matched_pulses[rake_offset[i]+int((1+unit_offset*i) * unit_time * sum_table[k] + (2 * jitter - 1) + 
                            unit_offset * unit_time * sum_table[k])]) * 
                            (matched_pulses[rake_offset[i]+int((1+unit_offset*i) * unit_time * sum_table[k] + (2 * jitter - 1) + 
                            unit_offset * unit_time * sum_table[k])]);
                }
            }                
        }
*/
        /*
   		if(noise_power == 0) {
			noise_power = current * current;
		} else { 
			noise_power = 0.00001 * current * current + (1-0.00001) * noise_power;
		}
        */
		
		//std::cout << "noise_power change into: " << noise_power << " at: " << sample_counter << std::endl;

        max_header_response = 0;
        if(matched_pulses.front() != 0 && start == 0 && sample_counter > (1/0.00001) * subsample_rate) {
        //if(matched_pulses.front() != 0 && start == 0) {
            //std::cout << "filling buffers done at " << sample_counter << std::endl;
            start = 1;
        }


    	//dump data out
	    //std::cout << long_matched_out[0] << std::endl;
			//noise_power = (all_pulse_energy[num_rake_filter-1] - long_matched_out[num_rake_filter-1])/(matched_pulses.size() - window_length[num_rake_filter-1]);

		//std::cout << noise_power << std::endl;



        for(int i = 0; i < num_rake_filter; i++){
            /*
            aggregated_header[i] = (long_matched_out[i] - 
                    (all_pulse_energy[i] * window_length[i])/matched_pulses[i].size())/
                    (all_pulse_energy[i] - long_matched_out[i]);
            */
            //aggregated_header[i] = (long_matched_out[i] - noise_power * window_length[i])/(noise_power * window_length[i]);
            //aggregated_header[i] = long_matched_out[i];
    //if(i == 0 && start == 1){
        //std::cout << "noise_power: " << noise_power << ", window_length: " << window_length[i] << std::endl;
        //std::cout << "background: " << noise_power * window_length[i] << std::endl;
        //std::cout << "sig: " << long_matched_out[i] << std::endl;
    //}	
            //aggregated_header[i] = (long_matched_out[i] - noise_power * window_length[i]);	
    

            if(long_matched_out[i] > max_header_response){
                max_header_response = long_matched_out[i];
                time_offset = i;
            }
                /*
            if(start == 1 && i < 40){
                //std::cout << aggregated_header[i] << ";";
                //std::cout << scores[i] << ";";
            }
        */
        }



			//std::cout << std::endl;
            //if(start == 0) {
				//std::cout << std::endl;
	    		//std::cout << 0 << ";" << sample_counter << std::endl;
            //}
           
	 
            //if(max_header_response > threshold) {
	            //std::cout << max_header_response << ", " << sample_counter << ", " << start << std::endl;
            //}

            if(sync == 1) {
                peak_distance++;
            }

            //if(sync == 0) {
            if(start == 1 && max_header_response > threshold) {
                if(sync == 0) {
                    if(max_header_response > last_max_response) {
                        last_max_response = max_header_response;
                        last_offset = time_offset;
                    } else {
                        time_offset = last_offset;
                        correct_offset = last_offset + 0.5;
                        /*
                       	std::cout << "find header, clock offset = " << time_offset << std::endl;
                        std::cout << "response = " << last_max_response << std::endl;
                        std::cout << "sample_counter = " << sample_counter << std::endl;
                        for(int x = 1; x<16; x++) {
                            std::cout << "matched_pulse[" << x << "] = " << pulse_queue[last_offset][x] << "; loc: " << 
                                pulse_loc_queue[last_offset][x] - pulse_loc_queue[last_offset][x-1] << std::endl;
                        }
                        std::cout << "header_data_distance: " << subsample_rate * header_data_distance * 
                                                (1 + unit_offset * correct_offset) << std::endl;
                        std::cout << "dis_one_to_zero: " << subsample_rate * dis_one_to_zero * 
                                                (1 + unit_offset * correct_offset) << std::endl;
                        std::cout << "dis_one_to_one: " << subsample_rate * dis_one_to_one *
                                                (1 + unit_offset * correct_offset) << std::endl;
                        std::cout << "dis_zero_to_zero: " << subsample_rate * dis_zero_to_zero * 
                                                (1 + unit_offset * correct_offset) << std::endl;
                        std::cout << "dis_zero_to_one: " << subsample_rate * dis_zero_to_one * 
                                                (1 + unit_offset * correct_offset) << std::endl;
                        */
                        sync = 1;
			            last_peak_response = last_max_response;
                        data_acquired = 0;
                        pos = 0;
                        target_zero_pos = header_data_distance * (1 + unit_offset * correct_offset) - data_window;
                        target_one_pos = (header_data_distance + dis_zero_to_zero * one_zero_rate/(1))*(1 + unit_offset * correct_offset)
                                            - data_window;
    	                data_pos_0[0] = target_zero_pos;
	                	data_pos_1[0] = target_one_pos;
                        peak_distance = 0;
                        //std::cout << "2.5" << std::endl;
                        /*
                        for(int i = 0; i < matched_pulses[time_offset].size(); i++) {
                            std::cout << matched_pulses[time_offset][i] << std::endl;
                        }
                        */
                    }
                } else if(sync == 1) {
                    if(max_header_response > last_max_response) {
                        last_max_response = max_header_response;
                        last_offset = time_offset;
                    } else {
                        if(peak_distance < (100000/subsample_rate) && last_max_response > last_peak_response + 0.2) {
                            //valid stronger peak!!!
                            time_offset = last_offset;
                            correct_offset = last_offset + 0.5;
                            /*
                            std::cout << "find new header, clock offset = " << time_offset << std::endl;
                            std::cout << "response = " << last_max_response << ", last peak = " << last_peak_response << std::endl;
                            std::cout << "sample_counter = " << sample_counter << std::endl;
                            for(int x = 1; x<16; x++) {
                                std::cout << "matched_pulse[" << x << "] = " << pulse_queue[last_offset][x] << "; loc: " << 
                                    pulse_loc_queue[last_offset][x] - pulse_loc_queue[last_offset][x-1] << std::endl;
                            }
                            std::cout << "header_data_distance: " << subsample_rate * header_data_distance * 
                                                    (1 + unit_offset * correct_offset) << std::endl;
                            std::cout << "dis_one_to_zero: " << subsample_rate * dis_one_to_zero * 
                                                    (1 + unit_offset * correct_offset) << std::endl;
                            std::cout << "dis_one_to_one: " << subsample_rate * dis_one_to_one *
                                                    (1 + unit_offset * correct_offset) << std::endl;
                            std::cout << "dis_zero_to_zero: " << subsample_rate * dis_zero_to_zero * 
                                                    (1 + unit_offset * correct_offset) << std::endl;
                            std::cout << "dis_zero_to_one: " << subsample_rate * dis_zero_to_one * 
                                                    (1 + unit_offset * correct_offset) << std::endl;
                            */
                            sync = 1;
			                last_peak_response = last_max_response;
                            data_acquired = 0;
                            pos = 0;
                            target_zero_pos = header_data_distance * (1 + unit_offset * correct_offset) - data_window;
                            target_one_pos = (header_data_distance + dis_zero_to_zero*one_zero_rate/(1))*(1 + unit_offset*correct_offset)
                                                - data_window;
    	                    data_pos_0[0] = target_zero_pos;
	                    	data_pos_1[0] = target_one_pos;
                            peak_distance = 0;
                            reset_data = 1;
                            //std::cout << "2.5" << std::endl;
                            /*
                            for(int i = 0; i < matched_pulses[time_offset].size(); i++) {
                                std::cout << matched_pulses[time_offset][i] << std::endl;
                            }
                            */          
                        } else {
                            reset_data = 0;
                        }
                    }
                } else {
                    std::cout << "unexist, ERROR!!!" << std::endl;
                }
            }
            

    
            //header identified, find the data
            if(sync == 1 && n < N) {
                if(reset_data == 1) {
    	            for(int i = 0; i < N; i++){
	            		data_energy_0[i] = 0;
	            		data_energy_1[i] = 0;
	            		data_pos_0[i] = 0;
	            		data_pos_1[i] = 0;
                        //std::cout << "reset pos to zero for new header" << std::endl;
	            	}
                    data_energy_out = 0;
		            n = 0;
		            demod_data.clear();
                    data_acquired = 0;
                    pos = 0;
                    target_zero_pos = header_data_distance * (1 + unit_offset * correct_offset) - data_window;
                    target_one_pos = (header_data_distance + dis_zero_to_zero * one_zero_rate/(1)) * (1 + unit_offset * correct_offset)
                                        - data_window;
	                data_pos_0[0] = target_zero_pos;
	            	data_pos_1[0] = target_one_pos;
                    reset_data = 0;
                }    
 
                //std::cout << current << std::endl;
		        //std::cout << "finding data" << std::endl;
                //std::cout << "pos: " << pos << "; n: " << n << "; data_acquired: " << data_acquired << std::endl;
                if(pos >= (target_zero_pos -
                    ((unit_offset * 0) * (data_distance * 15 + header_data_distance) + data_window * jitter)) && // window 7
                    pos <= (target_zero_pos +
                    ((unit_offset * 0) * (data_distance * 15 + header_data_distance) + data_window * jitter))){  //window 7.49

                    //data_energy_0[i] = current * current + data_energy_0[i];
                    if(data_energy_0[n] < current * current){
                        data_energy_0[n] = current * current;
                        data_pos_0[n] = pos;
                        //std::cout << "current pos zero" << n << " change to " << data_pos_0[n] << std::endl;
                    }
                   //std::cout << "pos: " << pos << " zero[" << n << "]: " << current << std::endl;
	            }
                if (pos>=int(target_one_pos -
                    ((unit_offset * 0) * (data_distance * 15 + header_data_distance) + data_window * jitter)) && //window 7.49
                    pos<=int(target_one_pos +
                    ((unit_offset * 0) * (data_distance * 15 + header_data_distance) + data_window * jitter))){ //window 7
                    if(data_energy_1[n] < current * current) {
                        data_energy_1[n] = current * current;
                        data_pos_1[n] = pos;
                        //std::cout << "current pos one" << n << " change to " << data_pos_1[n] << std::endl;
                    }
                    if (pos == int(target_one_pos + 
                                ((unit_offset * 0) * (data_distance * 15 + header_data_distance) + data_window * jitter))){ // window 7
                        data_acquired = 1;
                    }
                    //std::cout << "pos: " << pos << " one[" << n << "]: " << current << std::endl;
                }

                //std::cout << "current pos: " << data_pos_0[n] << "; " << data_pos_1[n] << std::endl;
                //compare and get data
                if((data_energy_0[n] >= data_energy_1[n]) && data_acquired){
                        demod_data.push_back(1);
                        n++;
                        target_one_pos = (dis_zero_to_one * (1 + unit_offset * correct_offset))  + data_pos_0[n - 1];
                        target_zero_pos = (dis_zero_to_zero * (1 + unit_offset * correct_offset)) + data_pos_0[n - 1];
                        //std::cout << "current: " << data_pos_0[n] << 
                        //        "; next_zero: " << target_zero_pos << "; next_one: " << target_one_pos << std::endl; 
                        data_acquired = 0;
                } else if((data_energy_0[n] <= data_energy_1[n]) && data_acquired) {
                        demod_data.push_back(0);
                        n++;
                        target_zero_pos = (dis_one_to_zero  * (1 + unit_offset * correct_offset)) + data_pos_1[n - 1];
                        target_one_pos = (dis_one_to_one  * (1 + unit_offset * correct_offset)) + data_pos_1[n - 1];
                        //std::cout << "current: " << data_pos_1[n] << 
                        //        "; next_zero: " << target_zero_pos << "; next_one: " << target_one_pos << std::endl; 
                        data_acquired = 0;
                }
                pos++;
                
            }


            //got all data
			if(n == N && peak_distance > (100000/subsample_rate)){                            // we've looked for all the data
				//Prepare outgoing packet for GATD
				//std::cout << "got all data" << std::endl;
				std::vector<uint8_t> demod_data_out;
				for(int ii=0; ii < d_gatd_id.size(); ii++)
					demod_data_out.push_back((uint8_t)d_gatd_id[ii]);
				for(int ii=0; ii < demod_data.size(); ii++)
					demod_data_out.push_back(demod_data[ii]);

				//Push message out with packet data
				pmt::pmt_t value = pmt::init_u8vector(demod_data_out.size(), (const uint8_t*)&demod_data_out[0]);
				pmt::pmt_t new_message = pmt::cons(pmt::PMT_NIL, value);
				message_port_pub(pmt::mp("frame_out"), new_message);
				float prf_length = roundf(mean(prf_vec)/sample_period);
				time_t current_time = time(0);
	            double seconds = difftime(current_time, last_time);
	            last_time = current_time;
				char* dt = std::ctime(&current_time);
					
				std::cout << "SENDING MESSAGE" << std::endl;
                std::cout << "offset: " << correct_offset << std::endl;
				std::cout << "@@@" << dt << ", " << seconds << " second." << " bitrate: " <<
                    (10000000/(subsample_rate * (1+correct_offset * unit_offset) * data_distance)) << std::endl;
				d_log_file << "SENDING MESSAGE" << std::endl;
			        d_log_file << "@@@" << dt << ", " << seconds << " second." << " bitrate: " << 
                        (10000000/(subsample_rate * (1+correct_offset * unit_offset) * data_distance)) << std::endl;
                    /*
				for(int ii=0; ii < demod_data.size(); ii++){
					std::cout << (int)(demod_data[ii]) << ", ";
					d_log_file << (int)(demod_data[ii]) << ", ";
				}*/
				for(int ii = demod_data.size() - 1; ii >= 0; ii--){
					std::cout << (int)(demod_data[ii]) << ", ";
					d_log_file << (int)(demod_data[ii]) << ", ";
				}
				d_log_file << std::endl;
				std::cout << std::endl;
				
			
				sync = 0;                        // start over looking for sync
				//prf_vec.clear();
				//pulse_vec.clear();
				demod_data.clear();
				//prf_count = 0;
				pulse_count = 0;
				n = 0;
				//sync_prf = 0;
				//sync_prf2 = 0;
				//sync_pulse = 0;
				valid_pulse = 0;
				valid_count = 0;
                last_max_response = 0;
                last_peak_response = 0;
                reset_data = 0;
                peak_distance = 0;
                last_offset = 0;
		        correct_offset = 0;
                data_acquired = 0;
                pos = 0;
                target_one_pos = 0;
                target_zero_pos = 0;
				for(int i = 0; i < N; i++){
	            	data_energy_0[i] = 0; 
	            	data_energy_1[i] = 0; 
	            	data_pos_0[i] = 0; 
	            	data_pos_1[i] = 0; 
                    //std::cout << "reset pos to zero for new header" << std::endl;
	          	}
                data_energy_out = 0;
                time_offset = 0;
                last_offset = 0;
				//std::cout << "found all bits, clear valid_counter: " << valid_count << std::endl;
			}
		}
	
		//Increment the data index pointer
		nn++;
	}
	// Tell runtime system how many output items we produced.
	return noutput_items;
}

} /* namespace nearfield */
} /* namespace gr */






