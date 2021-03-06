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

#ifndef INCLUDED_NEARFIELD_NEARFIELD_DEMOD_IMPL_H
#define INCLUDED_NEARFIELD_NEARFIELD_DEMOD_IMPL_H

#include <nearfield/nearfield_demod.h>
#include <fstream>
#include <queue>


namespace gr {
  namespace nearfield {
    class nearfield_demod_impl;
    struct object {
        nearfield_demod_impl* C;
        int start_num;
        int end_num;
        int thread_num;
    };



    class nearfield_demod_impl : public nearfield_demod
    {
     private:
	float N;
	float threshold;
	float SDR_sample_rate;
	float pulse_max;
	float pulse_min;
	float prf_max;
	float prf_min;
	float sample_period;
	float last_pulse_length;
	float last_pulse_distance;
	float last_data;
	float pulse_count;
	float prf_count;
	float sync;
	float valid_count;
	float error;
	float fuzz;
	float last_prf_count;
	float header;
	float sync_prf;
	float sync_prf2;
	float sync_pulse;
	float prf_win_cnt;
	float valid_pulse;
	float n;
	float d_pulse_len;
	float d_pulse_len_accuracy;
	float d_post_pulse_len_accuracy;
	float d_bitrate;
	float d_bitrate_accuracy;
	float d_post_bitrate_accuracy;
	float last_prf;
	float last_pulse;
	float max_sample;
    int last_offset;    
    int correct_offset;    
    long int sample_counter;
	std::queue<float> lastpulses;
	float threshold_sync;
    float unit_time;
    int time_offset;
    float window_size;
	//std::deque<float> matched_pulses;
	std::deque<float> data_queue[100];
	float energy;
	float all_pulse_energy[40];
	int sample_ctr;
	std::vector<float> pulse_vec;
	std::vector<float> prf_vec;
	std::vector<uint8_t> demod_data;
	float sav_pulse;
	std::ofstream d_log_file;
	std::string d_gatd_id;
	time_t last_time;

	//FIR filter necessities
	gr::filter::kernel::fir_filter_fff *d_fir;

    float aggregated_header[40];
	int scores[40];
	int sub_sample_counter;
    int distance_table[16];
    int seed_table[16];
	int window_length[40];
    int sum_table[16];
    std::deque<float> matched_pulses;
	float max_header_response;
    float last_max_response;
	int pos;	
	int jitter;
	int start;
	float last[40];
    int subsample_rate;
	float last_peak_response;
	float noise_power;
	float max_current;
	float avg_current;
	float data_energy_0[100];
	float data_energy_1[100];
    float reset_data;
	int delay;
    float peak_distance;
    float long_matched_out[40];
	float data_energy_out;
    int num_rake_filter;
    //std::deque<float> matched_pulses;
    float unit_offset;
    float max_offset;
    int rake_offset[40];
    int process_counter;

    pthread_t threads_0;
    pthread_t threads_1;
    pthread_t threads_2;
    pthread_t threads_3;

    pthread_t threads_4;
    pthread_t threads_5;
    pthread_t threads_6;
    pthread_t threads_7;

    pthread_mutex_t locks_0;
    pthread_mutex_t locks_1;
    pthread_mutex_t locks_2;
    pthread_mutex_t locks_3;

    pthread_mutex_t locks_4;
    pthread_mutex_t locks_5;
    pthread_mutex_t locks_6;
    pthread_mutex_t locks_7;

    pthread_mutex_t shared_lock;
    pthread_cond_t shared_cond;
    int count;
    int irets[8];
    object Objs_0;
    object Objs_1;
    object Objs_2;
    object Objs_3;
    object Objs_4;
    object Objs_5;
    object Objs_6;
    object Objs_7;


     public:
      nearfield_demod_impl(float sample_rate, float bitrate, float bitrate_accuracy, float post_bitrate_accuracy, float pulse_len, float pulse_len_accuracy, float post_pulse_len_accuracy, int packet_len, int header_len, const std::string gatd_id);
      ~nearfield_demod_impl();

      // Where all the action really happens
      int work(int noutput_items,
	       gr_vector_const_void_star &input_items,
	       gr_vector_void_star &output_items);

      float getLastObservedBitrate();
      float getLastObservedPulseLen();
      float getThreshold();
      void setPulseLen(float pulse_len_in);
      void setPulseLenAccuracy(float pulse_len_accuracy_in);
      void setPostPulseLenAccuracy(float post_pulse_len_accuracy_in);
      void setBitrate(float bitrate_in);
      void setBitrateAccuracy(float bitrate_accuracy_in);
      void setPostBitrateAccuracy(float post_bitrate_accuracy_in);
      void setSampleRate(float sample_rate_in);
      void setPacketLen(int packet_len_in);
      void setHeaderLen(int header_len_in);
      void rake_filter_process(int start_num, int end_num, int thread_num);


    };

    void* rake_filter_process_helper(void * obj);
  } // namespace nearfield
} // namespace gr

#endif /* INCLUDED_NEARFIELD_NEARFIELD_DEMOD_IMPL_H */

