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
    long int sample_counter;
	std::queue<float> lastpulses;
	float threshold_sync;
    float unit_time;
    int time_offset;
    float window_size;
	std::deque<float> matched_pulses;
	std::deque<float> data_queue[100];
	std::deque<float> lastsamples;
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
    float long_matched_out[40];
    float aggregated_header[40];
	int sub_sample_counter;
    int distance_table[16];
	int window_length[40];
    int sum_table[16];
	float max_header_response;
    float last_max_response;
	int pos;	
	int jitter;
	int start;
	float last[40];
    int subsample_rate;
	float noise_power;
	float max_current;
	float avg_current;
	float data_energy[100];
	float data_energy_out;
    int num_rake_filter;
    float unit_offset;
    float max_offset;
    int rake_offset[40];

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
      void* rake_filter_process(void* start_num);
    };

  } // namespace nearfield
} // namespace gr

#endif /* INCLUDED_NEARFIELD_NEARFIELD_DEMOD_IMPL_H */

