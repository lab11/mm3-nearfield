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

#include <gnuradio/io_signature.h>
#include "nearfield_demod_impl.h"
#include <numeric>
#include <iterator>
#include <queue>
#include <ctime>
#include <iostream>
namespace gr {
namespace nearfield {

float mean(std::vector<float> in_vec){
	float sum = std::accumulate(in_vec.begin(), in_vec.end(), 0.0);
	float m = sum / in_vec.size();
	return m;
}

nearfield_demod::sptr nearfield_demod::make(float sample_rate, float bitrate, float bitrate_accuracy, float post_bitrate_accuracy, float pulse_len, float pulse_len_accuracy, float post_pulse_len_accuracy, int packet_len, int header_len, const std::string gatd_id) {
		return gnuradio::get_initial_sptr
			(new nearfield_demod_impl(sample_rate, bitrate, bitrate_accuracy, post_bitrate_accuracy, pulse_len, pulse_len_accuracy, post_pulse_len_accuracy, packet_len, header_len, gatd_id));
	}

/*
 * The private constructor
 */
nearfield_demod_impl::nearfield_demod_impl(float sample_rate, float bitrate, float bitrate_accuracy, float post_bitrate_accuracy, float pulse_len, float pulse_len_accuracy, float post_pulse_len_accuracy, int packet_len, int header_len, const std::string gatd_id)
	: gr::sync_block("nearfield_demod",
			gr::io_signature::make(1, 1, sizeof(float)),
			gr::io_signature::make(0, 1, sizeof(float))),
	  d_log_file("nearfield_log.txt"), d_gatd_id(gatd_id) {

	// variables
	threshold = 0.7;            // threshold set after observing data
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
    subsample_rate = 5;
	
	threshold_sync = 0.7;
	max_header_response = 0;
    unit_time = 55000/(16 * subsample_rate);
    time_offset = 0;
	pos = 0;
    jitter = int(unit_time * 16 * 0.03);
    //std::cout << jitter << std::endl;
	//jitter = 1;
    num_rake_filter = 40;
    unit_offset = 0.0025;
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
        for (int k = 0; k < 100; k++ ) {
        	lastpulses.push(0);
	}
        for (int k = 0; k < 7; k++ ) {
        	lastsamples.push_front(0);
	}
	for (int k = 0; k < num_rake_filter; k++){

		//std::cout << "size of buffer[" << k << "] = " << (345 * unit_time * (1+unit_offset*k) + 2 * jitter + 345 * unit_offset * unit_time) << std::endl; 
 		if(k == num_rake_filter - 1){
            for(int j = 0; j < (345 * unit_time * (1+unit_offset*k) + 2 * jitter + 345 * unit_offset * unit_time); j++){
			    matched_pulses.push_front(0);
            }
        }
        rake_offset[k] = (345 * unit_time * (1+unit_offset*39) + 2 * jitter + 345 * unit_offset * unit_time) -
                            (345 * unit_time * (1+unit_offset*k) + 2 * jitter + 345 * unit_offset * unit_time);


		all_pulse_energy[k] = 0;

		//std::cout << "size of buffer[" << k << "] = " << matched_pulses[k].size() << std::endl; 
	}
	energy = 0;
        for(int i = 0; i < num_rake_filter; i++){
		long_matched_out[i] = 0;     
        aggregated_header[i] = 0;
	}
        for(int i = 0; i < 100; i++){
		data_energy[i] = 0; 
        /*
        for(int j = 0; j < unit_time * 16 * 0.1; j++){
            data_queue[i].push_front(0);
        }
        */
	}

    data_energy_out = 0;
        //init tables
        distance_table[0] = 0;
        distance_table[1] = 23;
        distance_table[2] = 16;
        distance_table[3] = 20;
        distance_table[4] = 25;
        distance_table[5] = 29;
        distance_table[6] = 27;
        distance_table[7] = 24;
        distance_table[8] = 22;
        distance_table[9] = 28;
        distance_table[10] = 17;
        distance_table[11] = 21;
        distance_table[12] = 30;
        distance_table[13] = 18;
        distance_table[14] = 26;
        distance_table[15] = 19;
       
        last_offset = 0;
	start = 0;
    last_max_response = 0;
        sum_table[0] = 0;
        for(int i = 1; i < 16; i++){
            sum_table[i] = sum_table[i-1] + distance_table[i];
		//std::cout << sum_table[i] << std::endl;
        }
    //std::cout << N << std::endl;
	for(int i = 0; i< num_rake_filter; i++){
		for(int j = 0; j < 16; j++) {
			if(j == 0) {
				window_length[i] = 2*jitter;
			} else {
				window_length[i] += int((1+unit_offset*i) * unit_time * sum_table[j] + (2 * jitter - 1) + 
						unit_offset * unit_time * sum_table[j]) - int((1+unit_offset*i) * unit_time * sum_table[j]);
			}
		}
	}

	message_port_register_out(pmt::mp("frame_out"));
}

/*
 * Our virtual destructor.
 */
nearfield_demod_impl::~nearfield_demod_impl() {
	d_log_file.close();
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
		//float past = lastpulses.front();
		//lastpulses.pop();
		//energy = energy + in[nn] * in[nn] - past * past;
		//lastpulses.push(in[nn]);
		/*
		for(int m=0; m < 7; m++){
			lastsamples[m] = lastsamples[m+1];
		}
		*/
		//lastsamples.pop_front();
		//lastsamples.push_back(in[nn]);
		/*
		for (std::deque<float>::iterator it = lastsamples.begin(); it!=lastsamples.end(); ++it)
    				std::cout << ' ' << *it;
		*/
		//std::cout << "in: " << in[nn-filt_count] << std::endl;	
		//float energy_template = 4.5211;
		//float energy_template = 0.00017858 + 0.00096827 + 0.011 + 0.0827 + 0.64595 + 2.27846 + 2.019532 + 0.12233;
		//float energy_template = 5.1611885;
		//current = (in[nn] * 0.5914 + in[nn-1] * 1.1921 + in[nn-2] * 1.2286 + 
		//in[nn-3] * 0.8965 + in[nn-4] * 0.5363 + in[nn-5] * 0.324 + 
		//in[nn-6] * 0.1764 + in[nn-7] * 0.1156)/sqrt(energy*4.5211);
		//std::cout << "energy: " << energy << std::endl;
		//std::cout << "current: " << current << std::endl;	
		//current = (in[nn] * 0.1156 + in[nn-1] * 0.1764 + in[nn-2] * 0.324 + 
		//	in[nn-3] * 0.5363 + in[nn-4] * 0.8965 + in[nn-5] * 1.2286 + 
		//	in[nn-6] * 1.1921 + in[nn-7] * 0.5914)/sqrt(energy*4.5211);
		
		//current = (in[nn] * 0.01336336 + in[nn-1] * 0.03111696 + in[nn-2] * 0.104976 + 
		//	in[nn-3] * 0.28761769 + in[nn-4] * 0.80371225 + in[nn-5] * 1.50945796 + 
		//	in[nn-6] * 1.42110241 + in[nn-7] * 0.34975num_rake_filter-16)/sqrt(energy * energy_template);
		
		//current = (lastsamples[0] * 0.01336336 + lastsamples[1] * 0.03111696 + lastsamples[2] * 0.104976 + 
		//	lastsamples[3] * 0.28761769 + lastsamples[4] * 0.80371225 + lastsamples[5] * 1.50945796 + 
		//	lastsamples[6] * 1.42110241 + lastsamples[7] * 0.34975num_rake_filter-16)/sqrt(energy * energy_template);
	    
        /*
		current = (lastsamples[0] * 0.1156 + lastsamples[1] * 0.1764 + lastsamples[2] * 0.324 + 
			lastsamples[3] * 0.5363 + lastsamples[4] * 0.8965 + lastsamples[5] * 1.2286 + 
			lastsamples[6] * 1.1921 + lastsamples[7] * 0.5914)/sqrt(energy*4.5211);
		*/

		sub_sample_counter++;
    	if(in[nn] > max_current){
			max_current = in[nn];
	    	avg_current += in[nn];
        }

		    //current = average_current;
            //insert into the deque
		    //std::cout << "pushing: " << current << ", poping: " << matched_pulses[0].back() << std::endl;
    	if(sub_sample_counter == subsample_rate){
    		sub_sample_counter = 0;
	    	current = max_current;
		    max_current = 0;
    		avg_current = 0;

		    //std::cout << "finding header" << std::endl;
		    float last[num_rake_filter];
			for(int i = 0; i < num_rake_filter; i++){
		    	last[i] = 0;
		    }

            for(int i = 0; i< num_rake_filter; i++){
	    		last[i] = matched_pulses[rake_offset[i]];
		    	all_pulse_energy[i] = all_pulse_energy[i] + current * current - last[i] * last[i];

    		}
		    matched_pulses.pop_front();
	        matched_pulses.push_back(current);

            for(int k = 0; k < 16; k++) {
                for(int i = 0; i < num_rake_filter; i++){
                //for(int j = 0; j < matched_pulses[i].size(); j++){
		        	if(k == 0) {
				    	/*
					    if(i == 0) {
					    //std::cout << "k = " << k << ", insert: 0" << ", delete: " << 2 * jitter + 1 << std::endl;
    					std::cout << "k = " << k << ", insert: [0]: " << matched_pulses[i][0] * matched_pulses[i][0] 
	    					<< ", delete: [" << 2 * jitter + 1 << "]: " <<  
							matched_pulses[i][2 * jitter + 1] * matched_pulses[i][2 * jitter + 1] << std::endl;
		    			}
			    		*/
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
		    //std::cout << std::endl;
		    
    		if(matched_pulses.front() != 0 && start == 0) {
    			//std::cout << "filling buffers done" << std::endl;
		    	start = 1;
			}


    		//dump data out
	    	//std::cout << "matched: " << long_matched_out[0] << "energy: " << sqrt(all_pulse_energy[0]) << std::endl;
			noise_power = (all_pulse_energy[num_rake_filter-1] - long_matched_out[num_rake_filter-1])/(matched_pulses.size() - window_length[num_rake_filter-1]);
            max_header_response = 0;


            for(int i = 0; i < num_rake_filter; i++){
				/*
                aggregated_header[i] = (long_matched_out[i] - 
						(all_pulse_energy[i] * window_length[i])/matched_pulses[i].size())/
	    				(all_pulse_energy[i] - long_matched_out[i]);
	    		*/
                aggregated_header[i] = (long_matched_out[i] - noise_power * window_length[i])/(noise_power * window_length[i]);	

                if(aggregated_header[i] > max_header_response){
			        max_header_response = aggregated_header[i];
		            time_offset = i;
                }
                /*
	    		if(start == 1 && i < 7){
		    		//std::cout << aggregated_header[i] << ";";
		        }
                */
    		}
            if(start == 1) {
				//std::cout << std::endl;
	    		//std::cout << max_header_response << std::endl;

            }
            /*
            if(max_header_response > threshold) {
	            std::cout << max_header_response << ", " << sample_counter << ", " << start << std::endl;
            }
            */


            if(sync == 0) {
                if(start == 1 && max_header_response > threshold) {
                    if(max_header_response > last_max_response) {
                        last_max_response = max_header_response;
                        last_offset = time_offset;
                    } else {
                        //time_offset = last_offset;
                        std::cout << "find header, clock offset = " << time_offset << std::endl;
                        std::cout << "sample_counter = " << sample_counter << std::endl;
                        sync = 1;
                        pos = 0;
                        //std::cout << "2.5" << std::endl;
                        /*
                        for(int i = 0; i < matched_pulses[time_offset].size(); i++) {
                            std::cout << matched_pulses[time_offset][i] << std::endl;
                        }
                        */
                    }
                }
            }
    
            //header identified, find the data
            if(sync == 1) {
                
                //std::cout << current << std::endl;
		        //std::cout << "finding data" << std::endl;
	            if(pos <= ( ((1+unit_offset*last_offset) * unit_time * 16.4 * (N)) + ((unit_offset*(last_offset))*unit_time*16.4*16+2*jitter) )){
	            	for(int i = 0; i < N; i++){
                        //std::cout << "last_offset: " << last_offset << ", unit_offset: " << unit_offset << std::endl;
                        //std::cout << "offset: " << unit_offset * last_offset << std::endl;
                        //std::cout << "center: " << int(((1+unit_offset*last_offset)*unit_time*16*(i+1))) << std::endl;
                        //std::cout << "span: " << int(((unit_offset*last_offset)*unit_time*16*16)+jitter) << std::endl;
                        //std::cout << "down: " << int(((1+unit_offset*last_offset)*unit_time*16*(i+1))-((unit_offset*last_offset)*unit_time*16*16)-jitter) << std::endl;
                        //std::cout << "up: " << int(((1+unit_offset*last_offset)*unit_time*16*(i+1))+((unit_offset*last_offset)*unit_time*16*16)+jitter) << std::endl;
	            		if(pos>=int(((1+unit_offset*last_offset)*unit_time*16.4*(i+1))-((unit_offset*(last_offset))*unit_time*16.4*16)-2*jitter)&&
	            		pos<=int(((1+unit_offset*last_offset)*unit_time*16.4*(i+1))+((unit_offset*(last_offset))*unit_time*16.4*16)+2*jitter)){
                            if(pos==1+int(((1+unit_offset*last_offset)*unit_time*16.4*(i+1))-((unit_offset*(last_offset))*unit_time*16.4*16)-2*jitter)){
                                //std::cout << "2" << std::endl;
                            }
	            			data_energy[i] = current * current + data_energy[i];
                            //data_queue[i].pop_front();
                            //data_queue[i].push_back(current);
                            //std::cout << current << std::endl;
	            		} else{
                            data_energy_out = current * current + data_energy_out;
                        }
	            	}
	            	pos++;
	            }
	            else {
                    //std::cout << "looking at data" << std::endl;
	            	for(int i = 0; i < N; i++) {
                        //float sig_power = data_energy[i] / (data_energy_out/N); 
                        float sig_power = (data_energy[i] - (noise_power * (((unit_offset*time_offset)*unit_time*16.4*16+2*jitter)*2)))/
	            			                (noise_power * (((unit_offset*time_offset)*unit_time*16.4*16+2*jitter)*2));
                        //std::cout << (data_energy[i] - (noise_power * (1+unit_offset*time_offset) * unit_time * 16 * 0.1))/
	            		//	(noise_power * (1+unit_offset*time_offset) * unit_time * 16 * 0.1) << std::endl;
                        std::cout << sig_power << std::endl;
	            		if(sig_power > threshold_sync){
								demod_data.push_back(1);
								n++;
	            		} else {
								demod_data.push_back(0);
								n++;
	            		}
	            	}
                    sync = 0;
                    //std::cout << "n: " << n << ", N: " << N << std::endl;
                }
            }


            //got all data
			if(n == N){                            // we've looked for all the data
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
				/*
				std::cout << "SENDING MESSAGE" << std::endl;
				std::cout << "@@@" << dt << ", " << seconds << " second." << " bitrate: " << (1/last_prf) << std::endl;
				d_log_file << "SENDING MESSAGE" << std::endl;
			        d_log_file << "@@@" << dt << ", " << seconds << " second." << " bitrate: " << (1/last_prf) << std::endl;
				for(int ii=0; ii < demod_data.size(); ii++){
					std::cout << (int)(demod_data[ii]) << ", ";
					d_log_file << (int)(demod_data[ii]) << ", ";
				}
				d_log_file << std::endl;
				std::cout << std::endl;
				*/
				sync = 0;                        // start over looking for sync
				prf_vec.clear();
				pulse_vec.clear();
				demod_data.clear();
				prf_count = 0;
				pulse_count = 0;
				n = 0;
				sync_prf = 0;
				sync_prf2 = 0;
				sync_pulse = 0;
				valid_pulse = 0;
				valid_count = 0;
                last_max_response = 0;
                last_offset = 0;
				for(int i = 0; i < N; i++){
	            	data_energy[i] = 0; 
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






