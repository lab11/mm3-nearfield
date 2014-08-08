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

namespace gr {
namespace nearfield {

float mean(std::vector<float> in_vec){
	float sum = std::accumulate(in_vec.begin(), in_vec.end(), 0.0);
	float m = sum / in_vec.size();
	return m;
}

nearfield_demod::sptr nearfield_demod::make(float sample_rate, float bitrate, float bitrate_accuracy, float post_bitrate_accuracy, float pulse_len, float pulse_len_accuracy, float post_pulse_len_accuracy, int packet_len, int header_len) {
		return gnuradio::get_initial_sptr
			(new nearfield_demod_impl(sample_rate, bitrate, bitrate_accuracy, post_bitrate_accuracy, pulse_len, pulse_len_accuracy, post_pulse_len_accuracy, packet_len, header_len));
	}

/*
 * The private constructor
 */
nearfield_demod_impl::nearfield_demod_impl(float sample_rate, float bitrate, float bitrate_accuracy, float post_bitrate_accuracy, float pulse_len, float pulse_len_accuracy, float post_pulse_len_accuracy, int packet_len, int header_len)
	: gr::sync_block("nearfield_demod",
			gr::io_signature::make(1, 1, sizeof(float)),
			gr::io_signature::make(0, 1, sizeof(float))),
	  d_log_file("nearfield_log.txt") {

	// variables
	threshold = 0.1;            // threshold set after observing data
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
	
	// Treating this like it's streaming data so I won't use some of the 
	// efficient Matlab functions that would speed this up.
	last_data = 0;      // previous received data point
	pulse_count = 0;    // length of current pulse (1s)
	prf_count = 0;      // length of current prf (0s)
	sync = 0;           // current consecutive sync pulses
	error = 0;
	fuzz = 0;
	last_prf_count = 0;
	valid_count = 0;
	sync_prf = 0;       // prf counter in sync part of loop
	sync_prf2 = 0;      // prf counter in sync part of loop (prf_window)
	sync_pulse = 0;     // pulse counter in sync part of loop
	prf_win_cnt = 0;    // counter for window of pulses between prf_min and prf_max in sync routine
	valid_pulse = 0;    // flag for pulse detection
	n = 0;              // counter for N
	max_sample = 0;

	pulse_vec.clear();
	prf_vec.clear();
	demod_data.clear();

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

	//int valid_count = 0;
	//std::cout << "reset at beginning" << std::endl;
	while(nn < noutput_items){
		//Local non-persistent variables
		float rx_data;
		float transition;

		//Update sample counter for AGC logic
		sample_ctr++;

		if(sample_ctr > 1e6){
		
			if(max_sample < threshold/1.2){
				threshold /= 1.1;
		                //threshold = 0.999*threshold + 0.001*in[nn];
				//std::cout << "threshold decrease: " << threshold << std::endl;
			}
		

		        //threshold = 0.999*threshold + 0.001*in[nn];
			sample_ctr = 0;
			max_sample = 0;

		}

		

		//threshold = 0.999*threshold + 0.001*in[nn];
	
		// STEP 1 --------------------------------------------------------------
/*
		if(in[nn] > 1) {
			std::cout << "start!" << std::endl;
			std::cout << "volt: " << in[nn] << ", threshold: " << threshold << std::endl;
		}
*/

		//std::cout << "volt: " << in[nn] << ", threshold: " << threshold << std::endl;
		if(in[nn] > threshold){
			max_sample = (in[nn] > max_sample) ? in[nn] : max_sample;
			rx_data = 1;
			sample_ctr = 0;
			//threshold = 0.999*threshold + 0.001*in[nn];             // increment threshold since we're probably picking up noise
		}else
			rx_data = 0;
	
		// STEP 2 --------------------------------------------------------------
		transition = rx_data - last_data;     // see if we are at an edge of a pulse
		// NOT SYNCHRONIZED YET
		if(sync < header) {

			if(transition > 0) {           // rising pulse
				
				//std::cout << "transition > 0, volt = " << in[nn] << std::endl; 
				//std::cout << "prf_count: " << prf_count << ", prf_val: " << prf_count*sample_period << std::endl;
				//std::cout << "prf_min: " << prf_min << "prf_max: " << prf_max << std::endl;
				threshold = 0.95*threshold + 0.05*in[nn];             // increment threshold since we're probably picking up noise
				//std::cout << "increase threshold to: " << threshold << std::endl; 
				// check for prf
				if(sync >= 1){
					float prf_val = prf_count*sample_period;
					if(prf_val > prf_min && prf_val < prf_max){
						prf_vec.push_back(prf_val);
						//std::cout << "sync valid " << sync << std::endl;
					} else if (fuzz == 0 &&  prf_val <= prf_min){
						fuzz += 1;
						threshold = 0.95*threshold + 0.05*in[nn];             // increment threshold since we're probably picking up noise
						//std::cout << "tolerant 1 fuzz " << std::endl;
						//std::cout << "threshold increase: " << threshold << std::endl;

					} else {
						sync = 0;                     // reset the sync counter
						valid_count = 0;
						prf_vec.clear();                 // reset prf vector
						pulse_vec.clear();               // reset pulse vector
						threshold *= 1.1;
						threshold = 0.95*threshold + 0.05*in[nn];             // increment threshold since we're probably picking up noise
						//std::cout << "sync abandoned " << std::endl;
						//std::cout << "threshold increase: " << threshold << std::endl;
						//std::cout << "clear, valid_coun: " << valid_count << ", sync: " << sync << std::endl; 
					}	
					
					//Update last_prf for every header bit
					last_prf = prf_val;
				}
				pulse_count = 1;     
				last_prf_count = prf_count;       
				prf_count = 0;  
			} else if(transition < 0) {       // falling edge of pulse
				//std::cout << "transition < 0, volt = " << in[nn] << std::endl;
				/* 
				if ( threshold > 0.8 && sync > 0) {
					pulse_count += 3;
				} else if (threshold > 0.7 && threshold <= 0.8 && sync > 0) {
					pulse_count += 2;
				} else if (threshold > 0.6 && threshold <= 0.7 && sync > 0) {
					pulse_count += 1;
				} else {
					pulse_count = pulse_count;
				}
				*/

				// determine if pulse is valid or not
				float pulse_val = pulse_count*sample_period;
		                float pulse_val_min = (pulse_count - 1) * sample_period; 
		                float pulse_val_max = (pulse_count + 1) * sample_period; 
				//Update last_pulse for every header bit
				if(pulse_count > 1)
					last_pulse = pulse_val;
				//std::cout << "pulse_count: " << pulse_count << ", pulse_val: " << pulse_val << std::endl;
				//std::cout << "pulse_min: " << pulse_min << ", pulse_max: " << pulse_max << std::endl;
				
				if((pulse_val_min > pulse_min && pulse_val_min < pulse_max) || 
						(pulse_val_max > pulse_min && pulse_val_max < pulse_max)) {
					// increment sync counter and pulse vector
					sync = sync+1;              // valid pulse
					error = 0;		//leave space for an error
					fuzz = 0;
					pulse_vec.push_back(pulse_val);
					pulse_count = 0;
					if (pulse_val > pulse_min && pulse_val < pulse_max) {
						valid_count = valid_count + 1;
						//std::cout << "accept, sync count: " << sync << ", valid_count: " << valid_count << std::endl;
					} else {
						//std::cout << "sudo accept, sync count: " << sync << ", valid_count: " << valid_count << std::endl;
					}
				} else {				
			 		if (error == 0 && sync >= 1 ) {	
						//std::cout << "tolerant 1 error!!!" << std::endl;
						error = + 1;
						prf_count = last_prf_count + pulse_count;
						pulse_count = 0;
					} else {
						pulse_count = 0;
						sync = 0;
						error = 0;
						valid_count = 0;
						//std::cout << "clear, valid_counter: " << valid_count << ", sync: " << sync << std::endl;
					}
				}
			} else {    // no transition
				if(rx_data == 1)
					pulse_count = pulse_count+1;
				prf_count = prf_count+1;
			}
			last_data = rx_data;
		}
		// SYNCHRONIZED
		if(sync == header && valid_count < header - 1) {
			valid_count = 0;
			sync = 0;
		}
		if(sync == header && valid_count >= (header - 1)) {
			error = 0;
			fuzz = 1;

			// find average pulse width and prf
			float pulse_length = roundf(mean(pulse_vec)/sample_period);
			float prf_length = roundf(mean(prf_vec)/sample_period);
			// using this, demodulate the next N bits
			// right after last sync pulse is a prf window
			float pulse_length_min = pulse_length-pulse_length*d_post_pulse_len_accuracy/1e2 - 1;
			float pulse_length_max = pulse_length+pulse_length*d_post_pulse_len_accuracy/1e2 + 1;
			float prf_length_min = roundf(prf_length - prf_length*d_post_bitrate_accuracy/1e2);
			float prf_length_max = roundf(prf_length + prf_length*d_post_bitrate_accuracy/1e2);
			float prf_window = prf_length_max - prf_length_min;
			// advance to the edge of the prf_min window and look for a pulse
			// anywhere in the prf_max-prf_min frame.
			// then return to prf_length and repeat.
			if(sync_prf < prf_length_min) {         // don't care until we get to the point where a pulse might come
				sync_prf = sync_prf + 1;
			} else {                                 // start looking for a pulse (we are in the prf window)
				sync_prf2 = sync_prf2 + 1;       // do this for prf window duration
			}
			if(rx_data == 1) {
				sync_pulse = sync_pulse + 1; // might be a pulse
			} else {
				/*
				if(threshold > 0.7) {
					sync_pulse += 2;
				} else if (threshold > 0.6 && threshold <= 0.7) {
					sync_pulse += 1;
				} else {
					sync_pulse = sync_pulse;
				}
				*/
				
				if(sync_pulse >= pulse_length_min && sync_pulse <= pulse_length_max) {     // valid pulse
					valid_pulse = 1;
					sav_pulse = sync_pulse;
				}
				//std::cout << "pulse too short, recognized as zero" << std::endl;
				sync_pulse = 0;              // reset
			}
			if(valid_pulse == 1) {              // we found a pulse
				demod_data.push_back(1);
				n = n + 1;
				sync_prf = sav_pulse;        // prf is defined as rising edge of pulse to the next rising edge of pulse
				sync_prf2 = 0;
				valid_pulse = 0;             // reset
			} else if(sync_prf2 == prf_window) {   // ran through whole prf window w/o finding pulse
				demod_data.push_back(0);
				n = n + 1;
				sync_prf = prf_length_max - prf_length;      // don't set to 0...reset to middle of window for timing
				sync_prf2 = 0;
			}
		}
		if(n == N){                            // we've looked for all the data
			//Push message out with packet data
			pmt::pmt_t new_message_dict = pmt::make_dict();
			pmt::pmt_t key = pmt::from_long((long)(0));
			pmt::pmt_t value = pmt::init_u8vector(demod_data.size(), (const uint8_t*)&demod_data[0]);
			new_message_dict = pmt::dict_add(new_message_dict, key, value);
			pmt::pmt_t new_message = pmt::cons(new_message_dict, pmt::PMT_NIL);
			message_port_pub(pmt::mp("frame_out"), new_message);

			std::cout << "SENDING MESSAGE" << std::endl;
			for(int ii=0; ii < demod_data.size(); ii++){
				std::cout << (int)(demod_data[ii]) << ", ";
				d_log_file << (int)(demod_data[ii]);
			}
			d_log_file << std::endl;
			std::cout << std::endl;

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
			//std::cout << "found all bits, clear valid_counter: " << valid_count << std::endl;
		}
	
		//Increment the data index pointer
		nn++;
	}

	// Tell runtime system how many output items we produced.
	return noutput_items;
}

} /* namespace nearfield */
} /* namespace gr */

