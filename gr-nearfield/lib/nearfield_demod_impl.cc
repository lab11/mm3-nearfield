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

namespace gr {
namespace nearfield {

float mean(std::vector<float> in_vec){
	float sum = std::accumulate(std::begin(in_vec), std::end(in_vec), 0.0);
	float m = sum / in_vec.size();
	return m;
}

nearfield_demod::sptr nearfield_demod::make() {
		return gnuradio::get_initial_sptr
			(new nearfield_demod_impl());
	}

/*
 * The private constructor
 */
nearfield_demod_impl::nearfield_demod_impl()
	: gr::sync_block("nearfield_demod",
			gr::io_signature::make(1, 1, sizeof(float)),
			gr::io_signature::make(0, 1, sizeof(float))) {

	// variables
	N = 5;                      % provided information
	threshold = 0.6;            % threshold set after observing data
	SDR_sample_rate = 12.5e6;   % sample rate of software defined radio receiver
	pulse_max = 550e-9;         % maximum pulse width
	pulse_min = 450e-9;         % minimum pulse width
	prf_max = 12e-3;            % maximum pulse repetition frequency
	prf_min = 10e-3;            % minimum pulse repetition freuqnecy
	
	// setup computations
	sample_period = 1/SDR_sample_rate;
	
	// Treating this like it's streaming data so I won't use some of the 
	// efficient Matlab functions that would speed this up.
	last_data = 0;      % previous received data point
	pulse_count = 0;    % length of current pulse (1s)
	prf_count = 0;      % length of current prf (0s)
	sync = 0;           % current consecutive sync pulses
	header = 4;         % # sync pulses needed
	sync_prf = 0;       % prf counter in sync part of loop
	sync_prf2 = 0;      % prf counter in sync part of loop (prf_window)
	sync_pulse = 0;     % pulse counter in sync part of loop
	prf_win_cnt = 0;    % counter for window of pulses between prf_min and prf_max in sync routine
	valid_pulse = 0;    % flag for pulse detection
	n = 0;              % counter for N

	pulse_vec.clear();
	prf_vec.clear();
	demod_data.clear();

	message_port_register_out(pmt::mp("frame_out"));
}

/*
 * Our virtual destructor.
 */
nearfield_demod_impl::~nearfield_demod_impl() {
}

int nearfield_demod_impl::work(int noutput_items,
		gr_vector_const_void_star &input_items,
		gr_vector_void_star &output_items) {

	const float *in = (const float *) input_items[0];
	float *out = (float *) output_items[0];

	int nn = 0;
	while(nn < noutput_items){
		//Local non-persistent variables
		float rx_data;
		float transition;
	
		// STEP 1 --------------------------------------------------------------
		if in[nn] > threshold
			rx_data = 1;
		else
			rx_data = 0;
	
		// STEP 2 --------------------------------------------------------------
		transition = rx_data - last_data;     // see if we are at an edge of a pulse
		// NOT SYNCHRONIZED YET
		if(sync < header) {
			if(transition > 0) {           // rising pulse
				prf_val = prf_count*sample_period;
				if(prf_val > prf_min && prf_val < prf_max){
					prf_vec.push_back(prf_val);
				} else {
					sync = 0;                     // reset the sync counter
					prf_vec.clear();                 // reset prf vector
					pulse_vec.clear();               // reset pulse vector
				}
				pulse_count = 1;     
				prf_count = 0;         
			} else if(transition < 0) {       // falling pulse
				pulse_val = pulse_count*sample_period;
				if(pulse_val > pulse_min && pulse_val < pulse_max) {
					sync = sync+1;              // valid pulse
					pulse_vec.push_back(pulse_val);
				} else {
					prf_count = pulse_count;    // it was noise...not pulse - count it towards the prf
				}
				pulse_count = 0;
				prf_count = prf_count+1;
			} else {    // no transition
				if(rx_data == 1)
					pulse_count = pulse_count+1;
				else
					prf_count = prf_count+1;
			}
			last_data = rx_data;
		}
		// SYNCHRONIZED
		if(sync == header) {
			// find average pulse width and prf
			pulse_length = roundf(mean(pulse_vec)/sample_period);
			prf_length = roundf(mean(prf_vec)/sample_period);
			// using this, demodulate the next N bits
			// right after last sync pulse is a prf window
			pulse_length_min = pulse_length - 1;
			pulse_length_max = pulse_length + 1;
			prf_length_min = roundf(prf_length - prf_length*0.1);
			prf_length_max = roundf(prf_length + prf_length*0.1);
			prf_window = prf_length_max - prf_length_min;
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
				if(sync_pulse >= pulse_length_min && sync_pulse <= pulse_length_max) {     // valid pulse
					valid_pulse = 1;
					sav_pulse = sync_pulse;
				}
				sync_pulse = 0;              // reset
			}
			if(valid_pulse == 1) {              // we found a pulse
				demod_data.push_back(1);
				n = n + 1;
				sync_prf = sav_pulse;        // prf is defined as rising edge of pulse to the next rising edge of pulse
				sync_prf2 = 0;
				valid_pulse = 0;             // reset
			} else if(sync_prf2 == prf_window) {   // ran through whole prf window w/o finding pulse
				demod_data.push_back(1);
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

			sync = 0;                        // start over looking for sync
			prf_vec.clear();
			pulse_vec.clear();
			prf_count = 0;
			pulse_count = 0;
			n = 0;
			sync_prf = 0;
			sync_prf2 = 0;
			sync_pulse = 0;
			valid_pulse = 0;
		}
	
		//Increment the data index pointer
		nn++;
	}

	// Tell runtime system how many output items we produced.
	return noutput_items;
}

} /* namespace nearfield */
} /* namespace gr */

