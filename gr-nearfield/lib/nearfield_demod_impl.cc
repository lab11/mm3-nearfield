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

nearfield_demod::sptr nearfield_demod::make() {
		return gnuradio::get_initial_sptr
			(new nearfield_demod_impl());
	}

/*
 * The private constructor
 */
nearfield_demod_impl::nearfield_demod_impl()
	: gr::sync_block("nearfield_demod",
			gr::io_signature::make(<+MIN_IN+>, <+MAX_IN+>, sizeof(<+ITYPE+>)),
			gr::io_signature::make(<+MIN_OUT+>, <+MAX_OUT+>, sizeof(<+OTYPE+>)))
{}

/*
 * Our virtual destructor.
 */
nearfield_demod_impl::~nearfield_demod_impl() {
}

int nearfield_demod_impl::work(int noutput_items,
		gr_vector_const_void_star &input_items,
		gr_vector_void_star &output_items) {

	const <+ITYPE+> *in = (const <+ITYPE+> *) input_items[0];
	<+OTYPE+> *out = (<+OTYPE+> *) output_items[0];

	// Do <+signal processing+>
	// STEP 1 --------------------------------------------------------------
	if data_pwr(x) > threshold
		rx_data = 1;
	else
		rx_data = 0;

	// STEP 2 --------------------------------------------------------------
	transition = rx_data - last_data;     // see if we are at an edge of a pulse
	// NOT SYNCHRONIZED YET
	if sync < header       
		if transition > 0           // rising pulse
			prf_val = prf_count*sample_period;
			disp('prf_count');
			disp(prf_val);
			if prf_val > prf_min && prf_val < prf_max
				prf_vec = [prf_vec prf_val];
			else
				sync = 0;                     // reset the sync counter
				prf_vec = [];                 // reset prf vector
				pulse_vec = [];               // reset pulse vector
			pulse_count = 1;     
			prf_count = 0;         
		elseif transition < 0       // falling pulse
			pulse_val = pulse_count*sample_period;
			disp('pulse_count');
			disp(pulse_val);
			if pulse_val > pulse_min && pulse_val < pulse_max
				sync = sync+1;              // valid pulse
				pulse_vec = [pulse_vec pulse_val];
				disp('valid pulse');
			else
				prf_count = pulse_count;    // it was noise...not pulse - count it towards the prf
			pulse_count = 0;
			prf_count = prf_count+1;
		else    // no transition
			if rx_data == 1
				pulse_count = pulse_count+1;
			else
				prf_count = prf_count+1;
		last_data = rx_data;
	// SYNCHRONIZED
	if sync == header
		disp('x');
		disp(x);
		// find average pulse width and prf
		pulse_length = round(mean(pulse_vec)/sample_period);
		prf_length = round(mean(prf_vec)/sample_period);
		// using this, demodulate the next N bits
		// right after last sync pulse is a prf window
		pulse_length_min = pulse_length - 1;
		pulse_length_max = pulse_length + 1;
		prf_length_min = round(prf_length - prf_length*0.1);
		prf_length_max = round(prf_length + prf_length*0.1);
		prf_window = prf_length_max - prf_length_min;
		// advance to the edge of the prf_min window and look for a pulse
		// anywhere in the prf_max-prf_min frame.
		// then return to prf_length and repeat.
		if sync_prf < prf_length_min         // don't care until we get to the point where a pulse might come
			disp('sync_prf');
			disp(sync_prf);
			sync_prf = sync_prf + 1;
		else                                 // start looking for a pulse (we are in the prf window)
			disp('sync_prf2');
			disp(sync_prf2);
			sync_prf2 = sync_prf2 + 1;       // do this for prf window duration
		if rx_data == 1
			sync_pulse = sync_pulse + 1; // might be a pulse
		else
			if sync_pulse >= pulse_length_min & sync_pulse <= pulse_length_max     // valid pulse
				valid_pulse = 1;
				sav_pulse = sync_pulse;
			sync_pulse = 0;              // reset
		if valid_pulse == 1              // we found a pulse
			demod_data = [demod_data 1];
			n = n + 1;
			sync_prf = sav_pulse;        // prf is defined as rising edge of pulse to the next rising edge of pulse
			sync_prf2 = 0;
			valid_pulse = 0;             // reset
		elseif sync_prf2 == prf_window   // ran through whole prf window w/o finding pulse
			demod_data = [demod_data 0];
			n = n + 1;
			sync_prf = prf_length_max - prf_length;      // don't set to 0...reset to middle of window for timing
			sync_prf2 = 0;
		end
	end
	if n == N                            // we've looked for all the data
		sync = 0;                        // start over looking for sync
		prf_vec = [];
		pulse_vec = [];
		prf_count = 0;
		pulse_count = 0;
		n = 0;
		sync_prf = 0;
		sync_prf2 = 0;
		sync_pulse = 0;
		valid_pulse = 0;
	//            break;

	// Tell runtime system how many output items we produced.
	return noutput_items;
}

} /* namespace nearfield */
} /* namespace gr */

