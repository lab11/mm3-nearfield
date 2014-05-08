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
	float last_data;
	float pulse_count;
	float prf_count;
	float sync;
	float header;
	float sync_prf;
	float sync_prf2;
	float sync_pulse;
	float prf_win_cnt;
	float valid_pulse;
	float n;
	std::vector<float> pulse_vec;
	std::vector<float> prf_vec;
	std::vector<uint8_t> demod_data;
	float sav_pulse;

     public:
      nearfield_demod_impl();
      ~nearfield_demod_impl();

      // Where all the action really happens
      int work(int noutput_items,
	       gr_vector_const_void_star &input_items,
	       gr_vector_void_star &output_items);
    };

  } // namespace nearfield
} // namespace gr

#endif /* INCLUDED_NEARFIELD_NEARFIELD_DEMOD_IMPL_H */

