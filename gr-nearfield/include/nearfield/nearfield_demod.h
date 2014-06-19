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


#ifndef INCLUDED_NEARFIELD_NEARFIELD_DEMOD_H
#define INCLUDED_NEARFIELD_NEARFIELD_DEMOD_H

#include <nearfield/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace nearfield {

    /*!
     * \brief <+description of block+>
     * \ingroup nearfield
     *
     */
    class NEARFIELD_API nearfield_demod : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<nearfield_demod> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of nearfield::nearfield_demod.
       *
       * To avoid accidental use of raw pointers, nearfield::nearfield_demod's
       * constructor is in a private implementation
       * class. nearfield::nearfield_demod::make is the public interface for
       * creating new instances.
       */
      static sptr make(float sample_rate, float bitrate, float bitrate_accuracy, float pulse_len, float pulse_len_accuracy, int packet_len, int header_len);

      virtual float getLastObservedBitrate() = 0;
      virtual float getLastObservedPulseLen() = 0;
      virtual float getThreshold() = 0;
      virtual void setPulseLen(float pulse_len_in) = 0;
      virtual void setPulseLenAccuracy(float pulse_len_accuracy_in) = 0;
      virtual void setPostPulseLenAccuracy(float post_pulse_len_accuracy_in) = 0;
      virtual void setBitrate(float bitrate_in) = 0;
      virtual void setBitrateAccuracy(float bitrate_accuracy_in) = 0;
      virtual void setPostBitrateAccuracy(float post_bitrate_accuracy_in) = 0;
      virtual void setSampleRate(float sample_rate_in) = 0;
      virtual void setPacketLen(int packet_len_in) = 0;
      virtual void setHeaderLen(int header_len_in) = 0;
    };

  } // namespace nearfield
} // namespace gr

#endif /* INCLUDED_NEARFIELD_NEARFIELD_DEMOD_H */

