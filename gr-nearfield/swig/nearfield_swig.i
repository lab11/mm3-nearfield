/* -*- c++ -*- */

#define NEARFIELD_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "nearfield_swig_doc.i"

%{
#include "nearfield/nearfield_demod.h"
%}


%include "nearfield/nearfield_demod.h"
GR_SWIG_BLOCK_MAGIC2(nearfield, nearfield_demod);
