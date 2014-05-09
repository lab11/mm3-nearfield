#!/usr/bin/env python

from gnuradio import gr, gru
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio.eng_option import eng_option
from optparse import OptionParser

import nearfield

import struct
import sys

class my_top_block(gr.top_block):
   def __init__(self, options):
	gr.top_block.__init__(self)

	self.source = blocks.file_source(gr.sizeof_gr_complex, "usrp_out_5mm_125e6_150e6.dat", True)
	#self.source = uhd.usrp_source(device_addr=options.args, stream_args=uhd.stream_args('fc32'))
	#self.u.set_samp_rate(12.5e6)
	#self.u.set_gain(100)

	self.mag = blocks.complex_to_mag_squared()

	self.tm_framer = nearfield.nearfield_demod()
	self.socket_pdu = blocks.socket_pdu("TCP_SERVER", "127.0.0.1", "12912", 10000)
	self.msg_connect(self.tm_framer, "frame_out", self.socket_pdu, "pdus")

	self.connect(self.source, self.mag, self.tm_framer)

def main():
	parser = OptionParser (option_class=eng_option, conflict_handler="resolve")
	expert_grp = parser.add_option_group("Expert")

        parser.add_option("-a", "--args", type="string", default="fpga=usrp_b200_dssdr.bin",
                          help="UHD device address args [default=%default]")

	(options, args) = parser.parse_args ()
	
	if len(args) != 0:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# build the graph
	tb = my_top_block(options)

	r = gr.enable_realtime_scheduling()
	if r != gr.RT_OK:
		print "Warning: Failed to enable realtime scheduling."

	tb.start()        # start flow graph
	tb.wait()         # wait for it to finish

if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		pass
