#!/usr/bin/env python

from gnuradio import gr, gru, uhd
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio.eng_option import eng_option
from optparse import OptionParser
import gnuradio.gr.gr_threading as _threading
from grc_gnuradio import wxgui as grc_wxgui
from gnuradio.wxgui import form, forms
import wx
import nearfield
import time

import struct
import sys

def make_normal(widget):
   font = widget.GetFont()
   font.SetWeight(wx.FONTWEIGHT_NORMAL)
   widget.SetFont(font)

class my_top_block(grc_wxgui.top_block_gui):
   def __init__(self, options):
	grc_wxgui.top_block_gui.__init__(self, title="mm3 Nearfield")

	self._down_gain = options.gain
	self._down_freq = options.freq
	self._down_bitrate = options.bitrate
	self._down_bitrate_accuracy = options.bitrate_accuracy
	self._pulse_len = options.pulse_len
	self._pulse_len_accuracy = options.pulse_len_accuracy
	self._header_len = options.header_len
	self._packet_len = options.packet_len
	self._sample_rate = options.sample_rate
	self._trigger_level = 1.0

	#Set up the GUI elements
	self.nb0 = wx.Notebook(self.GetWin(), style=wx.NB_TOP)
	self.nb0.AddPage(grc_wxgui.Panel(self.nb0), "RX")
	self.Add(self.nb0)
	self._carrier_text_box = forms.text_box(
		parent=self.nb0.GetPage(0).GetWin(),
		value=self._down_freq,
		callback=self.setDownFreq,
		label="Carrier Frequency",
		converter=forms.float_converter(),
	)
	self.nb0.GetPage(0).Add(self._carrier_text_box)
	self._gain_text_box = forms.text_box(
		parent=self.nb0.GetPage(0).GetWin(),
		value=self._down_gain,
		callback=self.setDownGain,
		label="Baseband Gain",
		converter=forms.float_converter(),
	)
	self.nb0.GetPage(0).Add(self._gain_text_box)
	self._bitrate_text_box = forms.text_box(
		parent=self.nb0.GetPage(0).GetWin(),
		value=self._down_bitrate,
		callback=self.setDownBitrate,
		label="Bitrate",
		converter=forms.float_converter(),
	)
	self.nb0.GetPage(0).Add(self._bitrate_text_box)
	self._bitrate_accuracy_text_box = forms.text_box(
		parent=self.nb0.GetPage(0).GetWin(),
		value=self._down_bitrate_accuracy,
		callback=self.setDownBitrateAccuracy,
		label="Bitrate Accuracy (%)",
		converter=forms.float_converter(),
	)
	self.nb0.GetPage(0).Add(self._bitrate_accuracy_text_box)
	self._pulse_len_text_box = forms.text_box(
		parent=self.nb0.GetPage(0).GetWin(),
		value=self._pulse_len,
		callback=self.setPulseLen,
		label="Pulse Length",
		converter=forms.float_converter(),
	)
	self.nb0.GetPage(0).Add(self._pulse_len_text_box)
	self._pulse_len_accuracy_text_box = forms.text_box(
		parent=self.nb0.GetPage(0).GetWin(),
		value=self._pulse_len_accuracy,
		callback=self.setPulseLenAccuracy,
		label="Pulse Length Accuracy (%)",
		converter=forms.float_converter(),
	)
	self.nb0.GetPage(0).Add(self._pulse_len_accuracy_text_box)
	self._header_len_text_box = forms.text_box(
		parent=self.nb0.GetPage(0).GetWin(),
		value=self._header_len,
		callback=self.setHeaderLen,
		label="Header Length (bits)",
		converter=forms.int_converter(),
	)
	self.nb0.GetPage(0).Add(self._header_len_text_box)
	self._packet_len_text_box = forms.text_box(
		parent=self.nb0.GetPage(0).GetWin(),
		value=self._packet_len,
		callback=self.setPacketLen,
		label="Packet Length (bits)",
		converter=forms.int_converter(),
	)
	self.nb0.GetPage(0).Add(self._packet_len_text_box)

	self._samplerate_text_box = forms.text_box(
		parent=self.nb0.GetPage(0).GetWin(),
		value=self._sample_rate,
		callback=self.setDownSampleRate,
		label="Sample Rate",
		converter=forms.float_converter(),
	)
	self.nb0.GetPage(0).Add(self._samplerate_text_box)
	self.prf_text = forms.static_text(
		parent=self.nb0.GetPage(0).GetWin(),
		value="",
		label="",
		converter=forms.str_converter(),
	)
	self.nb0.GetPage(0).Add(self.prf_text)
	self.pulse_len_text = forms.static_text(
		parent=self.nb0.GetPage(0).GetWin(),
		value="",
		label="",
		converter=forms.str_converter(),
	)
	self.nb0.GetPage(0).Add(self.pulse_len_text)
	self.threshold_text = forms.static_text(
		parent=self.nb0.GetPage(0).GetWin(),
		value="",
		label="",
		converter=forms.str_converter(),
		width=500
	)
	self.nb0.GetPage(0).Add(self.threshold_text)

	#self.source2 = blocks.file_source(gr.sizeof_gr_complex, "usrp_out_5mm_125e6_150e6.dat", True)
	#self.source = blocks.throttle(gr.sizeof_gr_complex,self._sample_rate)
	#self.connect(self.source2, self.source)
	self.source = uhd.usrp_source(device_addr=options.args, stream_args=uhd.stream_args('fc32'))
	self.source.set_center_freq(options.freq)
	self.source.set_samp_rate(self._sample_rate)
	self.source.set_gain(self._down_gain)

	self.mag = blocks.complex_to_mag_squared()

	self.tm_framer = nearfield.nearfield_demod(self._sample_rate, self._down_bitrate, self._down_bitrate_accuracy, self._pulse_len, self._pulse_len_accuracy, self._packet_len, self._header_len)
	self.socket_pdu = blocks.socket_pdu("TCP_SERVER", "127.0.0.1", "12912", 10000)
	self.msg_connect(self.tm_framer, "frame_out", self.socket_pdu, "pdus")

	self.connect(self.source, self.mag, self.tm_framer)

	_threading.Thread(target=self.watchFramer).start()

   def watchFramer(self):
	while True:
		time.sleep(1)
		obs_bitrate = self.tm_framer.getLastObservedBitrate()
		bit_fmt = '%.3e' % obs_bitrate
		if (obs_bitrate < self._down_bitrate*(1-self._down_bitrate_accuracy/100)):
			self.prf_text.set_value('Last Observed Bitrate: ' + bit_fmt + '  (TOO LOW)')
		elif (obs_bitrate > self._down_bitrate*(1+self._down_bitrate_accuracy/100)):
			self.prf_text.set_value('Last Observed Bitrate: ' + bit_fmt + '  (TOO HIGH)')
		else:
			self.prf_text.set_value('Last Observed Bitrate: ' + bit_fmt)

		obs_pl = self.tm_framer.getLastObservedPulseLen()
		pl_fmt = '%.3e' % obs_pl
		if (obs_pl < self._pulse_len*(1-self._pulse_len_accuracy/100)):
			self.pulse_len_text.set_value('Last Observed Pulse Length: ' + pl_fmt + '  (TOO LOW)')
		elif (obs_pl > self._pulse_len*(1+self._pulse_len_accuracy/100)):
			self.pulse_len_text.set_value('Last Observed Pulse Length: ' + pl_fmt + '  (TOO HIGH)')
		else:
			self.pulse_len_text.set_value('Last Observed Pulse Length: ' + pl_fmt)

		obs_thr = self.tm_framer.getThreshold()
		thr_fmt = '%.3e' % obs_thr
		if (obs_thr > 0.8):
			self.threshold_text.set_value('Current Threshold: ' + thr_fmt + '  (TOO HIGH)')
		elif (obs_thr < 0.1):
			self.threshold_text.set_value('Current Threshold: ' + thr_fmt + '  (TOO LOW)')
		else:
			self.threshold_text.set_value('Current Threshold: ' + thr_fmt)
		
   def setDownBitrate(self,arg):
	self._down_bitrate = float(arg)
	self.tm_framer.setBitrate(self._down_bitrate)

   def setDownFreq(self,arg):
	self._down_freq = float(arg)
	self.source.set_center_freq(self._down_freq)

   def setDownGain(self,arg):
	self._down_gain = float(arg)
	self.source.set_gain(self._down_gain)

   def setDownBitrateAccuracy(self,arg):
	self._down_bitrate_accuracy = float(arg)
	self.tm_framer.setBitrateAccuracy(self._down_bitrate_accuracy)

   def setPulseLen(self,arg):
	self._pulse_len = float(arg)
	self.tm_framer.setPulseLen(self._pulse_len)

   def setPulseLenAccuracy(self,arg):
	self._pulse_len_accuracy = float(arg)
	self.tm_framer.setPulseLenAccuracy(self._pulse_len_accuracy)

   def setHeaderLen(self,arg):
	self._header_len = int(arg)
	self.tm_framer.setHeaderLen(self._header_len)

   def setPacketLen(self,arg):
	self._packet_len = int(arg)
	self.tm_framer.setPacketLen(self._packet_len)

   def setDownSampleRate(self,arg):
	self._sample_rate = float(arg)
	self.source.set_samp_rate(self._sample_rate)
	self.tm_framer.setSampleRate(self._sample_rate)

def main():
	parser = OptionParser (option_class=eng_option, conflict_handler="resolve")
	expert_grp = parser.add_option_group("Expert")

        parser.add_option("-a", "--args", type="string", default="addr=192.168.10.2",
                          help="UHD device address args [default=%default]")
	parser.add_option("-f", "--freq", type="float", default=915e6, help="USRP carrier frequency [default=%default]")
	parser.add_option("-g", "--gain", type="float", default=100, help="USRP gain [default=%default]")
	parser.add_option("-b", "--bitrate", type="float", default=91, help="RX Bitrate [default=%default]")
	parser.add_option("-B", "--bitrate-accuracy", type="float", default=10, help="RX Bitrate Accuracy (%) [default=%default]")
	parser.add_option("-p", "--pulse-len", type="float", default=500e-9, help="Pulse Length [default=%default]")
	parser.add_option("-P", "--pulse-len-accuracy", type="float", default=80, help="Pulse Length Accuracy (%)")
	parser.add_option("-h", "--header-len", type="int", default=4, help="Header Length (bits)")
	parser.add_option("-n", "--packet-len", type="int", default=5, help="Packet Length (bits")
	parser.add_option("-s", "--sample-rate", type="float", default=12.5e6, help="RX Sample Rate [default=%default]")

	(options, args) = parser.parse_args ()
	
	if len(args) != 0:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# build the graph
	tb = my_top_block(options)

	r = gr.enable_realtime_scheduling()
	if r != gr.RT_OK:
		print "Warning: Failed to enable realtime scheduling."

	tb.Start(True)        # start flow graph
	tb.Wait()         # wait for it to finish

if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		pass
