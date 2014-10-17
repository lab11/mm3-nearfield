#!/usr/bin/env python

import socket
import struct
import random

IP = 'inductor.eecs.umich.edu'
PORT = 4001

o = False
def val():
	global o
	if not o:
		o = True
		return 1
	return 0

def rval():
	return random.random() > 0.5

s1 = struct.pack('!10s 28B', 'N9NaoNAJzi', *[val() for x in range(28)])
s2 = struct.pack('!10s 28B', 'N9NaoNAJzi', *[rval() for x in range(28)])

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock.sendto(s1, (IP, PORT))
sock.sendto(s2, (IP, PORT))
