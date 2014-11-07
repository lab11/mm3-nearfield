#!/usr/bin/env python

import socket
import struct
import random

IP = 'inductor.eecs.umich.edu'
#IP = 'localhost'
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

def r():
	return int(round(28 + random.random() * 4))

l = r()
print("sending packets with len {}".format(l))
s1 = struct.pack('!10s {}B'.format(l), 'N9NaoNAJzi', *[val() for x in range(l)])
s2 = struct.pack('!10s {}B'.format(l), 'N9NaoNAJzi', *[rval() for x in range(l)])

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock.sendto(s1, (IP, PORT))
sock.sendto(s2, (IP, PORT))
