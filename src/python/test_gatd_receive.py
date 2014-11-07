#!/usr/bin/env python

import socket
import struct
import random

HOST = 'localhost'
PORT = 4001

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock.bind((HOST, PORT))

print("Waiting for data...")
while True:
	data, addr = sock.recvfrom(1024)
	print("Received Message: ({} bytes)".format(len(data)))
	print("\traw: {}".format(data))
	print("\thex: {}".format(data.encode('hex')))
