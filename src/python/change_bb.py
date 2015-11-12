#!/usr/bin/python
# Python script to change bb-freq
# 11/2015 Gyouho Kim

import sys
import os

change_num = sys.argv[1]
cmd = "vim ppm_header_cali.py '+g/bb-freq/s/default=\A\+,/default=" + change_num + ",' '+wq'"
print(cmd)
os.system(cmd)

