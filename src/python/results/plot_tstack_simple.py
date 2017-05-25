#!/usr/bin/python
# Python script to parse Snail code radio output
# Looks for transmissions of Temp, VBAT, and Light measurement data
# 8/9/2016 Inhee - two 24-bit data
# 7/14/2016 Gyouho Kim
# 5/25/2017 Gyouho Kim - adding A0 as the data header

import sys
import csv
import matplotlib.pyplot as plt
import numpy as np

print()

text_file = open(sys.argv[1], "r")
lines = text_file.readlines()
text_file.close()
output_file = sys.argv[1].replace('.txt','')
output_file += '.csv'
output_fig = sys.argv[1].replace('.txt','.png')
#print lines
#print len(lines)

date = []
time = []

temp_int = []

i=0
j=0

for line in lines:
	if line[0]=='0' or line[0]=='1':
		data = line.strip()
		data = data.replace(',','')
		data = data.replace(' ','')
		data_header = (int(data,2) & 0xFF0000)>>16
		print(data_header)
		if data_header == 0xA0:
			temp_int.append(int(data,2) & 0x00FFFF)
			timestamp_split = lines[j-2].split()
			time.append(timestamp_split[3])
			date.append(timestamp_split[1] + " " + timestamp_split[2])
			i += 1 # data index
	j += 1 # line index

#for i in range(0,min(len(batt_int),len(lux_int))):
#	cmeas_cal.append(float(batt_int[i]+temp_int[i]-batt_int[i])*0.5/float(lux_int[i])*100000)

print "\nTEMP:\n"
print temp_int
print len(temp_int)


# Export to CSV
rows = zip(date,time,temp_int)

wr = csv.writer(open(output_file,'w'), delimiter=',', lineterminator='\n')
wr.writerow(['DATE','TIME','TEMP'])
for row in rows:
	wr.writerow(row)

# Generate Plots using MatLab

len_temp_int = len(temp_int)+1

fig = plt.figure()

plt.plot(range(1,len_temp_int),list(reversed(temp_int)),'+-')
plt.ylabel('Temperature')
plt.xlabel('Time (pts)')

fig.savefig(output_fig)
plt.show()

