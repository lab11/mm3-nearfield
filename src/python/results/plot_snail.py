#!/usr/bin/python
# Python script to parse Snail code radio output
# Looks for transmissions of Temp, VBAT, and Light measurement data
# 8/9/2016 Inhee - two 24-bit data
# 7/14/2016 Gyouho Kim

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
#print lines
#print len(lines)

date = []
time = []

batt_int = []
lux_int = []
temp_int = []

i=0
j=0

for line in lines:
	if line[0]=='0' or line[0]=='1':
		data = line.strip()
		data = data.replace(',','')
		data = data.replace(' ','')
		if (i%2 == 1):
			batt_int.append(int(data,2) >> 16)
			temp_int.append(int(data,2) & 0x00FFFF)
		elif (i%2 == 0):
			lux_int.append(int(data,2))
			
			timestamp_split = lines[j-2].split()
			time.append(timestamp_split[3])
			date.append(timestamp_split[1] + " " + timestamp_split[2])

		i += 1 # data index
	j += 1 # line index

#for i in range(0,min(len(batt_int),len(lux_int))):
#	cmeas_cal.append(float(batt_int[i]+temp_int[i]-batt_int[i])*0.5/float(lux_int[i])*100000)

print "\nBATT:\n"
print batt_int
print "\nLUX:\n"
print lux_int
print "\nTEMP:\n"
print temp_int
print len(batt_int)
print len(lux_int)
print len(temp_int)


# Export to CSV
rows = zip(date,time,batt_int,lux_int,temp_int)

wr = csv.writer(open(output_file,'w'), delimiter=',', lineterminator='\n')
wr.writerow(['DATE','TIME','BAT','LUX','TEMP'])
for row in rows:
	wr.writerow(row)

# Generate Plots using MatLab

len_batt_int = len(batt_int)+1
len_lux_int = len(lux_int)+1
len_temp_int = len(temp_int)+1

fig = plt.figure()

plt.subplot(3,1,1)
plt.title('Stack Testing Results for Snail Project')
plt.plot(range(1,len_batt_int),list(reversed(batt_int)),'+-')
plt.ylabel('Battery Voltage')
plt.xlabel('Time (pts)')

plt.subplot(3,1,2)
plt.plot(range(1,len_lux_int),list(reversed(lux_int)),'+-')
plt.ylabel('Light Intensity')
plt.xlabel('Time (pts)')


plt.subplot(3,1,3)
plt.plot(range(1,len_temp_int),list(reversed(temp_int)),'+-')
plt.ylabel('Temperature')
plt.xlabel('Time (pts)')

plt.show()

