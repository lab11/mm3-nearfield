#!/usr/bin/python
# Python script to parse Snail code radio output
# Looks for transmissions of Temp, VBAT, and Light measurement data
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
		if (i%3 == 1):
			lux_int.append(int(data,2))
		elif (i%3 == 2):
			temp_int.append(int(data,2))
		elif (i%3 == 0):
			batt_int.append(int(data,2))
			
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

fig1 = plt.figure()

plt.plot(range(1,len(batt_int)+1),list(reversed(batt_int)))
plt.ylabel('BATT')
plt.xlabel('Time (pts)')

plt.show()

