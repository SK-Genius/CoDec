import os
import time

for i in range(1, 8):
	name = 'LenaColor'
	fileIn  = open(name+'.dat', 'rb')
	fileOut = open(name+'_.dat', 'wb')
	fileOut.write(fileIn.read(i * 1000))
	fileIn.close()
	fileOut.close()
	os.system('main -d '+name+'_.dat '+name+'_.bmp')
	time.sleep(2)