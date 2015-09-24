import os
import time

n1 = 10
n2 = 10
for i in range(1, 38):
	name = 'LenaColor'
	fileIn  = open(name+'.dat', 'rb')
	fileOut = open(name+'_.dat', 'wb')
	size = n1 + n2
	n1 = n2
	n2 = size
	print(size)
	fileOut.write(fileIn.read(size))
	fileIn.close()
	fileOut.close()
	os.system('main -d '+name+'_.dat '+name+'_.bmp')
	time.sleep(2)
