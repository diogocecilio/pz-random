import matplotlib.pyplot as plt
import numpy as np
import os.path
datavec=[]

summ=0.

for i in range(1000):
	a = "fs";
	a+=str(i)
	a+=".dat"
	file_exists = os.path.exists(a)
	print(i)
	print(file_exists)
	if file_exists==True:
		f = open(a, "r")
		val=f.read(7)
		if val!='':
			summ+=float(val)
			mean=summ/(i+1)
			datavec.append(mean)

	else:
		
		print(file_exists)
		
for i in datavec:
 print(i)
 
plt.plot(datavec)
plt.ylabel('Factor of Safety')
plt.xlabel('Sample')
plt.show()
