import matplotlib.pyplot as plt
import numpy as np
import os.path
datafem=[1.05, 1.03, 1.05, 1.09, 1.13, 1.17, 1.26, 1.36, 1.57, 1.94, 2.18]
datamoreg=[0.90, 0.95, 1.00, 1.05, 1.11, 1.22, 1.38, 1.55, 1.72, 1.92, 2.20]
ord=[1., 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.]


 
# Plotting both the curves simultaneously
plt.plot(ord,datafem, marker = 'o',color='r', label='FEM')
plt.plot(ord,datamoreg,  marker = '*',color='g', label='Moregenstern 1963')
  
# Naming the x-axis, y-axis and the whole graph
plt.ylabel('FOS')
plt.xlabel('hw/h')
plt.title(" $ \ frac{c}{\gamma H} $   2:1")
  
# Adding legend, which helps us recognize the curve according to it's color
plt.legend()
  
# To load the display window
plt.show()
  

