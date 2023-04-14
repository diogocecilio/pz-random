import matplotlib.pyplot as plt
import numpy as np
import os.path
datafem=[1.53,1.58,0.0,1.94,2.39,3.00]
datamoreg=[1.30, 1.40, 1.60, 1.90, 2.40, 3.00]
ord=[1., 0.8, 0.6, 0.4, 0.2, 0.]


 
# Plotting both the curves simultaneously
plt.plot(ord,datafem, marker = 'o',color='r', label='FEM')
plt.plot(ord,datamoreg,  marker = '*',color='g', label='Michalowisck 1963')
  
# Naming the x-axis, y-axis and the whole graph
plt.ylabel('FOS')
plt.xlabel('hw/h')
plt.title(" $ \ frac{c}{\gamma H} $   2:1")
  
# Adding legend, which helps us recognize the curve according to it's color
plt.legend()
  
# To load the display window
plt.show()
  

