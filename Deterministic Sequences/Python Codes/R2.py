# Reference: http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/

from scipy.io import savemat
import numpy as np
import matplotlib.pyplot as plt

# Using the above nested radical formula for g=phi_d 
# or you could just hard-code it. 
# phi(1) = 1.6180339887498948482 
# phi(2) = 1.32471795724474602596 
def phi(d): 
  x=2.0000 
  for i in range(10): 
    x = pow(1+x,1/(d+1)) 
  return x



# Number of dimensions. 
d=2

# number of required points 
n=1024

g = phi(d) 
alpha = np.zeros(d) 
for j in range(d): 
  alpha[j] = pow(1/g,j+1) %1 
z = np.zeros((n, d)) 

# This number can be any real number. 
# Common default setting is typically seed=0
# But seed = 0.5 is generally better. 
for i in range(n): 
  z[i] = (0.5 + alpha*(i+1)) %1 
print(z)
print("Size of Z:", z.size)
print("Shape of Z:", z.shape)
mymat={'z':z}
savemat('R2_1.mat',mymat)


f = plt.figure()
	# Add space between Main Title and the subplot results
f.subplots_adjust(hspace=0.4, top=0.85)

# Set the figure size
plt.rcParams["figure.figsize"] = [7.50, 6.50]
plt.rcParams["figure.autolayout"] = True

# Scatter plot
ax1 = f.add_subplot(1,2,1)
plt.scatter(z[:, 0], z[:, 1])
#plt.show()

# Random data of 100×3 dimension
data = np.array(np.random.random((1024, 2)))
ax2 = f.add_subplot(1,2,2)
plt.scatter(data[:, 0], data[:, 1])

plt.show()

