import matplotlib.pyplot as plt
import numpy as np

# printaible domainn 196, 276 (in mm)
xx = 196.  # mm
yy = 276.  # mm

fig = plt.figure(figsize=(xx / 25.4, yy / 25.4))  # din A4 in inch
fig.add_axes([0., 0., 1., 1.]) 

N = 73  # every 5 degree
phi = np.linspace(0, 2 * np.pi, N)
x1 = np.array([ [500 * np.cos(a), 500 * np.sin(a)] for a in phi]) 
for i, x in enumerate(x1):
    if i % 9 == 0:
        plt.plot([0, x[0]], [0, x[1]], 'k')
    else:
        plt.plot([0, x[0]], [0, x[1]], 'b:')

N = 1000
phi = np.linspace(0, 2 * np.pi, N)
x2 = np.array([ [np.cos(a), np.sin(a)] for a in phi])

a = np.linspace(5, 200, 40)
for a_ in a:
    plt.plot(a_ * x2[:, 0], a_ * x2[:, 1], 'k', alpha=0.4)

plt.xlim([-xx / 2, xx / 2])  
plt.ylim([-yy / 2, yy / 2])
plt.show()
