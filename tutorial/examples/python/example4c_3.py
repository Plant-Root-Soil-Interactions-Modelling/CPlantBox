import py_rootbox as rb
import numpy as np
from rb_tools import *
import matplotlib.pyplot as plt

N = 4  # layers between a and b
a = -3  # cm
b = -7  # cm
scale_elongation = rb.EquidistantGrid1D(a, b, N + 1)
scales = range(0, N)
scale_elongation.data = a2v(scales)  # set proportionality factors

z_ = np.linspace(0, -10, 1000)
y_ = np.zeros((1000))

for i, z in enumerate(z_):
    y_[i] = scale_elongation.getValue(rb.Vector3d(0, 0, z))

plt.plot(z_, y_)

plt.plot(a, 0, "r*")  # interval borders
plt.plot(b, N - 1, "r*")

plt.xlabel("z (cm)")
plt.ylabel("y value (?)")

plt.show()

print("done.")
