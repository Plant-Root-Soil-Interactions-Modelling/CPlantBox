"""Visualize the carrot secondary-growth radius model r(a, z) at a fixed age."""

import matplotlib.pyplot as plt
import numpy as np


def carrot_radius(age, z, r0=0.05, rmax=2.5, a_s=30.0, k=6.4e-4, p=1.5, Ls=15, q=2):
    """Secondary-growth radius model r(a, z) = r0 + (rmax - r0) * g(a) * h(z).

    g(a) = 1 - exp(-k * (a - a_s)^p) ramps up secondary growth once age a exceeds a_s.
    h(z) = exp(-(z / Ls)^q) concentrates thickening near the root crown (z=0).

    age: root age [day]
    z:   distance from the root crown [cm]
    r0:  primary-growth (unthickened) radius [cm]
    rmax: fully thickened radius [cm]
    a_s: age at which secondary growth starts [day]
    """
    if age <= a_s:
        return np.full_like(z, r0, dtype=float)
    g = 1.0 - np.exp(-k * (age - a_s) ** p)
    h = np.exp(-((z / Ls) ** q))
    print(f"age={age}, z={z}, g={g}, h={h}")
    return r0 + (rmax - r0) * g * h


ages = range(0, 141, 20)  # root ages [day]

addParams = {"r0": 0.05, "rmax": 3, "a_s": 30.0, "k": 6.4e-4, "p": 1.5, "Ls": 20, "q": 3}

z = np.linspace(0, 60, 120)  # distance from the root crown [cm]
colors = plt.cm.viridis(np.linspace(0, 1, len(ages)))

fig, ax = plt.subplots()
for age, color in zip(ages, colors):
    radius = np.array([carrot_radius(age, zi, **addParams) for zi in z])
    ax.plot(radius, z, color=color, label=f"age={age} d")
ax.set_xlabel("radius [cm]")
ax.set_ylabel("distance from root crown, z [cm]")
ax.set_title("carrot_radius over time")
ax.invert_yaxis()  # z=0 (crown) at the top
ax.legend()
ax.grid(True)
plt.tight_layout()
plt.show()
