from ase.io import read
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Map of the Hg adsorption energy over one t-Se(10-10) surface unit cell (from the site scan),
# with the positions of the topmost Se atoms.

d = np.loadtxt('Hg_tSe_site_scan.dat')
slab = read('../02_tSe_slab/tSe_slab_relaxed.vasp')

nx, ny = int(d[:, 0].max()) + 1, int(d[:, 1].max()) + 1
x = d[:, 2].reshape(nx, ny)
y = d[:, 3].reshape(nx, ny)
e = d[:, 5].reshape(nx, ny) * 96.485   # kJ/mol

fig, ax = plt.subplots(figsize=(5.5, 5))
m = ax.pcolormesh(x, y, e, shading='nearest', cmap='viridis')
fig.colorbar(m, ax=ax, label='E_ads (kJ/mol), frozen slab')

# top-layer Se atoms near the scanned cell
z = slab.positions[:, 2]
top = slab.positions[z > z.max() - 1.8]
pad = 1.0
sel = (top[:, 0] > x.min() - pad) & (top[:, 0] < x.max() + pad) & \
      (top[:, 1] > y.min() - pad) & (top[:, 1] < y.max() + pad)
ax.scatter(top[sel, 0], top[sel, 1], s=120, facecolors='none', edgecolors='orange',
           linewidths=2, label='topmost Se')
i = np.argmin(d[:, 5])
ax.plot(d[i, 2], d[i, 3], 'r*', ms=14, label='best grid point')
ax.set_xlabel('x (A)  across the chains')
ax.set_ylabel('y (A)  along the chain axis')
ax.set_title('Hg on t-Se(10-10), MACE-MP-0+D3')
ax.legend(loc='upper right', fontsize=8)
ax.set_aspect('equal')
fig.tight_layout()
fig.savefig('Hg_tSe_site_scan.png', dpi=150)
print("saved Hg_tSe_site_scan.png")
