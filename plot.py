# -*- coding: utf-8 -*-

from scipy.io import FortranFile
import numpy as np
import matplotlib.pyplot as plt

f1 = FortranFile('outputs/r_lims', 'r')
f2 = FortranFile('outputs/theta_lims', 'r')
f3 = FortranFile('outputs/values', 'r')

a = f1.read_reals(dtype='float32')
b = f2.read_reals(dtype='float32')
c = f3.read_reals(dtype='float32')

num_values_a = len(a)//3

# Data came from the form (x, i+1)
# Data read in as linear array, varying by 2nd index first, so must transform
# x = 0 corresponds to R_a
# x = 1, 2 correspond to r_a_min, r_a_max
a = np.reshape(a, [num_values_a, 3]).T

# Plot R_a against r_a (min, max)
fig0, ax0 = plt.subplots()
ax0.scatter(a[0, :], a[1, :], label=r'$r_{\alpha, min}$', color='tab:blue', linewidth=1, marker='.')
ax0.scatter(a[0, :], a[2, :], label=r'$r_{\alpha, max}$', color='tab:red', linewidth=1, marker='.')
ax0.set_xlabel(r'$R_\alpha$')
ax0.set_ylabel(r'$r_\alpha$')
ax0.legend()

title = r'$M_\alpha$' + f'={c[0]:.0f}, ' + \
        r'$M_\beta$' + f'={c[1]:.0f}, ' + \
        r'$M_\gamma$' + f'={c[2]:.0f}, ' + \
        r'$R_\alpha$' + f'={c[3]:.0f} - {c[4]:.0f}, ' + \
        r'$R_\beta$' + f'={c[5]:.0f} - {c[6]:.0f}, ' + \
        r'$\gamma_{\alpha\beta}$' + f'={c[7]:.0f} - ' + r'$\pi$'

ax0.set_title(title)
# ax0.set_ylim([-0.5, 10.5])
ax0.xaxis.get_ticklocs(minor=True)
ax0.minorticks_on()
fig0.savefig('figures/R_a_r_a_limits.pdf', dpi=400)

num_values_b = int(np.sqrt(len(b)//4))
# Data in the form (x, i+1, j+1)
# x = 0, 1 correspond to R_a, r_a
# x = 2, 3 correspond to theta_a_min and theta_a_max
#  Second index changes R_a, third index changes r_a
b = np.reshape(b, [num_values_b, num_values_b, 4]).T

# Plot R_a against theta_a (min, max) for "fixed" r_a
fig1, ax1 = plt.subplots()
r_a_index = int(num_values_b*0.9)
ax1.set_xlabel(r'$R_\alpha$')
ax1.set_ylabel(r'$\theta_\alpha$')
ax1.scatter(b[0, :, r_a_index], b[2, :, r_a_index], label=r'$\theta_{\alpha, min}$', color='tab:blue', linewidth=1, marker='.')
ax1.scatter(b[0, :, r_a_index], b[3, :, r_a_index], label=r'$\theta_{\alpha, max}$', color='tab:red', linewidth=1, marker='.')
ax1.legend()
# ax1.set_ylim([-0.2, 3.5])
ax1.set_title(title)
fig1.savefig('figures/R_a_theta_a_limits.pdf', dpi=400)

# Plor r_a against theta_a (min, max) for fixed R_a
fig2, ax2 = plt.subplots()
R_a_index = int(num_values_b*0.9)
ax2.scatter(b[1, R_a_index, :], b[2, R_a_index, :], label=r'$\theta_{\alpha, min}$', color='tab:blue', linewidth=1, marker='.')
ax2.scatter(b[1, R_a_index, :], b[3, R_a_index, :], label=r'$\theta_{\alpha, max}$', color='tab:red', linewidth=1, marker='.')
ax2.set_xlabel(r'$r_\alpha$')
ax2.set_ylabel(r'$\theta_\alpha$')
# ax2.set_ylim([-0.2, 3.5])
ax2.legend()
ax2.set_title(title)
fig2.savefig('figures/small_r_a_theta_a_limits.pdf', dpi=400)

b = np.reshape(b, [4, num_values_b**2])

# Plot R_a against r_a and theta_a
fig3 = plt.figure()
ax3 = plt.axes(projection ='3d')
# ::n reduce points to make plot less dense
y = b[0, ::10]
x = b[1, ::10]
z1 = b[2, ::10]
z2 = b[3, ::10]
c1 = x + y
c2 = c1

ax3.scatter(x, y, z1, c='tab:blue', alpha=0.8, linewidth=0.5, marker='x', label=r'$\theta_{\alpha, min}$')
ax3.scatter(x, y, z2, c='tab:red', alpha=0.5, linewidth=0.1, marker='o', label=r'$\theta_{\alpha, max}$')
ax3.view_init(15, -115)
ax3.set_xlabel(r'$r_\alpha$')
ax3.set_ylabel(r'$R_\alpha$')
ax3.set_zlabel(r'$\theta_\alpha$')
# ax3.set_zlim([0, 3.5])
ax3.legend()
ax3.set_title(title)
fig3.savefig('figures/3d_limits.pdf', dpi=400)

f1.close()
f2.close()
f3.close()