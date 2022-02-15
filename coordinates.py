# -*- coding: utf-8 -*-

import math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

mass_a = 2
mass_b = 3
mass_c = 4

x_a = [0, 0, 0]
R_a = 2.97
r_a = 5.66
gamma_a = np.pi/1.1

# R_a starts at the centre of mass of B and C (start of R_a) and ends at A (0, 0)
R_a_x_1 = -R_a
R_a_x_2 = 0
R_a_y_1 = 0
R_a_y_2 = 0

# r_a starts at B and ends at C, meeting the start of R_a

# Length from midpoint to B and C
r_a_B_length = mass_c * r_a / (mass_c + mass_b)
r_a_C_length = mass_b * r_a / (mass_c + mass_b)

# Add sin/cos * length to start of R_a to get positions of B and C
r_a_x_1 = R_a_x_1 - r_a_B_length * math.cos(gamma_a)
r_a_x_2 = R_a_x_1 + r_a_C_length * math.cos(gamma_a)
r_a_y_1 = R_a_y_1 - r_a_B_length * math.sin(gamma_a)
r_a_y_2 = R_a_y_1 + r_a_C_length * math.sin(gamma_a)

# r_b starts at C (end of r_a) and ends at A (end of R_a)
r_b_x_1 = r_a_x_2
r_b_x_2 = R_a_x_2
r_b_y_1 = r_a_y_2
r_b_y_2 = R_a_y_2

# r_c starts at A (end of R_a) and ends at C (start of r_a)
r_c_x_1 = R_a_x_2
r_c_x_2 = r_a_x_1
r_c_y_1 = R_a_y_2
r_c_y_2 = r_a_y_1

# R_b starts at centre of mass of A and C and ends at B
# Lengths of r_b with x and y components
r_b_x_length = (r_b_x_2 - r_b_x_1)
r_b_y_length = (r_b_y_2 - r_b_y_1)
r_b = math.sqrt(r_b_x_length**2 + r_b_y_length**2)

# Length from C to centre of mass of A and C
r_b_C_length = mass_a * r_b / (mass_a + mass_c)

# Fraction along C -> A centre of mass is
r_b_C_frac = r_b_C_length / r_b

R_b_x_1 = r_b_x_1 + r_b_C_frac * r_b_x_length
R_b_x_2 = r_a_x_1
R_b_y_1 = r_b_y_1 + r_b_C_frac * r_b_y_length
R_b_y_2 = r_a_y_1

# R_c starts at centre of mass of A and B and ends at C
# Lengths of r_c with x and y components
r_c_x_length = (r_c_x_2 - r_c_x_1)
r_c_y_length = (r_c_y_2 - r_c_y_1)
r_c = math.sqrt(r_c_x_length**2 + r_c_y_length**2)

# Length from A to centre of mass of A and B
r_c_A_length = mass_b * r_c / (mass_a + mass_b)

# Fraction along B -> A centre of mass is
r_c_A_frac = r_c_A_length / r_c

R_c_x_1 = r_c_x_1 + r_c_A_frac * r_c_x_length
R_c_x_2 = r_a_x_2
R_c_y_1 = r_c_y_1 + r_c_A_frac * r_c_y_length
R_c_y_2 = r_a_y_2

fig, ax = plt.subplots()

# Line2D: [x_1, x_2], [y_1, y_2]

# R_a
line_1 = Line2D([R_a_x_1, R_a_x_2], [R_a_y_1, R_a_y_2], linewidth=1, linestyle = "-", color="green")
ax.annotate(r'$R_\alpha$', xy=(0,0), xytext=((R_a_x_1 + R_a_x_2) / 2, (R_a_y_1 + R_a_y_2) / 2))
ax.add_line(line_1)

# r_a
line_2 = Line2D([r_a_x_1, r_a_x_2], [r_a_y_1, r_a_y_2], linewidth=1, linestyle = "-", color="green")
ax.annotate(r'$r_\alpha$', xy=(0,0), xytext=((r_a_x_1 + r_a_x_2) / 2, (r_a_y_1 + r_a_y_2) / 2))
ax.add_line(line_2)

# r_b
line_3 = Line2D([r_b_x_1, r_b_x_2], [r_b_y_1, r_b_y_2], linewidth=1, linestyle = "-", color="green")
ax.annotate(r'$r_\beta$', xy=(0,0), xytext=((r_b_x_1 + r_b_x_2) / 2, (r_b_y_1 + r_b_y_2) / 2))
ax.add_line(line_3)

# r_c
line_4 = Line2D([r_c_x_1, r_c_x_2], [r_c_y_1, r_c_y_2], linewidth=1, linestyle = "-", color="green")
ax.annotate(r'$r_\gamma$', xy=(0,0), xytext=((r_c_x_1 + r_c_x_2) / 2, (r_c_y_1 + r_c_y_2) / 2))
ax.add_line(line_4)

# R_b
line_5 = Line2D([R_b_x_1, R_b_x_2], [R_b_y_1, R_b_y_2], linewidth=1, linestyle = "-", color="green")
ax.annotate(r'$R_\beta$', xy=(0,0), xytext=((R_b_x_1 + R_b_x_2) / 2, (R_b_y_1 + R_b_y_2) / 2))
ax.add_line(line_5)

# R_c
# line_6 = Line2D([R_c_x_1, R_c_x_2], [R_c_y_1, R_c_y_2], linewidth=1, linestyle = "-", color="green")
# ax.annotate(r'$R_\gamma$', xy=(0,0), xytext=((R_c_x_1 + R_c_x_2) / 2, (R_c_y_1 + R_c_y_2) / 2))
# ax.add_line(line_6)

# Masses A, B, C)
ax.annotate("A", xy=(R_a_x_2 + 0.1, R_a_y_2 + 0.1))
plt.scatter(R_a_x_2, R_a_y_2, color='r')

ax.annotate("B", xy=(r_a_x_1 + 0.1, r_a_y_1 + 0.1))
plt.scatter(r_a_x_1, r_a_y_1, color='r')

ax.annotate("C", xy=(r_a_x_2 + 0.1, r_a_y_2 + 0.1))
plt.scatter(r_a_x_2, r_a_y_2, color='r')

ax.set_xlim(-9, 2)
ax.set_ylim(-2, 2)
plt.draw()
plt.show()

fig.savefig('figures/coords.png', dpi=200)

# Lengths of R_b with x and y components
R_b_x_length = (R_b_x_2 - R_b_x_1)
R_b_y_length = (R_b_y_2 - R_b_y_1)
R_b = math.sqrt(R_b_x_length**2 + R_b_y_length**2)
print(f'R_b = {R_b}')

calc_gamma_a = (r_a_x_2 - r_a_x_1) * (R_a_x_2 - R_a_x_1) + (r_a_y_2 - r_a_y_1) * (R_a_y_2 - R_a_y_1)
calc_gamma_a /= (R_a * r_a)
calc_gamma_a = math.acos(calc_gamma_a)
frac_gamma_a = np.pi / gamma_a
frac_calc_gamma_a = np.pi / calc_gamma_a
print(f'Input gamma_a = {gamma_a} (~pi/{frac_gamma_a:.1f}) rad')
print(f'Output gamma_a = {calc_gamma_a} (~pi/{frac_calc_gamma_a:.1f}) rad')

gamma_ab = (R_b_x_2 - R_b_x_1) * (R_a_x_2 - R_a_x_1) + (R_b_y_2 - R_b_y_1) * (R_a_y_2 - R_a_y_1)
gamma_ab /= (R_a * R_b)
gamma_ab = math.acos(gamma_ab)
frac_gamma_ab = np.pi / gamma_ab
print(f'Output gamma_ab = {gamma_ab} (~pi/{frac_gamma_ab:.1f}) rad')