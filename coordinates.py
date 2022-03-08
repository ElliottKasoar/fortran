# -*- coding: utf-8 -*-

import math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from matplotlib.patches import Arc
from scipy.io import FortranFile

def draw_line(ax, x_1, x_2, y_1, y_2, label, linewidth=0.5, linestyle='-', color='black', offset=[0, 0]):

    # Line2D: [x_1, x_2], [y_1, y_2]
    line = Line2D([x_1, x_2], [y_1, y_2], linewidth=linewidth, linestyle=linestyle, color=color)
    ax.annotate(label, xy=(0,0), xytext=(offset[0] + (x_1 + x_2) / 2, offset[1] + (y_1 + y_2) / 2))
    ax.add_line(line)

    return

def draw_mass(ax, x, y, label, color='r', offset=[0.1, 0.1]):

    ax.annotate(label, xy=(x + offset[0], y + offset[1]))
    plt.scatter(x, y, color=color)

    return


def draw_angle(ax, theta_2, origin=[0, 0], len_x_axis=1, len_y_axis=1, offset=1, theta_1=0, angle=0, color='b', label='', label_offset=[0, 0]):

    # Convert to degrees
    theta_1 = math.degrees(theta_1)
    theta_2 = math.degrees(theta_2)
    angle = math.degrees(angle)

    angle_plot = Arc(origin, len_x_axis*offset, len_y_axis*offset, angle, theta_1, theta_2, color=color, label=label)
    ax.add_patch(angle_plot)
    ax.annotate(label, xy=(0, 0), xytext=(origin[0] + label_offset[0], origin[1] + label_offset[1]))

    return


def calc_single_coords(mass_a, mass_b, mass_c, R_a_mag, r_a_mag, gamma_a):
    # R_a starts at the centre of mass of B and C (start of R_a) and ends at A (0, 0)
    R_a_x_1 = -R_a_mag
    R_a_x_2 = 0
    R_a_y_1 = 0
    R_a_y_2 = 0

    # r_a starts at B and ends at C, meeting the start of R_a

    # Length from midpoint to B and C
    r_a_B_length = mass_c * r_a_mag / (mass_c + mass_b)
    r_a_C_length = mass_b * r_a_mag  / (mass_c + mass_b)

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

    R_a = [R_a_x_1, R_a_x_2, R_a_y_1, R_a_y_2]
    r_a = [r_a_x_1, r_a_x_2, r_a_y_1, r_a_y_2]
    r_b = [r_b_x_1, r_b_x_2, r_b_y_1, r_b_y_2]
    r_c = [r_c_x_1, r_c_x_2, r_c_y_1, r_c_y_2]
    R_b = [R_b_x_1, R_b_x_2, R_b_y_1, R_b_y_2]
    R_c = [R_c_x_1, R_c_x_2, R_c_y_1, R_c_y_2]

    coords = R_a, r_a, r_b, r_c, R_b, R_c

    return coords


def calc_mixed_coords(mass_a, mass_b, mass_c, R_a_mag, R_b_mag, gamma_ab):

    M = mass_a + mass_b + mass_c
    m_a = mass_b * mass_c / (mass_b + mass_c)
    m_b = mass_a * mass_c / (mass_a + mass_c)
    mu_a = mass_a * mass_b * mass_c / (m_a * M)

    r_a_mag = (R_a_mag / mass_c)**2 + (R_b_mag/m_b)**2 + 2 * R_a_mag * R_b_mag * np.cos(gamma_ab) / (mass_c * m_b)
    r_a_mag = np.sqrt(r_a_mag)
    r_a_mag *= mu_a

    gamma_a = (R_a_mag / mass_c) + R_b_mag * np.cos(gamma_ab) / m_b
    gamma_a *= mu_a * np.sign(-gamma_ab) / r_a_mag
    gamma_a = np.arccos(gamma_a)

    coords = calc_single_coords(mass_a, mass_b, mass_c, R_a_mag, r_a_mag, gamma_a)

    return coords


def calc_mag(x):

    x_length = x[1] - x[0]
    y_length = x[3] - x[2]
    mag = math.sqrt(x_length**2 + y_length**2)

    return mag


def calc_angle(x, y):

    x_mag = calc_mag(x)
    y_mag = calc_mag(y)
    angle = (x[1] - x[0]) * (y[1] - y[0]) + (x[3] - x[2]) * (y[3] - y[2])
    angle /= (x_mag * y_mag)
    angle = math.acos(angle)

    return angle


def calc_grad_intercept(x):

    if (x[1] < x[0]):
        temp = x[1]
        x[1] = x[0]
        x[0] = temp

        temp = x[3]
        x[3] = x[2]
        x[2] = temp

    m = (x[3] - x[2]) / (x[1] - x[0])
    c = x[2] - m * x[0]

    return m, c


def calc_intersect(x, y):

    m_1, c_1 = calc_grad_intercept(x)
    m_2, c_2 = calc_grad_intercept(y)

    intersect = np.zeros(2)
    intersect[0] = c_1 - c_2 / (m_2 - m_1)
    intersect[1] = c_1 * m_2 - c_2 * m_1 / (m_2 - m_1)

    return intersect


def plot_all(R_a, r_a, r_b, r_c, R_b, R_c, gamma_a, gamma_ab, R_a_mag, r_a_mag, R_b_mag, consts, single_jacobi=True, save_plots=False):

    # Calulate min and max values for x and y to auto-scale
    x_coords = np.array([R_a[:2], r_a[:2], r_b[:2], r_c[:2], R_b[:2]])
    x_min = np.min(x_coords)
    x_max = np.max(x_coords)

    y_coords = np.array([R_a[2:], r_a[2:], r_b[2:], r_c[2:], R_b[2:]])
    y_min = np.min(y_coords)
    y_max = np.max(y_coords)

    x_min = round(x_min - 1)
    x_max = round(x_max + 1)
    y_min = round(y_min - 1)
    y_max = round(y_max + 1)

    fig, ax = plt.subplots()

    # Draw lines R_a, r_a, r_b, r_c and R_b (R_c not displayed to reduce clutter)
    draw_line(ax, R_a[0], R_a[1], R_a[2], R_a[3], r'$R_\alpha$', offset=[0, 0.1])
    draw_line(ax, r_a[0], r_a[1], r_a[2], r_a[3], r'$r_\alpha$', offset=[-0.2, -0.1])
    draw_line(ax, r_b[0], r_b[1], r_b[2], r_b[3], r'$r_\beta$', offset=[0, 0.1])
    draw_line(ax, r_c[0], r_c[1], r_c[2], r_c[3], r'$r_\gamma$', offset=[0, -0.2])
    draw_line(ax, R_b[0], R_b[1], R_b[2], R_b[3], r'$R_\beta$', offset=[0.05, -0.1])
    # draw_line(ax, R_c[0], R_c[1], R_c[2], R_c[3], r'$R_\gamma$', offset=[0, 0])

    # Draw masses A, B, C
    draw_mass(ax, R_a[1], R_a[3], 'A')
    draw_mass(ax, r_a[0], r_a[2], 'B')
    draw_mass(ax, r_a[1], r_a[3], 'C')

    # Draw angles gamma_a and gamma_ab
    draw_angle(ax, gamma_a, [R_a[0], R_a[2]], x_max - x_min, y_max - y_min, 0.08, label=r'$\gamma_\alpha$', label_offset=[0.1, 0.4])
    intersect = calc_intersect(R_a, R_b)
    draw_angle(ax, gamma_ab, intersect, y_max - y_min, x_max - x_min, 0.1, angle=-gamma_ab, label=r'$\gamma_{\alpha\beta}$', label_offset=[0.1, -0.5])

    ax.set_xlim(-6.75, 1.25)
    ax.set_ylim(-3.5, 4)

    ax.axis('off')

    if single_jacobi:
        text_1 = r'$R_\alpha$' + f' = {R_a_mag:.2f}, ' + r'$r_\alpha$' + f' = {r_a_mag:.2f}, ' + \
            r'$\gamma_\alpha$' + f' = {gamma_a:.2f}, ' + r'$R_\beta$' + f' = {R_b_mag:.2f}, ' + \
            r'$\gamma_{\alpha\beta}$' + f' = {gamma_ab:.2f}'
    else:
        text_1 = r'$R_\alpha$' + f' = {R_a_mag:.2f}, ' + r'$R_\beta$' + f' = {R_b_mag:.2f}, ' + \
            r'$\gamma_{\alpha\beta}$' + f' = {gamma_ab:.2f}, ' + r'$r_\alpha$' + f' = {r_a_mag:.2f}, ' + \
            r'$\gamma_\alpha$' + f' = {gamma_a:.2f}'

    text_2 = r'$M_\alpha$' + f'={consts[0]:.0f}, ' + \
            r'$M_\beta$' + f'={consts[1]:.0f}, ' + \
            r'$M_\gamma$' + f'={consts[2]:.0f}, ' + \
            r'$R_\alpha$' + f'={consts[3]:.0f} - {consts[4]:.0f}, ' + \
            r'$R_\beta$' + f'={consts[5]:.0f} - {consts[6]:.0f}, ' + \
            r'$\gamma_{\alpha\beta}$' + f'={consts[7]:.0f} - ' + r'$\pi$'

    ax.text(-6, 4, text_1)
    ax.text(-6, 5, text_2)

    plt.draw()
    plt.show()

    if save_plots:
        if single_jacobi:
            fig.savefig(f'figures/single/coords_R_a-{R_a_mag:.2f}_r_a-{r_a_mag:.2f}_gamma_a-{gamma_a:.2f}.pdf', dpi=400, bbox_inches = "tight")
        else:
            fig.savefig(f'figures/mixed/coords_R_a-{R_a_mag:.2f}_R_b-{R_b_mag:.2f}_gamma_ab-{gamma_ab:.2f}.pdf', dpi=400, bbox_inches = "tight")
    return


def main():

    f = FortranFile('outputs/values', 'r')
    consts = f.read_reals(dtype='float32')
    f.close()

    mass_a = 1
    mass_b = 1
    mass_c = 1

    # Used for single and mixed Jacobi
    R_a_mag = 2.97

    # Used for single Jacobi
    r_a_mag = 5.66
    gamma_a = np.pi / 3

    # Used for mixed Jacobi
    R_b_mag = 4.5
    gamma_ab = 2.7

    # Control flags
    verbose = True
    save_plots = True
    single_jacobi = False

    if single_jacobi:
        R_a, r_a, r_b, r_c, R_b, R_c = calc_single_coords(mass_a, mass_b, mass_c, R_a_mag, r_a_mag, gamma_a)

        R_b_mag = calc_mag(R_b)
        if verbose:
            print(f'R_b = {R_b_mag}')
    else:
        R_a, r_a, r_b, r_c, R_b, R_c = calc_mixed_coords(mass_a, mass_b, mass_c, R_a_mag, R_b_mag, gamma_ab)

        r_a_mag = calc_mag(r_a)
        if verbose:
            print(f'r_a = {r_a_mag}')

    # Calculate gamma_a = r_a . R_a / r_a R_a
    calc_gamma_a = calc_angle(r_a, R_a)

    # Calculate gamma_ab = R_b . R_a / R_a R_b
    calc_gamma_ab = calc_angle(R_b, R_a)

    if single_jacobi:
        gamma_ab = calc_gamma_ab
    else:
        gamma_a = calc_gamma_a

    if verbose:

        # Input gamma_a or gamma_ab as fraction of pi
        frac_gamma_a = np.pi / gamma_a
        frac_gamma_ab = np.pi / gamma_ab

        frac_calc_gamma_a = np.pi / calc_gamma_a
        frac_calc_gamma_ab = np.pi / gamma_ab

        if single_jacobi:
            print(f'Input gamma_a = {gamma_a} (~pi/{frac_gamma_a:.1f}) rad')
            print(f'Output gamma_a = {calc_gamma_a} (~pi/{frac_calc_gamma_a:.1f}) rad')
            print(f'Output gamma_ab = {calc_gamma_ab} (~pi/{frac_calc_gamma_ab:.1f}) rad')
        else:
            print(f'Input gamma_ab = {gamma_ab} (~pi/{frac_gamma_ab:.1f}) rad')
            print(f'Output gamma_ab = {calc_gamma_ab} (~pi/{frac_calc_gamma_ab:.1f}) rad')
            print(f'Output gamma_a = {calc_gamma_a} (~pi/{frac_calc_gamma_a:.1f}) rad')


    plot_all(R_a, r_a, r_b, r_c, R_b, R_c, gamma_a, gamma_ab, R_a_mag, r_a_mag, R_b_mag, consts, single_jacobi, save_plots)

    return

main()