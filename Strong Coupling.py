import sys

sys.path.insert(0, '..')
import sif_reader
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import tkinter as tk
from tkinter.filedialog import askopenfilenames
from matplotlib.widgets import Slider, TextBox, Button

# set initial values and constants

h = 6.62607015 * 10 ** (-34)  # [Js]
hbar = h / (2 * np.pi)  # [Js]
c = 299792458  # [m/s]
me = 9.1093837 * 10 ** (-31)  # [kg]

LAMBDA_C = 838 * 10 ** (-9)  # [m]
EL = 1.602176634 * 10 ** (-19)  # [C]
EX0 = 1.484272 * EL  # [J]
EC0_1 = 1.4816 * EL  # [J]
EC0_2 = 1.4816 * EL  # [J]
EC0_3 = 1.49 * EL  # [J]
MC = 3.348 * 10 ** (-35)  # [kg]
RABI = 2.2 * 10 ** (-3) * EL / hbar  # [Hz]

# the (common) k-axis
NA = 0.42
idx_k0 = 494
idx_kmin = 294
idx_kmax = 692
k_max = 2 * np.pi * NA / LAMBDA_C  # [1/m]
# get proper step size for k-array
k_step = k_max / max(abs(idx_kmax - idx_k0), abs(idx_k0 - idx_kmin))
k_axis = np.linspace(-idx_k0 * k_step, (1023 - idx_k0) * k_step, 1024)
# the k's in the relevant range
k_relevant = np.linspace(k_axis[idx_kmin], k_axis[idx_kmax], 10000)


# load 3 files to be evaluated
def import_dat_data():
    # not the best style, but a global image storage will do the job for now
    global IMG_DATA1
    global IMG_DATA2
    global IMG_DATA3
    global IMG_ENERGIES1
    global IMG_ENERGIES2
    global IMG_ENERGIES3
    try:
        root = tk.Tk(className='Import sif-file')
        FILEPATHS = askopenfilenames()
        print(FILEPATHS)
        IMG_DATA1, IMG_INFO1 = sif_reader.np_open(FILEPATHS[0])
        IMG_DATA2, IMG_INFO2 = sif_reader.np_open(FILEPATHS[1])
        IMG_DATA3, IMG_INFO3 = sif_reader.np_open(FILEPATHS[2])
        IMG_DATA1 = np.asarray(IMG_DATA1, dtype=float)
        IMG_DATA2 = np.asarray(IMG_DATA2, dtype=float)
        IMG_DATA3 = np.asarray(IMG_DATA3, dtype=float)
        IMG_WAVELENGTHS1 = sif_reader.utils.extract_calibration(IMG_INFO1) * 10 ** (-9)  # [m]
        IMG_WAVELENGTHS2 = sif_reader.utils.extract_calibration(IMG_INFO2) * 10 ** (-9)  # [m]
        IMG_WAVELENGTHS3 = sif_reader.utils.extract_calibration(IMG_INFO3) * 10 ** (-9)  # [m]
        IMG_ENERGIES1 = h * c / IMG_WAVELENGTHS1 / EL
        IMG_ENERGIES2 = h * c / IMG_WAVELENGTHS2 / EL
        IMG_ENERGIES3 = h * c / IMG_WAVELENGTHS3 / EL
        root.destroy()
    except:
        print('Wrong file')
        root.destroy()


def E_up(k, ex0, ec0, mc, rabi):
    # return 0.5 * (ex0 + ec0 + (hbar * k) ** 2 / (2 * mc) + np.sqrt(ec0 - ex0 + (hbar * k) ** 2 / (2 * mc) +
    #                                                                4 * (hbar * rabi) ** 2))
    return ((ex0 + ec0 + (k ** 2) * ((hbar ** 2) / (2 * mc)) + np.sqrt(
        (ec0 - ex0 + (k ** 2) * ((hbar ** 2) / (2 * mc))) ** 2 + 4 * (hbar ** 2) * rabi ** 2)) * 0.5)


def E_down(k, ex0, ec0, mc, rabi):
    global EL
    global massfactor
    global rabifactor
    global hbar
    # return 0.5 * (ex0 + ec0 + (hbar * k) ** 2 / (2 * mc) - np.sqrt(ec0 - ex0 + (hbar * k) ** 2 / (2 * mc) +
    #                                                                4 * (hbar * rabi) ** 2))
    return ((ex0 + ec0 + (k ** 2) * ((hbar ** 2) / (2 * mc)) - np.sqrt(
        (ec0 - ex0 + (k ** 2) * ((hbar ** 2) / (2 * mc))) ** 2 + 4 * (hbar ** 2) * rabi ** 2)) * 0.5)


def hop_Cdown(k, ex0, ec0, mc, rabi):
    global hbar
    return - (hbar * rabi) / np.sqrt(
        (hbar * rabi) ** 2 + (E_down(k, ex0, ec0, mc, rabi) - (ec0 + (hbar * k) ** 2 / (2 * mc))) ** 2)


def hop_Xdown(k, ex0, ec0, mc, rabi):
    global hbar
    return - (ec0 + (hbar * k) ** 2 / (2 * mc) - E_down(k, ex0, ec0, mc, rabi)) / np.sqrt(
        (hbar * rabi) ** 2 + (E_down(k, ex0, ec0, mc, rabi) - (ec0 + (hbar * k) ** 2 / (2 * mc))) ** 2)


# --- handlers for all sliders

def on_ec1_changed(event):
    global EC0_1
    global EL
    # the slider is in units of eV
    EC0_1 = slider_ec01.val * EL
    redraw_1()
    # print(EC0_1)


def on_ec2_changed(event):
    global EC0_2
    global EL
    # the slider is in units of eV
    EC0_2 = slider_ec02.val * EL
    redraw_2()
    # print(EC0_2)


def on_ec3_changed(event):
    global EC0_3
    global EL
    # the slider is in units of eV
    EC0_3 = slider_ec03.val * EL
    redraw_3()
    # print(EC0_3)


def on_mc_changed(event):
    global MC
    global massfactor
    # the slider is in units of 10e-5 electron masses
    MC = slider_mc.val * massfactor
    redraw_all()
    # print(MC)


def on_rabi_changed(event):
    global RABI
    global rabifactor
    # the slider is in units of hbar*Omega in eV
    RABI = slider_rabi.val / rabifactor
    redraw_all()
    # print(RABI)


def on_ex0_changed(event):
    global EX0
    global EL
    # the slider is in units of eV
    EX0 = slider_ex0.val * EL
    redraw_all()
    # print(EX0)


def redraw_1():
    global EX0, EC0_1, MC, RABI
    # get rid of all lines but not the underlying image
    for l in range(len(ax_im1.lines)):
        ax_im1.lines.pop()
    for l in range(len(ax_hop1.lines)):
        ax_hop1.lines.pop()
    ax_im1.plot(k_relevant * 10 ** (-6), E_down(k_relevant, EX0, EC0_1, MC, RABI) / EL, color='red')
    ax_im1.plot(k_relevant * 10 ** (-6), E_up(k_relevant, EX0, EC0_1, MC, RABI) / EL, color='blue')
    ax_hop1.plot(k_axis * 10 ** (-6), hop_Cdown(k_axis, EX0, EC0_1, MC, RABI) ** 2,
                 label=r'$C_L^2(k_{||})$', color='magenta')
    ax_hop1.plot(k_axis * 10 ** (-6), hop_Xdown(k_axis, EX0, EC0_1, MC, RABI) ** 2,
                 label=r'$X_L^2(k_{||})$', color='forestgreen')
    ax_hop1.legend(markerscale=3, fontsize=17)
    fig_1.canvas.draw()


def redraw_2():
    global EX0, EC0_2, MC, RABI
    # get rid of all lines but not the underlying image
    for l in range(len(ax_im2.lines)):
        ax_im2.lines.pop()
    for l in range(len(ax_hop2.lines)):
        ax_hop2.lines.pop()
    ax_im2.plot(k_relevant * 10 ** (-6), E_down(k_relevant, EX0, EC0_2, MC, RABI) / EL, color='red')
    ax_im2.plot(k_relevant * 10 ** (-6), E_up(k_relevant, EX0, EC0_2, MC, RABI) / EL, color='blue')
    ax_hop2.plot(k_axis * 10 ** (-6), hop_Cdown(k_axis, EX0, EC0_2, MC, RABI) ** 2, label=r'$C_L^2(k_{||})$',
                 color='magenta')
    ax_hop2.plot(k_axis * 10 ** (-6), hop_Xdown(k_axis, EX0, EC0_2, MC, RABI) ** 2, label=r'$X_L^2(k_{||})$',
                 color='forestgreen')
    ax_hop2.legend(markerscale=3, fontsize=17)
    fig_2.canvas.draw()


def redraw_3():
    global EX0, EC0_3, MC, RABI
    # get rid of all lines but not the underlying image
    for l in range(len(ax_im3.lines)):
        ax_im3.lines.pop()
    for l in range(len(ax_hop3.lines)):
        ax_hop3.lines.pop()
    ax_im3.plot(k_relevant * 10 ** (-6), E_down(k_relevant, EX0, EC0_3, MC, RABI) / EL, color='red')
    ax_im3.plot(k_relevant * 10 ** (-6), E_up(k_relevant, EX0, EC0_3, MC, RABI) / EL, color='blue')
    ax_hop3.plot(k_axis * 10 ** (-6), hop_Cdown(k_axis, EX0, EC0_3, MC, RABI) ** 2, label=r'$C_L^2(k_{||})$',
                 color='magenta')
    ax_hop3.plot(k_axis * 10 ** (-6), hop_Xdown(k_axis, EX0, EC0_3, MC, RABI) ** 2, label=r'$X_L^2(k_{||})$',
                 color='forestgreen')
    ax_hop3.legend(markerscale=3, fontsize=17)
    fig_3.canvas.draw()


def redraw_all():
    redraw_1()
    redraw_2()
    redraw_3()


if __name__ == '__main__':

    # ------------- LOAD DATA ------------- #

    # get file
    import_dat_data()

    # ----------------------- PLOT EVERYTHING ------------------------- #
    # ----- PLOT 1
    fig_1, (ax_im1, ax_hop1) = plt.subplots(1, 2, figsize=(13, 8))
    ax_im1.set_ylabel("E [eV]", fontsize=18)
    ax_im1.set_xlabel(r'$k_{||}$ $[\frac{1}{\mu m}$]', fontsize=18)
    ax_im1.tick_params(axis='both', width=2, length=5, which='major', labelsize=15.5)
    # set limits to not overshoot later when fitting
    ax_im1.set_ylim([IMG_ENERGIES1[-1], IMG_ENERGIES1[0]])
    ax_hop1.set_ylabel("Hopfield coefficients squared", fontsize=18)
    ax_hop1.set_xlabel(r'$k_{||}$ $[\frac{1}{\mu m}$]', fontsize=18)
    ax_hop1.tick_params(axis='both', width=2, length=7, which='major', labelsize=15.5, direction='in')
    for a in ['top', 'bottom', 'left', 'right']:
        ax_hop1.spines[a].set_linewidth(2.6)

    # display the corresp. image
    colmap = cm.get_cmap('Greens').reversed()
    im1 = ax_im1.imshow(np.rot90(IMG_DATA1[0]), interpolation='none', origin='lower',
                        extent=[k_axis[0] * 10 ** (-6), k_axis[-1] * 10 ** (-6), IMG_ENERGIES1[-1], IMG_ENERGIES1[0]],
                        aspect="auto", cmap=colmap)
    # use up as much space as possible
    fig_1.subplots_adjust(left=0.07, right=0.98, bottom=0.10, top=0.93, wspace=0.19, hspace=0.0)

    # add pretty colorbar
    divider1 = make_axes_locatable(ax_im1)
    cax1 = divider1.append_axes('right', size='5%', pad=0.05)
    fig_1.colorbar(im1, cax=cax1, orientation='vertical')
    cax1.set_ylabel("Intensity (arb. unit)", rotation=270, fontsize=18, labelpad=17.5)
    # hide ticks and labels for colorbar
    cax1.xaxis.set_major_locator(ticker.NullLocator())
    cax1.yaxis.set_major_locator(ticker.NullLocator())

    # ----- PLOT 2
    fig_2, (ax_im2, ax_hop2) = plt.subplots(1, 2, figsize=(13, 8))
    ax_im2.set_ylabel("E [eV]", fontsize=18)
    ax_im2.set_xlabel(r'$k_{||}$ $[\frac{1}{\mu m}$]', fontsize=18)
    ax_im2.tick_params(axis='both', width=2, length=5, which='major', labelsize=15.5)
    # set limits to not overshoot later when fitting
    ax_im2.set_ylim([IMG_ENERGIES2[-1], IMG_ENERGIES2[0]])
    ax_hop2.set_ylabel("Hopfield coefficients squared", fontsize=18)
    ax_hop2.set_xlabel(r'$k_{||}$ $[\frac{1}{\mu m}$]', fontsize=18)
    ax_hop2.tick_params(axis='both', width=2, length=7, which='major', labelsize=15.5, direction='in')
    for a in ['top', 'bottom', 'left', 'right']:
        ax_hop2.spines[a].set_linewidth(2.6)

    # display the corresp. image
    colmap = cm.get_cmap('Greens').reversed()
    im2 = ax_im2.imshow(np.rot90(IMG_DATA2[0]), interpolation='none', origin='lower',
                        extent=[k_axis[0] * 10 ** (-6), k_axis[-1] * 10 ** (-6), IMG_ENERGIES2[-1], IMG_ENERGIES2[0]],
                        aspect="auto", cmap=colmap)
    # use up as much space as possible
    fig_2.subplots_adjust(left=0.07, right=0.98, bottom=0.10, top=0.93, wspace=0.19, hspace=0.0)

    # add pretty colorbar
    divider2 = make_axes_locatable(ax_im2)
    cax2 = divider2.append_axes('right', size='5%', pad=0.05)
    fig_2.colorbar(im2, cax=cax2, orientation='vertical')
    cax2.set_ylabel("Intensity (arb. unit)", rotation=270, fontsize=18, labelpad=17.5)
    # hide ticks and labels for colorbar
    cax2.xaxis.set_major_locator(ticker.NullLocator())
    cax2.yaxis.set_major_locator(ticker.NullLocator())

    # ----- PLOT 3
    fig_3, (ax_im3, ax_hop3) = plt.subplots(1, 2, figsize=(13, 8))
    ax_im3.set_ylabel("E [eV]", fontsize=18)
    ax_im3.set_xlabel(r'$k_{||}$ $[\frac{1}{\mu m}$]', fontsize=18)
    ax_im3.tick_params(axis='both', width=2, length=5, which='major', labelsize=15.5)
    # set limits to not overshoot later when fitting
    ax_im3.set_ylim([IMG_ENERGIES3[-1], IMG_ENERGIES3[0]])
    ax_hop3.set_ylabel("Hopfield coefficients squared", fontsize=18)
    ax_hop3.set_xlabel(r'$k_{||}$ $[\frac{1}{\mu m}$]', fontsize=18)
    ax_hop3.tick_params(axis='both', width=2, length=7, which='major', labelsize=15.5, direction='in')
    for a in ['top', 'bottom', 'left', 'right']:
        ax_hop3.spines[a].set_linewidth(2.6)

    # display the corresp. image
    colmap = cm.get_cmap('Greens').reversed()
    im3 = ax_im3.imshow(np.rot90(IMG_DATA3[0]), interpolation='none', origin='lower',
                        extent=[k_axis[0] * 10 ** (-6), k_axis[-1] * 10 ** (-6), IMG_ENERGIES3[-1], IMG_ENERGIES3[0]],
                        aspect="auto", cmap=colmap)
    # use up as much space as possible
    fig_3.subplots_adjust(left=0.07, right=0.98, bottom=0.10, top=0.93, wspace=0.19, hspace=0.0)

    # add pretty colorbar
    divider3 = make_axes_locatable(ax_im3)
    cax3 = divider3.append_axes('right', size='5%', pad=0.05)
    fig_3.colorbar(im3, cax=cax3, orientation='vertical')
    cax3.set_ylabel("Intensity (arb. unit)", rotation=270, fontsize=18, labelpad=17.5)
    # hide ticks and labels for colorbar
    cax3.xaxis.set_major_locator(ticker.NullLocator())
    cax3.yaxis.set_major_locator(ticker.NullLocator())
    # ---------------------------------------------------------------------------------------- #
    # ------------------ FIGURE WITH CONTROLS
    axis_color = 'lightgoldenrodyellow'
    fig_ctrl = plt.figure(figsize=(7, 7))
    # create controls and set starting value to the predefined initial guesses
    # note: energies printed in eV
    steps = 1000
    # EC(0) for all plots
    slider_ec01_ax = fig_ctrl.add_axes([0.1, 0.2, 0.3, 0.05], facecolor=axis_color)
    slider_ec01 = Slider(slider_ec01_ax, r'$E_{c,1}(0)$ [eV]', EC0_1 / EL * 0.97, EC0_1 / EL * 1.03,
                         valinit=EC0_1 / EL, valfmt='%1.7f', valstep=EC0_1 / EL / 6 / steps)
    slider_ec02_ax = fig_ctrl.add_axes([0.1, 0.5, 0.3, 0.05], facecolor=axis_color)
    slider_ec02 = Slider(slider_ec02_ax, r'$E_{c,2}(0)$ [eV]', EC0_2 / EL * 0.97, EC0_2 / EL * 1.03,
                         valinit=EC0_2 / EL, valfmt='%1.7f', valstep=EC0_2 / EL / 6 / steps)
    slider_ec03_ax = fig_ctrl.add_axes([0.1, 0.8, 0.3, 0.05], facecolor=axis_color)
    slider_ec03 = Slider(slider_ec03_ax, r'$E_{c,3}(0)$ [eV]', EC0_3 / EL * 0.97, EC0_3 / EL * 1.03,
                         valinit=EC0_3 / EL, valfmt='%1.7f', valstep=EC0_3 / EL / 6 / steps)

    # EX(0)
    slider_ex0_ax = fig_ctrl.add_axes([0.575, 0.8, 0.3, 0.05], facecolor=axis_color)
    slider_ex0 = Slider(slider_ex0_ax, r'$E_x(0)$ [eV]', EX0 / EL * 0.97, EX0 / EL * 1.03,
                        valinit=EX0 / EL, valfmt='%1.7f', valstep=EX0 / EL / 6 / steps)

    # MC
    # mass slider in units of 10*(-5) electron masses
    massfactor = me * 10 ** (-5)
    slider_mc_ax = fig_ctrl.add_axes([0.575, 0.5, 0.3, 0.05], facecolor=axis_color)
    slider_mc = Slider(slider_mc_ax, r'$m_c$ [$ 10^{-5} m_e$]', MC / massfactor * 0.20,
                       MC / massfactor * 1.80,
                       valinit=MC / massfactor, valfmt='%1.4f', valstep=MC / massfactor / 160 / (steps * 5))

    # RABI
    rabifactor = hbar / (EL * 10 ** (-3))
    slider_rabi_ax = fig_ctrl.add_axes([0.575, 0.2, 0.3, 0.05], facecolor=axis_color)
    slider_rabi = Slider(slider_rabi_ax, r'$\hbar \Omega$[meV]', RABI * rabifactor * 0.75, RABI * rabifactor * 1.25,
                         valinit=RABI * rabifactor, valfmt='%1.4f', valstep=RABI * rabifactor / 50 / (steps * 5))

    # connect sliders to their handlers
    slider_ec01.on_changed(on_ec1_changed)
    slider_ec02.on_changed(on_ec2_changed)
    slider_ec03.on_changed(on_ec3_changed)
    slider_ex0.on_changed(on_ex0_changed)
    slider_mc.on_changed(on_mc_changed)
    slider_rabi.on_changed(on_rabi_changed)

    # ------------- PLOTTING -------------- #
    # ---- 1 ----- #
    # plot polariton dispersions
    ax_im1.plot(k_relevant * 10 ** (-6), E_down(k_relevant, EX0, EC0_1, MC, RABI) / EL, color='red')
    ax_im1.plot(k_relevant * 10 ** (-6), E_up(k_relevant, EX0, EC0_1, MC, RABI) / EL, color='blue')

    # plot hopfield coefficients
    ax_hop1.plot(k_axis * 10 ** (-6), hop_Cdown(k_axis, EX0, EC0_1, MC, RABI) ** 2,
                 label=r'$C_L^2(k_{||})$', color='magenta')
    ax_hop1.plot(k_axis * 10 ** (-6), hop_Xdown(k_axis, EX0, EC0_1, MC, RABI) ** 2,
                 label=r'$X_L^2(k_{||})$', color='forestgreen')
    ax_hop1.legend(markerscale=3, fontsize=17)

    # ---- 2 ----- #
    # plot polariton dispersions
    ax_im2.plot(k_relevant * 10 ** (-6), E_down(k_relevant, EX0, EC0_2, MC, RABI) / EL, color='red')
    ax_im2.plot(k_relevant * 10 ** (-6), E_up(k_relevant, EX0, EC0_2, MC, RABI) / EL, color='blue')

    # plot hopfield coefficients
    ax_hop2.plot(k_axis * 10 ** (-6), hop_Cdown(k_axis, EX0, EC0_2, MC, RABI) ** 2, label=r'$C_L^2(k_{||})$', color='magenta')
    ax_hop2.plot(k_axis * 10 ** (-6), hop_Xdown(k_axis, EX0, EC0_2, MC, RABI) ** 2, label=r'$X_L^2(k_{||})$', color='forestgreen')
    ax_hop2.legend(markerscale=3, fontsize=17)

    # ---- 3 ----- #
    # plot polariton dispersions
    ax_im3.plot(k_relevant * 10 ** (-6), E_down(k_relevant, EX0, EC0_3, MC, RABI) / EL, color='red')
    ax_im3.plot(k_relevant * 10 ** (-6), E_up(k_relevant, EX0, EC0_3, MC, RABI) / EL, color='blue')

    # plot hopfield coefficients
    ax_hop3.plot(k_axis * 10 ** (-6), hop_Cdown(k_axis, EX0, EC0_3, MC, RABI) ** 2, label=r'$C_L^2(k_{||})$', color='magenta')
    ax_hop3.plot(k_axis * 10 ** (-6), hop_Xdown(k_axis, EX0, EC0_3, MC, RABI) ** 2, label=r'$X_L^2(k_{||})$', color='forestgreen')
    ax_hop3.legend(markerscale=3, fontsize=17)

    plt.show()
