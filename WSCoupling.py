import sys

sys.path.insert(0, '..')
import sif_reader
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import tkinter as tk
from tkinter.filedialog import askopenfilename
from matplotlib.widgets import Slider, TextBox, Button
from scipy.optimize import curve_fit
import numbers

# global variables
hbar = 6.582119569 * 10 ** (-16)  # [eV*s]
h = 4.135667696 * 10 ** (-15)  # [eV*s] .... 6.62607015×10−34 in [J*s]
c = 299792458  # [m/s]


# --------------------------- LOAD & STORE DATA ------------------------------- #
def import_dat_data():
    # not the best style, but a global image storage will do the job for now
    global IMG_DATA
    global IMG_INFO
    global IMG_WAVELENGTHS
    global IMG_ENERGIES
    global DAT_FILE_PATH
    try:
        root = tk.Tk(className='Import sif-file')
        DAT_FILE_PATH = askopenfilename()
        IMG_DATA, IMG_INFO = sif_reader.np_open(DAT_FILE_PATH)
        IMG_DATA = np.asarray(IMG_DATA, dtype=float)
        # num_img = np.shape(IMG_DATA)[0]
        # for i in range(num_img):
        #     IMG_DATA[i] = IMG_DATA[i]/np.max(IMG_DATA[i])
        IMG_WAVELENGTHS = sif_reader.utils.extract_calibration(IMG_INFO)
        IMG_ENERGIES = h * c * 10 ** 9 / IMG_WAVELENGTHS  # c*10**9 --> everything in nanometers
        root.destroy()
    except:
        print('Wrong file')
        root.destroy()


# --------------------------- FITTING RELATED ------------------------------- #
def gaussian(x, A, x0, sig):
    """

    :param x: input variable
    :param A: amplitude
    :param x0:
    :param sig: standard deviation
    :return:
    """
    return A / (sig * np.sqrt(2 * np.pi)) * np.exp(-(x - x0) ** 2 / (2 * sig ** 2))


def lorentz_single(x, A, x0, sig):
    return 2 * A / np.pi * sig / (4 * (x - x0) ** 2 + sig ** 2)


def lorentz_multi(x, *beta):
    """

    :param x:
    :param beta: fit parameters for variable number of peaks in the following form:
    [A1, A2,..., x01, x01, ..., sig1, sig2, ..., y0]
    :return:
    """
    num = int((len(beta) - 1) // 3)

    A = beta[0:num]
    x0 = beta[num:2 * num]
    sig = beta[2 * num:3 * num]
    y0 = beta[3 * num]

    result = np.zeros(x.shape[0])

    for o in range(num):
        result = result + lorentz_single(x, A[o], x0[o], sig[o])

    return result + y0


# ------------------ GUI FUNCTIONALITIES ----------------------- #
def on_button_singleLor_pressed(event):
    # global peak_legendline
    global leg_loc
    global lorentz_peaks
    # find array indices where to start and end the fitting procedure
    idx_left = min(range(len(IMG_WAVELENGTHS)), key=lambda i: abs(IMG_WAVELENGTHS[i] - slider_xmin.val))
    idx_right = min(range(len(IMG_WAVELENGTHS)), key=lambda i: abs(IMG_WAVELENGTHS[i] - slider_xmax.val))
    print("left index:  ", idx_left, "  right index:  ", idx_right)
    # read in input values fo lorentz fit function ...
    fit_init = [currentA1, currentx01, currentsigma1]
    print("FIT INIT: ", fit_init)
    # ... and perform lorentz fit and plot if they are valid
    if isinstance(fit_init[0], numbers.Number) and isinstance(fit_init[1], numbers.Number) and isinstance(fit_init[2],
                                                                                                          numbers.Number):
        pars, covs = curve_fit(lorentz_single, IMG_WAVELENGTHS[idx_left:idx_right],
                               IMG_DATA[0][CURRENTSLICE][idx_left:idx_right], p0=fit_init)

        # plot only within fitting range .. looks better this way
        col = cmap_lorentz(np.random.uniform())  # random color for gauss peak
        lorentz_lines.append(mlines.Line2D([], [], color=col, markersize=15, label='Fitted Lorentzian'))
        # plot over input range but with smaller step size to really see the Lorentz fit
        plotarray = np.linspace(IMG_WAVELENGTHS[idx_left], IMG_WAVELENGTHS[idx_right], 10000)
        ax_slice.plot(plotarray, lorentz_single(plotarray, *pars), '--', color=col,
                      label="lorentz fit", linewidth=2.1)

        # add corresponding text box and line
        ax_slice.axvline(linewidth=1, color=color_peaks, x=pars[1], ls='--', alpha=0.6)
        nums_lines_added.append(2)  # 1 for the plot, 1 for the dashed peak line

        # ax_slice.text(pars[1] + IMG_WAVELENGTHS[-1] * 0.005, lorentz_single(pars[1], *pars),
        #               r"$x_0={:.2f}\:$nm".format(pars[1], pars[2]), fontsize=8.5)
        # at last: add peak information to the predefined array
        lorentz_peaks.append([pars, covs])
    fig_2D.canvas.draw()
    fig_slice.canvas.draw()


def on_button_doubleLor_pressed(event):
    # global peak_legendline
    global leg_loc
    global lorentz_peaks
    # find array indices where to start and end the fitting procedure
    idx_left = min(range(len(IMG_WAVELENGTHS)), key=lambda i: abs(IMG_WAVELENGTHS[i] - slider_xmin.val))
    idx_right = min(range(len(IMG_WAVELENGTHS)), key=lambda i: abs(IMG_WAVELENGTHS[i] - slider_xmax.val))
    # read in input values fo lorentz fit function ...
    fit_init = [currentA1, currentA2, currentx01, currentx02, currentsigma1, currentsigma2, currenty0]
    # fit_init = [currentA1, currentA2, currentx01, currentx02, currentsigma1, currentsigma2]
    # print(fit_init)
    # ... and perform lorentz fit and plot if they are valid
    if all(isinstance(x, numbers.Number) for x in fit_init):
        pars, covs = curve_fit(lorentz_multi, IMG_WAVELENGTHS[idx_left:idx_right],
                               IMG_DATA[0][CURRENTSLICE][idx_left:idx_right], p0=fit_init)
        # pars, covs = curve_fit(lambda x, a1, x01, sig1, a2, x02, sig2:
        #                        lorentz_multi(x, a1, x01, sig1, a2, x02, sig2, currenty0),
        #                        IMG_WAVELENGTHS[idx_left:idx_right], IMG_DATA[0][CURRENTSLICE][idx_left:idx_right],
        #                        p0=fit_init)
        # pars = np.append(pars, currenty0)
        # print(pars)
        # plot only within fitting range .. looks better this way
        col = cmap_lorentz(np.random.uniform())  # random color for gauss peak
        lorentz_lines.append(mlines.Line2D([], [], color=col, markersize=15, label='Fitted Lorentzian'))
        # plot over input range but with smaller step size to really see the Lorentz fit
        plotarray = np.linspace(IMG_WAVELENGTHS[idx_left], IMG_WAVELENGTHS[idx_right], 10000)
        ax_slice.plot(plotarray, lorentz_multi(plotarray, *pars), '--', color=col,
                      label="lorentz fit", linewidth=2.1)
        # add text boxes and lines next to peaks TODO: CONTINUE HERE
        ax_slice.axvline(linewidth=1, color=color_peaks, x=pars[2], ls='--', alpha=0.6)
        ax_slice.axvline(linewidth=1, color=color_peaks, x=pars[3], ls='--', alpha=0.6)
        # keep track of number of lines and texts added
        nums_lines_added.append(3)  # 1 for the plot, 2 for the dashed peak lines
        # ax_slice.text(pars[1] + IMG_WAVELENGTHS[-1] * 0.005, lorentz_single(pars[1], *pars),
        #               r"$x_0={:.2f}\:$nm".format(pars[1], pars[2]), fontsize=8.5)
        # at last: add peak information to the predefined array
        lorentz_peaks.append([pars, covs])
        # for now only double peaks allowed
        print_doublepeaks_info(lorentz_peaks)
    fig_2D.canvas.draw()
    fig_slice.canvas.draw()


def on_del_pressed(event):
    # delete last element...
    if len(ax_slice.lines) > 1:  # make sure the original slice is never deleted
        lines_to_pop = nums_lines_added[-1]  # num of lines painted for the currently last fit
        for i in range(lines_to_pop):
            ax_slice.lines.pop()
        # now remove actual lorentz fit from legend ...
        lorentz_lines.pop()
        # ... and the peak information from the array storage ...
        lorentz_peaks.pop()
        # and the now deprecated line number information
        nums_lines_added.pop()
        # same thing for peak texts
        # if len(ax_slice.texts) > 0:
        #     texts_to_pop = nums_texts_added[-1]
        # for i in range(texts_to_pop):
        #     ax_slice.texts.pop()
        # nums_texts_added.pop()  #TODO: reactivate
        # for now only double peaks allowed
        print_doublepeaks_info(lorentz_peaks)
    fig_2D.canvas.draw()
    fig_slice.canvas.draw()


def on_sliders_changed(event):
    global lorentzrange
    # first reset background color, otherwise the new part gets painted on top
    lorentzrange.remove()
    lorentzrange = ax_slice.axvspan(slider_xmin.val, slider_xmax.val, facecolor='#2ca02c', alpha=0.2)
    fig_2D.canvas.draw()
    fig_slice.canvas.draw()


def on_slice_changed(event):
    global CURRENTSLICE
    CURRENTSLICE = int(slider_slice.val)
    # axis should always have at least one line since we add one immediately after initialization
    ax_2D.lines.pop()
    # ax_2D.axhline(linewidth=1, color=color_slice, y=CURRENTSLICE, ls='--', alpha=0.8)
    ax_2D.axvline(linewidth=1, color=color_slice, x=k_axis[CURRENTSLICE], ls='--', alpha=0.8)
    # plot current slice
    ax_slice.clear()
    ax_slice.plot(IMG_WAVELENGTHS, IMG_DATA[0][CURRENTSLICE], 'x')
    fig_2D.canvas.draw()
    fig_slice.canvas.draw()


def on_tboxA1_changed(expr):
    global currentA1
    currentA1 = float(expr)


def on_tboxx01_changed(expr):
    global currentx01
    currentx01 = float(expr)


def on_tboxsigma1_changed(expr):
    global currentsigma1
    currentsigma1 = float(expr)


def on_tboxA2_changed(expr):
    global currentA2
    currentA2 = float(expr)


def on_tboxx02_changed(expr):
    global currentx02
    currentx02 = float(expr)


def on_tboxsigma2_changed(expr):
    global currentsigma2
    currentsigma2 = float(expr)


def on_tboxy0_changed(expr):
    global currenty0
    currenty0 = float(expr)


def print_doublepeaks_info(peaks):
    strout = ""
    for entry in peaks:
        A1, A2, x1, x2, s1, s2, y0 = entry[0]
        covmat = entry[1]
        # errors are square roots of the diagonal elements
        x1_err = np.sqrt(covmat[2][2])
        s1_err = np.sqrt(covmat[4][4])
        x2_err = np.sqrt(covmat[3][3])
        s2_err = np.sqrt(covmat[5][5])
        # -------------- Q factor ---------------- #
        factor = 2 / np.sqrt(2)  # factor in Q-factor calculation
        # calc quality factor
        R1 = x1 / (factor * np.abs(s1))
        R2 = x2 / (factor * np.abs(s2))
        R1_err = R1 * np.sqrt((x1_err / x1) ** 2 + ((factor * s1_err) / (factor * s1)) ** 2)
        R2_err = R2 * np.sqrt((x2_err / x2) ** 2 + ((factor * s2_err) / (factor * s2)) ** 2)
        # ----------------- resonance energy -------------------- #
        E1 = h * c / (x1 * 10 ** (-9))
        # calc E-error by getting the maximum and minimum and then the bigger difference as the error
        E1_up = h * c / ((x1 - x1_err) * 10 ** (-9))
        E1_down = h * c / ((x1 + x1_err) * 10 ** (-9))
        del_E1 = max(abs(E1_up - E1), abs(E1 - E1_down))  # abs not necessary if we are sure what we do
        E2 = h * c / (x2 * 10 ** (-9))
        # calc E-error by getting the maximum and minimum and then the bigger difference as the error
        E2_up = h * c / ((x2 - x2_err) * 10 ** (-9))
        E2_down = h * c / ((x2 + x2_err) * 10 ** (-9))
        del_E2 = max(abs(E2_up - E2), abs(E2 - E2_down))  # abs not necessary if we are sure what we do

        # ---------------- print results ---------------- #
        strout += "\n ({:.5f},  {:.5f},  {:.5f},  {:.5f},  {:.5f},  {:.5f},  {:.5f})".format(A1, x1, s1, A2, x2, s2, y0)
        strout += "\npeak 1:"
        strout += "\nvalue: {}".format(lorentz_multi(np.array([x1]), A1, A2, x1, x2, s1, s2, y0))
        strout += "\nx01 = {} +- {}       sig01 = {} +- {}".format(x1, x1_err, s1, s1_err)
        strout += "\nE1= {} +- {}".format(E1, del_E1)
        strout += "\npeak 2 value: {}".format(lorentz_multi(np.array([x2]), A1, A2, x1, x2, s1, s2, y0))
        strout += "\nx02 = {} +- {}       sig02 = {} +- {}".format(x2, x2_err, s2, s2_err)
        strout += "\nE2= {} +- {}".format(E2, del_E2)
        strout += ("\n\n --->" + "R1 = ({} , {})").format(R1, R1_err)
        strout += ("\n --->" + "R2 = ({} , {})").format(R2, R2_err)
        print("Current Peaks:" + strout)


def print_peaks_info(peaks):
    strout = ""
    linecounter = 0
    for line in peaks:
        linecounter += 1
        num_peaks = int((len(line) - 1) // 3)
        # extract peak infos ...
        As = line[0:num_peaks]
        x0s = line[num_peaks:2 * num_peaks]
        sigs = line[2 * num_peaks:3 * num_peaks]
        y0s = line[3 * num_peaks]
        strout += "\n" + str(linecounter)
        for p in range(num_peaks):
            strout += "\n ({:.5f}  ,  {:.5f},  {:.5f},  {:.5f}),".format(As[num_peaks], x0s[num_peaks], sigs[num_peaks],
                                                                         y0s[num_peaks])
    print("Current Peaks:" + strout)


if __name__ == '__main__':
    # ------------- LOAD DATA ------------- #

    # get file
    import_dat_data()

    # ------------- init some variables ------------- #
    CURRENTSLICE = 494  # horizontal slice of current image

    currentA1 = " "
    currentx01 = " "
    currentsigma1 = " "

    currentA2 = " "
    currentx02 = " "
    currentsigma2 = " "

    currenty0 = " "

    cmap_lorentz = cm.get_cmap('hsv')
    color_peaks = 'chocolate'
    color_slice = 'red'

    # storing the lines for gauss fit legends
    lorentz_lines = []

    # storing the peak information for output
    lorentz_peaks = []

    # store in chronological order how many lines/texts were added during all performed fits
    # that way, when "Delete Last Fit" is pressed, we can track how many lines/texts have to be popped from
    # the 'Axes' object
    nums_lines_added = []
    nums_texts_added = []

    # further presets for plotting
    peak_legendline = mlines.Line2D([], [], color=color_peaks, markersize=15,
                                    label='Fitted Gaussian', linestyle='--')  # legend entry for Gauss peaks
    # TODO: set this according to current plot! ------------------- #
    leg_loc = 'upper center'
    # options: best, upper right, upper left, lower right, lower left, right, center left, center right, lower center
    #   upper center, center
    # -------------------------------------------------------------- #

    # -------------------------- GRAPHICS INIT ------------------------ #

    # ------------ figure with slice

    axis_color = 'lightgoldenrodyellow'
    fig_slice, ax_slice = plt.subplots(figsize=(18, 12))

    ax_slice.set_xlabel(r'$\lambda$[nm]', fontsize=26)
    ax_slice.set_ylabel("intensity (arb.unit)", fontsize=26)
    ax_slice.tick_params(width=2.8, length=7, top=True, right=True, direction='in', labelsize=22)
    for a in ['top', 'bottom', 'left', 'right']:
        ax_slice.spines[a].set_linewidth(3.3)

    # initial plot
    ax_slice.plot(IMG_WAVELENGTHS, IMG_DATA[0][CURRENTSLICE], '--')

    # change all spines
    for a in ['top', 'bottom', 'left', 'right']:
        ax_slice.spines[a].set_linewidth(2)

    # Adjust the subplots region to leave some space for the controls
    # fig_slice.subplots_adjust(bottom=0.2)
    fig_slice.subplots_adjust(left=0.08, right=0.96, bottom=0.17, top=0.98, wspace=0.2, hspace=0.0)

    # add buttons to the free area
    button_lor_single_ax = fig_slice.add_axes([0.04, 0.07, 0.1, 0.02], facecolor=axis_color)
    button_lor_double_ax = fig_slice.add_axes([0.04, 0.04, 0.1, 0.02], facecolor=axis_color)
    button_delete_ax = fig_slice.add_axes([0.04, 0.01, 0.1, 0.02], facecolor=axis_color)

    # ---- button to add single lorentz fit
    button_lor_single = Button(button_lor_single_ax, 'Fit Single Lorentz')
    button_lor_single.on_clicked(on_button_singleLor_pressed)

    # ---- button to add double lorentz fit
    button_lor_double = Button(button_lor_double_ax, 'Fit Double Lorentz')
    button_lor_double.on_clicked(on_button_doubleLor_pressed)

    # ---- button to delete last fit
    button_delete = Button(button_delete_ax, 'Delete Last Fit')
    button_delete.on_clicked(on_del_pressed)

    # --- sliders for fitting range
    slider_slice_ax = fig_slice.add_axes([0.2, 0.07, 0.4, 0.02], facecolor=axis_color)
    slider_slice = Slider(slider_slice_ax, 'slice', 0, np.shape(IMG_DATA)[2],
                          valinit=CURRENTSLICE, valfmt='%d', valstep=1)
    slider_xmin_ax = fig_slice.add_axes([0.2, 0.04, 0.4, 0.02], facecolor=axis_color)
    slider_xmin = Slider(slider_xmin_ax, '$x_{min}$', IMG_WAVELENGTHS[0], IMG_WAVELENGTHS[-1],
                         valinit=833, valfmt='%1.2f')
    slider_xmax_ax = fig_slice.add_axes([0.2, 0.01, 0.4, 0.02], facecolor=axis_color)
    slider_xmax = Slider(slider_xmax_ax, '$x_{max}$', IMG_WAVELENGTHS[0], IMG_WAVELENGTHS[-1],
                         valinit=839, valfmt='%1.2f')

    # connect sliders to handler
    slider_xmin.on_changed(on_sliders_changed)
    slider_xmax.on_changed(on_sliders_changed)
    slider_slice.on_changed(on_slice_changed)

    # initial fit-range coloring
    lorentzrange = ax_slice.axvspan(slider_xmin.val, slider_xmax.val, facecolor='#2ca02c', alpha=0.2)

    # --- text boxes for fitting parameters
    textbox_x01_ax = fig_slice.add_axes([0.7, 0.07, 0.05, 0.02], facecolor=axis_color)
    textbox_x01 = TextBox(textbox_x01_ax, r'$x_{0,1}$')
    textbox_A1_ax = fig_slice.add_axes([0.7, 0.04, 0.05, 0.02], facecolor=axis_color)
    textbox_A1 = TextBox(textbox_A1_ax, r"$A_1$")
    textbox_sigma1_ax = fig_slice.add_axes([0.7, 0.01, 0.05, 0.02], facecolor=axis_color)
    textbox_sigma1 = TextBox(textbox_sigma1_ax, r'$\sigma_1$')

    textbox_x02_ax = fig_slice.add_axes([0.8, 0.07, 0.05, 0.02], facecolor=axis_color)
    textbox_x02 = TextBox(textbox_x02_ax, r'$x_{0,2}$')
    textbox_A2_ax = fig_slice.add_axes([0.8, 0.04, 0.05, 0.02], facecolor=axis_color)
    textbox_A2 = TextBox(textbox_A2_ax, r"$A_2$")
    textbox_sigma2_ax = fig_slice.add_axes([0.8, 0.01, 0.05, 0.02], facecolor=axis_color)
    textbox_sigma2 = TextBox(textbox_sigma2_ax, r'$\sigma_2$')

    textbox_y0_ax = fig_slice.add_axes([0.9, 0.07, 0.05, 0.02], facecolor=axis_color)
    textbox_y0 = TextBox(textbox_y0_ax, r'$y_0$')

    # connect all boxes to handlers
    textbox_A1.on_submit(on_tboxA1_changed)
    textbox_x01.on_submit(on_tboxx01_changed)
    textbox_sigma1.on_submit(on_tboxsigma1_changed)
    textbox_A2.on_submit(on_tboxA2_changed)
    textbox_x02.on_submit(on_tboxx02_changed)
    textbox_sigma2.on_submit(on_tboxsigma2_changed)
    textbox_y0.on_submit(on_tboxy0_changed)

    # and change their value, triggering the handlers for the first time
    textbox_A1.set_val("1000")
    textbox_x01.set_val("835")
    textbox_sigma1.set_val("-1")
    textbox_A2.set_val("1000")
    textbox_x02.set_val("841")
    textbox_sigma2.set_val("-1")
    textbox_y0.set_val("50")

    # ------------ original 2D image
    # colorbar stuff taken from https://stackoverflow.com/questions/32462881/add-colorbar-to-existing-axis
    fig_2D, ax_2D = plt.subplots(figsize=(13, 8))
    # adjustments so it looks good
    ax_2D.set_ylabel("E [eV]", fontsize=24)
    ax_2D.set_xlabel(r'k [$\frac{1}{\mu m}$]', fontsize=24)
    fig_2D.subplots_adjust(left=0.10, right=0.92, bottom=0.11, top=0.98, wspace=0.2, hspace=0.0)
    ax_2D.tick_params(axis='both', width=2, length=5, which='major', labelsize=18)
    # interpolate to k-axis and plot
    lambda_c = 838 * 10 ** (-9)  # [m]
    NA = 0.42
    idx_k0 = 494
    idx_kmin = 294
    idx_kmax = 692
    k_max = 2 * np.pi * NA / (lambda_c * 10 ** 6)  # micrometers
    # get proper step size for k-array
    k_step = k_max / max(abs(idx_kmax - idx_k0), abs(idx_k0 - idx_kmin))
    k_axis = np.linspace(-idx_k0 * k_step, (1023 - idx_k0) * k_step, 1024)
    # im = ax_2D.imshow(np.rot90(IMG_DATA[0]), interpolation='none', origin='lower',
    #                   extent=[0, 1024, IMG_ENERGIES[-1], IMG_ENERGIES[0]], aspect="auto")
    im = ax_2D.imshow(np.rot90(IMG_DATA[0]), interpolation='none', origin='lower',
                      extent=[k_axis[0], k_axis[-1], IMG_ENERGIES[-1], IMG_ENERGIES[0]], aspect="auto")
    # add pretty color bar
    divider = make_axes_locatable(ax_2D)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig_2D.colorbar(im, cax=cax, orientation='vertical')
    cax.set_ylabel("Intensity (arb. unit)", rotation=270, fontsize=24, labelpad=26)
    # hide ticks and labels for colorbar
    cax.xaxis.set_major_locator(ticker.NullLocator())
    cax.yaxis.set_major_locator(ticker.NullLocator())

    # ax_2D.axhline(linewidth=1, color=color_slice, y=CURRENTSLICE, ls='--')
    ax_2D.axvline(linewidth=1, color=color_slice, x=k_axis[CURRENTSLICE], ls='--')

    plt.show()
