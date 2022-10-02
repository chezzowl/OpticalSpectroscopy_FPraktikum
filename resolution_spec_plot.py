import numpy as np
import matplotlib.pyplot as plt

# array containing with i-th entr given by: [slit_width_i , grating_i, resolution_i, resol_error_i]

data150 = np.array([
    [10, 4316.50761415362, 117.64671468348013],
    [25, 4000, 117],
    [50, 3000, 117],
    [100, 2000, 117],
    [250, 500, 117]
])

data300 = np.array([
    [10, 6800, 369],
    [25, 5863.287165995443, 369.8781341998863],
    [50, 4700, 369],
    [100, 4000, 369],
    [250, 3000, 369]
])

data1200 = np.array([
    [10, 11959.819938555414, 613.7875882603876],
    [25, 15207.513925721032, 287.2540219336781],
    [50, 11688.549398283856, 221.23777476682628],
    [100, 8000, 221],
    [250, 5000, 221]
])

# plot
fig, ax = plt.subplots()
# make plot pretty
plt.subplots_adjust(left=0.08, right=0.96, bottom=0.17, top=0.98, wspace=0.2, hspace=0.0)
ax.set_xlabel(r'slit width [$\mu$m]', fontsize=26)
ax.set_ylabel("resolution", fontsize=26)
ax.tick_params(width=2.8, length=10, top=True, right=True, direction='in', labelsize=22)
ax.tick_params(which='minor', top=True, right=True, width=1.8, length=5.5, direction='in')
for a in ['top', 'bottom', 'left', 'right']:
    ax.spines[a].set_linewidth(3.3)

ax.set_xscale("log")
ax.set_yscale("log")

# plot all slit width - resolution combis
ax.errorbar(data150[:, 0], data150[:, 1], data150[:, 2], ls='none', marker='x', capsize=3.3, ms=8, elinewidth=1.4,
            label='150 lines/mm')
ax.errorbar(data300[:, 0], data300[:, 1], data300[:, 2], ls='none', marker='x', capsize=3.3, ms=8, elinewidth=1.4,
            label='300 lines/mm')
ax.errorbar(data1200[:, 0], data1200[:, 1], data1200[:, 2], ls='none', marker='x', capsize=3.3, ms=8, elinewidth=1.4,
            label='1200 lines/mm')
ax.legend(prop={'size': 24})
plt.show()