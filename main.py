import sys

sys.path.insert(0, '..')

import matplotlib.pyplot as plt
import numpy as np
import sif_reader

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # filepath = "C:\\Users\\ivayl\\Documents\\Uni\\F Praktikum\\OS\\Angelov_Strunck\\B_150_10_2s.sif"
    filepath = "C:\\Users\\ivayl\\Documents\\Uni\\F Praktikum\\OS\\Angelov_Strunck\\Emission_largenegative.sif"
    data, info = sif_reader.np_open(filepath)
    data = np.asarray(data, dtype=float)
    wavelengths = sif_reader.utils.extract_calibration(info)
    # (data, info) = sif_reader.utils.parse(filepath)
    # place data into a pandas Series
    # df = pd.Series(data[:, 1], index=data[:, 0])

    print(data[0])
    print(data[0][1])
    print(np.shape(data))
    # print(len(data), "   ", len(data[0]))
    print(info)
    # testdata = [wavelengths, data[2]]

    plt.subplot()
    plt.imshow(data[0], interpolation='none', origin='lower',
               extent=[wavelengths[0], wavelengths[-1], 0, 1024], aspect="auto")
    # plt.imshow(data[0], interpolation='none', extent=[wavelengths[0], wavelengths[-1], 1024, 0], aspect="auto")
    # plt.imshow(testdata)

    plt.show()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
