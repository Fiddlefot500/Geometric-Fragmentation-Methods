import random as rng
import matplotlib.pyplot as plt
import numpy as np
import math


def find_nearest(array, value):
    # Edited version of a function found on stackoverflow. Original used lists, this uses numpy arrays. Can't find the
    # original post now, will add it here if I can manage to find it.
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]


def integer_method(sample, size, fracture_number, bin_num):
    # This method picks random integer numbers between 0 and the length of the bar as the fracture points. At large
    # enough values for the bar length, this method resembles the real method

    print("Selecting fracture points...")

    # Random break distance with set number of fractures
    method = []
    for x in range(sample):
        fracture_points4 = []
        for y in range(fracture_number):
            temp_val = rng.randint(0, size - 1)
            fracture_points4.append(temp_val)
        fracture_points4 = sorted(fracture_points4)
        method.extend([fracture_points4])

    print("Catching standard values...")

    # Getting a list of the number of masses in each run for every method, deleting any possible zero masses
    mass_num = []
    for x in range(len(method)):
        mass_num.append(len(method[x])+1)

    # Getting a list of every mass
    masses = []
    zero = []
    for x in range(sample):
        zero.append(0)
        for y in range(mass_num[x]):
            if y == mass_num[x] - 1:
                new_mass = (size - 1) - method[x][y - 1]
                if new_mass > 0:
                    masses.append(new_mass)
                else:
                    zero[x] = zero[x]+1
            elif y == 0:
                new_mass = method[x][0]
                if new_mass > 0:
                    masses.append(new_mass)
                else:
                    zero[x] = zero[x]+1
            else:
                new_mass = method[x][y] - method[x][y - 1]
                if new_mass > 0:
                    masses.append(new_mass)
                else:
                    zero[x] = zero[x]+1
    mass_num = np.subtract(mass_num, zero)

    avg_mass = np.divide((size - 1), mass_num)

    # normalizing the list of masses with the expected average mass
    norm_masses = []

    for x in range(sample):
        for y in range(mass_num[x]):
            norm_masses.append(np.divide(masses[y], avg_mass[x]))

    print("Shoving masses into bins...")

    # Creating bins to categorize every normalized mass to count how many masses are within a certain size
    bincats = []
    for x in range(bin_num):
        bincats.append(x * 10 / bin_num)
    bins = []
    for x in range(bin_num):
        bins.append(0)
    for x in range(len(bincats)):
        for y in range(len(norm_masses)):
            if x == 0:
                if norm_masses[y] <= bincats[x]:
                    bins[0] = bins[0] + 1
            elif bincats[x] >= norm_masses[y] >= bincats[x - 1]:
                bins[x] = bins[x] + 1
            elif x == range(len(bincats) - 1):
                if norm_masses[y] >= bincats[x]:
                    bins[x] = bins[x] + 1

    print("Plotting...")

    norm_bins = np.divide(bins, len(norm_masses))
    cumnormbins = np.cumsum(norm_bins)

    fig, counted_hist = plt.subplots()

    # Checking the mass adds up to 1000
    mass_check_bool = False
    mass_check = []
    start_val = 0
    for x in range(sample):
        temp_sum = 0
        start_val_temp = 0
        for y in range(mass_num[x]):
            temp_sum = temp_sum + masses[start_val + y]
            start_val_temp = start_val_temp + 1
        mass_check.append(temp_sum)
        start_val = start_val + start_val_temp
    for x in range(len(mass_check)):
        if not mass_check[x] == size - 1:
            mass_check_bool = True
    if mass_check_bool:
        print("Number 4 broke the law")
        print(mass_check)

    binspl = np.linspace(bincats[0], bincats[-1], bin_num)
    fit = np.polyfit(bincats, cumnormbins, 10)
    CDF_func = np.poly1d(fit)
    # counted_hist.plot(binspl,CDF_func(binspl),'b')

    PDF_func = np.polyder(CDF_func)

    # counted_hist.plot(binspl,PDF_func(binspl),'m',label='PDF using fit CDF data')

    x_axis = np.arange(0, 5, 0.1)
    CDF_lineau = 1 - np.exp(-x_axis)
    PDF_lineau = [0]
    for x in range(len(x_axis)):
        if x > 0:
            temp_deriv = (CDF_lineau[x] - CDF_lineau[x - 1]) / (x_axis[x] - (x_axis[x - 1]))
            PDF_lineau.append(temp_deriv)
    PDF_deriv = [0]
    for x in range(len(bincats)):
        if x > 0:
            temp_deriv = (cumnormbins[x] - cumnormbins[x - 1]) / (bincats[x] - (bincats[x - 1]))
            PDF_deriv.append(temp_deriv)

    PDF_avg = []
    for x in range(len(bincats)):
        temp_x = find_nearest(binspl, bincats[x])
        PDF_avg.append((PDF_func(temp_x) + PDF_deriv[x]) / 2)

    # counted_hist.plot(x_axis, CDF_lineau, 'k--',label='Theoretical CDF')
    counted_hist.plot(x_axis[1:-1], PDF_lineau[1:-1], 'k--', label='Theoretical PDF')
    # counted_hist.plot(bincats[1:-1],PDF_deriv[1:-1], 'g--',label='PDF of Derivative CDF Data')
    counted_hist.plot(bincats[1:-1], PDF_avg[1:-1], 'r', label='PDF using average of fit and derivative CDF data')
    plt.ylim(top=1.2)
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    plt.xlim(right=5)
    plt.title('Integer PDF and CDF')
    plt.legend()
    plt.xlabel('Normalized mass of fragments (m/m0)')
    plt.ylabel('Probability')

    avg_diff = []
    for x in range(len(bincats)):
        near = find_nearest(x_axis, bincats[x])
        temp_x = np.argwhere(x_axis == near)
        avg_diff.append(abs(PDF_avg[x] - PDF_lineau[temp_x[0][0]]))

    fig, dif = plt.subplots()
    dif.plot(bincats[1:-1], avg_diff[1:-1], label='Absolute difference between Theoretical PDF and Averaged PDF data')
    plt.xlabel('Normalized mass of fragments')
    plt.ylabel('Absolute difference in probability')
    plt.title('Integer Method Absolute Difference between Theoretical PDF and Averaged PDF data')
    plt.xlim(left=0)
    plt.xlim(right=5)


def real_method(sample, size, fracnum, bin_num):
    # This method picks random real numbers between 0 and the length of the bar as the fracture points.

    print("Selecting fracture points...")
    method = []
    for x in range(sample):
        frac_points = []
        for y in range(fracnum):
            temp_val = rng.random()*size
            frac_points.append(temp_val)
        frac_points = sorted(frac_points)
        method.extend([frac_points])

    print("Catching standard values...")

    # Getting a list of the number of masses
    mass_num = []
    for x in range(len(method)):
        mass_num.append(len(method[x])+1)

    # Getting a list of every mass
    masses = []
    zero = []
    for x in range(len(method)):
        zero.append(0)
    for x in range(sample):
        for y in range(mass_num[x]):
            if y == mass_num[x] - 1:
                new_mass = (size - 1) - method[x][y - 1]
                if new_mass > 0:
                    masses.append(new_mass)
                else:
                    zero[x] = zero[x] + 1
            elif y == 0:
                new_mass = method[x][0]
                if new_mass > 0:
                    masses.append(new_mass)
                else:
                    zero[x] = zero[x] + 1
            else:
                new_mass = method[x][y] - method[x][y - 1]
                if new_mass > 0:
                    masses.append(new_mass)
                else:
                    zero[x] = zero[x] + 1

    avg_mass = np.divide((size-1), mass_num)

    # Normalizing the masses
    norm_masses = []
    for x in range(sample):
        for y in range(mass_num[x]):
            norm_masses.append(np.divide(masses[y], avg_mass[x]))

    print("Shoving masses into bins...")

    # Creating bins to categorize the normalized masses by "weight"
    bincats = []
    for x in range(bin_num):
        bincats.append(x * 10/bin_num)
    bins = []
    for x in range(bin_num):
        bins.append(0)
    for x in range(len(bincats)):
        for y in range(len(norm_masses)):
            if x == 0:
                if norm_masses[y] <= bincats[x]:
                    bins[0] = bins[0] + 1
            elif bincats[x] >= norm_masses[y] >= bincats[x - 1]:
                bins[x] = bins[x] + 1
            elif x == range(len(bincats) - 1):
                if norm_masses[y] >= bincats[x]:
                    bins[x] = bins[x] + 1

    print("Plotting...")

    norm_bins = np.divide(bins, len(norm_masses))
    cumnormbins = np.cumsum(norm_bins)

    fig, counted_hist = plt.subplots()

    # Checking the mass adds up to 1000
    mass_check_bool = False
    mass_check = []
    start_val = 0
    for x in range(sample):
        temp_sum = 0
        start_val_temp = 0
        for y in range(mass_num[x]):
            temp_sum = temp_sum + masses[start_val + y]
            start_val_temp = start_val_temp + 1
        mass_check.append(temp_sum)
        start_val = start_val + start_val_temp
    for x in range(len(mass_check)):
        if not mass_check[x] == size - 1:
            mass_check_bool = True
    if mass_check_bool:
        print("Real numbers broke the law")
        print(mass_check)

    binspl = np.linspace(bincats[0], bincats[-1], bin_num)
    fit = np.polyfit(bincats, cumnormbins, 10)
    CDF_func = np.poly1d(fit)
    # counted_hist.plot(binspl, CDF_func(binspl), 'b')

    PDF_func = np.polyder(CDF_func)

    x_axis = np.arange(0, 5.1, 0.1)
    CDF_lineau = 1 - np.exp(-x_axis)
    PDF_lineau = [0]
    for x in range(len(x_axis)):
        if x > 0:
            temp_deriv = (CDF_lineau[x] - CDF_lineau[x - 1]) / (x_axis[x] - (x_axis[x - 1]))
            PDF_lineau.append(temp_deriv)
    PDF_deriv = [0]
    for x in range(len(bincats)):
        if x > 0:
            temp_deriv = (cumnormbins[x] - cumnormbins[x - 1]) / (bincats[x] - (bincats[x - 1]))
            PDF_deriv.append(temp_deriv)

    PDF_avg = []
    for x in range(len(bincats)):
        temp_x = find_nearest(binspl, bincats[x])
        PDF_avg.append((PDF_func(temp_x) + PDF_deriv[x])/2)

    avg_diff = []
    for x in range(len(bincats)):
        near = find_nearest(x_axis, bincats[x])
        temp_x = np.argwhere(x_axis == near)
        avg_diff.append(abs(PDF_avg[x] - PDF_lineau[temp_x[0][0]]))

    # counted_hist.plot(x_axis,CDF_lineau,'k',label='Theoretical CDF')
    counted_hist.plot(x_axis[1:-1], PDF_lineau[1:-1], 'k--', label='Theoretical PDF')

    # counted_hist.plot(binspl,PDF_func(binspl),'m',label='PDF using fit CDF data')
    # counted_hist.plot(bincats[1:-1],PDF_deriv[1:-1], 'g--',label='PDF using Derivative of CDF data')
    counted_hist.plot(bincats[1:-1], PDF_avg[1:-1], 'r', label='PDF using average of fit and derivative CDF data')
    plt.ylim(top=1.2)
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    plt.xlim(right=5)
    plt.legend()
    plt.title('Real Number PDF and CDF')
    plt.xlabel('Normalized mass of fragments (m/m0)')
    plt.ylabel('Probability')

    fig, dif = plt.subplots()
    dif.plot(bincats[1:-1], avg_diff[1:-1], label='Absolute difference between Theoretical PDF and Averaged PDF data')
    plt.xlabel('Normalized mass of fragments')
    plt.ylabel('Absolute difference in probability')
    plt.title('Real Method Absolute Difference between Theoretical PDF and Averaged PDF data')
    plt.xlim(left=0)
    plt.xlim(right=5)


sample_size = 1  # only really need this to be equal to 1, but you could have it larger I guess.
length = 10**9 + 1  # heavily recommended to be larger than the number of fractures
number_of_fractures = 10**4  # recommended to be between 10**3 and the effective length
number_of_bins = 50  # recommend to be smaller the larger the length and number of fractures

integer_method(sample_size, length, number_of_fractures, number_of_bins)
real_method(sample_size, length, number_of_fractures, number_of_bins)

plt.show()
