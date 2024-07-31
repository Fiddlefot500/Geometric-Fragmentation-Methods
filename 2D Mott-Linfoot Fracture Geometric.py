import random as rng
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
import time as t
import math

def find_nearest(array,value):
    # Edited version of a function found on stackoverflow. Original used lists, this uses numpy arrays. Can't find the
    # original post now, will add it here if I can manage to find it.
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def integer_method_2d_ml(long, tall, x_linenum, y_linenum, bin_num, show_xy, show_img):
    # This function uses the method Mott and Linfoot posed where fragments would be created by drawing lines of random
    # position through a 2D plane, though they must be either horizontal or vertical. This code uses a more "arithmatic"
    # method where only the borders of fragments used for calculating the area of each fragment

    fragnum = (x_linenum + 1) * (y_linenum + 1)

    print('Getting fracture lines...')
    # Choosing random lines across the 2D array
    x_lines = []  # lines parallel with the x-axis
    y_lines = []  # lines parallel with the y-axis
    frag_size = np.zeros(fragnum)
    frag_length = np.zeros(y_linenum+1)
    frag_height = np.zeros(x_linenum+1)
    for x in range(x_linenum):
        x_temp = rng.randint(1, tall - 1)
        while x_lines.__contains__(x_temp):
            x_temp = rng.randint(1, tall - 1)
        x_lines.append(x_temp)
    for x in range(y_linenum):
        y_temp = rng.randint(1, long - 1)
        while y_lines.__contains__(y_temp):
            y_temp = rng.randint(1, long - 1)
        y_lines.append(y_temp)

    x_lines.append(long)
    y_lines.append(tall)

    x_lines = sorted(x_lines)
    y_lines = sorted(y_lines)

    print('Drawing lines, fragmenting, and weighting...')

    for x in range(x_linenum + 1):
        if (x == 0):
            frag_height[x] = x_lines[x]
        else:
            frag_height[x] = x_lines[x] - x_lines[x - 1]
    for x in range(y_linenum + 1):
        if (x == 0):
            frag_length[x] = y_lines[x]
        else:
            frag_length[x] = y_lines[x] - y_lines[x - 1]
    frag_index = 0
    for x in range(x_linenum+1):
        print((x * 100) / (x_linenum+1), " % Complete...")
        for y in range(y_linenum+1):
            frag_size[frag_index] = frag_length[y] * frag_height[x]
            frag_index = frag_index + 1

    if not sum(frag_size) == (tall * long):
        print('You broke the law of conservation of mass. Congrats!')
        print(sum(frag_size))

    avg_mass = (tall * long) / fragnum
    avg_length = long / x_linenum
    avg_height = tall / y_linenum
    norm_frag_size = np.divide(frag_size, avg_mass)
    norm_frag_length = np.divide(frag_length, avg_length)
    norm_frag_height = np.divide(frag_height, avg_height)

    print('Shoving masses into bins...')
    bincats = []
    for x in range(bin_num):
        bincats.append(x * 5.5 / bin_num)
    bins = []
    for x in range(bin_num):
        bins.append(0)
    for x in range(len(bincats)):
        for y in range(len(norm_frag_size)):
            if x == 0:
                if norm_frag_size[y] <= bincats[x]:
                    bins[0] = bins[0] + 1
            elif bincats[x] >= norm_frag_size[y] > bincats[x - 1]:
                bins[x] = bins[x] + 1
            elif x == (len(bincats) - 1):
                if norm_frag_size[y] >= bincats[x]:
                    bins[x] = bins[x] + 1

    norm_bins = np.divide(bins, len(norm_frag_size))
    cumnormbins = np.cumsum(norm_bins)

    print("Plotting...")

    fig, df = plt.subplots()
    # df.plot(bincats,norm_bins,'ro')
    df.plot(bincats, cumnormbins, 'bx', label='CDF Data')

    binspl = np.linspace(bincats[0], bincats[-1], bin_num)
    fit = np.polyfit(bincats, cumnormbins, 10)
    CDF_func = np.poly1d(fit)
    df.plot(binspl, CDF_func(binspl), 'r')
    PDF_func = np.polyder(CDF_func)
    df.plot(binspl, PDF_func(binspl), 'm', label='PDF using CDF data')

    PDF_deriv = [0]
    for x in range(len(bincats)):
        if x > 0:
            temp_deriv = (cumnormbins[x] - cumnormbins[x - 1]) / (bincats[x] - (bincats[x - 1]))
            PDF_deriv.append(temp_deriv)

    a = np.arange(0, bincats[-1], 0.1)

    CDF_mott = 1 - sorted(2 * np.sqrt(a)) * scp.special.k1(sorted(2 * np.sqrt(a)))
    if CDF_mott[0] == np.NaN:
        np.put(CDF_mott, [0], [0])
    PDF_mott_deriv = [0]
    for x in range(len(a)):
        if x > 0:
            temp_deriv = (CDF_mott[x] - CDF_mott[x - 1]) / (a[x] - (a[x - 1]))
            PDF_mott_deriv.append(temp_deriv)
    df.plot(a, CDF_mott, 'k', label='Theoretical Mott-Linfoot CDF')
    df.plot(a, PDF_mott_deriv, 'k--', label='Theoretical Mott-Linfoot PDF')
    # df.plot(bincats[1:-1],PDF_deriv[1:-1], 'g--',label='PDF using Derivative CDF data')

    fit_diff = []
    deriv_diff = []
    for x in range(len(a)):
        fit_diff.append(abs(PDF_func(a[x]) - PDF_mott_deriv[x]))
    for x in range(len(bincats)):
        near = find_nearest(a,bincats[x])
        temp_x = np.argwhere(a == near)
        deriv_diff.append(abs(PDF_deriv[x] - PDF_mott_deriv[temp_x[0][0]]))

    plt.legend()
    plt.title('2D Mott-Linfoot PDF and CDF')
    plt.xlabel('Fragment Size')
    plt.ylabel('Probability Distribution')
    plt.xlim(right=5)
    plt.xlim(left=0)
    plt.ylim(top=1.2)
    plt.ylim(bottom=0)

    fig, dif = plt.subplots()
    dif.plot(a[1:-1], fit_diff[1:-1],'o--',label='Difference between Theoretical and Fit Polynomial')
    #dif.plot(bincats[2:-1],deriv_diff[2:-1],'rx-',label='Difference between Theoretical and Derivated CDF Data')
    plt.xlabel('Normalized mass of fragments')
    plt.ylabel('Absolute difference in probability')
    #plt.legend()
    plt.xlim(left=0)
    plt.xlim(right=5)
    plt.ylim(bottom=0)

    if show_img:
        fig, img = plt.subplots()
        for x in range(len(x_lines)):
            img.plot([0,length],[x_lines[x],x_lines[x]],'k')
        for x in range(len(y_lines)):
            img.plot([y_lines[x],y_lines[x]],[0,height],'k')
        plt.ylim(bottom=0)
        plt.ylim(top=height)
        plt.xlim(left=0)
        plt.xlim(right=length)



    if show_xy:

        ybincats = []
        for x in range(y_linenum):
            ybincats.append(x * 10 / y_linenum)
        ybins = []
        for x in range(y_linenum):
            ybins.append(0)
        for x in range(len(ybincats)):
            for y in range(len(norm_frag_height)):
                if x == 0:
                    if norm_frag_height[y] <= ybincats[x]:
                        ybins[0] = ybins[0] + 1
                elif ybincats[x] >= norm_frag_height[y] > ybincats[x - 1]:
                    ybins[x] = ybins[x] + 1
                elif x == (len(ybincats) - 1):
                    if norm_frag_height[y] >= ybincats[x]:
                        ybins[x] = ybins[x] + 1

        xbincats = []
        for x in range(x_linenum):
            xbincats.append(x * 10 / x_linenum)
        xbins = []
        for x in range(x_linenum):
            xbins.append(0)
        for x in range(len(xbincats)):
            for y in range(len(norm_frag_length)):
                if x == 0:
                    if norm_frag_length[y] <= xbincats[x]:
                        xbins[0] = xbins[0] + 1
                elif xbincats[x] >= norm_frag_length[y] > xbincats[x - 1]:
                    xbins[x] = xbins[x] + 1
                elif x == (len(xbincats) - 1):
                    if norm_frag_length[y] >= xbincats[x]:
                        xbins[x] = xbins[x] + 1

        norm_xbins = np.divide(xbins, len(norm_frag_length))
        norm_ybins = np.divide(ybins, len(norm_frag_height))
        cumnormxbins = np.cumsum(norm_xbins)
        cumnormybins = np.cumsum(norm_ybins)

        fig, xdf = plt.subplots()
        xdf.plot(xbincats, cumnormxbins, 'bx', label='CDF Data')

        xbinspl = np.linspace(xbincats[0], xbincats[-1], x_linenum)
        xfit = np.polyfit(xbincats, cumnormxbins, 10)
        xCDF_func = np.poly1d(xfit)
        xPDF_func = np.polyder(xCDF_func)
        xdf.plot(xbinspl, xCDF_func(xbinspl), 'r')
        xdf.plot(xbinspl, xPDF_func(xbinspl), 'mo--', label='PDF using CDF data')
        x_axis = np.arange(0, 5, 0.1)
        xCDF_lineau = 1 - np.exp(-x_axis)
        xPDF_lineau = [0]
        for x in range(len(x_axis)):
            if x > 0:
                temp_deriv = (xCDF_lineau[x] - xCDF_lineau[x - 1]) / (x_axis[x] - (x_axis[x - 1]))
                xPDF_lineau.append(temp_deriv)
        xdf.plot(x_axis, xCDF_lineau, 'k', label='Theoretical CDF')
        xdf.plot(x_axis[1:-1], xPDF_lineau[1:-1], 'k--', label='Theoretical PDF')

        plt.legend()
        plt.title('Length PDF and CDF')
        plt.xlabel('Fragment Length')
        plt.ylabel('Probability Distribution')
        plt.xlim(right=10)
        plt.xlim(left=0)
        plt.ylim(top=1.2)
        plt.ylim(bottom=0)

        fig, ydf = plt.subplots()
        ydf.plot(ybincats, cumnormybins, 'bx', label='CDF Data')

        ybinspl = np.linspace(ybincats[0], ybincats[-1], y_linenum)
        yfit = np.polyfit(ybincats, cumnormybins, 10)
        yCDF_func = np.poly1d(yfit)
        yPDF_func = np.polyder(yCDF_func)
        ydf.plot(ybinspl, yCDF_func(ybinspl), 'r')
        ydf.plot(ybinspl, yPDF_func(ybinspl), 'mo--', label='PDF using CDF data')
        ydf.plot(x_axis, xCDF_lineau, 'k', label='Theoretical CDF')
        ydf.plot(x_axis[1:-1], xPDF_lineau[1:-1], 'k--', label='Theoretical PDF')

        plt.legend()
        plt.title('Height PDF and CDF')
        plt.xlabel('Fragment Height')
        plt.ylabel('Probability Distribution')
        plt.xlim(right=10)
        plt.xlim(left=0)
        plt.ylim(top=1.2)
        plt.ylim(bottom=0)

    print("Number of fragments: ", fragnum)


t1 = t.time()
length = 10 ** 5
height = 10 ** 5
number_of_x_lines = 10**2
number_of_y_lines = 10**2
number_of_bins = 50
integer_method_2d_ml(length, height, number_of_x_lines, number_of_y_lines, number_of_bins, show_xy=False, show_img=True)
t2 = t.time()
print(t2 - t1, ' seconds')

plt.show()
