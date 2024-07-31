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


def integer_method_2d_cylseg(long, tall, x_linenum, y_linemax, bin_num, show_img):
    # This code is based on the method Mott posed where fragments would be created by drawing random horizontal lines
    # and then drawing random vertical lines between each horizontal line. This code uses a more "arithmatic" method
    # where only the borders of fragments used for calculating the area of each fragment

    print('Getting fracture lines')
    x_lines = []
    for x in range(x_linenum):
        print((x*100)/x_linenum,' % Complete...')
        x_temp = rng.randint(1, tall - 1)
        while x_lines.__contains__(x_temp):
            x_temp = rng.randint(1, tall - 1)
        x_lines.append(x_temp)
    x_lines = sorted(x_lines)
    x_lines.append(tall)
    frag_height = np.zeros(x_linenum+1)
    for x in range(x_linenum+1):
        if (x == 0):
            frag_height[x] = x_lines[x]
        else:
            frag_height[x] = x_lines[x] - x_lines[x - 1]
    y_lines=[]
    for x in range(len(x_lines)):
        print((x * 100) / len(x_lines), ' % Complete...')
        y_linenum = round(y_linemax)
        if(y_linenum < 1):
            y_linenum = 1
        temp_y_lines = []
        for y in range(y_linenum):
            y_temp = rng.randint(1, long - 1)
            while temp_y_lines.__contains__(y_temp):
                y_temp = rng.randint(1, long - 1)
            temp_y_lines.append(y_temp)
        temp_y_lines = sorted(temp_y_lines)
        temp_y_lines.append(long)
        y_lines.append([x_lines[x],temp_y_lines])

    print('Getting fragment sizes...')
    frag_length = []
    for x in range(len(y_lines)):
        for y in range(len(y_lines[x][1])):
            if(y==0):
                frag_length.append([frag_height[x],y_lines[x][1][y]])
            else:
                frag_length.append([frag_height[x],(y_lines[x][1][y] - y_lines[x][1][y-1])])

    frag_sizes = []
    for x in range(len(frag_length)):
        frag_sizes.append(frag_length[x][0]*frag_length[x][1])

    if(not sum(frag_sizes) == tall*long):
        print('You broke the law of conservation of mass. Congrats!')
        print(sum(frag_sizes))

    print('There are ',len(frag_sizes),' fragments')

    fragnum = len(frag_sizes)
    avg_mass = (tall*long)/fragnum
    norm_frag_sizes = np.divide(frag_sizes,avg_mass)

    print('Shoving masses into bins...')
    bincats = []
    for x in range(bin_num):
        bincats.append(x * 10 / bin_num)
    bins = []
    for x in range(bin_num):
        bins.append(0)
    for x in range(len(bincats)):
        for y in range(len(norm_frag_sizes)):
            if x == 0:
                if norm_frag_sizes[y] <= bincats[x]:
                    bins[0] = bins[0] + 1
            elif bincats[x] >= norm_frag_sizes[y] > bincats[x - 1]:
                bins[x] = bins[x] + 1
            elif x == (len(bincats) - 1):
                if norm_frag_sizes[y] >= bincats[x]:
                    bins[x] = bins[x] + 1

    norm_bins = np.divide(bins, len(norm_frag_sizes))
    cumnormbins = np.cumsum(norm_bins)


    print("Plotting...")

    fig, df = plt.subplots()
    #df.plot(bincats,norm_bins,'ro')
    df.plot(bincats, cumnormbins, 'bx', label='CDF Data')

    binspl = np.linspace(bincats[0], bincats[-1], bin_num)
    fit = np.polyfit(bincats, cumnormbins, 10)
    CDF_func = np.poly1d(fit)
    #df.plot(binspl, CDF_func(binspl), 'r')
    PDF_func = np.polyder(CDF_func)
    df.plot(binspl, PDF_func(binspl), 'm-', label='PDF using CDF data')
    a = np.arange(0, bincats[-1], 0.1)

    PDF_deriv = [0]
    for x in range(len(bincats)):
        if x > 0:
            temp_deriv = (cumnormbins[x] - cumnormbins[x - 1]) / (bincats[x] - (bincats[x - 1]))
            PDF_deriv.append(temp_deriv)

    # The equations used for displaying the theoretical behavior are technically not the associated equations found by
    # Mott. However, they are so similar that they are practically the same behavior.
    CDF_mott = 1 - sorted(2 * np.sqrt(a)) * scp.special.k1(sorted(2 * np.sqrt(a)))
    if CDF_mott[0] == np.NaN:
        np.put(CDF_mott, [0], [0])
    PDF_mott = [0]
    for x in range(len(a)):
        if x > 0:
            temp_deriv = (CDF_mott[x] - CDF_mott[x - 1]) / (a[x] - (a[x - 1]))
            PDF_mott.append(temp_deriv)
    df.plot(a, CDF_mott, 'k-', label='Theoretical Mott-Linfoot CDF')
    df.plot(a, PDF_mott, 'k--', label='Theoretical Mott-Linfoot PDF')
    df.plot(bincats[1:-1],PDF_deriv[1:-1], 'g--',label='PDF using Derivative CDF data')

    plt.legend()
    plt.title('2D Cylindrical Segmentation PDF and CDF')
    plt.xlabel('Fragment Size')
    plt.ylabel('Probability Distribution')
    plt.xlim(right=5)
    plt.xlim(left=0)
    plt.ylim(top=1.2)
    plt.ylim(bottom=0)

    fit_diff = []
    deriv_diff = []
    for x in range(len(a)):
        fit_diff.append(abs(PDF_func(a[x]) - PDF_mott[x]))
    for x in range(len(bincats)):
        near = find_nearest(a, bincats[x])
        temp_x = np.argwhere(a == near)
        deriv_diff.append(abs(PDF_deriv[x] - PDF_mott[temp_x[0][0]]))

    fig, dif = plt.subplots()
    dif.plot(a[1:-1], fit_diff[1:-1], 'o--', label='Difference between Theoretical and Fit Polynomial')
    dif.plot(bincats[2:-1], deriv_diff[2:-1], 'rx-', label='Difference between Theoretical and Derivated CDF Data')
    plt.xlabel('Normalized mass of fragments')
    plt.ylabel('Absolute difference in probability')
    plt.legend()
    plt.xlim(left=0)
    plt.xlim(right=5)

    if(show_img):
        fig, img = plt.subplots()
        for x in range(len(x_lines)):
            img.plot([0, length], [x_lines[x], x_lines[x]], 'k')
            for y in range(len(y_lines)):
                if(x_lines[x] < tall):
                    img.plot([y_lines[x][1],y_lines[x][1]],[x_lines[x+1],x_lines[x]],'k')
        plt.xlim(left=0)
        plt.xlim(right=length)
        plt.ylim(bottom=0)
        plt.ylim(top=height)

t1 = t.time()
length = 10**5
height = 10**5
number_of_x_lines = 10**1
max_number_of_y_lines = 10**1

number_of_bins = 50
integer_method_2d_cylseg(length, height, number_of_x_lines, max_number_of_y_lines, number_of_bins,show_img=True)
t2 = t.time()
print(t2 - t1, ' seconds')

plt.show()