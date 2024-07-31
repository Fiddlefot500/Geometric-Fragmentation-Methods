import random as rng
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
import math


def find_nearest(array, value):
    # Edited version of a function found on stackoverflow. Original used lists, this uses numpy arrays. Can't find the
    # original post now, will add it here if I can manage to find it.
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]


def integer_method_2D_voronoi(long, tall, fracnum, bin_num, show_img):
    # This function creates a 2D array, picks points in the 2D array at random, breaks the array apart by assigning a
    # fragment number to each point judged by which fragment point is closest, categorizes these fragments by their
    # "mass", and plots the CDF and PDF for these fragments

    print('Creating grid...')
    # Creating a discretized 2D array for the plane of fragments
    grid = np.zeros((tall, long))

    print('Getting fracture points...')
    # Choosing random points across the 2D array
    x_coords = []
    y_coords = []
    frag_size = []
    for x in range(fracnum):
        x_temp = rng.randint(0, tall-1)
        y_temp = rng.randint(0, long-1)
        while not grid[x_temp][y_temp] == 0:
            x_temp = rng.randint(0, tall-1)
            y_temp = rng.randint(0, long-1)
        grid[x_temp][y_temp] = x+1
        x_coords.append(x_temp)
        y_coords.append(y_temp)
        frag_size.append(1)

    print('Creating fragments...')
    # Assigning fragment numbers
    for x in range(tall):
        for y in range(long):
            percent = (((x*long)+y) / (long*tall))
            print(percent*100, ' % Complete...')
            if grid[x][y] == 0:
                short_dist = 0
                short_frag = 0
                for z in range(fracnum):
                    x_dist = x_coords[z] - x
                    y_dist = y_coords[z] - y
                    dist = np.sqrt((x_dist**2)+(y_dist**2))
                    if (not short_dist == 0) and (dist <= short_dist):
                        short_dist = dist
                        short_frag = z+1
                    elif short_dist == 0:
                        short_dist = dist
                        short_frag = z+1
                grid[x][y] = short_frag
                temp_index = int(grid[x][y])
                frag_size[temp_index-1] = frag_size[temp_index-1]+1
    if not sum(frag_size) == (long * tall):
        print("Uh oh, looks like you created mass! "
              "We don't particularly like that, so please either fix the code or run again")
        print("This was your final mass: ", sum(frag_size))
        print(x_coords)
        print(y_coords)
        print(frag_size)

    # These plots show the actual fragment image. Pretty cool, but unfortunately thats not all we're doing here
    if show_img:
        fig, frags = plt.subplots()
        fragment_heatmap = frags.imshow(grid)
        plt.colorbar(fragment_heatmap)

    avg_mass = (tall*long)/fracnum
    norm_frag_size = np.divide(frag_size, avg_mass)

    print('Shoving masses into bins...')
    bincats = []
    for x in range(bin_num):
        bincats.append(x * 10 / bin_num)
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

    PDF_deriv = [0]
    for x in range(len(bincats)):
        if x > 0:
            temp_deriv = (cumnormbins[x] - cumnormbins[x - 1]) / (bincats[x] - (bincats[x - 1]))
            PDF_deriv.append(temp_deriv)

    fig, df = plt.subplots()
    # df.plot(bincats,norm_bins,'ro')
    df.plot(bincats, cumnormbins, 'bx', label='CDF Data')
    # df.plot(bincats,PDF_deriv, 'g-', label='PDF using Derivative CDF data')

    binspl = np.linspace(bincats[0], bincats[-1], bin_num)
    fit = np.polyfit(bincats, cumnormbins, 10)
    CDF_func = np.poly1d(fit)
    # df.plot(binspl, CDF_func(binspl), 'k', label = 'CDF fit line')

    PDF_func = np.polyder(CDF_func)

    x_axis = np.arange(0, 5.1, 0.1)
    mu = 1
    c = 3.52440
    PDF_voronoi = (1/mu)*(4/(scp.special.gamma(4)))*(((4*x_axis)/mu)**3)*np.exp(-4*x_axis/mu)
    PDF_kiang = ((c**c)/scp.special.gamma(c))*(x_axis**(c-1))*np.exp(-c*x_axis)
    # Both equations are quite close, though Kiang's equation is slightly closer for larger numbers of fragments
    df.plot(x_axis, PDF_voronoi, 'r--', label='Theoretical PDF (Voronoi)')
    df.plot(x_axis, PDF_kiang, 'k--', label='Theoretical PDF (Kiang)')

    df.plot(binspl, PDF_func(binspl), 'm', label='PDF using CDF data')
    plt.legend()
    plt.title('Voronoi PDF and CDF')
    plt.xlabel('Normalized mass of fragments (m/m0)')
    plt.ylabel('Probability')
    plt.ylim(bottom=0)
    plt.ylim(top=1.1)
    plt.xlim(left=0)
    plt.xlim(right=5)

    kiang_dif = []
    voronoi_dif = []
    for x in range(len(x_axis)):
        kiang_dif.append(abs(PDF_kiang[x] - PDF_func(x_axis[x])))
        voronoi_dif.append(abs(PDF_voronoi[x] - PDF_func(x_axis[x])))

    fig, dif = plt.subplots()
    dif.plot(x_axis, kiang_dif, 'o--', label='Difference between Kiang Theoretical and Fit Polynomial')
    dif.plot(x_axis, voronoi_dif, 'rx-', label='Difference between Voronoi Theoretical and Fit Polynomial')
    plt.xlabel('Normalized mass of fragments')
    plt.ylabel('Absolute difference in probability')
    plt.legend()
    plt.xlim(left=0)
    plt.xlim(right=5)


length = 10**3
height = 10**3
number_of_fractures = 10**2
number_of_bins = 50
integer_method_2D_voronoi(length, height, number_of_fractures, number_of_bins, show_img=True)

plt.show()
