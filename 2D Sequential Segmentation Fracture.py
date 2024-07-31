import random as rng
import matplotlib.pyplot as plt
import numpy as np
import time as t
import math

# This code is far from perfect with trying to replicate Grady and Kipp's Sequential Segmentation method. It sometimes
# gives segments large panhandles, likely creating differences between the data and the theoretical curves. This code
# is also incredibly slow at determining the segments, which limits its capability a non-insignificant amount. The
# highest I would recommend the settings to go is 10^3 on the length, height, and number of fragments, and at that scale
# the plotted results do not seem granular enough to get usable data.


def find_nearest(array, value):
    # Edited version of a function found on stackoverflow. Original used lists, this uses numpy arrays. Can't find the
    # original post now, will add it here if I can manage to find it.
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]


def integer_method_2d_seqseg(long, tall, fragnum, bin_num, show_img):
    print('Creating grid...')
    # Creating a discretized 2D array for the plane of fragments

    grid = []
    frag_grid = []
    for x in range(tall):
        row = []
        frag_row = np.zeros(long)
        for y in range(long):
            el = [0]
            row.append(el)
        grid.extend([row])
        frag_grid.extend([frag_row])

    # Adding the outer points to the available start points for lines

    avail_points = []
    for x in range(tall):
        for y in range(long):
            if x == 0:
                avail_points.append([x, y, 'v', 'd'])
            elif y == 0:
                avail_points.append([x, y, 'h', 'r'])
    avail_points.remove([0, 0, 'v', 'd'])
    avail_points.remove([0, long-1, 'v', 'd'])

    current_fragnum = 1
    print('Getting fracture points...')
    while not current_fragnum == fragnum:
        print((current_fragnum*100)/fragnum, ' % Complete...')
        point = avail_points[rng.randint(0, len(avail_points)-1)]
        x_start = point[1]
        x_pos = point[1]
        y_start = point[0]
        y_pos = point[0]
        temp_frag_values = grid[y_pos][x_pos]
        hit_x_wall = False
        hit_y_wall = False
        avail_points.remove([y_pos, x_pos, point[2], point[3]])
        if point[2] == 'h':
            # Horizontal line
            if point[3] == 'l':
                # Go left
                while not hit_y_wall:
                    while np.array_equal(grid[y_pos][x_pos], temp_frag_values):
                        grid[y_pos][x_pos] = np.append(grid[y_pos][x_pos], current_fragnum)
                        if y_pos == y_start:
                            if avail_points.__contains__([y_pos, x_pos, 'v', 'd']):
                                avail_points.remove([y_pos, x_pos, 'v', 'd'])
                            if avail_points.__contains__([y_pos, x_pos, 'v', 'u']):
                                avail_points.remove([y_pos, x_pos, 'v', 'u'])
                            if avail_points.__contains__([y_pos, x_pos, 'h', 'r']):
                                avail_points.remove([y_pos, x_pos, 'h', 'r'])
                            if avail_points.__contains__([y_pos, x_pos, 'h', 'l']):
                                avail_points.remove([y_pos, x_pos, 'h', 'l'])
                            if y_pos+1 < tall:
                                if not (avail_points.__contains__(
                                        [y_pos+1, x_pos, 'v', 'd']) or avail_points.__contains__(
                                        [y_pos+1, x_pos, 'v', 'u']) or avail_points.__contains__(
                                        [y_pos+1, x_pos, 'h', 'r']) or avail_points.__contains__(
                                        [y_pos+1, x_pos, 'h', 'l'])):
                                    avail_points.append([y_pos+1, x_pos, 'v', 'd'])
                            if y_pos-1 > 0:
                                if not (avail_points.__contains__(
                                        [y_pos-1, x_pos, 'v', 'd']) or avail_points.__contains__(
                                        [y_pos-1, x_pos, 'v', 'u']) or avail_points.__contains__(
                                        [y_pos-1, x_pos, 'h', 'r']) or avail_points.__contains__(
                                        [y_pos-1, x_pos, 'h', 'l'])):
                                    avail_points.append([y_pos-1, x_pos, 'v', 'u'])
                        x_pos = x_pos - 1
                        if not (x_pos >= 0):
                            break
                    if y_pos < tall-1:
                        y_pos = y_pos + 1
                        x_pos = x_start
                    else:
                        hit_y_wall = True
                current_fragnum = current_fragnum + 1
            else:
                # Go right
                while not hit_y_wall:
                    while np.array_equal(grid[y_pos][x_pos], temp_frag_values):
                        grid[y_pos][x_pos] = np.append(grid[y_pos][x_pos], current_fragnum)
                        if y_pos == y_start:
                            if avail_points.__contains__([y_pos, x_pos, 'v', 'd']):
                                avail_points.remove([y_pos, x_pos, 'v', 'd'])
                            if avail_points.__contains__([y_pos, x_pos, 'v', 'u']):
                                avail_points.remove([y_pos, x_pos, 'v', 'u'])
                            if avail_points.__contains__([y_pos, x_pos, 'h', 'r']):
                                avail_points.remove([y_pos, x_pos, 'h', 'r'])
                            if avail_points.__contains__([y_pos, x_pos, 'h', 'l']):
                                avail_points.remove([y_pos, x_pos, 'h', 'l'])
                            if y_pos + 1 < tall:
                                if not (avail_points.__contains__(
                                        [y_pos, x_pos, 'v', 'd']) or avail_points.__contains__(
                                        [y_pos, x_pos, 'v', 'u']) or avail_points.__contains__(
                                        [y_pos, x_pos, 'h', 'r']) or avail_points.__contains__(
                                        [y_pos, x_pos, 'h', 'l'])):
                                    avail_points.append([y_pos, x_pos, 'v', 'd'])
                            if y_pos - 1 > 0:
                                if not (avail_points.__contains__(
                                        [y_pos-1, x_pos, 'v', 'd']) or avail_points.__contains__(
                                        [y_pos-1, x_pos, 'v', 'u']) or avail_points.__contains__(
                                        [y_pos-1, x_pos, 'h', 'r']) or avail_points.__contains__(
                                        [y_pos-1, x_pos, 'h', 'l'])):
                                    avail_points.append([y_pos - 1, x_pos, 'v', 'u'])
                        x_pos = x_pos + 1
                        if not (x_pos <= long - 1):
                            break

                    if y_pos < tall - 1:
                        y_pos = y_pos + 1
                        x_pos = x_start
                    else:
                        hit_y_wall = True
                current_fragnum = current_fragnum + 1
        else:
            # Vertical line
            if point[3] == 'u':
                # Go up
                while not hit_x_wall:
                    while np.array_equal(grid[y_pos][x_pos], temp_frag_values):
                        grid[y_pos][x_pos] = np.append(grid[y_pos][x_pos], current_fragnum)
                        if x_pos == x_start:
                            if avail_points.__contains__([y_pos, x_pos, 'v', 'd']):
                                avail_points.remove([y_pos, x_pos, 'v', 'd'])
                            if avail_points.__contains__([y_pos, x_pos, 'v', 'u']):
                                avail_points.remove([y_pos, x_pos, 'v', 'u'])
                            if avail_points.__contains__([y_pos, x_pos, 'h', 'r']):
                                avail_points.remove([y_pos, x_pos, 'h', 'r'])
                            if avail_points.__contains__([y_pos, x_pos, 'h', 'l']):
                                avail_points.remove([y_pos, x_pos, 'h', 'l'])
                            if x_pos + 1 < long:
                                if not (avail_points.__contains__(
                                        [y_pos, x_pos, 'v', 'd']) or avail_points.__contains__(
                                        [y_pos, x_pos, 'v', 'u']) or avail_points.__contains__(
                                        [y_pos, x_pos, 'h', 'r']) or avail_points.__contains__(
                                        [y_pos, x_pos, 'h', 'l'])):
                                    avail_points.append([y_pos, x_pos, 'h', 'r'])
                            if x_pos - 1 > 0:
                                if not (avail_points.__contains__(
                                        [y_pos, x_pos-1, 'v', 'd']) or avail_points.__contains__(
                                        [y_pos, x_pos-1, 'v', 'u']) or avail_points.__contains__(
                                        [y_pos, x_pos-1, 'h', 'r']) or avail_points.__contains__(
                                        [y_pos, x_pos-1, 'h', 'l'])):
                                    avail_points.append([y_pos, x_pos-1, 'h', 'l'])
                        y_pos = y_pos - 1
                        if not (y_pos >= 0):
                            break

                    if x_pos < long - 1:
                        x_pos = x_pos + 1
                        y_pos = y_start
                    else:
                        hit_x_wall = True
                current_fragnum = current_fragnum + 1
            else:
                # Go down
                while not hit_x_wall:
                    while np.array_equal(grid[y_pos][x_pos], temp_frag_values):
                        grid[y_pos][x_pos] = np.append(grid[y_pos][x_pos], current_fragnum)
                        if x_pos == x_start:
                            if avail_points.__contains__([y_pos, x_pos, 'v', 'd']):
                                avail_points.remove([y_pos, x_pos, 'v', 'd'])
                            if avail_points.__contains__([y_pos, x_pos, 'v', 'u']):
                                avail_points.remove([y_pos, x_pos, 'v', 'u'])
                            if avail_points.__contains__([y_pos, x_pos, 'h', 'r']):
                                avail_points.remove([y_pos, x_pos, 'h', 'r'])
                            if avail_points.__contains__([y_pos, x_pos, 'h', 'l']):
                                avail_points.remove([y_pos, x_pos, 'h', 'l'])
                            if x_pos + 1 < long:
                                if not (avail_points.__contains__(
                                        [y_pos, x_pos, 'v', 'd']) or avail_points.__contains__(
                                        [y_pos, x_pos, 'v', 'u']) or avail_points.__contains__(
                                        [y_pos, x_pos, 'h', 'r']) or avail_points.__contains__(
                                        [y_pos, x_pos, 'h', 'l'])):
                                    avail_points.append([y_pos, x_pos, 'h', 'r'])
                            if x_pos - 1 > 0:
                                if not (avail_points.__contains__(
                                        [y_pos, x_pos-1, 'v', 'd']) or avail_points.__contains__(
                                        [y_pos, x_pos-1, 'v', 'u']) or avail_points.__contains__(
                                        [y_pos, x_pos-1, 'h', 'r']) or avail_points.__contains__(
                                        [y_pos, x_pos-1, 'h', 'l'])):
                                    avail_points.append([y_pos, x_pos - 1, 'h', 'l'])
                        y_pos = y_pos + 1
                        if not (y_pos <= tall - 1):
                            break

                    if x_pos < long - 1:
                        x_pos = x_pos + 1
                        y_pos = y_start
                    else:
                        hit_x_wall = True
                current_fragnum = current_fragnum + 1

    for x in range(tall):
        for y in range(long):
            frag_grid[x][y] = grid[x][y][-1]

    frag_size = np.zeros(fragnum)
    for x in range(tall):
        for y in range(long):
            frag_size[int(frag_grid[x][y])] = frag_size[int(frag_grid[x][y])] + 1

    if not sum(frag_size) == tall * long:
        print('You broke the law of conservation of mass. Congrats!')
        print(sum(frag_size))

    if show_img:
        fig, frags = plt.subplots()
        fragment_heatmap = frags.imshow(frag_grid)
        plt.colorbar(fragment_heatmap)

    avg_mass = (tall * long) / fragnum
    norm_frag_sizes = np.divide(frag_size, avg_mass)

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
    # df.plot(bincats,norm_bins,'ro')
    df.plot(bincats, cumnormbins, 'bx', label='CDF Data')

    binspl = np.linspace(bincats[0], bincats[-1], bin_num)
    fit = np.polyfit(bincats, cumnormbins, 10)
    CDF_func = np.poly1d(fit)
    # df.plot(binspl, CDF_func(binspl), 'r')
    PDF_func = np.polyder(CDF_func)
    df.plot(binspl, PDF_func(binspl), 'm-', label='PDF using CDF data')

    PDF_deriv = [0]
    for x in range(len(bincats)):
        if x > 0:
            temp_deriv = (cumnormbins[x] - cumnormbins[x - 1]) / (bincats[x] - (bincats[x - 1]))
            PDF_deriv.append(temp_deriv)

    x_axis = np.arange(0, 5, 0.1)
    CDF_GK = 1 - np.exp(-x_axis)
    PDF_GK = [0]
    for x in range(len(x_axis)):
        if x > 0:
            temp_deriv = (CDF_GK[x] - CDF_GK[x - 1]) / (x_axis[x] - (x_axis[x - 1]))
            PDF_GK.append(temp_deriv)
    df.plot(x_axis, CDF_GK, 'k', label='Theoretical CDF')
    df.plot(x_axis[1:-1], PDF_GK[1:-1], 'k--', label='Theoretical PDF')
    df.plot(bincats[1:-1], PDF_deriv[1:-1], 'g--', label='PDF using Derivative CDF data')

    plt.legend()
    plt.title('2D Sequential Segmentation PDF and CDF')
    plt.xlabel('Fragment Size')
    plt.ylabel('Probability Distribution')
    plt.xlim(right=5)
    plt.xlim(left=0)
    plt.ylim(top=1.2)
    plt.ylim(bottom=0)

    fit_diff = []
    deriv_diff = []
    for x in range(len(x_axis)):
        fit_diff.append(abs(PDF_func(x_axis[x]) - PDF_GK[x]))
    for x in range(len(bincats)):
        near = find_nearest(x_axis, bincats[x])
        temp_x = np.argwhere(x_axis == near)
        deriv_diff.append(abs(PDF_deriv[x] - PDF_GK[temp_x[0][0]]))

    fig, dif = plt.subplots()
    dif.plot(x_axis[1:-1], fit_diff[1:-1], 'o--', label='Difference between Theoretical and Fit Polynomial')
    dif.plot(bincats[2:-1], deriv_diff[2:-1], 'rx-', label='Difference between Theoretical and Derivated CDF Data')
    plt.xlabel('Normalized mass of fragments')
    plt.ylabel('Absolute difference in probability')
    plt.legend()
    plt.xlim(left=0)
    plt.xlim(right=5)


t1 = t.time()
length = 10**3
height = 10**3
number_of_fragments = 10**1

number_of_bins = 50
integer_method_2d_seqseg(length, height, number_of_fragments, number_of_bins, show_img=True)
t2 = t.time()
print(t2 - t1, ' seconds')

plt.show()
