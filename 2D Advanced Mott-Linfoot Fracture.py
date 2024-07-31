import random as rng
import matplotlib.pyplot as plt
import numpy as np
import time as t
import math
import fractions


def find_nearest(array, value):
    # Edited version of a function found on stackoverflow. Original used lists, this uses numpy arrays. Can't find the
    # original post now, will add it here if I can manage to find it.
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]


def find_smallest_index(array):
    k = array.argmin()
    ncol = array.shape[1]
    return math.floor(k / ncol), k % ncol


def integer_method_2d_ml2(long, tall, linenum, bin_num, show_img):
    # This function uses the method Mott and Linfoot posed where fragments would be created by drawing lines of random
    # orientation and position through a 2D plane. This code uses a more "geometric" process using a 2D array to denote
    # which points are part of which fragment.

    max_fragnum = ((linenum**2) + linenum + 2)/2
    x_axis = np.arange(1, long+1)

    grid = np.ndarray.tolist(np.zeros((tall, long, 1)))
    grid2 = np.ndarray.tolist(np.zeros((tall, long)))

    pts_x = []
    pts_y = []
    print('Getting fracture lines...')
    # Choosing random lines across the 2D array
    for x in range(linenum):
        print((x*100)/linenum, ' % Complete...')
        # Choosing points
        end_pt1 = [0, 0]
        end_pt2 = [0, 0]
        face1 = str(rng.choices(['top', 'bottom', 'left', 'right'])[0])
        face2 = str(rng.choices(['top', 'bottom', 'left', 'right'])[0])
        while face1 == face2:
            face2 = str(rng.choices(['top', 'bottom', 'left', 'right'])[0])
        if face1.__eq__('top'):
            end_pt1[1] = 0
        if face1.__eq__('bottom'):
            end_pt1[1] = tall
        if face1.__eq__('left'):
            end_pt1[0] = 0
        if face1.__eq__('right'):
            end_pt1[0] = long
        if face2.__eq__('top'):
            end_pt2[1] = 0
        if face2.__eq__('bottom'):
            end_pt2[1] = tall
        if face2.__eq__('left'):
            end_pt2[0] = 0
        if face2.__eq__('right'):
            end_pt2[0] = long
        if face1.__eq__('top') or face1.__eq__('bottom'):
            end_pt1[0] = rng.randint(0, long-1)
        if face1.__eq__('left') or face1.__eq__('right'):
            end_pt1[1] = rng.randint(0, tall-1)
        if face2.__eq__('top') or face2.__eq__('bottom'):
            end_pt2[0] = rng.randint(0, long - 1)
        if face2.__eq__('left') or face2.__eq__('right'):
            end_pt2[1] = rng.randint(0, tall - 1)

        if end_pt2[0] < end_pt1[0]:
            # Going left to right, so renaming the points
            temp_pt = end_pt1
            end_pt1 = end_pt2
            end_pt2 = temp_pt
        # Finding linear relationship between the two points
        x_pts = [end_pt1[0], end_pt2[0]]
        y_pts = [end_pt1[1], end_pt2[1]]
        pts_x.append(x_pts)
        pts_y.append(y_pts)
        if x_pts[0] == x_pts[1]:
            # vertical line
            for y in range(tall):
                for z in range(x_pts[0], long):
                    if not (grid[y][z-1].__contains__(x + 1)):
                        grid[y][z-1].append(x + 1)
        else:
            slope = fractions.Fraction((end_pt2[1] - end_pt1[1]), (end_pt2[0] - end_pt1[0]))
            fit = np.polyfit(x_pts, y_pts, 1)
            lineq = np.poly1d(fit)
            y_vals = np.zeros(long)
            if slope.numerator > 0:
                for y in range(tall):
                    for z in range(0, end_pt1[0]):
                        if not (grid[y][z].__contains__(x + 1)):
                            grid[y][z].append(x + 1)
            if slope.numerator <= 0:
                for y in range(tall):
                    for z in range(end_pt2[0], long):
                        if not (grid[y][z].__contains__(x + 1)):
                            grid[y][z].append(x + 1)
            for y in range(long):
                y_vals[y] = np.round(lineq(x_axis[y]))
                if 0 <= y_vals[y] < tall:
                    if slope.numerator <= 0:
                        # Going (up) right:
                        for z in range(x_axis[y], long):
                            if not (grid[int(y_vals[y])][z].__contains__(x + 1)):
                                grid[int(y_vals[y])][z].append(x+1)
                # Going down
                for z in range(int(y_vals[y]), tall):
                    if 0 <= y_vals[y]:
                        if not (grid[z][x_axis[y] - 1].__contains__(x + 1)):
                            grid[z][x_axis[y] - 1].append(x + 1)

    print('Weighting fragments')
    frag_list = []
    frag_size = []
    frag_length = []
    frag_height = []
    for x in range(tall):
        print((x*100)/tall, ' % Complete...')
        for y in range(long):
            if frag_list.__contains__(grid[x][y]):
                temp_index = frag_list.index(grid[x][y])
                frag_size[temp_index] = frag_size[temp_index] + 1
                grid2[x][y] = temp_index + 1
                if (x - frag_height[temp_index][0]) > (frag_height[temp_index][1] - frag_height[temp_index][0]):
                    frag_height[temp_index][1] = x
                if (y - frag_length[temp_index][0]) > (frag_length[temp_index][1] - frag_length[temp_index][0]):
                    frag_length[temp_index][1] = y
                elif y < frag_length[temp_index][0]:
                    frag_length[temp_index][0] = y
            else:
                frag_list.append(grid[x][y])
                frag_length.append([y, y])
                frag_height.append([x, x])
                frag_size.append(1)
                grid2[x][y] = len(frag_size)

    if not sum(frag_size) == (tall * long):
        print('You broke the law of conservation of mass. Congrats!')
        print(sum(frag_size))

    frag_ar = []
    for x in range(len(frag_size)):
        frag_ar.append(((1+frag_length[x][1] - frag_length[x][0])/(1+frag_height[x][1] - frag_height[x][0])))

    frag_diam = np.zeros(len(frag_length))
    for x in range(len(frag_diam)):
        frag_diam[x] = 2*np.sqrt(frag_size[x]/np.pi)

    frag_ar = sorted(frag_ar)
    frag_diam = sorted(frag_diam)
    avg_diam = sum(frag_diam) / len(frag_diam)

    radsum_numerator = 0
    radsum_denominator = 0
    for x in range(len(frag_diam)):
        radsum_numerator = radsum_numerator + (frag_diam[x])**3
        radsum_denominator = radsum_denominator + (frag_diam[x])**2
    smd = radsum_numerator/radsum_denominator

    avg_ar = sum(frag_ar)/len(frag_ar)
    for x in range(len(frag_ar)):
        frag_ar[x] = frag_ar[x]/avg_ar
    eflength = (sum(frag_diam))/smd

    norm_frag_diam = np.divide(frag_diam, avg_diam)

    cutoff = norm_frag_diam[len(frag_diam)-math.floor(eflength)]

    print('Shoving masses into bins...')
    bincats = np.linspace(0, 5.5, bin_num)
    bins = []
    radii_bins = []
    for x in range(bin_num):
        bins.append(0)
        radii_bins.append(0)
    for x in range(len(bincats)):
        for y in range(len(frag_ar)):
            if x == 0:
                if frag_ar[y] <= bincats[x]:
                    bins[0] = bins[0] + 1
                if norm_frag_diam[y] <= bincats[x]:
                    if cutoff >= norm_frag_diam[y]:
                        radii_bins[0] = radii_bins[0] + 1  # norm_frag_diam[y]/cutoff
                    else:
                        radii_bins[0] = radii_bins[0] + 1
            if bincats[x] >= frag_ar[y] > bincats[x - 1]:
                bins[x] = bins[x] + 1
            if bincats[x] >= norm_frag_diam[y] > bincats[x-1]:
                radii_bins[x] = radii_bins[x] + 1
            if x == (len(bincats) - 1):
                if frag_ar[y] >= bincats[x]:
                    bins[x] = bins[x] + 1
                if norm_frag_diam[y] >= bincats[x]:
                    if cutoff >= norm_frag_diam[y]:
                        radii_bins[x] = radii_bins[x] + 1  # norm_frag_diam[y]/cutoff
                    else:
                        radii_bins[x] = radii_bins[x] + 1

    norm_bins = np.divide(bins, len(frag_ar))
    cumnormbins = np.cumsum(norm_bins)

    norm_diam_bins = np.divide(radii_bins, smd)
    cumnormdiambins = np.cumsum(norm_diam_bins)

    print("Plotting...")

    fig, bars = plt.subplots()
    bars.bar(bincats, radii_bins, width=0.075)
    bars.plot([1, 1], [0, len(frag_diam)], 'c', label='Average Diameter Size')
    bars.plot([smd / avg_diam, smd / avg_diam], [0, len(frag_diam)], 'c--', label='Sauter Mean Diameter Size')
    plt.legend()
    plt.ylim(top=max(radii_bins))

    fig, df = plt.subplots()
    df.plot(bincats, cumnormbins, 'bx', label='CDF Aspect Ratio Data')
    # df.plot(bincats,cumnormdiambins, 'ro', label='CDF Diameter Data')
    # df.plot(bincats,norm_diam_bins,'gx')

    PDF_deriv = [0]
    for x in range(len(bincats)):
        if x > 0:
            temp_deriv = (cumnormbins[x] - cumnormbins[x - 1]) / (bincats[x] - (bincats[x - 1]))
            PDF_deriv.append(temp_deriv)
    df.plot(bincats, PDF_deriv, 'g-', label='PDF using Derivative CDF data')

    binspl = np.linspace(bincats[0], bincats[-1], bin_num)
    fit = np.polyfit(bincats, cumnormbins, 15)
    CDF_func = np.poly1d(fit)
    df.plot(binspl, CDF_func(binspl), 'r', label='Fit line for CDF data')
    PDF_func = np.polyder(CDF_func)
    df.plot(binspl, PDF_func(binspl), 'm', label='PDF using CDF Aspect Ratio data')

    fit2 = np.polyfit(bincats, cumnormdiambins, 20)
    CDF_func2 = np.poly1d(fit2)
    PDF_func2 = np.polyder(CDF_func2)

    kiang_axis = np.arange(0, bincats[-1], 0.1)

    k = np.arange(1, 3, 0.1)
    l = np.arange(0.5, 1.5, 0.1)
    PDF_weib = np.zeros((len(k), len(l), len(kiang_axis)))
    CDF_weib = np.zeros((len(k), len(l), len(kiang_axis)))
    weib_dif = np.zeros((len(k), len(l), len(kiang_axis)))
    avg_weib_dif = np.zeros((len(k), len(l)))
    for x in range(len(k)):
        for y in range(len(l)):
            for z in range(len(kiang_axis)):
                PDF_weib[x, y, z] = ((k[x] / l[y]) * ((kiang_axis[z] / l[y]) ** (k[x] - 1)) *
                                     np.exp(-((kiang_axis[z] / l[y]) ** k[x])))
                CDF_weib[x, y, z] = 1 - np.exp(-((kiang_axis[z] / l[y]) ** k[x]))
                weib_dif[x, y, z] = abs(PDF_weib[x, y, z] - PDF_func(kiang_axis[z]))
            avg_weib_dif[x, y] = sum(weib_dif[x, y])/len(weib_dif[x, y])

    close_k_pos, close_l_pos = find_smallest_index(avg_weib_dif)
    close_k = k[close_k_pos]
    close_l = l[close_l_pos]
    print('The Weibull distribution (Aspect Ratio) for this size is using k=', close_k, ' and lambda =', close_l)
    df.plot(kiang_axis, PDF_weib[close_k_pos, close_l_pos], 'k--',
            label='Closest Weibull Distribution for Aspect Ratio')

    PDF_weib2 = np.zeros((len(k), len(l), len(kiang_axis)))
    CDF_weib2 = np.zeros((len(k), len(l), len(kiang_axis)))
    weib_dif2 = np.zeros((len(k), len(l), len(kiang_axis)))
    avg_weib_dif2 = np.zeros((len(k), len(l)))
    for x in range(len(k)):
        for y in range(len(l)):
            for z in range(len(kiang_axis)):
                PDF_weib2[x, y, z] = (k[x] / l[y]) * ((kiang_axis[z] / l[y]) ** (k[x] - 1)) * np.exp(-((kiang_axis[z] /
                                                                                                        l[y]) ** k[x]))
                CDF_weib2[x, y, z] = 1 - np.exp(-((kiang_axis[z] / l[y]) ** k[x]))
                weib_dif2[x, y, z] = abs(PDF_weib2[x, y, z] - PDF_func2(kiang_axis[z]))
            avg_weib_dif2[x, y] = sum(weib_dif2[x, y]) / len(weib_dif2[x, y])

    close_k_pos2, close_l_pos2 = find_smallest_index(avg_weib_dif2)
    close_k2 = k[close_k_pos2]
    close_l2 = l[close_l_pos2]
    print('The Weibull distribution (Diameter) for this size is using k=', close_k2, ' and lambda =', close_l2)

    max_density_diam = max(PDF_func2(binspl))
    max_density_ar = max(PDF_func(binspl))
    max_density_diam_size = binspl[np.where(PDF_func2(binspl) == max_density_diam)]
    max_density_ar_size = binspl[np.where(PDF_func(binspl) == max_density_ar)]
    print('The most dense diameter size is ', max_density_diam_size)
    print('The most dense aspect ratio size is ', max_density_ar_size)

    plt.legend()
    plt.title('2D Mott-Linfoot PDF and CDF')
    plt.xlabel('Fragment Size')
    plt.ylabel('Probability Distribution')
    plt.xlim(right=5)
    plt.xlim(left=0)
    plt.ylim(top=1.2)
    plt.ylim(bottom=0)

    print('The number of fragments is', len(frag_size), ', with an effective fragment number of ', eflength)
    print('Dust particles are considered to be ', cutoff, ' or smaller in size')
    print('Average diameter: ', avg_diam, ' Sauter mean diameter: ', smd)
    print('Normalized, the Sauter mean diamter is ', smd/avg_diam)

    if len(frag_size) > max_fragnum:
        print("You broke the Lazy caterer's sequence. Congrats!")

    # Visualization of fragments
    if show_img:
        randomized = np.random.permutation(len(frag_diam))
        for x in range(tall):
            for y in range(long):
                grid2[x][y] = randomized[grid2[x][y]-1]
        fig, frags = plt.subplots()
        fragment_heatmap = frags.imshow(grid2)
        plt.colorbar(fragment_heatmap)
        # Uncomment the following line to see all of the chosen starting and ending points for the fracture lines
        # frags.plot(pts_y,pts_x,'o','orange')


t1 = t.time()
length = 10**3
height = 10**3
number_of_lines = 150
number_of_bins = 50
integer_method_2d_ml2(length, height, number_of_lines, number_of_bins, show_img=True)
t2 = t.time()
print(t2 - t1, ' seconds')

plt.show()
