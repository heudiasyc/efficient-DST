#!/usr/bin/env python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import csv
import re
import os

families = ["random", "almost_consonant", "almost_bayesian"]

plot_FMT = True
# plot_FMT = False
support_size = 26
family = families[2]
mobius = True
# mobius = False

directory = "../DST_experiments/"
filelist = os.listdir(directory)
colors = ["tab:red", "tab:blue", "tab:green", "tab:orange", "tab:gray", "tab:yellow", "tab:black", "tab:brown"]
fig, ax1 = plt.subplots()

if mobius:
    plt.title("Möbius transform\nSupport size = "+ str(support_size) +"\n"+family+" sampling")
    row_offset = 1
    base_label = "Möbius with "
else:
    plt.title("Zeta transform\nSupport size = "+ str(support_size) +"\n"+family+" sampling")
    row_offset = 0
    base_label = "Zeta with "

ax1.set_xlabel('N')
ax1.set_ylabel('time (sec)')#, color=color)

ax2 = ax1.twinx()
ax2.set_ylabel("Average scheme support size")

color_index = 0
max_support_size = 0.
max_time = 0.
x = []
yF = []
for file in filelist:
    if re.search("^.*support-"+str(support_size)+"_family-"+family+".*.csv$", file) is not None:
        print(file)
        with open(directory + file, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter=',')
            x = []
            yE = []
            yF = []
            avg_support_sizes = []
            for row in plots:
                x.append(int(row[0]))
                yF.append(float(row[row_offset+1]))
                yE.append(float(row[row_offset+3]))
                avg_support_sizes.append(float(row[5]))
            max_support_size = max(max_support_size, max(avg_support_sizes))
            max_time = max(max_time, max(yE))
            if re.search("^.*_scheme-lattice.*.csv$", file) is not None:
                label = "lattice scheme"
            elif re.search("^.*_scheme-semilattice.*.csv$", file) is not None:
                label = "semilattice scheme"
            else:
                label = "direct scheme"
            ax1.plot(x,yE, label=base_label+label, color=colors[color_index])
            ax2.plot(x, avg_support_sizes, linestyle="--", label="Average scheme support size", color=colors[color_index])
            color_index += 1

if len(x) > 0:
    if plot_FMT:
        # max_time = max(max_time, max(yF))
        avg_support_sizes = []
        for i in range(len(x)):
            avg_support_sizes.append(pow(2, x[i]))
        ax1.plot(x,yF, label=base_label+"FMT", color=colors[color_index])
        ax2.plot(x, avg_support_sizes, linestyle="--", label="Average scheme support size", color=colors[color_index])

    ax1.set_ylim([0,max_time])
    ax2.set_ylim([0,int(max_support_size+0.5)])
    ax1.set_xlim([min(x), max(x)])
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    ax1.legend()
    ax2.legend(loc='center left')
    plt.show()
