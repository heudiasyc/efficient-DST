#!/usr/bin/env python
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import csv
import re
import os
import pandas as pd
import argparse
import math

directory = "./benchmark/"
filelist = os.listdir(directory)
# families = ["random", "almost_consonant", "almost_bayesian"]
# order_relations = ["superset", "subset"]
# proportions = [0.2, 0.4, 0.6, 0.8, 1]
families = set()
order_relations = set()
proportions = set()
schemes = set()

data = []
for file in filelist:
    N = int(re.search("(?<=N-)[0-9]+", file).group(0))
    p = float(re.search("(?<=prop-)[0-9.]+", file).group(0))
    o = re.search("(?<=order-)[a-z]+", file).group(0)
    f = re.search("(?<=family-).+(?=_prop)", file).group(0)
    s = re.search("(?<=scheme-).+(?=_order)", file).group(0)
    families.add(f)
    order_relations.add(o)
    proportions.add(p)
    schemes.add(s)
    print(file)
    with open(directory + file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        z = []
        m = []
        for row in plots:
            z.append(float(row[0]))
            m.append(float(row[1]))
        z = np.array(z, dtype=float)
        m = np.array(m, dtype=float)
        data.append([N, f, p, s, o, z.mean(), z.std(), m.mean(), m.std()])

proportions = list(proportions)
proportions.sort()
df = pd.DataFrame(data, columns = ['N', 'mass_family', 'proportion', 'scheme', 'order', 't_zeta', 't_zeta_std', 't_mobius', 't_mobius_std'])

argparser = argparse.ArgumentParser()
argparser.add_argument(
    '-n',
    default=None,
    type=int,
    help='N max')
argparser.add_argument(
    '-f',
    nargs='+',
    default=None,
    type=str,
    help='Mass families')
argparser.add_argument(
    '-o',
    nargs='+',
    default=None,
    type=str,
    help='Order relations')
argparser.add_argument(
    '-p',
    nargs='+',
    default=None,
    type=float,
    help='Proportions')
argparser.add_argument(
    '-s',
    nargs='+',
    default=None,
    type=str,
    help='Schemes')

args = argparser.parse_args()

color_palette = ["red", "blue", "green", "orange", "gray", "black", "brown", "yellow"]
colors = {}
for i, f in enumerate(schemes):
    colors[f] = color_palette[i]
# fig, ax1 = plt.subplots()
families = args.f if args.f is not None else families
order_relations = args.o if args.o is not None else order_relations
proportions = args.p if args.p is not None else proportions
schemes = args.s if args.s is not None else schemes
N_max = args.n if args.n is not None else df.N.max()

for f in families:
    print(f)
    for o in order_relations:
        print(o)
        plt.clf()
        plt.xlabel('N')
        plt.ylabel('time (sec)')#, color=color)
        for p in proportions:
            print(p)
            print("Intersection entre reduced_FMT et EMT_semilattice à environ N =", 3+(4.45 + 5 / (100*p)) * math.log(100*p))
            df_e = df[
                (df.mass_family == f) & \
                (df.proportion == p) & \
                (df.order == o) & \
                (df.N <= N_max)
            ]
            for s in df_e.scheme.unique():
                if s in schemes:
                    print(s)
                    df_s = df_e[df_e.scheme == s]
                    df_s = df_s.sort_values(by=['N'])
                    plt.plot(df_s.N, df_s["t_zeta"], label="t_zeta "+s, color=colors[s], marker='^')
                    plt.plot(df_s.N, df_s["t_mobius"], label="t_mobius "+s, color=colors[s], marker='v')
                    plt.fill_between(df_s.N, (df_s.t_zeta-df_s.t_zeta_std), (df_s.t_zeta+df_s.t_zeta_std), color=colors[s], alpha=.1)
                    plt.fill_between(df_s.N, (df_s.t_mobius-df_s.t_mobius_std), (df_s.t_mobius+df_s.t_mobius_std), color=colors[s], alpha=.1)
        plt.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.legend()
        plt.title(" ".join([f, o, str(p)]))
        plt.show()
