#!/usr/bin/python

import re
import os
import sys

color_mapping = {'ad': 'red', 'ik': 'blue', 'ik-i': 'deepskyblue', 'ad-i': 'firebrick', 'ad-r': 'brown', 'ik-r': 'cyan'}

def build_plot_string(dataset, color_mapping):
    if dataset['method'] in color_mapping:
        color_string = "lc rgb '%s'" % color_mapping[dataset['method']]
    else:
        color_string = "lc rgb 'black'"
    return "'%s' w lp t '%s' %s" % (dataset['filename'],dataset['method'],color_string)

file_re = re.compile('p3m-(.*?)-([0-9]*)-to-([0-9]*)-prec-(.*)-rcut-(.*)-dens-(.*)-q-(.*)\.dat');

params = ['method', 'lownum', 'highnum', 'prec', 'rcut', 'dens', 'q']

datasets = []
precs = []

for i in os.listdir(sys.argv[1]):
    match = file_re.match(i)
    if not match is None:
        seq = zip(params, match.groups())
        seq.append(('filename', i))
        datasets.append(dict(seq))
        if not datasets[-1]['prec'] in precs:
            precs.append(datasets[-1]['prec']);
            print datasets[-1]['prec']

plot_strings = []    

for i in datasets:
    plot_strings.append(build_plot_string(i, color_mapping))
    print i

f = open('all.gnuplot', 'w')


f.write("plot %s" %(", \\\n".join(plot_strings)))








