#!/usr/bin/python

import re
import os
import sys

file_re = re.compile('p3m-(.*?)-([0-9]*)-to-([0-9]*)-prec-(.*)-rcut-(.*)-dens-(.*)-q-(.*)\.dat');

params = ['method', 'lownum', 'highnum', 'prec', 'rcut', 'dens', 'q']

color_mapping = {'ad': 'red', 'ik': 'blue', 'ik-i': 'deepskyblue', 'ad-i': 'firebrick', 'ad-r': 'brown', 'ik-r': 'cyan'}

def build_plot_string(dataset, color_mapping):
    if dataset['method'] in color_mapping:
        color_string = "lc rgb '%s'" % color_mapping[dataset['method']]
    else:
        color_string = "lc rgb 'black'"
    return "'%s' w lp t '%s' %s" % (dataset['filename'],dataset['method'],color_string)

def write_plot_file(filename, plot_strings):
    f = open("%s.gnuplot" % filename, 'w')
    f.write("set key left top\n")
    f.write("set terminal postscript enhanced color\n")
    f.write("set output '%s.ps'\n" % filename) 
    f.write("plot %s" %(", \\\n".join(plot_strings)))
    f.close()

datasets = []
precs = {}

for i in os.listdir(sys.argv[1]):
    match = file_re.match(i)
    if not match is None:
        seq = zip(params, match.groups())
        seq.append(('filename', i))
        datasets.append(dict(seq))
        ds = datasets[-1]
        if not ds['prec'] in precs.keys():
            precs[ds['prec']] = {}
        if not ds['rcut'] in precs[ds['prec']].keys():
            precs[ds['prec']][ds['rcut']] = {}
        if not ds['dens'] in precs[ds['prec']][ds['rcut']].keys():
            precs[ds['prec']][ds['rcut']][ds['dens']] = []        
        precs[ds['prec']][ds['rcut']][ds['dens']].append(ds)
            

for prec,v in precs.items():
    for rcut,j in v.items():
        for dens,v in j.items():
            plot_strings = []
            for ds in v:
                plot_strings.append(build_plot_string(ds, color_mapping))
            write_plot_file("plot-prec-%s-rcut-%s-dens-%s" % (prec,rcut,dens), plot_strings)



















