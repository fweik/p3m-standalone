#!/usr/bin/python

import re
import os
import sys
import sqlite3

conn = sqlite3.connect('data.db')
c = conn.cursor()


file_re = re.compile('p3m-(.*?)-([0-9]*)-to-([0-9]*)-prec-(.*)-rcut-(.*)-dens-(.*)-q-(.*)\.dat');

method_ids = {'ik': 1, 'ik-r': 2, 'ik-i': 3, 'ad': 4, 'ad-r': 5, 'ad-i': 6}

params = ['method', 'lownum', 'highnum', 'prec', 'rcut', 'dens', 'q']

datasets = []

#CREATE TABLE data (type INTEGER, prec REAL, rcut REAL, dens REAL, q REAL, npart INTERGER, timing REAL, t_c REAL, t_g REAL, t_f REAL, cao INTEGER, mesh INTEGER, alpha REAL);

for i in os.listdir(sys.argv[1]):
    match = file_re.match(i)
    if not match is None:
        seq = zip(params, match.groups())
        seq.append(('filename', i))
        datasets.append(dict(seq))
        ds = datasets[-1]
        if ds['method'] in method_ids.keys():
            method_id = method_ids[ds['method']]
        else:
            method_id = 0

        with open(ds['filename']) as f:
            for line in f:
                if line.lstrip()[0] == '#':
                    continue
                values = map(float,  line.split())
                query = "INSERT INTO data VALUES (%d, %e, %e, %e, %e, %d, %e, %e, %e, %e, %d, %d, %e)" % (method_id, float(ds['prec']), float(ds['rcut']), float(ds['dens']), float(ds['q']), int(values[0]), values[5],values[6],values[7],values[8], int(values[2]), int(values[3]), values[4])
                c.execute(query)
                print query
        conn.commit()
        f.close()







    



















