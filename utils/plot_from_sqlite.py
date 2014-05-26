#!/usr/bin/python

import sqlite3
import os
import sys

argc = len(sys.argv)

method_ids = {'ik': 1, 'ik-r': 2, 'ik-i': 3, 'ad': 4, 'ad-r': 5, 'ad-i': 6}

if argc == 2:
    f = sys.stdout
elif argc == 3:
    f = open(sys.argv[2], "w")
else:
    print "Usage:\n"
    sys.exit(255)

conn = sqlite3.connect(sys.argv[1])
c = conn.cursor()

values  = c.execute("SELECT DISTINCT type FROM data").fetchall()

types = []

for i in values:
    types.append(i[0])

values  = c.execute("SELECT DISTINCT prec FROM data").fetchall()

precs = []

for i in values:
    precs.append(i[0])

values  = c.execute("SELECT DISTINCT rcut FROM data").fetchall()

rcuts = []

for i in values:
    rcuts.append(i[0])

values  = c.execute("SELECT DISTINCT dens FROM data").fetchall()

densities = []

pold = 0.0

for i in values:
    densities.append(i[0])

    t = 1
    dens = 1.0
    rcut = 3.0
    blocks = 0
    dens_tmp = []
    dens_tmp.append(1.0);
    for dens in dens_tmp:
        query = "SELECT npart,timing,prec,type,dens,rcut FROM data WHERE (dens = %e) AND (rcut = %e) AND ( (npart = 100000) OR (npart = 50000) OR (npart=1000) OR (npart=10000) OR (npart=100)) ORDER BY npart,type,prec" % (dens, rcut)
        values = c.execute(query).fetchall()
        for row in values:
            if row[3] != pold:
                f.write("\n\n")
                pold = row[3]
                blocks += 1
                f.write("%s\n" % method_ids.keys()[row[3]-1])
            f.write(" ".join(map(str, row)) + " %s\n" % method_ids.keys()[row[3]-1])

        f.write("\n\n")
                
print "Wrote %d blocks." % blocks









