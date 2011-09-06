#!/usr/bin/gnuplot

set xlabel "particles"
set ylabel "time [s]"

filename='./timings.d.weak.10k.dat'

plot filename using 1:4 title "ik" w lp, \
     filename using 1:7 title "ik-i" w lp, \
     filename using 1:10 title "ad" w lp, \
     filename using 1:13 title "ad-i" w lp, \
     filename using 1:3 title "ik (k)" w lp, \
     filename using 1:6 title "ik-i (k)" w lp, \
     filename using 1:9 title "ad (k)" w lp, \
     filename using 1:12 title "ad-i (k)" w lp

pause -1
