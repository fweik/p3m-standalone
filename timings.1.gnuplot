#!/usr/bin/gnuplot

set xlabel "particles"
set ylabel "time [s]"

plot './timings.1.dat' using 1:2 title "ik" w lp, \
     './timings.1.dat' using 1:3 title "ik-i" w lp, \
     './timings.1.dat' using 1:4 title "ad" w lp, \
     './timings.1.dat' using 1:5 title "ad-i" w lp

pause -1
