#!/usr/bin/gnuplot

set xlabel "particles"
set ylabel "time [s]"

plot './timings.d.dat' using 1:4 title "ik" w lp, \
     './timings.d.dat' using 1:7 title "ik-i" w lp, \
     './timings.d.dat' using 1:10 title "ad" w lp, \
     './timings.d.dat' using 1:13 title "ad-i" w lp, \
     './timings.d.dat' using 1:3 title "ik (k)" w lp, \
     './timings.d.dat' using 1:6 title "ik-i (k)" w lp, \
     './timings.d.dat' using 1:9 title "ad (k)" w lp, \
     './timings.d.dat' using 1:12 title "ad-i (k)" w lp

pause -1
