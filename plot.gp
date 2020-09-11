#set pm3d map
#set title 'Benney explicit solve with dt=2.73e-14, F=0, and initially h=1+0.001sin(2πx)' 
#set title 'Benney equation explicit solve with F=0 initial h=1+1e-2sin(2πx/64), dt=1e-5, 128 grid points'
#set xlabel 'x'
#set ylabel 't'
#splot 'out/explicit.out' u 2:1:3

set terminal gif animate delay 2
set output 'animation.gif'
stats 'out/explicit.out'

set xrange [0:64]
set xtics 0,8,64

set y2tics

do for [i=0:STATS_blank-1] {
plot 'out/explicit.out' every :::i::i u 2:3 w lines, 'out/explicit.out' every :::i::i u 2:4 w lines axes x1y2
}