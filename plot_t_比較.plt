set xtics 60
set ytics  1
set yrange [0:6]
set xlabel "t(min)"
set ylabel "c(x0.01ppm)"
set title "concentration profile at 10km downstream section"
  plot "result_nx500.txt" with lines linewidth 1
replot "result_nx1000.txt" with lines linewidth 2
replot "result_nx2000.txt" with lines linewidth 4
replot "result_nx4000.txt" with lines linewidth 6
