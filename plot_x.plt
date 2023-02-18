set xtics 1
set ytics 1
set yrange [0:10]
set xlabel "x(km)"
set ylabel "c(x0.01ppm)"
set title "concentration profile in x-direction"
  plot "result.0000" title   "0(min)" with lines linewidth 6
replot "result.0060" title  "60(min)" with lines linewidth 5
replot "result.0120" title "120(min)" with lines linewidth 4
replot "result.0180" title "180(min)" with lines linewidth 3
replot "result.0240" title "240(min)" with lines linewidth 2
