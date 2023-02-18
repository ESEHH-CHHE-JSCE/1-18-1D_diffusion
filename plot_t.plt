set xtics 60
set ytics  1
set yrange [0:6]
set xlabel "t(min)"
set ylabel "c(x0.01ppm)"
set title "concentration profile at 10km downstream section"
plot "result.txt" with lines linewidth 3
