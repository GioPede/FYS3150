# Definizioni Generali di Stile

set mxtics 5
set mytics 5
set key bmargin
set grid
set bars small
set terminal pngcairo enhanced font 'Arial, 20' background "#ffffff" fontscale 1.0 dashed size 1280,720

set style line 1  linetype 2 linecolor rgb "red"  linewidth 1 pointtype 2 pointinterval 0
set style line 2  linetype 2 linecolor rgb "green"  linewidth 1 pointtype 2 pointinterval 0

# Parametri Specifici Del Grafico
set output 'plot.png'
#set xformat "%.0f"
#set yformat "%.0f"
#set xrange [-1:]
#set yrange [0.01:20]
set logscale x 
set logscale y
set xlabel "Max N"
set ylabel "Discrepancy"

# Plot del modulo
plot 'data.txt' u ($1):($2) ls 1 ti 'Float', 'data.txt' u ($1):($3) ls 2 ti 'Double'


