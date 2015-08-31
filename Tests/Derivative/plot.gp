# Definizioni Generali di Stile

set mxtics 10
set mytics 10
set key bmargin
set grid
set bars small
set terminal pngcairo enhanced font 'Arial, 20' background "#ffffff" fontscale 1.0 dashed size 1280,1280

set style line 1  linetype 2 linecolor rgb "red"  linewidth 1 pointtype 2 pointinterval 0
set style line 2  linetype 2 linecolor rgb "green"  linewidth 1 pointtype 2 pointinterval 0

# Parametri Specifici Del Grafico
set output 'plot.png'
#set xformat "%.0f"
#set yformat "%.0f"
set xrange [10**(-14):1]
#set yrange [0.01:20]
set logscale x 
set logscale y
set xlabel "h-step"
set ylabel "Error"

# Plot del modulo
plot 'data_err.txt' u ($1):($4) ls 1 ti 'Total Error'


