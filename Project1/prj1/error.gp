# Definizioni Generali di Stile

set mxtics 10
set mytics 10
set key bmargin
set grid
set bars small
set terminal pngcairo enhanced font 'Arial, 20' background "#ffffff" fontscale 1.0 dashed size 1280,1280

set style line 1  linetype 1 linecolor rgb "red"  linewidth 1 
set style line 3  linetype 1 linecolor rgb "yellow"  linewidth 1 
set style line 4  linetype 1 linecolor rgb "green"  linewidth 1 
set style line 6  linetype 1 linecolor rgb "blue"  linewidth 1 
set style line 5  linetype 1 linecolor rgb "cyan"  linewidth 1 
set style line 2  linetype 1 linecolor rgb "orange"  linewidth 1 
set style line 7  linetype 1 linecolor rgb "black"  linewidth 4

# Parametri Specifici Del Grafico
set output 'error.png'
#set xformat "%.0f"
#set yformat "%.0f"
set logscale x 
set logscale y
set xlabel "h-step"
set ylabel "Error"

f(x) = 1 - (1 - exp(-10)) * x - exp(-10*x)

# Plot del modulo
plot  'Data/error.txt' ls 1  ti 'Relative Error'

