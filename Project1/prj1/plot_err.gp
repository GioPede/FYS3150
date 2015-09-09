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

# Parametri Specifici Del Grafico
set output 'plot_err.png'
#set xformat "%.0f"
#set yformat "%.0f"
set xrange [0:1]
#set yrange [0:1.5]
#set logscale x 
set logscale y
set xlabel "x"
set ylabel "u(x)"

# Plot del modulo
plot 'Data/10err.txt' ls 1 w linespoints ti '10', 'Data/100err.txt' ls 2  w linespoints ti '100', 'Data/1000err.txt' ls 3  w linespoints ti '1000', 'Data/10000err.txt' ls 4  w linespoints ti '10000', 'Data/100000err.txt' ls 5  w linespoints ti '100000' , 'Data/1000000err.txt' ls 6  w linespoints ti '1000000' 

#plot 'Data/1000err.txt' ls 3  w linespoints ti '1000', 'Data/10000err.txt' ls 4  w linespoints ti '10000', 'Data/100000err.txt' ls 5  w linespoints ti '100000' , 'Data/1000000err.txt' ls 6  w linespoints ti '1000000' 


set output 'plot_errLU.png'
plot  'Data/10errLU.txt' ls 1 w linespoints ti '10', 'Data/100errLU.txt' ls 2  w linespoints ti '100', 'Data/1000errLU.txt' ls 3  w linespoints ti '1000',

