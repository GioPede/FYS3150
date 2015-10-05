# Definizioni Generali di Stile

set mxtics 10
set mytics 10
set key top right
set grid
set bars small
set terminal pngcairo enhanced font 'Arial, 20' background "#ffffff" fontscale 1.0 dashed size 1280,1280

set style line 1  linetype 1 linecolor rgb "red"  linewidth 1 
set style line 3  linetype 1 linecolor rgb "yellow"  linewidth 1 
set style line 4  linetype 1 linecolor rgb "green"  linewidth 1 
set style line 6  linetype 1 linecolor rgb "blue"  linewidth 1 
set style line 5  linetype 1 linecolor rgb "cyan"  linewidth 1 
set style line 2  linetype 1 linecolor rgb "orange"  linewidth 1 
set style line 7  linetype 1 linecolor rgb "black"  linewidth 1

# Parametri Specifici Del Grafico
set output 'plot.png'
#set xformat "%.0f"
#set yformat "%.0f"
#set xrange [0:]
#set yrange [0:1.5]
#set logscale x 
#set logscale y
set xlabel "x"
set ylabel "u(x)"

f(x) = 1 - (1 - exp(-10)) * x - exp(-10*x)

# Plot del modulo
plot  'Data/0.01.txt' every ::1 ls 1  w linespoints ti sprintf("{/Symbol w} _{%s} = %s",  columnhead(1), columnhead(2)), \
	  'Data/0.5.txt'  every ::1 ls 2  w linespoints ti sprintf("{/Symbol w} _{%s} = %s",  columnhead(1), columnhead(2)), \
	  'Data/1.txt'    every ::1 ls 4  w linespoints ti sprintf("{/Symbol w} _{%s} = %s",  columnhead(1), columnhead(2)), \
      'Data/5.txt'    every ::1 ls 6  w linespoints ti sprintf("{/Symbol w} _{%s} = %s",  columnhead(1), columnhead(2)), 
