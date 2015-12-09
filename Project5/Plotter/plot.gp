# Definizioni Generali di Stile

set mxtics 10
set mytics 10
set key top right
set grid
set bars small
set terminal pngcairo enhanced font 'Arial, 20' background "#ffffff" fontscale 1.0 dashed size 1280,1280

set style line 1  linetype 1 linecolor rgb "red"  linewidth 2 
set style line 3  linetype 1 linecolor rgb "yellow"  linewidth 2 
set style line 4  linetype 1 linecolor rgb "green"  linewidth 2 
set style line 6  linetype 1 linecolor rgb "blue"  linewidth 2 
set style line 5  linetype 1 linecolor rgb "cyan"  linewidth 2 
set style line 2  linetype 1 linecolor rgb "orange"  linewidth 2 
set style line 7  linetype 1 linecolor rgb "black"  linewidth 2

# Parametri Specifici Del Grafico
set output 'plot_energy.png'
#set xformat "%.0f"
#set yformat "%.0f"
#set xrange [0:]
#set yrange [0:1.5]
#set logscale x 
#set logscale y
set xlabel "number of step"
set ylabel "Energy (eV)"

# Plot del modulo
plot  'data.txt' u ($1):($4) ls 1 w linespoints ti sprintf("E_k"), \
	  'data.txt' u ($1):($5) ls 2  w linespoints ti sprintf("U"), \
	  'data.txt' u ($1):($6) ls 3  w linespoints ti sprintf("E_{TOT}")



set output 'plot_momentum.png'
#set xformat "%.0f"
set format y "%e"
#set xrange [0:]
#set yrange [0:1.5]
#set logscale x 
#set logscale y
set xlabel "number of steps"
set ylabel "Momentum"

# Plot del modulo
plot  'data.txt' u ($1):($9) ls 1 w linespoints ti sprintf("P")


set output 'plot_temperature.png'
#set xformat "%.0f"
set format y "%.3f"
#set xrange [0:]
#set yrange [0:1.5]
#set logscale x 
#set logscale y
set xlabel "number of steps"
set ylabel "Temperature (K)"

# Plot del modulo
plot  'data.txt' u ($1):($3*119.8) ls 1 w linespoints ti sprintf("T")


set output 'plot_diffusion.png'
#set xformat "%.0f"
#set format y "%e"
#set xrange [0:]
#set yrange [0:1.5]
#set logscale x 
#set logscale y
set xlabel "number of steps"
set ylabel ""

# Plot del modulo
plot  'data.txt' u ($1):($7) every ::30 ls 1 w linespoints ti sprintf("D")

set output 'plot_pressure.png'
#set xformat "%.0f"
#set format y "%e"
#set xrange [0:]
#set yrange [0:1.5]
#set logscale x 
#set logscale y
set xlabel "x"
set ylabel "u(x)"

# Plot del modulo
plot  'data.txt' u ($1):($8) ls 1 w linespoints ti sprintf("P")



set terminal pngcairo enhanced font 'Arial, 20' background "#ffffff" fontscale 1.0 dashed size 1280,1280
set output 'plot_Tfinale.png'
#set format x "%.0e"
set format y "%f"
#set xrange [2e-17:]
set yrange [0:1]
#set size ratio -1
#set logscale x 
#set logscale y
set xlabel "T_{initial} (K)"
set ylabel "T_{final}/T_{initial}"
set key top left
unset size


# Plot del modulo
plot  'temperatura.txt' u ($1):($3) ls 1 ti sprintf("T_{final}/T_{initial}")


set terminal pngcairo enhanced font 'Arial, 20' background "#ffffff" fontscale 1.0 dashed size 1280,1280
set output 'plot_Tfinale_volumi.png'
#set format x "%.0e"
set format y "%f"
#set xrange [2e-17:]
set yrange [0:2]
#set size ratio -1
#set logscale x 
#set logscale y
set xlabel "Lattice Constant (A)"
set ylabel "T_{final}/T_{initial}"
set key top left
unset size

# Plot del modulo
plot  'volumi.txt' u ($1):($3 / 300) ls 1 ti sprintf("T_{final}/T_{initial}")




set terminal pngcairo enhanced font 'Arial, 20' background "#ffffff" fontscale 1.0 dashed size 1280,1280
set output 'plot_error.png'
set format x "%.0e"
set format y "%.0e"
set xrange [2e-17:]
#set yrange [0:1.5]
set size ratio -1
set logscale x 
set logscale y
set xlabel "{/Symbol D}t (s)"
set ylabel "{/Symbol s}_{E} (eV)"
set key top left

# Plot del modulo
plot  'verlet.txt'  ls 1 w linespoints ti sprintf("Verlet"),\
	  'euler.txt'   ls 2 w linespoints ti sprintf("Euler"),




