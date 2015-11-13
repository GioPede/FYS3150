# Definizioni Generali di Stile

set mxtics 10
set mytics 10
set key top right
set grid
set bars small
set terminal pngcairo enhanced font 'Arial, 20' background "#ffffff" fontscale 1.0 dashed size 1280,1280

set style line 1  linetype 1 linecolor rgb "red"     linewidth 2 
set style line 3  linetype 1 linecolor rgb "yellow"  linewidth 2 
set style line 4  linetype 1 linecolor rgb "green"   linewidth 2 
set style line 6  linetype 1 linecolor rgb "blue"    linewidth 2 
set style line 5  linetype 1 linecolor rgb "cyan"    linewidth 2 
set style line 2  linetype 1 linecolor rgb "orange"  linewidth 2 
set style line 7  linetype 1 linecolor rgb "black"   linewidth 2

# Parametri Specifici Del Grafico

#set xformat "%.0f"
#set yformat "%.0f"
#set xrange [0:]
#set yrange [0:1.5]
#set logscale x 
#set logscale y



# Plot Heat Capacity
set output 'Plots/Heat_Capacity.png'
set xlabel "K_BT"
set ylabel "C_V"

f(x) = abs(2.269 - x)**(-0.1)

plot  'Data/20.txt'  u ($1):($3) ls 1  w linespoints ti sprintf(" 20x20" ),\
      'Data/40.txt'  u ($1):($3) ls 2  w linespoints ti sprintf(" 40x40" ),\
      'Data/60.txt'  u ($1):($3) ls 3  w linespoints ti sprintf(" 60x60" ),\
      'Data/80.txt'  u ($1):($3) ls 4  w linespoints ti sprintf(" 80x80" ),\
      'Data/100.txt' u ($1):($3) ls 5  w linespoints ti sprintf("100x100"),\
      'Data/150.txt' u ($1):($3) ls 6  w linespoints ti sprintf("150x150"),\
      'Data/200.txt' u ($1):($3) ls 7  w linespoints ti sprintf("200x200")



# Plot Magnetization
set output 'Plots/Magnetization.png'
set xlabel "K_BT"
set ylabel "M"

plot  'Data/20.txt'  u ($1):($4) ls 1  w linespoints ti sprintf(" 20x20" ),\
      'Data/40.txt'  u ($1):($4) ls 2  w linespoints ti sprintf(" 40x40" ),\
      'Data/60.txt'  u ($1):($4) ls 3  w linespoints ti sprintf(" 60x60" ),\
      'Data/80.txt'  u ($1):($4) ls 4  w linespoints ti sprintf(" 80x80" ),\
      'Data/100.txt' u ($1):($4) ls 5  w linespoints ti sprintf("100x100"),\
      'Data/150.txt' u ($1):($4) ls 6  w linespoints ti sprintf("150x150"),\
      'Data/200.txt' u ($1):($4) ls 7  w linespoints ti sprintf("200x200")


# Plot Suceptibility
set output 'Plots/Suscettibility.png'
set xlabel "K_BT"
set ylabel "{/Symbol c}"

f(x) = abs(2.269 - x)**(-7/4)

plot  'Data/20.txt'  u ($1):($5) ls 1  w linespoints ti sprintf(" 20x20" ),\
      'Data/40.txt'  u ($1):($5) ls 2  w linespoints ti sprintf(" 40x40" ),\
      'Data/60.txt'  u ($1):($5) ls 3  w linespoints ti sprintf(" 60x60" ),\
      'Data/80.txt'  u ($1):($5) ls 4  w linespoints ti sprintf(" 80x80" ),\
      'Data/100.txt' u ($1):($5) ls 5  w linespoints ti sprintf("100x100"),\
      'Data/150.txt' u ($1):($5) ls 6  w linespoints ti sprintf("150x150"),\
      'Data/200.txt' u ($1):($5) ls 7  w linespoints ti sprintf("200x200")


# Plot energy
set output 'Plots/Energy.png'
set xlabel "K_BT"
set ylabel "E"
set key top left

plot  'Data/20.txt'  u ($1):($2) ls 1  w linespoints ti sprintf(" 20x20" ),\
      'Data/40.txt'  u ($1):($2) ls 2  w linespoints ti sprintf(" 40x40" ),\
      'Data/60.txt'  u ($1):($2) ls 3  w linespoints ti sprintf(" 60x60" ),\
      'Data/80.txt'  u ($1):($2) ls 4  w linespoints ti sprintf(" 80x80" ),\
      'Data/100.txt' u ($1):($2) ls 5  w linespoints ti sprintf("100x100"),\
      'Data/150.txt' u ($1):($2) ls 6  w linespoints ti sprintf("150x150"),\
      'Data/200.txt' u ($1):($2) ls 7  w linespoints ti sprintf("200x200")

