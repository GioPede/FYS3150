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

#set xformat "%.0f"
#set yformat "%.0f"
#set yrange [0:1.5]
#set logscale x 




# Plot flips
set output 'Flip_Counts.png'
set xlabel "K_BT"
set ylabel "Normalized Count of Accepted Moves"

plot  'flips.txt'  u ($1):($2/10000000) ls 1  w linespoints ti sprintf(" Accepted Moves " )

# Plot energy
set output 'Energy_Counts.png'
set xlabel "E"
set ylabel "P(E)"
set xrange [-800:0]

plot  '10.txt'  u ($1):($2/10000000) ls 1  w linespoints ti sprintf(" 1.0 " ),\
      '14.txt'  u ($1):($2/10000000) ls 2  w linespoints ti sprintf(" 1.4 " ),\
      '18.txt'  u ($1):($2/10000000) ls 3  w linespoints ti sprintf(" 1.8 " ),\
      '22.txt'  u ($1):($2/10000000) ls 4  w linespoints ti sprintf(" 2.2 " ),\
      '24.txt'  u ($1):($2/10000000) ls 5  w linespoints ti sprintf(" 2.4 " ),\
      '26.txt'  u ($1):($2/10000000) ls 6  w linespoints ti sprintf(" 2.6 " ),\
      '30.txt'  u ($1):($2/10000000) ls 7  w linespoints ti sprintf(" 3.0 " )



