set term postscript landscape monochrom  
#set term x11 
set out "graph_0.ps"
set title "Values at Nodes (/home/jmoran/Desktop/Complete_Thin_Cylinder/CTC_Dynamic_2rps.frd)"
set grid
set ylabel " ALL      "
set xlabel " Dataset "
plot "graph_0.out" using 2:5 title 'Node=19861' with linespoints
