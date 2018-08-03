set term postscript landscape monochrom  
#set term x11 
set out "graph_1.ps"
set title "Values at Nodes (/home/jmoran/Desktop/Complete_Thin_Cylinder/CTC_Dynamic_2rps.frd)"
set grid
set ylabel " ALL      "
set xlabel " Time "
plot "graph_1.out" using 3:5 title 'Node=19861' with linespoints
 