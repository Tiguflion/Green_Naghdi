 set title "Hauteur d'eau à T =    2.0000000000000000      "
 set terminal png
 set output "frame/plot.1869.png"
 plot "data/data1869.dat" using 1:2 title "h app"w l, "data/data1869.dat" using 1:3 title "h theo" w l
