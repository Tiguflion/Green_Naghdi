 set title "Hauteur d'eau à T =   0.22743544728602280      "
 set terminal png
 set output "frame/plot.0001.png"
 plot "data/data0001.dat" using 1:2 title "h app"w l, "data/data0001.dat" using 1:3 title "h theo" w l
