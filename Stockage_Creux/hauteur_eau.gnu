 set title "Hauteur d'eau à T =   0.29999999999999999      "
 set terminal png
 set output "frame/plot.1219.png"
 plot "data/data1219.dat" using 1:2 title "h app"w l, "data/data1219.dat" using 1:3 title "h theo" w l
