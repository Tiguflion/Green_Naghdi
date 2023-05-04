 set title "Hauteur d'eau à T =   0.29999999999999999      "
 set terminal png
 set output "plot.0299.png"
 plot "data/data0299.dat" using 1:2 title "h app"w l, "data/data0299.dat" using 1:3 title "h theo" w l
