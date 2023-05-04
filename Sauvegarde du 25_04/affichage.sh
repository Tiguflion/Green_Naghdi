#!/bin/bash 

function tout(){
       gnuplot -persist hauteur_eau.gnu
       gnuplot -persist debit.gnu
       gnuplot -persist debit_vert.gnu
       gnuplot -persist sigma.gnu
       #gnuplot -persist maree.gnu
       echo "Affichage terminé"
}      

function hauteur(){
	gnuplot -persist hauteur_eau.gnu
}

function maree(){
	gnuplot -persist maree.gnu
}

function debit(){
	gnuplot -persist debit.gnu
	gnuplot -persist debit_vert.gnu
}

tout

