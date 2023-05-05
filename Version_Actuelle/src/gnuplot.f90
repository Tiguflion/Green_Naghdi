MODULE mod_gnuplot

USE numerics

IMPLICIT NONE 

CONTAINS 

SUBROUTINE init_fichier_data_frame()
	CALL system("mkdir data" )
	CALL system("mkdir frame")
END SUBROUTINE init_fichier_data_frame

SUBROUTINE suppr_fichier_data_frame()
	CALL system("rm -rf data")
	CALL system("rm -rf frame")
END SUBROUTINE suppr_fichier_data_frame


SUBROUTINE script_sauv_png(title)
	CHARACTER(21), INTENT(IN) :: title
	CHARACTER(50) :: var_poubelle

	OPEN(UNIT = 111, FILE = "hauteur_eau.gnu", STATUS = 'UNKNOWN')
		WRITE(111,*) 'set title "',title,date, '"'
		WRITE(111,*) 'set terminal png'
		WRITE(111,*) 'set output "frame/',frame_name,'"'
		WRITE(111,*) 'plot "data/',data_name,'" using 1:2 title "h app"w l, "data/'&
		&,data_name, '" using 1:3 title "h theo" w l' 
	CLOSE(111)
	
	CALL system("gnuplot hauteur_eau.gnu")
	!CALL system("rm -rf hauteur_eau.gnu")
END SUBROUTINE script_sauv_png

SUBROUTINE Ecrits_data()
CHARACTER(19) :: nom_file 
REAL(RP) :: pos
INTEGER :: i

WRITE(nom_file,'(a,i4.4,a)') "data/data",niter,".dat"

open(unit=21,file=  nom_file,status = 'UNKNOWN')

	do i=1,Nx
        	write(21,*) X(i),W(1,i),h_temp(i)
	enddo
close(21)

END SUBROUTINE Ecrits_data


SUBROUTINE init_frame_data_names()
	WRITE(frame_name,'(a,i4.4,a)') 'plot.',niter,'.png'
	WRITE(data_name,'(a,i4.4,a)') 'data',niter,'.dat'
END SUBROUTINE init_frame_data_names


SUBROUTINE Animation()
	!CALL system("ffmpeg -i frame/plot%.png animation.avi")
	CALL system("ffmpeg -r 20 -i frame/plot.%04d.png -pix_fmt yuv420p out.mp4")
	CALL system("rm -f *.png")
	
END SUBROUTINE Animation

	
	
	



END MODULE mod_gnuplot
