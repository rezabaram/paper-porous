for i in *full.dat
do 

file=${i%.*}
 
echo $file
sed "s/#1/${file}/" lc_x.2.F90 > lc_x_${file}.F90
sed "s/#1/${file}/" condorSerial > condorSerial_${file}.cmd 

#condor_compile gfortran lc_x_${file}.F90 -o ${file}.run

done
