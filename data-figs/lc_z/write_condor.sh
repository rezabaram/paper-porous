for i in *full.dat
do 

file=${i%.*}
 
echo $file
sed "s/#1/${file}/" lc.F90 > lc_${file}.F90
sed "s/#1/${file}/" condorSerial > condorSerial_${file}.cmd 

#condor_compile gfortran lc_${file}.F90 -o ${file}.run

done
