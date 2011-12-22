L=256
dx=1
dt=1
tau=0.8
rho0=1000

for mth in mrt
do
(
    for dir in eta*_xi*
    do 	
	file=$(ls -1 ${dir}/profile-z_${dir}_${mth}_t00070000-*.asc)
	awk '(NR>8)' $file > tmp
	phi=$(awk '{sum+=$9} END {print 1-sum/256}' ${file})
	
	head -n1 tmp > tmpi
	tail -n1 tmp > tmpf
	
	rhoi=$(awk < tmpi '{ print $8 }')
	rhof=$(awk < tmpf '{ print $8 }')
	
	v=$(awk '{sum+=$7*$9} END { print sum/NR}' tmp)

cat > tmp.f <<- -eof1-
      program tmp
      double precision nu, A, v, cs2
      double precision dp, k, rhoi, rhof

      cs2 = $dx**2 / $dt**2 /3.
      nu =  cs2 * $dt / 2. * (2. * $tau - 1.)
      rhof = $rhof 
      rhoi = $rhoi

      A = $L * $L
      v = $v * $dx / $dt

      dp = $rho0 * ( rhof - rhoi ) * cs2 / 
     $  ( $L * $dx )

      k = - v * $rho0 * nu / dp

      write(*,*) dp, v, k/A
      end program
-eof1-
         gfortran tmp.f
	 output=$(./a.out)
	 
	 echo ${dir:3:3} ${dir:9:3} $phi $output 
	 if [ ${dir:9:3} == 2.5 ]
	 then
	     echo " "
	 fi
    done
) > data_${mth}.dat 
done

rm tmp tmpf tmpi tmp.f a.out
