#! /bin/bash

cat > temp.sh <<- -eof1-
#! /bin/bash

file=\$1

gnuplot temp.gp

latex \${file}.tex
dvips \${file}.dvi
#sed -i 's/^%%BoundingBox: 0 0 241 236.*$/%%BoundingBox: 0 10 241 189/' \${file}.ps
ps2eps \${file}.ps
epstopdf \${file}.eps
#ps2pdf \${file}.ps

for i in .tex .log .aux .dvi -inc.eps
# .ps
do
rm \${file}\${i}
done
-eof1-
chmod +x temp.sh

head2d="set term epslatex color standalone rounded size 8.5cm, 6.3cm lw 1.0 9.0"
head3d="set term epslatex color standalone rounded size 8.5cm, 8.3cm lw 1.0 9.0 header '\\usepackage{rotating}'"
margin2d="set tmargin 3; set bmargin 3.5; set lmargin 10.5; set rmargin 2"
margin3d="set lmargin 4; set rmargin 0"

plot=porosity
cat > temp.gp <<- -eof1-
$head3d
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%2.2f}"
set ztics 0.02
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\\\\begin{turn}{90}\$\\\\phi\$\\\\end{turn}" offset -1.5,0
unset title

splot \
 'z-permeability/data.dat' u 1:2:(1-\$3) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=permeability_z
cat > temp.gp <<- -eof1-
$head3d
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%2.1f}"
set ztics 0.5
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\\\\begin{turn}{90} \$\\\\kappa_{z}\\\\,[10^{6}]\$ \\\\end{turn}" offset -1.5,0

splot \
 'z-permeability/data.dat' u 1:2:(\$7*1e6) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=permeability_x
cat > temp.gp <<- -eof1-
$head3d
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%2.1f}"
set ztics 0.5
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\\\\begin{turn}{90} \$\\\\kappa_{x}\\\\,[10^{6}]\$ \\\\end{turn}" offset -1.5,0

splot \
 'x-permeability/data.dat' u 1:2:(\$7*1e6) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=permeability_ratio
cat > temp.gp <<- -eof1-
$head3d
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%2.1f}"
set ztics 0.2
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\\\\begin{turn}{90} \$\\\\kappa_{z}/\\\\kappa_{x}\$ \\\\end{turn}" offset -1.5,0

splot \
 '< paste z-permeability/data.dat x-permeability/data.dat' u 1:2:(\$7/\$14) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=lc_z
cat > temp.gp <<- -eof1-
$head3d
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%2.3f}"
set ztics 0.001
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\\\\begin{turn}{90} \$l_{c}\$ \\\\end{turn}" offset -1.5,0

splot \
 'lc_z/lc_512.asc' u 1:2:3 notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=lc_x
cat > temp.gp <<- -eof1-
$head3d
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%2.3f}"
set ztics 0.001
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\\\\begin{turn}{90} \$l_{c}\$ \\\\end{turn}" offset -1.5,0

splot \
 'lc_x/lc_512.2.asc' u 1:2:3 notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=k_lc_z
cat > temp.gp <<- -eof1-
$head3d
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%2.2f}"
set ztics 0.01
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\\\\begin{turn}{90} \$\\\\kappa/{l_{c}}^{2}\$ \\\\end{turn}" offset -1.5,0

splot \
 '< paste z-permeability/data.dat lc_z/lc_512.asc' u 1:2:(\$7/(\$10)**2) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=k_lc_x
cat > temp.gp <<- -eof1-
$head3d
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%2.2f}"
set ztics 0.01
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\$\\\\kappa/{l_{c}}^{2}\$" offset -1.5,0

splot \
 '< paste x-permeability/data.dat lc_x/lc_512.2.asc' u 1:2:(\$7/(\$10)**2) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=phi_k_lc_z
cat > temp.gp <<- -eof1-
$head2d
set out "${plot}.tex"

$margin2d

set xrange [0.24:0.36]
set format x "%2.2f"
set xtics 0.04

set yrange [0.010:0.040]
set format y "%2.3f"
set ytics 0.005

set xlabel "\$\\\\phi\$"
set ylabel "\$ \\\\kappa/{l_{c}}^{2}\$"

plot \
 '< paste z-permeability/data.dat lc_z/lc_512.asc' u (1-\$3):(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}


plot=zeta_k_lc_z
cat > temp.gp <<- -eof1-
$head2d
set out "${plot}.tex"

$margin2d

set xrange [1.0:2.5]
set format x "%2.1f"
set xtics 0.5

set yrange [0.010:0.040]
set format y "%2.3f"
set ytics 0.005

set xlabel "\$\\\\zeta\$"
set ylabel "\$ \\\\kappa/{l_{c}}^{2}\$"

plot \
 '< paste z-permeability/data.dat lc_z/lc_512.asc' u 1:(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}


plot=xi_k_lc_z
cat > temp.gp <<- -eof1-
$head2d
set out "${plot}.tex"

$margin2d

set xrange [1.0:2.5]
set format x "%2.1f"
set xtics 0.5

set yrange [0.010:0.040]
set format y "%2.3f"
set ytics 0.005

set xlabel "\$\\\\eta\$"
set ylabel "\$ \\\\kappa/{l_{c}}^{2}\$"

plot \
 '< paste z-permeability/data.dat lc_z/lc_512.asc' u 2:(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}


plot=phi_k_lc_x
cat > temp.gp <<- -eof1-
$head2d
set out "${plot}.tex"

$margin2d

# set xrange [0.24:0.36]
# set format x "%2.2f"
# set xtics 0.04

# set yrange [0.010:0.040]
# set format y "%2.3f"
# set ytics 0.005

set xlabel "\$\\\\phi\$"
set ylabel "\$ \\\\kappa/{l_{c}}^{2}\$"

plot \
 '< paste x-permeability/data.dat lc_x/lc_512.2.asc' u (1-\$3):(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}


plot=zeta_k_lc_x
cat > temp.gp <<- -eof1-
$head2d
set out "${plot}.tex"

$margin2d

set xrange [1.0:2.5]
set format x "%2.1f"
set xtics 0.5

set yrange [0.005:0.020]
set format y "%2.3f"
set ytics 0.005

set xlabel "\$\\\\zeta\$"
set ylabel "\$ \\\\kappa/{l_{c}}^{2}\$"

plot \
 '< paste x-permeability/data.dat lc_x/lc_512.2.asc' u 1:(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}


plot=xi_k_lc_x
cat > temp.gp <<- -eof1-
$head2d
set out "${plot}.tex"

$margin2d

set xrange [1.0:2.5]
set format x "%2.1f"
set xtics 0.5

set yrange [0.005:0.020]
set format y "%2.3f"
set ytics 0.005

set xlabel "\$\\\\eta\$"
set ylabel "\$ \\\\kappa/{l_{c}}^{2}\$"

plot \
 '< paste x-permeability/data.dat lc_x/lc_512.2.asc' u 2:(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}



plot=kozeny_carman
cat > temp.gp <<- -eof1-
$head3d
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%2.1f}"
set ztics 0.5
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\\\\begin{turn}{90} \$\\\\kappa\\\\,[10^{6}]\$ \\\\end{turn}" offset -1.5,0

D(x,y,z)=(x*y*z)**(1./3.)
sphericity(x,y,z)= D(x,y,z)**2 / ( ( (x*y)**p + (x*z)**p + (y*z)**p )/3.0 )**(1.0/p)
p=1.6075

splot \
 '< paste ellips.dat z-permeability/data.dat' u 1:2:(sphericity(\$3,\$4,\$5)**2*(2.0*D(\$3,\$4,\$5))**2/150.0*(1-\$8)**3/(\$8)**2*1e6) notitle w pm3d
-eof1-
./temp.sh ${plot}

rm temp.sh temp.gp

