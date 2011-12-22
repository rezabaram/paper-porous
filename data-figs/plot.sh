#! /bin/bash

cat > temp.sh <<- -eof1-
#! /bin/bash

file=\$1

gnuplot temp.gp

latex \${file}.tex
dvips \${file}.dvi
ps2pdf \${file}.ps

for i in .tex .log .aux .dvi -inc.eps .ps
do
rm \${file}\${i}
done
-eof1-
chmod +x temp.sh

head="set term epslatex color standalone rounded size 8.5cm, 6.3cm lw 1.0 9.0"
margin2d="set tmargin 3; set bmargin 3.5; set lmargin 10.5; set rmargin 2"
margin3d="set tmargin 3; set bmargin 3.5; set lmargin 6.5; set rmargin 2"

plot=porosity
cat > temp.gp <<- -eof1-
$head
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%3.2f}"
set ztics 0.02
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\$\\\\phi\$" offset -1.5,0
#set title 'Porosity'

splot \
 'original_orientation_full/data.dat' u 1:2:(1-\$3) notitle w pm3d
-eof1-
./temp.sh ${plot}

plot=old_porosity
cat > temp.gp <<- -eof1-
$head
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%3.2f}"
set ztics 0.02
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\$\\\\phi\$" offset -1.5,0
#set title 'Porosity'

splot \
 'original_orientation_full/density.dat' u 1:2:(1-\$3) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=permeability_z
cat > temp.gp <<- -eof1-
$head
set out "${plot}.tex"

$margin3d

set format x "\\\\tiny{%2.1f}"
set xtics 0.2 offset 0,-0.25

set format y "\\\\tiny{%2.1f}"
set ytics 0.2 offset 1.5,0.0

set format z "\\\\tiny{%3.2f}"
set ztics 0.5
unset colorbox
set ticslevel 0

set xlabel "\$\\\\zeta\$" offset 0,-0.5
set ylabel "\$\\\\eta\$" offset 0.5,0
set zlabel "\$\\\\kappa_{z}\\\\,[10^{6}]\$" offset -1.5,0
#set title 'Permeability in \$z\$-depth'

splot \
 'original_orientation_full/data.dat' u 1:2:(\$7*1e6) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=permeability_x
cat > temp.gp <<- -eof1-
$head
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
set zlabel "\$\\\\kappa_{x}\\\\,[10^{6}]\$" offset -1.5,0
#set title 'Permeability \$x\$-depth'

splot \
 'rotated_orientation_full/data.dat' u 1:2:(\$7*1e6) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=permeability_ratio
cat > temp.gp <<- -eof1-
$head
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
set zlabel "\$\\\\kappa_{z}/\\\\kappa_{x}\$" offset -1.5,0
#set title 'Permeability Ratio'

splot \
 '< paste original_orientation_full/data.dat rotated_orientation_full/data.dat' u 1:2:(\$7/\$14) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=lc_z
cat > temp.gp <<- -eof1-
$head
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
set zlabel "\$l_{c}\$" offset -1.5,0
#set title 'Carachteristic Length in \$z\$-depth'

splot \
 'lc_z/lc_512.asc' u 1:2:3 notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=lc_x
cat > temp.gp <<- -eof1-
$head
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
set zlabel "\$l_{c}\$" offset -1.5,0
#set title 'Characteristic Length in \$x\$-depth'

splot \
 'lc_x/lc_512.2.asc' u 1:2:3 notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=k_lc_z
cat > temp.gp <<- -eof1-
$head
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
#set title 'Permeability normalized by \$l_{c}\$ in \$z\$-depth'

splot \
 '< paste original_orientation_full/data.dat lc_z/lc_512.asc' u 1:2:(\$7/(\$10)**2) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=k_lc_x
cat > temp.gp <<- -eof1-
$head
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
#set title 'Permeability normalized  by \$l_{c}\$ in \$x\$-depth'

splot \
 '< paste rotated_orientation_full/data.dat lc_x/lc_512.2.asc' u 1:2:(\$7/(\$10)**2) notitle w pm3d
-eof1-
./temp.sh ${plot}


plot=phi_k_lc_z
cat > temp.gp <<- -eof1-
$head
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
#set title "\$\\\\kappa\$ normalized by \$l_{c}\$ in \$z\$-depth"

plot \
 '< paste original_orientation_full/data.dat lc_z/lc_512.asc' u (1-\$3):(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}


plot=eta_k_lc_z
cat > temp.gp <<- -eof1-
$head
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
#set title "\$\\\\kappa\$ normalized by \$l_{c}\$ in \$z\$-depth"

plot \
 '< paste original_orientation_full/data.dat lc_z/lc_512.asc' u 1:(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}


plot=xi_k_lc_z
cat > temp.gp <<- -eof1-
$head
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
#set title "\$\\\\kappa\$ normalized by \$l_{c}\$ in \$z\$-depth"

plot \
 '< paste original_orientation_full/data.dat lc_z/lc_512.asc' u 2:(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}


plot=phi_k_lc_x
cat > temp.gp <<- -eof1-
$head
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
#set title "\$\\\\kappa\$ normalized by \$l_{c}\$ in \$\x\$-depth"

plot \
 '< paste rotated_orientation_full/data.dat lc_x/lc_512.asc' u (1-\$3):(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}


plot=eta_k_lc_x
cat > temp.gp <<- -eof1-
$head
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
#set title "\$\\\\kappa\$ normalized by \$l_{c}\$ in \$x\$-depth"

plot \
 '< paste rotated_orientation_full/data.dat lc_x/lc_512.asc' u 1:(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}


plot=xi_k_lc_x
cat > temp.gp <<- -eof1-
$head
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
#set title "\$\\\\kappa\$ normalized by \$l_{c}\$ in \$x\$-depth"

plot \
 '< paste rotated_orientation_full/data.dat lc_x/lc_512.asc' u 2:(\$7/(\$10)**2) notitle pt 6 ps 2 lw 2
-eof1-
./temp.sh ${plot}


rm temp.sh temp.gp

