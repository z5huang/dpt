unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')
eval setdef('band','1')

output='node-vs-flux'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 4in,2in  #fontscale 0.6
  set output output.'.tex'
}

wf='wf.dat'
curv='curv.dat'
conn='conn.dat'

NBAND=3 # number of bands
NORM=100 # normalization of the connection vector
#band=1	# which band to plot

# column of a wavefunction component (basis) of band
wf_col(band,basis,nband) = 2 + 2*nband*(band - 1) + (2*basis - 1)
# column of the curvature of band
curv_col(band) = band+2
# column of the x- and y-components of the Berry connection of band
conn_x_col(band) = 1+band*2
conn_y_col(band) = 2+band*2

###
# do stuff
#set xrange [0:2]
#set yrange [0:2]
set xyplane 0.1
set xlabel '$k_x$'
set ylabel '$k_y$'

set multiplot layout 2,NBAND

set size ratio 1

#set cbrange [0:1]
do for [basis=1:NBAND] {
    set title sprintf('$|\langle %d|\psi^%d\rangle|', basis, band)
    p wf u ($1/pi):($2/pi):(column(wf_col(band,basis,NBAND))) w image notit
}
set title sprintf('Berry curvature density of $\psi^%d$',band)
set cbrange [*:*]
p curv u ($1/pi):($2/pi):(column(curv_col(band))) w image notit, conn u ($1/pi):($2/pi):(column(conn_x_col(band))/NORM):(column(conn_y_col(band))/NORM) w vector notit

unset multiplot
#
###

if (tex == 1){
  unset output
  set term wxt
  build = buildtex(output)
  print '"eval build" to build and preview'
} else {
  #print "press enter to continue"
  #pause -1
}
