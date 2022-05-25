unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')

output='qcoef-ppt'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 4in,2in  #fontscale 0.6
  set output output.'.tex'
}

fn='qcoef'
dat=fn.'.dat'
par=fn.'.par'
load par

###
# do stuff

#
###

set xrange [0:2]
set xlabel '$k_x/\pi$'
set xtics ( 0, '$\frac{2}{3}$' 2.0/3, '$\frac{4}{3}$' 4.0/3, 2 )
set ytics 0,0.5 format '%g'
set grid

p for [b=1:3] dat u ($1/pi):(column(2*b+2)**2) w l notit

if (tex == 1){
  unset output
  set term wxt
  build = buildtex(output)
  print '"eval build" to build and preview'
} else {
  #print "press enter to continue"
  #pause -1
}
