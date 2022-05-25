unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')

output='ssh3-overlap'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 3in,2in  #fontscale 0.6
  set output output.'.tex'
}

fn='ssh3-overlap'
dat=fn.'.dat'
par=fn.'.par'
load par

###
# do stuff

set xlabel '$k/\pi$'
set ylabel '$\left| \langle \phi^{(n)} | \psi\rangle \right|^2$'
set xrange [0:2]
set yrange [-0.05:1.05]
set ytics 0,0.5,1
set xtics 0,0.5,2

set grid xtics ytics
set key left center

p for [b=1:3] dat u ($1/pi):(column(1+b)) w linesp pi 10 tit sprintf( b == 1 ? '$n = %d$' : '$%d$', b)

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
