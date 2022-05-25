unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')

output='qcoef'
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
set xyplane 0.1
set hidden3d
set xlabel '$k_x$'
set ylabel '$k_y$'
set zlabel 'Overlap squared'
set xrange [0:2.0*pi/kxfrac]
set yrange [0:2.0*pi/kyfrac]
#sp for [b=1:q] dat u 1:2:(column(2+b*2)**2) w l tit sprintf('$|\left<\psi_%d | \psi \right>|^2$', b)

b=1
sp dat u 1:2:(column(2+b*2)**2) w l tit sprintf('$|\left<\psi_%d | \psi \right>|^2$', b)
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
