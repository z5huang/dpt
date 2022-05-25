unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')

output='wf-maxmin'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 4in,2in  #fontscale 0.6
  set output output.'.tex'
}

fn='wf_maxmin'
dat=fn.'.dat'
par=fn.'.par'
load par

iband=1
###
# do stuff
set xrange [0:4]
set yrange [-0.02:1.02]
set xtics 0.5
set ytics 0,0.2,1
set mxtics 5
set mytics 5

set xlabel '$t$'
set ylabel 'MaxMin'

bounds = "0.338 1.195 2.59 3.253 3.6"

setbound(b,yfrom,yto) = sprintf('set for [bb in "%s"] arrow from first bb, first %g to first bb, first %g nohead lt 0 front', b,yfrom,yto)

eval setbound(bounds,-0.02,1.02)

setlabel(txt,x,y,etc) = sprintf('set label "%s" at first %g, first %g front %s', txt, x,y,etc)
eval setlabel('\\tiny{[0,0,0] }', 0.01, 0.2, '')
eval setlabel('\\tiny{[0,1,-1]}', 0.6 , 0.2, '')
eval setlabel('\\tiny{[0,0,0 ]}', 1.7 , 0.2, '')
eval setlabel('\\tiny{[1,-1,0]}', 2.7 , 0.2, '')
eval setlabel('\\tiny{[1,-2,1]}', 3.23, 0.2, '')
eval setlabel('\\tiny{[2,-3,1]}', 3.61, 0.2, '')

eval setlabel('1', 0.1 , 0.25, '')
eval setlabel('2', 0.7 , 0.25, '')
eval setlabel('3', 1.8 , 0.25, '')
eval setlabel('4', 2.8 , 0.25, '')
eval setlabel('5', 3.35, 0.25, '')
eval setlabel('6', 3.7 , 0.25, '')

set key center bottom at first 2,first 0.7

#p for [iband=1:3] dat every ::iband u 1:(column(1+iband)) w linesp pi 3 pt 2*iband tit sprintf(iband == 1 ? '$n = %d$' : '$%d$',iband)
p for [iband=1:3] dat every ::iband u 1:(column(1+iband)) w linesp pi 3 pt 2*iband tit sprintf('$\psi^{(%d)}$', iband)

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
