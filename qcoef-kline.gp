unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')

output='qcoef-kline'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle preamble '\usepackage{amsmath}' size 2in,2in  #fontscale 0.6
  set output output.'.tex'
}

fn='qcoef_kline'
dat=fn.'.dat'
par=fn.'.par'
load par

###
# do stuff
set macros

set xlabel '$\boldsymbol{k}/\pi$'
set xrange [0:nk]
set xtics @tics

set ylabel '$\left| \langle \phi^{(a)} | \psi^{('.sprintf('%d', bidx+1).')}\rangle \right|^2$'
set yrange [-0.05:1.05]
set ytics 0,0.5,1
set mytics 5

set grid lt 0 lc rgb color='black'
#set key center center spacing 2
#set key out right
set key top left spacing 2
#set key center right spacing 2
#p for [b=1:q] dat u 0:(column(2+b*2)**2) w linesp pi 5 tit sprintf('$\left|\left<\phi_%d | \psi_%d\right>\right|^2$',b,bidx+1)
#p for [b=1:q] dat u 0:(column(2+b*2)**2) w linesp pi 5 tit (b == 1 ? sprintf('$|\langle \phi_a | \psi_%d \rangle |^2\ , \ a = %d$', bidx+1, b) : sprintf('$%d$',b))
p for [b=1:q] dat u 0:(column(2+b*2)**2) w linesp pi 5 tit (b == 1 ? sprintf('$a = %d$', b) : sprintf('$%d$',b))
#sp for [b=1:q] dat u 1:2:(column(2+b*2)**2) w l tit sprintf('$|\left<\psi_%d | \psi \right>|^2$', b)

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
