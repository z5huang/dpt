unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')

output='erg_vs_t'
if (tex == 1) {
  # set lmargin at screen 0
  # set rmargin at screen 1
  # set bmargin at screen 0
  # set tmargin at screen 1
  set term lua tikz standalone createstyle preamble '\usepackage{amsmath}' size 4in,6in  #fontscale 0.6
  set output output.'.tex'
}

fn = 'erg_vs_t'
dat = fn.'.dat'
par = fn.'.par'
load par

#tclist = "-0.268 -1 -3.73 0.268 1 3.73"
tclist = '-0.924 -0.707 -0.382 0 0.382 0.707 0.924'

###
# do stuff
#set xlabel '$t$'
set ylabel '$E$'
set multiplot layout 3,1

set yrange [-5:5]
set ytics -5,2,5
set xtics format ''
set for [tc in tclist] arrow nohead lt 0 back from first tc,first -5 to first tc,first 5

# k = (0,0)
#set label 1024 front  '$\boldsymbol{k} = (0,0)$' at first -0.8,first 3.5
set title '$\boldsymbol{k} = (0,0)$'
p for [b=1:q] dat every :::0::0 u 3:(column(3+b)) w l lt 2 notit

#set label 1024 front  '$\boldsymbol{k} = (0,\pi)$ and $(\pi/q,0)$' at first -0.8,first 3.5
set title '$\boldsymbol{k} = (0,\pi)$ and $(\pi/q,0)$'
p for [b=1:q] dat every :::1::1 u 3:(column(3+b)) w l lt 2 notit

#set title '$\boldsymbol{k} = (\pi/q,0)$'# and $(\pi/q,0)$'
#p for [b=1:q] dat every :::2::2 u 3:(column(3+b)) w l lt 2 notit

#set label 1024 front  '$\boldsymbol{k} = (\pi/q,\pi)$' at first -0.8,first 3.5
set title '$\boldsymbol{k} = (\pi/q,\pi)$'
set xlabel '$t$'
set xtics format '%g'
p for [b=1:q] dat every :::3::3 u 3:(column(3+b)) w l lt 2 notit

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
