unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')

output='t_vs_erg'
if (tex == 1) {
  # set lmargin at screen 0
  # set rmargin at screen 1
  # set bmargin at screen 0
  # set tmargin at screen 1
  set term lua tikz standalone createstyle preamble '\usepackage{amsmath}' size 6in,3in  #fontscale 0.6
  set output output.'.tex'
}

fn = 'erg_vs_t'
dat = fn.'.dat'
par = fn.'.par'
load par

## p/q = 1/4
#tclist = "-0.268 -1 -3.73 0.268 1 3.73"

## p/q = 1/3
#tclist = '-0.924 -0.707 -0.382 0 0.382 0.707 0.924'

## p/q = 2/5
#tclist = '-1.618 -1.023 -0.744 -0.618 -0.414 0.414 0.618 0.744 1.023 1.618'

## p/q = 1/5
#tclist = '-1.618 -0.618 -0.432 0.432 -0.213 0.213 0.618 1.618'
#tclist2 = '0.685 -0.685'

## p/q=1/3, triangular
tclist = '-0.5 0.5'

## p/q=1/4, triangular
#tclist = '0.382 -0.382 0.924 -0.924 0.708 -0.708 0'

xmax=10
dx = 2
ymax=10
###
# do stuff
#set xlabel '$t$'
set xlabel '$E$'
set multiplot layout 1,3

set lmargin 5
set rmargin 0

set xrange [-xmax:xmax]
set xtics -xmax,dx,xmax
set mxtics dx

#set xtics format ''
#set ylabel '$t$'
set yrange [-ymax:ymax]
set ytics nomirror -5,1,5
set ytics add ('$t\quad 0$' 0)
unset y2tics
set for [tc in tclist] arrow nohead lt 0 back from first -xmax,first tc to first xmax,first tc
#set for [tc in tclist2] arrow nohead lt 0 back from first -xmax,first tc to first xmax,first tc

set y2tics -5,0,-5
set for [tc in tclist] y2tics add ('' tc)
#set for [tc in tclist2] y2tics add ('' tc)

# k = (0,0)
#set label 1024 front  '$\boldsymbol{k} = (0,0)$' at first -0.8,first 3.5
set title '$\boldsymbol{k} = (0,0)$'
p for [b=1:q] dat every :::0::0 u (column(3+b)):3 w l lt 2 notit

set ytics auto
set ytics nomirror format ''
#set label 1024 front  '$\boldsymbol{k} = (0,\pi)$ and $(\pi/q,0)$' at first -0.8,first 3.5
set title '$\boldsymbol{k} = (0,\pi)$ and $(\pi/q,0)$'
p for [b=1:q] dat every :::1::1 u (column(3+b)):3 w l lt 2 notit

#set title '$\boldsymbol{k} = (\pi/q,0)$'# and $(\pi/q,0)$'
#p for [b=1:q] dat every :::2::2 u 3:(column(3+b)) w l lt 2 notit

#set label 1024 front  '$\boldsymbol{k} = (\pi/q,\pi)$' at first -0.8,first 3.5
set title '$\boldsymbol{k} = (\pi/q,\pi)$'
#set xlabel '$t$'
#set xtics format '%g'
#set for [tc in tclist] y2tics add (sprintf('\scriptsize{$%s$}',tc) tc)
set for [tc in tclist] label sprintf('\scriptsize{$%s$}',tc) at first xmax+0.4, first tc
#set for [tc in tclist2] label sprintf('\scriptsize{$%s$}',tc) at first xmax+0.4, first tc*1.1 #y2tics add (sprintf('\scriptsize{$\quad %s$}',tc) tc)

# show Chern numbers
set label '\scriptsize{$[-1,-1,4,-6,4]$}' at first 8, first -1.9
set label '\scriptsize{$[-1,-1,4,-6,4]$}' at first 8, first -1.2
set label '\scriptsize{$[-1,-1,-1,-1,4]$}' at first 8, first -0.7
set label '\scriptsize{$[-1,-1,-1,-1,4]$}' at first 8, first -0.5
set label '\scriptsize{$[-1,-1,-1,4,-1]$}' at first 8, first -0.3
set label '\scriptsize{$[-1,-1,4,-1,-1]$}' at first 8, first 0
set label '\scriptsize{$[4,-6,4,-1,-1]$}' at first 8, first 1.9
set label '\scriptsize{$[4,-6,4,-1,-1]$}' at first 8, first 1.2
set label '\scriptsize{$[4,-1,-1,-1,-1]$}' at first 8, first 0.7
set label '\scriptsize{$[4,-1,-1,-1,-1]$}' at first 8, first 0.5
set label '\scriptsize{$[-1,4,-1,-1,-1]$}' at first 8, first 0.3


p for [b=1:q] dat every :::3::3 u (column(3+b)):3 w l lt 2 notit



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
