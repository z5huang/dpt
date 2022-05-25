unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')

output='min_gap'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 4in,2in  #fontscale 0.6
  set output output.'.tex'
}

fn='min_gap'
dat=fn.'.dat'
par=fn.'.par'
load par

set xrange [tmin:tmax]
set yrange [mmin:mmax]
set xlabel '$t$'
set ylabel '$m$'
set xtics 0.5
set ytics 1
set mxtics
set mytics
set grid front xtics ytics mxtics mytics
min(x,y) = x < y ? x : y

set palette gray
set cbrange [-6:-2]
set cbtics 1
unset colorbox

setlabel(chern,x,y,etc) = sprintf('set label "%s" at first %g, first %g front %s', chern, x,y,etc)

eval setlabel('\\tiny{[-1,2,-1]}', 0.1 , 0.1 , 'rot')
eval setlabel('\\tiny{[2,-1,-1]}', 0.33, 0.4 , '')
eval setlabel('\\tiny{[2,-1,-1]}', 2   , 0.5 , '')
eval setlabel('\\tiny{[1,0,-1]}' , 0.4 , 1.2 , '')
eval setlabel('\\tiny{[0,1,-1]}' , 0.5 , 2.5 , '')
eval setlabel('\\tiny{[0,0,0]}'  , 0.5 , 5   , '')
eval setlabel('\\tiny{[0,0,0]}'  , 1.8 , 3.5 , '')
eval setlabel('\\tiny{[1,-1,0]}' , 2.5 , 2.5 , '')
eval setlabel('\\tiny{[2,-2,0]}' , 3   , 1.5 , '')
eval setlabel('\\tiny{[2,-3,1]}' , 3.55, 2.5 , '')
eval setlabel('\\tiny{[1,-2,1]}' , 3.5 , 3.6 , '')
eval setlabel('\\tiny{[1,0,-1]}' , 1.5 , 1.3 , '')
eval setlabel('\\tiny{[0,-1,1]}' , 3   , 5   , '')

eval setlabel('\\tiny{[1,0,-1]}' , 0.1 , 1.8 , 'boxed')
eval setlabel('\\tiny{[0,1,-1]}' , 1   , 5.5 , 'boxed')
eval setlabel('\\tiny{[2,-4,2]}' , 3.65, 1.8 , 'boxed')

set arrow from first 0.23, first 1.5 to first 0.28, first 0.6  front
set arrow from first 1.1 , first 5.2 to first 1.13, first 4.85 front
set arrow from first 3.8 , first 1.5 to first 3.9 , first 0.35  front 

idatmin = 0
idatmax = 9
# NB: the +1e-16 below is to ensure that log() returns a real number
#p dat u 1:2:(log(min($3,$5)+1e-16)) w image notit, for [i=idatmin:idatmax] fn.'_'.i.'.dat' u 1:2:(log(min($3,$5)+1e-16)) w image notit
p for [i=idatmin:idatmax] fn.'_'.i.'.dat' u 1:2:(log(min($3,$5)+1e-16)) w image notit

if (tex == 1){
  unset output
  set term wxt
  build = buildtex(output)
  print '"eval build" to build and preview'
} else {
  #print "press enter to continue"
  #pause -1
}
