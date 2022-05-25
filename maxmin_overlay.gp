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

min(x,y) = x < y ? x : y
phasedat = 'min_gap_'
ifrom = 0
ito = 9

iband=1
###
# do stuff
set palette gray
set cbrange [-6:-2]
unset colorbox

set xrange [0:4]
set yrange [0:6]
set y2range [0:1]
set ytics 1 nomirror
set y2tics
set mxtics 5
set mytics 5
set grid front xtics mxtics ytics mytics

set xlabel 't'
set ylabel 'm'
set y2label 'MaxMin'

set title sprintf('ti = %g, mi = %g, mf = %g', ti,mi,mf)

setlabel(chern,x,y,etc) = sprintf('set label "%s" at first %g, first %g front %s', chern, x,y,etc)

eval setlabel('{[-1,2,-1]}', 0.1 , 0.1 , 'rot')
eval setlabel('{[2,-1,-1]}', 0.33, 0.4 , '')
eval setlabel('{[2,-1,-1]}', 2   , 0.5 , '')
eval setlabel('{[1,0,-1]}' , 0.4 , 1.2 , '')
eval setlabel('{[0,1,-1]}' , 0.5 , 2.5 , '')
eval setlabel('{[0,0,0]}'  , 0.5 , 5   , '')
eval setlabel('{[0,0,0]}'  , 1.8 , 3.5 , '')
eval setlabel('{[1,-1,0]}' , 2.5 , 2.5 , '')
eval setlabel('{[2,-2,0]}' , 3   , 1.5 , '')
eval setlabel('{[2,-3,1]}' , 3.55, 2.5 , '')
eval setlabel('{[1,-2,1]}' , 3.5 , 3.6 , '')
eval setlabel('{[1,0,-1]}' , 1.5 , 1.3 , '')
eval setlabel('{[0,-1,1]}' , 3   , 5   , '')

eval setlabel('{[1,0,-1]}' , 0.1 , 1.8 , 'boxed')
eval setlabel('{[0,1,-1]}' , 1   , 5.5 , 'boxed')
eval setlabel('{[2,-4,2]}' , 3.65, 1.8 , 'boxed')

set arrow from first 0.23, first 1.5 to first 0.28, first 0.6  front
set arrow from first 1.1 , first 5.2 to first 1.13, first 4.85 front
set arrow from first 3.8 , first 1.5 to first 3.9 , first 0.35  front 


set object circle at first ti, first mi front fc rgbcolor "red" fs solid size 0.02

p for [i=ifrom:ito] phasedat.i.'.dat' u 1:2:(log(min($3,$5)+1e-16)) w image notit, mf w l lt -1 lc rgbcolor "red" tit 'mf', for [iband=1:3] dat u 1:(column(1+iband)) axis x1y2 w linesp tit sprintf('MaxMin(init band = %d)',iband)

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
