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

#set xrange [tmin:tmax]
#set yrange [mmin:mmax]
#set xtics 0.5
#set ytics 0.5
#set mxtics
#set mytics
#set grid front xtics ytics mxtics mytics
min(x,y) = x < y ? x : y

set palette gray
#set cbrange [-6:-2]

#setlabel(chern,x,y) = sprintf('set label "%s" at first %g, first %g front', chern, x,y)

#eval setlabel('[-1,2,-1]', 0.0, 0.4)
#eval setlabel('[2,-1,-1]', 0.4, 0.2)
#eval setlabel('[2,-1,-1]', 2  , 0.5)
#eval setlabel('[1,0,-1]' , 0.5, 1.2)
#eval setlabel('[0,1,-1]' , 0.5, 2.5)
#eval setlabel('[0,0,0]'  , 0.5, 5  )
#eval setlabel('[0,0,0]'  , 1.8, 3.5)
#eval setlabel('[1,-1,0]' , 2.5, 2.5)
#eval setlabel('[2,-2,0]' , 3  , 1.5)
#eval setlabel('[2,-3,1]' , 3.7, 2.4)
#eval setlabel('[1,-2,1]' , 3.5, 3.4)
#eval setlabel('[1,0,-1]' , 1.5, 1.2)
#eval setlabel('[0,-1,1]' , 3  , 5  )

idatmin = 2
idatmax = 9
#p dat u 1:2:(log(min(column(3),column(5)))) w image notit#, for [i=idatmin:idatmax] fn.'_'.i.'.dat' u 1:2:(log(min($3,$5))) w image notit
p dat u 1:2:5 w image notit

###
# do stuff

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
