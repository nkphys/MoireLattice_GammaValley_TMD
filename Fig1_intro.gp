set terminal epslatex size 5,4 standalone color colortext 8
set output 'Fig1_intro.tex'

set multiplot
##################################################################

##################################################################


set lmargin at screen 0.16
set rmargin at screen 0.5
set tmargin at screen 0.98
set bmargin at screen 0.08


#set xtics offset 0,2
set xlabel '\fontsize{16}{60}\selectfont$$'
set xlabel offset 17,1.5
set ylabel '\fontsize{16}{60}\selectfont$E(meV)$'
set ylabel offset 4.5,2.5

set key at 1.5,0.06
set key maxrows 8
set key spacing 1
set key width -1.7
set key samplen 1.0

#set label 2 '\fontsize{12}{60}\selectfont$$'
#set label 2 at 0.01,5

#set arrow 1 from 3.3,-1 to 3.3,15 nohead dt 4 lw 2 lc "black"
#set arrow 2 from 2.9,-1 to 2.9,20 nohead dt 4 lw 2 lc "black"


#set arrow 1 from 0.1,11.5 to 0.1,9.5
#set label 6 at 0.2,0.46
#set label 6 '\fontsize{12}{60}\selectfont $\textcolor{black}{(a)\; \textrm{3-orb model, L=16}} $' front


#set label 7 at 2.8,8.0
#set label 7 '\fontsize{16}{60}\selectfont{\color{red} $ ? $}'


#set label 8 at 6,18
#set label 8 '\fontsize{8}{60}\selectfont ZigZag Phase'


#set label 7 at 0.2,0.75
#set label 7 '\fontsize{7}{60}\selectfont$$'

#set xlabel offset 1.35,1.1
#set xlabel '\fontsize{8}{60}\selectfont$i$'
#set format x ""
#set format y "%1.0f"
#set xtics ('0.1' 0.1, '' 0.2, ''0.3, '' 0.4, '0.5' 0.5, '' 0.6, '' 0.7, '' 0.8, '' 0.9, '1.0' 1.0, '2.0' 2.0 , '' 3.0, '' 4.0, '5.0' 5.0, '' 6.0, '' 7.0, '' 8.0, '' 9.0, '10.0' 10.0)

#set xr [0.5:10]
#set yr [*:0.08]

#set logscale x

#set xr [0.2:10]


#set ytics 1.0
set xtics offset 0,-0.5
set ytics offset 0.4,0
#set ytics ( '\fontsize{12}{60}\selectfont$0.0$' 0.0 , '\fontsize{12}{60}\selectfont$0.1$' 0.1 ,'\fontsize{12}{60}\selectfont$0.2$' 0.2  ,'\fontsize{12}{60}\selectfont$0.3$' 0.3, '\fontsize{12}{60}\selectfont$0.4$' 0.4  ,'\fontsize{12}{60}\selectfont$0.5$' 0.5)
#set xtics ( '\fontsize{12}{60}\selectfont$0$' 0 , '\fontsize{12}{60}\selectfont$\pi/2$' 1.5707 ,'\fontsize{12}{60}\selectfont$\pi$' 3.14159  ,'\fontsize{12}{60}\selectfont$3\pi/2$' 4.7123850, '\fontsize{12}{60}\selectfont$2\pi$' 6.28318)
#set cbtics ( '\fontsize{12}{60}\selectfont$0$' 0 ,'\fontsize{12}{60}\selectfont$5$' 5  ,'\fontsize{12}{60}\selectfont$10$' 10, '\fontsize{12}{60}\selectfont$15$' 15  ,'\fontsize{12}{60}\selectfont$20$' 20  )

set ytics ( '\fontsize{16}{60}\selectfont-$80$' -80.0,  '\fontsize{16}{60}\selectfont$$' -70.0,'\fontsize{16}{60}\selectfont-$60$' -60.0,  '\fontsize{16}{60}\selectfont$$' -50.0, '\fontsize{16}{60}\selectfont-$40$' -40.0,  '\fontsize{16}{60}\selectfont$$' -30.0,'\fontsize{16}{60}\selectfont-$20$' -20.0,  '\fontsize{16}{60}\selectfont$$' -10.0, '\fontsize{16}{60}\selectfont$0$' 0.0,  '\fontsize{16}{60}\selectfont$$' 10.0,'\fontsize{16}{60}\selectfont$$' 20.0)
set xtics ( '\fontsize{16}{60}\selectfont$\Gamma$' 0.0, '\fontsize{16}{60}\selectfont$K$' 24, '\fontsize{16}{60}\selectfont$M$' 36, '\fontsize{16}{60}\selectfont$\Gamma$' 72 )
#set ytics offset 0.7

#set pm3d map
#set pm3d corners2color c1
set xr [0.0:72]
set yr [-80:4]
#set cbr [0:20]
#set palette define (0 "white", 0.05 "white", 0.15 "#C5E1EF", 0.25 "#9EC9E2", 0.38 "#3C93C2", 0.6 "#0D4A70", 1.0 "black")

#LzLzkw_OBC_lambda0.0_L36_eta0.02.txt  MzMzkw_OBC_lambda0.0_L36_eta0.02.txt
#sp "MzMzkw_OBC_lambda0.45_L16_eta0.05_3orb_Parzen.txt" u 1:2:3 w pm3d ti ''

#p "Bands_energy_MoS2_theta_1.5.txt" u 1:(($4-3)*0.17) w p pt 4 ps 1.5 lc "red" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:(($5-3)*0.17) w p pt 4 ps 1.5 lc "red" ti "",
p "Bands_energy_MoS2_theta_1.5.txt" u 1:6 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:7 w l lw 10 lc "black" ti "",  "Bands_energy_MoS2_theta_1.5.txt" u 1:8 w l lw 10 lc "red" ti "",  "Bands_energy_MoS2_theta_1.5.txt" u 1:9 w l lw 10 lc "red" ti "",  "Bands_energy_MoS2_theta_1.5.txt" u 1:10 w l lw 10 lc "red" ti "",  "Bands_energy_MoS2_theta_1.5.txt" u 1:11 w l lw 10 lc "red" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:12 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:13 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:14 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:15 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:16 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:17 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:18 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:19 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:20 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:21 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:22 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:23 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:24 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:25 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:26 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:27 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:28 w l lw 10 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:29 w l lw 10 lc "black" ti ""

unset label 7
unset label 8
unset label 6
unset arrow 1
unset arrow 2


###################################################3
set lmargin at screen 0.63
set rmargin at screen 0.95
set tmargin at screen 0.98
set bmargin at screen 0.6

#set xtics offset 0,2
set xlabel '\fontsize{16}{60}\selectfont$$'
set xlabel offset 17,1.5
set ylabel '\fontsize{16}{60}\selectfont$E(meV)$'
set ylabel offset 5.5,2.5

set key at 1.5,0.06
set key maxrows 8
set key spacing 1
set key width -1.7
set key samplen 1.0

#set label 2 '\fontsize{12}{60}\selectfont$$'
#set label 2 at 0.01,5

#set arrow 1 from 3.3,-1 to 3.3,15 nohead dt 4 lw 2 lc "black"
#set arrow 2 from 2.9,-1 to 2.9,20 nohead dt 4 lw 2 lc "black"


#set arrow 1 from 48,0.6 to 58,0.6 nohead lw 8 lc "black" front 

#set arrow 1 from 0.1,11.5 to 0.1,9.5
#set label 6 at 0.2,0.46
#set label 6 '\fontsize{12}{60}\selectfont $\textcolor{black}{(a)\; \textrm{3-orb model, L=16}} $' front


#set label 7 at 2.8,8.0
#set label 7 '\fontsize{16}{60}\selectfont{\color{red} $ ? $}'


#set label 8 at 6,18
#set label 8 '\fontsize{8}{60}\selectfont ZigZag Phase'


#set label 7 at 0.2,0.75
#set label 7 '\fontsize{7}{60}\selectfont$$'

#set xlabel offset 1.35,1.1
#set xlabel '\fontsize{8}{60}\selectfont$i$'
#set format x ""
#set format y "%1.0f"
#set xtics ('0.1' 0.1, '' 0.2, ''0.3, '' 0.4, '0.5' 0.5, '' 0.6, '' 0.7, '' 0.8, '' 0.9, '1.0' 1.0, '2.0' 2.0 , '' 3.0, '' 4.0, '5.0' 5.0, '' 6.0, '' 7.0, '' 8.0, '' 9.0, '10.0' 10.0)

#set xr [0.5:10]
#set yr [*:0.08]

#set logscale x

#set xr [0.2:10]


#set ytics 1.0
set xtics offset 0,-0.5
set ytics offset 0.4,0
#set ytics ( '\fontsize{12}{60}\selectfont$0.0$' 0.0 , '\fontsize{12}{60}\selectfont$0.1$' 0.1 ,'\fontsize{12}{60}\selectfont$0.2$' 0.2  ,'\fontsize{12}{60}\selectfont$0.3$' 0.3, '\fontsize{12}{60}\selectfont$0.4$' 0.4  ,'\fontsize{12}{60}\selectfont$0.5$' 0.5)
#set xtics ( '\fontsize{12}{60}\selectfont$0$' 0 , '\fontsize{12}{60}\selectfont$\pi/2$' 1.5707 ,'\fontsize{12}{60}\selectfont$\pi$' 3.14159  ,'\fontsize{12}{60}\selectfont$3\pi/2$' 4.7123850, '\fontsize{12}{60}\selectfont$2\pi$' 6.28318)
#set cbtics ( '\fontsize{12}{60}\selectfont$0$' 0 ,'\fontsize{12}{60}\selectfont$5$' 5  ,'\fontsize{12}{60}\selectfont$10$' 10, '\fontsize{12}{60}\selectfont$15$' 15  ,'\fontsize{12}{60}\selectfont$20$' 20  )

set ytics ( '\fontsize{16}{60}\selectfont$$' -1.5,  '\fontsize{16}{60}\selectfont-$1$' -1,'\fontsize{16}{60}\selectfont$$' -0.5,  '\fontsize{16}{60}\selectfont$0$' 0.0, '\fontsize{16}{60}\selectfont$$' 0.5,  '\fontsize{16}{60}\selectfont$1$' 1.0,'\fontsize{16}{60}\selectfont$$' 1.5)

set xtics ( '\fontsize{16}{60}\selectfont$\Gamma$' 0.0, '\fontsize{16}{60}\selectfont$K$' 32, '\fontsize{16}{60}\selectfont$M$' 48, '\fontsize{16}{60}\selectfont$\Gamma$' 96 )
#set ytics offset 0.7

#set pm3d map
#set pm3d corners2color c1
set xr [0.0:96]
set yr [-1.5:1.5]
#set cbr [0:20]
#set palette define (0 "white", 0.05 "white", 0.15 "#C5E1EF", 0.25 "#9EC9E2", 0.38 "#3C93C2", 0.6 "#0D4A70", 1.0 "black")

#LzLzkw_OBC_lambda0.0_L36_eta0.02.txt  MzMzkw_OBC_lambda0.0_L36_eta0.02.txt
#sp "MzMzkw_OBC_lambda0.45_L16_eta0.05_3orb_Parzen.txt" u 1:2:3 w pm3d ti ''

p "Bands_energy_MoS2_theta_1.5.txt" u 1:6 w l lw 16 lc "red" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:7 w l lw 16 lc "red" ti "",  "Bands_energy_MoS2_theta_1.5.txt" u 1:(($4-3)*0.17) w l dt 3 lw 8 lc "black" ti "", "Bands_energy_MoS2_theta_1.5.txt" u 1:(($5-3)*0.17) w l dt 3 lw 8 lc "black" ti ""


unset label 7
unset label 8
unset label 6
unset arrow 1
unset arrow 2


####################################################################################
set lmargin at screen 0.63
set rmargin at screen 0.95
set tmargin at screen 0.45
set bmargin at screen 0.10


set xlabel '\fontsize{16}{60}\selectfont$x(a_{m})$'
set xlabel offset 8,0.4
set ylabel '\fontsize{16}{60}\selectfont$y(a_{m})$'
set ylabel offset 5,4

set key at 0.26,29
set key maxrows 8
set key spacing 1.5
set key width -1.7
set key samplen 1.0

#set label 2 '\fontsize{12}{60}\selectfont$$'
#set label 2 at 0.01,5

#set arrow 1 from 3.3,-1 to 3.3,15 nohead dt 4 lw 2 lc "black"
#set arrow 2 from 2.9,-1 to 2.9,20 nohead dt 4 lw 2 lc "black"


#set arrow 1 from 0.1,11.5 to 0.1,9.5
#set label 6 at 0.0155555,27.5
#set label 6 '\fontsize{12}{60}\selectfont $(c)$' front


set label 7 at 1.5,3
set label 7 '\fontsize{16}{60}\selectfont{$\Delta(r)$}'

set arrow 146 from -2.5,2.0 to -2.30940107564,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 106 from -2.5,-2.0 to -2.30940107564,-2.0 nohead dt 1 lw 2 lc "black" front
set arrow 116 from -2.5,-1.0 to -2.30940107564,-1.0 nohead dt 1 lw 2 lc "black" front
set arrow 126 from -2.5,0 to -2.30940107564,0 nohead dt 1 lw 2 lc "black" front
set arrow 136 from -2.5,1.0 to -2.30940107564,1.0 nohead dt 1 lw 2 lc "black" front
set arrow 146 from -2.5,2.0 to -2.30940107564,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 367 from -2.02072594132,2.5 to -2.309401075910,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 537 from -2.02072594132,-2.5 to -2.309401075910,-2.0 nohead dt 1 lw 2 lc "black" front
set arrow 107 from -2.02072594132,-1.5 to -1.44337567214,-1.5 nohead dt 1 lw 2 lc "black" front
set arrow 327 from -2.02072594132,-1.5 to -2.309401075910,-2.0 nohead dt 1 lw 2 lc "black" front
set arrow 547 from -2.02072594132,-1.5 to -2.309401075910,-1.0 nohead dt 1 lw 2 lc "black" front
set arrow 117 from -2.02072594132,-.5 to -1.44337567214,-.5 nohead dt 1 lw 2 lc "black" front
set arrow 337 from -2.02072594132,-.5 to -2.309401075910,-1.0 nohead dt 1 lw 2 lc "black" front
set arrow 557 from -2.02072594132,-.5 to -2.309401075910,0 nohead dt 1 lw 2 lc "black" front
set arrow 127 from -2.02072594132,.5 to -1.44337567214,.5 nohead dt 1 lw 2 lc "black" front
set arrow 347 from -2.02072594132,.5 to -2.309401075910,0 nohead dt 1 lw 2 lc "black" front
set arrow 567 from -2.02072594132,.5 to -2.309401075910,1.0 nohead dt 1 lw 2 lc "black" front
set arrow 137 from -2.02072594132,1.5 to -1.44337567214,1.5 nohead dt 1 lw 2 lc "black" front
set arrow 357 from -2.02072594132,1.5 to -2.309401075910,1.0 nohead dt 1 lw 2 lc "black" front
set arrow 577 from -2.02072594132,1.5 to -2.309401075910,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 367 from -2.02072594132,2.5 to -2.309401075910,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 98 from -1.15470053782,-2.0 to -.57735026864,-2.0 nohead dt 1 lw 2 lc "black" front
set arrow 318 from -1.15470053782,-2.0 to -1.443375672410,-2.5 nohead dt 1 lw 2 lc "black" front
set arrow 538 from -1.15470053782,-2.0 to -1.443375672410,-1.5 nohead dt 1 lw 2 lc "black" front
set arrow 108 from -1.15470053782,-1.0 to -.57735026864,-1.0 nohead dt 1 lw 2 lc "black" front
set arrow 328 from -1.15470053782,-1.0 to -1.443375672410,-1.5 nohead dt 1 lw 2 lc "black" front
set arrow 548 from -1.15470053782,-1.0 to -1.443375672410,-.5 nohead dt 1 lw 2 lc "black" front
set arrow 118 from -1.15470053782,0 to -.57735026864,0 nohead dt 1 lw 2 lc "black" front
set arrow 338 from -1.15470053782,0 to -1.443375672410,-.5 nohead dt 1 lw 2 lc "black" front
set arrow 558 from -1.15470053782,0 to -1.443375672410,.5 nohead dt 1 lw 2 lc "black" front
set arrow 128 from -1.15470053782,1.0 to -.57735026864,1.0 nohead dt 1 lw 2 lc "black" front
set arrow 348 from -1.15470053782,1.0 to -1.443375672410,.5 nohead dt 1 lw 2 lc "black" front
set arrow 568 from -1.15470053782,1.0 to -1.443375672410,1.5 nohead dt 1 lw 2 lc "black" front
set arrow 138 from -1.15470053782,2.0 to -.57735026864,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 358 from -1.15470053782,2.0 to -1.443375672410,1.5 nohead dt 1 lw 2 lc "black" front
set arrow 578 from -1.15470053782,2.0 to -1.443375672410,2.5 nohead dt 1 lw 2 lc "black" front
set arrow 529 from -.28867513432,-2.5 to -.577350268910,-2.0 nohead dt 1 lw 2 lc "black" front
set arrow 99 from -.28867513432,-1.5 to .28867513486,-1.5 nohead dt 1 lw 2 lc "black" front
set arrow 319 from -.28867513432,-1.5 to -.577350268910,-2.0 nohead dt 1 lw 2 lc "black" front
set arrow 539 from -.28867513432,-1.5 to -.577350268910,-1.0 nohead dt 1 lw 2 lc "black" front
set arrow 109 from -.28867513432,-.5 to .28867513486,-.5 nohead dt 1 lw 2 lc "black" front
set arrow 329 from -.28867513432,-.5 to -.577350268910,-1.0 nohead dt 1 lw 2 lc "black" front
set arrow 549 from -.28867513432,-.5 to -.577350268910,0 nohead dt 1 lw 2 lc "black" front
set arrow 119 from -.28867513432,.5 to .28867513486,.5 nohead dt 1 lw 2 lc "black" front
set arrow 339 from -.28867513432,.5 to -.577350268910,0 nohead dt 1 lw 2 lc "black" front
set arrow 559 from -.28867513432,.5 to -.577350268910,1.0 nohead dt 1 lw 2 lc "black" front
set arrow 129 from -.28867513432,1.5 to .28867513486,1.5 nohead dt 1 lw 2 lc "black" front
set arrow 349 from -.28867513432,1.5 to -.577350268910,1.0 nohead dt 1 lw 2 lc "black" front
set arrow 569 from -.28867513432,1.5 to -.577350268910,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 359 from -.28867513432,2.5 to -.577350268910,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 90 from .57735026918,-2.0 to 1.15470053836,-2.0 nohead dt 1 lw 2 lc "black" front
set arrow 310 from .57735026918,-2.0 to .288675134590,-2.5 nohead dt 1 lw 2 lc "black" front
set arrow 530 from .57735026918,-2.0 to .288675134590,-1.5 nohead dt 1 lw 2 lc "black" front
set arrow 100 from .57735026918,-1.0 to 1.15470053836,-1.0 nohead dt 1 lw 2 lc "black" front
set arrow 320 from .57735026918,-1.0 to .288675134590,-1.5 nohead dt 1 lw 2 lc "black" front
set arrow 540 from .57735026918,-1.0 to .288675134590,-.5 nohead dt 1 lw 2 lc "black" front
set arrow 110 from .57735026918,0 to 1.15470053836,0 nohead dt 1 lw 2 lc "black" front
set arrow 330 from .57735026918,0 to .288675134590,-.5 nohead dt 1 lw 2 lc "black" front
set arrow 550 from .57735026918,0 to .288675134590,.5 nohead dt 1 lw 2 lc "black" front
set arrow 120 from .57735026918,1.0 to 1.15470053836,1.0 nohead dt 1 lw 2 lc "black" front
set arrow 340 from .57735026918,1.0 to .288675134590,.5 nohead dt 1 lw 2 lc "black" front
set arrow 560 from .57735026918,1.0 to .288675134590,1.5 nohead dt 1 lw 2 lc "black" front
set arrow 130 from .57735026918,2.0 to 1.15470053836,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 350 from .57735026918,2.0 to .288675134590,1.5 nohead dt 1 lw 2 lc "black" front
set arrow 570 from .57735026918,2.0 to .288675134590,2.5 nohead dt 1 lw 2 lc "black" front
set arrow 521 from 1.44337567268,-2.5 to 1.154700538090,-2.0 nohead dt 1 lw 2 lc "black" front
set arrow 91 from 1.44337567268,-1.5 to 2.02072594186,-1.5 nohead dt 1 lw 2 lc "black" front
set arrow 311 from 1.44337567268,-1.5 to 1.154700538090,-2.0 nohead dt 1 lw 2 lc "black" front
set arrow 531 from 1.44337567268,-1.5 to 1.154700538090,-1.0 nohead dt 1 lw 2 lc "black" front
set arrow 101 from 1.44337567268,-.5 to 2.02072594186,-.5 nohead dt 1 lw 2 lc "black" front
set arrow 321 from 1.44337567268,-.5 to 1.154700538090,-1.0 nohead dt 1 lw 2 lc "black" front
set arrow 541 from 1.44337567268,-.5 to 1.154700538090,0 nohead dt 1 lw 2 lc "black" front
set arrow 111 from 1.44337567268,.5 to 2.02072594186,.5 nohead dt 1 lw 2 lc "black" front
set arrow 331 from 1.44337567268,.5 to 1.154700538090,0 nohead dt 1 lw 2 lc "black" front
set arrow 551 from 1.44337567268,.5 to 1.154700538090,1.0 nohead dt 1 lw 2 lc "black" front
set arrow 121 from 1.44337567268,1.5 to 2.02072594186,1.5 nohead dt 1 lw 2 lc "black" front
set arrow 341 from 1.44337567268,1.5 to 1.154700538090,1.0 nohead dt 1 lw 2 lc "black" front
set arrow 561 from 1.44337567268,1.5 to 1.154700538090,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 351 from 1.44337567268,2.5 to 1.154700538090,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 82 from 2.30940107618,-2.0 to 2.5,-2.0 nohead dt 1 lw 2 lc "black" front
set arrow 302 from 2.30940107618,-2.0 to 2.020725941590,-2.5 nohead dt 1 lw 2 lc "black" front
set arrow 522 from 2.30940107618,-2.0 to 2.020725941590,-1.5 nohead dt 1 lw 2 lc "black" front
set arrow 92 from 2.30940107618,-1.0 to 2.5,-1.0 nohead dt 1 lw 2 lc "black" front
set arrow 312 from 2.30940107618,-1.0 to 2.020725941590,-1.5 nohead dt 1 lw 2 lc "black" front
set arrow 532 from 2.30940107618,-1.0 to 2.020725941590,-.5 nohead dt 1 lw 2 lc "black" front
set arrow 102 from 2.30940107618,0 to 2.5,0 nohead dt 1 lw 2 lc "black" front
set arrow 322 from 2.30940107618,0 to 2.020725941590,-.5 nohead dt 1 lw 2 lc "black" front
set arrow 542 from 2.30940107618,0 to 2.020725941590,.5 nohead dt 1 lw 2 lc "black" front
set arrow 112 from 2.30940107618,1.0 to 2.5,1.0 nohead dt 1 lw 2 lc "black" front
set arrow 332 from 2.30940107618,1.0 to 2.020725941590,.5 nohead dt 1 lw 2 lc "black" front
set arrow 552 from 2.30940107618,1.0 to 2.020725941590,1.5 nohead dt 1 lw 2 lc "black" front
set arrow 122 from 2.30940107618,2.0 to 2.5,2.0 nohead dt 1 lw 2 lc "black" front
set arrow 342 from 2.30940107618,2.0 to 2.020725941590,1.5 nohead dt 1 lw 2 lc "black" front
set arrow 562 from 2.30940107618,2.0 to 2.020725941590,2.5 nohead dt 1 lw 2 lc "black" front

#set label 8 at 6,18
#set label 8 '\fontsize{8}{60}\selectfont ZigZag Phase'


#set label 7 at 0.2,0.75
#set label 7 '\fontsize{7}{60}\selectfont$$'

#set xlabel offset 1.35,1.1
#set xlabel '\fontsize{8}{60}\selectfont$i$'
#set format x ""
#set format y "%1.0f"
#set xtics ('0.1' 0.1, '' 0.2, ''0.3, '' 0.4, '0.5' 0.5, '' 0.6, '' 0.7, '' 0.8, '' 0.9, '1.0' 1.0, '2.0' 2.0 , '' 3.0, '' 4.0, '5.0' 5.0, '' 6.0, '' 7.0, '' 8.0, '' 9.0, '10.0' 10.0)

#set xr [0.5:10]
#set yr [*:0.08]

#set logscale x

#set xr [0.2:10]

#set ytics 1.0
set xtics offset 0,0.5
set ytics offset 1.0,0


set ytics ( '\fontsize{12}{60}\selectfont-$2$' -2 , '\fontsize{12}{60}\selectfont-$1$'-1 ,'\fontsize{12}{60}\selectfont$0$' 0  ,'\fontsize{12}{60}\selectfont$1$' 1, '\fontsize{12}{60}\selectfont$2$' 2)
set xtics ( '\fontsize{12}{60}\selectfont-$2$' -2 , '\fontsize{12}{60}\selectfont-$1$'-1 ,'\fontsize{12}{60}\selectfont$0$' 0  ,'\fontsize{12}{60}\selectfont$1$' 1, '\fontsize{12}{60}\selectfont$2$' 2)
#set ytics offset 0.7

set xr [-2.5:2.5]
set yr [-2.5:2.5]

#LzLzkw_OBC_lambda0.0_L36_eta0.02.txt  MzMzkw_OBC_lambda0.0_L36_eta0.02.txt
#sp "SzSzkw_OBC_lambda0.0_L36_eta0.02.txt" u 1:2:3 w pm3d ti ''
#Cut_q0_3orb_lambda0.45_L16_eta0.05.txt   Cut_q0_SL_lambda0.45_L48_eta0.02.txt     Cut_qPi_3orb_lambda0.45_L16_eta0.05.txt  Cut_qPi_SL_lambda0.45_L48_eta0.02.txt
 
#p "Cut_q0_3orb_lambda0.45_L16_eta0.05.txt" u 2:3 w lp pt 4 ps 1 dt 4 lw 10 lc "#3C93C2" ti '\fontsize{12}{60}\selectfont $\textrm{3-orb},\; q=0$',  "Cut_qPi_3orb_lambda0.45_L16_eta0.05.txt" u 2:3 w lp pt 5 ps 1 dt 5 lw 10 lc "black" ti '\fontsize{12}{60}\selectfont $\textrm{3-orb},\; q=\pi$',  "Cut_q0_SL_lambda0.45_L48_eta0.02.txt" u 2:3 w lp pt 6 ps 1 dt 3 lw 10 lc "#31a354" ti '\fontsize{12}{60}\selectfont $\textrm{S-L},\; q=0$',  "Cut_qPi_SL_lambda0.45_L48_eta0.02.txt" u 2:3 w lp pt 9 ps 2 dt 3 lw 10 lc "#f03b20" ti '\fontsize{12}{60}\selectfont $\textrm{S-L},\; q=\pi$'


set pm3d map
set pm3d corners2color c1
set palette define (0.0 "#5082FF",  0.78 "white", 0.91 "#FF7050", 0.97 "#FFCF50",  1.0 "#FFF250")
#set palette define (0.0 "#5082FF",  0.78 "white", 0.91 "#FF7050", 0.97 "#FFCF50",  1.0 "#FFF250")

#set palette define (0.0 "white", 0.85 "white", 1.0 "blue")
#set arrow 1 from .577350269,0 to 1.15470053837,0 nohead dt 1 lw 1 lc "black" front
#set arrow 2 from .577350269,0 to .577350269,1 nohead dt 1 lw 1 lc "black" front
#set arrow 3 from .577350269,0 to 1.44337567278,0.5 nohead dt 1 lw 1 lc "black" front

set colorbox horizontal
set colorbox user origin 0.63,0.46 size 0.23,0.02
unset cbtics
set cbtics offset 0,2.3

set cbtics ('\fontsize{8}{60}\selectfont-$300$' -300,'\fontsize{8}{60}\selectfont-$200$' -200 ,'\fontsize{8}{60}\selectfont-$100$' -100,   '\fontsize{8}{60}\selectfont$0$' 0, '\fontsize{8}{60}\selectfont$100$' 100 )

set cbr [-350:100]
#sp "Moire_Potential_rough.txt" u ($1/125.86):($2/125.86):3 w pm3d ti ""

unset multiplot

