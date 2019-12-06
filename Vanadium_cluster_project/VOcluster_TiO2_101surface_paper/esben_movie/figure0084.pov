#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White}
camera {orthographic
  right -9.91*x up 12.31*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}
light_source {<  2.00,   3.00,  125.00> color White
  area_light <0.95, 0, 0>, <0, 0.80, 0>, 5, 4
  adaptive 1 jitter}

#declare simple = finish {phong 0.7}
#declare pale = finish {ambient .5 diffuse .85 roughness .001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.10 roughness 0.04 }
#declare vmd = finish {ambient .0 diffuse .65 phong 0.1 phong_size 40. specular 0.500 }
#declare jmol = finish {ambient .2 diffuse .6 specular 1 roughness .001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.70 roughness 0.04 reflection 0.15}
#declare mj_mol = finish {ambient .00 diffuse .55 phong 0.0 phong_size 0.specular .250 roughness 0.1 brilliance 0.8 reflection 0.0 }
#declare ase3 = finish {ambient .15 brilliance 2 diffuse .6 metallic specular 1. roughness .001 reflection .0}
#declare glass = finish {ambient .05 diffuse .3 specular 1. roughness .001}
#declare glass2 = finish {ambient .0 diffuse .3 specular 1. reflection .25 roughness .001}
#declare Rcell = 0.100;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      translate LOC}
#end

atom(< -4.97,  -4.33, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #0 
atom(<-10.14,  -6.25, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #1 
atom(<-12.03,  -6.25, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #2 
atom(< -6.86,  -4.33, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #3 
atom(< -3.92,  -4.33, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #4 
atom(< -9.09,  -6.25, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #5 
atom(< -3.58,  -6.25, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #6 
atom(< -8.76,  -4.33, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #7 
atom(< -5.82,  -4.33, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #8 
atom(<-10.99,  -6.25, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #9 
atom(< -7.20,  -6.25, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #10 
atom(<-12.37,  -4.33, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #11 
atom(< -3.63,  -4.34,  -9.19>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #12 
atom(< -8.81,  -6.25,  -9.20>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #13 
atom(<-10.67,  -6.25,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #14 
atom(< -5.50,  -4.33,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #15 
atom(<-12.95,  -4.34,  -7.63>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #16 
atom(< -7.76,  -6.26,  -7.66>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #17 
atom(<-12.52,  -6.25,  -7.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #18 
atom(< -7.34,  -4.33,  -7.73>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #19 
atom(< -4.38,  -4.33,  -6.79>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #20 
atom(< -9.53,  -6.23,  -6.81>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #21 
atom(< -5.85,  -6.25,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #22 
atom(<-11.04,  -4.33,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #23 
atom(<-12.45,  -4.35,  -5.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #24 
atom(< -7.25,  -6.24,  -5.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #25 
atom(< -9.13,  -6.25,  -4.78>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #26 
atom(< -3.96,  -4.33,  -4.73>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #27 
atom(<-11.48,  -4.34,  -4.19>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #28 
atom(< -6.32,  -6.25,  -4.18>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #29 
atom(<-10.96,  -6.22,  -3.96>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #30 
atom(< -5.81,  -4.33,  -3.95>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #31 
atom(<-13.09,  -4.31,  -3.28>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #32 
atom(< -7.94,  -6.26,  -3.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #33 
atom(< -4.32,  -6.24,  -4.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #34 
atom(< -9.49,  -4.33,  -4.75>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #35 
atom(< -4.97,  -0.51, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #36 
atom(<-10.14,  -2.42, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #37 
atom(<-12.03,  -2.42, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #38 
atom(< -6.86,  -0.51, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #39 
atom(< -3.92,  -0.51, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #40 
atom(< -9.09,  -2.42, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #41 
atom(< -3.58,  -2.42, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #42 
atom(< -8.76,  -0.51, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #43 
atom(< -5.82,  -0.51, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #44 
atom(<-10.99,  -2.42, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #45 
atom(< -7.20,  -2.42, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #46 
atom(<-12.37,  -0.51, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #47 
atom(< -3.62,  -0.51,  -9.16>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #48 
atom(< -8.81,  -2.42,  -9.20>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #49 
atom(<-10.67,  -2.42,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #50 
atom(< -5.44,  -0.51,  -8.38>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #51 
atom(<-12.93,  -0.51,  -7.55>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #52 
atom(< -7.75,  -2.41,  -7.66>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #53 
atom(<-12.52,  -2.42,  -7.70>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #54 
atom(< -7.34,  -0.51,  -7.72>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #55 
atom(< -4.38,  -0.51,  -6.77>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #56 
atom(< -9.53,  -2.43,  -6.81>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #57 
atom(< -5.85,  -2.42,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #58 
atom(<-11.03,  -0.51,  -8.29>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #59 
atom(<-12.45,  -0.50,  -5.62>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #60 
atom(< -7.25,  -2.43,  -5.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #61 
atom(< -9.14,  -2.42,  -4.79>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #62 
atom(< -3.91,  -0.51,  -4.70>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #63 
atom(<-11.44,  -0.51,  -3.99>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #64 
atom(< -6.32,  -2.41,  -4.18>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #65 
atom(<-10.97,  -2.45,  -3.99>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #66 
atom(< -5.81,  -0.51,  -3.95>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #67 
atom(<-13.12,  -0.52,  -3.22>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #68 
atom(< -7.95,  -2.38,  -3.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #69 
atom(< -4.33,  -2.42,  -4.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #70 
atom(< -9.52,  -0.51,  -4.77>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #71 
atom(<-12.51,  -2.31,  -0.61>, 0.80, rgb <0.19, 0.31, 0.97>, 0.0, pale) // #72 
atom(<-12.40,  -3.20,  -0.19>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #73 
atom(<-13.13,  -1.80,   0.00>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #74 
atom(<-12.98,  -2.45,  -1.47>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #75 
atom(<-10.48,  -0.59,  -1.96>, 0.80, rgb <0.19, 0.31, 0.97>, 0.0, pale) // #76 
atom(<-11.06,  -1.20,  -1.34>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #77 
atom(< -9.55,  -1.00,  -2.03>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #78 
atom(<-10.38,   0.32,  -1.54>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #79 
atom(< -4.97,   3.31, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #80 
atom(<-10.14,   1.40, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #81 
atom(<-12.03,   1.40, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #82 
atom(< -6.86,   3.31, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #83 
atom(< -3.92,   3.31, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #84 
atom(< -9.09,   1.40, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #85 
atom(< -3.58,   1.40, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #86 
atom(< -8.76,   3.31, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #87 
atom(< -5.82,   3.31, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #88 
atom(<-10.99,   1.40, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #89 
atom(< -7.20,   1.40, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #90 
atom(<-12.37,   3.31, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #91 
atom(< -3.63,   3.31,  -9.19>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #92 
atom(< -8.81,   1.40,  -9.20>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #93 
atom(<-10.67,   1.40,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #94 
atom(< -5.50,   3.31,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #95 
atom(<-12.95,   3.31,  -7.63>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #96 
atom(< -7.76,   1.39,  -7.66>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #97 
atom(<-12.52,   1.40,  -7.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #98 
atom(< -7.34,   3.31,  -7.73>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #99 
atom(< -4.38,   3.31,  -6.79>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #100 
atom(< -9.53,   1.41,  -6.81>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #101 
atom(< -5.85,   1.40,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #102 
atom(<-11.04,   3.31,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #103 
atom(<-12.45,   3.30,  -5.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #104 
atom(< -7.25,   1.40,  -5.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #105 
atom(< -9.13,   1.40,  -4.78>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #106 
atom(< -3.96,   3.32,  -4.73>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #107 
atom(<-11.48,   3.31,  -4.19>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #108 
atom(< -6.32,   1.39,  -4.18>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #109 
atom(<-10.96,   1.42,  -3.96>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #110 
atom(< -5.81,   3.32,  -3.95>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #111 
atom(<-13.09,   3.33,  -3.28>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #112 
atom(< -7.94,   1.38,  -3.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #113 
atom(< -4.32,   1.40,  -4.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #114 
atom(< -9.49,   3.31,  -4.75>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #115 
atom(< -4.97,   7.13, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #116 
atom(<-10.14,   5.22, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #117 
atom(<-12.03,   5.22, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #118 
atom(< -6.86,   7.13, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #119 
atom(< -3.92,   7.13, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #120 
atom(< -9.09,   5.22, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #121 
atom(< -3.58,   5.22, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #122 
atom(< -8.76,   7.13, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #123 
atom(< -5.82,   7.13, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #124 
atom(<-10.99,   5.22, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #125 
atom(< -7.20,   5.22, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #126 
atom(<-12.37,   7.13, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #127 
atom(< -3.62,   7.13,  -9.16>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #128 
atom(< -8.81,   5.22,  -9.20>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #129 
atom(<-10.67,   5.22,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #130 
atom(< -5.44,   7.13,  -8.38>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #131 
atom(<-12.93,   7.13,  -7.55>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #132 
atom(< -7.75,   5.23,  -7.66>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #133 
atom(<-12.52,   5.22,  -7.70>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #134 
atom(< -7.34,   7.13,  -7.72>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #135 
atom(< -4.38,   7.13,  -6.77>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #136 
atom(< -9.53,   5.21,  -6.81>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #137 
atom(< -5.85,   5.22,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #138 
atom(<-11.03,   7.13,  -8.29>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #139 
atom(<-12.45,   7.14,  -5.62>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #140 
atom(< -7.25,   5.22,  -5.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #141 
atom(< -9.14,   5.22,  -4.79>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #142 
atom(< -3.91,   7.14,  -4.70>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #143 
atom(<-11.44,   7.14,  -3.99>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #144 
atom(< -6.32,   5.24,  -4.18>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #145 
atom(<-10.97,   5.20,  -3.99>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #146 
atom(< -5.81,   7.14,  -3.95>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #147 
atom(<-13.12,   7.13,  -3.22>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #148 
atom(< -7.95,   5.26,  -3.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #149 
atom(< -4.33,   5.23,  -4.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #150 
atom(< -9.52,   7.14,  -4.77>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #151 
atom(<-12.51,   5.33,  -0.61>, 0.80, rgb <0.19, 0.31, 0.97>, 0.0, pale) // #152 
atom(<-12.40,   4.44,  -0.19>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #153 
atom(<-13.13,   5.84,   0.00>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #154 
atom(<-12.98,   5.19,  -1.47>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #155 
atom(<-10.48,   7.05,  -1.96>, 0.80, rgb <0.19, 0.31, 0.97>, 0.0, pale) // #156 
atom(<-11.06,   6.45,  -1.34>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #157 
atom(< -9.55,   6.64,  -2.03>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #158 
atom(<-10.38,   7.97,  -1.54>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #159 
atom(<  5.37,  -4.33, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #160 
atom(<  0.20,  -6.25, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #161 
atom(< -1.69,  -6.25, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #162 
atom(<  3.48,  -4.33, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #163 
atom(<  6.42,  -4.33, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #164 
atom(<  1.25,  -6.25, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #165 
atom(<  6.76,  -6.25, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #166 
atom(<  1.59,  -4.33, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #167 
atom(<  4.53,  -4.33, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #168 
atom(< -0.65,  -6.25, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #169 
atom(<  3.14,  -6.25, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #170 
atom(< -2.03,  -4.33, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #171 
atom(<  6.71,  -4.34,  -9.19>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #172 
atom(<  1.53,  -6.25,  -9.20>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #173 
atom(< -0.33,  -6.25,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #174 
atom(<  4.84,  -4.33,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #175 
atom(< -2.61,  -4.34,  -7.63>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #176 
atom(<  2.58,  -6.26,  -7.66>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #177 
atom(< -2.18,  -6.25,  -7.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #178 
atom(<  3.00,  -4.33,  -7.73>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #179 
atom(<  5.96,  -4.33,  -6.79>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #180 
atom(<  0.81,  -6.23,  -6.81>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #181 
atom(<  4.50,  -6.25,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #182 
atom(< -0.69,  -4.33,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #183 
atom(< -2.10,  -4.35,  -5.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #184 
atom(<  3.09,  -6.24,  -5.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #185 
atom(<  1.21,  -6.25,  -4.78>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #186 
atom(<  6.38,  -4.33,  -4.73>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #187 
atom(< -1.14,  -4.34,  -4.19>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #188 
atom(<  4.02,  -6.25,  -4.18>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #189 
atom(< -0.62,  -6.22,  -3.96>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #190 
atom(<  4.54,  -4.33,  -3.95>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #191 
atom(< -2.74,  -4.31,  -3.28>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #192 
atom(<  2.40,  -6.26,  -3.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #193 
atom(<  6.02,  -6.24,  -4.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #194 
atom(<  0.85,  -4.33,  -4.75>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #195 
atom(<  5.37,  -0.51, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #196 
atom(<  0.20,  -2.42, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #197 
atom(< -1.69,  -2.42, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #198 
atom(<  3.48,  -0.51, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #199 
atom(<  6.42,  -0.51, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #200 
atom(<  1.25,  -2.42, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #201 
atom(<  6.76,  -2.42, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #202 
atom(<  1.59,  -0.51, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #203 
atom(<  4.53,  -0.51, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #204 
atom(< -0.65,  -2.42, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #205 
atom(<  3.14,  -2.42, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #206 
atom(< -2.03,  -0.51, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #207 
atom(<  6.72,  -0.51,  -9.16>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #208 
atom(<  1.53,  -2.42,  -9.20>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #209 
atom(< -0.33,  -2.42,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #210 
atom(<  4.90,  -0.51,  -8.38>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #211 
atom(< -2.58,  -0.51,  -7.55>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #212 
atom(<  2.59,  -2.41,  -7.66>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #213 
atom(< -2.18,  -2.42,  -7.70>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #214 
atom(<  3.00,  -0.51,  -7.72>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #215 
atom(<  5.96,  -0.51,  -6.77>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #216 
atom(<  0.81,  -2.43,  -6.81>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #217 
atom(<  4.50,  -2.42,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #218 
atom(< -0.69,  -0.51,  -8.29>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #219 
atom(< -2.10,  -0.50,  -5.62>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #220 
atom(<  3.09,  -2.43,  -5.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #221 
atom(<  1.20,  -2.42,  -4.79>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #222 
atom(<  6.44,  -0.51,  -4.70>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #223 
atom(< -1.10,  -0.51,  -3.99>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #224 
atom(<  4.02,  -2.41,  -4.18>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #225 
atom(< -0.63,  -2.45,  -3.99>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #226 
atom(<  4.53,  -0.51,  -3.95>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #227 
atom(< -2.78,  -0.52,  -3.22>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #228 
atom(<  2.39,  -2.38,  -3.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #229 
atom(<  6.02,  -2.42,  -4.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #230 
atom(<  0.82,  -0.51,  -4.77>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #231 
atom(< -2.16,  -2.31,  -0.61>, 0.80, rgb <0.19, 0.31, 0.97>, 0.0, pale) // #232 
atom(< -2.06,  -3.20,  -0.19>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #233 
atom(< -2.79,  -1.80,   0.00>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #234 
atom(< -2.64,  -2.45,  -1.47>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #235 
atom(< -0.14,  -0.59,  -1.96>, 0.80, rgb <0.19, 0.31, 0.97>, 0.0, pale) // #236 
atom(< -0.72,  -1.20,  -1.34>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #237 
atom(<  0.79,  -1.00,  -2.03>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #238 
atom(< -0.04,   0.32,  -1.54>, 0.30, rgb <1.00, 1.00, 1.00>, 0.0, pale) // #239 
atom(<  5.37,   3.31, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #240 
atom(<  0.20,   1.40, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #241 
atom(< -1.69,   1.40, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #242 
atom(<  3.48,   3.31, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #243 
atom(<  6.42,   3.31, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #244 
atom(<  1.25,   1.40, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #245 
atom(<  6.76,   1.40, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #246 
atom(<  1.59,   3.31, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #247 
atom(<  4.53,   3.31, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #248 
atom(< -0.65,   1.40, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #249 
atom(<  3.14,   1.40, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #250 
atom(< -2.03,   3.31, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #251 
atom(<  6.71,   3.31,  -9.19>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #252 
atom(<  1.53,   1.40,  -9.20>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #253 
atom(< -0.33,   1.40,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #254 
atom(<  4.84,   3.31,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #255 
atom(< -2.61,   3.31,  -7.63>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #256 
atom(<  2.58,   1.39,  -7.66>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #257 
atom(< -2.18,   1.40,  -7.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #258 
atom(<  3.00,   3.31,  -7.73>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #259 
atom(<  5.96,   3.31,  -6.79>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #260 
atom(<  0.81,   1.41,  -6.81>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #261 
atom(<  4.50,   1.40,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #262 
atom(< -0.69,   3.31,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #263 
atom(< -2.10,   3.30,  -5.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #264 
atom(<  3.09,   1.40,  -5.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #265 
atom(<  1.21,   1.40,  -4.78>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #266 
atom(<  6.38,   3.32,  -4.73>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #267 
atom(< -1.14,   3.31,  -4.19>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #268 
atom(<  4.02,   1.39,  -4.18>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #269 
atom(< -0.62,   1.42,  -3.96>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #270 
atom(<  4.54,   3.32,  -3.95>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #271 
atom(< -2.74,   3.33,  -3.28>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #272 
atom(<  2.40,   1.38,  -3.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #273 
atom(<  6.02,   1.40,  -4.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #274 
atom(<  0.85,   3.31,  -4.75>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #275 
atom(<  5.37,   7.13, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #276 
atom(<  0.20,   5.22, -12.78>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #277 
atom(< -1.69,   5.22, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #278 
atom(<  3.48,   7.13, -12.02>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #279 
atom(<  6.42,   7.13, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #280 
atom(<  1.25,   5.22, -11.14>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #281 
atom(<  6.76,   5.22, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #282 
atom(<  1.59,   7.13, -11.27>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #283 
atom(<  4.53,   7.13, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #284 
atom(< -0.65,   5.22, -10.38>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #285 
atom(<  3.14,   5.22, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #286 
atom(< -2.03,   7.13, -11.89>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #287 
atom(<  6.72,   7.13,  -9.16>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #288 
atom(<  1.53,   5.22,  -9.20>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #289 
atom(< -0.33,   5.22,  -8.40>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #290 
atom(<  4.90,   7.13,  -8.38>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #291 
atom(< -2.58,   7.13,  -7.55>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #292 
atom(<  2.59,   5.23,  -7.66>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #293 
atom(< -2.18,   5.22,  -7.70>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #294 
atom(<  3.00,   7.13,  -7.72>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #295 
atom(<  5.96,   7.13,  -6.77>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #296 
atom(<  0.81,   5.21,  -6.81>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #297 
atom(<  4.50,   5.22,  -8.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #298 
atom(< -0.69,   7.13,  -8.29>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #299 
atom(< -2.10,   7.14,  -5.62>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #300 
atom(<  3.09,   5.22,  -5.71>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #301 
atom(<  1.20,   5.22,  -4.79>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #302 
atom(<  6.44,   7.14,  -4.70>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #303 
atom(< -1.10,   7.14,  -3.99>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #304 
atom(<  4.02,   5.24,  -4.18>, 1.30, rgb <0.75, 0.76, 0.78>, 0.0, pale) // #305 
atom(< -0.63,   5.20,  -3.99>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #306 
atom(<  4.53,   7.14,  -3.95>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #307 
atom(< -2.78,   7.13,  -3.22>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #308 
atom(<  2.39,   5.26,  -3.31>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #309 
atom(<  6.02,   5.23,  -4.69>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #310 
atom(<  0.82,   7.14,  -4.77>, 0.74, rgb <0.84, 0.32, 0.33>, 0.0, pale) // #311 
