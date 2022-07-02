--This file computes Betti tables for P^1 for d = 10 and b = 8
A := degreesRing 2
new HashTable from {
--tb stands for Total Betti numbers
"tb"=>new HashTable from {(7,0) => 240, (6,1) => 0, (7,1) => 0, (8,0) => 45, (8,1) => 0, (9,0) => 0, (9,1) => 1, (0,0) => 9, (0,1) => 0, (1,0) => 80, (2,0) => 315, (1,1) => 0, (3,0) => 720, (2,1) => 0, (4,0) => 1050, (3,1) => 0, (5,0) => 1008, (4,1) => 0, (5,1) => 0, (6,0) => 630},
--mb stands for Multigraded Betti numbers
"mb"=>new HashTable from {(7,0) => A_0^50*A_1^28+2*A_0^49*A_1^29+3*A_0^48*A_1^30+5*A_0^47*A_1^31+7*A_0^46*A_1^32+9*A_0^45*A_1^33+12*A_0^44*A_1^34+15*A_0^43*A_1^35+17*A_0^42*A_1^36+19*A_0^41*A_1^37+20*A_0^40*A_1^38+20*A_0^39*A_1^39+20*A_0^38*A_1^40+19*A_0^37*A_1^41+17*A_0^36*A_1^42+15*A_0^35*A_1^43+12*A_0^34*A_1^44+9*A_0^33*A_1^45+7*A_0^32*A_1^46+5*A_0^31*A_1^47+3*A_0^30*A_1^48+2*A_0^29*A_1^49+A_0^28*A_1^50, (6,1) => 0, (8,0) => A_0^52*A_1^36+A_0^51*A_1^37+2*A_0^50*A_1^38+2*A_0^49*A_1^39+3*A_0^48*A_1^40+3*A_0^47*A_1^41+4*A_0^46*A_1^42+4*A_0^45*A_1^43+5*A_0^44*A_1^44+4*A_0^43*A_1^45+4*A_0^42*A_1^46+3*A_0^41*A_1^47+3*A_0^40*A_1^48+2*A_0^39*A_1^49+2*A_0^38*A_1^50+A_0^37*A_1^51+A_0^36*A_1^52, (7,1) => 0, (9,0) => 0, (8,1) => 0, (9,1) => A_0^54*A_1^54, (0,0) => A_0^8+A_0^7*A_1+A_0^6*A_1^2+A_0^5*A_1^3+A_0^4*A_1^4+A_0^3*A_1^5+A_0^2*A_1^6+A_0*A_1^7+A_1^8, (0,1) => 0, (1,0) => A_0^17*A_1+2*A_0^16*A_1^2+3*A_0^15*A_1^3+4*A_0^14*A_1^4+5*A_0^13*A_1^5+6*A_0^12*A_1^6+7*A_0^11*A_1^7+8*A_0^10*A_1^8+8*A_0^9*A_1^9+8*A_0^8*A_1^10+7*A_0^7*A_1^11+6*A_0^6*A_1^12+5*A_0^5*A_1^13+4*A_0^4*A_1^14+3*A_0^3*A_1^15+2*A_0^2*A_1^16+A_0*A_1^17, (2,0) => A_0^25*A_1^3+2*A_0^24*A_1^4+4*A_0^23*A_1^5+6*A_0^22*A_1^6+9*A_0^21*A_1^7+12*A_0^20*A_1^8+16*A_0^19*A_1^9+19*A_0^18*A_1^10+23*A_0^17*A_1^11+25*A_0^16*A_1^12+27*A_0^15*A_1^13+27*A_0^14*A_1^14+27*A_0^13*A_1^15+25*A_0^12*A_1^16+23*A_0^11*A_1^17+19*A_0^10*A_1^18+16*A_0^9*A_1^19+12*A_0^8*A_1^20+9*A_0^7*A_1^21+6*A_0^6*A_1^22+4*A_0^5*A_1^23+2*A_0^4*A_1^24+A_0^3*A_1^25, (1,1) => 0, (2,1) => 0, (3,0) => A_0^32*A_1^6+2*A_0^31*A_1^7+4*A_0^30*A_1^8+7*A_0^29*A_1^9+11*A_0^28*A_1^10+16*A_0^27*A_1^11+22*A_0^26*A_1^12+29*A_0^25*A_1^13+36*A_0^24*A_1^14+43*A_0^23*A_1^15+49*A_0^22*A_1^16+54*A_0^21*A_1^17+57*A_0^20*A_1^18+58*A_0^19*A_1^19+57*A_0^18*A_1^20+54*A_0^17*A_1^21+49*A_0^16*A_1^22+43*A_0^15*A_1^23+36*A_0^14*A_1^24+29*A_0^13*A_1^25+22*A_0^12*A_1^26+16*A_0^11*A_1^27+11*A_0^10*A_1^28+7*A_0^9*A_1^29+4*A_0^8*A_1^30+2*A_0^7*A_1^31+A_0^6*A_1^32, (3,1) => 0, (4,0) => A_0^38*A_1^10+2*A_0^37*A_1^11+4*A_0^36*A_1^12+7*A_0^35*A_1^13+12*A_0^34*A_1^14+17*A_0^33*A_1^15+25*A_0^32*A_1^16+33*A_0^31*A_1^17+43*A_0^30*A_1^18+52*A_0^29*A_1^19+62*A_0^28*A_1^20+69*A_0^27*A_1^21+77*A_0^26*A_1^22+80*A_0^25*A_1^23+82*A_0^24*A_1^24+80*A_0^23*A_1^25+77*A_0^22*A_1^26+69*A_0^21*A_1^27+62*A_0^20*A_1^28+52*A_0^19*A_1^29+43*A_0^18*A_1^30+33*A_0^17*A_1^31+25*A_0^16*A_1^32+17*A_0^15*A_1^33+12*A_0^14*A_1^34+7*A_0^13*A_1^35+4*A_0^12*A_1^36+2*A_0^11*A_1^37+A_0^10*A_1^38, (4,1) => 0, (5,0) => A_0^43*A_1^15+2*A_0^42*A_1^16+4*A_0^41*A_1^17+7*A_0^40*A_1^18+11*A_0^39*A_1^19+17*A_0^38*A_1^20+24*A_0^37*A_1^21+32*A_0^36*A_1^22+41*A_0^35*A_1^23+50*A_0^34*A_1^24+59*A_0^33*A_1^25+67*A_0^32*A_1^26+73*A_0^31*A_1^27+77*A_0^30*A_1^28+78*A_0^29*A_1^29+77*A_0^28*A_1^30+73*A_0^27*A_1^31+67*A_0^26*A_1^32+59*A_0^25*A_1^33+50*A_0^24*A_1^34+41*A_0^23*A_1^35+32*A_0^22*A_1^36+24*A_0^21*A_1^37+17*A_0^20*A_1^38+11*A_0^19*A_1^39+7*A_0^18*A_1^40+4*A_0^17*A_1^41+2*A_0^16*A_1^42+A_0^15*A_1^43, (6,0) => A_0^47*A_1^21+2*A_0^46*A_1^22+4*A_0^45*A_1^23+6*A_0^44*A_1^24+10*A_0^43*A_1^25+14*A_0^42*A_1^26+20*A_0^41*A_1^27+25*A_0^40*A_1^28+32*A_0^39*A_1^29+37*A_0^38*A_1^30+43*A_0^37*A_1^31+46*A_0^36*A_1^32+50*A_0^35*A_1^33+50*A_0^34*A_1^34+50*A_0^33*A_1^35+46*A_0^32*A_1^36+43*A_0^31*A_1^37+37*A_0^30*A_1^38+32*A_0^29*A_1^39+25*A_0^28*A_1^40+20*A_0^27*A_1^41+14*A_0^26*A_1^42+10*A_0^25*A_1^43+6*A_0^24*A_1^44+4*A_0^23*A_1^45+2*A_0^22*A_1^46+A_0^21*A_1^47, (5,1) => 0},
--sb represents the betti numbers as sums of Schur functors
"sb"=>new HashTable from {(7,0) => {({50,28},1)}, (6,1) => {}, (7,1) => {}, (8,0) => {({52,36},1)}, (8,1) => {}, (9,0) => {}, (9,1) => {({54,54},1)}, (0,0) => {({8,0},1)}, (0,1) => {}, (1,0) => {({17,1},1)}, (2,0) => {({25,3},1)}, (1,1) => {}, (3,0) => {({32,6},1)}, (2,1) => {}, (4,0) => {({38,10},1)}, (3,1) => {}, (5,0) => {({43,15},1)}, (4,1) => {}, (5,1) => {}, (6,0) => {({47,21},1)}},
--dw encodes the dominant weights in each entry
"dw"=>new HashTable from {(7,0) => {{50,28}}, (6,1) => {}, (7,1) => {}, (8,0) => {{52,36}}, (8,1) => {}, (9,0) => {}, (9,1) => {{54,54}}, (0,0) => {{8,0}}, (0,1) => {}, (1,0) => {{17,1}}, (2,0) => {{25,3}}, (1,1) => {}, (3,0) => {{32,6}}, (2,1) => {}, (4,0) => {{38,10}}, (3,1) => {}, (5,0) => {{43,15}}, (4,1) => {}, (5,1) => {}, (6,0) => {{47,21}}},
--lw encodes the lex leading weight in each entry
"lw"=>new HashTable from {(7,0) => {50,28}, (6,1) => {}, (7,1) => {}, (8,0) => {52,36}, (8,1) => {}, (9,0) => {}, (9,1) => {54,54}, (0,0) => {8,0}, (0,1) => {}, (1,0) => {17,1}, (2,0) => {25,3}, (1,1) => {}, (3,0) => {32,6}, (2,1) => {}, (4,0) => {38,10}, (3,1) => {}, (5,0) => {43,15}, (4,1) => {}, (5,1) => {}, (6,0) => {47,21}},
--nr encodes the number of disctinct reprsentations in each entry
"nr"=>new HashTable from {(7,0) => 1, (6,1) => 0, (7,1) => 0, (8,0) => 1, (8,1) => 0, (9,0) => 0, (9,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (3,0) => 1, (2,1) => 0, (4,0) => 1, (3,1) => 0, (5,0) => 1, (4,1) => 0, (5,1) => 0, (6,0) => 1},
--nrm encodes the number of representations with multiplicity in each entry
"nrm"=>new HashTable from {(7,0) => 1, (6,1) => 0, (7,1) => 0, (8,0) => 1, (8,1) => 0, (9,0) => 0, (9,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (3,0) => 1, (2,1) => 0, (4,0) => 1, (3,1) => 0, (5,0) => 1, (4,1) => 0, (5,1) => 0, (6,0) => 1},
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
"er"=>new HashTable from {(7,0) => 240, (6,1) => 0, (7,1) => 0, (8,0) => 45, (8,1) => 0, (9,0) => 0, (9,1) => 1, (0,0) => 9, (0,1) => 0, (1,0) => 80, (2,0) => 315, (1,1) => 0, (3,0) => 720, (2,1) => 0, (4,0) => 1050, (3,1) => 0, (5,0) => 1008, (4,1) => 0, (5,1) => 0, (6,0) => 630},
--bs encodes the Boij-Soederberg coefficients each entry
"bs"=>{3628800/1},
}