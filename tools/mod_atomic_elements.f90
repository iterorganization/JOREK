!> Element names by atomic number.
!> Used to automatically generate filenames in some cases
!> Data taken from http://pastebin.com/raw/CKwm136x
!> Converted with `awk -F, '{printf("\'%-4s\', & ! %d, %s\n", $2, $3, $1)}'`
!> Remove extra space with :%s/^' /'/g in vim
module mod_atomic_elements
implicit none

private
public :: element_symbols, atomic_weights

!> NOTE: special solution to indicate isotopes of hydrogen -> Deuterium and Tritium

character(len=3), parameter :: element_symbols(-3:118) = [&
'T  ', & ! -3, Tritium
'D  ', & ! -2, Deuterium
'e  ', & ! -1, Electron
'n  ', & ! 0, Neutron
'H  ', & ! 1, Hydrogen
'He ', & ! 2, Helium
'Li ', & ! 3, Lithium
'Be ', & ! 4, Beryllium
'B  ', & ! 5, Boron
'C  ', & ! 6, Carbon
'N  ', & ! 7, Nitrogen
'O  ', & ! 8, Oxygen
'F  ', & ! 9, Fluorine
'Ne ', & ! 10, Neon
'Na ', & ! 11, Sodium
'Mg ', & ! 12, Magnesium
'Al ', & ! 13, Aluminum
'Si ', & ! 14, Silicon
'P  ', & ! 15, Phosphorus
'S  ', & ! 16, Sulfur
'Cl ', & ! 17, Chlorine
'Ar ', & ! 18, Argon
'K  ', & ! 19, Potassium
'Ca ', & ! 20, Calcium
'Sc ', & ! 21, Scandium
'Ti ', & ! 22, Titanium
'V  ', & ! 23, Vanadium
'Cr ', & ! 24, Chromium
'Mn ', & ! 25, Manganese
'Fe ', & ! 26, Iron
'Co ', & ! 27, Cobalt
'Ni ', & ! 28, Nickel
'Cu ', & ! 29, Copper
'Zn ', & ! 30, Zinc
'Ga ', & ! 31, Gallium
'Ge ', & ! 32, Germanium
'As ', & ! 33, Arsenic
'Se ', & ! 34, Selenium
'Br ', & ! 35, Bromine
'Kr ', & ! 36, Krypton
'Rb ', & ! 37, Rubidium
'Sr ', & ! 38, Strontium
'Y  ', & ! 39, Yttrium
'Zr ', & ! 40, Zirconium
'Nb ', & ! 41, Niobium
'Mo ', & ! 42, Molybdenum
'Tc ', & ! 43, Technetium
'Ru ', & ! 44, Ruthenium
'Rh ', & ! 45, Rhodium
'Pd ', & ! 46, Palladium
'Ag ', & ! 47, Silver
'Cd ', & ! 48, Cadmium
'In ', & ! 49, Indium
'Sn ', & ! 50, Tin
'Sb ', & ! 51, Antimony
'Te ', & ! 52, Tellurium
'I  ', & ! 53, Iodine
'Xe ', & ! 54, Xenon
'Cs ', & ! 55, Cesium
'Ba ', & ! 56, Barium
'La ', & ! 57, Lanthanum
'Ce ', & ! 58, Cerium
'Pr ', & ! 59, Praseodymium
'Nd ', & ! 60, Neodymium
'Pm ', & ! 61, Promethium
'Sm ', & ! 62, Samarium
'Eu ', & ! 63, Europium
'Gd ', & ! 64, Gadolinium
'Tb ', & ! 65, Terbium
'Dy ', & ! 66, Dysprosium
'Ho ', & ! 67, Holmium
'Er ', & ! 68, Erbium
'Tm ', & ! 69, Thulium
'Yb ', & ! 70, Ytterbium
'Lu ', & ! 71, Lutetium
'Hf ', & ! 72, Hafnium
'Ta ', & ! 73, Tantalum
'W  ', & ! 74, Tungsten
'Re ', & ! 75, Rhenium
'Os ', & ! 76, Osmium
'Ir ', & ! 77, Iridium
'Pt ', & ! 78, Platinum
'Au ', & ! 79, Gold
'Hg ', & ! 80, Mercury
'Tl ', & ! 81, Thallium
'Pb ', & ! 82, Lead
'Bi ', & ! 83, Bismuth
'Po ', & ! 84, Polonium
'At ', & ! 85, Astatine
'Rn ', & ! 86, Radon
'Fr ', & ! 87, Francium
'Ra ', & ! 88, Radium
'Ac ', & ! 89, Actinium
'Th ', & ! 90, Thorium
'Pa ', & ! 91, Protactinium
'U  ', & ! 92, Uranium
'Np ', & ! 93, Neptunium
'Pu ', & ! 94, Plutonium
'Am ', & ! 95, Americium
'Cm ', & ! 96, Curium
'Bk ', & ! 97, Berkelium
'Cf ', & ! 98, Californium
'Es ', & ! 99, Einsteinium
'Fm ', & ! 100, Fermium
'Md ', & ! 101, Mendelevium
'No ', & ! 102, Nobelium
'Lr ', & ! 103, Lawrencium
'Rf ', & ! 104, Rutherfordium
'Db ', & ! 105, Dubnium
'Sg ', & ! 106, Seaborgium
'Bh ', & ! 107, Bohrium
'Hs ', & ! 108, Hassium
'Mt ', & ! 109, Meitnerium
'Uun', & ! 110, Ununnilium
'Uuu', & ! 111, Unununium
'Uub', & ! 112, Ununbium
'Uut', & ! 113, Ununtrium
'Uuq', & ! 114, Ununquadium
'Uup', & ! 115, Ununpentium
'Uuh', & ! 116, Ununhexium
'Uus', & ! 117, Ununseptium
'Uuo']   ! 118, Ununoctium

! Atomic weights in atomic mass units (averaged over isotope when not specified per isotope)
real*4, parameter :: atomic_weights(-3:118) = [ &
3.0160492, & ! -3, Tritium
2.01410178, & ! -2, Deuterium
5.48579909e-4, & ! -1, Electron
1.0086649159, & ! 0, Neutron
1.007940, & ! 1, Hydrogen
4.002602, & ! 2, Helium
6.941000, & ! 3, Lithium
9.012180, & ! 4, Beryllium
10.811000, & ! 5, Boron
12.011000, & ! 6, Carbon
14.006740, & ! 7, Nitrogen
15.999400, & ! 8, Oxygen
18.998403, & ! 9, Fluorine
20.179700, & ! 10, Neon
22.989768, & ! 11, Sodium
24.305000, & ! 12, Magnesium
26.981539, & ! 13, Aluminum
28.085500, & ! 14, Silicon
30.973762, & ! 15, Phosphorus
32.066000, & ! 16, Sulfur
35.452700, & ! 17, Chlorine
39.948000, & ! 18, Argon
39.098300, & ! 19, Potassium
40.078000, & ! 20, Calcium
44.955910, & ! 21, Scandium
47.880000, & ! 22, Titanium
50.941500, & ! 23, Vanadium
51.996100, & ! 24, Chromium
54.938050, & ! 25, Manganese
55.847000, & ! 26, Iron
58.933200, & ! 27, Cobalt
58.693400, & ! 28, Nickel
63.546000, & ! 29, Copper
65.390000, & ! 30, Zinc
69.723000, & ! 31, Gallium
72.610000, & ! 32, Germanium
74.921590, & ! 33, Arsenic
78.960000, & ! 34, Selenium
79.904000, & ! 35, Bromine
83.800000, & ! 36, Krypton
85.467800, & ! 37, Rubidium
87.620000, & ! 38, Strontium
88.905850, & ! 39, Yttrium
91.224000, & ! 40, Zirconium
92.906380, & ! 41, Niobium
95.940000, & ! 42, Molybdenum
97.907200, & ! 43, Technetium
101.070000, & ! 44, Ruthenium
102.905500, & ! 45, Rhodium
106.420000, & ! 46, Palladium
107.868200, & ! 47, Silver
112.411000, & ! 48, Cadmium
114.818000, & ! 49, Indium
118.710000, & ! 50, Tin
121.760000, & ! 51, Antimony
127.600000, & ! 52, Tellurium
126.904470, & ! 53, Iodine
131.290000, & ! 54, Xenon
132.905430, & ! 55, Cesium
137.327000, & ! 56, Barium
138.905500, & ! 57, Lanthanum
140.115000, & ! 58, Cerium
140.907650, & ! 59, Praseodymium
144.240000, & ! 60, Neodymium
144.912700, & ! 61, Promethium
150.360000, & ! 62, Samarium
151.965000, & ! 63, Europium
157.250000, & ! 64, Gadolinium
158.925340, & ! 65, Terbium
162.500000, & ! 66, Dysprosium
164.930320, & ! 67, Holmium
167.260000, & ! 68, Erbium
168.934210, & ! 69, Thulium
173.040000, & ! 70, Ytterbium
174.967000, & ! 71, Lutetium
178.490000, & ! 72, Hafnium
180.947900, & ! 73, Tantalum
183.840000, & ! 74, Tungsten
186.207000, & ! 75, Rhenium
190.230000, & ! 76, Osmium
192.220000, & ! 77, Iridium
195.080000, & ! 78, Platinum
196.966540, & ! 79, Gold
200.590000, & ! 80, Mercury
204.383300, & ! 81, Thallium
207.200000, & ! 82, Lead
208.980370, & ! 83, Bismuth
208.982400, & ! 84, Polonium
209.987100, & ! 85, Astatine
222.017600, & ! 86, Radon
223.019700, & ! 87, Francium
226.025400, & ! 88, Radium
227.027800, & ! 89, Actinium
232.038100, & ! 90, Thorium
231.035880, & ! 91, Protactinium
238.028900, & ! 92, Uranium
237.048000, & ! 93, Neptunium
244.064200, & ! 94, Plutonium
243.061400, & ! 95, Americium
247.070300, & ! 96, Curium
247.070300, & ! 97, Berkelium
251.079600, & ! 98, Californium
252.083000, & ! 99, Einsteinium
257.095100, & ! 100, Fermium
258.100000, & ! 101, Mendelevium
259.100900, & ! 102, Nobelium
262.110000, & ! 103, Lawrencium
0.000000, & ! 104, Rutherfordium
0.000000, & ! 105, Dubnium
0.000000, & ! 106, Seaborgium
0.000000, & ! 107, Bohrium
0.000000, & ! 108, Hassium
0.000000, & ! 109, Meitnerium
0.000000, & ! 110, Ununnilium
0.000000, & ! 111, Unununium
0.000000, & ! 112, Ununbium
0.000000, & ! 113, Ununtrium
0.000000, & ! 114, Ununquadium
0.000000, & ! 115, Ununpentium
0.000000, & ! 116, Ununhexium
0.000000, & ! 117, Ununseptium
0.000000]   ! 118, Ununoctium

end module mod_atomic_elements
