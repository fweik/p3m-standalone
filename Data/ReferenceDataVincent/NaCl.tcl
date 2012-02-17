#
#    Simulation of a NaCl crystal or molten salt
#
#                                                       written by V. Ballenegger, dec. 2011
#unit of length: 1 angstrom
#unit of energy 1 kJ/mol
#unit of mass: 22.9898 u

#WARNING: I simulated a system where both Na and Cl ions have the SAME MASS!!!
#mass of Na: 22.9898 u
#mass of Cl: 35.453 u
#Note that this is not a problem, because equilibrium quantities do not depend on the masses.
#It would be better however to correct this script to use the correct masses.

#Note on converting units:
#1 u = atomic mass unit = 1.66053886 × 10-27 kg
#1 J = 0.001*N_A kJ/mol = 6.023E20 kJ/mol
#    ==> 1 kJ/mol = 1.66E-21 J

#The unit of time is
#   length unit*sqrt(mass unit/energy unit) = 10^-10*sqrt(22.9898*1.66E-27/1.66E-21) = 4.8x10^-13 s
# My time step (deltaT = 0.01) corresponds thus to 4.8 fs.
# The time between two "independent" snapshots is 1000*deltaT = 4.8 ps
# This is rather short, it might be better to use a longer time
# to be sure that the snapshots are independent.


set n_part 512

set epsilon 4.465
set sigma 2.32

setmd periodic 1 1 1
#==================== setup startup configuration ============================
set system 2

if {$system == 1} {
# =======================
# System 1: molten salt
# =======================
#Setup the system defined in Fig. 1 of
#J. Anwar, D. Frenkel, M. G. Noro, JCP 118 (2003): 728

#State point
set kT    [expr 2.0*$epsilon]
#Compute box length so that the reduced density (=rho*sigma^3) is equal to 0.383
set reduced_density 0.383
set box_l [expr pow($n_part/$reduced_density,1./3.)*$sigma]
setmd box_l $box_l $box_l $box_l

#
# Second test system
#
#Setup the system defined in Fig. 1 of
#  F. Lantelme, P. Turq, B. Quentrec and J.E.W. Lewis, Mol. Phys. 28 (1974):1537
#set n_part 432
#set kT 11.1456
#set box_l 24.64
#set bjerrum [expr 1389.3549/$kT]

#put particles in box
#type 1: Cl (charge -1)
#type 0: Na (charge +1)
set q 1;
set type 0
for {set i 0} { $i < $n_part } {incr i} {
  set posx [expr $box_l*[t_random]]
  set posy [expr $box_l*[t_random]]
  set posz [expr $box_l*[t_random]]
  
  set q    [expr -$q];
  set type [expr 1-$type]
  
  puts "$i $posx $posy $posz"
  part $i pos $posx $posy $posz q $q type $type
}
}

if {$system == 2} {
# =======================
# System 2: NaCl crystal
# =======================
puts "Setting up a NaCl crystal"
#true value at T=298K:
set unit_cell_length 5.682
   #Test: use the following value to get the same density as the melt
   #set unit_cell_length [expr 2*25.557108/8]
   #no, because the crystal melts...
set lattice_spacing [expr $unit_cell_length/2]
set nb_unit_cells 4
#nb of particles along a direction
set crystal_size [expr $nb_unit_cells*2]

set box_l [expr $nb_unit_cells*$unit_cell_length]
setmd box_l $box_l $box_l $box_l

   set part_no 0
   set q 1;
   set type 0
for {set i 0} { $i < $crystal_size } {incr i} {
for {set j 0} { $j < $crystal_size } {incr j} {
for {set k 0} { $k < $crystal_size } {incr k} {
  set posx [expr $i*$lattice_spacing]
  set posy [expr $j*$lattice_spacing]
  set posz [expr $k*$lattice_spacing]
  set q [expr pow(-1,$i+$j+$k)];
  if {$q > 0} {set type 0}
  if {$q < 0} {set type 1}
  #puts "$part_no $posx $posy $posz"
  part $part_no pos $posx $posy $posz q $q type $type
  incr part_no
}}}
#Crystal at room temperature:
#set T 298
set T 1000

set kT    [expr 8.31451*0.001*$T]
}

#Prefactor of coulombic forces [to have coulombic forces in kJ/(mol.A)]
set bjerrum [expr 1389.35485/$kT]


puts "----- system -------"
puts "Nb part: $n_part"
puts "box_l: $box_l"
puts "kT: $kT"
puts "--------------------"
  
#Langevin thermostat
setmd time_step 0.01; setmd skin 0.4

set gamma 1
thermostat langevin $kT $gamma

#Get # of degrees of freedom
   if { [regexp "ROTATION" [code_info]] } { 
     set deg_free 6
     puts "Rotations are turned on!"
   } else { set deg_free 3 }

#Generate a VTF file (for visualization in VMD)
set filePtr [open "conf_start.vtf" "w"]
writevsf $filePtr
writevcf $filePtr folded
close $filePtr

#==================== start of simulation ============================

#integration
set integ_steps 200

if {$system == 1} {
    #
	#Equilibrate system first using LJ forces
	#
	#Lennard-Jones interactions
	#1.12246 = 2^(1/6)
	set sig $sigma; set cut   [expr 1.12246*$sig]
	set eps $epsilon; set shift [expr 0.25*$eps]
	inter 0 0 lennard-jones $eps $sig $cut $shift 0
	inter 1 0 lennard-jones $eps $sig $cut $shift 0
	inter 1 1 lennard-jones $eps $sig $cut $shift 0
	
	#equilibration with capped forces
	puts "Equilibration: phase 1 (capped purely repulsive LJ forces)"
	for {set cap 20} {$cap < 200} {incr cap 20} {
	  set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
	  puts "t=[setmd time] E=[analyze energy total], T=$temp"
	  inter ljforcecap $cap;
	  integrate $integ_steps 
	}
	#remove LJ cap
	inter ljforcecap 0
	
	#Switch off LJ interactions
	inter 0 0 lennard-jones 0.0 $sig $cut 0.0 0
	inter 0 1 lennard-jones 0.0 $sig $cut 0.0 0
	inter 1 1 lennard-jones 0.0 $sig $cut 0.0 0
	
}

puts "Equilibration: phase 2 (Coulomb + Fumi-Tosi interactions)"

#Add Fumi-Tosi interactions
#
# warning: cutoff 5 is somewhat small. Should use 9 instead!
#
inter 0 0 bmhtf-nacl 25.4435  3.1546  101.1719  48.1771  2.34  9
inter 0 1 bmhtf-nacl 20.3548  3.1546  674.4793  837.0770 2.755  9
inter 1 1 bmhtf-nacl 15.2661  3.1546  6985.6786 14031.5785  3.170  9

#Add Coulomb interactions
inter coulomb $bjerrum p3m tunev2 accuracy 1e-3
#mesh 32

puts "Interactions:"
puts "[inter 0 0]"
for {set i 0} { $i < 10 } { incr i} {
  integrate [expr $integ_steps] 
  set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
  puts "t=[setmd time] E=[analyze energy total], T=$temp"
}

puts "Start of production run"
#INTEGRATION
setmd time 0.0
set integ_steps 1000

set nb_configs 10
for {set i 0} { $i < $nb_configs } { incr i} {
  set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
  puts "t=[setmd time] E=[analyze energy total], T=$temp"
  #integration
  integrate $integ_steps 

#Save config
#file format:
#   1: espresso  (needed for visualizing in VMD)
#   2: konfig    (VB format, for computing force accuracies)
set fileFormat 2

if {$fileFormat == 1} {
	set f [open "config_$i" "w"]
	blockfile $f write tclvariable {box_l}
	blockfile $f write variable box_l
	blockfile $f write particles {id pos type}
	close $f
}
if {$fileFormat == 2} {
	set f [open "konfig_$i" "w"]
	puts $f $n_part
	puts $f $box_l
	for {set j 0} { $j < $n_part } {incr j} {
	   puts -nonewline $f "[part $j print pos]"
	   puts $f " [part $j print q]"
	}
	close $f
}

#continue simulation loop
}

#Generate a VSF file (for visualization in VMD)
#this works by reading configs in fileFormat = 1, i.e. config_X files!!!
set filePtr [open "traj.vtf" "w"]
writevsf $filePtr
for {set i 0} { $i < $nb_configs } { incr i} {
   #read the config
   set f [open "config_$i" "r"]
   while { [blockfile $f read auto] != "eof" } {}
   close $f
   writevcf $filePtr folded
}
close $filePtr






