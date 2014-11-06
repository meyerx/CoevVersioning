#!/usr/bin/wish
#
# coev.tcl: Simple GUI for the coev program.
#
# Wim Hordijk   Last modified: 26 April 2013
#

#
# Application title
#
wm title . "Coev"

#
# Main frames
#
frame .frames -bd 3 -relief groove
frame .action
pack  .frames
pack  .action -side bottom

#
# Method
#
frame       .frames.mlabel -bd 3 -relief groove
pack        .frames.mlabel -fill x
label       .frames.mlabel.label -text "Method" -bg gray50 -anchor c
pack        .frames.mlabel.label -fill x

frame       .frames.method
pack        .frames.method -anchor w
label       .frames.method.label -text "Method:" -width 12 -anchor w
radiobutton .frames.method.rb1 -text "Bayesian    " \
            -variable method -value "bayes" \
            -command {.frames.sfreq.label configure -state normal
            .frames.sfreq.box configure -state normal
            .frames.pfreq.label configure -state normal
            .frames.pfreq.box configure -state normal
            .frames.it.label configure -state normal
            .frames.it.box configure -state normal
            .frames.burnin.label configure -state normal
            .frames.burnin.box configure -state normal}
radiobutton .frames.method.rb2 -text "ML" -variable method \
            -value "ml" -command {.frames.sfreq.label configure -state disabled
            .frames.sfreq.box configure -state disabled
            .frames.pfreq.label configure -state disabled
            .frames.pfreq.box configure -state disabled
            .frames.it.label configure -state disabled
            .frames.it.box configure -state disabled
            .frames.burnin.label configure -state disabled
            .frames.burnin.box configure -state disabled}
pack        .frames.method.label .frames.method.rb1 .frames.method.rb2 \
            -side left

frame       .frames.data
pack        .frames.data -anchor w

#
# Input & output
#
frame       .frames.input -bd 3 -relief groove
pack        .frames.input -fill x
label       .frames.input.label -text "Input/output" -bg gray50 -anchor c
pack        .frames.input.label -fill x

frame       .frames.tree
pack        .frames.tree -anchor w
label       .frames.tree.label -text "Tree file:" -width 12 -anchor w
entry       .frames.tree.entry -textvariable treeFile -width 25
button      .frames.tree.button -text "Browse..." \
            -command {set treeFile [tk_getOpenFile \
	    -title "Select a tree file..."]
	    .frames.tree.entry xview end}
pack        .frames.tree.label .frames.tree.entry .frames.tree.button -side left

frame       .frames.align
pack        .frames.align -anchor w
label       .frames.align.label -text "Alignment file:" -width 12 -anchor w
entry       .frames.align.entry -textvariable alignFile -width 25
button      .frames.align.button -text "Browse..." \
            -command {set alignFile [tk_getOpenFile \
	    -title "Select an alignment file..."]
	    .frames.align.entry xview end}
pack        .frames.align.label .frames.align.entry .frames.align.button \
            -side left

frame       .frames.cols
pack        .frames.cols -anchor w
label       .frames.cols.label1 -text "Position 1:" -width 12 -anchor w
spinbox     .frames.cols.box1 -textvariable col1 -width 6 -from 1 \
            -to 1000000 -increment 1
label       .frames.cols.label2 -text "    Position 2:" -width 12 -anchor w
spinbox     .frames.cols.box2 -textvariable col2 -width 6 -from 1 \
            -to 1000000 -increment 1
pack        .frames.cols.label1 .frames.cols.box1 .frames.cols.label2 \
            .frames.cols.box2 -side left

frame       .frames.log
pack        .frames.log -anchor w
label       .frames.log.label -text "Log file:" -width 12 -anchor w
entry       .frames.log.entry -textvariable logFile -width 25
button      .frames.log.button -text "Browse..." \
            -command {set logFile [tk_getOpenFile \
	    -title "Select a log file..."]
	    .frames.log.entry xview end}
pack        .frames.log.label .frames.log.entry .frames.log.button -side left

#
# Baysian parameters
#
frame       .frames.print -bd 3 -relief groove
pack        .frames.print -fill x
label       .frames.print.label -text "Bayesian parameters" -bg gray50 -anchor c
pack        .frames.print.label -fill x

frame       .frames.sfreq
pack        .frames.sfreq -anchor w
label       .frames.sfreq.label -text "Log frequency:" -width 12 -anchor w
spinbox     .frames.sfreq.box -textvariable sfreq -width 8 -from 0 \
            -to 1000000 -increment 1000
pack        .frames.sfreq.label .frames.sfreq.box -side left

frame       .frames.pfreq
pack        .frames.pfreq -anchor w
label       .frames.pfreq.label -text "Print frequency:" -width 12 -anchor w
spinbox     .frames.pfreq.box -textvariable pfreq -width 8 -from 0 \
            -to 1000000 -increment 1000
pack        .frames.pfreq.label .frames.pfreq.box -side left

frame       .frames.it
pack        .frames.it -anchor w
label       .frames.it.label -text "Nr. iterations:" -width 12 -anchor w
spinbox     .frames.it.box -textvariable IT -width 8 -from 0 \
            -to 1000000000 -increment 1000
pack        .frames.it.label .frames.it.box -side left

frame       .frames.burnin
pack        .frames.burnin -anchor w
label       .frames.burnin.label -text "Burn in:" -width 12 -anchor w
spinbox     .frames.burnin.box -textvariable burnin -width 8 -from 0 \
            -to 1000000000 -increment 1000
pack        .frames.burnin.label .frames.burnin.box -side left

#
# General parameters
#
frame       .frames.params -bd 3 -relief groove
pack        .frames.params -fill x
label       .frames.params.label -text "Instantaneous matrix paramters" -bg gray50 -anchor c
pack        .frames.params.label -fill x

frame       .frames.s
pack        .frames.s -anchor w
label       .frames.s.label -text "s: " -width 12 -anchor w
spinbox     .frames.s.box -width 8 -from 0.0 -to 1.0 -increment 0.001 \
            -textvariable s
pack        .frames.s.label .frames.s.box -side left

frame       .frames.d
pack        .frames.d -anchor w
label       .frames.d.label -text "d: " -width 12 -anchor w
spinbox     .frames.d.box -width 8 -from 0.0 -to 100.0 -increment 1.0 \
            -textvariable d
pack        .frames.d.label .frames.d.box -side left

frame       .frames.r1
pack        .frames.r1 -anchor w
label       .frames.r1.label -text "r1: " -width 12 -anchor w
spinbox     .frames.r1.box -width 8 -from 0.0 -to 10.0 -increment 0.1 \
            -textvariable r1
pack        .frames.r1.label .frames.r1.box -side left

frame       .frames.r2
pack        .frames.r2 -anchor w
label       .frames.r2.label -text "r2: " -width 12 -anchor w
spinbox     .frames.r2.box -width 8 -from 0.0 -to 10.0 -increment 0.1 \
            -textvariable r2
pack        .frames.r2.label .frames.r2.box -side left

#
# Action buttons
#
button .action.reset -text "Reset" -command {resetValues}
button .action.run   -text "Run"   -command {runCoev}
button .action.quit  -text "Quit"  -command {exit}
pack   .action.reset .action.run .action.quit -side left -padx 5

#
# Output window.
#
toplevel  .output
wm title  .output "Coev output"

text      .output.text -yscrollcommand {.output.scroll set}
scrollbar .output.scroll -command {.output.text yview}
pack      .output.text .output.scroll -side left -fill y

#
# resetValues: Reset all variables to their default values.
#
proc resetValues {} {
  global treeFile alignFile col1 col2 logFile IT sfreq pfreq burnin s d r1 r2 \
         method data

  set method "bayes"
  set data "nt"
  set treeFile "treeInput.txt"
  set alignFile "alignment.txt"
  set col1 1
  set col2 2
  set logFile "output.log"
  set IT 10000
  set sfreq 1000
  set pfreq 1000
  set burnin 0
  set s 0.001
  set d 100
  set r1 1
  set r2 1

  .frames.sfreq.label configure -state normal
  .frames.sfreq.box configure -state normal
  .frames.pfreq.label configure -state normal
  .frames.pfreq.box configure -state normal
  .frames.it.label configure -state normal
  .frames.it.box configure -state normal
  .frames.burnin.label configure -state normal
  .frames.burnin.box configure -state normal
}

#
# runCoev: Run the coev program with the given options.
#
proc runCoev {} {
  global treeFile alignFile col1 col2 logFile IT sfreq pfreq burnin s d r1 r2 \
         execDir method data

  #
  # Run the requested program and display the output.
  #
  .output.text delete 1.0 end
  set execFile [file join $execDir "coev"]
  set fp [open "|$execFile -method $method -data $data -IT $IT -sfreq $sfreq \
                 -pfreq $pfreq -burnin $burnin -s $s -d $d -r1 $r1 -r2 $r2 \
                 -tree $treeFile -align $alignFile -out $logFile -cols $col1 \
                 $col2" r]
  while {[gets $fp line] >= 0} {
    .output.text insert end $line
    .output.text insert end "\n"
    .output.text see end
    update idletasks
  }
  if {[catch {close $fp} err]} {
    tk_messageBox -icon error -title "Error" -type ok -message "The coev\
    program caused an error: $err"
  } else {
    tk_messageBox -icon info -title "Info" -type ok -message "The coev program\
    finished successfully."
  }
}

#
# Set default values and go...
#
set execDir [file dirname $argv0]
resetValues

#
# EOF: coev.tcl
#
