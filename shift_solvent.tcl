####################################################################################################################
# Uses VMD to shift coordinates and calculate PBC vectors
#
# This scripts is useful for fixing starting MD structure where the protein is positioned off-center in the z-axis
# This issue can make some visual analyses problematic
# The solvent including ions are shifted such that the protein lies in the middle. 
#
# -2020- ali.kusay@sydney.edu.au
####################################################################################################################

proc shift_solvent {args} {
    set default_protein_sel "protein"
    set default_solvent_sel "water or resname SOL WAT HOH TIP3 NA CL SOD CLA"
    if {[llength $args] == 1} {
        puts "Using protein selection $default_protein_sel and\nsolvent selection $default_solvent_sel"
    } elseif {[llength $args] == 2} {
        set default_protein_sel [lindex $args 1]
        puts "Using ADJUSTED protein selection $default_protein_sel and\nsolvent selection $default_solvent_sel"
    } elseif {[llength $args] == 3} {
        set default_protein_sel [lindex $args 1]
        set default_solvent_sel [lindex $args 2]
        puts "Using ADJUSTED protein selection $default_protein_sel and\nsolvent selection $default_solvent_sel"
    } else {
        error "1-3 arguments accepted: molid (required) protein_selection (optional) solvent_selection (optional)"
    }

    set molid [lindex $args 0]

    set protein_sel [atomselect $molid $default_protein_sel]
    set solvent_sel [atomselect $molid $default_solvent_sel]

    # get z center coordinate for protein and solvent
    set protein_z_mean [lindex [measure center $protein_sel] 2]
    set solvent_z_mean [lindex [measure center $solvent_sel] 2]

    # The 1.5 seems to work best for ensuring the same amount of solvent ends up on both sides of the protein
    set half_difference [expr ($protein_z_mean - $solvent_z_mean)/1.5]

    # get z min and max coordinate for solvent
    set sovent_min_max [measure minmax $solvent_sel]
    set solvent_z_min [lindex [lindex $sovent_min_max 0] 2]
    set solvent_z_max [lindex [lindex $sovent_min_max 1] 2]

    set shift_by [expr $solvent_z_max - $solvent_z_min]

    # Solvent below this will be shifted
    set threshold [expr $solvent_z_min + $half_difference]

    # same residue as ensures there are no residues straddling the cutoff
    set solvent_to_move [atomselect $molid "same residue as (($default_solvent_sel) and z < $threshold)"]

    set shift_coords [subst {0 0 $shift_by}]
    set no_to_shift [$solvent_to_move num]
    puts "Shifting $no_to_shift solvent atoms"
    $solvent_to_move moveby $shift_coords
}
