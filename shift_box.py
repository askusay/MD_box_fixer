'''
Uses Mdtraj to shift coordinates and calculate PBC vectors

This scripts is useful to fix starting MD structure where the protein is positioned off-center in the z-axis, making some visual analyses problematic, the solvent including ions is shifted such that the protein lies in the middle. 
In addition, the whole structure can be shifted such that the solvent box vertex aligns with a desired coordinate i.e. 0,0,0 which is necessary to avoid some visual PBC artefacts using some MD softwares
Lastly, the script enables calculation of PBC vectors which are appended to the pdb file

-2020- ali.kusay@sydney.edu.au
'''

import argparse
import sys
import logging

import numpy as np

import mdtraj as md

logging.basicConfig(stream=sys.stdout,
                    level=logging.DEBUG,
                    format='[%(asctime)s %(levelname)s] %(message)s',
                    datefmt='%Y/%m/%d %H:%M:%S')

parser = argparse.ArgumentParser()
parser.add_argument('structure',help='Input pdb file or equivalent')
parser.add_argument('output', help='Output file name')
parser.add_argument('--shift_solvent', action='store_true', dest='solv', help='Sift solvent atoms such that protein lies in the centre', default=False)
parser.add_argument('--shift_box', action='store', dest='box', help='Align box vertex (solvent) with desired position - default (0,0,0), entered as 0:0:0')
parser.add_argument('--find_pbc', action='store_true', dest='pbc', help='Calculate unitcell vectors and include in output structure', default=True)

args = parser.parse_args()

logging.info('Loading structure and finding solvent indicies...')
pdb = md.load(args.structure)

salt_list = {'CLA','SOD','POT','CAL','NA','CL','K','MG','CA'}
water_list = {'HOH', 'SPC','TIP3','T3P','TIP'}
solvent_atoms = np.array([atom.index for atom in pdb.topology.atoms if atom.residue.name in (salt_list | water_list)])

if args.solv:
    # findind indicies
    logging.info('---Shifting solvent to center protein---')
    protein_atoms = pdb.topology.select('protein')

    exclude_list = [atom.residue.name for atom in pdb.topology.atoms if atom.index not in solvent_atoms]
    exclude_names = list(dict.fromkeys(exclude_list))

    # finding atoms below cut-off
    positions = pdb.xyz[0]
    protein_coords = positions[protein_atoms]
    solvent_coords = positions[solvent_atoms]
    half_difference = (protein_coords.mean(axis=0)[-1] - solvent_coords.mean(axis=0)[-1])/2

    solvent_min, solvent_max = solvent_coords.min(axis=0)[-1], solvent_coords.max(axis=0)[-1]
    shift_by = solvent_max - solvent_min

    threshold = solvent_min + half_difference
    alter_pos = solvent_coords[:,2] < (threshold)
    rough_sele = solvent_atoms[alter_pos]
    logging.info(f'{len(rough_sele)} solvent atoms found below cut-off: {round(threshold, 3)} nm, checking for straddling residues...')

    # accounting for residues that straddle across the cut-off
    final_sele = list()

    for residue in pdb.topology.residues:
        if residue.name in exclude_names:
            continue
        indicies = {atom.index for atom in residue.atoms}
        if len(indicies & set(rough_sele)) > 0:
            final_sele.extend(indicies)

    logging.info(f'Found {len(final_sele)-len(rough_sele)} additional solvent atoms to be shifted')

    positions[final_sele] += [0,0,shift_by]

    logging.info(f'Solvent atoms shifted')

if args.box:
    logging.info(f'---Aligning box---')
    try:
        xyz = [float(i) for i in args.box.split(':')]
    except:
        logging.debug(f'{args.box} is not entered corrently, must be 3 digits seperated by ":" i.e. "0:0:0" - exiting.')
        sys.exit(1)

    positions = pdb.xyz[0]
    min_pos = positions[solvent_atoms].min(axis=0)

    positions += (xyz - min_pos)
    logging.info(f'All atoms shifted by {xyz - min_pos}')

if args.pbc:
    logging.info(f'---Finding PBC vectors and adding to stucture---')

    if pdb.unitcell_vectors is not None:
        logging.warning(f'Structure already contains PBC vectors')
    
    pbc_arr = np.zeros([1,3,3],dtype='float32')
    pbc_bool = np.array([[[True,False,False],[False,True,False],[False,False,True]]])

    positions = pdb.xyz[0]
    pbc_vectors = positions.max(axis=0) - positions.min(axis=0)
    pbc_arr[pbc_bool] = pbc_vectors

    pdb.unitcell_vectors = pbc_arr
    logging.info(f'{pbc_vectors} nm vectors added to the structure')

logging.info(f'Saving structure as {args.output}')

pdb.save(args.output)