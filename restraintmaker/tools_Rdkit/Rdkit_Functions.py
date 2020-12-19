"""
.. automodule::  RDKIT_Functions
    :members:
    TODO:write DOCU!

"""

import glob
import os
import typing as t

import numpy as np

# import matplotlib.pyplot as plt    #TODO: remove

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Draw
    from rdkit.Chem import rdFMCS
except ImportError:
    raise ImportError("Could not import Rdkit package!")

# sys.path.append(os.path.dirname(__file__))

np.set_printoptions(suppress=True)


def parse_pdb_blocks_to_rdkit(pdb_mols: str) -> t.List[Chem.Mol]:
    rdk_mols = []
    for ind, lig in enumerate(pdb_mols):
        mol = AllChem.MolFromPDBBlock(lig, removeHs=False)
        setattr(mol, "name", 'mol_' + str(ind))
        setattr(mol, "resi", ind)
        rdk_mols.append(mol)

    return rdk_mols


def parse_pdb_to_rdkit(
        path: str = "/home/bschroed/code/restraintmaker/test/test_files/ligand_system/single_ligs_good/*pdb"):
    mols = []
    for ind, lig in enumerate(glob.glob(path)):
        mol = AllChem.MolFromPDBFile(lig, removeHs=False)
        setattr(mol, "name", os.path.splitext(os.path.basename(lig))[0])
        setattr(mol, "resi", ind)

        ring_obj = mol.GetRingInfo()
        atomRings = ring_obj.AtomRings()
        bonds = mol.GetBonds()
        atoms = mol.GetAtoms()

        print("Ammount of rings: " + str(ring_obj.NumRings()))
        print("Atoms in rings: " + str(atomRings))
        print("Atom num: " + str(len(atoms)))
        collect = ""
        for atom in atoms:
            collect += str(atom.GetHybridization()) + " "
        print("Hybridizations: " + str(collect))

        mols.append(mol)
    Draw.MolsToGridImage(mols, molsPerRow=len(mols))
    # h_mols = [AllChem.AddHs(mol) for mol in mols] #add H-atoms
    # sh_mols = [AllChem.EmbedMolecule(mol, AllChem.ETKDG()) for mol in h_mols ]    #generate 3d Coords
    # Chem.Kekulize(mol) #check correct mol structure
    return mols


# TODO: DOES NOT WORK. THE PROBLEM IS NOT IN THE MOLEUCLUE ORDER> IT ALSO DOES NOT WORK IF I PRESELECT ONE MOLEUCLE
def ring_atom_filter(selected: list, mols: t.List[Chem.Mol], selected_mols: dict):
    def _get_mol_offset(selected_mol: Chem.Mol, mols: t.List[Chem.Mol]) -> int:
        # deal with offsets of ligands:
        mol_offset = 1
        for mol in sorted(mols, key=lambda x: x.resi):
            if (selected_mol.name == mol.name):
                print("FOUND_MOLNAME: ", mol.name)
                break
            else:
                mol_offset += len(mol.GetAtoms())
        else:
            raise ValueError('The molecule is not a member of the list of moleucles.')

        return mol_offset

    already_selected_atomIdx = []
    collect = []
    for ind, mol in enumerate(mols):
        # deal with offsets:
        mol_offset = _get_mol_offset(selected_mol=mol, mols=mols)

        # print("Mol " + str(mol_ind) +" "+str(mol.name)+" has offset: " + str(mol_offset))
        all_atoms = mol.GetAtoms()
        # print("All atoms:"+str([atom.GetIdx() for atom in all_atoms]))
        # get all ring atoms
        all_atoms_in_rings = list(filter(lambda x: x.IsInRing(), all_atoms))
        print("all_atoms in rings: " + str([atom.GetIdx() for atom in all_atoms_in_rings]))
        # print out all ring atoms
        test = ""
        for atom in all_atoms_in_rings:
            test += str(mol_offset + atom.GetIdx()) + " "
        print("all ring Atoms +offset: " + str(test))
        # reduce ring atoms by already selected atoms
        if len(already_selected_atomIdx) > 0:
            atoms_selected_in_rings = list(
                filter(lambda x: (x.GetIdx() + int(mol_offset)) in already_selected_atomIdx, all_atoms_in_rings))
        else:
            atoms_selected_in_rings = all_atoms_in_rings
        # collect filtered atoms
        # print("only preselected atoms in rings+off"+str([mol_offset + atom.GetIdx() for atom in atoms_selected_in_rings]))
        collect += [mol_offset + atom.GetIdx() for atom in atoms_selected_in_rings]

    print("Finally Selected atoms: ", collect)
    return collect


# get molecule selection
def mcs_selection(mols: t.List[Chem.Mol], min_MCS_size: int = 6, verbose: bool = True):
    '''
    :warning: results do not seem reasonable YET. TODO: Check what exactly is not going as expected
    :param mols:
    :type mols:
    :param min_MCS_size:
    :type min_MCS_size:
    :param verbose:
    :type verbose:
    :return:
    :rtype:
    '''

    def _get_mol_offset(selected_mol: Chem.Mol, mols: t.List[Chem.Mol]) -> int:
        # deal with offsets of ligands:
        mol_offset = 1
        for mol in sorted(mols, key=lambda x: x.resi):
            if (selected_mol.name == mol.name):
                # print("FOUND_MOLNAME: ", mol.name)
                break
            else:
                mol_offset += len(mol.GetAtoms())
        return mol_offset

    # generate MCS out of tools_Rdkit mols:
    ## most common substructres
    if verbose: print("Calculate MCS")
    mcs = [[] for mol in mols]
    for ind1, mol1 in enumerate(mols):
        if verbose: print(ind1)
        for ind2, mol2 in enumerate(mols):
            if verbose: print("\t", ind2)
            if (mol1 == mol2):
                continue
            else:
                tmp_mcs = rdFMCS.FindMCS([mol1, mol2], ringMatchesRingOnly=True)
                mcs[ind1].append(tmp_mcs)
                mcs[ind2].append(tmp_mcs)
    if verbose: print("MCS: got ", len(mcs), "\t aim ", len(mols) * len(mols))
    ## convert result to smart strings and make these to new tools_Rdkit mols:
    print(mcs)

    print(mcs[0])
    smart_mols = [Chem.MolFromSmarts(mcs_s.smartsString) for molecule_smarts in mcs for mcs_s in molecule_smarts]
    if verbose: print("SMARTS: got ", len(smart_mols), "\t aim ", len(mols) * len(mols))

    ## vis
    # names = [mol.name for mol2 in mols for mol in mols if(mol!=mol2)]
    # Draw.MolsToGridImage(smart_mols, molsPerRow=len(mols), legends=names)
    # Draw.MolsToGridImage(smart_mols, molsPerRow=len(mols), legends=names)

    # get atom nums of mcs_smart match for each row and col - Carefull there could be multiple matches!
    matching_atoms = []
    mol_matrix = [mol2 for mol2 in mols for mol1 in mols if (mol1 != mol2)]  # flattened mol matrix
    for mol, smart in zip(mol_matrix, smart_mols):
        substructre = mol.GetSubstructMatches(smart)
        print(substructre)
        matching_atoms.append(substructre)
    if verbose: print("Matching Atoms: got ", len(matching_atoms), "\t aim ", len(mols) * len(mols))

    if verbose: print("MATRIX LENGHT: ", len(matching_atoms))
    ##add offset to atom indices
    offset_added_mol_matches = []

    for ind, mol_matches in enumerate(matching_atoms):
        offset_added_matches = []
        if verbose: print("MOL: ", ind // (len(mols) - 1))
        for match in mol_matches:
            if (len(match) >= min_MCS_size):
                if verbose: print("\t", match)
                offset = _get_mol_offset(mols[ind // (len(mols) - 1)], mols)
                offset_added_matches.extend(tuple(map(lambda x: x + offset, match)))
            else:
                continue
        offset_added_mol_matches.extend(offset_added_matches)
    if verbose: print("MCS: got ", len(offset_added_matches), "\t aim ", len(mols) * len(mols))

    if verbose:  print('len:', str(len(offset_added_mol_matches)), ' -', offset_added_mol_matches)
    return offset_added_mol_matches

    ## reduce to a molecule pairwise selection:
    mol_pairwise_selection = []
    ###SOMETHING

    return mol_pairwise_selection


# PCA
# FUNCS
# Calc PCA
def _calc_pca(coords, dims, verbose: bool = True):
    if (verbose): print("PCA: ")
    # center data: (scale)
    geometric_mean = np.mean(coords, axis=0)
    centered_coords = coords - geometric_mean
    # get covariance_matrix
    covaricance_matrix = np.cov(centered_coords.T)
    if (verbose): print("\tcovariance: \n\t", "\n\t ".join(map(str, list(covaricance_matrix))))
    # eigenvalue decomp. => principle components

    eigen_values, eigen_vectors = np.linalg.eig(covaricance_matrix)
    if (verbose): print("\tPC: \n\t", "\n\t ".join(map(str, list(eigen_vectors))))
    if (verbose): print("\teigenvalues: ", eigen_values)

    # SORT for size
    sorted_eigen = list(reversed(sorted(zip(eigen_values, eigen_vectors), key=lambda x: x[0])))
    reduced_eigenvalues = [eigen[0] for eigen in sorted_eigen][:dims]
    if (verbose): print(sorted_eigen)

    # transform coordinates onto new coords.
    if (verbose): print("\tPC: \n\t", "\n\t ".join(map(str, list(eigen_vectors))))

    transformed_coords = eigen_vectors.T.dot(centered_coords.T)
    reduced_transformed_coords = np.array([transf_dim for transf_dim, eigen_val in zip(transformed_coords, eigen_values)
                                           if (eigen_val in reduced_eigenvalues)])
    if (verbose): print(reduced_transformed_coords)
    return reduced_transformed_coords, eigen_vectors, eigen_values, covaricance_matrix


def _calc_pca_without_scaling(coords, dims, verbose: bool = True):
    '''Analogous to _calc_pca, but it does not do any coordinate transfromations. (Rescaling etc. This should be done before calling this function.
    This will help save time, if the function is called many times'''
    # if (verbose): print("PCA: ")

    # get covariance_matrix
    # coords = list(list(c) for c in coords)

    coords = np.array(coords)
    covaricance_matrix = np.cov(coords.T)
    # if (verbose): print("\tcovariance: \n\t", "\n\t ".join(map(str, list(covaricance_matrix))))
    # eigenvalue decomp. => principle components

    eigen_values, eigen_vectors = np.linalg.eig(covaricance_matrix)
    # if (verbose): print("\tPC: \n\t", "\n\t ".join(map(str, list(eigen_vectors))))
    # if (verbose): print("\teigenvalues: ", eigen_values)

    # SORT for size
    sorted_eigen = list(reversed(sorted(zip(eigen_values, eigen_vectors), key=lambda x: x[0])))
    reduced_eigenvalues = [eigen[0] for eigen in sorted_eigen][:dims]
    # print(sorted_eigen)

    return reduced_eigenvalues


# Calculate polygon area!    #from plotly -
def PolygonSort(corners, title=""):
    n = len(corners)
    cx = float(sum(x for x, y in corners)) / n
    cy = float(sum(y for x, y in corners)) / n
    cornersWithAngles = []
    for x, y in corners:
        an = (np.arctan2(y - cy, x - cx) + 2.0 * np.pi) % (2.0 * np.pi)
        cornersWithAngles.append((x, y, an))
    cornersWithAngles.sort(key=lambda tup: tup[2])
    last = []
    for x, y, an in cornersWithAngles:
        last.append((x, y))
    return last


def PolyArea(corners, title=""):  # Shoelace formula#from plotly -
    print("small", corners)
    corners = PolygonSort(corners[:, :2])
    print("Corners", corners)
    n = len(corners)  # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    print("area", area)
    area = abs(area) / 2.0
    return area


# PLOT DATA
def _plot_data(coords, eigen_vectors, title):
    ##original data
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(coords[:, 0], coords[:, 1], c="k", label="data")
    ax.plot((-eigen_vectors[0, 0], eigen_vectors[0, 0]), (-eigen_vectors[0, 1], eigen_vectors[0, 1]), label="pc1")
    ax.plot((- eigen_vectors[1, 0], eigen_vectors[1, 0]), (-eigen_vectors[1, 1], eigen_vectors[1, 1]), label="pc2")
    ax.legend()
    ax.set_title("orig data - " + title)

    fig.set_tight_layout("thight")
    fig.show()


def _plot_pc_projection(projected_coords, title):
    ##porjection of data
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(projected_coords[0], projected_coords[1], c="k", label="data")
    ax.scatter(projected_coords[0], [0 for x in projected_coords[0]], c="r", label="prPC1")
    ax.scatter([0 for x in projected_coords[1]], projected_coords[1], c="b", label="prPC2")
    ax.legend()
    ax.set_title("pca_projection - " + title)

    fig.set_tight_layout("thight")
    fig.show()


def _plot_2D_PCA(coords, selected_arr, mol_name="", eig_vecs=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(coords[:, 0], coords[:, 1], c="k", alpha=0.4)
    ax.scatter(selected_arr[:, 0], selected_arr[:, 1], c="r", label="selected")

    if (not eig_vecs):
        for ind, vec in enumerate(eig_vecs):
            ax.plot((vec[0], vec[1]), label="pca" + str(ind))

    ax.set_xlabel('X ')
    ax.set_ylabel('Y ')
    ax.set_title("2DPCA " + mol_name + " molecule")
    ax.legend()

    fig.set_tight_layout("thight")
    fig.show()


def _plot_3D_PCA(coords, eig_vecs=None, title=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ##vis mol
    ax.scatter(list(coords[:, 0]), list(coords[:, 1]),
               list(coords[:, 2]), c="k", alpha=0.4)
    ##vis selected
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c="r", label="selected")

    ##vis pca vector
    if (type(eig_vecs) != None and len(eig_vecs) > 0):
        for ind, vec in enumerate(eig_vecs):
            ax.plot(xs=[0, 1 * vec[0]], ys=[0, 1 * vec[1]], zs=[0, 1 * vec[2]], label="pca" + str(ind))

    ##vis tuning
    ax.set_xlabel('X ')
    ax.set_ylabel('Y ')
    ax.set_zlabel('Z ')
    ax.set_title("3DPCA " + title + " molecule")
    # ax.legend()

    fig.set_tight_layout("thight")
    fig.show()
