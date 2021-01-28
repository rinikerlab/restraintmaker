"""
    This module is linking the restraintmaker functionality to rdkit
     - MCS search
     - molecule ring filter
     - etc.
"""

import numpy as np
import typing as t

from restraintmaker.utils.Utilities import print

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Draw
    from rdkit.Chem import rdFMCS
except ImportError:
    raise ImportError("Could not import Rdkit package!")

np.set_printoptions(suppress=True)


def parse_pdb_blocks_to_rdkit(pdb_mols: str) -> t.List[Chem.Mol]:
    """
        Convert a pdblock to and rdkit-Mol

    Parameters
    ----------
    pdb_mols: str

    Returns
    -------
    rdkit.Mol
        an rdkit mol obj
    """
    rdk_mols = []
    for ind, lig in enumerate(pdb_mols):
        print(lig, mv=5)
        mol = AllChem.MolFromPDBBlock(lig, removeHs=False)
        if(mol is None):
            raise ValueError("Rdkit Mol was None - "+str(mol)+"!\n got "+str(mol))

        setattr(mol, "name", 'mol_' + str(ind))
        setattr(mol, "resi", ind)
        rdk_mols.append(mol)

    return rdk_mols


def ring_atom_filter(selected: list, mols: t.List[Chem.Mol], selected_mols: dict) -> t.List[int]:
    """
        Filter molecules for rings and remove non-ring atoms from the selection

    Parameters
    ----------
    mols: t.List[Chem.Mol]
        molecules to select from

    Returns
    -------
    t.List[int]
        list of stil selected atom indices, inside a ring
    """

    def _get_mol_offset(selected_mol: Chem.Mol, mols: t.List[Chem.Mol]) -> int:
        """
            determines the molecule offset
        """
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
def mcs_selection(mols: t.List[Chem.Mol], min_MCS_size: int = 6) -> t.List:
    """
        The selection is getting only the MCS of pairwise molecules

        Warnings:  results do not seem reasonable YET.
                TODO: Check what exactly is not going as expected

    Parameters
    ----------
    mols: t.List[Chem.Mol]
        molecules
    min_MCS_size :
        what is the minimal MCS size
    verbose : bool
        loud and noisy

    Returns
    -------
    List

    """

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
    print("Calculate MCS", mv=0)
    mcs = [[] for mol in mols]
    for ind1, mol1 in enumerate(mols):
        print(ind1)
        for ind2, mol2 in enumerate(mols):
            print("\t", ind2)
            if (mol1 == mol2):
                continue
            else:
                tmp_mcs = rdFMCS.FindMCS([mol1, mol2], ringMatchesRingOnly=True)
                mcs[ind1].append(tmp_mcs)
                mcs[ind2].append(tmp_mcs)
    print("MCS: got ", len(mcs), "\t aim ", len(mols) * len(mols), mv=1)

    ## convert result to smart strings and make these to new tools_Rdkit mols:
    smart_mols = [Chem.MolFromSmarts(mcs_s.smartsString) for molecule_smarts in mcs for mcs_s in molecule_smarts]
    print("SMARTS: got ", len(smart_mols), "\t aim ", len(mols) * len(mols), mv=1)

    # get atom nums of mcs_smart match for each row and col - Carefull there could be multiple matches!
    matching_atoms = []
    mol_matrix = [mol2 for mol2 in mols for mol1 in mols if (mol1 != mol2)]  # flattened mol matrix
    for mol, smart in zip(mol_matrix, smart_mols):
        substructre = mol.GetSubstructMatches(smart)
        print(substructre)
        matching_atoms.append(substructre)
    print("Matching Atoms: got ", len(matching_atoms), "\t aim ", len(mols) * len(mols), mv=1)
    print("MATRIX LENGHT: ", len(matching_atoms), mv=1)

    ##add offset to atom indices
    offset_added_mol_matches = []

    for ind, mol_matches in enumerate(matching_atoms):
        offset_added_matches = []
        print("MOL: ", ind // (len(mols) - 1), mv=1)
        for match in mol_matches:
            if (len(match) >= min_MCS_size):
                print("\t", match)
                offset = _get_mol_offset(mols[ind // (len(mols) - 1)], mols)
                offset_added_matches.extend(tuple(map(lambda x: x + offset, match)))
            else:
                continue
        offset_added_mol_matches.extend(offset_added_matches)
    print("MCS: got ", len(offset_added_matches), "\t aim ", len(mols) * len(mols), mv=1)
    print('len:', str(len(offset_added_mol_matches)), ' -', offset_added_mol_matches, mv=1)
    return offset_added_mol_matches


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
