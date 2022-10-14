from rdkit import Chem
from rdkit.Chem import AllChem
from MDAnalysis.converters.RDKit import (
    _infer_bo_and_charges, _standardize_patterns
)
from pdb_split import split_into_residues
from utils import get_ccd_slib


ccd_slib = get_ccd_slib()
# TODO: fix matching without bond orders and charges
sparam = Chem.SubstructMatchParameters()
sparam.aromaticMatchesConjugated = True


def infer_mol(base_mol) -> Chem.Mol:
    residues = split_into_residues(base_mol)
    mols = []
    for res in residues:
        mol = infer_from_ccd(res)
        # TODO: uncomment try and catch the right exception(s)
        # try:
        #     mol = infer_from_ccd(res)
        # except:
        #     mol = unknown_res_to_mol(res)
        mols.append(mol)
    mol = reassemble_fragments(mols)
    return mol


def unknown_res_to_mol(pdb_block: str) -> Chem.Mol:
    mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    _infer_bo_and_charges(mol)
    mol = _standardize_patterns(mol)
    return mol


def infer_from_ccd(mol):
    matches = ccd_slib.GetMatches(mol)
    template = next(iter(matches))
    return AllChem.AssignBondOrdersFromTemplate(template, mol)


def reassemble_fragments(mols):
    # TODO
    return mols[0]
