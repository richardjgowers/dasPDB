from rdkit import Chem


def split_into_residues(m) -> list[Chem.Mol]:
    """Split along any bond between different residues
    
    Residues are defined by a tuple of (chain, resnum, icode, resname)
    """
    def mi_hash(atom):
        monomerinfo = atom.GetMonomerInfo()

        return (
            monomerinfo.GetChainId(),
            monomerinfo.GetResidueNumber(),
            monomerinfo.GetInsertionCode(),
            monomerinfo.GetResidueName(),
        )

    bonds_to_break = []
    for b in m.GetBonds():
        if mi_hash(b.GetBeginAtom()) != mi_hash(b.GetEndAtom()):
            bonds_to_break.append(b.GetIdx())


    broken_m = Chem.rdmolops.FragmentOnBonds(
        m, bonds_to_break, addDummies=True
    )

    residues = Chem.rdmolops.GetMolFrags(broken_m, asMols=True, sanitizeFrags=False)
    
    return residues
