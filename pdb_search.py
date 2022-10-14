import gemmi
from rdkit import Chem


def block_to_mol(b: gemmi.cif.Block, exclude_leaving_atoms=False) -> Chem.Mol:
    """Convert gemmi Block to RDKit Mol"""
    orders = {
        'SING': Chem.rdchem.BondType.SINGLE,
        'DOUB': Chem.rdchem.BondType.DOUBLE,
        'TRIP': Chem.rdchem.BondType.TRIPLE,
    }

    atom_chiralities = {
        'R': Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
        'S': Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
        'N': Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
    }

    mol = Chem.EditableMol(Chem.Mol())

    # populate atoms
    atoms = b.get_mmcif_category('_chem_comp_atom.')
    if not atoms:  # e.g. UNL
        return Chem.Mol()
    
    leaving_atom_indices = set()
    for i, (elem, c, a, s, l) in enumerate(zip(
                                atoms['type_symbol'],
                                atoms['charge'],
                                atoms['pdbx_aromatic_flag'],
                                atoms['pdbx_stereo_config'],
                                atoms['pdbx_leaving_atom_flag'])):
        isotope = 0
        if exclude_leaving_atoms and l == 'Y':
            elem = '*'
            a = 'N'
            c = 0
            s = 'N'
            leaving_atom_indices.add(i)
        else:
            if elem == 'X': # e.g. ASX
                elem = '*'
            elif elem == 'D':  # e.g. D3O
                elem = 'H'
                isotope = 3
            if c is None:  # e.g. 3CD
                c = 0
        atom = Chem.Atom(elem.capitalize())
        atom.SetIsotope(isotope)
        atom.SetIsAromatic(a == 'Y')
        atom.SetFormalCharge(int(c))
        # atom.SetNoImplicit(True)
        atom.SetChiralTag(atom_chiralities[s])
        mol.AddAtom(atom)

    # bonds with orders
    bonds = b.get_mmcif_category('_chem_comp_bond.')    
    if bonds:
        id_to_index = {nm: i for i, nm in enumerate(b.get_mmcif_category('_chem_comp_atom.')['atom_id'])}

        for x, y, o in zip(bonds['atom_id_1'], bonds['atom_id_2'], bonds['value_order']):
            i, j = id_to_index[x], id_to_index[y]
            if not (i in leaving_atom_indices and j in leaving_atom_indices):
                mol.AddBond(i, j, order=orders[o])
            
    if exclude_leaving_atoms:
        to_remove = []
        # remove any leaving atoms that now have no bonds
        # e.g. ARG removes -OH leaving -*.* with an orphan *
        m = mol.GetMol()
        for i in leaving_atom_indices:
            if not m.GetAtomWithIdx(i).GetBonds():
                to_remove.append(i)
        # remove in descending order to avoid grief
        for i in sorted(to_remove, reverse=True):
            mol.RemoveAtom(i)

    return mol.GetMol()