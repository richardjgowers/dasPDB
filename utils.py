import gzip
from pathlib import Path
from tempfile import NamedTemporaryFile

import gemmi
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import rdSubstructLibrary

from pdb_search import block_to_mol


HERE = Path(__file__).parent

def gen_ccd_slib():
    print("Downloading components library")
    cif_url = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"
    with requests.get(cif_url, stream=True) as resp:
        resp.raise_for_status()
        with NamedTemporaryFile("wb") as temp:
            for chunk in resp.iter_content(chunk_size=8192):
                temp.write(chunk)
            temp.seek(0)
            with gzip.open(temp.name, "rb") as fh:
                doc = gemmi.cif.read_string(fh.read().decode())
    slib = rdSubstructLibrary.SubstructLibrary()
    print("Building fragments library")
    RDLogger.DisableLog("rdApp.*")
    for b in doc:
        mol = block_to_mol(b, exclude_leaving_atoms=True)
        mol.UpdatePropertyCache(strict=False)
        try:
            m = Chem.RemoveHs(mol)
        except (Chem.AtomValenceException, Chem.AtomKekulizeException):
            m = Chem.RemoveHs(mol, sanitize=False)
        slib.AddMol(m)
    RDLogger.EnableLog("rdApp.*")
    print("Saving library for future use")
    with open(HERE / "ccd_slib.bin", "wb") as fh:
        fh.write(slib.Serialize())
    return slib


def get_ccd_slib():
    ccd_slib_path = Path(HERE / "ccd_slib.bin")
    if ccd_slib_path.is_file():
        with open(ccd_slib_path, "rb") as fh:
            ccd_slib = rdSubstructLibrary.SubstructLibrary(fh.read())
    else:
        ccd_slib = gen_ccd_slib()
    return ccd_slib