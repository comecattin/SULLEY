import symmetry_poltype as symm
from rdkit import Chem
from rdkit.Chem import AllChem

class Poltype():
        def __init__(
                self, mol,
                molstructfname = 'mol.mol',
                peditinfile = 'local_frame.txt',
                paramfile = 'polarize.prm',
            ) -> None:

            self.usepoleditframes = False
            self.usesymtypes = True
            self.totalcharge = 0
            self.prmstartidx = 0
            self.indextotypefile = None
            self.peditinfile = peditinfile
            self.molstructfname = molstructfname
            self.indextompoleframefile = None
            self.forcefield = 'TEST'
            self.updatedsmallmoleculepolarizeprmlib = paramfile



            self.canonicallabel = [ 0 ] * mol.NumAtoms()
            self.localframe1 = [ 0 ] * mol.NumAtoms()
            self.localframe2 = [ 0 ] * mol.NumAtoms()
            
            (
                self.idxtosymclass,
                self.symmetryclass
            ) = symm.gen_canonicallabels(self,mol,None,self.usesymtypes,True)
        
        def WriteToLog(self, message):
            print(message)

        def CheckInputCharge(self, molecule):
            return molecule, 0
        
        def smiles_to_mol_file(self, smiles, file_name):
            # Create a molecule from the SMILES
            mol = Chem.MolFromSmiles(smiles)
            
            if mol is None:
                print("Invalid SMILES string")
                return
            
            # Add hydrogens (optional, but often useful)
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            Chem.AllChem.Compute3DCoords(mol)

            # Write the molecule to a MOL file
            with open(file_name, 'w') as mol_file:
                mol_block = Chem.MolToMolBlock(mol)
                mol_file.write(mol_block)
            print(f"MOL file '{file_name}' created successfully.")