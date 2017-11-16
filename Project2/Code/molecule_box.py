import mbuild as mb

from molecules import Alcohol, Alkane


class MoleculeBox(mb.Compound):
    """ A box of molecules """
    def __init__(self, molecule, n_molecules, box_length):
        """Initialize MoleculeBox

        Parameters
        ----------
        molecule : mb.Compound
            Molecule to replicate within the box
        n_molecules : int
            Number of molecules to add to the box
        box_length : float
            Length of the box sides (in nm)

        """
        super(MoleculeBox, self).__init__()

        box = mb.Box(lengths=[box_length, box_length, box_length])
        box_of_molecules = mb.fill_box(compound=molecule,
                                       n_compounds=n_molecules,
                                       box=box)
        self.add(box_of_molecules)
