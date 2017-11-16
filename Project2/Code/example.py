from molecules import Alcohol, Alkane
from molecule_box import MoleculeBox

molecule = Alkane(chainlength=5)
box = MoleculeBox(molecule=molecule,
                  n_molecules=10,
                  box_length=2)
box.save('box.hoomdxml', forcefield_name='trappe-ua', overwrite=True)
