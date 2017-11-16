import mbuild as mb

class CH2(mb.Compound):
    """Defines a united-atom CH2 particle """
    def __init__(self):
        super(CH2, self).__init__()
        self.add(mb.Particle(name='_CH2'))

        self.add(mb.Port(anchor=self[0], orientation=[0, 1, 0],
                         separation=0.075), 'up')
        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], 
                         separation=0.075), 'down')

class CH3(mb.Compound):
    """Defines a united-atom CH3 particle """
    def __init__(self):
        super(CH3, self).__init__()
        self.add(mb.Particle(name='_CH3'))

        self.add(mb.Port(anchor=self[0], separation=0.075), 'up')

class CH4(mb.Compound):
    """Defines a united-atom methane """
    def __init__(self):
        super(CH4, self).__init__()
        self.add(mb.Particle(name='_CH4'))

class Hydroxyl(mb.Compound):
    """Defines a hydroxyl group """
    def __init__(self):
        super(Hydroxyl, self).__init__()
        oxygen = mb.Particle(name='O', pos=[0.0, 0.0, 0.0])
        hydrogen = mb.Particle(name='H', pos=[0.0945, 0.0, 0.0])
        self.add((oxygen, hydrogen))
        self.add_bond((oxygen, hydrogen))
        self.add(mb.Port(anchor=oxygen, orientation=[-1, 0, 0], separation=0.075),
                 label='up')

class Alkane(mb.Compound):
    """A united-atom alkane chain of variable length """
    def __init__(self, chainlength):
        """Initialize a united-atom alkane chain

        Parameters
        ----------
        chainlength : int
            Number of particles in the chain backbone

        """
        super(Alkane, self).__init__()

        if chainlength < 1:
            raise ValueError('Chain length must be greater than one!')
        elif chainlength == 1:
            self.add(CH4())
        elif chainlength == 2:
            ua1 = CH3()
            ua2 = CH3()
            self.add(ua1)
            self.add(ua2)

            mb.force_overlap(ua2, ua2['up'], ua1['up'])
        else:
            ua1 = CH3()
            ua2 = CH3()
            self.add(ua1)
            self.add(ua2)

            chain_middle = mb.Polymer(CH2(), n=chainlength-2,
                                      port_labels=('up', 'down'))
            self.add(chain_middle)
            mb.force_overlap(chain_middle, chain_middle['up'], ua1['up'])
            mb.force_overlap(ua2, ua2['up'], chain_middle['down'])

class Alcohol(mb.Compound):
    """A united-atom alcohol chain of variable length """
    def __init__(self, chainlength):
        """Initialize a united-atom alcohol chain

        Parameters
        ----------
        chainlength : int
            Number of *carbons* in the chain backbone. This does not include
            the OH group

        """
        super(Alcohol, self).__init__()

        hydroxyl = Hydroxyl()
        ch3 = CH3()
        self.add(hydroxyl)
        self.add(ch3)

        if chainlength < 1:
            raise ValueError('Chain length must be greater than one!')
        elif chainlength == 1:
            mb.force_overlap(ch3, ch3['up'], hydroxyl['up'])
        else:
            chain_middle = mb.Polymer(CH2(), n=chainlength-1,
                                      port_labels=('up', 'down'))
            self.add(chain_middle)
            mb.force_overlap(chain_middle, chain_middle['up'], hydroxyl['up'])
            mb.force_overlap(ch3, ch3['up'], chain_middle['down'])
