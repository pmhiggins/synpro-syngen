import sys
import os
from reaktoro import *
import math

this_dir = os.path.dirname(__file__)


class BioMolecule:
    """

    class for storing amino acid properties, or any other useful biological molecule really

    """
    name = ''
    Mr = 0. # molecular mass
    conc_mol_per_l = None	# mol/l
    conc_mol_per_cell = None # mol/cell

    gamma = 1.		# activity coefficient

    std_formation_gibbs = None		# J/mol
    std_formation_R = None #J/mol

    #booleans of state
    thermo = True # whether we have the thermodynamic data available

    copies_per_cell = None
    frequency = None
    probrange = None

    #generic counter to use how you please
    counter = 0

    """

    INITIALISATION METHODS

    """

    def __init__(self, name, Mr, N, conc_mol_per_cell=None, conc_mol_per_l=None, thermo=True, gamma=1., T=298.15, P=101325., pH=7., charge=0., I=0. ):
        self.name = name
        self.Mr = Mr

        if thermo: #pass Thermo as False to update thermochemical parameters yourself
            self.GetThermoParams(T, P, N, pH=pH, I=I, z=charge)

        self.thermo = thermo
        self.conc_mol_per_cell = conc_mol_per_cell
        self.conc_mol_per_l = conc_mol_per_l
        self.gamma = gamma

    def GetThermoParams(self, T, P, N, pH=7., I=0., z=0.):

        db = Database(this_dir+"/supcrt07-organics_AABB.xml")
        thermodynamics = Thermo(db)
        self.std_formation_gibbs = thermodynamics.standardPartialMolarGibbsEnergy(T, P, self.name).val *0.001
        AABB = thermodynamics.standardPartialMolarGibbsEnergy(T, P, 'AABB').val *0.001
		# self.std_formation_gibbs = thermodynamics.standardPartialMolarGibbsEnergy(T=T, P=P, species=self.name).val
		# AABB = thermodynamics.standardPartialMolarGibbsEnergy(T=T, P=P, species='AABB').val
        self.std_formation_R = self.std_formation_gibbs - AABB

        self.stdbio_formation_gibbs = BioMolecule.biogibbs(
          self.std_formation_gibbs, T, pH, I, z, N)

        stdbio_AABB = BioMolecule.biogibbs(
          self.std_formation_gibbs, T, pH, I, z, N)

        self.stdbio_formation_R = self.stdbio_formation_gibbs - stdbio_AABB

    def biogibbs(stdGibbs, T, pH, N, I=0., z=0):
        if I == 0.:
            return stdGibbs - (N*0.00831*T*math.log(10**(-pH)))
        else:
            return stdGibbs - (2.91482*z*z*(I**0.5)/(1+(B*(I**0.5)))) - N*(
              - 2.91482*(I**0.5)/(1+(B*(I**0.5)))
              + 0.00831*T*math.log(10**(-pH)))
