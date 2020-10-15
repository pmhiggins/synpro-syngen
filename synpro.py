#----------0.0 To initialise and define values---------

import sys,os, math, statistics
#import pandas as pd
import numpy as np
from Bio import SeqIO
from io import StringIO
import matplotlib.pyplot as plt

from BioMolecule import BioMolecule

#dGfAABB= -82.38 # kcal mol - Amend J. (2000)
dGfAABB= -344.6779 #KJ mol
#PBB= -20.15 #kcal mol - Amend J. (2000)
PBB= -84.3076 #KJ mol
dGfH2O= -237.3 #kJ mol−1 - McCollom and Amend (2005)
R= 0.008314472 #Gas constant KJ mol K
T= 298.15 # Standard temperature
#------------------ A.-1 To get the information ---------
print('Please type in the name of a fasta file containing a protein sequence')
fname = input('fasta file name: ')
print('Choose one of the three model organisms corresponding to the protein sequence')
celltype= input('Is this a bactetria, yeast or mammalian cell? [bacteria/yeast/mammalian]')

#------------------ A.0 To open a file without biopython ---------

#inFile = open(fname,'r')

#headerList = []
#seqList = []
#currentSeq = ''
#for line in inFile:
#   if line[0] == ">":
#      headerList.append(line[1:].strip())
#      if currentSeq != '':
#         seqList.append(currentSeq)

#      currentSeq = ''
 #  else:
#      currentSeq += line.strip()

#seqList.append(currentSeq)
#-------------------------A.1 To open a file with biopython-----------
protname = []
proteins = []

fasta = open(fname,'r')
for record in SeqIO.parse(fasta,'fasta'):
   protname.append(record.id)
   proteins.append(str(record.seq))
print('there are', len(proteins), 'sequences in this proteome')
#-------------------B. Define the dictionary--------------------------------
#Glycine is considered a 'Special case'
#The dictionary contains:
#                       0= AA 3 letter code
#                       1= dGf for the R group #(KJ mol) - Amend & Helgeson (2000)
#                       2= Molecular weight #of the AA (Da)
#                       3= Concentration #(mol/L)E. Coli - Bennett B. et al 2009 & Park J. et al 2016
#                       4= Concentration #(mol/L)E. Mammalian - Bennett B. et al 2009 & Park J. et al 2016
#                       5= Concentration #(mol/L)E. Yeast - Bennett B. et al 2009 & Park J. et al 2016
#                       6 = AA dGf (KJ mol) - McCollom & Amend (2005) - [AA] (Neutral)
#                       7 = R group dGf (KJ mol) - Equilibrator - [R] Cell
#                       8 = AA dGf (KJ mol) - Equilibrator - [AA] Cell (Charged -1,-2 or 1 and ionic strength of 0.1M)

ReactionGibbs_neu = []
ReactionGibbs_chr = []


#real_ReactionGibbs = []

TK = np.linspace(300,400, num=5)


for T in TK:

    _AABB = BioMolecule('AABB', 74.05866, 4, T=T)
    _PBB = BioMolecule('PBB', 56.04918, 2, T=T)
    _Water = BioMolecule('H2O(l)', 18, 2, T=T)

    dGfAABB= _AABB.std_formation_gibbs
    PBB = _PBB.std_formation_gibbs
    dGfH2O = _Water.std_formation_gibbs


    _ALA = BioMolecule('Alanine(aq)', 89, 7, T=T)
    _ARG = BioMolecule('Arginine(aq)', 174, 14, T=T)
    _ASN = BioMolecule('Asparagine(aq)', 132, 8, T=T)
    _ASP = BioMolecule('Aspartic-Acid(aq)', 133, 7, T=T)
    _CYS = BioMolecule('Cysteine(aq)', 121, 7, T=T)
    _GLN = BioMolecule('Glutamic-Acid(aq)', 148, 9, T=T)
    _GLU = BioMolecule('Glutamine(aq)', 147, 10, T=T)
    _HIS = BioMolecule('Histidine(aq)', 155, 9, T=T)
    _ILE = BioMolecule('Isoleucine(aq)', 131, 13, T=T)
    _LEU = BioMolecule('Leucine(aq)', 131, 13, T=T)
    _LYS = BioMolecule('Lysine(aq)', 146, 14, T=T)
    _MET = BioMolecule('Methionine(aq)', 149, 11, T=T)
    _PHE = BioMolecule('Phenylalanine(aq)', 165, 11, T=T)
    _PRO = BioMolecule('Proline(aq)', 115, 9, T=T)
    _SER = BioMolecule('Serine(aq)', 105, 7, T=T)
    _THR = BioMolecule('Threonine(aq)', 119, 9, T=T)
    _TRP = BioMolecule('Tryptophan(aq)', 181, 12, T=T)
    _TYR = BioMolecule('Tyrosine(aq)', 181, 11, T=T)
    _VAL = BioMolecule('Valine(aq)', 117, 11, T=T)
    _GLY = BioMolecule('Glycine(aq)', 75, 5, T=T)


    AAdict = {
      'A': ['ALA', _ALA.std_formation_R, 89, 2.6e-03, 6.98e-3, 2.23e-2,
        _ALA.std_formation_gibbs,
        257.47,
        -87.2],
      'R': ['ARG', _ARG.std_formation_R, 174, 5.7e-4, 2.55e-4, 2.18e-2,
        _ARG.std_formation_gibbs,
        683.17,
        338.5],
      'N': ['ASN', _ASN.std_formation_R, 132, 5.1e-4, 2.15e-4, 5.69e-3,
        _ASN.std_formation_gibbs,
        143.077,
        -201.6],
      'D': ['ASP', _ASP.std_formation_R, 133, 4.2e-3, 1.49e-2, 6.29e-3,
        _ASP.std_formation_gibbs,
        -97.92,
        -442.6],
      'C': ['CYS', _CYS.std_formation_R, 121, 3.7e-4, 8.40E-5, 3.7e-4,
        _CYS.std_formation_gibbs,
        289.67,
        -55.0],
      'Q': ['GLN', -_GLN.std_formation_R, 146, 3.8e-3, 1.72e-2, 3.55e-2,
        _GLN.std_formation_gibbs,
        224.37,
        -120.3],
      'E': ['GLU', _GLU.std_formation_R, 147, 9.6e-2, 6.38e-2, 3.91e-2,
        _GLU.std_formation_gibbs,
        0.12,
        -344.8],
      'H': ['HIS', _HIS.std_formation_R, 155, 6.8e-5, 4.10e-4, 6.76e-5,
        _HIS.std_formation_gibbs,
        529.77,
        185.1],
      'I': ['ILE', _ILE.std_formation_R, 131, 1.5e-4, 1.76e-3, 3.53e-4,
        _ILE.std_formation_gibbs,
        527.87,
        183.2],
      'L': ['LEU', _LEU.std_formation_R, 131, 1.5e-4, 1.76e-3, 3.53e-4,
        _LEU.std_formation_gibbs,
        519.57,
        174.9],
      'K': ['LYS', _LYS.std_formation_R, 146, 4.1e-4, 5.06e-4, 5.16e-3,
        _LYS.std_formation_gibbs,
        608.57,
        263.9],
      'M': ['MET', _MET.std_formation_R, 149, 1.5e-4, 6.39e-4, 1.91e-4,
        _MET.std_formation_gibbs,
        287.77,
        -56.9],
      'F': ['PHE', _PHE.std_formation_R, 165, 1.8e-5, 8.40e-4, 2.73e-4,
        _PHE.std_formation_gibbs,
        583.57,
        238.9],
      'P': ['PRO', _PRO.std_formation_R, 115, 3.9e-4, 1.23e-3, 1.36e-3,
        _PRO.std_formation_gibbs,
        423.97,
        79.3],
      'S': ['SER', _SER.std_formation_R, 105, 1.13e-3, 4.86e-3, 3.87e-3,
        _SER.std_formation_gibbs,
        106.47,
        -238.2],
      'T': ['THR', _THR.std_formation_R, 119, 1.26e-3, 6.69e-3, 6.69e-3,
        _THR.std_formation_gibbs,
        180.67,
        -164.0],
      'W': ['TRP', _TRP.std_formation_R, 204, 1.2e-5, 1.80e-4, 5.55e-5,
        _TRP.std_formation_gibbs,
        719.27,
        374.6],
      'Y': ['TYR', _TYR.std_formation_R, 181, 2.9e-5, 9.38e-4, 2.48e-4,
        _TYR.std_formation_gibbs,
        420.0779,
        75.4],
      'V': ['VAL', _VAL.std_formation_R, 117, 4e-3, 1.51e-3, 2.50e-3, -357.2,
        _VAL.std_formation_gibbs,
        431.97,
        87.3],
      'G': ['GLY', _GLY.std_formation_R, 75, 5.82e-4, 3.71E-3, 3.71E-3,
        _GLY.std_formation_gibbs,
        167.47,
        -177.2]}

    #---------------C. To get the count of AA in each sequence-------------------
    #To append the values of each protein, define a counter outside the loop
    sumavgaa=0
    protlen=0
    totmw=0
    proteinGs_chr=[]
    proteinGs_neu=[]
    dGfs_chr=[]
    dGfs_neu=[]

    #---- First count the AA in the protein, this block will identify one or/and three letter AA codes
    for protein in proteins:
        CAla=(protein.count('A') + protein.count('ALA'))
        CArg=(protein.count('R') + protein.count('ARG'))
        CAsn=(protein.count('N') + protein.count('ASN'))
        CAsp=(protein.count('D') + protein.count('ASP'))
        CCys=(protein.count('C') + protein.count('CYS'))
        CGln=(protein.count('Q') + protein.count('GLN'))
        CGlu=(protein.count('E') + protein.count('GLU'))
        CGly=(protein.count('G') + protein.count('GLY'))
        CHis=(protein.count('H') + protein.count('HIS'))
        CIle=(protein.count('I') + protein.count('ILE'))
        CLeu=(protein.count('L') + protein.count('LEU'))
        CLys=(protein.count('K') + protein.count('LYS'))
        CMet=(protein.count('M') + protein.count('MET'))
        CPhe=(protein.count('F') + protein.count('PHE'))
        CPro=(protein.count('P') + protein.count('PRO'))
        CSer=(protein.count('S') + protein.count('SER'))
        CThr=(protein.count('T') + protein.count('THR'))
        CTrp=(protein.count('W') + protein.count('TRP'))
        CTyr=(protein.count('Y') + protein.count('TYR'))
        CVal=(protein.count('V') + protein.count('VAL'))

        #Number of AA in the protein #(gets rid of the special characters and helps identifying specific aa)
        nAA= (CAla  + CArg + CAsn + CAsp
        + CCys + CGln + CGlu + CGly
        + CHis + CIle + CLeu + CLys
        + CMet + CPhe + CPro + CSer
        + CThr + CTrp + CTyr + CVal)
        #Number Glycines in the protein
        nGly= CGly

        #To get the molecular weight of each protein
        protmw = (
        (CAla * AAdict['A'][2]) + (CArg * AAdict['R'][2]) + (CAsn * AAdict['N'][2]) + (CAsp * AAdict['D'][2])
        + (CCys * AAdict['C'][2]) + (CGln * AAdict['Q'][2]) + (CGlu * AAdict['E'][2]) + (CGly * AAdict['G'][2])
        + (CHis * AAdict['H'][2]) + (CIle * AAdict['I'][2]) + (CLeu * AAdict['L'][2]) + (CLys * AAdict['K'][2])
        + (CMet * AAdict['M'][2]) + (CPhe * AAdict['F'][2]) + (CPro * AAdict['P'][2]) + (CSer * AAdict['S'][2])
        + (CThr * AAdict['T'][2]) + (CTrp * AAdict['W'][2]) + (CTyr * AAdict['Y'][2]) + (CVal * AAdict['V'][2])- (nAA*18))

        """Take away nAA H2O weights"""

        #For the average molecular weight of all the amino acids in each protein
        avmw= protmw/nAA
        #print('**avmw** is the average molecular weight of this protein is=', avmw, 'Da')

    #Append the values for:
        sumavgaa+= avmw #1. For the average molecular weight of the amino acids for all the proteins
        protlen+= nAA #2. For the proteins lenght
        totmw+= protmw #3. For the total weight of the proteins

        #D ---To obtain the Gibbs free energy of formation for each protein (dGf[P]) at different conditions ---
        #dGf[P] = dGf[PBB] + ((n - ngly - 1) * dGf[AABB]) + summatory of the dGf[Ri] of each amino acid + dGf[GLY]
        dGfPBB = (nAA - nGly - 1) * PBB
        dGfPBB_dGfAABB = dGfPBB + dGfAABB

        #-----------------------------------Charged---------------------------------

        #sigmAAR_chr will get the dGf of the R group of each amino acid contained in the protein with a net charge of -2, -1 or 1 depending on the AA
        sigmAAR_chr = (
        (CAla * AAdict['A'][7]) + (CArg * AAdict['R'][7]) + (CAsn * AAdict['N'][7]) + (CAsp * AAdict['D'][7])
        + (CCys * AAdict['C'][7]) + (CGln * AAdict['Q'][7]) + (CGlu * AAdict['E'][7]) + (CGly * AAdict['G'][7])
        + (CHis * AAdict['H'][7]) + (CIle * AAdict['I'][7]) + (CLeu * AAdict['L'][7]) + (CLys * AAdict['K'][7])
        + (CMet * AAdict['M'][7]) + (CPhe * AAdict['F'][7]) + (CPro * AAdict['P'][7]) + (CSer * AAdict['S'][7])
        + (CThr * AAdict['T'][7]) + (CTrp * AAdict['W'][7]) + (CTyr * AAdict['Y'][7]) + (CVal * AAdict['V'][7]))

        #sigmAA wil ge the dGf of each amino acid
        sigmAA_chr = (
        (CAla * AAdict['A'][8]) + (CArg * AAdict['R'][8]) + (CAsn * AAdict['N'][8]) + (CAsp * AAdict['D'][8])
        + (CCys * AAdict['C'][8]) + (CGln * AAdict['Q'][8]) + (CGlu * AAdict['E'][8]) + (CGly * AAdict['G'][8])
        + (CHis * AAdict['H'][8]) + (CIle * AAdict['I'][8]) + (CLeu * AAdict['L'][8]) + (CLys * AAdict['K'][8])
        + (CMet * AAdict['M'][8]) + (CPhe * AAdict['F'][8]) + (CPro * AAdict['P'][8]) + (CSer * AAdict['S'][8])
        + (CThr * AAdict['T'][8]) + (CTrp * AAdict['W'][8]) + (CTyr * AAdict['Y'][8]) + (CVal * AAdict['V'][8]))

        #For the entire equation dGf[P]
        dGfP_chr= dGfPBB_dGfAABB + sigmAAR_chr
        dGfs_chr.append(dGfP_chr)
        #print('The standard free energy of formation for this protein is',dGfP, 'kJ mol-1')

        #---------E. The Gibbs free energy of reaction dGr[P]-----------------------
        dGr_chr= (((nAA-1)*dGfH2O) + dGfP_chr) - sigmAA_chr #KJ mol
        proteinGs_chr.append(dGr_chr)

        #---------------------------------Neutral charge----------------------------
        #sigmAAR_neu will get the dGf of the R group of each amino acid contained in the protein with a net charge of 0 on the AA
        sigmAAR_neu = (
        (CAla * AAdict['A'][1]) + (CArg * AAdict['R'][1]) + (CAsn * AAdict['N'][1]) + (CAsp * AAdict['D'][1])
        + (CCys * AAdict['C'][1]) + (CGln * AAdict['Q'][1]) + (CGlu * AAdict['E'][1]) + (CGly * AAdict['G'][1])
        + (CHis * AAdict['H'][1]) + (CIle * AAdict['I'][1]) + (CLeu * AAdict['L'][1]) + (CLys * AAdict['K'][1])
        + (CMet * AAdict['M'][1]) + (CPhe * AAdict['F'][1]) + (CPro * AAdict['P'][1]) + (CSer * AAdict['S'][1])
        + (CThr * AAdict['T'][1]) + (CTrp * AAdict['W'][1]) + (CTyr * AAdict['Y'][1]) + (CVal * AAdict['V'][1]))

        #sigmAA wil ge the dGf of each amino acid
        # was 8, now 6
        sigmAA_neu = (
        (CAla * AAdict['A'][6]) + (CArg * AAdict['R'][6]) + (CAsn * AAdict['N'][6]) + (CAsp * AAdict['D'][6])
        + (CCys * AAdict['C'][6]) + (CGln * AAdict['Q'][6]) + (CGlu * AAdict['E'][6]) + (CGly * AAdict['G'][6])
        + (CHis * AAdict['H'][6]) + (CIle * AAdict['I'][6]) + (CLeu * AAdict['L'][6]) + (CLys * AAdict['K'][6])
        + (CMet * AAdict['M'][6]) + (CPhe * AAdict['F'][6]) + (CPro * AAdict['P'][6]) + (CSer * AAdict['S'][6])
        + (CThr * AAdict['T'][6]) + (CTrp * AAdict['W'][6]) + (CTyr * AAdict['Y'][6]) + (CVal * AAdict['V'][6]))

        #For the entire equation dGf[P]
        dGfP_neu= dGfPBB_dGfAABB + sigmAAR_neu
        dGfs_neu.append(dGfP_neu)

        #print('The standard free energy of formation for this protein is',dGfP, 'kJ mol-1')
        #-----------------------------------E. The Gibbs free energy of reaction dGr[P]-------------------------------------
        dGr_neu= ((((nAA-1)*dGfH2O) + dGfP_neu) - sigmAA_neu) #KJ mol
        proteinGs_neu.append(dGr_neu)

    #----------------------------------F. Number of proteins per type of cell (estimation)----------------------------
    #As proposed by Milo R. Insights & Perspectives (2013)
    # Protein density = cell density (1.1 g/mL) * Water content * protein fraction of dry mass
    # p𝑐= 𝐶𝑑 ∗ 𝑤 ∗ 𝑃𝑓
    #For E. coli / bacteria
    pdbacteria= (1.1)*(1-0.7)*(0.55)
    #For Yeast
    pdyeast= (1.1)*(1-0.6)*(0.40)
    #For human (Albe K. J.  theor. Biol. (1990) 143, 163-195)
    pdhuman= (1.1)*(1-0.65)*(0.52)
    # Protein mass = protein density+ Avogadro's number (6.022e23) + Cell volume (10e-12 mL/um3)
    # pm= 𝑐 ∗ 𝐴𝑁 ∗ 𝐶𝑣
    #For E. coli / bacteria
    pmbacteria= ((pdbacteria) * (6.022*10**23) * (10**-12))
    #For Yeast
    pmyeast= ((pdyeast) * (6.022*10**23) * (10**-12))
    #For human
    pmhuman= ((pdhuman) * (6.022*10**23) * (10**-12))
    #Average protein size of this cell
    avP= totmw/ len(proteins) #average size in Da
    av_lenght= protlen/ len(proteins) #average lenght

    #------------------------------------- Protein concentration per type of cell-----------------------------------------------------
    #Protein number = ((protein mass / (Average lenght, in amino acids, of the proteins * Average amino acid molecular weight in Da)) x (characteristic volume in each type of cell)
    if celltype == "bacteria":
        host_volume= 1 #um3
        host_mass= 2.85e-13 #dry grams
        #---------------------------------------------------------
        """Protein concentration in E. coli"""
        moles = pmbacteria/(avP*6.022e23) # estimated moles of protein per um3 in E. coli
        mol_L = moles/(1e-15) #estimated moles per litre
        protein_conc = (pmbacteria/avP) #Protein concentration in 1um3 (proteins/um3)
        number_prot = (protein_conc * host_volume)#number of proteins (an E. coli cell is 1 um3)
        proteome_weight = number_prot * avP #Total weight of the proteome in Da

    elif celltype == "mammalian":
        host_volume=3000 #um3
        host_mass= 1.23e-09 #dry grams
        #---------------------------------------------------------
        """Protein concentration in a mammalian cell"""
        moles = pmmammalian/(avP*6.022e23) # estimated moles of protein per um3 in E. coli
        mol_L = moles/(1e-15) #estimated moles per litre
        protein_conc= (pmmammalian/avP) #Protein concentration in 1um3 (proteins/um3)
        number_prot= (protein_conc * host_volume) #number of proteins in one mammalian cell
        proteome_weight= number_prot * avP #Total weight of the proteome in Da

    elif celltype == "yeast":
        host_volume= 30 #um3
        host_mass= 3.16e-11 #dry grams
        #---------------------------------------------------------
        """Protein concentration in yeast"""
        moles = pmyeast/(avP*6.022e23) # estimated moles of protein per um3 in E. coli
        mol_L = moles/(1e-15) #estimated moles per litre
        protein_conc= (pmyeast/avP) #Protein concentration in 1um3 (proteins/um3)
        number_prot= (protein_conc * host_volume) #number of proteins in one yeast cell
        proteome_weight = number_prot * avP #Total weight of the proteome in Da

    #---------------------------------------------------------
    """Real values"""
    #real_mol_L= len(proteins)/(6.022e23*1e-15)

    #--------------------------------------------------------H. For the molal Gibbs Free Energy (Neutral Charge)-----------------------------------
    PC = mol_L
    #real_PC = real_mol_L

    #make a list of the free energies for each protein
    reactionGs_neu = []
    #real_reactionGs = []
    for protein, stdG in zip(proteins, proteinGs_neu):
        lnQ = math.log(PC/len(proteins))
        for a in protein:
            if a in AAdict.keys():
                if celltype == "bacteria":
                    lnQ -= math.log(AAdict[a][3])
                elif celltype == "mammalian":
                    lnQ -= math.log(AAdict[a][4])
                elif celltype == "yeast":
                    lnQ -= math.log(AAdict[a][5])
        reactionGs_neu.append((stdG+(R*(T)*lnQ)))

    #take the average, which is in kJ/mol
    avgReactionG_kJmol_neu = statistics.mean(reactionGs_neu)
    # convert to kJ / dry g
    # 1e-15*PC/len(proteins) is the mean number of moles of protein per cell
    # 2.8e-13 converts to dry mass

    avgReactionG_kJdryg_neu = avgReactionG_kJmol_neu * 1/(avP)
    ReactionGibbs_neu.append(avgReactionG_kJdryg_neu)


    #--------------------------------------------------------H. For the molal Gibbs Free Energy (Charged)-----------------------------------

    #make a list of the free energies for each protein
    reactionGs_chr = []
    #real_reactionGs = []
    for protein, stdG in zip(proteins, proteinGs_chr):
        lnQ = math.log(PC/len(proteins))
        for a in protein:
            if a in AAdict.keys():
                if celltype == "bacteria":
                    lnQ -= math.log(AAdict[a][3])
                elif celltype == "mammalian":
                    lnQ -= math.log(AAdict[a][4])
                elif celltype == "yeast":
                    lnQ -= math.log(AAdict[a][5])
        reactionGs_chr.append((stdG+(R*(T)*lnQ)))

    #take the average, which is in kJ/mol
    avgReactionG_kJmol_chr = statistics.mean(reactionGs_chr)
    # convert to kJ / dry g
    # 1e-15*PC/len(proteins) is the mean number of moles of protein per cell
    # 2.8e-13 converts to dry mass

    avgReactionG_kJdryg_chr = avgReactionG_kJmol_chr * 1/(avP)
    ReactionGibbs_chr.append(avgReactionG_kJdryg_chr)

print('--------These are the theoretical values:--------------------------------')
"""
print('The total theoretical size of the proteome according to the input sequence is ', '{:.2e}'.format(number_prot), 'proteins in the cell, with a diversity of', len(proteins),'proteins')
print('There are', moles,'moles of proteins in the proteome with a concentration of',mol_L,' mol per litre')
print('The total weight of this proteome is', proteome_weight, ' Da')
print('The estimated energy to build this proteome at 298.15 K with neutral AA is:',ReactionGibbs_neu[45],'KJ mol and with charged AA is:',ReactionGibbs_chr[45],'KJ mol')
print('The average dGf for each protein with neutral AAs is:', statistics.mean(dGfs_neu),'and', statistics.mean(dGfs_chr), 'KJ mol with charged AAs')
print('The average dGr for each protein with neutral AAs is:', avgReactionG_kJdryg_neu, 'and', avgReactionG_kJdryg_chr, 'KJ mol with charged AAs')
"""
#-----------------------------------------------------I. To plot in different temperatures---------------------------------
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy calculated at two conditions for a bacteria cell')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy calculated at two conditions for a mammalian cell')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy calculated at two conditions for a yeast cell')

ax.plot(TK, ReactionGibbs_neu, label='Proteome energy with neutrual AA', c='g', linewidth=3)
ax.plot(TK, ReactionGibbs_chr, label='Proteome energy with AA charged', c='r', linewidth=3)
#ax.plot(TK, real_ReactionGibbs, label='Input proteins synthesis', c='r', linewidth=3)
ax.set_ylabel(r'Energetic cost [kJ per dry g]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)


ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()
