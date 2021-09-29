"""Syngen. version G.04.07.21"""
"""-0-"""
import sys,os, math, statistics, re
import numpy as np
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.Seq import transcribe
import matplotlib.pyplot as plt
from BioMolecule import BioMolecule
#--------------------------------
from time import sleep
import time
import progressbar
from tqdm import tqdm, trange
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    AdaptiveETA, FileTransferSpeed, FormatLabel, Percentage, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer, UnknownLength
"""Progress bars"""
widgets0=['Parsing through',' ', progressbar.Counter(),' ', 'DNA sequence(s)',
]
widgets1=['Generating the complementary strand', progressbar.Bar('§'), progressbar.Counter(), '(', progressbar.Percentage(), ' complete)',
    ' [', progressbar.Timer(), '] '
    ' (', progressbar.ETA(), ') ',
]
widgets2=['Transcribing the DNA sequence', progressbar.Bar('|'), progressbar.Counter(), '(', progressbar.Percentage(), ' complete)',
    ' [', progressbar.Timer(), '] '
    ' (', progressbar.ETA(), ') ',
]
widgets3=['Analysing sequences', progressbar.Bar('|'), progressbar.Counter(), '(', progressbar.Percentage(), ' complete)',
    ' [', progressbar.Timer(), '] '
    ' (', progressbar.ETA(), ') ',
]
widgets4=['dsDNA', progressbar.Bar('|'), progressbar.Counter(), '(', progressbar.Percentage(), ' complete)',
    ' [', progressbar.Timer(), '] '
    ' (', progressbar.ETA(), ') ',
]
widgets5=['mRNAs', progressbar.Bar('|'), progressbar.Counter(), '(', progressbar.Percentage(), ' complete)',
    ' [', progressbar.Timer(), '] '
    ' (', progressbar.ETA(), ') ',
]


"""-A-"""
#--------Transcription------
def transcription(ORF):
    mRNAlst = [] # save the fragments in a list
    mRNA= ''
    foundStart = False
    foundEnd = False
    for i in range(0, len(ORF), 3):
        codon= "".join(ORF[i:i+3]) # Define a codon
        if codon == 'ATG' and not foundStart:
            foundStart = True
            foundEnd = False # start recording until it finds TAG.
        if foundStart and not foundEnd:
            cc=transcribe(codon)
            mRNA = mRNA + cc
        if codon == 'TAA': # First possible stop codon
            foundEnd = True
        elif codon == 'TGA': # Second possible stop codon
            foundEnd = True
        elif codon == 'TAG': # Third possible stop codon
            foundEnd = True
            foundStart = False # Start looking again
            mRNAlst.append(mRNA) # save what we have to the list
            mRNA='' # reset for starting again

    # if we get to the end with no TAG, add everything saved to mRNA so far
    if mRNA != '':
        mRNAlst.append(mRNA)
    return(mRNAlst)


"""-B-"""
#----Input sequence---------
print('Enter the name of a DNA fasta file (e.g. Ecoli_genome.fasta)')
fname= input('File name: ')
file_name= "%s.fasta" % fname
print('---------------------------------Model cell-----------------------------')
celltype= input('Is this a bactetria, yeast or mammalian cell? [bacteria/yeast/mammalian]')
if celltype == 'mammalian':
    TK = np.linspace(275,400, num=2)
else:
    TK = np.linspace(275,400, num=5)
R= 0.008314472 #Gas constant KJ mol K

temp= input('Do you want to calculate the energy at any specific temperature? [yes/no] ')
if temp == 'yes':
    stemperature= int(input('Type the temperature (in kelvin) you would like to print results for: '))
    TK=np.append(TK,stemperature)
    TK=np.sort(TK)
    selected_temp= np.nonzero(TK == stemperature)[0][0]


print(' ')
"""-C-"""
#---Parsing with biopython----
Seqname = []
genelist = []
records= []
raw_mRNAs=[]
mRNAs=[]
transcripts = []
noncoding_RNAs= []

fasta = open(fname,'r')
for record in progressbar.progressbar(SeqIO.parse(fasta,'fasta'), widgets=widgets0):
    time.sleep(0.01)
    Seqname.append(record.id) #Name of the gene
    genelist.append(str(record.seq)) #Genes in a list of strings
    records.append((record.seq)) # Genes in a list of int to be used by biopyhton

print('There are', len(records), 'DNA sequences in this file.')

"""-D-"""
#---Complement and transcribe------
#"""Complement"""
print(' ')
print('DNA double chain')
for record in progressbar.progressbar(records, widgets=widgets1):
    time.sleep(0.01)
    complementary=record.complement()
    genelist.append(str(complementary))

#"""Transcribe"""
print(' ')
print('Transcription')
for gene in progressbar.progressbar(records, widgets=widgets2):
    time.sleep(0.01)
    mRNA = transcription(gene)
    if mRNA != '':
        raw_mRNAs.extend(mRNA) #  Add a list on the end to be cleared
mRNAs= [x for x in raw_mRNAs if x]

#"""Clean the transcripts"""
#To clean the list by getting rid of the false positives, zeroes, empty lists, short and very long mRNAs
for rna in mRNAs:
    if len(rna) > 18 and len(rna) < 750: #Shea J. Andrews and Joseph A. Rothnagel. Nature Reviews  2014. "The smallest translated coding sORF described so far is six codons long."
        transcripts.append(rna)
    else:
        noncoding_RNAs.append(rna)
print('There are',len(transcripts), 'sequences fulfilling the characteristics for a mRNA, whilst',len(noncoding_RNAs), 'non-coding short or long RNA(s) that will not be considered for the calculations')


print(' ')
if celltype == "bacteria":
    print('The program will calculate the energy of the', len(records),'sequences in the file + the complementary strand (',len(genelist),'single strands in total) and',len(transcripts),'transcripts in this bacteria cell')

elif celltype == "mammalian":
     print('The program will calculate the energy of the', len(records),'sequences in the file + the complementary strand (',len(genelist),'single strands in total) and',len(transcripts),'transcripts in this mammalian cell')

elif celltype == "yeast":
     print('The program will calculate the energy of the', len(records),'sequences in the file + the complementary strand (',len(genelist),'single strands in total) and',len(transcripts),'transcripts in this yeast cell')
print(' ')
"""-E-"""
#----Nucleotides library---
#Start the for loop to perform calculations at different temperatures
print(' ')
print('Energy calculation at different temperatures')
#DNA chain method Energy
DNA_dGf_av=[]
DNA_dGr_av=[]
DNA_ReactionGibbs=[]
DNA_ReactionGibbs_t=[]
#DNA chain method Energy - Oxic
DNA_oxic_dGf_av=[]
DNA_oxic_dGr_av=[]
DNA_oxic_ReactionGibbs=[]
#DNA block method 1 energy
DNA_block_dGf_av=[]
DNA_block_dGr_av=[]
DNA_block_ReactionGibbs=[]
#DNA block method 2 energy
DNA_NS_dGf_av=[]
DNA_NS_dGr_av=[]
DNA_NS_ReactionGibbs=[]
#RNA chain method nergy
RNA_dGf_av=[]
RNA_dGr_av=[]
RNA_ReactionGibbs = []
RNA_ReactionGibbs_t=[]
#RNA block method 1 energy
RNA_block_dGf_av=[]
RNA_block_dGr_av=[]
RNA_block_ReactionGibbs = []
#DNA block method 2 energy
RNA_NS_dGf_av=[]
RNA_NS_dGr_av=[]
RNA_NS_ReactionGibbs=[]
#Membrane
membrane_dGf=[]
membrane_dGr=[]



for T in TK:
    print('Formation, reaction and molar Gibbs energy at',T,'K')
    # Gibbs energy of formation
    #Inorganic compounds (KJ/mol) #McCollom & Amend 2005
    dGf_H=0
    dGf_e=0

    #miscelaneous
    _Water = BioMolecule('H2O(l)', 18, 2, T=T)
    _H2 = BioMolecule('H2(g)', 2.01, 2, T=T)
    _NH4 = BioMolecule('NH4+', 18.03, 4, T=T, z=1)


    # Building blocks for DNA and RNA
    '''A'''
    # BioMolecules are declared with (name, Mr, no. H, T, charge)
    #Base
    _A = BioMolecule('Adenine(aq)', 135, 5, T=T)
    #DNA Nucleoside
    _dAS = BioMolecule('Deoxyadenosine(aq)', 251.24192, 13, T=T)
    #RNA nucleoside
    _AS = BioMolecule('Adenosine(aq)', 267.24132, 13, T=T)
    #DNA Nucleotide
    _dAMP2 = BioMolecule('d+H2AMP-(aq)', 300.24806, 14, T=T)
    _dAMP = BioMolecule('dHAMP-', 330.213881, 13, T=T, z=-1)
    #ion
    _dAMPion = BioMolecule('dAMP2-', 329.205941, 12, T=T, z=-2)
    #RNA nucleotide
    #_AMP = BioMolecule('+H2AMP-(aq)', 347.221221, 14, T=T)
    _AMP = BioMolecule('HAMP-', 346.205341, 13, T=T, z=-1)
    #ion
    _AMPion = BioMolecule('AMP2-', 345.205341, 12, T=T, z=-2)

    '''C'''
    # BioMolecules are declared with (name, Mr, no. H, T, charge)
    #Base
    _C = BioMolecule('Cytosine(aq)', 111.102, 5, T=T)
    #DNA Nucleoside
    _dCS = BioMolecule('Deoxycytidine(aq)', 227.21722, 13, T=T)
    #RNA nucleoside
    _CS = BioMolecule('Cytidine(aq)', 243.21662, 13, T=T)
    #DNA Nucleotide
    #_dCMP = BioMolecule('d+H2CMP-(aq)', 307.197121, 14, T=T)
    _dCMP = BioMolecule('dHCMP-', 306.189181, 13, T=T, z=-1)
    #ion
    _dCMPion = BioMolecule('dCMP2-', 305.181241, 12, T=T, z=-2)
    #RNA nucleotide
    #_CMP = BioMolecule('+H2CMP-(aq)', 323.196521, 14, T=T)
    _CMP = BioMolecule('HCMP-', 322.188581, 13, T=T, z=-1)
    #ion
    _CMPion = BioMolecule('CMP2-', 321.180641, 12, T=T, z=-2)

    '''G'''
    # BioMolecules are declared with (name, Mr, no. H, T, charge)
    #Base
    _G = BioMolecule('Guanine(aq)', 151.1261, 5, T=T)
    #DNA nucleoside
    _dGS = BioMolecule('Deoxyguanosine(aq)', 267.24132, 13, T=T)
    #RNA nucleoside
    _GS = BioMolecule('Guanosine(aq)', 283.24072, 13, T=T)
    #DNA nucleotide
    #_dGMP = BioMolecule('d+H2GMP-(aq)', 316.24746, 14, T=T)
    _dGMP = BioMolecule('dHGMP-', 346.213281, 13, T=T, z=-1)
    #ion
    _dGMPion = BioMolecule('dGMP2-', 345.205341, 12, T=T, z=-2)
    #RNA nucleotide
    #_GMP = BioMolecule('+H2GMP-(aq)', 363.220621, 14, T=T)
    _GMP = BioMolecule('HGMP-', 362.212681, 13, T=T, z=-1)
    #ion
    _GMPion = BioMolecule('GMP2-', 361.204741, 12, T=T, z=-2)

    '''T'''
    # BioMolecules are declared with (name, Mr, no. H, T, charge)
    #Base
    _T = BioMolecule('Thymine(aq)', 126.11334, 6, T=T)
    #DNA nucleoside
    _dTS = BioMolecule('Deoxythymidine(aq)', 242.22856, 14, T=T)
    #DNA nucleotide
    #_dTMP = BioMolecule('d+H2TMP-(aq)', 322.208461, 15, T=T)
    _dTMP =BioMolecule('dHTMP-', 321.200521, 14, T=T, z=-1)
    #ion
    _dTMPion =BioMolecule('dTMP2-', 320.192581, 13, T=T, z=-2)

    '''U'''
    # BioMolecules are declared with (name, Mr, no. H, T, charge)
    #Base
    _U = BioMolecule('Uracil(aq)', 112.086, 4, T=T)
    #RNA nucleoside
    _US = BioMolecule('Uridine(aq)', 244.20138, 12, T=T)
    #RNA nucleotide
    #_UMP = BioMolecule('+H2UMP-(aq)', 324.181281, 13, T=T)
    _UMP = BioMolecule('HUMP-', 323.173341, 12, T=T, z=-1)
    #ion
    _UMPion = BioMolecule('UMP2-', 322.165401, 11, T=T, z=-2)

    # sugars
    _rib = BioMolecule('Ribose(aq)', 150.1299, 10, T=T)
    _drib = BioMolecule('Deoxyribose(aq)', 134.1305, 10, T=T)
    _R5P = BioMolecule('Ribose-5-Phosphate-2', 228.093921, 9, T=T, z=-2)

    # Phosphates
    _H2PO4 = BioMolecule('H2PO4-', 96.987241, 2, T=T, z=-1)
    _H3PO4 = BioMolecule('H3PO4(aq)', 97.995181, 3, T=T)
    _HPO4 = BioMolecule('HPO4--', 95.979301, 1, T=T, z=-2)
    _PO4 = BioMolecule('PO4---', 94.971361, 0, T=T, z=-3)
    _P2O7 = BioMolecule('P2O7----', 173.943322, 0, T=T, z=-4)
    # _ribphos = BioMolecule('Ribose-5-Phosphate-2', 150.1299, 10, T=T, z=-2)

    #Amino acids
    _C5H10N2O3 = BioMolecule('Glutamine(aq)', 146.1445, 10, T=T)
    _C3H7NO3 = BioMolecule('Serine(aq)', 105.09258, 7, T=T)

    # Other metabolites
    _C6H12O6 = BioMolecule('Glucose(aq)', 180.15588, 12, T=T)
    _C3H2O4 = BioMolecule('Malonate--', 102.04558, 2, T=T, z=-2)
    _C10H16N5O13P3 = BioMolecule('+H4ATP-(aq)', 507.181023, 16, T=T)
    _C10H15N5O10P2 = BioMolecule('+H3ADP1-(aq)', 427.201122, 15, T=T)
    _O2 = BioMolecule('O2(aq)', 31.9988, 0, T=T)
    _CO2 = BioMolecule('CO2(aq)', 44.0095, 0, T=T)



    '''Values lab'''
    # formation energies
    glucose= _C6H12O6.stdbio_formation_gibbs
    ATP= _C10H16N5O13P3.stdbio_formation_gibbs
    ADP= _C10H15N5O10P2.stdbio_formation_gibbs
    glutamine= _C5H10N2O3.stdbio_formation_gibbs
    serine= _C3H7NO3.stdbio_formation_gibbs
    O2= _O2.stdbio_formation_gibbs
    CO2= _CO2.stdbio_formation_gibbs
    H2O = _Water.stdbio_formation_gibbs
    PO4= _PO4.stdbio_formation_gibbs
    H2PO4= _H2PO4.stdbio_formation_gibbs
    Phosphate= _HPO4.stdbio_formation_gibbs
    malonate= _C3H2O4.std_formation_gibbs
    dGf_R5P= _R5P.std_formation_gibbs

    #Group contribution
    #Molecules or functional groups
    # estimate that the -OH group on a sugar is similar to the difference
    # between ribose and deoxyribose (one OH bond on the ring).
    dGf_OH = _rib.stdbio_formation_gibbs - _drib.stdbio_formation_gibbs
    oxygen_m= _H2.stdbio_formation_gibbs - _Water.stdbio_formation_gibbs
    hydrogen= dGf_OH - oxygen_m

    # Bonds (KJ/mol)
    #ester_bond= -22.175 # Phosphodiester bond At 25 °C, pH 7 (Hydrolysis) Dickson K. 2000
    dGf_PhdB= 22.175 # (Reversed) Phosphodiester bond At 25 °C, pH 7  Dickson K. 2000
    #dGf_PhdB= 103.76 # Phosphodiester bond At 65 °C, pH 8 (Hydrolysis) Molina R. 2015
    ester_bond= _P2O7.stdbio_formation_gibbs - ((2*_HPO4.stdbio_formation_gibbs) - H2O)
    e_bond= (_rib.stdbio_formation_gibbs + Phosphate) - (dGf_R5P - H2O)
    glyco_bond= _dAS.stdbio_formation_gibbs - ((_A.stdbio_formation_gibbs + _drib.stdbio_formation_gibbs)- H2O)
    hydro_bond= dGf_OH - _Water.stdbio_formation_gibbs
    Hbond=  (dGf_OH + hydrogen) - _Water.stdbio_formation_gibbs


    #Nucleotides or nucleosides
    damp= (_HPO4.stdbio_formation_gibbs + ester_bond + _drib.stdbio_formation_gibbs + glyco_bond + _A.stdbio_formation_gibbs) - (2*H2O)
    NuSi= (_HPO4.stdbio_formation_gibbs + ester_bond + _dAS.stdbio_formation_gibbs) - (H2O)
    _dAMP_dGf= _dAMP.stdbio_formation_gibbs


    #Backbones
    #DNA phosphate backbone = Phosphate + ester bond + Deoxyribose
    phosba= (
    Phosphate + ester_bond + _drib.stdbio_formation_gibbs - H2O)
    #RNA phosphate backbone = Phosphate + ester bond + Ribose
    RNA_phosba= (
    Phosphate + ester_bond + _rib.stdbio_formation_gibbs - H2O)
    #Phosohate + ester bond
    phosphoester= (
    Phosphate + ester_bond - H2O)

    #Metabolites energies
    glucose_6_p= ((ATP - ADP) + glucose)
    G3P= (((glucose*0.5) + (ATP*1)) - ((ADP*1)))
    choline = (((serine*3) + (ATP*7.4) + (H2O*0.2)) - ((ADP*7.8) + (PO4*6.6) + (dGf_H*8.8)))
    glycerol = (((glucose_6_p*0.5) + (ATP*2.5) + (O2*1.5)) - ((ADP*2.5) + (PO4*3)))
    pyruvate= (((G3P*1) + (ADP*2) + (PO4*1)) - ((ATP*2) + (H2O*1)))
    palmitate= (((pyruvate*8) + (malonate*2)) - ((CO2*14) + (H2O*2)))
    oleate= (((pyruvate*9) + (malonate*2)) - ((CO2*15) + (H2O*3)))

    #Phospholipids
    POPC= (((glucose*5) + (serine*1) + (ATP*1) + (pyruvate*10) + (malonate*1)) - ((ADP*1) + (H2O*14) + (CO2*24)))

    damp= (_HPO4.stdbio_formation_gibbs + ester_bond + _drib.stdbio_formation_gibbs + glyco_bond + _A.stdbio_formation_gibbs) - (2*H2O)
    NuSi= (_HPO4.stdbio_formation_gibbs + ester_bond + _dAS.stdbio_formation_gibbs) - (H2O)
    print('damp by block method (H12)',damp,'KJ/mol')
    print('damp by lock method 2 (H12)',NuSi,'KJ/mol')
    print('dAMPion (database) (H12)',(_dAMPion.stdbio_formation_gibbs),'KJ/mol')
    print('dAMP (database) (H13)',(_dAMP.stdbio_formation_gibbs),'KJ/mol')
    print('dAMP2 (database)(H14)',(_dAMP2.stdbio_formation_gibbs),'KJ/mol')


    print('Phosphate backbone')
    print('_PO4',_PO4.stdbio_formation_gibbs,'KJ/mol')
    print('Deoxyribose',_drib.stdbio_formation_gibbs,'KJ/mol')
    print('ribose',_rib.stdbio_formation_gibbs,'KJ/mol')
    print(' ')

    print('Bases')
    print('Adenine',_A.stdbio_formation_gibbs,'KJ/mol')
    print('Cytosine',_C.stdbio_formation_gibbs,'KJ/mol')
    print('Guanine',_G.stdbio_formation_gibbs,'KJ/mol')
    print('Thymine',_T.stdbio_formation_gibbs,'KJ/mol')
    print('Uracil',_U.stdbio_formation_gibbs,'KJ/mol')
    print(' ')


    print('DNA nucleosides')
    print('Deoxyadenosine',_dAS.stdbio_formation_gibbs,'KJ/mol')
    print('Deoxycytidine',_dCS.stdbio_formation_gibbs,'KJ/mol')
    print('Deoxyguanosine',_dGS.stdbio_formation_gibbs,'KJ/mol')
    print('Deoxythymidine',_dTS.stdbio_formation_gibbs,'KJ/mol')
    print(' ')

    print('RNA nucleosides')
    print('adenosine',_AS.stdbio_formation_gibbs,'KJ/mol')
    print('cytidine',_CS.stdbio_formation_gibbs,'KJ/mol')
    print('guanosine',_GS.stdbio_formation_gibbs,'KJ/mol')
    print('Uradine',_US.stdbio_formation_gibbs,'KJ/mol')
    print(' ')

    print('DNA nucleotides')
    print('dAMP',_dAMP.stdbio_formation_gibbs,'KJ/mol')
    print('dCMP',_dCMP.stdbio_formation_gibbs,'KJ/mol')
    print('dGMP',_dGMP.stdbio_formation_gibbs,'KJ/mol')
    print('dTMP',_dTMP.stdbio_formation_gibbs,'KJ/mol')
    print(' ')

    print('RNA nucleotides')
    print('AMP',_AMP.stdbio_formation_gibbs,'KJ/mol')
    print('CMP',_CMP.stdbio_formation_gibbs,'KJ/mol')
    print('GMP',_GMP.stdbio_formation_gibbs,'KJ/mol')
    print('UMP',_UMP.stdbio_formation_gibbs,'KJ/mol')
    print(' ')

    print('Metabolites')
    print('Glucose',glucose,'KJ/mol')
    print('ATP',ATP,'KJ/mol')
    print('ADP',ADP,'KJ/mol')
    print('Glutamine',glutamine,'KJ/mol')
    print('O2',O2,'KJ/mol')
    print('CO2',CO2,'KJ/mol')
    print('glucose_6_p',glucose_6_p,'KJ/mol')
    print('P from ATP',(ATP - ADP), 'KJ/mol')





    # NB for this we CAN'T change it with temperature - so this brings uncertainty

    """-F-"""
    #---Dictionaries----
    #The dictionaries for both DNA  and RNA contain (In the same order):
    #                       0= Nucleotide name
    #                       1= Base name
    #                       2= Nucleotide's molecular weight (Da)
    #                       3= Nucleotide's absolute intracellular concentration (mol/L - E. coli)
    #                       4= Nucleotide's absolute intracellular concentration (mol/L - Yeast)
    #                       5= Nucleotide's absolute intracellular concentration (mol/L - Mammalian)
    #                       6= dGf - base (database)
    #                       7= dGf - Nucleoside (database)
    #                       8= dGf - Nucleotide (database)
    #                       9= dGf - Nucleotide ion

    #Dictionary for DNA
    DNAdict = {
      'A': ['dAMP','A',_dAMP.Mr,8.8e-6, 8.8e-6, 1.68e-5,
        _A.stdbio_formation_gibbs, _dAS.stdbio_formation_gibbs,
        _dAMP.stdbio_formation_gibbs, _dAMPion.stdbio_formation_gibbs],
      'C': ['dCMP','C',_dCMP.Mr,3.71e-05, 3.71e-05, 3.71e-05,
        _C.stdbio_formation_gibbs, _dCS.stdbio_formation_gibbs,
        _dCMP.stdbio_formation_gibbs,_dCMPion.stdbio_formation_gibbs],
      'T': ['dTMP','T', _dTMP.Mr,1.18e-5,1.18e-5, 1.18e-5,
        _T.stdbio_formation_gibbs, _dTS.stdbio_formation_gibbs,
        _dTMP.stdbio_formation_gibbs, _dTMPion.stdbio_formation_gibbs],
      'G': ['dGMP','G',_dGMP.Mr,5.1e-5,5.1e-5,5.1e-5,
        _G.stdbio_formation_gibbs, _dGS.stdbio_formation_gibbs,
        _dGMP.stdbio_formation_gibbs, _dGMPion.stdbio_formation_gibbs]}

    #Dictionary for RNA
    RNAdict = {
      'A': ['AMP','A',_AMP.Mr,2.8e-4,8.12e-5,4.23e-5,
        _A.stdbio_formation_gibbs, _AS.stdbio_formation_gibbs,
        _AMP.stdbio_formation_gibbs, _AMPion.stdbio_formation_gibbs],
      'C': ['CMP','C',_CMP.Mr,3.6e-4,5.18e-6, 1.18e-5,
        _C.stdbio_formation_gibbs, _CS.stdbio_formation_gibbs,
        _CMP.stdbio_formation_gibbs, _CMPion.stdbio_formation_gibbs],
      'U': ['UMP','U',_UMP.Mr,1.45E-05,1.45E-05,1.45E-05,
        _U.stdbio_formation_gibbs, _US.stdbio_formation_gibbs,
        _UMP.stdbio_formation_gibbs, _UMPion.stdbio_formation_gibbs],
      'G': ['GMP','G',_GMP.Mr,2.37e-5,1.02e-5, 1.81e-5,
        _G.stdbio_formation_gibbs, _GS.stdbio_formation_gibbs,
        _GMP.stdbio_formation_gibbs, _GMPion.stdbio_formation_gibbs]}

    #Dictionary for the building block concentrations:

    #       Base  | Phosphate | Deoxyribose | Ribose | DNA nucleoside | RNA nucleoside
    # E coli   0        1           2           3           4           5
    # S cerevisiae 6    7           8           9           10          11
    # Mammalian  12     13          14          15          16          17
    BB_conc = {
      'A': [1.47E-06, 2.39E-02, 3.03E-04, 1.52E-04, 2.82E-06, 1.31E-07,
            1.47E-06, 4.93E-02, 3.03E-04, 1.52E-04, 2.82E-06, 1.31E-07,
            1.47E-06, 5.83E-03, 3.03E-04, 7.83E-05, 2.82E-06, 1.31E-07],
      'C': [2.59E-06, 2.39E-02, 3.03E-04, 1.52E-04, 2.82E-06, 1.41E-05,
            2.59E-06, 4.93E-02, 3.03E-04, 1.52E-04, 2.82E-06, 1.41E-05,
            2.59E-06, 5.83E-03, 3.03E-04, 7.83E-05, 2.82E-06, 1.41E-05],
      'G': [1.88E-04, 2.39E-02, 3.03E-04, 1.52E-04, 5.22E-07, 1.62E-06,
            1.88E-04, 4.93E-02, 3.03E-04, 1.52E-04, 5.22E-07, 1.62E-06,
            1.88E-04, 5.83E-03, 3.03E-04, 7.83E-05, 5.22E-07, 1.35E-06],
      'T': [1.47E-06, 2.39E-02, 3.03E-04, 1.52E-04, 5.22E-07, 0,
            1.47E-06, 4.93E-02, 3.03E-04, 1.52E-04, 5.22E-07, 0,
            1.47E-06, 5.83E-03, 3.03E-04, 7.83E-05, 5.22E-07, 0],
      'U': [2.10E-03, 2.39E-02, 3.03E-04, 1.52E-04, 0, 2.09E-03,
            2.10E-03, 4.93E-02, 3.03E-04, 1.52E-04, 0, 2.09E-03,
            2.10E-03, 5.83E-03, 3.03E-04, 7.83E-05, 0, 2.09E-03]}

    #         G-6-P  |  ATP  |   ADP   |  Glutamite | CO2  |
    # E coli    0       1          2        3           4
    # S cerevisiae 5    6          7        8           9
    # Mammalian  10     11         12       13          14
    metabolites_conc = {
      'A': [7.88E-03, 9.63E-03, 5.55E-04, 3.81E-03, 7.52E-05,
            5.31E-03, 1.93E-03, 4.88E-04, 3.55E-02, 8.16E-05,
            6.75E-04, 4.67E-03, 5.69E-04, 1.62E-02, 7.63E-03],
      'C': [7.88E-03, 9.63E-03, 5.55E-04, 3.81E-03, 7.52E-05,
            5.31E-03, 1.93E-03, 4.88E-04, 3.55E-02, 8.16E-05,
            6.75E-04, 4.67E-03, 5.69E-04, 1.62E-02, 7.63E-03],
      'G': [7.88E-03, 9.63E-03, 5.55E-04, 3.81E-03, 7.52E-05,
            5.31E-03, 1.93E-03, 4.88E-04, 3.55E-02, 8.16E-05,
            6.75E-04, 4.67E-03, 5.69E-04, 1.62E-02, 7.63E-03],
      'T': [7.88E-03, 9.63E-03, 5.55E-04, 3.81E-03, 7.52E-05,
            5.31E-03, 1.93E-03, 4.88E-04, 3.55E-02, 8.16E-05,
            6.75E-04, 4.67E-03, 5.69E-04, 1.62E-02, 7.63E-03],
      'U': [7.88E-03, 9.63E-03, 5.55E-04, 3.81E-03, 7.52E-05,
            5.31E-03, 1.93E-03, 4.88E-04, 3.55E-02, 8.16E-05,
            6.75E-04, 4.67E-03, 5.69E-04, 1.62E-02, 7.63E-03]}




    """-G-"""
    #---DNA---
    """ G.0 - Initialise values"""
    #--- Number of nucleotides on each strand and total weight---
    tot_N=0
    totweight=0
    #---Energy of each strand (Chain method)---
    DNA_dGf=[]
    DNA_dGr=[]
    #---Energy of each strand (Chain method - Oxic)---
    DNA_oxic_dGf=[]
    DNA_oxic_dGr=[]
    #---Energy of each strand  (Block method 1)---
    DNA_block_dGf=[]
    DNA_block_dGr=[]
    #---Energy of each strand  (Block method 2)---
    DNA_NS_dGf=[]
    DNA_NS_dGr=[]
    #Molar Gibbs free energies for all the genes
    DNA_reactionGs = []
    DNA_oxic_reactionGs = []
    DNA_block_reactionGs= []
    DNA_NS_reactionGs = []


    """G.1 - Nucleotides count"""
    #---Iterate through the DNA sequences----
    for sequence in progressbar.progressbar(genelist, widgets=widgets4):
        #Count for each nucleotide
        CA=sequence.count('A')
        CT=sequence.count('T')
        CC=sequence.count('C')
        CG=sequence.count('G')
        #Total number of nucleotides
        Nn= len(sequence)
        tot_N+= Nn

        """G.2 - Weight of the genome"""
        #The product of each nucleotide by its corresponding mass
        gene_weight = (
        (CA * DNAdict['A'][2]) + (CT * DNAdict['T'][2]) + (CC * DNAdict['C'][2]) + (CG * DNAdict['G'][2]) - ((Nn-1)*18))
        totweight+= gene_weight

        """G.3 - Energy of the nucleotides"""
        # At standard biological conditions, with pH 7 and ionic strength 0.

        #sigmaB will get the sum of the bases formation energy in the sequence
        sigmaB = (
        (CA * DNAdict['A'][6]) + (CT * DNAdict['T'][6]) + (CC * DNAdict['C'][6]) + (CG * DNAdict['G'][6]))

        #sigmaB will get the sum of the bases formation energy in the sequence
        sigmaNS = (
        (CA * DNAdict['A'][7]) + (CT * DNAdict['T'][7]) + (CC * DNAdict['C'][7]) + (CG * DNAdict['G'][7]))

        #sigmaN will get the sum of the nucleotides formation energy in the sequence
        sigmaN = (
            (CA * DNAdict['A'][8]) + (CT * DNAdict['T'][8]) + (CC * DNAdict['C'][8]) + (CG * DNAdict['G'][8]))

        #sigmaNion will get the sum of energy of formation of the 2- charged  nucleotides in the sequence (e.g. 0 OH group on the phosphate)
        sigmaNion = (
        (CA * DNAdict['A'][9]) + (CT * DNAdict['T'][9]) + (CC * DNAdict['C'][9]) + (CG * DNAdict['G'][9]))

        #sigmaH will get the total number of hydrogen bonds
        sigmaH = (
        (CA * hydro_bond) + (CT * hydro_bond) + (CC * (hydro_bond*1.5)) + (CG * (hydro_bond*1.5)))




        """G.4 - Formation Energy - Block method 1"""
        #By block - Linking ester bond + phosphate backbone + glycosidic bond + base - H2O)
        """Gibbs Formation energy"""
        #Nion_block= ((Nn*_HPO4.stdbio_formation_gibbs) + (Nn*ester_bond) + (Nn*_drib.stdbio_formation_gibbs) + (Nn*glyco_bond) + sigmaB) - (Nn*(2*H2O))
        Nion_block= (Nn*(phosba + glyco_bond - H2O)) + sigmaB
        DNA_dGf_byblock = ((Nn-1)*(ester_bond - dGf_OH)) + Nion_block
        DNA_block_dGf.append(DNA_dGf_byblock)
        """Gibbs Reaction Energy"""
        DNA_dGr_byblock = ((((Nn-1)*H2O) + DNA_dGf_byblock) - sigmaN)
        DNA_block_dGr.append(DNA_dGr_byblock)

        """G.5 - Formation Energy - Block method 2"""
        #By block method 2 - Linking ester bond + phosphoester + nucleoside
        """Gibbs Formation energy"""
        DNA_dGf_byblock2 = ((Nn-1)*(ester_bond - dGf_OH)) + (Nn*phosphoester) + sigmaNS
        DNA_NS_dGf.append(DNA_dGf_byblock2)
        #print('damp2 (block method)',nucleoside,'KJ/mol')
        #print('dAMPion (database)',((Nn)*_dAMPion.stdbio_formation_gibbs),'KJ/mol')
        """Gibbs Reaction Energy"""
        DNA_dGr_byblock2 = (((Nn-1)*H2O) + DNA_dGf_byblock2) - sigmaN
        DNA_NS_dGr.append(DNA_dGr_byblock2)
        #DNA_block_totdGr += DNA_block_dGr

        """G.6 - Formation Energy - Chain method"""
        """Gibbs Formation energy"""
        #Chain method - Linking ester bond + Nucleotide ion
        # sum of the nucleotide 2- ions, with n-1 Phosphodiester bonds made and an OH- removed from the sugar.
        #DNA_dGf_chain = ((Nn-1)*(dGf_PhdB - dGf_OH)) + sigmaNion
        DNA_dGf_chain = ((Nn-1)*(ester_bond - dGf_OH)) + sigmaNion
        DNA_dGf.append(DNA_dGf_chain)
        """Gibbs Reaction Energy"""
        DNA_dGr_chain = ((((Nn-1)*H2O) + DNA_dGf_chain) - sigmaN)
        DNA_dGr.append(DNA_dGr_chain)

        """G.7 - Oxic method"""
        """Gibbs Reaction Energy"""
        #DNA_dGf_oxic = ((Nn-1)*(ester_bond - dGf_OH)) + ((glucose_6_p*0.5)+(ATP*0.5)+(glutamine*2.5)+(O2*6.5))
        DNA_dGf_oxic = ((Nn-1)*(ester_bond - dGf_OH)) + sigmaNion
        DNA_oxic_dGf.append(DNA_dGf_chain)

        DNA_dGr_oxic = (((Nn-1)*H2O) + DNA_dGf_oxic + (Nn*(((glucose_6_p*0.5)+(ATP*0.5)+(glutamine*2.5)+(O2*6.5))-((ADP*0.5)+(H2O*9)+(CO2*5.5))))) - sigmaN
        #DNA_dGr_oxic = (((Nn-1)*H2O) + DNA_dGf_oxic + ((Nn)*((H2O*9)+(ADP*0.5)+(CO2*5.5)))) - sigmaN
        DNA_oxic_dGr.append(DNA_dGr_oxic)





    """-H-"""
    #Average values for molecular weight and gene lenght
    #DNA oncentration. for the reaction quotient
    nucleotide_avmw= totweight/tot_N #Da
    gene_meanlength= tot_N/len(genelist) #Average lenght of the strands
    avG= totweight/len(genelist) #average gene size in Daltons
    bp= tot_N/2
    #host_volume= 1 #um3
    # To get the moles
    #if celltype == "bacteria":  #15 mg/mL with avG (15mg/ml from bionumbers [Elowitz MB, Surette MG, Wolf PE, Stock JB, Leibler S. Protein mobility in the cytoplasm of Escherichia coli. J Bacteriol. 1999 Jan181(1):197-203])
    #    conc_Da_ml = 8.55/(1.66054e-21) # 1.66054e-21 is the conversion Da-mg
    #elif celltype == "mammalian":
    #    conc_Da_ml = 1.88/(1.66054e-21)
    #elif celltype == "yeast":
    #    conc_Da_ml = 0.76/(1.66054e-21)

    #conc_molecules_ml = conc_Da_ml / avG #number molecules per ml
    #mol_ml = conc_molecules_ml/6.022e23
    #DNA_mol_L = mol_ml*1000

    '''By grams'''
    if celltype == "bacteria":
        DNA_grams= 8.54e-15 #g
        cell_moles = DNA_grams/nucleotide_avmw #mol
        cell_volume= 1e-15 #L
    elif celltype == "yeast":
        DNA_grams= 3.25E-14 #g
        cell_moles = DNA_grams/nucleotide_avmw #mol
        cell_volume= 1.1e-13 #L
    elif celltype == "mammalian":
        DNA_grams= 6.01e-12 #g
        cell_moles = DNA_grams/nucleotide_avmw #mol
        cell_volume= 3E-12 #L

    DNA_mol_Lit= cell_moles/cell_volume


    """-I-"""
    #Reaction quotient for molar Gibbs energy of each gene with the chain method
    PC = DNA_mol_Lit

    for gene, stdG in zip(genelist, DNA_dGr):
        DNA_lnQ = math.log(PC/len(genelist))
        for base in gene:
            if base in DNAdict.keys():
                                if celltype == "bacteria":
                                    DNA_lnQ -= math.log(DNAdict[base][3])
                                elif celltype == "yeast":
                                    DNA_lnQ -= math.log(DNAdict[base][4])
                                elif celltype == "mammalian":
                                    DNA_lnQ -= math.log(DNAdict[base][5])

        DNA_reactionGs.append((stdG+(R*(T)*DNA_lnQ)))

    #--------Oxic-----------

    for gene_o, o_stdG in zip(genelist, DNA_oxic_dGr):
        DNA_o_lnQ = (math.log(PC/len(genelist)) + math.log(metabolites_conc[base][2]) + math.log(metabolites_conc[base][4]))
        for base in gene_o:
            if base in metabolites_conc.keys():
                                if celltype == "bacteria":
                                    DNA_o_lnQ -= math.log(metabolites_conc[base][0]) +  math.log(metabolites_conc[base][1]) + math.log(metabolites_conc[base][3])
                                elif celltype == "yeast":
                                    DNA_o_lnQ -= math.log(metabolites_conc[base][5]) +  math.log(metabolites_conc[base][6]) + math.log(metabolites_conc[base][8])
                                elif celltype == "mammalian":
                                    DNA_o_lnQ -= math.log(metabolites_conc[base][10]) +  math.log(metabolites_conc[base][11]) + math.log(metabolites_conc[base][13])

        DNA_oxic_reactionGs.append((o_stdG+(R*(T)*DNA_o_lnQ)))



    #--------By blocks-----------
    #Reaction quotient for molar Gibbs energy of each gene with the Block 1 method

    for geneblock, B_stdG in zip(genelist, DNA_block_dGr):
        DNA_B_lnQ = math.log(PC/len(genelist))
        for base in geneblock:
            if base in BB_conc.keys():
                                if celltype == "bacteria":
                                    DNA_B_lnQ -= math.log(BB_conc[base][0]) +  math.log(BB_conc[base][1]) + math.log(BB_conc[base][2])
                                elif celltype == "yeast":
                                    DNA_B_lnQ -= math.log(BB_conc[base][6]) +  math.log(BB_conc[base][7]) + math.log(BB_conc[base][8])
                                elif celltype == "mammalian":
                                    DNA_B_lnQ -= math.log(BB_conc[base][12]) +  math.log(BB_conc[base][13]) + math.log(BB_conc[base][14])

        DNA_block_reactionGs.append((B_stdG+(R*(T)*DNA_B_lnQ)))

    #Reaction quotient for molar Gibbs energy of each gene with the Block 2 method

    for geneNS, NS_stdG in zip(genelist, DNA_NS_dGr):
        DNA_NS_lnQ = math.log(PC/len(genelist))
        for base in geneNS:
            if base in BB_conc.keys():
                                if celltype == "bacteria":
                                    DNA_NS_lnQ -= math.log(BB_conc[base][1]) + math.log(BB_conc[base][4])
                                elif celltype == "yeast":
                                    DNA_NS_lnQ -= math.log(BB_conc[base][7]) + math.log(BB_conc[base][10])
                                elif celltype == "mammalian":
                                    DNA_NS_lnQ -= math.log(BB_conc[base][13]) + math.log(BB_conc[base][16])

        DNA_NS_reactionGs.append((NS_stdG+(R*(T)*DNA_NS_lnQ)))



    """-J-"""
    #Chain method
    #take the average, which is in kJ/mol
    DNA_avgReactionG_kJmol = statistics.mean(DNA_reactionGs)
    # convert to kJ / dry g
    #avG is the average gene size
    DNA_avgReactionG_kJdryg = DNA_avgReactionG_kJmol * 1/(avG)
    DNA_ReactionGibbs.append(DNA_avgReactionG_kJdryg)
    #Sum all the energies for the average

    #---oxic---
    DNA_oxic_avgReactionG_kJmol = statistics.mean(DNA_oxic_reactionGs)
    DNA_oxic_avgReactionG_kJdryg = DNA_oxic_avgReactionG_kJmol * 1/(avG)
    DNA_oxic_ReactionGibbs.append(DNA_oxic_avgReactionG_kJdryg)


    #---Block method 1---
    DNA_block_avgReactionG_kJmol = statistics.mean(DNA_block_reactionGs)
    DNA_block_avgReactionG_kJdryg = DNA_block_avgReactionG_kJmol * 1/(avG)
    DNA_block_ReactionGibbs.append(DNA_block_avgReactionG_kJdryg)
    #Sum all the energies for the average

    #---Block method 2---
    DNA_NS_avgReactionG_kJmol = statistics.mean(DNA_NS_reactionGs)
    DNA_NS_avgReactionG_kJdryg = DNA_NS_avgReactionG_kJmol * 1/(avG)
    DNA_NS_ReactionGibbs.append(DNA_NS_avgReactionG_kJdryg)
    #Sum all the energies for the average

    #single cell Values
    DNA_avgReactionG_kJmol_t = statistics.mean(DNA_reactionGs)
    DNA_avgReactionG_kJdryg_t = DNA_avgReactionG_kJmol_t * (1/(avG) * DNA_grams)
    DNA_ReactionGibbs_t.append(DNA_avgReactionG_kJdryg_t)






    """-K-"""
    #---RNA---
    """ K.0 - Initialise values"""
    #--- Number of nucleotides on each strand and total weight---
    transcriptome_lenght=0
    transcriptome_totweight=0
    #---Energy of each strand (Chain method)---
    RNA_dGf=[]
    RNA_dGr=[]
    #----RNA energy of each strand  (Block method 1)---
    RNA_block_dGf=[]
    RNA_block_dGr=[]
    #---RNA energy of each strand  (Block method 2)---
    RNA_NS_dGf=[]
    RNA_NS_dGr=[]
    #Molar Gibbs free energies of all the genes
    RNA_reactionGs = []
    RNA_block_reactionGs= []
    RNA_NS_reactionGs = []

    """K.1 - Nucleotides count"""
    #---Iterate through the DNA sequences----
    for m in progressbar.progressbar(transcripts, widgets=widgets5):
        time.sleep(0.01)
        #Total number of RNA nucleotides
        if len(m) != 0: #To get rid of blank lists
            RNA_Nn= len(m)
            transcriptome_lenght+= RNA_Nn
            #Number of each nucleotide
            RA=m.count('A')
            RU=m.count('U')
            RC=m.count('C')
            RG=m.count('G')

            """K.2 - Weight of the genome"""
            #To get the product of each nucleotides weight by the number of nucleotides
            RNA_weight = (
            (RA * RNAdict['A'][2]) + (RU * RNAdict['U'][2]) + (RC * RNAdict['C'][2]) + (RG * RNAdict['G'][2]) - ((RNA_Nn-1)*18))
            """We take away the H2O weight for every nucleotide"""
            transcriptome_totweight+= RNA_weight

            """K.3 - Energy of the nucleotides"""
            # we use standard biological conditions, with pH 7 and ionic strength 0.

            #RNA_sigmaB will get the energy of formation of the bases in the sequence at cellular conditions
            RNA_sigmaB = (
            (RA * RNAdict['A'][6]) + (RU * RNAdict['U'][6]) + (RC * RNAdict['C'][6]) + (RG * RNAdict['G'][6]))

            #sigmaB will get the sum of the bases formation energy in the sequence
            RNA_sigmaNS = (
            (CA * DNAdict['A'][7]) + (CT * DNAdict['T'][7]) + (CC * DNAdict['C'][7]) + (CG * DNAdict['G'][7]))

            #RNA_sigmaN will get the energy of formation of the 1- Charged nucleotides in the sequence at cellular conditions
            RNA_sigmaN = (
            (RA * RNAdict['A'][8]) + (RU * RNAdict['U'][8]) + (RC * RNAdict['C'][8]) + (RG * RNAdict['G'][8]))

            #RNA_sigmaNion will get the sum of energy of formation of the 2- charged  nucleotides in the sequence (e.g. 0 OH group on the phosphate, for group contribution)
            RNA_sigmaNion = (
            (RA * RNAdict['A'][9]) + (RU * RNAdict['U'][9]) + (RC * RNAdict['C'][9]) + (RG * RNAdict['G'][9]))


            """K.4 - RNA formation Energy - Block method 1"""
            #By block - Linking ester bond + phosphate backbone + glycosidic bond + base - H2O)
            """Gibbs Formation energy"""
            #Nion_block= ((Nn*_HPO4.stdbio_formation_gibbs) + (Nn*ester_bond) + (Nn*_drib.stdbio_formation_gibbs) + (Nn*glyco_bond) + sigmaB) - (Nn*(2*H2O))
            RNA_Nion_block= (RNA_Nn*(RNA_phosba + glyco_bond - H2O)) + RNA_sigmaB
            RNA_dGf_byblock = ((RNA_Nn-1)*(ester_bond - dGf_OH)) + RNA_Nion_block
            RNA_block_dGf.append(RNA_dGf_byblock)
            """Gibbs Reaction Energy"""
            RNA_dGr_byblock = ((((RNA_Nn-1)*H2O) + RNA_dGf_byblock) - RNA_sigmaN)
            RNA_block_dGr.append(RNA_dGr_byblock)

            """G.5 - RNA formation Energy - Block method 2"""
            #By block method 2 - Linking ester bond + phosphoester + nucleoside
            """Gibbs Formation energy"""
            #RNA_dGf_byblock2 = ((RNA_Nn-1)*(ester_bond - dGf_OH)) + (RNA_Nn*phosphoester) + RNA_sigmaNS
            RNA_dGf_byblock2 = ((RNA_Nn-1)*(ester_bond - dGf_OH)) + (RNA_Nn*((_HPO4.stdbio_formation_gibbs + ester_bond) - H2O)) + RNA_sigmaNS
            RNA_NS_dGf.append(RNA_dGf_byblock2)
            #print('damp2 (block method)',nucleoside,'KJ/mol')
            #print('dAMPion (database)',((Nn)*_dAMPion.stdbio_formation_gibbs),'KJ/mol')
            """Gibbs Reaction Energy"""
            RNA_dGr_byblock2 = (((RNA_Nn-1)*H2O) + RNA_dGf_byblock2) - RNA_sigmaN
            RNA_NS_dGr.append(RNA_dGr_byblock2)
            #DNA_block_totdGr += DNA_block_dGr

            """K.6 - RNA Formation Energy - Chain method"""
            """Gibbs Formation energy"""
            #Chain method - Linking ester bond + Nucleotide ion
            # sum of the nucleotide 2- ions, with n-1 Phosphodiester bonds made and an OH- removed from the sugar.
            #DNA_dGf_chain = ((Nn-1)*(dGf_PhdB - dGf_OH)) + sigmaNion
            RNA_dGf_chain = ((RNA_Nn-1)*(ester_bond - dGf_OH)) + RNA_sigmaNion
            RNA_dGf.append(RNA_dGf_chain)
            """Gibbs Reaction Energy"""
            RNA_dGr_chain = ((((RNA_Nn-1)*H2O) + RNA_dGf_chain) - RNA_sigmaN)
            RNA_dGr.append(RNA_dGr_chain)


    """-L-"""
    #get average values for molecular weight and lenght
    RNA_nucleotide_avmw= transcriptome_totweight/transcriptome_lenght #Da
    transcript_meanlength= transcriptome_lenght/len(transcripts) #Average lenght of the transcript
    transcript_meanweight= transcriptome_totweight/len(transcripts) #average transcript size in Daltons

    '''By grams'''
    if celltype == "bacteria":
        RNA_grams= 5.84e-14
        cell_RNA_moles = RNA_grams/RNA_nucleotide_avmw #mol
    elif celltype == "yeast":
        RNA_grams= 1.46E-12
        cell_RNA_moles = RNA_grams/RNA_nucleotide_avmw #mol
    elif celltype == "mammalian":
        RNA_grams= 2.10E-11
        cell_RNA_moles = RNA_grams/RNA_nucleotide_avmw #mol

    RNA_mol_Lit= cell_RNA_moles/cell_volume


    """-M-"""
    #Reaction quotient for the entire transcriptome
    RNA_PC = RNA_mol_Lit

    #Chain method
    for transcript, RNA_stdG in zip(transcripts, RNA_dGr):
        RNA_lnQ = math.log(RNA_PC/len(transcripts))
        for base in transcript:
            if base in RNAdict.keys():
                                if celltype == "bacteria":
                                    RNA_lnQ -= math.log(RNAdict[base][3])
                                elif celltype == "yeast":
                                    RNA_lnQ -= math.log(RNAdict[base][4])
                                elif celltype == "mammalian":
                                    RNA_lnQ -= math.log(RNAdict[base][5])

        RNA_reactionGs.append((RNA_stdG+(R*(T)*RNA_lnQ)))

    #--------By block method 1----------
    #Reaction quotient for the Gibbs energy by building block

    for b_transcript, RNA_B_stdG in zip(transcripts, RNA_block_dGr):
        RNA_B_lnQ = math.log(PC/len(transcripts))
        for base in b_transcript:
            if base in BB_conc.keys():
                                if celltype == "bacteria":
                                    RNA_B_lnQ -= math.log(BB_conc[base][0] * BB_conc[base][1] * BB_conc[base][3])
                                elif celltype == "yeast":
                                    RNA_B_lnQ -= math.log(BB_conc[base][6] * BB_conc[base][7] * BB_conc[base][9])
                                elif celltype == "mammalian":
                                    RNA_B_lnQ -= math.log(BB_conc[base][12] * BB_conc[base][13] * BB_conc[base][15])

        RNA_block_reactionGs.append((RNA_B_stdG+(R*(T)*RNA_B_lnQ)))

    #--------By block method 2----------
    #Reaction quotient for the Gibbs energy by building block

    for NS_transcript, RNA_NS_stdG in zip(transcripts, RNA_NS_dGr):
        RNA_NS_lnQ = math.log(PC/len(transcripts))
        for base in NS_transcript:
            if base in BB_conc.keys():
                                if celltype == "bacteria":
                                    RNA_NS_lnQ -= math.log(BB_conc[base][1] * BB_conc[base][15])
                                elif celltype == "yeast":
                                    RNA_NS_lnQ -= math.log(BB_conc[base][7] * BB_conc[base][11])
                                elif celltype == "mammalian":
                                    RNA_NS_lnQ -= math.log(BB_conc[base][13] * BB_conc[base][17])

        RNA_NS_reactionGs.append((RNA_NS_stdG+(R*(T)*RNA_NS_lnQ)))




    """-N-"""
    #Molar Gibbs Energy
    #--Chain method
    #take the average, which is in kJ/mol
    RNA_avgReactionG_kJmol = statistics.mean(RNA_reactionGs)
    # convert to kJ / dry g
    #avG is the average gene size
    #Bucket
    RNA_avgReactionG_kJdryg = RNA_avgReactionG_kJmol * 1/(transcript_meanweight)
    RNA_ReactionGibbs.append(RNA_avgReactionG_kJdryg)

    #---By block method 1---
    RNA_avgReactionG_block_kJmol = statistics.mean(RNA_block_reactionGs)
    RNA_avgReactionG_block_kJdryg = RNA_avgReactionG_block_kJmol * 1/(transcript_meanweight)
    RNA_block_ReactionGibbs.append(RNA_avgReactionG_block_kJdryg)

    #---By block method 2---
    RNA_avgReactionG_NS_kJmol = statistics.mean(RNA_NS_reactionGs)
    RNA_avgReactionG_NS_kJdryg = RNA_avgReactionG_NS_kJmol * 1/(transcript_meanweight)
    RNA_NS_ReactionGibbs.append(RNA_avgReactionG_NS_kJdryg)

    #Single cell values
    RNA_avgReactionG_kJmol_t = statistics.mean(RNA_reactionGs)
    RNA_avgReactionG_kJdryg_t = RNA_avgReactionG_kJmol_t * (1/(transcript_meanweight)*RNA_grams)
    RNA_ReactionGibbs_t.append(RNA_avgReactionG_kJdryg_t)



    """Averages"""
    #DNA chain energy
    DNA_dGf_av.append(statistics.mean(DNA_dGf)) #Formation Gibbs Energy
    DNA_dGr_av.append(statistics.mean(DNA_dGr)) #Reaction Gibbs Energy
    #DNA block method 1
    DNA_block_dGf_av.append(statistics.mean(DNA_block_dGf)) #Formation Gibbs Energy
    DNA_block_dGr_av.append(statistics.mean(DNA_block_dGr)) #Reaction Gibbs Energy
    #DNA block method 2
    DNA_NS_dGf_av.append(statistics.mean(DNA_NS_dGf)) #Formation Gibbs Energy
    DNA_NS_dGr_av.append(statistics.mean(DNA_NS_dGr)) #Reaction Gibbs Energy
    #RNA chain
    RNA_dGf_av.append(statistics.mean(RNA_dGf)) #Formation Gibbs Energy
    RNA_dGr_av.append(statistics.mean(RNA_dGr)) #Reaction Gibbs Energy
    #RNA block method 1
    RNA_block_dGf_av.append(statistics.mean(RNA_block_dGf)) #Formation Gibbs Energy
    RNA_block_dGr_av.append(statistics.mean(RNA_block_dGr)) #Reaction Gibbs Energy
    #RNA block method 2
    RNA_NS_dGf_av.append(statistics.mean(RNA_NS_dGf)) #Formation Gibbs Energy
    RNA_NS_dGr_av.append(statistics.mean(RNA_NS_dGr)) #Reaction Gibbs Energy



    """------------------------------"""
    """Membrane and cell walls"""

    POPC_w= 760.091 #g/mol

    if celltype == "bacteria":
        PL_grams= 2.59E-14
        cell_volume= 1e-15 #L
        PL_moles = PL_grams/POPC_w #mol
        membrane_number= PL_moles * 6.0221e23
    elif celltype == "yeast":
        PL_grams= 1.14E-12
        cell_volume= 1.1e-13 #L
        PL_moles = PL_grams/POPC_w #mol
        membrane_number= PL_moles * 6.0221e23
    elif celltype == "mammalian":
        PL_grams= 6.83E-11
        cell_volume= 13E-12 #L
        PL_moles = PL_grams/POPC_w #mol
        membrane_number= PL_moles * 6.0221e23


    membrane_mol_Lit= PL_moles/cell_volume


    """X - Membrane Energy """
    """Gibbs Formation energy"""
    membrane = membrane_number*((oleate+palmitate+choline+glycerol+H2PO4)-(4*H2O))
    membrane_dGf.append(membrane)
    """Gibbs Reaction Energy"""
    membrane_dGr_chain = (membrane+(membrane_number*(ADP+(H2O*14)+(CO2*24)))) - (membrane_number*((glucose*5)+serine+ATP+(pyruvate*10)+malonate))
    membrane_dGr.append(membrane_dGr_chain)



    M_PC = membrane_mol_Lit

    membrane_lnQ = math.log(M_PC/membrane_number)
    if celltype == "bacteria":
        membrane_lnQ -= math.log(7.88E-03) +  math.log(1.13E-03) + math.log(9.63E-03) + math.log(3.66E-03) + math.log(7.26E-05)
    elif celltype == "yeast":
        membrane_lnQ -= math.log(5.31E-03) +  math.log(3.87E-03) + math.log(1.93E-03) + math.log(9.40E-03) + math.log(7.26E-05)
    elif celltype == "mammalian":
        membrane_lnQ -= math.log(6.75E-04) +  math.log(4.86E-03) + math.log(4.67E-03) + math.log(5.88E-03) + math.log(7.26E-05)

    membrane_reactionGs=(membrane_dGr+(R*(T)*membrane_lnQ))
    membrane_ReactionGibbs = membrane_reactionGs * (1/(POPC_w*membrane_number))


"""-O-"""
#---Printing values---
if temp == 'yes':
    """DNA"""
    print('The energy required for this DNA sequence and the transcripts at',stemperature, 'K is:')
    print(' ')
    print('Genome')
    print('Formation energy')
    print('Average formation Gibbs energy of the genes using chain method:', "%.3f" % DNA_dGf_av[selected_temp], 'KJ/mol')
    print('Reaction energy')
    print('Average reaction Gibbs energy of the genes using chain method:', "%.3f" % DNA_dGr_av[selected_temp], 'KJ/mol')
    print('Molar energies')
    print('Average molar energy of the genome (chain method):', "%.3f" % DNA_ReactionGibbs[selected_temp], 'KJ/g')
    print('Average molar energy of the genome by block:', "%.3f" % DNA_block_ReactionGibbs[selected_temp], 'KJ/g')
    print('Average molar energy of the genome by block 2:', "%.3f" % DNA_NS_ReactionGibbs[selected_temp], 'KJ/g')
    print('Average molar energy oxic conditions', "%.3f" % DNA_oxic_ReactionGibbs[selected_temp], 'KJ/g')
    print('Molar energy of a single cell genome:', DNA_ReactionGibbs_t[selected_temp], 'KJ/g')


    print(' ')
    """RNA"""
    print('Transcriptome')
    print('Formation energy')
    print('Average formation Gibbs energy of the transcripts using chain method:', "%.3f" % RNA_dGf_av[selected_temp], 'KJ/mol')
    print('Reaction energy')
    print('Average reaction Gibbs energy of the transcripts using chain method:', "%.3f" % RNA_dGr_av[selected_temp], 'KJ/mol')
    print('Molar energies')
    print('Molar energy of the transcriptome (chain method):', "%.3f" % RNA_ReactionGibbs[selected_temp], 'KJ/g')
    print('Molar energy of the transcriptome by block:', "%.3f" % RNA_block_ReactionGibbs[selected_temp], 'KJ/g')
    print('Molar energy of the transcriptome by block 2:', "%.3f" % RNA_NS_ReactionGibbs[selected_temp], 'KJ/g')
    print('Molar energy of a single cell genome:', RNA_ReactionGibbs_t[selected_temp], 'KJ/g')

    print(' ')
    """Membrane"""
    print('Membrane')
    print('Formation energy')
    print('Average formation Gibbs energy of the membrane and cell walls:', "%.3f" % membrane_dGf[selected_temp], 'KJ/mol')
    print('Reaction energy')
    print('Average reaction Gibbs energy of the tmembrane and cell walls:', "%.3f" %membrane_dGr[selected_temp], 'KJ/mol')
    print('Molar energy')
    print('Molar energy of the membrane and cell walls:', "%.3f" % membrane_ReactionGibbs[selected_temp], 'KJ/g')



"""-P-"""
#Single cell
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy for a bacteria genome')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy for a mammalian genome')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy for a yeast genome')

ax.plot(TK,DNA_ReactionGibbs_t, label='SC Genome Chain Energy', color='orange', marker='^', linewidth=2)
ax.plot(TK,RNA_ReactionGibbs_t, label='SC Transcriptome Chain Energy', color='green', marker='D', linewidth=2)

ax.set_ylabel(r'Energetic cost [KJ per cell]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

#ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()

#--- 3 methods DNA ----
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy for a bacteria genome with three methods')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy for a mammalian genome with three methods')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy for a yeast genome with three methods')

ax.plot(TK,DNA_ReactionGibbs, label='Genome Chain Energy', color='crimson', marker='1', linewidth=2)
ax.plot(TK,DNA_block_ReactionGibbs, label='Genome block 1 energy', color='deepskyblue', marker='d', linewidth=2)
ax.plot(TK,DNA_NS_ReactionGibbs, label='Genome block 2 energy', color='darkviolet', marker='o', linewidth=2)

ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

#ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()

#--- 3 methods RNA ----
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy for a bacteria transcriptome with three methods')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy for a mammalian transcriptome with three methods')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy for a yeast transcriptome with three methods')

ax.plot(TK,RNA_ReactionGibbs, label='Transcriptome Chain Energy', linestyle='dotted', c='r', linewidth=2)
ax.plot(TK,RNA_block_ReactionGibbs, label='Transcriptome block 1 energy', linestyle='dotted', c='b', linewidth=2)
#ax.plot(TK,RNA_NS_ReactionGibbs, label='Transcriptome block 2 energy', linestyle='dotted', c='y', linewidth=2)

ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

#ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()

#---Chain energy DNA + RNA ----
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy of the genome and transcriptome for a bacteria cell')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy of the genome and transcriptome for a mammalian cell')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy of the genome and transcriptome for a yeast cell')

ax.plot(TK,DNA_ReactionGibbs, label='Genome Chain Energy', color='crimson', marker='1', linewidth=2)
ax.plot(TK,RNA_ReactionGibbs, label='Transcriptome Chain Energy', linestyle='dotted', color='crimson', marker='1')



ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

#ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()

#---DNA all----
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy of the genome and transcriptome for a bacteria cell')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy of the genome and transcriptome for a mammalian cell')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy of the genome and transcriptome for a yeast cell')

ax.plot(TK,DNA_oxic_ReactionGibbs, label='Genome Chain Energy from metabolites', color='gold', marker='1', linewidth=2)
ax.plot(TK,DNA_block_ReactionGibbs, label='Genome block 1 energy', color='deepskyblue', marker='d', linewidth=2)
ax.plot(TK,DNA_NS_ReactionGibbs, label='Genome block 2 energy', color='darkviolet', marker='o', linewidth=2)
ax.plot(TK,DNA_ReactionGibbs, label='Genome Chain Energy', color='crimson', marker='1', linewidth=2)

ax.plot(TK,DNA_ReactionGibbs_t, label='SC Genome Chain Energy', color='orange', marker='^', linewidth=2)

ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

#ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()

#--- Membrane ----
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy for a bacteria genome with three methods')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy for a mammalian genome with three methods')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy for a yeast genome with three methods')

ax.plot(TK,membrane_ReactionGibbs, label='membrane energy', color='red', marker='v', linewidth=2)

ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()
