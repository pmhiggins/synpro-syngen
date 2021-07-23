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
#DNA Energy
DNA_dGf_av=[]
DNA_dGr_av=[]
DNA_ReactionGibbs=[]
#DNA Energy by block
DNA_block_dGf_av=[]
DNA_block_dGr_av=[]
DNA_block_ReactionGibbs=[]
#RNA Energy
RNA_dGf_av=[]
RNA_dGr_av=[]
RNA_ReactionGibbs = []
#RNA Energy by block
RNA_block_dGf_av=[]
RNA_block_dGr_av=[]
RNA_block_ReactionGibbs = []

for T in TK:
    print('Formation, reaction and molar Gibbs energy at',T,'K')
    # Gibbs energy of formation
    #Inorganic compounds (KJ/mol) #McCollom & Amend 2005
    dGf_H=0
    dGf_e=0

    _Water = BioMolecule('H2O(l)', 18, 2, T=T)
    dGf_H2O = _Water.stdbio_formation_gibbs

    # bases and nucleotides
    # BioMolecules are declared with (name, Mr, no. H, T, charge)
    _A = BioMolecule('Adenine(aq)', 135, 5, T=T)
    #_AMP = BioMolecule('HAMP-', 346.205341, 13, T=T, z=-1)
    _AMPion = BioMolecule('AMP2-', 345.205341, 12, T=T, z=-2)
    #_dAMP = BioMolecule('dHAMP-', 330.213881, 13, T=T, z=-1)
    _dAMPion = BioMolecule('dAMP2-', 329.205941, 12, T=T, z=-2)
    #--
    _AMP = BioMolecule('+H2AMP-(aq)', 347.221221, 14, T=T)
    _dAMP = BioMolecule('d+H2AMP-(aq)', 300.24806, 14, T=T)

    _C = BioMolecule('Cytosine(aq)', 111.102, 5, T=T)
    #_CMP = BioMolecule('HCMP-', 322.188581, 13, T=T, z=-1)
    _CMPion = BioMolecule('CMP2-', 321.180641, 12, T=T, z=-2)
    #_dCMP = BioMolecule('dHCMP-', 306.189181, 13, T=T, z=-1)
    _dCMPion = BioMolecule('dCMP2-', 305.181241, 12, T=T, z=-2)
    #--
    _CMP = BioMolecule('+H2CMP-(aq)', 323.196521, 14, T=T)
    _dCMP = BioMolecule('d+H2CMP-(aq)', 307.197121, 14, T=T)

    _G = BioMolecule('Guanine(aq)', 151.1261, 5, T=T)
    #_GMP = BioMolecule('HGMP-', 362.212681, 13, T=T, z=-1)
    _GMPion = BioMolecule('GMP2-', 361.204741, 12, T=T, z=-2)
    #_dGMP = BioMolecule('dHGMP-', 346.213281, 13, T=T, z=-1)
    _dGMPion = BioMolecule('dGMP2-', 345.205341, 12, T=T, z=-2)
    #--
    _dGMP = BioMolecule('d+H2GMP-(aq)', 316.24746, 14, T=T)
    _GMP = BioMolecule('+H2GMP-(aq)', 363.220621, 14, T=T)

    _T = BioMolecule('Thymine(aq)', 126.11334, 6, T=T)
    #_dTMP =BioMolecule('dHTMP-', 321.200521, 14, T=T, z=-1)
    _dTMPion =BioMolecule('dTMP2-', 320.192581, 13, T=T, z=-2)
    #--
    _dTMP = BioMolecule('d+H2TMP-(aq)', 322.208461, 15, T=T)

    _U = BioMolecule('Uracil(aq)', 112.086, 4, T=T)
    #_UMP = BioMolecule('HUMP-', 323.173341, 12, T=T, z=-1)
    _UMPion = BioMolecule('UMP2-', 322.165401, 11, T=T, z=-2)
    #--
    _UMP = BioMolecule('+H2UMP-(aq)', 324.181281, 13, T=T)


    # sugars
    _rib = BioMolecule('Ribose(aq)', 150.1299, 10, T=T)
    _drib = BioMolecule('Deoxyribose(aq)', 134.1305, 10, T=T)
    # estimate that the -OH group on a sugar is similar to the difference
    # between ribose and deoxyribose (one OH bond on the ring).
    dGf_OH = _rib.stdbio_formation_gibbs - _drib.stdbio_formation_gibbs

    # other useful molecules - remove comments as needed!
    _H2PO4 = BioMolecule('H2PO4-', 96.987241, 2, T=T, z=-1)
    _H3PO4 = BioMolecule('H3PO4(aq)', 97.995181, 3, T=T)
    # _HPO4 = BioMolecule('HPO4--', 95.979301, 1, T=T, z=-2)
    # _PO4 = BioMolecule('PO4---', 94.971361, 0, T=T, z=-3)
    # _ribphos = BioMolecule('Ribose-5-Phosphate-2', 150.1299, 10, T=T, z=-2)
    # _AS = BioMolecule('Adenosine(aq)', 267.24132, 13, T=T)
    # _CS = BioMolecule('Cytidine(aq)', 243.21662, 13, T=T)
    # _GS = BioMolecule('Guanosine(aq)', 283.24072, 13, T=T)
    # _US = BioMolecule('Uridine(aq)', 244.20138, 12, T=T)
    # _dAS = BioMolecule('Deoxyadenosine(aq)', 251.24192, 13, T=T)
    # _dCS = BioMolecule('Deoxycytidine(aq)', 227.21722, 13, T=T)
    # _dGS = BioMolecule('Deoxyguanosine(aq)', 267.24132, 13, T=T)
    # _dTS = BioMolecule('Deoxythymidine(aq)', 242.22856, 14, T=T)

    # Bonds (KJ/mol)
    #dGf_PhdB= -22.175 # Phosphodiester bond At 25 °C, pH 7 (Repair) Dickson K. 2000
    dGf_PhdB= 22.175 # (Reversed) Phosphodiester bond At 25 °C, pH 7 (Repair) Dickson K. 2000
    #dGf_PhdB= 103.76 # Phosphodiester bond At 65 °C, pH 8 (Hydrolysis) Molina R. 2015

    # NB for this we CAN'T change it with temperature - so this brings uncertainty

    """-F-"""
    #---Dictionaries----
    #The dictionaries for both DNA  and RNA contain (In the same order):
    #                       0= Nucleotide name
    #                       1= Base name
    #                       2= Molecular weight (Da)
    #                       3= Absolute intracellular concentration mol/L - E. coli
    #                       4= Absolute intracellular concentration mol/L - Mammalian
    #                       5= Absolute intracellular concentration mol/L - Yeast
    #                       6= dGf base - Cell conditions (Noor E. 2005 - eQuilibrator) (KJ mol)
    #                       7= dGf base - Environmental conditions (Noor E. 2005 - eQuilibrator) (KJ mol)
    #                       8= dGf - Nucleotide at cell conditions (Noor E. 2005 eQuilibrator) (KJ mol)
    #                       9= dGf - Nucleotide at environmental conditions (McCollom & Amend 2005) (KJ mol)
    #                        10= dGf - base (database)
    #                        11= dGf - Nucleotide (database)
    #                        12= dGf - Nucleotide ion

    #Dictionary for DNA
    DNAdict = {
      'A': ['dAMP','A',_dAMP.Mr,8.8e-6, 1.68e-5, 8.8e-6, 335.5, -41.3, -365.2, -864.9,
        _A.stdbio_formation_gibbs, _dAMP.stdbio_formation_gibbs,
        _dAMPion.stdbio_formation_gibbs],
      'C': ['dCMP','C',_dCMP.Mr,3.71e-5, 3.71e-5,3.71e-5,188.1,33.3,-692.8,-1212.9,
        _C.stdbio_formation_gibbs, _dCMP.stdbio_formation_gibbs,
        _dCMPion.stdbio_formation_gibbs],
      'T': ['dTMP','T', _dTMP.Mr,1.18e-5,1.18e-5, 1.18e-5, 8.4,-163.8,6.87e-3,-842.6,
        _T.stdbio_formation_gibbs, _dTMP.stdbio_formation_gibbs,
        _dTMPion.stdbio_formation_gibbs],
      'G': ['dGMP','G',_dGMP.Mr,5.1e-5,5.1e-5,5.1e-5,293.7,56.1,-596.4,-1115.8,
        _G.stdbio_formation_gibbs, _dGMP.stdbio_formation_gibbs,
        _dGMPion.stdbio_formation_gibbs]}

    #Dictionary for RNA
    RNAdict = {
      'A': ['AMP','A',_AMP.Mr,2.8e-4,4.23e-5,8.12e-5,335.5,-41.3,-529.3,-1020.5,
        _A.stdbio_formation_gibbs, _AMP.stdbio_formation_gibbs,
        _AMPion.stdbio_formation_gibbs],
      'C': ['CMP','C',_CMP.Mr,3.6e-4,1.18e-5,5.18e-6,188.1,33.3,-862.5,-1368.5,
        _C.stdbio_formation_gibbs, _CMP.stdbio_formation_gibbs,
        _CMPion.stdbio_formation_gibbs],
      'U': ['UMP','U',_UMP.Mr,1.45e-5,1.45e-5,1.45e-5,-78.5,-92.4,-1091.9,-1583.1,
        _U.stdbio_formation_gibbs, _UMP.stdbio_formation_gibbs,
        _UMPion.stdbio_formation_gibbs],
      'G': ['GMP','G',_GMP.Mr,2.4e-5,1.81e-5,1.02e-5,293.7,56.1,-764.8,-1271.4,
        _G.stdbio_formation_gibbs, _GMP.stdbio_formation_gibbs,
        _GMPion.stdbio_formation_gibbs]}

    """-G-"""
    #---DNA---
    """ G.0 - Initialise values"""
    #--- Number of nucleotides on each strand and total weight---
    genome_length=0
    totweight=0
    #---Energy of each strand---
    DNA_dGf=[]
    DNA_dGr=[]
    #---By block---
    DNA_block_dGf=[]
    DNA_block_dGr=[]
    #Molar Gibbs free energies for all the genes
    DNA_reactionGs = []
    DNA_block_reactionGs= []


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
        genome_length+= Nn

        """G.2 - Weight of the genome"""
        #The product of each nucleotide by its corresponding mass
        gene_weight = (
        (CA * DNAdict['A'][2]) + (CT * DNAdict['T'][2]) + (CC * DNAdict['C'][2]) + (CG * DNAdict['G'][2]) - ((Nn-1)*18))
        totweight+= gene_weight

        """G.3 - Energy of the nucleotides"""
        # At standard biological conditions, with pH 7 and ionic strength 0.

        #sigmaB will get the sum of energy of formation of the bases in the sequence
        sigmaB = (
        (CA * DNAdict['A'][10]) + (CT * DNAdict['T'][10]) + (CC * DNAdict['C'][10]) + (CG * DNAdict['G'][10]))

        #sigmaN will get the sum of energy of formation of the 1- charged  nucleotides in the sequence (e.g. 1 OH group on the phosphate)
        sigmaN = (
            (CA * DNAdict['A'][11]) + (CT * DNAdict['T'][11]) + (CC * DNAdict['C'][11]) + (CG * DNAdict['G'][11]))

        #sigmaNion will get the sum of energy of formation of the 2- charged  nucleotides in the sequence (e.g. 0 OH group on the phosphate, for group contribution)
        sigmaNion = (
        (CA * DNAdict['A'][12]) + (CT * DNAdict['T'][12]) + (CC * DNAdict['C'][12]) + (CG * DNAdict['G'][12]))

        #Phosphate backbone = sugar + phosphate
        phosba= (
        Nn * (_drib.stdbio_formation_gibbs + _H2PO4.stdbio_formation_gibbs))

        """G.4 - Formation and reaction energy of the genome"""

        """Gibbs Formation energy"""
        # sum of the nucleotide 2- ions, with n-1 Phosphodiester bonds made and an OH- removed from the sugar.
        DNA_dGf_chain = sigmaNion + (Nn-1)*(dGf_PhdB - dGf_OH)
        DNA_dGf.append(DNA_dGf_chain)
        """Gibbs Reaction Energy"""
        DNA_dGr_ss = ((((Nn-1)*dGf_H2O) + DNA_dGf_chain) - sigmaN)
        DNA_dGr.append(DNA_dGr_ss)


        """G.5 - Formation and reaction energy of the genome by blocks"""
        """Gibbs Formation energy"""
        DNA_block_dGf_chain = sigmaB + (phosba - dGf_OH)
        DNA_block_dGf.append(DNA_block_dGf_chain)
        """Gibbs Reaction Energy"""
        DNA_block_dGr_ss = ((Nn*dGf_H2O) + DNA_block_dGf_chain) - sigmaN
        DNA_block_dGr.append(DNA_block_dGr_ss)
        #DNA_block_totdGr += DNA_block_dGr


    """-H-"""
    #Average values for molecular weight and gene lenght
    #DNA oncentration. for the reaction quotient
    nucleotide_avmw= 327 #Da
    gene_meanlength= genome_length/len(genelist) #Average lenght of the strands
    avG= totweight/len(genelist) #average gene size in Daltons
    #host_volume= 1 #um3
    # To get the moles
    if celltype == "bacteria":  #15 mg/mL with avG (15mg/ml from bionumbers [Elowitz MB, Surette MG, Wolf PE, Stock JB, Leibler S. Protein mobility in the cytoplasm of Escherichia coli. J Bacteriol. 1999 Jan181(1):197-203])
        conc_Da_ml = 15/(1.66054e-21) # 1.66054e-21 is the conversion Da-mg
    elif celltype == "mammalian":
        conc_Da_ml = 3/(1.66054e-21)
    elif celltype == "yeast":
        conc_Da_ml = 0.36425/(1.66054e-21)

    conc_molecules_ml = conc_Da_ml / avG #number molecules per ml
    mol_ml = conc_molecules_ml/6.022e23
    DNA_mol_L = mol_ml*1000


    """-I-"""
    #Reaction quotient for this gene
    PC = DNA_mol_L

    for gene, stdG in zip(genelist, DNA_dGr):
        DNA_lnQ = math.log(PC/len(genelist))
        for base in gene:
            if base in DNAdict.keys():
                                if celltype == "bacteria":
                                    DNA_lnQ -= math.log(DNAdict[base][3])
                                elif celltype == "mammalian":
                                    DNA_lnQ -= math.log(DNAdict[base][4])
                                elif celltype == "yeast":
                                    DNA_lnQ -= math.log(DNAdict[base][5])

        DNA_reactionGs.append((stdG+(R*(T)*DNA_lnQ)))
    #--------By blocks-----------
    #Reaction quotient for the Gibbs energy by building block

    for geneblock, B_stdG in zip(genelist, DNA_block_dGr):
        DNA_B_lnQ = math.log(PC/len(genelist))
        for base in geneblock:
            if base in DNAdict.keys():
                                if celltype == "bacteria":
                                    DNA_B_lnQ -= math.log(DNAdict[base][3])
                                elif celltype == "mammalian":
                                    DNA_B_lnQ -= math.log(DNAdict[base][4])
                                elif celltype == "yeast":
                                    DNA_B_lnQ -= math.log(DNAdict[base][5])

        DNA_block_reactionGs.append((B_stdG+(R*(T)*DNA_B_lnQ)))


    """-J-"""
    #Molar Gibbs Energy
    #take the average, which is in kJ/mol
    DNA_avgReactionG_kJmol = statistics.mean(DNA_reactionGs)
    # convert to kJ / dry g
    #avG is the average gene size
    DNA_avgReactionG_kJdryg = DNA_avgReactionG_kJmol * 1/(avG)
    DNA_ReactionGibbs.append(DNA_avgReactionG_kJdryg)
    #Sum all the energies for the average


    #---By block---
    DNA_block_avgReactionG_kJmol = statistics.mean(DNA_block_reactionGs)
    DNA_block_avgReactionG_kJdryg = DNA_block_avgReactionG_kJmol * 1/(avG)
    DNA_block_ReactionGibbs.append(DNA_block_avgReactionG_kJdryg)
    #Sum all the energies for the average



    """-K-"""
    #---RNA---
    """ K.0 - Initialise values"""
    #--- Number of nucleotides on each strand and total weight---
    transcriptome_lenght=0
    transcriptome_totweight=0
    #---Energy of each mRNA---
    RNA_dGf=[]
    RNA_dGr=[]
    #---By block---
    RNA_block_dGf=[]
    RNA_block_dGr=[]
    #Molar Gibbs free energies for all the genes
    RNA_reactionGs = []
    RNA_block_reactionGs= []

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
            (RA * RNAdict['A'][10]) + (RU * RNAdict['U'][10]) + (RC * RNAdict['C'][10]) + (RG * RNAdict['G'][10]))

            #RNA_sigmaN will get the energy of formation of the 1- Charged nucleotides in the sequence at cellular conditions
            RNA_sigmaN = (
            (RA * RNAdict['A'][11]) + (RU * RNAdict['U'][11]) + (RC * RNAdict['C'][11]) + (RG * RNAdict['G'][11]))

            #RNA_sigmaNion will get the sum of energy of formation of the 2- charged  nucleotides in the sequence (e.g. 0 OH group on the phosphate, for group contribution)
            RNA_sigmaNion = (
            (RA * RNAdict['A'][12]) + (RU * RNAdict['U'][12]) + (RC * RNAdict['C'][12]) + (RG * RNAdict['G'][12]))

            RNA_phosba= (
            RNA_Nn * (_rib.stdbio_formation_gibbs + _H2PO4.stdbio_formation_gibbs))


            """K.4 - Formation and reaction energy of the transcriptome"""
            #This method will calculate the energy based on the dGf of the phosphate
            #backbone and the dGf of the base without considering intramolecular interactions

            """Gibbs Formation energy"""
            # sum of the nucleotide 2- ions, with n-1 Phosphodiester bonds made and an OH- removed.
            RNA_dGf_chain = RNA_sigmaNion + (RNA_Nn-1)*(dGf_PhdB - dGf_OH)
            RNA_dGf.append(RNA_dGf_chain)

            """Gibbs reaction energy"""
            RNA_dGr_chain= ((((RNA_Nn-1)*dGf_H2O) + RNA_dGf_chain) - RNA_sigmaN)
            RNA_dGr.append(RNA_dGr_chain)

            """G.5 - Formation and reaction energy of the transcriptome by blocks"""

            """Gibbs Formation energy"""
            RNA_dGf_byblock = RNA_sigmaB + (RNA_phosba - dGf_OH)
            RNA_block_dGf.append(RNA_dGf_byblock)

            """Gibbs reaction energy"""
            RNA_dGr_byblock= ((Nn*dGf_H2O) + RNA_dGf_byblock) - RNA_sigmaN
            RNA_block_dGr.append(RNA_dGr_byblock)


    """-L-"""
    #get average values for molecular weight and lenght
    RNA_nucleotide_avmw= 339.5 #Da
    messenger_meanlenght= transcriptome_lenght/len(transcripts) #Average lenght of the strands
    avm= transcriptome_totweight/len(transcripts) #average gene weight in Daltons
    # same technique for concentrations as DNA
    host_volume= 1 #um3
    # To get the moles
    if celltype == "bacteria":
        RNA_conc_Da_ml = 100/(1.66054e-21) # 1.66054e-21 is the conversion Da-mg
    elif celltype == "mammalian":
        RNA_conc_Da_ml = 10/(1.66054e-21)
    elif celltype == "yeast":
        RNA_conc_Da_ml = 13/(1.66054e-21)

    RNA_conc_molecules_ml = RNA_conc_Da_ml / avm #number molecules per ml
    RNA_mol_ml = RNA_conc_molecules_ml/6.022e23
    RNA_mol_L = RNA_mol_ml*1000


    """-M-"""
    #Reaction quotient for the entire transcriptome
    RNA_PC = RNA_mol_L

    #Bucket method at cellular conditions
    for transcript, RNA_stdG in zip(transcripts, RNA_dGr):
        RNA_lnQ = math.log(RNA_PC/len(transcripts))
        for base in transcript:
            if base in RNAdict.keys():
                                if celltype == "bacteria":
                                    RNA_lnQ -= math.log(RNAdict[base][3])
                                elif celltype == "mammalian":
                                    RNA_lnQ -= math.log(RNAdict[base][4])
                                elif celltype == "yeast":
                                    RNA_lnQ -= math.log(RNAdict[base][5])

        RNA_reactionGs.append((RNA_stdG+(R*(T)*RNA_lnQ)))

    #---By block---
    for b_transcript, RNA_B_stdG in zip(transcripts, RNA_block_dGr):
        RNA_B_lnQ = math.log(RNA_PC/len(transcripts))
        for base in b_transcript:
            if base in RNAdict.keys():
                                if celltype == "bacteria":
                                    RNA_B_lnQ -= math.log(RNAdict[base][3])
                                elif celltype == "mammalian":
                                    RNA_B_lnQ -= math.log(RNAdict[base][4])
                                elif celltype == "yeast":
                                    RNA_B_lnQ -= math.log(RNAdict[base][5])

        RNA_block_reactionGs.append((RNA_B_stdG+(R*(T)*RNA_B_lnQ)))


    """-N-"""
    #Molar Gibbs Energy
    #take the average, which is in kJ/mol

    RNA_avgReactionG_kJmol = statistics.mean(RNA_reactionGs)
    # convert to kJ / dry g
    #avG is the average gene size
    #Bucket
    RNA_avgReactionG_kJdryg = RNA_avgReactionG_kJmol * 1/(avm)
    RNA_ReactionGibbs.append(RNA_avgReactionG_kJdryg)

    #---By block---
    RNA_avgReactionG_block_kJmol = statistics.mean(RNA_block_reactionGs)
    RNA_avgReactionG_block_kJdryg = RNA_avgReactionG_block_kJmol * 1/(avm)
    RNA_block_ReactionGibbs.append(RNA_avgReactionG_block_kJdryg)
    print(' ')


    """Averages"""
    #DNA energy
    DNA_dGf_av.append(statistics.mean(DNA_dGf)) #Formation Gibbs Energy
    DNA_dGr_av.append(statistics.mean(DNA_dGr)) #Reaction Gibbs Energy
    #DNA energy by block
    DNA_block_dGf_av.append(statistics.mean(DNA_block_dGf)) #Formation Gibbs Energy
    DNA_block_dGr_av.append(statistics.mean(DNA_block_dGr)) #Reaction Gibbs Energy
    #RNA Energy
    RNA_dGf_av.append(statistics.mean(RNA_dGf)) #Formation Gibbs Energy
    RNA_dGr_av.append(statistics.mean(RNA_dGr)) #Reaction Gibbs Energy
    #RNA energy by block
    RNA_block_dGf_av.append(statistics.mean(RNA_block_dGf)) #Formation Gibbs Energy
    RNA_block_dGr_av.append(statistics.mean(RNA_block_dGr)) #Reaction Gibbs Energy




"""-O-"""
#---Printing values---
if temp == 'yes':
    """DNA"""
    print('The energy required for this DNA sequence and the transcripts at',stemperature, 'K is:')
    print(' ')
    print('Genome')
    print('Average formation Gibbs energy of the genes:', "%.3f" % DNA_dGf_av[selected_temp], 'KJ/mol')
    print('Average reaction Gibbs energy of the genes:', "%.2f" % DNA_dGr_av[selected_temp], 'KJ/mol')
    print('Average molar energy of the genome:', "%.2f" % DNA_ReactionGibbs[selected_temp], 'KJ/g')
    print('Average molar energy of the genome by block:', "%.2f" % DNA_block_ReactionGibbs[selected_temp], 'KJ/g')
    print(' ')
    """RNA"""
    print('Transcriptome')
    print('Average formation Gibbs energy of the transcripts:', "%.2f" % RNA_dGf_av[selected_temp], 'KJ/mol')
    print('Average reaction Gibbs energy of the transcripts:', "%.2f" % RNA_dGr_av[selected_temp], 'KJ/mol')
    print('Molar energy of the transcriptome:', "%.2f" % RNA_ReactionGibbs[selected_temp], 'KJ/g')
    print('Molar energy of the transcriptome by block:', "%.2f" % RNA_block_ReactionGibbs[selected_temp], 'KJ/g')

"""-P-"""
#---Ploting values---
#--- Both ----
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy for a bacteria cell with two methods')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy for a mammalian cell with two methods')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy for a yeast cell with two methods')

ax.plot(TK,DNA_ReactionGibbs, label='Genome Energy', c='r', linewidth=2)
ax.plot(TK,DNA_block_ReactionGibbs, label='Genome energy by block', linestyle='dotted', c='b', linewidth=2)
#---RNA
ax.plot(TK,RNA_ReactionGibbs, label='Transcriptome Energy', linestyle='solid', c='b', linewidth=2)
ax.plot(TK,RNA_block_ReactionGibbs, label='Transcriptome Energy by block', linestyle='dotted', c='g', linewidth=2)



ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

#ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()

#---Energy ----
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy of the transcriptome for a bacteria cell')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy of the transcriptome for a mammalian cell')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy of the transcriptome for a yeast cell')

ax.plot(TK,DNA_ReactionGibbs, label='Genome Energy', c='r', linewidth=2)
ax.plot(TK,RNA_ReactionGibbs, label='Transcriptome Energy', linestyle='dashed', c='b', linewidth=2)
#ax.plot(TK,RNA_block_ReactionGibbs, label='RNA Energy by block', linestyle='dotted', c='g', linewidth=2)
#---DNA
#ax.plot(TK,DNA_block_ReactionGibbs, label='DNA Energy by block', linestyle='dotted', c='b', linewidth=2)

ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

#ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()

#---By block ----
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy of the transcriptome for a bacteria cell by block')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy of the transcriptome for a mammalian cell by block')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy of the transcriptome for a yeast cell by block')

#ax.plot(TK,RNA_ReactionGibbs, label='Energy to synthesise this transcriptome', linestyle='dashed', c='b', linewidth=2)
ax.plot(TK,RNA_block_ReactionGibbs, label='RNA Energy by block', linestyle='dotted', c='g', linewidth=2)
# ax.plot(TK, RNA_ReactionGibbs_env, label='Transcriptime energy at environment conditions', marker='o', c='b', linewidth=3)
#---DNA
ax.plot(TK,DNA_block_ReactionGibbs, label='DNA Energy by block', linestyle='dotted', c='b', linewidth=2)

ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

#ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()

print(' ')
if temp == 'yes':
    byblockresults=input('Do you want to print the results for the by block calculaitons? [y/n]')
    if byblockresults == 'y':
        """DNA"""
        print(' ')
        print('The energy required for this DNA sequence and the transcripts at',stemperature, 'K by block is:')
        print(' ')
        print('Genome')
        print('Average formation Gibbs energy of the genes:', "%.3f" % DNA_block_dGf_av[selected_temp], 'KJ/mol')
        print('Average reaction Gibbs energy of the genes:', "%.2f" % DNA_block_dGr_av[selected_temp], 'KJ/mol')
        print('Average molar energy of the genome by block:', "%.2f" % DNA_block_ReactionGibbs[selected_temp], 'KJ/g')
        print(' ')
        """RNA"""
        print('Transcriptome')
        print('Average formation Gibbs energy of the transcripts:', "%.2f" % RNA_block_dGf_av[selected_temp], 'KJ/mol')
        print('Average reaction Gibbs energy of the transcripts:', "%.2f" % RNA_block_dGr_av[selected_temp], 'KJ/mol')
        print('Molar energy of the transcriptome by block:', "%.2f" % RNA_block_ReactionGibbs[selected_temp], 'KJ/g')
        print(' ')



        # #------------------------G.4.2 Environmental conditions-------------------
        # #The environmental conditions are considered here as the standard conditions
        # #without accounting for pH, ionic strength or any other cellular factors
        # #(aqueous solution abd with a net charge of -2/-3 per molecule)
        #
        # #sigmaB_env will get the energy of formation of the bases in the sequence at environental conditions
        # RNA_sigmaB_env = (
        # (RA * RNAdict['A'][7]) + (RU * RNAdict['U'][7]) + (RC * RNAdict['C'][7]) + (RG * RNAdict['G'][7]))
        #
        # #sigmaN_env will get the energy of formation of the nucleotides in the sequence at environental conditions
        # RNA_sigmaN_env = (
        # (RA * RNAdict['A'][9]) + (RU * RNAdict['U'][9]) + (RC * RNAdict['C'][9]) + (RG * RNAdict['G'][9]))
        #
        # """----------Energy calculation at environmental conditions--------------"""
        # #This method will calculate the energy based on the dGf of the phosphate
        # #backbone and the dGf of the base without considering interactions
        #
        # """Gibbs Formation energy"""
        # # dGf using the equation for the the bucket method at environmental conditions
        # RNA_dGf_chain_env= (Nn*RNA_dGf_PhBB_env) + sigmaB_env
        # RNA_totdGf_env+= RNA_dGf_chain_env
        #
        # """Gibbs Reaction Energy"""
        # RNA_dGr_env= ((((Nn-1)*dGf_H2O) + RNA_dGf_chain_env) - sigmaN_env)
        # RNA_geneGs_env.append(RNA_dGr_env)
        # RNA_totdGr_env+= RNA_dGr_env
        #
        # """Reation quotient at environental conditions"""
        # PC = RNA_mol_L
        # RNA_ReactionGibbs_env = []
        # TK = np.linspace(253.15,400)
        #
        # for T in TK:
        #     #make a list of the free energies for each gene
        #     RNA_reactionGs_env = []
        #     #Bucket method at environmental conditions
        #     for gene, RNA_stdG in zip(genelist, RNA_geneGs_env):
        #         RNA_lnQ_env = math.log(PC/len(genelist))
        #         for base in gene:
        #             if base in RNAdict.keys():
        #                                 if celltype == "bacteria":
        #                                     RNA_lnQ_env -= math.log(RNAdict[base][3])
        #                                 elif celltype == "mammalian":
        #                                     RNA_lnQ_env -= math.log(RNAdict[base][4])
        #                                 elif celltype == "yeast":
        #                                     RNA_lnQ_env -= math.log(RNAdict[base][5])
        #         RNA_reactionGs_env.append((RNA_stdG+(R*(T)*RNA_lnQ_env)))
        #
        #
        #     #take the average, which is in kJ/mol
        #     RNA_avgReactionG_kJmol_env = statistics.mean(RNA_reactionGs_env)
        #     # convert to kJ / dry g
        #     #avG is the average gene size
        #     #Bucket
        #     RNA_avgReactionG_kJdryg_env = RNA_avgReactionG_kJmol_env * 1/(avG)
        #     RNA_ReactionGibbs_env.append(RNA_avgReactionG_kJdryg_env)
