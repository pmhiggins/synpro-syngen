import sys,os, math, statistics, re
import numpy as np
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.Seq import transcribe
import matplotlib.pyplot as plt

from BioMolecule import BioMolecule

#----------For Extracting the transcripts from the RNA sequence input file-------
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
            foundStart = False # this tells it to start looking again
            mRNAlst.append(mRNA) # save what we have to the list
            mRNA='' # reset for starting again

    # if we got to the end with no TAG, add everything saved to mRNA so far
    if mRNA != '':
        mRNAlst.append(mRNA)
    return(mRNAlst)


R= 0.008314472 #Gas constant KJ mol K
RNA_ReactionGibbs = []
DNA_ReactionGibbs = []
TK = np.linspace(275,400, num=6)


#------------B. Obtaining the sequences and definiting the conditions-----------
print('Please type the name of a RNA fasta file, the file should be contained in the same folder')
fname= input('fasta file name: ')
print('---------------------------------Environments----------------------------')
print('This program can calculate the energy to build a genome and the respective transcriptome for three model organisms in two conditions:')
print('-->Environmental conditions: Here denominated as the standard conditions without accounting for pH, ionic strength or any other cellular factors (aqueous solution abd with a net charge of -2/-3 per molecule)')
print('-->Cell conditions: Standard conditions plus a pH of 7 and a default ionic strength of 0.1 M. (in an aqueous solution and with a net charge of 0)')
print('---------------------------------Model cell-----------------------------')
celltype= input('Do you want to calculate the Gibbs energy for a bactetria, yeast or mammalian cell? [bacteria/yeast/mammalian]')

#-----------------------------C. Opening with biopython-------------------------
Seqname = []
genelist = []
records= []
mRNAs=[]

fasta = open(fname,'r')
for record in SeqIO.parse(fasta,'fasta'):
    Seqname.append(record.id) #Name of the gene
    genelist.append(str(record.seq)) #Genes in a list of strings
    records.append((record.seq)) # Genes in a list of int to be used by biopyhton

#for record, name in zip(records, Seqname):
#    print(record, name)

#-------------------------D. Obtaining the complementary strand-----------------

for record in records:
    complementary=record.complement()
    genelist.append(str(complementary))

print('there are', len(genelist)/2, 'sequences in this file, threfore, the calculations will be done for',len(genelist),'strands to get the values of the double strand DNA')

for gene in records:
    mRNA = transcription(gene)
    mRNAs.extend(mRNA) #  Add a list on the end



for T in TK:
    #-----------------------------A. Initial values---------------------------------

    # Gibbs energy of formation
    #Inorganic compounds (KJ/mol) #McCollom & Amend 2005
    dGf_H=0
    dGf_e=0

    _Water = BioMolecule('H2O(l)', 18, 2, T=T)

    dGf_H2O = _Water.stdbio_formation_gibbs

    # bases and nucleotides
    # BioMolecules are declared with (name, Mr, no. H, T, charge)
    _A = BioMolecule('Adenine(aq)', 135, 5, T=T)
    _AMP = BioMolecule('HAMP-', 346.205341, 13, T=T, z=-1)
    _AMPion = BioMolecule('AMP2-', 345.205341, 12, T=T, z=-2)
    _dAMP = BioMolecule('dHAMP-', 330.213881, 13, T=T, z=-1)
    _dAMPion = BioMolecule('dAMP2-', 329.205941, 12, T=T, z=-2)

    _C = BioMolecule('Cytosine(aq)', 111.102, 5, T=T)
    _CMP = BioMolecule('HCMP-', 322.188581, 13, T=T, z=-1)
    _CMPion = BioMolecule('CMP2-', 321.180641, 12, T=T, z=-2)
    _dCMP = BioMolecule('dHCMP-', 306.189181, 13, T=T, z=-1)
    _dCMPion = BioMolecule('dCMP2-', 305.181241, 12, T=T, z=-2)

    _G = BioMolecule('Guanine(aq)', 151.1261, 5, T=T)
    _GMP = BioMolecule('HGMP-', 362.212681, 13, T=T, z=-1)
    _GMPion = BioMolecule('GMP2-', 361.204741, 12, T=T, z=-2)
    _dGMP = BioMolecule('dHGMP-', 346.213281, 13, T=T, z=-1)
    _dGMPion = BioMolecule('dGMP2-', 345.205341, 12, T=T, z=-2)

    _T = BioMolecule('Thymine(aq)', 126.11334, 6, T=T)
    _dTMP =BioMolecule('dHTMP-', 321.200521, 14, T=T, z=-1)
    _dTMPion =BioMolecule('dTMP2-', 320.192581, 13, T=T, z=-2)

    _U = BioMolecule('Uracil(aq)', 112.086, 4, T=T)
    _UMP = BioMolecule('HUMP-', 323.173341, 12, T=T, z=-1)
    _UMPion = BioMolecule('UMP2-', 322.165401, 11, T=T, z=-2)


    # sugars
    _rib = BioMolecule('Ribose(aq)', 150.1299, 10, T=T)
    _drib = BioMolecule('Deoxyribose(aq)', 134.1305, 10, T=T)
    # estimate that the -OH group on a sugar is similar to the difference
    # between ribose and deoxyribose (one OH bond on the ring).
    dGf_OH = _rib.stdbio_formation_gibbs - _drib.stdbio_formation_gibbs


    # other useful molecules - remove comments as needed!

    # _H2PO4 = BioMolecule('H2PO4-, 96.987241, 2, T=T, z=-1)
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
    dGf_PhdB = 22.175 # Phosphodiester bond (hidrolysis) Dickson K. 2000

    # NB for this we CAN'T change it with temperature - so this brings uncertainty

    """Dictionaries"""
    #-----------------------------E. Dictionaries-----------------------------------
    #The dictionaries for both DNA  and RNA contain (In the same order):
    #                       0= monophosphate name
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
    #Dictionary for DNA
    DNAdict = {
      'A': ['dAMP','A',_dAMP.Mr,8.8e-6, 1.68e-5, 8.8e-6,335.5,-41.3,-365.2,-864.9,
        _A.stdbio_formation_gibbs, _dAMP.stdbio_formation_gibbs,
        _dAMPion.stdbio_formation_gibbs],
      'C': ['dCMP','C',_dCMP.Mr,8.02e-3, 3.71e-5,8.02e-3,188.1,33.3,-692.8,-1212.9,
        _C.stdbio_formation_gibbs, _dCMP.stdbio_formation_gibbs,
        _dCMPion.stdbio_formation_gibbs],
      'T': ['dTMP','T', _dTMP.Mr,6.87e-3, 1.18e-5,-1434.2,8.4,-163.8,6.87e-3,-842.6,
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
      'U': ['UMP','U',_UMP.Mr,6.6e-4,1.45e-5,1.45e-5,-78.5,-92.4,-1091.9,-1583.1,
        _U.stdbio_formation_gibbs, _UMP.stdbio_formation_gibbs,
        _UMPion.stdbio_formation_gibbs],
      'G': ['GMP','G',_GMP.Mr,2.4e-5,1.81e-5,1.02e-5,293.7,56.1,-764.8,-1271.4,
        _G.stdbio_formation_gibbs, _GMP.stdbio_formation_gibbs,
        _GMPion.stdbio_formation_gibbs]}

    """----------------------------G. DNA----------------------------------------"""
    #--------------------------G.1 Initialise values--------------------------------
    genome_lenght=0
    totweight=0
    #---------------------For the Energy of each sequence---------------------------

    DNA_totdGf=0
    DNA_geneGs=[]
    DNA_totdGr= 0

    #make a list of the free energies for each gene
    DNA_reactionGs = []

    #------------------------G.2 Iterate through the DNA sequences------------------
    for sequence in genelist:
        #Total number of nucleotides
        Nn= len(sequence)
        genome_lenght+= Nn
        #Count for each nucleotide
        CA=sequence.count('A')
        CT=sequence.count('T')
        CC=sequence.count('C')
        CG=sequence.count('G')

        #-------------------------G.3 Weight of the genome--------------------------
        #To get the product of each nucleotides weight by the number of nucleotides
        gene_weight = (
        (CA * DNAdict['A'][2]) + (CT * DNAdict['T'][2]) + (CC * DNAdict['C'][2]) + (CG * DNAdict['G'][2]) - ((Nn-1)*18))
        totweight+= gene_weight

        #---------------------G.4 Miscellaneous values of the genome--------------------


        #--------G.5 Energy of the genome -------

        # we use standard biological conditions, with pH 7 and ionic strength 0.

        #sigmaB will get the sum of energy of formation of the bases in the sequence
        sigmaB = (
        (CA * DNAdict['A'][10]) + (CT * DNAdict['T'][10]) + (CC * DNAdict['C'][10]) + (CG * DNAdict['G'][10]))

        #sigmaN will get the sum of energy of formation of the 1- charged  nucleotides in the sequence (e.g. 1 OH group on the phosphate)
        sigmaN = (
            (CA * DNAdict['A'][11]) + (CT * DNAdict['T'][11]) + (CC * DNAdict['C'][11]) + (CG * DNAdict['G'][11]))

        #sigmaNion will get the sum of energy of formation of the 2- charged  nucleotides in the sequence (e.g. 0 OH group on the phosphate, for group contribution)
        sigmaNion = (
        (CA * DNAdict['A'][12]) + (CT * DNAdict['T'][12]) + (CC * DNAdict['C'][12]) + (CG * DNAdict['G'][12]))

        """---------- Energy calculation from nucleotides --------------"""
        #The bucket method will calculate the energy based on the dGf of the phosphate
        #backbone and the dGf of the base without considering interactions

        """Gibbs Formation energy"""
        #GCA for dGf of the chain
        # sum of the nucleotide 2- ions, with n-1 Phosphodiester bonds made and an OH- removed from the sugar.
        DNA_dGf_chain = sigmaNion + (Nn-1)*(dGf_PhdB - dGf_OH)

        DNA_totdGf += DNA_dGf_chain

        """Gibbs Reaction Energy"""
        DNA_dGr = ((((Nn-1)*dGf_H2O) + DNA_dGf_chain) - sigmaN)

        DNA_geneGs.append(DNA_dGr)
        DNA_totdGr += DNA_dGr



    """Concentrations"""
    #get average values for molecular weight and gene lenght
    nucleotide_avmw= 327 #Da
    gene_meanlenght= genome_lenght/len(genelist) #Average lenght of the strands
    avG= totweight/len(genelist) #average gene size in Daltons

    host_volume= 1 #um3
    # To get the moles
    # use 15 mg/mL with avG (15mg/ml from bionumbers)
    conc_Da_ml = 15/(1.66054e-21) # 1.66054e-21 is the conversion Da-mg
    conc_molecules_ml = conc_Da_ml / avG #number molecules per ml
    mol_ml = conc_molecules_ml/6.022e23
    DNA_mol_L = mol_ml*1000

    ### not sure how this worked
    # moles= genome_lenght/6.022e23
    # To get mol per litre
    # print(genome_lenght/(6.022e23 * host_volume * 1e-15))

    """Reaction quotient for this gene"""
    PC = DNA_mol_L

    for gene, stdG in zip(genelist, DNA_geneGs):
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




    #take the average, which is in kJ/mol
    DNA_avgReactionG_kJmol = statistics.mean(DNA_reactionGs)
    # convert to kJ / dry g
    #avG is the average gene size
    #Bucket
    DNA_avgReactionG_kJdryg = DNA_avgReactionG_kJmol * 1/(avG)
    DNA_ReactionGibbs.append(DNA_avgReactionG_kJdryg)


    #------------------------G.4.2 Environmental conditions-------------------



    """------------------------------J. RNA--------------------------------------"""
    #-------------------------J.1 Initialise values---------------------------------
    transcriptome_lenght=0
    transcriptome_totweight=0
    #---------------------For the Energy of each sequence---------------------------

    RNA_totdGf=0
    RNA_geneGs=[]
    RNA_totdGr= 0

    #-------------------------J.2 Iterate through the RNA sequences-----------------
    for m in mRNAs:
        #Total number of RNA nucleotides
        if len(m) != 0: #To get rid of blank lists
            RNA_Nn= len(m)
            transcriptome_lenght+= RNA_Nn

            #Number of each nucleotide
            RA=m.count('A')
            RU=m.count('U')
            RC=m.count('C')
            RG=m.count('G')

            #-----------------------J.3 Weight of the transcriptome-----------------
            #To get the product of each nucleotides weight by the number of nucleotides
            RNA_weight = (
            (RA * RNAdict['A'][2]) + (RU * RNAdict['U'][2]) + (RC * RNAdict['C'][2]) + (RG * RNAdict['G'][2]) - ((RNA_Nn-1)*18))
            """We take away the H2O weight for every nucleotide"""
            transcriptome_totweight+= RNA_weight

            #--------------J.5 Miscellaneous values for the transcriptome-------------------


            #--------J.5 Energy of the proteome -------

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

            """---------- Energy calculation --------------"""
            #The bucket method will calculate the energy based on the dGf of the phosphate
            #backbone and the dGf of the base without considering interactions

            """Gibbs Formation energy at cellular conditions"""
            # sum of the nucleotide 2- ions, with n-1 Phosphodiester bonds made and an OH- removed.
            RNA_dGf_chain = RNA_sigmaNion + (Nn-1)*(dGf_PhdB - dGf_OH)
            RNA_totdGf+= RNA_dGf_chain

            """Gibbs Reaction Energy at ceullar conditions"""
            RNA_dGr= ((((Nn-1)*dGf_H2O) + RNA_dGf_chain) - RNA_sigmaN)
            RNA_geneGs.append(RNA_dGr)
            RNA_totdGr+= RNA_dGr

    """Concentrations"""
    #get average values for molecular weight and lenght
    RNA_nucleotide_avmw= 339.5 #Da
    messenger_meanlenght= transcriptome_lenght/len(mRNAs) #Average lenght of the strands
    avm= transcriptome_totweight/len(mRNAs) #average gene size in Daltons

    # same technique for concentrations as DNA
    host_volume= 1 #um3
    # To get the moles
    # use 100 mg/mL with avG (100mg/ml from bionumbers)
    RNA_conc_Da_ml = 100/(1.66054e-21) # 1.66054e-21 is the conversion Da-mg
    RNA_conc_molecules_ml = RNA_conc_Da_ml / avm #number molecules per ml
    RNA_mol_ml = RNA_conc_molecules_ml/6.022e23
    RNA_mol_L = RNA_mol_ml*1000




    """Reation quotient at cellular conditions"""
    RNA_PC = RNA_mol_L

    #make a list of the free energies for each gene
    RNA_reactionGs = []
    #Bucket method at cellular conditions
    for gene, RNA_stdG in zip(mRNAs, RNA_geneGs):
        RNA_lnQ = math.log(RNA_PC/len(genelist))
        for base in gene:
            if base in RNAdict.keys():
                                if celltype == "bacteria":
                                    RNA_lnQ -= math.log(RNAdict[base][3])
                                elif celltype == "mammalian":
                                    RNA_lnQ -= math.log(RNAdict[base][4])
                                elif celltype == "yeast":
                                    RNA_lnQ -= math.log(RNAdict[base][5])

        RNA_reactionGs.append((RNA_stdG+(R*(T)*RNA_lnQ)))


    #take the average, which is in kJ/mol
    RNA_avgReactionG_kJmol = statistics.mean(RNA_reactionGs)
    # convert to kJ / dry g
    #avG is the average gene size
    #Bucket
    RNA_avgReactionG_kJdryg = RNA_avgReactionG_kJmol * 1/(avm)
    RNA_ReactionGibbs.append(RNA_avgReactionG_kJdryg)








#------------------------I. Printing values of the genome-----------------------
print('-------------------------------------------Genome information-------------------------------------------------------------------')
# if celltype == "bacteria":
#     print('These are the values of the genome in the bacterial cell at cellular and environental conditions')
#
# elif celltype == "mammalian":
#     print('These are the values of the genome in the mammalian cell at cellular and environental conditions')
#
# elif celltype == "yeast":
#     print('These are the values of the genome in the yeast cell at cellular and environental conditions')

# print('--> Energy at Cellular conditions:')
# print('Average Gibbs Free energy of Formation [dGf]:', DNA_totdGf/(len(genelist)), 'KJ per sequence')
# print('Average Gibbs Free energy of reaction [dGr]:', DNA_totdGr/(len(genelist)), 'KJ per sequence')
# print('The energy to build the genome adjusted to 1 gram at 298.15K:')
# print(DNA_ReactionGibbs[45],'KJ/g considering the absolute intracellular metabolite concentrations for the selected cell type and the sequence lenght')


# print('--> Energy at Environmental conditions:')
# print('Average Gibbs Free energy of Formation [dGf]:', DNA_totdGf_env/(len(genelist)), 'KJ per sequence')
# print('Average Gibbs Free energy of reaction [dGr]:', DNA_totdGr_env/(len(genelist)), 'KJ per sequence')
# print('The energy to build the genome adjusted to 1 gram at 298.15K:')
# print(DNA_ReactionGibbs_env[45],'KJ/g considering the absolute intracellular metabolite concentrations for the selected cell type and the sequence lenght')


print('Average values per gene sequence (Wheight and lenght):')
print('The weight of the double strand genome, based on the input file is:', totweight,'Da')
print('The average weight of the genes or genome of the input file is' ,avG,'Da with', gene_meanlenght,'nucleotides in average')
print('The concentration is',DNA_mol_L,' moles of DNA per litre')



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

#------------------------I. Printing values of the transcriptome-----------------------
print('-------------------------------------------Transcriptome information-------------------------------------------------------------------')
# if celltype == "bacteria":
#     print('These are the values of the transcriptome in the bacterial cell at cellular and environental conditions')
#
# elif celltype == "mammalian":
#     print('These are the values of the transcriptome in the mammalian cell at cellular and environental conditions')
#
# elif celltype == "yeast":
#     print('These are the values of the transcriptome in the yeast cell at cellular and environental conditions')
#
#
# print('--> Energy at Cellular conditions:')
# print('Average Gibbs Free energy of Formation [dGf]:', RNA_totdGf/(len(genelist)), 'KJ per sequence')
# print('Average Gibbs Free energy of reaction [dGr]:', RNA_totdGr/(len(genelist)), 'KJ per sequence')
# print('The energy to build the transcriptome adjusted to 1 gram at 298.15K:')
# print(RNA_ReactionGibbs[45],'KJ/g considering the absolute intracellular metabolite concentrations for the selected cell type and the sequence lenght')
#
# print('--> Energy at Environmental conditions:')
# print('Average Gibbs Free energy of Formation [dGf]:', RNA_totdGf_env/(len(genelist)), 'KJ per sequence')
# print('Average Gibbs Free energy of reaction [dGr]:', RNA_totdGr_env/(len(genelist)), 'KJ per sequence')
# print('The energy to build the transcriptome adjusted to 1 gram at 298.15K:')
# print(RNA_ReactionGibbs_env[45],'KJ/g considering the absolute intracellular metabolite concentrations for the selected cell type and the sequence lenght')

print('Average values per gene sequence (Wheight and lenght):')
print('The weight of the mRNA, based on the input file is:', totweight,'Da')
print('The average weight of the mRNAs transcribed from DNA the input file is' ,avG,'Da with', gene_meanlenght,'nucleotides in average')
print('The concentration is',RNA_mol_L,' moles of RNA per litre')


#----------------------Plotting for the DNA ------------------------------------
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy of the genome calculated in two conditions for a bacteria cell')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy of the genome calculated in two conditions for a mammalian cell')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy of the genome calculated in two conditions for a yeast cell')

ax.plot(TK,DNA_ReactionGibbs, label='Genome energy at cell conditions', linestyle='dashed', c='r', linewidth=3)
# ax.plot(TK, DNA_ReactionGibbs_env, label='Genome energy at environment conditions', marker='o', c='b', linewidth=3)

ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()

#----------------------Plotting for the RNA ------------------------------------
fig = plt.figure(figsize = (7,5))
ax= fig.add_subplot(111)

if celltype == "bacteria":
    plt.title('Molar Gibbs Energy of the transcriptome calculated in two conditions for a bacteria cell')
elif celltype == "mammalian":
    plt.title('Molar Gibbs Energy of the transcriptome calculated in two conditions for a mammalian cell')
elif celltype == "yeast":
    plt.title('Molar Gibbs Energy of the transcriptome calculated in two conditions for a yeast cell')

ax.plot(TK,RNA_ReactionGibbs, label='Transcriptome energy at cell conditions', linestyle='dashed', c='r', linewidth=3)
# ax.plot(TK, RNA_ReactionGibbs_env, label='Transcriptime energy at environment conditions', marker='o', c='b', linewidth=3)

ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()
