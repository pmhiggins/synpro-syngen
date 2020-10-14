import sys,os, math, statistics, re
import numpy as np
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.Seq import transcribe
import matplotlib.pyplot as plt

#-----------------------------A. Initial values---------------------------------
# Gibbs energy of formation
#Inorganic compounds (KJ/mol) #McCollom & Amend 2005
dGf_H2O= -237.3
dGf_H=0
dGf_e=0
dGf_HCO3= -587.2
dGf_NH4= -79.5
dGf_H2= 17.7
dGf_H2S= -27.9
dGf_HPO4= -1089.7
dGf_NO3= -111.0
dGf_NO2= -32.2
dGf_SO4_2= -744.8

#Phosphate backbones
"""environmental conditions"""
DNA_dGf_PhBB_env = -1293.70 #(KJ mol)
RNA_dGf_PhBB_env = -1461.50 #(KJ mol)

"""cell conditions"""
DNA_dGf_PhBB_cell = -966.70 #(KJ mol)
RNA_dGf_PhBB_cell = -1122.50 #(KJ mol)

#Bonds (KJ/mol)
dGf_PhdB =-22.175 #Phosphodiester bond (hidrolysis) Dickson K. 2000

#Standard values
R= 0.008314472 #Gas constant KJ mol K
T= 298.15 #Standard temperature in Kelvin

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

"""Functions"""
#----------E. Extracting the transcripts from the RNA sequence input file-------
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

for gene in records:
    mRNA = transcription(gene)
    mRNAs.extend(mRNA) #  Add a list on the end
#print(mRNAs)

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
#Dictionary for DNA
DNAdict = {'A': ['dAMP','A',331.2,8.8e-6, 1.68e-5, 8.8e-6,335.5,-41.3,-365.2,-864.9], 'C': ['dCMP','C',307.2,8.02e-3, 3.71e-5,8.02e-3,188.1,33.3,-692.8,-1212.9],
           'T': ['dTMP','T',322.2,6.87e-3, 1.18e-5,-1434.2,8.4,-163.8,6.87e-3,-842.6], 'G': ['dGMP','G',347.2,5.1e-5,5.1e-5,5.1e-5,293.7,56.1,-596.4,-1115.8]}

#Dictionary for RNA
RNAdict = {'A': ['AMP','A',347.2,2.8e-4,4.23e-5,8.12e-5,335.5,-41.3,-529.3,-1020.5], 'C': ['CMP','C',323.2,3.6e-4,1.18e-5,5.18e-6,188.1,33.3,-862.5,-1368.5],
           'U': ['UMP','U',324.2,6.6e-4,1.45e-5,1.45e-5,-78.5,-92.4,-1091.9,-1583.1], 'G': ['GMP','G',363.2,2.4e-5,1.81e-5,1.02e-5,293.7,56.1,-764.8,-1271.4]}

"""----------------------------G. DNA----------------------------------------"""
#--------------------------G.1 Initialise values--------------------------------
genome_lenght=0
totweight=0
#---------------------For the Energy of each sequence---------------------------
"""Environmental conditions (McCollom & Amend 2005)"""
#With the bucket method
DNA_totdGf_env=0
DNA_geneGs_env=[]
DNA_totdGr_env= 0

"""Cell conditions (E Noor 2013 - eQuilibrator)"""
#With the bucket method
DNA_totdGf_cell=0
DNA_geneGs_cell=[]
DNA_totdGr_cell= 0

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
    (CA * DNAdict['A'][2]) + (CT * DNAdict['T'][2]) + (CC * DNAdict['C'][2]) + (CG * DNAdict['G'][2]) - (Nn*18))
    totweight+= gene_weight

    #---------------------G.4 Miscellaneous values of the genome--------------------
    #get average values for molecular weight and gene lenght
    nucleotide_avmw= 327 #Da
    gene_meanlenght= genome_lenght/len(genelist) #Average lenght of the strands
    avG= totweight/len(genelist) #average gene size in Daltons
    """Concentrations"""
    host_volume= 1 #um3
    #To get the moles
    moles= genome_lenght/6.022e23
    #To get mol per litre
    DNA_mol_L= genome_lenght/(6.022e23 * host_volume * 1e-15)

    #--------G.5 Energy of the genome at different conditions and methods-------
    #------------------------G.4.1 Cellular conditions-------------------
    #Cell conditions are considered here as the standard conditions plus a
    #pH of 7 and a default ionic strength of 0.1 M. (in an aqueous solution and with a net charge of 0)
    #Adittionally, the formation energies where calculated as if the reactants had a concentration of 1M

    #sigmaB_cell will get the energy of formation of the bases in the sequence at cellular conditions
    sigmaB_cell = (
    (CA * DNAdict['A'][6]) + (CT * DNAdict['T'][6]) + (CC * DNAdict['C'][6]) + (CG * DNAdict['G'][6]))

    #sigmaN_cell will get the energy of formation of the nucleotides in the sequence at cellular conditions
    sigmaN_cell = (
        (CA * DNAdict['A'][8]) + (CT * DNAdict['T'][8]) + (CC * DNAdict['C'][8]) + (CG * DNAdict['G'][8]))

    """----------Energy calculation at cellular conditions--------------"""
    #The bucket method will calculate the energy based on the dGf of the phosphate
    #backbone and the dGf of the base without considering interactions

    """Gibbs Formation energy at cellular conditions"""
    # dGf using the equation for the at cellular conditions
    DNA_dGf_chain_cell= (Nn*DNA_dGf_PhBB_cell) + sigmaB_cell
    DNA_totdGf_cell+= DNA_dGf_chain_cell

    """Gibbs Reaction Energy at ceullar conditions"""
    DNA_dGr_cell= ((((Nn-1)*dGf_H2O) + DNA_dGf_chain_cell) - sigmaN_cell)
    DNA_geneGs_cell.append(DNA_dGr_cell)
    DNA_totdGr_cell+= DNA_dGr_cell

    """Reation quotient at cellular conditions"""
    PC = DNA_mol_L
    DNA_ReactionGibbs_cell = []
    TK = np.linspace(253.15,400)

    for T in TK:
        #make a list of the free energies for each gene
        DNA_reactionGs_cell = []
        #Bucket method at cellular conditions
        for gene, DNA_stdG in zip(genelist, DNA_geneGs_cell):
            DNA_lnQ_cell = math.log(PC/len(genelist))
            for base in gene:
                if base in DNAdict.keys():
                                    if celltype == "bacteria":
                                        DNA_lnQ_cell -= math.log(DNAdict[base][3])
                                    elif celltype == "mammalian":
                                        DNA_lnQ_cell -= math.log(DNAdict[base][4])
                                    elif celltype == "yeast":
                                        DNA_lnQ_cell -= math.log(DNAdict[base][5])
            DNA_reactionGs_cell.append((DNA_stdG+(R*(T)*DNA_lnQ_cell)))


        #take the average, which is in kJ/mol
        DNA_avgReactionG_kJmol_cell = statistics.mean(DNA_reactionGs_cell)
        # convert to kJ / dry g
        #avG is the average gene size
        #Bucket
        DNA_avgReactionG_kJdryg_cell = DNA_avgReactionG_kJmol_cell * 1/(avG)
        DNA_ReactionGibbs_cell.append(DNA_avgReactionG_kJdryg_cell)


    #------------------------G.4.2 Environmental conditions-------------------
    #The environmental conditions are considered here as the standard conditions
    #without accounting for pH, ionic strength or any other cellular factors
    #(aqueous solution abd with a net charge of -2/-3 per molecule)

    #sigmaB_env will get the energy of formation of the bases in the sequence at environental conditions
    sigmaB_env = (
    (CA * DNAdict['A'][7]) + (CT * DNAdict['T'][7]) + (CC * DNAdict['C'][7]) + (CG * DNAdict['G'][7]))

    #sigmaN_env will get the energy of formation of the nucleotides in the sequence at environental conditions
    sigmaN_env = (
        (CA * DNAdict['A'][9]) + (CT * DNAdict['T'][9]) + (CC * DNAdict['C'][9]) + (CG * DNAdict['G'][9]))

    """----------Energy calculation at environmental conditions--------------"""
    #This method will calculate the energy based on the dGf of the phosphate
    #backbone and the dGf of the base without considering interactions

    """Gibbs Formation energy"""
    # dGf using the equation for the the bucket method at environmental conditions
    DNA_dGf_chain_env= (Nn*DNA_dGf_PhBB_env) + sigmaB_env
    DNA_totdGf_env+= DNA_dGf_chain_env

    """Gibbs Reaction Energy"""
    DNA_dGr_env= ((((Nn-1)*dGf_H2O) + DNA_dGf_chain_env) - sigmaN_env)
    DNA_geneGs_env.append(DNA_dGr_env)
    DNA_totdGr_env+= DNA_dGr_env

    """Reation quotient at environental conditions"""
    PC = DNA_mol_L
    DNA_ReactionGibbs_env = []
    TK = np.linspace(253.15,400)

    for T in TK:
        #make a list of the free energies for each gene
        DNA_reactionGs_env = []
        #Bucket method at environmental conditions
        for gene, DNA_stdG in zip(genelist, DNA_geneGs_env):
            DNA_lnQ_env = math.log(PC/len(genelist))
            for base in gene:
                if base in DNAdict.keys():
                                    if celltype == "bacteria":
                                        DNA_lnQ_env -= math.log(DNAdict[base][3])
                                    elif celltype == "mammalian":
                                        DNA_lnQ_env -= math.log(DNAdict[base][4])
                                    elif celltype == "yeast":
                                        DNA_lnQ_env -= math.log(DNAdict[base][5])
            DNA_reactionGs_env.append((DNA_stdG+(R*(T)*DNA_lnQ_env)))


        #take the average, which is in kJ/mol
        DNA_avgReactionG_kJmol_env = statistics.mean(DNA_reactionGs_env)
        # convert to kJ / dry g
        #avG is the average gene size
        #Bucket
        DNA_avgReactionG_kJdryg_env = DNA_avgReactionG_kJmol_env * 1/(avG)
        DNA_ReactionGibbs_env.append(DNA_avgReactionG_kJdryg_env)

#------------------------I. Printing values of the genome-----------------------
print('-------------------------------------------Genome information-------------------------------------------------------------------')
if celltype == "bacteria":
    print('These are the values of the genome in the bacterial cell at cellular and environental conditions')

elif celltype == "mammalian":
    print('These are the values of the genome in the mammalian cell at cellular and environental conditions')

elif celltype == "yeast":
    print('These are the values of the genome in the yeast cell at cellular and environental conditions')

print('--> Energy at Cellular conditions:')
print('Average Gibbs Free energy of Formation [dGf]:', DNA_totdGf_cell/(len(genelist)), 'KJ per sequence')
print('Average Gibbs Free energy of reaction [dGr]:', DNA_totdGr_cell/(len(genelist)), 'KJ per sequence')
print('The energy to build the genome adjusted to 1 gram at 298.15K:')
print(DNA_ReactionGibbs_cell[45],'KJ/g considering the absolute intracellular metabolite concentrations for the selected cell type and the sequence lenght')


print('--> Energy at Environmental conditions:')
print('Average Gibbs Free energy of Formation [dGf]:', DNA_totdGf_env/(len(genelist)), 'KJ per sequence')
print('Average Gibbs Free energy of reaction [dGr]:', DNA_totdGr_env/(len(genelist)), 'KJ per sequence')
print('The energy to build the genome adjusted to 1 gram at 298.15K:')
print(DNA_ReactionGibbs_env[45],'KJ/g considering the absolute intracellular metabolite concentrations for the selected cell type and the sequence lenght')



print('Average values per gene sequence (Wheight and lenght):')
print('The weight of the double strand genome, based on the input file is:', totweight,'Da')
print('The average weight of the genes or genome of the input file is' ,avG,'Da with', gene_meanlenght,'nucleotides in average')
print('The concentration is',DNA_mol_L,' moles of DNA per litre')


"""------------------------------J. RNA--------------------------------------"""
#-------------------------J.1 Initialise values---------------------------------
transcriptome_lenght=0
transcriptome_totweight=0
#---------------------For the Energy of each sequence---------------------------
"""Environmental conditions (McCollom & Amend 2005)"""
RNA_totdGf_env=0
RNA_geneGs_env=[]
RNA_totdGr_env= 0

"""Cell conditions (E Noor 2013 - eQuilibrator)"""
RNA_totdGf_cell=0
RNA_geneGs_cell=[]
RNA_totdGr_cell= 0
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
        (RA * RNAdict['A'][2]) + (RU * RNAdict['U'][2]) + (RC * RNAdict['C'][2]) + (RG * RNAdict['G'][2]) - (RNA_Nn*18))
        """We take away the H2O weight for every nucleotide"""
        transcriptome_totweight+= RNA_weight

        #--------------J.5 Miscellaneous values for the transcriptome-------------------
        #get average values for molecular weight and lenght
        RNA_nucleotide_avmw= 339.5 #Da
        messenger_meanlenght= transcriptome_lenght/len(mRNAs) #Average lenght of the strands
        avm= transcriptome_totweight/len(mRNAs) #average gene size in Daltons
        """Concentrations"""
        host_volume= 1 #um3
        #To get the moles
        RNA_moles= transcriptome_lenght/6.022e23
        #To get mol per litre
        RNA_mol_L= transcriptome_lenght/(6.022e23 * host_volume * 1e-15)


        #--------G.4 Energy of the proteome at different conditions and methods-------

        #------------------------G.4.1 Cellular conditions-------------------
        #Cell conditions are considered here as the standard conditions plus a
        #pH of 7 and a default ionic strength of 0.1 M. (in an aqueous solution and with a net charge of 0)
        #Adittionally, the formation energies where calculated as if the reactants had a concentration of 1M

        #sigmaB_cell will get the energy of formation of the bases in the sequence at cellular conditions
        RNA_sigmaB_cell = (
        (RA * RNAdict['A'][6]) + (RU * RNAdict['U'][6]) + (RC * RNAdict['C'][6]) + (RG * RNAdict['G'][6]))

        #sigmaN_cell will get the energy of formation of the nucleotides in the sequence at cellular conditions
        RNA_sigmaN_cell = (
        (RA * RNAdict['A'][8]) + (RU * RNAdict['U'][8]) + (RC * RNAdict['C'][8]) + (RG * RNAdict['G'][8]))

        """----------Energy calculation at cellular conditions--------------"""
        #The bucket method will calculate the energy based on the dGf of the phosphate
        #backbone and the dGf of the base without considering interactions

        """Gibbs Formation energy at cellular conditions"""
        # dGf using the equation for the at cellular conditions
        RNA_dGf_chain_cell= (Nn*RNA_dGf_PhBB_cell) + sigmaB_cell
        RNA_totdGf_cell+= RNA_dGf_chain_cell

        """Gibbs Reaction Energy at ceullar conditions"""
        RNA_dGr_cell= ((((Nn-1)*dGf_H2O) + RNA_dGf_chain_cell) - sigmaN_cell)
        RNA_geneGs_cell.append(RNA_dGr_cell)
        RNA_totdGr_cell+= RNA_dGr_cell

        """Reation quotient at cellular conditions"""
        PC = RNA_mol_L
        RNA_ReactionGibbs_cell = []
        TK = np.linspace(253.15,400)

        for T in TK:
            #make a list of the free energies for each gene
            RNA_reactionGs_cell = []
            #Bucket method at cellular conditions
            for gene, RNA_stdG in zip(genelist, RNA_geneGs_cell):
                RNA_lnQ_cell = math.log(PC/len(genelist))
                for base in gene:
                    if base in RNAdict.keys():
                                        if celltype == "bacteria":
                                            RNA_lnQ_cell -= math.log(RNAdict[base][3])
                                        elif celltype == "mammalian":
                                            RNA_lnQ_cell -= math.log(RNAdict[base][4])
                                        elif celltype == "yeast":
                                            RNA_lnQ_cell -= math.log(RNAdict[base][5])
                RNA_reactionGs_cell.append((RNA_stdG+(R*(T)*RNA_lnQ_cell)))


            #take the average, which is in kJ/mol
            RNA_avgReactionG_kJmol_cell = statistics.mean(RNA_reactionGs_cell)
            # convert to kJ / dry g
            #avG is the average gene size
            #Bucket
            RNA_avgReactionG_kJdryg_cell = RNA_avgReactionG_kJmol_cell * 1/(avG)
            RNA_ReactionGibbs_cell.append(RNA_avgReactionG_kJdryg_cell)


        #------------------------G.4.2 Environmental conditions-------------------
        #The environmental conditions are considered here as the standard conditions
        #without accounting for pH, ionic strength or any other cellular factors
        #(aqueous solution abd with a net charge of -2/-3 per molecule)

        #sigmaB_env will get the energy of formation of the bases in the sequence at environental conditions
        RNA_sigmaB_env = (
        (RA * RNAdict['A'][7]) + (RU * RNAdict['U'][7]) + (RC * RNAdict['C'][7]) + (RG * RNAdict['G'][7]))

        #sigmaN_env will get the energy of formation of the nucleotides in the sequence at environental conditions
        RNA_sigmaN_env = (
        (RA * RNAdict['A'][9]) + (RU * RNAdict['U'][9]) + (RC * RNAdict['C'][9]) + (RG * RNAdict['G'][9]))

        """----------Energy calculation at environmental conditions--------------"""
        #This method will calculate the energy based on the dGf of the phosphate
        #backbone and the dGf of the base without considering interactions

        """Gibbs Formation energy"""
        # dGf using the equation for the the bucket method at environmental conditions
        RNA_dGf_chain_env= (Nn*RNA_dGf_PhBB_env) + sigmaB_env
        RNA_totdGf_env+= RNA_dGf_chain_env

        """Gibbs Reaction Energy"""
        RNA_dGr_env= ((((Nn-1)*dGf_H2O) + RNA_dGf_chain_env) - sigmaN_env)
        RNA_geneGs_env.append(RNA_dGr_env)
        RNA_totdGr_env+= RNA_dGr_env

        """Reation quotient at environental conditions"""
        PC = RNA_mol_L
        RNA_ReactionGibbs_env = []
        TK = np.linspace(253.15,400)

        for T in TK:
            #make a list of the free energies for each gene
            RNA_reactionGs_env = []
            #Bucket method at environmental conditions
            for gene, RNA_stdG in zip(genelist, RNA_geneGs_env):
                RNA_lnQ_env = math.log(PC/len(genelist))
                for base in gene:
                    if base in RNAdict.keys():
                                        if celltype == "bacteria":
                                            RNA_lnQ_env -= math.log(RNAdict[base][3])
                                        elif celltype == "mammalian":
                                            RNA_lnQ_env -= math.log(RNAdict[base][4])
                                        elif celltype == "yeast":
                                            RNA_lnQ_env -= math.log(RNAdict[base][5])
                RNA_reactionGs_env.append((RNA_stdG+(R*(T)*RNA_lnQ_env)))


            #take the average, which is in kJ/mol
            RNA_avgReactionG_kJmol_env = statistics.mean(RNA_reactionGs_env)
            # convert to kJ / dry g
            #avG is the average gene size
            #Bucket
            RNA_avgReactionG_kJdryg_env = RNA_avgReactionG_kJmol_env * 1/(avG)
            RNA_ReactionGibbs_env.append(RNA_avgReactionG_kJdryg_env)

#------------------------I. Printing values of the transcriptome-----------------------
print('-------------------------------------------Transcriptome information-------------------------------------------------------------------')
if celltype == "bacteria":
    print('These are the values of the transcriptome in the bacterial cell at cellular and environental conditions')

elif celltype == "mammalian":
    print('These are the values of the transcriptome in the mammalian cell at cellular and environental conditions')

elif celltype == "yeast":
    print('These are the values of the transcriptome in the yeast cell at cellular and environental conditions')


print('--> Energy at Cellular conditions:')
print('Average Gibbs Free energy of Formation [dGf]:', RNA_totdGf_cell/(len(genelist)), 'KJ per sequence')
print('Average Gibbs Free energy of reaction [dGr]:', RNA_totdGr_cell/(len(genelist)), 'KJ per sequence')
print('The energy to build the transcriptome adjusted to 1 gram at 298.15K:')
print(RNA_ReactionGibbs_cell[45],'KJ/g considering the absolute intracellular metabolite concentrations for the selected cell type and the sequence lenght')

print('--> Energy at Environmental conditions:')
print('Average Gibbs Free energy of Formation [dGf]:', RNA_totdGf_env/(len(genelist)), 'KJ per sequence')
print('Average Gibbs Free energy of reaction [dGr]:', RNA_totdGr_env/(len(genelist)), 'KJ per sequence')
print('The energy to build the transcriptome adjusted to 1 gram at 298.15K:')
print(RNA_ReactionGibbs_env[45],'KJ/g considering the absolute intracellular metabolite concentrations for the selected cell type and the sequence lenght')

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

ax.plot(TK,DNA_ReactionGibbs_cell, label='Genome energy at cell conditions', linestyle='dashed', c='r', linewidth=3)
ax.plot(TK, DNA_ReactionGibbs_env, label='Genome energy at environment conditions', marker='o', c='b', linewidth=3)

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

ax.plot(TK,RNA_ReactionGibbs_cell, label='Transcriptome energy at cell conditions', linestyle='dashed', c='r', linewidth=3)
ax.plot(TK, RNA_ReactionGibbs_env, label='Transcriptime energy at environment conditions', marker='o', c='b', linewidth=3)

ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
ax.set_xlabel('Temperature [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xlim(270, 400)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()
