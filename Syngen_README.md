# synpro-syngen

**Function**
The program is designed to read a DNA sequence contained in a single file to obatain the energy to synthesise the nucleotides contained in the sequence plus the complementary chain.
Furthermore, the program transcribes the lead strand reading the first possible series of Open Reading Frames (ORF) and calculates the energy to synthesise the nucleotides of the transcriptome.

**Structure**
The program is divided in 17 blocks:

0. Initialise: Define the initial values and import the necessary functions, methods, libraries etc.

A. Transcription: Here is defined the transcription method to be called later in the program.
   The program will look for the first start codon (ATG) available in the DNA sequence and transcribe it until the next stop codon 'TAA', 'TGA' or 'TAG' in that order of appearance.
   Once, the program has transcribed a gene it will begin to look for the next start codon to transcribe the next gene and so on until parsing the whole sequence.
   All the transcripts will be stored in an independent list (mRNA).

B. Input sequence: The user provides the DNA sequence to be used by the program in a FASTA file.
   The user can choose between the three different model cells according to the input sequence

C. Parsing with biopython: This block helps the program reading and obtaining the information contained in the FASTA file.

D. Complement and transcribe: A biopython block that will obtain the complementary chain and transcribe the input DNA.

E. Nucleotides library: Different values for each nucleotide and nucleotide components at different temperatures

F. Dictionaries: Here is defined the dictionary with the values to be used of each nucleotide in DNA and RNA.

--DNA--

G. DNA: This block will do calculations for DNA in four sections:  G.1: will ge the count of each nucleotide from the input sequence
                                                                   G.2: Calculates the weight of the genome by multiplying the count of each nucleotide by its respective molecular weight.
                                                                   G.3: Gets the total energy value of each building block involved in the formation of the DNA
                                                                   G.4: Obtains the formation Gibbs energy, using a group contribution algorithm to finally obtain the reaction Gibbs energy of the the chain.
                                                                   G.5: Computes the formation and reaction Gibbs energy by splitting the equation into building blocks

H. DNA concentration: This block will generate some miscelaneous values of interest such as average length of the gene in case there are more than one sequence, the total weight of the genetic material and the DNA concentration in mol/L for the reaction quotient

I. Reaction quotient: It calculates the reaction quotient based on the DNA concentration and the building blocks concentration reported in the three different model cells

J. Molar Gibbs energy of the chain: Using the equation: ∆G= ∆G°+RTlnQ this block will use the reaction energy calculated in G.4, a different temperature from 275K to 400K and the reaction quotient.

--RNA--

K. RNA: This block will do calculations for the RNA messengers in four sections:  K.1: will ge the count of each nucleotide from the input sequence
                                                                                  K.2: Calculates the weight of the proteome by multiplying the count of each nucleotide by its respective molecular weight.
                                                                                  K.3: Gets the total energy value of each building block involved in the formation of the RNAm
                                                                                  K.4: Obtains the formation Gibbs energy, using a group contribution algorithm to finally obtain the reaction Gibbs energy of the the chain.

L. RNA concentration: This block will generate some miscelaneous values of interest such as average lenght of the gene in case there are more than one sequence, the total weight of the genetic material and the DNA concentration in mol/L for the reaction quotient

M. Reaction quotient: It calculates the reaction quotient based on the RNA concentration and the building blocks concentration reported in the three different model cells

N. Molar Gibbs energy of the chain: Using the equation: ∆G= ∆G°+RTlnQ this block will use the reaction energy calculated in K.4, a different temperature from 275K to 400K and the reaction quotient.

O. Printing results

P. Plotting results


**Use**
To run this program type "python syngen.py" and follow the instructions, it will ask for a fasta file.
**Note**
This program needs an environment called reaktoro to run properly.
You will first need to install miniconda or conda to use different environments and then install reaktoro.
Once having conda installed type:
conda install reaktoro
Or go to the reaktoro page for further information


**Installing Biopython:**
If Biopython is already installed you can run the code normally. To install biopython before running the program type: pip install Biopython

**Limitations:**
1. It can only read fasta files.
2. It reads the first Open Reading Frame (ORF) available, excluding other potential ORFs in the sequence.
3. It is no considering other molecules at the moment.
4. Most of the values are adjusted to E. coli cells
5. It needs biopython and a specific environment to run.
