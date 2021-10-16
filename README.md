# synpro-syngen (also known as syncell)

### Welcome to syncell, a collection of python code to estimate the energy required to build biomacromolecules and organisms.

This repository is a work-in-progress, but check back for big updates soon!

---------------------------------------------------------------------------------------------------------------------------
Getting ready to use Syncell

To run Syngen and Synpro contained in Syncell is necessary to install the reaktoro environment and the biopython, matplotlib, progressbar and tqdm packages in case you don't have them already.
To do this simply follow the next steps:

Reaktoro

  0. In the terminal choose the work carpet where the programs and your files are contained
  1. Instal miniconda from the following link https://docs.conda.io/en/latest/miniconda.html
  2. Download the 64-bit bash file
  3. In the terminal window type 'bash' and drag the downloaded file next to bash in the terminal window
  4. Follow the instructions to install miniconda and accept the terms and conditions
  4. Next, add the conda-forge channels by typing in the terminal: conda config --append channels conda-forge
  5. Once done this, install reaktoro by typing: conda install reaktoro=1
  6. To create the environment for reaktoro type: conda create -n rkt reaktoro
  7. Before running Syncell is necessary to activate reaktoro by typing: conda activate rkt or conda activate (rkt path)
Biopython
  8. In the terminal window type: pip install biopython
Matplotlib
  9. In the terminal type: python -m pip install -U matplotlib
ProgressBar
  10. Similarly, type: pip install progressbar2
ProgressBar
  11. Finally, type: pip install tqdm


Using Syncell
  1. Type: pyhton syncell.py
  2. type the name of the fasta files you want to use, one for proteins and one for nucleic acids
