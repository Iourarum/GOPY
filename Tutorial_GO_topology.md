The tutorial is complete, although some parts may be hard to read as it's still "raw". Contact me if you get stuck - this way I will also know what needs improvement.
(I think this is the tutorial I have used when I first needed to make such topologies: https://ambermd.org/tutorials/advanced/tutorial1_orig/ . It might be useful in case I missed anything here, as in explanation-wise for topology generation using xleap/tleap.)

###### Please read the bolded paragraph below:
##### Link to GOPY paper: https://doi.org/10.1016/j.softx.2020.100586 . We appreciate all citations! This tutorial is here to help those struggling with creating GROMACS or AMBER topology and coordinate files for Graphene Oxide (GO). In the following steps, GO is assumed to be functionalized with COOH, -O- and OH groups. We'll take you step by step, as if doing actual research, as we did it (to the best of our memory).




### Tutorial on Creating Graphene Oxide Topology Files (AMBER or GROMACS - to be used with Forcefields from the AMBER-family or OPLSAA)

###### Pre-requisites: Python (prefferably 3.9 and not 3.8), AmberTools (from source or through "conda install -c conda-forge ambertools"). We will make use of tleap or xleap and parmed from AmberTools. In tleap / xleap we create .lib files for defining a residue and use .frcmod files to store parameters not included in a forcefield by default. Thus, please make sure you have these pre-requisites installed. Sometimes, when installing AmberTools through conda we have encountered some issues in the past, which seemed to occur randomly as they did not appear on all the machines we performed the installation on, such as: parmed might not be installed together with AmberTools, you can try to write to me regarding these too - I will do my best to help.

We assume you are familiar with the structure of a PDB file (you may see it quickly in the GOPY paper) and somewhat familiar with the files of GROMACS, .top (topology) and .gro (coordinates). You don't need to be familiar with .lib and .frcmod files to understand the workflow, though it helps if you are.

###### What's important is that when we made GOPY, we initially planned for one CX atom to represent one full GGG residue (pristine graphene), C, O, O, H to represent a C1A residue (carboxyl, 4 atoms only), O to represent one E1A residue (epoxy, 1 atom only) and H, O to represent one H1A residue (hydroxyl 2 atoms). As it is today, GOPY names the carbon atoms on the graphene plane to which functional groups are connected as "CY" instead of "CX" (or one "CY" and one "CZ" in the case of epoxy). Essentially, CY and CZ are CX atoms. Long story short, in the current form of GOPY, carboxyl residues now contain the CY atom and the C, O, O and H (total 5 atoms), epoxy residues contain the O and the two carbon atoms on the graphene plane (total 3 atoms), called CY and CZ, hydroxyl residues contain the H, O and the CY atom on the graphene plane (total 3 atoms).

###### The next issue one has to be aware of regards the choice of parameters for the four different residues: GGG - pristine graphene, C1A - carboxyl, E1A - epoxy, H1A - hydroxyl. In the study in which we used the GOPY script ( https://doi.org/10.3390/coatings10030289 ) we cited these two studies (similar authors) for parameters for pristine graphene and functional groups:

https://doi.org/10.1088/0022-3727/47/50/505401

https://doi.org/10.1088/0022-3727/48/27/275402

(Note: Renaming some of the CX atoms to CY and CZ, as explained in the earlier paragraph, should make it easier for one to use the parameters present in this third study, here https://doi.org/10.1002/chem.201701733 . This will not be part of this tutorial, unless we decide to expand it slightly.)

#### Creating a GO PDB file using GOPY
##### Please download GOPY. We are going to create a small, 5 nm x 5 nm GO molecule using the formula C20(OH)2(-O-)2(COOH)1. This formula was described as a suitable model in the first two studies we have provided the DOIs for. 
0. Navigate to the folder containing GOPY.py (if you encounter any error using GOPY, such as the encoding one, please see the readme file on github, use "python GOPY.py help" or "python3 GOPY.py help" for instructions (preferably use Python3.9 and not Python3.8, should work with other versions sometimes too though.)
1. Create a 5 nm x 5 nm PG PDB file using: "python GOPY.py generate_PG 50 50 PG.pdb" or "python3 ..."
2. VMD detects 1009 atoms when visualising PG.pdb, meaning we need 50 COOH (1009/20) and 101 OH and -O- groups (1009/10) to add. Use: "python GOPY.py generate_GO path/to/file.pdb 50 101 101 GO.pdb". Now you should have a GO PDB file as described. If you encounter any errors here make sure everything you typed makes sense. We do not expect any errors so far regarding the use of GOPY.

#### Parameters for PG, -COOH, -OH, -O- ... we need these, follow the instructions below.
###### Essentially, the authors state regarding functional groups: "The parameters of hydroxyl, carboxyl and epoxy groups were taken from the AMBER99SB force field for serine, glutamic acid and dialkyl ether, respectively.". I had to look in GROMACS force field files to figure the parameters. How one interprets their sentence depends, I guess. For example, if you go in your GROMACS installation /share/top/amber99sb.ff, open aminoacids.rtp, find "GLH" - that's where I picked the parameters for carboxyl (meaning partial charges, atom types etc). Well, following this line of thought I eventually ended up with the parameters I think are right for all three functional groups. 
###### Check the GO_tutorial folder. You will find all lib files, frcmod file etc. over there.
###### I have uploaded here library files (.lib) for the COOH functional group, but not COO. One can pick parameters for a COO functional group similarly as above. Have a look through the existing parameters. For the carbon atoms making up the pristine graphene, the authors state all parameters. All one needs to do is to convert them to GROMACS units (if I remember right). This aspect can be searched and understood easily.

0. Open the earlier generated GO.pdb file in an editor of choice (gedit, nano etc.)
1. Single out one CX/GGG atom, one carboxyl, one hydroxyl and one epoxy functional groups and save these as independent PDB files, you may use those present in the GO_tutorial folder.
2. To be able to generate AMBER coordinates and topology files, we will need to follow some steps in xleap.
2. 1 Open xleap through typing "xleap" in terminal (xleap is installed as part of AmberTools). You may also use tleap if you preffer.
2. 2 Import the AMBER99SB default parameters by typing 
###### source oldff/leaprc.ff99SB. 
Attention: Xleap has its own set of quirks you need to take into account. Keep your cursor inside the typing area on xleap and keep NumLock off) - assuming you have AmberTools installed in a conda environment called "AmberTools20", you may look for other parameter sets at: ../anaconda3/envs/AmberTools20/dat/leap/cmd/ (other than those in oldff/leaprc.ff99SB ; as we have been involved in nucleic acid simulations we used parmbsc1 or "source leaprc.DNA.bsc1")
4. 3 Moving on, open all PDB files of interest, using:
###### mol_GGG = loadpdb /path/to/GGG.pdb
###### mol_C1A = loadpdb /path/to/C1A.pdb
###### mol_H1A = loadpdb /path/to/H1A.pdb
###### mol_E1A = loadpdb /path/to/E1A.pdb
2. 4 Edit parameters for each new molecule
You should be familiar with this from other AMBER tutorials.
Or you may use the lib files we have provided in the GO_tutorial folder - these are identical to the ones we used for this publication: https://doi.org/10.3390/coatings10030289 .
If you prefer to do it manually, do:
###### saveoff mol_GGG GGG.lib
###### saveoff mol_C1A C1A.lib
###### saveoff mol_E1A E1A.lib
###### saveoff mol_H1A H1A.lib
Now you will have your own lib files.
2. 5 The frcmod file present in GO_tutorial will look like a mess because the CX, CY and CZ atoms are, technically, all sharing the same parameters for this tutorial (or method described here). You can maybe improve them, possibly by using the parameters present in this study: https://doi.org/10.1002/chem.201701733 ... but I can't certify it's better or anything, I just figured you would be able to make use of the different names of these atoms to even further customize their parameters. You'll know what to do if you are experienced with MD.

#### Creating the actual GO topology
###### This is the section you actually came for (and the files, maybe).

Thus, open xleap and write (or tleap if you don't need GUI):
###### source leaprc.DNA.bsc1 (this loads forcefield parameters or use source oldff/leaprc.ff99SB depending on which forcefield you would like to use)
###### loadoff /path/to/library.lib (load all four: GGG, C1A, E1A, H1A)
###### loadamberparams /path/to/file.frcmod (load the frcmod file with the parameters)
###### mol = loadpdb /path/to/GO.pdb (load your GO molecule)
###### edit mol (opens the molecule structure in a graphical format)
##### Now, the way this was designed you should have all parameters imported well and no editing of parameters should be necessary. This saves so much time...
###### bondbydistance mol (this will create all bonds based on distance, you should see it in the previously opened window if using xleap)
###### saveamberparm mol mol.prmtop mol.inpcrd (prmtop is the topology file, inpcrd is the coordinate file, they can have any name really before the dot)

Then edit and use the python script attached (ambertogro.py) to obtain the conversion gromacs topologies.

- You will need AmberTools installed, either from source or through conda, both are fine from what I remember. From AmberTools you will mostly need xleap (GUI) or tleap (just terminal), xleap is more useful cause you can graphically check your structure (you will see). Parmed is used, I think, in the conversion script (this script already exists). Note: As you will see/know about xleap,... it may act "weird", make sure your cursor is placed in the window/box you want to type in, check your caps lock etc. (if anything behaves weird).
- The second thing you will need are library files and frcmod files. As you can see in GOPY, I chose to consider each CX atom as a full residue, so one GGG residue has only one atom. As I mentioned previously, initially I considered each COOH residue to contain just a C, O, O and a H atom, but later I added the "CY" atoms and considered the functional group residues to also include the atom it binds to on the graphene plane. (Note: epoxy contains one CY and one CZ). Both CY and CZ were CX atoms and were just renamed. So one will need one library file for graphene atoms, one for carboxyl, one for epoxy and one for hydroxyl - a total of four. Parameters that cannot be found in a standard set of parameters will need to be specified in the frcmod file - one file only. Thus I have attached 7 library files, 3 contain the functional groups without the CY/CZ atoms, 3 contain the functional groups with the CY/CZ atoms and one is for the pristine graphene residue GGG.lib.

(I think this is the tutorial I have used when I first needed to make such topologies: https://ambermd.org/tutorials/advanced/tutorial1_orig/ . It might be useful in case I miss anything.)

FORCE FIELD: Because I was interested in DNA simulations, I used parmbsc1. This is not available by default in GROMACS, but its easy to "install" if needed. Look here https://mmb.irbbarcelona.org/ParmBSC1/help.php?id=download and here http://www.gromacs.org/Downloads/User_contributions/Force_fields. It's preferable you will load a forcefield in xleap, maybe look here $AMBERHOME/dat/leap/cmd/ . In our case it's less relevant as the parameters will mostly be common for the available forcefields (except for what we add in the frcmods).

LIB FILES: Here we go, there are tutorials explaining the way to create lib file. If you need help, I can also try to help with anything new. I attached the seven library files as mentioned earlier. These lib files define the atoms making up a residue.
How do you do a new lib file? There must be some tutorials online. However, for example, for a COO residue, I would create a PDB file with only three atoms, COO. Then, with AmberTools activated, type "xleap" in your terminal. Then do:
mol = loadpdb /home/path/to/file.pdb
edit mol (new window apears, left click and rectangle select all atoms until they are purple/violet > then click on Edit menu > Edit selected atoms and a new table appears)
Edit the table accordingly, then click on Table menu > Save (Whats important: Type and Charge, you can write "true" in unused I think. I added a screenshot of the table for C1A/COOH as an example)
Then, click back on the main window of xleap and do:
saveoff mol mol.lib
One more step, open the mol.lib file with gedit or similar and find and replace all mol with the desired name of the residue and rename the file with that name .lib . I don't remember why this mattered.

FRCMOD FILE: This will look like a mess. Initially it looked better/more organised when I only had CX to worry about. Then I added CY. Then I added CZ today (not sure if I added all needed CZ lines yet, hope I did).

Anyway, the process to get the topologies is this:
open xleap and write
source leaprc.DNA.bsc1 (this loads forcefield parameters)
loadoff /path/to/library.lib (load all four)
loadamberparams /path/to/frcmod.frcmod (load the frcmod file and parameters)
mol = loadpdb /path/to/GO.pdb
edit mol (opens the molecule structure in a graphical format)
bondbydistance mol (this will create all bonds based on distance, you should see it in the previously opened window)
saveamberparm mol mol.prmtop mol.inpcrd (prmtop is the topology file, inpcrd is the coordinate file, they can have any name really before the dot)

Then edit and use the python script attached (ambertogro.py) to obtain the conversion gromacs topologies.

Essentially this is it. Good luck with research, hope for the best!
Don't hesitate to contact me if anyting is unclear or you think I can be of help with something.

CX           6      12.01    0.0000  A   3.39967e-01  3.59824e-01


