##### Update: I did a quick update recently. Personally I've been using GOPY for producing PDB files which I later imported into xleap / AmberTools. When it comes to graphene oxide (GO), an epoxy group contains two carbon atoms and both were given the same name in GOPY, "CY". I recently found that giving the second a different name worked better, so now the second is named "CZ". If you would like both of them to be "CY", either change the code or open a text editor and change all "CZ" to "CY". I will be perfectting this later (checking for errors etc). I did this because when I used xleap, I defined a lib file for the graphene atom, COOH residue, epoxy residue and hydroxyl residue and the "bondbydistance" command worked better this way. If you need help, please write to me at sebmuraru@gmail.com . I can actually help sometimes!

# GOPY: A tool for building 2D graphene-based computational models
### Paper available at: https://doi.org/10.1016/j.softx.2020.100586

GOPY is a free and open-source Python tool written in order to automate the generation of 2D graphene-based molecular models such as pristine graphene (PG), together with different forms of graphene oxide (GO/ rGO/ GO-COOH/ GO-OH/ rGO-PEG-NH2 / N-doped graphene). Key advantages to using GOPY instead of manually building the molecular models are: significantly speeding up the process, reducing potential bias due to the manual placing of functional groups and facilitating the generation of much larger and more complex models than are usually built manually. Each model is outputted in the PDB format which is easily convertible to a wide array of other molecular formats.

##### Requires Python 3, numpy and scipy.

### Quick Start
 You can use GOPY in the following manner to generate graphene-based 2D PDB models:

##### python GOPY.py generate_PG X Y file_to_save  

- Used to generate a pristine graphene layer. The X and Y dimensions are required in Å. The file_to_save represents the name under which the new PDB file should be saved.
E.g. 'python GOPY.py generate_PG 30 30' will generate a 3 nm x 3 nm PG layer. 

##### python GOPY.py generate_GO path_to_file X Y Z file_to_save  

- Used to generate a graphene oxide layer. The 'path_to_file' points to an existing pristine graphene PDB file. You may first generate a PG layer using GOPY. X corresponds to the desired number of carboxyl functional groups, Y corresponds to the number of epoxy groups and Z corresponds to the number of hydroxyl groups. The file_to_save represents the name under which the new PDB file should be saved.
E.g. 'python GOPY.py /path/to/PG.pdb 30 60 60' will generate a GO layer, attempting to place
30 carboxyl groups, 60 epoxy groups and 60 hydroxyl groups.

##### python GOPY.py generate_rGO_PEG_NH2 path_to_file X Y Z file_to_save     

- Used to generate a rGO-PEG-NH2 layer. The 'path_to_file' points to an existing graphene oxide (GO) PDB file. X Y and Z represent the percentages of carboxyl, epoxy, respectively hydroxyl groups to be removed. All removed hydroxyl and epoxy groups that were removed will be replaced by the PEG-NH2 chains. The file_to_save represents the name under which the new PDB file should be saved.
E.g. 'python GOPY.py /path/to/GO.pdb 0.6 1 1 PEG-NH2.pdb' removes 60% of carboxyl groups
and 100% of epoxy and hydroxyl groups.

##### python GOPY.py generate_hole path_to_file N R1 R2 ARG1 ARG2 C file_to_save  

- Used to generate holes in a PG layer. N represents the number of holes to be created. 
R1, R2 represent the range expressed as a list such as "[x, y]", ARG1 should be either
"u" (uni-directional) or "m" (multi-directional), C should be "a" if cleanup should
be performed and file_to_save represents the name under which to save the new file.
E.g. 'python GOPY.py /path/to/PG.pdb 10 11 20 m i a hole.pdb' reads in the PG.pdb file and
attempts to create 10 holes with a size between 11 and 20 atoms in a multi-directional manner,
not allowing the holes to touch edges. Cleanup is performed at the end and the file is saved 
as hole.pdb.

##### python GOPY.py generate_N_doped path_to_file 10 9 8 file_to_save 

- Used to generate an N-doped graphene layer. X Y and Z represent the number of N-graphitic,
N-pyridinic and N-pyrrolic atoms respectively.
E.g. 'python GOPY.py /path/to/PG.pdb 10 10 10 Ndoped.pdb' generates an N-doped molecule of the 
PG layer given as input with 10 N-graphitic atoms, 9 N-pyridinic atoms and 8 N-pyrrolic atoms.  

### Illustrative Examples
![Illustrative Examples](wast2.png)

#### Acknowledgements
This work was supported by a grant of the Ministry of Research and Innovation, Operational Program
Competitiveness Axis 1—Section E, Program co-financed from European Regional Development Fund under the
project number 154/25.11.2016, P_37_221/2015, “A novel graphene biosensor testing osteogenic potency; capturing best stem cell performance for regenerative medicine” (GRABTOP). 

#### Please cite us:
https://doi.org/10.1016/j.softx.2020.100586
