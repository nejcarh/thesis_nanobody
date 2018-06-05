## Description:
After applying crystal symmetry transformations, this script searches for atoms within 4 Å of any cobalt ion. If at least two atoms come from chains created with different symmetry operations, the cobalt ion is presumed to lie between two chains on a crystal contact. All such structures are written in an output file along with the names of atoms and their parent residues and with predicted secondary structure (with the help of DSSP). Manual overview was required to check for the types of metal contacts found.


## Usage:
Place the **script** in a folder alongside the **DSSP executable** and all the **protein structure files** you wish to check. Only .PDB format is currently supported, though PDBx/mmCIF support could probably be added easily using Biopython's MMCIFParser module. I recommend installing one of the conda distributions since they contain Python with all the necessary libraries.

## Quick modification:
`maxdist = 4` ... change to modify the maxiumum search distance between around metal atoms in   
`if a_original.get_id() == 'CO':` ... change to modify type of atom to look for in metal contacts

## Output:
- **log_full**: lists all structures checked, number of models found, and all atoms surrounding cobalt on crystal contacts
- **log_hits**: only lists the names of structures, where a cobalt atom was found on crystal contacts
- **log_hits_details**: only lists structures in log_hits along with the corresponding data from log_full

## Notes on interpretation:
I recommend observing log_hits for a list of structures with cobalt on crystal contacts and log_hits_details for the nature of said contacts. The full log was mainly used for debugging less cooperative structures. Metal atoms on crystal contacts are often duplicated with different symmetry operations - in this case they will show up a very small distance away (0-0.5 Å) with a "...CLASH?" warning to the right. Non-protein (DNA) structures will correctly recognize cobalt on crystal structures and display "...file error?" message because of an exception in DSSP. Secondary structures of parent residues are given in "struct: " column in standard [DSSP notation](https://en.wikipedia.org/wiki/Protein_secondary_structure#DSSP_classification). Since the script does not further classify the number or nature of metal contacts you are required to check the structures manually. Chimera, for example, allows for easy visualisation of "Higher order structures". To immediately display contacts of interest, use command:  
`show Co z<4; select Co z<4; focus #0@Co; clip off; bondrepr stick; color byhet`

## Possible future improvements:
The output files require manual overview which can be tedious and inexact when large amounts candidate structures are found. An automatic scoring/sorting system and better presentation system would be a priority in an updated version of this script. Ease of use and portability could be enhanced with an external parameters file. Other than that, I believe that the core concept works well to locate specific atoms on crystal contacts. Support for multi-atom ligands could also be added.

## Notes on licensing:
The [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/index.html) executable is required for secondary structure prediction and **is not my work**. Please observe its license and credit/cite the original authors if you decide to use it yourself.
