## Description:
This one identifies pairs of residues from different chains in the unit cells that are less than 12 Å (or any other distance) appart. Lists of all such pairs along with their Cα-Cα distance and Cα-Cβ angle give an approximate sense of the size and nature of the contacts found. The script proved to be of limited use when applied to the majority of β-sheet containing structures since the output had to be checked by hand and any kind of comparison in terms of bonding tendency was pretty much impossible. When used on a smaller set of nanobody structures though, it did clearly identify "5IMO" as containing a sizable β-β contact surface. 

**script_beta_ANY** looks for any pairs of residues on crystal contacts and **script_beta_E** is an example of a script that only looks for pairs of glutamines.

## Usage:
Place the **script** in a folder alongside the **DSSP executable** and all the **protein structure files** you wish to check. Only .PDB format is currently supported, though PDBx/mmCIF support could probably be added easily using Biopython's MMCIFParser module. I recommend installing one of the conda distributions since they contain Python with all the necessary libraries.

## Quick modification:
`DISTANT_ATOM_CONDITION = "E"` ... [change](https://en.wikipedia.org/wiki/Protein_secondary_structure#DSSP_classification) to modify the secondary structure of contact residues   
`RES_NAME = "GLU"` ... change to modify the residue looked for in crystal contacts. Only applies to the beta_E example, beta_ANY looks for all residues. Nearby parameters decide which atoms should be used for residue direction determination and what constitutes nearby residues.

## Output:
- **log_full**: lists all structures checked, number of models found and all specified residue pairs on crystal contacts
- **log_hits**: only lists the names of structures, where specified nearby residue pairs were found on crystal contacts
- **log_hits_details**: only lists structures in log_hits along with the corresponding data from log_full

## Notes on interpretation:
I recommend observing log_hits a list of structures with specified residue pairs on crystal contacts and log_hits_details for the nature of said contacts. The full log was mainly used for debugging less cooperative structures. The number of contacts found does not fully approximate the contact surface size (even in "any residue" setting) it might just be a single residue contact between multiple different chains (e.g. A.60/B.60, C.60/D.60...) or duplicates of a single contact (e.g. A.60/B.60, B.60/A.60). Secondary structures of parent residues are given in "S: " column in standard [DSSP notation](https://en.wikipedia.org/wiki/Protein_secondary_structure#DSSP_classification). Distances between C_alpha atoms in Å are given in "d: " column. Angles between residues in degrees are given in "a: " column. A high number of pairs with angles > 90° vaguely point to a sizeable crystal contact. Since the script does not further classify or nature of crystal contacts you are required to check the structures manually. Chimera, for example, allows for easy visualisation of "Higher order structures".

## Possible future improvements:
The output files require manual overview which can be tedious and inexact when large amounts candidate structures are found. An automatic scoring/sorting system and better presentation system would be a priority in an updated version of this script. Ease of use and portability would definitely be enhanced with an external parameters file. Single/any residue search behaviour could be trivially unified in one script. Other summarizing features such as total contact area should be looked into instead of the raw residue list.

## Notes on licensing:
The [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/index.html) executable is required for secondary structure prediction and **is not my work**. Please observe its license and credit/cite the original authors if you decide to use it yourself.
