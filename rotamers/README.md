## Description:
This one looks for sets of mutations that enable coordination of a central dummy atom in the artificial nanobody dimer. Rotamers of His, Glu or Asp that could conceivably coordinate the dummy atom in a certain position are added to a list, out of which groups forming octahedral or tetrahedral geometries are extracted. Each possible combination is saved as a separate variant structure in a PDB file. Variants are grouped in candidate sets based on which mutations (e.g. K19E F79E) they require. Additionally, two more scripts were written in order to make sense of the tens of thousands of files generated.
**Sorter_distribution** simply counts the number of times mutations were used to form a variant and adds them to a table. A mutation that was used many times is presumed to be flexible in its (possibly unexpected) ability to form metal contacts.
**Sorter_leaderboards sorts** all candidate sets by multiple criteria thought to indicate whether or not their variants are viable. Categories include the number of corners filled and the number of steric clashes. Both distribution and leaderboards are easily viewable in Excel or similar software, preferably with some kind of conditional formatting.


## Usage:
Place the **script** in a folder alongside the **DSSP executable**, the input folder and all the **protein structure files** you wish to check. `input/structure/dimer.pdb` should contain the initial dimer structure. For now the search volume coordinates are hardcoded but one might be able to substitute their own structure by aligning the searchable centers in-between the chains. The main script (rotamers) can take hours to finish its job and it should deposit all newly discovered structures as it goes. After that, run the sorter scripts to generate "|" separated text files, viewable in excel.

Only .PDB format is currently supported, though PDBx/mmCIF support could probably be added easily using Biopython's MMCIFParser module. I recommend installing one of the conda distributions since they contain Python with all the necessary libraries.

## Some quick modifications (script_rotamers):
`search_volume_x,y,z = ...` ... search space dimensions in Å  
`search_spacing = 0.5` ... search step in Å (how much the atom moves in-between steps)  
`coord_radius = 2.1` ... the ideal coordination bond length in Å  
`coord_tolerance = 0.5` ... allowed deviation from ideal coordination bond length in Å  
`clash_dist = 2.5` ... the maximum distance in Å still considered a clash in steric checks  
`CA_CO_max = 8` ... maximum distance between test atom and Cα	atom, for a residue to be considered for coordination  
`candidate_ids = [17, 19, 21, 23, 68, 70, 72, 73, 74, 75, 77, 79, 81, 83, 84]` ... residues to consider for mutation
`candidate_chains = ["A", "B"]` ... chains to test for mutations, not tested for more than two  
... and more, some in code and some in outside files. Before delving deeper read the following:

**A note on repurposing script_rotamers:** Writing it was a learning experience. A priority was making it do the one thing required with little knowledge or time for good programming practice. Attempts to change anything might break stuff in fascinating ways. If anyone for any reason ever attempts to do so, feel free to contact me and/or look for improved versions somewhere on my GitHub.

## Output:
The main script will create structures of possible coordination complexes using different sets of mutated residues on given chains. File structure in /output directory should look something like this:
```
output/
├── HIS # structures containing histidines
    ├── 1,2... # folders containing structures with a certain number of histidines
        ├── MUTATION_SET (NUM) # folders containing variants requiring the same mutations, number of variants is in brackets
            ├── VARIANT.PDB # structure (.pdb) and evaluation (.txt) files for all candidate variants
├── NUM  # all structures, sorted by how filled the coordination geometries are
    ├── O,T  # octahedral/tetrahedral geometry
        ├── 1,2,3,4... # folders containing sets/variants that use varying numbers of vertices in their coordination geometry  
            ├── MUTATION_SET (N)  # folders containing variants requiring the same mutations, number of variants is in brackets
                ├──VARIANT.PDB  # structure (.pdb) and evaluation (.txt) files for all candidate variants
```
Both sorter scripts will generate "|" separated text files, viewable in Excel.


## Notes on interpretation:
The main scripts generates tens of thousands of  variant structure files. Displaying sorter_leaderboards output in a spreadsheet viewer is the only way to gain quick an overview of all results and spot the outliers in mutation set quality. Still, a manual review is required to verify structural plausability of predicted variants. It is important to keep the following properties, quirks and flaws of the main script when interpreting results:
- Chains are always treated as different and distinct. A set of mutations might enable two metal contacts around a symmetry axis and both of these may be found as separate variants, but will not be integrated into a single two-contact structure.
- As a consequence, equivalent contacts around a symmetry axis are treated *and counted* as different variants. This might make sense in terms of entropy and number of possible states, but is decieving in terms of different ways form a metal contact. In completely symmetric systems we can therefore expect a 2x overestimate of the number of *different* contacts.
- Clashes may or may not mean that the structure is impossible or even unlikely. Sometimes the offending residue has many other plausible rotamers to occupy. Whether or not that is realistic is another question. As a guide, clashes with atoms on ends of long residues might not be a problem while clashes with Cα or the backbone.
- The script will add "water" in the form of calcium atoms (because they're small and bright green in when viewed in Chimera). These "waters" will also be a part of the clash check to make sure all vertices can at least in theory be filled. However, the water-placing system is dodgy and sometimes adds too many waters or adds them to odd places. Manual review of interesting variants is recommended.
- Coordinating residues might sometimes end up too close to each other without triggering a clash. Manual review is recommended.
- The script treats a "mutation" as "an aligned rotamer succesfuly coordinated the metal ion" and does not check if a change in amino acid was even required. Therefore "K19E F79E" and "K19E D72D F79" will be treated as different mutations sets, even though the second one uses a residue that is already present.
- The way in which Glu/Asp coordination is detected might not be optimal and sometimes leads to odd angles.
- The entire process is sensitive to fairly slight variations in sidechain orientation. The script starts off with aligning all rotamers to the backbone and the Cα-Cβ bond. The latter pointing in a slightly different direction (~20-30°) can mean ~1 Å differences in coordinating atom locations case of aminoacids with longer side chains like glutamate. This seems to be a fundamental problem with this static approach. A way around it might be to take snapshots of structures at different points in a short MD simulations as starting points of parallel searches with the script.


## Possible future improvements:
The aim is to make the main script fit for for use on any set of protein chains and perhaps even expand it to other kinds of interactions (H-bonds, stacked hydrophobic, ionic, combined,...). Better candidate evaluation with tiered (more and more demanding) checks for more promising candidates will have to be added (basic force field function, hydrogens,...). Due to my inexperience starting out on this project, some of the code seems to be follow [this] design pholosophy (https://www.se.rit.edu/~tabeec/RIT_441/Resources_files/How%20To%20Write%20Unmaintainable%20Code.pdf), so a rewrite using more sane programming practices will have be the first priority. Any improved or rewritten versions I happen to finish will end up somewhere on my GitHub.

## Notes on licensing:
The [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/index.html) executable is required for secondary structure prediction and **is not my work**. Please observe its license and credit/cite the original authors if you decide to use it yourself.                        
