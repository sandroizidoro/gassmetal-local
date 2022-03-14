# GASS-Metal

GASS-Metal is a method based on a genetic algorithm to search for similar metal-binding sites in proteins. In addition to finding similar metal-binding sites, the method can find inter-domain sites and perform not exact matches using a substitution matrix (conservative mutations).

In this new version, GASS-Metal uses parallel genetic algorithms to create an initial population (seeds), improve accuracy and decrease processing time. 

## Installation


GASS-Metal does not require a installation, however the user's system needs to be able to compile c++ source code (g++) and some python packages to execute the main script. This packages can be installed through pip, using:

```shell
pip install requests json
```


## Usage

The GASS-Metal binary can be executed directly with the arguments:
```shell
./thread_GASS <path_to_target*> <path_to_template*> <centroid_x> <centroid_y> <centroid_z>
```
Where: 

+ path_to_target = path to a folder containing a protein.dat file of the target protein.
+ path_to_template = path to a template folder containing the reference .dat file, a Template.txt file, a GA_Conf.txt file and a SubstitutionMatrix.txt file.
+ centroid(optional) = values of the protein center in 3D space to split the data in 4 quadrants.

But the recommended way to execute is through the python script: 

```shell
python rungass.py <reference_protein*> <template_site*> <mutations> <target_protein*> <reference_atom*>
```
Where: 

+ reference_protein = path to a pdb file, pdb code or UniProt(AlphaFold) code of a protein to be used as template.
+ template_site = the residues that form the template site on the format: "Residue name, position on PDB, chain;". Example: **HIS,57,E; ASP,102,E; SER,195,E;**
+ mutations(optional) = swap of residues that may occur on this specific template site on the format "Residue name, residue name;". Example: **HIS,CYS;GLU,ASP;SER,HIS;**
+ target_protein = path to a pdb file, pdb code or UniProt(AlphaFold) code of the protein where the site will be searched.
+ reference_atom = the atom to be used as reference. The values may be only "**CA**" or "**LHA**".

Example:
```shell
python rungass.py 3nos "CYS,94,B;CYS,99,A;CYS,99,B;CYS,94,A;" 3nos CA
```
The python script then will automatically download the pdb files from the [RCSB PDB](https://www.rcsb.org/) or [AlphaFold](https://alphafold.com/) site, extract the reference atoms position from the structure files, calculate the centroid of the target protein, run the GASS-Metal compiled binary (compiling it if necessary) and show the output on the terminal.

The script however is designed to be executed on a **Linux machine**, the execution on a Windows machine is not granted by now.
## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)