# GASS-Metal

GASS-Metal is a method based on a genetic algorithm to search for similar metal-binding sites in proteins. In addition to finding similar metal-binding sites, the method can find inter-domain sites and perform not exact matches using a substitution matrix (conservative mutations).

In this new version, GASS-Metal uses parallel genetic algorithms to create an initial population (seeds), improve accuracy and decrease processing time. 

## Installation


GASS-Metal does not require a installation, however the user needs to compile the source cpp files with any cpp compiler and be able to execute python scripts necessary for the pdb files preprocessing.

A simple preparation to run GASS-Metal is:
```shell
g++ thread_GASS.cpp -o thread_GASS
g++ filtroC_dados.cpp -o filtroC_dados
```


## Usage

The GASS-Metal binary can be executed directly with the arguments:
```shell
./thread_GASS <path_to_target*> <path_to_template*> <centroid_x> <centroid_y> <centroid_z>
```
Where: 

+ path_to_target = path to .dat file of the target protein.
+ path_to_template = path to a template folder containing ...TODO
+ centroid = values of the protein center in 3D space to split the data in 4 quadrants.

But the recommended way to execute is through the python script: 

```shell
./thread_GASS <path_to_target*> <path_to_template*> <reference_atom> 
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)