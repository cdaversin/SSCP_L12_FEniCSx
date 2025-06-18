# SSCP_L12_FEniCSx
Additional material for SSCP L12 Lecture - FEniCS(x) introduction

This repository provides FEniCSx code for the examples used in the L12 lecture introducing the Finite Element Method and FEniCS(x).
These scripts can be run within an environment or container with FEniCSx installed :

## Conda environment
Conda can be installed through [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions).
The environment to be used for running the example scripts can be created and activated with the following commands: 
```
conda env create environment.yml
conda activate fenicsx_sscp_2025
```
The environment can be deactivated with:
```
conda deactivate
```

## Docker container
Docker can be installed in various ways as shown in the [Docker documentation](https://docs.docker.com/).
Docker also provides a [convenience script](https://docs.docker.com/engine/install/ubuntu/#install-using-the-convenience-script).

The stable FEniCSx Docker container can be run using the following command.
Note : the option `-v` is used to share the current directory `($pwd)` of the host machine with the Docker container. Therefore, we recommand to run the command within the `SSCP_L12_FEniCSx` directory to give the Docker container access to the example scripts.

```
docker run -it --name sscp_2025 -v $(pwd):/home dolfinx/dolfinx:stable
```
