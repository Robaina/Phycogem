# Phycogem: Reconstructing phycospheric metabolic networks

Updating CarveME's universe curation pipeline in [this fork](https://github.com/Robaina/carveme_expanded_universe)


# :bulb: About

This repo contains Python code and Jupyter Notebooks aimed at facilitating the reconstruction and analysis of community genome-scale models reconstructed from genome data.

## :wrench: Installation

To install the package, first clone the repo and create the conda environment. We recomend using [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) instead of conda to speed up the process:

```bash
git clone https://github.com/Robaina/Phycogem.git
cd Phycogem
mamba env create -f envs/phycogem-dev.yml
conda activate phycogem-dev
```

Then build and install the package:

```bash
poetry build
pip install dist/phycogem*.whl
```

## :rocket: Usage

To use the package, activate the Conda environment:

```bash
conda activate phycogem-dev
```

This environment includes a Jupyter Notebook ipykernel, so the package can be used in a notebook as well.

## :notebook_with_decorative_cover: Notebooks

This repo contains several Notebooks to exemply the use of the package. They can be found in the `notebooks` folder.