# Casboundary

Casboundary is the first method able to automatically define CRISPR cassette boundaries. In addition, it contains a Cas type predictive model used to assign each gene located in the region defined by a cassette’s boundaries a Cas label from a set of pre-defined Cas types. Furthermore, the proposed method can detect potentially new cas genes and decompose a cassette into its modules.

Casboundary can be easily integrated with [CRISPRcasIdentifier](https://github.com/BackofenLab/CRISPRcasIdentifier), a recent tool for the classification of CRISPR cassettes. Casboundary outputs a set of Fasta files containing the identified cassettes, which can be given as input to CRISPRcasIdentifier. By integrating these tools, the users have a complete CRISPR detection and classification pipeline.

## Citation

If you use Casboundary, please cite our paper: [Casboundary: Automated definition of integral Cas cassettes](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa984/5998667) Victor A. Padilha, Omer S. Alkhnbashi, Van Dinh Tran, Shiraz A. Shah, André C. P. L. F. de Carvalho, Rolf Backofen, Bioinformatics, 2020, DOI: [10.1093/bioinformatics/btaa984](https://doi.org/10.1093/bioinformatics/btaa984).

## Installation and requirements

Casboundary has been tested with Python 3.7.7. To run it, we recommend installing the same library versions we used. Since we exported our classifiers following the [model persistence guideline from scikit-learn](https://scikit-learn.org/stable/modules/model_persistence.html), it is not guaranteed that they will work properly if loaded using other Python and/or library versions. For such, we recommend the use of our docker image or a conda virtual environment. They make it easy to install the correct Python and library dependencies without affecting the whole operating system (see below).

### First step: download the last version of the tool and extract it

```
wget https://github.com/BackofenLab/Casboundary/archive/v1.0.0.tar.gz
tar -xzf v1.1.0.tar.gz
```

### Second step: download the Hidden Markov (HMM) and Machine Learning (ML) models

Due to GitHub's file size constraints, we made our HMM and ML models available in Google Drive. You can download them from the following links:

* [Machine Learning Models](https://drive.google.com/file/d/1SkyS03hQG0P7bO7KvaQXgQ_J2Vfng85K/view?usp=sharing)
* [General HMM Models](https://drive.google.com/file/d/1yZ3rl0LPIk-LDLKl2KHb3ra5O2h9_jpA/view?usp=sharing)
* [Signature HMM Models](https://drive.google.com/file/d/1itiqV8djmrfwgwsByZLGe9olv7Uri5hV/view?usp=sharing)
* [Cas HMM Models](https://drive.google.com/file/d/1zRJNlgqAC6A8BEDXPNKmtN4AfnVhjJIR/view?usp=sharing)

Save all tar.gz files inside Casboundary's folder. It is not necessary to extract them, since the tool will do that the first time it is run.

Next, you can choose which third step to follow: either [using a docker container](#third-step-docker) or [using a conda environment](#third-step-conda).

### Third step (docker)

First, you need to install docker (please refer to its [installation guideline](https://docs.docker.com/get-docker/) for details). Next, you can either [pull our image from DockerHub](#pull-the-image-from-dockerhub) or [build the image manually](#build-the-image).

#### Pull the image from DockerHub

```
docker pull padilha/casboundary:1.0.0
```

Inside Casboundary's folder, run the docker image.

```
docker run --rm -v "$(pwd):/home/crispr/Casboundary" -it padilha/casboundary:1.0.0 /bin/bash
```

Since we are using the volume option (-v), Casboundary's folder will be shared between the host machine and the docker container. Thus, there is no need to move files from one to the other.

After this step, everything should be set up and you can skip to [How to use](#how-to-use).

#### Build the image

Inside Casboundary's folder, build the image from the Dockerfile.

```
docker build -t casboundary .
```

Inside Casboundary's folder, run the docker image.

```
docker run --rm -v "$(pwd):/home/crispr/Casboundary" -it casboundary:latest /bin/bash
```

Since we are using the volume option (-v), Casboundary's folder will be shared between the host machine and the docker container. Thus, there is no need to move files from one to the other.

After this step, everything should be set up and you can skip to [How to use](#how-to-use).

### Third step (conda)

Another way to install the correct python version and its dependencies to run Casboundary is by using [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Install Miniconda.

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

Create and activate environment for Casboundary.

```
conda env create -f casboundary-env.yml -n casboundary-env
conda activate casboundary-env
```

After using Casboundary you can deactivate the environment.

```
conda deactivate
```

After this step, everything should be set up and you can skip to [How to use](#how-to-use).

## How to use

To list the available command line arguments type

    python Casboundary.py -h

The available options are:

* `-h` : displays the help message.

* `-f path/to/file.fa` : input DNA fasta file path.

* `-c` : sequence completeness. Available options: `complete` or `partial` (default: `complete`).

* `-o` : output directory path. (default: `./casboundary_output`).

* `-n` : number of CPUs to use for hmmsearch (default: 1).

* `-g` : maximum number of contiguous gaps allowed in a cassette (default: 1).

* `-m` : which ML models to use. Available options: `ERT` or `DNN` (default: `ERT`).

* `-ho` : Hmmsearch output directory path (default: `./hmmsearch_output`).

* `-d` : draws the found cassettes.

## Examples

We provide two simple examples in the `examples` folder:

* `NC_000961.fasta`, which must be run as: `python Casboundary.py -f examples/NC_000961.fasta -c complete`.

* `NC_017104.fasta`, which must be run as: `python Casboundary.py -f examples/NC_017104.fasta -c complete`

* The fasta files generated in the output path in the `predictions` directory can be given as input to [CRISPRcasIdentifier](https://github.com/BackofenLab/CRISPRcasIdentifier). For example:
`python CRISPRcasIdentifier.py -f ./casboundary_output/predictions/cassette_1.fasta`

## License (GPLv3)

    Casboundary
    Copyright (C) 2020 Victor Alexandre Padilha <victorpadilha@usp.br>,
                       Omer Salem Alkhnbashi <alkhanbo@informatik.uni-freiburg.de>,
                       Van Dinh Tran <dinh@informatik.uni-freiburg.de>,
                       Shiraz Ali Shah <shiraz.shah@dbac.dk>,
                       André Carlos Ponce de Leon Ferreira de Carvalho <andre@icmc.usp.br>,
                       Rolf Backofen <backofen@informatik.uni-freiburg.de>
    
    This file is part of Casboundary.
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
