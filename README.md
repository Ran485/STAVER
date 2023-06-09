<p align="center">
    <br>
    <img src="https://github.com/Ran485/STAVER/blob/main/docs/_static/STAVER_logo.svg" width="400"/>
    <br>
    <h2 align="center">
    STAVER: A Standardized Dataset-Based Algorithm for Efficient Variation Reduction
    </h2>
<p>


<div align="center">
  <a href="#">
  <a href="https://github.com/Ran485/STAVER/stargazers">
  <img alt="Downloads" src="https://img.shields.io/github/stars/Ran485/STAVER?logo=GitHub&color=red">
  </a>
  <img src="https://img.shields.io/badge/Python-3.7+-blue">
  </a>
  <a href="https://github.com/dwyl/esta/issues">
  <img alt="PRs welcome" src="https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat">
  </a>
  <a href="https://opensource.org/licenses/MIT">
  <img alt="DOI" src="https://img.shields.io/badge/License-MIT-yellow.svg">
  </a>
</div>

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Getting Started](#getting-started)
- [Documentation](#documentation)
- [How to Contribute](#how-to-contribute)
- [Contact Us](#contact-us)
- [License](#license)


## Introduction

STAVER is Python library that presents a standardized dataset-based algorithm designed to reduce variation in large-scale data-independent acquisition (DIA) mass spectrometry data. By employing a reference dataset to standardize mass spectrometry signals, STAVER effectively reduces noise and enhances protein quantification accuracy, especially in the context of multi-library search. The effectiveness of STAVER is demonstrated in several large-scale DIA datasets, showing improved identification and quantification of thousands of proteins. STAVER, featuring a modular design, provides flexible compatibility with existing DIA MS data analysis pipelines. The project aims to promote the adoption of multi-library search and improve the quality of DIA proteomics data through the open-source STAVER software package. A comprehensive overview of the research workflow and STAVER algorithm architecture are summarized in the following figure: ![alt text](https://github.com/Ran485/STAVER/blob/main/docs/_static/STAVER_pipeline.png)

## Installation

You can install ``staver`` package from [PyPI](https://pypi.org/project/dia-staver/) by calling the following command: 

``` shell
pip install dia-staver
```
You may install from source by cloning the STAVER repo, navigating to the root directory and using one of the following commands ``pip install .``, or ``pip install -e .`` to install in editable mode:

``` shell
# clone the source repo
git clone https://github.com/Ran485/STAVER.git

# install the package in editable mode
pip install .

# or using the following command
pip install -e .
```
You may install additional environmental dependencies:

``` shell
pip install -r requirements_dev.txt
pip install -r requirements.txt
```

## Getting Started

For example code and an introduction to the library, see the Jupyter notebooks in
[tutorials](https://pypi.org/project/dia-staver/), and the guided walkthrough
[here](https://pypi.org/project/dia-staver/). A straightforward command-line demonstration for a quick start can be discovered in the following block.

```shell
python  ./staver_pipeline.py \
        --thread_numbers < The CPU worker numbers, Default to [nmax-2] > \
        --input < The DIA data input directory > \
        --output_peptide < The processed DIA peptide data output directory > \
        --output_protein < The processed DIA protein data output directory > \
        --count_cutoff_same_libs < Default to 1 > \
        --count_cutoff_diff_libs < Default to 2 > \
        --proteins_cv_thresh < Default to 0.3 > \
        --na_threshold < Default to 0.3 > \
        --top_precursor_ions < Default to 3 > \
        --file_suffix < Default to "_F1_R1" >  \
```
Run the `test-data` in the following block
```shell
python  ./staver/staver_pipeline.py \
        --thread_numbers 16 \
        --input ./staver/data/ \
        --reference_dataset_path ./data/ \
        --output_peptide ./staver/results/peptides/ \
        --output_protein ./staver/results/proteins/ \
        --count_cutoff_same_libs 1 \
        --count_cutoff_diff_libs 2 \
        --peptides_cv_thresh 0.3 \
        --proteins_cv_thresh 0.3 \
        --na_threshold 0.3 \
        --top_precursor_ions 5 \
        --file_suffix _F1_R1 \
```

## Documentation
To gain a comprehensive understanding of STAVER's application and to thoroughly appreciate the function and purpose of each parameter, we highly recommend perusing the all-encompassing STAVER [documentation](https://pypi.org/project/dia-staver/). This resource provides detailed, step-by-step instructions, accompanied by illustrative examples and clear explanations, equipping users with the knowledge to skillfully navigate and exploit the software's complete potential.

## How to Contribute
We welcome the contribution from the open-source community to improve the library!

To add a new explanation method/feature into the library, please follow the template and steps demonstrated in this 
[documentation](https://pypi.org/project/dia-staver/).

## Contact Us
If you have any questions, comments or suggestions, please do not hesitate to contact us at 21112030023@m.fudan.edu.cn

## License
The STAVER project licensed under the [MIT License](https://opensource.org/licenses/MIT), granting users open access and the freedom to employ, adapt, and share the software as needed, while preserving the original copyright and license acknowledgements.
