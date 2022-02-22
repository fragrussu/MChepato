# MChepato
This repository contains the code and the synthetic data used for: **"Diffusion MRI signal cumulants and hepatocyte microstructure at fixed diffusion time: Insights from simulations, 9.4T imaging, and histology"**. *Grussu F, Bernatowicz K, Casanova-Salas I, Castro N, Nuciforo P, Mateo J, Barba I, Perez-Lopez R*; Magnetic Resonance in Medicine 2022 (early view ahead of print), [doi: 10.1002/mrm.29174](https://doi.org/10.1002/mrm.29174).

**FG has received funding from the postdoctoral fellowships programme Beatriu de Pinós (2020 BP 00117), funded by the Secretary of Universities and Research (Government of Catalonia).**

![GenCatFund](https://github.com/fragrussu/MChepato/blob/main/funder.png)

# Dependencies
This repository contains python and bash shell scripts. The code relies on the following dependencies:

* a python 3 distribution with:
   * [Matplotlib](https://matplotlib.org/stable/index.html)
   * [Numpy](https://numpy.org)
   * [statsmodels](https://www.statsmodels.org/stable/index.html)
   * [scikit-learn](https://scikit-learn.org/stable/)
   * [Scipy](https://www.scipy.org/)
   * [DiPy](https://dipy.org/)
   * [MP-PCA python code](https://github.com/NYU-DiffusionMRI/mppca_denoise/blob/master/mpdenoise.py) by New York University
   * [Gibbs unringing python code](https://github.com/RafaelNH/gibbs-removal/blob/master/gibbs_removal.py) by Henriques R
   * [MRItools](https://github.com/fragrussu/MRItools)
   
* [MCDC diffusion simulator](https://github.com/jonhrafe/MCDC_Simulator_public) by Patiño R et al
* [FMRIB Software Library](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
* [NiftyReg](http://cmictig.cs.ucl.ac.uk/wiki/index.php/NiftyReg)
* [QuPath](https://qupath.github.io/)

# Content
This repository contains three sub-folders, each with its own README file within it. The sub-folders are:

* [`simulations`](https://github.com/fragrussu/MChepato/tree/main/simulations): folder containing the code used to perform simulations. README file [here](https://github.com/fragrussu/MChepato/blob/main/simulations/README.md);
* [`data`](https://github.com/fragrussu/MChepato/tree/main/data): folder containing syntehtic data processed by the code stored in [`simulations`](https://github.com/fragrussu/MChepato/tree/main/simulations). README file [here](https://github.com/fragrussu/MChepato/blob/main/data/README.md);
* [`exvivo`](https://github.com/fragrussu/MChepato/tree/main/exvivo): folder containing the code written to analyse 9.4T _ex vivo_ MRI scans of fixed mouse livers and their co-localised histological images. README file [here](https://github.com/fragrussu/MChepato/blob/main/exvivo/README.md).

The folder [`perturbed`](https://github.com/fragrussu/MChepato/tree/main/data/perturbed) within [`data`](https://github.com/fragrussu/MChepato/tree/main/data) contains synthetic cell meshes in .ply format used to perform Monte Carlo simulations. The meshes are illustrated below (15 cell shapes for each unique cell size).

<img src="https://github.com/fragrussu/MChepato/blob/main/cellmesh.png" width="700"> 

Researchers interested in accessing the _ex vivo_ mouse liver MRI and histology data can contact FG at [<fgrussu@vhio.net>](mailto:fgrussu@vhio.net) to stipulate relevant research and data transfer agreements.

# License
This repository is distributed under the Attribution-ShareAlike 4.0 International license ([CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)). Copyright (c) 2022, Fundació Privada Institut d’Investigació Oncològica de Vall d’Hebron. All rights reserved. Link to license [here](https://github.com/fragrussu/MChepato/blob/main/LICENSE.txt). 

**The use of MChepato MUST also comply with the individual licenses of all of its dependencies.**

# Funding
**FG has received funding from the postdoctoral fellowships programme Beatriu de Pinós (2020 BP 00117), funded by the Secretary of Universities and Research (Government of Catalonia).** This project was also supported by Fundació La Caixa and by the investigator-initiated PREdICT study at the Vall d'Hebron Institute of Oncology (Barcelona), funded by AstraZeneca and supporting FG. KB is funded by a Beatriu de Pinós post-doctoral grant (2019 BP 00182). RPL is supported by a CRIS Foundation Talent Award (TALENT19-05), the Instituto de Salud Carlos III-Investigación en Salud (PI18/01395), Spanish Ministry for Science, Innovation and Universities (RTI2018-095209-B-C21, FIS-G64384969), Prostate Cancer Foundation Young Investigator Award and Fero Foundation. ICS is supported by a fellowship from Fundació ”la Caixa” (ID 100010434) and the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 847648, fellowship code LCF/BQ/PI20/1176003.

# Disclosures
FG was supported by PREdICT, a study funded by AstraZeneca at the Vall d’Hebron Institute of Oncology (Barcelona, Spain). There are no conflicts of interest: AstraZeneca was not involved in any aspect concerning this study; it has not influenced the analysis of the data, the interpretation of the results, and the decision to submit the manuscript for publication.
