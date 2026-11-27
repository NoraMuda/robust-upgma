# Data

Copies of the published data sets used in the paper, kept here so that the results can be
reproduced offline. They are **not** covered by the MIT licence of the code.

## `fcps/` - Fundamental Clustering Problem Suite

Atom, Chainlink, EngyTime, Hepta, Lsun3D, Target, Tetra, TwoDiamonds and WingNut (benchmarks)
and GolfBall (data without cluster structure, Table 7), as `.rda` files from the FCPS R package,
<https://github.com/Mthrun/FCPS> (GPL-3).

Thrun, M. C. & Ultsch, A. Clustering benchmark datasets exploiting the fundamental clustering
problems. *Data in Brief* 30, 105501 (2020).

## `pmlb/` - labelled real data

Iris, Wine, Breast Cancer (Wisconsin diagnostic), E. coli, Penguins and MAGIC gamma telescope, as
distributed by the Penn Machine Learning Benchmarks collection, <https://github.com/EpistasisLab/pmlb>,
which redistributes the UCI Machine Learning Repository versions (see each data set's original
licence there).

Romano, J. D. et al. PMLB v1.0: an open-source dataset collection for benchmarking machine
learning methods. *Bioinformatics* 38, 878-880 (2022).

Preprocessing (in `python/experiments/common.py::load_dataset` and `R/data_loader.R`): features
z-scored, constant features removed, Penguins restricted to the four morphometric variables
(complete cases). EngyTime and MAGIC are class-stratified subsampled to 2,000 observations in
every benchmark realisation.
