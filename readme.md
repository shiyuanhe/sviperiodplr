# Stochastic variational inference of periods and period-luminosity relations for Mira variable stars


## Data

There are three datasets used in the paper. These datasets are packed in the folder `./data/`.

- The first dataset consists of the 5000 synthetic Mira O-rich samples from Yuan et al. (2018). This dataset is used in our first simulation comparision in Section 4. 
- The second dataset is the actuall M33 candidate Mira observations  from Yuan et al. (2018). This dataset is for the real data analysis in Section 5.
- Our third dataset corresponds to the simulation in the Supplementary Materials. It consists of 90 sub-datasets. Each has 1,000 Mira observations. 


The table in the text file `./data/id_map_o.dat` has eight columns. It contains the ground truth for the first simulation dataset (Yuan et al. 2018) in Section 4.

| Column Name | Description |
|-------|--------------|
| simID | Simulation Sample ID |
| ogleID | The OGLE sample ID. It indicates which OGLE observation is used as template to generate this synthetic observation.|
| m33ID | The M33 observation ID. It indicates which M33 observation cadence is used to generate this synthetic observation. |
| type | The subtype of Miras. All of the synthetic observations are O-rich Miras. |
| period| The true period of the synthetic observation. |
| mJ | The average J band magnitude. |
| mH | The average H band magnitude. |
| mK | The average K band magnitude. |

The table in the text file `./data/simu_PtP1P2.dat` is the estimation result of Yuan et al. (2018). It has six columns. Our paper will employ four columns as decribed below. In particular, the third column named P1 is used for accuracy comparision in our simulation study.

| Column Name | Description |
|-------|--------------|
| simID | Simulation Sample ID |
| type | The subtype of Miras. All of the synthetic observations are O-rich Miras. |
| P1 | The period of the primary peak of MSP |
| P2 | The period of the secondary peak of MSP |


The folder `./data/simu1Raw/` contains the 5,000 files. Each file corresponds to one synthetic Mira observation in Section 4. Each file is also named by its simulation sample ID. According to the ID in the file name, its ground truth information can be found in the table  `./data/id_map_o.dat`. Each file in the folder `./data/simu1Raw/` has six columns, as described below.

| Column Number | Description |
|-------|--------------|
| 1 | The time (a quantity related to the Julian date) of the observational light curve point. |
| 2 | The band of the observed light curve point. |
| 3 | The observation survey. |
| 4 |  The magnitude of the observed light curve point.  |
| 5 |  The photometric uncertainty of the observed light curve point.  |
| 6 |  The transformation uncertainty of the observed light curve point.  |


The text file `./data/yuan2018.dat` is from Yuan et al. (2018). It contains the information for the M33 candidate Miras. It has 26 columns. Our paper uses several columns as described below.

| Column Name | Description |
|-------|--------------|
| id | The M33 Observation ID |
| P | The estimated period by Yuan et al. (2018).|
| m33ID | The M33 observation ID. It indicates which M33 observation cadence is used to generate this synthetic observation. |
| period| The true period of the synthetic observation. |
| J | The estimated average J band magnitude by Yuan et al. (2018). |
| H | The estimated average H band magnitude by Yuan et al. (2018). |
| K | The estimated  average K band magnitude by Yuan et al. (2018). |
| cl | The classified subtype of Mira. The symbols N, O, C indicate Unknown, O-rich and C-rich, respectively. |

The folder `./data/elcs3` contains the actual observation for M33 candidate Miras. Each file corresponds to one Mira observation, and has six columns as described below.

| Column Number | Description |
|-------|--------------|
| 1 | The time (a quantity related to the Julian date) of the observational light curve point. |
| 2 | The band of the observed light curve point. |
| 3 | The observation survey. |
| 4 |  The magnitude of the observed light curve point.  |
| 5 |  The photometric uncertainty of the observed light curve point.  |
| 6 |  The transformation uncertainty of the observed light curve point.  |

The file `./data/simu2/output.tar.gz` is the compressed datasets for the second simulation in the Supplementary Materials.  Exctract the file to get the folder `./data/simu2/output/`. It has three subfolders: `pttn1` `pttn2` and `pttn3`, corresponding to different sampling cadence. In each of them also exits three subfolders: `noise1` `noise2` and `noise3`, corresponding to various noise levels. Also, each of them contains another level of subfolders `size1`, `size2`, ... `size10`, corresponding to distinct number of sampling points. Suppose $n_K$ is the number of points for the K band light curves, and $n_I$ is the number of points for the I band light curves. The deepest level folder name maps to the number of points as the table below. In addition, at the deepest folder level, there is a folder called `pars`. It contains the ground truth for the synthetic dataset.


| Folder | Description |
|-------|--------------|
| size1 | $(n_K, n_I) = (5,5)$|
| size2 | $(n_K, n_I) = (5,10)$|
| size3 | $(n_K, n_I) = (5,20)$|
| size4 | $(n_K, n_I) = (5,30)$|
| size5 | $(n_K, n_I) = (10,10)$|
| size6 | $(n_K, n_I) = (10,20)$|
| size7 | $(n_K, n_I) = (10,30)$|
| size8 | $(n_K, n_I) = (20,20)$|
| size9 | $(n_K, n_I) = (20,30)$|
| size10 | $(n_K, n_I) = (30,30)$|
| pars | Each file contains the true frequency, the true average K band and I band magnitude. |




## Code

In the simulation experiments, we have compared our proposed method with the generalized Lomb-Scargle method (GLS, Zechmeister and Kurster, 2009) and its multi-band version (MGLS, Vanderplas andIvezic, 2015), the single band semi-parametric model (SP, He et al., 2016a) and its multi-band extension (MSP, Yuan et al., 2018). The codes contain the R and bash scripts to reproduce the simulation results in Section 4, the real data analysis results in Section 5, and the additional simulation resutls in the Supplementary Materials. The large synthetic dataset can also be reproduced by the code for the second simulation in the Supplementary Material.


The code depends on the following R and C++ packages. 

- microbenchmark (Version 1.4.7 [https://cran.r-project.org/package=microbenchmark](https://cran.r-project.org/package=microbenchmark))
- snowfall  (Version 1.84.6.1, [https://cran.r-project.org/package=snowfall](https://cran.r-project.org/package=snowfall)).
- readr (Version 1.3.1, [https://CRAN.R-project.org/package=readr](https://CRAN.R-project.org/package=readr)).
- scales (Version 1.1.0, [https://cran.r-project.org/package=scales](https://cran.r-project.org/package=scales)).
- tidyverse (Version 1.3.0, [https://CRAN.R-project.org/package=tidyverse](https://CRAN.R-project.org/package=tidyverse)).
- Rcpp (Version 1.0.3, [https://CRAN.R-project.org/package=Rcpp](https://CRAN.R-project.org/package=Rcpp)).
- RcppArmadillo (Version 0.9.850.1.0, [https://CRAN.R-project.org/package=RcppArmadillo](https://CRAN.R-project.org/package=RcppArmadillo)).
- RcppProgress (Version 0.4.1, [https://CRAN.R-project.org/package=RcppProgress](https://CRAN.R-project.org/package=RcppProgress)).
- multiband (Version 0.1.0, [https://cran.r-project.org/package=multiband](https://cran.r-project.org/package=multiband)).
- the varStar package available at GitHub repo ([https://github.com/shiyuanhe/varStar](https://github.com/shiyuanhe/varStar)).
- The C++ package armadillo (Version 9.870.2, [http://arma.sourceforge.net/](http://arma.sourceforge.net/)).
- The TRNG c++ library(Version 4.19. [https://www.numbercrunch.de/trng/](https://www.numbercrunch.de/trng/).


To run the multi-band semi-parametric model (MSP, Yuan et al., 2018), compile the C++ file `./MSP/mgp.cpp` by the following command. It will generate an executable called `mGP` in the folder `./MSP`. Before compiling, make sure the C++ package ([armadillo](http://arma.sourceforge.net/)) has been installed. 

```
clang++ ./MSP/mgp.cpp -o ./MSP/mGP -Wall -O2 -std=c++11 -larmadillo
```

See the documentation for the installation and compilation details of the  C++ package ([armadillo](http://arma.sourceforge.net/)). The source file can also be compiled without the runtime library. But the user need to provide the path to the armadillo headers.

```
clang++ ./MSP/mgp.cpp -o ./MSP/mGP -Wall -O2 -std=c++11 -I/PATH/TO/ARMADILLO/include
```

Our proposed method is packed in the R package `./variationalMira_0.1.0.tar.gz`. Install the package by the following command. Before installing the package, make sure its dependent packages have been installed.

```
R CMD INSTALL variationalMira_0.1.0.tar.gz
```

## Steps to Reproduce


The reproducing procedures are divided into two steps, as some of the compared methods take a long time to finish. In the first step, each individual method is run. Their estimation results are stored as R data files (“.RData”). Based on these cached estimation result files, an R markdown file will be able to generate the figures and tables in the paper and the Supplmentary Material. 

- The first step is to convert the file format of Yuan et al. (2018) to the format required by our code. The conversion is done by the following command 

```bash
Rscript ./code/simu1/simu1_prepare.R
Rscript ./code/realdata/RealDataPrepare.R
```

- The estimation results of the first simulation can be reproduced by running the following commands from the terminal. The computed results will be saved into the folder `./result/simu1`. For the MSP methods on the first simulation, we directly use the estimation result in the text file `./data/simu_PtP1P2.dat` provided by the authors of Yuan et al. (2018).

```bash
 Rscript ./code/simu1/simu1SP.R
 Rscript ./code/simu1/simu1GLS.R
 Rscript ./code/simu1/simu1MGLS.R
 Rscript ./code/simu1/simu1SVI.R
```

- The estimation results of the second simulation in the Supplemtary Materials can be reproduced by running the following bash scripts from the terminal. The results of all methods will be saved into the folder `./result/simu2`. Before running the code, make sure to uncompress the file `./data/simu2/output.tar.gz` into the folder `./data/simu2/output/`.

```bash
 bash simu2SVI.sh
 bash simu2gls.sh
 bash simu2MGLS.sh
 bash simu2SP.sh
 bash simu2MSP.sh
```

- The results of the real data analysis can be reproduced by running the follwoing command from the linux bash.

```bash
Rscript ./code/realdata/realData.R 
```

- Finally, all the resulting figures and tables in the paper can be generated by running the R markdown file `NumericalResults.Rmd`. It will generate a “NumericalResults.html” file with all the tables and figures. Its computation is based on the cached estimation results in the above steps.


