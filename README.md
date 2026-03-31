
# CRISPRhitchhike

The goal of this package is to run the simulations corresponding to the
paper “The benefits of CRISPR–Cas immunity to discriminate good from
better infections” by Domingues et al.

## Installation

You can install the package as follows:

``` r
library("devtools")
install_github("bramkuijper\CRISPRhitchhike")
```

## Usage

Use the function `CRISPRhh()` to run a single simulation, data of which
will be returned in a `data.frame`:

``` r
library("CRISPRhitchhike")
resulting_data <- CRISPRhh()
```
