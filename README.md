# PBD

[![Build Status](https://travis-ci.org/richelbilderbeek/PBD.svg?branch=master)](https://travis-ci.org/richelbilderbeek/PBD)
[![codecov.io](https://codecov.io/github/richelbilderbeek/PBD/coverage.svg?branch=master)](https://codecov.io/github/richelbilderbeek/PBD?branch=master)

Protracted Birth-Death model in R

This is a development version before the official release on CRAN.

## Installing PBD

The PBD package has a stable version on CRAN and 
a development version on GitHub.

### From CRAN

From within R, do:

```
install.packages("PBD")
```

### From GitHub

Because the PBD package is located in the folder `PBD`, do:

```
devtools::install_github("rsetienne/PBD/PBD")
```

## Using PBD as a package dependency

### From CRAN

To your DESCRIPTION file, add `PBD` as any normal package.

If your package directly uses `PBD`:

```
Imports:
  PBD
```

If your package uses `PBD` in its perepherals (e.g. vignettes and tests):

```
Suggests:
  PBD
```

### From GitHub

Because the PBD package is located in the folder `PBD`, do:

```
Removes:
  rsetienne/PBD/PBD
```
