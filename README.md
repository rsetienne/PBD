# PBD

Branch|[![Travis CI logo](TravisCI.png)](https://travis-ci.org)|[![Codecov logo](Codecov.png)](https://www.codecov.io)
---|---|---
master|[![Build Status](https://travis-ci.org/rsetienne/PBD.svg?branch=master)](https://travis-ci.org/rsetienne/PBD) | [![codecov.io](https://codecov.io/github/rsetienne/PBD/coverage.svg?branch=master)](https://codecov.io/github/rsetienne/PBD?branch=master)
master|[![Build Status](https://travis-ci.org/richelbilderbeek/PBD.svg?branch=master)](https://travis-ci.org/richelbilderbeek/PBD) | [![codecov.io](https://codecov.io/github/richelbilderbeek/PBD/coverage.svg?branch=master)](https://codecov.io/github/richelbilderbeek/PBD?branch=master)
develop|[![Build Status](https://travis-ci.org/richelbilderbeek/PBD.svg?branch=develop)](https://travis-ci.org/richelbilderbeek/PBD) | [![codecov.io](https://codecov.io/github/richelbilderbeek/PBD/coverage.svg?branch=develop)](https://codecov.io/github/richelbilderbeek/PBD?branch=develop)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/PBD)](https://cran.r-project.org/package=PBD)
[![](http://cranlogs.r-pkg.org/badges/grand-total/PBD)]( https://CRAN.R-project.org/package=PBD)
[![](http://cranlogs.r-pkg.org/badges/PBD)](https://CRAN.R-project.org/package=PBD)

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
Remotes:
  rsetienne/PBD/PBD
```
