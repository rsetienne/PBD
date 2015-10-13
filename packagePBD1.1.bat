d:
cd "d:\data\ms\PBD"
R CMD build PBD
R CMD INSTALL --build PBD_1.1.tar.gz
pause
R CMD check --timings --as-cran PBD_1.1.tar.gz
pause