d:
cd "d:\data\ms\PBD"
R CMD build PBD
R CMD INSTALL --build PBD_1.3.tar.gz
pause
R CMD check --timings --as-cran PBD_1.3.tar.gz
pause