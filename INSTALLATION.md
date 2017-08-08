# exploBATCH

Explore batch effect (exploBATCH) is a package for discovering and correcting for batch effect using an approach in Nyamundanda et al (2017). Detailed discription on how to run exploBATCH is given in form of a vignette in folder named vignettes.

## INSTALLING exploBATCH.
To install exploBATCH straight from syspremed github account you require devtools package. Otherwise, below are two lines of R codes to install exploBATCH.

require(devtools)

install_github("syspremed/exploBATCH") 

Most problems associated with this installation is due to Rcpp package. Make sure you have up to date gfortran binary for mac, otherwise you can get a “-lgfortran” error or a “-lquadmath” error since few files are missing from their expected locations.  and also update your Rcpp package. For updating gfortran or solving this problem follow this post here
http://www.thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/

