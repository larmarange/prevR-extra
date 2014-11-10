prevR-extra
===========

Some examples of code using the prevR package (http://cran.r-project.org/package=prevR) and some additional functions.

## Code examples

- Infant mortality in Burkina Faso

## Additional functions

### Noptim()

Compute an optimal value of N based on http://cybergeo.revues.org/24606

*Argument:*

- object: a prevR object

### prevR()

A automated function computing a kernel density estimation with bandwiths based on a fixed number of observations (N parameter) and with 

*Arguments:*

`prevR(object, N=Noptim(object), nb.cells=100, cell.size=NULL, weighted=NULL, plot.results=TRUE, return.results=FALSE, legend.title="%", cex=0.7, progression=TRUE)`

- object: a prevR object
- N: N parameter (see ?rings and ?kde). An integer or a list of integer. By default, use an optimal value computed with Noptim().
- nb.cells: number of cells (see ?kde)
- cell.size: size of each cell (see ?kde)
- weighted: use weighted data (see ?kde), by default use weighted data if available in object
- plot.results: plot the obtained map?
- return.results: return computed surfaces as a list of SpatialPixelsDataFrame objects
- legend.title: title of the legend
- cex: size of radius labels
- progression: show progression bars?

### prevR.comp()

This function could be used to compare evolution between two DHSs for the same country, by using the same bandwidths for the two estimated surfaces. Typically, we will have two DHSs conducted in the same country, one in year A and one in year B. A prevR object will be created for each year. The prevalence surface for A will be computed using the usual prevR approach (KDE with adaptive bandwiths of same number of observations). The prevalence surface for B will be computed using the bandwiths computed for A, i.e. the two surfaces will have the same bandwiths.

`prevR.comp(A, B, N = Noptim(A), weighted=NULL, nb.cells=100, cell.size=NULL, plot.results=TRUE, return.results=FALSE, labA="A", labB="B", legend.title="%", cex=0.7, progression=TRUE)`


