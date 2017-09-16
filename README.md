# scaleboot

What is scaleboot?

scaleboot is an add-on package of R. This is for calculating
approximately unbiased (AU) p-values from a set of multiscale
bootstrap probabilities for a hypothesis. Scaling is equivalent to
changing the sample size of data set in bootstrap resampling. We
calculate bootstrap probabilities at several scales, from which a very
accurate p-value is calculated. This multiscale bootstrap method has
been implemented in CONSEL software and pvclust package. The thrust of
scaleboot package is to calculate an improved version of AU p-values
which are justified even for hypotheses with nonsmooth boundaries by
taking care of the singularity.

scaleboot package includes an interface to pvclust package of R for
bootstrapping hierarchical clustering. We use pvclust to calculate
multiscale bootstrap probabilities, from which we calculate an
improved version of AU p-values using scaleboot.

scaleboot has a front end for phylogenetic inference, and it can
replace CONSEL software for testing phylogenetic trees. Currently,
scaleboot does not have a method for file conversion of several
phylogenetic software, and so we must use CONSEL for this purpose
before applying scaleboot to calculate an improved version of AU
p-values for trees and edges.

The package vignette "Multiscale Bootstrap Using Scaleboot Package"
(usesb.pdf) explains the methodology. It includes a simple example for
illustration. It also includes real applications in hierarchical
clustering and phylogenetic inference. Further description is given in
Shimodaira (2008). For the use of scaleboot, Shimodaira (2008) may be
referenced.

The official website of scaleboot is
http://stat.sys.i.kyoto-u.ac.jp/prog/scaleboot/

The CRAN site is
http://cran.r-project.org/web/packages/scaleboot/

The package is now developped at the github site
http://github.com/shimo-lab/scaleboot
