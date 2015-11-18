consensusSeekeR : Detection of consensus regions inside a group of experiments using genomic positions and genomic ranges
=====================

[![Build Status](https://travis-ci.org/ArnaudDroitLab/consensusSeekeR.svg?branch=master)](https://travis-ci.org/ArnaudDroitLab/consensusSeekeR)
[![codecov.io](https://codecov.io/github/ArnaudDroitLab/consensusSeekeR/coverage.svg?branch=master)](https://codecov.io/github/ArnaudDroitLab/consensusSeekeR?branch=master)

This R package compares multiple narrowPeak data from different experiments to extract common peak regions. 
The size of the analyzed region is adjustable, as well
as the number of experiences in which a peak must be present to tag a 
potential region as a consensus region. If needed, the consensus regions can be extended to cover the entire regions of enclosed peaks.

## Bioconductor Package ##

[![Bioconductor Time](http://bioconductor.org/shields/years-in-bioc/consensusSeekeR.svg)](http://bioconductor.org/packages/devel/bioc/html/consensusSeekeR.html "Bioconductor status")

consensusSeekeR is now an official package of [Bioconductor](http://bioconductor.org/). The current development release can be directly downloaded from their website:
[Current devel release](http://bioconductor.org/packages/devel/bioc/html/consensusSeekeR.html )


## Authors ##

[Astrid Desch&ecirc;nes](http://ca.linkedin.com/in/astriddeschenes 
"Astrid Desch&ecirc;nes"), 
[Fabien Claude Lamaze](http://ca.linkedin.com/in/fabienlamaze/en 
"Fabien Claude Lamaze"), 
[Pascal Belleau](http://ca.linkedin.com/in/pascalbelleau 
"Pascal Belleau") 
and [Arnaud Droit](http://ca.linkedin.com/in/drarnaud 
"Arnaud Droit").

See [Arnaud Droit Lab](http://bioinformatique.ulaval.ca/home/ 
"Arnaud Droit Lab") website.

## Maintainer ##

[Astrid Desch&ecirc;nes](http://ca.linkedin.com/in/astriddeschenes 
"Astrid Desch&ecirc;nes")

## License ##

This package and the underlying consensusSeekeR code are distributed under the 
Artistic license 2.0. You are free to use and redistribute this software. 

For more information on Artistic 2.0 License see
[http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)

## Bugs/Feature requests ##

If you have any bugs or feature requests,
[let us know](https://github.com/ArnaudDroitLab/consensusSeekeR/issues). Thanks!
