# TephraFits

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3590642.svg)](https://doi.org/10.5281/zenodo.3590642)

*TephraFits* is a Matlab function for describing the geometry of tephra fallout deposits. Using field data, it allows to:
- Fit contour maps using using exponential, power-law and Weibull methods (Fierstein and Nathenson 1992; Bonadonna and Houghton 2005; Bonadonna and Costa 2012);
- Calculate the deposit volume using isopach;
- Calculate the deposit mass using isomass;
- Estimate the decay of clast diameter using isopleth;
- Calculate transects as a distance from the vent;
- Classify eruptions using Pyle (1989) and Bonadonna and Costa (2013);
- Perform uncertainty assessment using Monte-Carlo simulations (Biass et al. 2014).

The associated paper has been accepted in Journal of Applied Volcanology.

## Documentation
For a detail description of the usage of the function, navigate to the root of *TephraFits* in Matlab and type:
~~~
>> help tephraFits
~~~
Additionally, the file `demo.m` contains examples and additional updates are described at  [https://e5k.github.io](https://e5k.github.io).

## Citation
*TephraFits* was published as a Methodology paper in *Journal of Applied Volcanology* available [here](https://www.researchgate.net/publication/330399543_A_step-by-step_evaluation_of_empirical_methods_to_quantify_eruption_source_parameters_from_tephra-fall_deposits) and [here](https://appliedvolc.biomedcentral.com/articles/10.1186/s13617-018-0081-1). Please cite as:
> Biass S, Bonadonna C, Houghton BF (2019) A step-by-step evaluation of empirical methods to quantify eruption source parameters from tephra-fall deposits. J Appl Volcanol 8:1–16

> Biass S, Bonadonna C, Houghton BF (2019). TephraFits. [doi:10.5281/zenodo.3590642](https://zenodo.org/record/3590642#.XgBi1y2p2Vk)

## License
*TephraFits* is released under a GPL3 license, which means that everybody should 
feel free to contribute, comment, suggest and modify the code for as long as any 
new update remains open-source. A copy of the complete license comes with this function.

The GPL3 License covers tephraFits. The authors of dependencies used in this function
retain the copyright to their works. In particular, I am grateful to:
- J. Landsey for the function [bplot](https://au.mathworks.com/matlabcentral/fileexchange/42470-box-and-whiskers-plot-without-statistics-toolbox)
- J. D'Errico for the function [fminsearchbnd](https://au.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon)
- J. Wells for the function [rsquare](https://au.mathworks.com/matlabcentral/fileexchange/34492-r-square-the-coefficient-of-determination)
- A. Horchler for the function [waittext](https://au.mathworks.com/matlabcentral/fileexchange/56424-waittext)

Don't hesitate to contact me by email should you have a suggestion or find a bug.

Hope this code will help!

Copyright Seb Biass (sbiasse AT ntu DOT edu DOT sg) - 2018

## References
- Biass, S., Bagheri, G., Aeberhard, W., Bonadonna, C., 2014. TError: towards a better quantification of the uncertainty propagated during the characterization of tephra deposits. Stat. Volcanol. 1, 1–27. doi:10.5038/2163-338X.1.2
- Bonadonna, C., Costa, A., 2012. Estimating the volume of tephra deposits: A new simple strategy. Geology 40, 415–418.
- Bonadonna, C., Houghton, B., 2005. Total grain-size distribution and volume of tephra-fall deposits. Bull Volcanol 67, 441–456.
- Bonadonna, C., Cioni, R., Pistolesi, M., Connor, C., Scollo, S., Pioli, L., Rosi, M., 2013. Determination of the largest clast sizes of tephra deposits for the characterization of explosive eruptions: a study of the IAVCEI commission on tephra hazard modelling. Bull. Volcanol. 75, 1–15. doi:10.1007/s00445-012-0680-3
- Fierstein, J., Nathenson, M., 1992. Another look at the calculation of fallout tephra volumes. Bull Volcanol 54, 156–167.
- Pyle, D., 1989. The thickness, volume and grainsize of tephra fall deposits. Bull Volcanol 51, 1–15.