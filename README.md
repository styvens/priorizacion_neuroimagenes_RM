# Combined mathematical model for the classification of neuroimaging based on measures of similarity between brain hemispheres.


The healthy brain is highly symmetric, therefore, we proposed to use similarity measures in neuroimaging data to distinguish between a healthy and a pathological brain. Three measures were used, the first one Normal Cross Correlation (NCC), the second one is a test based on permutations of values ​​of intensities between hemispheres (HCP) and the third one is a measurement of entropy of values local exchanged between hemispheres (HPLE). Furthermore, we also used the maps derived from BIANCA as an input for the similarity measures

HPLE(with 5% significance)  provided better results compared to the HCP measure both obtained accuracy values ​​equal to 66%, however, the ability to detect pathological brains was greater for HPLE. However, HPLE is not good enough for clinical applications. 

The use of NCC correlation map derived from FLAIR images to measure similarity between hemispheres allowed the classification between healthy and pathological brains by generating a subspace with the Eigenfaces method and applying classification models to the data projected to the proper space found. The reproducibility process included in the generation of own spaces and the adjustment of the classification models, allowed to attribute the results to the method and not to the choice of a specific subspace. Results were observed for this combination of methods with an accuracy of 81% and predictive capacity of 79%.

The main contribution of this project is the combination of similarity measures, methods for the construction of subspaces, and classification models. Specifically, the NCC was used as a measure of similarity, which was projected into the own space of decomposition of singular values following the Eigenfaces methodology, finally, classification models were fitted.

## Prerequisites

It is required to have the images in nifty format. The processing of the images in this work was done with FSL (it is possible to use others). Python 3.5 was used for data analysis.


## Contributing

The following people have contributed to the project with their valuable knowledge and experience:

* [Gabriel Castrillón Guzmán](https://github.com/gabocas)
* [Henry Laniado Rodas](http://www.eafit.edu.co/docentes-investigadores/Paginas/henry-laniado-rodas.aspx)
* [Catalina Bustamante](https://github.com/catalinabustam)
* [Luisa Sanchez]()
* [Aura Puche](https://github.com/acpuche)


## Authors

* [Danny Styvens Cardona Pineda](https://scienti.minciencias.gov.co/cvlac/visualizador/generarCurriculoCv.do?cod_rh=0000079706)