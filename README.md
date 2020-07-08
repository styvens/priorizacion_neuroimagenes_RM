# Automatic Prioritization of Neuroimaging in Magnetic Resonance


## Objective 
Predict the clinical priority of neuroimaging on magnetic resonance imaging using statistical models that can be implemented in a production environment.


## Introduction
According to the hypothesis of symmetry between brain hemispheres, similarity measures are proposed for the detection of healthy and pathological brain, with neuroimaging. Three measures are explored, the first is a cross-correlation index (NCC), the second is a test based on permutations between hemispheres (HCP) and the third is the measurement of entropy of local values permuted between hemispheres (HPLE) . Probability maps with Bianca were also calculated to use as an input of similarity measures. In addition, the Eigenfaces technique was applied on the NCC map and the probability map, with this method subspaces were obtained with the main characteristics of the image set. Finally, classification models were adjusted for the prediction of the healthy/pathological categories.
 
 
## Results
The HCP was unable to capture enough information from the images to distinguish between healthy and pathological, with low detection of pathological and healthy (61% and 70% respectively). With HPLE, better results were observed compared to HCP, with an f1-score equal to 64% and a recall equal to 78% with a significance level of 5%, however, it is not enough yet to suggest that it can be used in clinical practice unless improvements are made.
 
With the Eigenfaces method, the best results were obtained, especially when working with the NCC correlation map instead of the probability map. The model with the best performance and predictive power was the SVM with 98% predictive capacity.
 
 
## Conclusions
The use of the NCC correlation map obtained from a FLAIR sequence to measure the similarity between hemispheres is useful for the prediction of healthy/pathological images, after generating a subspace with the Eigenfaces method and applying a vector support machine classification model. These methods applied together allow adequate prioritization with 98% predictive capacity.


## Prerequisites

It is required to have the images in nifti format. The processing of the images in this work was done with FSL (it is possible to use others). Python 3.5 was used for data analysis.


## Contributing

The following people have contributed to the project with their valuable knowledge and experience:

* [Gabriel Castrillón Guzmán](https://github.com/gabocas)
* [Henry Laniado Rodas](http://www.eafit.edu.co/docentes-investigadores/Paginas/henry-laniado-rodas.aspx)
* [Catalina Bustamante](https://github.com/catalinabustam)
* [Luisa Sanchez Marin]()
* [Aura Puche](https://github.com/acpuche)


## Authors

* **Danny Cardona Pineda**
