# Scattering Coefficients

## Related works:
* The scattering coefficients of superconducting microwave resonators: I. Transfer-matrix approacch [(arXiv:2109.07762)](https://arxiv.org/abs/2109.07762)
* The scattering coefficients of superconducting microwave resonators: II. System-bath approach [(arXiv:2109.07766)](https://arxiv.org/abs/2109.07766)

## Description:
This repository contains all the codes for generating and analyzing the data of the works mentioned above. Given a circuit, it generates the scattering coefficients. Given the scattering coefficients, it corrects distortions and provides a full characterization of the resonator. 

* **Toolkit_Simulation.py** contains both the numerical and analytical methods for calculating the scattering coefficients.

* **HangerType.py**, **NecklaceType.py**, and **CrossType.py** contain examples for three types of resonators. The corresponding results are **Fig_hanger.pdf**, **Fig_necklace.pdf**, **Fig_cross.pdf**.

* **HangerType_Chain.py** and **NecklaceType_Chain.py** contain examples for two types of resonator crystals. The corresponding results are **Fig_hanger_chain.pdf**, **Fig_necklace_chain.pdf**, **Fig_cross.pdf**.

* **Toolkit_Characterization.py** contains the step by step recipe for characterizing a general resonator. Here, we assume that we are measuring the reflection response of a symmetric necklace-type &lambda;/2 resonator. To run the code, one needs to load the measured complex-valued scattering response to the array "s" with the corresponding measurement frequencies to the array "f".  

* **Data_20190926_15mmgap_f3.npy** is the experimental data with S11 measurement of a necklace-type &lambda;/2 resonator.

* **Data_necklace_lamb2_asymmetry.npy** is fake data of a necklace-type &lambda;/2 resonator, which is numerically generated for crosschecking the results. The characterization result is **Fig_characterization.pdf**.

## Acknowledgement:
We are delighted if you find this projects helpful to your own study. Feel free to contact us if you have questions, suggestions, criticisms, etc. You are free to copy, share, and build on this project without notifying the authors. Citing our publications are appreciated but not required. 
