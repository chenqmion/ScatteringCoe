# Scattering Coe(ficients

## Related works:
* Scattering coefficients of superconducting microwave resonators: I. Transfer-matrix approacch [[Phys. Rev. B 106, 214505 (2022)]](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.214505)
* Scattering coefficients of superconducting microwave resonators: II. System-bath approach [[Phys. Rev. B 106, 214506 (2022)]](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.214506)

## Description:
This repository contains all the codes for generating and analyzing the data of the works mentioned above. Given a circuit, it generates the scattering coefficients. Given the scattering coefficients, it corrects distortions and provides a full characterization of the resonator. 

* **Toolkit_Simulation.py** contains both the numerical and analytical methods for calculating the scattering coefficients.

* **HangerType.py**, **NecklaceType.py**, and **CrossType.py** contain examples for three types of resonators. The corresponding results are **Fig_hanger.pdf**, **Fig_necklace.pdf**, **Fig_cross.pdf**.

* **HangerType_Chain.py** and **NecklaceType_Chain.py** contain examples for two types of resonator crystals. The corresponding results are **Fig_hanger_chain.pdf** and **Fig_necklace_chain.pdf**.

* **Toolkit_Characterization.py** contains the step by step recipe for characterizing a general resonator. Here, we assume that we are measuring the reflection response of a symmetric necklace-type &lambda;/2 resonator. To run the code, one needs to load the measured complex-valued scattering response to the array "s" with the corresponding measurement frequencies to the array "f".  

* **Data_20190926_15mmgap_f3.npy** is the experimental data with S11 measurement of a necklace-type &lambda;/2 resonator.

* **Data_necklace_lamb2_asymmetry.npy** is fake data of a necklace-type &lambda;/2 resonator, which is numerically generated for crosschecking the results. The characterization result is **Fig_characterization.pdf**.

## Acknowledgement:
We are delighted if you find this project helpful to your own study. Feel free to contact us if you have questions, suggestions, criticisms, etc. You are free to copy, share, and build on this project without notifying the authors. Citing our publications is appreciated but not required. 
