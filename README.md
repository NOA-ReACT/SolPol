# SolPol
Instrument manual and data processing algorithm for the NOA/UH - Solar Polarimeter (SolPol)

This repository contains everything related to the operation, characterization and data manipulation of the SolPol - Solar polarimeter, as operated by the NOA-ReACT team. The instrument was conferred to the National Observatory of Athens by the University of Hertfordshire in the framework of the D-TECT/ERC project.

**SolPol**:
- is a passive polarimetric ground-based instrument that targets directly the Sun
- has measuring capabilities of a few parts per million (ppms) in linear polarization and of the order of $10^{-7}$ in circular polarization
- is based on the existing design of the PlanetPol high sensitivity polarimeter as described in Hough et al., 2006. 
- usage under laboratory conditions is described by Martin et al., 2010.


## Index

In this repository you can find:
- `solpol_data.py`: a python 3.7-based code for processing the raw data retrieved from the instrument
- `drivers`: the drivers for each peripheral component required for the operation of SolPol
- `manual`: a complete guide for the instrument set-up, initialization and measurement sequence

## Acknowledgements
This research was supported by D-TECT (Grant Agreement 725698) funded by the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme.

## Blame

- **Team lead**: Lilly Daskalopoulou <vdaskalop@noa.gr> @ModusElectrificus
- **Research consultant**: Alexandra Tsekeri <atsekeri@noa.gr>
- **Instrument operator**: Panagiotis-Ioannis Raptis <piraptis@meteo.noa.gr>
- **Instrument PI**: Vassilis Amiridis <vamoir@noa.gr>
- **Consulting**: Thanasis Georgiou <ageorgiou@noa.gr> @thgeorgiou


Period of operation: extensive measurements and tests between August 2018 to November 2022.
Operation locations: i. the PANGEA observatory in Antikythera (primarily), ii. the National Observatory of Athens (partially), iii. the Cyprus Institute, Nicosia, Cyprus (during 2019 campaign) and iv. the OSCM Institute, Mindelo, Cabo Verde (during ASKOS 2022 campaign).
