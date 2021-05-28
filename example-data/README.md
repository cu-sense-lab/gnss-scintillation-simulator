## Description of the data

The data was obtained from Septentrio PolaRxS Pro GNSS receiver and converted to MAT-file.

### Variables

```
LOG/
  -PRN
  -GPSTime
  -SatElevation
  -SatAzimuth
  -{L1, L2, L5} PSR (1-Hz)
  -{L1, L2, L5} ADR (100-Hz)
  -{L1, L2, L5} PWR (50-Hz)
  -{L1, L2, L5} CN0 (1-Hz)

eph/
  - (GPS ephemeris parameters)

origin_llh/
  - (longitude, latitude, height of the receiver)

```

## Online Repository link

* [DataRepository](https://github.com/cu-sense-lab/gnss-scintillation-simulator/tree/master/example-data) - Link to the data repository.


## Papers

This Hong Kong example dataset is used in the following publications:

* Sun, K., Chang, H., Lee, J., Seo, J., Morton, Y. J., & Pullen, S. (2020). *Performance Benefit from Dual-Frequency GNSS-based Aviation Applications under Ionospheric Scintillation: A New Approach to Fading Process Modeling*. Paper presented at 2020 International Technical Meeting of The Institute of Navigation, San Diego, CA. https://doi.org/10.33012/2020.17184
* Sun, A. K., Chang, H., Pullen, S., Kil, H., Seo, J., Morton, Y. J., & Lee, J. (2021). Markov Chain-based Stochastic Modeling of Deep Signal Fading: Availability Assessment of Dual-frequency GNSS-based Aviation under Ionospheric Scintillation. *Space Weather*, (Submitted)


## Acknowledgments

* The data was collected by a system built at the Satellite Navigation and Sensing Lab at the University of Colorado Boulder

* The receiver that collected the Hong Kong dataset is hosted by Prof. Zhizhao Liu from Hong Kong Polytechnic University.
