# SenSAgri Classification Chain
The SenSAgri Classification Chain performs supervised classification based on Sentinel-1 and Sentinel-2 time series
to output a Crop Mask and a Crop Type products each provided with their respective confidence map.

## Dependencies
The chain needs the following dependancies that are commonly installed on Unix/GNU-Linux systems designs for scientific calculations:
- gcc 4.8.5
- cmake 3.9.1
- python 2.7.5 with numpy, scypy, matplotlib and ogr/gdal libraries.
- mpirun 3.0.4
- pdflatex

Slighlty lower versions of these dependencies might still work but it has not been tested.

## Installing
Details about installation and usage can be found on the Delivrable D5.4 of the SenSAgri project.
(see documents/D5.4_Software_of_the_seasonal_crop_type_mapping_prototype_OBSOLETE.pdf )
Note: This document is very obsolete.

## Versioning
Version 0.1 (2019)

## Authors
* **CesBIO** - *Initial work* - [sensagricesbio](https://github.com/sensagricesbio)

## License
This is free software under the GNU Affero General Public License v3.0. The licence is shared between:
- Université Paul Sabatier - Toulouse III
- CesBIO
- SenSAgri project officier

See http://www.gnu.org/licenses/agpl.html for more details.
