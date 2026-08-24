\## Code supporting the manuscript "**Rain Check for Climate Models**"



Natural Gamma Radiation data used for this paper are freely available through the National Offshore Petroleum Information Management System (NOPIMS, Geoscience Australia) (https://www.ga.gov.au/nopims). The core-based NGR data for the IODP Sites are available through the IODP repository (https://web.iodp.tamu.edu/LORE/?appl=LORE\&reportName=ngr). The scientific downhole logging data is available at https://mlp.ldeo.columbia.edu/logdb/. 



NGR data for all 105 sites is stored in .csv files in two folders: "Sites Data\_Depth-NGR" and "Sites Data\_Age-NGR"



The folder "DTW\_RScripts" in the repository contains the R scripts for the DTW calculations performed and presented for the paper. The separate R script `Custom Step Pattern.R` has the custom step pattern asymmetricP1.1 required for the DTW correlation technique, which must be run first. 



The folder "Time\_Slice\_Map" contains the R scripts `Deep-Site\_NGR.R` and `Shallow-Site\_NGR.R` for plotting Figure 2 and Extended Figure 1, respectively. The other PNG files are required to run the R script `Time-Slice-Map.R`. The required "Sites" and "CENOGRID" data are available as .csv files in the RScripts\&Data folder. The Reef shapefiles for the timeslice maps are provided in the "Reef Shp Files" folder.     



The folder "NetCDF" contains all the NetCDF files for HadCM3 modelled P-E and runoff simulations. The masks for the Australian paleocoastlines are also provided in this folder. The R Scripts and outputs of HadCM3 modelled P-E and runoff are provided in their respective folders. R Scripts required to plot Figure 4 are provided in the "Latitudinal Distribution" folder. 


Due to file size constraints on GitHub, the bathymetry file (`gebco_australia_bathymetry.tif`, ~329 MB) used in `Time_Slice_Map/Time-Slice-Map.R` 
is not included in this repository. This file is publicly available from GEBCO: https://www.gebco.net/data_and_products/gridded_bathymetry_data/



\## License



Shield: \[!\[CC BY 4.0]\[cc-by-shield]]\[cc-by]



This work is licensed under a

\[Creative Commons Attribution 4.0 International License]\[cc-by].



\[!\[CC BY 4.0]\[cc-by-image]]\[cc-by]



\[cc-by]: http://creativecommons.org/licenses/by/4.0/

\[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png

\[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg



\## Dependencies



* Giorgino T (2009). “Computing and Visualizing Dynamic Time Warping Alignments in

&#x20; R: The dtw Package.” \_Journal of Statistical Software\_, \*31\*(7), 1-24.

&#x20; doi:10.18637/jss.v031.i07 <https://doi.org/10.18637/jss.v031.i07>.

&#x20; 

* Meyers, S.R. (2014). Astrochron: An R Package for Astrochronology.

&#x20; https://cran.r-project.org/package=astrochron



* Signorell A (2024). \_DescTools: Tools for Descriptive Statistics\_. R package

&#x20; version 0.99.54, <https://CRAN.R-project.org/package=DescTools>.



* GEBCO Compilation Group (2023) GEBCO 2023 Grid. https://doi.org/10.5285/f98b053b-0cbc-6c23-e053-6c86abc0af7b.



* Kocsis, Á. T., Raja, N. B., Williams, S. \& Dowding, E. M. rgplates: R interface for the GPlates Web Service and Desktop Application. (2024)



* Merdith, A. S. et al. Extending full-plate tectonic models into deep time: Linking the Neoproterozoic and the Phanerozoic. Earth-Science Reviews 214, 103477 (2021).



* Samant, R., Bialik, O. M. \& De Vleeschouwer, D. Industrial Gamma Radiation Signals: An Untapped Resource to Reconstruct Neogene Paleoclimate Off Northwest Australia. Paleoceanography and Paleoclimatology 40, e2025PA005157 (2025)



* Samant, R., Giorgino, T., Jarochowska, E. \& De Vleeschouwer, D. Enhancing the Accuracy of Dynamic Time Warping by Integrating Stratigraphic Constraints Into the Automated Correlation of Sedimentary Sequences. Paleoceanography and Paleoclimatology 40, e2024PA005082 (2025).



* Meyers, S. R. Astrochron: An R package for astrochronology. (2014).



* Thronberens, S., Back, S., Bourget, J., Allan, T. \& Reuning, L. 3-D seismic chronostratigraphy of reefs and drifts in the Browse Basin, NW Australia. GSA Bulletin 134, 3155–3175 (2022).



* McCaffrey, J. C., Wallace, M. W. \& Gallagher, S. J. A Cenozoic Great Barrier Reef on Australia’s North West shelf. Global and Planetary Change 184, 103048 (2020).



* McCaffrey, J. C. et al. The Rowley Shoals atolls: Remnants of a Miocene great barrier reef on the north-west Australian margin. Global and Planetary Change 245, 104688 (2025).



* Valdes, P. J. et al. The BRIDGE HadCM3 family of climate models: HadCM3@Bristol v1.0. Geoscientific Model Development 10, 3715–3743 (2017).



