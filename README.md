# Adult Chinook salmon diets in the Salish Sea

Data and code for Greentree et al. Adult Chinook salmon diets identify distinct forage communities in the Salish Sea.

If using these data for analyses of Chinook salmon diet, see `data/out/Salish-Sea-Chinook-salmon-diet.RDS`. Additional data are available on request (please email wesleygreentree@uvic.ca, uvicsalmondiet@gmail.com).

**Order of scripts:**

1.  `data-processing.R`: takes raw catch and stomach contents data from Adult Salmon Diet Program database and processes it into proportional diet composition data. The input data are not provided to protect angler confidentiality, meaning that the first six sections cannot be run (up to section 6: Save data). Sections 7 and 8 can be run with the output of section 6.
2.  `select-grid.R`: determines the optimal grid given the sampling coverage.
3.  `summer-grid-analysis.R` and `winter-grid-analysis.R`: runs cluster analysis and supporting analyses/visualizations in summer and winter.
4.  `grid-position-sensitivity.R`: runs sensitivity analyses for the optimal number of clusters and if the grid selection procedure is robust.
5. `size-sensitivity-analysis.R`: runs sensitivity analyses related to regionally-varying minimum size regulations.
6. `species-accumulation-curves.R`: assesses if the minimum sample size, 10 stomachs per grid cell, is appropriate.
7. `study-area-map.R`: produces study area maps for main text and supplemental material.

