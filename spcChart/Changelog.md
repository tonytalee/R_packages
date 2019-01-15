# Changelog of spcChart

## [0.2.1] - 2018-12-18

### Added

- custom color set for spc plot

### Fixed

- arguments of plot width and height didn't pass to plot when color_var existed

### Changed

- release the constrains of group_var and color_var must be in info_var
- spc chart return a list contains spc chart, df with ooc labels and df with original info

## [0.2.0] - 2018-11-10

### Changed

- use cl_xxx() functions to count control limits
- options in xxx_chart(..., R_chart= TRUE) to choose NOT including R chart in Xbar-R for example.

### Added

- functions to count control limits and within-group standard deviation
