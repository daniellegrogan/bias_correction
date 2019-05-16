# bias_correction
R code for climate model bias correction, following Hempel et al. 2013

The following steps are done for temperature, and precipitation:
1. Calculate bias correction factors based on the observational data and historical model output in the reference period.
a.	Temperature correction factors: Tparam_calc.R
i.	C, monthly mean correction
ii.	B, transfer function slope 

b.	Precipitation correction factors: Pparam_calc.R
i.	C, monthly mean correction
ii.	A, transfer function intercept
iii.	B, transfer function slope
iv.	Epsilon_m, monthly mean threshold for “wet” definition
v.	Epsolon_d, daily threshold for “wet” definition

B.	Apply the correction factors to the application period and the historical period:
Apply_T_Param.R
Apply_P_Param.R


