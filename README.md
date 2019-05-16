# bias_correction
R code for climate model bias correction, following Hempel et al. 2013

There are two main steps to this bias correction process.

The following steps are done for temperature, and precipitation:
1. Calculate bias correction factors based on the observational data and historical model output in the reference period.

Tparam_calc.R: calculate semperature correction factors: 
- C, monthly mean correction
- B, transfer function slope 

Pparam_calc.R: calculate precipitation correction factors: 
- C, monthly mean correction
-	A, transfer function intercept
-	B, transfer function slope
-	Epsilon_m, monthly mean threshold for “wet” definition
-	Epsilon_d, daily threshold for “wet” definition

2.	Apply the correction factors to the application period and the historical period:

Apply_T_Param.R: apply the temperature correction factors

Apply_P_Param.R: apply the precipitation correction factors


