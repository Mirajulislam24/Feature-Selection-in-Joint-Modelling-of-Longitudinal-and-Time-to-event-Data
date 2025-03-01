This repository includes the code for a Bayesian feature selection approach that models multiple longitudinal (3 continuous) outcomes together with a time-to-event assuming different association structures (value, slope, area and threshold).

Specifically:

"DataPrep": includes the preparation of the data (like ARIC)
"ModelJags": includes the joint model for jags and its implementation in R.
"SimData_1R2F: includes data simulation for the case where one risk factor with two features is important for the hazard.

How does it work:

Download all files and place them in one folder.
Set as working directory in R this folder.
Run the code in "ModelJags" for fitting the joint model.
