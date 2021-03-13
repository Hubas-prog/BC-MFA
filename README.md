# CUSTOM FUNCTION THAT PERFORMS A BETWEEN CLASS MULTIPLE FACTOR ANALYSIS (BC-MFA)

# Packages

library(ade4),
library(vegan),
library(scales),
library(factoextra)

# Description

The folowing function performs a specific supervised analysis. It performs a Multiple Factor Analysis by using the dudi.pca function of package ade4. The generated objects of class pca and dudi. are then used to perform a Between Class Analysis. The final result is a supervised MFA called BC-MFA

# Usage

bc.mfa(df,bloc,fac,cos2)

# Arguments

note: The user must identify several groups of variables within the data frame in odert to build argument bloc. Please modify col.pal and row.pal to change scores and loadings plot label colors

df => a data frame with n rows (individuals) and p columns (numeric variables).

bloc => a vector or factor object giving the groups for the corresponding groups of variable of df.

fac => an external factor used for the supervised analysis (BCA). 

cos2 => A numerical value giving the cos2 by which variables text and symbols should be magnified in the variable plot. by default = 0.

# Note

The function returns a plot of the Multiple factor analysis ordiantion (topleft pannel). The function returns also a plot of the supervised MFA (i.e. BC-MFA) in the bottom left pannel. Total Inertia Explained (T.I.E) by the chosen factor is given in the second plot. The function returns also a plot of the loadings (i.e. variables) of the BC-MFA. percentages of inertia explained by each axis is also given in axis titles. The functions also returns the result of the generic function randtest() whch performs a Monte-Carlo test of the BC-MFA
