# CUSTOM FUNCTION THAT PERFORMS A BETWEEN CLASS MULTIPLE FACTOR ANALYSIS (BC-MFA)

# Packages

library(ade4),
library(vegan),
library(scales)

# Description

The folowing function performs a specific supervised analysis. It performs a Multiple Factor Analysis (MFA) as described by Escofier and Pages [1] by using the dudi.pca function of package ade4. The generated objects of class pca and dudi. are then used to perform a Between Class Analysis as described by Dolédec and Chessel [2]. The final result is a supervised MFA called BC-MFA

# Usage

bc.mfa(df,bloc,fac,spcos=0,X=1,Y=2,...)

# Arguments

note: The user must identify several groups of variables within the data frame in odert to build argument bloc. Please modify col.pal and row.pal to change scores and loadings plot label colors

df => a data frame with n rows (individuals) and p columns (numeric variables).

bloc => a vector or factor object giving the groups for the corresponding groups of variable of df.

fac => an external factor used for the supervised analysis (BCA). 

spcos => A numerical value giving the cos2 by which variables text and symbols should be magnified in the variable plot. by default = 0.

X and Y => the dimension of the principal components to be plotted. Please note that the BCA allows the extraction of a number k-1 of components which will be a function of the number of modalities (k) of the factor. If the fac factor has less than 3 modalities, the analysis wont be possible. 

# Note

The function returns a plot of the Multiple factor analysis ordiantion (topleft pannel). The function returns also a plot of the supervised MFA (i.e. BC-MFA) in the bottom left pannel. Total Inertia Explained (T.I.E) by the chosen factor is given in the second plot. The function returns also a plot of the loadings (i.e. variables) of the BC-MFA. percentages of inertia explained by each axis is also given in axis titles. The functions also returns the result of the generic function randtest() whch performs a Monte-Carlo test of the BC-MFA

# References

[1] Escofier B. and Pagès J. (1994) Multiple factor analysis (AFMULT package). Computational Statistics & Data Analysis. 18(1):121-140. https://doi.org/10.1016/0167-9473(94)90135-X
[2] Dolédec, S. and Chessel, D. (1987) Rythmes saisonniers et composantes stationnelles en milieu aquatique I- Description d'un plan d'observations complet par projection de variables. Acta Oecologica, Oecologia Generalis, 8, 3, 403–426.
