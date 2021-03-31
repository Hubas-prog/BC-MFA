########################### SCRIPT ################################
# by Cédric Hubas 
########################### SCRIPT ################################

#####################
# PACKAGES
#####################
library(ade4)
library(ggplot2)
library(scales)
library(cowplot)
#library(devtools)

#####################
# AESTHETICS
#####################
col.pal <- colorRampPalette(c("red3","orange","green3","royalblue"))
row.pal <- colorRampPalette(c("#FFF587","#FF8C64","#FF665A","#7D6B7D","#A3A1A8"))
# please modify col.pal and row.pal to change scores and loadings plot label colors

#####################
# Description
#####################

# The folowing function performs a specific supervised analysis.
# It performs a Multiple Factor Analysis as described by Escofier and Pages [1] by using the dudi.pca function of package ade4.
# The generated objects of class pca and dudi. are then used to perform a Between Class Analysis as described by Dolédec and Chessel [2].
# The final result is a supervised MFA called BC-MFA

#####################
# Usage
#####################

#bc.mfa(df,bloc,fac,spcos=0,X=1,Y=2)

#####################
# Arguments
#####################
# note: The used must identify several groups of variables within the data frame in odert to build argument bloc.

# df => a data frame with n rows (individuals) and p columns (numeric variables).
# bloc => a vector or factor object giving the groups for the corresponding groups of variable of df.
# fac => an external factor used for the supervised analysis (BCA). 
# spcos => a numerical value giving the cos2 by which variables text and symbols should be magnified in the variable plot.
# X and Y => the dimension of the principal components to be plotted.

#####################
# Note
#####################

# The function returns a plot of the Multiple factor analysis ordination (top-left pannel)
# The function returns also a plot of the supervised MFA (i.e. BC-MFA) in the bottom-left pannel
# Total Inertia Explained (T.I.E) by the chosen factor is given in the second plot
# The function returns also a plot of the loadings (i.e. variables) of the BC-MFA
# percentages of inertia explained by each axis is also given in axis titles
# The functions also returns the result of the generic function randtest() whch performs a Monte-Carlo test of the BC-MFA
# Please note that the BCA allows the extraction of a number k-1 of components which will be a function of the number of modalities (k) of the factor.
# If the fac factor has less than 3 modalities, the analysis will not be possible. 

#####################
# geom_convexhull function
# by : Charles Martin https://github.com/cmartin/ggConvexHull
#####################

# devtools::install_github("cmartin/ggConvexHull",force=T)

StatConvexHull <- ggplot2::ggproto(
  "StatConvexHull",
  ggplot2::Stat,
  required_aes = c("x", "y"),
  compute_group = function(self, data, scales, params) {
    data[chull(data$x, data$y), ]
  }
)

stat_convexhull <- function(mapping = NULL, data = NULL, geom = "polygon",
                            position = "identity", show.legend = NA,
                            inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatConvexHull,
    data = data, mapping = mapping, geom = geom, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes, params = list(...)
  )
}

geom_convexhull <- function (mapping = NULL, data = NULL, stat = "convex_hull", position = "identity",
                             ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat, geom = ggplot2::GeomPolygon,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#####################
# bc.mfa function
# by Cédric Hubas : https://github.com/Hubas-prog/BC-MFA
#####################

bc.mfa<-function(df,bloc,fac,spcos=0,X=1,Y=2,...){
  v<-NULL
  for(i in 1:length(bloc)){
    v[i]<-rep(paste("Group",i,sep=""))
  }
  splitfac<-factor(rep(v,bloc))
  
  LI<-NULL
  for(i in 1:length(bloc)){
    LI[i]<-list(df[,splitfac==paste("Group",i,sep="")])
  }
  names(LI)<-v
  
# MFA
eig<-unlist(lapply(LI,function(x){dudi.pca(x,scannf=F,nf=2)$eig[1]}))
res.mfa<-dudi.pca(df,col.w=rep(1/eig,bloc),scannf=F,nf=ncol(df))
varexp1<-res.mfa$eig*100/sum(res.mfa$eig)
  
# BC-MFA
k<-length(levels(factor(fac)))
res.bcmfa<-bca(res.mfa,fac,scannf=F,nf=k-1)
varexp2<-res.bcmfa$eig*100/sum(res.bcmfa$eig)

#plots

mfaX<-res.mfa$li[,X]
mfaY<-res.mfa$li[,Y]

MFA.ind<-ggplot(res.mfa$li,aes(x=res.mfa$li[,X],y=res.mfa$li[,Y],col=fac))+
  geom_point()+
  geom_convexhull(alpha = 0.3,aes(fill = fac))+
  scale_fill_manual(values=row.pal(length(levels(fac))))+
  scale_colour_manual(values=row.pal(length(levels(fac))))+
  ggtitle("MFA scores")+
  xlab(paste("Axis ",X," : ",round(varexp1[X],2),"%"))+
  ylab(paste("Axis ",Y," : ",round(varexp1[Y],2),"%"))+
  labs(fill="Grouping variable",col="Grouping variable")+
  theme_bw()

BCMFA.ind<-ggplot(res.bcmfa$ls,aes(x=res.bcmfa$ls[,X],y=res.bcmfa$ls[,Y],col=fac))+
  geom_point()+
  scale_fill_manual(values=row.pal(length(levels(fac))))+
  scale_colour_manual(values=row.pal(length(levels(fac))))+
  geom_convexhull(alpha = 0.3,aes(fill = fac))+
  ggtitle(paste("BC-MFA scores"," | T.I.E=",
                round(res.bcmfa$ratio,2)*100,
                "%",
                sep=""))+
  xlab(paste("Axis ",X," : ",round(varexp2[X],2),"%"))+
  ylab(paste("Axis ",Y," : ",round(varexp2[Y],2),"%"))+
  labs(fill="Explanatory variable",col="Explanatory variable")+
  theme_bw()

cos2 <- as.matrix(res.bcmfa$co[,c(X,Y)])*as.matrix(res.bcmfa$co[,c(X,Y)])
var.group<-factor(rep(v,bloc),levels=v)
res.bcmfa$co$col<-col.pal(length(bloc))[var.group]
col.filter<-cos2[,1]>spcos | cos2[,2]>spcos

if(spcos==0) {
  BCMFA.var<-ggplot(res.bcmfa$co,aes(x=res.bcmfa$co[,X],
                                     y=res.bcmfa$co[,Y],col=var.group))+
    scale_colour_manual(values=col.pal(length(levels(var.group))))+
    geom_segment(aes(x = 0, y = 0,
                     xend = res.bcmfa$co[,X],
                     yend = res.bcmfa$co[,Y]),
                 arrow = arrow(length = unit(0.5, "cm")),
                 col="lightgrey")+
    annotate(geom="text",
             x=res.bcmfa$co[,X],
             y=res.bcmfa$co[,Y],
             label=rownames(res.mfa$co),
             size=4,
             col=res.bcmfa$co$col)+
    ggtitle("BC-MFA variables")+
    xlab(paste("Axis ",X," : ",round(varexp2[X],2),"%"))+
    ylab(paste("Axis ",Y," : ",round(varexp2[Y],2),"%"))+
    xlim(c(-1,1))+
    ylim(c(-1,1))+
    theme_bw()
}else {
  BCMFA.var<-ggplot(res.bcmfa$co,aes(x=res.bcmfa$co[,X],
                                     y=res.bcmfa$co[,Y],col=var.group))+
    geom_segment(aes(x = 0, y = 0,
                     xend = res.bcmfa$co[,X],
                     yend = res.bcmfa$co[,Y]),
                 arrow = arrow(length = unit(0.5, "cm")),
                 col="lightgrey")+
    annotate(geom="text",
             x=res.bcmfa$co[,X],
             y=res.bcmfa$co[,Y],
             label=rownames(res.mfa$co),
             size=4,
             col=alpha(res.bcmfa$co$col,c(0.15,1)[factor(col.filter)]))+
    ggtitle("BC-MFA variables")+
    xlab(paste("Axis ",X," : ",round(varexp2[X],2),"%"))+
    ylab(paste("Axis ",Y," : ",round(varexp2[Y],2),"%"))+
    xlim(c(-1,1))+
    ylim(c(-1,1))+
    theme_bw()
}

left<-plot_grid(MFA.ind,BCMFA.ind,labels=c("a","b"),ncol=1)
PLOT<-plot_grid(left,BCMFA.var,labels=c("","c"),ncol=2)
return(list(PLOT,perm.test=randtest(res.bcmfa),
                 BCMFAcos2=cos2,
                 BCMFAco=res.bcmfa$co,
                 BCMFAind=res.bcmfa$ls))
}

# col=alpha(res.bcmfa$co$col,c(0.15,1)[factor(col.filter)])


#####################
# Example1
#####################

data(meaudret)
fauna <- meaudret$spe
enviro <- meaudret$env
ISvariable1 <- meaudret$design$season
ISvariable2 <- meaudret$design$site
dataset <- data.frame(fauna,enviro)
bloc <- c(dim(fauna)[2],dim(enviro)[2])

bc.mfa(df=dataset,
       bloc=bloc,
       fac=factor(ISvariable1),
       X=1,
       Y=2)

#####################
# Example2
#####################
# From the following article
# Michelet C, Zeppilli D, Hubas C, Baldrighi E, Cuny P, Dirberg G, Militon C, Walcker R, Lamy D, Jézéquel R, Receveur J, Gilbert F, Houssainy AE, Dufour A, Heimbürger-Boavida L-E, Bihannic I, Sylvi L, Vivier B, Michaud E.
# First Assessment of the Benthic Meiofauna Sensitivity to Low Human-Impacted Mangroves in French Guiana.
# Forests. 2021; 12(3):338.
# https://doi.org/10.3390/f12030338
# see repository : https://github.com/Hubas-prog/Script-meiofauna-sensitivity
# Data available at : https://doi.org/10.5281/zenodo.4592299

conta<-read.csv("https://zenodo.org/record/4592300/files/Contaminants.csv?download=1",h=T,sep=";",dec=",")
pigments<-read.csv("https://zenodo.org/record/4592300/files/pigments.csv?download=1",h=T,sep=";",dec=",")
CHN<-read.csv("https://zenodo.org/record/4592300/files/CHN.csv?download=1",h=T,sep=";",dec=",")
enviro<-read.csv("https://zenodo.org/record/4592300/files/dataENV.csv?download=1",h=T,sep=";",dec=",")
DATA<-cbind(conta[,2:32],pigments[,2:19],enviro[,3:7],CHN[,2:3])
DATA$sites<-substr(conta$groupe,1,2)
DATA$layer<-paste("L",substr(conta$groupe,4,4),sep="")
DATA$core<-substr(conta$groupe,3,3)

bloc<-c(dim(conta[,2:32])[2],
        dim(pigments[,2:19])[2],
        dim(enviro[,3:7])[2],
        dim(CHN[,2:3])[2])

bc.mfa(df=DATA[,1:56],
       bloc=bloc,
     fac=factor(DATA$sites),
       spcos=0)


#####################
# References
#####################

#[1] Escofier B. and Pagès J. (1994) Multiple factor analysis (AFMULT package). Computational Statistics & Data Analysis. 18(1):121-140. https://doi.org/10.1016/0167-9473(94)90135-X
#[2] Dolédec, S. and Chessel, D. (1987) Rythmes saisonniers et composantes stationnelles en milieu aquatique I- Description d'un plan d'observations complet par projection de variables. Acta Oecologica, Oecologia Generalis, 8, 3, 403–426.
