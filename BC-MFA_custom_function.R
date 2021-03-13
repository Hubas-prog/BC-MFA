########################### SCRIPT ################################
# by Cédric Hubas 
########################### SCRIPT ################################

#####################
# PACKAGES
#####################
library(ade4)
library(vegan)
library(scales)

#####################
# AESTHETICS
#####################
row.pal <- colorRampPalette(c("red","green","blue"))
col.pal <- colorRampPalette(c("red3","orange","green3","royalblue","purple"))

#####################
# Description
#####################

# The folowing function performs a specific supervised analysis.
# It performs a Multiple Factor Analysis by using the dudi.pca function of package ade4.
# The generated objects of class pca and dudi. are then used to perform a Between Class Analysis.
# The fianl result is a supervised MFA called BC-MFA

#####################
# Usage
#####################

#bc.mfa(df,bloc,fac,cos2)

#####################
# Arguments
#####################
#note: #The used must identify several groups of variables within the data frame in odert to build argument bloc.

#df => a data frame with n rows (individuals) and p columns (numeric variables).
# bloc => a vector or factor object giving the groups for the corresponding groups of variable of df.
# fac => an external factor used for the supervised analysis (BCA). 
# spcos => A numerical value giving the cos2 by which variables text and symbols should be magnified in the variable plot.

#####################
# Note
#####################

# The function returns a plot of the Multiple factor analysis ordination (top-left pannel)
# The function returns also a plot of the supervised MFA (i.e. BC-MFA) in the bottom-left pannel
# Total Inertia Explained (T.I.E) by the chosen factor is given in the second plot
# The function returns also a plot of the loadings (i.e. variables) of the BC-MFA
# percentages of inertia explained by each axis is also given in axis titles
# The functions also returns the result of the generic function randtest() whch performs a Monte-Carlo test of the BC-MFA

#####################
# Function
#####################

bc.mfa<-function(df,bloc,fac,spcos=0,...){
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
  res.mfa<-dudi.pca(df,col.w=rep(1/eig,bloc),scannf=F,nf=2)
  varexp1<-res.mfa$eig*100/sum(res.mfa$eig)
  
  # BC-MFA
  res.bcmfa<-bca(res.mfa,fac,scannf=F,nf=2)
  varexp2<-res.bcmfa$eig*100/sum(res.bcmfa$eig)
  
  #plots
  par(layout(matrix(c(1,2,3,3),nrow=2)),mar=c(4,4,4,4))
  
  plot(res.mfa$li[,2]~res.mfa$li[,1],
       pch=21,cex=2,col="white",
       bg=row.pal(length(levels(fac)))[fac],
       xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),
       ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),
       main="MFA scores")
  
  ordihull(res.mfa$li,fac,lab=T)
  
  plot(res.bcmfa$ls[,2]~res.bcmfa$ls[,1],
       pch=21,cex=2,col="white",
       bg=row.pal(length(levels(fac)))[fac],
       xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),
       ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),
       main=paste("BC-MFA scores","| T.I.E=",round(res.bcmfa$ratio,2)*100,"%",sep=""))
  
  ordihull(res.bcmfa$ls,fac,lab=T)
  
  plot(res.bcmfa$co[,2]~res.bcmfa$co[,1],type="n",
       xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),
       ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),
       ylim=c(-1,1),xlim=c(-1,1),
       main="BC-MFA loadings")
  
  arrows(x0=0,y0=0,x1=res.bcmfa$co[,1],y1=res.bcmfa$co[,2],col="lightgrey",length=0.1)
  
  cos2 <- as.matrix(res.bcmfa$co[,1:2])*as.matrix(res.bcmfa$co[,1:2])
  
  var.group<-factor(rep(v,bloc),
                    levels=v)
  res.bcmfa$co$col<-col.pal(length(bloc))[var.group]
  
  newco <- res.bcmfa$co[cos2[,1]>spcos | cos2[,2]>spcos,] 
  oldco <- res.bcmfa$co[cos2[,1]<spcos & cos2[,2]<spcos,]
  
  if(spcos==0){
    text(newco[,1], newco[,2],rownames(newco),col=newco$col,cex=1.2)
  }else{
    text(newco[,1], newco[,2],rownames(newco),col=newco$col,cex=1.2)
    text(oldco[,1], oldco[,2],rownames(oldco),col=alpha(oldco$col,0.3),cex=0.8)
  }
  
  abline(h=0,lty="dashed")
  abline(v=0,lty="dashed")
  
  legend("topleft",
         v,
         text.col=col.pal(length(bloc)),cex=0.8,box.lty=0)

  return(randtest(res.bcmfa))
}

#####################
# Examples
#####################
# From the following article
# First assessment of the benthic meiofauna sensitivity to low human-impacted mangroves in French Guiana
# by : Claire Michelet, Daniela Zeppilli, Cédric Hubas, Elisa Baldrighi, Philippe Cuny, Guillaume Dirberg, Cécile Militon, Romain Walcker, Dominique Lamy, Ronan Jézéquel, Justine Receveur, Franck Gilbert, Amonda El Houssainy, Aurélie Dufour, Lars-Eric Heimbürger-Boavida, Isabelle Bihannic, Léa Sylvi, Baptiste Vivier, Emma Michaud
# see repository : https://github.com/Hubas-prog/Script-meiofauna-sensitivity

conta<-read.csv("Contaminants.csv",h=T,sep=";",dec=",")
pigments<-read.csv("pigments.csv",h=T,sep=";",dec=",")
CHN<-read.csv("CHN.csv",h=T,sep=";",dec=",")
bact<-read.csv("Communaute-Bacterienne.csv",h=T,sep=";",dec=",")
enviro<-read.csv("dataENV.csv",h=T,sep=";",dec=",")
DATA<-cbind(conta[,2:32],pigments[,2:19],bact[,3:4],enviro[,3:7],CHN[,2:3])
DATA$sites<-substr(conta$groupe,1,2)
DATA$layer<-paste("L",substr(conta$groupe,4,4),sep="")
DATA$core<-substr(conta$groupe,3,3)

bloc<-c(dim(conta[,2:32])[2],
        dim(pigments[,2:19])[2],
        dim(bact[,3:4])[2],
        dim(enviro[,3:7])[2],
        dim(CHN[,2:3])[2])

bc.mfa(df=DATA[,1:58],bloc=bloc,fac=fact<-factor(DATA$sites),spcos=0.40)
