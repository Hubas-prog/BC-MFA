#################################
# Hubas et al.
# Identification of microbial exopolymer producers in sandy and muddy intertidal sediments by compound-specific isotope analysis
# Script version: V2.O
#################################

#################################
# Packages
#################################
library(ggplot2)
library(cowplot)
library(agricolae)
library(car)
library(RVAideMemoire)

#################################
# Aesthetics
#################################
my.palette <- colorRampPalette(c("#1AC0A4","#DAF7A6","#FFC300","#FF5733","#C70039","#900C3F"))
Clabiso<- expression(paste(delta^{13}, "C (\u2030)"))
Nlabiso<-expression(paste(delta^{15}, "N (\u2030)"))

#####################
# M&M
#####################

CHLA<-read.csv("https://zenodo.org/record/8116986/files/chla.csv?download=1",
               sep=";",
               h=T)
BRA<-read.csv("https://zenodo.org/record/8116986/files/bra.csv?download=1",
              sep=";",
              h=T)

CHLA.plot<-ggplot(CHLA,aes(y= CHLA,x= site)) +
  geom_boxplot()+
  ylab(expression(paste("Chlorophyll a (µg.g sediment dry weight (SDW",")"^-1,")")))+
  theme_bw(base_size = 14)+
  theme(axis.title.x=element_blank(),
        legend.position="none",
        axis.text.x=element_text(angle=90))

BRA.plot<-ggplot(BRA,aes(y= BRA,x= site)) +
  geom_boxplot()+
  ylab(expression(paste("Branched fatty acids (% of total)")))+
  theme_bw(base_size = 14)+
  theme(axis.title.x=element_blank(),
        legend.position="none",
        axis.text.x=element_text(angle=90))

Cratio.plot<-ggplot(BRA,aes(y= C16ratio,x= site)) +
  geom_boxplot()+
  ylab(expression(paste("16:0/16:1w7 ratio")))+
  theme_bw(base_size = 14)+
  theme(axis.title.x=element_blank(),
        legend.position="none",
        axis.text.x=element_text(angle=90))

plot_grid(Cratio.plot,CHLA.plot,BRA.plot,labels = c("a","b","c"),ncol=3)

#####################
# STATISTICS M&M
#####################

tapply(BRA$BRA,BRA$site,shapiro.test)
var.test(BRA$BRA~BRA$site)
wilcox.test(BRA$BRA~BRA$site)

tapply(BRA$C16ratio,BRA$site,shapiro.test)
var.test(BRA$C16ratio~BRA$site)
wilcox.test(BRA$C16ratio~BRA$site)

tapply(CHLA$CHLA,CHLA$site,shapiro.test)
var.test(CHLA$CHLA~CHLA$site)
t.test(CHLA$CHLA~CHLA$site,var.equal=F)

#####################
# FIGURE 2
#####################

EPS.data<-read.csv("https://zenodo.org/record/7351531/files/EPS_SI_ratio.csv?download=1",
                   sep=";",
                   h=T)

d13CNat<-ggplot(EPS.data,aes(y= d13C,x= groups,col=type)) +
  geom_boxplot()+ 
  ylab(Clabiso)+
  theme_bw(base_size = 14)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90))+
  labs(col="EPS type")+
  scale_color_manual(values=c("#e76f51","#264653"))

d15NNat<-ggplot(EPS.data,aes(y= d15N,x= groups,col=type)) +
  geom_boxplot()+ 
  ylab(Nlabiso)+
  theme_bw(base_size = 14)+
  labs(col="EPS type")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90))+
  scale_color_manual(values=c("#e76f51","#264653"))

Ccontent<-ggplot(EPS.data,aes(y= C.content,x= groups,col=type)) +
  geom_boxplot()+
  ylab(expression(paste("C content (µg.",mg^-1,")")))+
  theme_bw(base_size = 14)+
  theme(axis.title.x=element_blank(),
        legend.position="none",
        axis.text.x=element_text(angle=90))+
  scale_color_manual(values=c("#e76f51","#264653"))

Ncontent<-ggplot(EPS.data,aes(y= N.content,x= groups,col=type)) +
  geom_boxplot()+
  ylab(expression(paste("N content (µg.",mg^-1,")")))+
  theme_bw(base_size = 14)+
  theme(axis.title.x=element_blank(),
        legend.position="none",
        axis.text.x=element_text(angle=90))+
  scale_color_manual(values=c("#e76f51","#264653"))

plot_grid(Ccontent,d13CNat,Ncontent,d15NNat,labels = c("a","b","c","d"),ncol=2)

#####################
# STATISTICS Fig2
#####################

# creation of a cross factor (to test EPS type x sampling occasion)
EPS.data$f<-paste(EPS.data$type, EPS.data$groups)

# Van der Waerden test for Fig.2-C.content
res.waer1<-waerden.test(EPS.data$C.content,
                        factor(EPS.data$f),
                        group = TRUE)
res.waer1

# Van der Waerden test for Fig.2-N.content
res.waer2<-waerden.test(EPS.data$N.content,
                        factor(EPS.data$f),
                        group = TRUE)
res.waer2

# Van der Waerden test for Fig.2-d13C
res.waer3<-waerden.test(EPS.data$d13C,
                        factor(EPS.data$f),
                        group = TRUE)
res.waer3

# Van der Waerden test for Fig.2-d15N
res.waer4<-waerden.test(EPS.data$d15N,
                        factor(EPS.data$f),
                        group = TRUE)
res.waer4

# creation of a cross factor (to test EPS type x sediment type)
EPS.data$f2<-paste(EPS.data$type, EPS.data$site)

# Van der Waerden test for Fig.2-C.content
res.waer1.f2<-waerden.test(EPS.data$C.content,
                           factor(EPS.data$f2),
                           group = TRUE)
res.waer1.f2

# Van der Waerden test for Fig.2-N.content
res.waer2.f2<-waerden.test(EPS.data$N.content,
                           factor(EPS.data$f2),
                           group = TRUE)
res.waer2.f2
# Van der Waerden test for Fig.2-d13C
res.waer3.f2<-waerden.test(EPS.data$d13C,
                           factor(EPS.data$f2),
                           group = TRUE)
res.waer3.f2
# Van der Waerden test for Fig.2-d15N
res.waer4.f2<-waerden.test(EPS.data$d15N,
                           factor(EPS.data$f2),
                           group = TRUE)
res.waer4.f2

#####################
# FIGURE 3
#####################

FA.data<-read.csv("https://zenodo.org/record/7351531/files/Fatty_acids_SI_ratio.csv?download=1",
                  sep=";",
                  h=T)

boxFA<-ggplot(FA.data[FA.data$fatty.name!="24:0",],aes(y=d13Ccor,x=family))+
  geom_boxplot(aes(col=family))+
  theme_bw(base_size = 14)+
  ylab(Clabiso)+
  scale_color_manual(values=my.palette(4))+
  theme(legend.position="none",axis.title.x=element_blank())+
  labs(col="Fatty acid family")+
  facet_grid(rows=vars(site))+
  ylim(c(-40,-5))

plotC<-ggplot(EPS.data,aes(y= d13C,x=type)) +
  geom_boxplot(aes(col=type))+
  ylab(Clabiso) +
  scale_color_manual(values=c("#e76f51","#264653"))+
  theme_bw(base_size = 14)+
  theme(legend.position="none",axis.title.x=element_blank())+
  facet_grid(rows=vars(site))+
  labs(col="EPS type")+
  ylim(c(-40,-5))

plot_grid(plotC, boxFA,labels=c("a","b"),nrow=1)

#####################
# STATISTICS Fig3
#####################

# Stats about SFA 24:0
summary(FA.data[FA.data$fatty.name=="24:0",]$d13C)
summary(FA.data[FA.data$fatty.name=="24:0"
                &FA.data$site=="Mud",]$d13C)
summary(FA.data[FA.data$fatty.name=="24:0"
                &FA.data$site=="Sand",]$d13C)

# Stats Fig3a EPS (Mud and Sand separately)
site.EPS<-split(EPS.data,
                EPS.data$site)

var.test(site.EPS$Mud$d13C~site.EPS$Mud$type)
tapply(site.EPS$Mud$d13C,
       site.EPS$Mud$type,
       shapiro.test)
perm.t.test(site.EPS$Mud$d13C~site.EPS$Mud$type)

var.test(site.EPS$Sand$d13C~site.EPS$Sand$type)
tapply(site.EPS$Sand $d13C,
       site.EPS$Sand $type,
       shapiro.test)
t.test(site.EPS$Sand $d13C~site.EPS$Sand $type,
       var.equal=T)

# Stats Fig.3b FA (Mud and Sand separately)
site.FA<-split(FA.data[FA.data$fatty.name!="24:0",],
               FA.data[FA.data$fatty.name!="24:0",]$site)

leveneTest(d13Ccor ~ family,data = site.FA$Mud)
shapiro.test(aov(site.FA$Mud$d13Ccor ~ site.FA$Mud$family)$residuals)

res<-oneway.test(d13Ccor ~ family,
                 data = site.FA$Mud,var.equal=F)$statistic
F.perm=NULL
for(i in 1:9999){
  F.perm[i]<-oneway.test(d13Ccor ~ sample(site.FA$Mud$family),
                         data = site.FA$Mud,var.equal=F)
  
} # Permutation of welsh anova
perm.F<-unlist(F.perm)
(1+length(which(perm.F>res)))/9999 #p.value
hist(perm.F)
TukeyHSD(aov(site.FA$Mud$d13Ccor ~ site.FA$Mud$family))

## Stats Fig.3b Sand
leveneTest(d13Ccor ~ family,data = site.FA$Sand)
shapiro.test(aov(site.FA$Sand$d13Ccor ~ site.FA$Sand$family)$residuals)
boxplot(d13Ccor ~ family,data = site.FA$Sand)
oneway.test(d13Ccor ~ family,data = site.FA$Sand,var.equal=T)
TukeyHSD(aov(site.FA$Sand$d13Ccor ~ site.FA$Sand$family))

#####################
# FIGURE 4
#####################

density.SFA<-ggplot(FA.data[
  FA.data$fatty.name=="14:0"|
    FA.data$fatty.name=="15:0"|
    FA.data$fatty.name=="16:0"|
    FA.data$fatty.name=="17:0"|
    FA.data$fatty.name=="18:0"|
    FA.data$fatty.name=="22:0"|
    FA.data$fatty.name=="24:0",],aes(x=d13Ccor))+
  geom_density(aes(lty=fatty.name),
               col=my.palette(4)[2],
               fill="lightgrey",
               bw=2,
               alpha=0.2)+
  theme_bw(base_size = 14)+
  xlab(Clabiso)+
  facet_grid(~site)+
  theme(
    legend.position="right",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0,0.5,0,0.5), "lines"))+
  xlim(-30,-5)

density.BFA<-ggplot(FA.data[
  FA.data$fatty.name=="15:0anteiso"|
    FA.data$fatty.name=="15:0iso"|
    FA.data$fatty.name=="17:0iso",],aes(x=d13Ccor))+
  geom_density(aes(lty=fatty.name),
               col=my.palette(4)[1],
               fill="lightgrey",
               bw=2,
               alpha=0.2)+
  theme_bw(base_size = 14)+
  xlab(Clabiso)+
  facet_grid(~site)+
  theme(
    legend.position="right",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0,0.5,0,0.5), "lines"))+
  xlim(-30,-5)

density.MUFA<-ggplot(FA.data[
  FA.data$fatty.name=="16:1n-5"|
    FA.data$fatty.name=="16:1n-7"|
    FA.data$fatty.name=="17:1n-5"|
    FA.data$fatty.name=="17:1n-7"|
    FA.data$fatty.name=="18:1n-7"|
    FA.data$fatty.name=="18:1n-9+18:3n-3",],aes(x=d13Ccor))+
  geom_density(aes(lty=fatty.name),
               col=my.palette(4)[3],
               fill="lightgrey",
               bw=2,
               alpha=0.2)+
  theme_bw(base_size = 14)+
  xlab(Clabiso)+
  facet_grid(~site)+
  theme(
    legend.position="right",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0,0.5,0,0.5), "lines"))+
  xlim(-30,-5)

density.PUFA<-ggplot(FA.data[
  FA.data$fatty.name=="16:4w3+16:3n-4"|
    FA.data$fatty.name=="16:4n-1"|
    FA.data$fatty.name=="18:2n-6"|
    FA.data$fatty.name=="18:4n-3"|
    FA.data$fatty.name=="18:4n-6"|
    FA.data$fatty.name=="20:4n-6"|
    FA.data$fatty.name=="20:5n-3"|
    FA.data$fatty.name=="20:4n-6"|
    FA.data$fatty.name=="22:5n-3"|
    FA.data$fatty.name=="22:4n-6"|
    FA.data$fatty.name=="22:6n-3",],aes(x=d13Ccor))+
  geom_density(aes(lty=fatty.name),
               col=my.palette(4)[4],
               fill="lightgrey",
               bw=2,
               alpha=0.2)+
  theme_bw(base_size = 14)+
  xlab(Clabiso)+
  facet_grid(~site)+
  theme(
    legend.position="right",
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0,0.5,0,0.5), "lines"))+
  xlim(-30,-5)

EPS.data$cat<-factor(paste(EPS.data$site,EPS.data$type))
centroids <- aggregate(cbind(d15N,d13C)~cat,data=EPS.data,mean)
EPS.data2 <- merge(EPS.data,centroids,by="cat",
                   suffixes=c("",".centroid"))

biplots<-ggplot(EPS.data2,aes(y= d15N,x=d13C)) +
  geom_point(aes(col=type),
             size=3,
             show.legend = F)+
  geom_segment(aes(x=d13C.centroid,
                   y=d15N.centroid,
                   xend=d13C,
                   yend=d15N),
               col=rep(c("#e76f51","#264653"),2)[EPS.data2$cat])+
  facet_grid(~site)+
  ylab(Nlabiso) +
  xlab(Clabiso) +
  scale_color_manual(values=c("#e76f51","#264653"))+
  theme_bw(base_size = 14)+
  xlim(-30,-5)

# Final plot
densities<-plot_grid(density.BFA,
                     density.SFA,
                     density.MUFA,
                     density.PUFA,
                     biplots,
                     labels=c("a","b","c","d","e"),
                     ncol=1,
                     align = "v",
                     axis="r")
densities
#####################
# FIGURE SF1
#####################

EPS_colorimetry<-read.csv2("https://zenodo.org/record/7351531/files/EPS_colorimetry.csv?download=1",
                           h=T)

ggplot(EPS_colorimetry,aes(y=as.numeric(EPS),x= groups,col=type)) +
  geom_boxplot()+
  facet_grid(~compound)+
  ylab(expression(paste("EPS concentration (mg eqGlucose or BSA.",ml^-1,")")))+
  theme_bw(base_size = 14)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90))+
  labs(col="EPS type")+
  scale_color_manual(values=c("#e76f51","#264653"))


