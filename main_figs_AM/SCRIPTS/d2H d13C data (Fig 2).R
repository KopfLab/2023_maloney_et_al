## Code for fattyacid/water d2H fractionation and biomass/substrate d13C-fractionation for DOI: 10.1073/pnas.2310771121 by A.Maloney
library(tidyverse)
library("DataCombine")
library("readxl")
#data----
Hdata = read_excel("DATA/GCIRMS H All Culture Summary.xlsx",sheet = "supp table no units") 
galactose2H = 40
glucose2H = 8
glycerol2H = -45
Ldata = read_excel("DATA/GCFID SUMMARY lipid.xlsx",sheet = "summary") 
Cdata = read_excel("DATA/EAIRMS C N All Culture Summary.xlsx",sheet = "summary") 

#FIG2 set up data parts----
FY4glucose_16<- filter(Hdata, Strain == "FY4" & type == "c" & `Lipid16:0` == "C16:0" & Media == "0.75% Glucose")
FY4glucose_batch<- filter(Hdata, Strain == "FY4" & type == "b" & `Lipid16:0` == "C16:0" & Media == "0.75% Glucose")
FY4galactose_16<- filter(Hdata, Strain == "FY4" & type == "c" & `Lipid16:0` == "C16:0" & Media == "0.75% Galactose")
FY4galactose_batch<- filter(Hdata, Strain == "FY4" & type == "b" & `Lipid16:0` == "C16:0" & Media == "0.75% Galactose")
EXFglucose_16<- filter(Hdata, Strain == "EXF-4126" & type == "c" & `Lipid16:0` == "C16:0" & Media == "2% Glucose")
EXFglucose_batch<- filter(Hdata, Strain == "EXF-4126" & type == "b" & `Lipid16:0` == "C16:0" & Media == "2% Glucose")
EXFglycerol_16<- filter(Hdata, Strain == "EXF-4126" & type == "c" & `Lipid16:0` == "C16:0" & Media == "2% Glycerol")
EXFglycerol_batch<- filter(Hdata, Strain == "EXF-4126" & type == "b" & `Lipid16:0` == "C16:0" & Media == "2% Glycerol")
FY4glucose_bigbatch<- filter(Hdata, Strain == "FY4" & type == "bigbatch" & `Lipid16:0` == "C16:0" & Media == "0.75% Glucose")
FY4galactose_bigbatch<- filter(Hdata, Strain == "FY4" & type == "bigbatch" & `Lipid16:0` == "C16:0" & Media == "0.75% Galactose")
EXFglucose_bigbatch<- filter(Hdata, Strain == "EXF-4126" & type == "bigbatch" & `Lipid16:0` == "C16:0" & Media == "2% Glucose")
EXFglycerol_bigbatch<- filter(Hdata, Strain == "EXF-4126" & type == "bigbatch" & `Lipid16:0` == "C16:0" & Media == "2% Glycerol")

FY4glucose_L<- filter(Ldata, Strain == "FY4" & type == "chemostat" & Substrate == "glucose")
FY4galactose_L<- filter(Ldata, Strain == "FY4" & type == "chemostat" & Substrate == "galactose")
EXFglucose_L<- filter(Ldata, Strain == "EXF4126" & type == "chemostat" & Substrate == "glucose")
EXFglycerol_L<- filter(Ldata, Strain == "EXF4126" & type == "chemostat" & Substrate == "glycerol")

FY4glucose_C<- filter(Cdata, Strain == "FY4" & type == "chemostat" & Substrate == "glucose")
FY4galactose_C<- filter(Cdata, Strain == "FY4" & type == "chemostat" & Substrate == "galactose")
EXFglucose_C<- filter(Cdata, Strain == "EXF4126" & type == "chemostat" & Substrate == "glucose")
EXFglycerol_C<- filter(Cdata, Strain == "EXF4126" & type == "chemostat" & Substrate == "glycerol")
FY4glucose_Cbatch<- filter(Cdata, Strain == "FY4" & type == "batch" & Substrate == "glucose")
FY4galactose_Cbatch<- filter(Cdata, Strain == "FY4" & type == "batch" & Substrate == "galactose")
EXFglucose_Cbatch<- filter(Cdata, Strain == "EXF4126" & type == "batch" & Substrate == "glucose")
EXFglycerol_Cbatch<- filter(Cdata, Strain == "EXF4126" & type == "batch" & Substrate == "glycerol")
FY4glucose_Cbigbatch<- filter(Cdata, Strain == "FY4" & type == "bigbatch" & Substrate == "glucose")
FY4galactose_Cbigbatch<- filter(Cdata, Strain == "FY4" & type == "bigbatch" & Substrate == "galactose")
EXFglucose_Cbigbatch<- filter(Cdata, Strain == "EXF4126" & type == "bigbatch" & Substrate == "glucose")
EXFglycerol_Cbigbatch<- filter(Cdata, Strain == "EXF4126" & type == "bigbatch" & Substrate == "glycerol")

FY4galactose_16slow <- filter(FY4galactose_16, D ==1.34)
FY4galactose_16med <- filter(FY4galactose_16, D ==3.59)
FY4galactose_16fast <- filter(FY4galactose_16, D > 5)
FY4glucose_16slow <- filter(FY4glucose_16, D ==1.18)
FY4glucose_16med <- filter(FY4glucose_16, D ==3.36)
FY4glucose_16fast <- filter(FY4glucose_16, D > 5)
EXFglucose_16slow <- filter(EXFglucose_16, D ==0.55)
EXFglucose_16med <- filter(EXFglucose_16, D ==1.16)
EXFglucose_16fast <- filter(EXFglucose_16, D > 2)
EXFglycerol_16slow <- filter(EXFglycerol_16, D ==0.65)
EXFglycerol_16med <- filter(EXFglycerol_16, D ==1.35)
EXFglycerol_16fast <- filter(EXFglycerol_16, D > 2)

FY4galactose_Lslow <- filter(FY4galactose_L, D ==1.34)
FY4galactose_Lmed <- filter(FY4galactose_L, D ==3.59)
FY4galactose_Lfast <- filter(FY4galactose_L, D > 5)
FY4glucose_Lslow <- filter(FY4glucose_L, D ==1.18)
FY4glucose_Lmed <- filter(FY4glucose_L, D ==3.36)
FY4glucose_Lfast <- filter(FY4glucose_L, D > 5)
EXFglucose_Lslow <- filter(EXFglucose_L, D ==0.55)
EXFglucose_Lmed <- filter(EXFglucose_L, D ==1.16)
EXFglucose_Lfast <- filter(EXFglucose_L, D > 2)
EXFglycerol_Lslow <- filter(EXFglycerol_L, D ==0.65)
EXFglycerol_Lmed <- filter(EXFglycerol_L, D ==1.35)
EXFglycerol_Lfast <- filter(EXFglycerol_L, D > 2)

FY4galactose_Cslow <- filter(FY4galactose_C, D ==1.34)
FY4galactose_Cmed <- filter(FY4galactose_C, D ==3.59)
FY4galactose_Cfast <- filter(FY4galactose_C, D > 5)
FY4glucose_Cslow <- filter(FY4glucose_C, D ==1.18)
FY4glucose_Cmed <- filter(FY4glucose_C, D ==3.36)
FY4glucose_Cfast <- filter(FY4glucose_C, D > 5)
EXFglucose_Cslow <- filter(EXFglucose_C, D ==0.55)
EXFglucose_Cmed <- filter(EXFglucose_C, D ==1.16)
EXFglucose_Cfast <- filter(EXFglucose_C, D > 2)
EXFglycerol_Cslow <- filter(EXFglycerol_C, D ==0.65)
EXFglycerol_Cmed <- filter(EXFglycerol_C, D ==1.35)
EXFglycerol_Cfast <- filter(EXFglycerol_C, D > 2)




#Fig 2 dD and d13C----
pdf("FIGURES/d2H d13C FIG 2.pdf", width=4, height=5.5,encoding="MacRoman")
par(mar = c(0,.2,0,0), oma=c(6.1,4.5,3,1))
par(mfrow=c(2,4))
#dD----
#EXF glycerol
y=c(mean(EXFglycerol_Lslow$WeightAvg),mean(EXFglycerol_Lmed$WeightAvg),mean(EXFglycerol_Lfast$WeightAvg))
x=c(mean(EXFglycerol_Lslow$D),mean(EXFglycerol_Lmed$D),mean(EXFglycerol_Lfast$D))
plot(NULL, xlim=c(0,3.5), ylim=c(-250,400), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)
abline(h=c(0))
legend(-2.5,520,"",title="A",bty="n",x.intersp=.1,y.intersp=.1,text.font=2,xpd=NA,cex=2)

axis(side=2,mgp=c(2,.5,0),las=1)
mtext(side=2,expression(paste("   "^"2", epsilon,"","" ["fattyacid/water"], " (\u2030)")),line=2.5,cex=1)
rect(x-0.3, 0, x+0.3, y, col="skyblue")

ySize=c(mean(EXFglycerol_Lslow$`percentC160`),mean(EXFglycerol_Lmed$`percentC160`),mean(EXFglycerol_Lfast$`percentC160`))
y=c(mean(EXFglycerol_Lslow$`ε16:0`),mean(EXFglycerol_Lmed$`ε16:0`),mean(EXFglycerol_Lfast$`ε16:0`))
points(x,y, cex=ySize*4,pch=21,bg="black")
ySize=c(mean(EXFglycerol_Lslow$`percentC161`),mean(EXFglycerol_Lmed$`percentC161`),mean(EXFglycerol_Lfast$`percentC161`))
y=c(mean(EXFglycerol_Lslow$`ε16:1`),mean(EXFglycerol_Lmed$`ε16:1`),mean(EXFglycerol_Lfast$`ε16:1`))
points(x,y, cex=ySize*4,pch=21,col="black")
ySize=c(mean(EXFglycerol_Lslow$`percentC181`),mean(EXFglycerol_Lmed$`percentC181`),mean(EXFglycerol_Lfast$`percentC181`))
y=c(mean(EXFglycerol_Lslow$`ε18:1`),mean(EXFglycerol_Lmed$`ε18:1`),mean(EXFglycerol_Lfast$`ε18:1`))
rect(0, 410, 7, 500, col=adjustcolor("gray", alpha.f=0.4),xpd=NA, border=NA)
text(3.5,480,"strain",cex=1,col="black",xpd=NA)
text(3.5,440,"EXF-4126",cex=1,col="black",xpd=NA)

arrows(0.65,mean(EXFglycerol_Lslow$WeightAvg)-sd(EXFglycerol_Lslow$WeightAvg),0.65,mean(EXFglycerol_Lslow$WeightAvg)+sd(EXFglycerol_Lslow$WeightAvg),code=3,length=0)
arrows(1.35,mean(EXFglycerol_Lmed$WeightAvg)-sd(EXFglycerol_Lmed$WeightAvg),1.35,mean(EXFglycerol_Lmed$WeightAvg)+sd(EXFglycerol_Lmed$WeightAvg),code=3,length=0)
arrows(2.69,mean(EXFglycerol_Lfast$WeightAvg)-sd(EXFglycerol_Lfast$WeightAvg),2.69,mean(EXFglycerol_Lfast$WeightAvg)+sd(EXFglycerol_Lfast$WeightAvg),code=3,length=0)

#EXF glucose
y=c(mean(EXFglucose_Lslow$WeightAvg),mean(EXFglucose_Lmed$WeightAvg),mean(EXFglucose_Lfast$WeightAvg))
x=c(mean(EXFglucose_Lslow$D),mean(EXFglucose_Lmed$D),mean(EXFglucose_Lfast$D))
plot(NULL, xlim=c(0,3.5), ylim=c(-250,400), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)
abline(h=c(0))
rect(x-0.3, 0, x+0.3, y, col="orange3")
ySize=c(mean(EXFglucose_Lslow$`percentC160`),mean(EXFglucose_Lmed$`percentC160`),mean(EXFglucose_Lfast$`percentC160`))
y=c(mean(EXFglucose_Lslow$`ε16:0`),mean(EXFglucose_Lmed$`ε16:0`),mean(EXFglucose_Lfast$`ε16:0`))
points(x,y, cex=ySize*4,pch=21,bg="black")
ySize=c(mean(EXFglucose_Lslow$`percentC161`),mean(EXFglucose_Lmed$`percentC161`),mean(EXFglucose_Lfast$`percentC161`))
y=c(mean(EXFglucose_Lslow$`ε16:1`),mean(EXFglucose_Lmed$`ε16:1`),mean(EXFglucose_Lfast$`ε16:1`))
points(x,y, cex=ySize*4,pch=21,col="black",xpd=NA)
arrows(0.55,mean(EXFglucose_Lslow$WeightAvg)-sd(EXFglucose_Lslow$WeightAvg),0.55,mean(EXFglucose_Lslow$WeightAvg)+sd(EXFglucose_Lslow$WeightAvg),code=3,length=0)
arrows(1.16,mean(EXFglucose_Lmed$WeightAvg)-sd(EXFglucose_Lmed$WeightAvg),1.16,mean(EXFglucose_Lmed$WeightAvg)+sd(EXFglucose_Lmed$WeightAvg),code=3,length=0)
arrows(2.37,mean(EXFglucose_Lfast$WeightAvg)-sd(EXFglucose_Lfast$WeightAvg),2.37,mean(EXFglucose_Lfast$WeightAvg)+sd(EXFglucose_Lfast$WeightAvg),code=3,length=0)

abline(v=3.5, lwd=1.5, lty=3)

#FY4 glucose
y=c(mean(FY4glucose_Lslow$WeightAvg,na.rm=T),mean(FY4glucose_Lmed$WeightAvg,na.rm=T),mean(FY4glucose_Lfast$WeightAvg,na.rm=T))
x=c(mean(FY4glucose_Lslow$D),mean(FY4glucose_Lmed$D),mean(FY4glucose_Lfast$D))
plot(NULL, xlim=c(0,7), ylim=c(-250,400), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)
abline(h=c(0))
rect(x-0.6, 0, x+0.6, y, col="orange")
ySize=c(mean(FY4glucose_Lslow$`percentC160`,na.rm=T),mean(FY4glucose_Lmed$`percentC160`,na.rm=T),mean(FY4glucose_Lfast$`percentC160`,na.rm=T))
y=c(mean(FY4glucose_Lslow$`ε16:0`,na.rm=T),mean(FY4glucose_Lmed$`ε16:0`,na.rm=T),mean(FY4glucose_Lfast$`ε16:0`,na.rm=T))
points(x,y, cex=ySize*4,pch=21,bg="black")
ySize=c(mean(FY4glucose_Lslow$`percentC161`,na.rm=T),mean(FY4glucose_Lmed$`percentC161`,na.rm=T),mean(FY4glucose_Lfast$`percentC161`,na.rm=T))
y=c(mean(FY4glucose_Lslow$`ε16:1`,na.rm=T),mean(FY4glucose_Lmed$`ε16:1`,na.rm=T),mean(FY4glucose_Lfast$`ε16:1`,na.rm=T))
points(x,y, cex=ySize*4,pch=21,col="black")

rect(0, 410, 15, 500, col=adjustcolor("gray", alpha.f=0.4),xpd=NA, border=NA)
text(7,480,"strain",cex=1,col="black",xpd=NA)
text(7,440,"FY4",cex=1,col="black",xpd=NA)

arrows(1.18,mean(FY4glucose_Lslow$WeightAvg)-sd(FY4glucose_Lslow$WeightAvg),1.18,mean(FY4glucose_Lslow$WeightAvg)+sd(FY4glucose_Lslow$WeightAvg),code=3,length=0)
arrows(3.36,mean(FY4glucose_Lmed$WeightAvg)-sd(FY4glucose_Lmed$WeightAvg),3.36,mean(FY4glucose_Lmed$WeightAvg)+sd(FY4glucose_Lmed$WeightAvg),code=3,length=0)
arrows(5.68333,mean(FY4glucose_Lfast$WeightAvg)-sd(FY4glucose_Lfast$WeightAvg),5.68333,mean(FY4glucose_Lfast$WeightAvg)+sd(FY4glucose_Lfast$WeightAvg),code=3,length=0)

#FY4 galactose
y=c(mean(FY4galactose_Lslow$WeightAvg,na.rm=T),mean(FY4galactose_Lmed$WeightAvg,na.rm=T),mean(FY4galactose_Lfast$WeightAvg,na.rm=T))
x=c(mean(FY4galactose_Lslow$D),mean(FY4galactose_Lmed$D),mean(FY4galactose_Lfast$D))
plot(NULL, xlim=c(0,7), ylim=c(-250,400), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)
abline(h=c(0))
rect(x-0.6, 0, x+0.6, y, col="gold",xpd=NA)
ySize=c(mean(FY4galactose_Lslow$`percentC160`,na.rm=T),mean(FY4galactose_Lmed$`percentC160`,na.rm=T),mean(FY4galactose_Lfast$`percentC160`,na.rm=T))
y=c(mean(FY4galactose_Lslow$`ε16:0`,na.rm=T),mean(FY4galactose_Lmed$`ε16:0`,na.rm=T),mean(FY4galactose_Lfast$`ε16:0`,na.rm=T))
points(x,y, cex=ySize*4,pch=21,bg="black")
ySize=c(mean(FY4galactose_Lslow$`percentC161`,na.rm=T),mean(FY4galactose_Lmed$`percentC161`,na.rm=T),mean(FY4galactose_Lfast$`percentC161`,na.rm=T))
y=c(mean(FY4galactose_Lslow$`ε16:1`,na.rm=T),mean(FY4galactose_Lmed$`ε16:1`,na.rm=T),mean(FY4galactose_Lfast$`ε16:1`,na.rm=T))
points(x,y, cex=ySize*4,pch=21,col="black",xpd=NA)

arrows(1.34,mean(FY4galactose_Lslow$WeightAvg)-sd(FY4galactose_Lslow$WeightAvg),1.34,mean(FY4galactose_Lslow$WeightAvg)+sd(FY4galactose_Lslow$WeightAvg),code=3,length=0)
arrows(3.59,mean(FY4galactose_Lmed$WeightAvg)-sd(FY4galactose_Lmed$WeightAvg),3.59,mean(FY4galactose_Lmed$WeightAvg)+sd(FY4galactose_Lmed$WeightAvg),code=3,length=0)
arrows(6.51,mean(FY4galactose_Lfast$WeightAvg)-sd(FY4galactose_Lfast$WeightAvg),6.51,mean(FY4galactose_Lfast$WeightAvg)+sd(FY4galactose_Lfast$WeightAvg),code=3,length=0)

points(c(-2.2,.2),c(300,300),pch=21,col=c("black","black"),bg=c("black","white"),xpd=NA,cex=c(.2*4))
points(c(-2.2,.2),c(250,250),pch=21,col=c("black","black"),bg=c("black","white"),xpd=NA,cex=c(.4*4))
points(c(-2.2,.2),c(200,200),pch=21,col=c("black","black"),bg=c("black","white"),xpd=NA,cex=c(.6*4))
text(-3.2,345,expression(paste("",italic(""),"C" ["16:0"],italic("")," C" ["16:1"])),xpd=NA,cex=1,adj=0)
text(-4,300,"20%",xpd=NA)
text(-4,250,"40%",xpd=NA)
text(-4,200,"60%",xpd=NA)

#13C----
#EXF glycerol
y=c(mean(EXFglycerol_Cslow$C_epsilon),mean(EXFglycerol_Cmed$C_epsilon),mean(EXFglycerol_Cfast$C_epsilon))
x=c(mean(EXFglycerol_Cslow$D),mean(EXFglycerol_Cmed$D),mean(EXFglycerol_Cfast$D))
plot(NULL, xlim=c(0,3.5), ylim=c(0,7.5), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)
legend(-2.5,7.8,"",title="B",bty="n",x.intersp=.1,y.intersp=.1,text.font=2,xpd=NA,cex=2)

axis(side=2,mgp=c(2,.5,0),las=1)
mtext(side=2,expression(paste(""^"13", epsilon,"","" ["biomass/substrate"], " (\u2030)         ")),line=2.5,cex=1)
rect(x-0.3, 0, x+0.3, y, col="skyblue")
abline(h=c(0),lwd=2)
axis(side=1, at=c(0,1,2,3),cex=.8,mgp=c(2,.5,0),labels=NA)
axis(side=1, at=c(0,1,2,3),cex=.8,mgp=c(2,.5,0))
text(1.5,-.75,expression(paste(" d"^"-1")),cex=.9,xpd=NA)
text(1.5,-1.2,"2% Glycerol",cex=1,col="black",xpd=NA)
rect(0, -1.5, 3, -2, col=adjustcolor("skyblue", alpha.f=1),xpd=NA, border=NA)
text(1.5,-1.75,"Respiration",cex=1,col="black",xpd=NA)

arrows(0.65,mean(EXFglycerol_Cslow$C_epsilon)-sd(EXFglycerol_Cslow$C_epsilon),0.65,mean(EXFglycerol_Cslow$C_epsilon)+sd(EXFglycerol_Cslow$C_epsilon),code=3,length=0)
arrows(1.35,mean(EXFglycerol_Cmed$C_epsilon)-sd(EXFglycerol_Cmed$C_epsilon),1.35,mean(EXFglycerol_Cmed$C_epsilon)+sd(EXFglycerol_Cmed$C_epsilon),code=3,length=0)
arrows(2.69,mean(EXFglycerol_Cfast$C_epsilon)-sd(EXFglycerol_Cfast$C_epsilon),2.69,mean(EXFglycerol_Cfast$C_epsilon)+sd(EXFglycerol_Cfast$C_epsilon),code=3,length=0)

#EXF glucose
y=c(mean(EXFglucose_Cslow$C_epsilon),mean(EXFglucose_Cmed$C_epsilon),mean(EXFglucose_Cfast$C_epsilon))
x=c(mean(EXFglucose_Cslow$D),mean(EXFglucose_Cmed$D),mean(EXFglucose_Cfast$D))
plot(NULL, xlim=c(0,3.5), ylim=c(0,7.5), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)
abline(h=c(0))
rect(x-0.3, 0, x+0.3, y, col="orange3")
axis(side=1, at=c(0,1,2,3),cex=.8,mgp=c(2,.5,0),labels=NA)
axis(side=1, at=c(0,1,2,3),cex=.8,mgp=c(2,.5,0))
text(1.5,-.75,expression(paste(" d"^"-1")),cex=.9,xpd=NA)
text(1.5,-1.2,"2% Glucose",cex=1,col="black",xpd=NA)
rect(0, -1.5, 11, -2, col=adjustcolor("orange3", alpha.f=1),border=NA,xpd=NA)
abline(h=c(0),lwd=2)
arrows(0.55,mean(EXFglucose_Cslow$C_epsilon)-sd(EXFglucose_Cslow$C_epsilon),0.55,mean(EXFglucose_Cslow$C_epsilon)+sd(EXFglucose_Cslow$C_epsilon),code=3,length=0)
arrows(1.16,mean(EXFglucose_Cmed$C_epsilon)-sd(EXFglucose_Cmed$C_epsilon),1.16,mean(EXFglucose_Cmed$C_epsilon)+sd(EXFglucose_Cmed$C_epsilon),code=3,length=0)
arrows(2.37,mean(EXFglucose_Cfast$C_epsilon)-sd(EXFglucose_Cfast$C_epsilon),2.37,mean(EXFglucose_Cfast$C_epsilon)+sd(EXFglucose_Cfast$C_epsilon),code=3,length=0)

abline(v=3.5, lwd=1.5, lty=3)
text(3.5, -2.5, "Growth Rate", cex=1.5,xpd=NA)
#FY4 glucose
y=c(mean(FY4glucose_Cslow$C_epsilon),mean(FY4glucose_Cmed$C_epsilon),mean(FY4glucose_Cfast$C_epsilon))
x=c(mean(FY4glucose_Cslow$D),mean(FY4glucose_Cmed$D),mean(FY4glucose_Cfast$D))
plot(NULL, xlim=c(0,7), ylim=c(0,7.5), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)
abline(h=c(0))
rect(x-0.6, 0, x+0.6, y, col="orange")
abline(h=c(0),lwd=2)
axis(side=1, at=c(0,1,2,3,4,5,6),cex=.8,mgp=c(2,.5,0),labels=NA)
axis(side=1, at=c(0,3,6),cex=.8,mgp=c(2,.5,0))
text(4,-.75,expression(paste(" d"^"-1")),cex=.9,xpd=NA)
text(3.25,-1.2,"  .75% Glucose",cex=1,col="black",xpd=NA)

text(3.5,-1.75,"Fermentation",cex=1,col="black",xpd=NA)
arrows(1.18,mean(FY4glucose_Cslow$C_epsilon)-sd(FY4glucose_Cslow$C_epsilon),1.18,mean(FY4glucose_Cslow$C_epsilon)+sd(FY4glucose_Cslow$C_epsilon),code=3,length=0)
arrows(3.36,mean(FY4glucose_Cmed$C_epsilon)-sd(FY4glucose_Cmed$C_epsilon),3.36,mean(FY4glucose_Cmed$C_epsilon)+sd(FY4glucose_Cmed$C_epsilon),code=3,length=0)
arrows(5.683333,mean(FY4glucose_Cfast$C_epsilon)-sd(FY4glucose_Cfast$C_epsilon),5.6833331,mean(FY4glucose_Cfast$C_epsilon)+sd(FY4glucose_Cfast$C_epsilon),code=3,length=0)

#FY4 galactose
y=c(mean(FY4galactose_Cslow$C_epsilon),mean(FY4galactose_Cmed$C_epsilon),mean(FY4galactose_Cfast$C_epsilon))
x=c(mean(FY4galactose_Cslow$D),mean(FY4galactose_Cmed$D),mean(FY4galactose_Cfast$D))
plot(NULL, xlim=c(0,7), ylim=c(0,7.5), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)
abline(h=c(0))
rect(x-0.6, 0, x+0.6, y, col="gold",xpd=NA)
abline(h=c(0),lwd=2)
axis(side=1, at=c(0,1,2,3,4,5,6,7),cex=.8,mgp=c(2,.5,0),labels=NA)
axis(side=1, at=c(0,3,6),cex=.8,mgp=c(2,.5,0))
text(4,-.75,expression(paste(" d"^"-1")),cex=.9,xpd=NA)
text(3.7,-1.2,".75% Galactose  ",cex=1,col="black",xpd=NA)
arrows(1.34,mean(FY4galactose_Cslow$C_epsilon)-sd(FY4galactose_Cslow$C_epsilon),1.34,mean(FY4galactose_Cslow$C_epsilon)+sd(FY4galactose_Cslow$C_epsilon),code=3,length=0)
arrows(3.59,mean(FY4galactose_Cmed$C_epsilon)-sd(FY4galactose_Cmed$C_epsilon),3.59,mean(FY4galactose_Cmed$C_epsilon)+sd(FY4galactose_Cmed$C_epsilon),code=3,length=0)
arrows(6.51,mean(FY4galactose_Cfast$C_epsilon)-sd(FY4galactose_Cfast$C_epsilon),6.51,mean(FY4galactose_Cfast$C_epsilon)+sd(FY4galactose_Cfast$C_epsilon),code=3,length=0)
#off----
dev.off()

#OLD STUFF----
#C/N----
pdf("FIGURES/CN ratio (extra, not used).pdf", width=3.9, height=3.2,encoding="MacRoman")
par(mar = c(2,.2,2,0), oma=c(3.3,3.7,1,3.7))
par(mfrow=c(1,4))
#EXF glycerol
y=c(mean(EXFglycerol_Cslow$CN),mean(EXFglycerol_Cmed$CN),mean(EXFglycerol_Cfast$CN))
x=c(mean(EXFglycerol_Cslow$D),mean(EXFglycerol_Cmed$D),mean(EXFglycerol_Cfast$D))
plot(NULL, xlim=c(0,3.5), ylim=c(0,25), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)

axis(side=2,mgp=c(2,.5,0),las=1)
mtext(side=2,"C/N",line=1.5,cex=1)
rect(x-0.3, 0, x+0.3, y, col="skyblue")
abline(h=c(0),lwd=2)
axis(side=1, at=c(0,1,2,3),cex=.8,mgp=c(2,.5,0),labels=NA)
axis(side=1, at=c(0,1,2,3),cex=.8,mgp=c(2,.5,0))
text(1.5,-3,expression(paste(" d"^"-1")),cex=.9,xpd=NA)
text(1.5,-4.5,"2% Glycerol",cex=1,col="black",xpd=NA)
rect(0, -5.5, 3, -6.5, col=adjustcolor("skyblue", alpha.f=0.3),xpd=NA, border=NA)
text(1.5,-6,"Respiration",cex=1,col="black",xpd=NA)

arrows(0.65,mean(EXFglycerol_Cslow$CN)-sd(EXFglycerol_Cslow$CN),0.65,mean(EXFglycerol_Cslow$CN)+sd(EXFglycerol_Cslow$CN),code=3,length=0)
arrows(1.35,mean(EXFglycerol_Cmed$CN)-sd(EXFglycerol_Cmed$CN),1.35,mean(EXFglycerol_Cmed$CN)+sd(EXFglycerol_Cmed$CN),code=3,length=0)
arrows(2.69,mean(EXFglycerol_Cfast$CN)-sd(EXFglycerol_Cfast$CN),2.69,mean(EXFglycerol_Cfast$CN)+sd(EXFglycerol_Cfast$CN),code=3,length=0)
rect(0, 22.5, 7, 25.5, col=adjustcolor("gray", alpha.f=0.3),xpd=NA, border=NA)
text(3.5,25,"strain",cex=1,col="black",xpd=NA)
text(3.5,23.5,"EXF-4126",cex=1,col="black",xpd=NA)

#EXF glucose
y=c(mean(EXFglucose_Cslow$CN),mean(EXFglucose_Cmed$CN),mean(EXFglucose_Cfast$CN))
x=c(mean(EXFglucose_Cslow$D),mean(EXFglucose_Cmed$D),mean(EXFglucose_Cfast$D))
plot(NULL, xlim=c(0,3.5), ylim=c(0,25), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)
abline(h=c(0))
rect(x-0.3, 0, x+0.3, y, col="orange3")
axis(side=1, at=c(0,1,2,3),cex=.8,mgp=c(2,.5,0),labels=NA)
axis(side=1, at=c(0,1,2,3),cex=.8,mgp=c(2,.5,0))
text(1.5,-3,expression(paste(" d"^"-1")),cex=.9,xpd=NA)
text(1.5,-4.5,"2% Glucose",cex=1,col="black",xpd=NA)

rect(0, -5.5, 11, -6.5, col=adjustcolor("orange", alpha.f=0.3),border=NA,xpd=NA)
abline(h=c(0),lwd=2)
arrows(0.55,mean(EXFglucose_Cslow$CN)-sd(EXFglucose_Cslow$CN),0.55,mean(EXFglucose_Cslow$CN)+sd(EXFglucose_Cslow$CN),code=3,length=0)
arrows(1.16,mean(EXFglucose_Cmed$CN)-sd(EXFglucose_Cmed$CN),1.16,mean(EXFglucose_Cmed$CN)+sd(EXFglucose_Cmed$CN),code=3,length=0)
arrows(2.37,mean(EXFglucose_Cfast$CN)-sd(EXFglucose_Cfast$CN),2.37,mean(EXFglucose_Cfast$CN)+sd(EXFglucose_Cfast$CN),code=3,length=0)
abline(v=3.5, lwd=1.5, lty=3)
#FY4 glucose
y=c(mean(FY4glucose_Cslow$CN),mean(FY4glucose_Cmed$CN),mean(FY4glucose_Cfast$CN))
x=c(mean(FY4glucose_Cslow$D),mean(FY4glucose_Cmed$D),mean(FY4glucose_Cfast$D))
plot(NULL, xlim=c(0,7), ylim=c(0,25), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)
abline(h=c(0))
rect(x-0.6, 0, x+0.6, y, col="orange")
abline(h=c(0),lwd=2)
axis(side=1, at=c(0,1,2,3,4,5,6),cex=.8,mgp=c(2,.5,0),labels=NA)
axis(side=1, at=c(0,3,6),cex=.8,mgp=c(2,.5,0))
text(4,-3,expression(paste(" d"^"-1")),cex=.9,xpd=NA)
text(3.25,-4.5,"  .75% Glucose",cex=1,col="black",xpd=NA)

text(3.5,-6,"Fermentation",cex=1,col="black",xpd=NA)
arrows(1.18,mean(FY4glucose_Cslow$CN)-sd(FY4glucose_Cslow$CN),1.18,mean(FY4glucose_Cslow$CN)+sd(FY4glucose_Cslow$CN),code=3,length=0)
arrows(3.36,mean(FY4glucose_Cmed$CN)-sd(FY4glucose_Cmed$CN),3.36,mean(FY4glucose_Cmed$CN)+sd(FY4glucose_Cmed$CN),code=3,length=0)
arrows(5.683333,mean(FY4glucose_Cfast$CN)-sd(FY4glucose_Cfast$CN),5.6833331,mean(FY4glucose_Cfast$CN)+sd(FY4glucose_Cfast$CN),code=3,length=0)

rect(0, 22.2, 15, 25.5, col=adjustcolor("gray", alpha.f=0.3),xpd=NA, border=NA)
text(7,25,"strain",cex=1,col="black",xpd=NA)
text(7,23.5,"FY4",cex=1,col="black",xpd=NA)
#FY4 galactose
y=c(mean(FY4galactose_Cslow$CN),mean(FY4galactose_Cmed$CN),mean(FY4galactose_Cfast$CN))
x=c(mean(FY4galactose_Cslow$D),mean(FY4galactose_Cmed$D),mean(FY4galactose_Cfast$D))
plot(NULL, xlim=c(0,7), ylim=c(0,25), bty="n",
     xaxt="n", xaxs="i", yaxs="i", xlab="", ylab="",axes=F)
abline(h=c(0))
rect(x-0.6, 0, x+0.6, y, col="gold",xpd=NA)
abline(h=c(0),lwd=2)
axis(side=1, at=c(0,1,2,3,4,5,6,7),cex=.8,mgp=c(2,.5,0),labels=NA)
axis(side=1, at=c(0,3,6),cex=.8,mgp=c(2,.5,0))
text(4,-3,expression(paste(" d"^"-1")),cex=.9,xpd=NA)
text(3.7,-4.5,".75% Galactose  ",cex=1,col="black",xpd=NA)
arrows(1.34,mean(FY4galactose_Cslow$CN)-sd(FY4galactose_Cslow$CN),1.34,mean(FY4galactose_Cslow$CN)+sd(FY4galactose_Cslow$CN),code=3,length=0)
arrows(3.59,mean(FY4galactose_Cmed$CN)-sd(FY4galactose_Cmed$CN),3.59,mean(FY4galactose_Cmed$CN)+sd(FY4galactose_Cmed$CN),code=3,length=0)
arrows(6.51,mean(FY4galactose_Cfast$CN)-sd(FY4galactose_Cfast$CN),6.51,mean(FY4galactose_Cfast$CN)+sd(FY4galactose_Cfast$CN),code=3,length=0)

dev.off()

