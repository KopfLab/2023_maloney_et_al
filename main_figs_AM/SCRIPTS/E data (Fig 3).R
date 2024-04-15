## Code for Enzyme plots for DOI: 10.1073/pnas.2310771121 by A.Maloney
library(tidyverse)
library("DataCombine")
library("readxl")
#data----
Edata = read_excel("DATA/ENZYME ASSAY SUMMARY.xlsx",sheet = "avg") 
allEdata = read_excel("DATA/ENZYME ASSAY SUMMARY.xlsx",sheet = "all") 
# load data back in from table_S1 excel sheet
table_s1_path <- "DATA/DATA_S1_all_data.xlsx"
growth_rates <- readxl::read_excel(table_s1_path, sheet = "growth_rates")
lipids_all <- readxl::read_excel(table_s1_path, sheet = "lipids_all")
enzymes_all <- readxl::read_excel(table_s1_path, sheet = "enzymes_all")
bulk_all <- readxl::read_excel(table_s1_path, sheet = "bulk_all_chemostat")

AEM04=filter(allEdata, Culture=="AEM04" )
AEM03=filter(allEdata, Culture=="AEM03" )
#FIG 3 MAIN ----
pdf("FIGURES/ENZYME Fig3 IDH, ALDH, oxPPP.pdf", width=4.6, height=4.4,encoding="MacRoman")
par(mar = c(1.1,.5,1.5,0), oma=c(2,5.5,1,10))
layout(matrix(c(1,2,3,4,5,5),nrow=3,ncol=2,byrow=T))
#g6p----
AEM03G=filter(AEM03, enzyme=="g6p")
G6pdh03=lm(Edata$G6PDH03*1000~Edata$D03); summary(G6pdh03)
plot(AEM03G$D,AEM03G$`V/BCA`*1000,axes=F ,ylim=c(0,.25*1000),xlim=c(0,3),cex=.2,cex.main=1,col="skyblue"); abline(G6pdh03,col="skyblue")
#axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1)); box()
text(1.5,.305*1000,"Glucose-6-phosphate ",xpd=NA)
text(1.5,.275*1000,"dehydrogenase (G6PDH)",xpd=NA)
axis(side=2,mgp=c(2,.5,0),at=c(0,100,200),las=1);box()
axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1),labels=NA)
#mtext(side=1,expression(paste("d"^"-1")),line=2)
points(Edata$D03,Edata$G6PDH03*1000,bg=adjustcolor("skyblue",alpha.f=.7) ,pch=21,cex=2); #abline(G6pdh03,col="skyblue",lty=3)
legend("topleft","",title="A",bty="n",x.intersp=.1,y.intersp=.1,text.font=2)

AEM04G=filter(AEM04, enzyme=="g6p")
G6pdh04=lm(Edata$G6PDH04*1000~Edata$D04); summary(G6pdh04)
points(AEM04G$D,AEM04G$`V/BCA`*1000,col="orange3",cex=.2 ); abline(G6pdh04,col="orange3")
points(Edata$D04,Edata$G6PDH04*1000,bg=adjustcolor("orange3",alpha.f=.7),pch=22 ,cex=2);# abline(G6pdh04,col="orange")

rp = vector('expression',2) #thanks https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(summary(G6pdh03)$adj.r.squared,digits=2)))[2]
rp[2] = substitute(expression(italic(R)^2 == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(summary(G6pdh04)$adj.r.squared,digits=2)))[2]
legend('bottomright', cex=.85,legend = rp, bty = 'n',text.col=c("skyblue","orange3"))

#6pg----
AEM03p=filter(AEM03, enzyme=="6pg")
pdh03=lm(Edata$p6pg03*1000~Edata$D03); summary(pdh03)
plot(AEM03p$D,AEM03p$`V/BCA`*1000,axes=F ,ylim=c(0,.25*1000),xlim=c(0,3),cex=.2,cex.main=1,col="skyblue"); abline(pdh03,col="skyblue")
text(1.5,.305*1000,"6-Phosphogluconate ",xpd=NA)
text(1.5,.275*1000,"dehydrogenase (6PGDH)",xpd=NA)
axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1),labels=NA)
points(Edata$D03,Edata$p6pg03*1000,bg=adjustcolor("skyblue",alpha.f=.7) ,pch=21,cex=2); #abline(pdh03,col="skyblue",lty=3)
legend("topleft","",title="B",bty="n",x.intersp=.1,y.intersp=.1,text.font=2)


AEM04p=filter(AEM04, enzyme=="6pg")
pdh04=lm(Edata$p6pg04*1000~Edata$D04); summary(pdh04)
points(AEM04p$D,AEM04p$`V/BCA`*1000,col="orange3" ,cex=.2); abline(pdh04,col="orange3")
points(Edata$D04,Edata$p6pg04*1000,bg=adjustcolor("orange3",alpha.f=.7) ,pch=22,cex=2); #abline(pdh04,col="red")
box()

rp = vector('expression',2) 
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(summary(pdh03)$adj.r.squared,digits=2)))[2]
rp[2] = substitute(expression(italic(R)^2 == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(summary(pdh04)$adj.r.squared,digits=2)))[2]
legend('topright', cex=.85,legend = rp, bty = 'n',text.col=c("skyblue","orange3"))


legend(3.2,.25*1000,title="strain EXF-4126",legend=c("Respiration","(2% Glycerol)", "Fermentation", "(2% Glucose)"), text.col=c("black","skyblue","black","orange3"),
       pch=c(21,NA,22,NA),pt.bg=c(adjustcolor("skyblue",alpha.f=.7),NA,adjustcolor("orange3",alpha.f=.7),NA),pt.cex=2,bty="n",xpd=NA)


#ALDH----
AEM03ALDH=filter(AEM03, enzyme=="ALDH 5mM" | enzyme=="ALDH 10mM")
ALDH03=lm(Edata$ALDH03*1000~Edata$D03); summary(ALDH03)
plot(AEM03ALDH$D,AEM03ALDH$`V/BCA`*1000,axes=F ,ylim=c(0,.25*1000),xlim=c(0,3),cex=.2,cex.main=1,col="skyblue"); abline(ALDH03,col="skyblue")
text(1.5,.305*1000,"Aldehyde ",xpd=NA)
text(1.5,.275*1000,"dehydrogenase (ALDH)",xpd=NA)
axis(side=2,at=c(0,100,200),mgp=c(2,.5,0),las=1);box()
axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1)); box()
#mtext(side=1,expression(paste("d"^"-1")),cex=.9,line=2)
points(Edata$D03,Edata$ALDH03*1000 ,pch=21,bg=adjustcolor("skyblue",alpha.f=.7),cex=2); #abline(ALDH03)
legend("topleft","",title="C",bty="n",x.intersp=.1,y.intersp=.1,text.font=2)
mtext(side=1,expression(paste("d"^"-1","")),xpd=NA,cex=.7, line=1.1)
AEM04ALDH=filter(AEM04, enzyme=="ALDH 5mM" | enzyme=="ALDH 10mM")
ALDH04=lm(Edata$ALDH04*1000~Edata$D04); summary(ALDH04)
points(AEM04ALDH$D,AEM04ALDH$`V/BCA`*1000,col="orange3",cex=.2); abline(ALDH04,col="orange3")
points(Edata$D04,Edata$ALDH04*1000 ,pch=22,bg=adjustcolor("orange3",alpha.f=.7),cex=2); #abline(ALDH04)

mtext(side=2,"Enzyme Activity",xpd=NA,srt=90,cex=1,line=4)
mtext(side=2,expression(paste("(","mol NADPH ","mg total protein"^"-1"," min"^"-1",")")),xpd=NA,srt=90,cex=.8,line=2.5)

rp = vector('expression',2) 
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(summary(ALDH03)$adj.r.squared,digits=2)))[2]
rp[2] = substitute(expression(italic(R)^2 == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(summary(ALDH04)$adj.r.squared,digits=2)))[2]
legend('topright', cex=.85,legend = rp, bty = 'n',text.col=c("skyblue","orange3"))
#IDH----
AEM03IDH=filter(AEM03, enzyme=="IDH")
IDH03=lm(Edata$IDH03*1000~Edata$D03); summary(IDH03)
plot(AEM03IDH$D,AEM03IDH$`V/BCA`*1000,axes=F ,ylim=c(0,.25*1000),xlim=c(0,3),cex=.2,cex.main=1,col="skyblue");  abline(IDH03,col="skyblue")
text(1.5,.305*1000,"Isocitrate ",xpd=NA)
text(1.5,.275*1000,"dehydrogenase (IDH)",xpd=NA)
axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1))
mtext(side=1,expression(paste("d"^"-1","")),xpd=NA,cex=.7, line=1.1)
box()
#mtext(side=1,expression(paste("d"^"-1")),line=2,cex=.9);box()
points(AEM03IDH$D,AEM03IDH$`V/BCA`*1000,col="skyblue" ,cex=.2); #abline(IDH03,col="green",lty=3)
points(Edata$D03,Edata$IDH03*1000,bg=adjustcolor("skyblue",alpha.f=.7),pch=21,cex=2 ); #abline(IDH03,col="green")
legend("topleft","",title="D",bty="n",x.intersp=.1,y.intersp=.1,text.font=2)
text(3.2,.14*1000, "cytosolic",col="skyblue",xpd=NA,adj=0)
text(3.2,.11*1000, "+ mitochondrial IDH",col="skyblue",xpd=NA,adj=0)
text(3.2,.05*1000, "mitochondrial IDH only",col="orange3",xpd=NA,adj=0)


AEM04IDH=filter(AEM04, enzyme=="IDH")
IDH04=lm(Edata$IDH04*1000~Edata$D04); summary(IDH04)
points(AEM04IDH$D,AEM04IDH$`V/BCA`*1000,col="orange3" ,cex=.2); abline(IDH04,col="orange3")
points(Edata$D04,Edata$IDH04*1000,bg=adjustcolor("orange3",alpha.f=.7),pch=22,cex=2 ); 
rp = vector('expression',2) 
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(summary(IDH03)$adj.r.squared,digits=2)))[2]
rp[2] = substitute(expression(italic(R)^2 == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(summary(IDH04)$adj.r.squared,digits=2)))[2]
legend('topright', cex=.85,legend = rp, bty = 'n',text.col=c("skyblue","orange3"))
#total----
table_s1_path <- "DATA/DATA_S1_all_data.xlsx"
table_s2_path <- "DATA/DATA_S2_data_summary.xlsx"
# get growth rates from S1
growth_rates <- readxl::read_excel(table_s1_path, sheet = "growth_rates")
# get the enzyme summary from S2
enzymes_sum <- readxl::read_excel(table_s2_path, sheet = "enzymes_strainEXF4126_summary")
# combine
enzymes_w_growth_rates <- enzymes_sum |>
  left_join(growth_rates, by = c("strain","dataset", "growth"))
allenzymes_w_growth_rates <- filter(enzymes_w_growth_rates,enzyme == "all")
all03 <- filter(allenzymes_w_growth_rates,dataset == "glycerol")
all04 <- filter(allenzymes_w_growth_rates,dataset == "glucose")

plot(all03$growth_rate.1_d,all03$rate_mean ,axes=F,ylim=c(0,500),xlim=c(0,3),cex=.2,cex.main=1,col="skyblue"); 
axis(side=2,mgp=c(2,.5,0),las=1);box()
axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1)); box()
lmALL=lm(all03$rate_mean ~ all03$growth_rate.1_d); summary(lmALL)
abline(lmALL,col="skyblue")
arrows(all03$growth_rate.1_d,all03$rate_mean-all03$rate_sd,all03$growth_rate.1_d,all03$rate_mean+all03$rate_sd,code=3,length=0)
points(all03$growth_rate.1_d,all03$rate_mean ,pch=21,bg=adjustcolor("skyblue",alpha.f=.7),cex=2); 

lmALL04=lm(all04$rate_mean ~ all04$growth_rate.1_d); summary(lmALL04)
abline(lmALL04,col="orange3")
arrows(all04$growth_rate.1_d,all04$rate_mean-all04$rate_sd,all04$growth_rate.1_d,all04$rate_mean+all04$rate_sd,code=3,length=0)
points(all04$growth_rate.1_d,all04$rate_mean ,pch=22,bg=adjustcolor("orange3",alpha.f=.7),cex=2); 


rp = vector('expression',2) #thanks https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(summary(lmALL)$adj.r.squared,digits=2)))[2]
rp[2] = substitute(expression(italic(R)^2 == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(summary(lmALL04)$adj.r.squared,digits=2)))[2]
legend('bottomright', cex=.85,legend = rp, bty = 'n',text.col=c("skyblue","orange3"))

text(1.5,.175*1000,"Total NADPH production ",xpd=NA)
text(1.5,.10*1000,"(G6PDH+6PGDH+ALDH+IDH)",xpd=NA)
legend("topleft","",title="E",bty="n",x.intersp=.1,y.intersp=.1,text.font=2)
mtext(side=1,expression(paste("Growth Rate (d"^"-1",")")),xpd=NA,cex=1,line=2); box()

dev.off()




#FIG SI enzyme (percent of total)----
pdf("FIGURES/ENZYME SI  IDH, ALDH, oxPPP.pdf", width=6, height=4.4,encoding="MacRoman")
par(mar = c(1,.5,1,0), oma=c(2.2,5.5,1,10))
layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T))
table_s1_path <- "DATA/DATA_S1_all_data.xlsx"
table_s2_path <- "DATA/Data_S2_data_summary.xlsx"
# get growth rates from S1
growth_rates <- readxl::read_excel(table_s1_path, sheet = "growth_rates")
# get the enzyme summary from S2
enzymes_sum <- readxl::read_excel(table_s2_path, sheet = "enzymes_strainEXF4126_summary")
# combine
enzymes_w_growth_rates <- enzymes_sum |>
  left_join(growth_rates, by = c("strain","dataset", "growth"))

allenzymes_w_growth_rates <- filter(enzymes_w_growth_rates,enzyme == "all")
all03 <- filter(allenzymes_w_growth_rates,dataset == "glycerol")
all04 <- filter(allenzymes_w_growth_rates,dataset == "glucose")

#g6p----
g6p_enzymes_w_growth_rates <- filter(enzymes_w_growth_rates,enzyme == "g6p")
percent_g6p=g6p_enzymes_w_growth_rates$rate_mean/allenzymes_w_growth_rates$rate_mean

g6p_enzymes_w_per <- cbind(g6p_enzymes_w_growth_rates,percent_g6p)
AEM03G <- filter(g6p_enzymes_w_per,dataset == "glycerol")
AEM04G <- filter(g6p_enzymes_w_per,dataset == "glucose")

G6pdh03=lm(AEM03G$percent_g6p~AEM03G$growth_rate.1_d); summary(G6pdh03)
plot(AEM03G$growth_rate.1_d,AEM03G$percent_g6p,axes=F ,ylim=c(0,1),xlim=c(0,3),cex=.2,cex.main=1,col="skyblue"); 
abline(G6pdh03,col="skyblue")
#axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1)); box()
#text(1.5,1.2,"Glucose-6-phosphate ",xpd=NA)
#text(1.5,1.1,"dehydrogenase (G6PDH)",xpd=NA)
axis(side=2,mgp=c(2,.5,0),at=c(0,.5,1),las=1);box()
axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1),labels=NA)
#mtext(side=1,expression(paste("d"^"-1")),line=2)
points(AEM03G$growth_rate.1_d,AEM03G$percent_g6p,bg=adjustcolor("skyblue",alpha.f=.8) ,pch=21,cex=2); #abline(G6pdh03,col="skyblue",lty=3)
legend("topleft","",title="A. G6PDH",bty="n",x.intersp=.1,y.intersp=.1,text.font=2)

G6pdh04=lm(AEM04G$percent_g6p~AEM04G$growth_rate.1_d); summary(G6pdh04)
points(AEM04G$growth_rate.1_d,AEM04G$percent_g6p,col="orange3",cex=.2 ); 
abline(G6pdh04,col="orange3")
points(AEM04G$growth_rate.1_d,AEM04G$percent_g6p,bg=adjustcolor("orange3",alpha.f=.8),pch=22 ,cex=2);# abline(G6pdh04,col="orange")
legend("bottomright", cex=.85,legend=c((paste("R2=",format(summary(G6pdh03)$adj.r.squared,digits=2))),
                                       (paste("R2=",format(summary(G6pdh04)$adj.r.squared,digits=2)))),text.col=c("skyblue","orange3"),bty="n")
#6pg----
enzymes_w_growth_rates_6pg <- filter(enzymes_w_growth_rates,enzyme == "6pg")
percent_6pg=enzymes_w_growth_rates_6pg$rate_mean/allenzymes_w_growth_rates$rate_mean

enzymes_w_per_6pg <- cbind(enzymes_w_growth_rates_6pg,percent_6pg)
AEM03p <- filter(enzymes_w_per_6pg,dataset == "glycerol")
AEM04p <- filter(enzymes_w_per_6pg,dataset == "glucose")

pdh03=lm(AEM03p$percent_6pg~AEM03p$growth_rate.1_d); summary(pdh03)
plot(AEM03p$growth_rate.1_d,AEM03p$percent_6pg,axes=F ,ylim=c(0,1),xlim=c(0,3),cex=.2,cex.main=1,col="skyblue"); 
     abline(pdh03,col="skyblue")
#text(1.5,1.2,"6-Phosphogluconate ",xpd=NA)
#text(1.5,1.1,"dehydrogenase (6PGDH)",xpd=NA)
axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1),labels=NA)
points(AEM03p$growth_rate.1_d,AEM03p$percent_6pg,bg=adjustcolor("skyblue",alpha.f=.8) ,pch=21,cex=2); #abline(pdh03,col="skyblue",lty=3)
legend("topleft","",title="B. 6PGDH",bty="n",x.intersp=.1,y.intersp=.1,text.font=2)

pdh04=lm(AEM04p$percent_6pg~AEM04p$growth_rate.1_d); summary(pdh04)
points(AEM04p$growth_rate.1_d,AEM04p$percent_6pg,col="orange3" ,cex=.2); abline(pdh04,col="orange3")
points(AEM04p$growth_rate.1_d,AEM04p$percent_6pg,bg=adjustcolor("orange3",alpha.f=.8) ,pch=22,cex=2); #abline(pdh04,col="red")
box()
legend("topright",cex=.85, legend=c((paste("R2=",format(summary(pdh03)$adj.r.squared,digits=2))),
                                    (paste("R2=",format(summary(pdh04)$adj.r.squared,digits=2)))),text.col=c("skyblue","orange3"),bty="n")
legend(3.2,1,title="strain EXF-4126",legend=c("Respiration","(2% Glycerol)", "Fermentation", "(2% Glucose)"), text.col=c("black","skyblue","black","orange3"),
       pch=c(21,NA,22,NA),pt.bg=c(adjustcolor("skyblue",alpha.f=.8),NA,adjustcolor("orange3",alpha.f=.8),NA),pt.cex=2,bty="n",xpd=NA)
#ALDH----
enzymes_w_growth_rates_aldh <- filter(enzymes_w_growth_rates,enzyme == "ALDH")
percent_aldh=enzymes_w_growth_rates_aldh$rate_mean/allenzymes_w_growth_rates$rate_mean

enzymes_w_per_aldh <- cbind(enzymes_w_growth_rates_aldh,percent_aldh)
AEM03ALDH <- filter(enzymes_w_per_aldh,dataset == "glycerol")
AEM04ALDH <- filter(enzymes_w_per_aldh,dataset == "glucose")

ALDH03=lm(AEM03ALDH$percent_aldh~AEM03ALDH$growth_rate.1_d); summary(ALDH03)
plot(AEM03ALDH$growth_rate.1_d,AEM03ALDH$percent_aldh,axes=F ,ylim=c(0,1),xlim=c(0,3),cex=.2,cex.main=1,col="skyblue"); abline(ALDH03,col="skyblue")
#text(1.5,1.2,"Aldehyde ",xpd=NA)
#text(1.5,1.1,"dehydrogenase (ALDH)",xpd=NA)
axis(side=2,at=c(0,.5,1),mgp=c(2,.5,0),las=1);box()
axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1)); box()
#mtext(side=1,expression(paste("d"^"-1")),cex=.9,line=2)
points(AEM03ALDH$growth_rate.1_d,AEM03ALDH$percent_aldh ,pch=21,bg=adjustcolor("skyblue",alpha.f=.8),cex=2); #abline(ALDH03)
legend("topleft","",title="C. ALDH",bty="n",x.intersp=.1,y.intersp=.1,text.font=2)
mtext(side=1,expression(paste("d"^"-1","")),xpd=NA,cex=.7, line=1.1)


ALDH04=lm(AEM04ALDH$percent_aldh~AEM04ALDH$growth_rate.1_d); summary(ALDH04)
points(AEM04ALDH$growth_rate.1_d,AEM04ALDH$percent_aldh,col="orange3",cex=.2); abline(ALDH04,col="orange3")
points(AEM04ALDH$growth_rate.1_d,AEM04ALDH$percent_aldh,pch=22,bg=adjustcolor("orange3",alpha.f=.8),cex=2); #abline(ALDH04)

text(-1.4,1.25,"Enzyme Activity",xpd=NA,srt=90,cex=1.5)
text(-1,1.25,expression(paste("(% total activity measured)")),xpd=NA,srt=90,cex=1)

legend("topright",cex=.85, legend=c((paste("R2=",format(summary(ALDH03)$adj.r.squared,digits=2))),
                                    (paste("R2=",format(summary(ALDH04)$adj.r.squared,digits=2)))),text.col=c("skyblue","orange3"),bty="n")
#IDH----
enzymes_w_growth_rates_idh<- filter(enzymes_w_growth_rates,enzyme == "IDH")
percent_idh=enzymes_w_growth_rates_idh$rate_mean/allenzymes_w_growth_rates$rate_mean

enzymes_w_per_idh <- cbind(enzymes_w_growth_rates_idh,percent_idh)
AEM03IDH <- filter(enzymes_w_per_idh,dataset == "glycerol")
AEM04IDH <- filter(enzymes_w_per_idh,dataset == "glucose")

IDH03=lm(AEM03IDH$percent_idh~AEM03IDH$growth_rate.1_d); summary(IDH03)
plot(AEM03IDH$growth_rate.1_d,AEM03IDH$percent_idh,axes=F ,ylim=c(0,1),xlim=c(0,3),cex=.2,cex.main=1,col="skyblue");  abline(IDH03,col="skyblue")
#text(1.5,1.2,"Isocitrate ",xpd=NA)
#text(1.5,1.1,"dehydrogenase (IDH)",xpd=NA)
axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1))
mtext(side=1,expression(paste("d"^"-1","")),xpd=NA,cex=.7, line=1.1)
box()
#mtext(side=1,expression(paste("d"^"-1")),line=2,cex=.9);box()
points(AEM03IDH$growth_rate.1_d,AEM03IDH$percent_idh,col="skyblue" ,cex=.2); #abline(IDH03,col="green",lty=3)
points(AEM03IDH$growth_rate.1_d,AEM03IDH$percent_idh,bg=adjustcolor("skyblue",alpha.f=.8),pch=21,cex=2 ); #abline(IDH03,col="green")
legend("topleft","",title="D. IDH",bty="n",x.intersp=.1,y.intersp=.1,text.font=2)


text(3.2,.5, "cytosolic IDH",col="skyblue",xpd=NA,adj=0)
text(3.2,.4, "+ mitochondrial IDH",col="skyblue",xpd=NA,adj=0)
text(3.2,.05, "mitochondrial IDH only",col="orange3",xpd=NA,adj=0)

arrows(3.2,.45,2,.35,code=2,length=.1,col="skyblue",xpd=NA)
arrows(3.2,0.05,2,.15,code=2,length=.1,col="orange3",xpd=NA)

IDH04=lm(AEM04IDH$percent_idh~AEM04IDH$growth_rate.1_d); summary(IDH04)
points(AEM04IDH$growth_rate.1_d,AEM04IDH$percent_idh,col="orange3" ,cex=.2); abline(IDH04,col="orange3")
points(AEM04IDH$growth_rate.1_d,AEM04IDH$percent_idh,bg=adjustcolor("orange3",alpha.f=.8),pch=22,cex=2 ); 
legend("topright", cex=.85,legend=c((paste("R2=",format(summary(IDH03)$adj.r.squared,digits=2))),
                                    (paste("R2=",format(summary(IDH04)$adj.r.squared,digits=2)))),text.col=c("skyblue","orange3"),bty="n")
text(-0.2,-.3,expression(paste("Growth Rate (d"^"-1",")")),xpd=NA,cex=1.1); box()
dev.off()
#FIG SI RELATIVE Enzyme activity glycerol/glucose----
pdf("FIGURES/ENZYME SI relaive enzyme.pdf", width=4, height=3,encoding="MacRoman")
par(mar = c(1.5,0,1,.6), oma=c(2,4,1,1))
par(mfrow=c(1,2))

Dboth=data.frame(Edata$D03,Edata$D04)
plot(rowMeans(Dboth),Edata$ALDH03/Edata$ALDH04,axes=F ,ylim=c(0,6),xlim=c(0,3),cex=.2,cex.main=.7);
mtext(side=1,expression(paste("Average dilution rate d"^"-1")),line=1.5)
mtext(side=2,"Relative Enzyme Activity",line=2.5,cex=.9)
mtext(side=2," (Glycerol/Glucose)",line=1.5,cex=.9)
axis(side=1,mgp=c(2,.5,0),at=seq(0,3,by=1))
axis(side=2,mgp=c(2,.5,0));box()

PE_aldh=(Edata$SD_aldh03/Edata$ALDH03 + Edata$SD_aldh04/Edata$ALDH04)*Edata$ALDH03/Edata$ALDH04
PE_idh=(Edata$SD_idh03/Edata$IDH03 + Edata$SD_idh04/Edata$IDH04)*Edata$IDH03/Edata$IDH04
PE_g6pdh=(Edata$SD_g03/Edata$G6PDH03 + Edata$SD_g04/Edata$G6PDH04)*Edata$G6PDH03/Edata$G6PDH04
PE_p6=(Edata$SD_603/Edata$p6pg03 + Edata$SD_604/Edata$p6pg04)*Edata$p6pg03/Edata$p6pg04


arrows(rowMeans(Dboth),(Edata$ALDH03/Edata$ALDH04)-PE_aldh,rowMeans(Dboth),(Edata$ALDH03/Edata$ALDH04)+PE_aldh,code=3,length=0)
arrows(rowMeans(Dboth),(Edata$IDH03/Edata$IDH04)-PE_idh,rowMeans(Dboth),(Edata$IDH03/Edata$IDH04)+PE_idh,code=3,length=0, col="black")
arrows(rowMeans(Dboth),(Edata$G6PDH03/Edata$G6PDH04)-PE_g6pdh,rowMeans(Dboth),(Edata$G6PDH03/Edata$G6PDH04)+PE_g6pdh,code=3,length=0, col="gray60")
arrows(rowMeans(Dboth),(Edata$p6pg03/Edata$p6pg04)-PE_p6,rowMeans(Dboth),(Edata$p6pg03/Edata$p6pg04)+PE_p6,code=3,length=0, col="gray80")

points(rowMeans(Dboth),Edata$ALDH03/Edata$ALDH04,bg="white" ,pch=21,cex=1); #abline(ALDH03)
points(rowMeans(Dboth),Edata$IDH03/Edata$IDH04,bg="black",pch=21 ,cex=1); #abline(IDH03,col="green")
points(rowMeans(Dboth),Edata$G6PDH03/Edata$G6PDH04,bg="gray60" ,pch=21,cex=1); #abline(G6pdh03,col="orange")
points(rowMeans(Dboth),Edata$p6pg03/Edata$p6pg04,bg="gray80" ,pch=21,cex=1); #abline(pdh03,col="red")
legend(4,5, pch=21, legend=c("G6PDH","6PGDH","ALDH","IDH"),pt.bg=c("gray60","gray80",NA,"black"),
       cex=.7,xpd=NA)

dev.off()
