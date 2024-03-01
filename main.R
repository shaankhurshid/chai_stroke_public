# Depends
library(data.table)
library(stringr)
library(plyr)
library(pROC)
library(survival)
library(prodlim)
library(epiR)
library(sktools)

# Load stroke discrimination dataset
load(file='stroke_discrim_complete.RData')

# Remove people previously excluded for being in training set (N=1431 - 102 = 1329)
final_set <- fread(file='chai_stroke_set_121523.csv')
excluded <- stroke_discrim_complete[!(EMPI %in% final_set$EMPI)]
stroke_discrim_complete <- stroke_discrim_complete[EMPI %in% final_set$EMPI]

# Final output
## Load adjudications
adjudications <- fread(file='no_ecg_adjudicated_122523.csv')
setkey(adjudications,MGH_MRN); setkey(stroke_discrim_complete,MGH_MRN)
stroke_discrim_complete[adjudications,':='(stroke_date2 = as.Date(i.stroke_date,format='%Y-%m-%d'),
                                           changed_stroke_date = i.stroke_date_wrong,
                                           no_ecg_sk = i.no_ecg_sk)]
stroke_discrim_complete[,stroke_date := as.Date(ifelse(!is.na(stroke_date2),stroke_date2,stroke_date),origin='1970-01-01')]
write.csv(stroke_discrim_complete,file='stroke_set_011024.csv',row.names=F)

# Add BWH IDs
c3po_linker <- fread(file='mgh_bwh_mrn_empi_linker.txt')
bwh_linker <- c3po_linker[Hospital=='BWH']

# Add BWH IDs
setkey(stroke_discrim_complete,EMPI); setkey(bwh_linker,EMPI)
stroke_discrim_complete[bwh_linker,bwh_mrn := i.MRN]

# Recalculate CHARGE-AF to accommodate changes to stroke dates
stroke_discrim_complete[,':='(stroke_age = (as.numeric(stroke_date) - as.numeric(Date_of_Birth))/365.25)]
stroke_discrim_complete[,':='(charge_age = stroke_age/5)]
stroke_discrim_complete[,charge := 0.508*charge_age + charge_race*0.465 + 0.248*charge_height 
                        + 0.115*charge_weight + 0.197*charge_sbp + (-0.101)*charge_dbp + 0.359*charge_tobacco 
                        + 0.349*htn_med + 0.237*DM2 + 0.701*heart_failure + 0.496*mi]

# Load AF inferences
infer <- fread(file='chai_stroke_ecg_survival_curve_af_mgh_inference_v2024_01_12.tsv')
infer_bwh <- fread(file='chai_stroke_ecg_survival_curve_af_bwh_inference_v2023_12_20.tsv')

# Format and merge
### MGH
infer[,':='(infer_ecg_date = as.Date(datetime,format='%Y-%m-%d'),
            stroke_date = as.Date(stroke_date,format='%Y-%m-%d'))]
infer[,abs_diff_days := abs(infer_ecg_date-stroke_date)]
infer <- infer[af==FALSE]
setkey(infer,abs_diff_days)
infer <- infer[,.SD[which.min(abs_diff_days)],by='MGH_MRN']

setkey(infer,MGH_MRN); setkey(stroke_discrim_complete,MGH_MRN)
stroke_discrim_complete[infer,':='(survival_curve_af_prediction = i.survival_curve_af_prediction,
                                   infer_ecg_date = i.infer_ecg_date)]

### BWH
infer_bwh[,':='(infer_ecg_date = as.Date(datetime_y,format='%Y-%m-%d'),
                stroke_date = as.Date(datetime_x,format='%Y-%m-%d'))]
infer_bwh[,abs_diff_days := abs(infer_ecg_date-stroke_date)]
infer_bwh <- infer_bwh[af==FALSE]
setkey(infer_bwh,abs_diff_days)
infer_bwh <- infer_bwh[,.SD[which.min(abs_diff_days)],by='bwh_mrn']

setkey(infer_bwh,bwh_mrn); setkey(stroke_discrim_complete,bwh_mrn)
stroke_discrim_complete[infer_bwh,':='(survival_curve_af_prediction_bwh = i.survival_curve_af_prediction,
                                       infer_ecg_date_bwh = i.infer_ecg_date)]

# Combine
stroke_discrim_complete[,':='(infer_ecg_date_combined = ifelse(!is.na(infer_ecg_date),infer_ecg_date,infer_ecg_date_bwh),
                              survival_curve_af_prediction_combined = ifelse(!is.na(survival_curve_af_prediction),survival_curve_af_prediction,
                                                                             survival_curve_af_prediction_bwh))]

# Calculate components
stroke_discrim_complete[,':='(ecg_logit = log(survival_curve_af_prediction_combined/(1-survival_curve_af_prediction_combined)))]
stroke_discrim_complete[,':='(chai = 0.35655*charge + 0.44266*ecg_logit)]

# Remove missing ECG (1069 - 105 = 964)
stroke_discrim_complete_ecg <- stroke_discrim_complete[!is.na(ecg_logit)]

#################### PRIMARY
# Standardized variables
stroke_discrim_complete_ecg[,':='(charge_std = (charge - mean(charge))/sd(charge),
                                  ecg_logit_std = (ecg_logit - mean(ecg_logit))/sd(ecg_logit),
                                  chai_std = (chai - mean(chai))/sd(chai))]

# Metrics
### HR
charge_hr <- glm(def_cardioembolic~charge_std,data=stroke_discrim_complete_ecg,family = 'binomial')
charge_or <- c(exp(charge_hr$coefficients[2]),exp(confint(charge_hr)[2,1]),exp(confint(charge_hr)[2,2]))
ecg_hr <- glm(def_cardioembolic~ecg_logit_std,data=stroke_discrim_complete_ecg,family = 'binomial')
ecg_or <- c(exp(ecg_hr$coefficients[2]),exp(confint(ecg_hr)[2,1]),exp(confint(ecg_hr)[2,2]))
chai_hr <- glm(def_cardioembolic~chai_std,data=stroke_discrim_complete_ecg,family = 'binomial')
chai_or <- c(exp(chai_hr$coefficients[2]),exp(confint(chai_hr)[2,1]),exp(confint(chai_hr)[2,2]))

# OR PLOT
pdf('or_plot.pdf',height=5,width=7,
    pointsize=3)
par(oma=c(1,1,0,1))
par(mar=c(5,10,0,1))
plot(x=c(chai_or[1],ecg_or[1],charge_or[1]),
     y=seq(2.5,0.5,-1),xlim=c(0.9,4),ylim=c(0,3),log='x',
     xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col=c('#8da0cb','#66c2a5','#fc8d62'),cex=4,bty='n')
axis(1,cex.axis=2.5,at=seq(1,4,1),pos=-0.02)
axis(2,at=seq(2.5,0.5,-1),labels=FALSE,cex=2.5,pos=0.95)
mtext('Odds ratio for CE stroke',side=1,line=2.5,cex=2.5)
segments(c(chai_or[2],ecg_or[2],charge_or[2]),seq(2.5,0.5,-1),
         c(chai_or[3],ecg_or[3],charge_or[3]),seq(2.5,0.5,-1),col=c('#8da0cb','#66c2a5','#fc8d62'),
         lwd=2.2)
segments(1,0,1,2.5,col='black',lty=5)

plot_names <- c('CH-AI','ECG-AI','CHARGE-AF')
for (i in seq(2.5,0.5,-1)){
  n <- 3-i
  mtext(paste0(plot_names[i+0.5]),
        side=2,line=-2,las=2,cex=2.5,at=n)
}
dev.off()

### AUC
charge_roc <- roc(response=stroke_discrim_complete_ecg$def_cardioembolic,
                  predictor=stroke_discrim_complete_ecg$charge,plot=TRUE)
c(ci.auc(charge_roc)[1],ci.auc(charge_roc)[2],ci.auc(charge_roc)[3])
ecg_roc <- roc(response=stroke_discrim_complete_ecg$def_cardioembolic,
               predictor=stroke_discrim_complete_ecg$ecg_logit,plot=TRUE)
c(ci.auc(ecg_roc)[1],ci.auc(ecg_roc)[2],ci.auc(ecg_roc)[3])
chai_roc <- roc(response=stroke_discrim_complete_ecg$def_cardioembolic,
                predictor=stroke_discrim_complete_ecg$chai,plot=TRUE)
c(ci.auc(chai_roc)[1],ci.auc(chai_roc)[2],ci.auc(chai_roc)[3])

## Plotting
pdf(file='roc.pdf',height=5,width=5,
    pointsize=3)

# Plot settings
par(oma=c(1,1,1,1))
par(mar=c(4,4.8,1,1))
par(xpd=TRUE)

# Specifics
## Plot 1
plot(charge_roc,col='#fc8d62',lwd=1.4,axes=F,xlab='',ylab='')
plot(ecg_roc,add=T,col='#66c2a5',lwd=1.2,xlab='',ylab='')
plot(chai_roc,add=T,col='#8da0cb',lwd=1.2,xlab='',ylab='')

## Axes
axis(2,at=seq(1,0,-0.2),labels=c('1.0','0.8','0.6','0.4','0.2','0.0'),cex.axis=1.6,las=2)
axis(1,at=seq(0,1,0.2),cex.axis=1.6,las=1)

## Labels
title(xlab='1 - Specificity',line=2.6,cex.lab=1.8)
title(ylab='Sensitivity',line=3.4,cex.lab=1.8)

## Legend
legend(0.5,0.2,legend=c('CH-AI (0.736)','ECG-AI (0.721)','CHARGE-AF (0.662)'),
       col=c('#8da0cb','#66c2a5','#fc8d62'),
       lty=1,lwd=1,pch=1,bty='n',cex=1.5)

## Stop
dev.off()

# Density plots
# CHARGE stratified by CE
x <- list(v1=stroke_discrim_complete_ecg[def_cardioembolic==1]$charge_std,
          v2=stroke_discrim_complete_ecg[def_cardioembolic==0]$charge_std)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(-4,4,1),expand=c(0,0),limits=c(-4,4)) +
  scale_y_continuous(breaks=seq(0,0.6,0.1),expand=c(0,0),limits=c(0,0.6)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('CE','No CE')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=30,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=30,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=30),legend.text=element_text(size=30)) +
  labs(x='CHARGE-AF',y='Density') 
ggsave('density_charge.pdf',
       height=2,width=2.5,units='in',scale=4)

# ECG-AI stratified by CE
x <- list(v1=stroke_discrim_complete_ecg[def_cardioembolic==1]$ecg_logit_std,
          v2=stroke_discrim_complete_ecg[def_cardioembolic==0]$ecg_logit_std)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(-3,5,1),expand=c(0,0),limits=c(-3,5)) +
  scale_y_continuous(breaks=seq(0,0.6,0.1),expand=c(0,0),limits=c(0,0.6)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('CE','No CE')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=30,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=30,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=30),legend.text=element_text(size=30)) +
  labs(x='ECG-AI',y='Density') 
ggsave('density_ecgai.pdf',
       height=2,width=2.5,units='in',scale=4)

# CH-AI stratified by CE
x <- list(v1=stroke_discrim_complete_ecg[def_cardioembolic==1]$chai_std,
          v2=stroke_discrim_complete_ecg[def_cardioembolic==0]$chai_std)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(-4,4,1),expand=c(0,0),limits=c(-4,4)) +
  scale_y_continuous(breaks=seq(0,0.6,0.1),expand=c(0,0),limits=c(0,0.6)) +
  scale_fill_manual(values=c("#f03b20","#2b8cbe"),name='',labels=c('CE','No CE')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=30,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=30,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=30),legend.text=element_text(size=30)) +
  labs(x='CH-AI',y='Density')
ggsave('density_chai.pdf',
       height=2,width=2.5,units='in',scale=4)

##################### SECONDARY NO AF
# Add AF outcome and other ages
stroke_af <- fread(file='af_def_012224.csv')
death_last_enc <- fread(file='death_last_enc_011924.csv')
setkey(stroke_af,MGH_MRN); setkey(stroke_discrim_complete_ecg,MGH_MRN); setkey(death_last_enc,rpdr_mrn)
stroke_discrim_complete_ecg[stroke_af,age_at_af_mgh := i.age_at_af]
stroke_discrim_complete_ecg[death_last_enc,':='(death_age = as.numeric(i.Date_Of_Death_Age),
                                                last_enc = as.numeric(i.last_enc_age))]
stroke_discrim_complete_ecg[,last_enc := as.numeric(last_enc,origin='1970-01-01')]
setkey(stroke_af,bwh_mrn); setkey(stroke_discrim_complete_ecg,bwh_mrn)
stroke_discrim_complete_ecg[stroke_af,age_at_af_bwh := i.age_at_af]
stroke_discrim_complete_ecg[,age_at_af_combined := ifelse(!is.na(age_at_af_mgh),age_at_af_mgh,age_at_af_bwh)]
stroke_discrim_complete_ecg[,stroke_age := (as.numeric(stroke_date) - as.numeric(Date_of_Birth))/365.25]
stroke_discrim_complete_ecg[,stroke_age_blanked := ((as.numeric(stroke_date) - as.numeric(Date_of_Birth))+30)/365.25]
stroke_discrim_complete_ecg[,':='(has_af = ifelse(!is.na(age_at_af_combined),1,0))]
stroke_discrim_complete_ecg[,':='(time_to_af = ifelse(has_af==1,(age_at_af_combined - stroke_age_blanked),
                                                      pmin((last_enc - stroke_age_blanked),(death_age - stroke_age_blanked),na.rm=T)))]
stroke_discrim_complete_ecg[,incd_af := ifelse(has_af==1 & (time_to_af > 0),1,0)]
stroke_discrim_complete_ecg[,prev_af := ifelse(has_af==1 & (time_to_af <= 0),1,0)]

# Remove people with prevalent AF, AF during stroke admission, or AF within 30 days of stroke (N= 964 - 231 = 733)
stroke_discrim_complete_ecg_noaf <- stroke_discrim_complete_ecg[!is.na(subsequent_af) & subsequent_af==0 & prev_af==0]

# Restandardize
stroke_discrim_complete_ecg_noaf[,':='(charge_std = (charge - mean(charge))/sd(charge),
                                       ecg_logit_std = (ecg_logit - mean(ecg_logit))/sd(ecg_logit),
                                       chai_std = (chai - mean(chai))/sd(chai))]

# Metrics
### HR
charge_hr <- glm(def_cardioembolic~charge_std,data=stroke_discrim_complete_ecg_noaf,family = 'binomial')
c(exp(charge_hr$coefficients[2]),exp(confint(charge_hr)[2,1]),exp(confint(charge_hr)[2,2]))
ecg_hr <- glm(def_cardioembolic~ecg_logit_std,data=stroke_discrim_complete_ecg_noaf,family = 'binomial')
c(exp(ecg_hr$coefficients[2]),exp(confint(ecg_hr)[2,1]),exp(confint(ecg_hr)[2,2]))
chai_hr <- glm(def_cardioembolic~chai_std,data=stroke_discrim_complete_ecg_noaf,family = 'binomial')
c(exp(chai_hr$coefficients[2]),exp(confint(chai_hr)[2,1]),exp(confint(chai_hr)[2,2]))

### AUC
charge_roc <- roc(response=stroke_discrim_complete_ecg_noaf$def_cardioembolic,
                  predictor=stroke_discrim_complete_ecg_noaf$charge)
c(ci.auc(charge_roc)[1],ci.auc(charge_roc)[2],ci.auc(charge_roc)[3])
ecg_roc <- roc(response=stroke_discrim_complete_ecg_noaf$def_cardioembolic,
               predictor=stroke_discrim_complete_ecg_noaf$ecg_logit)
c(ci.auc(ecg_roc)[1],ci.auc(ecg_roc)[2],ci.auc(ecg_roc)[3])
chai_roc <- roc(response=stroke_discrim_complete_ecg_noaf$def_cardioembolic,
                predictor=stroke_discrim_complete_ecg_noaf$chai)
c(ci.auc(chai_roc)[1],ci.auc(chai_roc)[2],ci.auc(chai_roc)[3])

# Youden using ECG ROC
best_ecgai <- coords(ecg_roc, "best", best.method="youden")

##################### INCIDENT AF AMONG NON-CE
incident_af_set <- stroke_discrim_complete_ecg_noaf[!is.na(time_to_af) & (time_to_af > 0)]

incident_af_set[,':='(charge_std = (charge - mean(charge))/sd(charge),
                      ecg_logit_std = (ecg_logit - mean(ecg_logit))/sd(ecg_logit),
                      chai_std = (chai - mean(chai))/sd(chai))]

mod_af_ecgai <- coxph(Surv(time_to_af,incd_af) ~ ecg_logit_std,data=incident_af_set)

incident_af_set[,high_ecgai := ifelse(ecg_logit_std <= quantile(ecg_logit_std,0.9),1,0)]
incident_af_set[,high_ecgai_youden := ifelse(ecg_logit > best_ecgai$threshold,1,0)]

incident_af_set[,low_ecgai := ifelse(ecg_logit_std < quantile(ecg_logit_std,0.1),1,0)]
incident_af_set[,low_ecgai_youden := ifelse(ecg_logit <= best_ecgai$threshold,1,0)]

mod_af_ecgai_bin <- coxph(Surv(time_to_af,incd_af) ~ high_ecgai,data=incident_af_set)
mod_af_ecgai_bin_y <- coxph(Surv(time_to_af,incd_af) ~ high_ecgai_youden,data=incident_af_set)

mod_af_ecgai_bin_low <- coxph(Surv(time_to_af,incd_af) ~ low_ecgai,data=incident_af_set)
mod_af_ecgai_bin_low_y <- coxph(Surv(time_to_af,incd_af) ~ low_ecgai_youden,data=incident_af_set)

### High risk
prodlim_ecgai <- prodlim(Hist(time_to_af,incd_af)~high_ecgai,data=incident_af_set)
surv_10y_high <- sktools::survivor(data=incident_af_set,risk_data='high_ecgai',event='incd_af',time='time_to_af',breakpoint=10)

# KM
pdf(file='km_af.pdf',height=3.5,width=3.8,
    pointsize=3)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
par(xpd=TRUE)
plot(prodlim_ecgai,'cuminc',ylim=c(0,0.7),xlim=c(0,10),
     axis2.at=seq(0,0.7,0.1),axis2.las=2,lwd=1.4,background=F,
     axis1.at=seq(0,10,5),axis1.labels=as.character(seq(0,10,5)),
     atrisk.times=seq(0,10,5),col=c("#f03b20",'darkgray'),atrisk.col='black',confint=FALSE,
     legend.x=-0.5,legend.y=0.85,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,legend.title='',
     legend.cex=2,legend.legend=c("Top 10%","Bottom 90%"),legend.y.intersp=1,
     atrisk.title=("                     "),atrisk.pos=-0.8,atrisk.line=c(7,8.5),
     atrisk.cex=1.8,atrisk.interspace=0.8,xlab='',ylab='')
mtext("Cumulative risk of AF (%)",side=2,line=-1.2,at=0.3,cex=2.5)
mtext("Years",side=1, line=-2.5,cex=2.5)
mtext('Top 10%',side=1, line=-0,cex=1.8,at=-3)
mtext('Bottom 90%',side=1, line=1.5,cex=1.8,at=-3)
dev.off()

### Low risk
prodlim_ecgai_low <- prodlim(Hist(time_to_af,incd_af)~low_ecgai,data=incident_af_set)
surv_10y_low <- sktools::survivor(data=incident_af_set,risk_data='low_ecgai',event='incd_af',time='time_to_af',breakpoint=10)

# KM
pdf(file='km_af_low.pdf',height=3.5,width=3.8,
    pointsize=3)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
par(xpd=TRUE)
plot(prodlim_ecgai_low,'cuminc',ylim=c(0,0.7),xlim=c(0,10),
     axis2.at=seq(0,0.7,0.1),axis2.las=2,lwd=1.4,background=F,
     axis1.at=seq(0,10,5),axis1.labels=as.character(seq(0,10,5)),
     atrisk.times=seq(0,10,5),col=c("darkgray",'#66c2a5'),atrisk.col='black',confint=FALSE,
     legend.x=-0.5,legend.y=0.85,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,legend.title='',
     legend.cex=2,legend.legend=c("Top 90%","Bottom 10%"),legend.y.intersp=1,
     atrisk.title=("                     "),atrisk.pos=-0.8,atrisk.line=c(7,8.5),
     atrisk.cex=1.8,atrisk.interspace=0.8,xlab='',ylab='')
mtext("Cumulative risk of AF (%)",side=2,line=-1.2,at=0.3,cex=2.5)
mtext("Years",side=1, line=-2.5,cex=2.5)
mtext('Top 90%',side=1, line=-0,cex=1.8,at=-3)
mtext('Bottom 10%',side=1, line=1.5,cex=1.8,at=-3)
dev.off()

############# Test chars
# Function to make a 2x2 table
table2x2 <- function(data,disease,test,key){
  true_pos <- data[data[,disease]==1,]
  true_neg <- data[data[,disease]==0,]
  test_pos <- data[data[,test]==1,]
  test_neg <- data[data[,test]==0,]
  a <- nrow(true_pos[true_pos[,key] %in% test_pos[,key],])
  c <- nrow(true_pos[true_pos[,key] %in% test_neg[,key],])
  b <- nrow(true_neg[true_neg[,key] %in% test_pos[,key],])
  d <- nrow(true_neg[true_neg[,key] %in% test_neg[,key],])
  table <- matrix(c(a,c,b,d),ncol=2,nrow=2)
}
stroke_discrim_complete_ecg[,id := 1:nrow(stroke_discrim_complete_ecg)]
stroke_discrim_complete_ecg_noaf[,id := 1:nrow(stroke_discrim_complete_ecg_noaf)]

########## Overall
## CHARGE-AF
charge_roc <- roc(response=stroke_discrim_complete_ecg$def_cardioembolic,
                  predictor=stroke_discrim_complete_ecg$charge_std)
charge_coord70 <- coords(charge_roc, 0.10, 'threshold') 
charge_coord80 <- coords(charge_roc, -0.13, 'threshold') 
charge_coord90 <- coords(charge_roc, -0.63, 'threshold') 
charge_coord95 <- coords(charge_roc, -1.24, 'threshold') 
stroke_discrim_complete_ecg[,':='(charge_sens70 = ifelse(charge_std > 0.10,1,0),
                                  charge_sens80 = ifelse(charge_std > -0.13,1,0),
                                  charge_sens90 = ifelse(charge_std > -0.63,1,0),
                                  charge_sens95 = ifelse(charge_std > -1.24,1,0))]
charge_coord_spec_70 <- coords(charge_roc, 0.50, 'threshold') 
charge_coord_spec_80 <- coords(charge_roc, 0.77, 'threshold') 
charge_coord_spec_90 <- coords(charge_roc, 1.11, 'threshold') 
charge_coord_spec_95 <- coords(charge_roc, 1.38, 'threshold') 
stroke_discrim_complete_ecg[,':='(charge_spec70 = ifelse(charge_std > 0.50,1,0),
                                  charge_spec80 = ifelse(charge_std > 0.77,1,0),
                                  charge_spec90 = ifelse(charge_std > 1.11,1,0),
                                  charge_spec95 = ifelse(charge_std > 1.38,1,0))]
stroke_discrim_complete_ecg_df <- as.data.frame(stroke_discrim_complete_ecg)
sens70_charge <- table2x2(data=stroke_discrim_complete_ecg_df,
                          disease='def_cardioembolic',test='charge_sens70',key='id')
epi.tests(sens70_charge)
sens80_charge <- table2x2(data=stroke_discrim_complete_ecg_df,
                          disease='def_cardioembolic',test='charge_sens80',key='id')
epi.tests(sens80_charge)
sens90_charge <- table2x2(data=stroke_discrim_complete_ecg_df,
                          disease='def_cardioembolic',test='charge_sens90',key='id')
epi.tests(sens90_charge)
sens95_charge <- table2x2(data=stroke_discrim_complete_ecg_df,
                          disease='def_cardioembolic',test='charge_sens95',key='id')
epi.tests(sens95_charge)
spec70_charge <- table2x2(data=stroke_discrim_complete_ecg_df,
                          disease='def_cardioembolic',test='charge_spec70',key='id')
epi.tests(spec70_charge)
spec80_charge <- table2x2(data=stroke_discrim_complete_ecg_df,
                          disease='def_cardioembolic',test='charge_spec80',key='id')
epi.tests(spec80_charge)
spec90_charge <- table2x2(data=stroke_discrim_complete_ecg_df,
                          disease='def_cardioembolic',test='charge_spec90',key='id')
epi.tests(spec90_charge)
spec95_charge <- table2x2(data=stroke_discrim_complete_ecg_df,
                          disease='def_cardioembolic',test='charge_spec95',key='id')
epi.tests(spec95_charge)

## ECG-AI
ecgai_roc <- roc(response=stroke_discrim_complete_ecg$def_cardioembolic,
                 predictor=stroke_discrim_complete_ecg$ecg_logit_std)
ecgai_coord70 <- coords(ecgai_roc, 0.03, 'threshold') 
ecgai_coord80 <- coords(ecgai_roc, -0.25, 'threshold') 
ecgai_coord90 <- coords(ecgai_roc, -0.58, 'threshold') 
ecgai_coord95 <- coords(ecgai_roc, -0.86, 'threshold') 
stroke_discrim_complete_ecg[,':='(ecgai_sens70 = ifelse(ecg_logit_std > 0.03,1,0),
                                  ecgai_sens80 = ifelse(ecg_logit_std > -0.25,1,0),
                                  ecgai_sens90 = ifelse(ecg_logit_std > -0.58,1,0),
                                  ecgai_sens95 = ifelse(ecg_logit_std > -0.86,1,0))]
ecgai_coord_spec_70 <- coords(ecgai_roc, 0.15, 'threshold') 
ecgai_coord_spec_80 <- coords(ecgai_roc, 0.45, 'threshold') 
ecgai_coord_spec_90 <- coords(ecgai_roc, 0.83, 'threshold') 
ecgai_coord_spec_95 <- coords(ecgai_roc, 1.24, 'threshold') 
stroke_discrim_complete_ecg[,':='(ecgai_spec70 = ifelse(ecg_logit_std > 0.15,1,0),
                                  ecgai_spec80 = ifelse(ecg_logit_std > 0.45,1,0),
                                  ecgai_spec90 = ifelse(ecg_logit_std > 0.83,1,0),
                                  ecgai_spec95 = ifelse(ecg_logit_std > 1.24,1,0))]
stroke_discrim_complete_ecg_df <- as.data.frame(stroke_discrim_complete_ecg)
sens70_ecgai <- table2x2(data=stroke_discrim_complete_ecg_df,
                         disease='def_cardioembolic',test='ecgai_sens70',key='id')
epi.tests(sens70_ecgai)
sens80_ecgai <- table2x2(data=stroke_discrim_complete_ecg_df,
                         disease='def_cardioembolic',test='ecgai_sens80',key='id')
epi.tests(sens80_ecgai)
sens90_ecgai <- table2x2(data=stroke_discrim_complete_ecg_df,
                         disease='def_cardioembolic',test='ecgai_sens90',key='id')
epi.tests(sens90_ecgai)
sens95_ecgai <- table2x2(data=stroke_discrim_complete_ecg_df,
                         disease='def_cardioembolic',test='ecgai_sens95',key='id')
epi.tests(sens95_ecgai)
spec70_ecgai <- table2x2(data=stroke_discrim_complete_ecg_df,
                         disease='def_cardioembolic',test='ecgai_spec70',key='id')
epi.tests(spec70_ecgai)
spec80_ecgai <- table2x2(data=stroke_discrim_complete_ecg_df,
                         disease='def_cardioembolic',test='ecgai_spec80',key='id')
epi.tests(spec80_ecgai)
spec90_ecgai <- table2x2(data=stroke_discrim_complete_ecg_df,
                         disease='def_cardioembolic',test='ecgai_spec90',key='id')
epi.tests(spec90_ecgai)
spec95_ecgai <- table2x2(data=stroke_discrim_complete_ecg_df,
                         disease='def_cardioembolic',test='ecgai_spec95',key='id')
epi.tests(spec95_ecgai)

## CHAI
chai_roc <- roc(response=stroke_discrim_complete_ecg$def_cardioembolic,
                predictor=stroke_discrim_complete_ecg$chai_std)
chai_coord70 <- coords(chai_roc, 0.25, 'threshold') 
chai_coord80 <- coords(chai_roc, 0.02, 'threshold') 
chai_coord90 <- coords(chai_roc, -0.51, 'threshold') 
chai_coord95 <- coords(chai_roc, -0.88, 'threshold') 
stroke_discrim_complete_ecg[,':='(chai_sens70 = ifelse(chai_std > 0.25,1,0),
                                  chai_sens80 = ifelse(chai_std > 0.02,1,0),
                                  chai_sens90 = ifelse(chai_std > -0.51,1,0),
                                  chai_sens95 = ifelse(chai_std > -0.88,1,0))]
chai_coord_spec_70 <- coords(chai_roc, 0.3, 'threshold') 
chai_coord_spec_80 <- coords(chai_roc, 0.59, 'threshold') 
chai_coord_spec_90 <- coords(chai_roc, 0.97, 'threshold') 
chai_coord_spec_95 <- coords(chai_roc, 1.3, 'threshold') 
stroke_discrim_complete_ecg[,':='(chai_spec70 = ifelse(chai_std > 0.3,1,0),
                                  chai_spec80 = ifelse(chai_std > 0.59,1,0),
                                  chai_spec90 = ifelse(chai_std > 0.97,1,0),
                                  chai_spec95 = ifelse(chai_std > 1.3,1,0))]
stroke_discrim_complete_ecg_df <- as.data.frame(stroke_discrim_complete_ecg)
sens70_chai <- table2x2(data=stroke_discrim_complete_ecg_df,
                        disease='def_cardioembolic',test='chai_sens70',key='id')
epi.tests(sens70_chai)
sens80_chai <- table2x2(data=stroke_discrim_complete_ecg_df,
                        disease='def_cardioembolic',test='chai_sens80',key='id')
epi.tests(sens80_chai)
sens90_chai <- table2x2(data=stroke_discrim_complete_ecg_df,
                        disease='def_cardioembolic',test='chai_sens90',key='id')
epi.tests(sens90_chai)
sens95_chai <- table2x2(data=stroke_discrim_complete_ecg_df,
                        disease='def_cardioembolic',test='chai_sens95',key='id')
epi.tests(sens95_chai)
spec70_chai <- table2x2(data=stroke_discrim_complete_ecg_df,
                        disease='def_cardioembolic',test='chai_spec70',key='id')
epi.tests(spec70_chai)
spec80_chai <- table2x2(data=stroke_discrim_complete_ecg_df,
                        disease='def_cardioembolic',test='chai_spec80',key='id')
epi.tests(spec80_chai)
spec90_chai <- table2x2(data=stroke_discrim_complete_ecg_df,
                        disease='def_cardioembolic',test='chai_spec90',key='id')
epi.tests(spec90_chai)
spec95_chai <- table2x2(data=stroke_discrim_complete_ecg_df,
                        disease='def_cardioembolic',test='chai_spec95',key='id')
epi.tests(spec95_chai)

########## No AF subset
## CHARGE-AF
charge_roc <- roc(response=stroke_discrim_complete_ecg_noaf$def_cardioembolic,
                  predictor=stroke_discrim_complete_ecg_noaf$charge_std)
charge_coord70 <- coords(charge_roc, -0.04, 'threshold') 
charge_coord80 <- coords(charge_roc, -0.81, 'threshold') 
charge_coord90 <- coords(charge_roc, -1.19, 'threshold') 
charge_coord95 <- coords(charge_roc, -1.25, 'threshold') 
stroke_discrim_complete_ecg_noaf[,':='(charge_sens70 = ifelse(charge_std > -0.04,1,0),
                                       charge_sens80 = ifelse(charge_std > -0.81,1,0),
                                       charge_sens90 = ifelse(charge_std > -1.19,1,0),
                                       charge_sens95 = ifelse(charge_std > -1.25,1,0))]
charge_coord_spec_70 <- coords(charge_roc, 0.56, 'threshold') 
charge_coord_spec_80 <- coords(charge_roc, 0.88, 'threshold') 
charge_coord_spec_90 <- coords(charge_roc, 1.18, 'threshold') 
charge_coord_spec_95 <- coords(charge_roc, 1.43, 'threshold') 
stroke_discrim_complete_ecg_noaf[,':='(charge_spec70 = ifelse(charge_std > 0.56,1,0),
                                       charge_spec80 = ifelse(charge_std > 0.88,1,0),
                                       charge_spec90 = ifelse(charge_std > 1.18,1,0),
                                       charge_spec95 = ifelse(charge_std > 1.43,1,0))]
stroke_discrim_complete_ecg_noaf_df <- as.data.frame(stroke_discrim_complete_ecg_noaf)
sens70_charge <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                          disease='def_cardioembolic',test='charge_sens70',key='id')
epi.tests(sens70_charge)
sens80_charge <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                          disease='def_cardioembolic',test='charge_sens80',key='id')
epi.tests(sens80_charge)
sens90_charge <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                          disease='def_cardioembolic',test='charge_sens90',key='id')
epi.tests(sens90_charge)
sens95_charge <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                          disease='def_cardioembolic',test='charge_sens95',key='id')
epi.tests(sens95_charge)
spec70_charge <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                          disease='def_cardioembolic',test='charge_spec70',key='id')
epi.tests(spec70_charge)
spec80_charge <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                          disease='def_cardioembolic',test='charge_spec80',key='id')
epi.tests(spec80_charge)
spec90_charge <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                          disease='def_cardioembolic',test='charge_spec90',key='id')
epi.tests(spec90_charge)
spec95_charge <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                          disease='def_cardioembolic',test='charge_spec95',key='id')
epi.tests(spec95_charge)

## ECG-AI
ecgai_roc <- roc(response=stroke_discrim_complete_ecg_noaf$def_cardioembolic,
                 predictor=stroke_discrim_complete_ecg_noaf$ecg_logit_std)
ecgai_coord70 <- coords(ecgai_roc, -0.07, 'threshold') 
ecgai_coord80 <- coords(ecgai_roc, -0.38, 'threshold') 
ecgai_coord90 <- coords(ecgai_roc, -0.58, 'threshold') 
ecgai_coord95 <- coords(ecgai_roc, -1.49, 'threshold') 
stroke_discrim_complete_ecg_noaf[,':='(ecgai_sens70 = ifelse(ecg_logit_std > -0.07,1,0),
                                       ecgai_sens80 = ifelse(ecg_logit_std > -0.38,1,0),
                                       ecgai_sens90 = ifelse(ecg_logit_std > -0.58,1,0),
                                       ecgai_sens95 = ifelse(ecg_logit_std > -1.49,1,0))]
ecgai_coord_spec_70 <- coords(ecgai_roc, 0.42, 'threshold') 
ecgai_coord_spec_80 <- coords(ecgai_roc, 0.72, 'threshold') 
ecgai_coord_spec_90 <- coords(ecgai_roc, 1.2, 'threshold') 
ecgai_coord_spec_95 <- coords(ecgai_roc, 1.65, 'threshold') 
stroke_discrim_complete_ecg_noaf[,':='(ecgai_spec70 = ifelse(ecg_logit_std > 0.42,1,0),
                                       ecgai_spec80 = ifelse(ecg_logit_std > 0.72,1,0),
                                       ecgai_spec90 = ifelse(ecg_logit_std > 1.2,1,0),
                                       ecgai_spec95 = ifelse(ecg_logit_std > 1.65,1,0))]
stroke_discrim_complete_ecg_noaf_df <- as.data.frame(stroke_discrim_complete_ecg_noaf)
sens70_ecgai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                         disease='def_cardioembolic',test='ecgai_sens70',key='id')
epi.tests(sens70_ecgai)
sens80_ecgai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                         disease='def_cardioembolic',test='ecgai_sens80',key='id')
epi.tests(sens80_ecgai)
sens90_ecgai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                         disease='def_cardioembolic',test='ecgai_sens90',key='id')
epi.tests(sens90_ecgai)
sens95_ecgai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                         disease='def_cardioembolic',test='ecgai_sens95',key='id')
epi.tests(sens95_ecgai)
spec70_ecgai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                         disease='def_cardioembolic',test='ecgai_spec70',key='id')
epi.tests(spec70_ecgai)
spec80_ecgai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                         disease='def_cardioembolic',test='ecgai_spec80',key='id')
epi.tests(spec80_ecgai)
spec90_ecgai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                         disease='def_cardioembolic',test='ecgai_spec90',key='id')
epi.tests(spec90_ecgai)
spec95_ecgai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                         disease='def_cardioembolic',test='ecgai_spec95',key='id')
epi.tests(spec95_ecgai)

## CHAI
chai_roc <- roc(response=stroke_discrim_complete_ecg_noaf$def_cardioembolic,
                predictor=stroke_discrim_complete_ecg_noaf$chai_std)
chai_coord70 <- coords(chai_roc, 0.07, 'threshold') 
chai_coord80 <- coords(chai_roc, -0.32, 'threshold') 
chai_coord90 <- coords(chai_roc, -0.93, 'threshold') 
chai_coord95 <- coords(chai_roc, -1.58, 'threshold') 
stroke_discrim_complete_ecg_noaf[,':='(chai_sens70 = ifelse(chai_std > 0.07,1,0),
                                       chai_sens80 = ifelse(chai_std > -0.32,1,0),
                                       chai_sens90 = ifelse(chai_std > -0.93,1,0),
                                       chai_sens95 = ifelse(chai_std > -1.58,1,0))]
chai_coord_spec_70 <- coords(chai_roc, 0.53, 'threshold') 
chai_coord_spec_80 <- coords(chai_roc, 0.83, 'threshold') 
chai_coord_spec_90 <- coords(chai_roc, 1.22, 'threshold') 
chai_coord_spec_95 <- coords(chai_roc, 1.53, 'threshold') 
stroke_discrim_complete_ecg_noaf[,':='(chai_spec70 = ifelse(chai_std > 0.53,1,0),
                                       chai_spec80 = ifelse(chai_std > 0.83,1,0),
                                       chai_spec90 = ifelse(chai_std > 1.22,1,0),
                                       chai_spec95 = ifelse(chai_std > 1.53,1,0))]
stroke_discrim_complete_ecg_noaf_df <- as.data.frame(stroke_discrim_complete_ecg_noaf)
sens70_chai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                        disease='def_cardioembolic',test='chai_sens70',key='id')
epi.tests(sens70_chai)
sens80_chai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                        disease='def_cardioembolic',test='chai_sens80',key='id')
epi.tests(sens80_chai)
sens90_chai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                        disease='def_cardioembolic',test='chai_sens90',key='id')
epi.tests(sens90_chai)
sens95_chai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                        disease='def_cardioembolic',test='chai_sens95',key='id')
epi.tests(sens95_chai)
spec70_chai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                        disease='def_cardioembolic',test='chai_spec70',key='id')
epi.tests(spec70_chai)
spec80_chai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                        disease='def_cardioembolic',test='chai_spec80',key='id')
epi.tests(spec80_chai)
spec90_chai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                        disease='def_cardioembolic',test='chai_spec90',key='id')
epi.tests(spec90_chai)
spec95_chai <- table2x2(data=stroke_discrim_complete_ecg_noaf_df,
                        disease='def_cardioembolic',test='chai_spec95',key='id')
epi.tests(spec95_chai)