library(dplyr)
library(ggplot2)
library(stringr)
library(survival)
library(survminer)
library(reshape2)
library(ggforce)
library(ggrepel)

rm(list=ls())

GeneralTheme <- theme(axis.title = element_text(color = "black", size = 20),
                      axis.text = element_text(color = "black", size = 18),
                      strip.text = element_text(size = 18, color = "black"),
                      axis.ticks = element_line(color = "black"),
                      legend.title = element_text(color = "black", size = 20),
                      legend.text = element_text(color = "black", size = 18),
                      legend.position = "top")

# Table 1 -----------------------------------------------------------------

for(j in c("sex", "smk_na", "Anemia", "Thrombocytopenia", "prev_cancer_yn", "het_yn")){
  for(i in c(1, 2)){
    
    print(j)
    print(i)
    print(table(PhenoDF$AnyLCHIP_init_filt, PhenoDF[,(which(colnames(PhenoDF)==j))])[i,])
    print(prop.table(table(PhenoDF$AnyLCHIP_init_filt, PhenoDF[,(which(colnames(PhenoDF)==j))])[i,])*100)
    print(fisher.test(PhenoDF$AnyLCHIP_init_filt, PhenoDF[,(which(colnames(PhenoDF)==j))])[[1]])
    
  }
}

by(PhenoDF$age, PhenoDF$AnyLCHIP_init_filt, quantile)
wilcox.test(PhenoDF$age ~ PhenoDF$AnyLCHIP_init_filt)


# Table 2 -----------------------------------------------------------------

for(j in c("sex", "smk_na", "Anemia", "Thrombocytopenia", "prev_cancer_yn", "het_yn")){
  for(i in c(1, 2, 3, 4)){
    
    print(j)
    print(i)
    print(table(PhenoDF$LMCHIP_Comp, PhenoDF[,(which(colnames(PhenoDF)==j))])[i,])
    print(prop.table(table(PhenoDF$LMCHIP_Comp, PhenoDF[,(which(colnames(PhenoDF)==j))])[i,])*100)
    
  }
}

int_df <- PhenoDF%>%
  select(c("LMCHIP_Comp", "sex", "smk_na", "Anemia", "Thrombocytopenia", "prev_cancer_yn", "het_yn", "age"))
int_df1 <- int_df[int_df$LMCHIP_Comp %in% c("OnlyLCHIP", "Neither"),]
int_df2 <- int_df[int_df$LMCHIP_Comp %in% c("OnlyMCHIP", "Neither"),]
int_df3 <- int_df[int_df$LMCHIP_Comp %in% c("OnlyLCHIP", "OnlyMCHIP"),]

for(j in c("sex", "smk_na", "Anemia", "Thrombocytopenia", "prev_cancer_yn", "het_yn")){
  
  print(j)
  print(fisher.test(as.character(int_df1$LMCHIP_Comp), int_df1[,(which(colnames(int_df1)==j))])[[1]])
  print(fisher.test(as.character(int_df2$LMCHIP_Comp), int_df2[,(which(colnames(int_df2)==j))])[[1]])
  print(fisher.test(as.character(int_df3$LMCHIP_Comp), int_df3[,(which(colnames(int_df3)==j))])[[1]])
  
}

by(PhenoDF$age, PhenoDF$LMCHIP_Comp, quantile)
pairwise.wilcox.test(int_df$age, int_df$LMCHIP_Comp, p.adjust.method = "none")

# Figure 1 ----------------------------------------------------------------

UKB_Freq_Tab <- LCHIP_vars%>%
  distinct(SampID, Gene.refGene)%>%
  pull(Gene.refGene)%>%
  table()%>%
  data.frame()%>%
  rename("Var1" = ".")%>%
  arrange(-Freq)

#Figure 1A
UKB_Freq_Tab%>%
  head(n=30)%>%
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq))+
  geom_col()+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black", size = 18, angle = 45, hjust = 1))+
  GeneralTheme+
  xlab("Gene")+
  ylab("Number of individuals")

#Figure 1B
table(LCHIP_vars$SampID)%>%
  data.frame()%>%
  mutate(FreqN = ifelse(Freq>3, 4, Freq))%>%
  pull(FreqN)%>%
  table()%>%
  data.frame()%>%
  rename("Var1" = ".")%>%
  mutate(Var1 = as.numeric(as.character(Var1)))%>%
  ggplot(aes(x = Freq, y = reorder(Var1, -Var1)))+
  geom_col()+
  theme_classic()+
  GeneralTheme+
  ylab("Number of L-CHIP mutations")+
  xlab("Number of individuals")

#Figure 1C
AgeDFPrevLCHIP <- table(cut(PhenoDF$age, c(35, 45, 50, 55, 60, 65, 100)), PhenoDF$AnyLCHIP_init_filt)%>%
  data.frame()%>%
  group_by(Var1)%>%
  mutate(s = sum(Freq))%>%
  ungroup()%>%
  mutate(P = Freq/s*100)%>%
  filter(Var2 == TRUE)

AgeDFPrevLCHIP%>%
  mutate(Var1 = factor(Var1, levels = c("(35,45]", "(45,50]", "(50,55]", "(55,60]", "(60,65]", "(65,100]"),
                       labels = c("<45", "46-50", "51-55", "56-60", "61-65", ">65")))%>%
  ggplot(aes(x = Var1, y = P, group = Var2))+
  geom_line(color = "#E64B35")+
  geom_point(size = 3, color = "#E64B35")+
  theme_classic()+
  GeneralTheme+
  theme(legend.position = "none",
        strip.background = element_blank())+
  scale_y_continuous(breaks = c(2.25, 2.50, 2.75, 3.00),
                     labels = paste0(c(2.25, 2.50, 2.75, 3.00), "%"))+
  ylab("L-CHIP Prevalence")+
  xlab("Age group")



# Figure 2 ----------------------------------------------------------------

table(PhenoDF$CLL_CENSOR)
PhenoDF%>%
  filter(CLL_CENSOR==1)%>%
  pull(CLL_follow)%>%
  quantile()

int_df <- PhenoDF

surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)

fit.surv <- survfit(surv_object ~ het_yn, data = int_df)
#Figure2A
ggsurvplot(fit.surv, data = int_df, 
           #pval = TRUE, 
           censor = FALSE,
           fun = "event",
           palette= c("#595959FF", "#00A087FF"),
           risk.table.col="strata",
           risk.table.y.text=FALSE,
           #surv.median.line = "hv",
           break.time.by=5,
           surv.scale="percent",
           ylab="% Individuals",
           xlab="Years",
           risk.table = TRUE,
           #ggtheme=theme_classic(axis.title = element_text(color = "black")),
           ggtheme=theme(axis.text = element_text(color = "black", size = 18),
                         axis.title = element_text(color = "black", size = 20),
                         axis.line = element_line(color = "black"),
                         axis.ticks = element_blank(),
                         panel.background = element_blank(),
                         panel.grid = element_blank())
           #ylim = c(0, 1)
)

fit.coxph <- coxph(surv_object ~ as.numeric(het_yn) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

#Figure2B

fit.coxph <- coxph(surv_object ~ count_het + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

fit.coxph <- coxph(surv_object ~ mMSS + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

fit.coxph <- coxph(surv_object ~ mMSS+count_het + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

F2BTP <- readxl::read_xlsx("SupplementaryDataFile.xlsx", sheet = "Figure 2", skip = 5)

F2BTP%>%
  mutate(VarGrp = ifelse(grepl("Not adjusted", `Additional adjustment`), "Not adjusted", "Mutually adjusted"))%>%
  ggplot(aes(x = `exp(coef)`, y = paste(`Independent variable`, VarGrp, sep = "\n")))+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_point()+
  geom_errorbar(aes(xmin = `lower .95`, xmax = `upper .95`), width = 0)+
  theme_classic()+
  GeneralTheme+
  scale_x_continuous(breaks = c(1, 2, 4, 6))+
  facet_grid(`Independent variable`~., scales = "free_y")


# Figure 3; Supplementary Table 8 ----------------------------------------------------------------

int_df <- PhenoDF

surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)

fit.coxph <- coxph(surv_object ~ as.numeric(AnyLCHIP_init_filt)+as.numeric(AnyCHIP_init_filt) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

fit.coxph <- coxph(surv_object ~ as.numeric(AnyLCHIP_init_filt) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

int_df <- PhenoDF

int_df$GOI <- ifelse(int_df$AnyLCHIP_init_filt_VAF10, 1,
                     ifelse(int_df$AnyLCHIP_init_filt==FALSE, 0, NA))

surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)

fit.coxph <- coxph(surv_object ~ GOI + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)



int_df <- PhenoDF

int_df$GOI <- ifelse(int_df$NOM_LCHIP>1, 1,
                     ifelse(int_df$AnyLCHIP_init_filt==FALSE, 0, NA))

surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)

fit.coxph <- coxph(surv_object ~ GOI + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)


int_df <- PhenoDF

int_df$GOI <- ifelse(int_df$LCHIP_CLL_HR, 1,
                     ifelse(int_df$AnyLCHIP_init_filt==FALSE, 0, NA))

surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)

fit.coxph <- coxph(surv_object ~ GOI + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

LCHIPgvect <- names(sort(table(LCHIP_vars$Gene.refGene[LCHIP_vars$SampID %in% (PhenoDF$id[(PhenoDF$CLL_CENSOR == 1)])])))


int_df <- PhenoDF%>%
  filter(AnyLCHIP_init_filt==TRUE)

Tabfindf <- data.frame()
for(i in LCHIPgvect){
  
  int_df$GOI <- int_df$SampID %in% (LCHIP_vars$SampID[LCHIP_vars$Gene.refGene == i])
  
  
  surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)
  
  fit.coxph <- coxph(surv_object ~ as.numeric(GOI) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                     data = int_df)
  
  Tabtmpdf <- as.data.frame(summary(fit.coxph)$coefficients)[1,]
  Tabtmpdf$GOI <- i
  
  Tabtmpdf <- cbind(Tabtmpdf, as.data.frame(exp(confint(fit.coxph)))[1,])
  Tabfindf <- rbind(Tabfindf, Tabtmpdf)
  
  print(i)  
}


Tabfindf <- (LCHIP_vars[LCHIP_vars$SampID %in% (PhenoDF$SampID[(PhenoDF$CLL_CENSOR == 1)]),])%>%
  distinct(SampID, Gene.refGene)%>%
  pull(Gene.refGene)%>%
  table()%>%
  data.frame()%>%
  rename("Var1" = ".")%>%
  left_join(Tabfindf, ., by = c("GOI" = "Var1"))

#Supplementary Table 8
Tabfindf <- (LCHIP_vars)%>%
  distinct(SampID, Gene.refGene)%>%
  pull(Gene.refGene)%>%
  table()%>%
  data.frame()%>%
  rename("Var1" = ".")%>%
  rename("FreqTot" = "Freq")%>%
  left_join(Tabfindf, ., by = c("GOI" = "Var1"))









#A
LCRCLL <- readxl::read_xlsx("SupplementaryDataFile.xlsx", sheet = "Figure 3", range = "A2:H6")

LCRCLL%>%
  mutate(`Independent variable` = str_remove(`Independent variable`, " vs.*"))%>%
  mutate(Variable = factor(`Independent variable`, levels = rev(c(`Independent variable`))))%>%
  ggplot(aes(x = log10(`exp(coef)`), y = Variable))+
  geom_point()+
  geom_vline(xintercept = 0, linetype = 2)+
  geom_errorbar(aes(xmin = log10(`lower .95`), xmax = log10(`upper .95`)), width = 0)+
  theme_classic()+
  GeneralTheme+
  scale_x_continuous(breaks = c(log10(1), log10(10), log10(100)),
                     labels = c(1, 10, 100))+
  xlab("HR 95% CI")+
  ylab("")

#B
Tabfindf%>%
  filter(Freq>1)%>%
  ggplot(aes(x = log10(`exp(coef)`), y = reorder(GOI, `exp(coef)`)))+
  geom_vline(xintercept = log10(1), linetype = 2)+
  geom_point()+
  geom_errorbar(aes(xmin = log10(`2.5 %`), xmax = log10(`97.5 %`)), width = 0)+
  theme_classic()+
  GeneralTheme+
  scale_x_continuous(breaks = c(log10(1), log10(10), log10(100), log10(1000)),
                     labels = c(1, 10, 100, 1000))+
  xlab("HR 95% CI")+
  ylab("")

# Figure 4 ----------------------------------------------------------------

# Figure 4B; Supplementary Table 10 -------------------------------

FreqGenesToAssess <- c("MYD88")

CHIPfiltch <- c()
Genech <- c()
hvarch <- c()
ORch <- c()
SEch <- c()
HetCHch <- c()
TotCHch <- c()
PercHet <- c()
pch <- c()



for(i in c("AnyLCHIP_init_filt", "AnyLCHIP_init_filt_Small10", "AnyLCHIP_init_filt_VAF10", "NOMsingle", "NOMmulti", "HRG", "nonHRG")){
  
  int_df <- PhenoDF
  
  int_df$HCov0 <- int_df$count_het>0
  
  int_CHIP <- LCHIP_vars
  
  int_CHIP$Gene <- int_CHIP$Gene.refGene
  
  
  if(grepl("VAF10", i)){
    int_CHIP <- int_CHIP[int_CHIP$VAF>=10,]
    int_df <- int_df[(int_df$AnyLCHIP_init_filt_VAF10 == TRUE)|(int_df$AnyLCHIP_init_filt == FALSE),]
  }
  
  
  if(i == "AnyLCHIP_init_filt_Small10"){
    int_CHIP <- int_CHIP[int_CHIP$VAF<10,]
    int_df <- int_df[(int_df$AnyLCHIP_init_filt_Small10 == TRUE)|(int_df$AnyLCHIP_init_filt == FALSE),]
  }
  
  if(i == "NOMmulti"){
    int_df <- int_df[(int_df$NOM_LCHIP > 1)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  if(i == "NOMsingle"){
    int_df <- int_df[(int_df$NOM_LCHIP == 1)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  if(i == "HRG"){
    int_df <- int_df[(int_df$LCHIP_CLL_HR==TRUE)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  if(i == "nonHRG"){
    int_df <- int_df[(int_df$LCHIP_CLL_HR==FALSE)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  

  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID[int_CHIP$Gene == g]
    
    
    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID
    }
    
    
    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]
      
      
      ForLOGDF <- int_df%>%filter(count_het==0 | het_VOI)%>%filter(AnyLCHIP_init_filt == 0| IsGOI)
      
      
      si_ct <- with(ForLOGDF, table(IsGOI, het_VOI))%>%
        data.frame()%>%
        filter(IsGOI == TRUE)%>%
        mutate(s = sum(Freq))%>%
        filter(het_VOI == TRUE)%>%
        mutate(P = Freq/s*100)
      
      
      
      Logistic_multi <- glm(het_VOI ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
      si_logistic_sum <- summary(Logistic_multi)
      
      CHIPfiltch <- c(CHIPfiltch, i)
      Genech <- c(Genech, g)
      hvarch <- c(hvarch, hvar)
      
      if(length(si_ct$IsGOI)==1){
        
        HetCHch <- c(HetCHch, si_ct$Freq)
        TotCHch <- c(TotCHch, si_ct$s)
        PercHet <- c(PercHet, si_ct$P)
        
      }else{
        
        HetCHch <- c(HetCHch, NA)
        TotCHch <- c(TotCHch, NA)
        PercHet <- c(PercHet, NA)
        
      }
      
      ORch <- c(ORch, (exp(Logistic_multi$coefficients)[[2]]))
      SEch <- c(SEch, (exp(si_logistic_sum$coefficients)[2,2]))
      pch <- c(pch, (si_logistic_sum$coefficients[2,4]))
      
      rm(Logistic_multi)
      rm(si_logistic_sum)
      
    }
  }
}


si_ct <- with(PhenoDF, table((AnyLCHIP_init_filt==FALSE), (count_het>0)))%>%
  data.frame()%>%
  filter(Var1 == TRUE)%>%
  mutate(s = sum(Freq))%>%
  filter(Var2 == TRUE)%>%
  mutate(P = Freq/s*100)




CHIPfiltch <- c(CHIPfiltch, "AnyLCHIP_init_filt_Small")
Genech <- c(Genech, "No CHIP")
hvarch <- c(hvarch, "HCov0")
HetCHch <- c(HetCHch, si_ct$Freq)
TotCHch <- c(TotCHch, si_ct$s)
PercHet <- c(PercHet, si_ct$P)
ORch <- c(ORch, NA)
SEch <- c(SEch, NA)
pch <- c(pch, NA)




FreqGeneCHIPhet_df <- data.frame(CHIPfiltch, Genech, hvarch, ORch, SEch, HetCHch, TotCHch, PercHet, pch)

FreqGeneCHIPhet_df <- FreqGeneCHIPhet_df%>%
  mutate(PlotGF = paste(CHIPfiltch, Genech))%>%
  mutate(PlotGF = factor(PlotGF, levels = rev(c("AnyLCHIP_init_filt_Small No CHIP", "AnyLCHIP_init_filt any", "AnyLCHIP_init_filt_Small10 any", "AnyLCHIP_init_filt_VAF10 any", "NOMsingle any", "NOMmulti any", "nonHRG any", "HRG any")),
                         labels = rev(c("No CHIP", "CHIP", "2%<VAF<10%", "VAF>10%", "Single mutation", "Multiple mutations", "No high-risk gene", "High-risk gene"))))%>%
  filter(!(is.na(PlotGF)))

FreqGeneCHIPhet_df$SplVar <- case_when(FreqGeneCHIPhet_df$PlotGF %in% c("No CHIP", "CHIP") ~ "wSC",
                                       FreqGeneCHIPhet_df$PlotGF %in% c("2%<VAF<10%", "VAF>10%") ~ "xSC",
                                       FreqGeneCHIPhet_df$PlotGF %in% c("Single mutation", "Multiple mutations") ~ "ySC",
                                       TRUE ~ "zG")


FreqGeneCHIPhet_df$Significant <- ifelse(FreqGeneCHIPhet_df$PlotGF %in% c("No CHIP"), "Not significant", "Significant")

FreqGeneCHIPhet_df%>%
  ggplot(aes(y = PlotGF, x = PercHet, fill = Significant))+
  geom_col(position = position_dodge())+
  theme_classic()+
  GeneralTheme+
  facet_grid(SplVar~., scales = "free_y", space = "free_y")+
  scale_fill_manual(values = c("#595959FF", "#E64B35FF"))+
  scale_x_continuous(breaks = 0:10*10,
                     labels = paste0(0:10*10, "%"))+
  geom_text(aes(x = 1, y = PlotGF, label = HetCHch), color = "white", hjust = 0, size = 6)+
  ylab("Gene")+
  xlab("Heteroplasmy prevalence")


ForLOGDF <- PhenoDF

Logistic_multi <- glm(het_yn ~ as.numeric(AnyLCHIP_init_filt) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)

table(ForLOGDF$AnyLCHIP_init_filt, ForLOGDF$het_yn)
prop.table(table(ForLOGDF$AnyLCHIP_init_filt, ForLOGDF$het_yn)[1,])*100
prop.table(table(ForLOGDF$AnyLCHIP_init_filt, ForLOGDF$het_yn)[2,])*100

ForLOGDF <- PhenoDF%>%
  filter(AnyLCHIP_init_filt == TRUE)

Logistic_multi <- glm(het_yn ~ as.numeric(AnyLCHIP_init_filt_VAF10) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$AnyLCHIP_init_filt_VAF10, ForLOGDF$het_yn)
prop.table(table(ForLOGDF$AnyLCHIP_init_filt_VAF10, ForLOGDF$het_yn)[1,])*100
prop.table(table(ForLOGDF$AnyLCHIP_init_filt_VAF10, ForLOGDF$het_yn)[2,])*100


Logistic_multi <- glm(het_yn ~ as.numeric(NOM_LCHIP>1) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$NOM_LCHIP>1, ForLOGDF$het_yn)
prop.table(table(ForLOGDF$NOM_LCHIP>1, ForLOGDF$het_yn)[1,])*100
prop.table(table(ForLOGDF$NOM_LCHIP>1, ForLOGDF$het_yn)[2,])*100

Logistic_multi <- glm(het_yn ~ as.numeric(LCHIP_CLL_HR) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$LCHIP_CLL_HR, ForLOGDF$het_yn)
prop.table(table(ForLOGDF$LCHIP_CLL_HR, ForLOGDF$het_yn)[1,])*100
prop.table(table(ForLOGDF$LCHIP_CLL_HR, ForLOGDF$het_yn)[2,])*100


Logistic_multi <- glm(het_yn ~ as.numeric(LCHIP_CLL_HR)+as.numeric(AnyLCHIP_init_filt_VAF10) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$LCHIP_CLL_HR, ForLOGDF$het_yn)
prop.table(table(ForLOGDF$LCHIP_CLL_HR, ForLOGDF$het_yn)[1,])*100
prop.table(table(ForLOGDF$LCHIP_CLL_HR, ForLOGDF$het_yn)[2,])*100

#Supplementary Table 10

ForLOGDF <- PhenoDF%>%
  filter(AnyCHIP_init_filt==FALSE)

Logistic_multi <- glm(het_yn ~ as.numeric(AnyLCHIP_init_filt) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$AnyLCHIP_init_filt)
table(ForLOGDF$AnyLCHIP_init_filt, ForLOGDF$het_yn)
prop.table(table(ForLOGDF$AnyLCHIP_init_filt, ForLOGDF$het_yn)[1,])*100
prop.table(table(ForLOGDF$AnyLCHIP_init_filt, ForLOGDF$het_yn)[2,])*100

ForLOGDF <- PhenoDF%>%
  filter(AnyLCHIP_init_filt == TRUE)%>%
  filter(AnyCHIP_init_filt==FALSE)

Logistic_multi <- glm(het_yn ~ as.numeric(AnyLCHIP_init_filt_VAF10) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$AnyLCHIP_init_filt_VAF10)
table(ForLOGDF$AnyLCHIP_init_filt_VAF10, ForLOGDF$het_yn)
prop.table(table(ForLOGDF$AnyLCHIP_init_filt_VAF10, ForLOGDF$het_yn)[1,])*100
prop.table(table(ForLOGDF$AnyLCHIP_init_filt_VAF10, ForLOGDF$het_yn)[2,])*100


Logistic_multi <- glm(het_yn ~ as.numeric(NOM_LCHIP>1) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$NOM_LCHIP>1)
table(ForLOGDF$NOM_LCHIP>1, ForLOGDF$het_yn)
prop.table(table(ForLOGDF$NOM_LCHIP>1, ForLOGDF$het_yn)[1,])*100
prop.table(table(ForLOGDF$NOM_LCHIP>1, ForLOGDF$het_yn)[2,])*100

Logistic_multi <- glm(het_yn ~ as.numeric(LCHIP_CLL_HR) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$LCHIP_CLL_HR)
table(ForLOGDF$LCHIP_CLL_HR, ForLOGDF$het_yn)
prop.table(table(ForLOGDF$LCHIP_CLL_HR, ForLOGDF$het_yn)[1,])*100
prop.table(table(ForLOGDF$LCHIP_CLL_HR, ForLOGDF$het_yn)[2,])*100


# Figure 4C; Supplementary Table 10 -------------------------------

FreqGenesToAssess <- c("MYD88")

CHIPfiltch <- c()
Genech <- c()
hvarch <- c()
ORch <- c()
SEch <- c()
HetCHch <- c()
TotCHch <- c()
PercHet <- c()
pch <- c()

for(i in c("AnyLCHIP_init_filt", "AnyLCHIP_init_filt_Small10", "AnyLCHIP_init_filt_VAF10", "NOMsingle", "NOMmulti", "nonHRG", "HRG")){
  
  int_df <- PhenoDF%>%
    filter(count_het>0)
  
  int_df$HCov0 <- int_df$count_het>1
  int_df$MSSov0 <- int_df$MSS>0
  int_df$MSSov010 <- int_df$MSS>0.10
  int_df$MSSov050 <- int_df$MSS>0.50
  
  int_CHIP <- LCHIP_vars

  int_CHIP$Gene <- int_CHIP$Gene.refGene
  
  
  if(grepl("VAF10", i)){
    int_CHIP <- int_CHIP[int_CHIP$VAF>=10,]
    int_df <- int_df[(int_df$AnyLCHIP_init_filt_VAF10 == TRUE)|(int_df$AnyLCHIP_init_filt == FALSE),]
  }
  
  
  if(i == "AnyLCHIP_init_filt_Small10"){
    int_CHIP <- int_CHIP[int_CHIP$VAF<10,]
    int_df <- int_df[(int_df$AnyLCHIP_init_filt_Small10 == TRUE)|(int_df$AnyLCHIP_init_filt == FALSE),]
  }
  
  if(i == "NOMmulti"){
    int_df <- int_df[(int_df$NOM_LCHIP > 1)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  if(i == "NOMsingle"){
    int_df <- int_df[(int_df$NOM_LCHIP == 1)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  if(i == "HRG"){
    int_df <- int_df[(int_df$LCHIP_CLL_HR==TRUE)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  if(i == "nonHRG"){
    int_df <- int_df[(int_df$LCHIP_CLL_HR==FALSE)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID[int_CHIP$Gene == g]
    
    
    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID
    }
    table(int_df$IsGOI, int_df$AnyLCHIP_init_filt_Small10)
    
    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]
      
      ForLOGDF <- int_df%>%filter(count_het==1 | het_VOI)%>%filter(AnyLCHIP_init_filt == 0| IsGOI)
      
      
      si_ct <- with(ForLOGDF, table(IsGOI, het_VOI))%>%
        data.frame()%>%
        filter(IsGOI == TRUE)%>%
        mutate(s = sum(Freq))%>%
        filter(het_VOI == TRUE)%>%
        mutate(P = Freq/s*100)
      
      
      
      Logistic_multi <- glm(het_VOI ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
      si_logistic_sum <- summary(Logistic_multi)
      
      CHIPfiltch <- c(CHIPfiltch, i)
      Genech <- c(Genech, g)
      hvarch <- c(hvarch, hvar)
      
      if(length(si_ct$IsGOI)==1){
        
        HetCHch <- c(HetCHch, si_ct$Freq)
        TotCHch <- c(TotCHch, si_ct$s)
        PercHet <- c(PercHet, si_ct$P)
        
      }else{
        
        HetCHch <- c(HetCHch, NA)
        TotCHch <- c(TotCHch, NA)
        PercHet <- c(PercHet, NA)
        
      }
      
      ORch <- c(ORch, (exp(Logistic_multi$coefficients)[[2]]))
      SEch <- c(SEch, (exp(si_logistic_sum$coefficients)[2,2]))
      pch <- c(pch, (si_logistic_sum$coefficients[2,4]))
      
      rm(Logistic_multi)
      rm(si_logistic_sum)
      
    }
  }
}


si_ct <- with(PhenoDF%>%filter(count_het>0), table((AnyLCHIP_init_filt==FALSE), (count_het>1)))%>%
  data.frame()%>%
  filter(Var1 == TRUE)%>%
  mutate(s = sum(Freq))%>%
  filter(Var2 == TRUE)%>%
  mutate(P = Freq/s*100)




CHIPfiltch <- c(CHIPfiltch, "AnyLCHIP_init_filt_Small")
Genech <- c(Genech, "No CHIP")
hvarch <- c(hvarch, "HCov0")
HetCHch <- c(HetCHch, si_ct$Freq)
TotCHch <- c(TotCHch, si_ct$s)
PercHet <- c(PercHet, si_ct$P)
ORch <- c(ORch, NA)
SEch <- c(SEch, NA)
pch <- c(pch, NA)



FreqGeneCHIPhet_df <- data.frame(CHIPfiltch, Genech, hvarch, ORch, SEch, HetCHch, TotCHch, PercHet, pch)


FreqGeneCHIPhet_df <- FreqGeneCHIPhet_df%>%
  mutate(PlotGF = paste(CHIPfiltch, Genech))%>%
  mutate(PlotGF = factor(PlotGF, levels = rev(c("AnyLCHIP_init_filt_Small No CHIP", "AnyLCHIP_init_filt any", "AnyLCHIP_init_filt_Small10 any", "AnyLCHIP_init_filt_VAF10 any", "NOMsingle any", "NOMmulti any", "nonHRG any", "HRG any")),
                         labels = rev(c("No CHIP", "CHIP", "2%<VAF<10%", "VAF>10%", "Single mutation", "Multiple mutations", "No high-risk gene", "High-risk gene"))))%>%
  filter(!(is.na(PlotGF)))

FreqGeneCHIPhet_df$SplVar <- case_when(FreqGeneCHIPhet_df$PlotGF %in% c("No CHIP", "CHIP") ~ "wSC",
                                       FreqGeneCHIPhet_df$PlotGF %in% c("2%<VAF<10%", "VAF>10%") ~ "xSC",
                                       FreqGeneCHIPhet_df$PlotGF %in% c("Single mutation", "Multiple mutations") ~ "ySC",
                                       TRUE ~ "zG")

FreqGeneCHIPhet_df$Significant <- ifelse(FreqGeneCHIPhet_df$PlotGF %in% c("No CHIP"), "Not significant", "Significant")

FreqGeneCHIPhet_df%>%
  ggplot(aes(y = PlotGF, x = PercHet, fill = Significant))+
  geom_col(position = position_dodge())+
  theme_classic()+
  GeneralTheme+
  facet_grid(SplVar~., scales = "free_y", space = "free_y")+
  scale_fill_manual(values = c("#595959FF", "#E64B35FF"))+
  scale_x_continuous(breaks = 0:10*10,
                     labels = paste0(0:10*10, "%"))+
  geom_text(aes(x = 1, y = PlotGF, label = HetCHch), color = "white", hjust = 0, size = 6)+
  ylab("Gene")+
  xlab("Prevalence of multiple heteroplasmies")








ForLOGDF <- PhenoDF%>%
  filter(count_het>0)



Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(AnyLCHIP_init_filt) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt)
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt)[,1])*100
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt)[,2])*100

ForLOGDF <- PhenoDF%>%
  filter(AnyLCHIP_init_filt == TRUE & count_het>0)



Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(AnyLCHIP_init_filt_VAF10) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt_VAF10)
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt_VAF10)[,1])*100
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt_VAF10)[,2])*100

Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(NOM_LCHIP>1) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$count_het>1, ForLOGDF$NOM_LCHIP>1)
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$NOM_LCHIP>1)[,1])*100
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$NOM_LCHIP>1)[,2])*100

Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(LCHIP_CLL_HR) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$count_het>1, ForLOGDF$LCHIP_CLL_HR)
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$LCHIP_CLL_HR)[,1])*100
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$LCHIP_CLL_HR)[,2])*100

#Supplementary Table 10
ForLOGDF <- PhenoDF%>%
  filter(count_het>0)%>%
  filter(AnyCHIP_init_filt==FALSE)

Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(AnyLCHIP_init_filt) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt)
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt)[,1])*100
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt)[,2])*100

ForLOGDF <- PhenoDF%>%
  filter(AnyLCHIP_init_filt == TRUE & count_het>0)%>%
  filter(AnyCHIP_init_filt==FALSE)

Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(AnyLCHIP_init_filt_VAF10) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt_VAF10)
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt_VAF10)[,1])*100
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$AnyLCHIP_init_filt_VAF10)[,2])*100

Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(NOM_LCHIP>1) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$count_het>1, ForLOGDF$NOM_LCHIP>1)
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$NOM_LCHIP>1)[,1])*100
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$NOM_LCHIP>1)[,2])*100

Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(LCHIP_CLL_HR) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)
table(ForLOGDF$count_het>1, ForLOGDF$LCHIP_CLL_HR)
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$LCHIP_CLL_HR)[,1])*100
prop.table(table(ForLOGDF$count_het>1, ForLOGDF$LCHIP_CLL_HR)[,2])*100

# Figure 4D; Supplementary Table 9; Supplementary Table 10 ------------------------------------------

FreqGenesToAssess <- c("MYD88")

UKB_LIN_DF <- data.frame()

for(i in c("AnyLCHIP_init_filt", "AnyLCHIP_init_filt_Small10", "AnyLCHIP_init_filt_VAF10", "NOMsingle", "NOMmulti", "nonHRG", "HRG")){
  
  int_df <- PhenoDF
  int_df$LogmMSS <- log10(int_df$mMSS+1)
  int_df$HCov0 <- int_df$count_het>0
  
  int_CHIP <- LCHIP_vars
  int_CHIP$VAF <- int_CHIP$newVAF
  int_CHIP$Gene <- int_CHIP$Gene.refGene
  
  
  if(grepl("VAF10", i)){
    int_CHIP <- int_CHIP[int_CHIP$VAF>=10,]
    int_df <- int_df[(int_df$AnyLCHIP_init_filt_VAF10 == TRUE)|(int_df$AnyLCHIP_init_filt == FALSE),]
  }
  
  
  if(i == "AnyLCHIP_init_filt_Small10"){
    int_CHIP <- int_CHIP[int_CHIP$VAF<10,]
    int_df <- int_df[(int_df$AnyLCHIP_init_filt_Small10 == TRUE)|(int_df$AnyLCHIP_init_filt == FALSE),]
  }
  
  if(i == "NOMmulti"){
    int_df <- int_df[(int_df$NOM_LCHIP > 1)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  if(i == "NOMsingle"){
    int_df <- int_df[(int_df$NOM_LCHIP == 1)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  if(i == "HRG"){
    int_df <- int_df[(int_df$LCHIP_CLL_HR==TRUE)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  if(i == "nonHRG"){
    int_df <- int_df[(int_df$LCHIP_CLL_HR==FALSE)|(int_df$AnyLCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }

  

  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID[int_CHIP$Gene == g]

    
    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID
    }
    
    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]
      
      
      
      ForLOGDF <- int_df%>%filter(het_VOI)%>%filter(AnyLCHIP_init_filt == 0| IsGOI)
      
      
      ToBindDf <- ForLOGDF%>%
        filter(IsGOI)%>%
        select(mMSS)
      
      if(length(ToBindDf$mMSS)==0){
        next
      }
      
      Linear_multi <- lm(LogmMSS ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, data = ForLOGDF)
      si_linear_sum <- summary(Linear_multi)
      
      ToBindDf$CHIPfiltch <- i
      ToBindDf$Genech <- g
      ToBindDf$hvarch <- hvar
      ToBindDf$pch <- (si_linear_sum$coefficients[2,4])
      
      UKB_LIN_DF <- rbind(UKB_LIN_DF, ToBindDf)
      
      rm(Linear_multi)
      rm(si_linear_sum)
      
    }
  }
}



UKB_LIN_DF <- PhenoDF%>%
  filter(AnyLCHIP_init_filt==FALSE & count_het>0)%>%
  select(mMSS)%>%
  mutate(CHIPfiltch = "AnyLCHIP_init_filt_Small", Genech = "No CHIP", hvarch = "HCov0", pch = NA)%>%
  rbind(UKB_LIN_DF, .)





UKB_LIN_DF <- UKB_LIN_DF%>%
  mutate(PlotGF = paste(CHIPfiltch, Genech))%>%
  mutate(PlotGF = factor(PlotGF, levels = rev(c("AnyLCHIP_init_filt_Small No CHIP", "AnyLCHIP_init_filt any", "AnyLCHIP_init_filt_Small10 any", "AnyLCHIP_init_filt_VAF10 any", "NOMsingle any", "NOMmulti any", "nonHRG any", "HRG any")),
                         labels = rev(c("No CHIP", "CHIP", "2%<VAF<10%", "VAF>10%", "Single mutation", "Multiple mutations", "No high-risk gene", "High-risk gene"))))%>%
  filter(!(is.na(PlotGF)))

UKB_LIN_DF$SplVar <- case_when(UKB_LIN_DF$PlotGF %in% c("No CHIP", "CHIP") ~ "wSC",
                               UKB_LIN_DF$PlotGF %in% c("2%<VAF<10%", "VAF>10%") ~ "xSC",
                               UKB_LIN_DF$PlotGF %in% c("Single mutation", "Multiple mutations") ~ "ySC",
                               TRUE ~ "zG")


UKB_LIN_DF$Significant <- ifelse(UKB_LIN_DF$PlotGF %in% c("No CHIP"), "Not significant", "Significant")





UKB_LIN_DF%>%
  ggplot(aes(x = PlotGF, y = mMSS, fill = Significant))+
  geom_violin(scale = "width")+
  geom_boxplot(alpha = 0)+
  theme_classic()+
  GeneralTheme+
  facet_grid(SplVar~., scales = "free_y", space = "free_y")+
  scale_fill_manual(values = c("#595959FF", "#E64B35FF"))+
  ylab("mMSS")+
  xlab("Gene")+
  coord_flip(ylim = c(0, 1))


by(UKB_LIN_DF$mMSS, UKB_LIN_DF$PlotGF, quantile)

#Supplementary Table 9
ForLOGDF <- PhenoDF%>%filter(count_het>0)

ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)
Linear_multi <- lm(LogmMSS ~ as.numeric(AnyLCHIP_init_filt) + rms::rcs(age, df = 4)+count_het + sex + smk_na + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$AnyLCHIP_init_filt, quantile)


ForLOGDF <- PhenoDF%>%filter(count_het>0)%>%filter(AnyLCHIP_init_filt == TRUE)
ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)

Linear_multi <- lm(LogmMSS ~ as.numeric(AnyLCHIP_init_filt_VAF10) + rms::rcs(age, df = 4)+count_het + sex + smk_na + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$AnyLCHIP_init_filt_VAF10, quantile)

Linear_multi <- lm(LogmMSS ~ as.numeric(NOM_LCHIP>1) + rms::rcs(age, df = 4)+count_het + sex + smk_na + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$NOM_LCHIP>1, quantile)

Linear_multi <- lm(LogmMSS ~ as.numeric(LCHIP_CLL_HR) + rms::rcs(age, df = 4)+count_het + sex + smk_na + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$LCHIP_CLL_HR, quantile)


#Supplementary Table 10

ForLOGDF <- PhenoDF%>%filter(count_het>0)%>%
  filter(AnyCHIP_init_filt==FALSE)

ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)
Linear_multi <- lm(LogmMSS ~ as.numeric(AnyLCHIP_init_filt) + rms::rcs(age, df = 4)+count_het + sex + smk_na + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$AnyLCHIP_init_filt, quantile)

ForLOGDF <- PhenoDF%>%filter(count_het>0)%>%filter(AnyLCHIP_init_filt == TRUE)%>%
  filter(AnyCHIP_init_filt==FALSE)
ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)

Linear_multi <- lm(LogmMSS ~ as.numeric(AnyLCHIP_init_filt_VAF10) + rms::rcs(age, df = 4)+count_het + sex + smk_na + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$AnyLCHIP_init_filt_VAF10, quantile)

Linear_multi <- lm(LogmMSS ~ as.numeric(NOM_LCHIP>1) + rms::rcs(age, df = 4)+count_het + sex + smk_na + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$NOM_LCHIP>1, quantile)

Linear_multi <- lm(LogmMSS ~ as.numeric(LCHIP_CLL_HR) + rms::rcs(age, df = 4)+count_het + sex + smk_na + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$LCHIP_CLL_HR, quantile)


# Figure 5 ----------------------------------------------------------------

int_df <- PhenoDF

j <- "count_het"
m <- "AnyLCHIP_init_filt"

int_df$BindingHetVAR <- case_when(j == "count_het" ~ (int_df$count_het>0))

int_df <- int_df[((int_df$AnyLCHIP_init_filt==0)|(int_df[,which(colnames(int_df) == m)])),]
int_df <- int_df[((int_df$BindingHetVAR)|(int_df$count_het==0)),]

int_df$AnyCHIP_VAR <- (int_df[,which(colnames(int_df) == m)])
int_df$CHIP_VAR <- paste0((int_df[,which(colnames(int_df) == m)]),
                          (int_df$BindingHetVAR))

int_df$CHIP_VAR <- factor(int_df$CHIP_VAR, levels = c("FALSEFALSE", "FALSETRUE", "TRUEFALSE", "TRUETRUE"),
                          labels = c("NO_CHIP_NO_Het", "NO_CHIP_YES_Het", "YES_CHIP_NO_Het", "YES_CHIP_YES_Het"))

surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)

fit.surv <- survfit(surv_object ~ CHIP_VAR, data = int_df)

#Figure 5B
ggsurvplot(fit.surv, data = int_df, 
           censor = FALSE,
           fun = "event",
           palette= c("#595959FF", "#00A087FF", "#4DBBD5FF", "#E64B35FF"),
           risk.table.col="strata",
           risk.table.y.text=FALSE,
           #surv.median.line = "hv",
           break.time.by=5,
           surv.scale="percent",
           ylab="% Individuals",
           xlab="Years",
           risk.table = TRUE,
           ggtheme=theme(axis.text = element_text(color = "black", size = 18),
                         axis.title = element_text(color = "black", size = 20),
                         axis.line = element_line(color = "black"),
                         axis.ticks = element_blank(),
                         panel.background = element_blank(),
                         panel.grid = element_blank())
)

fit.coxph <- coxph(surv_object ~ CHIP_VAR + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)

exp(fit.coxph[[1]])
exp(confint(fit.coxph))
summary(fit.coxph)



#Figure 5C
int_df <- PhenoDF%>%
  filter(AnyLCHIP_init_filt==TRUE)

int_df$MultiLCHIP <- int_df$NOM_LCHIP>1

surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)


fit.coxph <- coxph(surv_object ~ count_het+AnyLCHIP_init_filt_VAF10+MultiLCHIP+LCHIP_CLL_HR+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

fit.coxph <- coxph(surv_object ~ mMSS+AnyLCHIP_init_filt_VAF10+MultiLCHIP+LCHIP_CLL_HR+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

fit.coxph <- coxph(surv_object ~ mMSS+count_het+AnyLCHIP_init_filt_VAF10+MultiLCHIP+LCHIP_CLL_HR+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)






SF3BTP <- readxl::read_xlsx("SupplementaryDataFile.xlsx", sheet = "Figure 5", range = "A14:I18")

SF3BTP%>%
  mutate(VarGrp = ifelse(grepl("Not adjusted", `Additional adjustment`), "Not adjusted", "Mutually adjusted"))%>%
  ggplot(aes(x = `exp(coef)`, y = paste(`Independent variable`, VarGrp, sep = "\n")))+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_point()+
  geom_errorbar(aes(xmin = `lower .95`, xmax = `upper .95`), width = 0)+
  theme_classic()+
  GeneralTheme+
  scale_x_continuous(breaks = c(1, 2, 4, 6))+
  facet_grid(`Independent variable`~., scales = "free_y")


#Figure 5D

int_df <- PhenoDF%>%
  filter(AnyLCHIP_init_filt==FALSE)


surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)

fit.coxph <- coxph(surv_object ~ count_het+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

fit.coxph <- coxph(surv_object ~ mMSS+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

fit.coxph <- coxph(surv_object ~ mMSS+count_het+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)



SF3BTP <- readxl::read_xlsx("SupplementaryDataFile.xlsx", sheet = "Figure 5", range = "A21:I25")

SF3BTP%>%
  mutate(VarGrp = ifelse(grepl("Not adjusted", `Additional adjustment`), "Not adjusted", "Mutually adjusted"))%>%
  ggplot(aes(x = `exp(coef)`, y = paste(`Independent variable`, VarGrp, sep = "\n")))+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_point()+
  geom_errorbar(aes(xmin = `lower .95`, xmax = `upper .95`), width = 0)+
  theme_classic()+
  GeneralTheme+
  scale_x_continuous(breaks = c(1, 2, 4, 6))+
  facet_grid(`Independent variable`~., scales = "free_y")

# Supplementary Table 3 ----------------------------------------

table(PhenoDF$het_yn)
for(j in c("sex", "smk_na", "Anemia", "Thrombocytopenia", "prev_cancer_yn", "AnyLCHIP_init_filt")){
  for(i in c(1, 2)){
    
    print(j)
    print(i)
    print(table(PhenoDF$het_yn, PhenoDF[,(which(colnames(PhenoDF)==j))])[i,])
    print(prop.table(table(PhenoDF$het_yn, PhenoDF[,(which(colnames(PhenoDF)==j))])[i,])*100)
    print(fisher.test(PhenoDF$het_yn, PhenoDF[,(which(colnames(PhenoDF)==j))])[[1]])
    
  }
}

by(PhenoDF$age, PhenoDF$het_yn, quantile)
wilcox.test(PhenoDF$age ~ PhenoDF$het_yn)

# Supplementary Table 5; Supplementary Figure 3 ---------------------------------------------------

int_df <- PhenoDF

ComplexNameCh <- c()
TotalComplexCh <- c()
CLLComplexCh <- c()
MNComplexCh <- c()

for(j in which(grepl("hetcount_complex_", colnames(int_df)))){
  
  ComplexNameCh <- c(ComplexNameCh, (colnames(int_df)[j]))
  TotalComplexCh <- c(TotalComplexCh, sum(((int_df[,j])>0)))
  CLLComplexCh <- c(CLLComplexCh, sum(((int_df[,j])>0)&(int_df$CLL_CENSOR>0)))
  MNComplexCh <- c(MNComplexCh, sum(((int_df[,j])>0)&(int_df$HEME_CENSOR>0)))
  
}

data.frame(ComplexNameCh, TotalComplexCh, CLLComplexCh, MNComplexCh)

surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)
fit.coxph <- coxph(surv_object ~ mMSS_complex_I+mMSS_complex_III+mMSS_complex_IV+mMSS_complex_V+mMSS_complex_DLOOP+mMSS_complex_RRNA+mMSS_complex_TRNA + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)

si_tmp <- summary(fit.coxph)
si_tmp <- si_tmp$coefficients
si_tmp <- as.data.frame(si_tmp)
si_tmp$VOI <- rownames(si_tmp)
si_tmp$GOI <- "CLL"
si_tmp <- cbind(si_tmp, as.data.frame(exp(confint(fit.coxph))))
ROHD <- si_tmp[grepl("mMSS_complex_", si_tmp$VOI),]

surv_object <- Surv(time = int_df$first_incidence_fin, event = int_df$HEME_CENSOR)


fit.coxph <- coxph(surv_object ~ mMSS_complex_I+mMSS_complex_III+mMSS_complex_IV+mMSS_complex_V+mMSS_complex_DLOOP+mMSS_complex_RRNA+mMSS_complex_TRNA + rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)

si_tmp <- summary(fit.coxph)

si_tmp <- si_tmp$coefficients
si_tmp <- as.data.frame(si_tmp)
si_tmp$VOI <- rownames(si_tmp)
si_tmp$GOI <- "MN"
si_tmp <- cbind(si_tmp, as.data.frame(exp(confint(fit.coxph))))

ROHD <- rbind(ROHD, si_tmp[grepl("mMSS_complex_", si_tmp$VOI),])

ROHD$VOI <- factor(ROHD$VOI, levels = rev(c("mMSS_complex_I", "mMSS_complex_III", "mMSS_complex_IV", "mMSS_complex_V",
                                            "mMSS_complex_DLOOP","mMSS_complex_RRNA","mMSS_complex_TRNA")))

ROHD$GOI <- factor(ROHD$GOI, levels = c("MN", "CLL"))

ROHD$SplVar <- (ROHD$VOI %in% c("mMSS_complex_DLOOP","mMSS_complex_RRNA","mMSS_complex_TRNA"))

# Supplementary Figure 3
ROHD%>%
  ggplot(aes(x = log10(`exp(coef)`), y = VOI, color = GOI))+
  geom_vline(xintercept = 0, linetype = 2)+
  geom_point(position = position_dodge(width = 0.75))+
  geom_errorbar(aes(xmin = log10(`2.5 %`), xmax = log10(`97.5 %`)), position = position_dodge(width = 0.75), width = 0)+
  theme_classic()+
  GeneralTheme+
  scale_color_manual(values = c("#4DBBD5", "#E64B35"))+
  coord_cartesian(xlim = c(log10(0.1), log10(1000)))+
  scale_x_continuous(breaks = c(log10(0.1), log10(1), log10(10), log10(100), log10(1000)),
                     labels = c(0.1, 1, 10, 100, 1000))+
  facet_grid(SplVar~., scales = "free_y", space = "free_y")


Genech <- c()
Complexch1 <- c()
Complexch2 <- c()
Zscorech <- c()
pvalch <- c()


for(i in as.character(unique(ROHD$VOI))){
  
  ROHD_tmp <- ROHD[ROHD$VOI %in% i,]  
  
  k <- which(ROHD_tmp$GOI=="CLL")
  l <- which(ROHD_tmp$GOI=="MN")
  
  zscr <- (ROHD_tmp$coef[k]-ROHD_tmp$coef[l])/(sqrt((ROHD_tmp$`se(coef)`[k])^2+(ROHD_tmp$`se(coef)`[l])^2))
  
  Genech <- c(Genech, i)
  Complexch1 <- c(Complexch1, ROHD_tmp$VOI[k])
  Complexch2 <- c(Complexch2, ROHD_tmp$VOI[l])
  Zscorech <- c(Zscorech, zscr)
  if(zscr>=0){
    pvalch <- c(pvalch, 2*pnorm(zscr, lower.tail=FALSE))
  }else{
    pvalch <- c(pvalch, 2*pnorm(zscr, lower.tail=TRUE))
  }
  
}


ROHTOUT <- data.frame(Genech, Complexch1, Complexch2, Zscorech, pvalch)


# Supplementary Table 6 ---------------------------------------------------

int_df <- PhenoDF

surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)

fit.coxph <- coxph(surv_object ~ HC_VAF03+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)
fit.coxph <- coxph(surv_object ~ count_het+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)
fit.coxph <- coxph(surv_object ~ HC_VAF10+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

fit.coxph <- coxph(surv_object ~ mMSS_VAF03+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)
fit.coxph <- coxph(surv_object ~ mMSS+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)
fit.coxph <- coxph(surv_object ~ mMSS_VAF10+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)


# Supplementary Table 7 ---------------------------------------------------

int_df <- PhenoDF

surv_object <- Surv(time = int_df$CLL_follow, event = int_df$CLL_CENSOR)

fit.coxph <- coxph(surv_object ~ count_het+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn, 
                   data = int_df)
summary(fit.coxph)

fit.coxph <- coxph(surv_object ~ count_het+log(CN)+ rms::rcs(age, df = 4) + sex + smk_na + prev_cancer_yn+
                     neutro+lymph+rbc+plt+mono, 
                   data = int_df)
summary(fit.coxph)

# Supplementary Figure 1 ------------------------------------

HCCOMP_DF <- (PhenoDF[,c(1, which(grepl("^hetcount_complex_", colnames(PhenoDF))))])%>%
  melt(., id.vars = "id")

#A
table(HCCOMP_DF$variable, HCCOMP_DF$value>0)%>%
  data.frame()%>%
  filter(Var2 != FALSE)%>%
  mutate(Var1 = str_remove(Var1, "hetcount_complex_"))%>%
  mutate(Var1 = factor(Var1, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                       labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  ggplot(aes(x = Freq, y = Var1))+
  geom_col()+
  theme_classic()+
  GeneralTheme+
  ylab("Complex")+
  xlab("Number of individuals")



KIV <- HCCOMP_DF%>%
  filter(value>0)%>%
  mutate(CombVar = paste0(id, "_", variable))%>%
  pull(CombVar)%>%
  as.character()%>%
  str_replace(., "hetcount", "mMSS")

MSSCOMP_DF <- (PhenoDF[,c(1, which(grepl("^mMSS_complex_", colnames(PhenoDF))))])%>%
  melt(., id.vars = "id")

MSSCOMP_DF <- MSSCOMP_DF[paste0(MSSCOMP_DF$id, "_", MSSCOMP_DF$variable) %in% KIV,]


#B
MSSCOMP_DF%>%
  mutate(variable = str_remove(variable, "mMSS_complex_"))%>%
  mutate(variable = factor(variable, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                           labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  ggplot(aes(y = variable, x = value))+
  geom_violin(scale = "width")+
  geom_boxplot(alpha = 0, width = 0.15)+
  theme_classic()+
  GeneralTheme+
  xlab("mMSS")+
  ylab("Complex")+
  coord_cartesian(xlim = c(0, 2))+
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2))


#C
table(PhenoDF$count_het)%>%
  data.frame()%>%
  mutate(Var1 = as.numeric(as.character(Var1)))%>%
  filter(Var1>0)%>%
  ggplot(aes(x = Freq, y = reorder(Var1, -Var1)))+
  geom_col()+
  theme_classic()+
  GeneralTheme+
  ylab("Number of heteroplasmies")+
  xlab("Number of individuals")



# Supplementary Figure 5 --------------------------------------------------

FreqGenesToAssess <- LCHIP_vars%>%
  distinct(SampID, Gene.refGene)%>%
  pull(Gene.refGene)%>%
  unique()%>%
  as.character()

CHIPfiltch <- c()
Genech <- c()
hvarch <- c()
HetCHch <- c()
TotCHch <- c()
PercHet <- c()
pch <- c()

Medch <- c()


for(i in c("AnyLCHIP_init_filt")){
  
  int_df <- PhenoDF
  
  int_df$HCov0 <- int_df$count_het>0
  
  int_CHIP <- LCHIP_vars
  int_CHIP$Gene <- int_CHIP$Gene.refGene
  
  
  
  
  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID[int_CHIP$Gene == g]
    
    
    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID
    }
    
    
    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]
      
      
      ForLOGDF <- int_df%>%filter(count_het==0 | het_VOI)%>%filter(AnyLCHIP_init_filt == 0| IsGOI)
      
      
      si_ct <- with(ForLOGDF, table(IsGOI, het_VOI))%>%
        data.frame()%>%
        filter(IsGOI == TRUE)%>%
        mutate(s = sum(Freq))%>%
        filter(het_VOI == TRUE)%>%
        mutate(P = Freq/s*100)
      
      
      
      CHIPfiltch <- c(CHIPfiltch, i)
      Genech <- c(Genech, g)
      hvarch <- c(hvarch, hvar)
      
      if(length(si_ct$IsGOI)==1){
        
        HetCHch <- c(HetCHch, si_ct$Freq)
        TotCHch <- c(TotCHch, si_ct$s)
        PercHet <- c(PercHet, si_ct$P)
        
        if(si_ct$Freq == 0){
          Medch <- c(Medch, 0)
        }else{
          Medch <- c(Medch, (ForLOGDF%>%
                               filter(IsGOI & het_VOI)%>%
                               pull(mMSS)%>%
                               median()))
          
        }
        
      }else{
        
        HetCHch <- c(HetCHch, NA)
        TotCHch <- c(TotCHch, NA)
        PercHet <- c(PercHet, NA)
        
      }
      
      
      
      
      rm(Logistic_multi)
      
      
    }
  }
}

DF_General_Gene1 <- table(PhenoDF$het_yn[(PhenoDF$AnyLCHIP_init_filt==FALSE)])%>%
  data.frame()%>%
  mutate(s = sum(Freq))%>%
  mutate(P = Freq/s*100)%>%
  filter(Var1==1)
DF_General_Gene1$Var1 <- "No L-CHIP"
DF_General_Gene1$m <- median(PhenoDF$mMSS[((PhenoDF$AnyLCHIP_init_filt==FALSE)&(PhenoDF$het_yn==1))])


DF_General_Gene2 <- table(PhenoDF$het_yn[(PhenoDF$AnyLCHIP_init_filt==TRUE)])%>%
  data.frame()%>%
  mutate(s = sum(Freq))%>%
  mutate(P = Freq/s*100)%>%
  filter(Var1==1)
DF_General_Gene2$Var1 <- "L-CHIP"
DF_General_Gene2$m <- median(PhenoDF$mMSS[((PhenoDF$AnyLCHIP_init_filt==TRUE)&(PhenoDF$het_yn==1))])


DF_General_Gene3 <- table(PhenoDF$het_yn[((PhenoDF$AnyCHIP_init_filt==FALSE)&(PhenoDF$AnyLCHIP_init_filt==FALSE))])%>%
  data.frame()%>%
  mutate(s = sum(Freq))%>%
  mutate(P = Freq/s*100)%>%
  filter(Var1==1)
DF_General_Gene3$Var1 <- "Neither L-CHIP nor M-CHIP"
DF_General_Gene3$m <- median(PhenoDF$mMSS[((PhenoDF$AnyCHIP_init_filt==FALSE)&(PhenoDF$AnyLCHIP_init_filt==FALSE)&(PhenoDF$het_yn==1))])


PlotDFGene <- LCHIP_vars%>%
  distinct(SampID, Gene.refGene)%>%
  left_join(., (PhenoDF%>%select(id, het_yn, mMSS)), by = c("SampID" = "id"))

HetPrev_Gene_DF <- table(PlotDFGene$Gene.refGene, PlotDFGene$het_yn)%>%
  data.frame()%>%
  group_by(Var1)%>%
  mutate(s = sum(Freq))%>%
  ungroup()%>%
  filter(Var2 == 1)%>%
  mutate(P = Freq/s*100)

MedmMSS_Gene_DF <- PlotDFGene%>%
  filter(het_yn == 1)%>%
  group_by(Gene.refGene)%>%
  summarise(m = median(mMSS))%>%
  ungroup()

table(MedmMSS_Gene_DF$Gene.refGene %in% HetPrev_Gene_DF$Var1)
table(HetPrev_Gene_DF$Var1 %in% MedmMSS_Gene_DF$Gene.refGene)

Tog_Gene_DF <- left_join(HetPrev_Gene_DF, MedmMSS_Gene_DF, by = c("Var1" = "Gene.refGene"))
Tog_Gene_DF$m <- ifelse(is.na(Tog_Gene_DF$m), 0, Tog_Gene_DF$m)
Tog_Gene_DF <- Tog_Gene_DF%>%
  select(-Var2)

Tog_Gene_DF <- rbind(DF_General_Gene1, DF_General_Gene2, DF_General_Gene3, Tog_Gene_DF)

Tog_Gene_DF$SzCat <- cut(Tog_Gene_DF$s, c(0, 10, 100, 1000, 1000000))

sum(is.na(Tog_Gene_DF$SzCat))

UKB_Freq_Tab <- LCHIP_vars%>%
  distinct(SampID, Gene.refGene)%>%
  pull(Gene.refGene)%>%
  table()%>%
  data.frame()%>%
  rename("Var1" = ".")%>%
  arrange(-Freq)


FreqGeneTabvect <- UKB_Freq_Tab%>%
  head(n=30)%>%
  pull(Var1)%>%
  as.character()

options(ggrepel.max.overlaps = Inf)
Tog_Gene_DF$HRG <- (Tog_Gene_DF$Var1 %in% c("MYD88", "IGLL5", "NOTCH1", "XPO1", "LTB", "TCL1A", "TAF1", "DUSP2"))
Tog_Gene_DF$FreqGene <- (Tog_Gene_DF$Var1 %in% FreqGeneTabvect)

with(Tog_Gene_DF%>%filter(s<10000), hist(s))
sort(Tog_Gene_DF$s)

Tog_Gene_DF%>%
  filter(FreqGene|HRG|s>1000)%>%
  ggplot(aes(x = P, y = m))+
  geom_point(aes(size = SzCat, color = HRG))+
  geom_text_repel(aes(label = Var1), size = 3)+
  theme_classic()+
  GeneralTheme+
  scale_x_continuous(breaks = seq(0, 100, by = 10),
                     labels = paste0(seq(0, 100, by = 10), "%"))+
  xlab("Heteroplasmy prevalence")+
  ylab("Median mMSS")+
  scale_color_manual(values = c("#595959", "#E64B35"))





