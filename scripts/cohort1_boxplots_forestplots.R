### boxplots forestplots with Trp, CRP
rm(list = ls())
setwd("~/Projects/trpcid_pub/")
source("./scripts/init_dat.R")
library(lmerTest)
library(ggplot2)
library(MASS)
library(car)
library(rstatix)
library(ggpubr)

#### Set up colors ####
disease.colours <- c(Control = "lightgrey",
                     GI = "#0a9396",
                     MSK = "#94d2bd",
                     Skin ="#bb3e03")

hauptdiagnose.colours <- merge(clindat[!duplicated(Hauptdiagnose), .(Hauptdiagnose, Disease.cohort)], 
                               as.data.table(cbind(disease.colours), keep.rownames = "Disease.cohort"))

### Cohort overview ####
nrow(clindat)
clindat[!duplicated(Patient.no), .N, Hauptdiagnose]
clindat[!duplicated(Patient.no) & Hauptdiagnose != "Control", .N]

##### Trp by disease cohort #### 

lmm.all.mod <- lmer("log10(Trp) ~ Hauptdiagnose + sex + (1|Patient.no)", clindat)
plot(fitted(lmm.all.mod), resid(lmm.all.mod))
plot(resid(lmm.all.mod))
lmm.all <- summary(lmm.all.mod)
lmm.all
y_shift <- list() 
for (disease in levels(clindat$Hauptdiagnose)) {
  y_shift[[disease]] <- clindat[Hauptdiagnose %in% disease, quantile(Trp, probs = c(0.35, 0.6), na.rm = TRUE, names = FALSE)]
    # The significance stars are to be plotted at the 0.6 quantile
}

lmm.all.conf <- confint(lmm.all.mod)
hauptdiagnose.coef <- as.data.table(cbind(Hauptdiagnose = c(levels(clindat$Hauptdiagnose), "sex"),
                                          signif(lmm.all$coefficients[,-5], 2),
                                          `CI 2.5%` = signif(lmm.all.conf[-c(1:2), "2.5 %"], 2),
                                          `CI 97.5%` = signif(lmm.all.conf[-c(1:2), "97.5 %"], 2),
                                          lower.tail.p = c(NA,signif(as.numeric(pt(q = lmm.all$coefficients[c(-1,-15),4], lmm.all$coefficients[c(-1,-15),3], lower.tail = TRUE)), 2), NA),
                                          fdr = c(NA, signif(p.adjust(as.numeric(pt(q = lmm.all$coefficients[c(-1,-15),4], lmm.all$coefficients[c(-1,-15),3], lower.tail = TRUE)), method = "fdr"), 2), NA),
                                          # calculate lower tail p-value from t value and df
                                          y.sig = c(sapply(y_shift, "[[", 2), NA), y.no = c(sapply(y_shift, "[[", 1), NA))) # add y shift and sig stars

hauptdiagnose.coef[, `:=`(y.sig = as.numeric(y.sig), y.no = as.numeric(y.no), fdr = as.numeric(fdr))] # y gives plotting location for significance stars 
hauptdiagnose.coef[!is.na(fdr), Sig := pasterics(fdr)]
hauptdiagnose.coef[fdr > 0.05, Sig := paste0("italic('",round(fdr, 2),"')")]
hauptdiagnose.coef <- hauptdiagnose.coef[!Hauptdiagnose %in% "sex"] # do not include sex in output
y_iqr_ctrl <- quantile(clindat[Hauptdiagnose == "Control", Trp], probs = c(.25, .75))   # to graph grey IQR bar
hauptdiagnose.coef <- merge(hauptdiagnose.coef, clindat[!duplicated(Patient.no), .N, Hauptdiagnose], by = "Hauptdiagnose")

write.csv(hauptdiagnose.coef, "./out/ehealth_lmm_by_disease_entity_05062023.csv", row.names = FALSE)
hauptdiagnose.coef[, Disease.cohort := Hauptdiagnose]
disease.entity.box <- ggplot(clindat, aes(Hauptdiagnose, Trp, fill = Disease.cohort)) +
  annotate('ribbon', x = c(-Inf, Inf), ymin = y_iqr_ctrl[1], ymax = y_iqr_ctrl[2],
           alpha = 1, fill = 'lightgrey') +
  geom_boxplot() +
  scale_fill_manual(values = disease.colours) +
  labs(y = "Trp [µM]", x = "") +
  ylim(c(0,100))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        #panel.grid = element_blank(),
        panel.border =  element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  geom_text(data = hauptdiagnose.coef[fdr <= 0.05], aes(Hauptdiagnose, y = y.sig, label = Sig)) +
  geom_text(data = hauptdiagnose.coef[fdr > 0.05], aes(Hauptdiagnose, y = y.sig, label = Sig), size = 3, parse = TRUE)# +
  #geom_text(data = hauptdiagnose.coef, aes(Hauptdiagnose, y = y.no, label = N), size = 2.5)



disease.entity.box
svg("./out/boxplot_disease_entities_15122023.svg", width = 6, height = 3.33)
disease.entity.box +
  # theme(legend.position = c(0.15, 0.92),
  #       legend.background = element_blank(),
  #       legend.text = element_text(size = 8),
  #       legend.title = element_text(size = 8),
  #       legend.key.size = unit(3.5, "mm"))+
  theme(legend.position = "none") + 
  guides(fill = guide_legend(nrow = 2, title = "Affected organ system", size = 12))
dev.off()


### Inactive Trp by disease ####
clindat[CRP < 5, crp.inactive := "Inactive"]
clindat[Hauptdiagnose == "Control", crp.inactive := "Control"]

lmm.crp.inactive.mod <- lmer("log10(Trp) ~ Hauptdiagnose + sex + (1|Patient.no)",clindat[!is.na(crp.inactive)])
plot(fitted(lmm.crp.inactive.mod), resid(lmm.crp.inactive.mod))
plot(resid(lmm.crp.inactive.mod))
plot(qqnorm(resid(lmm.crp.inactive.mod)))
lmm.crp.inactive <- summary(lmm.crp.inactive.mod)
lmm.crp.inactive
lmm.crp.inactive.conf <-confint(lmm.all.mod)
y_shift_crp <- list() 
for (disease in levels(clindat$Hauptdiagnose)) {
  y_shift_crp[[disease]] <- clindat[!is.na(crp.inactive) & Hauptdiagnose %in% disease, quantile(Trp, probs = c(0.35, 0.6), na.rm = TRUE, names = FALSE)]
  # The significance stars are to be plotted at the 0.6 quantile
}


hauptdiagnose.coef.crp <- as.data.table(cbind(Hauptdiagnose = c(levels(clindat$Hauptdiagnose), "sex"),
                                          signif(lmm.crp.inactive$coefficients, 2),
                                          `CI 2.5%` = signif(lmm.crp.inactive.conf[-c(1:2), "2.5 %"], 2),
                                          `CI 97.5%` = signif(lmm.crp.inactive.conf[-c(1:2), "97.5 %"], 2),
                                          fdr = c(NA, signif(p.adjust(lmm.crp.inactive$coefficients[c(-1, -15), 5], method = "fdr"), 2), NA), #NOT LOWER TAILED
                                          y.sig = c(sapply(y_shift_crp, "[[", 2), NA), y.no = c(sapply(y_shift_crp, "[[", 1), NA))) # add y shift and sig stars

hauptdiagnose.coef.crp[, `:=`(y.sig = as.numeric(y.sig), y.no = as.numeric(y.no), fdr = as.numeric(fdr))] # y gives plotting location for significance stars
hauptdiagnose.coef.crp[fdr <= 0.05, Sig := pasterics(as.numeric(fdr))]
hauptdiagnose.coef.crp[fdr > 0.05, Sig := paste0("italic('", round(fdr, 2), "')")]

hauptdiagnose.coef.crp[Hauptdiagnose == "CD", pasterics(as.numeric(fdr))]
hauptdiagnose.coef.crp <- hauptdiagnose.coef.crp[!Hauptdiagnose %in% "sex"] # do not include sex in output
hauptdiagnose.coef.crp <- merge(hauptdiagnose.coef.crp, clindat[!duplicated(Patient.no) & !is.na(crp.inactive), .N, Hauptdiagnose], by = "Hauptdiagnose")
hauptdiagnose.coef.crp[, round(.SD, 4), .SDcols = c("Estimate"#, "Std. Error", "df", "t value", "Pr(>|t|)"
                                                    )]

write.csv(hauptdiagnose.coef.crp, "./out/ehealth_lmm_by_disease_entity_inactive_crp_28062023.csv", row.names = FALSE)
hauptdiagnose.coef.crp[, Disease.cohort := Hauptdiagnose]
disease.entity.inactive.crp.box <- ggplot(clindat[!is.na(crp.inactive)], aes(Hauptdiagnose, Trp, fill = Disease.cohort)) +
  annotate('ribbon', x = c(-Inf, Inf), ymin = y_iqr_ctrl[1], ymax = y_iqr_ctrl[2],
           alpha = 1, fill = 'lightgrey') +
  geom_boxplot(#outlier.shape = NA
               ) +
  scale_fill_manual(values = disease.colours) +
  labs(y = "Trp [µM]", x = "") +
  ylim(c(0,100))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        #panel.grid = element_blank(),
        panel.border =  element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  geom_text(data = hauptdiagnose.coef.crp[fdr <= 0.05], aes(Hauptdiagnose, y = y.sig, label = Sig)) +
  geom_text(data = hauptdiagnose.coef.crp[fdr > 0.05], aes(Hauptdiagnose, y = y.sig, label = Sig), size = 3, parse = TRUE)# +
  #geom_text(data = hauptdiagnose.coef.crp, aes(Hauptdiagnose, y = y.no, label = N), size = 2)

disease.entity.inactive.crp.box
hauptdiagnose.coef.crp
svg("./out/boxplot_disease_entities_inactive_crp_28062023.svg", width = 6, height = 3.33)
disease.entity.inactive.crp.box +
  # theme(legend.position = c(0.15, 0.92),
  #       legend.background = element_blank(),
  #       legend.text = element_text(size = 8),
  #       legend.title = element_text(size = 8),
  #       legend.key.size = unit(3.5, "mm"))+
  guides(fill = guide_legend(nrow = 2, title = "Affected organ system", size = 8))
dev.off()

##### Trp and CRP #####
diagnoses <- clindat[Hauptdiagnose != "Control", unique(Hauptdiagnose)]
trp_lmm_by_disease <- list()
for (x in diagnoses) {
    trp_lmm_by_disease[[x]]<- list(sum.mod = summary(lmer("log10(CRP) ~ log10(Trp) + sex + (1|Patient.no)",clindat[Hauptdiagnose == x])),
                                   conf.mod = confint(lmer("log10(CRP) ~ log10(Trp) + sex + (1|Patient.no)",clindat[Hauptdiagnose == x])))
    
    plot(resid(trp_lmm_by_disease[[x]]$sum.mod), main = x)
} 


coef_list <- as.data.table(cbind(t(sapply(trp_lmm_by_disease, function(x) x$sum.mod$coefficients[2,])), t(sapply(trp_lmm_by_disease, function(x) x$conf.mod["log10(Trp)",]))), keep.rownames = "Hauptdiagnose")

coef_list[, `:=`(#min.est = Estimate-`Std. Error`, max.est = Estimate+`Std. Error`,
  fdr = p.adjust(`Pr(>|t|)`, method = "fdr"))]
coef_list[fdr <= 0.05, p.asterix := pasterics(fdr)]
coef_list[fdr > 0.05, p.round := pasterics(fdr)]
coef_list <- merge(coef_list, hauptdiagnose.colours, by = "Hauptdiagnose", all.x = T)


write.csv(coef_list, "./out/ehealth_lmm_trp_and_crp_14032023.csv")

coef_list[, Hauptdiagnose := factor(Hauptdiagnose, levels = rev(levels(clindat$Hauptdiagnose)))]
crp_forest <- 
  ggplot(coef_list, mapping = aes(x = Hauptdiagnose, y = Estimate, ymin =`2.5 %`, ymax=`97.5 %`)) + 
  geom_errorbar(aes(color = Disease.cohort), size = 1, width = 0)+
  geom_point(aes(color = Disease.cohort), shape = 18, size = 3) +
  geom_hline(yintercept=0, lty=2, color = "darkred", linewidth = 0.7, alpha = 0.7) +  # add a dotted line at x=1 after flip
  coord_flip()+   # flip coordinates (puts labels on y axis)
  xlab(element_blank()) + 
  ylab("Estimate (± 95% CI)") +
  scale_color_manual(values = disease.colours,
                     name = "Affected organ\nsystem") +
    
  geom_text(label = coef_list$p.asterix, vjust = +0.1, color = "#005f73") +
  geom_text(label = coef_list$p.round, vjust = -0.7, color = "#005f73", size = 3) +
  theme_classic() + 
  theme(axis.text = element_text(color = "black"))

  

cairo_pdf("./out/forest_plot_trp_crp_14032023.pdf", width = 4, height = 3.6)
crp_forest +
  theme(legend.position = "none")
dev.off()

### By disease activity score ####
clindat[!is.na(BASDAI), range(BASDAI)]
clindat[!is.na(DAS28), range(DAS28)]
clindat[!is.na(PASI), range(PASI)]


lmm.basdai <- lmer(BASDAI ~ log10(Trp) + sex+ 
                     (1|Patient.no), clindat)
lmm.das28.neg <- lmer(log(DAS28) ~ log10(Trp) + sex + 
                        (1|Patient.no), clindat[Hauptdiagnose == "RAneg"])
lmm.das28.pos <- lmer(log(DAS28) ~ log10(Trp) + sex + 
                        (1|Patient.no), clindat[Hauptdiagnose == "RApos"])
lmm.cdai <- lmer(cdai ~ log10(Trp) + sex +
                   (1|Patient.no), clindat)
lmm.cmayo <- lmer(complete.mayo ~ log10(Trp) + sex + 
                    (1|Patient.no), clindat)
lm.pasi <- glm(PASI ~ log10(Trp), clindat, family = "gaussian") # only one female, one observation per patient#


df.residual(lm.pasi) #degrees of freedom
plot(fitted(lmm.basdai), resid(lmm.basdai))
plot(fitted(lmm.das28.neg), resid(lmm.das28.neg))
plot(fitted(lmm.das28.pos), resid(lmm.das28.pos))
plot(residuals(lm.pasi), fitted(lm.pasi))
summary(lmm.das28.neg)
summary(lmm.das28.pos)
summary(lmm.basdai)
summary(lm.pasi)
summary(lmm.cdai)
summary(lmm.cmayo)
range(clindat$DAS28, na.rm = T)

count_dai_patients <- clindat[!is.na(BASDAI) | !is.na(PASI) | !is.na(DAS28)]
count_dai_patients[!duplicated(Patient.no) & !is.na(BASDAI), .N, Hauptdiagnose]
count_dai_patients[!duplicated(Patient.no) & !is.na(DAS28), .N, .(Hauptdiagnose, sex)]
count_dai_patients[!duplicated(Patient.no) & !is.na(PASI), .N, Hauptdiagnose]

summary(lmm.basdai)
clindat[BASDAI <= 2.1 | Hauptdiagnose == "Control", Trp]


clindat[, dai.inactive := fifelse(BASDAI < 2.1 | DAS28 < 2.6 | cdai < 150 | complete.mayo <= 4, "Inactive", NA_character_, NA_character_)]
clindat[!is.na(dai.inactive), dai.inactive := paste0(dai.inactive, ".", Hauptdiagnose)]
clindat[Hauptdiagnose == "Control", dai.inactive := "Control"]
clindat[, .N, dai.inactive]


basdai.point <- ggplot(clindat[!is.na(BASDAI)], aes(y = BASDAI, x = Trp)) + 
  geom_point(color = "#94d2bd") + 
  geom_smooth(method = "lm", color = "black") +
  scale_x_continuous(limits = c(0,100)) + 
  scale_y_continuous(limits = c(0,8)) + 
  labs(x = "Trp [µM]") +
  annotate("text", y = 8*0.95, x = 50, label = paste0("p = 0.010; LMM"), size = 4) +
  theme_bw()
basdai.point


das28.point.neg <- ggplot(clindat[Hauptdiagnose == "RAneg"], aes(y = DAS28, x = Trp)) + 
  geom_point(color = "#94d2bd") + 
  geom_smooth(method = "lm", color = "black") +
  scale_x_continuous(limits = c(0,100)) + 
  scale_y_continuous(limits = c(0,8)) +
  annotate("text", y = 8*0.95, x = 50, label = paste0("p = 0.0029; LMM"), size = 4) +
  labs(x = "Trp [µM]", y = "DAS28  - RAneg") +
  theme_bw()


das28.point.pos <- ggplot(clindat[Hauptdiagnose == "RApos"], aes(y = DAS28, x = Trp)) + 
  geom_point(color = "#94d2bd") + 
  geom_smooth(method = "lm", color = "black") +
  scale_x_continuous(limits = c(0,100)) + 
  scale_y_continuous(limits = c(0,8)) + 
  annotate("text", y = 8*0.95, x = 50, label = paste0("p = 0.58; LMM"), size = 4) +
  labs(x = "Trp [µM]", y = "DAS28  - RApos") +
  theme_bw()
das28.point.pos


pasi.point <- ggplot(clindat, aes(x = Trp, y = PASI)) + 
  geom_point(color = "#bb3e03") + 
  geom_smooth(method = "lm", color = "black") +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0,45)) +
  annotate("text", y = 45*0.95, x = 50, label = paste0("p = 0.95; linear regression"), size = 4) +
  labs(x = "Trp [µM]") +
  theme_bw()
pasi.point

cdai.point <- ggplot(clindat, aes(x = Trp, y = cdai)) +
  geom_point(color = "#0a9396") +
  geom_smooth(method = "lm", color = "black") +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0,650)) +
  annotate("text", y = 650*0.95, x = 50, label = paste0("p = 0.00021; LMM"), size = 4) +
  labs(x = "Trp [µM]", y = "CDAI") +
  theme_bw()
cdai.point
#log10(Trp)   -143.69      38.42  377.00  -3.740 0.000213 ***

cmayo.point <- ggplot(clindat, aes(x = Trp, y = complete.mayo)) +
  geom_point(color = "#0a9396") +
  geom_smooth(method = "lm", color = "black") +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0,15)) +
  annotate("text", y = 15*0.95, x = 50, label = paste0("p < 0.0001; LMM"), size = 4) +
  labs(x = "Trp [µM]", y = "Total Mayo") +
  theme_bw()
cmayo.point
#log10(Trp)   -9.3816     1.7001 164.3449  -5.518 1.31e-07 ***


dai_combi <- ggarrange(cdai.point, cmayo.point, basdai.point, das28.point.neg, das28.point.pos, pasi.point, ncol = 3, nrow = 2)
dai_combi

cairo_pdf("./out/dai_point_plots_04042023.pdf", 12,8)
dai_combi
dev.off()



dai_coef <- rbind(PASI = c(summary(lm.pasi)$coef[2,1:2], df = df.residual(lm.pasi), summary(lm.pasi)$coef[2,3:4], confint(lm.pasi)["log10(Trp)",]),
                  DAS28_RAneg = c(summary(lmm.das28.neg)$coef[2,], confint(lmm.das28.neg)["log10(Trp)",]),
                  DAS28_RApos = c(summary(lmm.das28.pos)$coef[2,], confint(lmm.das28.pos)["log10(Trp)",]),
                  BASDAI = c(summary(lmm.basdai)$coef[2,], confint(lmm.basdai)["log10(Trp)",]),
                  CDAI = c(summary(lmm.cdai)$coef[2,], confint(lmm.cdai)["log10(Trp)",]),
                  cMayo = c(summary(lmm.cmayo)$coef[2,], confint(lmm.cmayo)["log10(Trp)",]))
dai_coef <- as.data.table(dai_coef, keep.rownames = "DAI")
write.csv(dai_coef, "./out/dai_lmm_lm_estimate_tables.csv", row.names = FALSE)


#### Variance by patient ####
clindat[no.observations >= 3, `:=`(trp_sd = sd(Trp, na.rm = TRUE), crp_sd = sd(CRP, na.rm = TRUE)), Patient.no] 
clindat[!duplicated(Patient.no) & no.observations >= 3 & Hauptdiagnose != "Control", .N]
clindat[!duplicated(Patient.no) & Hauptdiagnose != "Control", .(mean.obs = mean(no.observations), sd = sd(no.observations))]

plot(density(clindat[!duplicated(Patient.no) & !is.na(crp_sd), crp_sd]))
plot(density(clindat[!duplicated(Patient.no) & !is.na(trp_sd), trp_sd]))
wcxt_trp <- wilcox.test(trp_sd ~ biologicum.status, data = clindat[!duplicated(Patient.no)])
wilcox.test(trp_sd ~ biologicum.status, data = clindat[!duplicated(Patient.no) & sex == "Female"])
wilcox.test(trp_sd ~ biologicum.status, data = clindat[!duplicated(Patient.no) & sex == "Male"])

wcxt_crp <- wilcox.test(crp_sd ~ biologicum.status, data = clindat[!duplicated(Patient.no)])
wilcox.test(crp_sd ~ biologicum.status, data = clindat[!duplicated(Patient.no) & sex == "Female"])
wilcox.test(crp_sd ~ biologicum.status, data = clindat[!duplicated(Patient.no) & sex == "Male"])


varplot_trp <- 
  ggplot(clindat[!duplicated(Patient.no) & !is.na(biologicum.status)], aes(x = biologicum.status, y = trp_sd)) +
  geom_boxplot(aes(fill = biologicum.status)) + 
  labs(y = "Trp standard deviation", x = "Therapy escalation") + 
  scale_fill_manual(values = c("#ee9b00ff", "#8a6623ff")) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        panel.border =  element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"))+
  scale_y_continuous(trans = "log10", labels = ~(format(.x, scientific = FALSE)), limits = c(0.01,300)) +
  geom_line(data = data.table(x = c(1, 2), y = c(70, 70)), 
            aes(x = x, y = y), inherit.aes = FALSE) +
  geom_text(data = data.table(x = c(1.5), y = c(100)), 
            aes(x = x, y = y), 
            label = paste0("p = ", signif(wcxt_trp$p.value, 2)), inherit.aes = FALSE, size = 4)

pdf("out/ehealth_trp_variance_biologicum_barplots_16032023.pdf", width = 3, height = 3)
varplot_trp
dev.off()



varplot_crp <- 
  ggplot(clindat[!duplicated(Patient.no) &  !is.na(biologicum.status) & !is.na(crp_sd)], aes(x = biologicum.status, y = crp_sd)) +
  geom_boxplot(aes(fill = biologicum.status)) +
  labs(y = "CRP Standard Deviation", x = "Therapy escalation") + 
  scale_fill_manual(values = c("#ee9b00ff", "#8a6623ff")) +
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        panel.border =  element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  scale_y_continuous(trans = "log10", labels = ~(format(.x, scientific = FALSE)), limits = c(0.01, 300)) +
  geom_line(data = data.table(x = c(1, 2), y = c(175, 175)),
            aes(x = x, y = y), inherit.aes = FALSE) +
  geom_text(data = data.table(x = c(1.5), y = c(250)),
            aes(x = x, y = y),
            label = paste0("p = ", signif(wcxt_crp$p.value, 2)), inherit.aes = FALSE, size = 4)

pdf("out/ehealth_crp_variance_biologicum_barplots_16032023.pdf", width = 3, height = 3)
varplot_crp
p == dev.off()



individualplot <- ggplot(clindat[!is.na(Trp) & Patient.no %in% c(383, 30) & visit.no <= 20,  ],
                              aes(x = visit.no, y = Trp, group = Patient.no, color = biologicum.status)) +
  facet_wrap(~biologicum.status, nrow = 1) +
  geom_line()+
  geom_point() +
  theme_bw() +
  ylim(c(0,60)) +
  ylab("Serum Tryptophan [µM]") +
  xlab("Visit number") +
  labs(title = "Therapy escalation") + 
  theme(strip.background =element_rect(fill="#f0f0f0"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(colour = 'black'), legend.position = "none") +
  scale_color_manual(values = c("#ee9b00ff", "#8a6623ff"))



cairo_pdf("./out/individual_trp_course_examples_wide_20102022.pdf", height = 2, width = 3)
individualplot
dev.off()

variance_combi <- ggarrange(individualplot, varplot_trp, varplot_crp, ncol = 3)
pdf("./out/cohort1_trp_crp_variance_therapy_escalation_10082023.pdf", width = 9, height = 3)
variance_combi
dev.off()



#### Cohort overview ####
library(tableone)
overview_table <- CreateTableOne(vars = "Trp", data = clindat, strata = "Hauptdiagnose")
overview_table_by_patient <- CreateTableOne(vars = c("sex", "age", "biologicum.status"), data = clindat[!duplicated(Patient.no)], strata = "Hauptdiagnose")

overview_table_csv <- print(overview_table)
overview_table_by_patient_csv<- print(overview_table_by_patient, cramVars = c("sex", "biologicum.status"))
write.csv(overview_table_csv, file = "./out/ehealth_cohort_overview_all_obs_11082023.csv")
write.csv(overview_table_by_patient_csv, file = "./out/ehealth_cohort_overview_by_patient_11082023.csv")

tmp <- clindat[!duplicated(Patient.no) & Hauptdiagnose != "Control", .N, Hauptdiagnose]
sum(tmp$N)
clindat[, .N] #total observations
mean(clindat[Hauptdiagnose != "Control", .N, Patient.no]$N) #mean number of visits
sd(clindat[Hauptdiagnose != "Control", .N, Patient.no]$N)
clindat[!is.na(cdai), .N]
clindat[!is.na(complete.mayo), .N]

### Low CRP and low Trp with disease activity
clindat[, low.trp := fifelse(Trp < 45, "low Trp", "normal or high", NA_character_)]
clindat[, low.crp := fifelse(CRP < 5, "low CRP", "high CRP", NA_character_)]

cor.test(clindat[low.crp == "low CRP", Trp], clindat[low.crp == "low CRP", complete.mayo], use = "complete.obs", method = "pearson")
cor.test(clindat[low.crp == "low CRP", Trp], clindat[low.crp == "low CRP", cdai], use = "complete.obs", method = "pearson")
cor.test(clindat[low.crp == "low CRP",CRP], clindat[low.crp == "low CRP", complete.mayo], use = "complete.obs", method = "pearson")
cor.test(clindat[low.crp == "low CRP", CRP], clindat[low.crp == "low CRP", cdai], use = "complete.obs", method = "pearson")
cor.test(clindat[low.crp == "low CRP", CRP], clindat[low.crp == "low CRP", BASDAI], use = "complete.obs", method = "pearson")
cor.test(clindat[low.crp == "low CRP", CRP], clindat[low.crp == "low CRP", Trp], use = "complete.obs", method = "pearson")


png("./tmp/lowcrp_lowtrp_completemayo.png")
ggplot(clindat[!is.na(low.trp) & low.crp == "low CRP"], aes(x = low.trp, y = complete.mayo)) +
  geom_boxplot()
dev.off()

png("./tmp/lowcrp_lowtrp_cdai.png")
ggplot(clindat[!is.na(low.trp)& low.crp == "low CRP"], aes(x = low.trp, y = cdai)) +
  geom_boxplot()
dev.off()

png("./tmp/lowcrp_lowtrp_basdai.png")
ggplot(clindat[!is.na(low.trp)& low.crp == "low CRP"], aes(x = low.trp, y = BASDAI)) +
  geom_boxplot()
dev.off()
png("./tmp/lowcrp_lowtrp_das28.png")
ggplot(clindat[!is.na(low.trp)& low.crp == "low CRP"], aes(x = low.trp, y = DAS28)) +
  geom_boxplot()
dev.off()



#### Association with DAI for patients with low CRP ####
summary(lmer(cdai ~ low.trp + (1|Patient.no), clindat[low.crp == "low CRP"]))
summary(lmer(Trp ~ CRP + (1|Patient.no), clindat[low.crp == "low CRP"]))
summary(lmer(Trp ~ CRP + (1|Patient.no), clindat[low.crp == "low CRP" & Hauptdiagnose %in% c("UC")]))
summary(lmer(Trp ~ CRP + (1|Patient.no), clindat[low.crp == "low CRP" & Hauptdiagnose %in% c("RAneg")]))
summary(lmer(Trp ~ CRP + (1|Patient.no), clindat[low.crp == "high CRP"]))

summary(lmer(complete.mayo ~ low.trp + (1|Patient.no), clindat[low.crp == "low CRP"]))
summary(glm(BASDAI ~ low.trp, clindat[low.crp == "low CRP"], family = "gaussian"))
summary(glm(DAS28 ~ low.trp, clindat[low.crp == "low CRP"], family = "gaussian"))
clindat[!is.na(DAS28), Hauptdiagnose]
clindat[low.crp == "low CRP" & !is.na(cdai), .N, low.trp]
clindat[low.crp == "low CRP" & !is.na(complete.mayo), .N, low.trp]

