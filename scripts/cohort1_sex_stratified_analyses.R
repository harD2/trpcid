### Sex-based differences in Trp metabolism  ###
rm(list = ls())
setwd("~/Projects/trpcid_pub/")
source("./scripts/init_dat.R")
library(ggpubr)
library(lmerTest)
clindat[, CID.status := fifelse(Hauptdiagnose == "Control", "Control", "CID")]
male.colors <- c(Control = "#C9B5D5", CID = "#9874AF")
female.colors <- c(Control = "#9BC6CA", CID = "#4daeb7")

### 1. Directly plot sex differences in Trp
summary(lmer(log10(Trp) ~ sex + (1|Patient.no), clindat))
trp.box <- ggplot(clindat[!is.na(sex)], aes(x = sex, y = Trp, fill= sex)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        #legend.position = "none",
        panel.border =  element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  labs(x = element_blank(), y = "Trp [µM]", fill = "Sex") +
  scale_fill_manual(values = c("#4daeb7", "#9874AF")) +
  geom_line(data = data.table(x = c(1.1, 1.9), y = c(155, 155)), 
            aes(x = x, y = y), inherit.aes = FALSE) +
  geom_text(data = data.table(x = c(1.5), y = c(161)), 
            aes(x = x, y = y), 
            label = "p < 0.0001", inherit.aes = FALSE, size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

### 2. CRP plotted
# summary(lmer(log10(CRP) ~ sex + (1|Patient.no), clindat))
# confint(lmer(log10(CRP) ~ sex + (1|Patient.no), clindat))
# sexMale     -0.003997701 0.08826662
crp.box <- ggplot(clindat[!is.na(sex)], aes(x = sex, y = CRP, fill= sex)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        panel.border =  element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  labs(x = element_blank(), y = "CRP [mg/L]") +
  scale_fill_manual(values = c("#4daeb7", "#9874AF")) +
  geom_line(data = data.table(x = c(1.1, 1.9), y = c(450, 450)), 
            aes(x = x, y = y), inherit.aes = FALSE) +
  geom_text(data = data.table(x = c(1.5), y = c(700)), 
            aes(x = x, y = y), 
            label = "p = 0.074", inherit.aes = FALSE, size = 4) + 
  scale_y_continuous(trans = "log10", labels = ~(format(.x, scientific = FALSE)), limits = c(0.01,1000)) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

### 3a. Trp by disease: female
lmm.females <- lmer(log10(Trp) ~ Hauptdiagnose + (1|Patient.no), clindat[sex == "Female"])

plot(fitted(lmm.females), resid(lmm.females))
plot(resid(lmm.females))
lmm.females.sum <- summary(lmm.females)
lmm.females.sum
y_shift_female <- list() 
for (disease in levels(clindat$Hauptdiagnose)) {
  y_shift_female[[disease]] <- clindat[Hauptdiagnose %in% disease & sex == "Female", quantile(Trp, probs = c(0.35, 0.6), na.rm = TRUE, names = FALSE)]
  # The significance stars are to be plotted at the 0.6 quantile
}
signif(lmm.females.sum$coefficients[,-5], 2)
lmm.females.conf <- confint(lmm.females)
hauptdiagnose.coef.female <- as.data.table(cbind(Hauptdiagnose = c(levels(clindat$Hauptdiagnose)),
                                          signif(lmm.females.sum$coefficients[,-5], 2),
                                          `CI 2.5%` = signif(lmm.females.conf[-c(1:2), "2.5 %"], 2),
                                          `CI 97.5%` = signif(lmm.females.conf[-c(1:2), "97.5 %"], 2),
                                          lower.tail.p = c(NA,as.numeric(pt(q = lmm.females.sum$coefficients[c(-1,-15),4], lmm.females.sum$coefficients[c(-1,-15),3], lower.tail = TRUE))),
                                          fdr = c(NA, p.adjust(as.numeric(pt(q = lmm.females.sum$coefficients[c(-1,-15),4], lmm.females.sum$coefficients[c(-1,-15),3], lower.tail = TRUE)), method = "fdr")),
                                          # calculate lower tail p-value from t value and df
                                          y.sig = c(sapply(y_shift_female, "[[", 2)), y.no = c(sapply(y_shift_female, "[[", 1)))) # add y shift and sig stars

hauptdiagnose.coef.female[, `:=`(y.sig = as.numeric(y.sig), y.no = as.numeric(y.no), fdr = as.numeric(fdr))] # y gives plotting location for significance stars 
hauptdiagnose.coef.female[!is.na(fdr), Sig := pasterics(fdr)]
hauptdiagnose.coef.female[fdr > 0.05, Sig := paste0("italic('",signif(fdr, 2),"')")]
hauptdiagnose.coef.female[, CID.status := fifelse(Hauptdiagnose == "Control", "Control", "CID")]

y_iqr_ctrl <- quantile(clindat[Hauptdiagnose == "Control" & sex == "Female", Trp], probs = c(.25, .75))   # to graph grey IQR bar
hauptdiagnose.coef.female <- merge(hauptdiagnose.coef.female, clindat[!duplicated(Patient.no), .N, Hauptdiagnose], by = "Hauptdiagnose")

write.csv(hauptdiagnose.coef.female, "./out/ehealth_lmm_by_disease_entity_females_03012024.csv", row.names = FALSE)
hauptdiagnose.coef.female[, Disease.cohort := Hauptdiagnose]
disease.entity.female <- ggplot(clindat[sex == "Female"], aes(Hauptdiagnose, Trp, fill = CID.status)) +
  annotate('ribbon', x = c(-Inf, Inf), ymin = y_iqr_ctrl[1], ymax = y_iqr_ctrl[2],
           alpha = 1, fill = "#9BC6CA") +
  geom_boxplot() +
  scale_fill_manual(values = female.colors) +
  labs(y = "Trp [µM]", x = element_blank()) +
  ylim(c(0,100))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        #panel.grid = element_blank(),
        panel.border =  element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  geom_text(data = hauptdiagnose.coef.female[fdr <= 0.05], aes(Hauptdiagnose, y = y.sig, label = Sig)) +
  geom_text(data = hauptdiagnose.coef.female[fdr > 0.05], aes(Hauptdiagnose, y = y.sig, label = Sig), size = 3, parse = TRUE)# +
#geom_text(data = hauptdiagnose.coef.female, aes(Hauptdiagnose, y = y.no, label = N), size = 2.5)
disease.entity.female



### 3b. Trp by disease: male
lmm.males <- lmer(log10(Trp) ~ Hauptdiagnose + (1|Patient.no), clindat[sex == "Male"])

plot(fitted(lmm.males), resid(lmm.males))
plot(resid(lmm.males))
lmm.males.sum <- summary(lmm.males)
lmm.males.sum
y_shift_male <- list() 
for (disease in levels(clindat$Hauptdiagnose)) {
  y_shift_male[[disease]] <- clindat[Hauptdiagnose %in% disease & sex == "Male", quantile(Trp, probs = c(0.35, 0.6), na.rm = TRUE, names = FALSE)]
  # The significance stars are to be plotted at the 0.6 quantile
}
signif(lmm.males.sum$coefficients[,-5], 2)
lmm.males.conf <- confint(lmm.males)
hauptdiagnose.coef.male <- as.data.table(cbind(Hauptdiagnose = c(levels(clindat$Hauptdiagnose)),
                                          signif(lmm.males.sum$coefficients[,-5], 2),
                                          `CI 2.5%` = signif(lmm.males.conf[-c(1:2), "2.5 %"], 2),
                                          `CI 97.5%` = signif(lmm.males.conf[-c(1:2), "97.5 %"], 2),
                                          lower.tail.p = c(NA,as.numeric(pt(q = lmm.males.sum$coefficients[c(-1,-15),4], lmm.males.sum$coefficients[c(-1,-15),3], lower.tail = TRUE))),
                                          fdr = c(NA, p.adjust(as.numeric(pt(q = lmm.males.sum$coefficients[c(-1,-15),4], lmm.males.sum$coefficients[c(-1,-15),3], lower.tail = TRUE)), method = "fdr")),
                                          # calculate lower tail p-value from t value and df
                                          y.sig = c(sapply(y_shift_male, "[[", 2)), y.no = c(sapply(y_shift_male, "[[", 1)))) # add y shift and sig stars

hauptdiagnose.coef.male[, `:=`(y.sig = as.numeric(y.sig), y.no = as.numeric(y.no), fdr = as.numeric(fdr))] # y gives plotting location for significance stars 
hauptdiagnose.coef.male[!is.na(fdr), Sig := pasterics(fdr)]
hauptdiagnose.coef.male[fdr > 0.05, Sig := paste0("italic('",signif(fdr, 2),"')")]
hauptdiagnose.coef.male[, CID.status := fifelse(Hauptdiagnose == "Control", "Control", "CID")]

y_iqr_ctrl <- quantile(clindat[Hauptdiagnose == "Control" & sex == "Male", Trp], probs = c(.25, .75))   # to graph grey IQR bar
hauptdiagnose.coef.male <- merge(hauptdiagnose.coef.male, clindat[!duplicated(Patient.no), .N, Hauptdiagnose], by = "Hauptdiagnose")

write.csv(hauptdiagnose.coef.male, "./out/ehealth_lmm_by_disease_entity_males_03012024.csv", row.names = FALSE)



disease.entity.male <- 
  ggplot(clindat[sex == "Male"], aes(Hauptdiagnose, Trp, fill = CID.status)) +
  annotate('ribbon', x = c(-Inf, Inf), ymin = y_iqr_ctrl[1], ymax = y_iqr_ctrl[2],
           alpha = 1, fill = "#C9B5D5") +
  geom_boxplot() +
  scale_fill_manual(values = male.colors) + 
  labs(y = "Trp [µM]", x = element_blank()) +
  ylim(c(0,100))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        #panel.grid = element_blank(),
        panel.border =  element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  geom_text(data = hauptdiagnose.coef.male[fdr <= 0.05], aes(Hauptdiagnose, y = y.sig, label = Sig)) +
  geom_text(data = hauptdiagnose.coef.male[fdr > 0.05], aes(Hauptdiagnose, y = y.sig, label = Sig), size = 3, parse = TRUE)# +

disease.entity.male

### 4a. Trp and CRP female
diagnoses <- clindat[Hauptdiagnose != "Control", unique(Hauptdiagnose)]
trp_lmm_by_disease_female <- list()
for (x in diagnoses) {
  trp_lmm_by_disease_female[[x]]<- list(sum.mod = summary(lmer("log10(CRP) ~ log10(Trp) +  (1|Patient.no)",clindat[Hauptdiagnose == x & sex == "Female"])),
                                 conf.mod = confint(lmer("log10(CRP) ~ log10(Trp) + (1|Patient.no)",clindat[Hauptdiagnose == x  & sex == "Female"])))
  
  plot(resid(trp_lmm_by_disease_female[[x]]$sum.mod), main = x)
} 


coef_list_female <- as.data.table(cbind(t(sapply(trp_lmm_by_disease_female, function(x) x$sum.mod$coefficients["log10(Trp)",])), t(sapply(trp_lmm_by_disease_female, function(x) x$conf.mod["log10(Trp)",]))), keep.rownames = "Hauptdiagnose")

coef_list_female[, `:=`(#min.est = Estimate-`Std. Error`, max.est = Estimate+`Std. Error`,
  fdr = p.adjust(`Pr(>|t|)`, method = "fdr"))]
coef_list_female[fdr <= 0.05, p.asterix := pasterics(fdr)]
coef_list_female[fdr > 0.05, p.round := pasterics(fdr)]
coef_list_female[, CID.status := fifelse(Hauptdiagnose == "Control", "Control", "CID")]

write.csv(coef_list_female, "./out/ehealth_lmm_trp_and_crp_female_03012024.csv")

coef_list_female[, Hauptdiagnose := factor(Hauptdiagnose, levels = rev(levels(clindat$Hauptdiagnose)))]
crp_forest_female <- 
  ggplot(coef_list_female, mapping = aes(x = Hauptdiagnose, y = Estimate, ymin =`2.5 %`, ymax=`97.5 %`)) + 
  geom_hline(yintercept=0, lty=2, color = "#471A19", linewidth = 0.7, alpha = 0.7) + 
  geom_errorbar(aes(color = CID.status), size = 1, width = 0)+
  geom_point(aes(color = CID.status), shape = 18, size = 3) +
  coord_flip()+   # flip coordinates (puts labels on y axis)
  xlab(element_blank()) + 
  ylab("Estimate (± 95% CI)") +
  scale_color_manual(values = female.colors) +
  geom_text(label = coef_list_female$p.asterix, vjust = +0.1, color = "#04404C") +
  geom_text(label = coef_list_female$p.round, vjust = -0.7, color = "#04404C", size = 3) +
  theme_classic() + 
  theme(axis.text = element_text(color = "black"))

cairo_pdf("./out/forest_plot_trp_crp_14032023.pdf", width = 4, height = 3.6)
crp_forest_female +
  theme(legend.position = "none")
dev.off()

### 4b. Trp and CRP male
trp_lmm_by_disease_male <- list()
for (x in diagnoses) {
  trp_lmm_by_disease_male[[x]]<- list(sum.mod = summary(lmer("log10(CRP) ~ log10(Trp) +  (1|Patient.no)",clindat[Hauptdiagnose == x & sex == "Male"])),
                                        conf.mod = confint(lmer("log10(CRP) ~ log10(Trp) + (1|Patient.no)",clindat[Hauptdiagnose == x  & sex == "Male"])))
  
  plot(resid(trp_lmm_by_disease_male[[x]]$sum.mod), main = x)
} 


coef_list_male <- as.data.table(cbind(t(sapply(trp_lmm_by_disease_male, function(x) x$sum.mod$coefficients["log10(Trp)",])), t(sapply(trp_lmm_by_disease_male, function(x) x$conf.mod["log10(Trp)",]))), keep.rownames = "Hauptdiagnose")

coef_list_male[, `:=`(#min.est = Estimate-`Std. Error`, max.est = Estimate+`Std. Error`,
  fdr = p.adjust(`Pr(>|t|)`, method = "fdr"))]
coef_list_male[fdr <= 0.05, p.asterix := pasterics(fdr)]
coef_list_male[fdr > 0.05, p.round := pasterics(fdr)]
coef_list_male[, CID.status := fifelse(Hauptdiagnose == "Control", "Control", "CID")]


write.csv(coef_list_male, "./out/ehealth_lmm_trp_and_crp_male_03012024.csv")

coef_list_male[, Hauptdiagnose := factor(Hauptdiagnose, levels = rev(levels(clindat$Hauptdiagnose)))]
crp_forest_male <- 
  ggplot(coef_list_male, mapping = aes(x = Hauptdiagnose, y = Estimate, ymin =`2.5 %`, ymax=`97.5 %`)) + 
  geom_hline(yintercept=0, lty=2, color = "#471A19", linewidth = 0.7, alpha = 0.7) +  # add a dotted line at x=1 after flip 
  geom_errorbar(aes(color = CID.status), size = 1, width = 0)+
  geom_point(aes(color = CID.status), shape = 18, size = 3) +
  coord_flip()+   # flip coordinates (puts labels on y axis)
  xlab(element_blank()) + 
  ylab("Estimate (± 95% CI)") +
  scale_color_manual(values = male.colors,
                     name = "Affected organ\nsystem") +
  geom_text(label = coef_list_male$p.asterix, vjust = +0.1, color = "#04404C") +
  geom_text(label = coef_list_male$p.round, vjust = -0.7, color = "#04404C", size = 3) +
  theme_classic() + 
  theme(axis.text = element_text(color = "black"))

cairo_pdf("./out/forest_plot_trp_crp_14032023.pdf", width = 4, height = 3.6)
crp_forest_male +
  theme(legend.position = "none")
dev.off()

### 5. Trp and disease entities
# DAI: female
lmm.basdai.female <- lmer(BASDAI ~ log10(Trp) +
                     (1|Patient.no), clindat[sex == "Female"])
lmm.das28.neg.female <- lmer(log(DAS28) ~ log10(Trp) + 
                        (1|Patient.no), clindat[Hauptdiagnose == "RAneg" & sex == "Female"])
lmm.das28.pos.female <- lmer(log(DAS28) ~ log10(Trp) + 
                        (1|Patient.no), clindat[Hauptdiagnose == "RApos" & sex == "Female"])
lmm.cdai.female <- lmer(cdai ~ log10(Trp) + 
                   (1|Patient.no), clindat[sex == "Female"])
lmm.cmayo.female <- lmer(complete.mayo ~ log10(Trp) + 
                    (1|Patient.no), clindat[sex == "Female"])

# DAI: male
lmm.basdai.male <- lmer(BASDAI ~ log10(Trp) +
                            (1|Patient.no), clindat[sex == "Male"])
lmm.das28.neg.male <- lmer(log(DAS28) ~ log10(Trp) + 
                               (1|Patient.no), clindat[Hauptdiagnose == "RAneg" & sex == "Male"])
lmm.das28.pos.male <- lmer(log(DAS28) ~ log10(Trp) + 
                               (1|Patient.no), clindat[Hauptdiagnose == "RApos" & sex == "Male"])
lmm.cdai.male <- lmer(cdai ~ log10(Trp) + 
                          (1|Patient.no), clindat[sex == "Male"])
lmm.cmayo.male <- lmer(complete.mayo ~ log10(Trp) + 
                           (1|Patient.no), clindat[sex == "Male"])

paste0("p = ", signif(summary(lmm.basdai.female)$coef["log10(Trp)", "Pr(>|t|)"], 2), " (females)\n p = ", signif(summary(lmm.basdai.male)$coef["log10(Trp)", "Pr(>|t|)"], 2), " (males)")

### Graphs
basdai.point <- 
  ggplot(clindat[!is.na(BASDAI)], aes(y = BASDAI, x = Trp, color = sex, fill = sex)) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "grey22") +
  scale_x_continuous(limits = c(0,100)) + 
  scale_y_continuous(limits = c(0,8)) + 
  labs(x = "Trp [µM]", fill = "Sex", color = "Sex") +
  annotate("text", y = 8*0.95, x = 50, label = paste0("p = ", signif(summary(lmm.basdai.female)$coef["log10(Trp)", "Pr(>|t|)"], 2), 
                                                      " (Female)\np = ", signif(summary(lmm.basdai.male)$coef["log10(Trp)", "Pr(>|t|)"], 2), " (Male)"), size = 4) +
  scale_fill_manual(values =  c("#30c5d2", "#6C18A0")) +
  scale_color_manual(values =  c("#30c5d2", "#6C18A0")) +
  theme_bw()
basdai.point



das28.point.neg <- ggplot(clindat[Hauptdiagnose == "RAneg"], aes(y = DAS28, x = Trp, color = sex, fill = sex)) + 
  geom_point() +
  scale_fill_manual(values =  c("#30c5d2", "#6C18A0")) +
  scale_color_manual(values =  c("#30c5d2", "#6C18A0")) + 
  geom_smooth(method = "lm", color = "grey22") +
  scale_x_continuous(limits = c(0,100)) + 
  scale_y_continuous(limits = c(0,8)) +
  annotate("text", y = 8*0.95, x = 50, label = paste0("p = ", signif(summary(lmm.das28.neg.female)$coef["log10(Trp)", "Pr(>|t|)"], 2), 
                                                      " (Female)\np = ", signif(summary(lmm.das28.neg.male)$coef["log10(Trp)", "Pr(>|t|)"], 2), " (Male)"), size = 4) +
  labs(x = "Trp [µM]", y = "DAS28  - RAneg", fill = "Sex", color = "Sex") +
  theme_bw()
das28.point.neg


das28.point.pos <- ggplot(clindat[Hauptdiagnose == "RApos"], aes(y = DAS28, x = Trp, color = sex, fill = sex)) + 
  geom_point()+
  scale_fill_manual(values =  c("#30c5d2", "#6C18A0")) +
  scale_color_manual(values =  c("#30c5d2", "#6C18A0")) + 
  geom_smooth(method = "lm", color = "grey22") +
  scale_x_continuous(limits = c(0,100)) + 
  scale_y_continuous(limits = c(0,8)) + 
  annotate("text", y = 8*0.95, x = 50, label = paste0("p = ", signif(summary(lmm.das28.pos.female)$coef["log10(Trp)", "Pr(>|t|)"], 2), 
                                                      " (Female)\np = ", signif(summary(lmm.das28.pos.male)$coef["log10(Trp)", "Pr(>|t|)"], 2), " (Male)"), size = 4) +
  labs(x = "Trp [µM]", y = "DAS28  - RApos", fill = "Sex", color = "Sex") +
  theme_bw()
das28.point.pos


cdai.point <- ggplot(clindat[!is.na(sex)], aes(x = Trp, y = cdai, color = sex, fill = sex)) +
  geom_point() +
  scale_fill_manual(values =  c("#30c5d2", "#6C18A0")) +
  scale_color_manual(values =  c("#30c5d2", "#6C18A0")) + 
  geom_smooth(method = "lm", color = "grey22") +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0,650)) +
  annotate("text", y = 650*0.95, x = 50, label = paste0("p = ", signif(summary(lmm.cdai.female)$coef["log10(Trp)", "Pr(>|t|)"], 2), 
                                                      " (Female)\np = 0.0003 (Male)"), size = 4) +
  labs(x = "Trp [µM]", y = "CDAI", fill = "Sex", color = "Sex") +
  theme_bw()
cdai.point

cmayo.point <- ggplot(clindat[!is.na(sex)], aes(x = Trp, y = complete.mayo, color = sex, fill = sex)) +
  geom_point() +
  scale_fill_manual(values =  c("#30c5d2", "#6C18A0")) +
  scale_color_manual(values =  c("#30c5d2", "#6C18A0")) + 
  geom_smooth(method = "lm", color = "grey22") +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0,15)) +
  annotate("text", y = 15*0.95, x = 50, label = paste0("p < 0.0001 (Female)\np = ", signif(summary(lmm.cmayo.male)$coef["log10(Trp)", "Pr(>|t|)"], 2), " (Male)"), size = 4) +
  labs(x = "Trp [µM]", y = "Total Mayo", fill = "Sex", color = "Sex") +
  theme_bw()
cmayo.point


pdf("./out/sex_stratified_analysis_boxplots_03012024.pdf", width = 10, height = 12)
ggarrange(ggarrange(trp.box, disease.entity.female, crp.box, disease.entity.male, widths = c(0.5, 1), common.legend = TRUE),
          ggarrange(crp_forest_female + theme(legend.position = "none"), crp_forest_male + theme(legend.position = "none")), ncol = 1, heights = c(2, 1))
dev.off()

pdf("./out/sex_stratified_analysis_points_03012024.pdf", width = 10, height = 8)
ggarrange(basdai.point, das28.point.pos, das28.point.neg, cdai.point, cmayo.point, common.legend = TRUE, legend = "top")
dev.off()
