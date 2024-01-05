# Question: How do Trp catabolites differ in CID patients with high/low serum Trp and control patients?
rm(list = ls())

library(data.table)
library(tableone)
library(lmerTest)
library(ggplot2)
library(data.table)
library(ppcor)

setwd('~/Projects/trpcid_pub/')
source("./scripts/init_dat.R")

trp.derivs.tested = c("IPA", "IAA", "IndSO4", "Kynurenine", "Serotonin", "Kyn.Trp")
mets[, paste0(trp.derivs.tested, ".log") := log10(.SD), .SDcols = c("IPA", "IAA", "IndSO4", "Kynurenine", "Serotonin", "Kyn.Trp")]
mets[, Disease.entity := factor(Disease.entity, c("Control", "CD","IC", "UC","RAneg","Pso","SjogrenSyndrome","Polymyositis","SLE","SSc","RApos","ConnTisUnsp","SpA", "RAunsp"))]

lmm.highlow <- lmm.by.status(dat = mets, metnames = trp.derivs.tested, fixed.effect.tested = "Disease.cohort", fixed.effect.controlled = "Sex")
  # uses log transformed data
high.low.coef.table <- rbindlist(lapply(lmm.highlow$coefs, rbindlist, idcol = "Metabolite"), idcol = "Trp.status")
high.low.coef.table <- high.low.coef.table[grep("^D", rn)] # Only the disease rows, no intercepts
high.low.coef.table[, `:=`(Sig = pasterics(p.adjust(`Pr(>|t|)`, method = "fdr")), FDR = p.adjust(`Pr(>|t|)`, method = "fdr"))]
high.low.coef.table[Sig == "0.2", Sig := "0.20"]
high.low.coef.table[grep("^D", rn),.N]
plot.trp.derivs(plot.mets = trp.derivs.tested, directory = "./out/trp_deriv_boxplots_04012024") 


high.vs.low.lmer <- list()
for (i in trp.derivs.tested) {
  print(i)
  print(summary(lmer(log10(get(i)) ~ Trp.status + Disease.cohort + Sex + (1|Patient.no), mets[Trp.status != "Control"]))$coef)
  high.vs.low.lmer[[i]] <- merge(as.data.table(summary(lmer(log10(get(i)) ~ Trp.status + Disease.cohort + Sex + (1|Patient.no), 
                                                      mets[Trp.status != "Control"]))$coef, keep.rownames = TRUE),
                                 as.data.table(confint(lmer(log10(get(i)) ~ Trp.status + Disease.cohort + Sex + (1|Patient.no), 
                                                            mets[Trp.status != "Control"])), keep.rownames = TRUE), by = "rn")
}

high.vs.low.lmer <- rbindlist(high.vs.low.lmer, idcol = "Metabolite")
high.vs.low.lmer[rn == "Trp.statusLow", FDR := p.adjust(`Pr(>|t|)`, method = "fdr")]

### amino acid differences
amino_acids <- names(mets)[c(16:31, 33:34)]
aa.by.trp.status <- list()
for (aa in amino_acids) {
  print(aa)
  aa.by.trp.status[[aa]] <- merge(as.data.table(summary(lmer(log10(get(aa)) ~ Trp.status + Disease.cohort + Sex + (1|Patient.no), mets))$coef, keep.rownames = TRUE),
                                  as.data.table(confint(lmer(log10(get(aa)) ~ Trp.status + Disease.cohort + Sex + (1|Patient.no), mets)), keep.rownames = TRUE), by = "rn")
}


aa.by.trp.table <- rbindlist(aa.by.trp.status, idcol = "Metabolite")
aa.by.trp.table[Sig == "0.4", Sig := 0.40]
aa.by.trp.table[rn %in% c("Trp.statusLow", "Trp.statusHigh"), FDR := p.adjust(`Pr(>|t|)`, method = "fdr")]
aa.by.trp.table[!is.na(FDR), Sig := pasterics(FDR)]
aa.by.trp.table[Sig == "", Sig := "n.s."]
aa.by.trp.table.low <- aa.by.trp.table[rn == "Trp.statusLow"]
aa.by.trp.table.low[FDR >= 0.05 | `t value` > 0]
aa.by.trp.table.low[FDR < 0.05 & `t value` < 0]
aa.by.trp.table.low[FDR < 0.05 & `t value` > 0]

aa.by.trp.table.high <- aa.by.trp.table[rn == "Trp.statusHigh"]
aa.by.trp.table.high[FDR >= 0.05| `t value` < 0]
aa.by.trp.table.high[FDR < 0.05 & `t value` >  0]

plot_aa <- list()
for (aa in amino_acids) {
  place_sig <- quantile(mets[Trp.status == "High", get(aa)], probs = 1)[1]
  plot_aa <- ggplot(mets, aes(x = Trp.status, y = get(aa), fill = Disease.cohort)) +
    annotate('ribbon', x = c(-Inf, Inf),
             ymin = quantile(mets[Disease.entity == "Control", get(aa)])["25%"],
             ymax = quantile(mets[Disease.entity == "Control", get(aa)])["75%"],
             alpha = 0.5, fill = 'lightgrey') +
    geom_boxplot() + 
    scale_fill_manual(values = c(Control = "lightgrey", GI = "#0a9396",MSK = "#94d2bd", Skin ="#bb3e03")) +
    labs(y = paste(aa, " [μM]"), x = "Trp status", fill = "Affected organ") +
    geom_text(data = data.table(x = 2, y = place_sig), 
              aes(x = x, y = y), 
              label = paste0(aa.by.trp.table[Metabolite == aa & rn == "Trp.statusHigh", Sig]), inherit.aes = FALSE, size = 5) +
    geom_text(data = data.table(x = 3, y = place_sig), 
              aes(x = x, y = y), 
              label =  paste0(aa.by.trp.table[Metabolite == aa & rn == "Trp.statusLow", Sig]), inherit.aes = FALSE, size = 5) +
    theme_bw()
  
  cairo_pdf(paste0("./out/aa_boxplots_21032023/", aa, "_high_vs_low_boxplots_21032023.pdf"), width = 4, height = 3)
  print(plot_aa)
  dev.off()

  print(plot_aa)
}



#### Cohort overview ####
names(mets)
nrow(mets)
nrow(mets[!duplicated(Patient.no)])
mets[!duplicated(Patient.no), .N, Disease.entity]
mets[!duplicated(Patient.no), .N, Disease.cohort]
cohort_overview <- CreateTableOne(vars = c("Trp", "Sex", "Age.visit", "Trp.status"), 
                                  data = mets[!duplicated(Patient.no)],
                                  strata = "Disease.cohort")
cohort_overview_csv <- print(cohort_overview, cramVar = "Sex")
write.csv(cohort_overview_csv, "out/cohort2_cohort_overview_02032023.csv")

cohort_overview_by_observation <- CreateTableOne(vars = c("Trp", "Sex", "Age.visit", "Trp.status"), 
                                                 data = mets,
                                                 strata = "Disease.cohort")
cohort_overview_by_observation_csv <- print(cohort_overview_by_observation, cramVar = "Sex")
write.csv(cohort_overview_by_observation_csv, "out/cohort2__overview_by_observation_19052023.csv")
mets[!duplicated(Patient.no) & Disease.cohort != "Control", .N]
mets[!duplicated(Patient.no), .N, Disease.entity]


trp.bp <- 
  ggplot(mets[!is.na(Trp) & !is.na(Disease.cohort)],
         aes(x = Disease.cohort, y = Trp, fill = Trp.status)) +
  annotate('ribbon', x = c(-Inf, Inf),
           ymin = quantile(mets[Disease.entity == "Control", Trp])["25%"],
           ymax = quantile(mets[Disease.entity == "Control", Trp])["75%"],
           alpha = 0.5, fill = 'lightgrey') +
  geom_boxplot(position = position_dodge(width =1), outlier.shape = NA) +
  geom_point(pch=21, size=1.5, position=position_jitterdodge(jitter.width=0.1, dodge.width=1), alpha = 0.5) +
  scale_fill_manual(values = c("#DBA8AC", "#26547C", "#79A9D1")) +
  labs(x = "Affected organ group", y = "Trp [μM]", fill = "Trp Status")+
  scale_y_log10() +
  theme_bw()
trp.bp

cairo_pdf("./out/trp_deriv_boxplots_29062023/trp_boxplot_diseasecohort_trpstatus_21032023.pdf", 3, 3)
trp.bp + theme(legend.position = "none")
dev.off()



##### Partial correlations ######
#### Feature set-up ####
mets_cor <- mets[!duplicated(Patient.no)]
mets_cor[, Sex := as.numeric(as.factor(Sex))]

trp_high <- mets_cor[Trp.status == "High"]
trp_low <- mets_cor[Trp.status == "Low"]
trp_control <- mets_cor[Trp.status == "Control"]
#trp_control <- trp_control[1:79]
nrow(trp_high)
nrow(trp_low)
nrow(trp_control)
#### Partial correlations using pcorr_mets_cor ####
length(names(mets_cor)[16:34])
names(mets_cor)
which(names(mets_cor) == "Ala") 
which(names(mets_cor) == "Val") 
metabolite <- names(mets_cor)[16:34]
metabolite <- metabolite[metabolite != "Trp"]

pcorr_high  <- pcorr_mets(dat = trp_high, metabolite = metabolite)
pcorr_low  <- pcorr_mets(dat = trp_low, metabolite = metabolite)
pcorr_control <- pcorr_mets(dat = trp_control, metabolite = metabolite) 


cor_coef_table <- as.data.table(cbind(cor.coef.control = sapply(pcorr_control, function(x) (rbind(x$estimate[1, "Trp"]))), 
                                      cor.coef.high = sapply(pcorr_high, function(x) (rbind(x$estimate[1, "Trp"]))), 
                                      cor.coef.low = sapply(pcorr_low, function(x) (rbind(x$estimate[1, "Trp"]))), 
                                      rawp.control = sapply(pcorr_control, function(x) (rbind(x$p.value[1, "Trp"]))), 
                                      rawp.high = sapply(pcorr_high, function(x) (rbind(x$p.value[1, "Trp"]))), 
                                      rawp.low = sapply(pcorr_low, function(x) (rbind(x$p.value[1, "Trp"])))), keep.rownames = "metabolite")

cor_coef_table[, `:=`(fdr.control = p.adjust(rawp.control, method = "fdr"),
                      fdr.high = p.adjust(rawp.high, method = "fdr"),
                      fdr.low = p.adjust(rawp.low, method = "fdr"))]

nrow(cor_coef_table[fdr.control <= 0.05])
nrow(cor_coef_table[fdr.high <= 0.05])
nrow(cor_coef_table[fdr.low <= 0.05])


group.corrs <- cor_coef_table[, 1:4]
setnames(group.corrs, c("metabolite", "cor.coef.control", "cor.coef.high", "cor.coef.low"),
         c("Metabolite", "Control", "High", "Low"))
group.corrs <- as.data.table(melt(group.corrs, variable.name = "Class", value.name = "Corr"))
group.fdr <- cor_coef_table[, c(1,8,9,10)]
setnames(group.fdr, c("metabolite", "fdr.control", "fdr.high", "fdr.low"),
         c("Metabolite", "Control", "High", "Low"))
group.fdr <- as.data.table(melt(group.fdr, variable.name = "Class", value.name = "FDR"))
group.fdr[, negLogFDR := -log10(FDR)]
group.corrs <- merge.data.table(group.corrs, group.fdr, by = c("Metabolite", "Class"))
group.corrs[, Class := factor(Class, c("Low", "High", "Control"))]
group.corrs[FDR < 0.05, negLogFDRsig := negLogFDR]

corr.plot <- ggplot(group.corrs, aes(Metabolite, Class, fill = Corr)) +
  geom_point(aes(size = abs(negLogFDR), color = Corr), shape=19) + 
  geom_point(aes(size = abs(negLogFDRsig)), shape=21) + 
  labs(x = NULL, y = NULL, fill = "Spearman's ρ", title="Serum Amino Acid Correlations with Trp") + 
  theme_classic() +
  scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1), guide = FALSE) +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  theme(text=element_text(family="Roboto")) +
  scale_size(range=c(1,10), guide=NULL)


cairo_pdf("./out/partial_correlation_grouped_plot_02032023.pdf", width = 7, height = 2.5)
corr.plot
dev.off()

### Creatinine associations
pcor(mets_cor[Disease.cohort != "Control", .(Creatinine, `Ind-SO4`, Sex)], method = "spearman")
pcor(mets_cor[Disease.cohort == "Control", .(Creatinine, `Ind-SO4`, Sex)], method = "spearman")
pcor(mets_cor[Disease.cohort != "Control", .(Creatinine, Trp, Sex)], method = "spearman")
pcor(mets_cor[Disease.cohort == "Control", .(Creatinine, Trp, Sex)], method = "spearman")
pcor(mets_cor[Disease.cohort != "Control", .(`Ind-SO4`, Trp, Sex)], method = "spearman")

##### Stool Trp in high/low serum Trp groups ####
t.test(log10(Trp) ~ Trp.status, stool)

stool_trp_violin <- ggplot(stool, aes(x=Trp.status, y=Trp, fill = Trp.status)) + 
  geom_violin(draw_quantiles = 0.5) +
  scale_y_continuous(trans = 'log10') + 
  scale_fill_manual(values = c("#a1cca5ff","#709775ff")) +
  scale_color_manual(values = c("grey90", "grey50", "grey10")) + 
  geom_jitter(height = 0, width = 0.1, alpha = 0.5, aes(col = disease_class)) +
  theme_bw() +
  labs(x = "Serum Trp status", y = "Faecal Trp [µM]", col = "Affected organ", fill = "Trp status") 

pdf("./out/stool_trp_violin_03032023.pdf", 4.5, 3)
stool_trp_violin
dev.off()

