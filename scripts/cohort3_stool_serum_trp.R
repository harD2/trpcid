setwd("~/Projects/trpcid_pub/")
library(data.table)
library(lmerTest)
stool <- fread("./data/clean/cohort3_export_pub_stool_04012024.csv")
## Calculate difference in dates between stool and serum
stool <- stool[!is.na(Trp.clin)]
stool[, .(mean_diff = mean(diff_days_stool_serum, na.rm = TRUE),
          sd_diff = sd(diff_days_stool_serum, na.rm = TRUE),
          range_diff = range(diff_days_stool_serum, na.rm = TRUE))]

## Create LMM for the relationship between stool and serum Trp
stool.lmm <- lmer(log10(Trp) ~ log10(Trp.clin) + Sex_factor + Disease_cohort + (1|Patient_no), stool)
summary(lmer(log10(Trp) ~ log10(Trp.clin) + Sex_factor + (1|Patient_no), stool[Disease_cohort == "GI"]))
confint(lmer(log10(Trp) ~ log10(Trp.clin) + Sex_factor + (1|Patient_no), stool[Disease_cohort == "GI"]))

summary(lmer(log10(Trp) ~ log10(Trp.clin) + Sex_factor + (1|Patient_no), stool[Disease_cohort == "MSK"]))
confint(lmer(log10(Trp) ~ log10(Trp.clin) + Sex_factor + (1|Patient_no), stool[Disease_cohort == "MSK"]))

stool.conf <- confint(stool.lmm)
stool.lmm.overview <- cbind(summary(stool.lmm)$coef, stool.conf[3:6,])

stool$Trp.clin
## Visualize relationship 
disease.colours <- c(Control = "lightgrey",
                     GI = "#0a9396",
                     MSK = "#94d2bd",
                     Skin ="#bb3e03")

ggplot(stool, aes(y = Trp, x = Trp.clin, color = Disease_cohort)) + 
  geom_point() +
  scale_color_manual(values = c("#0a9396", "#94d2bd")) + 
  scale_x_continuous(limits = c(0,100)) +
  labs(y = "Stool Trp [µM]", x = "Serum Trp [µM]", color = "Affected organ system") +
  #annotate("text", y = 500*0.75, x = 50, label = paste0("B. est: -0.18, 95% CI: -1.21 - 0.82\n p = 0.728"), size = 4) +
  theme_bw()

cairo_pdf("./out/cohort3_stool_serum_trp_04012024.pdf", 6,4)
ggplot(stool, aes(y = Trp, x = Trp.clin, color = Disease_cohort)) + 
  geom_point() +
  scale_color_manual(values = c("#0a9396", "#94d2bd")) + 
  scale_x_continuous(limits = c(0,100)) +
  labs(y = "Stool Trp [µM]", x = "Serum Trp [µM]", color = "Affected organ system") +
  #annotate("text", y = 500*0.75, x = 50, label = paste0("B. est: -0.18, 95% CI: -1.21 - 0.82\n p = 0.728"), size = 4) +
  theme_bw()
dev.off()


unique(stool$Treatment_1)
stool$ICD_10
### Cohort overview ###
#### Cohort overview ####
library(tableone)

stool[!duplicated(Patient_no), .N, ICD_10]

overview_table_by_patient <- CreateTableOne(vars = c("Sex_factor", "Age", "Trp.clin", "Treatment_1"), strata = "Disease_cohort", data = stool[!duplicated(Patient_no)])
overview_table_by_patient_csv <- print(overview_table_by_patient, cramVars = "Sex_factor")

overview_table_by_week <- CreateTableOne(vars = "Time_weeks", strata = "Disease_cohort", factorVars = "Time_weeks", stool)
overview_table_by_week_csv <- print(overview_table_by_week)

write.csv(overview_table_by_patient_csv, file = "./out/cohort3_stool_sample_overview_by_patient_04012024.csv")
write.csv(overview_table_by_week_csv, file = "./out/cohort3_stool_sample_overview_by_week_04012024.csv")

