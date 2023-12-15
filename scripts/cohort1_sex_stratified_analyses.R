## Sex specific analysis of cohort 1
crp.mod <- summary(lmer(log10(CRP) ~ sex + (1|Patient.no), clindat))
plot(resid(crp.mod))
ggplot(clindat[!is.na(sex)], aes(x = sex, y = Trp)) +
  geom_boxplot()
plot(density(log10(clindat$CRP), na.rm = T))

summary(lmer(log10(Trp) ~ sex + (1|Patient.no), clindat))
ggplot(clindat[!is.na(sex)], aes(x = sex, y = Trp)) +
  geom_boxplot()
