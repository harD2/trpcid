library(data.table)
clindat <- fread("./data/clean/cohort1_public_08082023.csv")
clindat[, `:=`(Hauptdiagnose = factor(Hauptdiagnose, 
                                  levels = c("Control", 
                                             "CD", "UC", 
                                             "AxSpA","PsA", "RAneg","RApos",# joints
                                             "M35.9", "Sicca","SLE","SSc",  # connective tissue
                                             "GCA", "PMR", # blood vessels
                                             "Pso"),
                                  labels = c("Control", 
                                             "CD", "UC", 
                                             "AxSpA","PsA", "RAneg","RApos",# joints
                                             "CTD", "SjS","SLE","SSc",  # connective tissue
                                             "GCA", "PMR", # blood vessels
                                             "Pso")),
               Disease.cohort = factor(Disease.cohort),
               biologicum.status = factor(biologicum.status, c("never", "change"), 
                                          c("No", "Yes")))]
clindat[DAS28 > 9.4, DAS28 := NA] # outside of actual DAS28 range 0 - 9.4
clindat[!is.na(Trp), `:=`(visit.no = 1:.N, no.observations = .N), Patient.no]
mets <- fread("./data/clean/cohort2_serum_biocrates_public_05062023.csv")
setnames(mets, c("3-IAA", "3-IPA", "Ind-SO4"), c("IAA", "IPA", "IndSO4"))
stool <- fread("./data/clean/cohort2_stool_biocrates_public_05062023.csv")

#### Functions ####
pasterics <- function(p) {
  out <- rep("", length(p))
  for(i in 1:length(p)) {
    if(p[i] < 0.1)
      out[i] <- "°"
    if(p[i] < 0.05)
      out[i] <- "*"
    if(p[i] < 0.01)
      out[i] <- "**"
    if(p[i] < 0.001)
      out[i] <- "***"
  }
  return(out)
}


lmm.by.status <- function(dat = mets, metnames = NULL, 
                          fixed.effect.controlled = NULL, # parameters to be included in all LMM (to negate the effect of e.g. sex)
                          fixed.effect.tested = NULL,
                          random.effect = "(1| Patient.no)",
                          directory = "tmp") {
  lmm.low = list()
  lmm.high = list()
  # Split data into Trp high/Low and perform LMM to check for differences between disease cohort and healthy control for metabolites in metnames
  # Output coefficients and p-values to lmm.ligh/low lists
  for (i in metnames) {
    message("Metabolite: ", i)
    H1.high <- lmer(paste0(i, ".log ~ +", fixed.effect.tested, "+", fixed.effect.controlled, "+", random.effect),
                    data = dat[Trp.status %in% c("Control", "High")], REML = FALSE)
    
    H1.low <- lmer(paste0(i, ".log ~ +", fixed.effect.tested, "+", fixed.effect.controlled, "+", random.effect),
                   data = dat[Trp.status %in% c("Control", "Low")], REML = FALSE)
    
    
    lmm.high[[i]] <- as.data.table(summary(H1.high)$coef, keep.rownames = TRUE)
    lmm.low[[i]] <- as.data.table(summary(H1.low)$coef, keep.rownames = TRUE)
    
  }
  return(list(coefs = list(High = lmm.high, Low = lmm.low), mods = list(High.mod = H1.high, Low.mod = H1.low)))
}


plot.trp.derivs <- function (dat = mets, plot.mets = NULL, sig.table = high.low.coef.table, directory = "tmp") {
  
  for (i in plot.mets){
    
    plotms <- ggplot(dat[!is.na(Disease.cohort) & !is.na(Trp.status)], aes(x = Disease.cohort, y = get(i),
                                                                           fill = Trp.status)) +
      annotate('ribbon', x = c(-Inf, Inf),
               ymin = quantile(dat[Disease.entity == "Control", get(i)])["25%"],
               ymax = quantile(dat[Disease.entity == "Control", get(i)])["75%"],
               alpha = 0.5, fill = 'lightgrey') +
      geom_boxplot(position = position_dodge(width =1), outlier.shape = NA) +
      geom_point(pch=21, size=1.5, position=position_jitterdodge(jitter.width=0.1, dodge.width=1), alpha = 0.5) +
      theme_bw() +
      labs(x = "Affected organ group",
           y = if (grepl("\\.", i)) {
             gsub("\\.", ":", i) }
               else {
             paste0(i, " [μM]")
             }, fill = "Trp Status") +
      scale_fill_manual(values = c("#DBA8AC", "#26547C", "#79A9D1")) +
      stat_summary(fun=mean, geom="point", shape=23, size=3, color="black",
                   position = position_dodge2(width = 1,   
                                              preserve = "single")) + 
      geom_text(data = data.frame(x= 1.75, y= quantile(dat[, get(i)], probs = 0.999)),
                aes(x=x, y=y),
                inherit.aes = FALSE, label = sig.table[Trp.status == "High" & Metabolite == i & rn == "Disease.cohortGI", Sig]) +
      geom_text(data = data.frame(x= 2.25, y= quantile(dat[, get(i)], probs = 0.999)),
                aes(x=x, y=y),
                inherit.aes = FALSE, label = sig.table[Trp.status == "Low" & Metabolite == i &rn == "Disease.cohortGI", Sig]) +
      geom_text(data = data.frame(x= 2.75, y= quantile(dat[, get(i)], probs = 0.999)),
                aes(x=x, y=y),
                inherit.aes = FALSE, label = sig.table[Trp.status == "High" & Metabolite == i & rn == "Disease.cohortMSK", Sig]) +
      geom_text(data = data.frame(x= 3.25, y= quantile(dat[, get(i)], probs = 0.999)),
                aes(x=x, y=y),
                inherit.aes = FALSE, label = sig.table[Trp.status == "Low" & Metabolite == i & rn == "Disease.cohortMSK", Sig]) +
      geom_text(data = data.frame(x= 3.75, y= quantile(dat[, get(i)], probs = 0.999)),
                aes(x=x, y=y),
                inherit.aes = FALSE, label = sig.table[Trp.status == "High" & Metabolite == i & rn == "Disease.cohortSkin", Sig]) +
      geom_text(data = data.frame(x= 4.25, y= quantile(dat[, get(i)], probs = 0.999)),
                aes(x=x, y=y),
                inherit.aes = FALSE, label = sig.table[Trp.status == "Low" & Metabolite == i & rn == "Disease.cohortSkin", Sig]) +
      scale_y_log10()
    
    if (dir.exists(directory)) { } # create folder to save the plots to
    else {dir.create(directory)}
    
    cairo_pdf(filename = paste0("./", directory, "/", # save each as .pdf
                                i, "_boxplot_diseasecohort_trpstatus_21032023.pdf"),
              width = 3, height = 3)
    
    print(plotms + theme(legend.position = "none")) # need to call to print or it doesn't work
    dev.off()
    
    print(plotms) # print each plot to screen
    
    
  }
  
}

pcorr_mets <- function(dat = NULL, metabolite = "C0", main_cor_parameter = "Trp", other_cor_parameters = c("Sex")) {
  # Take metabolite and do partial correlation against another metabolite of interest (main_cor_parameter), 
  # controlling for any other parameters (other_cor_parameters). Must be numeric, NAs not tolerated 
  cor_list = list()
  for(i in metabolite){
    cor_list[[i]] <- pcor(dat[, c(..i, ..main_cor_parameter, ..other_cor_parameters)], method = "spearman")
  }
  return(cor_list)
}

mets[, Disease.entity.copy := Disease.entity]
for (i in unique(mets$Disease.entity)) {
  if(mets[Disease.entity == i, .N, Disease.entity][1,2] <= 5) {
    mets[Disease.entity == i, Disease.entity.copy := "Other"]
  }
  
}
