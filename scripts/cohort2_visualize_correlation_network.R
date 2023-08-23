### Different correlation network visualization approaches

library(data.table)
library(psych)
library(ggplot2)
library(ggraph)
library(patchwork)
setwd('~/Projects/trpcid_pub/')
source("./scripts/init_dat.R")

mets$Sex <- as.numeric(as.factor(mets$Sex))
mets <- mets[!duplicated(Patient.no)]

mets[, .N, Trp.status] # Lowest number of observations within a group is 82

pcor_FDR_filter <- function(dat, trp.status) {
  pcor_full <- list()
  pcor_pruned <- list()
  pcor_pval <- list()
  pcor_pruned_dt <- list()
  
  for (i in trp.status) {
    message("Partial Spearman correlation controlling for sex in Trp status: ", i)  
    # Partial Spearman correlation
    tmp <- dat[Trp.status == i]
    # tmp <- tmp[1:82] # only plot the first 82 observations per group to match n of Trp.status Low
    pcor_full[[i]] <- partial.r(data = tmp, x = 16:34, y = "Sex", method = "spearman")
    
    # Save results for filtering by FDR, take only top half of matrix
    nrow_col = nrow(pcor_full[[i]])
    pcor_pruned[[i]] <- pcor_full[[i]]
    
    # Get p-values and adjust for FDR - only upper half has adjusted value
    pcor_pval[[i]] <- corr.p(pcor_full[[i]], n= dat[Trp.status == i, .N], adjust = "fdr")
    # All correlations >= FDR 0.05 are set to 999
    pcor_pruned[[i]][which(pcor_pval[[i]]$p >= 0.05)] = NA
    pcor_pruned[[i]][lower.tri(pcor_full[[i]])] = 999
    pcor_pruned_dt[[i]] <- as.data.table(pcor_pruned[[i]], keep.rownames = "node1")
    pcor_pruned_dt[[i]] <- melt(pcor_pruned_dt[[i]], variable.name = "node2",
                                value.name = "correlation")
    
    # Remove the second instance of the correlation and the diagonal
    pcor_pruned_dt[[i]] <- pcor_pruned_dt[[i]][!correlation %in% 999 & node1 != node2]
    
  }
  
  return(list(pcor_full = pcor_full, pcor_pval = pcor_pval, graph_input = pcor_pruned_dt))
  
}

prepped_partial_cors <- pcor_FDR_filter(dat = mets, trp.status = c("Control", "High", "Low"))


# Create the correlation network plots
pControl <- ggraph(graph =prepped_partial_cors$graph_input$Control, layout = "stress") +
  geom_edge_link(aes(edge_width = abs(correlation), color = correlation)) +
  geom_node_label(aes(label = name), family = "Arial") +
  scale_edge_width(range = c(0.5, 3), name = "Correlation (abs)", limits = c(0,1)) +
  scale_edge_color_gradient2(low = "darkred", mid = "white", high = "blue", midpoint = 0, 
                             name = "Correlation", limits = c(-0.5,1))+
  theme_graph(base_family = "Arial") +
  labs(title = "Control") +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 10)))

pHigh <- ggraph(graph =prepped_partial_cors$graph_input$High, layout = "stress") +
  geom_edge_link(aes(edge_width = abs(correlation), color = correlation)) +
  geom_node_label(aes(label = name), family = "Arial") +
  scale_edge_width(range = c(0.5, 3), name = "Correlation (abs)", limits = c(0,1)) +
  scale_edge_color_gradient2(low = "darkred", mid = "white", high = "blue", midpoint = 0, 
                             name = "Correlation", limits = c(-0.5,1))+
  theme_graph(base_family = "Arial") +
  labs(title = "High serum Trp") +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 10)))

pLow <- ggraph(graph =prepped_partial_cors$graph_input$Low, layout = "stress") +
  geom_edge_link(aes(edge_width = abs(correlation), color = correlation)) +
  geom_node_label(aes(label = name), family = "Arial") +
  scale_edge_width(range = c(0.5, 3), name = "Correlation (abs)", limits = c(0,1)) +
  scale_edge_color_gradient2(low = "darkred", mid = "white", high = "blue", midpoint = 0, 
                             name = "Correlation", limits = c(-0.5,1))+
  theme_graph(base_family = "Arial") +
  labs(title = "Low serum Trp") +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 10)))

prepped_partial_cors$graph_input$Control[node1 == "Trp"]

pcor_networks <- pControl + pHigh +pLow + plot_layout(ncol = 3) + plot_layout(guides = "collect")

cairo_pdf("./out/trpcid_pcor_networks_22062023.pdf", 16, 5)
pcor_networks
dev.off()
