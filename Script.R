##### Required Packages #####
library(dplyr)
library(plyr)
library(stats)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(data.table)
library(pheatmap)
library(igraph)
library(stringdist)
library(readxl)
library(ggraph)
library(tidygraph)
library(circlize)



##### GLIPH differential motif analsyis #####
# Reading GLIPH output
gliph_output <- read_xlsx("GLIPH_Output.csv") %>% filter(Fisher_score<0.05)


# Calculating basic statistics
## freq_summed - sum of the frequency column (normalized by the number of productive sequences)
## contribution - contribution of each TCR in a GLIPH group (pattern)
## TCR contribution by group (e.g. Naive subsets cotribution in )
## fold change based on TCR contribution from each group (e.g. TCM vs Naive)
## summed group contribution for poisson statistics (round off to the nearest integer)
## perform poisson statistics
gliph_output <- gliph_output %>% group_by(pattern) %>% mutate(freq_summed=sum(Freq), contribution=(Freq/freq_summed)*100, 
                                                              TCM=sum(.data$contribution[grep("_TCM",.data$Sample)])/length(grep("_TCM",.data$Sample)),
                                                              TEM=sum(.data$contribution[grep("_TEM",.data$Sample)])/length(grep("_TEM",.data$Sample)),
                                                              Naive=sum(.data$contribution[grep("_Naive",.data$Sample)])/length(grep("_Naive",.data$Sample)),
                                                              TCM_by_Naive=TCM/Naive,
                                                              TEM_by_Naive=TEM/Naive,
                                                              log2FC_TCM_by_Naive=log2(TCM_by_Naive),
                                                              log2FC_TEM_by_Naive=log2(TEM_by_Naive),
                                                              summed_TCM_contribution=round(sum(.data$contribution[grep("_TCM",.data$Sample)])),
                                                              summed_TEM_contribution=round(sum(.data$contribution[grep("_TEM",.data$Sample)])),
                                                              summed_Naive_contribution=round(sum(.data$contribution[grep("_Naive",.data$Sample)])))

write.csv(gliph_output, "Stats_GLIPH2_Naive_TEM_TCM.csv")



# Poisson Statistics
TCM_vs_Naive <- gliph_output %>% filter(!(TCM_by_Naive)=="NaN")
gliph_TCM_vs_Naive_stats <- NULL

for(i in unique(TCM_vs_Naive$pattern)){
  print(paste0("Calculating statistics for: ",i))
  df <- gliph_output %>% filter(pattern==i)
  gliph_TCM_vs_Naive_stats[[i]] <- poisson.test(c(unique(df$summed_TCM_contribution), unique(df$summed_Naive_contribution)), c(length(grep("_TCM",df$Sample)),length(grep("_Naive",df$Sample))))
  
}

TCM_vs_Naive_stats <- rbindlist(gliph_TCM_vs_Naive_stats, idcol = T) %>% distinct(.id, p.value,estimate, parameter)


TCM_vs_Naive_stats_summarized <- merge(TCM_vs_Naive, TCM_vs_Naive_stats, by.x="pattern", by.y=1) %>% distinct(pattern,Fisher_score,number_unique_cdr3,final_score, hla_score, vb_score, expansion_score, length_score,
                                                                                                                    type,freq_summed,TCM, Naive, log2FC_TCM_by_Naive, summed_TCM_contribution, p.value)
write.csv(TCM_vs_Naive_stats_summarized, "GLIPH2_TCM_vs_Naive_stats_summarized.csv")



TEM_vs_Naive <- gliph_output %>% filter(!(TEM_by_Naive)=="NaN")
gliph_TEM_vs_Naive_stats <- NULL

for(i in unique(TEM_vs_Naive$pattern)){
  print(paste0("Calculating statistics for: ",i))
  df <- gliph_output %>% filter(pattern==i)
  gliph_TEM_vs_Naive_stats[[i]] <- poisson.test(c(unique(df$summed_TEM_contribution), unique(df$summed_Naive_contribution)), c(length(grep("_TEM",df$Sample)),length(grep("_Naive",df$Sample))))
  
}

TEM_vs_Naive_stats <- rbindlist(gliph_TEM_vs_Naive_stats, idcol = T) %>% distinct(.id, p.value,estimate, parameter)

TEM_vs_Naive_stats_summarized <- merge(TEM_vs_Naive, TEM_vs_Naive_stats, by.x="pattern", by.y=1) %>% distinct(pattern,Fisher_score,number_unique_cdr3,final_score, hla_score, vb_score, expansion_score, length_score,
                                                                                                              type, freq_summed, TEM, Naive, log2FC_TEM_by_Naive, summed_TEM_contribution, p.value)

write.csv(TEM_vs_Naive_stats_summarized, "GLIPH2_TEM_vs_Naive_stats_summarized.csv")



# Volcano Plot TCM vs Naive
TCM_vs_Naive_stats_summarized <- read.csv("GLIPH2_TCM_vs_Naive_stats_summarized.csv", row.names = 1)
TEM_vs_Naive_stats_summarized <- read.csv("GLIPH2_TEM_vs_Naive_stats_summarized.csv", row.names = 1)


mycolors <- c("blue", "red", "black")
names(mycolors) <- c("Up in Naive", "Up in TCM", "NA")

TCM_vs_Naive_stats_summarized$diffexpressed <- "NA"

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
TCM_vs_Naive_stats_summarized$diffexpressed[TCM_vs_Naive_stats_summarized$log2FC_TCM_by_Naive > 0 & TCM_vs_Naive_stats_summarized$p.value < 0.05] <- "Up in TCM"

# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
TCM_vs_Naive_stats_summarized$diffexpressed[TCM_vs_Naive_stats_summarized$log2FC_TCM_by_Naive < 0 & TCM_vs_Naive_stats_summarized$p.value < 0.05] <- "Up in Naive"


pointSize <- 14
lineWidth <- 1 / 2.835

png("TCM_vs_Naive_VolcanoPlot.png", width =10, height = 10, res = 600, units = "in",bg = "white")
TCM_vs_Naive_stats_summarized %>% filter(expansion_score<0.05) %>% ggplot(aes(x=log2FC_TCM_by_Naive, y=-log10(p.value), label=pattern, col=diffexpressed)) + geom_point(aes(size = as.numeric(expansion_score))) + theme_minimal()  + geom_text_repel(aes(size = expansion_score)) +
  geom_vline(xintercept=c(-2,2), col="red") +
  geom_hline(yintercept=-log10(0.00001), col="red") + scale_color_manual(values = mycolors)+
  scale_size(trans = 'reverse')+
  labs(x = "log2 fold change", y = "-log10(p-value)") + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+ scale_x_continuous(breaks = scales::pretty_breaks(n =8), limits = )+
  guides(colour = guide_legend(override.aes = list(size=4))) + scale_size(name="expansion_score", breaks = c(0.001,0.01,0.05),
                                                                          range=c(2.5, 1), limits = c(0.001,0.05), labels = c("p \u2264 0.001", "p \u2264 0.01", "p \u2264 0.05")) +
  theme(
    text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(size = lineWidth, colour = "black"),
    plot.title  = element_text(size = pointSize, colour = "black",face = "bold.italic"),
    axis.title  = element_text(size = pointSize, colour = "black", face = "bold.italic"),
    axis.text.x  = element_text(size = pointSize * 0.8, colour = "black",vjust = 0.5, hjust=1, face ="bold.italic"),
    axis.text.y  = element_text(size = pointSize * 0.8, colour = "black", hjust = 1, face = "bold.italic"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = pointSize * 0.8, colour = "black", face = "bold.italic"),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.3, "cm"),
    axis.line = element_line(size = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
  ) 
dev.off()



# Volcano Plot TEM vs Naive
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("Up in Naive", "Up in TEM", "NA")

TEM_vs_Naive_stats_summarized$diffexpressed <- "NA"

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
TEM_vs_Naive_stats_summarized$diffexpressed[TEM_vs_Naive_stats_summarized$log2FC_TEM_by_Naive > 0 & TEM_vs_Naive_stats_summarized$p.value < 0.05] <- "Up in TEM"

# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
TEM_vs_Naive_stats_summarized$diffexpressed[TEM_vs_Naive_stats_summarized$log2FC_TEM_by_Naive < 0 & TEM_vs_Naive_stats_summarized$p.value < 0.05] <- "Up in Naive"


pointSize <- 14
lineWidth <- 1 / 2.835

png("TEM_vs_Naive_VolcanoPlot.png", width =10, height = 10, res = 600, units = "in",bg = "white")
TEM_vs_Naive_stats_summarized %>% filter(expansion_score<0.05) %>% ggplot(aes(x=log2FC_TEM_by_Naive, y=-log10(p.value), label=pattern, col=diffexpressed))  + geom_point(aes(size = as.numeric(expansion_score))) + theme_minimal()  + geom_text_repel(aes(size = expansion_score)) +
  geom_vline(xintercept=c(-2,2), col="red") +
  geom_hline(yintercept=-log10(0.00001), col="red") + scale_color_manual(values = mycolors)+
  scale_size(trans = 'reverse')+
  labs(x = "log2 fold change", y = "-log10(p-value)") + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+ scale_x_continuous(breaks = scales::pretty_breaks(n =8), limits = )+
  guides(colour = guide_legend(override.aes = list(size=4))) + scale_size(name="expansion_score", breaks = c(0.001,0.01,0.05),
                                                                          range=c(2.5, 1), limits = c(0.001,0.05), labels = c("p \u2264 0.001", "p \u2264 0.01", "p \u2264 0.05")) +
  theme(
    text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(size = lineWidth, colour = "black"),
    plot.title  = element_text(size = pointSize, colour = "black",face = "bold.italic"),
    axis.title  = element_text(size = pointSize, colour = "black", face = "bold.italic"),
    axis.text.x  = element_text(size = pointSize * 0.8, colour = "black",vjust = 0.5, hjust=1, face ="bold.italic"),
    axis.text.y  = element_text(size = pointSize * 0.8, colour = "black", hjust = 1, face = "bold.italic"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = pointSize * 0.8, colour = "black", face = "bold.italic"),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.3, "cm"),
    axis.line = element_line(size = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
  ) 
dev.off()



#### TCR Sample Proportion ####
# Reading GLIPH output
sig_patterns <- read_xlsx("Significant_Motifs.xlsx") %>% dplyr::select(pattern)


png("Sample_Proportion_byDonor.png", width = 7, height = 4, units = "in", bg = "white", res = 400)
ggplot(gliph_subset, aes(Donor, fill = Type)) +
  geom_bar(position = "fill") + ylab("Percentage %")+
  scale_y_continuous(labels = scales::percent)
dev.off()



##### V and V/J analysis by donor #####
# Frequency Tables
sheets <- excel_sheets(path = "Enriched_Motifs_byDonors.xlsx")

for(i in sheets){
  v_vj_usage <- read_xlsx("Enriched_Motifs_byDonors.xlsx", sheet =i) 
  v_vj_usage$Type <- str_split(v_vj_usage$Sample, ":", simplify = T)[,2] %>% str_replace_all(pattern = "_", replacement = " ")
  v_vj_usage$Donor <- str_split(v_vj_usage$Sample, "/", simplify = T)[,1]
  v_vj_usage$`V Gene` <- str_split(v_vj_usage$`V`, "\\*", simplify = T)[,1]
  v_vj_usage$`J Gene` <- str_split(v_vj_usage$`J`, "\\*", simplify = T)[,1]
  v_vj_usage %>% group_by(Donor, Type, `V Gene`) %>% dplyr::summarise(frequency=length(`V Gene`)) %>% pivot_wider(names_from = "Donor", values_from ="frequency") %>% write.csv(.,paste0(i,"_Donor_V_Frequency_Table.csv"))
  v_vj_pair <- table(v_vj_usage$`V Gene`, v_vj_usage$`J Gene`) %>% as.data.frame()
  names(v_vj_pair) <- c("V", "J", "Frequency")
  write.csv(v_vj_pair, paste0(i,"_Donor_VJ_Frequency_Table.csv"))
}


# Heatmap
Donor1_VJ <- read.csv("Donor1_Enriched_Motifs.csv")
png("VJ_Usage_Overlap.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
pheatmap(table(Donor1_VJ$V, Donor1_VJ$J), color=colorRampPalette(c("blue","yellow", "red"))((100)), name = "Frequency",display_numbers = T, number_color = "black", main = "VJ Usage Overlap Frequency")                              
dev.off()



##### Circos Plots #####
donor <- read_xlsx("Donor1_Enriched_Motifs.xlsx") %>% mutate(Donor= str_split(sample_name, "_", simplify = T)[,1])


# Name  manipulation if required
donor$Sample <- case_when(
  donor$sample_name %like%  "TCM" ~ "TCM",
  donor$sample_name %like%  "TEM" ~ "TEM",
  donor$sample_name %like%  "Naive" ~ "Naive",
  TRUE ~ NA
)



# Inner join or create an adjacency matrix
x <- donor %>% inner_join(donor, by="amino_acid") %>% filter(!(Sample.x=="TCM"&Sample.y=="TCM"))

#df_adj <- get.adjacency(graph_from_data_frame(data.frame(x$Sample.x, x$Sample.y), directed=FALSE)) %>% as.matrix()


grid.col = structure( c("red", "black", "blue", "orange"),
                      names = c("ApoB AIM +", "TCM", "TEM", "Naive"))

x$Sample.y <- factor(x$Sample.y, levels=c("ApoB AIM +","TCM", "TEM", "Naive" ))

transparency <- ifelse(x$Sample.x==x$Sample.y,1,0.5)

png("RibbonPlot_Donor1.png", width = 5, height=5, units = "in", res = 400,bg = "white")
chordDiagramFromDataFrame(data.frame(FROM=x$Sample.x, TO=x$Sample.y), scale =T, directional = 0, grid.col = grid.col, 
                          link.visible = x$Sample.x %in% "TCM")
dev.off()



#####

##### Chord Plots #####
donor1 <- read_xlsx("donor1_filtered.xlsx") %>% mutate(`V Gene`=str_split(v_resolved, "\\*", simplify = T)[,1], `J Gene`=str_split(j_resolved, "\\*", simplify = T)[,1])%>% filter(`V Gene` %in% target_v_genes) 
donor2 <- read_xlsx("donor2_filtered.xlsx")%>% mutate(`V Gene`=str_split(v_resolved, "\\*", simplify = T)[,1], `J Gene`=str_split(j_resolved, "\\*", simplify = T)[,1])%>% filter(`V Gene` %in% target_v_genes) 
donor3 <- read_xlsx("donor3_filtered.xlsx")%>% mutate(`V Gene`=str_split(v_resolved, "\\*", simplify = T)[,1], `J Gene`=str_split(j_resolved, "\\*", simplify = T)[,1])%>% filter(`V Gene` %in% target_v_genes) 
donor4 <- read_xlsx("donor4_filtered.xlsx")%>% mutate(`V Gene`=str_split(v_resolved, "\\*", simplify = T)[,1], `J Gene`=str_split(j_resolved, "\\*", simplify = T)[,1])%>% filter(`V Gene` %in% target_v_genes) 
donor5 <- read_xlsx("donor5_filtered.xlsx")%>% mutate(`V Gene`=str_split(v_resolved, "\\*", simplify = T)[,1], `J Gene`=str_split(j_resolved, "\\*", simplify = T)[,1])%>% filter(`V Gene` %in% target_v_genes) 
donor6 <- read_xlsx("donor6_filtered.xlsx")%>% mutate(`V Gene`=str_split(v_resolved, "\\*", simplify = T)[,1], `J Gene`=str_split(j_resolved, "\\*", simplify = T)[,1])%>% filter(`V Gene` %in% target_v_genes) 



target_v_genes <- c("TRBV6-5", "TRBV7-2", "TRBV19", "TRBV11-2", "TRBV12-3", "TRBV13-1", "TRBV5-7", "TRBV4-1") # Preferred V-genes across donors


grid.col = structure( c("#fe7f00", "#3BCEAC", "black","darkgreen","darkblue", "red","yellow","#540D6E", rep("grey85",length(unique(donor1$`J Gene`)))),
                      names = c(target_v_genes, unique(donor1$`J Gene`)))

png("Donor1_vj_pairs.png", height=12, width=9, units = "in",bg = "white", res = 400)
chordDiagram(data.frame(v_gene=donor1$`V Gene`, j_gene =donor1$`J Gene`), grid.col = grid.col, annotationTrack ="grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(data.frame(v_gene=donor1$`V Gene`, j_gene =donor1$`J Gene`)))))))


circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))}, bg.border = NA)
dev.off()




#####

