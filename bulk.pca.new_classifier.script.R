library(edgeR)
library(stringr)
library(ggplot2)
library(scattermore)
library(cowplot)
library(tidyr)
library(data.table)

# load("/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/testrun.m.rda")
# meta.df<-fread("/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/samples.tsv",data.table=F)
# genes.df <- rtracklayer::readGFF("/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/gencode.vM25.chr_patch_hapl_scaff.annotation.gff3", filter = list(type = "gene"))
# bulk<-read.table("/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/recount3 mouse data/bulk and sc samples/bulk.samples")
# scRNA<-read.table("/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/recount3 mouse data/bulk and sc samples/sc.samples")
# source("/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/recount3 mouse data/load.studies.R")

genes.df <- rtracklayer::readGFF("/users/kchen1/recount3/mouse_expression/gencode.vM25.chr_patch_hapl_scaff.annotation.gff3", filter = list(type = "gene"))
meta.df<-fread("/users/kchen1/recount3/mouse_expression/recount3 mouse data/mouse sample metadata/samples.tsv",data.table=F)
data.m<-readRDS("/users/kchen1/recount3/mouse_expression/protein.coding.genes.updated.Rds")
bulk<-read.table("/users/kchen1/recount3/mouse_expression/recount3 mouse data/bulk and sc samples/bulk.samples")
scRNA<-read.table("/users/kchen1/recount3/mouse_expression/recount3 mouse data/bulk and sc samples/sc.samples")
source("/users/kchen1/recount3/mouse_expression/recount3 mouse data/load.studies.R")

colnames(bulk)[1]<-"Run"
colnames(scRNA)[1]<-"Run"
meta.df$tissue<-NA
meta.df$sample.type<-NA  #
meta.df$cell.type<-NA
meta.df$cell.line<-NA
#source("/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/recount3 mouse data/bulk.labeling.R")
source("/users/kchen1/recount3/mouse_expression/recount3 mouse data/bulk.labeling.R")

# these are sample types that were pre-labeled by Chris
meta.df$sample.type[which(meta.df$run_acc %in% bulk$Run)]<-"bulk"      
meta.df$sample.type[which(meta.df$run_acc %in% scRNA$Run)]<-"scRNA"

meta.df<-meta.df[match(colnames(data.m),meta.df$run_acc),]
genes.df<-genes.df[match(rownames(data.m),genes.df$ID),]
dge.o <- DGEList(counts = data.m, genes = genes.df, samples = meta.df)

rm(list = c("genes.df", "data.m","meta.df"))
gc()

#filter out the lowly expressed genes
filter.idx <- rowMeans(dge.o$counts) > 0.5
dge.o <- dge.o[filter.idx, ]

zero.pct.v <-  colMeans(dge.o$counts < 1)

dge.o$samples$zero.pct<-zero.pct.v

bulk.idx<-which(dge.o$samples$sample.type=="bulk")
scRNA.idx<-which(dge.o$samples$sample.type=="scRNA")


# compute the probability of each sample
bulk.dens<-density(dge.o$samples$zero.pct[bulk.idx])
scRNA.dens<-density(dge.o$samples$zero.pct[scRNA.idx])
bulk.f<-approxfun(bulk.dens)
scRNA.f<-approxfun(scRNA.dens)
bulk.bayes<-sapply(dge.o$samples$zero.pct,function(x){
  f<-integrate(bulk.f,lower=x-0.001,upper=x+0.001)$value
  g<-integrate(scRNA.f,lower=x-0.001,upper=x+0.001)$value
  f/(f + g)
})

dge.o$samples$bulk.bayes<-bulk.bayes

dge.o$samples$pred.type<-NA
dge.o$samples$pred.type[dge.o$samples$bulk.bayes>=0.825]<-"bulk"
dge.o$samples$pred.type[dge.o$samples$bulk.bayes<=0.1]<-"scRNA"


library(MASS)
pred.type<-dge.o$samples$pred.type[!is.na(dge.o$samples$pred.type)]
real.type<-dge.o$samples$sample.type[!is.na(dge.o$samples$pred.type)]


confusion<-table(pred.type,real.type)
correctness<-(confusion[1,1]+confusion[2,2])/sum(confusion)

# write.table(as.matrix(confusion), "/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/Sept:18:2020_test _file/confusion_table",quote=F)
# write(correctness,file="confusion_table",append=T)
# write(paste0("percentage of all samples kept: ",sum(confusion)/nrow(dge.o$samples)),file="confusion_table",append=T)
# write(paste0("percentage of bulk samples kept: ",sum(pred.type=="bulk")/nrow(dge.o$samples)),file="confusion_table",append=T)
# write(paste0("correctness among bulk samples kept: ",confusion[1,1]/(confusion[1,1]+confusion[1,2])),file="confusion_table",append=T)
# write(paste0("percentage of scRNA samples kept: ", sum(pred.type=="scRNA")/nrow(dge.o$samples)),file="confusion_table",append=T)
# write(paste0("correctness among scRNA samples kept: ", confusion[2,2]/(confusion[2,2]+confusion[2,1])),file="confusion_table",append=T)

write.table(as.matrix(confusion), "/users/kchen1/recount3/mouse_expression/sept_18_2020/confusion_table",quote=F)
write(correctness,file="/users/kchen1/recount3/mouse_expression/sept_18_2020/confusion_table",append=T)
write(paste0("percentage of all samples kept: ",sum(confusion)/nrow(dge.o$samples)),file="/users/kchen1/recount3/mouse_expression/sept_18_2020/confusion_table",append=T)
write(paste0("percentage of bulk samples kept: ",sum(pred.type=="bulk")/nrow(dge.o$samples)),file="/users/kchen1/recount3/mouse_expression/sept_18_2020/confusion_table",append=T)
write(paste0("correctness among bulk samples kept: ",confusion[1,1]/(confusion[1,1]+confusion[1,2])),file="/users/kchen1/recount3/mouse_expression/sept_18_2020/confusion_table",append=T)
write(paste0("percentage of scRNA samples kept: ", sum(pred.type=="scRNA")/nrow(dge.o$samples)),file="/users/kchen1/recount3/mouse_expression/sept_18_2020/confusion_table",append=T)
write(paste0("correctness among scRNA samples kept: ", confusion[2,2]/(confusion[2,2]+confusion[2,1])),file="/users/kchen1/recount3/mouse_expression/sept_18_2020/confusion_table",append=T)

#pdf(file="/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/sept_18_2020/prob_dist.pdf")
pdf(file="/users/kchen1/recount3/mouse_expression/sept_18_2020/prob_dist.pdf")
hist(dge.o$samples$bulk.bayes, breaks=500, main="histogram of the probability of the samples")
dev.off()

# classify only the samples with probability >= 0.85 as bulk
idx <- which(dge.o$samples$bulk.bayes>=0.825)
dge.o <- dge.o[, idx]

sample.idx <- which(colSums(dge.o$counts) > 10^6)    # get rid of small sample sizes

dge.o <- dge.o[, sample.idx, keep.lib.sizes = TRUE]

keep.idx <- (aveLogCPM(dge.o, prior.count = 2) > 1) & (dge.o$genes$seqid %in% str_c("chr", 1:19))
length(keep.idx)   
sum(keep.idx) 
dge.o <- dge.o[keep.idx, , keep.lib.sizes=FALSE]

cpm.m <- cpm(dge.o, normalized.lib.sizes = TRUE, log = FALSE)
log2cpm.m <- log2(cpm.m + 1)

rm(cpm.m)
gc()

library(rsvd)
set.seed(1)
scale.m <- t(scale(t(log2cpm.m), center = TRUE, scale = FALSE))
rpca.o <- rpca(A = scale.m, k = 30, center = FALSE, scale = FALSE)
expVar.v <- summary(rpca.o )[3, ]

save(rpca.o,file="/users/kchen1/recount3/mouse_expression/sept_18_2020/rpca.o.Rda")

rm(scale.m)
gc()

# figure for variance explained by first 30 PCs
expVar.df <- data.frame(fraction = expVar.v, id = seq_along(expVar.v))

fvar.point.p <- ggplot(expVar.df, aes(x = id, y = fraction)) + 
  geom_line(size = 0.2, color = "#86BBD8", linetype = "dashed") + 
  geom_point( shape = 18, size = 1, color = "#86BBD8") + 
  scale_x_continuous(name = "Principal component", breaks = seq(1, length(expVar.v), 5), expand = c(0, 0), limits = c(0.5, 30.5)) + 
  labs(title = str_c("PCA on protein coding genes (n=", nrow(rpca.o$rotation),")"), y = "Fraction of explained variance")  + 
  guides(color = F) + 
  # panel_border(colour="black", linetype = 1, size= 0.5) +
  theme_bw(base_size = 8) + 
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.25),
        axis.line.x = element_line(colour="black", linetype = 1, size= 0.25),
        axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
        axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        plot.title = element_text(face = "plain", size = 8, hjust = 0.85),
        plot.margin = unit(c(5, 5, 6, 6), "pt"),
        legend.background = element_rect(color= "black", linetype = 1, size= 0.2, fill = NA),
        legend.position = c(0.3, 0.5),
        legend.justification = c(0, 0),
        legend.key.size = unit(8, "pt"),
        legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text.align = 0.5,
        legend.text = element_text(size = 5),
        strip.text = element_text(color = "white"),
        strip.background = element_rect(fill = "black", color = "black"))


#save_plot("/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/final docs/bulk.fvar.pdf", fvar.point.p,base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 1)
save_plot("/users/kchen1/recount3/mouse_expression/sept_18_2020/bulk.fvar.pdf", fvar.point.p,base_height = 2, base_width = 2*1.2, nrow = 1, ncol = 1)

tmp.df <- data.frame(pc1 = rpca.o$rotation[, 1], pc2 = rpca.o$rotation[, 2], 
                     pc3 = rpca.o$rotation[, 3], pc4 = rpca.o$rotation[, 4],
                     libsize = log2(dge.o$samples$lib.size),
                     Study=dge.o$samples$study_acc, cell.type=dge.o$samples$cell.type,
                     cell.line=dge.o$samples$cell.line,
                     tissue = dge.o$samples$tissue, zero.pct=dge.o$samples$zero.pct)

tmp.df$tissue[(tmp.df$tissue %in% names(table(tmp.df$tissue)[which(table(tmp.df$tissue) < 10)])) ] <- NA
tmp.df$tissue<-factor(tmp.df$tissue)

save(tmp.df,file="/users/kchen1/recount3/mouse_expression/sept_18_2020/tmp.df.Rda")

pc12_libsize.point.p <- ggplot(tmp.df, aes(x = pc1, y = pc2, color = libsize)) +
  geom_scattermore(pointsize = 0.5, shape = 1, alpha = 0.7) +
  #coord_cartesian(xlim=c(-0.0075,0.005),ylim=c(-0.01,0.004)) +
  scale_color_gradient(low = "grey", high = "#F2712B", name = expression(paste(log[2], "(lib.size)"))) +
  labs( y = str_c("PC2: ", sprintf(expVar.v[2]*100, fmt = '%#.1f'), "%"), x = str_c("PC1: ", sprintf(expVar.v[1]*100, fmt = '%#.1f'), "%"),
        title = str_c("Bulk samples (n=", sum(!is.na(tmp.df$libsize)),")")) +
  theme_bw(base_size = 8) + 
  theme(panel.border = element_rect(colour="black", linetype = 1, size= 0.25),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.25),
        axis.line.x = element_line(colour="black", linetype = 1, size= 0.25),
        axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
        axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(size = 7),
        plot.title = element_text(face = "plain", size = 8, hjust = 0.5),
        plot.margin = unit(c(5, 5, 6, 6), "pt"),
        legend.background = element_rect(color= "black", linetype = 1, size= 0.2, fill = NA),
        legend.position = "right",
        legend.justification = c(0.5, 0.5),
        legend.key.size = unit(8, "pt"),
        legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
        legend.title = element_text(size = 5),
        legend.title.align = 0.5,
        legend.text.align = 0.5,
        legend.text = element_text(size = 5),
        strip.text = element_text(color = "white"),
        strip.background = element_rect(fill = "black", color = "black"))

pc34_libsize.point.p <- ggplot(tmp.df, aes(x = pc3, y = pc4, color = libsize)) +
  geom_scattermore(pointsize = 0.5, shape = 1, alpha = 0.7) +
  scale_color_gradient(low = "grey", high = "#F2712B", name = expression(paste(log[2], "(lib.size)"))) +
  labs( y = str_c("PC4: ", sprintf(expVar.v[4]*100, fmt = '%#.1f'), "%"), x = str_c("PC3: ", sprintf(expVar.v[3]*100, fmt = '%#.1f'), "%"),
        title = str_c("Bulk samples (n=", sum(!is.na(tmp.df$libsize)),")")) +
  theme_bw(base_size = 8) + 
  theme(panel.border = element_rect(colour="black", linetype = 1, size= 0.25),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.25),
        axis.line.x = element_line(colour="black", linetype = 1, size= 0.25),
        axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
        axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(size = 7),
        plot.title = element_text(face = "plain", size = 8, hjust = 0.5),
        plot.margin = unit(c(5, 5, 6, 6), "pt"),
        legend.background = element_rect(color= "black", linetype = 1, size= 0.2, fill = NA),
        legend.position = "right",
        legend.justification = c(0.5, 0.5),
        legend.key.size = unit(8, "pt"),
        legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
        legend.title = element_text(size = 5),
        legend.title.align = 0.5,
        legend.text.align = 0.5,
        legend.text = element_text(size = 5),
        strip.text = element_text(color = "white"),
        strip.background = element_rect(fill = "black", color = "black"))


colors.llv <- readRDS("/users/kchen1/recount3/mouse_expression/colors.llv.rds")
pc12_tissue.point.p <- ggplot(tmp.df, aes(x = pc1, y = pc2, color = tissue)) +
  geom_scattermore(pointsize = 0.6, shape = 1,alpha=0.7) +
  #coord_cartesian(xlim=c(-0.0075,0.005),ylim=c(-0.01,0.004)) +
  scale_color_manual(name = "Bulk tissue", values = colors.llv$color50[[2]][seq_len(nlevels(tmp.df$tissue))]) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1, shape = 19))) + 
  labs( y = str_c("PC2: ", sprintf(expVar.v[2]*100, fmt = '%#.1f'), "%"), x = str_c("PC1: ", sprintf(expVar.v[1]*100, fmt = '%#.1f'), "%"),
        title = str_c("Bulk tissue (n=", sum(!is.na(tmp.df$tissue)),")")) +
  theme_bw(base_size = 8) + 
  theme(panel.border = element_rect(colour="black", linetype = 1, size= 0.25),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.25),
        axis.line.x = element_line(colour="black", linetype = 1, size= 0.25),
        axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
        axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(size = 7),
        plot.title = element_text(face = "plain", size = 8, hjust = 0.5),
        plot.margin = unit(c(5, 5, 6, 6), "pt"),
        legend.background = element_rect(color= "black", linetype = 1, size= 0.2, fill = NA),
        legend.position = "right",
        legend.justification = c(0.5, 0.5),
        legend.key.size = unit(6, "pt"),
        legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
        legend.title = element_text(size = 4),
        legend.title.align = 0.5,
        legend.text.align = 0.5,
        legend.text = element_text(size = 3),
        strip.text = element_text(color = "white"),
        strip.background = element_rect(fill = "black", color = "black"))

set.seed(101)
pc34_tissue.point.p <- ggplot(tmp.df, aes(x = pc3, y = pc4, color = tissue)) +
  geom_scattermore(pointsize = 0.6, shape = 1) +
  scale_color_manual(name = "Bulk tissue", values = colors.llv$color50[[2]][seq_len(nlevels(tmp.df$tissue))]) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 1, shape = 19))) + 
  labs( y = str_c("PC4: ", sprintf(expVar.v[4]*100, fmt = '%#.1f'), "%"), x = str_c("PC3: ", sprintf(expVar.v[3]*100, fmt = '%#.1f'), "%"),
        title = str_c("Bulk tissue (n=", sum(!is.na(tmp.df$tissue)),")")) +
  theme_bw(base_size = 8) + 
  theme(panel.border = element_rect(colour="black", linetype = 1, size= 0.25),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.25),
        axis.line.x = element_line(colour="black", linetype = 1, size= 0.25),
        axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
        axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(size = 7),
        plot.title = element_text(face = "plain", size = 8, hjust = 0.5),
        plot.margin = unit(c(5, 5, 6, 6), "pt"),
        legend.background = element_rect(color= "black", linetype = 1, size= 0.2, fill = NA),
        legend.position = "right",
        legend.justification = c(0.5, 0.5),
        legend.key.size = unit(6, "pt"),
        legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
        legend.title = element_text(size = 4),
        legend.title.align = 0.5,
        legend.text.align = 0.5,
        legend.text = element_text(size = 3),
        strip.text = element_text(color = "white"),
        strip.background = element_rect(fill = "black", color = "black"))

pca.mp.1 <- plot_grid(pc12_libsize.point.p, pc12_tissue.point.p, nrow = 2, ncol = 1, label_size = 10, labels = NULL, align = "hv", axis = "tblr")
pca.mp.2 <- plot_grid(pc34_libsize.point.p, pc34_tissue.point.p, nrow = 2, ncol = 1, label_size = 10, labels = NULL, align = "hv", axis = "tblr")
#save_plot("/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/final docs/bulk.pca.pdf", pca.mp,base_height = 2, base_width = 2*1.8, nrow = 2, ncol = 2)
save_plot("/users/kchen1/recount3/mouse_expression/sept_18_2020/bulk.pca.1.pdf", pca.mp.1,base_height = 3, base_width = 2*2.8, nrow = 2, ncol = 1)
save_plot("/users/kchen1/recount3/mouse_expression/sept_18_2020/bulk.pca.2.pdf", pca.mp.2,base_height = 3, base_width = 2*2.8, nrow = 2, ncol = 1)

pc12_tissue.point.lp <- lapply(levels(tmp.df$tissue), function(tis) {
  plot.df <- tmp.df
  plot.df$color <- NA_character_
  plot.df$color[which(plot.df$tissue == tis)] <- as.character(plot.df$tissue[which(plot.df$tissue == tis)])
  
  p <- ggplot(data = dplyr::filter(plot.df, is.na(plot.df$color)), aes(x = pc1, y = pc2, color = color)) +
    geom_scattermore(pointsize = 0.6, shape = 1, alpha = 0.8) +
    geom_scattermore(data = dplyr::filter(plot.df, color == tis), pointsize = 0.6, shape = 1, alpha = 1) +
    scale_color_manual(name = "", na.value = "grey", values = colors.llv$color50[[2]][match(tis, levels(tmp.df$tissue))], guide = FALSE) +
    
    labs( y = str_c("PC2: ", sprintf(expVar.v[2]*100, fmt = '%#.1f'), "%"), x = str_c("PC1: ", sprintf(expVar.v[1]*100, fmt = '%#.1f'), "%"),
          title = str_c(tis, " (n=", sum(tmp.df$tissue == tis, na.rm = TRUE),")")) +
    ylim(range(tmp.df$pc2)) + xlim(range(tmp.df$pc1)) +
    theme_bw(base_size = 8) + 
    theme(panel.border = element_rect(colour="black", linetype = 1, size= 0.25),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.25),
          axis.line.x = element_line(colour="black", linetype = 1, size= 0.25),
          axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
          axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 7),
          plot.title = element_text(face = "plain", size = 8, hjust = 0.5),
          plot.margin = unit(c(5, 5, 6, 6), "pt"),
          legend.background = element_rect(color= "black", linetype = 1, size= 0.2, fill = NA),
          legend.position = "right",
          legend.justification = c(0.5, 0.5),
          legend.key.size = unit(6, "pt"),
          legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
          legend.title = element_text(size = 4),
          legend.title.align = 0.5,
          legend.text.align = 0.5,
          legend.text = element_text(size = 3),
          strip.text = element_text(color = "white"),
          strip.background = element_rect(fill = "black", color = "black"))
  return(p)
})

tissue.mp <- plot_grid(plotlist = pc12_tissue.point.lp,
                       nrow = 6, ncol = 8, label_size = 10, labels = NULL, align = "hv", axis = "tblr")
#save_plot("/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/final docs/bulk.tissue.pca.facets.pdf", tissue.mp,base_height = 2, base_width = 2*1.4, nrow = 6, ncol = 8)
save_plot("/users/kchen1/recount3/mouse_expression/sept_18_2020/bulk.tissue.pca.facets.pdf", tissue.mp,base_height = 2, base_width = 2*1.4, nrow = 6, ncol = 8)

pc12_tissue.point.studies.lp <- lapply(levels(tmp.df$tissue), function(tis) {
  plot.df <- tmp.df
  plot.df$color <- NA_character_
  plot.df$color[which(plot.df$tissue == tis)] <- as.character(plot.df$Study[which(plot.df$tissue == tis)])
  tmp.levels<-unique(as.character(dplyr::filter(plot.df,!is.na(plot.df$color))$Study))
  p <- ggplot(data = dplyr::filter(plot.df, is.na(plot.df$color)), aes(x = pc1, y = pc2, color = color)) +
    geom_scattermore(pointsize = 0.6, shape = 1, alpha = 0.8) +
    geom_scattermore(data = dplyr::filter(plot.df, color %in% tmp.levels), pointsize = 0.6, shape = 1, alpha = 1) +
    scale_color_manual(name = "Studies", na.value = "grey", values = colors.llv$color50[[2]][seq_len(length(tmp.levels))]) +
    
    labs( y = str_c("PC2: ", sprintf(expVar.v[2]*100, fmt = '%#.1f'), "%"), x = str_c("PC1: ", sprintf(expVar.v[1]*100, fmt = '%#.1f'), "%"),
          title = str_c(tis, " (n=", sum(tmp.df$tissue == tis, na.rm = TRUE),")")) +
    ylim(range(tmp.df$pc2)) + xlim(range(tmp.df$pc1)) +
    theme_bw(base_size = 8) + 
    theme(panel.border = element_rect(colour="black", linetype = 1, size= 0.25),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.25),
          axis.line.x = element_line(colour="black", linetype = 1, size= 0.25),
          axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
          axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 7),
          plot.title = element_text(face = "plain", size = 8, hjust = 0.5),
          plot.margin = unit(c(5, 5, 6, 6), "pt"),
          legend.background = element_rect(color= "black", linetype = 1, size= 0.2, fill = NA),
          legend.position = "right",
          legend.justification = c(0.5, 0.5),
          legend.key.size = unit(6, "pt"),
          legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
          legend.title = element_text(size = 4),
          legend.title.align = 0.5,
          legend.text.align = 0.5,
          legend.text = element_text(size = 3),
          strip.text = element_text(color = "white"),
          strip.background = element_rect(fill = "black", color = "black"))
  return(p)
})

tissue.studies.mp <- plot_grid(plotlist = pc12_tissue.point.studies.lp,
                               nrow = 6, ncol = 8, label_size = 10, labels = NULL, align = "hv", axis = "tblr")
#save_plot("/Users/kevinchen/Documents/MS Bioinformatic program courses and supplement courses/Hansen Lab/Recount_3/Big Mouse Studies/final docs/bulk.tissue.pca.facets.studies.pdf", tissue.studies.mp,base_height = 2, base_width = 2*1.4, nrow = 6, ncol = 8)
save_plot("/users/kchen1/recount3/mouse_expression/sept_18_2020/bulk.tissue.pca.facets.studies.pdf", tissue.studies.mp,base_height = 2, base_width = 2*1.4, nrow = 6, ncol = 8)


