library(ggplot2)
library(ggfortify)
library(mixOmics)
#install.packages('heatmaply')

#badNames_full <- list.files('C:/Users/Sid/OneDrive/Research/Multiomics/drought_root/diff_genes_subtypes_50_with_Pheno_kegg', full.names = T)
#badNames_full
#setwd("C:/Users/ssz74/OneDrive/Research/Multiomics/TCGA_new/28_nov_2018/second_run")

#upload mRNA maize_drought
mRNA_drought <- read.csv('diff_genes_field_rna_seq.txt',header=TRUE,sep='\t',row.names=1, check.names=FALSE)
mRNA_drought[1:10,1:10]

##filtering
mRNA_drought_large <- data.frame(mRNA_drought)
dim(mRNA_drought_large) #45736    18
mRNA_drought_large[1:10,1:10]
keep <- rowSums(mRNA_drought_large[1:3838,]) != 0
mRNA_drought_filt <- mRNA_drought_large[keep,]
dim(mRNA_drought_filt) ##33632    18
##
mRNA_drought_dt <- data.frame(t(mRNA_drought_filt))
mRNA_drought_dt[1:10,1:10]

##upload protein drought
#prot_BRCA <- read.csv('BRCA_protein_small.txt',header=TRUE,sep='\t',row.names=1, check.names=FALSE)
#prot_BRCA_dt <- data.frame(prot_BRCA)
#prot_BRCA_dt[1:10,1:10]
#swap_2 <- data.frame(t(prot_BRCA_dt))
#prot_BRCA_dt <- data.frame(t(swap_2))
#prot_BRCA_dt[1:10,1:10]
#rm(swap_2)

##upload metabolite
metab_drought <- read.csv('field_metab_kegg.txt',header=TRUE,sep='\t',row.names=1, check.names=FALSE)
metab_drought_t <- data.frame(metab_drought[1:18])
metab_drought_t[1:10,1:10]
metab_drought_dt <- data.frame(t(metab_drought_t))
metab_drought_dt[1:10,1:10]

##upload Lab Phenotype info
pheno_drought <- read.csv('upd_lab_phenotype_data.txt',header=TRUE,sep='\t',row.names=1, check.names=FALSE)
pheno_drought_t <- data.frame(pheno_drought)
pheno_drought_t[1:10,1:10]
pheno_drought_dt <- data.frame(t(pheno_drought_t))
pheno_drought_dt[1:10,1:10]


##info files
info_drought <- read.csv('drought_root_info.txt',header=TRUE,sep='\t',row.names=1, check.names=FALSE)
info_drought_dt <- data.frame(info_drought)
info_drought_dt[1:10,]

##
dim(mRNA_drought_dt)
#mRNA_BRCA_dt[1:10,1:10]
#dim(prot_drought_dt)
#prot_BRCA_dt[1:10,1:10]
dim(metab_drought_dt)
dim(pheno_drought_dt)
#miRNA_BRCA_dt[1:10,1:10]

#data_brca = list(mRNA = mRNA_BRCA_dt,
#                 miRNA = miRNA_BRCA_dt,
#                 proteomics = prot_BRCA_dt)
data_drought = list(mRNA = mRNA_drought_dt,
                 metab = metab_drought_dt,
                 pheno = pheno_drought_dt)



lapply(data_drought, dim)

Y_drought = info_drought_dt$tip_region
#Y_read = info_READ_dt$Sample_Type
summary(Y_drought)
Y_drought

##pca mRNA
Sample_Type <- info_drought_dt$tip_region
dim(mRNA_drought_dt)
mRNA_drought_dtinfo <- cbind(mRNA_drought_dt,Sample_Type)
dim(mRNA_drought_dtinfo)
#jpeg('pca_BRCA_mRNA.jpg',margin=c(8,2))
autoplot(prcomp(mRNA_drought_dtinfo[1:3836]), data = mRNA_drought_dtinfo,colour = "Sample_Type", label = FALSE, label.size = 1,frame = TRUE) ##<<< THIS WORKS
#dev.off()

##pca metab
Sample_Type <- info_drought_dt$tip_region
dim(metab_drought_dt)
metab_drought_dtinfo <- cbind(metab_drought_dt,Sample_Type)
dim(metab_drought_dtinfo)
#jpeg('pca_BRCA_miRNA.jpg',margin=c(8,2))
autoplot(prcomp(metab_drought_dtinfo[1:569]), data = metab_drought_dtinfo,colour = "Sample_Type", label = FALSE, label.size = 4,frame = TRUE) ##<<< THIS WORKS
#dev.off()

##pca pheno
Sample_Type <- info_drought_dt$tip_region
dim(pheno_drought_dt)
pheno_drought_dtinfo <- cbind(pheno_drought_dt,Sample_Type)
dim(pheno_drought_dtinfo)
#jpeg('pca_BRCA_miRNA.jpg',margin=c(8,2))
autoplot(prcomp(pheno_drought_dtinfo[1:10]), data = pheno_drought_dtinfo,colour = "Sample_Type", label = FALSE, label.size = 4,frame = TRUE) ##<<< THIS WORKS
dev.off()

##pca protein
#Sample_Type <- info_BRCA_dt$Type
#dim(prot_BRCA_dt)
#prot_BRCA_dtinfo <- cbind(prot_BRCA_dt,Sample_Type)
#dim(prot_BRCA_dtinfo)
#autoplot(prcomp(prot_BRCA_dtinfo[1:224]), data = prot_BRCA_dtinfo,colour = "Sample_Type", label = FALSE, label.size = 4) ##<<< THIS does not WORK

##matrix design

design = matrix(0.1, ncol = length(data_drought), nrow = length(data_drought), 
                dimnames = list(names(data_drought), names(data_drought)))  
diag(design) = 0

design 

##tuning
sgccda.res = block.splsda(X = data_drought, Y = Y_drought, ncomp = 5, 
                          design = design)

sgccda.res
#In cor(A[[k]], variates.A[[k]]) : the standard deviation is zero

set.seed(123)

perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 3, nrepeat = 10)

##plot error rates - multiple curves
#jpeg('error_rate_BRCA.jpg',margin=c(8,2))
plot(perf.diablo)
#dev.off()

perf.diablo$choice.ncomp$WeightedVote

ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
#ncomp = 3
ncomp

##### THIS WILL TAKE TIME
set.seed(123) # for reproducibility, only when the `cpus' argument is not used
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,50,5)),
                   metab = c(5:9, seq(10, 18, 2), seq(20,50,5)),
                   pheno = c(5:9, seq(10, 18, 2), seq(20,50,5)))

#test.keepX = list (mRNA = c(seq(2, 100, 2)),
#                   miRNA = c(seq(2, 100, 2)),
#                   proteomics = c(seq(2, 100, 2)))

t1 = proc.time()
tune.drought = tune.block.splsda(X = data_drought, Y = Y_drought, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 3, nrepeat = 5,
                              cpus = 4, dist = "centroids.dist")
t2 = proc.time()
running_time = t2 - t1; running_time

list.keepX = tune.drought$choice.keepX
list.keepX


##save tuning
#save(tune.TCGA,list.keepX, file = 'BRCA_analysis_tune.RData')
#save.image(file = "drought_root_50_work_space.RData")
#load('BRCA_analysis_tune.RData')
#tune.TCGA$choice.keepX
#tune.TCGA$choice.keepX$mRNA

## ------------------------------------------------------------------------
sgccda.res = block.splsda(X = data_drought, Y = Y_drought, ncomp = ncomp, 
                          keepX = list.keepX, design = design)
#sgccda.res   # list the different functions of interest related to that object

## ------------------------------------------------------------------------
sgccda.res$design


#save(tune.drought,list.keepX, file = 'BRCA_small_subtype_tune.RData')
save.image(file = "diff_genes_lab_50_KEGG_Pheno_work_space.RData")
## ------------------------------------------------------------------------
# mrna variables selected on component 1

mRNA_result <- data.frame(selectVar(sgccda.res, block = 'mRNA', comp = 1)$mRNA$name)
mRNA_result
#write.csv(mRNA_result,"mra_1.csv",row.names = FALSE)
#write.csv(mRNA_result,"mra_2.csv")

metab_result <- data.frame(selectVar(sgccda.res, block = 'metab', comp = 1)$metab$name)
metab_result
#write.csv(metab_result,"metab_2.csv")

pheno_result <- data.frame(selectVar(sgccda.res, block = 'pheno', comp = 1)$pheno$name)
pheno_result
#write.csv(pheno_result,"pheno_2.csv")
####comp2####
mRNA_result <- data.frame(selectVar(sgccda.res, block = 'mRNA', comp = 2)$mRNA$name)
mRNA_result
#write.csv(mRNA_result,"mra_1.csv",row.names = FALSE)

metab_result <- data.frame(selectVar(sgccda.res, block = 'metab', comp = 2)$metab$name)
metab_result

pheno_result <- data.frame(selectVar(sgccda.res, block = 'pheno', comp = 2)$pheno$name)
pheno_result
#####


#selectVar(sgccda.res, block = 'proteomics', comp = 1)$proteomics$name

## ------------------------------------------------------------------------
#jpeg('plotdiablo_comp1_BRCA.jpg',margin=c(8,2))
#plotDiablo(sgccda.res, ncomp = 2,legend.ncol = 2)
#dev.off()

## ------------------------------------------------------------------------
#jpeg('sgccda_BRCA.jpg',margin=c(8,2))
#plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
#dev.off()
#plotIndiv(sgccda.res, ind.names = TRUE, legend = TRUE, title = 'DIABLO')

## ------------------------------------------------------------------------
#jpeg('arrowPlot_BRCA.jpg',margin=c(8,2))
#plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
#dev.off()

## ------------------------------------------------------------------------
#jpeg('varplot_BRCA.jpg',margin=c(8,2))
plotVar(sgccda.res, var.names = TRUE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17, 15), cex = c(1,1.2,1), col = c('darkgreen','brown1','blue')) #, 'lightgreen darkorchid'
?plotVar

circleplot_1 <- plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
                        pch = c(16, 17, 15), cex = c(2,1,2), col = c('darkgreen','brown1','lightgreen'))

circleplot_1 ## has X & y coordinates

#write.table(circleplot_1,file="circle_plot_XY.txt",sep = "\t",row.names = T,col.names = T)

####correlation for all 3 levels####
#plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
#        pch = c(16, 17, 15), cex = c(2,2,2), col = c('darkgreen', 'brown1', 'lightgreen')) #, 'lightgreen darkorchid'
#dev.off()

## ------------------------------------------------------------------------
#jpeg('circos_DE_gene_maize_root.jpg')
##comp_1
circos_comp1 <-circosPlot(sgccda.res,comp=1 ,cutoff = 0.9, line = TRUE, 
            color.blocks= c('darkorchid', 'brown1','green'), #, 'lightgreen'
            color.cor = c("chocolate3","grey20"), size.labels = 0.5,size.variables = 0.7,ncol.legend = 1)

par(mfrow=c(1,1))
circosPlot(sgccda.res,comp=1,cutoff = 0.9, line = TRUE,ncol.legend = 1,showIntraLinks = TRUE)

circos_comp1 <- circosPlot(sgccda.res,comp=1 ,cutoff = 0.0, line = TRUE, 
                           color.blocks= c('darkorchid', 'brown1','green'), #, 'lightgreen'
                           color.cor = c("chocolate3","grey20"), size.labels = 0.5,size.variables = 0.7,ncol.legend = 1)

?circosPlot

library(corrplot)
library(RColorBrewer)

#M <-cor(mtcars)
corrplot(circos_comp1, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

circosPlot(sgccda.res,comp=2 ,cutoff = 0.7, line = TRUE, 
           color.blocks= c('darkorchid', 'brown1','green'), #, 'lightgreen'
           color.cor = c("chocolate3","grey20"), size.labels = 0.5,size.variables = 0.7,ncol.legend = 1)

circos_comp2 <- circosPlot(sgccda.res,comp=2 ,cutoff = 0.0, line = TRUE, 
                           color.blocks= c('darkorchid', 'brown1','green'), #, 'lightgreen'
                           color.cor = c("chocolate3","grey20"), size.labels = 0.5,size.variables = 0.7,ncol.legend = 1)

corrplot(circos_comp2, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))


write.table(circos_comp1,file = "circos_comp1_matrix.txt",sep = "\t",col.names = T, row.names = T)
write.table(circos_comp2,file = "circos_comp2_matrix.txt",sep = "\t",col.names = T, row.names = T)
#dev.off()

## ---- eval = TRUE--------------------------------------------------------
#jpeg('network_BRCA.jpg')
# my.network = network(sgccda.res, blocks = c(1,2,3), #,3
#         color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.7, #, 'lightgreen'
#         save = 'jpeg',
#         name.save = 'network_de_gene_maize_subtype')
# dev.off()

library(igraph)
my.network = network(sgccda.res, blocks = c(1,2,3), #,3
                     color.node = c('darkorchid', 'brown1','lightgreen'), cutoff = 0.7) #, 'lightgreen'
my.network
write.graph(my.network$gR, file = "DE_gene_maize_50_07.gml", format = "gml")
dev.off()
#jpeg('loading_1.jpg',width = 1280, height = 800)
jpeg('loading_2.jpg',width = 1280, height = 800)
plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median',size.name = 1,size.legend = 1.5,border = TRUE)
loading_1<-plotLoadings(sgccda.res,block = "mRNA" ,comp = 1, contrib = 'max', method = 'median',size.name = 1,size.legend = 1.5,border = TRUE)
write.table(loading_1,file = "loading_1_metab.txt",sep="\t")
dev.off()

jpeg('heatmap_diff_gene.jpg')
sgccda.res
cimDiablo(sgccda.res,margin=c(10,4),legend = FALSE)
#dev.off()
##custom heatmap
heatmap_1 <- cimDiablo(sgccda.res,margin=c(10,20),legend = FALSE)
write.table(heatmap_1,file = "heatmap_values.txt",sep="\t")

library("RColorBrewer")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap(heatmap_1, scale = "none", col =  col)

library("pheatmap")

pheatmap(t(heatmap_1), cutree_rows = 3,cutree_cols = 3,fontsize_row = 9,fontsize_col= 14, border_color = "grey60")
?pheatmap()
dev.off()

row.names(t(heatmap_1))

library(heatmaply)
#heatmaply(x=t(heatmap_1),color = plasma,Colv = NULL,hclust_method = 'complete', k_row = NA,file = c('int_plot_pheno.html', 'int_plot_pheno.png'))
heatmaply(x=t(heatmap_1),color = col,Colv = NULL,hclust_method = 'complete', k_row = NA,dendrogram = 'both',fontsize_row = 5)
?heatmaply


rbind.data.frame(metab_drought[1,1:3],metab_drought[1,4:6])
