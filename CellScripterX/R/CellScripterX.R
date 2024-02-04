# 导入所需的包
library_list <- c("dplyr", "tidyverse", "patchwork", "reshape2", 
                  "ggstatsplot", "ggcorrplot", "corrplot", "ggplot2", 
                  "pheatmap", "RColorBrewer", "cols4all", "Seurat", 
                  "decontX", "DoubletFinder", "AnnotationHub", 
                  "org.Hs.eg.db", "clusterProfiler", "Rgraphviz", "SingleR", 
                  "harmony", "figpatch","reticulate",'Matrix')

lapply(library_list, library, character.only = TRUE)

# 全局设置
options(connectionObserver = NULL)


celltype_major = list(
  Epi=c('EPCAM', 'KRT5','KRT18'),
  Myeloid=c('CSF1R','CD68','LYZ'),
  Fib=c('COL1A1','LUM','DCN'),
  Endo=c('PECAM1','CLDN5','RAMP2'),
  T=c('CD3D','CD3E','TRAC'),
  B=c('CD19','CD79A','MS4A1'),
  Plasma=c('MZB1','IGHG1','IGHG2'),
  Mast=c('TPSB2','CPA3','TPSAB1'),
  Cycling=c('MKI67','PCNA','TOP2A')
)

T_minor = list(
  CD4T=c('CD4'),
  CD8T=c('CD8A'),
  Treg=c('FOXP3','IL32','TNFRSF18','TNFRSF4'),
  Th17=c('IL17A','IL17F','CD40LG'),
  NaiveT=c('TCF7','SELL','LEF1','CCR7'),
  MemoryT=c('CD8','CD44','CD127','TNFRSF7','CD28'),
  NK=c('NKG7','KLRF1','KLRD1','NCAM1')
)

Myeloid_minor = list(
  CD14_mono=c('CD14'),
  CD16_mono=c('FCGR3A'),
  Macro=c('C1QA','APOE','SEPP1','RNASE1'),
  Neu=c('CSF3R','CXCL8','SOD2','NAMPT'),
  mDC=c('CD1C','HLA-DRA','HLA-DPB1','CST3','HLA-DPA'),
  pDC=c('LILRA4','PTGDS','SOX4','GZMB','IRF7'),
  HSC=c('CYTL1','GATA2')
)

Nervous = list(
  excitatory = c("SLC17A6"),
  inhibitory = c('GAD2'), 
  GABAergic = c('GAD2','GRIK1'), 
  dopaminergic_neurons = c('TH'),
  astrocytes = c("AQP4", "ADGRV1", "GPC5", "RYR3"),
  endothelial = c("CLDN5", "ABCB1", "EBF1"),
  excitatory = c("CAMK2A", "CBLN2", "LDB2"),
  inhibitory = c("GAD1", "LHFPL3", "PCDH15"), 
  microglia = c("C3", "LRMDA", "DOCK8"),
  oligodendrocytes = c("MBP", "PLP1", "ST18"),
  OPC=c('Tnr,Igsf21,Neu4,Gpr17'),
  Ependymal=c('Cfap126,Fam183b,Tmem212,pifo,Tekt1,Dnah12'),
  pericyte=c(  'DCN', 'LUM',  'GSN' ,'FGF7','MME', 'ACTA2','RGS5')
)

cluster_col = c4a("poly.wright25")
#"#1F78C8" "#FF0000" "#33A02C" "#6A33C2" "#FF7F00" "#565656" "#FFD700" "#A6CEE3" "#FB6496" "#B2DF8A" "#CAB2D6" "#FDBF6F" "#999999" "#EEE685" "#C8308C" "#FF83FA" "#C814FA" "#0000FF" "#36648B" "#00E2E5" "#00FF00" "#778B00" "#BEBE00" "#8B3B00" "#A52A3C"
celltype_col = cluster_col
sample_col = c4a("tableau.20")
#"#1F77B4" "#FF7F0E" "#2CA02C" "#D62728" "#9467BD" "#8C564B" "#E377C2" "#BCBD22" "#17BECF"
#"#4E79A7" "#A0CBE8" "#F28E2B" "#FFBE7D" "#59A14F" "#8CD17D" "#B6992D" "#F1CE63" "#499894" "#86BCB6" "#E15759" "#FF9D9A" "#79706E" "#D37295" "#FABFD2" "#B07AA1" "#D4A6C8" "#9D7660" "#D7B5A6"
group_col = c4a('light')
#"#77AADD" "#EE8866" "#EEDD88" "#FFAABB" "#99DDFF" "#44BB99" "#BBCC33" "#AAAA00"


set_directory <- function(dir_name) {
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  setwd(dir_name)
}

save_plots <- function(plot, file_name, width, height) {
  ggsave(filename = paste0(file_name, ".pdf"), plot = plot, device = 'pdf', dpi = 300, width = width, height = height)
  ggsave(filename = paste0(file_name, ".png"), plot = plot, device = 'png', dpi = 300, width = width, height = height)
}




pipe_01qcX = function(mydata){
  ### 01_QC ###-------------------------------------------------------------------------------------------------------------------------------
  set_directory("01_QC")

  mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "^MT-")
  mydata[["percent.ribo"]] <- PercentageFeatureSet(mydata, pattern = "^RP[SL]")
  mydata[["percent.hb"]] <- PercentageFeatureSet(mydata, pattern = "^HB[^(P)]")

  vlnplot=VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"), 
                  ncol = 5, cols = sample_col,pt.size = 0)

  ggsave(filename = "qc_vlnplot.pdf", plot = vlnplot,device = 'pdf',dpi = 300, width = 23, height = 8)
  ggsave(filename = "qc_vlnplot.png", plot = vlnplot,device = 'png',dpi = 300, width = 23, height = 8)

  print('01_QC process finished')
  setwd("..")

  return(mydata)
}


pipe_02clusteringX = function(mydata, hvg=2000, pca=15, res=0.1){
  ### 02_cluster ###--------------------------------------------------------------------------------------------------------------------------
  set_directory("02_clustering")
  
  mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = hvg, verbose = F)
  mydata <- ScaleData(mydata, features = rownames(mydata), verbose = F)
  mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata), verbose = F)
  mydata <- FindNeighbors(mydata, dims = 1:pca, verbose = F)
  mydata <- FindClusters(mydata, resolution = res, verbose = F)
  mydata <- RunUMAP(mydata, dims = 1:pca, verbose = F)
  gc(verbose = F)
  
  umap_plot <- advanced_dimplot(mydata, color = cluster_col)
  save_plots(umap_plot, "clusters_umap", 9, 8)
  umap_plot <- advanced_dimplot(mydata, group.by = 'orig.ident', color = cluster_col)
  save_plots(umap_plot, "samples_umap", 9, 8)

  cell_clustering=as.data.frame(cbind(colnames(mydata),mydata@meta.data$seurat_clusters))
  colnames(cell_clustering) = c('barcodes','seurat_clusters')
  write.csv(cell_clustering,file = "cell_clustering.csv")

  cell_umap = as.data.frame(mydata@reductions$umap@cell.embeddings)
  write.csv(cell_umap,file = "cell_umap.csv")

  # pheatmap
  av=AverageExpression(mydata, group.by = "seurat_clusters", assays = "RNA", verbose = F)
  av=av[[1]]
  cg=names(tail(sort(apply(av,1,sd)),1000))
  cluster_cor = round(cor(as.matrix(av[cg,]), method = "pearson"),3)
  cluster_corp = round(cor_pmat(as.matrix(av[cg,]), method = "pearson"),3)
  pdf(file = "clusters_cor_heatmap.pdf", width = 9, height = 8)
  corrplot(cluster_cor, method = "square",  addgrid.col = "darkgray", col.lim = c(0,1),
                  #p.mat = cluster_corp, insig = "label_sig", sig.level = c(.05),pch.cex = 1,pch = 1,
           order="hclust", addrect = 4, rect.col = "black", rect.lwd = 5,
           cl.pos = "b", tl.col = "indianred4", tl.cex = 1.5, cl.cex = 1.5,
           addCoef.col = "white", number.digits = 2, number.cex = 0.75,
           col = colorRampPalette(c("khaki4","palegreen","midnightblue","white","darkred"))(100))
  dev.off()
  png(file = "clusters_cor_heatmap.png", width = 900, height = 800)
  corrplot(cluster_cor, method = "square",  addgrid.col = "darkgray", col.lim = c(0,1),
           #p.mat = cluster_corp, insig = "label_sig", sig.level = c(.05),pch.cex = 1,pch = 1,
           order="hclust", addrect = 4, rect.col = "black", rect.lwd = 5,
           cl.pos = "b", tl.col = "indianred4", tl.cex = 1.5, cl.cex = 1.5,
           addCoef.col = "white", number.digits = 2, number.cex = 0.75,
           col = colorRampPalette(c("khaki4","palegreen","midnightblue","white","darkred"))(100))
  dev.off()
  #ggsave(filename = "clusters_cor_heatmap.pdf",plot = replayPlot(corp), device = 'pdf',dpi = 300, width = 9, height = 8)
  #ggsave(filename = "clusters_cor_heatmap.png",plot = replayPlot(corp), device = 'png',dpi = 300, width = 9, height = 8)

  print('02_clustering process finished')
  setwd("..")

  return(mydata)
}

pipe_03clustering_qcfilterX = function(mydata, mt=TRUE, mt_cutoff = 20, contam = TRUE,contam_cutoff = 0.25, singlet=TRUE, hvg=2000, pca=15, res=0.1){
  ### 03_clustering_qcfilter ###--------------------------------------------------------------------------------------------------------------------------
  set_directory("03_clustering_qcfilter")
  
  meta = mydata@meta.data
  meta$cellQC = 'Pass'
  meta$idents <- 1
  if(isTRUE(mt)){
    meta[meta$percent.mt>mt_cutoff,]$cellQC = 'Fail-MT'
  }
  if(isTRUE(contam)){
    meta[meta$Contam_scores>contam_cutoff,]$cellQC = 'Fail-Contam'
  }
  if(isTRUE(singlet)){
    meta[meta$doublet_scrublet=='TRUE',]$cellQC = 'Fail-Doublet'
  }
  write.csv(meta,file = "cellQC_stats.csv")
  stats_plot <- ggbarstats(data = meta, 
                   x = cellQC, y = orig.ident,
                   package = 'ggsci',
                   palette = 'category20c_d3',
                   colors = c('#1BA624','#7B539E','#EB852F','#B39B7C'),
                   ggstatsplot.layer = FALSE,
                   results.subtitle = FALSE,
                   bf.message = FALSE,
                   proportion.test = FALSE,
                   label.args = list(size = 5, 
                                     fill = 'black', 
                                     label.size = 0,
                                     alpha = 0,
                                     fontface = 'bold'),
                   perc.k = 2,
                   title = '',
                   xlab = '',
                   legend.title = 'cellQC',
                   ggtheme = ggpubr::theme_pubclean()) +
    theme(axis.ticks.x = element_line(color = 'black', lineend = 'round'),
      axis.ticks.y = element_line(color = 'black', lineend = 'round'),
      axis.title.x=element_text(size=18, angle = 90),
      legend.position = 'right',
      axis.text.x = element_text(size = 18, color = 'black',angle = 0),
      axis.text.y = element_text(size = 18, color = 'black'),
      legend.text = element_text(size = 18, color = 'black'),
      legend.title = element_text(size = 18, color = 'black'))+
    scale_fill_manual(values=c('#1BA624','#7B539E','#EB852F','#B39B7C'))
  save_plots(stats_plot, "samples_cellQC", 11, 8)

  if(isTRUE(mt)){
    cellnum1 = ncol(mydata)
    mydata <- subset(mydata, subset = percent.mt < mt_cutoff)
    cellnum2 = ncol(mydata)
    cellnum = cellnum1-cellnum2
    print(paste('During the QC process of mitochondria, ',cellnum,' cells were screened out.', sep = ''))
  }else{
    cellnum = 0
    print(paste('During the QC process of mitochondria, ',cellnum,' cells were screened out.', sep = ''))
  }

  if(isTRUE(contam)){
    cellnum1 = ncol(mydata)
    mydata <- subset(mydata, subset = Contam_scores < contam_cutoff)
    cellnum2 = ncol(mydata)
    cellnum = cellnum1-cellnum2
    print(paste('During the QC process of contamination, ',cellnum,' cells were screened out.', sep = ''))
  }else{
    cellnum = 0
    print(paste('During the QC process of contamination, ',cellnum,' cells were screened out.', sep = ''))
  }

  if(isTRUE(singlet)){
    cellnum1 = ncol(mydata)
    mydata <- subset(mydata, subset = doublet_scrublet == 'FALSE')
    cellnum2 = ncol(mydata)
    cellnum = cellnum1-cellnum2
    print(paste('During the QC process of Doublet, ',cellnum,' cells were screened out.', sep = ''))
  }else{
    cellnum = 0
    print(paste('During the QC process of Doublet, ',cellnum,' cells were screened out.', sep = ''))
  }

  mydata <- pipe_02clusteringX(mydata)

  print('03_clustering_qcfilter process finished')
  setwd("..")

  return(mydata)
  gc(verbose = F)
}

pipe_doubletfinder =function(mydata){
  ###### DoubletFinder ######
  mpK <- 10
  annotations <- mydata@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.04*ncol(mydata@assays$RNA@data))
  seurat_filterDouble <- doubletFinder_v3(mydata, PCs = 1:15, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  cn<-colnames(seurat_filterDouble@meta.data)
  cn[length(cn)] <- "Doublet"
  colnames(seurat_filterDouble@meta.data)<-cn
  mydata$doublet_finder=seurat_filterDouble$Doublet
  
  cell_clustering=as.data.frame(cbind(colnames(mydata),mydata@meta.data$doublet_finder))
  colnames(cell_clustering) = c('barcodes','doublet_finder')
  write.csv(cell_clustering,file = "cell_doubletfinder.csv")  
  
  Doubletplot <- advanced_dimplot(mydata, group.by='doublet_finder')
  save_plots(Doubletplot, "doubletfinder_umap", 9, 8)
  
  return(mydata)
  gc(verbose = F)
}

pipe_scrublet =function(mydata){
  ###### Scrublet ######
  scrublet <- import('scrublet')
  scr = scrublet$Scrublet(t(Matrix(mydata@assays$RNA@counts, sparse = TRUE)))
  doublet_scores = scr$scrub_doublets()
  mydata$doublet_scrublet_scores <- doublet_scores[[1]]
  mydata$doublet_scrublet <- doublet_scores[[2]]
  
  cell_clustering=as.data.frame(cbind(colnames(mydata),mydata@meta.data$doublet_scrublet_scores,mydata@meta.data$doublet_scrublet))
  colnames(cell_clustering) = c('barcodes','doublet_scrublet_scores','doublet_scrublet')
  write.csv(cell_clustering,file = "cell_scrublet.csv")
  
  Contamplot <- FeaturePlot(mydata,features = 'doublet_scrublet_scores', pt.size = .1, reduction = 'umap')
  save_plots(Contamplot, "scrublet_featureplot", 9, 8)
  Doubletplot <- advanced_dimplot(mydata, group.by='doublet_scrublet')
  save_plots(Doubletplot, "scrublet_umap", 9, 8)
  
  return(mydata)
  gc(verbose = F)
}

pipe_contam = function(mydata){
  ###### Contamination ######
  counts <- mydata@assays$RNA@counts
  decontX_results <- decontX(counts)
  mydata$Contam_scores = decontX_results$contamination
  mydata$Contam = ifelse(decontX_results$contamination<0.25,'not_contaminated','contaminated')
  
  cell_clustering=as.data.frame(cbind(colnames(mydata),mydata@meta.data$Contam_scores,mydata@meta.data$Contam))
  colnames(cell_clustering) = c('barcodes','Contam_scores','Contam')
  write.csv(cell_clustering,file = "cell_Contamination.csv")
  
  Contamplot <- FeaturePlot(mydata,features = 'Contam_scores', pt.size = .1, reduction = 'umap')
  save_plots(Contamplot, "contamination_featureplot", 9, 8)
  Contamplot <- advanced_dimplot(mydata, group.by='Contam')
  save_plots(Contamplot, "contamination_umap", 9, 8)
  
  
  return(mydata)
  gc(verbose = F)
}

pipe_cycling = function(mydata){
  ###### cellcycling ######
  mydata<- CellCycleScoring(mydata, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
  
  cell_clustering=as.data.frame(cbind(colnames(mydata),mydata@meta.data$Phase))
  colnames(cell_clustering) = c('barcodes','Phase')
  write.csv(cell_clustering,file = "cell_Phase.csv")
  
  cyclingplot <- advanced_dimplot(mydata, group.by='Phase',color =  c('#9AC9DBFF','#F8AC8CFF','#2878B5FF'))
  save_plots(cyclingplot, "cellcycling_umap", 9, 8)
  
  return(mydata)
  gc(verbose = F)
}

pipe_03cellstatusX = function(mydata, contam=TRUE, doublet=TRUE, cycing=TRUE){
  ### 03_check abnormal cell ###------------------------------------------------------------------------------------------------------------
  set_directory("03_cellstatus")
  
  ###### 03_sub contaminated RNA ######
  if(isTRUE(contam)){
    print('RUNNING contamination detection')
    mydata = pipe_contam(mydata)
  }else{
    print('PASS contamination detection')
  }

  ###### 03_sub Scrublet ######
  if(isTRUE(doublet)){
    print('RUNNING doublet detection')
    mydata = pipe_scrublet(mydata)
  }else{
    print('PASS doublet detection')
  }

  ###### 03_sub cell cycling ######
  if(isTRUE(doublet)){
    print('RUNNING cellcycing detection')
    mydata = pipe_cycling(mydata)
  }else{
    print('PASS cellcycing detection')
  }
  

  print('03_cellstatus process finished')
  setwd("..")

  return(mydata)
  gc(verbose = F)
}


pipe_04DEGsX = function(mydata, by = 'seurat_clusters', minpct=0.5, logfc=1, minp=0.01){
  ### 04_DEG ###--------------------------------------------------------------------------------------------------------------------------
  set_directory("04_DEGs")

  ###### 04_sub findmarkers ######
  mydata@active.ident <- as.factor(mydata$seurat_clusters)
  pbmc.markers <- FindAllMarkers(mydata, only.pos = TRUE, min.pct = minpct, logfc.threshold = logfc, return.thresh = minp, verbose = FALSE)
  Top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  degs_heatmap<-DoHeatmap(mydata, features = Top10$gene, angle = -50, hjust=0.8, raster = FALSE,assay = 'RNA',
                          label = T ,draw.lines = TRUE,size = 8,
                          group.colors=c(cluster_col),
                          group.bar.height = 0.03) +
    theme(axis.text.y = element_text(size = 9))+
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =5, name = "RdBu")))
  save_plots(degs_heatmap, "TopDEGs_heatmap", 11, 12)
  write.csv(pbmc.markers, file='clusters_DEGs.csv')

  av=AverageExpression(mydata, group.by = "seurat_clusters", assays = "RNA", verbose = F, return.seurat = T)
  av$seurat_clusters = levels(factor(mydata$seurat_clusters))
  degs_heatmap<-DoHeatmap(av, group.by = "seurat_clusters",features = Top10$gene,assay = 'RNA',
                          angle = -50, hjust=0.8, raster = FALSE,
                          label = T ,draw.lines = F,size = 8,
                          group.colors=c(cluster_col),
                          group.bar.height = 0.03) +
    theme(axis.text.y = element_text(size = 9))+
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =5, name = "RdBu")))
  save_plots(degs_heatmap, "TopDEGs_heatmap_avg", 11, 12)


  ###### 04_sub featureplot ######
  Top20 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  for (i in levels(Top20$cluster)){
    name <- paste("cluster",i,sep="")
    dir.create(name)
    sub <- subset(Top20, cluster %in% i)
    for(j in sub$gene){
      p<-advanced_featureplot(mydata,genelist=j)
      save_plots(p, paste(name,"/",j,sep=""), 9, 8)
    }
  }

  print('04_DEGs process finished')
  setwd("..")

  return(mydata)
  gc(verbose = F)
}


pipe_05pathwayX = function(mydata, by = 'seurat_clusters', minpct=0.5, logfc=1, minp=0.01){
  ### 05_GOKEGG ###--------------------------------------------------------------------------------------------------------------------------
  set_directory("05_GOKEGG")

  if (file.exists("../04_DEGs/clusters_DEGs.csv") == T){
    group_g <- read.csv('../04_DEGs/clusters_DEGs.csv', header = T, row.names = 1)
  }else {
    mydata@active.ident <- as.factor(mydata$seurat_clusters)
    group_g <- FindAllMarkers(mydata, only.pos = TRUE, min.pct = minpct, logfc.threshold = logfc, return.thresh = minp, verbose = FALSE)
  }
  tmp <- bitr(group_g$gene, fromType="SYMBOL",
              toType = c("ENSEMBL", "ENTREZID"),
              OrgDb="org.Hs.eg.db")
  de_gene_clusters=merge(tmp,group_g,by.x='SYMBOL',by.y='gene')

  ###### 05_sub GO enrichment ######
  formula_res <- compareCluster(
    ENTREZID~cluster,
    data=de_gene_clusters,
    fun="enrichGO",
    OrgDb="org.Hs.eg.db",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)

  t<-formula_res@compareClusterResult
  write.csv(t,file="formula_res_go.csv")

  lineage1_ego <- clusterProfiler::simplify(
    formula_res,
    cutoff=0.1,
    by="p.adjust",
    select_fun=min
  )

  formula_res %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = p.adjust) -> formula_res2

  go_plot = enrichplot::dotplot(formula_res2, showCategory=3)
  save_plots(go_plot, "enrichplot_go", 12, 11)


  ###### 05_sub KEGG enrichment ######
  options(clusterProfiler.download.method = "wget")
  formula_res <- compareCluster(
    ENTREZID~cluster,
    data=de_gene_clusters,
    fun="enrichKEGG",
    organism="hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  formula_res %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = p.adjust) -> formula_res2

  t<-formula_res@compareClusterResult
  write.csv(t,file="formula_res_kegg.csv")

  kegg_plot = enrichplot::dotplot(formula_res2, showCategory=3)
  save_plots(kegg_plot, "enrichplot_kegg", 12, 11)

  ###### 05_sub forloop enrichment ######
  # options(clusterProfiler.download.method = "wget")
  # for (i in levels(mydata$seurat_clusters)){
  #   dir_name=paste("cluster",i,sep="")
  #   dir.create(dir_name)
  #   markers <- subset(group_g, cluster==i)$gene
  #   genelist = bitr(markers, fromType="SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
  #   ego <- enrichGO(genelist$ENTREZID,OrgDb = "org.Hs.eg.db", ont = "BP",readable = T,pvalueCutoff = 0.5, qvalueCutoff = 1)
  #   df_ego <- summary(ego)
  #   if(dim(df_ego)[1]!=0){
  #     p1=barplot(ego,title="enrichGO")
  #     ggsave(filename = paste(dir_name,"/DEG_GO_cluster",i,".pdf",sep=""), plot = p1,device = 'pdf',dpi = 300, width = 12, height = 11)
  #     ggsave(filename = paste(dir_name,"/DEG_GO_cluster",i,".png",sep=""), plot = p1,device = 'png',dpi = 300, width = 12, height = 11)
  #     t1=ego@result
  #     write.csv(t1,file=paste(dir_name,"/DEG_GO_cluster",i,".csv",sep=""))
  #   }
  #   kk <- enrichKEGG(gene = genelist$ENTREZID, organism ="hsa", pvalueCutoff = 0.5, qvalueCutoff = 1, minGSSize = 1, use_internal_data =FALSE)
  #   df_kk <- summary(kk)
  #   if(dim(df_kk)[1]!=0){
  #     p2=barplot(kk,title="KEGG")
  #     ggsave(filename = paste(dir_name,"/DEG_KEGG_cluster",i,".pdf",sep=""), plot = p2,device = 'pdf',dpi = 300, width = 12, height = 11)
  #     ggsave(filename = paste(dir_name,"/DEG_KEGG_cluster",i,".png",sep=""), plot = p2,device = 'png',dpi = 300, width = 12, height = 11)
  #     t2=kk@result
  #     write.csv(t2,file=paste(dir_name,"/DEG_KEGG_cluster",i,".csv",sep=""))
  #   }
  # }
  print('05_GOKEGG process finished')
  setwd("..")

  return(mydata)
  gc(verbose = F)

}



pipe_06annoX = function(mydata){
  ### 06_Annotation ###----------------------------------------------------------------------------------------------------------------------
  set_directory("06_Anno")

  ###### 06_sub singleR ######
  refdata <- celldex::HumanPrimaryCellAtlasData()
  testdata <- GetAssayData(mydata, slot="data")
  clusters <- mydata$seurat_clusters
  cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main,
                      method = "cluster", clusters = clusters,
                      assay.type.test = "logcounts", assay.type.ref = "logcounts")
  celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

  mydata$SingleR = "NA"
  for(i in 1:nrow(celltype)){
    mydata@meta.data[which(mydata$seurat_clusters == celltype$ClusterID[i]),'SingleR'] <- celltype$celltype[i]
  }

  anno_plot<-advanced_dimplot(mydata, group.by='SingleR')
  save_plots(anno_plot, "SingleR_anno_umap", 9, 8)

  ###### 06_sub cellMarkers dotplot ######

  anno_plot = DotPlot(mydata, assay = "RNA", features = celltype_major, group.by = 'seurat_clusters') +
    theme(axis.text.x = element_text(angle = 45,  vjust = 0.9, hjust=0.9)) +
    scale_colour_gradient2(low = "steelblue", mid = "lightgrey", high = "#DD0000") +
    RotatedAxis()
  save_plots(anno_plot, "cellMarkers_dotplot", 9, 8)


  ###### 06_sub cellMarkers featureplot ######
  anno_plot = advanced_featureplot(mydata)
  save_plots(anno_plot, "cellMarkers_featureplot", 9, 23)

  ###### 06_sub cellMarkers violinplot ######
  anno_plot = advanced_violinplot(mydata, color = cluster_col)
  save_plots(anno_plot, "cellMarkers_vlnplot", 23, 9)


  print('06_Anno process finished')
  setwd("..")

  return(mydata)
  gc(verbose = F)
}


pipe_07saveX = function(mydata, name){
  ### 07_out ###----------------------------------------------------------------------------------------------------------------------
  set_directory("07_out")

  saveRDS(mydata,file=paste(name,'.rds', sep = ''))
  meta<-mydata@meta.data
  write.csv(meta,file=paste(name,'_metadata.csv', sep = ''))

  print('07_save process finished')
  setwd("..")
  gc()

  return(mydata)
}

advanced_featureplot = function(mydata, genelist=NULL, reduction='umap',ncol=5){
  if(is.null(genelist)) {
    genelist = unique(unlist(celltype_major))
    genelist = intersect(genelist, rownames(mydata))
  }

  if(length(genelist) == 1){
    FeaturePlot(mydata, features = genelist,pt.size = .1, reduction = 'umap')+
      scale_color_gradient(low = 'lightgrey',high = '#DD0000',name = 'Expr')+
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
      theme(
        legend.position = "right",
        legend.title = element_blank(), #去掉legend.title
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm') )
  } else {
    GeneExp <- FetchData(mydata,vars = genelist)
    pc <- Embeddings(mydata,reduction = reduction) %>% data.frame()
    colnames(pc) <- c('Dim_1','Dim_2')
    gbid <- cbind(pc,GeneExp)
    gbidlong <- melt(gbid,id.vars = c('Dim_1','Dim_2'),value.name = 'exp',variable.name = 'gene')

    ggplot(gbidlong,aes(x = Dim_1,y = Dim_2,color = exp)) +
      geom_point(size = 0,show.legend = T) +
      scale_color_gradient(low = 'lightgrey',high = '#DD0000',name = 'Expr') +
      theme_bw(base_size = 16) +
      theme(panel.grid = element_blank(),
            axis.ticks = element_blank(),
            aspect.ratio = 1,
            strip.background = element_rect(colour = NA,fill = NA),
            axis.text = element_blank(),
            plot.title = element_text(size=16,hjust = 0.5)) +
      facet_wrap(~gene,ncol = ncol)
  }

}

advanced_violinplot = function(mydata, genelist=NULL, group.by=NULL,color=NULL){
  if(is.null(genelist)){
    genelist = unique(unlist(celltype_major))
    genelist = intersect(genelist, rownames(mydata))
  }else{}

  dd = as.data.frame(t(as.matrix(mydata@assays$RNA@data[genelist,])))
  all(rownames(dd)==rownames(mydata@meta.data))

  if(is.null(group.by)){
    dd$group = mydata$seurat_clusters
  }else{
    dd$group = mydata@meta.data[,group.by]
  }

  dd = melt(dd)
  colnames(dd) = c('cluster','genes', 'expr')

  if(is.null(color)){
    mycol<- c4a('20')
  }else{
    mycol <- color
  }

  p = ggplot(data = dd,aes(x = expr, y = cluster, fill = cluster)) +
    geom_violin(scale = 'width',
                draw_quantiles= c(0.25, 0.5, 0.75),
                color= 'black',
                size= 0.45,
                alpha= 0.8) +
    facet_grid(cols = vars(genes), scales = 'free_x')+
    scale_fill_manual(values = mycol) + #填充色修改
    scale_x_continuous(breaks = seq(0, 8, by = 4)) + #x轴刻度显示调整
    theme_bw()+
    theme(
      panel.grid = element_blank(), #移除背景网格线
      axis.text.x = element_blank(), #x轴标签大小调整
      axis.text.y = element_text(size = 16), #y轴标签大小调整
      axis.title.x = element_text(size = 16), #x轴标题大小调整
      axis.title.y = element_blank(), #移除y轴标题
      axis.ticks.x = element_blank(),
      strip.background = element_blank(), #移除分面外围边框
      strip.text.x = element_text(size = 16, angle = 60), #分面文本标签倾斜60°
      legend.title = element_text(size = 16), #图例标题大小调整
      legend.text = element_text(size = 15) #图例标签大小调整
    ) +
    labs(x = 'Log Normalized Expression')

  return(p)
}


advanced_dimplot = function(mydata, group.by=NULL, reduction=NULL,color=NULL){
  if(is.null(group.by)){
    group.by = 'seurat_clusters'
  }else{}

  if(is.null(reduction)){
    reduction = 'umap'
  }else{}

  reduct = unlist(mydata@reductions[reduction][[1]])
  umap = reduct@cell.embeddings %>%
    as.data.frame() %>%
    cbind(cell_type = mydata@meta.data[,group.by])
  colnames(umap) = c('DIM_1','DIM_2','clusters')

  if(is.null(color)){
    mycol<- c4a("poly.wright25")
  }else{
    mycol <- color
  }

  ggplot(umap,aes(x= DIM_1 , y = DIM_2 ,color = clusters)) +
    geom_point(size = 0.1 , alpha =0.8 ) +
    scale_color_manual(values = mycol)+
    #scale_color_discrete_c4a_cat("carto.safe")+
    theme(panel.grid.major = element_blank(), #主网格线
          panel.grid.minor = element_blank(), #次网格线
          panel.border = element_blank(), #边框
          axis.title = element_blank(),  #轴标题
          axis.text = element_blank(), # 文本
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'white'), #背景色
          plot.background=element_rect(fill="white")) +
    theme(
      legend.title = element_blank(), #去掉legend.title
      legend.key=element_rect(fill='white'), #
      legend.text = element_text(size=20), #设置legend标签的大小
      legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
    guides(color = guide_legend(override.aes = list(size=5))) + #设置legend中 点的大小
    geom_segment(aes(x = min(DIM_1) , y = min(DIM_2) ,
                     xend = min(DIM_1) +3, yend = min(DIM_2) ),
                 colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+
    geom_segment(aes(x = min(DIM_1)  , y = min(DIM_2)  ,
                     xend = min(DIM_1) , yend = min(DIM_2) + 3),
                 colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
    annotate("text", x = min(umap$DIM_1) +1.5, y = min(umap$DIM_2) -1, label = "DIM_1",
             color="black",size = 6 ) +
    annotate("text", x = min(umap$DIM_1) -1, y = min(umap$DIM_2) + 1.5, label = "DIM_2",
             color="black",size = 6, angle=90) +
    theme(legend.position = "right",
          legend.text = element_text(size=23))
}

pipe_dataMerge = function(indir){
  indir = indir
  indirX = paste(indir, list.files(indir),sep = '/')
  names(indirX) = list.files(indir)
  mydata <- list()
  for(i in 1:length(list.files(indir))){
    counts <- Read10X(data.dir = indirX[i])
    mydata[[i]] <- CreateSeuratObject(counts, min.cells = 3, min.features = 200)
  }
  mydata <- merge(mydata[[1]], y=c(mydata[2:length(list.files(indir))]))
  saveRDS(mydata, 'seuratObject_raw.rds')
  return(mydata)
  gc(verbose = F)
}

pipe_harmony = function(mydata, by=NULL,hvg=2000,pca=15, res=0.1){
  ### 03_clustering_harmonyInteg ###--------------------------------------------------------------------------------------------------------------------------
  set_directory("03_clustering_harmonyInteg")
  if(is.null(by)){
    by = 'orig.ident'
  }else{
    by = by
  }
  
  mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = hvg, verbose = F)
  mydata <- ScaleData(mydata, features = rownames(mydata), verbose = F)
  mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata), verbose = F)
  mydata <- mydata %>%
    RunHarmony(by, plot_convergence = FALSE)
  mydata <- FindNeighbors(mydata, dims = 1:pca,reduction = "harmony", verbose = F)
  mydata <- FindClusters(mydata, resolution = res, verbose = F)
  mydata <- RunUMAP(mydata, dims = 1:pca,reduction = "harmony", verbose = F)
  mydata$clusters_harmonyInteg = mydata$seurat_clusters
  gc(verbose = F)
  mydata$clusters_harmonyInteg = mydata$seurat_clusters

  umap=advanced_dimplot(mydata,color = cluster_col)
  save_plots(umap, "clusters_umap", 9, 8)
  umap=advanced_dimplot(mydata,group.by = 'orig.ident', color = sample_col)
  save_plots(umap, "samples_umap", 9, 8)
  
  cell_clustering=as.data.frame(cbind(colnames(mydata),mydata@meta.data$seurat_clusters))
  colnames(cell_clustering) = c('barcodes','seurat_clusters')
  write.csv(cell_clustering,file = "cell_clustering.csv")
  
  cell_umap = as.data.frame(mydata@reductions$umap@cell.embeddings)
  write.csv(cell_umap,file = "cell_umap.csv")
  
  # pheatmap
  av=AverageExpression(mydata, group.by = "seurat_clusters", assays = "RNA", verbose = F)
  av=av[[1]]
  cg=names(tail(sort(apply(av,1,sd)),1000))
  cluster_cor = round(cor(as.matrix(av[cg,]), method = "pearson"),3)
  cluster_corp = round(cor_pmat(as.matrix(av[cg,]), method = "pearson"),3)
  pdf(file = "clusters_cor_heatmap.pdf", width = 9, height = 8)
  corrplot(cluster_cor, method = "square",  addgrid.col = "darkgray", col.lim = c(0,1),
           #p.mat = cluster_corp, insig = "label_sig", sig.level = c(.05),pch.cex = 1,pch = 1,
           order="hclust", addrect = 4, rect.col = "black", rect.lwd = 5,
           cl.pos = "b", tl.col = "indianred4", tl.cex = 1.5, cl.cex = 1.5,
           addCoef.col = "white", number.digits = 2, number.cex = 0.75,
           col = colorRampPalette(c("khaki4","palegreen","midnightblue","white","darkred"))(100))
  dev.off()
  png(file = "clusters_cor_heatmap.png", width = 900, height = 800)
  corrplot(cluster_cor, method = "square",  addgrid.col = "darkgray", col.lim = c(0,1),
           #p.mat = cluster_corp, insig = "label_sig", sig.level = c(.05),pch.cex = 1,pch = 1,
           order="hclust", addrect = 4, rect.col = "black", rect.lwd = 5,
           cl.pos = "b", tl.col = "indianred4", tl.cex = 1.5, cl.cex = 1.5,
           addCoef.col = "white", number.digits = 2, number.cex = 0.75,
           col = colorRampPalette(c("khaki4","palegreen","midnightblue","white","darkred"))(100))
  dev.off()
  #ggsave(filename = "clusters_cor_heatmap.pdf",plot = replayPlot(corp), device = 'pdf',dpi = 300, width = 9, height = 8)
  #ggsave(filename = "clusters_cor_heatmap.png",plot = replayPlot(corp), device = 'png',dpi = 300, width = 9, height = 8)
  
  print('03_clustering_harmonyInteg process finished')
  setwd("..")
  
  return(mydata)
  gc(verbose = F)
}

pipe_seuratInteg = function(mydata, by=NULL, hvg=2000, pca=15, res=0.1){
  ### 03_clustering_seuratInteg ###--------------------------------------------------------------------------------------------------------------------------
  set_directory("03_clustering_seuratInteg")
  if(is.null(by)){
    by = 'orig.ident'
  }else{
    by = by
  }
  DefaultAssay(mydata) <- "RNA"
  split_seurat <- SplitObject(mydata, split.by = by)
  split_seurat <- lapply(X = split_seurat, FUN = function(x) {
    x <- NormalizeData(x, verbose = F)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg, verbose = F)
  })
  features <- SelectIntegrationFeatures(object.list = split_seurat,verbose = F)
  seurat.anchors <- FindIntegrationAnchors(object.list = split_seurat, anchor.features = features,
                                           dims = 1:pca, verbose = F)
  mydata <- IntegrateData(anchorset = seurat.anchors,verbose = F)
  DefaultAssay(mydata) <- "integrated"
  mydata <- ScaleData(mydata, features = rownames(mydata), verbose = F)
  mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata), verbose = F)
  mydata <- FindNeighbors(mydata, dims = 1:pca, verbose = F)
  mydata <- FindClusters(mydata, resolution = res, verbose = F)
  mydata$clusters_seuratInteg = mydata$seurat_clusters
  mydata <- RunUMAP(mydata, dims = 1:pca, verbose = F)
  
  umap=advanced_dimplot(mydata,color = cluster_col)
  save_plots(umap, "clusters_umap", 9, 8)
  umap=advanced_dimplot(mydata,group.by = 'orig.ident', color = sample_col)
  save_plots(umap, "samples_umap", 9, 8)
  
  cell_clustering=as.data.frame(cbind(colnames(mydata),mydata@meta.data$seurat_clusters))
  colnames(cell_clustering) = c('barcodes','seurat_clusters')
  write.csv(cell_clustering,file = "cell_clustering.csv")
  
  cell_umap = as.data.frame(mydata@reductions$umap@cell.embeddings)
  write.csv(cell_umap,file = "cell_umap.csv")
  
  # pheatmap
  av=AverageExpression(mydata, group.by = "seurat_clusters", assays = "RNA", verbose = F)
  av=av[[1]]
  cg=names(tail(sort(apply(av,1,sd)),1000))
  cluster_cor = round(cor(as.matrix(av[cg,]), method = "pearson"),3)
  cluster_corp = round(cor_pmat(as.matrix(av[cg,]), method = "pearson"),3)
  pdf(file = "clusters_cor_heatmap.pdf", width = 9, height = 8)
  corrplot(cluster_cor, method = "square",  addgrid.col = "darkgray", col.lim = c(0,1),
           #p.mat = cluster_corp, insig = "label_sig", sig.level = c(.05),pch.cex = 1,pch = 1,
           order="hclust", addrect = 4, rect.col = "black", rect.lwd = 5,
           cl.pos = "b", tl.col = "indianred4", tl.cex = 1.5, cl.cex = 1.5,
           addCoef.col = "white", number.digits = 2, number.cex = 0.75,
           col = colorRampPalette(c("khaki4","palegreen","midnightblue","white","darkred"))(100))
  dev.off()
  png(file = "clusters_cor_heatmap.png", width = 900, height = 800)
  corrplot(cluster_cor, method = "square",  addgrid.col = "darkgray", col.lim = c(0,1),
           #p.mat = cluster_corp, insig = "label_sig", sig.level = c(.05),pch.cex = 1,pch = 1,
           order="hclust", addrect = 4, rect.col = "black", rect.lwd = 5,
           cl.pos = "b", tl.col = "indianred4", tl.cex = 1.5, cl.cex = 1.5,
           addCoef.col = "white", number.digits = 2, number.cex = 0.75,
           col = colorRampPalette(c("khaki4","palegreen","midnightblue","white","darkred"))(100))
  dev.off()
  #ggsave(filename = "clusters_cor_heatmap.pdf",plot = replayPlot(corp), device = 'pdf',dpi = 300, width = 9, height = 8)
  #ggsave(filename = "clusters_cor_heatmap.png",plot = replayPlot(corp), device = 'png',dpi = 300, width = 9, height = 8)
  
  print('03_clustering_seuratInteg process finished')
  setwd("..")
  
  DefaultAssay(mydata) <- "RNA"
  return(mydata)
  gc(verbose = F)
}


pipe_integLoop_optimized <- function(mydata, outdir, by, hvgs, pcas, reses, methods) {
  # 创建主目录
  main_dir <- file.path(outdir, "03_clustering_loop")
  if (!dir.exists(main_dir)) {
    dir.create(main_dir)
  }
  setwd(main_dir)
  
  # 循环不同的参数组合
  for (method in methods) {
    for (hvg in hvgs) {
      for (pca in pcas) {
        for (res in reses) {
          # 处理每个参数组合
          process_combination(mydata, method, hvg, pca, res, main_dir, by)
        }
      }
    }
  }
}

process_combination <- function(mydata, method, hvg, pca, res, main_dir, by) {
  # 创建特定参数组合的目录
  specific_dir <- file.path(main_dir, paste(method, hvg, pca, res, sep = '_'))
  if (!dir.exists(specific_dir)) {
    dir.create(specific_dir)
  }
  setwd(specific_dir)
  mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 10000,verbose=F)
  mydata <- ScaleData(mydata, features = rownames(mydata),verbose=F)
  # 根据方法处理数据
  if (method == 'seurat') {
    split_seurat <- SplitObject(mydata, split.by = by)
    split_seurat <- lapply(X = split_seurat, FUN = function(x) {
      x <- NormalizeData(x,verbose = F)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg, verbose = F)
    })
    features <- SelectIntegrationFeatures(object.list = split_seurat,verbose = F)
    seurat.anchors <- FindIntegrationAnchors(object.list = split_seurat, anchor.features = features,
                                             dims = 1:pca, verbose = F)
    gc(verbose = F)
    mydata <- IntegrateData(anchorset = seurat.anchors,verbose = F)
    gc(verbose = F)
    DefaultAssay(mydata) <- "integrated"
    mydata <- ScaleData(mydata, features = rownames(mydata), verbose = F)
    mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata), verbose = F)
    mydata <- FindNeighbors(mydata, dims = 1:pca, verbose = F)
    mydata <- FindClusters(mydata, resolution = res, verbose = F)
    mydata <- RunUMAP(mydata, dims = 1:pca, verbose = F)
    mydata$clusters_seuratInteg = mydata$seurat_clusters
    gc(verbose = F)
    
  } else if (method == 'harmony') {
    DefaultAssay(mydata) <- "RNA"
    mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = hvg, verbose = F)
    mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata), verbose = F)
    mydata <- mydata %>%
      RunHarmony(by, plot_convergence = FALSE)
    mydata <- FindNeighbors(mydata, dims = 1:pca,reduction = "harmony", verbose = F)
    mydata <- FindClusters(mydata, resolution = res, verbose = F)
    mydata <- RunUMAP(mydata, dims = 1:pca,reduction = "harmony", verbose = F)
    mydata$clusters_harmonyInteg = mydata$seurat_clusters
    gc(verbose = F)
  } else {
    DefaultAssay(mydata) <- "RNA"
    mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = hvg, verbose = F)
    mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata), verbose = F)
    mydata <- FindNeighbors(mydata, dims = 1:pca,reduction = "pca", verbose = F)
    mydata <- FindClusters(mydata, resolution = res, verbose = F)
    mydata <- RunUMAP(mydata, dims = 1:pca,reduction = "pca", verbose = F)
    mydata$clusters_harmonyInteg = mydata$seurat_clusters
    gc(verbose = F)
    
  }
  
  umap=advanced_dimplot(mydata,color = cluster_col)
  save_plots(umap, "clusters_umap", 9, 8)
  umap=advanced_dimplot(mydata,group.by = 'orig.ident', color = sample_col)
  save_plots(umap, "samples_umap", 9, 8)
  
  cell_clustering=as.data.frame(cbind(colnames(mydata),mydata@meta.data$seurat_clusters))
  colnames(cell_clustering) = c('barcodes','seurat_clusters')
  write.csv(cell_clustering,file = "cell_clustering.csv")
  cell_umap = as.data.frame(mydata@reductions$umap@cell.embeddings)
  write.csv(cell_umap,file = "cell_umap.csv")
  
  p1 <- advanced_dimplot(mydata, group.by = 'seurat_clusters',color = cluster_col)+ ggtitle("seurat_clusters") 
  p2 <- advanced_dimplot(mydata, group.by = 'orig.ident',color = sample_col)+ ggtitle("orig.ident")
  p <- p1|p2
  save_plots(p, "umap_overview", 15, 8)
  
  
  p = DotPlot(mydata, assay = "RNA", features = celltype_major, group.by = 'seurat_clusters') +
    theme(axis.text.x = element_text(angle = 45,  vjust = 0.9, hjust=0.9)) +
    scale_colour_gradient2(low = "steelblue", mid = "lightgrey", high = "#DD0000") +
    RotatedAxis()
  save_plots(p, "cellMarkers_dotplot", 11, 9)
  
  p = advanced_featureplot(mydata)
  save_plots(p, "cellMarkers_featureplot", 9, 23)
  
  p = advanced_violinplot(mydata, color = cluster_col)
  save_plots(p, "cellMarkers_vlnplot", 23, 9)
}

pipe_samplerunX = function(indir,outdir,samplename=NULL,deep_qc = TRUE,save=TRUE, midsave=FALSE, integ=TRUE, integ_method=c('seurat','harmony')){

  ###### 10X matrix read ######
  setwd(outdir)
  if(!is.null(samplename)){
    dir.create(samplename)
    setwd(samplename)
  }else{
    samplename =strsplit(indir, '/')[[1]][length(strsplit(indir, '/')[[1]])]
    dir.create(samplename)
    setwd(samplename)
  }
  if(isTRUE(length(list.files(indir))>length(sample_col))){
    sample_col = c(c4a("tableau.20"), c4a("tableau.classic10"), 
                   c4a("brewer.set1"),c4a("brewer.set2"),c4a("brewer.set3"),
                   c4a("poly.palette36"),c4a("dark24"))
  }
  mydata = pipe_dataMerge(indir = indir)
  gc(verbose = F)
  mydata = pipe_01qcX(mydata)
  gc(verbose = F)
  mydata = pipe_02clusteringX(mydata)
  gc(verbose = F)
  mydata = pipe_03cellstatusX(mydata)
  gc(verbose = F)
  if(deep_qc==TRUE){
    mydata = pipe_03clustering_qcfilterX(mydata)
    gc(verbose = F)
  }
  if(integ==TRUE){
    if('seurat' %in% integ_method){
      mydata = pipe_seuratInteg(mydata)
      gc(verbose = F)
    }
    if('harmony' %in% integ_method){
      mydata = pipe_harmony(mydata)
      gc(verbose = F)
    }
  }
  mydata = pipe_04DEGsX(mydata)
  gc(verbose = F)
  mydata = pipe_05pathwayX(mydata)
  gc(verbose = F)
  mydata = pipe_06annoX(mydata)
  gc(verbose = F)
  pipe_07saveX(mydata,name = samplename)

  return(mydata)
}

pipe_layout = function(dir, deep_qc=TRUE){
  library(figpatch)
  library(ggplot2)
  setwd(dir)
  p1 = fig("./01_QC/qc_vlnplot.png")
  p2 = fig("./02_clustering/clusters_umap.pdf")
  p3 = fig("./02_clustering/clusters_cor_heatmap.pdf")
  p4 = fig("./03_cellstatus/cellcycling_umap.pdf")
  p5 = fig("./03_cellstatus/contamination_umap.pdf")
  p6 = fig("./03_cellstatus/doublet_umap.pdf")

  cow_A <- cowplot::plot_grid(p1, nrow = 1)
  cow_BC <- cowplot::plot_grid(p2, p3,nrow = 1)
  cow_DEF = cowplot::plot_grid(p4,p5,p6, nrow = 1)
  fig1 = cowplot::plot_grid(cow_A,cow_BC,cow_DEF,nrow = 3)
  ggsave(filename = "fig1.pdf", plot = fig1,device = 'pdf',dpi = 600, width = 15, height = 23)
  ggsave(filename = "fig1.png", plot = fig1,device = 'png',dpi = 600, width = 15, height = 23)

  if(isTRUE(deep_qc)){
    p1 = fig("./03_clustering_qcfilter/clusters_umap.pdf")
    p2 = fig("./03_clustering_qcfilter/clusters_cor_heatmap.pdf")
  }else{
    p1 = fig("./02_clustering/clusters_umap.pdf")
    p2 = fig("./02_clustering/clusters_cor_heatmap.pdf")
  }
  p3 = fig("./04_DEGs/TopDEGs_heatmap.pdf")
  p4 = fig("./04_DEGs/TopDEGs_heatmap_avg.pdf")
  p5 = fig("./05_GOKEGG/enrichplot_go.pdf")
  p6 = fig("./05_GOKEGG/enrichplot_kegg.pdf")

  cow_AB <- cowplot::plot_grid(p1, p2, nrow = 1)
  cow_CD <- cowplot::plot_grid(p3,p4,nrow = 1)
  cow_EF = cowplot::plot_grid(p5,p6, nrow = 1)
  fig2 = cowplot::plot_grid(cow_AB,cow_CD,cow_EF, nrow = 3)
  ggsave(filename = "fig2.pdf", plot = fig2,device = 'pdf',dpi = 600, width = 15, height = 23)
  ggsave(filename = "fig2.png", plot = fig2,device = 'png',dpi = 600, width = 15, height = 23)

  p1 = fig("./06_Anno/SingleR_anno_umap.pdf")
  p2 = fig("./06_Anno/cellMarkers_dotplot.pdf")
  p3 = fig("./06_Anno/cellMarkers_vlnplot.pdf")
  p4 = fig("./06_Anno/cellMarkers_featureplot.pdf")

  cow_AB <- cowplot::plot_grid(p1,p2, nrow = 2)
  cow_ABC <- cowplot::plot_grid(cow_AB,p4,ncol = 2)
  fig3 = cowplot::plot_grid(cow_ABC,p3, ncol = 1)
  ggsave(filename = "fig3.pdf", plot = fig3,device = 'pdf',dpi = 600, width = 15, height = 23)
  ggsave(filename = "fig3.png", plot = fig3,device = 'png',dpi = 600, width = 15, height = 23)
}

load_fig <- function(path) {
  if (file.exists(path)) {
    return(figpatch::fig(path))
  } else {
    warning("File not found: ", path)
    return(NULL)
  }
}

save_fig_layout <- function(layout, filename, width, height, dpi = 600) {
  ggsave(filename = paste0(filename, ".pdf"), plot = layout, device = 'pdf', dpi = dpi, width = width, height = height)
  ggsave(filename = paste0(filename, ".png"), plot = layout, device = 'png', dpi = dpi, width = width, height = height)
}

create_and_save_layout <- function(dir, fig_paths, layout_matrix, filename, width = 15, height = 23, dpi = 600) {
  figs <- lapply(fig_paths, load_fig)
  layout <- cowplot::plot_grid(plotlist = figs, layout_matrix = layout_matrix)
  save_fig_layout(layout, filename, width, height, dpi)
}

pipe_layout <- function(dir, layouts_info) {
  for (layout_info in layouts_info) {
    create_and_save_layout(dir, layout_info$fig_paths, layout_info$layout_matrix, layout_info$filename, layout_info$width, layout_info$height, layout_info$dpi)
  }
}

dir <- "D:/tmp/lung/"
layouts_info <- list(
  list(
    fig_paths = c("./01_QC/qc_vlnplot.png", "./02_clustering/clusters_umap.pdf","./02_clustering/clusters_cor_heatmap.pdf",
                  "./03_cellstatus/cellcycling_umap.pdf","./03_cellstatus/contamination_umap.pdf","./03_cellstatus/doublet_umap.pdf"),
    layout_matrix = matrix(c(1, 2, 3,4,5,6), nrow = 3),
    filename = "layout1",
    width = 15,
    height = 10
  ),
  list(
    fig_paths = c("./03_clustering_qcfilter/clusters_umap.pdf", "./03_clustering_qcfilter/clusters_cor_heatmap.pdf", 
                  "./04_DEGs/TopDEGs_heatmap.pdf", "./04_DEGs/TopDEGs_heatmap_avg.pdf",
                  "./05_GOKEGG/enrichplot_go.pdf", "./05_GOKEGG/enrichplot_kegg.pdf"),
    layout_matrix = matrix(c(1, 2, 3), nrow = 1),
    filename = "layout2",
    width = 15,
    height = 10
  )
  
)



#library(autoscRNA)
#indir = 'D:/tmp/lung/'
#outdir = 'D:/tmp_out'
#start_time <- Sys.time() # 记录初始时间
#mydata = pipe_samplerunX(indir = indir,outdir = outdir,integ = TRUE, integ_method = c('seurat,'harmony'))
#set_directory(outdir)
#samplename =strsplit(indir, '/')[[1]][length(strsplit(indir, '/')[[1]])]
#setwd(samplename)
#mydata = pipe_dataMerge(indir = indir)
#mydata = pipe_01qcX(mydata)
#mydata = pipe_02clusteringX(mydata, hvg = 2000, pca=20, res = 0.1)
#mydata = pipe_03cellstatusX(mydata, contam = T, doublet = T, cycing = T)
#mydata = pipe_03clustering_qcfilterX(mydata)
#mydata = pipe_seuratInteg(mydata)
#mydata = pipe_harmony(mydata)
#pipe_integLoop_optimized(mydata, "D:/tmp_out/lung/", "orig.ident", c(2000, 3000), c(20,30), c(0.1,0.2), "harmony")
#mydata = pipe_04DEGsX(mydata)
#mydata = pipe_05pathwayX(mydata)
#mydata = pipe_06annoX(mydata)
#pipe_07saveX(mydata,name = samplename)
#end_time <- Sys.time() # 记录终止时间
#end_time - start_time # 计算时间差
#pipe_layout('D:/tmp_out/sample_hg17/')



