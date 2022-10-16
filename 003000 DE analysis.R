###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




##################################################
#
#
#  clustering 
#  differential expression analysis
#  GSEA
#
#
##################################################




#  Vermeulen et al. The Lancet Oncology (2009) list of prognostic gene set
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/draw_highlights.R')


#  [run once] create list of prognostic genes by manually copying and pasting the Table from the SI
#             get the corresponding circRNA/mRNA expression in the neuroblastoma tumors
#{{{
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/load_circ_gene_expression.R')


#  load reference to resolve gene_ids
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_name', 'gene_type')]


verm<-data.frame(gene_name=c('AHCY', 'AKR1C1', 'ARHGEF7', 'BIRC5', 'CADM1', 'CAMTA1', 'CAMTA2', 'CD44', 'CDCA5', 'CDKN3', 'CHD5', 'CLSTN1', 'CPSG3', 'DDC', 'DPYSL3', 'ECEL1', 'ELAVL4', 'EPB41L3', 'EPHA5', 'EPN2', 'FYN', 'GNB1', 'HIVEP2', 'INPP1', 'MAP2K4', 'MAP7', 'MAPT', 'MCM2', 'MRPL3', 'MTSS1', 'MYCN', 'NHLH2', 'NME1', 'NRCAM', 'NTRK1', 'ODC1', 'PAICS', 'PDE4DIP', 'PIK3R1', 'PLAGL1', 'PLAT', 'PMP22', 'PRAME', 'PRDM2', 'PRKACB', 'PRKCZ', 'PTN', 'PTPRF', 'PTPRH', 'PTPRN2', 'QPCT', 'SCG2', 'SLC25A5', 'SLC6A8', 'SNAPC1', 'TNFRSF25', 'TYMS', 'ULK2', 'WSB1'), risk_group_up=c('HR', 'LR', 'LR', 'HR', 'LR', 'LR', 'LR', 'LR', 'HR', 'HR', 'LR', 'LR', 'HR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'HR', 'HR', 'LR', 'HR', 'HR', 'HR', 'LR', 'LR', 'HR', 'HR', 'LR', 'LR', 'LR', 'LR', 'LR', 'HR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'HR', 'HR', 'HR', 'LR', 'HR', 'LR', 'LR'))


#  CPSG3 (NM_004386) is actually NCAN
setdiff(verm$gene_name, hsa$gene_name)
verm$gene_name[ verm$gene_name %in% 'CPSG3' ]<-'NCAN'


#  add gene_ids
verm$gene_id<-hsa$gene_id[ match(verm$gene_name, hsa$gene_name) ]
verm<-verm[, c('gene_id', 'gene_name', 'risk_group_up')]


#  load unified set of circRNAs
#  limit genes to those of the Vermeulen geneset, the gene list is defined irrespective of whether a circRNA exists in the unified group of circRNAs
#  load mRNA expression in TPMs and compute circRNA CPMs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
GENES<-DataFrame( verm[, c('gene_id', 'gene_name')] )
CIRCS<-CIRCS[ CIRCS$gene_id %in% verm$gene_id ]
length(CIRCS)                  #  45 circRNAs
length(unique(CIRCS$gene_id))  #  18 genes
l<-load_circ_gene_expression(CIRCS, GENES, META=nb.meta, nb.tumors.only=T, nb.only=T, gene.tpm=T)
nb.meta<-l$meta
nb.meta<-nb.meta[, risk_group:=factor(risk_group, levels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))][ order(risk_group), ][, risk_group:=as.character(risk_group)]
nb.circs<-l$circs
nb.genes<-l$genes
stopifnot( length(setdiff( GENES$gene_id, rownames(l$counts) ))==0 )
nb.counts<-l$counts[ GENES$gene_id, ]   #  (genes x samples) matrix order (do not convert gene_id to gene_name because there are duplicates!)
rm(list=c('l','CIRCS.all','GENES',ls(pattern='vt.|hb.')))


#  save
save(verm, nb.meta, nb.circs, nb.genes, nb.counts, CIRCS, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/Vermeulen et al_prognostic_geneset.RData')

#}}}


#  load the gene set and the circRNA/mRNA expression
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/Vermeulen et al_prognostic_geneset.RData')


#  [genes, log10-transformed TPMs, all risk groups] boxplot of LR-up and HR-up genes
#{{{

#  split by risk group
nb.<-split(nb.genes[, ! colnames(nb.genes) %in% c('risk_group', 'cell_model', 'tissue')],  nb.genes$risk_group)
l.LR<-lapply(nb., function(x){ log10(1+unname(unlist(c(x[, colnames(x) %in% verm$gene_name[ verm$risk_group_up %in% 'LR' ]])))) })
l.HR<-lapply(nb., function(x){ log10(1+unname(unlist(c(x[, colnames(x) %in% verm$gene_name[ verm$risk_group_up %in% 'HR' ]])))) })
l.LR<-l.LR[ unique(nb.meta$risk_group) ]
l.HR<-l.HR[ unique(nb.meta$risk_group) ]
L<-append(l.LR, l.HR)


#  boxplot
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.0, 7.0, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
options(scipen=0)
YTICK<-pretty(c(0, max(ceiling(sapply(L, max, na.rm=T)))), 4)
plot(0:1, 0:1, xlim=c(0, length(L))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(L, col=unique(nb.meta$col), ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=unique(nb.meta$col), range=0, add=T)
mtext(expression(log[10](1+TPM)), side=2, line=2, padj=-0.5, las=0, cex=2.4)
mtext(text=c('LR-up', 'HR-up'), side=1, line=1, at=seq(median(seq_along(l.LR)), length(L), length(l.LR)), las=1, adj=0.5, cex=2.4, col='black')
mtext(text=c(paste0('(', table(verm$risk_group_up)['LR'], ')'), paste0('(', table(verm$risk_group_up)['HR'], ')')), side=1, line=3, at=seq(median(seq_along(l.LR)), length(L), length(l.LR)), las=1, adj=0.5, cex=2.4, col='black')
draw_highlights(L=length(L), STEP=length(l.LR), YMAX=max(YTICK))
legend('topleft', legend=names(l.LR), col=unique(nb.meta$col), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, xpd=T, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_Vermeulen et al_prognostic_geneset_gene_TPMs.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  clean up
dev.off()

#}}}


#  [genes] hierarchical clustering and PCA 
#{{{

#  round up gene counts of all genes to integers and variance-stabilize them
stopifnot( length(setdiff(colnames(nb.counts), nb.meta$bid))==0 )
dds<-DESeqDataSetFromMatrix(countData=as.matrix(ceiling(nb.counts[, nb.meta$bid])), colData=data.frame(Risk=factor(nb.meta$risk_group), row.names=nb.meta$bid), design=~Risk) 
nb.sf<-sizeFactors(estimateSizeFactors(dds, type='poscounts'))
nb.vs<-assay(varianceStabilizingTransformation(dds, fitType='local', blind=T))
rm(dds,nb.sf)


#  select the genes of interest and convert gene_id to gene_name
stopifnot( length(setdiff(verm$gene_id, rownames(nb.vs)))==0 )
nb.vs<-nb.vs[ verm$gene_id, ]
rownames(nb.vs)<-verm$gene_name[ match(rownames(nb.vs), verm$gene_id) ]


#  recycle
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
ex<-data.frame(Risk=factor(nb.meta[ match(colnames(nb.vs), bid), risk_group], exclude=F), row.names=colnames(nb.vs))
cl<-setNames(list(setNames( nb.meta[, unique(col)], nb.meta[, unique(risk_group)] )), colnames(ex))


#  [log10-transformed TPMs] clustering based on Spearman correlations of complete pairs
x<-log10(1+t(nb.genes[, ! colnames(nb.genes) %in% c('risk_group', 'cell_model', 'tissue')]))
stopifnot( length(setdiff(colnames(x), colnames(nb.vs)))==0 )
x<-x[, colnames(nb.vs)]
d<-cor(x, method='spearman', use='pairwise.complete.obs')
hc<-hclust(as.dist(1-d), method='ward.D2')
ph<-pheatmap(as.matrix(d), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(20)), border_color=NA, scale='none', 
        breaks=seq(0, 1, length.out=21),  #  reduced range
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_Vermeulen et al_prognostic_geneset_genes_tpm_cor.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x)


#  [variance-stabilized counts] clustering based on Euclidean distance 
d<-dist(t(nb.vs), method='euclidean')
hc<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_Vermeulen et al_prognostic_geneset_genes_vsc.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [variance-stabilized counts] PCA centered but not scaled 
p<-prcomp(t(nb.vs), center=T, scale.=F)    #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-nb.meta[ match( rownames(pca), bid ), risk_group ]
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
   #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl$Risk) + 
    scale_color_manual(name='risk group', values=cl$Risk) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=34, color='black'), axis.title.x=element_text(size=34), 
        axis.text.y=element_text(size=34, color='black'), axis.title.y=element_text(size=34), 
        legend.text=element_text(size=22), legend.title=element_text(size=22, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/plotPCA_Vermeulen et al_prognostic_geneset_genes_vsc.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}


#  [circRNAs, log10-transformed CPMs, all risk groups] boxplot of LR-up and HR-up genes
#{{{

#  split by risk group
x<-as.data.table(as.data.frame(mcols(nb.circs)[, c('bid', 'circ_name', 'cpm')]))
x<-dcast(x, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)
y<-as.data.frame(x[, -1])
rownames(y)<-x[, bid]
stopifnot( length(setdiff(rownames(y), rownames(nb.genes)))==0 )
x<-y[ rownames(nb.genes), ]
LRcircs<-colnames(x)[grep(paste0(verm$gene_name[ verm$risk_group_up %in% 'LR' ], collapse='|'), colnames(x))]
HRcircs<-colnames(x)[grep(paste0(verm$gene_name[ verm$risk_group_up %in% 'HR' ], collapse='|'), colnames(x))]
stopifnot( length(LRcircs)+ length(HRcircs)==ncol(x) )
nb.<-split(x,  nb.genes$risk_group)
l.LR<-lapply(nb., function(x){ log10(1+unname(unlist(c(x[, LRcircs])))) })
l.HR<-lapply(nb., function(x){ log10(1+unname(unlist(c(x[, HRcircs])))) })
l.LR<-l.LR[ unique(nb.meta$risk_group) ]
l.HR<-l.HR[ unique(nb.meta$risk_group) ]
L<-append(l.LR, l.HR)


#  boxplot
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.0, 10.0, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
options(scipen=0)
#YTICK<-pretty(c(0, max(ceiling(sapply(L, max, na.rm=T)))), 4)
YTICK<-pretty(c(0, 0.2), 4)
plot(0:1, 0:1, xlim=c(0, length(L))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(L, col=unique(nb.meta$col), ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=unique(nb.meta$col), range=0, add=T)
mtext(expression(log[10](1+CPM)), side=2, line=5, padj=-0.5, las=0, cex=2.4)
mtext(text=c('LR-up', 'HR-up'), side=1, line=1, at=seq(median(seq_along(l.LR)), length(L), length(l.LR)), las=1, adj=0.5, cex=2.4, col='black')
mtext(text=c(paste0('(', table(verm$risk_group_up)['LR'], ')'), paste0('(', table(verm$risk_group_up)['HR'], ')')), side=1, line=3, at=seq(median(seq_along(l.LR)), length(L), length(l.LR)), las=1, adj=0.5, cex=2.4, col='black')
draw_highlights(L=length(L), STEP=length(l.LR), YMAX=max(YTICK), YMIN=0.0)
legend('topleft', legend=names(l.LR), col=unique(nb.meta$col), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, xpd=T, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
#dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_Vermeulen et al_prognostic_geneset_circRNAs_CPMs.pdf', width=18, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_Vermeulen et al_prognostic_geneset_circRNAs_CPMs.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/Vermeulen et al_prognostic_geneset.RData



#  [run once] Prepare circRNAs and genes for the analysis. We compute:
#
#                 gene TPMs, raw counts and variance-stabilized log2-transformed counts
#                 circRNA raw counts and variance-stabilized log2-transformed counts based on size factors computed from gene counts
#
#                 PCA for genes and circRNAs is done THROUGHOUT TUMORS AND CELL LINES using centered but not scaled variance-stabilized 
#                 (and log2-transformed) counts
#
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(vsn)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(pheatmap)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/load_circ_gene_expression.R')


#  load reference 
#  discard chrM/chrY genes from the analysis
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene']
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]
seqlevels(hsa, pruning.mode='coarse')<-seqlevels(hsa)[ grep('chrM|chrY', seqlevels(hsa), invert=T) ]


#  remove failed samples and Pilot samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta<-meta[ !(failed) & !grepl('CBPilote', bid) ]


#  add sex metadata to the tumor samples
#load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/sex_determination.RData')
#meta<-sex[ meta, on='bid'][, -c(1:3)]
#rm(sex)
#
#  N.B. latest clinical annotation includes sex besides for CB3009 which we find it to be male
meta[, sex:=SEX][, SEX:=NULL]
meta[ PAT_ID_BERLIN %in% 'CB3009', sex:='M' ]


#  load our cohort of circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
rm(list=ls(pattern='\\.(genes|circs|meta|prefailed|all)|GENES'))


#  circRNA counts and gene counts (for all genes) + TPMs (for the selected genes)
#  make sure the correct samples are kept throughout
nb<-load_circ_gene_expression(CIRCS=CIRCS, GENES=mcols(hsa)[, c('gene_id', 'gene_name')], META=meta, nb.tumors.only=F, nb.only=T, gene.tpm=T)
circ<-nb$circs[ nb$circs$bid %in% meta$bid ]
gns.tpm<-t(nb$genes[ meta$bid, grep('risk_group|cell_model|tissue', colnames(nb$genes), invert=T)])
stopifnot( length(setdiff(hsa$gene_id, rownames(nb$counts)))==0 )
gns.cnt<-ceiling(nb$counts[hsa$gene_id, meta$bid])  #  remove chrM, chrY counts, the resulting matrix form is: (gene_id X sample)
rm(nb)


#  remove unexpressed genes 
gns.cnt<-gns.cnt[rowSums(gns.cnt)!=0,  ]
gns.tpm<-gns.tpm[rowSums(gns.tpm)!=0,  ]


#  compute size factors across all samples (tumors and cell lines)
#  compute variance-stabilized gene counts
dds<-DESeqDataSetFromMatrix(countData=gns.cnt, colData=data.frame(risk_group=factor(meta$risk_group), cell_model=factor(meta$cell_model), row.names=meta$bid), design=~1)
gns.sf<-sizeFactors(estimateSizeFactors(dds, type='poscounts'))
gns.vs<-assay(varianceStabilizingTransformation(dds, fitType='local'))
rm(dds)


#  PCA on variance-stabilized (and glog2-transformed) gene counts centered but not scaled
gns.pca<-prcomp(t(gns.vs), center=T, scale.=F)
gns.ve<-round(1000 * gns.pca$sdev^2/sum(gns.pca$sdev^2))/10 


#  CHECK: compute size factors across tumors
#         compute variance-stabilized gene counts for the tumors
#
#{{{
#m<-meta[ !(failed) & !is.na(risk_group) ]
#g<-gns.cnt[, m$bid]
#dds<-DESeqDataSetFromMatrix(countData=g, colData=data.frame(risk_group=factor(m$risk_group), row.names=m$bid), design=~1)
#gns.sf.tumors<-sizeFactors(estimateSizeFactors(dds, type='poscounts'))
#gns.vs.tumors<-assay(varianceStabilizingTransformation(dds, fitType='local'))
#rm(dds)


#  how different are the tumor sizefactors between the methods?
#x<-gns.sf[ names(gns.sf.tumors) ]
#y<-gns.sf.tumors
#summary(x-y)
#
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.07566 -0.01798 -0.00850 -0.00619  0.00501  0.04903 


#  PCA on variance-stabilized (and glog2-transformed) gene counts centered but not scaled
#gns.pca.tumors<-prcomp(t(gns.vs.tumors), center=T, scale.=F)
#gns.ve.tumors<-round(1000 * gns.pca.tumors$sdev^2/sum(gns.pca$sdev^2))/10 

#}}}


#  summarize circRNA counts at the gene level 
#  force size factors to be those based from gene counts
#  variance-stabilize
#  PCA on variance-stabilized (and glog2-transformed) counts centered by not scaled
crs.cnt<-data.table(data.frame(mcols(circ)[, c('gene_id', 'jc_count', 'bid')]))[, .(jc_count=sum(jc_count)), by=.(bid, gene_id)]
crs.cnt<-dcast(crs.cnt, bid ~ gene_id, value.var='jc_count', fun.aggregate=sum)
crs.cnt<-t(as.matrix(data.frame(crs.cnt[, -1], row.names=crs.cnt$bid, check.names=F)))[, meta$bid]
dds<-DESeqDataSetFromMatrix(countData=crs.cnt, colData=data.frame(risk_group=factor(meta$risk_group), cell_model=factor(meta$cell_model), row.names=meta$bid), design=~1)
sizeFactors(dds)<-gns.sf
crs.vs<-assay(varianceStabilizingTransformation(dds, fitType='local'))
crs.pca<-prcomp(t(crs.vs), center=T, scale.=F)
crs.ve<-round(1000 * crs.pca$sdev^2/sum(crs.pca$sdev^2))/10 
rm(dds,circ)


#  save all including the reference for convenience
save(list=c('hsa', ls(pattern='(gns|crs)\\.'), 'meta', 'CIRCS'), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData



#  [tumors + cell lines] clustering 
#                        DESeq2 analysis
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(vsn)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  functions
#{{{

gene_results<-function(res=NULL, annot=NULL){
    m<-res[ intersect(annot$gene_id, rownames(res)), c('baseMean', 'log2FoldChange')]
    m$gene_name<-annot$gene_name[ match( rownames(m), annot$gene_id ) ]
    m$col<-annot$col[ match( rownames(m), annot$gene_id ) ]
    return(m)
}


do_deseq2<-function(DDS, CT, Type='', condition='Risk', lfcThreshold=0.0, lfcShrinkType='apeglm', include.annot=F, names.trim, ...){
    #             DDS = DESeqDataSet object of samples belonging to two conditions
    #              CT = metadata data.frame of only samples belonging to the two conditions
    #            Type = 'genes' or 'circRNAs' to add to figure/data names when saved
    #       condition = CT metadata column to use that defines the two conditions so we can look up the corresponding colors
    #    lfcThreshold = lfcThreshold to use in DESeq2
    #   lfcShrinkType = which method to use to estimate shrunken MAP log2FC?
    #   include.annot = shall we annotate special genes in the MAplot?
    #      names.trim = regular expression to use to trim sample names when plotting, e.g. '-11-R01$', pass '' for no trimming
    #                   N.B. you NEED to define this if you pass down unnamed arguments with ...
    #             ... = list of plotting parameters for par() and for the y-axis mtext to pass down when doing the MA-plot


    #  save open graphics devices by start
    #DEV_START<-dev.list()


    #  variance stabilizing transformation including normalization by library size factors glog2-transformed back to counts
    VSC<-assay(varianceStabilizingTransformation(DDS, fitType='local', blind=T))

    
    #  remove genes that have identical expression throughout the samples
    VSC<-VSC[ apply(VSC, 1, function(x){ any(x!=x[1]) }), ]


    #  PCA on variance-stabilized (and glog2-transformed) counts centered but not scaled
    PCA<-prcomp(t(VSC), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
    VE<-round(1000 * PCA$sdev^2/sum(PCA$sdev^2))/10 
    stopifnot( all(rownames(PCA$x)==colnames(VSC)) )
    PCA$x<-cbind(PCA$x, CT)


    #  differential expression analysis
    #  use Cook's distance to flag outliers but do not replace their values
    DDS<-DESeq(DDS, fitType='local', minReplicatesForReplace=Inf, betaPrior=F)
    RES<-results(DDS, alpha=0.05, lfcThreshold=lfcThreshold, altHypothesis='greaterAbs', cooksCutoff=T)
    RES<-lfcShrink(dds=DDS, coef=tail(resultsNames(DDS), 1), res=RES, type=lfcShrinkType)  
    RES<-RES[ order(RES$padj, decreasing=F), ]
    CND<-paste0(rev(levels(colData(DDS)[, condition])), collapse=' vs ')
    x11()
    plotDispEsts(DDS)  #  just to see
    dev.off()


    #  add gene_names
    RES$gene_name<-hsa$gene_name[ match(rownames(RES), hsa$gene_id) ]


    # MA-plot
    x11(width=16, height=16, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
    YLIM<-pretty(range(RES$log2FoldChange, na.rm=T))
    YLIM<-c(YLIM[1], tail(YLIM, 1))
    XLIM<-pretty(range(RES$baseMean, na.rm=T), 5)
    XLIM<-c(0.1, tail(XLIM, 1))
    if(...length()>0){
        dots<-list(...)[[1]]  #  it is already a list
        par(dots$par)
        YLAB.LINE<-dots$ylab.line
    } else {
        par(mar=c(5.0, 7.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2)
        YLAB.LINE<-3
    }
    options(scipen=+20)
    plotMA(RES, xlab='', ylab='', main='', ylim=YLIM, xlim=XLIM, cex=1.2, colSig='red3')
    mtext(CND, side=3, line=0, padj=+1.5, cex=1.4)
    mtext(expression(log[2]('fold change')), side=2, line=YLAB.LINE, padj=-0.2, cex=2.4, las=3)
    mtext('Mean expression', side=1, line=4, padj=-0.3, cex=2.4, las=1)
    if (include.annot){
        r<-gene_results(RES, annot)
        points(r$baseMean , r$log2FoldChange, pch=21, lwd=6, col='black', bg=r$col, cex=1.4)
        legend('bottomright', legend=r$gene_name, bty='n', lty=0, lwd=0, pch=21, col='black', pt.bg=r$col, pt.cex=1.8, pt.lwd=4, cex=1.2, x.intersp=-0.4, y.intersp=0.4)
    }
    if(lfcThreshold>0){
        abline(h=c(-lfcThreshold, lfcThreshold), lty=1, lwd=4, col='cyan4')
    }
    #identify(RES$baseMean, RES$log2FoldChange, labels=RES$gene_name, cex=0.7, offset=0.2, xpd=T)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/plotMA_', Type, '_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    rm(YLIM,XLIM)
    dev.off()


    #  mean-sd plots to see if variance strongly depends on mean
    X11(width=12, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
    par(mar=c(5,4,0.1,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.4, cex.axis=1.4)
    meanSdPlot(VSC)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/meanSD_', Type, '_', sub(' vs ', '_', CND), '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    dev.off()


    #  heatmap based on Euclidean distances of variance stabilized (and glog2-transformed) normalized counts
    x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
    v<-VSC
    ex<-as.data.frame(colData(DDS)[, condition, drop=F])
    stopifnot( all.equal( rownames(ex), rownames(CT) ) )
    cl<-setNames(list(setNames( unique(CT$col), unique(CT[, condition, drop=T]) )), colnames(ex))
    if(names.trim!=''){
        colnames(v)<-sub(names.trim, '', colnames(v))
        rownames(ex)<-sub(names.trim, '', rownames(ex))
    }
    d<-dist(t(v), method='euclidean')
    hc<-hclust(d, method='ward.D2')
    ph<-pheatmap(as.matrix(d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
            cluster_rows=hc,
            cluster_cols=hc,
            #cutree_row=4, cutree_col=4,
            annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
            annotation_col=ex, annotation_row=ex, annotation_colors=cl,
            drop_levels=F, show_rownames=T, show_colnames=T, 
            display_numbers=T, number_format='%.1f', number_color='grey39',
            fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_', Type, '_vsc_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    rm(v,d,cl,ex,hc,ph)
    dev.off()


    #  heatmap based on Spearman correlations of raw counts
    x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
    a<-assay(DDS)
    ex<-as.data.frame(colData(DDS)[, condition, drop=F])
    stopifnot( all.equal( rownames(ex), rownames(CT) ) )
    cl<-setNames(list(setNames( unique(CT$col), unique(CT[, condition, drop=T]) )), colnames(ex))
    if(names.trim!=''){
        colnames(a)<-sub(names.trim, '', colnames(a))
        rownames(ex)<-sub(names.trim, '', rownames(ex))
    }
    d<-cor(a, method='spearman', use='pairwise.complete.obs')
    hc<-hclust(as.dist(1-d), method='ward.D2')
    ph<-pheatmap(as.matrix(d), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(20)), border_color=NA, scale='none', 
            breaks=seq(-1, 1, length.out=21),
            cluster_rows=hc,
            cluster_cols=hc,
            #cutree_row=4, cutree_col=4,
            annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
            annotation_col=ex, annotation_row=ex, annotation_colors=cl,
            drop_levels=F, show_rownames=T, show_colnames=T, 
            #display_numbers=T, number_format='%.1f', number_color='grey39',
            fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_', Type, '_cor_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    rm(a,d,cl,ex,hc,ph)
    dev.off()


    #  PCA based on variance-stabilized, glog2-transformed, centered but not scaled counts
    XLIM<-pretty(range(PCA$x[, 'PC1']), 5)
    XLIM<-c(XLIM[1], XLIM[length(XLIM)])
    YLIM<-pretty(range(PCA$x[, 'PC2']), 5)
    YLIM<-c(YLIM[1], YLIM[length(YLIM)])
    x11(width=12, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
    x<-PCA$x[, c('PC1', 'PC2', condition, 'col')]
    x[, condition]<-factor(x[, condition], levels=unique(x[, condition]))
    if(names.trim!=''){
        rownames(x)<-sub(names.trim, '', rownames(x))
    }
    print(
        ggplot(x, aes(PC1, PC2, color=get(condition))) + 
        geom_point(size=10) + 
        #geom_text_repel(aes(label=rownames(x)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
        scale_shape_manual(values=18) + 
        scale_fill_manual(name='sample', values=unique(x$col)) + 
        scale_color_manual(name='sample', values=unique(x$col)) + 
        theme(text = element_text(family='Arial'), axis.text.x=element_text(size=34), axis.title.x=element_text(size=34), 
            axis.title.y=element_text(size=34), axis.text.y=element_text(size=34), 
            legend.text=element_text(size=22), legend.title=element_text(size=22, face='plain'), aspect.ratio=1, panel.background = element_blank(),
            axis.line.x=element_line(size=0.2, colour = 'black', linetype='solid'), 
            axis.line.y=element_line(size=0.2, colour = 'black', linetype='solid'),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
        scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
        xlab(paste0('PC1: ', VE[1], '% variance')) + ylab(paste0('PC2: ', VE[2], '% variance')) 
    )
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/plotPCA_', Type, '_vsc_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    rm(XLIM,YLIM)
    dev.off()


    #  save DESeq results along with gene counts
    options(scipen=0)
    save(DDS,CND,RES,VSC,PCA,VE, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.RData'), compress=T)


    #  convert to table 
    x<-RES
    x$gene_id<-rownames(x)
    x<-x[, c('gene_id', 'gene_name', 'baseMean', 'log2FoldChange', 'pvalue', 'padj')]
    rownames(x)<-NULL
    x<-as.data.frame(x)
    write.table(x, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.tsv'), quote=F, sep='\t', row.names=F, col.names=T) 
    write.xlsx(x, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.xlsx'), row.names=F, col.names=T, sheetName=gsub(' ', '_', CND)) 
    rm(x)


   #readline("Hit ENTER to close all figures: ") 
   #for (n in setdiff(dev.list(), DEV_START)){ dev.off(n) }

   return(list(dds=DDS, res=RES, cnd=CND))
}

#}}}


#  load pre-prepared counts etc. for all samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')


#  [tumors, circRNAs + genes] clustering and PCA
#{{{

#  [tumors, circRNAs] heatmap based on Euclidean distance of variance stabilized (and glog2-transformed) counts 
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-meta[!is.na(risk_group), ]
x<-crs.vs[, m$bid]
colnames(x)<-sub('-11-R01', '', colnames(x))
d<-dist(t(x), method='euclidean')
hc<-hclust(d, method='ward.D2')
ex<-data.frame(Type=factor(m$risk_group, exclude=F), row.names=colnames(x))
cl<-setNames(list(setNames( unique(m$col), unique(m$risk_group) )), colnames(ex))
ph<-pheatmap(as.matrix(d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_circRNAs_vsc_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [tumors, genes] heatmap based on Euclidean distance of variance stabilized (and glog2-transformed) counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-meta[!is.na(risk_group), ]
x<-gns.vs[, m$bid]
colnames(x)<-sub('-11-R01', '', colnames(x))
d<-dist(t(x), method='euclidean')
hc<-hclust(d, method='ward.D2')
ex<-data.frame(Type=factor(m$risk_group, exclude=F), row.names=colnames(x))
cl<-setNames(list(setNames( unique(m$col), unique(m$risk_group) )), colnames(ex))
ph<-pheatmap(as.matrix(d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_genes_vsc_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [tumors, circRNAs] heatmaps based on Spearman correlations of raw counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-meta[!is.na(risk_group), ]
x<-crs.cnt[, m$bid]
colnames(x)<-sub('-11-R01', '', colnames(x))
d<-cor(x, method='spearman', use='pairwise.complete.obs')
hc<-hclust(as.dist(1-d), method='ward.D2')
ex<-data.frame(Type=factor(m$risk_group, exclude=F), row.names=colnames(x))
cl<-setNames(list(setNames( unique(m$col), unique(m$risk_group) )), colnames(ex))
ph<-pheatmap(as.matrix(d), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(20)), border_color=NA, scale='none', 
        breaks=seq(0, 1, length.out=21),
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_circRNAs_cor_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [tumors, genes] heatmaps based on Spearman correlations of raw counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-meta[!is.na(risk_group), ]
x<-gns.cnt[, m$bid]
colnames(x)<-sub('-11-R01', '', colnames(x))
d<-cor(x, method='spearman', use='pairwise.complete.obs')
hc<-hclust(as.dist(1-d), method='ward.D2')
ex<-data.frame(Type=factor(m$risk_group, exclude=F), row.names=colnames(x))
cl<-setNames(list(setNames( unique(m$col), unique(m$risk_group) )), colnames(ex))
ph<-pheatmap(as.matrix(d), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(20)), border_color=NA, scale='none', 
        breaks=seq(0.7, 1, length.out=21),  #  reduce range since samples are highly correlated
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_genes_cor_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x,m,d,ex,hc,ph,cl)
dev.off()


#  [tumors, circRNAs] PCA based on variance stabilized (and log2-transformed) counts centered but not scaled
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-meta[!is.na(risk_group), ]
x<-crs.vs[, m$bid]
colnames(x)<-sub('-11-R01', '', colnames(x))
p<-prcomp(t(x), center=T, scale.=F)
v<-round(1000 * p$sdev^2/sum(p$sdev^2))/10 
n<-cbind( data.frame( p$x[, c('PC1', 'PC2')] ), Type=m$risk_group, bid=colnames(x))
cl<-list('Type'=setNames( unique(m$col), unique(m$risk_group) ))
XLIM<-pretty(range(n[, 'PC1']), 2)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(n[, 'PC2']), 2)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(n[, c('PC1', 'PC2', 'Type', 'bid')], aes(PC1, PC2, color=Type)) + 
    geom_point(size=12) + 
    scale_shape_manual(values=18) + 
    scale_fill_manual(name='Type', values=cl$Type) + 
    scale_color_manual(name='Type', values=cl$Type) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=34), axis.title.x=element_text(size=34), 
        axis.title.y=element_text(size=34), axis.text.y=element_text(size=34), 
        legend.text=element_text(size=22), legend.title=element_text(size=22, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.2, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.2, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', v[1], '% variance')) + ylab(paste0('PC2: ', v[2], '% variance'))
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/plotPCA_circRNAs_vsc_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,p,v,n,cl,XLIM,YLIM)
dev.off()


#  [tumors, genes] PCA based on variance stabilized glog2-transformed counts centered but not scaled
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-meta[!is.na(risk_group), ]
x<-gns.vs[, m$bid]
colnames(x)<-sub('-11-R01', '', colnames(x))
p<-prcomp(t(x), center=T, scale.=F)
v<-round(1000 * p$sdev^2/sum(p$sdev^2))/10 
n<-cbind( data.frame( p$x[, c('PC1', 'PC2')] ), Type=m$risk_group, bid=colnames(x))
cl<-list('Type'=setNames( unique(m$col), unique(m$risk_group) ))
#XLIM<-pretty(range(n[, 'PC1']), 2)
#XLIM<-c(XLIM[1], XLIM[length(XLIM)])
#YLIM<-pretty(range(n[, 'PC2']), 2)
#YLIM<-c(YLIM[1], YLIM[length(YLIM)])
XLIM<-c(-150, 250)  #  force limits
YLIM<-c(-150, 100)  #  force limits
ggplot(n[, c('PC1', 'PC2', 'Type', 'bid')], aes(PC1, PC2, color=Type)) + 
    geom_point(size=12) + 
    scale_shape_manual(values=18) + 
    scale_fill_manual(name='Type', values=cl$Type) + 
    scale_color_manual(name='Type', values=cl$Type) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=34), axis.title.x=element_text(size=34), 
        axis.title.y=element_text(size=34), axis.text.y=element_text(size=34), 
        legend.text=element_text(size=22), legend.title=element_text(size=22, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.2, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.2, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=seq(XLIM[1], XLIM[2], length.out=5)) + 
    scale_y_continuous(lim=YLIM, breaks=seq(YLIM[1], YLIM[2], length.out=6)) + 
    xlab(paste0('PC1: ', v[1], '% variance')) + ylab(paste0('PC2: ', v[2], '% variance'))
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/plotPCA_genes_vsc_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,p,v,n,cl,XLIM,YLIM)
dev.off()

#}}}


#  [tumors, circRNAs] all pairwise DE analyses
#{{{

#  [MNA vs HR_nMNA] DESeq2 forcing library sizes from the gene counts and including sex in the design
m<-meta[risk_group %in% c('MNA', 'HR_nMNA') ]
M<-crs.cnt[, m$bid ]
CT<-data.frame(m[, c('sex', 'risk_group', 'col'), with=F], row.names=m$bid)
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(risk_group=factor(CT$risk_group), sex=factor(CT$sex), row.names=rownames(CT)), design=~sex+risk_group) 
DDS$risk_group<-relevel(DDS$risk_group, 'HR_nMNA')  #  define first level so that the comparison log2(level 2/level 1) is done
sizeFactors(DDS)<-gns.sf[ colnames(DDS) ]
d<-do_deseq2(DDS=DDS, CT=CT, Type='circRNAs', condition='risk_group', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='-11-R01$', list(par=list(mar=c(5.0, 7.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=3))
rm(m,M,CT,DDS)


#  [HR_nMNA vs IMR] DESeq2 forcing library sizes from the gene counts and including sex in the design
m<-meta[risk_group %in% c('IMR', 'HR_nMNA') ]
M<-crs.cnt[, m$bid ]
CT<-data.frame(m[, c('sex', 'risk_group', 'col'), with=F], row.names=m$bid)
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(risk_group=factor(CT$risk_group), sex=factor(CT$sex), row.names=rownames(CT)), design=~sex+risk_group) 
DDS$risk_group<-relevel(DDS$risk_group, 'IMR')  #  define first level so that the comparison log2(level 2/level 1) is done
sizeFactors(DDS)<-gns.sf[ colnames(DDS) ]
d<-do_deseq2(DDS=DDS, CT=CT, Type='circRNAs', condition='risk_group', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='-11-R01$', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)


#  [HR_nMNA vs LR] DESeq2 forcing library sizes from the gene counts and including sex in the design
m<-meta[risk_group %in% c('LR', 'HR_nMNA') ]
M<-crs.cnt[, m$bid ]
CT<-data.frame(m[, c('sex', 'risk_group', 'col'), with=F], row.names=m$bid)
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(risk_group=factor(CT$risk_group), sex=factor(CT$sex), row.names=rownames(CT)), design=~sex+risk_group) 
DDS$risk_group<-relevel(DDS$risk_group, 'LR')  #  define first level so that the comparison log2(level 2/level 1) is done
sizeFactors(DDS)<-gns.sf[ colnames(DDS) ]
d<-do_deseq2(DDS=DDS, CT=CT, Type='circRNAs', condition='risk_group', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='-11-R01$', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)


#  [IMR vs LR] DESeq2 forcing library sizes from the gene counts and including sex in the design
m<-meta[risk_group %in% c('LR', 'IMR') ]
M<-crs.cnt[, m$bid ]
CT<-data.frame(m[, c('sex', 'risk_group', 'col'), with=F], row.names=m$bid)
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(risk_group=factor(CT$risk_group), sex=factor(CT$sex), row.names=rownames(CT)), design=~sex+risk_group) 
DDS$risk_group<-relevel(DDS$risk_group, 'LR')  #  define first level so that the comparison log2(level 2/level 1) is done
sizeFactors(DDS)<-gns.sf[ colnames(DDS) ]
d<-do_deseq2(DDS=DDS, CT=CT, Type='circRNAs', condition='risk_group', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='-11-R01$', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)


#  [LR vs ST4S] DESeq2 forcing library sizes from the gene counts and including sex in the design
m<-meta[risk_group %in% c('ST4S', 'LR') ]
M<-crs.cnt[, m$bid ]
CT<-data.frame(m[, c('sex', 'risk_group', 'col'), with=F], row.names=m$bid)
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(risk_group=factor(CT$risk_group), sex=factor(CT$sex), row.names=rownames(CT)), design=~sex+risk_group) 
DDS$risk_group<-relevel(DDS$risk_group, 'LR')  #  define first level so that the comparison log2(level 2/level 1) is done
sizeFactors(DDS)<-gns.sf[ colnames(DDS) ]
d<-do_deseq2(DDS=DDS, CT=CT, Type='circRNAs', condition='risk_group', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='-11-R01$', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)

#}}}


#  [tumors, genes] all pairwise DE analyses
#{{{

#  [MNA vs HR_nMNA] DESeq2 including sex in the design
m<-meta[risk_group %in% c('MNA', 'HR_nMNA') ]
M<-gns.cnt[, m$bid ]
CT<-data.frame(m[, c('sex', 'risk_group', 'col'), with=F], row.names=m$bid)
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(risk_group=factor(CT$risk_group), sex=factor(CT$sex), row.names=rownames(CT)), design=~sex+risk_group) 
DDS$risk_group<-relevel(DDS$risk_group, 'HR_nMNA')  #  define first level so that the comparison log2(level 2/level 1) is done
d<-do_deseq2(DDS=DDS, CT=CT, Type='genes', condition='risk_group', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='-11-R01$', list(par=list(mar=c(5.0, 7.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=3))
rm(m,M,CT,DDS)


#  [HR_nMNA vs IMR] DESeq2 including sex in the design
m<-meta[risk_group %in% c('IMR', 'HR_nMNA') ]
M<-gns.cnt[, m$bid ]
CT<-data.frame(m[, c('sex', 'risk_group', 'col'), with=F], row.names=m$bid)
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(risk_group=factor(CT$risk_group), sex=factor(CT$sex), row.names=rownames(CT)), design=~sex+risk_group) 
DDS$risk_group<-relevel(DDS$risk_group, 'IMR')  #  define first level so that the comparison log2(level 2/level 1) is done
d<-do_deseq2(DDS=DDS, CT=CT, Type='genes', condition='risk_group', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='-11-R01$', list(par=list(mar=c(5.0, 7.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=3))
rm(m,M,CT,DDS)


#  [HR_nMNA vs LR] DESeq2 including sex in the design
m<-meta[risk_group %in% c('LR', 'HR_nMNA') ]
M<-gns.cnt[, m$bid ]
CT<-data.frame(m[, c('sex', 'risk_group', 'col'), with=F], row.names=m$bid)
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(risk_group=factor(CT$risk_group), sex=factor(CT$sex), row.names=rownames(CT)), design=~sex+risk_group) 
DDS$risk_group<-relevel(DDS$risk_group, 'LR')  #  define first level so that the comparison log2(level 2/level 1) is done
d<-do_deseq2(DDS=DDS, CT=CT, Type='genes', condition='risk_group', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='-11-R01$', list(par=list(mar=c(5.0, 7.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=3))
rm(m,M,CT,DDS)


#  [IMR vs LR] DESeq2 including sex in the design
m<-meta[risk_group %in% c('LR', 'IMR') ]
M<-gns.cnt[, m$bid ]
CT<-data.frame(m[, c('sex', 'risk_group', 'col'), with=F], row.names=m$bid)
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(risk_group=factor(CT$risk_group), sex=factor(CT$sex), row.names=rownames(CT)), design=~sex+risk_group) 
DDS$risk_group<-relevel(DDS$risk_group, 'LR')  #  define first level so that the comparison log2(level 2/level 1) is done
d<-do_deseq2(DDS=DDS, CT=CT, Type='genes', condition='risk_group', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='-11-R01$', list(par=list(mar=c(5.0, 7.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=3))
rm(m,M,CT,DDS)


#  [ST4S vs LR] DESeq2 including sex in the design
m<-meta[risk_group %in% c('ST4S', 'LR') ]
M<-gns.cnt[, m$bid ]
CT<-data.frame(m[, c('sex', 'risk_group', 'col'), with=F], row.names=m$bid)
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(risk_group=factor(CT$risk_group), sex=factor(CT$sex), row.names=rownames(CT)), design=~sex+risk_group) 
DDS$risk_group<-relevel(DDS$risk_group, 'LR')  #  define first level so that the comparison log2(level 2/level 1) is done
d<-do_deseq2(DDS=DDS, CT=CT, Type='genes', condition='risk_group', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='-11-R01$', list(par=list(mar=c(5.0, 7.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=3))
rm(m,M,CT,DDS)

#}}}


#  [cell models, circRNAs] DE analyses
#{{{

#  [MYCN +Tet 4h vs ETOH 4h] DESeq2 forcing library sizes from the gene counts 
m<-meta[grepl('CB-SKNAS-TR-MYCN', bid) & treatment %in% c('ETOH 4h', '+Tet 4h') ]
M<-crs.cnt[, m$bid ]
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=m$bid)
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, 'ETOH 4h')  #  define first level so that the comparison log2(level 2/level 1) is done
sizeFactors(DDS)<-gns.sf[ colnames(DDS) ]
d<-do_deseq2(DDS=DDS, CT=CT, Type='circRNAs', condition='Treatment', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)


#  [MYCN +Tet 48h vs ETOH 48h] DESeq2 forcing library sizes from the gene counts 
m<-meta[grepl('CB-SKNAS-TR-MYCN', bid) & treatment %in% c('ETOH 48h', '+Tet 48h') ]
M<-crs.cnt[, m$bid ]
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=m$bid)
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, 'ETOH 48h')  #  define first level so that the comparison log2(level 2/level 1) is done
sizeFactors(DDS)<-gns.sf[ colnames(DDS) ]
d<-do_deseq2(DDS=DDS, CT=CT, Type='circRNAs', condition='Treatment', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)


#  [MYCN +Tet 48h vs +Tet 4h ] DESeq2 forcing library sizes from the gene counts 
m<-meta[grepl('CB-SKNAS-TR-MYCN', bid) & treatment %in% c('+Tet 4h', '+Tet 48h') ]
M<-crs.cnt[, m$bid ]
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=m$bid)
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, '+Tet 4h')  #  define first level so that the comparison log2(level 2/level 1) is done
sizeFactors(DDS)<-gns.sf[ colnames(DDS) ]
d<-do_deseq2(DDS=DDS, CT=CT, Type='circRNAs', condition='Treatment', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)

#}}}


#  [cell models, genes] DE analyses 
#{{{

#  [MYCN +Tet 4h vs ETOH 4h] DESeq2 
m<-meta[grepl('CB-SKNAS-TR-MYCN', bid) & treatment %in% c('ETOH 4h', '+Tet 4h') ]
M<-gns.cnt[, m$bid ]
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=m$bid)
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, 'ETOH 4h')  #  define first level so that the comparison log2(level 2/level 1) is done
d<-do_deseq2(DDS=DDS, CT=CT, Type='genes', condition='Treatment', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)


#  [MYCN +Tet 48h vs ETOH 48h] DESeq2 
m<-meta[grepl('CB-SKNAS-TR-MYCN', bid) & treatment %in% c('ETOH 48h', '+Tet 48h') ]
M<-gns.cnt[, m$bid ]
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=m$bid)
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, 'ETOH 48h')  #  define first level so that the comparison log2(level 2/level 1) is done
d<-do_deseq2(DDS=DDS, CT=CT, Type='genes', condition='Treatment', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)


#  [MYCN +Tet 48h vs +Tet 4h ] DESeq2
m<-meta[grepl('CB-SKNAS-TR-MYCN', bid) & treatment %in% c('+Tet 4h', '+Tet 48h') ]
M<-gns.cnt[, m$bid ]
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=m$bid)
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=M, colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, '+Tet 4h')  #  define first level so that the comparison log2(level 2/level 1) is done
d<-do_deseq2(DDS=DDS, CT=CT, Type='genes', condition='Treatment', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)

#}}}

#}}}



#  [tumors + cell lines] process of all the DESeq2 results into HTML-table-ready objects
#                        do MSigDB C2 enrichment analysis
#                        do MSigDB GSEA analysis
#
#                        combine the MNA vs HR_nMNA, MYCN +-Tet 4h, MYCN +- 48h, MYCN +- 120h results
#
#{{{

#  [MNA vs HR_nMNA, circRNAs+genes] enrichment in MYCN targets
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)


#  circRNAs
#{{{

#  load circRNAs DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_MNA_HR_nMNA.RData')
circ<-subset(RES, baseMean>0)
rownames(circ)<-sub('\\.[0-9]*$', '', rownames(circ))
rm(CND, DDS, PCA, VE, VSC, RES)


#  MYCN targets 
#  separate to induced and repressed
load('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/MYCN_targets.RData')
IND<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'induced', gene_id])
REP<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'repressed', gene_id])
stopifnot( length(IND)+length(REP)==nrow(MYCN) )
rm(MYCN)


#  [induced targets, significantly up-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(circ, log2FoldChange>0 & padj<0.05) )
rst<-setdiff( rownames(circ), grp )    #  not significantly up-DE which will include significantly down-DE for example
stopifnot( length(grp) + length(rst)==nrow(circ) )
fisher.test(data.frame('in'=c(length(intersect(grp, IND)), 
                              length(setdiff(grp, IND))), 
                       'out'=c(length(intersect(rst, IND)), 
                               length(setdiff(rst, IND)))), alternative='greater')$p.value
#
#  => 0.2315609377


#  induced targets in the significantly down-DE genes
grp<-rownames( subset(circ, log2FoldChange<0 & padj<0.05) )
intersect(grp, IND)  #  ENSG00000144741


#  [repressed targets, significantly down-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(circ, log2FoldChange<0 & padj<0.05) )
rst<-setdiff( rownames(circ), grp )    #  not significantly down-DE which will include significantly up-DE for example
stopifnot( length(grp) + length(rst)==nrow(circ) )
fisher.test(data.frame('in'=c(length(intersect(grp, REP)), 
                              length(setdiff(grp, REP))), 
                       'out'=c(length(intersect(rst, REP)), 
                               length(setdiff(rst, REP)))), alternative='greater')$p.value
#
#  => 0.005987559572


#  no repressed targets in the significantly up-DE genes
grp<-rownames( subset(circ, log2FoldChange>0 & padj<0.05) )
intersect(grp, REP)

#}}}


#  genes
#{{{

#  load genes DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_MNA_HR_nMNA.RData')
gns<-subset(RES, baseMean>0)
rownames(gns)<-sub('\\.[0-9]*$', '', rownames(gns))
rm(CND, DDS, PCA, VE, VSC, RES)


#  MYCN targets 
#  separate to induced and repressed
load('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/MYCN_targets.RData')
IND<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'induced', gene_id])
REP<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'repressed', gene_id])
stopifnot( length(IND)+length(REP)==nrow(MYCN) )
rm(MYCN)


#  [induced targets, significantly up-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(gns, log2FoldChange>0 & padj<0.05) )
rst<-setdiff( rownames(gns), grp )    #  not significantly up-DE which will include significantly down-DE for example
stopifnot( length(grp) + length(rst)==nrow(gns) )
fisher.test(data.frame('in'=c(length(intersect(grp, IND)), 
                              length(setdiff(grp, IND))), 
                       'out'=c(length(intersect(rst, IND)), 
                               length(setdiff(rst, IND)))), alternative='greater')$p.value
#
#  => 3.929419516e-131


#  induced targets in the significantly down-DE genes
grp<-rownames( subset(gns, log2FoldChange<0 & padj<0.05) )
intersect(grp, IND)


#  [repressed targets, significantly down-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(gns, log2FoldChange<0 & padj<0.05) )
rst<-setdiff( rownames(gns), grp )    #  not significantly down-DE which will include significantly up-DE for example
stopifnot( length(grp) + length(rst)==nrow(gns) )
fisher.test(data.frame('in'=c(length(intersect(grp, REP)), 
                              length(setdiff(grp, REP))), 
                       'out'=c(length(intersect(rst, REP)), 
                               length(setdiff(rst, REP)))), alternative='greater')$p.value
#
#  => 5.558197207e-23


#  no repressed targets in the significantly up-DE genes
grp<-rownames( subset(gns, log2FoldChange>0 & padj<0.05) )
intersect(grp, REP)

#}}}

#}}}


#  [MYCN +Tet 4h vs ETOH 4h, circRNAs+genes] enrichment in MYCN targets
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)


#  circRNAs
#{{{

#  load circRNAs DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_+Tet 4h_ETOH 4h.RData')
circ<-subset(RES, baseMean>0)
rownames(circ)<-sub('\\.[0-9]*$', '', rownames(circ))
rm(CND, DDS, PCA, VE, VSC, RES)


#  MYCN targets 
#  separate to induced and repressed
load('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/MYCN_targets.RData')
IND<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'induced', gene_id])
REP<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'repressed', gene_id])
stopifnot( length(IND)+length(REP)==nrow(MYCN) )
rm(MYCN)


#  [induced targets, significantly up-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(circ, log2FoldChange>0 & padj<0.05) )
rst<-setdiff( rownames(circ), grp )    #  not significantly up-DE which will include significantly down-DE for example
stopifnot( length(grp) + length(rst)==nrow(circ) )
fisher.test(data.frame('in'=c(length(intersect(grp, IND)), 
                              length(setdiff(grp, IND))), 
                       'out'=c(length(intersect(rst, IND)), 
                               length(setdiff(rst, IND)))), alternative='greater')$p.value
#
#  => 1


#  induced targets in the significantly down-DE genes
grp<-rownames( subset(circ, log2FoldChange<0 & padj<0.05) )
intersect(grp, IND)


#  [repressed targets, significantly down-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(circ, log2FoldChange<0 & padj<0.05) )
rst<-setdiff( rownames(circ), grp )    #  not significantly down-DE which will include significantly up-DE for example
stopifnot( length(grp) + length(rst)==nrow(circ) )
fisher.test(data.frame('in'=c(length(intersect(grp, REP)), 
                              length(setdiff(grp, REP))), 
                       'out'=c(length(intersect(rst, REP)), 
                               length(setdiff(rst, REP)))), alternative='greater')$p.value
#
#  => 1


#  no repressed targets in the significantly up-DE genes
grp<-rownames( subset(circ, log2FoldChange>0 & padj<0.05) )
intersect(grp, REP)

#}}}


#  genes
#{{{

#  load genes DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_+Tet 4h_ETOH 4h.RData')
gns<-subset(RES, baseMean>0)
rownames(gns)<-sub('\\.[0-9]*$', '', rownames(gns))
rm(CND, DDS, PCA, VE, VSC, RES)


#  MYCN targets 
#  separate to induced and repressed
load('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/MYCN_targets.RData')
IND<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'induced', gene_id])
REP<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'repressed', gene_id])
stopifnot( length(IND)+length(REP)==nrow(MYCN) )
rm(MYCN)


#  [induced targets, significantly up-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(gns, log2FoldChange>0 & padj<0.05) )
rst<-setdiff( rownames(gns), grp )    #  not significantly up-DE which will include significantly down-DE for example
stopifnot( length(grp) + length(rst)==nrow(gns) )
fisher.test(data.frame('in'=c(length(intersect(grp, IND)), 
                              length(setdiff(grp, IND))), 
                       'out'=c(length(intersect(rst, IND)), 
                               length(setdiff(rst, IND)))), alternative='greater')$p.value
#
#  => 1.556497422e-46


#  induced targets in the significantly down-DE genes
grp<-rownames( subset(gns, log2FoldChange<0 & padj<0.05) )
intersect(grp, IND)  #  ENSG00000125835 ENSG00000178035 ENSG00000166482


#  [repressed targets, significantly down-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(gns, log2FoldChange<0 & padj<0.05) )
rst<-setdiff( rownames(gns), grp )    #  not significantly down-DE which will include significantly up-DE for example
stopifnot( length(grp) + length(rst)==nrow(gns) )
fisher.test(data.frame('in'=c(length(intersect(grp, REP)), 
                              length(setdiff(grp, REP))), 
                       'out'=c(length(intersect(rst, REP)), 
                               length(setdiff(rst, REP)))), alternative='greater')$p.value
#
#  => 0.01605289623


#  no repressed targets in the significantly up-DE genes
grp<-rownames( subset(gns, log2FoldChange>0 & padj<0.05) )
intersect(grp, REP)  #  ENSG00000078018, ENSG00000068366

#}}}

#}}}


#  [MYCN +Tet 48h vs ETOH 48h, circRNAs+genes] enrichment in MYCN targets
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)


#  circRNAs
#{{{

#  load circRNAs DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_+Tet 48h_ETOH 48h.RData')
circ<-subset(RES, baseMean>0)
rownames(circ)<-sub('\\.[0-9]*$', '', rownames(circ))
rm(CND, DDS, PCA, VE, VSC, RES)


#  MYCN targets 
#  separate to induced and repressed
load('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/MYCN_targets.RData')
IND<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'induced', gene_id])
REP<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'repressed', gene_id])
stopifnot( length(IND)+length(REP)==nrow(MYCN) )
rm(MYCN)


#  [induced targets, significantly up-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(circ, log2FoldChange>0 & padj<0.05) )
rst<-setdiff( rownames(circ), grp )    #  not significantly up-DE which will include significantly down-DE for example
stopifnot( length(grp) + length(rst)==nrow(circ) )
fisher.test(data.frame('in'=c(length(intersect(grp, IND)), 
                              length(setdiff(grp, IND))), 
                       'out'=c(length(intersect(rst, IND)), 
                               length(setdiff(rst, IND)))), alternative='greater')$p.value
#
#  => 1


#  induced targets in the significantly down-DE genes
grp<-rownames( subset(circ, log2FoldChange<0 & padj<0.05) )
intersect(grp, IND)


#  [repressed targets, significantly down-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(circ, log2FoldChange<0 & padj<0.05) )
rst<-setdiff( rownames(circ), grp )    #  not significantly down-DE which will include significantly up-DE for example
stopifnot( length(grp) + length(rst)==nrow(circ) )
fisher.test(data.frame('in'=c(length(intersect(grp, REP)), 
                              length(setdiff(grp, REP))), 
                       'out'=c(length(intersect(rst, REP)), 
                               length(setdiff(rst, REP)))), alternative='greater')$p.value
#
#  => 1


#  no repressed targets in the significantly up-DE genes
grp<-rownames( subset(circ, log2FoldChange>0 & padj<0.05) )
intersect(grp, REP)

#}}}


#  genes
#{{{

#  load genes DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_+Tet 48h_ETOH 48h.RData')
gns<-subset(RES, baseMean>0)
rownames(gns)<-sub('\\.[0-9]*$', '', rownames(gns))
rm(CND, DDS, PCA, VE, VSC, RES)


#  MYCN targets 
#  separate to induced and repressed
load('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/MYCN_targets.RData')
IND<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'induced', gene_id])
REP<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'repressed', gene_id])
stopifnot( length(IND)+length(REP)==nrow(MYCN) )
rm(MYCN)


#  [induced targets, significantly up-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(gns, log2FoldChange>0 & padj<0.05) )
rst<-setdiff( rownames(gns), grp )    #  not significantly up-DE which will include significantly down-DE for example
stopifnot( length(grp) + length(rst)==nrow(gns) )
fisher.test(data.frame('in'=c(length(intersect(grp, IND)), 
                              length(setdiff(grp, IND))), 
                       'out'=c(length(intersect(rst, IND)), 
                               length(setdiff(rst, IND)))), alternative='greater')$p.value
#
#  => 5.587745546e-44


#  induced targets in the significantly down-DE genes
grp<-rownames( subset(gns, log2FoldChange<0 & padj<0.05) )
intersect(grp, IND)  
#  ENSG00000184117, ENSG00000166482, ENSG00000135446, ENSG00000278615, ENSG00000125835, ENSG00000178035, ENSG00000172403, ENSG00000100226,
#  ENSG00000119977


#  [repressed targets, significantly down-DE genes] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#      MYCN target  |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |      x2          |          y2          | 
grp<-rownames( subset(gns, log2FoldChange<0 & padj<0.05) )
rst<-setdiff( rownames(gns), grp )    #  not significantly down-DE which will include significantly up-DE for example
stopifnot( length(grp) + length(rst)==nrow(gns) )
fisher.test(data.frame('in'=c(length(intersect(grp, REP)), 
                              length(setdiff(grp, REP))), 
                       'out'=c(length(intersect(rst, REP)), 
                               length(setdiff(rst, REP)))), alternative='greater')$p.value
#
#  => 2.468565789e-10


#  no repressed targets in the significantly up-DE genes
grp<-rownames( subset(gns, log2FoldChange>0 & padj<0.05) )
intersect(grp, REP)  #  ENSG00000078018, ENSG00000164604, ENSG00000068366, ENSG00000080823

#}}}

#}}}


#  [MYCN +Tet 120h vs ETOH 120h, circRNAs+genes] enrichment in MYCN targets: 
#
#      008000 MYCN Tet-induction 120h.R


#  [MNA vs HR_nMNA, circRNAs+genes] create html-table-ready objects to be loaded in the reports
#                                   do MSigDB C2 enrichment analysis
#                                   do MSigDB GSEA analysis using baseMean>=10 genes and sum(log2FC) for genes with identical name
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(clusterProfiler)


#  load unified circRNAs and metadata
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
meta<-nb.meta


#  [circRNAs] load DE results
#             separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_MNA_HR_nMNA.RData')
all.circ<-RES
up.circ<-subset(RES, log2FoldChange>0 & padj<0.05)
down.circ<-subset(RES, log2FoldChange<0 & padj<0.05)
bid.circ<-colnames(DDS)
meta.circ<-meta[ match( bid.circ, bid ), ]
rm(CND,DDS,PCA,RES,VE,VSC)


#  [genes] load DE results
#          separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_MNA_HR_nMNA.RData')
all.gns<-RES
up.gns<-subset(RES, log2FoldChange>0 & padj<0.05)
down.gns<-subset(RES, log2FoldChange<0 & padj<0.05)
bid.gns<-colnames(DDS)
meta.gns<-meta[ match( bid.gns, bid ), ]
rm(CND,DDS,PCA,RES,VE,VSC,meta)


#  check enrichment of up-/down-regulated circRNAs in up-/down-regulated genes
#{{{

#  [significantly up-regulated circRNAs] Fisher's exact test:
#
#                        |   genes sig up   |   genes not sig up   |
#  ----------------------|------------------|----------------------|
#     circRNAs sig up    |        x1        |          y1          |
#  ----------------------|------------------|----------------------|
#   circRNAs not sig up  |        x2        |          y2          | 
#
fisher.test(data.frame('in'=c(length(intersect(rownames(up.circ), rownames(up.gns))), 
                              length(intersect(setdiff(rownames(all.circ), rownames(up.circ)), rownames(up.gns)))), 
                       'out'=c(length(intersect(rownames(up.circ), setdiff(rownames(all.gns), rownames(up.gns)))), 
                               length(intersect(setdiff(rownames(all.gns), rownames(up.gns)), setdiff(rownames(all.circ), rownames(up.circ)))))))$p.value
#
#  => 2.906822976e-12



#  [significantly down-regulated circRNAs] Fisher's exact test:
#
#                          |   genes sig down   |   genes not sig down   |
#  ------------------------|--------------------|------------------------|
#     circRNAs sig down    |        x1          |          y1            |
#  ------------------------|--------------------|------------------------|
#   circRNAs not sig down  |        x2          |          y2            | 
#
fisher.test(data.frame('in'=c(length(intersect(rownames(down.circ), rownames(down.gns))), 
                              length(intersect(setdiff(rownames(all.circ), rownames(down.circ)), rownames(down.gns)))), 
                       'out'=c(length(intersect(rownames(down.circ), setdiff(rownames(all.gns), rownames(down.gns)))), 
                               length(intersect(setdiff(rownames(all.gns), rownames(down.gns)), setdiff(rownames(all.circ), rownames(down.circ)))))))$p.value
#
#  => 5.547784901e-62



#  [all up-regulated circRNAs (same for all down-regulated circRNAs)] Fisher's exact test:
#
#                    |   genes up   |     genes down   |
#  ------------------|--------------|------------------|
#     circRNAs up    |      x1      |        y1        |
#  ------------------|--------------|------------------|
#    circRNAS down   |      x2      |        y2        | 
#
all.gns.up<-rownames(subset(all.gns, log2FoldChange>0))
all.gns.down<-rownames(subset(all.gns, log2FoldChange<0))
all.circ.up<-rownames(subset(all.circ, log2FoldChange>0))
all.circ.down<-rownames(subset(all.circ, log2FoldChange<0))
fisher.test(data.frame('in'=c(length(intersect(all.circ.up, all.gns.up)),
                              length(intersect(all.circ.down, all.gns.up))),
                       'out'=c(length(intersect(all.circ.up, all.gns.down)),
                               length(intersect(all.circ.down, all.gns.down)))))$p.value
rm(all.gns.up, all.gns.down, all.circ.up, all.circ.down)
#
#  => 1.589789435e-59

#}}}


#  C2 MSigDB enrichment analysis
#  C2 MSigDB GSEA using baseMean>=10 genes and sum(log2FC) for genes with identical name
#{{{

#  load the C2 MSigDB gene sets
c2<-read.gmt('/fast/groups/ag_schulte/work/reference/annotation/MSigDB/c2.all.v7.0.symbols.gmt')


#  [genes] C2 MSigDB enrichment analysis
up.gns.c2<-enricher(up.gns$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.gns$gene_name, TERM2GENE=c2)  #  to access the main result: up.gns.c2@result
down.gns.c2<-enricher(down.gns$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.gns$gene_name, TERM2GENE=c2)


#  [circRNAs] C2 MSigDB enrichment analysis
up.circ.c2<-enricher(up.circ$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.circ$gene_name, TERM2GENE=c2)  #  to access the main result: up.circ.c2@result
down.circ.c2<-enricher(down.circ$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.circ$gene_name, TERM2GENE=c2)


#  [genes, baseMean>=10, add log2FC of genes with identical names] C2 MSigDB GSEA analysis
g<-data.table(data.frame(subset(all.gns, baseMean>=10)))[, .(log2FoldChange=sum(log2FoldChange)), by=.(gene_name)]
g<-sort( setNames( g$log2FoldChange , g$gene_name), decreasing=T ) 
gsea.gns.c2<-GSEA(g, pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, TERM2GENE=c2)
rm(g)


#  [circRNAs, baseMean>=1]
g<-subset(all.circ, baseMean>=1)
g<-sort( setNames( g$log2FoldChange , g$gene_name), decreasing=T ) 
gsea.circ.c2<-GSEA(g, pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, TERM2GENE=c2)
rm(g)

#}}}


#  [circRNAs] rounded mean circular junction coverage across samples for each circRNA isoform
#             build genomic locations as well (CIRI2 uses 1-based coordinate system)
circ<-data.table(data.frame(CIRCS.all[ CIRCS.all$bid %in% bid.circ ]))[, .(jc_count=ceiling(mean(jc_count)), gene_id=unique(gene_id), locus=paste0('[', sub('chr', '', seqnames), ':', start, '-', end, '](http://www.ensembl.org/Homo_sapiens/Location/View?r=', sub('chr', '', seqnames), ':', start, '-', end, ')')), by=.(gene_name, seqnames, start, end, strand)]


#  [circRNAs] collect circRNA isoform mean circular junction coverage and their genomic positions
circ<-circ[, .(n_circs=.N, jc_count=list(jc_count), gene_id=unique(gene_id), locus=list(locus)), by=.(gene_name)]
circ$gene_name<-paste0('[', circ$gene_name, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', circ$gene_name, ')')


#  [circRNAs] form the up-regulated and down-regulated data.frames, ready for kable()
up.circ<-cbind( as.data.frame(circ[ match(rownames(up.circ), circ$gene_id), ]), up.circ[, c(1:2, 6)] ) 
rownames(up.circ)<-up.circ$gene_id
up.circ<-as.data.frame(up.circ)[, c('gene_name', 'n_circs', 'jc_count', 'locus', 'baseMean', 'log2FoldChange', 'padj')]
down.circ<-cbind( as.data.frame(circ[ match(rownames(down.circ), circ$gene_id), ]), down.circ[, c(1:2, 6)] ) 
rownames(down.circ)<-down.circ$gene_id
down.circ<-as.data.frame(down.circ)[, c('gene_name', 'n_circs', 'jc_count', 'locus', 'baseMean', 'log2FoldChange', 'padj')]


#  [genes] form the up-regulated and down-regulated data.frames, ready for kable()
up.gns$gene_name<-paste0('[', up.gns$gene_name, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', up.gns$gene_name, ')')
up.gns<-as.data.frame(up.gns)[, c('gene_name', 'baseMean', 'log2FoldChange', 'padj')]
down.gns$gene_name<-paste0('[', down.gns$gene_name, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', down.gns$gene_name, ')')
down.gns<-as.data.frame(down.gns)[, c('gene_name', 'baseMean', 'log2FoldChange', 'padj')]


#  [circRNA] add host gene baseMean expression to the table
up.circ$host_gene_baseMean<-all.gns[ rownames(up.circ), 'baseMean']
down.circ$host_gene_baseMean<-all.gns[ rownames(down.circ), 'baseMean']


#  save them all
save(up.circ, down.circ, up.gns, down.gns, up.gns.c2, down.gns.c2, gsea.gns.c2, up.circ.c2, down.circ.c2, gsea.circ.c2, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_kable-ready_genes+circRNAs_MNA_HR_nMNA.RData')

#}}}


#  [MYCN +Tet 4h vs ETOH 4h, circRNAs+genes] create html-table-ready objects to be loaded in the reports
#                                            do MSigDB C2 enrichment analysis
#                                            do MSigDB GSEA analysis using baseMean>=10 genes and sum(log2FC) for genes with identical name
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(clusterProfiler)


#  load unified cohort of circRNAs identified in the cell models
#  isolate the current cell model
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_cell_models.RData')
meta<-meta[ cell_model %in% 'MYCN' & treatment %in% c('+Tet 4h', 'ETOH 4h') ]
circ<-CIRCS.all[ CIRCS.all$bid %in% meta$bid ]
rm(CIRCS.all)


#  [circRNAs] load DE results
#             separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_+Tet 4h_ETOH 4h.RData')
all.circ<-RES
up.circ<-subset(RES, log2FoldChange>0 & padj<0.05)
down.circ<-subset(RES, log2FoldChange<0 & padj<0.05)
bid.circ<-colnames(DDS)
meta.circ<-meta[ match( bid.circ, bid ), ]
rm(CND,DDS,PCA,RES,VE,VSC)


#  [genes] load DE results
#          separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_+Tet 4h_ETOH 4h.RData')
all.gns<-RES
up.gns<-subset(RES, log2FoldChange>0 & padj<0.05)
down.gns<-subset(RES, log2FoldChange<0 & padj<0.05)
bid.gns<-colnames(DDS)
meta.gns<-meta[ match( bid.gns, bid ), ]
rm(CND,DDS,PCA,RES,VE,VSC,meta)


#  C2 MSigDB enrichment analysis
#  C2 MSigDB GSEA using baseMean>=10 genes and sum(log2FC) for genes with identical name
#{{{

#  load the C2 MSigDB gene sets
c2<-read.gmt('/fast/groups/ag_schulte/work/reference/annotation/MSigDB/c2.all.v7.0.symbols.gmt')


#  [genes] C2 MSigDB enrichment analysis
up.gns.c2<-enricher(up.gns$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.gns$gene_name, TERM2GENE=c2)  #  to access the main result: up.gns.c2@result
down.gns.c2<-enricher(down.gns$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.gns$gene_name, TERM2GENE=c2)


#  [circRNAs] C2 MSigDB enrichment analysis
up.circ.c2<-enricher(up.circ$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.circ$gene_name, TERM2GENE=c2)  #  to access the main result: up.circ.c2@result
down.circ.c2<-enricher(down.circ$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.circ$gene_name, TERM2GENE=c2)


#  [genes, baseMean>=10, add log2FC of genes with identical names] C2 MSigDB GSEA analysis
g<-data.table(data.frame(subset(all.gns, baseMean>=10)))[, .(log2FoldChange=sum(log2FoldChange)), by=.(gene_name)]
g<-sort( setNames( g$log2FoldChange , g$gene_name), decreasing=T ) 
gsea.gns.c2<-GSEA(g, pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, TERM2GENE=c2)
rm(g)


#  [circRNAs, baseMean>=1]
g<-subset(all.circ, baseMean>=1)
g<-sort( setNames( g$log2FoldChange , g$gene_name), decreasing=T ) 
gsea.circ.c2<-GSEA(g, pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, TERM2GENE=c2)
rm(g)

#}}}


#  [circRNAs] rounded mean circular junction coverage across samples for each circRNA isoform
#             build genomic locations as well (CIRI2 uses 1-based coordinate system)
circ<-data.table(data.frame(circ[ circ$bid %in% bid.circ ]))[, .(jc_count=ceiling(mean(jc_count)), gene_id=unique(gene_id), locus=paste0('[', sub('chr', '', seqnames), ':', start, '-', end, '](http://www.ensembl.org/Homo_sapiens/Location/View?r=', sub('chr', '', seqnames), ':', start, '-', end, ')')), by=.(gene_name, seqnames, start, end, strand)]


#  [circRNAs] collect circRNA isoform mean circular junction coverage and their genomic positions
circ<-circ[, .(n_circs=.N, jc_count=list(jc_count), gene_id=unique(gene_id), locus=list(locus)), by=.(gene_name)]
circ$gene_name<-paste0('[', circ$gene_name, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', circ$gene_name, ')')


#  [circRNAs] form the up-regulated and down-regulated data.frames, ready for kable()
up.circ<-cbind( as.data.frame(circ[ match(rownames(up.circ), circ$gene_id), ]), up.circ[, c(1:2, 6)] ) 
rownames(up.circ)<-up.circ$gene_id
up.circ<-as.data.frame(up.circ)[, c('gene_name', 'n_circs', 'jc_count', 'locus', 'baseMean', 'log2FoldChange', 'padj')]
down.circ<-cbind( as.data.frame(circ[ match(rownames(down.circ), circ$gene_id), ]), down.circ[, c(1:2, 6)] ) 
rownames(down.circ)<-down.circ$gene_id
down.circ<-as.data.frame(down.circ)[, c('gene_name', 'n_circs', 'jc_count', 'locus', 'baseMean', 'log2FoldChange', 'padj')]


#  [genes] form the up-regulated and down-regulated data.frames, ready for kable()
up.gns$gene_name<-paste0('[', up.gns$gene_name, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', up.gns$gene_name, ')')
up.gns<-as.data.frame(up.gns)[, c('gene_name', 'baseMean', 'log2FoldChange', 'padj')]
down.gns$gene_name<-paste0('[', down.gns$gene_name, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', down.gns$gene_name, ')')
down.gns<-as.data.frame(down.gns)[, c('gene_name', 'baseMean', 'log2FoldChange', 'padj')]


#  [circRNA] add host gene baseMean expression to the table
up.circ$host_gene_baseMean<-all.gns[ rownames(up.circ), 'baseMean']
down.circ$host_gene_baseMean<-all.gns[ rownames(down.circ), 'baseMean']


#  save them all
save(up.circ, down.circ, up.gns, down.gns, up.gns.c2, down.gns.c2, gsea.gns.c2, up.circ.c2, down.circ.c2, gsea.circ.c2, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_kable-ready_genes+circRNAs_+Tet 4h_ETOH 4h.RData')

#}}}


#  [MYCN +Tet 48h vs ETOH 48h, circRNAs+genes] create html-table-ready objects to be loaded in the reports
#                                              do MSigDB C2 enrichment analysis
#                                              do MSigDB GSEA analysis using baseMean>=10 genes and sum(log2FC) for genes with identical name
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(clusterProfiler)


#  load unified cohort of circRNAs identified in the cell models
#  isolate the current cell model
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_cell_models.RData')
meta<-meta[ cell_model %in% 'MYCN' & treatment %in% c('+Tet 48h', 'ETOH 48h') ]
circ<-CIRCS.all[ CIRCS.all$bid %in% meta$bid ]
rm(CIRCS.all)


#  [circRNAs] load DE results
#             separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_+Tet 48h_ETOH 48h.RData')
all.circ<-RES
up.circ<-subset(RES, log2FoldChange>0 & padj<0.05)
down.circ<-subset(RES, log2FoldChange<0 & padj<0.05)
bid.circ<-colnames(DDS)
meta.circ<-meta[ match( bid.circ, bid ), ]
rm(CND,DDS,PCA,RES,VE,VSC)


#  [genes] load DE results
#          separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_+Tet 48h_ETOH 48h.RData')
all.gns<-RES
up.gns<-subset(RES, log2FoldChange>0 & padj<0.05)
down.gns<-subset(RES, log2FoldChange<0 & padj<0.05)
bid.gns<-colnames(DDS)
meta.gns<-meta[ match( bid.gns, bid ), ]
rm(CND,DDS,PCA,RES,VE,VSC,meta)


#  C2 MSigDB enrichment analysis
#  C2 MSigDB GSEA using baseMean>=10 genes and sum(log2FC) for genes with identical name
#{{{

#  load the C2 MSigDB gene sets
c2<-read.gmt('/fast/groups/ag_schulte/work/reference/annotation/MSigDB/c2.all.v7.0.symbols.gmt')


#  [genes] C2 MSigDB enrichment analysis
up.gns.c2<-enricher(up.gns$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.gns$gene_name, TERM2GENE=c2)  #  to access the main result: up.gns.c2@result
down.gns.c2<-enricher(down.gns$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.gns$gene_name, TERM2GENE=c2)


#  [circRNAs] C2 MSigDB enrichment analysis
up.circ.c2<-enricher(up.circ$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.circ$gene_name, TERM2GENE=c2)  #  to access the main result: up.circ.c2@result
down.circ.c2<-enricher(down.circ$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=all.circ$gene_name, TERM2GENE=c2)


#  [genes, baseMean>=10, add log2FC of genes with identical names] C2 MSigDB GSEA analysis
g<-data.table(data.frame(subset(all.gns, baseMean>=10)))[, .(log2FoldChange=sum(log2FoldChange)), by=.(gene_name)]
g<-sort( setNames( g$log2FoldChange , g$gene_name), decreasing=T ) 
gsea.gns.c2<-GSEA(g, pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, TERM2GENE=c2)
rm(g)


#  [circRNAs, baseMean>=1]
g<-subset(all.circ, baseMean>=1)
g<-sort( setNames( g$log2FoldChange , g$gene_name), decreasing=T ) 
gsea.circ.c2<-GSEA(g, pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, TERM2GENE=c2)
rm(g)

#}}}


#  [circRNAs] rounded mean circular junction coverage across samples for each circRNA isoform
#             build genomic locations as well (CIRI2 uses 1-based coordinate system)
circ<-data.table(data.frame(circ[ circ$bid %in% bid.circ ]))[, .(jc_count=ceiling(mean(jc_count)), gene_id=unique(gene_id), locus=paste0('[', sub('chr', '', seqnames), ':', start, '-', end, '](http://www.ensembl.org/Homo_sapiens/Location/View?r=', sub('chr', '', seqnames), ':', start, '-', end, ')')), by=.(gene_name, seqnames, start, end, strand)]


#  [circRNAs] collect circRNA isoform mean circular junction coverage and their genomic positions
circ<-circ[, .(n_circs=.N, jc_count=list(jc_count), gene_id=unique(gene_id), locus=list(locus)), by=.(gene_name)]
circ$gene_name<-paste0('[', circ$gene_name, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', circ$gene_name, ')')


#  [circRNAs] form the up-regulated and down-regulated data.frames, ready for kable()
up.circ<-cbind( as.data.frame(circ[ match(rownames(up.circ), circ$gene_id), ]), up.circ[, c(1:2, 6)] ) 
rownames(up.circ)<-up.circ$gene_id
up.circ<-as.data.frame(up.circ)[, c('gene_name', 'n_circs', 'jc_count', 'locus', 'baseMean', 'log2FoldChange', 'padj')]
down.circ<-cbind( as.data.frame(circ[ match(rownames(down.circ), circ$gene_id), ]), down.circ[, c(1:2, 6)] ) 
rownames(down.circ)<-down.circ$gene_id
down.circ<-as.data.frame(down.circ)[, c('gene_name', 'n_circs', 'jc_count', 'locus', 'baseMean', 'log2FoldChange', 'padj')]


#  [genes] form the up-regulated and down-regulated data.frames, ready for kable()
up.gns$gene_name<-paste0('[', up.gns$gene_name, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', up.gns$gene_name, ')')
up.gns<-as.data.frame(up.gns)[, c('gene_name', 'baseMean', 'log2FoldChange', 'padj')]
down.gns$gene_name<-paste0('[', down.gns$gene_name, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', down.gns$gene_name, ')')
down.gns<-as.data.frame(down.gns)[, c('gene_name', 'baseMean', 'log2FoldChange', 'padj')]


#  [circRNA] add host gene baseMean expression to the table
up.circ$host_gene_baseMean<-all.gns[ rownames(up.circ), 'baseMean']
down.circ$host_gene_baseMean<-all.gns[ rownames(down.circ), 'baseMean']


#  save them all
save(up.circ, down.circ, up.gns, down.gns, up.gns.c2, down.gns.c2, gsea.gns.c2, up.circ.c2, down.circ.c2, gsea.circ.c2, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_kable-ready_genes+circRNAs_+Tet 48h_ETOH 48h.RData')

#}}}


#  [MYCN +Tet 4h vs ETOH 4h, +Tet 48h vs ETOH 48h, MNA vs HR_nMNA] commonly DE genes and circRNAs
#                                                                  create html-table-ready objects
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)


#  load MYCN +Tet 4h vs ETOH 4h
#{{{

#  load gene DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_+Tet 4h_ETOH 4h.RData')
gns.mycn.4h<-subset(RES, padj<0.05)
gns.all.mycn.4h<-RES
rm(CND, DDS, PCA, VE, VSC, RES)


#  load circRNAs DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_+Tet 4h_ETOH 4h.RData')
circ.mycn.4h<-subset(RES, padj<0.05)
circ.all.mycn.4h<-RES
rm(CND, DDS, PCA, VE, VSC, RES)


#  load kable-ready up- and down-regulated circRNAs and genes
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_kable-ready_genes+circRNAs_+Tet 4h_ETOH 4h.RData')
gns.html.mycn.4h<-rbind(up.gns, down.gns)
circ.html.mycn.4h<-rbind(up.circ, down.circ)
rm(up.gns, down.gns, up.circ, down.circ)

#}}}


#  load MYCN +Tet 48h vs ETOH 48h
#{{{

#  load gene DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_+Tet 48h_ETOH 48h.RData')
gns.mycn.48h<-subset(RES, padj<0.05)
gns.all.mycn.48h<-RES
rm(CND, DDS, PCA, VE, VSC, RES)


#  load circRNAs DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_+Tet 48h_ETOH 48h.RData')
circ.mycn.48h<-subset(RES, padj<0.05)
circ.all.mycn.48h<-RES
rm(CND, DDS, PCA, VE, VSC, RES)


#  load kable-ready up- and down-regulated circRNAs and genes
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_kable-ready_genes+circRNAs_+Tet 48h_ETOH 48h.RData')
gns.html.mycn.48h<-rbind(up.gns, down.gns)
circ.html.mycn.48h<-rbind(up.circ, down.circ)
rm(up.gns, down.gns, up.circ, down.circ)

#}}}


#  load MYCN +Tet 120h vs ETOH 120h
#{{{

#  load gene DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_genes_+Tet 120h_ETOH 120h.RData')
gns.mycn.120h<-subset(RES, padj<0.05)
gns.all.mycn.120h<-RES
rm(CND, DDS, PCA, VE, VSC, RES)


#  load circRNAs DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_circRNAs_+Tet 120h_ETOH 120h.RData')
circ.mycn.120h<-subset(RES, padj<0.05)
circ.all.mycn.120h<-RES
rm(CND, DDS, PCA, VE, VSC, RES)


#  load kable-ready up- and down-regulated circRNAs and genes
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw//BATCH_202006/DESeq2_kable-ready_genes+circRNAs_+Tet 120h_ETOH 120h.RData')
gns.html.mycn.120h<-rbind(up.gns, down.gns)
circ.html.mycn.120h<-rbind(up.circ, down.circ)
rm(up.gns, down.gns, up.circ, down.circ)

#}}}


#  load MNA vs HR_nMNA results
#{{{

#  load gene DE between MNA vs HR_nMNA
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_MNA_HR_nMNA.RData')
gns.mna<-subset(RES, padj<0.05)
gns.all.mna<-RES
rm(CND, DDS, PCA, VE, VSC, RES)


#  [MNA vs HR_nMNA] load circRNAs DE between MNA vs HR_nMNA
#                   remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_MNA_HR_nMNA.RData')
circ.mna<-subset(RES, padj<0.05)
circ.all.mna<-RES
rm(CND, DDS, PCA, VE, VSC, RES)


#  load kable-ready up- and down-regulated circRNAs and genes
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_kable-ready_genes+circRNAs_MNA_HR_nMNA.RData')
gns.html.mna<-rbind(up.gns, down.gns)
circ.html.mna<-rbind(up.circ, down.circ)
rm(up.gns, down.gns, up.circ, down.circ)

#}}}


#  [MYCN +Tet 4h vs ETOH 4h and MNA vs HR_nMNA] common DE genes
#                                               common DE circRNAs
common<-intersect( rownames(gns.mycn.4h), rownames(gns.mna) )
common.gns.4h<-cbind( gns.html.mycn.4h[ common, -1], gns.html.mna[common, c(2:4, 1)] )
colnames(common.gns.4h)<-c('MYCN_4h_baseMean', 'MYCN_4h_log2FC', 'MYCN_4h_padj', 'MNA_baseMean', 'MNA_log2FC', 'MNA_padj', 'gene_name')
common<-intersect( rownames(circ.mycn.4h), rownames(circ.mna) )
common.circ.4h<-cbind( circ.html.mycn.4h[ common, c('baseMean', 'log2FoldChange', 'padj', 'host_gene_baseMean')], 
                       circ.html.mna[common, c('baseMean', 'log2FoldChange', 'padj', 'host_gene_baseMean', 'gene_name')] )
colnames(common.circ.4h)<-c('MYCN_4h_baseMean', 'MYCN_4h_log2FC', 'MYCN_4h_padj', 'MYCN_4h_host_gene_baseMean', 
                            'MNA_baseMean', 'MNA_log2FC', 'MNA_padj', 'MNA_host_gene_baseMean', 'gene_name')
rm(common)


#  [MYCN +Tet 48h vs ETOH 48h and MNA vs HR_nMNA] common DE genes
#                                                 common DE circRNAs
common<-intersect( rownames(gns.mycn.48h), rownames(gns.mna) )
common.gns.48h<-cbind( gns.html.mycn.48h[ common, -1], gns.html.mna[common, c(2:4, 1)] )
colnames(common.gns.48h)<-c('MYCN_48h_baseMean', 'MYCN_48h_log2FC', 'MYCN_48h_padj', 'MNA_baseMean', 'MNA_log2FC', 'MNA_padj', 'gene_name')
common<-intersect( rownames(circ.mycn.48h), rownames(circ.mna) )
common.circ.48h<-cbind( circ.html.mycn.48h[ common, c('baseMean', 'log2FoldChange', 'padj', 'host_gene_baseMean')], 
                       circ.html.mna[common, c('baseMean', 'log2FoldChange', 'padj', 'host_gene_baseMean', 'gene_name')] )
colnames(common.circ.48h)<-c('MYCN_48h_baseMean', 'MYCN_48h_log2FC', 'MYCN_48h_padj', 'MYCN_48h_host_gene_baseMean', 
                            'MNA_baseMean', 'MNA_log2FC', 'MNA_padj', 'MNA_host_gene_baseMean', 'gene_name')
rm(common)


#  [MYCN +Tet 120h vs ETOH 120h and MNA vs HR_nMNA] common DE genes
#                                                   common DE circRNAs
common<-intersect( rownames(gns.mycn.120h), rownames(gns.mna) )
common.gns.120h<-cbind( gns.html.mycn.120h[ common, -1], gns.html.mna[common, c(2:4, 1)] )
colnames(common.gns.120h)<-c('MYCN_120h_baseMean', 'MYCN_120h_log2FC', 'MYCN_120h_padj', 'MNA_baseMean', 'MNA_log2FC', 'MNA_padj', 'gene_name')
common<-intersect( rownames(circ.mycn.120h), rownames(circ.mna) )
common.circ.120h<-cbind( circ.html.mycn.120h[ common, c('baseMean', 'log2FoldChange', 'padj', 'host_gene_baseMean')], 
                       circ.html.mna[common, c('baseMean', 'log2FoldChange', 'padj', 'host_gene_baseMean', 'gene_name')] )
colnames(common.circ.120h)<-c('MYCN_120h_baseMean', 'MYCN_120h_log2FC', 'MYCN_120h_padj', 'MYCN_120h_host_gene_baseMean', 
                            'MNA_baseMean', 'MNA_log2FC', 'MNA_padj', 'MNA_host_gene_baseMean', 'gene_name')
rm(common)


#  common between all the comparisons 
common<-Reduce(intersect, list(rownames(common.gns.4h), rownames(common.gns.48h), rownames(common.gns.120h)))
common.gns.4h.48h.120h<-cbind( common.gns.4h[ common, 1:3], common.gns.48h[common, 1:3], common.gns.120h[common, ] )
common<-Reduce(intersect, list(rownames(common.circ.4h), rownames(common.circ.48h), rownames(common.circ.120h)))
common.circ.4h.48h.120h<-cbind( common.circ.4h[ common, 1:3], common.circ.48h[common, 1:3], common.circ.120h[common, ] )
rm(common)


#  common between 48h, 120h and MNA only 
common<-Reduce(intersect, list(rownames(common.gns.48h), rownames(common.gns.120h)))
common.gns.48h.120h<-cbind( common.gns.48h[common, 1:3], common.gns.120h[common, ] )
common<-Reduce(intersect, list(rownames(common.circ.48h), rownames(common.circ.120h)))
common.circ.48h.120h<-cbind(common.circ.48h[common, 1:3], common.circ.120h[common, ] )
rm(common)


#  save kable-ready results
save(list=ls(pattern='common\\.'), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_kable-ready_genes+circRNAs_MYCN+MNA.RData')

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_kable-ready_genes+circRNAs_MNA_HR_nMNA.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_kable-ready_genes+circRNAs_+Tet 4h_ETOH 4h.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_kable-ready_genes+circRNAs_+Tet 48h_ETOH 48h.RData
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_kable-ready_genes+circRNAs_MYCN+MNA.RData



#  [TMM-normalized circRNA DE analysis] is global downregulation of circRNAs in MNA vs HR_nMNA an effect of the global mRNA upregulation?
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(vsn)
library(limma)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(pheatmap)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_plotMD.R')


#  remove failed samples and Pilot samples
#  keep only MNA + HR_nMNA primary tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta<-meta[ !(failed) & !grepl('CBPilote', bid) & risk_group %in% c('MNA', 'HR_nMNA') ]


#  load reference 
#  discard chrM genes
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene']
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]
seqlevels(hsa, pruning.mode='coarse')<-seqlevels(hsa)[ grep('chrM', seqlevels(hsa), invert=T) ]


#  load our cohort of circRNAs
#  keep only MNA + HR_nMNA primary tumor predictions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
circs<-CIRCS.all[ CIRCS.all$bid %in% meta$bid ]


#  load featureCounts 
#  discard chrM genes
#  compute RAW count DGE to use for DE analyses (should not be normalized since normalization will be done there)
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')
gns<-unlist(totalrna[ meta$bid ])
gns$bid<-sub('\\.[0-9]*$', '', rownames(gns))
rownames(gns)<-NULL
gns<-gns[ gns$gene_id %in% hsa$gene_id, ]
gns<-data.table(gns[, c('bid', 'gene_id', 'length', 'counts')])
gns.cnt<-dcast(gns, bid ~ gene_id, value.var='counts', fun.aggregate=sum)
gns.cnt<-t(data.frame(gns.cnt[, -1], row.names=gns.cnt[, bid], check.names=F))[,  meta$bid]  #  RAW unnormalized counts


#  compute the TMM-based normalization factor for the counts 
#  is there a systematic difference in the normalization factors between MNA and HR_nMNA?
gns.DGE<-DGEList(counts=gns.cnt)
gns.DGE<-calcNormFactors(gns.DGE, method='TMM')
stopifnot( all.equal( rownames(gns.DGE$samples), meta$bid) )
gns.nf.rg<-split(setNames(gns.DGE$samples$norm.factors, rownames(gns.DGE$samples)), meta$risk_group)
wilcox.test(x=gns.nf.rg[['MNA']], y=gns.nf.rg[['HR_nMNA']], alternative='less')$p.value  #  0.004045581733


#  construct DGEs for the rounded up normalized counts and the TPMs
N<-data.table(bid=rownames(gns.DGE$samples), NF=gns.DGE$samples$norm.factors)
gns<-N[ gns, on='bid']
gns<-gns[, counts:=counts/NF][, NF:=NULL]
gns<-gns[, .(gene_id=gene_id, counts=counts, tpm=1e6*counts/length/sum(counts/length), cpm=1e6*counts/sum(counts)), by=.(bid)]
gns.tpm<-dcast(gns, bid ~ gene_id, value.var='tpm', fun.aggregate=sum)
gns.tpm<-t(data.frame(gns.tpm[, -1], row.names=gns.tpm[, bid], check.names=F))[,  meta$bid]
gns.norm<-dcast(gns, bid ~ gene_id, value.var='counts', fun.aggregate=sum)
gns.norm<-t(data.frame(gns.norm[, -1], row.names=gns.norm[, bid], check.names=F))[,  meta$bid]  #  TMM-normalized counts


#  add normalization factor to circRNAs to normalize raw counts at will
#  compute CPMs based on the gene-level effective library sizes
N<-setNames(gns.DGE$samples$norm.factors, rownames(gns.DGE$samples))
circs$norm.factor<-N[ circs$bid ]
N<-setNames(gns.DGE$samples$lib.size*gns.DGE$samples$norm.factors, rownames(gns.DGE$samples))
circs$nreads<-N[ circs$bid ]  #  effective library sizes
circs$cpm<-circs$jc_count/circs$nreads*1e6
rm(N, totalrna, gns)


#  remove unexpressed genes 
gns.cnt<-gns.cnt[rowSums(gns.cnt)!=0,  ]
gns.norm<-gns.norm[rowSums(gns.norm)!=0,  ]
gns.tpm<-gns.tpm[rowSums(gns.tpm)!=0,  ]


#  compute variance-stabilized gene counts based on the normalized count matrix
#  PCA on variance-stabilized (and glog2-transformed) gene counts centered but not scaled
gns.vs<-varianceStabilizingTransformation(ceiling(gns.norm), fitType='local')
gns.pca<-prcomp(t(gns.vs), center=T, scale.=F)
gns.ve<-round(1000 * gns.pca$sdev^2/sum(gns.pca$sdev^2))/10 


#  summarize circRNA RAW counts at the gene level 
#  summarize circRNA NORMALIZED counts at the gene level and compute the variance-stabilized version of them
#  PCA on variance-stabilized (and glog2-transformed) counts centered by not scaled
crs.cnt<-data.table(data.frame(mcols(circs)[, c('gene_id', 'jc_count', 'bid')]))[, .(jc_count=sum(jc_count)), by=.(bid, gene_id)]
crs.cnt<-dcast(crs.cnt, bid ~ gene_id, value.var='jc_count', fun.aggregate=sum)
crs.cnt<-t(as.matrix(data.frame(crs.cnt[, -1], row.names=crs.cnt$bid, check.names=F)))[, meta$bid]
crs.norm<-data.table(data.frame(mcols(circs)[, c('gene_id', 'jc_count', 'norm.factor', 'bid')]))[, .(jc_count=ceiling(sum(jc_count/norm.factor))), by=.(bid, gene_id)]
crs.norm<-dcast(crs.norm, bid ~ gene_id, value.var='jc_count', fun.aggregate=sum)
crs.norm<-t(as.matrix(data.frame(crs.norm[, -1], row.names=crs.norm$bid, check.names=F)))[, meta$bid]
crs.vs<-varianceStabilizingTransformation(crs.norm, fitType='local')
crs.pca<-prcomp(t(crs.vs), center=T, scale.=F)
crs.ve<-round(1000 * crs.pca$sdev^2/sum(crs.pca$sdev^2))/10 


#  prepare the annotations
hsa<-data.frame(mcols(hsa))
gns.annot<-data.frame(hsa[ match(rownames(gns.cnt), hsa$gene_id), ], row.names=hsa$gene_id[match(rownames(gns.cnt), hsa$gene_id)], check.names=F)
stopifnot( all.equal( rownames(gns.cnt), rownames(gns.annot) ) )
crs.annot<-data.frame(hsa[ match(rownames(crs.cnt), hsa$gene_id), ], row.names=hsa$gene_id[match(rownames(crs.cnt), hsa$gene_id)], check.names=F)
stopifnot( all.equal( rownames(crs.cnt), rownames(crs.annot) ) )


#  create the metadata matrix and choose the baseline levels 
CT<-data.frame(meta[, c('col', 'risk_group'), with=F], row.names=meta$bid)
colnames(CT)<-c('col', 'Risk')
CT$Risk<-factor(CT$Risk, levels=c('HR_nMNA', 'MNA'))


#  create the design matrix as a single factor which is straightforward to understand
group<-CT$Risk
DESIGN<-model.matrix(~0+group)
colnames(DESIGN)<-sub('^group', '', colnames(DESIGN))
rm(group)


#  create the DGEList using as library sizes the raw unnormalized gene count-based ones
#
#  N.B. we recompute the TMM normalization factors now using the circRNA raw counts and the gene-based library sizes
#      in order to mitigate the global downregulation effect. In the end, we still find the circRNAs globally downregulated.
#  
#  N.B. if we use the gene-based TMM normalization factors then the effect remains as pronounced as in the DESeq2 results.
#
DGE<-DGEList(counts=crs.cnt, genes=crs.annot, group=CT$Risk)
keep<-filterByExpr(DGE, min.count=2)
table(keep)
# keep
# FALSE  TRUE 
#   674  1628 
# 
DGE<-DGE[keep, , keep.lib.sizes=F]
stopifnot( all.equal( rownames(DGE$samples), rownames(gns.DGE$samples) ) )
DGE$samples$lib.size<-gns.DGE$samples$lib.size  #  set library sizes to be the library sizes based on raw unnormalized gene counts
DGE<-calcNormFactors(DGE, method='TMM')         #  recompute TMM normalization factors to mitigate the global downregulation effect
#DGE$samples$norm.factors<-gns.DGE$samples$norm.factors   #  use the gene-based TMM normalization factors to see similar results to DESeq2


#  is there a systematic difference in the normalization factors between MNA and HR_nMNA?
crs.nf.rg<-split(setNames(DGE$samples$norm.factors, rownames(DGE$samples)), meta$risk_group)
wilcox.test(x=crs.nf.rg[['MNA']], y=crs.nf.rg[['HR_nMNA']], alternative='less')$p.value  #  0.007596766998


#  voom transformation
DGE.voom<-voom(DGE, design=DESIGN, plot=T)  #  log2-CPMs after this point!!!!
dev.off()
rm(keep)


#  construct the contrasts of interest
FCC<-log2(1.0)
CONT.MAT<-makeContrasts(MNA_HR_nMNA='MNA-HR_nMNA', levels=DESIGN)
FIT<-lmFit(DGE.voom, design=DESIGN)
FIT<-contrasts.fit(FIT, CONT.MAT)
FIT<-treat(FIT, lfc=FCC, trend=T, robust=T)
colnames(FIT$coefficients)


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel')
COEF<-colnames(FIT$coefficients)
RES<-topTreat(FIT, coef=COEF, number=nrow(FIT), adjust.method='BH', sort.by='p', lfc=FCC)
colnames(RES)<-c('gene_id', 'gene_type', 'gene_name', 'log2FoldChange', 'AveExpr', 't', 'p.value', 'padj')
RES<-RES[, c('gene_name', 'gene_id', 'gene_type', 'log2FoldChange', 'AveExpr', 't', 'p.value', 'padj')]
my_plotMD(FIT, COEF, USE.ANNOT='gene_name', PV=0.1, LFC=FCC, TITLE='', TITLE.CEX=2.0,
        FIGDIR='', FIGROOT='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/plotMA_circRNAs_TMM-normalized_',
        #FIGDIR='', FIGROOT='',
        with.identify=F,
        list(par=list(mar=c(5.0, 7.0, 1.0, 4.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4),
        xline=4, xline.padj=-0.3, yline=3, yline.padj=-0.2, height=14, width=20))


#  save for later analysis 
save(list=c('circs', 'DGE', 'DGE.voom', 'COEF', 'RES', 'FIT', 'CT', 'DESIGN', 'FCC', 'CT', 'CONT.MAT', 'hsa', ls(pattern='gns\\.|crs\\.')), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/edgeR_circRNAs_TMM-normalized_MNA_vs_HR_nMNA.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/edgeR_circRNAs_TMM-normalized_MNA_vs_HR_nMNA.RData



#  [publication-ready DE plots]
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(vsn)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()



#  functions
#{{{

gene_results<-function(res=NULL, annot=NULL){
    m<-res[ intersect(annot$gene_id, rownames(res)), c('baseMean', 'log2FoldChange')]
    m$gene_name<-annot$gene_name[ match( rownames(m), annot$gene_id ) ]
    m$col<-annot$col[ match( rownames(m), annot$gene_id ) ]
    return(m)
}

#}}}



#  recycle
x11(width=16, height=16, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  union of SF
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/SpliceAid-F+GO.RData')
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/MSigDB_c2_v7.0_splicing.RData')
SFU<-unique(rbind(data.frame(db[, c('gene_name', 'gene_id')]), msigdbSF[, c('gene_name', 'gene_id')]))
SFU$col<-'khaki'
rm(db, msigdbSF)


#  MYCN targets 
load('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/MYCN_targets.RData')
MYCN.IND<-MYCN[ type %in% 'induced', gene_id]
MYCN.REP<-MYCN[ type %in% 'repressed', gene_id]
stopifnot( length(MYCN.IND)+length(MYCN.REP)==nrow(MYCN) )
rm(MYCN)


#  [genes] MNA vs HR_nMNA
#{{{

#  load DE results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_MNA_HR_nMNA.RData')


#  single genes to annotate specifically
ANNOT<-data.frame(gene_id=c('ENSG00000134323.12', 'ENSG00000233718.7', 'ENSG00000079785.15'), gene_name=c('MYCN', 'MYCNOS', 'DDX1'), col='red3')


#  MYCN targets 
IND<-intersect(MYCN.IND, rownames(RES))
REP<-intersect(MYCN.REP, rownames(RES))
IND<-data.frame(gene_id=IND, gene_name=RES[ IND, 'gene_name'], col='royalblue3', cex=0.2)
REP<-data.frame(gene_id=REP, gene_name=RES[ REP, 'gene_name'], col='springgreen1', cex=0.2)


#  splice factors
SF<-SFU[ SFU$gene_id %in% intersect(SFU$gene_id, rownames(RES)), , drop=F]


# MA-plot with annotations
par(mar=c(5.0, 7.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4)
YLIM<-pretty(range(RES$log2FoldChange, na.rm=T))
YLIM<-c(YLIM[1], tail(YLIM, 1))
XLIM<-pretty(range(RES$baseMean, na.rm=T), 4)
XLIM<-c(0.1, tail(XLIM, 1))
options(scipen=+10)
plot(subset(RES, padj>=0.05)$baseMean, subset(RES, padj>=0.05)$log2FoldChange, type='p', log='x', xaxt='n', pch=21, xlab='', ylab='', main='', ylim=YLIM, xlim=XLIM, cex=0.4, lwd=0, col=adjustcolor('grey60', alpha.f=0.2), bg=adjustcolor('grey60', alpha.f=0.2))
axis(1, at=axTicks(1, usr=log10(XLIM), log=T))
points(subset(RES, padj<0.05)$baseMean, subset(RES, padj<0.05)$log2FoldChange, pch=21, lwd=1, cex=0.8, col=adjustcolor('red3', alpha.f=0.1), bg=adjustcolor('red3', alpha.f=0.5))
abline(h=0, lty=1, lwd=4, col='grey50')
#mtext(CND, side=3, line=0, padj=+1.5, cex=1.4)
mtext(expression(log[2]('fold change')), side=2, line=3, padj=-0.2, cex=2.4, las=3)
mtext('Mean expression', side=1, line=4, padj=-0.3, cex=2.4, las=1)
#  
#  annotate single genes
#
r<-gene_results(RES, ANNOT)
points(r$baseMean , r$log2FoldChange, pch=21, lwd=1, col=r$col, bg=r$col, cex=1.2)
text(r$baseMean , r$log2FoldChange, labels=r$gene_name, adj=c(-0.15, -0.15), cex=1.2, col='black')
#
#  annotate group of genes
#
r<-gene_results(RES, IND)
points(r$baseMean , r$log2FoldChange, pch=21, lwd=1, col='black', bg=r$col, cex=0.9)
r<-gene_results(RES, REP)
points(r$baseMean , r$log2FoldChange, pch=21, lwd=1, col='black', bg=r$col, cex=0.9)
r<-gene_results(RES, SF)
points(r$baseMean , r$log2FoldChange, pch=21, lwd=1, col='black', bg=r$col, cex=0.9)
legend('topleft', legend=c('MYCN induced targets', 'MYCN repressed targets', 'Splicing factors'), bty='n', lty=0, lwd=0, pch=21, col='black', pt.bg=c(unique(IND$col), unique(REP$col), unique(SF$col)), pt.cex=1.8, pt.lwd=2, cex=1.2, x.intersp=-0.4, y.intersp=0.8)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/plotMA_genes_', sub(' vs ', '_', CND), '.svg'), width=14.5, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}



#  [genes] MYCN Tet-inducible 120h
#{{{

#  load DE results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_genes_+Tet 120h_ETOH 120h.RData')


#  single genes to annotate specifically
ANNOT<-data.frame(gene_id=c('ENSG00000134323.12'), gene_name=c('MYCN'), col='red3')


#  MYCN targets 
IND<-intersect(MYCN.IND, rownames(RES))
REP<-intersect(MYCN.REP, rownames(RES))
IND<-data.frame(gene_id=IND, gene_name=RES[ IND, 'gene_name'], col='royalblue3', cex=0.2)
REP<-data.frame(gene_id=REP, gene_name=RES[ REP, 'gene_name'], col='springgreen1', cex=0.2)


#  splice factors
SF<-SFU[ SFU$gene_id %in% intersect(SFU$gene_id, rownames(RES)), , drop=F]


# MA-plot with annotations
par(mar=c(5.0, 7.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4)
YLIM<-pretty(range(RES$log2FoldChange, na.rm=T))
YLIM<-c(YLIM[1], tail(YLIM, 1))
XLIM<-pretty(range(RES$baseMean, na.rm=T), 4)
XLIM<-c(0.1, tail(XLIM, 1))
options(scipen=+10)
plot(subset(RES, padj>=0.05)$baseMean, subset(RES, padj>=0.05)$log2FoldChange, type='p', log='x', xaxt='n', pch=21, xlab='', ylab='', main='', ylim=YLIM, xlim=XLIM, cex=0.4, lwd=0, col=adjustcolor('grey60', alpha.f=0.2), bg=adjustcolor('grey60', alpha.f=0.2))
axis(1, at=axTicks(1, usr=log10(XLIM), log=T))
points(subset(RES, padj<0.05)$baseMean, subset(RES, padj<0.05)$log2FoldChange, pch=21, lwd=1, cex=0.8, col=adjustcolor('red3', alpha.f=0.1), bg=adjustcolor('red3', alpha.f=0.5))
abline(h=0, lty=1, lwd=4, col='grey50')
#mtext(CND, side=3, line=0, padj=+1.5, cex=1.4)
mtext(expression(log[2]('fold change')), side=2, line=3, padj=-0.2, cex=2.4, las=3)
mtext('Mean expression', side=1, line=4, padj=-0.3, cex=2.4, las=1)
#  
#  annotate single genes
#
r<-gene_results(RES, ANNOT)
points(r$baseMean , r$log2FoldChange, pch=21, lwd=1, col=r$col, bg=r$col, cex=1.2)
text(r$baseMean , r$log2FoldChange, labels=r$gene_name, adj=c(-0.15, -0.15), cex=1.2, col='black')
#
#  annotate group of genes
#
r<-gene_results(RES, IND)
points(r$baseMean , r$log2FoldChange, pch=21, lwd=1, col='black', bg=r$col, cex=0.9)
r<-gene_results(RES, REP)
points(r$baseMean , r$log2FoldChange, pch=21, lwd=1, col='black', bg=r$col, cex=0.9)
r<-gene_results(RES, SF)
points(r$baseMean , r$log2FoldChange, pch=21, lwd=1, col='black', bg=r$col, cex=0.9)
legend('topleft', legend=c('MYCN induced targets', 'MYCN repressed targets', 'Splicing factors'), bty='n', lty=0, lwd=0, pch=21, col='black', pt.bg=c(unique(IND$col), unique(REP$col), unique(SF$col)), pt.cex=1.8, pt.lwd=2, cex=1.2, x.intersp=-0.4, y.intersp=0.8)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/plotMA_genes_', sub(' vs ', '_', CND), '.svg'), width=15.5, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#}}}





##########################################
#
#
#  housekeeping genes and the MYCN problem
#
#
##########################################




#  [GSE80153, SHEP-21N time course after MYCN induction shutdown] identify stable housekeeping genes using ERCC normalized FPKMs
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  fix the ggplot2 theme
#{{{
theme_set(theme_bw(base_size=35) + theme(panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_blank(),
    plot.margin=margin(t=1.0, b=0.1, l=0.5, r=1.0, unit='cm'),
    axis.line=element_line(color='black', size=0.5),
    axis.text=element_text(color='black', size=40, family='Arial', face='plain'), 
    axis.text.x=element_text(margin=margin(t=0.2, b=0, unit='cm')),
    axis.title=element_text(color='black', size=45, family='Arial', face='plain'), 
    axis.title.x=element_text(color='black', margin=margin(t=20, b=1)),
    axis.title.y=element_text(color='black', margin=margin(r=20, l=1)),
    axis.ticks=element_line(color='black', size=0.5), 
    axis.ticks.length=unit(0.5,  'cm'),
    legend.background=element_blank(),
    legend.justification=c(1, 1), 
    legend.position=c(0.15, 1.10), 
    legend.margin=margin(),
    legend.key=element_blank(),
    legend.title=element_text(size=24, family='Arial', face='plain', margin=margin(b=0.5, unit='cm')),
    legend.text=element_text(size=24, family='Arial', face='plain')
))
#}}}


#  load housekeeping genes
#  downloaded from: https://doi.org/10.1016/j.tig.2013.05.010
#  lifted over to GENCODE v30 with the help of HUGO as well.
hk<-fread('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/housekeeping.genes.gencode.v30.tsv', header=F, data.table=F)[, 1]


#  download ERCC normalized FPKMs:
#
#      wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE80153&format=file&file=GSE80153%5FSHEP21%5Fall%5Ffpkm%5Fexprs%5Fnorm%2Etxt%2Egz' -O - | gzip -d > GSE80153_SHEP21_all_fpkm_exprs_norm.txt


#  load them up (this is a data.frame with row.names, fread() will complain)
f.all<-f<-read.table('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/GSE80153/GSE80153_SHEP21_all_fpkm_exprs_norm.txt', sep='\t', header=T, row.names=1)
colnames(f)<-sub('SHEP21_', '', colnames(f))


#  percentage of housekeeping genes among all genes
round(length(intersect(hk, rownames(f)))/nrow(f)*100, 2)  #  15.13


#  keep only housekeeping genes
#  log10-transform
f<-log10(f[ rownames(f) %in% hk, ])
rm(hk)


#  remove 0h time point with only two replicates
f<-f[, grep('0H', colnames(f), invert=T)]
colnames(f)<-sub('_rep.*$', '', colnames(f))
timePoints<-unique(colnames(f))


#  do pairwise t-tests among successive time points
#  do not adjust for multiple hypothesis testing to be even more conservative since we are looking now at the non-DE
tt<-list()
for(n in seq(1, length(timePoints), 2)){
    p1<-colnames(f) %in% timePoints[n]
    p2<-colnames(f) %in% timePoints[n+1]
    #tt[[paste(timePoints[c(n, n+1)], sep='', collapse='_')]]<-p.adjust(apply(f, 1, function(x){ t.test(x=x[ p1 ], y=x[ p2 ], alternative='two.sided')$p.value }), 'fdr')
    tt[[paste(timePoints[c(n, n+1)], sep='', collapse='_')]]<-apply(f, 1, function(x){ t.test(x=x[ p1 ], y=x[ p2 ], alternative='two.sided')$p.value })
}
tt<-do.call(cbind, tt)
rm(n, p1, p2)


#  identify very well-expressed and not DE housekeeping genes
stopifnot( all.equal( rownames(tt), rownames(f) ) ) 
hk<-f[ rowMeans(f)>=2.7 & apply(tt, 1, function(x){ all(x>0.1) }), ]


#  save mean expression across all time points/replicates
wb<-createWorkbook()
addWorksheet(wb, 'stable housekeeping genes')
writeDataTable(wb, 'stable housekeeping genes', data.frame(gene_name=rownames(hk), mean=rowMeans(f.all[rownames(hk), ])))
setColWidths(wb, sheet=1, cols=seq_len(ncol(f.all)), widths='auto')
saveWorkbook(wb, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/GSE80153/top_stable_housekeepers.xlsx', overwrite=T)


#  recycle
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  stripchart (legend might look like shit on device but ok on SVG)
B<-data.table(gene_name=rownames(hk), hk)
B<-data.frame(melt(B, id.vars='gene_name', variable.name='time', value.name='normalized_expression'))
B<-B[ order(B$normalized_expression, decreasing=T), ]
B$time<-factor(B$time, levels=unique(colnames(hk)))
B$gene_name<-factor(B$gene_name, levels=unique(B$gene_name))
S.COL<-setNames(colorRampPalette(c('#00008B', '#006400', '#104E8B', '#2F4F4F', '#458B74', '#556B2F', '#8B1A1A', '#E9967A', '#FF8C00'))(length(levels(B$time))), levels(B$time))
#YTICK<-pretty(range(B$normalized_expression, na.rm=T), 5)
YTICK<-pretty(c(2.4, 4.0), 4)
XHIGH<-seq(2, length(unique(B$gene_name)), 2)
ggplot(B, aes(x=gene_name, y=normalized_expression, color=time, group=gene_name)) + 
    geom_jitter(height=0, width=0.2, size=6) +
    geom_rect(data=data.frame(xmin=XHIGH-0.5, xmax=XHIGH+0.5), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), 
              color='white', fill='grey39', alpha=0.2, show.legend=F, inherit.aes=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    labs(x='', y=expression(log[10]('Norm.Expr')), color='') + 
    scale_color_manual(values=S.COL) +
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_text(angle=90, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=0.5, hjust=1),
          plot.margin=margin(t=1.0, b=-2.0, l=0.5, r=1.0, unit='cm'),
          legend.justification=c(1, 1), 
          legend.position=c(1.01, 1.10)
    )
ggsave('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/GSE80153/stripchart_top_stable_housekeepers.svg', device=svg, width=26, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/GSE80153/top_stable_housekeepers.xlsx



#  [tumors] is there a global shift in housekeeping genes in MNA tumors compared to HR_nMNA?
#           is there a significant increase in TPMs of the genes producing circRNAs in the MNA tumors compared to HR_nMNA?
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(vsn)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_ecdfs.R')


#  load pre-prepared counts etc. for all samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')


#  load DE results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_MNA_HR_nMNA.RData')


#  load housekeeping genes
#  downloaded from: https://doi.org/10.1016/j.tig.2013.05.010
#  lifted over to GENCODE v30 with the help of HUGO as well.
hk<-fread('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/housekeeping.genes.gencode.v30.tsv', header=F, data.table=F)[, 1]


#  percentage of housekeeping genes among all genes
round(length(hk)/nrow(gns.tpm)*100, 2)  #  6.49


#  [housekeeping genes] pick them from the DE results
hk.r<-subset(RES, gene_name %in% hk)


#  [circRNA producing genes] pick their TPMs
gc.tpm<-gns.tpm[ rownames(gns.tpm) %in% hsa$gene_name[ match(rownames(crs.cnt), hsa$gene_id) ], ]


#  recycle
x11(width=16, height=16, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  housekeeping genes
#{{{

#  MA plot annotating the housekeeping genes
YLIM<-pretty(range(RES$log2FoldChange, na.rm=T))
YLIM<-c(YLIM[1], tail(YLIM, 1))
XLIM<-pretty(range(RES$baseMean, na.rm=T), 5)
XLIM<-c(0.1, tail(XLIM, 1))
par(mar=c(5.0, 7.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2)
options(scipen=+20)
DESeq2::plotMA(RES, xlab='', ylab='', main='', ylim=YLIM, xlim=XLIM, cex=1.2, colSig='red3')
points(hk.r$baseMean, hk.r$log2FoldChange, pch=19, col='springgreen', cex=0.5)
abline(h=median(hk.r$log2FoldChange), lty=1, lwd=4, col='springgreen')
mtext(CND, side=3, line=0, padj=+1.5, cex=1.4)
mtext(expression(log[2]('fold change')), side=2, line=3, padj=-0.2, cex=2.4, las=3)
mtext('Mean expression', side=1, line=4, padj=-0.3, cex=2.4, las=1)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/plotMA_genes_annotating_housekeeping_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(YLIM,XLIM)


#  boxplot of log10(1+TPMs)
par(mar=c(2.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-log10(1+gns.tpm[ rownames(gns.tpm) %in% hk.r$gene_name, meta[ risk_group %in% c('HR_nMNA', 'MNA'), bid] ])
colnames(B)<-meta[ match(colnames(B), bid), risk_group ]
B<-lapply(split(split(B, col(B)), colnames(B)), unlist)
#
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='two.sided')$p.value  #  4.213897874e-18
#
B.cl<-setNames(meta[ risk_group %in% c('HR_nMNA', 'MNA'), unique(col)], meta[ risk_group %in% c('HR_nMNA', 'MNA'), unique(risk_group)])
YTICK<-pretty(range(sapply(B, range)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
mtext(text='housekeeping genes', side=3, line=-1, las=1, padj=-0.1, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_genes_housekeeping_MNA_vs_HR_nMNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  ecdf
#
my_ecdfs(B, B.cl, XLIM=c(0.0, 3.0), XLAB=expression(log[10](1+'TPM')), MAIN='housekeeping genes', LTY=1, LWD=12, LEGEND='topleft', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_genes_housekeeping_MNA_vs_HR_nMNA.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)

#}}}


#  [circRNA-producing genes] boxplot of log10(1+TPMs)
par(mar=c(2.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-log10(1+gc.tpm[, meta[ risk_group %in% c('HR_nMNA', 'MNA'), bid] ])
colnames(B)<-meta[ match(colnames(B), bid), risk_group ]
B<-lapply(split(split(B, col(B)), colnames(B)), unlist)
#
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value  #  7.804481886e-77
#
B.cl<-setNames(meta[ risk_group %in% c('HR_nMNA', 'MNA'), unique(col)], meta[ risk_group %in% c('HR_nMNA', 'MNA'), unique(risk_group)])
YTICK<-pretty(range(sapply(B, range)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
mtext(text='host genes of circRNAs', side=3, line=-1, las=1, padj=-0.1, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_genes_of_circRNAs_MNA_vs_HR_nMNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  ecdf
#
my_ecdfs(B, B.cl, XLIM=c(0.0, 3.0), XLAB=expression(log[10](1+'TPM')), MAIN='host genes of circRNAs', LTY=1, LWD=12, LEGEND='topleft', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_genes_of_circRNAs_MNA_vs_HR_nMNA.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)

#}}}




#################
#
#
#  fooling around
#  obsolete code
#  scratch code
#
#
#################




#  PAGODA?
#{{{
library(org.Hs.eg.db)
library(GO.db)
library(Seurat)
library(scde)   #  install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")


#  turn NAs to zeros (there is no difference in the results when NAs are kept)
n.<-dcast(nb.circ, ssid ~ gene_name, value.var='jc_count', fun.aggregate=sum)    #  CDR1 is associated with two different gene_id...we aggregate
n.[is.na(n.)]<-0


#  convert to matrix and transpose
n.<-t(as.matrix(data.frame(n.[, -1], row.names=n.$ssid, check.names=F)))


#  min.lib.size : 50 circRNAs per sample
#     min.reads : 10 reads per circRNA (anyway they are filtered at this level)
#  min.detected : circRNAs in at least 5 samples (~8% frequency)
dge<-clean.counts(n., min.lib.size=50, min.reads=10, min.detected=5)


#  k nearest neighbor error model
#
#                    k : 5 number of nearest neighbors used for fitting (similar to frequency of circRNA occurrence)
#  min.count.threshold : cutoff for dropout events, in our case all circRNAs with counts>=10 are kept so it does not matter
#        min.nonfailed : 5 samples should support the circRNA in order to be considered
#     min.size.entries : 100 minimum number of circRNAs to use for model fitting
knn<-knn.error.models( dge, k=5, n.cores=20, min.count.threshold=1, min.nonfailed=5, max.model.plots=10, min.size.entries=100)
#
#          corr.a : slope of the correlated component fit 
#          corr.b : intersept of the correlated component fit
#  conc.a, conc.b : concomitant or associated fits
#      corr.theta : NB over-dispersion
#          fail.r : background Poisson rate (fixed)



#  normalize variance relative to circRNA-wide expectation with outlier trimming
varinfo<-pagoda.varnorm(knn, counts=dge, trim=3/ncol(dge), max.adj.var=5, n.cores=20, plot=T)
varinfo<-pagoda.subtract.aspect(varinfo, colSums(dge[, rownames(knn)]>0))


#  translate gene names to ids
ids<-unlist(lapply(mget(rownames(dge), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids<-names(ids)
names(rids) <- ids
go.env<-eapply(org.Hs.egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
go.env<-clean.gos(go.env, min.size=5)
desc<-unlist(lapply(mget(names(go.env), GOTERM, ifnotfound = NA), function(x) if(is.logical(x)) { return('') } else { slot(x, 'Term')}))
names(go.env)<-paste(names(go.env), desc) 
go.env<-list2env(go.env)



#  get overdispersed GO gene sets 
pwpca<-pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=20, n.internal.shuffles=10)


#  de novo gene sets
clpca<-pagoda.gene.clusters(varinfo, trim=7.1/ncol(varinfo$mat), n.clusters=150, n.cores=20, plot=T)


#  top overdispersed GO gene sets
tam<-pagoda.top.aspects(pwpca, z.score = 1.96)


#  cluster samples based on top overdispersed gene sets
hc<-pagoda.cluster.cells(tam, varinfo)


#  reduce GO term redundancy
tamr<-pagoda.reduce.loading.redundancy(tam, pwpca, n.cores=20)
tamr2<-pagoda.reduce.redundancy(tamr, distance.threshold=0.5, plot=T, cell.clustering=hc, labRow=NA, labCol=NA, box=T, margins=c(0.5, 0.5), trim=0)


#  heatmap
pagoda.view.aspects(tamr2, cell.clustering=hc, box=T, labCol=NA, margins=c(0.5, 60), col.cols=nb.meta[ match(hc$labels, ssid), col])




#}}}



#  ZINB-WaVE?
#{{{
library(zinbwave)
library(magrittr)


#  CDR1 is associated with two different gene_id so we aggregate which also turns all NA to zero
n.<-dcast(nb.circ, ssid ~ gene_name, value.var='jc_count', fun.aggregate=sum)


#  convert to matrix and transpose so that samples are columns and genes are rows
n.<-t(as.matrix(data.frame(n.[, -1], row.names=n.$ssid, check.names=F)))


#  keep 100 most variable circRNAs
# n. %>% log1p %>% rowVars -> vars
# names(vars)<-rownames(n.)
# vars<-sort(vars, decreasing=T)
# n.<-n.[ names(vars[1:100]), ]


#  keep only circRNAs with frequency of occurrence > 0.6
n.<-n.[ apply(n., 1, function(x){ sum(x>0) })/ncol(n.) > 0.6, ]


#  convert to SummarizedExperiment object
se<-SummarizedExperiment( assays=list(counts=n.),
                          colData=DataFrame(risk=nb.meta[ match(colnames(n.), ssid), risk_group], 
                                            nreads=nb.meta[ match(colnames(n.), ssid), nreads], 
                                            gc=nb.meta[ match(colnames(n.), ssid), gc], row.names=colnames(n.)))


#  fit the zero-inflated negative binomial model:
# 
#      n samples, J features with Y_{ij} the count of feature j for sample i modelled with a zero-inflated negative binomial distribution
#      with parameters \mu_{ij}, \theta_{ij} and \pi_{ij} regressed as follows:
#
#             ln( \mu_{ij} ) = (X\beta_\mu + (V\gamma_\mu)^T + W\alpha_\mu + O_\mu )_{ij}
#          logit( \pi_{ij} ) = (X\beta_\pi + (V\gamma_\pi)^T + W\alpha_\pi + O_\pi )_{ij}
#            ln(\theta_{ij}) = \zeta_j
#
#      X is a known nxM matrix corresponding to M cell-level covariates and \beta=(\beta_\mu,\beta_\pi) its associated MxJ matrices of 
#      regression parameters. X can typically include covariates that induce variation of interest, such as cell types, or covariates that 
#      induce unwanted variation, such as batch or quality control (QC) measures. By default, it includes only a constant column of ones to 
#      account for gene-specific intercepts.
#      
#      V is a known JxL matrix corresponding to J gene-level covariates, such as gene length or GC-content, and \gamma=(\gamma_\mu,\gamma_\pi) its 
#      associated Lxn matrices of regression parameters. By default, V only includes a constant column of ones to account for cell-specific intercepts, 
#      such as size factors representing differences in library sizes.
#      
#      W is an unobserved nxK matrix corresponding to K unknown cell-level covariates, which could be of unwanted variation or of interest 
#      (such as cell type), and \alpha=(\alpha_\mu,\alpha_\pi) its associated KxJ matrices of regression parameters. 
#
#      O_\mu and O_\pi are known nxJ matrices of offsets. 
#
zinb<-zinbFit(se, K=2, X='~gc+nreads+risk', epsilon=1000, verbose=T)    #  K=2 infers two latent variables after known variation is accounted for


#  get the estimated unobserved variation
W<-getW(zinb)
colnames(W)<-paste0('W', seq_len(ncol(W)))
data.frame( W, risk=colData(se)$risk, gc=colData(se)$gc, nreads=colData(se)$nreads ) %>% ggplot( aes(W1, W2, colour=risk)) + geom_point(size=10) + scale_color_brewer(type='qual', palette='Set1') + theme_classic() 

#}}}



#  tSNE 
#{{{
library(Rtsne)


#  functions
#{{{

source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/read.cef.R')

#}}}


#  [circRNAs] tSNE based on expression
#{{{

#  variance-stabilized and normalized circRNA counts
#n.<-apply(sweep(nb.circ, 2, nb.sf, '/'), 2, ceiling)
#n.<-nb.circ.vs


#  variance-stabilized and normalized circRNA counts of circRNAs expressed in at least 65% of samples
n.<-nb.circ.vs[ apply(nb.circ, 1, function(z){ sum(z!=0) } )/ncol(nb.circ)>=0.65, ]


#  PCA
p<-prcomp(t(n.), center=T, scale.=F)    #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
rm(p)


#  gene raw counts normalized by library-sizes and rounded up
# n.<-apply(sweep(nb.cnt, 2, nb.sf, '/'), 2, ceiling)
# n.<-nb.cnt.vs


#  recycle
ex<-data.frame(Risk=factor(nb.meta[ match(colnames(n.), ssid), risk_group], exclude=F), row.names=colnames(n.))
cl<-setNames( nb.meta[, unique(col)], nb.meta[, unique(risk_group)] )
ex$col<-cl[ match( ex$Risk, names(cl) ) ]


#  run t-SNE with given complexity (column dimension reduced, row dimension preserved)
tsne.pca<-Rtsne(pca, check_duplicates=F, pca=F, pca_center=F, pca_scale=F, verbose=T, theta=0.0, dims=2, max_iter=3e3, perplexity=12)
tsne.circs<-Rtsne(n., check_duplicates=F, pca=T, pca_center=T, pca_scale=F, verbose=T, theta=0.0, dims=2, max_iter=3e3, perplexity=12)
tsne<-Rtsne(t(n.), check_duplicates=F, pca=T, pca_center=T, pca_scale=F, verbose=T, theta=0.0, dims=2, max_iter=3e3, perplexity=12)


#  save all
save(tsne, tsne.circs, tsne.pca, pca, ve, n., file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/tSNE_circRNAs_CIRI2_geq_65p_samples_vsc_12.RData')


#  load data back
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/tSNE_circRNAs_CIRI2_geq_65p_samples_vsc_12.RData')


#  plot tSNE based on counts reducing samples
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.0,6.5,0.5,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.8, cex.axis=1.6)
XLIM<-range(pretty(range(tsne.circs$Y[, 1]), 4))
YLIM<-range(pretty(range(tsne.circs$Y[, 2]), 4))
plot(x=tsne.circs$Y[, 1], y=tsne.circs$Y[, 2], pch=19, col='darkgrey', xlim=XLIM, ylim=YLIM, xlab='', ylab='', cex=2.2)
mtext('t-SNE y-axis', side=2, line=4, padj=-0.5, cex=1.8, las=3)
mtext('t-SNE x-axis', side=1, line=3, padj=+0.5, cex=1.8, las=1)



#  plot tSNE based on PCA reducing PCs
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.0,6.5,0.5,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.8, cex.axis=1.6)
XLIM<-range(pretty(range(tsne.pca$Y[, 1]), 3))
YLIM<-range(pretty(range(tsne.pca$Y[, 2]), 3))
plot(x=tsne.pca$Y[, 1], y=tsne.pca$Y[, 2], pch=19, col=ex$col, xlim=XLIM, ylim=YLIM, xlab='', ylab='', cex=2.2)
mtext('t-SNE y-axis', side=2, line=4, padj=-0.5, cex=1.8, las=3)
mtext('t-SNE x-axis', side=1, line=3, padj=+0.5, cex=1.8, las=1)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/tSNE_circRNAs_CIRI2_geq_65p_samples_vsc_pca_tumors.pdf', width=14, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/tSNE_circRNAs_CIRI2_geq_65p_samples_vsc_pca_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  plot tSNE based on the counts reducing circRNAs
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.0,6.5,0.5,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.8, cex.axis=1.6)
XLIM<-range(pretty(range(tsne$Y[, 1]), 4))
YLIM<-range(pretty(range(tsne$Y[, 2]), 4))
plot(x=tsne$Y[, 1], y=tsne$Y[, 2], pch=19, col=ex$col, xlim=XLIM, ylim=YLIM, xlab='', ylab='', cex=2.2)
mtext('t-SNE y-axis', side=2, line=4, padj=-0.5, cex=1.8, las=3)
mtext('t-SNE x-axis', side=1, line=3, padj=+0.5, cex=1.8, las=1)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/tSNE_circRNAs_CIRI2_geq_65p_samples_vsc_tumors.pdf', width=14, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/tSNE_circRNAs_CIRI2_geq_65p_samples_vsc_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  common clustered MNA samples with backSPIN?
bs<-read.cef('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_geq_65p_samples_biclustered_all_2.cef')
bs<-names(bs$annot_samples[ bs$annot_samples==3 ])
te<-as.data.frame(tsne$Y)
rownames(te)<-colnames(n.)
colnames(te)<-c('x', 'y')
te$risk_group<-nb.meta[ match(rownames(te), ssid), risk_group ]
b<-te[ identify(te$x, te$y, labels=rownames(te), cex=1.2), ]    #  identify the x, y boarders of the cluster
l<-te[ te$x <= max(b$x) & te$y<=max(b$y), ]                     #  identify the cluster
text(l$x, l$y, label='*', cex=2.0)                              #  annotate the points with stars just to be sure
te<-rownames(l[ l$risk_group %in% 'MNA', ])
length(intersect( bs, te ))                                     #  commonly clustered MNA samples in both backSPIN and tSNE
length(setdiff( bs, te ))                                       #  only in backsPIN
length(setdiff( te, bs ))                                       #  only in tSNE
rm(bs,te,b,l)

#}}}


#  [genes] tSNE based on expression
#{{{

#  variance-stabilized gene expression
#n.<-nb.cnt.vs


#  variance-stabilized and normalized gene counts of genes expressed in at least 65% of samples
n.<-nb.cnt.vs[ apply(nb.cnt, 1, function(z){ sum(z!=0) } )/ncol(nb.cnt)>=0.65, ]


#  gene raw counts normalized by library-sizes and rounded up
# n.<-apply(sweep(nb.cnt, 2, nb.sf, '/'), 2, ceiling)
# n.<-nb.cnt.vs


#  recycle
ex<-data.frame(Risk=factor(nb.meta[ match(colnames(n.), ssid), risk_group], exclude=F), row.names=colnames(n.))
cl<-setNames( nb.meta[, unique(col)], nb.meta[, unique(risk_group)] )
ex$col<-cl[ match( ex$Risk, names(cl) ) ]


#  run t-SNE with given complexity
tsne<-Rtsne(t(n.), check_duplicates=F, pca=T, pca_center=T, pca_scale=F, verbose=T, theta=0.0, dims=2, max_iter=3e3, perplexity=12)


#  save
save(tsne, n., file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/tSNE_genes_geq_65p_samples_vsc_12.RData')


#  load the results back
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/tSNE_genes_geq_65p_samples_vsc_12.RData')


#  plot
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.0,6.5,0.5,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.8, cex.axis=1.6)
XLIM<-range(pretty(range(tsne$Y[, 1]), 5))
YLIM<-range(pretty(range(tsne$Y[, 2]), 5))
plot(x=tsne$Y[, 1], y=tsne$Y[, 2], pch=19, col=ex$col, xlim=XLIM, ylim=YLIM, xlab='', ylab='', cex=2.2)
mtext('t-SNE y-axis', side=2, line=4, padj=-0.5, cex=1.8, las=3)
mtext('t-SNE x-axis', side=1, line=3, padj=+0.5, cex=1.8, las=1)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/tSNE_genes_geq_65p_samples_vsc_tumors.pdf', width=14, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/tSNE_genes_geq_65p_samples_vsc_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  common clustered MNA samples with backSPIN?
bs<-read.cef('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_geq_65p_samples_biclustered_all_2.cef')
bs<-names(bs$annot_samples[ bs$annot_samples==3 ])
te<-as.data.frame(tsne$Y)
rownames(te)<-colnames(n.)
colnames(te)<-c('x', 'y')
te$risk_group<-nb.meta[ match(rownames(te), ssid), risk_group ]
b<-te[ identify(te$x, te$y, labels=rownames(te), cex=1.2), ]    #  identify the x, y boarders of the cluster
l<-te[ te$x >= min(b$x) & te$y>=min(b$y), ]                     #  identify the cluster (should be at the top left)
text(l$x, l$y, label='*', cex=2.0)                              #  annotate the points with stars just to be sure
te<-rownames(l[ l$risk_group %in% 'MNA', ])
length(intersect( bs, te ))                                     #  commonly clustered MNA samples in both backSPIN and tSNE
length(setdiff( bs, te ))                                       #  only in backsPIN
length(setdiff( te, bs ))                                       #  only in tSNE
rm(bs,te,b,l)

#}}}

#}}}



#  gene expression vs circRNA expression comparison
#{{{

#  identify the common set of expressed genes with circRNAs
x<-nb.circ
y<-nb.cnt
rownames(y)<-hsa$gene_name[ match( rownames(y) , hsa$gene_id ) ]
common<-intersect( rownames(x), rownames(y) )
x<-x[ match(common, rownames(x)), ]
y<-y[ match(common, rownames(y)), ]
stopifnot( all.equal( rownames(x), rownames(y) ) )
rm(common)


#  R^2 sample by sample for genes vs circRNAs (the result when common zeros are excluded is almost identical)
r2<-diag(cor(x, y, method='pearson', use='complete.obs'))^2


#  take the mean across tumors and log10-transform and do linear regression
x<-log10(rowMeans(x))
y<-log10(rowMeans(y))
l<-lm( y ~ x )


#  scatterplot of mean counts across all samples
x11(width=12, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
XLIM<-pretty(range(x), 5)
YLIM<-pretty(range(y), 5)
par(mar=c(4.5,5,0.1,2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.8, cex.axis=1.6)
den<-col2rgb(densCols(x, y, colramp=colorRampPalette(c('black', 'white'))))[1, ]+1L
cls<-colorRampPalette(c('#000099', '#00FEFF', '#45FE4F', '#FCFF00', '#FF9400', '#FF3100'))(256)[den]
plot(x[ order(den) ], y[order(den)], type='p', pch=19, cex=0.8, col=cls[order(den)], xlim=c(XLIM[1], tail(XLIM,1)), ylim=c(YLIM[1], tail(YLIM,1)), xlab='', ylab='')
abline(l, lty=1, lwd=4, col='black')
mtext(bquote(R^2 == .(format(summary(l)$r.squared, digits=2))), side=3, line=-1, padj=+0.5, cex=1.4, xpd=NA)
mtext(expression(log[10]('mean counts')~'[genes]'), side=2, line=2, padj=-0.3, cex=1.8, las=3)
mtext(expression(log[10]('mean counts')~'[circRNAs]'), side=1, line=3, padj=+0.2, cex=1.8, las=1)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/scatterplot_counts_vsc_genes_circRNAs_CIRI2_tumors.pdf', width=14, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/scatterplot_counts_vsc_genes_circRNAs_CIRI2_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  barplot of R^2
x11(width=25, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(4.5, 5.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.4, cex.axis=1.4)
stopifnot( all.equal( names(r2), nb.meta$ssid ) )
B<-r2
B.cl<-setNames(nb.meta$col, nb.meta$risk_group)
YTICK<-pretty(c(0, 0.2), 4)
bp<-barplot(B, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, tail(YTICK,1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col='white', axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=1.4, cex.axis=1.4, add=T)
for(i in seq_along(B)){ 
    b<-B
    b[-i]<-NA    #  remove all values but the current 
    barplot(b, border='white', col=B.cl[i], axes=F, axisnames=F, beside=F, yaxt='n', add=T)
}
axis(2, at=YTICK, line=0, cex.axis=1.2)
mtext(text=names(B), side=1, line=0, at=bp, las=2, adj=1, cex=0.9, col=B.cl)
mtext(expression(R^2), side=2, line=3, padj=-0.6, las=0, cex=1.4)
legend(x=par('usr')[2]*0.001, y=par('usr')[4]*1.001, legend=unique(names(B.cl)), col=unique(B.cl), bty='n', lty=1, lwd=10, pch=NA, cex=1.2, xpd=T)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_R2_counts_vsc_genes_circRNAs_CIRI2.pdf', width=25, height=12, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_R2_counts_vsc_genes_circRNAs_CIRI2.svg', width=25, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}



#  backSPIN 
#
#  N.B.: it should be run on raw counts normalized by library-sizes to resemble single-cell sequencing requirements
#
#{{{

#  functions
#{{{

source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/plot_cv.R')
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/read.cef.R')
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/write.cef.R')


backspin_process<-function(CEF='', log2.dge=NULL, scale.dge=NULL, plot.only=NULL){
    bs<-read.cef(CEF)

    if(missing(log2.dge)){ log2.dge<-F }

    if(missing(scale.dge)){ scale.dge<-F }

    if(!is.null(plot.only)){
        bs$annot_samples<-bs$annot_samples[ bs$annot_samples==plot.only ]
        bs$annot_genes<-bs$annot_genes[ bs$annot_genes==plot.only ]
        bs$dge<-bs$dge[ rownames(bs$dge) %in% names(bs$annot_genes),  colnames(bs$dge) %in% names(bs$annot_samples) ]
    }


    if(log2.dge){
        n.<-log2(1+bs$dge)
    } else if (scale.dge){
        n.<-t(scale(t(bs$dge), center=T, scale=T))  #  scale genes across samples
    } else {
        n.<-bs$dge
    }


    an.col<-data.frame(Risk=factor(nb.meta[ match(colnames(n.), ssid), risk_group], exclude=F), backSPIN=factor(bs$annot_samples), row.names=colnames(n.))
    an.row<-data.frame(backSPIN=factor(bs$annot_genes), row.names=rownames(n.))
    cl.risk<-c(ST4S='seagreen4', LR='darkgreen', IMR='cornflowerblue', HR_nMNA='chocolate1', MNA='coral4')
    cl.backspin<-setNames(colorRampPalette(c('black', 'navy', 'dodgerblue4', 'steelblue', 'steelblue1'))(length(levels(an.row$backSPIN))), levels(an.row$backSPIN))
    #ph<-pheatmap(as.matrix(n.), color=colorRampPalette(brewer.pal(9,'Reds'))(20), border_color=NA, scale='none', 
    ph<-pheatmap(as.matrix(n.), color=viridis(20), border_color=NA, scale='none', 
            cluster_rows=F, cluster_cols=F,
            gaps_row=as.integer(na.omit(match(levels(an.row$backSPIN)[-1], an.row$backSPIN)-1)),
            gaps_col=as.integer(na.omit(match(levels(an.row$backSPIN)[-1], an.col$backSPIN)-1)),
            annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
            annotation_col=an.col, annotation_row=an.row, annotation_colors=list(Risk=cl.risk, backSPIN=cl.backspin),
            drop_levels=F, show_rownames=T, show_colnames=T, 
            #display_numbers=T, number_format='%.1f', number_color='grey39',
            fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
    

    return(list(dge=bs$dge, annot_samples=bs$annot_samples, annot_genes=bs$annot_genes, ph=ph))
}

#}}}


#  open and recycle
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#
#  circRNAs
#


#  biclustering of circRNAs 
#{{{

#  circRNAs raw counts normalized by library-sizes and rounded up
n.<-apply(sweep(nb.circ, 2, nb.sf, '/'), 2, ceiling)


#  [only for count data] plot log2(CV) versus log2(Mean) and the non-linear fit
#g<-hsa$gene_name[ match(rownames(nb), hsa$gene_id) ]
cl<-setNames(rep('darkgrey', nrow(n.)), rownames(n.))
#cl[ g ]<-'black'
plot_cv(n., with.r=T, pt.cex=0.5, col.genes=cl)


#  save as CEF
write.cef(n., hd=c('DGE', 'neuroblastoma circRNAs'), fl='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_all.cef')


#  [top varying circRNAs] run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_all.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_all_biclustered_500_2.cef -f 500 -v -d 2 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_all_biclustered_500_2.cef', scale.dge=T)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_circRNAs_CIRI2_all_biclustered_500_2.pdf', width=16, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_circRNAs_CIRI2_all_biclustered_500_2.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [all circRNAs] run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
#{{{

system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_all.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_all_biclustered_all_2.cef -v -d 2 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_all_biclustered_all_2.cef', scale.dge=T)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_circRNAs_CIRI2_all_biclustered_all_2.pdf', width=16, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_circRNAs_CIRI2_all_biclustered_all_2.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [do not run but keep around for reference: tightest MNA cluster] topGO enrichment analysis using our GO annotations
#{{{

#  load backSPIN clustering results using all circRNAs
circ.bs<-read.cef('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_all_biclustered_all_2.cef')


#  group of circRNAs that determine the MNA cluster
grp<-hsa$gene_id[ hsa$gene_name %in% names(circ.bs$annot_genes[ circ.bs$annot_genes==2 ]) ]


#  universe is all the circRNAs used for biclustering
uni<-hsa$gene_id[ hsa$gene_name %in% names(circ.bs$annot_genes) ]


#  define factor that identifies the genes of interest in the universe
allG<-setNames(rep(0, length(uni)), uni)
allG[ grp  ]<-1
stopifnot(  table(allG)[2] == length(grp) )
allG<-as.factor(allG)


#  molecular function GO enrichment
mf.topgo<-new('topGOdata', description='enrichment test', ontology='MF', allGenes=allG, annot=annFUN.gene2GO, gene2GO=MF.gid2go, nodeSize=10)
mf.test<-runTest(mf.topgo, algorithm='weight01', statistic='fisher')
mf<-data.table(GenTable(mf.topgo, p.value=mf.test, orderBy='p.value', topNodes=geneData(mf.test)['SigTerms'], numChar=120))
mf$p.value<-as.numeric(mf$p.value)
mf[ p.value<0.05 ]


#  biological process GO enrichment
bp.topgo<-new('topGOdata', description='enrichment test', ontology='BP', allGenes=allG, annot=annFUN.gene2GO, gene2GO=BP.gid2go, nodeSize=10)
bp.test<-runTest(bp.topgo, algorithm='weight01', statistic='fisher')
bp<-data.table(GenTable(bp.topgo, p.value=bp.test, orderBy='p.value', topNodes=geneData(bp.test)['SigTerms'], numChar=120))
bp$p.value<-as.numeric(bp$p.value)
bp[ p.value<0.05 ]


#  save
save(grp, uni, allG, mf, bp, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/topGO_circRNAs_CIRI2_all_biclustered_all_2_MNA_cluster.RData')


#  plot with p-value cutoff < 0.05
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/topGO_circRNAs_CIRI2_all_biclustered_all_2_MNA_cluster.RData')
x11(width=30, height=18, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(4.0,49.5,0.5,1.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.4, cex.axis=1.4)
b<-bp[ p.value<0.05, c('Term', 'Significant', 'p.value')][ order(-Significant) ]
b$col<-'grey39'
b$col[ grep('chrom|telomere', b$Term) ]<-'deepskyblue3'
MAX<-tail(pretty(b$Significant, 5), 1)
br<-barplot(b$Significant, border='white', col=b$col, horiz=T, axisnames=F, xlab='', ylab='', las=1, yaxt='n', xlim=c(0, MAX))
axis(2, at=br, labels=b$Term, line=0, tick=F, cex.axis=1.8)
mtext('Count', side=1, line=2, padj=+0.8, las=0, cex=1.8)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/topGO_circRNAs_CIRI2_all_biclustered_all_2_MNA_cluster.pdf', width=30, height=18, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/topGO_circRNAs_CIRI2_all_biclustered_all_2_MNA_cluster.svg', width=30, height=18, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b,MAX,br)


dev.off()

#}}}

#}}}


#  [circRNAs supported by at least 65% of samples] run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
#{{{

#  pick circRNAs supported by at least 65% of samples 
#  normalize by library-sizes and round up
n.<-nb.circ[ apply(nb.circ, 1, function(z){ sum(z!=0) } )/ncol(nb.circ)>=0.65, ]
n.<-apply(sweep(n., 2, nb.sf, '/'), 2, ceiling)
write.cef(n., hd=c('DGE', 'neuroblastoma circRNAs supported by at least 65% of samples'), fl='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_geq_65p_samples.cef')


#  [only for count data] plot log2(CV) versus log2(Mean) and the non-linear fit
#g<-hsa$gene_name[ match(rownames(nb), hsa$gene_id) ]
cl<-setNames(rep('darkgrey', nrow(n.)), rownames(n.))
#cl[ g ]<-'black'
plot_cv(n., with.r=T, pt.cex=0.5, col.genes=cl)


#  run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_geq_65p_samples.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_geq_65p_samples_biclustered_all_2.cef -v -d 2 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_geq_65p_samples_biclustered_all_2.cef', scale.dge=T)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_circRNAs_CIRI2_geq_65p_samples_biclustered_all_2.pdf', width=16, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_circRNAs_CIRI2_geq_65p_samples_biclustered_all_2.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [MNA cluster] topGO enrichment analysis using our GO annotations
#{{{

#  load backSPIN clustering results
circ.bs<-read.cef('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_geq_65p_samples_biclustered_all_2.cef')


#  group of circRNAs that determine the MNA cluster
grp<-hsa$gene_id[ hsa$gene_name %in% names(circ.bs$annot_genes[ circ.bs$annot_genes==3 ]) ]


#  universe is all the circRNAs used for biclustering
uni<-hsa$gene_id[ hsa$gene_name %in% names(circ.bs$annot_genes) ]


#  define factor that identifies the genes of interest in the universe
allG<-setNames(rep(0, length(uni)), uni)
allG[ grp  ]<-1
stopifnot(  table(allG)[2] == length(grp) )
allG<-as.factor(allG)


#  molecular function GO enrichment
mf.topgo<-new('topGOdata', description='enrichment test', ontology='MF', allGenes=allG, annot=annFUN.gene2GO, gene2GO=MF.gid2go, nodeSize=10)
mf.test<-runTest(mf.topgo, algorithm='weight01', statistic='fisher')
mf<-data.table(GenTable(mf.topgo, p.value=mf.test, orderBy='p.value', topNodes=geneData(mf.test)['SigTerms'], numChar=120))
mf$p.value<-as.numeric(mf$p.value)
mf[ p.value<0.05 ]


#  biological process GO enrichment
bp.topgo<-new('topGOdata', description='enrichment test', ontology='BP', allGenes=allG, annot=annFUN.gene2GO, gene2GO=BP.gid2go, nodeSize=10)
bp.test<-runTest(bp.topgo, algorithm='weight01', statistic='fisher')
bp<-data.table(GenTable(bp.topgo, p.value=bp.test, orderBy='p.value', topNodes=geneData(bp.test)['SigTerms'], numChar=120))
bp$p.value<-as.numeric(bp$p.value)
bp[ p.value<0.05 ]


#  save
save(grp, uni, allG, mf, bp, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/topGO_circRNAs_CIRI2_geq_65p_samples_biclustered_all_2_MNA_cluster.RData')


#  plot with p-value cutoff < 0.05
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/topGO_circRNAs_CIRI2_geq_65p_samples_biclustered_all_2_MNA_cluster.RData')
x11(width=30, height=18, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(4.0,43.5,0.5,1.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.4, cex.axis=1.4)
b<-bp[ p.value<0.05, c('Term', 'Significant', 'p.value')][ order(-Significant) ]
b$col<-'grey39'
b$col[ grep('chrom|telomere', b$Term) ]<-'deepskyblue3'
MAX<-tail(pretty(b$Significant, 5), 1)
br<-barplot(b$Significant, border='white', col=b$col, horiz=T, axisnames=F, xlab='', ylab='', las=1, yaxt='n', xlim=c(0, MAX))
axis(2, at=br, labels=b$Term, line=0, tick=F, cex.axis=1.8)
mtext('Count', side=1, line=2, padj=+0.8, las=0, cex=1.8)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/topGO_circRNAs_CIRI2_geq_65p_samples_biclustered_all_2_MNA_cluster.pdf', width=30, height=18, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/topGO_circRNAs_CIRI2_geq_65p_samples_biclustered_all_2_MNA_cluster.svg', width=30, height=18, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b,MAX,br)


dev.off()

#}}}


#  [MNA cluster] exploratory heatmap and PCA
#{{{

#  load backSPIN clustering results
circ.bs<-read.cef('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_geq_65p_samples_biclustered_all_2.cef')


#  isolate the DGE and the strong MNA cluster
n.<-circ.bs$dge
n.<-n.[ rownames(n.) %in% names(circ.bs$annot_genes[ circ.bs$annot_genes %in% 3 ]), ]


#  isolate circRNAs associated with the two groups of the last split
#n.<-n.[ rownames(n.) %in% names(circ.bs$annot_genes[ circ.bs$annot_genes %in% c(2, 3) ]), ]


#  recycle
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
ex<-data.frame(Risk=factor(nb.meta[ match(colnames(n.), ssid), risk_group], exclude=F), row.names=colnames(n.))
cl<-setNames(list(setNames( nb.meta[, unique(col)], nb.meta[, unique(risk_group)] )), colnames(ex))


#  Pearson correlations heatmap
d<-cor(n., method='pearson', use='pairwise.complete.obs')
ph<-pheatmap(as.matrix(d), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0, 1, length.out=11),  #  reduce the range to positives
        cluster_rows=F,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)


#  heatmap of distance matrix of Pearson correlations 
d<-as.dist(1-cor(n., method='pearson', use='pairwise.complete.obs'))
ph<-pheatmap(as.matrix(d), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0, 2, length.out=11),
        cluster_rows=F,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)


#  heatmap of distance matrix based on log2-transformed expression
d<-dist(t(log2(1+n.)), method='euclidean')
ph<-pheatmap(as.matrix(d), color=colorRampPalette(brewer.pal(9,'GnBu'))(10), border_color=NA, scale='none', 
        cluster_rows=F,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)


#  PCA of centered but not scaled variance-stabilized and log2-transformed counts
n.<-log2(nb.circ.vs[ match( rownames(n.), rownames(nb.circ.vs)), colnames(n.)])
p<-prcomp(t(n.), center=T, scale.=F)    #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$risk<-nb.meta[ match( rownames(pca), ssid ), risk_group ]
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk')], aes(PC1, PC2, color=risk)) + 
    geom_point(size=10) + 
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk', values=cl$Risk) + 
    scale_color_manual(name='risk', values=cl$Risk) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=34), axis.title.x=element_text(size=34), 
        axis.title.y=element_text(size=34), axis.text.y=element_text(size=34), 
        legend.text=element_text(size=22), legend.title=element_text(size=22, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.2, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.2, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 

#}}}


#  [MNA cluster] enrichment in MYCN targets
#{{{

#  load backSPIN clustering results
circ.bs<-read.cef('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_geq_65p_samples_biclustered_all_2.cef')


#  group of circRNAs that determine the MNA cluster and the rest
grp<-names(circ.bs$annot_genes[ circ.bs$annot_genes==3 ])
rst<-names(circ.bs$annot_genes[ circ.bs$annot_genes!=3 ])
stopifnot( length(grp)+length(rst)==length(circ.bs$annot_genes) )


#  MYCN targets 
#  separate to induced and repressed
load('/data/annotation/GRCh38/MYCN_targets.RData')
IND<-MYCN[ type %in% 'induced', gene_name]
REP<-MYCN[ type %in% 'repressed', gene_name]
stopifnot( length(IND)+length(REP)==nrow(MYCN) )
rm(MYCN)


#  [induced targets] Fisher's exact test:
#
#                   |  in-MNA-cluster  |  out-of-MNA-cluster  |
#  -----------------|------------------|----------------------|
#      MYCN target  |          x1      |      y1              |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |          x2      |      y2              | 
fisher.test(data.frame('in'=c(length(intersect(grp, IND)), 
                              length(setdiff(grp, IND))), 
                       'out'=c(length(intersect(rst, IND)), 
                               length(setdiff(rst, IND)))), alternative='greater')$p.value
#
#  => 0.0063


#  there are no repressed MYCN targets in the MNA group
intersect(grp, REP)

#}}}

#}}}

#}}}


#  [do not run but keep around for reference] biclustering of MSigDB neuroblastoma-specific circRNAs
#{{{

#  circRNAs raw counts normalized by library-sizes and rounded up
g<-hsa$gene_name[ match(rownames(nb), hsa$gene_id) ]
n.<-apply(sweep(nb.circ[ rownames(nb.circ) %in% g, ], 2, nb.sf, '/'), 2, ceiling)
rm(g)


#  [only for count data] plot log2(CV) versus log2(Mean) and the non-linear fit
g<-hsa$gene_name[ match(rownames(nb), hsa$gene_id) ]
cl<-setNames(rep('darkgrey', nrow(n.)), rownames(n.))
cl[ g ]<-'black'
plot_cv(n., with.r=T, pt.cex=0.5, col.genes=cl)
rm(g)


#  save as CEF
write.cef(n., hd=c('DGE', 'neuroblastoma circRNAs selected by MSigDB analysis'), fl='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_MSigDB.cef')


#  run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_MSigDB.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_MSigDB_biclustered_all_2.cef -v -d 2 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_MSigDB_biclustered_all_2.cef', scale.dge=T)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_circRNAs_CIRI2_MSigDB_biclustered_all_2.pdf', width=16, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_circRNAs_CIRI2_MSigDB_biclustered_all_2.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [do not run but keep around for reference: MNA vs HR_nMNA DE results] biclustering based on specific circRNAs associated with the DE analysis
#{{{

#  load gene DE between MNA vs HR_nMNA
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/DESeq2_genes_MNA_HR_nMNA.RData')
x<-subset(RES, padj<0.05)
rm(CND, DDS, PCA, VE, VSC, RES)


#  load circRNAs DE between MNA vs HR_nMNA
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/DESeq2_circRNAs_CIRI2_MNA_HR_nMNA.RData')
y<-subset(RES, padj<0.05)
rm(CND, DDS, PCA, VE, VSC, RES)


#  genes common between them
common<-intersect( rownames(x), rownames(y) )
common<-data.frame( gene=x[common, 'baseMean'], gene_log2FC=x[common, 'log2FoldChange'], gene_padj=x[common, 'padj'],
                    circ=y[common, 'baseMean'], circ_log2FC=y[common, 'log2FoldChange'], circ_padj=y[common, 'padj'],
                    gene_name=y[common, 'gene_name'], row.names=common)


#  circRNAs raw counts normalized by library-sizes and rounded up
#n.<-apply(sweep(nb.circ[ common$gene_name, ], 2, nb.sf, '/'), 2, ceiling)                      #  common DE genes + circRNAs
n.<-apply(sweep(nb.circ[ rownames(nb.circ) %in% x$gene_name, ], 2, nb.sf, '/'), 2, ceiling)    #  circRNAs of DE genes (best results)
#n.<-apply(sweep(nb.circ[ rownames(nb.circ) %in% y$gene_name, ], 2, nb.sf, '/'), 2, ceiling)    #  DE circRNAs 


#  [only for count data] plot log2(CV) versus log2(Mean) and the non-linear fit
cl<-setNames(rep('darkgrey', nrow(n.)), rownames(n.))
cl[ names(cl) %in% rownames(nb) ]<-'black'
plot_cv(n., with.r=F, pt.cex=0.5, col.genes=cl)


#  save as CEF
write.cef(n., hd=c('DGE', 'MNA vs HS_nMNA specific circRNAs'), fl='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_specific_MNA_HR_nMNA.cef')


#  run backSPIN with maximum 1 split since we are dealing with DE results between two conditions
system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_specific_MNA_HR_nMNA.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_specific_MNA_HR_nMNA_biclustered_all_1.cef -v -d 1 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/circRNAs_CIRI2_specific_MNA_HR_nMNA_biclustered_all_1.cef', scale.dge=T)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_circRNAs_CIRI2_specific_MNA_HR_nMNA_biclustered_all_1.pdf', width=16, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_circRNAs_CIRI2_specific_MNA_HR_nMNA_biclustered_all_1.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#
#  genes
#


#  biclustering of genes
#{{{

#  gene raw counts normalized by library-sizes and rounded up
n.<-apply(sweep(nb.cnt, 2, nb.sf, '/'), 2, ceiling)


#  [only for count data] plot log2(CV) versus log2(Mean) and the non-linear fit
cl<-setNames(rep('darkgrey', nrow(n.)), rownames(n.))
#cl[ names(cl) %in% rownames(nb) ]<-'black'
plot_cv(n., with.r=T, pt.cex=0.5, col.genes=cl)


#  save as CEF
write.cef(n., hd=c('DGE', 'neuroblastoma genes'), fl='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_all.cef')


#  [top 5000 varying genes] run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_all.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_all_biclustered_5000_2.cef -f 5000 -v -d 2 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_all_biclustered_5000_2.cef', scale.dge=T)
TOP<-5000    #  save
dev.print(device=pdf, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_genes_all_biclustered_', TOP, '_2.pdf'), width=16, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_genes_all_biclustered_', TOP,'_2.svg'), width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [top 1000 varying genes] run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_all.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_all_biclustered_1000_2.cef -f 1000 -v -d 2 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_all_biclustered_1000_2.cef', scale.dge=T)
TOP<-1000    #  save
dev.print(device=pdf, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_genes_all_biclustered_', TOP, '_2.pdf'), width=16, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_genes_all_biclustered_', TOP,'_2.svg'), width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [all genes, takes ~ 4 days to finish] run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
#{{{

system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_all.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_all_biclustered_all_2.cef -v -d 2 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_all_biclustered_all_2.cef', scale.dge=T)
#  save only PDF, SVG of 1GB size needs >130GB of memory to be converted to PNG by inkscape
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_genes_all_biclustered_all_2.pdf', width=16, height=14, bg='white', colormodel='srgb', pointsize=20, useDingbats=F, family='Arial')  #  sRGB color model


#  [tightest MNA cluster] topGO enrichment analysis using our GO annotations
#{{{

#  load backSPIN clustering results 
gns.bs<-read.cef('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_all_biclustered_all_2.cef')


#  group of genes that determine the MNA cluster
grp<-names(gns.bs$annot_genes[ gns.bs$annot_genes==2 ])


#  universe is all the genes used for biclustering
uni<-names(gns.bs$annot_genes)


#  define factor that identifies the genes of interest in the universe
allG<-setNames(rep(0, length(uni)), uni)
allG[ grp  ]<-1
stopifnot(  table(allG)[2] == length(grp) )
allG<-as.factor(allG)


#  molecular function GO enrichment
mf.topgo<-new('topGOdata', description='enrichment test', ontology='MF', allGenes=allG, annot=annFUN.gene2GO, gene2GO=MF.gid2go, nodeSize=10)
mf.test<-runTest(mf.topgo, algorithm='weight01', statistic='fisher')
mf<-data.table(GenTable(mf.topgo, p.value=mf.test, orderBy='p.value', topNodes=geneData(mf.test)['SigTerms'], numChar=120))
mf$p.value<-as.numeric(mf$p.value)
mf[ p.value<0.05 ]


#  biological process GO enrichment
bp.topgo<-new('topGOdata', description='enrichment test', ontology='BP', allGenes=allG, annot=annFUN.gene2GO, gene2GO=BP.gid2go, nodeSize=10)
bp.test<-runTest(bp.topgo, algorithm='weight01', statistic='fisher')
bp<-data.table(GenTable(bp.topgo, p.value=bp.test, orderBy='p.value', topNodes=geneData(bp.test)['SigTerms'], numChar=120))
bp$p.value<-as.numeric(bp$p.value)
bp[ p.value<0.05 ]


#  save
save(grp, uni, allG, mf, bp, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/topGO_genes_all_biclustered_all_2_MNA_cluster.RData')


#  plot with p-value cutoff < 1e-5
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/topGO_genes_all_biclustered_all_2_MNA_cluster.RData')
x11(width=30, height=30, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(4.0, 45.5, 0.5,1.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.4, cex.axis=1.4)
b<-bp[ p.value<1e-5, c('Term', 'Significant', 'p.value')][ order(-Significant) ]
b$col<-'grey39'
b$col[ grep('chrom|telomere', b$Term) ]<-'deepskyblue3'
MAX<-tail(pretty(b$Significant, 5), 1)
br<-barplot(b$Significant, border='white', col=b$col, horiz=T, axisnames=F, xlab='', ylab='', las=1, yaxt='n', xlim=c(0, MAX))
axis(2, at=br, labels=b$Term, line=0, tick=F, cex.axis=1.8)
mtext('Count', side=1, line=2, padj=+0.8, las=0, cex=1.8)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/topGO_genes_all_biclustered_all_2_MNA_cluster.pdf', width=30, height=30, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/topGO_genes_all_biclustered_all_2_MNA_cluster.svg', width=30, height=30, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b,MAX,br)


dev.off()

#}}}

#}}}


#  [genes expressed in at least 65% of samples, takes ~ 2 days to finish] run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
#{{{

#  pick genes that are expressed in at least 65% of the samples 
#  normalize by library-sizes and round up
n.<-nb.cnt[ apply(nb.cnt, 1, function(z){ sum(z!=0) } )/ncol(nb.cnt)>=0.65, ]
n.<-apply(sweep(n., 2, nb.sf, '/'), 2, ceiling)
write.cef(n., hd=c('DGE', 'neuroblastoma genes expressed in at least 65% of the samples'), fl='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_geq_65p_samples.cef')


#  run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_geq_65p_samples.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_geq_65p_samples_biclustered_all_2.cef -v -d 2 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_geq_65p_samples_biclustered_all_2.cef', scale.dge=T)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_genes_geq_65p_samples_biclustered_all_2.pdf', width=16, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')  #  only PDF


#  [MNA cluster] topGO enrichment analysis using our GO annotations
#{{{

#  load backSPIN clustering results 
gns.bs<-read.cef('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_geq_65p_samples_biclustered_all_2.cef')


#  group of genes that determine the main MNA cluster
grp<-names(gns.bs$annot_genes[ gns.bs$annot_genes==2 ])


#  universe is all the genes used for biclustering
uni<-names(gns.bs$annot_genes)


#  define factor that identifies the genes of interest in the universe
allG<-setNames(rep(0, length(uni)), uni)
allG[ grp  ]<-1
stopifnot(  table(allG)[2] == length(grp) )
allG<-as.factor(allG)


#  molecular function GO enrichment
mf.topgo<-new('topGOdata', description='enrichment test', ontology='MF', allGenes=allG, annot=annFUN.gene2GO, gene2GO=MF.gid2go, nodeSize=10)
mf.test<-runTest(mf.topgo, algorithm='weight01', statistic='fisher')
mf<-data.table(GenTable(mf.topgo, p.value=mf.test, orderBy='p.value', topNodes=geneData(mf.test)['SigTerms'], numChar=120))
mf$p.value<-as.numeric(mf$p.value)
mf[ p.value<0.05 ]


#  biological process GO enrichment
bp.topgo<-new('topGOdata', description='enrichment test', ontology='BP', allGenes=allG, annot=annFUN.gene2GO, gene2GO=BP.gid2go, nodeSize=10)
bp.test<-runTest(bp.topgo, algorithm='weight01', statistic='fisher')
bp<-data.table(GenTable(bp.topgo, p.value=bp.test, orderBy='p.value', topNodes=geneData(bp.test)['SigTerms'], numChar=120))
bp$p.value<-as.numeric(bp$p.value)
bp[ p.value<0.05 ]


#  save
save(grp, uni, allG, mf, bp, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/topGO_genes_geq_65p_samples_biclustered_all_2_MNA_cluster.RData')


#  plot with p-value cutoff < 0.001
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/topGO_genes_geq_65p_samples_biclustered_all_2_MNA_cluster.RData')
x11(width=30, height=30, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(4.0,51.5,0.5,1.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.4, cex.axis=1.4)
b<-bp[ p.value<0.001, c('Term', 'Significant', 'p.value')][ order(-Significant) ]
b$col<-'grey39'
b$col[ grep('chrom|telomere', b$Term) ]<-'deepskyblue3'
MAX<-tail(pretty(b$Significant, 5), 1)
br<-barplot(b$Significant, border='white', col=b$col, horiz=T, axisnames=F, xlab='', ylab='', las=1, yaxt='n', xlim=c(0, MAX))
axis(2, at=br, labels=b$Term, line=0, tick=F, cex.axis=1.8)
mtext('Count', side=1, line=2, padj=+0.8, las=0, cex=1.8)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/topGO_genes_geq_65p_samples_biclustered_all_2_MNA_cluster.pdf', width=30, height=22, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/topGO_genes_geq_65p_samples_biclustered_all_2_MNA_cluster.svg', width=30, height=22, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b,MAX,br)


dev.off()

#}}}


#  [primary MNA cluster] enrichment in MYCN targets
#{{{

#  load backSPIN clustering results
gns.bs<-read.cef('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_geq_65p_samples_biclustered_all_2.cef')


#  group of genes that determine the MNA cluster and the rest
grp<-names(gns.bs$annot_genes[ gns.bs$annot_genes==2 ])
rst<-names(gns.bs$annot_genes[ gns.bs$annot_genes!=2 ])
stopifnot( length(grp)+length(rst)==length(gns.bs$annot_genes) )


#  MYCN targets 
#  separate to induced and repressed
load('/data/annotation/GRCh38/MYCN_targets.RData')
IND<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'induced', gene_id])
REP<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'repressed', gene_id])
stopifnot( length(IND)+length(REP)==nrow(MYCN) )
rm(MYCN)


#  [induced targets] Fisher's exact test:
#
#                   |  in-MNA-cluster  |  out-of-MNA-cluster  |
#  -----------------|------------------|----------------------|
#      MYCN target  |          x1      |      y1              |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |          x2      |      y2              | 
fisher.test(data.frame('in'=c(length(intersect(grp, IND)), 
                              length(setdiff(grp, IND))), 
                       'out'=c(length(intersect(rst, IND)), 
                               length(setdiff(rst, IND)))), alternative='greater')$p.value
#
#  => 3.8e-77


#  there are no repressed MYCN targets in the MNA group
intersect(grp, REP)  #  ENSG00000159409

#}}}


#  [secondary MNA cluster] enrichment in MYCN targets
#{{{

#  load backSPIN clustering results
gns.bs<-read.cef('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_geq_65p_samples_biclustered_all_2.cef')


#  group of genes that determine the last MNA cluster and the rest
grp<-names(gns.bs$annot_genes[ gns.bs$annot_genes==3 ])
rst<-names(gns.bs$annot_genes[ gns.bs$annot_genes!=3 ])
stopifnot( length(grp)+length(rst)==length(gns.bs$annot_genes) )


#  MYCN targets 
#  separate to induced and repressed
load('/data/annotation/GRCh38/MYCN_targets.RData')
IND<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'induced', gene_id])
REP<-sub('\\.[0-9]*$', '', MYCN[ type %in% 'repressed', gene_id])
stopifnot( length(IND)+length(REP)==nrow(MYCN) )
rm(MYCN)


#  [induced targets] Fisher's exact test:
#
#                   |  in-MNA-cluster  |  out-of-MNA-cluster  |
#  -----------------|------------------|----------------------|
#      MYCN target  |          x1      |      y1              |
#  -----------------|------------------|----------------------|
#  non-MYCN target  |          x2      |      y2              | 
fisher.test(data.frame('in'=c(length(intersect(grp, IND)), 
                              length(setdiff(grp, IND))), 
                       'out'=c(length(intersect(rst, IND)), 
                               length(setdiff(rst, IND)))), alternative='greater')$p.value
#
#  => 1.0


#  there are no repressed MYCN targets in the MNA group
intersect(grp, REP)  #  ENSG00000108830 ENSG00000108309 ENSG00000128266 ENSG00000175906 ENSG00000183828

#}}}

#}}}

#}}}


#  [do not run but keep around for reference] biclustering of MSigDB neuroblastoma-specific genes
#{{{

#  gene raw counts normalized by library-sizes and rounded up
n.<-apply(sweep(nb.cnt[ rownames(nb.cnt) %in% rownames(nb), ], 2, nb.sf, '/'), 2, ceiling)


#  [only for count data] plot log2(CV) versus log2(Mean) and the non-linear fit
cl<-setNames(rep('darkgrey', nrow(n.)), rownames(n.))
cl[ names(cl) %in% rownames(nb) ]<-'black'
plot_cv(n., with.r=F, pt.cex=0.5, col.genes=cl)


#  save as CEF
write.cef(n., hd=c('DGE', 'neuroblastoma genes based on MSigDB analysis'), fl='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_MSigDB.cef')


#  run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_MSigDB.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_MSigDB_biclustered_all_2.cef -v -d 2 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_MSigDB_biclustered_all_2.cef', scale.dge=T)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_genes_MSigDB_biclustered_all_2.pdf', width=16, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_genes_MSigDB_biclustered_all_2.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [do not run but keep around for reference] biclustering based on histone genes
#{{{

#  identify histone genes and pseudogenes
h<-hsa$gene_id[ grep('^HIST|^H2A|^H2B', hsa$gene_name) ]
h<-hsa$gene_id[ grep('^H2A|^H2B', hsa$gene_name) ]


#  gene raw counts normalized by library-sizes and rounded up
n.<-apply(sweep(nb.cnt[ rownames(nb.cnt) %in% h, ], 2, nb.sf, '/'), 2, ceiling)


#  [only for count data] plot log2(CV) versus log2(Mean) and the non-linear fit
cl<-setNames(rep('darkgrey', nrow(n.)), rownames(n.))
cl[ names(cl) %in% rownames(nb) ]<-'black'
plot_cv(n., with.r=T, pt.cex=0.5, col.genes=cl)


#  save as CEF
write.cef(n., hd=c('DGE', 'histone genes'), fl='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/histone_genes.cef')


#  run backSPIN with maximum 2 splits (2^2 maximum number of clusters)
system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/histone_genes.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/histone_genes_biclustered_all_2.cef -v -d 2 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/histone_genes_biclustered_all_2.cef', scale.dge=T)

#}}}


#  [do not run but keep around for reference: MNA vs HR_nMNA DE results] biclustering based on specific genes associated with the DE analysis
#{{{

#  load genes DE between MNA vs HR_nMNA
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/DESeq2_genes_MNA_HR_nMNA.RData')
x<-subset(RES, padj<0.05)
rownames(x)<-sub('\\.[0-9]*$', '', rownames(x))
rm(CND, DDS, PCA, VE, VSC, RES)


#  load circRNAs DE between MNA vs HR_nMNA
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unique/DESeq2_circRNAs_CIRI2_MNA_HR_nMNA.RData')
y<-subset(RES, padj<0.05)
rownames(y)<-sub('\\.[0-9]*$', '', rownames(y))
rm(CND, DDS, PCA, VE, VSC, RES)


#  genes common between them
common<-intersect( rownames(x), rownames(y) )
common<-data.frame( gene=x[common, 'baseMean'], gene_log2FC=x[common, 'log2FoldChange'], gene_padj=x[common, 'padj'],
                    circ=y[common, 'baseMean'], circ_log2FC=y[common, 'log2FoldChange'], circ_padj=y[common, 'padj'],
                    gene_name=y[common, 'gene_name'], row.names=common)


#  gene raw counts normalized by library-sizes and rounded up
#n.<-apply(sweep(nb.cnt[ rownames(common), ], 2, nb.sf, '/'), 2, ceiling)     #  common DE genes + circRNAs
#n.<-apply(sweep(nb.cnt[ rownames(y), ], 2, nb.sf, '/'), 2, ceiling)          #  genes of DE circRNAs 
n.<-apply(sweep(nb.cnt[ rownames(x), ], 2, nb.sf, '/'), 2, ceiling)          #  DE genes (best results)


#  [only for count data] plot log2(CV) versus log2(Mean) and the non-linear fit
cl<-setNames(rep('darkgrey', nrow(n.)), rownames(n.))
cl[ names(cl) %in% rownames(nb) ]<-'black'
plot_cv(n., with.r=F, pt.cex=0.5, col.genes=cl)


#  save as CEF
write.cef(n., hd=c('DGE', 'MNA vs HR_nMNA specific genes'), fl='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_specific_MNA_HR_nMNA.cef')


#  run backSPIN with maximum 1 split (2 maximum number of clusters) since we have DE results
system('backspin -i /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_specific_MNA_HR_nMNA.cef -o /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_specific_MNA_HR_nMNA_biclustered_all_1.cef -v -d 1 -s 0.1 -S 0.1')
bs<-backspin_process('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/genes_specific_MNA_HR_nMNA_biclustered_all_1.cef', scale.dge=T)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_genes_specific_MNA_HR_nMNA_biclustered_all_1.pdf', width=16, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/heatmap_genes_specific_MNA_HR_nMNA_biclustered_all_1.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#
#  process results for genes + circRNAs
#


#  overlap of GO enriched terms associated with genes and circRNAs biclustering MNA group
#{{{

#  load GO enrichment analysis for circRNAs biclustering MNA cluster
#  keep significantly enriched terms only
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/topGO_circRNAs_CIRI2_geq_65p_samples_biclustered_all_2_MNA_cluster.RData')
circ.bp<-bp[ p.value<0.05 ]
circ.mf<-mf[ p.value<0.05 ]
rm(bp,mf,grp,uni,allG)


#  load GO enrichment analysis for genes biclustering MNA cluster
#  keep significantly enriched terms only
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/topGO_genes_geq_65p_samples_biclustered_all_2_MNA_cluster.RData')
gns.bp<-bp[ p.value<0.05 ]
gns.mf<-mf[ p.value<0.05 ]
rm(bp,mf,grp,uni,allG)


#  common terms
common<-intersect( gns.bp$GO.ID , circ.bp$GO.ID )
circ.bp<-circ.bp[ GO.ID %in% common ]
gns.bp<-gns.bp[ GO.ID %in% common ]
common<-intersect( gns.mf$GO.ID , circ.mf$GO.ID )
circ.mf<-circ.mf[ GO.ID %in% common ]
gns.mf<-gns.mf[ GO.ID %in% common ]
rm(common)


#  plot with p-value cutoff < 0.05
x11(width=30, height=18, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(4.0,27.0,0.5,1.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.4, cex.axis=1.4)
b<-gns.bp[, c('Term', 'Significant', 'p.value')][ order(-Significant) ]
b$col<-'grey39'
b$col[ grep('chrom|telomere', b$Term) ]<-'deepskyblue3'
MAX<-tail(pretty(b$Significant, 5), 1)
br<-barplot(b$Significant, border='white', col=b$col, horiz=T, axisnames=F, xlab='', ylab='', las=1, yaxt='n', xlim=c(0, MAX))
axis(2, at=br, labels=b$Term, line=0, tick=F, cex.axis=1.8)
mtext('Count', side=1, line=2, padj=+0.8, las=0, cex=1.8)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/topGO_genes+circRNAs_CIRI2_geq_65p_samples_biclustered_all_2_MNA_cluster.pdf', width=30, height=18, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/backspin/figures/topGO_genes+circRNAs_CIRI2_geq_65p_samples_biclustered_all_2_MNA_cluster.svg', width=30, height=18, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b,MAX,br)


dev.off()

#}}}

#}}}



#  APClustering
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(vsn)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(pheatmap)
library(openxlsx)
library(apcluster)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load pre-prepared counts etc. for all samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')


#  cluster and plot
m<-meta[!is.na(risk_group), ]
x<-gns.vs[, m$bid]
#x<-crs.vs[, m$bid]
colnames(x)<-sub('-11-R01', '', colnames(x))
d<-dist(t(x), method='euclidean')
hc<-hclust(d, method='ward.D2')
ex<-data.frame(Type=factor(m$risk_group, exclude=F), row.names=colnames(x))
cl<-setNames(list(setNames( unique(m$col), unique(m$risk_group) )), colnames(ex))

#  affinity propagating clustering
apc<-apcluster(-as.matrix(d), details=T, includeSim=T, nonoise=T, q=0.1)
agg<-aggExCluster(x=apc)





#  hierarchical clustering showing the cut at the affinity propagating number of clusters length as well
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
ex$hierarchical<-as.character(cutree(hc, k=length(apc@exemplars)))
cl$hierarchical<-setNames(colorRampPalette(brewer.pal(8, 'Dark2'))(length(apc@exemplars)), unique(ex$hierarchical))
p1<-pheatmap(as.matrix(d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=hc,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0, silent=F)


#  affinity propagating clustering
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
ex$apc<-setNames(as.character(rep(seq_along(apc@clusters), lengths(apc@clusters))), names(unlist(apc@clusters)))[ rownames(ex) ]
cl$apc<-setNames(colorRampPalette(brewer.pal(8, 'Set2'))(length(apc@exemplars)), seq_along(apc@clusters))
p2<-pheatmap(as.matrix(d)[unlist(apc@clusters), unlist(apc@clusters)], color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=F,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0, silent=T)
dev.set(dev.next())
print(p2)


#  compare the results
summary(silhouette(as.integer(ex$hierarchical), dmatrix=as.matrix(d))) 
summary(silhouette(as.integer(ex$apc), dmatrix=as.matrix(d))) 

#}}}



#  [MYCN +-Tet] check circRNAs expression for Steffen
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(vsn)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load unified cohort of circRNAs identified in the cell models
#  isolate the current cell model
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_cell_models.RData')
meta<-meta[ cell_model %in% 'MYCN' & treatment %in% c('+Tet 4h', 'ETOH 4h') ]
circ<-CIRCS.all[ CIRCS.all$bid %in% meta$bid ]
rm(CIRCS.all)


#  [circRNAs, 4h] load DE results 
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_+Tet 4h_ETOH 4h.RData')
res.4h<-RES
dds.4h<-DDS
vsc.4h<-VSC
rm(list=setdiff(ls(), c(l, ls(pattern='\\.4h'))))


#  [circRNAs, 48h] load DE results 
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_+Tet 48h_ETOH 48h.RData')
res.48h<-RES
dds.48h<-DDS
vsc.48h<-VSC
rm(list=setdiff(ls(), c(l, ls(pattern='\\.48h'))))


#  validated circRNAs (hg38):
#
#     ARID1A  chr1:26729651-26732792   hsa_circ_0008494
#       EYA1  chr8:71299047-71334174   hsa_circ_0084764 (also: chr8:71269740-71334174, chr8:71317552-71334174)
#      SETD3  chr14:99458279-99465813  hsa_circ_0000567
#      HIPK3  chr11:33286413-33287511  hsa_circ_0000284
#    SMARCA5  chr4:143543509-143543972 hsa_circ_0001445
#     ZNF124  chr1:247159006-247159813
#     NFATC3  chr16:68121987-68126610  (also: chr16:68121987-68123121)
#      ASAP1  chr8:130152736-130180880 (also: chr8:130358017-130361771, chr8:130152736-130169067)
#       TET2  chr4:105233897-105237351 hsa_circ_0070562
#      MAML3  chr4:139889357-139890967 
#     AKAP12  chr6:151348711-151353752
#    CDR1-AS  chrX:140783175-140784659 hsa_circ_0001946
#
#  check if they are downregulated (even not significantly so):
check<-c('ARID1A', 'EYA1', 'SETD3', 'HIPK3', 'SMARCA5', 'ZNF124', 'NFATC3', 'ASAP1', 'TET2', 'MAML3', 'AKAP12', 'LINC00632')
check.4h<-res.4h[match(check, res.4h$gene_name), ]
check.48h<-res.48h[match(check, res.48h$gene_name), ]


#  [4h] show only downregulated ones
down.4h<-subset(check.4h, log2FoldChange<0)
#                             baseMean      log2FoldChange              lfcSE                stat            pvalue              padj   gene_name
#                            <numeric>           <numeric>          <numeric>           <numeric>         <numeric>         <numeric> <character>
# ENSG00000183576.13  74.4802602025483  -0.075158209689051   0.17949476284511  -0.418817797319143 0.675349296959083 0.998279639555248       SETD3
# ENSG00000110422.12  172.685516401783  -0.187935507535934  0.133540452179712   -1.40712306655049 0.159390899337655 0.998279639555248       HIPK3
# ENSG00000196418.12   11.276836300825  -0.102620057580469  0.144115337710206  -0.726452386590776 0.467561486205406 0.998279639555248      ZNF124
# ENSG00000196782.12  4.81743202711883 -0.0153755760156152  0.182185644976078 -0.0843609394257204 0.932769462362703 0.998279639555248       MAML3
# ENSG00000131016.17  5.98187740012636  -0.158403581040284  0.193751180272878  -0.817829975415627 0.413454279298537 0.998279639555248      AKAP12
# ENSG00000203930.11 0.865471752068695 -0.0144660753836937 0.0582273132370359  -0.252599959060137 0.800577361658186 0.998279639555248   LINC00632



#  [4h] show only downregulated ones
down.48h<-subset(check.48h, log2FoldChange<0)
#                            baseMean       log2FoldChange             lfcSE                stat            pvalue              padj   gene_name
#                           <numeric>            <numeric>         <numeric>           <numeric>         <numeric>         <numeric> <character>
# ENSG00000117713.20 69.0520392046263   -0.118533624349954 0.194290037846314  -0.610085316162143 0.541805293105961 0.998685166329935      ARID1A
# ENSG00000104313.19 33.5153118826472  -0.0184263695142833 0.233818058706799   -0.07874380960964  0.93723639861339 0.998685166329935        EYA1
# ENSG00000183576.13 87.2907306154451   -0.185735200966403  0.18608764644133  -0.997388720852996 0.318575864020125 0.998685166329935       SETD3
# ENSG00000110422.12 174.265048017466   -0.039264794681416 0.156146749091993  -0.251437053941074 0.801476222739776 0.998685166329935       HIPK3
# ENSG00000153147.6  214.044196128637 -0.00861794416715385 0.176449949942198 -0.0488339457551621 0.961051629635781 0.998685166329935     SMARCA5
# ENSG00000072736.19 97.7688744135811  -0.0353516850112922 0.184951608740055  -0.191133120438533 0.848421299587909 0.998685166329935      NFATC3
# ENSG00000168769.13 22.5900953600536   -0.287973056965198 0.246949690853176   -1.16290813827775 0.244866775153578 0.998685166329935        TET2
# ENSG00000196782.12 7.68139769544028  -0.0173869662527712 0.250803994374179 -0.0693421239740907 0.944717296328914 0.998685166329935       MAML3


#  common at 4h and 48h
intersect(down.4h$gene_name, down.48h$gene_name)
#
# SETD3 HIPK3 MAML3


#  make a matix for easy copy-pase
x<-data.frame(round(as.matrix(down.4h[, c('baseMean', 'log2FoldChange', 'pvalue', 'padj')]), 3), row.names=down.4h$gene_name)
x<-data.frame(round(as.matrix(down.48h[, c('baseMean', 'log2FoldChange', 'pvalue', 'padj')]), 3), row.names=down.48h$gene_name)

#}}}



#  [genes of interest] TPM expression stratified by risk-group
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_ecdfs.R')
source('~/bio/lib/grouplist2boxplot.R')


#  fix the ggplot2 theme
#{{{
theme_set(theme_bw(base_size=35) + theme(panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_blank(),
    plot.margin=margin(t=1.0, b=0.1, l=0.5, r=1.0, unit='cm'),
    axis.line=element_line(color='black', size=0.5),
    axis.text=element_text(color='black', size=40, family='Arial', face='plain'), 
    axis.text.x=element_text(margin=margin(t=0.2, b=0, unit='cm')),
    axis.title=element_text(color='black', size=45, family='Arial', face='plain'), 
    axis.title.x=element_text(color='black', margin=margin(t=20, b=1)),
    axis.title.y=element_text(color='black', margin=margin(r=20, l=1)),
    axis.ticks=element_line(color='black', size=0.5), 
    axis.ticks.length=unit(0.5,  'cm'),
    legend.background=element_blank(),
    legend.justification=c(1, 1), 
    legend.position=c(0.15, 1.10), 
    legend.margin=margin(),
    legend.key=element_blank(),
    legend.title=element_text(size=24, family='Arial', face='plain', margin=margin(b=0.5, unit='cm')),
    legend.text=element_text(size=24, family='Arial', face='plain')
))
#}}}


#  load gene expression (TPMs) and keep only tumors
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')
meta<-meta[ !is.na(risk_group) ]
stopifnot( length(setdiff(meta[, bid], colnames(gns.tpm)))==0 )
gns.tpm<-gns.tpm[, meta[, bid]]
gns.cnt<-gns.cnt[, meta[, bid]]
rm(list=setdiff(ls(), c(l, 'meta', 'hsa', 'gns.tpm', 'gns.cnt')))


#  [boxplot list object] we convert the gene log10(1+TPM) matrix to a list of matrices separating the samples by risk group 
L<-log10(1+gns.tpm)
stopifnot( all.equal( colnames(L), meta[, bid] ) )
colnames(L)<-meta[, risk_group]
L<-setNames(lapply(unique(colnames(L)), function(x){ L[, colnames(L) %in% x, drop=F] }), unique(colnames(L)))
L.COL<-setNames(meta[, unique(col)], meta[, unique(risk_group)])


#  [stripchart data.frame object] we melt the gene log10(1+TPM) matrix to (gene_name, tpm, risk_group) data.frame
S<-data.frame(log10(1+gns.tpm), check.names=F)
stopifnot( all.equal( colnames(S), meta[, bid] ) )
colnames(S)<-meta[, risk_group]
S<-data.table(t(S), risk_group=colnames(S))
S<-data.frame(melt(S, id.vars=c('risk_group'), variable.name='gene_name', value.name='tpm'), check.names=F)
S$risk_group<-factor(S$risk_group, levels=unique(S$risk_group))
S.COL<-L.COL


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  MYCN, MYC, MYCL
#{{{

#  pick the genes of interest
GENES<-c('MYCN', 'MYC', 'MYCL')


#  boxplots across tumors
grouplist2boxplot(L=lapply(L, function(x){ x[ match(GENES, rownames(x)), ,drop=F] }), L.COL=L.COL, YMAX=3, XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_across_tumors.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  stripchart
B<-S[ S$gene_name %in% GENES, , drop=F]
B$gene_name<-factor(B$gene_name, levels=GENES)
YTICK<-pretty(c(0, max(B$tpm, na.rm=T)), 10)
XHIGH<-seq(2, length(unique(B$gene_name)), 2)
ggplot(B, aes(x=gene_name, y=tpm, color=risk_group, group=gene_name)) + 
    geom_jitter(height=0, width=0.2, size=6) +
    geom_rect(data=data.frame(xmin=XHIGH-0.5, xmax=XHIGH+0.5), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), 
              color='white', fill='grey39', alpha=0.2, show.legend=F, inherit.aes=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    labs(x='', y=expression(log[10]('1+TPM')), color='') + 
    scale_color_manual(values=S.COL) +
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5),
          plot.margin=margin(t=1.0, b=-2.0, l=0.5, r=1.0, unit='cm'),
          legend.justification=c(1, 1), 
          legend.position=c(1.02, 1.10)
    )
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/stripchart_', paste0(GENES, collapse='+'), '_across_tumors.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  KHSRP
#{{{

#  pick the genes of interest
GENES<-c('KHSRP')


#  boxplots across tumors
grouplist2boxplot(L=lapply(L, function(x){ x[ match(GENES, rownames(x)), ,drop=F ] }), L.COL=L.COL, YMAX=2.5, XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_across_tumors.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  stripchart
B<-S[ S$gene_name %in% GENES, , drop=F]
B$gene_name<-factor(B$gene_name, levels=GENES)
YTICK<-pretty(range(B$tpm, na.rm=T), 4)
#XHIGH<-seq(2, length(unique(B$gene_name)), 2)
ggplot(B, aes(x=gene_name, y=tpm, color=risk_group, group=gene_name)) + 
    geom_jitter(height=0, width=0.2, size=6) +
    #geom_rect(data=data.frame(xmin=XHIGH-0.5, xmax=XHIGH+0.5), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), 
    #          color='white', fill='grey39', alpha=0.2, show.legend=F, inherit.aes=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    labs(x='', y=expression(log[10]('1+TPM')), color='') + 
    scale_color_manual(values=S.COL) +
    coord_cartesian(clip='off') +
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5),
          plot.margin=margin(t=1.0, b=-2.0, l=0.5, r=1.0, unit='cm'),
          legend.justification=c(1, 1), 
          legend.position=c(0.22, 1.10)
    )
ggsave(paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/stripchart_', paste0(GENES, collapse='+'), '_across_tumors.svg'), device=svg, width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  DHX9
#{{{

#  pick the genes of interest
GENES<-c('DHX9')


#  boxplots across tumors
grouplist2boxplot(L=lapply(L, function(x){ x[ match(GENES, rownames(x)), ,drop=F ] }), L.COL=L.COL, YMAX=3, XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_across_tumors.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  stripchart
B<-S[ S$gene_name %in% GENES, , drop=F]
B$gene_name<-factor(B$gene_name, levels=GENES)
YTICK<-pretty(range(B$tpm, na.rm=T), 4)
#XHIGH<-seq(2, length(unique(B$gene_name)), 2)
ggplot(B, aes(x=gene_name, y=tpm, color=risk_group, group=gene_name)) + 
    geom_jitter(height=0, width=0.2, size=6) +
    #geom_rect(data=data.frame(xmin=XHIGH-0.5, xmax=XHIGH+0.5), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), 
    #          color='white', fill='grey39', alpha=0.2, show.legend=F, inherit.aes=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    labs(x='', y=expression(log[10]('1+TPM')), color='') + 
    scale_color_manual(values=S.COL) +
    coord_cartesian(clip='off') +
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5),
          plot.margin=margin(t=1.0, b=-2.0, l=0.5, r=1.0, unit='cm'),
          legend.justification=c(1, 1), 
          legend.position=c(0.22, 1.10)
    )
ggsave(paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/stripchart_', paste0(GENES, collapse='+'), '_across_tumors.svg'), device=svg, width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  YBX1, HNRNPF
#{{{

#  pick the genes of interest
GENES<-c('YBX1', 'HNRNPF')


#  boxplots across tumors
grouplist2boxplot(L=lapply(L, function(x){ x[ match(GENES, rownames(x)), ,drop=F ] }), L.COL=L.COL, YMAX=3, XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, LEGEND='topright', DRAW.HIGHLIGHT=F, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_across_tumors.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  stripchart
B<-S[ S$gene_name %in% GENES, , drop=F]
B$gene_name<-factor(B$gene_name, levels=GENES)
YTICK<-pretty(range(B$tpm, na.rm=T), 4)
XHIGH<-seq(2, length(unique(B$gene_name)), 2)
ggplot(B, aes(x=gene_name, y=tpm, color=risk_group, group=gene_name)) + 
    geom_jitter(height=0, width=0.2, size=6) +
    geom_rect(data=data.frame(xmin=XHIGH-0.5, xmax=XHIGH+0.5), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), 
              color='white', fill='grey39', alpha=0.2, show.legend=F, inherit.aes=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    labs(x='', y=expression(log[10]('1+TPM')), color='') + 
    scale_color_manual(values=S.COL) +
    coord_cartesian(clip='off') +
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5),
          plot.margin=margin(t=1.0, b=-2.0, l=0.5, r=1.0, unit='cm'),
          legend.justification=c(1, 1), 
          legend.position=c(0.22, 1.10)
    )
ggsave(paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/stripchart_', paste0(GENES, collapse='+'), '_across_tumors.svg'), device=svg, width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [MNA, HR_nMNA only] DHX9, YBX1, HNRNPF
#{{{

#  pick the genes of interest
GENES<-c('DHX9', 'YBX1', 'HNRNPF')


#  boxplots across tumors
grouplist2boxplot(L=lapply(L[c('HR_nMNA', 'MNA')], function(x){ x[ match(GENES, rownames(x)), ,drop=F ] }), L.COL=L.COL[c('HR_nMNA', 'MNA')], YMAX=3, XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_HR_nMNA+MNA.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  stripchart
B<-S[ S$gene_name %in% GENES & S$risk_group %in% c('HR_nMNA', 'MNA'), , drop=F]
B$gene_name<-factor(B$gene_name, levels=GENES)
YTICK<-pretty(range(B$tpm, na.rm=T), 4)
XHIGH<-seq(2, length(unique(B$gene_name)), 2)
ggplot(B, aes(x=gene_name, y=tpm, color=risk_group, group=gene_name)) + 
    geom_jitter(height=0, width=0.2, size=6) +
    geom_rect(data=data.frame(xmin=XHIGH-0.5, xmax=XHIGH+0.5), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), 
              color='white', fill='grey39', alpha=0.2, show.legend=F, inherit.aes=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    labs(x='', y=expression(log[10]('1+TPM')), color='') + 
    scale_color_manual(values=S.COL) +
    coord_cartesian(clip='off') +
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5),
          plot.margin=margin(t=1.0, b=-2.0, l=0.5, r=1.0, unit='cm'),
          legend.justification=c(1, 1), 
          legend.position=c(0.22, 1.10)
    )
ggsave(paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/stripchart_', paste0(GENES, collapse='+'), '_HR_nMNA+MNA.svg'), device=svg, width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}



#  [MYCN +-Tet] TPM expression of genes of interest
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_ecdfs.R')
source('~/bio/lib/grouplist2boxplot.R')


#  fix the ggplot2 theme
#{{{
theme_set(theme_bw(base_size=35) + theme(panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_blank(),
    plot.margin=margin(t=1.0, b=0.1, l=0.5, r=1.0, unit='cm'),
    axis.line=element_line(color='black', size=0.5),
    axis.text=element_text(color='black', size=40, family='Arial', face='plain'), 
    axis.text.x=element_text(margin=margin(t=0.2, b=0, unit='cm')),
    axis.title=element_text(color='black', size=45, family='Arial', face='plain'), 
    axis.title.x=element_text(color='black', margin=margin(t=20, b=1)),
    axis.title.y=element_text(color='black', margin=margin(r=20, l=1)),
    axis.ticks=element_line(color='black', size=0.5), 
    axis.ticks.length=unit(0.5,  'cm'),
    legend.background=element_blank(),
    legend.justification=c(1, 1), 
    legend.position=c(0.15, 1.10), 
    legend.margin=margin(),
    legend.key=element_blank(),
    legend.title=element_text(size=24, family='Arial', face='plain', margin=margin(b=0.5, unit='cm')),
    legend.text=element_text(size=24, family='Arial', face='plain')
))
#}}}


#  [4h, 48h] load gene expression (TPMs) and split by condition
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')
meta<-meta[grepl('CB-SKNAS-TR-MYCN', bid) & treatment %in% c('ETOH 4h', '+Tet 4h', 'ETOH 48h', '+Tet 48h') ]
gns.tpm<-gns.tpm[, meta[, bid]]
gns.4h<-gns.tpm[, meta[ grep('4h', treatment), bid] ]
gns.48h<-gns.tpm[, meta[ grep('48h', treatment), bid] ]
rm(list=setdiff(ls(), c(l, 'hsa', 'meta', 'gns.4h', 'gns.48h')))


#  [120h] load gene expression (TPMs) and split by condition
#         replace gene_ids with gene_names in the rows
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_circRNAs+genes_all_libraries.RData')
gns.120h<-gns.tpm[, grep('120h', colnames(gns.tpm))]
rownames(gns.120h)<-hsa$gene_name[ match( rownames(gns.120h), hsa$gene_id ) ]
rm(list=setdiff(ls(), c(l, 'gns.120h')))


#  organize log10(1+TPMs) by condition in a list
L<-List()
L[['ETOH 4h']]<-log10(1+gns.4h[, grep('ETOH', colnames(gns.4h))])
L[['Tet 4h']]<-log10(1+gns.4h[, grep('ETOH', colnames(gns.4h), invert=T)])
L[['ETOH 48h']]<-log10(1+gns.48h[, grep('ETOH', colnames(gns.48h))])
L[['Tet 48h']]<-log10(1+gns.48h[, grep('ETOH', colnames(gns.48h), invert=T)])
L[['ETOH 120h']]<-log10(1+gns.120h[, grep('ETOH', colnames(gns.120h))])
L[['Tet 120h']]<-log10(1+gns.120h[, grep('ETOH', colnames(gns.120h), invert=T)])
L.COL<-setNames(c('cadetblue4', adjustcolor('cadetblue4', offset=c(0.2, 0.2, 0.2, 0.0)), 
                   'darkcyan', adjustcolor('darkcyan', offset=c(0.2, 0.2, 0.2, 0.0)),
                   'darkolivegreen', adjustcolor('darkolivegreen', offset=c(0.2, 0.2, 0.2, 0.0))), names(L))
rm(gns.4h, gns.48h, gns.120h)


#  melt all log10(1+TPMs) in a data.frame of (gene_name, tpm, condition) columns
S<-data.table(t(L[['ETOH 4h']]), condition='ETOH 4h')
S<-rbind(S, data.table(t(L[['Tet 4h']]), condition='Tet 4h'))
S<-data.frame(melt(S, id.vars=c('condition'), variable.name='gene_name', value.name='tpm'), check.names=F)
x<-data.table(t(L[['ETOH 48h']]), condition='ETOH 48h')
x<-rbind(x, data.table(t(L[['Tet 48h']]), condition='Tet 48h'))
x<-data.frame(melt(x, id.vars=c('condition'), variable.name='gene_name', value.name='tpm'), check.names=F)
S<-rbind(S, x)
x<-data.table(t(L[['ETOH 120h']]), condition='ETOH 120h')
x<-rbind(x, data.table(t(L[['Tet 120h']]), condition='Tet 120h'))
x<-data.frame(melt(x, id.vars=c('condition'), variable.name='gene_name', value.name='tpm'), check.names=F)
S<-rbind(S, x)
S$condition<-factor(S$condition, levels=c('ETOH 4h', 'Tet 4h', 'ETOH 48h', 'Tet 48h', 'ETOH 120h', 'Tet 120h'))
S.COL<-setNames(c('cadetblue4', adjustcolor('cadetblue4', offset=c(0.2, 0.2, 0.2, 0.0)), 
                   'darkcyan', adjustcolor('darkcyan', offset=c(0.2, 0.2, 0.2, 0.0)),
                   'darkolivegreen', adjustcolor('darkolivegreen', offset=c(0.2, 0.2, 0.2, 0.0))), levels(S$condition))
rm(x)


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  DHX9 
#{{{

#  pick the genes of interest
GENES<-'DHX9'


#  boxplots 
grouplist2boxplot(L=lapply(L, function(x){ x[ match(GENES, rownames(x)), ,drop=F] }), L.COL=L.COL, YMIN=1.9, YMAX=2.4, XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, LEGEND='topleft', mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_MYCN_Tet-inducible.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  stripchart
B<-S[ S$gene_name %in% GENES, , drop=F]
B$gene_name<-factor(B$gene_name, levels=GENES)
#YTICK<-pretty(c(0, max(B$tpm, na.rm=T)), 10)
YTICK<-pretty(c(1.9, 2.4), 4)
#XHIGH<-seq(2, length(unique(B$gene_name)), 2)
ggplot(B, aes(x=gene_name, y=tpm, color=condition, group=gene_name)) + 
    geom_jitter(height=0, width=0.2, size=10) +
    #geom_rect(data=data.frame(xmin=XHIGH-0.5, xmax=XHIGH+0.5), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), 
    #          color='white', fill='grey39', alpha=0.2, show.legend=F, inherit.aes=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    labs(x='', y=expression(log[10]('1+TPM')), color='') + 
    scale_color_manual(values=S.COL) +
    coord_cartesian(clip='off') +
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5),
          plot.margin=margin(t=1.0, b=-2.0, l=0.5, r=1.0, unit='cm'),
          legend.justification=c(1, 1), 
          legend.position=c(1.02, 1.10)
    )
ggsave(file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/stripchart_', paste0(GENES, collapse='+'), '_MYCN_Tet-inducible.svg'), device=svg, width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  YBX1, HNRNPF
#{{{

#  pick the genes of interest
GENES<-c('YBX1', 'HNRNPF')


#  boxplots 
grouplist2boxplot(L=lapply(L, function(x){ x[ match(GENES, rownames(x)), ,drop=F] }), L.COL=L.COL, YMIN=2.0, YMAX=3.0, XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, LEGEND='topright', DRAW.HIGHLIGHT=F, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_MYCN_Tet-inducible.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  stripchart
B<-S[ S$gene_name %in% GENES, , drop=F]
B$gene_name<-factor(B$gene_name, levels=GENES)
YTICK<-pretty(range(B$tpm, na.rm=T), 4)
#YTICK<-pretty(c(1.9, 2.4), 4)
XHIGH<-seq(2, length(unique(B$gene_name)), 2)
ggplot(B, aes(x=gene_name, y=tpm, color=condition, group=gene_name)) + 
    geom_jitter(height=0, width=0.2, size=10) +
    #geom_rect(data=data.frame(xmin=XHIGH-0.5, xmax=XHIGH+0.5), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), 
    #          color='white', fill='grey39', alpha=0.2, show.legend=F, inherit.aes=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    labs(x='', y=expression(log[10]('1+TPM')), color='') + 
    scale_color_manual(values=S.COL) +
    coord_cartesian(clip='off') +
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5),
          plot.margin=margin(t=1.0, b=-2.0, l=0.5, r=1.0, unit='cm'),
          legend.justification=c(1, 1), 
          legend.position=c(1.02, 1.10)
    )
ggsave(file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/stripchart_', paste0(GENES, collapse='+'), '_MYCN_Tet-inducible.svg'), device=svg, width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [48h, 120h only] DHX9, YBX1, HNRNPF
#{{{

#  pick the genes of interest
GENES<-c('DHX9', 'YBX1', 'HNRNPF')


#  boxplots 
grouplist2boxplot(L=lapply(L[grep('4h', names(L), invert=T)], function(x){ x[ match(GENES, rownames(x)), ,drop=F] }), L.COL=L.COL[grep('4h', names(L.COL), invert=T)], YMIN=1.6, YMAX=3.0, XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, LEGEND='topleft', DRAW.HIGHLIGHT=T, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_MYCN_Tet-inducible_48h+120h.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  stripchart
B<-S[ S$gene_name %in% GENES & ! grepl('4h', S$condition), , drop=F]
B$gene_name<-factor(B$gene_name, levels=GENES)
YTICK<-pretty(range(B$tpm, na.rm=T), 4)
#YTICK<-pretty(c(1.9, 2.4), 4)
XHIGH<-seq(2, length(unique(B$gene_name)), 2)
ggplot(B, aes(x=gene_name, y=tpm, color=condition, group=gene_name)) + 
    geom_jitter(height=0, width=0.2, size=10) +
    geom_rect(data=data.frame(xmin=XHIGH-0.5, xmax=XHIGH+0.5), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf), 
              color='white', fill='grey39', alpha=0.2, show.legend=F, inherit.aes=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    labs(x='', y=expression(log[10]('1+TPM')), color='') + 
    scale_color_manual(values=S.COL) +
    coord_cartesian(clip='off') +
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5),
          plot.margin=margin(t=1.0, b=-2.0, l=0.5, r=1.0, unit='cm'),
          legend.justification=c(1, 1), 
          legend.position=c(0.22, 1.10)
    )
ggsave(file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/stripchart_', paste0(GENES, collapse='+'), '_MYCN_Tet-inducible_48h+120h.svg'), device=svg, width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}



#  CHECK: PCA on the Berlin and Peifer cohort separately
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(vsn)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load pre-prepared counts etc. for all samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  Berlin cohort
#{{{

#  isolate Berlin cohort
ber<-meta[ !is.na(risk_group) & grepl('^CB2', PAT_ID_BERLIN) ]
ber.cnt<-gns.cnt[, ber$bid ]
colnames(ber.cnt)<-sub('-11-R0[0-9]*', '', colnames(ber.cnt))


#  variance-stabilized gene counts
dds<-DESeqDataSetFromMatrix(countData=ber.cnt, colData=data.frame(risk_group=factor(ber$risk_group), row.names=ber$PAT_ID_BERLIN), design=~1)
ber.sf<-sizeFactors(estimateSizeFactors(dds, type='poscounts'))
ber.vs<-assay(varianceStabilizingTransformation(dds, fitType='local'))


#  PCA on variance-stabilized (and glog2-transformed) gene counts centered but not scaled
ber.pca<-prcomp(t(ber.vs), center=T, scale.=F)
#x<-ber.vs[ apply(ber.vs, 1, function(x){ sd(x)!=0 }), ]
#ber.pca<-prcomp(t(x), center=T, scale.=T)
o<-order(apply(ber.vs, 1, function(x){ var(x)/mean(x) }), decreasing=T)
ber.pca<-prcomp(t(ber.vs[o[1:1000], ]), center=T, scale.=T)
ber.ve<-round(1000 * ber.pca$sdev^2/sum(ber.pca$sdev^2))/10 
all.equal( ber$PAT_ID_BERLIN, rownames(ber.pca$x) )


#  PCA 
n<-cbind( data.frame( ber.pca$x[, c('PC1', 'PC2')] ), Type=ber$risk_group, bid=ber$PAT_ID_BERLIN)
cl<-list('Type'=setNames( unique(ber$col), unique(ber$risk_group) ))
XLIM<-range(pretty(range(n$PC1), 2))
YLIM<-range(pretty(range(n$PC2), 2))
ggplot(n[, c('PC1', 'PC2', 'Type', 'bid')], aes(PC1, PC2, color=Type)) + 
    geom_point(size=12) + 
    scale_shape_manual(values=18) + 
    scale_fill_manual(name='Type', values=cl$Type) + 
    scale_color_manual(name='Type', values=cl$Type) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=34), axis.title.x=element_text(size=34), 
        axis.title.y=element_text(size=34), axis.text.y=element_text(size=34), 
        legend.text=element_text(size=22), legend.title=element_text(size=22, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.2, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.2, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=seq(XLIM[1], XLIM[2], length.out=5)) + 
    scale_y_continuous(lim=YLIM, breaks=seq(YLIM[1], YLIM[2], length.out=6)) + 
    xlab(paste0('PC1: ', ber.ve[1], '% variance')) + ylab(paste0('PC2: ', ber.ve[2], '% variance'))
dev.print(device=svg, file='~/Downloads/plotPCA_berlin.svg', width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  rest cohort
#{{{

#  isolate Berlin cohort
rest<-meta[ !is.na(risk_group) & !grepl('^CB2', PAT_ID_BERLIN) ]
rest.cnt<-gns.cnt[, rest$bid ]
colnames(rest.cnt)<-sub('-11-R0[0-9]*', '', colnames(rest.cnt))


#  variance-stabilized gene counts
dds<-DESeqDataSetFromMatrix(countData=rest.cnt, colData=data.frame(risk_group=factor(rest$risk_group), row.names=rest$PAT_ID_BERLIN), design=~1)
rest.sf<-sizeFactors(estimateSizeFactors(dds, type='poscounts'))
rest.vs<-assay(varianceStabilizingTransformation(dds, fitType='local'))


#  PCA on variance-stabilized (and glog2-transformed) gene counts centered but not scaled
#x<-rest.vs[ apply(rest.vs, 1, function(x){ sd(x)!=0 }), ]
#x<-rest.vs[ apply(rest.vs, 1, function(x){ sd(x)!=0 }), ]
#rest.pca<-prcomp(t(x), center=T, scale.=T)
rest.pca<-prcomp(t(rest.vs), center=T, scale.=F)
rest.ve<-round(1000 * rest.pca$sdev^2/sum(rest.pca$sdev^2))/10 
all.equal( rest$PAT_ID_BERLIN, rownames(rest.pca$x) )


#  PCA 
n<-cbind( data.frame( rest.pca$x[, c('PC1', 'PC2')] ), Type=rest$risk_group, bid=rest$PAT_ID_BERLIN)
cl<-list('Type'=setNames( unique(rest$col), unique(rest$risk_group) ))
XLIM<-range(pretty(range(n$PC1), 2))
YLIM<-range(pretty(range(n$PC2), 2))
ggplot(n[, c('PC1', 'PC2', 'Type', 'bid')], aes(PC1, PC2, color=Type)) + 
    geom_point(size=12) + 
    scale_shape_manual(values=18) + 
    scale_fill_manual(name='Type', values=cl$Type) + 
    scale_color_manual(name='Type', values=cl$Type) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=34), axis.title.x=element_text(size=34), 
        axis.title.y=element_text(size=34), axis.text.y=element_text(size=34), 
        legend.text=element_text(size=22), legend.title=element_text(size=22, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.2, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.2, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=seq(XLIM[1], XLIM[2], length.out=5)) + 
    scale_y_continuous(lim=YLIM, breaks=seq(YLIM[1], YLIM[2], length.out=6)) + 
    xlab(paste0('PC1: ', rest.ve[1], '% variance')) + ylab(paste0('PC2: ', rest.ve[2], '% variance'))

#}}}

#}}}






