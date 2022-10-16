###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




##########################################################################################
#
#
#  feature selection based on random forest classification of our cohort of neuroblastomas
#
#
##########################################################################################




#  [run once] Prepare circRNAs and genes for the analysis. We do the following:
#
#                 we remove chrY, chrM genes and unexpressed genes
#                 we sum isoform circRNA counts and gene counts of different gene_ids at the gene_name level
#                 we compute the size factors for each library based on the non-removed gene counts 
#                 we normalize the gene and circRNA counts per library by the corresponding size factor
#                 
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(vsn)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/load_circ_gene_expression.R')


#  load reference 
#  discard chrM|chrY genes from the analysis
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene']
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]
seqlevels(hsa, pruning.mode='coarse')<-seqlevels(hsa)[ grep('chrM|chrY', seqlevels(hsa), invert=T) ]


#  keep only non-failed tumor samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-meta.tum[ !(failed) ]


#  load our cohort of circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
rm(list=ls(pattern='\\.(genes|circs|meta|prefailed|all)|GENES'))


#  gene counts of all genes and circRNA counts
nb<-load_circ_gene_expression(CIRCS=CIRCS, GENES=mcols(hsa)[, c('gene_id', 'gene_name')], META=meta, nb.tumors.only=T, nb.only=T)
gns.cnt<-ceiling(nb$counts[hsa$gene_id, meta$bid]) 
gns.cnt<-data.table(gns.cnt, gene_id=rownames(gns.cnt), gene_name=hsa$gene_name[ match(rownames(gns.cnt), hsa$gene_id) ])
gns.cnt<-melt(gns.cnt, id.vars=c('gene_name', 'gene_id'), variable.name='bid', value.name='counts')[, .(counts=sum(counts)), by=.(bid, gene_name)]
gns.cnt<-dcast(gns.cnt, bid ~ gene_name, value.var='counts', fun.aggregate=sum)
gns.cnt<-t(data.frame(gns.cnt[, -1], row.names=gns.cnt[, bid], check.names=F))
crs.cnt<-data.table(data.frame(mcols(nb$circs)[, c('gene_id', 'gene_name', 'jc_count', 'bid')]))
crs.cnt<-dcast(crs.cnt, bid ~ gene_name, value.var='jc_count', fun.aggregate=sum)
crs.cnt<-t(data.frame(crs.cnt[, -1], row.names=crs.cnt[, bid], check.names=F))
rm(nb)


#  remove unexpressed genes 
gns.cnt<-gns.cnt[rowSums(gns.cnt)!=0,  ]


#  make sure metadata order and sample order is identical
gns.cnt<-gns.cnt[, meta$bid]
crs.cnt<-crs.cnt[, meta$bid]


#  compute size factors based on gene counts
#  replace gene count matrix with size-factor normalized counts
dds<-DESeqDataSetFromMatrix(countData=gns.cnt, colData=data.frame(risk_group=factor(meta$risk_group), row.names=meta$bid), design=~1)
dds<-estimateSizeFactors(dds, type='poscounts')
gns.sf<-sizeFactors(dds)
gns.cnt<-counts(dds, normalized=T)
rm(dds)


#  normalize circRNA counts as well based on the size factors of gene counts
dds<-DESeqDataSetFromMatrix(countData=crs.cnt, colData=data.frame(risk_group=factor(meta$risk_group), row.names=meta$bid), design=~1)
sizeFactors(dds)<-gns.sf
crs.cnt<-counts(dds, normalized=T)
rm(dds, gns.sf)


#  save all including the reference for convenience
save(hsa, meta, CIRCS, gns.cnt, crs.cnt, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/circRNAs+genes_normalized-counts_tumors.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/circRNAs+genes_normalized-counts_tumors.RData



#  [run once] filter genes by:
#
#                 above a given frequency of zeros across samples (normalized counts<1 are set to zero)
#                 above a given quantile of mean expression 
#                 below minimum log2(CV) value 
#                 above maximum log2(CV) value 
#
#             and add to features all unified circRNAs irrespective of gene criteria.
#{{{
rm(list=ls())
library(data.table)
library(GenomicFeatures)
library(rtracklayer)
library(matrixStats)  #  rowSds()
library(extrafont)    #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load normalized counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/circRNAs+genes_normalized-counts_tumors.RData')


#  set normalized counts<1 to zero
gns.cnt[ gns.cnt<1 ]<-0.0


#  remove features with frequency of zeros above threshold
ZEROF<-0.1 
table(pick<-apply(gns.cnt, 1, function(x){ sum(x==0)/length(x)<=ZEROF }))
# 
# FALSE  TRUE 
# 28327 25018 
#
gns.cnt<-gns.cnt[ pick, ]
rm(pick)


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')  #  recycle


#  any circRNAs/genes with zero standard deviation?
any(i<-rowSds(as.matrix(gns.cnt))==0)  #  FALSE
any(i<-rowSds(as.matrix(crs.cnt))==0)  #  FALSE


#  [genes] pick features 
#{{{

#  selection thresholds
QPICK<-0.10   #  quantile of mean expression above which we pick features
CVDIFF<-0.0   #  optional minimum log2(CV) difference from the fit model (0 turns this off)
CVMIN<- -5.0  #  minimum log2(CV) value to have irrespective of CVDIFF (-5 turns this off)
CVMAX<-5.0    #  maximum log2(CV) value to have irrespective of CVDIFF (5 turns this off)


#  fit a CV = Mean^a + b non-linear regression model
mn<-rowMeans(gns.cnt)
cv<-rowSds(as.matrix(gns.cnt))/mn 
summary(nl<-nls( cv ~ mn^a + b , start=list(a=-1.0, b=1.0), control=list(maxiter=500)))  #  non-linear regression
x<-log2(mn)
y<-log2(cv)


#  distribution of CV
par(mar=c(5.5, 8.0, 0.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
h<-hist(y, breaks=100, plot=F)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
XTICK<-pretty(h$mids, 5)
h<-hist(y, breaks=h$breaks, col='grey39', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(XTICK), add=F)
abline(v=CVMIN, lty=2, lwd=4, col='red3', xpd=F)
abline(v=CVMAX, lty=2, lwd=4, col='red3', xpd=F)
mtext(expression(log[2]('CV')), side=1, line=3, padj=+0.5, cex=2.4)
mtext('Frequency', side=2, line=6, padj=+0.2, cex=2.4, las=0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/histogram_cv_genes.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  pick and scatterplot the CV of features versus mean
y.<-log2(predict(nl))
table( pick<-( x>log2(quantile(mn, QPICK)) & y>=CVMIN & y<=CVMAX & abs(y-y.)>=CVDIFF ) )   #  picking features here...
lim1<-pretty(range(x), 5)
lim2<-pretty(range(y), 5)
par(mar=c(5.0, 5.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(x, y, type='p', pch=19, col='darkgrey', cex=0.4, xlab='', ylab='', xlim=range(lim1), ylim=range(lim2))
mtext(expression(log[2]('Mean')), side=1, line=3, padj=+0.3, cex=2.4)
mtext(expression(log[2]('CV')), side=2, line=3, padj=+0.2, cex=2.4, las=0)
n<-order(x, decreasing=F)
points(x[pick], y[pick], pch=19, cex=0.4, col='deepskyblue4')
lines(x[n], y.[n], col='red3', lty=1, lwd=8)
abline(h=0.0, lty=2, lwd=2)
abline(h=CVMIN, lty=3, lwd=6, col='deepskyblue4', xpd=F)
abline(h=CVMAX, lty=3, lwd=6, col='deepskyblue4', xpd=F)
#lines(seq(min(lim1), max(lim1), length.out=100), -0.5*seq(min(lim1), max(lim1), length.out=100), lty=1, lwd=6, col='darkgreen')  #  Poisson noise
mtext('Genes', side=3, line=0, padj=+0.5, cex=1.8, las=0)
legend('topleft', legend=c('fitted model', paste0('selected features (', sum(pick), ')')),  col=c('red3', 'deepskyblue4'), bty='n', lty=c(1, NA), pch=c(NA, 19), pt.cex=c(0, 2.0), lwd=c(15, 0), cex=1.4, y.intersp=0.8, x.intersp=0.5, seg.len=0.5, xpd=NA)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/scatterplot_cv_vs_mean_genes.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  pick the selected genes
gns.cnt<-gns.cnt[ pick, , drop=F]

#}}}


#  [circRNAs] we keep them all since they have been carefully selected
#             append 'circ' prefix to their names
x<-crs.cnt
stopifnot( !any(grepl('^circ', rownames(x))) )  #  sanity check that no gene names starts with 'circ'
rownames(x)<-sub('^', 'circ', rownames(x))


#  finalize the list of features
stopifnot( all.equal( colnames(x), colnames(gns.cnt) ) )
feat<-rbind(gns.cnt, x)
rm(x, gns.cnt, crs.cnt)


#  define the risk classes
stopifnot( all.equal( colnames(feat), meta$bid ) )
risk<-factor(setNames(meta$risk_group, meta$bid), levels=meta[, unique(risk_group)])


#  save 
save(meta, feat, risk, hsa, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/circRNAs+genes_normalized-counts_tumors_filtered.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/circRNAs+genes_normalized-counts_tumors_filtered.RData



#  [genes+circRNAs] self-contained machine learning algorithm that does the feature selection
#{{{
rm(list=ls())
library(rtracklayer)
library(data.table)
library(foreach)
library(doParallel)
library(iterators)
library(cluster)      #  silhouette()
library(randomForest)
source('~/bio/lib/enrichment_analysis.R')
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/feature_selection_functions.R')


#  global parameters
NTREE<-200               #  number of trees to grow in all random forest runs (we need good sampling of the cell lines)
LEAVE_OUT<-25            #  number of samples to leave out for cross-validation
MIN_IMPORTANCE<-3.0      #  minimum importance (MeanDecreaseAccuracy) to ask for selecting one-by-one the important features
MIN_OOB_ERROR<-0.60      #  minimum OOB class error to ask for selecting one-by-one the important features 
MAX_OOB_ERROR<-0.25      #  max OOB error to ask for at least one class during the second round of selection
FREQ_CUTOFF<-0.05        #  minimum frequency that an important feature survived trimmed random forest in order to be considered in the end 



#  load the filtered features
#  log10-transform normalized counts preserving zeros
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/circRNAs+genes_normalized-counts_tumors_filtered.RData')
feat<-log10(1+feat)



#  [step one] random forest on one-by-one all the features using all samples and the OOB error 
#{{{

#  register the cores
registerDoParallel(cores=detectCores())


#  run in parallel
Y<-risk
X<-feat[, names(Y)]
stopifnot(all.equal( names(Y), colnames(X) ) )
#ITRAIN<-sort(sample(ncol(X), ncol(X)-LEAVE_OUT, replace=F))        #  define training set once and for all 
#while (length(setdiff(levels(Y), levels(factor(Y[ITRAIN]))))!=0){  #  make sure all classes are included in the training set
#    ITRAIN<-sort(sample(ncol(X), ncol(X)-LEAVE_OUT, replace=F))   
#}
ITRAIN<-seq_len(ncol(X))
ONE_BY_ONE<-foreach_one_by_one_randomForest(t(X), Y, ITRAIN, NTREE)


#  unregister the cores
registerDoSEQ()


#  save the results
save(ONE_BY_ONE, ITRAIN, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one.RData')

#}}}



#  [step two] identify important classifiers 
#             cluster them by Spearman correlation distance 
#             within each cluster do repeated trimmed random forest collecting errors on randomly drawn test sets
#             identify trimmed sets with up to maximum test classification errors and within them features of at least a minimum frequency of occurrence
#             do repeated trimmed random forest on the identified features
# 
#             MSigDB C2, GO BP/MF enrichment analysis for the final list of unique trimmed features against all features
#{{{

#  explore the one-by-one positive importances to get an idea if the cutoffs above make sense
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one.RData')
imp<-ONE_BY_ONE[ imp>0 ]
hist(imp$imp, breaks=100)
summary(imp$imp)
#
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#  0.000000  3.173723  6.381798  7.435246 10.575484 41.833504 


#  MIN_IMPORTANCE of 3.0 removes approximately the bottom 25% which makes sense
imp<-imp[ imp>=MIN_IMPORTANCE ]
hist(imp$test_err, breaks=100)
dev.off()
summary(imp$test_err)
#
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000000 0.5000000 0.5833333 0.5517803 0.6206897 0.8214286 


#  MIN_OOB_ERROR of 0.6 removes top 25%
imp<-imp[ test_err<=MIN_OOB_ERROR ][ order(-imp, +test_err )]


#  cluster based on Spearman correlations cutting at 2 major clusters
NCLUST<-2
Y<-risk
X<-feat[ rownames(feat) %in% imp[, feature], ]
CO<-cor(t(X), method='spearman')
CO.D<-as.dist(1-CO)
H<-hclust(CO.D, method='ward.D2')  
imp.w<-cutree(H, k=NCLUST)  
print(summary(silhouette(imp.w, dist=CO.D)))
W.check<-cutree(H, k=NCLUST+1)
print(summary(silhouette(W.check, dist=CO.D)))
#cor_heatmap(CO, imp.w, H, new.fig=F)  #  optionally visualize the clusters 
rm(W.check)


#  do trimmed random forest within each cluster
#  identify the features appearing with at least a threshold frequency
imp.kept<-list()
for(i in seq_len(NCLUST)){
    X<-feat[ rownames(feat) %in% names(imp.w)[ imp.w==i ], , drop=F]
    NRUNS<-ceiling(nrow(X)/2) - ceiling(nrow(X)/2) %% 50 + 50
    cat('\n\nstarting the', NRUNS, 'runs...\n')

    #  OOB errors are good enough for this purpose do not leave anything out at this stage
    er<-collect_errors_and_features(NRUNS, X=X, Y=Y, RETURN.FITS=F, FUN=select_features_randomForest, REPEATS=200, NTREE=NTREE, LEAVE.OUT=NULL)

    f<-sort(table(unlist(er$features[ apply(er$errors, 1, function(e){ any(e<=MAX_OOB_ERROR) }) ]))/NRUNS, decreasing=T)
    imp.kept[[i]]<-f[ f>=FREQ_CUTOFF ]
    cat('\n...done!\n')
}


#  use all features that survived the filtering for a last round of repeated trimmed random forest with cross-validation
X<-feat[ unique(names(unlist(imp.kept))), , drop=F]
imp.errors<-collect_errors_and_features(nrow(X), X=X, Y=Y, RETURN.FITS=T, FUN=select_features_randomForest, REPEATS=200, NTREE=NTREE, LEAVE.OUT=LEAVE_OUT) 


#  MSigDB:C2, GO:BP, GO:MF enrichment analysis for the genes (not circRNAs) associated with the important features
#{{{

#  load reference and remove subversions from gene_ids
hsa<-import('/data/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_name', 'gene_id', 'gene_type')]
hsa$gene_id<-sub('\\.[0-9]*$', '', hsa$gene_id)


#  GO/MSigDB C2 enrichment analyses with universe being the whole set of non-circRNA features 
U<-rownames(feat[ grep('^circ', rownames(feat), invert=T), , drop=F])
U<-setNames(hsa$gene_id[ match(U, hsa$gene_name) ], U)
G<-unique(unlist(imp.errors$features))
G<-G[ grep('^circ', G, invert=T) ]
G<-setNames(hsa$gene_id[ match(G, hsa$gene_name) ], G)
imp.enrich<-enrichment_analysis(G, U, MS=setNames('/data/annotation/MSigDB/c2.all.v7.0.symbols.gmt', 'C2'))
rm(hsa, U, G)

#}}}


#  save
save(imp, imp.w, imp.kept, imp.errors, imp.enrich, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one_important.RData')

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one_important.RData



#  [circRNAs] self-contained machine learning algorithm that does the feature selection
#{{{
rm(list=ls())
library(rtracklayer)
library(data.table)
library(foreach)
library(doParallel)
library(iterators)
library(cluster)      #  silhouette()
library(randomForest)
source('~/bio/lib/enrichment_analysis.R')
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/feature_selection_functions.R')


#  global parameters
NTREE<-200               #  number of trees to grow in all random forest runs (we need good sampling of the cell lines)
LEAVE_OUT<-25            #  number of samples to leave out for cross-validation
MIN_IMPORTANCE<-3.0      #  minimum importance (MeanDecreaseAccuracy) to ask for selecting one-by-one the important features
MIN_OOB_ERROR<-0.60      #  minimum OOB class error to ask for selecting one-by-one the important features 
MAX_OOB_ERROR<-0.25      #  max OOB error to ask for at least one class during the second round of selection
FREQ_CUTOFF<-0.05        #  minimum frequency that an important feature survived trimmed random forest in order to be considered in the end 



#  load the filtered features
#  log10-transform normalized counts preserving zeros
#  keep circRNAs only
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/circRNAs+genes_normalized-counts_tumors_filtered.RData')
feat<-log10(1+feat)
feat<-feat[ grep('^circ', rownames(feat)), , drop=T]


#  [step one] random forest on one-by-one all the features using all samples and the OOB error 
#{{{

#  register the cores
registerDoParallel(cores=detectCores())


#  run in parallel
Y<-risk
X<-feat[, names(Y)]
stopifnot(all.equal( names(Y), colnames(X) ) )
#ITRAIN<-sort(sample(ncol(X), ncol(X)-LEAVE_OUT, replace=F))        #  define training set once and for all 
#while (length(setdiff(levels(Y), levels(factor(Y[ITRAIN]))))!=0){  #  make sure all classes are included in the training set
#    ITRAIN<-sort(sample(ncol(X), ncol(X)-LEAVE_OUT, replace=F))   
#}
ITRAIN<-seq_len(ncol(X))
ONE_BY_ONE<-foreach_one_by_one_randomForest(t(X), Y, ITRAIN, NTREE)


#  unregister the cores
registerDoSEQ()


#  save the results
save(ONE_BY_ONE, ITRAIN, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one_circRNAs.RData')

#}}}



#  [step two] identify important classifiers 
#             cluster them by Spearman correlation distance 
#             within each cluster do repeated trimmed random forest collecting errors on randomly drawn test sets
#             identify trimmed sets with up to maximum test classification errors and within them features of at least a minimum frequency of occurrence
#             do repeated trimmed random forest on the identified features
# 
#             MSigDB C2, GO BP/MF enrichment analysis for the final list of unique trimmed features against all features
#{{{

#  explore the one-by-one positive importances to get an idea if the cutoffs above make sense
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one_circRNAs.RData')
imp<-ONE_BY_ONE[ imp>0 ]
hist(imp$imp, breaks=100)
summary(imp$imp)
#
#        Min.     1st Qu.      Median        Mean     3rd Qu.        Max. 
#  0.02275671  3.14147573  6.26025493  6.92226787  9.75941908 25.91795380 


#  MIN_IMPORTANCE of 3.0 removes approximately the bottom 25% which makes sense
imp<-imp[ imp>=MIN_IMPORTANCE ]
hist(imp$test_err, breaks=100)
dev.off()
summary(imp$test_err)
#
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000000 0.2083333 0.3793103 0.3862250 0.5517241 0.8461538 


#  MIN_OOB_ERROR of 0.6 removes top 15%
imp<-imp[ test_err<=MIN_OOB_ERROR ][ order(-imp, +test_err )]


#  cluster based on Spearman correlations (no need to cut at different clusters)
NCLUST<-1
Y<-risk
X<-feat[ rownames(feat) %in% imp[, feature], ]
CO<-cor(t(X), method='spearman')
CO.D<-as.dist(1-CO)
H<-hclust(CO.D, method='ward.D2')  
imp.w<-cutree(H, k=NCLUST)  
print(summary(silhouette(imp.w, dist=CO.D)))
W.check<-cutree(H, k=NCLUST+1)
print(summary(silhouette(W.check, dist=CO.D)))
cor_heatmap(CO, imp.w, H, new.fig=F)  #  optionally visualize the clusters 
rm(W.check)


#  do trimmed random forest within each cluster
#  identify the features appearing with at least a threshold frequency
imp.kept<-list()
for(i in seq_len(NCLUST)){
    X<-feat[ rownames(feat) %in% names(imp.w)[ imp.w==i ], , drop=F]
    NRUNS<-ceiling(nrow(X)/2) - ceiling(nrow(X)/2) %% 50 + 50
    cat('\n\nstarting the', NRUNS, 'runs...\n')

    #  OOB errors are good enough for this purpose do not leave anything out at this stage
    er<-collect_errors_and_features(NRUNS, X=X, Y=Y, RETURN.FITS=F, FUN=select_features_randomForest, REPEATS=200, NTREE=NTREE, LEAVE.OUT=NULL)

    f<-sort(table(unlist(er$features[ apply(er$errors, 1, function(e){ any(e<=MAX_OOB_ERROR) }) ]))/NRUNS, decreasing=T)
    imp.kept[[i]]<-f[ f>=FREQ_CUTOFF ]
    cat('\n...done!\n')
}


#  use all features that survived the filtering for a last round of repeated trimmed random forest with cross-validation
X<-feat[ unique(names(unlist(imp.kept))), , drop=F]
imp.errors<-collect_errors_and_features(nrow(X), X=X, Y=Y, RETURN.FITS=T, FUN=select_features_randomForest, REPEATS=200, NTREE=NTREE, LEAVE.OUT=LEAVE_OUT) 


#  MSigDB:C2, GO:BP, GO:MF enrichment analysis for the genes (not circRNAs) associated with the important features
#{{{

#  load reference and remove subversions from gene_ids
hsa<-import('/data/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_name', 'gene_id', 'gene_type')]
hsa$gene_id<-sub('\\.[0-9]*$', '', hsa$gene_id)


#  GO/MSigDB C2 enrichment analyses with universe being the whole set of non-circRNA features 
U<-sub('^circ', '', rownames(feat))
U<-setNames(hsa$gene_id[ match(U, hsa$gene_name) ], U)
G<-sub('^circ', '', unique(unlist(imp.errors$features)))
G<-setNames(hsa$gene_id[ match(G, hsa$gene_name) ], G)
imp.enrich<-enrichment_analysis(G, U, MS=setNames('/data/annotation/MSigDB/c2.all.v7.0.symbols.gmt', 'C2'))
rm(hsa, U, G)

#}}}


#  save
save(imp, imp.w, imp.kept, imp.errors, imp.enrich, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one_important_circRNAs.RData')

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one_circRNAs.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one_important_circRNAs.RData




#########################################
#
# 
#  analysis of the classification results  
#
# 
#########################################




#  [genes+circRNAs] identify and save the consistently up-/downregulated important features per risk group
#                    do all sorts of visual representations of the results...
#{{{
rm(list=ls())
library(rtracklayer)
library(data.table)
library(pheatmap)
library(gplots)
library(ggplot2)
library(ggrepel)
library(randomForest)
library(cluster)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/draw_highlights.R')
source('~/bio/lib/my_ecdfs.R')
source('~/bio/lib/grouplist2boxplot.R')
source('~/bio/lib/enrichment_analysis.R')
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/feature_selection_functions.R')


#  load the filtered features
#  log10-transform normalized counts preserving zeros
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/circRNAs+genes_normalized-counts_tumors_filtered.RData')
feat<-log10(1+feat)


#  load the one-by-one feature selection results and the final list of important features
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one_important.RData')


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  PCA and hierarchical clustering based on MOST VARIABLE AND HIGHER EXPRESSED features
#{{{

#  pick features based on CV>median(CV) and mean expression>median(expression)
cv<-apply(feat, 1, function(x){ sd(x)/mean(x) })
mn<-rowMeans(feat)


#  PCA 
p<-prcomp(t(feat[ cv>quantile(cv, 0.5) & mn>quantile(mn, 0.5), ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)] )
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
    #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl) +
    scale_color_manual(name='risk group', values=cl) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/plotPCA_most-variable-expressed_features.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(p, ve, pca, cl, XLIM, YLIM)

#}}}


#  PCA, hierarchical and kmeans clustering based on ALL KEPT features
#{{{

#  PCA on the log10-transformed kept features centered but not scaled 
p<-prcomp(t(feat[ unique(names(unlist(imp.kept))), ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)] )
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
    #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl) +
    scale_color_manual(name='risk group', values=cl) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/plotPCA_log10_normalized_counts_kept_features.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(p, ve, pca, cl, XLIM, YLIM)


#  [log10-transformed kept features] hierarchical clustering of samples and features 
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
x<-feat[ unique(names(unlist(imp.kept))), ,drop=F]
d<-dist(t(x), method='euclidean')
hc.c<-hclust(d, method='ward.D2')
d<-dist(x, method='euclidean')
hc.r<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(x), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0, 5, length.out=11),
        cluster_rows=hc.r,
        cluster_cols=hc.c,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
        annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/heatmap_hclust_log10_normalized_counts_kept_features.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc.r, hc.c, ph)


#  [log10-transformed kept features] hierarchical clustering of samples based on Spearman correlations
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
x<-feat[ unique(names(unlist(imp.kept))), ,drop=F]
d<-cor(x, method='spearman', use='pairwise.complete.obs')
hc<-hclust(as.dist(1-d), method='ward.D2')
ph<-pheatmap(as.matrix(d), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0.6, 1, length.out=11),
        cluster_rows=hc,
        cluster_cols=hc,
        annotation_legend=T, annotation_names_row=F, annotation_names_col=F, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/heatmap_hclust_cor_log10_normalized_counts_kept_features.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc, ph)


#  [log10-transformed kept features] kmeans clustering of samples and hierarchical clustering of features
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
x<-feat[ unique(names(unlist(imp.kept))), ,drop=F]
hc.c<-kmeans(t(x), length(cl[['risk_group']]))
x<-x[, order(hc.c$cluster)]
stopifnot( all.equal( rownames(ex), names(hc.c$cluster) ) ) 
ex$cluster<-hc.c$cluster
cl[['cluster']]<-setNames(colorRampPalette(setdiff(brewer.pal(9, 'Greys'), c('#FFFFFF', '#F0F0F0')))(length(unique(hc.c$cluster))), unique(sort(hc.c$cluster)))
d<-dist(x, method='euclidean')
hc.r<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(x), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0, 5, length.out=11),
        cluster_rows=hc.r,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=F, annotation_names_col=F, 
        annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/heatmap_kmeans_log10_normalized_counts_kept_features.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, hc.c, hc.r, d, ph)

#}}}


#  [run once] identify and plot features that consistently classify a risk group
#             MSigDB C2, GO MF/BP enrichment analysis 
#             collect performance statistics
#             save
#{{{

#  for each risk group identify the consistently high performing feature sets 
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)] )
imp.up<-imp.down<-list()
for( r in names(cl)){
    #  select the features
    f<-imp.errors$errors[, r]<=0.05
    f<-data.frame(table(unlist(imp.errors$features[ imp.errors$errors[, r]<=0.05 ]))/sum(f))
    #f<-as.character(f[ f$Freq>quantile(f$Freq, 0.10), 1])
    f<-as.character(f[, 1])


    #  boxplot of expression distributions for those features with distributions that stand out for the corresponding risk group
    #
    #  p-value cutoff = 0.05
    #
    pv<-expression_boxplots(X=feat[f, , drop=F], Y=risk, main='', trim.chr=T, 
               p.values=list(x=r, y=setdiff(names(cl), r)), pv.cutoff=0.05, pv.as.groups=F, 
               legend.place='topright', yline=3, xline=0, xlas=2, xadj=1,
               figure.name=NULL, svg.width=16, colorize=NULL, cluster.colors=cl, 
               mar=c(11.5, 6.0, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
    imp.up[[r]]<-pv$greater[ apply(pv$greater, 1, function(p){ all(p<0.05) }), , drop=F]
    imp.down[[r]]<-pv$less[ apply(pv$less, 1, function(p){ all(p<0.05) }), , drop=F]
    readline('Press ENTER:')


    #  PCA on ALL selected features
    f<-c(rownames(imp.up[[r]]), rownames(imp.down[[r]]))
    if(length(f)>0){
        p<-prcomp(t(feat[ f, ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when samples are the rows
        ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
        pca<-as.data.frame(p$x)
        pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
        XLIM<-pretty(range(p$x[, 'PC1']), 5)
        XLIM<-c(XLIM[1], XLIM[length(XLIM)])
        YLIM<-pretty(range(p$x[, 'PC2']), 5)
        YLIM<-c(YLIM[1], YLIM[length(YLIM)])
        g<-ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
            geom_point(size=10) + 
            #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
            scale_shape_manual(values=c(18:21)) + 
            scale_fill_manual(name='risk group', values=cl) +
            scale_color_manual(name='risk group', values=cl) + 
            theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
                axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
                legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
                aspect.ratio=1, panel.background = element_blank(),
                axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
                axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
            scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
            xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
        print(g)
        readline('Press ENTER:')
    }
}


#  MSigDB C2 and GO BP/MF enrichment analysis for the identified features per risk group
h<-mcols(hsa)[, c('gene_id', 'gene_name')]
h$gene_id<-sub('\\.[0-9]*$', '', h$gene_id)
U<-rownames(feat[ grep('^circ', rownames(feat), invert=T), , drop=F])
U<-setNames(h$gene_id[ match(U, h$gene_name) ], U)
imp.up.enrich<-lapply(imp.up, function(r){ 
    G<-rownames(r)
    G<-G[ grep('^circ', G, invert=T) ]
    G<-setNames(h$gene_id[ match(G, h$gene_name) ], G)
    enrichment_analysis(G, U, MS=setNames('/data/annotation/MSigDB/c2.all.v7.0.symbols.gmt', 'C2'))
})
imp.down.enrich<-lapply(imp.down, function(r){ 
    G<-rownames(r)
    G<-G[ grep('^circ', G, invert=T) ]
    G<-setNames(h$gene_id[ match(G, h$gene_name) ], G)
    enrichment_analysis(G, U, MS=setNames('/data/annotation/MSigDB/c2.all.v7.0.symbols.gmt', 'C2'))
})
rm(U, h)


#  collect performance statistics
#{{{

#  constants
NTREE<-200
NRUNS<-200
LEAVE_OUT<-25


#  containers that will hold the results for all risk groups
imp.up.performance<-imp.down.performance<-list()


#  collect the statistics
for(n in seq_len(NRUNS)){
    cat('Starting iteration:', n)

    #  define same training set for all risk groups
    #  make sure to include all classes
    ITRAIN<-sort(sample(length(risk), length(risk)-LEAVE_OUT, replace=F))
    while (length(setdiff(levels(risk), levels(factor(risk[ITRAIN]))))!=0){
        ITRAIN<-sort(sample(length(risk), length(risk)-LEAVE_OUT, replace=F))   
    }


    #  go over the risk groups and do a single run
    for(r in names(imp.up)){
        X<-feat[ rownames(imp.up[[r]]), , drop=F]
        imp.up.performance[[r]][[n]]<-collect_performance(X=X, Y=risk, NTREE=NTREE, ITRAIN=ITRAIN)
        X<-feat[ rownames(imp.down[[r]]), , drop=F]
        imp.down.performance[[r]][[n]]<-collect_performance(X=X, Y=risk, NTREE=NTREE, ITRAIN=ITRAIN)
    }
    cat('                             \r')
}
cat('Done!\n')

#}}}


#
#  save
#
save(imp.up, imp.up.enrich, imp.down, imp.down.enrich, imp.up.performance, imp.down.performance, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_identified_per_risk-group.RData')

#}}}


#  load consistently up-/downregulated identified features and their performances
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_identified_per_risk-group.RData')


#  PCA and hierarchical clustering based on ONLY HIGHER EXPRESSED features in a given risk group
#{{{

#  pool 
length(keep<-unique(unlist(lapply(imp.up, rownames))))  #  194


#  PCA
p<-prcomp(t(feat[ keep, ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)] )
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
    #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl) +
    scale_color_manual(name='risk group', values=cl) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/plotPCA_log10_normalized_counts_identified_features_up.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(p, ve, pca, cl, XLIM, YLIM)


#  hierarchical clustering of samples and features
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
x<-feat[ keep, ,drop=F]
d<-dist(t(x), method='euclidean')
hc.c<-hclust(d, method='ward.D2')
d<-dist(x, method='euclidean')
hc.r<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(x), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0, 5, length.out=11),
        cluster_rows=hc.r,
        cluster_cols=hc.c,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
        annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/heatmap_hclust_log10_normalized_counts_identified_features_up.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc.r, hc.c, ph)


#  expression distribution boxplots
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
for(r in names(cl$risk_group)){
    pv<-expression_boxplots(X=feat[rownames(imp.up[[r]]), , drop=F], Y=risk, main='', trim.chr=T, 
                   p.values=list(x=r, y=setdiff(names(cl$risk_group), r)), pv.cutoff=0.05, pv.as.groups=F, 
                   legend.place='topright', yline=3, xline=0, xlas=2, xadj=1, xcex=1.6, ycex=2.4, 
                   figure.name=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_', r, '_up.svg'), svg.width=min(65, max(16, nrow(imp.up[[r]]))), colorize=NULL, cluster.colors=cl$risk_group, 
                   mar=c(10.0, 6.5, 0.5, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
    readline('Press ENTER:')


    ph<-pheatmap(as.matrix(feat[ rownames(imp.up[[r]]), , drop=F]), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
            cluster_rows=F,
            cluster_cols=F,
            annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
            annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
            drop_levels=F, show_rownames=F, show_colnames=F, 
            fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
    readline('Press ENTER:')


    #  PCA
    p<-prcomp(t(feat[ rownames(imp.up[[r]]), ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
    ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
    pca<-as.data.frame(p$x)
    pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
    XLIM<-pretty(range(p$x[, 'PC1']), 5)
    XLIM<-c(XLIM[1], XLIM[length(XLIM)])
    YLIM<-pretty(range(p$x[, 'PC2']), 5)
    YLIM<-c(YLIM[1], YLIM[length(YLIM)])
    g<-ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
        geom_point(size=10) + 
        scale_shape_manual(values=c(18:21)) + 
        scale_fill_manual(name='risk group', values=cl$risk_group) +
        scale_color_manual(name='risk group', values=cl$risk_group) + 
        theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
            axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
            legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
            aspect.ratio=1, panel.background = element_blank(),
            axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
            axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
        scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
        xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
    print(g)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/plotPCA_identified_features_', r, '_up.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    readline('Press ENTER:')
}

#}}}


#  PCA and hierarchical clustering based on ONLY LOWER EXPRESSED features in a given risk group
#{{{

#  pool 
length(keep<-unique(unlist(lapply(imp.down, rownames))))  #  279


#  PCA
p<-prcomp(t(feat[ keep, ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)] )
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
    #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl) +
    scale_color_manual(name='risk group', values=cl) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/plotPCA_log10_normalized_counts_identified_features_down.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(p, ve, pca, cl, XLIM, YLIM)


#  hierarchical clustering of samples and features
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
x<-feat[ keep, ,drop=F]
d<-dist(t(x), method='euclidean')
hc.c<-hclust(d, method='ward.D2')
d<-dist(x, method='euclidean')
hc.r<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(x), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0, 5, length.out=11),
        cluster_rows=hc.r,
        cluster_cols=hc.c,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
        annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/heatmap_hclust_log10_normalized_counts_identified_features_down.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc.r, hc.c, ph)


#  expression distribution boxplots
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
for(r in names(cl$risk_group)){
    pv<-expression_boxplots(X=feat[rownames(imp.down[[r]]), , drop=F], Y=risk, main='', trim.chr=T, 
                   p.values=list(x=r, y=setdiff(names(cl$risk_group), r)), pv.cutoff=0.05, pv.as.groups=F, 
                   legend.place='topright', yline=3, xline=0, xlas=2, xadj=1, xcex=1.6, ycex=2.4, 
                   figure.name=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_', r, '_down.svg'), svg.width=min(65, max(16, nrow(imp.down[[r]]))), colorize=NULL, cluster.colors=cl$risk_group, 
                   mar=c(10.0, 6.5, 0.5, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
    readline('Press ENTER:')


    ph<-pheatmap(as.matrix(feat[ rownames(imp.down[[r]]), , drop=F]), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
            cluster_rows=F,
            cluster_cols=F,
            annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
            annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
            drop_levels=F, show_rownames=F, show_colnames=F, 
            fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
    readline('Press ENTER:')


    #  PCA
    p<-prcomp(t(feat[ rownames(imp.down[[r]]), ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
    ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
    pca<-as.data.frame(p$x)
    pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
    XLIM<-pretty(range(p$x[, 'PC1']), 5)
    XLIM<-c(XLIM[1], XLIM[length(XLIM)])
    YLIM<-pretty(range(p$x[, 'PC2']), 5)
    YLIM<-c(YLIM[1], YLIM[length(YLIM)])
    g<-ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
        geom_point(size=10) + 
        scale_shape_manual(values=c(18:21)) + 
        scale_fill_manual(name='risk group', values=cl$risk_group) +
        scale_color_manual(name='risk group', values=cl$risk_group) + 
        theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
            axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
            legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
            aspect.ratio=1, panel.background = element_blank(),
            axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
            axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
        scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
        xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
    print(g)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/plotPCA_identified_features_', r, '_down.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    readline('Press ENTER:')
}

#}}}


#  [kept features] distributions of cross-validation errors
#{{{

B.cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)])
B<-imp.errors$errors
YTICK<-pretty(c(0, max(B)), 5)
par(mar=c(10.0, 8.0, 0.5, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(1, ncol(B))+c(-0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Prediction error (%)', side=2, line=5, padj=-0.2, las=0, cex=2.4)
mtext(text=colnames(B), at=seq_len(ncol(B)), side=1, line=0, las=2, padj=0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_cross-validation_errors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B, B.cl, YTICK, bp)

#}}}


#  [risk group-specific features] accuracy and precision plots 
#{{{

#  [up] accuracy = (TP+TN)/(P+N)
x<-lapply(imp.up.performance, function(r){ 100*do.call(rbind, lapply(r, function(m){ (m[, 'tp']+m[, 'tn'])/rowSums(m) })) })
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)])
L<-lapply(x, t)
grouplist2boxplot(L, L.COL=cl, YLAB='Accuracy (%)', XTEXT.LINE=0, XTEXT.ADJ=1, LEGEND='bottomleft', mar=c(8.0, 5.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_up_accuracy.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x, L, cl)


#  [down] accuracy = (TP+TN)/(P+N)
x<-lapply(imp.down.performance, function(r){ 100*do.call(rbind, lapply(r, function(m){ (m[, 'tp']+m[, 'tn'])/rowSums(m) })) })
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)])
L<-lapply(x, t)
grouplist2boxplot(L, L.COL=cl, YLAB='Accuracy (%)', XTEXT.LINE=0, XTEXT.ADJ=1, LEGEND='bottomleft', mar=c(8.0, 5.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_down_accuracy.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x, L, cl)


#  [up] precision = 1 - FDR = TP/(TP+FP)
#  N.B. NaN means 0/0, i.e. all cases have been classified as negative, if FN is also zero then precision = 1, otherwise precision = 0
x<-lapply(imp.up.performance, function(r){ 100*do.call(rbind, lapply(r, function(m){
    a<-setNames(m[, 'tp']/(m[, 'tp']+m[, 'fp']), rownames(m))
    a[is.na(a) & m[, 'fn']==0]<-1
    a[is.na(a) & m[, 'fn']!=0]<-0
    return(a)
}))})
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)])
L<-lapply(x, t)
grouplist2boxplot(L, L.COL=cl, YLAB='Precision (%)', XTEXT.LINE=0, XTEXT.ADJ=1, LEGEND='bottomright', mar=c(8.0, 5.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_up_precision.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x, L, cl)


#  [down] precision = 1 - FDR = TP/(TP+FP)
#  N.B. NaN means 0/0, i.e. all cases have been classified as negative, if FN is also zero then precision = 1, otherwise precision = 0
x<-lapply(imp.down.performance, function(r){ 100*do.call(rbind, lapply(r, function(m){
    a<-setNames(m[, 'tp']/(m[, 'tp']+m[, 'fp']), rownames(m))
    a[is.na(a) & m[, 'fn']==0]<-1
    a[is.na(a) & m[, 'fn']!=0]<-0
    return(a)
}))})
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)])
L<-lapply(x, t)
grouplist2boxplot(L, L.COL=cl, YLAB='Precision (%)', XTEXT.LINE=0, XTEXT.ADJ=1, LEGEND='bottomright', mar=c(8.0, 5.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_down_precision.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x, L, cl)

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_identified_per_risk-group.RData



#  [circRNAs] identify and save the consistently up-/downregulated important features per risk group
#             do all sorts of visual representations of the results...
#{{{
rm(list=ls())
library(rtracklayer)
library(data.table)
library(pheatmap)
library(gplots)
library(ggplot2)
library(ggrepel)
library(randomForest)
library(cluster)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/draw_highlights.R')
source('~/bio/lib/my_ecdfs.R')
source('~/bio/lib/grouplist2boxplot.R')
source('~/bio/lib/enrichment_analysis.R')
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/feature_selection_functions.R')


#  load the filtered features
#  log10-transform normalized counts preserving zeros
#  keep circRNAs only
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/circRNAs+genes_normalized-counts_tumors_filtered.RData')
feat<-log10(1+feat)
feat<-feat[ grep('^circ', rownames(feat)), , drop=T]


#  load the one-by-one feature selection results and the final list of important features
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one_circRNAs.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one_important_circRNAs.RData')


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  PCA, hierarchical and kmeans clustering based on all kept features
#{{{

#  PCA on the log10-transformed kept features centered but not scaled 
p<-prcomp(t(feat[ unique(names(unlist(imp.kept))), ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)] )
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
    #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl) +
    scale_color_manual(name='risk group', values=cl) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/plotPCA_log10_normalized_counts_kept_features_circRNAs.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(p, ve, pca, cl, XLIM, YLIM)


#  [log10-transformed kept features] hierarchical clustering of samples and features 
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
x<-feat[ unique(names(unlist(imp.kept))), ,drop=F]
d<-dist(t(x), method='euclidean')
hc.c<-hclust(d, method='ward.D2')
d<-dist(x, method='euclidean')
hc.r<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(x), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0, 4, length.out=11),
        cluster_rows=hc.r,
        cluster_cols=hc.c,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
        annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/heatmap_hclust_log10_normalized_counts_kept_features_circRNAs.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc.r, hc.c, ph)


#  [log10-transformed kept features] hierarchical clustering of samples based on Spearman correlations
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
x<-feat[ unique(names(unlist(imp.kept))), ,drop=F]
d<-cor(x, method='spearman', use='pairwise.complete.obs')
hc<-hclust(as.dist(1-d), method='ward.D2')
ph<-pheatmap(as.matrix(d), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0.0, 1, length.out=11),
        cluster_rows=hc,
        cluster_cols=hc,
        annotation_legend=T, annotation_names_row=F, annotation_names_col=F, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/heatmap_hclust_cor_log10_normalized_counts_kept_features_circRNAs.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc, ph)


#  [log10-transformed kept features] kmeans clustering of samples and hierarchical clustering of features
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
x<-feat[ unique(names(unlist(imp.kept))), ,drop=F]
hc.c<-kmeans(t(x), length(cl[['risk_group']]))
x<-x[, order(hc.c$cluster)]
stopifnot( all.equal( rownames(ex), names(hc.c$cluster) ) ) 
ex$cluster<-hc.c$cluster
cl[['cluster']]<-setNames(colorRampPalette(setdiff(brewer.pal(9, 'Greys'), c('#FFFFFF', '#F0F0F0')))(length(unique(hc.c$cluster))), unique(sort(hc.c$cluster)))
d<-dist(x, method='euclidean')
hc.r<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(x), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0, 4, length.out=11),
        cluster_rows=hc.r,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=F, annotation_names_col=F, 
        annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/heatmap_kmeans_log10_normalized_counts_kept_features_circRNAs.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, hc.c, hc.r, d, ph)

#}}}


#  [run once] identify and plot features that consistently classify a risk group
#             MSigDB C2, GO MF/BP enrichment analysis 
#             collect performance statistics
#             save
#{{{

#  for each risk group identify the consistently high performing feature sets 
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)] )
imp.up<-imp.down<-list()
for( r in names(cl)){
    #  select the features
    f<-imp.errors$errors[, r]<=0.05
    f<-data.frame(table(unlist(imp.errors$features[ imp.errors$errors[, r]<=0.05 ]))/sum(f))
    #f<-as.character(f[ f$Freq>quantile(f$Freq, 0.10), 1])
    f<-as.character(f[, 1])


    #  boxplot of expression distributions for those features with distributions that stand out for the corresponding risk group
    #
    #  p-value cutoff = 0.1
    #
    pv<-expression_boxplots(X=feat[f, , drop=F], Y=risk, main='', trim.chr=T, 
               p.values=list(x=r, y=setdiff(names(cl), r)), pv.cutoff=0.1, pv.as.groups=F, 
               legend.place='topright', yline=5, xline=0, xlas=2, xadj=1, xcex=1.6, ycex=2.4, 
               figure.name=NULL, svg.width=16, colorize=NULL, cluster.colors=cl, 
               mar=c(10.0, 8.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
    imp.up[[r]]<-pv$greater[ apply(pv$greater, 1, function(p){ all(p<0.1) }), , drop=F]
    imp.down[[r]]<-pv$less[ apply(pv$less, 1, function(p){ all(p<0.1) }), , drop=F]
    readline('Press ENTER:')


    #  PCA on ALL selected features
    f<-c(rownames(imp.up[[r]]), rownames(imp.down[[r]]))
    if(length(f)>1){
        p<-prcomp(t(feat[ f, ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when samples are the rows
        ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
        pca<-as.data.frame(p$x)
        pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
        XLIM<-pretty(range(p$x[, 'PC1']), 5)
        XLIM<-c(XLIM[1], XLIM[length(XLIM)])
        YLIM<-pretty(range(p$x[, 'PC2']), 5)
        YLIM<-c(YLIM[1], YLIM[length(YLIM)])
        g<-ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
            geom_point(size=10) + 
            #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
            scale_shape_manual(values=c(18:21)) + 
            scale_fill_manual(name='risk group', values=cl) +
            scale_color_manual(name='risk group', values=cl) + 
            theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
                axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
                legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
                aspect.ratio=1, panel.background = element_blank(),
                axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
                axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
            scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
            xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
        print(g)
        readline('Press ENTER:')
    }
}


#  MSigDB C2 and GO BP/MF enrichment analysis for the identified features per risk group
h<-mcols(hsa)[, c('gene_id', 'gene_name')]
h$gene_id<-sub('\\.[0-9]*$', '', h$gene_id)
U<-sub('^circ', '', rownames(feat))
U<-setNames(h$gene_id[ match(U, h$gene_name) ], U)
imp.up.enrich<-lapply(imp.up, function(r){ 
    G<-sub('^circ', '', rownames(r))
    G<-setNames(h$gene_id[ match(G, h$gene_name) ], G)
    enrichment_analysis(G, U, MS=setNames('/data/annotation/MSigDB/c2.all.v7.0.symbols.gmt', 'C2'))
})
imp.down.enrich<-lapply(imp.down, function(r){ 
    G<-sub('^circ', '', rownames(r))
    G<-setNames(h$gene_id[ match(G, h$gene_name) ], G)
    enrichment_analysis(G, U, MS=setNames('/data/annotation/MSigDB/c2.all.v7.0.symbols.gmt', 'C2'))
})
rm(U, h)


#  collect performance statistics
#{{{

#  constants
NTREE<-200
NRUNS<-200
LEAVE_OUT<-25


#  containers that will hold the results for all risk groups
imp.up.performance<-imp.down.performance<-list()


#  collect the statistics
for(n in seq_len(NRUNS)){
    cat('Starting iteration:', n)

    #  define same training set for all risk groups
    #  make sure to include all classes
    ITRAIN<-sort(sample(length(risk), length(risk)-LEAVE_OUT, replace=F))
    while (length(setdiff(levels(risk), levels(factor(risk[ITRAIN]))))!=0){
        ITRAIN<-sort(sample(length(risk), length(risk)-LEAVE_OUT, replace=F))   
    }


    #  go over the risk groups and do a single run
    for(r in names(imp.up)){
        X<-feat[ rownames(imp.up[[r]]), , drop=F]
        if(nrow(X)>0){
            imp.up.performance[[r]][[n]]<-collect_performance(X=X, Y=risk, NTREE=NTREE, ITRAIN=ITRAIN)
        }
        X<-feat[ rownames(imp.down[[r]]), , drop=F]
        if(nrow(X)>0){
            imp.down.performance[[r]][[n]]<-collect_performance(X=X, Y=risk, NTREE=NTREE, ITRAIN=ITRAIN)
        }
    }
    cat('                             \r')
}
cat('Done!\n')

#}}}


#
#  save
#
save(imp.up, imp.up.enrich, imp.down, imp.down.enrich, imp.up.performance, imp.down.performance, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_identified_per_risk-group_circRNAs.RData')

#}}}


#  load consistently up-/downregulated identified features and their performances
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_identified_per_risk-group_circRNAs.RData')


#  PCA and hierarchical clustering based on ONLY HIGHER EXPRESSED features in a given risk group
#{{{

#  pool 
length(keep<-unique(unlist(lapply(imp.up, rownames))))  #  42


#  PCA
p<-prcomp(t(feat[ keep, ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)] )
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
    #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl) +
    scale_color_manual(name='risk group', values=cl) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/plotPCA_log10_normalized_counts_identified_features_up_circRNAs.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(p, ve, pca, cl, XLIM, YLIM)


#  hierarchical clustering of samples and features
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
x<-feat[ keep, ,drop=F]
d<-dist(t(x), method='euclidean')
hc.c<-hclust(d, method='ward.D2')
d<-dist(x, method='euclidean')
hc.r<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(x), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0, 3.5, length.out=11),
        cluster_rows=hc.r,
        cluster_cols=hc.c,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
        annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/heatmap_hclust_log10_normalized_counts_identified_features_up_circRNAs.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc.r, hc.c, ph)


#  expression distribution boxplots
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
for(r in names(cl$risk_group)){
    f<-rownames(imp.up[[r]])


    if(length(f)>1){
        pv<-expression_boxplots(X=feat[f, , drop=F], Y=risk, main='', trim.chr=T, 
                       p.values=list(x=r, y=setdiff(names(cl$risk_group), r)), pv.cutoff=0.1, pv.as.groups=F, 
                       legend.place='topright', yline=5, xline=0, xlas=2, xadj=1, xcex=1.6, ycex=2.4, 
                       figure.name=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_', r, '_up_circRNAs.svg'), svg.width=min(65, max(16, nrow(imp.up[[r]]))), colorize=NULL, cluster.colors=cl$risk_group, 
                       mar=c(10.0, 8.5, 0.5, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
        readline('Press ENTER:')


        ph<-pheatmap(as.matrix(feat[f, , drop=F]), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
                cluster_rows=F,
                cluster_cols=F,
                annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
                annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
                drop_levels=F, show_rownames=F, show_colnames=F, 
                fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
        readline('Press ENTER:')


        #  PCA
        p<-prcomp(t(feat[f, ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
        ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
        pca<-as.data.frame(p$x)
        pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
        XLIM<-pretty(range(p$x[, 'PC1']), 5)
        XLIM<-c(XLIM[1], XLIM[length(XLIM)])
        YLIM<-pretty(range(p$x[, 'PC2']), 5)
        YLIM<-c(YLIM[1], YLIM[length(YLIM)])
        g<-ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
            geom_point(size=10) + 
            scale_shape_manual(values=c(18:21)) + 
            scale_fill_manual(name='risk group', values=cl$risk_group) +
            scale_color_manual(name='risk group', values=cl$risk_group) + 
            theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
                axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
                legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
                aspect.ratio=1, panel.background = element_blank(),
                axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
                axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
            scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
            xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
        print(g)
        dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/plotPCA_identified_features_', r, '_up_circRNAs.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
        readline('Press ENTER:')
    }
}

#}}}


#  PCA and hierarchical clustering based on ONLY LOWER EXPRESSED features in a given risk group
#{{{

#  pool 
length(keep<-unique(unlist(lapply(imp.down, rownames))))  #  84


#  PCA
p<-prcomp(t(feat[ keep, ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)] )
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
    #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl) +
    scale_color_manual(name='risk group', values=cl) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/plotPCA_log10_normalized_counts_identified_features_down_circRNAs.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(p, ve, pca, cl, XLIM, YLIM)


#  hierarchical clustering of samples and features
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
x<-feat[ keep, ,drop=F]
d<-dist(t(x), method='euclidean')
hc.c<-hclust(d, method='ward.D2')
d<-dist(x, method='euclidean')
hc.r<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(x), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0, 3.5, length.out=11),
        cluster_rows=hc.r,
        cluster_cols=hc.c,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
        annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/heatmap_hclust_log10_normalized_counts_identified_features_down_circRNAs.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc.r, hc.c, ph)


#  expression distribution boxplots
ex<-data.frame(risk_group=factor(meta[ match(colnames(feat), bid), risk_group], exclude=F), row.names=colnames(feat))
cl<-setNames(list(setNames(meta[, unique(col)], meta[, unique(risk_group)] )), 'risk_group')
for(r in names(cl$risk_group)){
    f<-rownames(imp.down[[r]])


    if(length(f)>1){
        pv<-expression_boxplots(X=feat[f, , drop=F], Y=risk, main='', trim.chr=T, 
                       p.values=list(x=r, y=setdiff(names(cl$risk_group), r)), pv.cutoff=0.1, pv.as.groups=F, 
                       legend.place='topright', yline=5, xline=0, xlas=2, xadj=1, xcex=1.6, ycex=2.4, 
                       figure.name=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_', r, '_down_circRNAs.svg'), svg.width=min(65, max(16, nrow(imp.down[[r]]))), colorize=NULL, cluster.colors=cl$risk_group, 
                       mar=c(10.0, 8.5, 0.5, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
        readline('Press ENTER:')


        ph<-pheatmap(as.matrix(feat[f, , drop=F]), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
                cluster_rows=F,
                cluster_cols=F,
                annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
                annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
                drop_levels=F, show_rownames=F, show_colnames=F, 
                fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
        readline('Press ENTER:')


        #  PCA
        p<-prcomp(t(feat[f, ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
        ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
        pca<-as.data.frame(p$x)
        pca$'risk group'<-factor(meta[ match( rownames(pca), bid ), risk_group ], levels=meta[, unique(risk_group)])
        XLIM<-pretty(range(p$x[, 'PC1']), 5)
        XLIM<-c(XLIM[1], XLIM[length(XLIM)])
        YLIM<-pretty(range(p$x[, 'PC2']), 5)
        YLIM<-c(YLIM[1], YLIM[length(YLIM)])
        g<-ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
            geom_point(size=10) + 
            scale_shape_manual(values=c(18:21)) + 
            scale_fill_manual(name='risk group', values=cl$risk_group) +
            scale_color_manual(name='risk group', values=cl$risk_group) + 
            theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
                axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
                legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
                aspect.ratio=1, panel.background = element_blank(),
                axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
                axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
            scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
            xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
        print(g)
        dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/plotPCA_identified_features_', r, '_down_circRNAs.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
        readline('Press ENTER:')
    }
}

#}}}


#  [kept features] distributions of cross-validation errors
#{{{

B.cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)])
B<-imp.errors$errors
YTICK<-pretty(c(0, max(B)), 5)
par(mar=c(10.0, 8.0, 0.5, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(1, ncol(B))+c(-0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Prediction error (%)', side=2, line=5, padj=-0.2, las=0, cex=2.4)
mtext(text=colnames(B), at=seq_len(ncol(B)), side=1, line=0, las=2, padj=0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_cross-validation_errors_circRNAs.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B, B.cl, YTICK, bp)

#}}}


#  [risk group-specific features] accuracy and precision plots 
#{{{

#  [up] accuracy = (TP+TN)/(P+N)
x<-lapply(imp.up.performance, function(r){ 100*do.call(rbind, lapply(r, function(m){ (m[, 'tp']+m[, 'tn'])/rowSums(m) })) })
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)])
for(r in setdiff(names(cl), names(x))){  #  if we have missing risk groups we need to add them with NA values
    x[[r]]<-matrix(rep(NA, length(cl)), ncol=length(cl), dimnames=list('', names(cl)))
}
x<-x[ names(cl) ]
L<-lapply(x, t)
grouplist2boxplot(L, L.COL=cl, YLAB='Accuracy (%)', XTEXT.LINE=0, XTEXT.ADJ=1, LEGEND='bottomleft', mar=c(8.0, 5.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_up_accuracy_circRNAs.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x, L, cl)


#  [down] accuracy = (TP+TN)/(P+N)
x<-lapply(imp.down.performance, function(r){ 100*do.call(rbind, lapply(r, function(m){ (m[, 'tp']+m[, 'tn'])/rowSums(m) })) })
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)])
for(r in setdiff(names(cl), names(x))){  #  if we have missing risk groups we need to add them with NA values
    x[[r]]<-matrix(rep(NA, length(cl)), ncol=length(cl), dimnames=list('', names(cl)))
}
x<-x[ names(cl) ]
L<-lapply(x, t)
grouplist2boxplot(L, L.COL=cl, YLAB='Accuracy (%)', XTEXT.LINE=0, XTEXT.ADJ=1, LEGEND='bottomleft', mar=c(8.0, 5.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_down_accuracy_circRNAs.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x, L, cl)


#  [up] precision = 1 - FDR = TP/(TP+FP)
#  N.B. NaN means 0/0, i.e. all cases have been classified as negative, if FN is also zero then precision = 1, otherwise precision = 0
x<-lapply(imp.up.performance, function(r){ 100*do.call(rbind, lapply(r, function(m){
    a<-setNames(m[, 'tp']/(m[, 'tp']+m[, 'fp']), rownames(m))
    a[is.na(a) & m[, 'fn']==0]<-1
    a[is.na(a) & m[, 'fn']!=0]<-0
    return(a)
}))})
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)])
L<-lapply(x, t)
grouplist2boxplot(L, L.COL=cl, YLAB='Precision (%)', XTEXT.LINE=0, XTEXT.ADJ=1, LEGEND='bottomright', mar=c(8.0, 5.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_up_precision_circRNAs.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x, L, cl)


#  [down] precision = 1 - FDR = TP/(TP+FP)
#  N.B. NaN means 0/0, i.e. all cases have been classified as negative, if FN is also zero then precision = 1, otherwise precision = 0
x<-lapply(imp.down.performance, function(r){ 100*do.call(rbind, lapply(r, function(m){
    a<-setNames(m[, 'tp']/(m[, 'tp']+m[, 'fp']), rownames(m))
    a[is.na(a) & m[, 'fn']==0]<-1
    a[is.na(a) & m[, 'fn']!=0]<-0
    return(a)
}))})
cl<-setNames(meta[, unique(col)], meta[, unique(risk_group)])
for(r in setdiff(names(cl), names(x))){  #  if we have missing risk groups we need to add them with NA values
    x[[r]]<-matrix(rep(NA, length(cl)), ncol=length(cl), dimnames=list('', names(cl)))
}
x<-x[ names(cl) ]
L<-lapply(x, t)
grouplist2boxplot(L, L.COL=cl, YLAB='Precision (%)', XTEXT.LINE=0, XTEXT.ADJ=1, LEGEND='topleft', mar=c(8.0, 5.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/boxplot_identified_features_down_precision_circRNAs.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x, L, cl)

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_identified_per_risk-group_circRNAs.RData




###########################################################################################
#
# 
#  classification of other neuroblastoma cohorts based on features identified in our cohort
#
# 
###########################################################################################




#  Zhang and Hertwig 2015 SEQC neuroblastomas
#{{{
rm(list=ls())
library(Biobase)
library(rtracklayer)
library(data.table)
library(pheatmap)
library(gplots)
library(ggplot2)
library(ggrepel)
library(randomForest)
library(cluster)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/draw_highlights.R')
source('~/bio/lib/grouplist2boxplot.R')
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/feature_selection_functions.R')


#  [run once] process the metadata and the counts
#             define the risk groups according to the scheme:
#             
#                    ST4S : stage 4s + mycn_status = 0 (or NA)
#                      LR : (stage 1 OR stage 2 OR (stage 3 AND age <547 days)) AND mycn_status = 0 (or NA)
#                     IMR : ((stage 3 AND age >= 547) OR (stage 4 AND age <547 days)) AND mycn_status = 0 (or NA)
#                 HR_nMNA : stage 4 AND age >= 547 days AND mycn_status = 0 (or NA)
#                     MNA : mycn_status = 1 (any stage)
#{{{

#  load the metadata
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/Zhang_Hertwig2015_SEQC_metadata.RData')
zhang.meta<-data.table(sid=rownames(pd), pd)
rm(pd)


#  define the risk groups
zhang.meta[ (inss_stage %in% '4s') & (mycn_status %in% 0 | is.na(mycn_status)), risk_group:='ST4S' ]
zhang.meta[ ((inss_stage %in% c('1', '2')) | (inss_stage %in% '3' & age_at_diagnosis<547)) & (mycn_status %in% 0 | is.na(mycn_status)), risk_group:='LR' ]
zhang.meta[ ((inss_stage %in% '3' & age_at_diagnosis>=547) | (inss_stage %in% '4' & age_at_diagnosis<547)) & (mycn_status %in% 0 | is.na(mycn_status)), risk_group:='IMR']
zhang.meta[ (inss_stage %in% '4' & age_at_diagnosis>=547) & (mycn_status %in% 0 | is.na(mycn_status)), risk_group:='HR_nMNA']
zhang.meta[ mycn_status %in% 1, risk_group:='MNA']


#  order them by risk
zhang.meta<-zhang.meta[, risk_group:=factor(risk_group, levels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))][ order(risk_group, sid) ][, risk_group:=as.character(risk_group)]


#  add colors
zhang.meta[ risk_group %in% 'ST4S', col:='seagreen4' ] 
zhang.meta[ risk_group %in% 'LR', col:='darkgreen' ] 
zhang.meta[ risk_group %in% 'IMR', col:='cornflowerblue' ] 
zhang.meta[ risk_group %in% 'HR_nMNA', col:='chocolate1' ] 
zhang.meta[ risk_group %in% 'MNA', col:='coral4' ] 


#  load library-size normalized and log2-transformed counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/Zhang_Hertwig2015_SEQC_counts.RData')
zhang.feat<-exprs(normX)
stopifnot( length(setdiff(colnames(zhang.feat), zhang.meta$sid))==0 )
zhang.feat<-zhang.feat[, zhang.meta$sid]
stopifnot( all.equal( colnames(zhang.feat), zhang.meta$sid) ) 
rm(normX)


#  save the metadata and the normalized counts
save(zhang.meta, zhang.feat, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/Zhang_Hertwig2015_SEQC_read_for_ML.RData')

#}}}


#  load the metadata and the normalized counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/Zhang_Hertwig2015_SEQC_read_for_ML.RData')


#  [our cohort] load the final list of important features, and the consistently up-/downregulated identified features and their performances
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_one_by_one_important.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/classification_ML_tumors_identified_per_risk-group.RData')


#  some of our identified features are not present in the count matrix
KEPT<-unique(names(unlist(imp.kept)))
length(intersect(KEPT, rownames(zhang.feat)))  #  466
KEPT<-intersect(KEPT, rownames(zhang.feat))


#  go over all our risk group-specific up-/downregulated features and keep only those found in the current data set
KEPT.UP<-lapply(imp.up, function(x){ x[ intersect(rownames(x), rownames(zhang.feat)), , drop=F] })
KEPT.DOWN<-lapply(imp.down, function(x){ x[ intersect(rownames(x), rownames(zhang.feat)), , drop=F] })


#  define risk factor 
RISK<-factor(setNames(zhang.meta$risk_group, zhang.meta$sid), levels=zhang.meta[, unique(risk_group)])


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  PCA and hierarchical clustering based on MOST VARIABLE AND HIGHER EXPRESSED features
#{{{

#  pick features based on CV>median(CV) and mean expression>median(expression)
cv<-apply(zhang.feat, 1, function(x){ sd(x)/mean(x) })
mn<-rowMeans(zhang.feat)


#  PCA 
p<-prcomp(t(zhang.feat[ cv>quantile(cv, 0.5) & mn>quantile(mn, 0.5), ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-factor(zhang.meta[ match( rownames(pca), sid ), risk_group ], levels=zhang.meta[, unique(risk_group)])
cl<-setNames(zhang.meta[, unique(col)], zhang.meta[, unique(risk_group)] )
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
    #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl) +
    scale_color_manual(name='risk group', values=cl) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/plotPCA_most-variable-expressed_features.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(p, ve, pca, cl, XLIM, YLIM)

#}}}


#  PCA and hierarchical clustering based on ALL KEPT features
#{{{

#  PCA 
p<-prcomp(t(zhang.feat[ KEPT, ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-factor(zhang.meta[ match( rownames(pca), sid ), risk_group ], levels=zhang.meta[, unique(risk_group)])
cl<-setNames(zhang.meta[, unique(col)], zhang.meta[, unique(risk_group)] )
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
    #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl) +
    scale_color_manual(name='risk group', values=cl) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/plotPCA_kept_features.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(p, ve, pca, cl, XLIM, YLIM)


#  hierarchical clustering of samples and features 
ex<-data.frame(risk_group=factor(zhang.meta[ match(colnames(zhang.feat), sid), risk_group], exclude=F), row.names=colnames(zhang.feat))
cl<-setNames(list(setNames(zhang.meta[, unique(col)], zhang.meta[, unique(risk_group)] )), 'risk_group')
x<-zhang.feat[ KEPT, ,drop=F]
d<-dist(t(x), method='euclidean')
hc.c<-hclust(d, method='ward.D2')
d<-dist(x, method='euclidean')
hc.r<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(x), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(2, 22, length.out=11),
        cluster_rows=hc.r,
        cluster_cols=hc.c,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
        annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/heatmap_hclust_kept_features.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc.r, hc.c, ph)


#  hierarchical clustering of samples based on Spearman correlations
ex<-data.frame(risk_group=factor(zhang.meta[ match(colnames(zhang.feat), sid), risk_group], exclude=F), row.names=colnames(zhang.feat))
cl<-setNames(list(setNames(zhang.meta[, unique(col)], zhang.meta[, unique(risk_group)] )), 'risk_group')
x<-zhang.feat[ KEPT, ,drop=F]
d<-cor(x, method='spearman', use='pairwise.complete.obs')
hc<-hclust(as.dist(1-d), method='ward.D2')
ph<-pheatmap(as.matrix(d), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(0.6, 1, length.out=11),
        cluster_rows=hc,
        cluster_cols=hc,
        annotation_legend=T, annotation_names_row=F, annotation_names_col=F, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/heatmap_hclust_cor_kept_features.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc, ph)

#}}}


#  PCA and hierarchical clustering based on ONLY HIGHER EXPRESSED KEPT features in a given risk group
#{{{

#  pool 
length(keep<-unique(unlist(lapply(KEPT.UP, rownames))))  #  151


#  PCA
p<-prcomp(t(zhang.feat[ keep, ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-factor(zhang.meta[ match( rownames(pca), sid ), risk_group ], levels=zhang.meta[, unique(risk_group)])
cl<-setNames(zhang.meta[, unique(col)], zhang.meta[, unique(risk_group)] )
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
    #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl) +
    scale_color_manual(name='risk group', values=cl) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/plotPCA_identified_features_up.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(p, ve, pca, cl, XLIM, YLIM)


#  hierarchical clustering of samples and features
ex<-data.frame(risk_group=factor(zhang.meta[ match(colnames(zhang.feat), sid), risk_group], exclude=F), row.names=colnames(zhang.feat))
cl<-setNames(list(setNames(zhang.meta[, unique(col)], zhang.meta[, unique(risk_group)] )), 'risk_group')
x<-zhang.feat[ keep, ,drop=F]
d<-dist(t(x), method='euclidean')
hc.c<-hclust(d, method='ward.D2')
d<-dist(x, method='euclidean')
hc.r<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(x), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(2, 22, length.out=11),
        cluster_rows=hc.r,
        cluster_cols=hc.c,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
        annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/heatmap_hclust_identified_features_up.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc.r, hc.c, ph)


#  expression distribution boxplots
ex<-data.frame(risk_group=factor(zhang.meta[ match(colnames(zhang.feat), sid), risk_group], exclude=F), row.names=colnames(zhang.feat))
cl<-setNames(list(setNames(zhang.meta[, unique(col)], zhang.meta[, unique(risk_group)] )), 'risk_group')
for(r in names(cl$risk_group)){
    pv<-expression_boxplots(X=zhang.feat[rownames(KEPT.UP[[r]]), , drop=F], Y=RISK, main='', trim.chr=T, 
                   legend.place='topright', yline=4, xline=0, xlas=2, xadj=1, xcex=1.6, ycex=2.4, 
                   ylab=expression(log[2]('Expression')), 
                   figure.name=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/boxplot_identified_features_', r, '_up.svg'), svg.width=min(65, max(16, nrow(KEPT.UP[[r]]))), colorize=NULL, cluster.colors=cl$risk_group, 
                   mar=c(10.0, 7.0, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
    readline('Press ENTER:')


    ph<-pheatmap(as.matrix(zhang.feat[ rownames(KEPT.UP[[r]]), , drop=F]), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), 
            border_color=NA, scale='none', 
            cluster_rows=F,
            cluster_cols=F,
            annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
            annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
            drop_levels=F, show_rownames=F, show_colnames=F, 
            fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
    readline('Press ENTER:')


    #  PCA
    p<-prcomp(t(zhang.feat[ rownames(KEPT.UP[[r]]), ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
    ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
    pca<-as.data.frame(p$x)
    pca$'risk group'<-factor(zhang.meta[ match( rownames(pca), sid ), risk_group ], levels=zhang.meta[, unique(risk_group)])
    XLIM<-pretty(range(p$x[, 'PC1']), 5)
    XLIM<-c(XLIM[1], XLIM[length(XLIM)])
    YLIM<-pretty(range(p$x[, 'PC2']), 5)
    YLIM<-c(YLIM[1], YLIM[length(YLIM)])
    g<-ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
        geom_point(size=10) + 
        scale_shape_manual(values=c(18:21)) + 
        scale_fill_manual(name='risk group', values=cl$risk_group) +
        scale_color_manual(name='risk group', values=cl$risk_group) + 
        theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
            axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
            legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
            aspect.ratio=1, panel.background = element_blank(),
            axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
            axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
        scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
        xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
    print(g)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/plotPCA_identified_features_', r, '_up.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    readline('Press ENTER:')
}

#}}}


#  PCA and hierarchical clustering based on ONLY LOWER EXPRESSED KEPT features in a given risk group
#{{{

#  pool 
length(keep<-unique(unlist(lapply(KEPT.DOWN, rownames))))  #  216


#  PCA
p<-prcomp(t(zhang.feat[ keep, ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
pca<-as.data.frame(p$x)
pca$'risk group'<-factor(zhang.meta[ match( rownames(pca), sid ), risk_group ], levels=zhang.meta[, unique(risk_group)])
cl<-setNames(zhang.meta[, unique(col)], zhang.meta[, unique(risk_group)] )
XLIM<-pretty(range(p$x[, 'PC1']), 5)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(p$x[, 'PC2']), 5)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
    geom_point(size=10) + 
    #geom_text_repel(aes(label=rownames(pca)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
    scale_shape_manual(values=c(18:21)) + 
    scale_fill_manual(name='risk group', values=cl) +
    scale_color_manual(name='risk group', values=cl) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/plotPCA_identified_features_down.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(p, ve, pca, cl, XLIM, YLIM)


#  hierarchical clustering of samples and features
ex<-data.frame(risk_group=factor(zhang.meta[ match(colnames(zhang.feat), sid), risk_group], exclude=F), row.names=colnames(zhang.feat))
cl<-setNames(list(setNames(zhang.meta[, unique(col)], zhang.meta[, unique(risk_group)] )), 'risk_group')
x<-zhang.feat[ keep, ,drop=F]
d<-dist(t(x), method='euclidean')
hc.c<-hclust(d, method='ward.D2')
d<-dist(x, method='euclidean')
hc.r<-hclust(d, method='ward.D2')
ph<-pheatmap(as.matrix(x), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), border_color=NA, scale='none', 
        breaks=seq(2, 22, length.out=11),
        cluster_rows=hc.r,
        cluster_cols=hc.c,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
        annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/heatmap_hclust_identified_features_down.svg', width=20, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(ex, cl, x, d, hc.r, hc.c, ph)


#  expression distribution boxplots
ex<-data.frame(risk_group=factor(zhang.meta[ match(colnames(zhang.feat), sid), risk_group], exclude=F), row.names=colnames(zhang.feat))
cl<-setNames(list(setNames(zhang.meta[, unique(col)], zhang.meta[, unique(risk_group)] )), 'risk_group')
for(r in names(cl$risk_group)){
    pv<-expression_boxplots(X=zhang.feat[rownames(KEPT.UP[[r]]), , drop=F], Y=RISK, main='', trim.chr=T, 
                   legend.place='topright', yline=4, xline=0, xlas=2, xadj=1, xcex=1.6, ycex=2.4, 
                   ylab=expression(log[2]('Expression')), 
                   figure.name=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/boxplot_identified_features_', r, '_down.svg'), svg.width=min(65, max(16, nrow(KEPT.UP[[r]]))), colorize=NULL, cluster.colors=cl$risk_group, 
                   mar=c(10.0, 7.0, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
    readline('Press ENTER:')


    ph<-pheatmap(as.matrix(zhang.feat[ rownames(KEPT.UP[[r]]), , drop=F]), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(10)), 
            border_color=NA, scale='none', 
            cluster_rows=F,
            cluster_cols=F,
            annotation_legend=T, annotation_names_row=T, annotation_names_col=F, 
            annotation_col=ex, annotation_row=NULL, annotation_colors=cl,
            drop_levels=F, show_rownames=F, show_colnames=F, 
            fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
    readline('Press ENTER:')


    #  PCA
    p<-prcomp(t(zhang.feat[ rownames(KEPT.UP[[r]]), ,drop=F]), center=T, scale.=F)  #  PCA on the samples is done when the samples are the rows
    ve<-round(100*p$sdev^2/sum(p$sdev^2), 1)
    pca<-as.data.frame(p$x)
    pca$'risk group'<-factor(zhang.meta[ match( rownames(pca), sid ), risk_group ], levels=zhang.meta[, unique(risk_group)])
    XLIM<-pretty(range(p$x[, 'PC1']), 5)
    XLIM<-c(XLIM[1], XLIM[length(XLIM)])
    YLIM<-pretty(range(p$x[, 'PC2']), 5)
    YLIM<-c(YLIM[1], YLIM[length(YLIM)])
    g<-ggplot(pca[, c('PC1', 'PC2', 'risk group')], aes(PC1, PC2, color=`risk group`)) + 
        geom_point(size=10) + 
        scale_shape_manual(values=c(18:21)) + 
        scale_fill_manual(name='risk group', values=cl$risk_group) +
        scale_color_manual(name='risk group', values=cl$risk_group) + 
        theme(text=element_text(family='Arial'), axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
            axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
            legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
            aspect.ratio=1, panel.background = element_blank(),
            axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
            axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
        scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
        xlab(paste0('PC1: ', ve[1], '% variance')) + ylab(paste0('PC2: ', ve[2], '% variance')) 
    print(g)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/figures/Zhang_Hertwig2015_SEQC/plotPCA_identified_features_', r, '_down.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    readline('Press ENTER:')
}

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/Zhang_Hertwig2015_SEQC_read_for_ML.RData




####################################################################################################
#
#
#  feature selection based on random forest classification of the Zhang and Hertwig 2015 SEQC cohort
#
#
####################################################################################################




#  self-contained machine learning algorithm that does the feature selection
#{{{
rm(list=ls())
library(rtracklayer)
library(data.table)
library(foreach)
library(doParallel)
library(iterators)
library(cluster)      #  silhouette()
library(randomForest)
source('~/bio/lib/enrichment_analysis.R')
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/feature_selection_functions.R')


#  global parameters
NTREE<-200               #  number of trees to grow in all random forest runs (we need good sampling of the cell lines)
LEAVE_OUT<-25            #  number of samples to leave out for cross-validation
MIN_IMPORTANCE<-3.0      #  minimum importance (MeanDecreaseAccuracy) to ask for selecting one-by-one the important features
MIN_OOB_ERROR<-0.60      #  minimum OOB class error to ask for selecting one-by-one the important features 
MAX_OOB_ERROR<-0.25      #  max OOB error to ask for at least one class during the second round of selection
FREQ_CUTOFF<-0.05        #  minimum frequency that an important feature survived trimmed random forest in order to be considered in the end 



#  load the metadata and the normalized counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/Zhang_Hertwig2015_SEQC_read_for_ML.RData')
feat<-zhang.feat
meta<-zhang.meta
rm(zhang.feat, zhang.meta)


#  define the risk groups
risk<-factor(setNames(meta$risk_group, meta$sid), levels=meta[, unique(risk_group)])


#  [step one] random forest on one-by-one all the features using all samples and the OOB error 
#{{{

#  register the cores
registerDoParallel(cores=detectCores())


#  run in parallel
Y<-risk
X<-feat[, names(Y)]
stopifnot(all.equal( names(Y), colnames(X) ) )
ITRAIN<-seq_len(ncol(X))
ONE_BY_ONE<-foreach_one_by_one_randomForest(t(X), Y, ITRAIN, NTREE)


#  unregister the cores
registerDoSEQ()


#  save the results
save(ONE_BY_ONE, ITRAIN, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/Zhang_Hertwig2015_SEQC_classification_ML_tumors_one_by_one.RData')

#}}}



#  [step two] identify important classifiers 
#             cluster them by Spearman correlation distance 
#             within each cluster do repeated trimmed random forest collecting errors on randomly drawn test sets
#             identify trimmed sets with up to maximum test classification errors and within them features of at least a minimum frequency of occurrence
#             do repeated trimmed random forest on the identified features
# 
#             MSigDB C2, GO BP/MF enrichment analysis for the final list of unique trimmed features against all features
#{{{

#  explore the one-by-one positive importances to get an idea if the cutoffs above make sense
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/Zhang_Hertwig2015_SEQC_classification_ML_tumors_one_by_one.RData')
imp<-ONE_BY_ONE[ imp>0 ]
hist(imp$imp, breaks=100)
summary(imp$imp)
#
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#  0.000000  6.770998 13.332815 16.662976 23.298565 85.976865 


#  MIN_IMPORTANCE of 3.0 removes approximately the bottom 10% 
imp<-imp[ imp>=MIN_IMPORTANCE ]
hist(imp$test_err, breaks=100)
dev.off()
summary(imp$test_err)
#
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0326087 0.4649123 0.5000000 0.4904960 0.5263158 0.6666667 


#  MIN_OOB_ERROR of 0.6 removes top 0.2%
imp<-imp[ test_err<=MIN_OOB_ERROR ][ order(-imp, +test_err )]


#  cluster based on Spearman correlations cutting at 2 major clusters
NCLUST<-2
Y<-risk
X<-feat[ rownames(feat) %in% imp[, feature], ]
CO<-cor(t(X), method='spearman')
CO.D<-as.dist(1-CO)
H<-hclust(CO.D, method='ward.D2')  
imp.w<-cutree(H, k=NCLUST)  
print(summary(silhouette(imp.w, dist=CO.D)))
W.check<-cutree(H, k=NCLUST+1)
print(summary(silhouette(W.check, dist=CO.D)))
#cor_heatmap(CO, imp.w, H, new.fig=F)  #  optionally visualize the clusters 
rm(W.check)


#  do trimmed random forest within each cluster
#  identify the features appearing with at least a threshold frequency
imp.kept<-list()
for(i in seq_len(NCLUST)){
    X<-feat[ rownames(feat) %in% names(imp.w)[ imp.w==i ], , drop=F]
    NRUNS<-ceiling(nrow(X)/2) - ceiling(nrow(X)/2) %% 50 + 50
    cat('\n\nstarting the', NRUNS, 'runs...\n')

    #  OOB errors are good enough for this purpose do not leave anything out at this stage
    er<-collect_errors_and_features(NRUNS, X=X, Y=Y, RETURN.FITS=F, FUN=select_features_randomForest, REPEATS=200, NTREE=NTREE, LEAVE.OUT=NULL)

    f<-sort(table(unlist(er$features[ apply(er$errors, 1, function(e){ any(e<=MAX_OOB_ERROR) }) ]))/NRUNS, decreasing=T)
    imp.kept[[i]]<-f[ f>=FREQ_CUTOFF ]
    cat('\n...done!\n')
}


#  use all features that survived the filtering for a last round of repeated trimmed random forest with cross-validation
X<-feat[ unique(names(unlist(imp.kept))), , drop=F]
imp.errors<-collect_errors_and_features(nrow(X), X=X, Y=Y, RETURN.FITS=T, FUN=select_features_randomForest, REPEATS=200, NTREE=NTREE, LEAVE.OUT=LEAVE_OUT) 


#  MSigDB:C2, GO:BP, GO:MF enrichment analysis for the genes (not circRNAs) associated with the important features
#{{{

#  load reference and remove subversions from gene_ids
hsa<-import('/data/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_name', 'gene_id', 'gene_type')]
hsa$gene_id<-sub('\\.[0-9]*$', '', hsa$gene_id)


#  GO/MSigDB C2 enrichment analyses with universe being the whole set of non-circRNA features 
U<-rownames(feat[ grep('^circ', rownames(feat), invert=T), , drop=F])
U<-setNames(hsa$gene_id[ match(U, hsa$gene_name) ], U)
G<-unique(unlist(imp.errors$features))
G<-G[ grep('^circ', G, invert=T) ]
G<-setNames(hsa$gene_id[ match(G, hsa$gene_name) ], G)
imp.enrich<-enrichment_analysis(G, U, MS=setNames('/data/annotation/MSigDB/c2.all.v7.0.symbols.gmt', 'C2'))
rm(hsa, U, G)

#}}}


#  save
save(imp, imp.w, imp.kept, imp.errors, imp.enrich, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/Zhang_Hertwig2015_SEQC_classification_ML_tumors_one_by_one_important.RData')

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/Zhang_Hertwig2015_SEQC_classification_ML_tumors_one_by_one.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/feature_selection/Zhang_Hertwig2015_SEQC_classification_ML_tumors_one_by_one_important.RData




##########
#
#
#  scratch
#
#
##########



#  fucking boxplot shit
#{{{

X<-feat[f, , drop=F]
Y<-risk
o<-names(sort(rowMeans(X), decreasing=T))
B<-unlist(lapply(split(names(Y), Y), function(x){ x<-X[o, x, drop=F]; unlist(apply(x, 1, list), recursive=F) }), recursive=F)
N<-nrow(X)
B<-B[c(sapply(seq_len(N), function(n){ seq(n, length(B), N) }))]  #  reorder so that clusters per feature are besides each other
o.cl<-setNames(rep('black', length(o)), o)
B.cl<-cl$risk_group
YTICK<-pretty(range(unlist(B)), 4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=rep(B.cl, N), ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=rep(B.cl, N), range=0, add=T)
ylab<-expression(log[10]('Expression'))
mtext(ylab, side=2, line=yline, padj=+0.1, las=0, cex=ycex)
mtext(text=o, side=1, line=xline, at=seq(2, length(B), length(levels(Y))), las=xlas, adj=1, padj=NA, cex=xcex, col=o.cl)

#}}}






