###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




##############################
#
#
#  circRNAs RBP motif analysis
#
#
##############################




#  [all RBPs] neuroblastoma-DE RBPs mRNA expression compared across tissues (human brain, various tumors)
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/load_gene_expression.R')
source('~/bio/lib/draw_highlights.R')


#  load all RBPs
#  inflate per motif the subgroups based on PWMs into one group of unique genes (motif-inflating needs to happen alone, different 1:many map than genes)
#  keep 6mers or longer
#  collect per RBP gene all its motifs
#  define GENES data.frame
load('/fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')
RBP<-db[, .(gene_name=unique(unlist(gene_name)), gene_id=unique(unlist(gene_id)), motif=motif), by=.(id)][, .(motif=unlist(motif)), by=.(gene_name, gene_id)][ sapply(motif, nchar)>=6, ][, .(motif=list(unique(motif))), by=.(gene_name, gene_id)]
GENES<-data.frame(RBP[, c('gene_name', 'gene_id')])
rm(db)


#  get CPMs across tissues/tumors
#  Mann-Whitney two-sided tests for neuroblastoma vs human brain and neuroblastoma vs various tumors for those genes found expressed in 
#  at least 50% of samples with at least 30% of them having a CPM>=1
#  FDR-correct p-values
#  define one p-value to be the sum of the FDR-corrected p-values if both are significant (<0.05)
#  prepare for plotting
#  save
#{{{

#  get their CPMs across tissues/tumors
x<-load_gene_expression(GENES, nb.tumors.only=T, vt.tumors.only=T)
nb.genes<-x$nb$genes
nb.meta<-x$nb$meta
hb.genes<-x$hb$genes
hb.meta<-x$hb$meta
vt.genes<-x$vt$genes
vt.meta<-x$vt$meta
rm(x,load_gene_expression)


#  [neuroblastoma] keep only tumors
#                  keep only CPMs and transpose
nb.genes<-t(nb.genes[ !is.na(nb.genes$risk_group), -tail(seq_len(ncol(nb.genes)),3)])


#  [human brain] keep only CPMs and transpose 
hb.genes<-t(hb.genes[, -tail(seq_len(ncol(hb.genes)),1) ])


#  [various tumors] keep only CPMs and transpose 
vt.genes<-t(vt.genes[, -tail(seq_len(ncol(vt.genes)),1) ])


#  go over each gene and do Mann-Whitney tests for the CPMs in NB vs brain tissue and NB vs various tumors
GENES$pv.hb<-GENES$pv.vt<-NA
for (n in seq_along(GENES$gene_name)){

    x<-nb.genes[ n, ]
    y<-hb.genes[ n, ]
    z<-vt.genes[ n, ]

    #  at least 50% of samples in the whole cohort should express it and at least 30% of those should have CPM>=1
    if(quantile(x[x!=0], 0.7)<1 | sum(x!=0)/length(x)<0.5){ next }  

    #  Mann-Whitney test that NB are higher expressed
    GENES$pv.hb[n]<-wilcox.test(x=x, y=y, alternative='two.sided')$p.value
    GENES$pv.vt[n]<-wilcox.test(x=x, y=z, alternative='two.sided')$p.value
}
rm(n,x,y,z)


#  FDR-adjust p-values
GENES$pv.hb<-p.adjust(GENES$pv.hb, method='BH')
GENES$pv.vt<-p.adjust(GENES$pv.vt, method='BH')


#  define a p-value which is the sum of p-values provided both FDR-corrected p-values are significant (<0.05)
GENES$pvalue<-GENES$pv.hb + GENES$pv.vt 
GENES$pvalue[ ! ( GENES$pv.hb<0.05 & GENES$pv.vt<0.05) ]<-NA


#  order by p-value
GENES<-GENES[ order(GENES$pvalue, decreasing=F), ]
rownames(GENES)<-NULL
nb.genes<-nb.genes[ GENES$gene_name, ]
hb.genes<-hb.genes[ GENES$gene_name, ]
vt.genes<-vt.genes[ GENES$gene_name, ]
stopifnot( all.equal( rownames(nb.genes), GENES$gene_name ) )
stopifnot( all.equal( rownames(hb.genes), GENES$gene_name ) )
stopifnot( all.equal( rownames(vt.genes), GENES$gene_name ) )


#  save for easy use
save(RBP, GENES, nb.genes, nb.meta, hb.genes, hb.meta, vt.genes, vt.meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues.RData')

#}}}


#  load above results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues.RData')


#  keep only significantly DE throughout
GENES<-subset(GENES, pvalue<0.05)
nb.genes<-nb.genes[ GENES$gene_name, ]
hb.genes<-hb.genes[ GENES$gene_name, ]
vt.genes<-vt.genes[ GENES$gene_name, ]


#  order by NB expression
n<-order(rowMeans(nb.genes), decreasing=T)
GENES<-GENES[n, ]
nb.genes<-nb.genes[n, ]
hb.genes<-hb.genes[n, ]
vt.genes<-vt.genes[n, ]
stopifnot( all.equal( rownames(nb.genes), GENES$gene_name ) )
stopifnot( all.equal( rownames(nb.genes), rownames(hb.genes)) )
stopifnot( all.equal( rownames(nb.genes), rownames(vt.genes)) )
rm(n)


#  identify those classified as groupI in Gerstberger et al. 2014 which are down-regulated upon brain development
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/validations/rbp/brainspan_census_of_human_RBPs_nrg3813/nrg3813-s7.RData')
GENES$groupI<-F
GENES$groupI[ GENES$gene_name %in% rbp[ groupI==T, gene_name] ]<-T


#  recycle
x11(width=25, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  heatmap on the log2-transformed (best contrast, compared to z-scores of raw CPMs, or z-scores of log2-transformed CPMs)
#  clustered by Euclidean distance in NB samples
X<-cbind(log2(1+nb.genes), log2(1+hb.genes), log2(1+vt.genes))  #  combine in one matrix
ex<-data.frame(tissue=factor(rep(c('neuroblastoma', 'brain tissue', 'various tumors'), c(ncol(nb.genes), ncol(hb.genes), ncol(vt.genes))), exclude=F), row.names=colnames(X))
cl<-setNames(list(setNames( c('seagreen4', 'cornflowerblue', 'coral4'), c('neuroblastoma', 'brain tissue', 'various tumors') )), colnames(ex))
hc<-hclust(dist(log2(1+nb.genes), method='euclidean'), method='ward.D2')  #  cluster based on distance between rows for NB only samples 
#hc<-hclust(dist(X, method='euclidean'), method='ward.D2')
ph<-pheatmap(X, color=viridis(10), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=NA, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=F, 
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_RBPs_across_tissues.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(X,hc,ex,cl,ph)


#  boxplots of log2-transformed CPMs
B<-c(split(log2(1+nb.genes), row(nb.genes)), split(log2(1+hb.genes), row(hb.genes)), split(log2(1+vt.genes), row(vt.genes)))
n<-c(rbind(rbind(1:nrow(nb.genes), nrow(nb.genes)+(1:nrow(nb.genes)), 2*nrow(nb.genes)+(1:nrow(nb.genes)))))  #  order them in triplets
B<-B[n]
B.cl<-rep(c('seagreen4', 'cornflowerblue', 'coral4'), nrow(nb.genes))
groupI.cl<-setNames(rep('black', nrow(nb.genes)), rownames(nb.genes))
#groupI.cl[ GENES$gene_name[ GENES$groupI ] ]<-'darkorange'
par(mar=c(6.0, 7.0, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(-2, 2), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2](1+'CPM')), side=2, line=3, padj=-0.2, las=0, cex=2.4)
mtext(text=rownames(nb.genes), side=1, line=0, at=seq(2, length(B), 3), las=2, adj=1, cex=1.2, col=groupI.cl)
draw_highlights(L=length(B), STEP=3, YMAX=max(YTICK))
legend('topleft', legend=c('neuroblastoma', 'brain tissue', 'various tumors'), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_RBPs_across_tissues.svg', width=35, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(YTICK,B,n,B.cl,bp)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues.RData



#  [MNA vs HR_nMNA] intersection of DE splicing factors with neuroblastoma-DE RBPs
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load the neuroblastoma-DE RBPs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues.RData')


#  load DE results for genes and split to up-/downregulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_MNA_HR_nMNA.RData')
up.gns<-subset(RES, log2FoldChange>0 & padj<0.05)
down.gns<-subset(RES, log2FoldChange<0 & padj<0.05)
rm(RES, CND, DDS, PCA, VE, VSC)


#  load union of splice-factors and identify the significantly up-/downregulated
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/SpliceAid-F+GO.RData')
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/MSigDB_c2_v7.0_splicing.RData')
db<-unique(rbind(data.frame(db[, c('gene_name', 'gene_id')]), msigdbSF[, c('gene_name', 'gene_id')]))
up.sf<-up.gns[ rownames(up.gns) %in% db$gene_id, , drop=F]
down.sf<-down.gns[ rownames(down.gns) %in% db$gene_id, , drop=F]
rm(db,msigdbSF)


#  upregulated splicing factors (basically all of them)
RBP[ gene_id %in% rownames(up.sf) ]
#     gene_name            gene_id                                                              motif
#  1:     SRSF3 ENSG00000112081.17          AACTTTAT,ACATTCAT,ATCATCAT,ATCTTCAC,ATCTTCAT,CACATCAT,...
#  2:      YBX1 ENSG00000065978.19                AGCGAGC,CGAGCGG,GAGCGAG,GAGCGGA,GCGAGCG,CCCTGCG,...
#  3:    DAZAP1 ENSG00000071626.16                    AAAAAAA,AATTTA,AGATAT,AGTAGG,GGGGGGG,GTAACG,...
#  4:       FUS ENSG00000089280.18                   AAAAAAA,GGGGGGG,TTTTTTT,CGGTGA,CGGTGG,GGGTGA,...
#  5:   HNRNPA1 ENSG00000135486.17                   AAAAAAA,AATTTA,AGATAT,AGTAGG,TTTTTTT,CCCCCCC,...
#  6:    HNRNPD ENSG00000138668.19                  AAAAAAA,AATTTA,AGATAT,AGTAGG,ATTTATTTA,TTAGAG,...
#  7:    HNRNPK ENSG00000165119.21                  AAAAAAA,GGGGGGG,TTTTTTT,CCCCCCC,ACCCAA,ACCCAT,...
#  8:    HNRNPU ENSG00000153187.20                                    AAAAAAA,GGGGGGG,TTTTTTT,CCCCCCC
#  9:     PCBP1  ENSG00000169564.6                  AAAAAAA,GGGGGGG,TTTTTTT,CCCCCCC,TTAGAG,TTAGGA,...
# 10:   SYNCRIP ENSG00000135316.17                                                    AAAAAAA,TTTTTTT
# 11:   HNRNPA0  ENSG00000177733.6                                               AATTTA,AGATAT,AGTAGG
# 12:    HNRNPC ENSG00000092199.17                 GGGGGGG,TTTTTTT,GGATAC,ATTTTTG,TTTTTTG,CTTTTTG,...
# 13:    HNRNPF ENSG00000169813.16             GGGGGGG,AGGGAT,AAGGTG,GGAGGA,AGGGAAGGGA,AGGGGAGGGG,...
# 14:   HNRNPH1 ENSG00000169045.17         GGGGGGG,AAGGTG,GGAGGA,AGGGAAGGGA,AGGGGAGGGG,CGGGGGGGGC,...
# 15:   HNRNPH3 ENSG00000096746.17         GGGGGGG,AAGGTG,GGAGGA,AGGGAAGGGA,AGGGGAGGGG,CGGGGGGGGC,...
# 16:    HNRNPM ENSG00000099783.12                                                    GGGGGGG,GAAGGAA
# 17:     PCBP2 ENSG00000197111.15                   GGGGGGG,TTTTTTT,CCCCCCC,TTAGAG,TTAGGA,TTAGGG,...
# 18:     PTBP1 ENSG00000011304.20                   GGGGGGG,TTTTTTT,CCCCCCC,ATCTTC,CTCTTA,CTCTTC,...
# 19:    ELAVL1 ENSG00000066044.15         TTTTTTT,ATTTATTTATTT,CCCCCCC,TTATTTATT,TTATTTT,TTGTTTT,...
# 20:     U2AF2 ENSG00000063244.12                                            TTTTTTT,TTTTTCC,TTTTTTC
# 21:     TRA2B ENSG00000136527.18               GAAGGA,GAAGAA,AAGAAG,AAGAAGAA,AAGAAGAAGAA,AAGAAC,...
# 22:     ESRP2 ENSG00000103067.13                    TGGGAAA,TGGGGAA,TGGGAAT,TGGGGAT,TGGGAAG,TGGGGAG
# 23:     SRSF1 ENSG00000136450.13                      AAGGTG,GGATAC,GCATAC,GGATAT,GGATTC,GGGTAC,...
# 24:     SRSF2 ENSG00000161547.16    CTAGACTAGA,GAGGAG,GGAGGA,GTAAGTACGC,AAAAGAGAAG,AGAGGAAGGCGA,...
# 25:   HNRNPA3 ENSG00000170144.20                                                        GCCAAGGAGCC
# 26: HNRNPA1L2 ENSG00000139675.12                    ATAGGGA,TTAGGGA,GTAGGGA,ATAGGGT,TTAGGGT,GTAGGGT
# 27:    HNRNPL ENSG00000104824.17                ACACAAA,ACACGAA,ACACAAC,ACACGAC,ACACAAG,ACACGAG,...
# 28:     RBM25 ENSG00000119707.14                                                    ATCGGGCA,CGGGCA
# 29:      RBMX ENSG00000147274.14                      ACCAAA,ATCAAA,ATCCCA,ATCCCC,TAAGAC,TCAAAA,...
# 30:     SNRPA  ENSG00000077312.9   AGGAGAT,ATTGCAC,ATTGCACC,GAGCAGTAGGC,GAGCAGTAGTC,GAGCAGTAGGG,...
# 31:     SRSF7 ENSG00000115875.19 AGAGGAAGGCGA,GAAGAAGAA,CTCTTCAC,AAAGGACAAA,ACGAATGAT,ACGAGACTA,...
# 32:     SRSF6 ENSG00000124193.15                 GAAGAAGA,GAGGAAGAA,ACCGGG,ACCGTC,AGCGGA,ATCGTA,...
#     gene_name            gene_id                                                              motif


#  downregulated SF (basically all of them)
RBP[ gene_id %in% rownames(down.sf) ]
#    gene_name            gene_id                                                                 motif
# 1:    TARDBP ENSG00000120948.17          TGTGTG,TGTGTGTG,TGTGTGTGTG,GTGAATGA,GTTGTGC,TGTGTGTGTGTG,...
# 2:     CELF6 ENSG00000140488.16                       TGTGAGG,TGTGTGG,TGTGGGG,TGTGATG,TGTGTTG,TGTGGTG
# 3:    ELAVL4 ENSG00000162374.17 AAAAAAA,TTTTTTT,ATTTATTTATTT,TTTATTTATTT,TTTATTTATTTA,TTTTTTATTTT,...
# 4:      PPIE ENSG00000084072.16                         TTTTTT,AAAAAA,AATAAA,TAAAAA,ATAAAA,TTAAAA,...
# 5:    RBFOX1 ENSG00000078328.20                     TGCATG,AGCATG,TGCATGT,AGCATGA,TGCATGA,AGCATGC,...
# 6:     SRSF4 ENSG00000116350.17                                                    GAAGAAGA,GAGGAAGAA
#    gene_name            gene_id                                                                 motif

#}}}



#  [circRNAs] DREME : 6mer, 7mer, 8mer motif enrichment
#               AME : ATtRACT human RBP motif enrichment for 6mers and above
#              FIMO : count individual motif occurrences for 6mers and above 
#{{{
rm(list=ls())
library(data.table)
library(GenomicAlignments)
library(rtracklayer)
library(XML)


#  create the subdirectory and run in separate screen sessions the different analyses
#{{{

system2('mkdir', args='-p /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme')
#
#
#  run MEME (de novo 7mer-50mer motif enrichment):
#  N.B. when full sequence lengths are used this becomes unbelievably slow, anyway there is no need to run it, AME is what we need
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results_meme/
#      fasta-get-markov -m 0 -rna -norc ../../circRNAs_sequences.fa background.model
#      meme -rna -mod anr -oc results_meme -minw 7 -maxw 20 -markov_order 0 -bfile background.model -nmotifs 50 -allw -neg ../../circRNAs_controls_sequences.fa -objfun de -test mhg -searchsize 0 ../../circRNAs_sequences.fa 
#
#  run DREME (6mer, 7mer, 8mer motif enrichment):
#  N.B. this is also unbelievably slow and anyway I never used these results, AME is what we need
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results_dreme/
#      dreme-py3 -p ../../circRNAs_sequences.fa -n ../../circRNAs_controls_sequences.fa -rna -norc -e 0.05 -mink 6 -maxk 8 -oc results_dreme -eps -png
#
#  run AME (known motif enrichment using the ATtRACT Human motifs converted to the MEME motif format):
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results_ame
#      fasta-get-markov -m 0 -rna -norc ../../circRNAs_sequences.fa background.model
#      ame --oc ./ --method fisher --scoring avg --evalue-report-threshold 0.05 --bfile background.model --control ../../circRNAs_controls_sequences.fa ../../circRNAs_sequences.fa /fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme
#
#  run FIMO (counting motif occurrences):
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results_fimo
#      fasta-get-markov -m 0 -rna -norc ../../circRNAs_controls_sequences.fa control_background.model
#      fasta-get-markov -m 0 -rna -norc ../../circRNAs_sequences.fa signal_background.model
#
#      fimo --oc ./ --norc --no-qvalue --thresh 1e-4 --bfile signal_background.model -max-stored-scores 1000000 /data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme ../../circRNAs_sequences.fa 
#      fimo --oc ./control --norc --no-qvalue --thresh 1e-4 --bfile control_background.model -max-stored-scores 1000000 /data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme ../../circRNAs_controls_sequences.fa

#}}}


#  load reference to add gene_id to the results
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]


#  AME results
#  convert p-values to numerics
#  add FDR = FP/(FP+TP) column 
ame<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results_ame/ame.tsv', header=T, sep='\t', strip.white=T, blank.lines.skip=F, drop=c('motif_DB', 'FASTA_max', 'PWM_min'))[, fdr:=round(FP/(FP+TP), 3)]  #  comments at the bottom of the file produce warning
colnames(ame)<-c('rank', 'motif_id', 'motif_alt_id', 'consensus', 'pvalue', 'padj', 'evalue', 'tests', 'pos', 'neg', 'tp', '%tp', 'fp', '%fp', 'fdr')
ame<-ame[, c('pvalue', 'padj', 'evalue'):=list(as.numeric(pvalue), as.numeric(padj), as.numeric(evalue))]


#  [do not run] DREME results
#{{{

#  import the XML version of the results keeping the consensus and all equivalent motifs
#  bind them to data.table converting to integer/numeric when appropriate and U->T
#
#      N : ACGT
#      V : ACG
#      H : ACT
#      D : AGT
#      B : CGT
#      M : AC
#      R : AG
#      W : AT
#      S : CG
#      Y : CT
#      K : GT
#
dreme<-do.call(rbind, lapply(xmlToList(xmlParse('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results_dreme/dreme.xml'))[[2]], function(x){
    a<-data.table(data.frame(t(x$.attrs), all_seq=unname(sapply(x[ grep('match', names(x)) ], '[', 1)), row.names=NULL))[, .(length=as.integer(length[1]), 
       nsites=as.integer(nsites[1]), 
       p=as.integer(p[1]), 
       n=as.integer(n[1]), 
       pvalue=as.numeric(pvalue[1]), 
       evalue=as.numeric(evalue[1]), 
       unerased_evalue=as.numeric(unerased_evalue[1]), 
       all_seq=list(all_seq)), by=.(id, alt, seq)]; 
    return(a)}))
dreme[, c('seq', 'all_seq'):=list(gsub('U', 'T', seq), lapply(all_seq, function(s){ gsub('U', 'T', s) }))]

#}}}


#  FIMO results
#{{{

#  FIMO cuts off motif_ids at 100 characters. We need to fix that.
load('/fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')
rm(pwm)

 
#  load FIMO results (warning about last line skipping is normal, it contains comments about FIMO command run)
fimo<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results_fimo/fimo.tsv', header=T, strip.white=T, blank.lines.skip=T, select=c('motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'score', 'p-value', 'q-value', 'matched_sequence'))
colnames(fimo)<-c('motif_id', 'motif_alt_id', 'circ_name', 'start', 'end', 'score', 'pvalue', 'qvalue', 'seq')
fimo.control<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results_fimo/control/fimo.tsv', header=T, strip.white=T, blank.lines.skip=T, select=c('motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'score', 'p-value', 'q-value', 'matched_sequence'))
colnames(fimo.control)<-c('motif_id', 'motif_alt_id', 'circ_name', 'start', 'end', 'score', 'pvalue', 'qvalue', 'seq')


#  fix cut motif_ids
#{{{

#  identify the motif_ids 100 characters long, most likely they were cut, and the corresponding alternative motif_ids
#  construct the full-length motif_ids 
#  check that each cut motif_id matches one full length motif_id
m<-fimo[nchar(motif_id)==100, unique(motif_id)]
a<-fimo[ motif_id %in% m, unique(motif_alt_id)]
d<-db[ motif %in% a, paste(sapply(gene_name, paste, sep='', collapse=','), id, sep='_') ]
for(i in m){
    k<-grep(i, d, fixed=T)
    if (length(k)==0){
        stop(paste('cut motif:', i, 'does not have a home!'))
    } else if (length(k)>1){
        stop(paste('cut motif:', i, 'has many homes!'))
    }
    fimo[ motif_id %in% i, motif_id:=d[k] ]
}


#  do the same for the controls
m<-fimo.control[nchar(motif_id)==100, unique(motif_id)]
a<-fimo.control[ motif_id %in% m, unique(motif_alt_id)]
d<-db[ motif %in% a, paste(sapply(gene_name, paste, sep='', collapse=','), id, sep='_') ]
for(i in m){
    k<-grep(i, d, fixed=T)
    if (length(k)==0){
        stop(paste('cut motif:', i, 'does not have a home!'))
    } else if (length(k)>1){
        stop(paste('cut motif:', i, 'has many homes!'))
    }
    fimo.control[ motif_id %in% i, motif_id:=d[k] ]
}

#}}}


#  identify the RBP groups for each motif
fimo<-fimo[, rbp_group:=lapply(strsplit(sub('_.*$', '', motif_id), ','), unique)]
fimo.control<-fimo.control[, rbp_group:=lapply(strsplit(sub('_.*$', '', motif_id), ','), unique)]


#  identify the overlapping motifs per sequence
fimo[, o_group:={g<-seq_len(.N);
                 o<-mcols(reduce(IRanges(start=start, end=end), with.revmap=T))$revmap  #  God bless for reverse-map!!!!!
                 invisible(sapply(seq_along(o), function(n){ g[ o[[n]] ]<<-n }))        #  assign same group number to all overlapping motifs
                 g
                }, by=.(circ_name)]
fimo.control[, o_group:={g<-seq_len(.N);
                         o<-mcols(reduce(IRanges(start=start, end=end), with.revmap=T))$revmap  #  God bless for reverse-map!!!!!
                         invisible(sapply(seq_along(o), function(n){ g[ o[[n]] ]<<-n }))        #  assign same group number to all overlapping motifs
                         g
                        }, by=.(circ_name)]


#  group overlapping motifs per sequence
fimo<-fimo[, .(motif_id=list(unique(motif_id)), motif_alt_id=list(unique(motif_alt_id)), start=list(start), end=list(end), score=list(score), pvalue=list(pvalue), qvalue=list(qvalue), seq=list(seq), rbp_group=list(unique(unlist(rbp_group)))), by=.(circ_name, o_group)][, o_group:=NULL]
fimo.control<-fimo.control[, .(motif_id=list(unique(motif_id)), motif_alt_id=list(unique(motif_alt_id)), start=list(start), end=list(end), score=list(score), pvalue=list(pvalue), qvalue=list(qvalue), seq=list(seq), rbp_group=list(unique(unlist(rbp_group)))), by=.(circ_name, o_group)][, o_group:=NULL]


#  keep on the side all motifs grouped by overlaps
fimo.all<-copy(fimo)
fimo.control.all<-copy(fimo.control)


#  identify number of non-overlapping motifs per sequence collapsing motif_ids, motif_alt_ids, motif sequences and the RBPs involved
fimo<-fimo[, .(ncount=.N, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(circ_name)]
fimo.control<-fimo.control[, .(ncount=.N, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(circ_name)]


#  add gene_ids to the RBPs involved per sequence
fimo[, gene_ids:=lapply(rbp_group, function(r){ hsa$gene_id[ match(r, hsa$gene_name) ] })]
fimo.control[, gene_ids:=lapply(rbp_group, function(r){ hsa$gene_id[ match(r, hsa$gene_name) ] })]


#  include lengths of the putative circRNA sequences 
#  do the same for the controls 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')
l<-setNames(width(CIRCS.exons.seqs), names(CIRCS.exons.seqs))
fimo[, length:=l[ circ_name ]]
fimo.all[, length:=l[ circ_name ]]
l<-setNames(width(CIRCS.controls.seqs.resized), names(CIRCS.controls.seqs.resized))  #  N.B. THE RESIZED SEQUENCES ARE USED AS CONTROLS!!!
fimo.control[, length:=l[ circ_name ]]
fimo.control.all[, length:=l[ circ_name ]]


#  compute count densities per 1K of sequence
fimo[, density:=ncount/length*1e3]
fimo.control[, density:=ncount/length*1e3]

#}}}


#  save
save(ame, fimo, fimo.all, fimo.control, fimo.control.all, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData



#  [circRNA introns] DREME : 6mer, 7mer, 8mer motif enrichment 
#                      AME : ATtRACT human RBP motif enrichment for 7mers and above
#                     FIMO : count individual motif occurrences for 7mers and above
#{{{
rm(list=ls())
library(data.table)
library(rtracklayer)
library(XML)


#  create the subdirectory and run in separate screen sessions the different analyses
#{{{

system2('mkdir', args='-p /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns')
#
#
#  combine upstream and downsteam intron sequences in one FASTA file BUT MAKE SURE TO ANNOTE WHICH IS WHICH:
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns
#      cat <(sed '/^>/{s/$/_up/}' ../../circRNAs_introns_up.fa) <(sed '/^>/{s/$/_down/}' ../../circRNAs_introns_down.fa) > introns.fa 
#
#
#  run DREME (6mer, 7mer, 8mer motif enrichment):
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results_dreme/
#      dreme-py3 -p ../introns.fa -n ../../../circRNAs_introns_controls.fa -rna -norc -e 0.05 -mink 6 -maxk 8 -oc results_dreme -eps
#
#  run AME (known motif enrichment using the ATtRACT Human motifs converted to the MEME motif format):
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results_ame/
#      fasta-get-markov -m 0 -rna -norc ../introns.fa background.model
#      ame --oc ./ --method fisher --scoring avg --evalue-report-threshold 0.05 --bfile background.model --control ../../../circRNAs_introns_controls.fa ../introns.fa /fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme
#
#  run FIMO (counting motif occurrences):
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results_fimo
#      fasta-get-markov -m 0 -rna -norc ../../../circRNAs_introns_controls.fa control_background.model
#      fasta-get-markov -m 0 -rna -norc ../introns.fa signal_background.model
#
#      fimo --oc ./ --norc --no-qvalue --thresh 1e-4 --bfile signal_background.model -max-stored-scores 1000000 /data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme ../introns.fa 
#      fimo --oc ./control --norc --no-qvalue --thresh 1e-4 --bfile control_background.model -max-stored-scores 1000000 /data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme ../../../circRNAs_introns_controls.fa

#}}}


#  load reference to add gene_id to the results
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]


#  AME results
#  convert p-values to numerics
#  add FDR = FP/(FP+TP) column 
ame<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results_ame/ame.tsv', header=T, sep='\t', strip.white=T, blank.lines.skip=F, drop=c('motif_DB', 'FASTA_max', 'PWM_min'))[, fdr:=round(FP/(FP+TP), 3)]  #  comments at the bottom of the file produce warning
colnames(ame)<-c('rank', 'motif_id', 'motif_alt_id', 'consensus', 'pvalue', 'padj', 'evalue', 'tests', 'pos', 'neg', 'tp', '%tp', 'fp', '%fp', 'fdr')
ame<-ame[, c('pvalue', 'padj', 'evalue'):=list(as.numeric(pvalue), as.numeric(padj), as.numeric(evalue))]


#  [do not run] DREME results
#{{{

#  import the XML version of the results keeping the consensus and all equivalent motifs
#  bind them to data.table converting to integer/numeric when appropriate and U->T
#
#      N : ACGT
#      V : ACG
#      H : ACT
#      D : AGT
#      B : CGT
#      M : AC
#      R : AG
#      W : AT
#      S : CG
#      Y : CT
#      K : GT
#
dreme<-do.call(rbind, lapply(xmlToList(xmlParse('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results_dreme/dreme.xml'))[[2]], function(x){
    a<-data.table(data.frame(t(x$.attrs), all_seq=unname(sapply(x[ grep('match', names(x)) ], '[', 1)), row.names=NULL))[, .(length=as.integer(length[1]), 
       nsites=as.integer(nsites[1]), 
       p=as.integer(p[1]), 
       n=as.integer(n[1]), 
       pvalue=as.numeric(pvalue[1]), 
       evalue=as.numeric(evalue[1]), 
       unerased_evalue=as.numeric(unerased_evalue[1]), 
       all_seq=list(all_seq)), by=.(id, alt, seq)]; 
    return(a)}))
dreme[, c('seq', 'all_seq'):=list(gsub('U', 'T', seq), lapply(all_seq, function(s){ gsub('U', 'T', s) }))]

#}}}


#  FIMO results
#{{{

#  FIMO cuts off motif_ids at 100 characters. We need to fix that.
load('/fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')
rm(pwm)

 
#  load FIMO results (warning about last line skipping is normal, it contains comments about FIMO command run)
fimo<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results_fimo/fimo.tsv', header=T, strip.white=T, blank.lines.skip=T, select=c('motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'score', 'p-value', 'q-value', 'matched_sequence'))
colnames(fimo)<-c('motif_id', 'motif_alt_id', 'intron_name', 'start', 'end', 'score', 'pvalue', 'qvalue', 'seq')
fimo.control<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results_fimo/control/fimo.tsv', header=T, strip.white=T, blank.lines.skip=T, select=c('motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'score', 'p-value', 'q-value', 'matched_sequence'))
colnames(fimo.control)<-c('motif_id', 'motif_alt_id', 'intron_name', 'start', 'end', 'score', 'pvalue', 'qvalue', 'seq')


#  fix cut motif_ids
#{{{

#  identify the motif_ids 100 characters long, most likely they were cut, and the corresponding alternative motif_ids
#  construct the full-length motif_ids 
#  check that each cut motif_id matches one full length motif_id
m<-fimo[nchar(motif_id)==100, unique(motif_id)]
a<-fimo[ motif_id %in% m, unique(motif_alt_id)]
d<-db[ motif %in% a, paste(sapply(gene_name, paste, sep='', collapse=','), id, sep='_') ]
for(i in m){
    k<-grep(i, d, fixed=T)
    if (length(k)==0){
        stop(paste('cut motif:', i, 'does not have a home!'))
    } else if (length(k)>1){
        stop(paste('cut motif:', i, 'has many homes!'))
    }
    fimo[ motif_id %in% i, motif_id:=d[k] ]
}


#  do the same for the controls
m<-fimo.control[nchar(motif_id)==100, unique(motif_id)]
a<-fimo.control[ motif_id %in% m, unique(motif_alt_id)]
d<-db[ motif %in% a, paste(sapply(gene_name, paste, sep='', collapse=','), id, sep='_') ]
for(i in m){
    k<-grep(i, d, fixed=T)
    if (length(k)==0){
        stop(paste('cut motif:', i, 'does not have a home!'))
    } else if (length(k)>1){
        stop(paste('cut motif:', i, 'has many homes!'))
    }
    fimo.control[ motif_id %in% i, motif_id:=d[k] ]
}

#}}}


#  identify the RBP groups for each motif
fimo<-fimo[, rbp_group:=lapply(strsplit(sub('_.*$', '', motif_id), ','), unique)]
fimo.control<-fimo.control[, rbp_group:=lapply(strsplit(sub('_.*$', '', motif_id), ','), unique)]


#  identify the overlapping motifs per sequence
fimo[, o_group:={g<-seq_len(.N);
                 o<-mcols(reduce(IRanges(start=start, end=end), with.revmap=T))$revmap  #  God bless for reverse-map!!!!!
                 invisible(sapply(seq_along(o), function(n){ g[ o[[n]] ]<<-n }))        #  assign same group number to all overlapping motifs
                 g
                }, by=.(intron_name)]
fimo.control[, o_group:={g<-seq_len(.N);
                         o<-mcols(reduce(IRanges(start=start, end=end), with.revmap=T))$revmap  #  God bless for reverse-map!!!!!
                         invisible(sapply(seq_along(o), function(n){ g[ o[[n]] ]<<-n }))        #  assign same group number to all overlapping motifs
                         g
                        }, by=.(intron_name)]


#  group overlapping motifs per sequence
fimo<-fimo[, .(motif_id=list(unique(motif_id)), motif_alt_id=list(unique(motif_alt_id)), start=list(start), end=list(end), score=list(score), pvalue=list(pvalue), qvalue=list(qvalue), seq=list(seq), rbp_group=list(unique(unlist(rbp_group)))), by=.(intron_name, o_group)][, o_group:=NULL]
fimo.control<-fimo.control[, .(motif_id=list(unique(motif_id)), motif_alt_id=list(unique(motif_alt_id)), start=list(start), end=list(end), score=list(score), pvalue=list(pvalue), qvalue=list(qvalue), seq=list(seq), rbp_group=list(unique(unlist(rbp_group)))), by=.(intron_name, o_group)][, o_group:=NULL]


#  keep on the side all motifs grouped by overlaps
fimo.all<-copy(fimo)
fimo.control.all<-copy(fimo.control)


#  identify number of non-overlapping motifs per sequence collapsing motif_ids, motif sequences and the RBPs involved
fimo<-fimo[, .(ncount=.N, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(intron_name)]
fimo.control<-fimo.control[, .(ncount=.N, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(intron_name)]


#  add gene_ids to the RBPs involved per circRNA
fimo[, gene_ids:=lapply(rbp_group, function(r){ hsa$gene_id[ match(r, hsa$gene_name) ] })]
fimo.control[, gene_ids:=lapply(rbp_group, function(r){ hsa$gene_id[ match(r, hsa$gene_name) ] })]


#  define circ_name out of intron_name
fimo[, circ_name:=sub('_up$|_down$', '', intron_name)]
fimo.all[, circ_name:=sub('_up$|_down$', '', intron_name)]
fimo.control[, circ_name:=sub('_[0-9]*$', '', intron_name)]
fimo.control.all[, circ_name:=sub('_[0-9]*$', '', intron_name)]


#  load introns
#  resolve upstream and downstream intron sequences and add the corresponding length to the FIMO objects
#  add corresponding lengths to the control introns as well
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_introns.RData')
fimo[ grepl('_up', intron_name), c('up', 'length'):=list(T , width(CIRCS.introns.up[ match(circ_name, CIRCS.introns.up$circ_name)]))]
fimo[ grepl('_down', intron_name), c('down', 'length'):=list(T , width(CIRCS.introns.down[ match(circ_name, CIRCS.introns.down$circ_name)]))]
fimo.control[, length:=width(CIRCS.controls.introns[ match(intron_name, CIRCS.controls.introns$intron_name) ])]


#  compute count densities per 1K of sequence 
#  N.B. densities of the upstream/downstream introns of the circRNAs remain separated
fimo[, density:=ncount/length*1e3]
fimo.control[, density:=ncount/length*1e3]

#}}}


#  save
save(ame, fimo, fimo.all, fimo.control, fimo.control.all, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results.RData



#  [circRNAs] analysis of the AME ATtRACT RBP motif enrichment results 
#{{{
rm(list=ls())
library(data.table)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(grid)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/draw_highlights.R')


#  load neuroblastoma-specific mRNA of RBPs results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues.RData')


#  load AME human RBP motif enrichment result
#  identify gene_names from the motif_id
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
ame[, gene_name:=list(strsplit(sub('_.*$', '', motif_id), ',', fixed=T))]


#  RBPs in neuroblastoma with significant enrichment of their motifs in circRNAs 
sig<-intersect(GENES$gene_name, unlist(ame$gene_name))
GENES<-GENES[ GENES$gene_name %in% sig, ]
rownames(GENES)<-NULL
nb.genes<-nb.genes[ GENES$gene_name, ]
hb.genes<-hb.genes[ GENES$gene_name, ]
vt.genes<-vt.genes[ GENES$gene_name, ]
ame<-ame[ sapply(gene_name, function(g){ any(g %in% sig) }), ]


#  order by NB expression
n<-order(rowMeans(nb.genes), decreasing=T)
nb.genes<-nb.genes[n, ]
hb.genes<-hb.genes[n, ]
vt.genes<-vt.genes[n, ]
stopifnot( all.equal( rownames(nb.genes), rownames(hb.genes) ) )
stopifnot( all.equal( rownames(nb.genes), rownames(vt.genes) ) )
rm(n)


#  identify those classified as groupI in Gerstberger et al. 2014 which are down-regulated upon brain development
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/validations/rbp/brainspan_census_of_human_RBPs_nrg3813/nrg3813-s7.RData')
GENES$groupI<-F
GENES$groupI[ GENES$gene_name %in% rbp[ groupI==T, gene_name] ]<-T


#  save for easy access
save(GENES, nb.genes, hb.genes, vt.genes, ame, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues_with_motifs_significantly_enriched_in_circRNAs.RData')


#  load back
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues_with_motifs_significantly_enriched_in_circRNAs.RData')


#  recycle
x11(width=25, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  heatmap on the log2-transformed (best contrast, compared to z-scores of raw CPMs, or z-scores of log2-transformed CPMs)
#  clustered by Euclidean distance in NB samples
X<-cbind(log2(1+nb.genes), log2(1+hb.genes), log2(1+vt.genes))  #  combine in one matrix
ex<-data.frame(tissue=factor(rep(c('neuroblastoma', 'brain tissue', 'various tumors'), c(ncol(nb.genes), ncol(hb.genes), ncol(vt.genes))), exclude=F), row.names=colnames(X))
cl<-setNames(list(setNames( c('seagreen4', 'cornflowerblue', 'coral4'), c('neuroblastoma', 'brain tissue', 'various tumors') )), colnames(ex))
hc<-hclust(dist(log2(1+nb.genes), method='euclidean'), method='ward.D2')  #  cluster based on distance between rows for NB only samples 
#hc<-hclust(dist(X, method='euclidean'), method='ward.D2')
grid.newpage()
ph<-pheatmap(X, color=viridis(10), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=F, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=F, 
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
grid.force()
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_RBPs_across_tissues_with_motifs_significantly_enriched_in_circRNAs.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(X,hc,ex,cl,ph)


#  boxplots of log2-transformed CPMs
B<-c(split(log2(1+nb.genes), row(nb.genes)), split(log2(1+hb.genes), row(hb.genes)), split(log2(1+vt.genes), row(vt.genes)))
n<-c(rbind(rbind(1:nrow(nb.genes), nrow(nb.genes)+(1:nrow(nb.genes)), 2*nrow(nb.genes)+(1:nrow(nb.genes)))))  #  order them in triplets
B<-B[n]
B.cl<-rep(c('seagreen4', 'cornflowerblue', 'coral4'), nrow(nb.genes))
groupI.cl<-setNames(rep('black', nrow(nb.genes)), rownames(nb.genes))
#groupI.cl[ GENES$gene_name[ GENES$groupI ] ]<-'darkorange'
par(mar=c(8.5, 7.0, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(-2, 2), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2](1+'CPM')), side=2, line=3, padj=-0.2, las=0, cex=2.4)
mtext(text=rownames(nb.genes), side=1, line=0, at= seq(2, length(B), 3), las=3, adj=1, cex=1.6, col=groupI.cl)
draw_highlights(L=length(B), STEP=3, YMAX=max(YTICK))
legend('topright', legend=c('neuroblastoma', 'brain tissue', 'various tumors'), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, xpd=T, y.intersp=0.60, x.intersp=0.2, seg.len=0.2)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_RBPs_across_tissues_with_motifs_significantly_enriched_in_circRNAs.svg', width=30, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(YTICK,B,n,B.cl,bp)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues_with_motifs_significantly_enriched_in_circRNAs.RData



#  [circRNA introns] analysis of the AME ATtRACT RBP motif enrichment results 
#{{{
rm(list=ls())
library(data.table)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/draw_highlights.R')


#  load neuroblastoma-specific mRNA of RBPs results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues.RData')


#  load AME human RBP motif enrichment result
#  identify gene_names from the motif_id
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results.RData')
ame[, gene_name:=list(strsplit(sub('_.*$', '', motif_id), ',', fixed=T))]


#  RBPs in neuroblastoma with also significant enrichment of their motifs in introns
sig<-intersect(GENES$gene_name, unlist(ame$gene_name))
GENES<-GENES[ GENES$gene_name %in% sig, ]
rownames(GENES)<-NULL
nb.genes<-nb.genes[ GENES$gene_name, ]
hb.genes<-hb.genes[ GENES$gene_name, ]
vt.genes<-vt.genes[ GENES$gene_name, ]
ame<-ame[ sapply(gene_name, function(g){ any(g %in% sig) }), ]


#  order by NB expression
n<-order(rowMeans(nb.genes), decreasing=T)
nb.genes<-nb.genes[n, ]
hb.genes<-hb.genes[n, ]
vt.genes<-vt.genes[n, ]
stopifnot( all.equal( rownames(nb.genes), rownames(hb.genes) ) )
stopifnot( all.equal( rownames(nb.genes), rownames(vt.genes) ) )
rm(n)


#  identify those classified as groupI in Gerstberger et al. 2014 which are down-regulated upon brain development
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/validations/rbp/brainspan_census_of_human_RBPs_nrg3813/nrg3813-s7.RData')
GENES$groupI<-F
GENES$groupI[ GENES$gene_name %in% rbp[ groupI==T, gene_name] ]<-T


#  save for easy access
save(GENES, nb.genes, hb.genes, vt.genes, ame, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues_with_motifs_significantly_enriched_in_introns.RData')


#  load back
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues_with_motifs_significantly_enriched_in_introns.RData')


#  recycle
x11(width=25, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  heatmap on the log2-transformed (best contrast, compared to z-scores of raw CPMs, or z-scores of log2-transformed CPMs)
#  clustered by Euclidean distance in NB samples
X<-cbind(log2(1+nb.genes), log2(1+hb.genes), log2(1+vt.genes))  #  combine in one matrix
ex<-data.frame(tissue=factor(rep(c('neuroblastoma', 'brain tissue', 'various tumors'), c(ncol(nb.genes), ncol(hb.genes), ncol(vt.genes))), exclude=F), row.names=colnames(X))
cl<-setNames(list(setNames( c('seagreen4', 'cornflowerblue', 'coral4'), c('neuroblastoma', 'brain tissue', 'various tumors') )), colnames(ex))
hc<-hclust(dist(log2(1+nb.genes), method='euclidean'), method='ward.D2')  #  cluster based on distance between rows for NB only samples 
#hc<-hclust(dist(X, method='euclidean'), method='ward.D2')
grid.newpage()
ph<-pheatmap(X, color=viridis(10), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=F, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=F, 
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
grid.force()
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_RBPs_across_tissues_with_motifs_significantly_enriched_in_introns.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(X,hc,ex,cl,ph)


#  boxplots of log2-transformed CPMs
B<-c(split(log2(1+nb.genes), row(nb.genes)), split(log2(1+hb.genes), row(hb.genes)), split(log2(1+vt.genes), row(vt.genes)))
n<-c(rbind(rbind(1:nrow(nb.genes), nrow(nb.genes)+(1:nrow(nb.genes)), 2*nrow(nb.genes)+(1:nrow(nb.genes)))))  #  order them in triplets
B<-B[n]
B.cl<-rep(c('seagreen4', 'cornflowerblue', 'coral4'), nrow(nb.genes))
groupI.cl<-setNames(rep('black', nrow(nb.genes)), rownames(nb.genes))
#groupI.cl[ GENES$gene_name[ GENES$groupI ] ]<-'darkorange'
par(mar=c(8.5, 7.0, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(-2, 2), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2](1+'CPM')), side=2, line=3, padj=-0.2, las=0, cex=2.4)
mtext(text=rownames(nb.genes), side=1, line=0, at= seq(2, length(B), 3), las=3, adj=1, cex=1.6, col=groupI.cl)
draw_highlights(L=length(B), STEP=3, YMAX=max(YTICK))
legend('topright', legend=c('neuroblastoma', 'brain tissue', 'various tumors'), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, xpd=T, y.intersp=0.60, x.intersp=0.2, seg.len=0.2)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_RBPs_across_tissues_with_motifs_significantly_enriched_in_introns.svg', width=40, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(YTICK,B,n,B.cl,bp)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues_with_motifs_significantly_enriched_in_introns.RData



#  [circRNAs+introns, FIMO results] distributions of non-overlapping RBP motif counts and count densities per 1kb of cumulative sequence
#                                   correlation of non-overlapping RBP motif counts and circRNA lengths
#                                   distribution comparison of non-overlapping RBP motif count densities per 1kb of cumulative sequence
#                                   enrichment of repeat elements in introns
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(pheatmap)
library(viridis)
library(grid)
library(RColorBrewer)
library(topGO)
library(org.Hs.eg.db)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/draw_highlights.R')


#  load FIMO results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
circs.fimo<-fimo
circs.fimo.all<-fimo.all
circs.fimo.control<-fimo.control
circs.fimo.control.all<-fimo.control.all
rm(fimo, fimo.all, fimo.control, fimo.control.all)


#  load FIMO human non-overlapping RBP motif counts 
#  recompute intron motif densities for cumulative intron sequences and counts (upstream + downstream for circRNAs, all control introns for controls)
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results.RData')
introns.fimo<-fimo
introns.fimo.all<-fimo.all
introns.fimo.control<-fimo.control
introns.fimo.control.all<-fimo.control.all
introns.fimo<-introns.fimo[, .(rbp_group=list(unique(unlist(rbp_group))), ncount=sum(ncount, na.rm=T), length=sum(length, na.rm=T)), by=.(circ_name)][, density:=ncount/length*1e3]
introns.fimo.control<-introns.fimo.control[, .(rbp_group=list(unique(unlist(rbp_group))), ncount=sum(ncount, na.rm=T), length=sum(length, na.rm=T)), by=.(circ_name)][, density:=ncount/length*1e3]
rm(fimo, fimo.all, fimo.control, fimo.control.all)
gc()


#  load circRNA exon and intron GRanges
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')
circs.exons<-CIRCS.exons
circs.exons.control<-CIRCS.controls
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_introns.RData')
circs.introns<-c(CIRCS.introns.up, CIRCS.introns.down)
circs.introns.controls<-CIRCS.controls.introns
rm(list=setdiff(ls(), c(l, 'circs.exons', 'circs.exons.control', 'circs.introns', 'circs.introns.controls')))


#  import RepeatMasker annotations
#
#  TRF : Tandem Repeat Finder
#
#  N.B. annotations overlap but we are going to keep all hits since they are anyway valid even for identical types, e.g.
#       if a TRF simple repeat overlaps with another within an intron range then this should count as two hits, not one.
#
#       This is not the same as in the RBP motif analysis where the RBP footprint ought to be taken into consideration rendering 
#       overlapping motifs practically impossible to be simultaneously bound.
#
rpmk<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/rpmk+simple_repeats.gtf')
seqlevels(rpmk, pruning.mode='coarse')<-seqlevels(unlist(circs.exons))


#  circRNAs basic analysis
#{{{

#  correlation of number of non-overlapping RBP motifs and circRNA length
circs.fimo[, cor(ncount, length, method='pearson')]          #  0.9817346438
circs.fimo.control[, cor(ncount, length, method='pearson')]  #  0.9731943631


#  explore the circRNAs-RBP pairs with no motifs found in the corresponding control sequences 
#{{{

x<-sapply(circs.fimo$circ_name, function(x){ ! any( circs.fimo.control[ circ_name %in% x, unlist(rbp_group)] %in% circs.fimo[ circ_name %in% x, unlist(rbp_group) ] ) }) 
x<-circs.fimo[ circ_name %in% names(x[x]) ][ order(-ncount) ]
x[, gene_id:=sub('\\.[0-9]*\\|.*$', '', circ_name)]


#  prepare for BP and MP enrichment analysis
grp<-x[, unique(gene_id)]
uni<-unique(sub('\\.[0-9]*\\|.*$', '', circs.fimo$circ_name)) 
allG<-setNames(rep(0, length(uni)), uni)
allG[ grp  ]<-1
stopifnot(  table(allG)[2] == length(grp) )
allG<-as.factor(allG)


#  molecular function GO enrichment
mf.topgo<-new('topGOdata', description='enrichment test', ontology='MF', allGenes=allG, annot=annFUN.org, mapping='org.Hs.eg.db', ID='Ensembl', nodeSize=1)
mf.test<-runTest(mf.topgo, algorithm='weight01', statistic='fisher')
mf<-data.table(GenTable(mf.topgo, p.value=mf.test, orderBy='p.value', topNodes=geneData(mf.test)['SigTerms'], numChar=120))
mf$p.value<-as.numeric(mf$p.value)
mf<-mf[ p.value<0.05 ]


#  biological process GO enrichment
bp.topgo<-new('topGOdata', description='enrichment test', ontology='BP', allGenes=allG, annot=annFUN.org, mapping='org.Hs.eg.db', ID='Ensembl', nodeSize=1)
bp.test<-runTest(bp.topgo, algorithm='weight01', statistic='fisher')
bp<-data.table(GenTable(bp.topgo, p.value=bp.test, orderBy='p.value', topNodes=geneData(bp.test)['SigTerms'], numChar=120))
bp$p.value<-as.numeric(bp$p.value)
bp<-bp[ p.value<0.05 ]

#}}}


#  recycle
x11(width=16, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  distribution of non-overlapping RBP motif counts
par(mar=c(5.0, 8.0, 1.0, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-circs.fimo[, ncount]
h<-hist(x, breaks=100, plot=F)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
h<-hist(x, breaks=h$breaks, col='grey39', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(pretty(x, 4)), add=F)
mtext('Frequency', side=2, line=6, padj=+0.1, las=0, cex=2.4)
mtext('Number of motifs', side=1, line=3, padj=+0.1, las=0, cex=2.4)
#dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_number_non-overlapping_RBP_motifs_circRNAs.svg', width=25, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x,h,YMAX)


#  distributions of non-overlapping RBP motif densities 
par(mar=c(5.0, 8.0, 1.0, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
i<-circs.fimo[, density]
j<-circs.fimo.control[, density]
#
wilcox.test(i, j, alternative='greater')$p.value  #  0.0516396029
#
#  histogram
#
hi<-hist(i, breaks=seq(min(i,j), max(i,j), length.out=100), plot=F)
hj<-hist(j, breaks=hi$breaks, plot=F)
YMAX<-max(pretty(c(1, max(hi$counts, hj$counts)), 4))
hj<-hist(j, breaks=hj$breaks, col='darkorange', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(pretty(c(i,j), 4)), add=F)
hi<-hist(i, breaks=hi$breaks, col=adjustcolor('blue4', alpha.f=0.65), border='white', add=T)
mtext('Frequency', side=2, line=6, padj=+0.1, las=0, cex=2.4)
mtext('Motif density', side=1, line=3, padj=+0.1, las=0, cex=2.4)
legend('topright', legend=c('circRNAs', 'controls'), col=c('blue4', 'darkorange'), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, y.intersp=0.8, x.intersp=0.5, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_density_non-overlapping_RBP_motifs_circRNAs.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  ecdf 
#
par(mar=c(5.0, 7.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
XLIM<-range(pretty(c(mean(i)-2*sd(i), mean(i)+2*sd(i)), 4))
hi<-curve(ecdf(i)(x), from=XLIM[1], to=XLIM[2], n=min(length(i), 100), ylab='', xlab='', pch=NA, col='blue4', lty=1, lwd=12, main='', xlim=c(XLIM[1], XLIM[2]), ylim=c(0, 1), xaxt='n')
hj<-curve(ecdf(j)(x), from=XLIM[1], to=XLIM[2], n=min(length(j), 100), ylab='', xlab='', pch=NA, col='darkorange', lty=1, lwd=12, add=T)
axis(1, at=pretty(range(hi$x)), labels=pretty(range(hi$x)))
mtext('Probability', side=2, line=5, padj=+0.1, las=0, cex=2.4)
mtext('Motif density', side=1, line=3, padj=+0.1, las=0, cex=2.4)
legend('topleft', legend=c('circRNAs', 'controls'), col=c('blue4', 'darkorange'), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, y.intersp=0.8, x.intersp=0.5, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_density_non-overlapping_RBP_motifs_circRNAs.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(i,j,hi,hj,YMAX,XLIM)

#}}}



#  introns basic analysis
#{{{

#  explore the intron-RBP pairs with no motifs found in the corresponding control sequences 
#{{{

x<-sapply(introns.fimo$circ_name, function(x){ ! any( introns.fimo.control[ circ_name %in% x, unlist(rbp_group)] %in% introns.fimo[ circ_name %in% x, unlist(rbp_group) ] ) }) 
x<-introns.fimo[ circ_name %in% names(x[x]) ][ order(-ncount) ]
x[, gene_id:=sub('\\.[0-9]*\\|.*$', '', circ_name)]


#  prepare for BP and MP enrichment analysis
grp<-x[, unique(gene_id)]
uni<-unique(sub('\\.[0-9]*\\|.*$', '', introns.fimo$circ_name)) 
allG<-setNames(rep(0, length(uni)), uni)
allG[ grp  ]<-1
stopifnot(  table(allG)[2] == length(grp) )
allG<-as.factor(allG)


#  molecular function GO enrichment
mf.topgo<-new('topGOdata', description='enrichment test', ontology='MF', allGenes=allG, annot=annFUN.org, mapping='org.Hs.eg.db', ID='Ensembl', nodeSize=1)
mf.test<-runTest(mf.topgo, algorithm='weight01', statistic='fisher')
mf<-data.table(GenTable(mf.topgo, p.value=mf.test, orderBy='p.value', topNodes=geneData(mf.test)['SigTerms'], numChar=120))
mf$p.value<-as.numeric(mf$p.value)
mf<-mf[ p.value<0.05 ]


#  biological process GO enrichment
bp.topgo<-new('topGOdata', description='enrichment test', ontology='BP', allGenes=allG, annot=annFUN.org, mapping='org.Hs.eg.db', ID='Ensembl', nodeSize=1)
bp.test<-runTest(bp.topgo, algorithm='weight01', statistic='fisher')
bp<-data.table(GenTable(bp.topgo, p.value=bp.test, orderBy='p.value', topNodes=geneData(bp.test)['SigTerms'], numChar=120))
bp$p.value<-as.numeric(bp$p.value)
bp<-bp[ p.value<0.05 ]

#}}}


#  recycle
x11(width=16, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  distribution of non-overlapping RBP motif counts
par(mar=c(5.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-introns.fimo[, ncount]
h<-hist(x, breaks=100, plot=F)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
h<-hist(x, breaks=h$breaks, col='grey39', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(pretty(x, 4)), add=F)
mtext('Frequency', side=2, line=6, padj=+0.1, las=0, cex=2.4)
mtext('Number of motifs', side=1, line=3, padj=+0.1, las=0, cex=2.4)
#dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_number_non-overlapping_RBP_motifs_introns.svg', width=25, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x,h,YMAX)


#  distributions of non-overlapping RBP motif densities 
par(mar=c(5.0, 8.0, 1.0, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
i<-introns.fimo[, density]
j<-introns.fimo.control[, density]
#
wilcox.test(i, j, alternative='less')$p.value  #  8.335119573e-06
#
#  histogram
#
hi<-hist(i, breaks=seq(min(i,j), max(i,j), length.out=100), plot=F)
hj<-hist(j, breaks=hi$breaks, plot=F)
YMAX<-max(pretty(c(1, max(hi$counts, hj$counts)), 4))
hj<-hist(j, breaks=hj$breaks, col='darkorange', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(pretty(c(i,j), 4)), add=F)
hi<-hist(i, breaks=hi$breaks, col=adjustcolor('blue4', alpha.f=0.65), border='white', add=T)
mtext('Frequency', side=2, line=6, padj=+0.1, las=0, cex=2.4)
mtext('Motif density', side=1, line=3, padj=+0.1, las=0, cex=2.4)
legend('topright', legend=c('circRNA introns', 'control introns'), col=c('blue4', 'darkorange'), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, y.intersp=0.8, x.intersp=0.5, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_density_non-overlapping_RBP_motifs_introns.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  ecdf 
#
par(mar=c(5.0, 7.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
XLIM<-range(pretty(c(mean(i)-2*sd(i), mean(i)+2*sd(i)), 4))
hi<-curve(ecdf(i)(x), from=XLIM[1], to=XLIM[2], n=min(length(i), 100), ylab='', xlab='', pch=NA, col='blue4', lty=1, lwd=12, main='', xlim=c(XLIM[1], XLIM[2]), ylim=c(0, 1), xaxt='n')
hj<-curve(ecdf(j)(x), from=XLIM[1], to=XLIM[2], n=min(length(j), 100), ylab='', xlab='', pch=NA, col='darkorange', lty=1, lwd=12, add=T)
axis(1, at=pretty(range(hi$x)), labels=pretty(range(hi$x)))
mtext('Probability', side=2, line=5, padj=+0.1, las=0, cex=2.4)
mtext('Motif density', side=1, line=3, padj=+0.1, las=0, cex=2.4)
legend('topleft', legend=c('circRNA introns', 'control introns'), col=c('blue4', 'darkorange'), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, y.intersp=0.8, x.intersp=0.5, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_density_non-overlapping_RBP_motifs_introns.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(i,j,hi,hj,YMAX,XLIM)

#}}}



introns_only<-setdiff(introns.fimo$circ_name, circs.fimo$circ_name)  #  ENSG00000164494.12|PDSS2_chr6-107193822-107197898
#
#  only the circRNA ENSG00000164494.12|PDSS2_chr6-107193822-107197898, which is 91nts long, does not have any significant RBP motif counts.
#
#  N.B. it is not the only small circRNA in the cohort


#  recycle
x11(width=16, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  [introns] enrichment of repeat elements 
#{{{

#  [~3mins] add all overlaps to the metadata
o.i<-findOverlaps(circs.introns, rpmk, select='all', type='any', minoverlap=6, ignore.strand=F)
o.i<-data.table(data.frame(rpmk[ subjectHits(o.i) ])[, c('gene_id'), drop=F])[, i:=queryHits(o.i)][, .(rpmk=list(unlist(gene_id))), by=.(i)]
m<-mcols(circs.introns)
m$rpmk<-vector('list', length(circs.introns))
for(n in o.i$i){ m$rpmk[n]<-o.i[ i==n, rpmk] }
mcols(circs.introns)<-m
o.i<-findOverlaps(circs.introns.controls, rpmk, select='all', type='any', minoverlap=6, ignore.strand=F)
o.i<-data.table(data.frame(rpmk[ subjectHits(o.i) ])[, c('gene_id'), drop=F])[, i:=queryHits(o.i)][, .(rpmk=list(unlist(gene_id))), by=.(i)]
m<-mcols(circs.introns.controls)
m$rpmk<-vector('list', length(circs.introns.controls))
for(n in o.i$i){ m$rpmk[n]<-o.i[ i==n, rpmk] }
mcols(circs.introns.controls)<-m
rm(m, o.i, n)


#  distributions of hits per 1Kb of intron sequence for signal and control
par(mar=c(5.0, 8.0, 1.0, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
i<-lengths(circs.introns$rpmk)/width(circs.introns)*1e3
j<-lengths(circs.introns.controls$rpmk)/width(circs.introns.controls)*1e3
#
wilcox.test(x=i, y=j, alternative='greater')$p.value  #  3.576e-91
#
#  histogram
#
hi<-hist(i, breaks=seq(min(i,j), max(i,j), length.out=100), plot=F)
hj<-hist(j, breaks=hi$breaks, plot=F)
YMAX<-max(pretty(c(1, max(hi$counts, hj$counts)), 4))
#XLIM<-range(pretty(c(i,j), 4))
XLIM<-c(0, 5)
hj<-hist(j, breaks=hj$breaks, col='darkorange', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=XLIM, add=F)
hi<-hist(i, breaks=hi$breaks, col=adjustcolor('blue4', alpha.f=0.65), border='white', add=T)
mtext('Frequency', side=2, line=6, padj=+0.1, las=0, cex=2.4)
mtext('Repeat elements per 1K', side=1, line=3, padj=+0.1, las=0, cex=2.4)
legend('topright', legend=c('circRNA introns', 'control introns'), col=c('blue4', 'darkorange'), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, y.intersp=0.8, x.intersp=0.5, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_rpmk_density_of_hits_introns.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  ecdf 
#
par(mar=c(5.0, 7.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#XLIM<-range(pretty(c(mean(i)-2*sd(i), mean(i)+2*sd(i)), 4))
XLIM<-c(-0.001, 3)
hi<-curve(ecdf(i)(x), from=XLIM[1], to=XLIM[2], n=min(length(i), 100), ylab='', xlab='', pch=NA, col='blue4', lty=1, lwd=12, main='', xlim=c(XLIM[1], XLIM[2]), ylim=c(0, 1), xaxt='n')
hj<-curve(ecdf(j)(x), from=XLIM[1], to=XLIM[2], n=min(length(j), 100), ylab='', xlab='', pch=NA, col='darkorange', lty=1, lwd=12, add=T)
#axis(1, at=pretty(range(hi$x)), labels=pretty(range(hi$x)))
axis(1, at=pretty(c(0, max(hi$x))), labels=pretty(c(0, max(hi$x))))
mtext('Probability', side=2, line=5, padj=+0.1, las=0, cex=2.4)
mtext('Repeat elements per 1K', side=1, line=3, padj=+0.1, las=0, cex=2.4)
legend('bottomright', legend=c('circRNA introns', 'control introns'), col=c('blue4', 'darkorange'), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, y.intersp=0.8, x.intersp=0.5, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_rpmk_density_of_hits_introns.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(i,j,hi,hj,YMAX,XLIM)

#}}}

#}}}



#  [ARID1A circRNA] RBP motif counts and corresponding mRNA expression across NB, human brain tissue, various tumors
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(grid)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/load_gene_expression.R')
source('~/bio/lib/draw_highlights.R')


#  run once
#{{{

#  load reference to add gene_ids to the gene_names
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]


#  non-overlapping RBP motif counts for circARID1A
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
rbp.motifs<-fimo.all[ grep('ARID1A', circ_name) ][, c('pvalue', 'qvalue', 'score', 'length'):=NULL]
rm(ame, fimo.all, fimo, fimo.control.all, fimo.control)
gc()


#  add gene_ids
#  identify the unique RBPs involved
#  count the number of times the unique RBPs appear in the non-overlapping motifs
rbp.motifs[, gene_ids:=lapply(rbp_group, function(r){ hsa$gene_id[ match(r, hsa$gene_name) ] })]
rbp<-unique(rbp.motifs[, c('rbp_group', 'gene_ids')][, .(gene_name=unlist(rbp_group), gene_id=unlist(gene_ids))])
rbp[, ncount:=rbp.motifs[, table(unlist(rbp_group))][ gene_name ]]
rbp<-rbp[ order(-ncount) ]


#  get CPMs across tissues/tumors
#  Mann-Whitney two-sided tests for neuroblastoma vs human brain and neuroblastoma vs various tumors asking for at least 50% of samples to 
#  express it and at least 30% of them to have CPM>=1
#  FDR-correct p-values
#  define one p-value to be the sum of the FDR-corrected p-values if both are significant (<0.05)
#{{{

#  get their CPMs across tissues/tumors
x<-load_gene_expression(data.frame(rbp), nb.tumors.only=T, vt.tumors.only=T)
nb.genes<-x$nb$genes
nb.meta<-x$nb$meta
hb.genes<-x$hb$genes
hb.meta<-x$hb$meta
vt.genes<-x$vt$genes
vt.meta<-x$vt$meta
rm(x,load_gene_expression)


#  [neuroblastoma] keep only tumors
#                  keep only CPMs and transpose
nb.genes<-t(nb.genes[ !is.na(nb.genes$risk_group), -tail(seq_len(ncol(nb.genes)),3)])


#  [human brain] keep only CPMs and transpose 
hb.genes<-t(hb.genes[, -tail(seq_len(ncol(hb.genes)),1) ])


#  [various tumors] keep only CPMs and transpose 
vt.genes<-t(vt.genes[, -tail(seq_len(ncol(vt.genes)),1) ])


#  go over each gene and do Mann-Whitney tests for the CPMs in NB vs brain tissue and NB vs various tumors
rbp$pv.hb<-rbp$pv.vt<-NA
for (n in seq_along(rbp$gene_name)){

    x<-nb.genes[ n, ]
    y<-hb.genes[ n, ]
    z<-vt.genes[ n, ]

    #  at least 50% of samples in the whole cohort should express it and at least 30% of those should have CPM>=1
    if(quantile(x, 0.7)<1 | sum(x!=0)/length(x)<0.5){ next }  

    #  Mann-Whitney test that NB are differentially expressed
    rbp$pv.hb[n]<-wilcox.test(x=x, y=y, alternative='two.sided')$p.value
    rbp$pv.vt[n]<-wilcox.test(x=x, y=z, alternative='two.sided')$p.value
}
rm(n,x,y,z)


#  FDR-adjust p-values
rbp$pv.hb<-p.adjust(rbp$pv.hb, method='BH')
rbp$pv.vt<-p.adjust(rbp$pv.vt, method='BH')


#  define a p-value which is the sum of p-values 
rbp$pvalue<-rbp$pv.hb + rbp$pv.vt 


#  order by p-value
rbp<-rbp[ order(pvalue) ]
nb.genes<-nb.genes[ rbp$gene_name, ]
hb.genes<-hb.genes[ rbp$gene_name, ]
vt.genes<-vt.genes[ rbp$gene_name, ]
stopifnot( all.equal( rownames(nb.genes), rbp$gene_name ) )
stopifnot( all.equal( rownames(hb.genes), rbp$gene_name ) )
stopifnot( all.equal( rownames(vt.genes), rbp$gene_name ) )

#}}}


#  save 
save(rbp.motifs, rbp, nb.genes, nb.meta, hb.genes, hb.meta, vt.genes, vt.meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_RBPs_across_tissues.RData')

#}}}


#  load above results:
#  
#    rbp.motifs : all non-overlapping RBP motifs found significantly enriched
#           rbp : all involved RBPs and their FDR-adjusted p-values for neuroblastoma-specific differential expression
#      nb.genes : CPMs in neuroblastoma tumors
#      hb.genes : CPMs in human brain tissue
#      vt.genes : CPMs in various tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_RBPs_across_tissues.RData')


#  recycle
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  order by decreasing counts
n<-order(rbp$ncount, decreasing=T)
rbp<-rbp[n, ]
nb.genes<-nb.genes[n, ]
hb.genes<-hb.genes[n, ]
vt.genes<-vt.genes[n, ]
rm(n)


#  heatmap on the log2-transformed (best contrast, compared to z-scores of raw CPMs, or z-scores of log2-transformed CPMs)
#  clustered by Euclidean distance in NB samples
X<-cbind(log2(1+nb.genes), log2(1+hb.genes), log2(1+vt.genes))  #  combine in one matrix
rownames(X)<-paste0(rbp$gene_name, ' (', rbp$ncount, ')')
col.an<-data.frame(tissue=factor(rep(c('neuroblastoma', 'brain tissue', 'various tumors'), c(ncol(nb.genes), ncol(hb.genes), ncol(vt.genes))), exclude=F), row.names=colnames(X))
row.an<-data.frame(neuroblastoma=factor(! (rbp$pvalue>=0.05 | is.na(rbp$pvalue)), levels=c(F, T), labels=c('not DE', 'DE')), row.names=rownames(X))
cl<-setNames(list(setNames( c('seagreen4', 'cornflowerblue', 'coral4'), c('neuroblastoma', 'brain tissue', 'various tumors') ),
                  setNames( c('steelblue3', 'darkgrey'), c('DE', 'not DE') )), 
             c(colnames(col.an), colnames(row.an)))
hc<-hclust(dist(log2(1+nb.genes), method='euclidean'), method='ward.D2')  #  cluster based on distance between rows for NB only samples 
#hc<-hclust(dist(X, method='euclidean'), method='ward.D2')
grid.newpage()
ph<-pheatmap(X, color=viridis(10), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=F, annotation_names_col=F, 
        annotation_col=col.an, annotation_row=row.an, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=F, 
        fontsize = 14, fontsize_row=14, fontsize_col=14, fontsize_number=5.0)
grid.force()
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_ARID1A_RBPs_across_tissues.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(X,hc,col.an,row.an,cl,ph)


#  boxplots of log2-transformed CPMs
B<-c(split(log2(1+nb.genes), row(nb.genes)), split(log2(1+hb.genes), row(hb.genes)), split(log2(1+vt.genes), row(vt.genes)))
n<-c(rbind(rbind(1:nrow(nb.genes), nrow(nb.genes)+(1:nrow(nb.genes)), 2*nrow(nb.genes)+(1:nrow(nb.genes)))))  #  order them in triplets
B<-B[n]
B.cl<-rep(c('seagreen4', 'cornflowerblue', 'coral4'), nrow(nb.genes))
par(mar=c(6.0, 7.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(-1, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2](1+'CPM')), side=2, line=3, padj=-0.2, las=0, cex=2.4)
mtext(text=paste0(rbp$gene_name, ' (', rbp$ncount, ')'), side=1, line=-1, at=seq(2, length(B), 3), las=2, adj=1.09, cex=1.2, col=ifelse( !( rbp$pvalue>=0.05 | is.na(rbp$pvalue) ), 'steelblue3', 'darkgrey'))
draw_highlights(L=length(B), STEP=3, YMAX=max(YTICK))
legend(x=par('usr')[1]*0.5, y=par('usr')[4]*1.07, legend=c('neuroblastoma', 'brain tissue', 'various tumors'), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
legend('topright', legend=c('DE', 'not DE'), col=c('steelblue3', 'darkgrey'), title='mRNA', bty='n', lty=1, lwd=15, pch=NA, cex=1.5, xpd=T, y.intersp=0.50, x.intersp=0.3, seg.len=0.3)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_ARID1A_RBPs_across_tissues.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(YTICK,B,n,B.cl,bp)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_RBPs_across_tissues.RData



#  [ARID1A introns] RBP motif counts and corresponding mRNA expression across NB, human brain tissue, various tumors
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(grid)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/load_gene_expression.R')
source('~/bio/lib/draw_highlights.R')


#  run once
#{{{

#  load reference to add gene_ids to the gene_names
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]


#  load intron non-overlapping RBP motif counts for circARID1A
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results.RData')
rbp.motifs<-fimo.all[ grep('ARID1A', circ_name) ][, c('pvalue', 'qvalue', 'score'):=NULL]
rm(ame, fimo.all, fimo, fimo.control.all, fimo.control)
gc()


#  add gene_ids
#  identify the unique RBPs involved
#  count the number of times the unique RBPs appear in the non-overlapping motifs
rbp.motifs[, gene_ids:=lapply(rbp_group, function(r){ hsa$gene_id[ match(r, hsa$gene_name) ] })]
rbp<-unique(rbp.motifs[, c('rbp_group', 'gene_ids')][, .(gene_name=unlist(rbp_group), gene_id=unlist(gene_ids))])
rbp[, ncount:=rbp.motifs[, table(unlist(rbp_group))][ gene_name ]]
rbp<-rbp[ order(-ncount) ]


#  get CPMs across tissues/tumors
#  Mann-Whitney two-sided tests for neuroblastoma vs human brain and neuroblastoma vs various tumors asking for at least 50% of samples to 
#  express it and at least 30% of them to have CPM>=1
#  FDR-correct p-values
#  define one p-value to be the sum of the FDR-corrected p-values if both are significant (<0.05)
#  save
#{{{

#  get their CPMs across tissues/tumors
x<-load_gene_expression(data.frame(rbp), nb.tumors.only=T, vt.tumors.only=T)
nb.genes<-x$nb$genes
nb.meta<-x$nb$meta
hb.genes<-x$hb$genes
hb.meta<-x$hb$meta
vt.genes<-x$vt$genes
vt.meta<-x$vt$meta
rm(x,load_gene_expression)


#  [neuroblastoma] keep only tumors
#                  keep only CPMs and transpose
nb.genes<-t(nb.genes[ !is.na(nb.genes$risk_group), -tail(seq_len(ncol(nb.genes)),3)])


#  [human brain] keep only CPMs and transpose 
hb.genes<-t(hb.genes[, -tail(seq_len(ncol(hb.genes)),1) ])


#  [various tumors] keep only CPMs and transpose 
vt.genes<-t(vt.genes[, -tail(seq_len(ncol(vt.genes)),1) ])


#  go over each gene and do Mann-Whitney tests for the CPMs in NB vs brain tissue and NB vs various tumors
rbp$pv.hb<-rbp$pv.vt<-NA
for (n in seq_along(rbp$gene_name)){

    x<-nb.genes[ n, ]
    y<-hb.genes[ n, ]
    z<-vt.genes[ n, ]

    #  at least 50% of samples in the whole cohort should express it and at least 30% of those should have CPM>=1
    if(quantile(x, 0.7)<1 | sum(x!=0)/length(x)<0.5){ next }  

    #  Mann-Whitney test that NB are differentially expressed
    rbp$pv.hb[n]<-wilcox.test(x=x, y=y, alternative='two.sided')$p.value
    rbp$pv.vt[n]<-wilcox.test(x=x, y=z, alternative='two.sided')$p.value
}
rm(n,x,y,z)


#  FDR-adjust p-values
rbp$pv.hb<-p.adjust(rbp$pv.hb, method='BH')
rbp$pv.vt<-p.adjust(rbp$pv.vt, method='BH')


#  define a p-value which is the sum of p-values 
rbp$pvalue<-rbp$pv.hb + rbp$pv.vt 


#  order by p-value
rbp<-rbp[ order(pvalue) ]
nb.genes<-nb.genes[ rbp$gene_name, ]
hb.genes<-hb.genes[ rbp$gene_name, ]
vt.genes<-vt.genes[ rbp$gene_name, ]
stopifnot( all.equal( rownames(nb.genes), rbp$gene_name ) )
stopifnot( all.equal( rownames(hb.genes), rbp$gene_name ) )
stopifnot( all.equal( rownames(vt.genes), rbp$gene_name ) )

#}}}


#  save for easy use
save(rbp.motifs, rbp, nb.genes, nb.meta, hb.genes, hb.meta, vt.genes, vt.meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_introns_RBPs_across_tissues.RData')

#}}}


#  load above results:
#  
#    rbp.motifs : all non-overlapping RBP motifs found significantly enriched
#           rbp : all involved RBPs and their FDR-adjusted p-values for neuroblastoma-specific differential expression
#      nb.genes : CPMs in neuroblastoma tumors
#      hb.genes : CPMs in human brain tissue
#      vt.genes : CPMs in various tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_introns_RBPs_across_tissues.RData')


#  recycle
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  order by decreasing counts
n<-order(rbp$ncount, decreasing=T)
rbp<-rbp[n, ]
nb.genes<-nb.genes[n, ]
hb.genes<-hb.genes[n, ]
vt.genes<-vt.genes[n, ]
rm(n)


#  heatmap on the log2-transformed (best contrast, compared to z-scores of raw CPMs, or z-scores of log2-transformed CPMs)
#  clustered by Euclidean distance in NB samples
X<-cbind(log2(1+nb.genes), log2(1+hb.genes), log2(1+vt.genes))  #  combine in one matrix
rownames(X)<-paste0(rbp$gene_name, ' (', rbp$ncount, ')')
col.an<-data.frame(tissue=factor(rep(c('neuroblastoma', 'brain tissue', 'various tumors'), c(ncol(nb.genes), ncol(hb.genes), ncol(vt.genes))), exclude=F), row.names=colnames(X))
row.an<-data.frame(neuroblastoma=factor(! (rbp$pvalue>=0.05 | is.na(rbp$pvalue)), levels=c(F, T), labels=c('not DE', 'DE')), row.names=rownames(X))
cl<-setNames(list(setNames( c('seagreen4', 'cornflowerblue', 'coral4'), c('neuroblastoma', 'brain tissue', 'various tumors') ),
                  setNames( c('steelblue3', 'darkgrey'), c('DE', 'not DE') )), 
             c(colnames(col.an), colnames(row.an)))
hc<-hclust(dist(log2(1+nb.genes), method='euclidean'), method='ward.D2')  #  cluster based on distance between rows for NB only samples 
#hc<-hclust(dist(X, method='euclidean'), method='ward.D2')
grid.newpage()
ph<-pheatmap(X, color=viridis(10), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=F, annotation_names_col=F, 
        annotation_col=col.an, annotation_row=row.an, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=F, 
        fontsize = 14, fontsize_row=14, fontsize_col=14, fontsize_number=5.0)
grid.force()
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_ARID1A_introns_RBPs_across_tissues.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(X,hc,col.an,row.an,cl,ph)


#  boxplots of log2-transformed CPMs
B<-c(split(log2(1+nb.genes), row(nb.genes)), split(log2(1+hb.genes), row(hb.genes)), split(log2(1+vt.genes), row(vt.genes)))
n<-c(rbind(rbind(1:nrow(nb.genes), nrow(nb.genes)+(1:nrow(nb.genes)), 2*nrow(nb.genes)+(1:nrow(nb.genes)))))  #  order them in triplets
B<-B[n]
B.cl<-rep(c('seagreen4', 'cornflowerblue', 'coral4'), nrow(nb.genes))
par(mar=c(8.5, 7.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(-1, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2](1+'CPM')), side=2, line=3, padj=-0.2, las=0, cex=2.4)
mtext(text=paste0(rbp$gene_name, ' (', rbp$ncount, ')'), side=1, line=0, at=seq(2, length(B), 3), las=2, adj=1, cex=1.2, col=ifelse( !( rbp$pvalue>=0.05 | is.na(rbp$pvalue) ), 'steelblue3', 'darkgrey'))
draw_highlights(L=length(B), STEP=3, YMAX=max(YTICK))
legend(x=par('usr')[1]*0.5, y=par('usr')[4]*1.07, legend=c('neuroblastoma', 'brain tissue', 'various tumors'), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=1.8, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
legend('topright', legend=c('DE', 'not DE'), col=c('steelblue3', 'darkgrey'), title='mRNA', bty='n', lty=1, lwd=15, pch=NA, cex=1.5, xpd=T, y.intersp=0.50, x.intersp=0.3, seg.len=0.3)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_ARID1A_introns_RBPs_across_tissues.svg', width=32, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(YTICK,B,n,B.cl,bp)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_introns_RBPs_across_tissues.RData



#  identify all RBPs with 90% motif similarity to the backspliced junction loci of all circRNAs
#{{{
rm(list=ls())
library(data.table)
library(GenomicFeatures)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)


#  [run once] create the circular junctions 6+6nts long
#             look for RBP hits for RBP motifs at least 7nts long
#{{{

#  [~5mins] create the circular junctions 6+6nts long
#           look for RBP hits for RBP motifs at least 7nts long
#{{{

#  load circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')


#  extract all exons per gene
#  exclude chrM exons
txdb<-loadDb('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/TxDb_GRCh38.gencode.v30.RData')
EXONS<-unlist(exonsBy(txdb, 'gene'))
EXONS$gene_id<-names(EXONS)
names(EXONS)<-NULL
seqlevels(EXONS, pruning.mode='coarse')<-seqlevels(CIRCS)


#  define exons upstream the 3' and downstream the 5' sites of the circular junction allowing for +1nt discrepancies and taking 
#  the exon with the longest part upstream of the 3' and downstream of the 5' sites 
#
#  a typical plus strand circRNA is defined as follows:
#
#            exon1               exon2           exon3
#      5'aaaaaaaaaaa-----------bbbbbbbb-------cccccccccccc3'
#
#  which would result in the circular junction:
#
#      cccccccccccc3'5'aaaaaaaaaaa
#      ==========================>
#
#  a typical minus strand circRNA is defined as follows:
#
#            exon3               exon2           exon1
#      3'ccccccccccc-----------bbbbbbbb-------aaaaaaaaaaaa5'
#
#  which would result in the circular junction:
#
#      aaaaaaaaaaaa5'3'ccccccccccc
#      <==========================
#
#  plus strand 3' exon (upstream):
P<-CIRCS[ strand(CIRCS) %in% '+' ]
o<-findOverlaps(resize(resize(P, width(P)+2, fix='center'), 2, fix='end'), EXONS, type='any', select='all')  #  allow +1nt discrepancies
o<-o[ P[ queryHits(o) ]$gene_id==EXONS[ subjectHits(o) ]$gene_id ]  #  make sure you stay within the gene of the corresponding circRNA
P<-P[ unique(queryHits(o)) ]  #  one isoform of CDR1 drops out
o<-findOverlaps(resize(resize(P, width(P)+2, fix='center'), 2, fix='end'), EXONS, type='any', select='all') 
o<-o[ P[ queryHits(o) ]$gene_id==EXONS[ subjectHits(o) ]$gene_id ]
o<-split(subjectHits(o), queryHits(o))
l<-GRangesList()
for(n in seq_along(o)){
    p<-P[ as.integer(names(o)[n]) ]
    ex<-EXONS[ o[[n]] ]
    end(ex)<-end(p)  #  fix the end of all exons to the end of the circRNA
    l[[n]]<-ex[which.max(width(ex))]  #  pick the first longest
}
stopifnot( all(lengths(l)==1) )  #  expecting exactly one exon assigned 
P.UP<-unlist(l)
stopifnot( all.equal( end(P.UP), end(P) ) )
#
#  plus strand 5' exon (downstream):
o<-findOverlaps(resize(resize(P, width(P)+2, fix='center'), 2, fix='start'), EXONS, type='any', select='all') 
o<-o[ P[ queryHits(o) ]$gene_id==EXONS[ subjectHits(o) ]$gene_id ]
stopifnot( length( setdiff( seq_along(P), queryHits(o) ) )==0 )  #  all circRNAs should overlap with at least one start exon as well
o<-split(subjectHits(o), queryHits(o))
l<-GRangesList()
for(n in seq_along(o)){
    p<-P[ as.integer(names(o)[n]) ]
    ex<-EXONS[ o[[n]] ]
    start(ex)<-start(p)  #  fix the start of all exons to the start of the circRNA
    l[[n]]<-ex[which.max(width(ex))]  #  pick the first longest
}
stopifnot( all(lengths(l)==1) )  #  expecting exactly one exon assigned 
P.DOWN<-unlist(l)
stopifnot( all.equal( start(P.DOWN), start(P) ) )
#
#  minus strand 3' exon (upstream) (N.B. fix='end' is transcription end, i.e. GRanges 'start'):
M<-CIRCS[ strand(CIRCS) %in% '-' ]
o<-findOverlaps(resize(resize(M, width(M)+2, fix='center'), 2, fix='end'), EXONS, type='any', select='all')
o<-o[ M[ queryHits(o) ]$gene_id==EXONS[ subjectHits(o) ]$gene_id ]  #  make sure you stay within the gene of the corresponding circRNA
M<-M[ unique(queryHits(o)) ]
o<-findOverlaps(resize(resize(M, width(M)+2, fix='center'), 2, fix='end'), EXONS, type='any', select='all') 
o<-o[ M[ queryHits(o) ]$gene_id==EXONS[ subjectHits(o) ]$gene_id ]
o<-split(subjectHits(o), queryHits(o))
l<-GRangesList()
for(n in seq_along(o)){
    m<-M[ as.integer(names(o)[n]) ]
    ex<-EXONS[ o[[n]] ]
    start(ex)<-start(m)  #  fix exon end positions to the circRNA end position (i.e. the GRanges 'start' position for minus strand features)
    l[[n]]<-ex[which.max(width(ex))]  #  pick the first longest
}
stopifnot( all(lengths(l)==1) )  #  expecting exactly one exon assigned 
M.UP<-unlist(l)
stopifnot( all.equal( start(M.UP), start(M) ) )
#
#  minus strand 5' exon (downstream) (N.B. fix='start' is transcription start, i.e. GRanges 'end'):
o<-findOverlaps(resize(resize(M, width(M)+2, fix='center'), 2, fix='start'), EXONS, type='any', select='all') 
o<-o[ M[ queryHits(o) ]$gene_id==EXONS[ subjectHits(o) ]$gene_id ]
stopifnot( length( setdiff( seq_along(M), queryHits(o) ) )==0 )  #  all circRNAs should overlap with at least one start exon as well
o<-split(subjectHits(o), queryHits(o))
l<-GRangesList()
for(n in seq_along(o)){
    m<-M[ as.integer(names(o)[n]) ]
    ex<-EXONS[ o[[n]] ]
    end(ex)<-end(m)  #  fix exon start positions to the circRNA start position (i.e. the GRanges 'end' position for minus strand features)
    l[[n]]<-ex[which.max(width(ex))]  #  pick the first longest
}
stopifnot( all(lengths(l)==1) )  #  expecting exactly one exon assigned 
M.DOWN<-unlist(l)
stopifnot( all.equal( end(M.DOWN), end(M) ) )
rm(o,n,m,p,ex,l)


#  add widths
P.UP$w<-width(P.UP)
P.DOWN$w<-width(P.DOWN)
M.UP$w<-width(M.UP)
M.DOWN$w<-width(M.DOWN)


#  keep only circular junctions of at least 12nts long
keep<-P.UP$w+P.DOWN$w>=12
cat('Number of plus strand circular junctions of at least 12nts = ', sum(keep), ' (number of dropouts = ', sum(!keep), ')\n', sep='')
P<-P[ keep ]
P.UP<-P.UP[ keep ]
P.DOWN<-P.DOWN[ keep ]
keep<-M.UP$w+M.DOWN$w>=12
cat('Number of minus strand circular junctions of at least 12nts = ', sum(keep), ' (number of dropouts = ', sum(!keep), ')\n', sep='')
M<-M[ keep ]
M.UP<-M.UP[ keep ]
M.DOWN<-M.DOWN[ keep ]
rm(keep)


#  trim exons to appropriate lengths so that 12nts circular junctions are created
x<-data.table(up=P.UP$w, down=P.DOWN$w)
stopifnot( nrow(x[ up<6 & down<6])==0 )  #  if one side is short of nts, the other should always compensate
x[ up<6, c('up.w', 'down.w'):=list(as.integer(up), as.integer(12-up)) ]
x[ down<6, c('up.w', 'down.w'):=list(as.integer(12-down), as.integer(down)) ]
x[ up>=6 & down>=6, c('up.w', 'down.w'):=list(as.integer(6), as.integer(6)) ]
P.UP<-resize(P.UP, x$up.w, fix='end')
mcols(P.UP)<-mcols(P.UP)[, 'gene_id', drop=F]
P.DOWN<-resize(P.DOWN, x$down.w, fix='start')
mcols(P.DOWN)<-mcols(P.DOWN)[, 'gene_id', drop=F]
stopifnot( all.equal(end(P.UP), end(P)) )
stopifnot( all.equal(start(P.DOWN), start(P)) )
x<-data.table(up=M.UP$w, down=M.DOWN$w)
stopifnot( nrow(x[ up<6 & down<6])==0 )  #  if one side is short of nts, the other should always compensate
x[ up<6, c('up.w', 'down.w'):=list(as.integer(up), as.integer(12-up)) ]
x[ down<6, c('up.w', 'down.w'):=list(as.integer(12-down), as.integer(down)) ]
x[ up>=6 & down>=6, c('up.w', 'down.w'):=list(as.integer(6), as.integer(6)) ]
M.UP<-resize(M.UP, x$up.w, fix='end')  #  transcription end is fixed, i.e. GRanges 'start'
mcols(M.UP)<-mcols(M.UP)[, 'gene_id', drop=F]
M.DOWN<-resize(M.DOWN, x$down.w, fix='start')
mcols(M.DOWN)<-mcols(M.DOWN)[, 'gene_id', drop=F]
stopifnot( all.equal(end(M.DOWN), end(M)) )
stopifnot( all.equal(start(M.UP), start(M)) )
rm(x)


#  reunite the strands
CIRCS<-c(P, M)
UP<-c(P.UP, M.UP)
DOWN<-c(P.DOWN, M.DOWN)
rm(P,M,P.UP,P.DOWN,M.UP,M.DOWN)


#  construct the circular junctions 
#
#  ----->!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!<-----
#  ----->!!!!!                                                                                     !!!!!<-----
#  ----->!!!!!  extractTranscriptSeqs() pastes together exons sequences IN THE ORDER THEY APPEAR   !!!!!<-----
#  ----->!!!!!  exonsBy(TxDb, 'tx') orders exons by rank which for minus-strand transcripts is     !!!!!<-----
#  ----->!!!!!  from right to left so the pasting would result in the correct transcript sequence  !!!!!<-----
#  ----->!!!!!                                                                                     !!!!!<-----
#  ----->!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!<-----
#
x<-c(UP, DOWN)
x$n<-rep(seq_along(UP), 2)
CIRCS.junctions<-split(x, x$n)  #  we checked that the order of the exons is preserved
names(CIRCS.junctions)<-paste0(CIRCS$gene_id, '|', CIRCS$gene_name,'_', as.character(seqnames(CIRCS)),as.character(strand(CIRCS)), start(CIRCS), '-', end(CIRCS))
CIRCS.junctions.seqs<-extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, CIRCS.junctions)
CIRCS$circ_name<-names(CIRCS.junctions)
rm(x,UP,DOWN)


#  save
save(CIRCS, CIRCS.junctions, CIRCS.junctions.seqs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_12nts_backspliced_junctions.RData')
writeXStringSet(CIRCS.junctions.seqs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_12nts_backspliced_junctions.fa', format='fasta')

#}}}


#  [~4h] look for hits for RBP motifs of at least 7nts long
#{{{

#  load the backspliced junctions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_12nts_backspliced_junctions.RData')


#  load the ATtRACT db
#  restrict it to motifs of at least 7nts long
#  convert to rows for nucleotide probabilities and columns for motif nucleotides
load('/fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')
pwm<-pwm[sapply(pwm, function(x){ nrow(x)>=7 })]
db<-db[id %in% names(pwm)][ match(id, names(pwm)), ]
pwm<-lapply(pwm, function(p){ p<-t(p);  rownames(p)<-c('A', 'C', 'G', 'T');  p})


#  go over each circRNA and each PWM and look for 90% similarities (at 7nts long motifs each hit will span the junction)
hits<-list()
for(i in seq_along(CIRCS.junctions.seqs)){
    hits[[names(CIRCS.junctions.seqs)[i]]]<-list()
    for(p in seq_along(pwm)){
        h<-unlist(as(matchPWM(pwm[[p]], CIRCS.junctions.seqs[[i]], min.score='90%'), 'CompressedIRangesList'))
        if(length(h)>0){
            mcols(h)<-extractAt(CIRCS.junctions.seqs[[i]], h)
            m<-DataFrame(data.frame(db[ id %in% names(pwm)[p], c('id', 'gene_name') ]))
            m$gene_name<-List(m$gene_name)
            mcols(h)<-cbind(mcols(h), m)
            colnames(mcols(h))<-c('motif', 'id', 'gene_name')
            hits[[ names(CIRCS.junctions.seqs)[i] ]]<-h
        }
        hits[[names(CIRCS.junctions.seqs)[i]]]<-unname(unlist(List(hits[[names(CIRCS.junctions.seqs)[i]]])))
    }
    cat('Finished backspliced junction', i, 'out of', length(CIRCS.junctions.seqs), '\n')
}
hits<-hits[ lengths(hits)>0 ]


#  save 
save(hits, db, pwm, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_12nts_backspliced_junctions_RBP_motif_hits.RData')

#}}}

#}}}


#  load results for further exploration?
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_12nts_backspliced_junctions.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_12nts_backspliced_junctions_RBP_motif_hits.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_12nts_backspliced_junctions.{RData,fa}
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_12nts_backspliced_junctions_RBP_motif_hits.RData



#  [neuroblastoma-specific circRNAs] RBP motif enrichments
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(grid)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/load_gene_expression.R')
source('~/bio/lib/draw_highlights.R')
source('~/bio/lib/enrichment_analysis.R')


#  functions
#{{{
markers_boxplot<-function(PICK, 
                          COLS=setNames(c('seagreen4','darkgreen','cornflowerblue','chocolate1','coral4'), c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA')), 
                          fileroot=''){
    #  internal function that depends on many global variables and plots a nice boxplot across risk groups

    x<-RBP.tpm[ RBP$gene_id[ match(PICK, RBP$gene_name) ], , drop=F]
    colnames(x)<-nb.meta[ match(colnames(x), bid), risk_group ]
    B<-setNames(lapply(c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'), function(r){ c(x[, colnames(x) %in% r]) }), c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))


    B.cl<-COLS[ names(B) ] 
    YTICK<-pretty(c(0.0, sapply(B, max)), 5)
    plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
    bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
    mtext('TPM', side=2, line=5, padj=-0.1, las=0, cex=2.4)
    mtext(text=paste0(names(B), ' (', lengths(B)/nrow(x), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
    mtext(text=paste(PICK, sep='', collapse=','), side=3, line=-1, padj=+0.2, cex=1.8)

    if (fileroot!=''){
    dev.print(device=svg, file=paste0(fileroot, paste(PICK, sep='', collapse=','), '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    }
}

#}}}


#  run once
#{{{

#  load reference to add gene_ids to the gene_names
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]


#  load the neuroblastoma-specific circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_neuroblastoma-specific.RData')
nb.circs<-CIRCS
rm(CIRCS, CIRCS.hg19, hb.circs, vt.circs, nb.genes, hb.genes, vt.genes, hb.meta, vt.meta)


#  load all RBPs and compute their TPMs across tumors
#{{{

env<-new.env()
local(load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/RBPs_across_tissues.RData'), envir=env)
RBP<-get('GENES', envir=env)
rm(env)


#  load featureCounts of all tumors
#  compute TPMs
#  convert to matrix
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')
RBP.tpm<-unlist(totalrna[ nb.meta[, bid] ])
RBP.tpm$bid<-sub('\\.[0-9]*$', '', rownames(RBP.tpm))
rownames(RBP.tpm)<-NULL
RBP.tpm<-data.table(RBP.tpm)[, .(gene_id=gene_id, length=length, counts=counts, tpm=counts/length/sum(counts/length)*1e6), by=.(bid)]
RBP.tpm<-dcast(RBP.tpm, bid ~ gene_id, value.var='tpm', fun.aggregate=sum)
RBP.tpm<-t(data.frame(RBP.tpm[, -1], row.names=RBP.tpm[, bid], check.names=F))
RBP.tpm<-RBP.tpm[ rownames(RBP.tpm) %in% RBP$gene_id, ]
rm(totalrna)

#}}}


#  identify the non-overlapping RBP motif counts for the neuroblastoma-specific circRNAs
#  add gene_ids to the RBP gene_names
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
rbp.motifs<-fimo.all[ circ_name %in% nb.circs$circ_name ][, c('pvalue', 'qvalue', 'score', 'length'):=NULL]
rbp.motifs[, gene_ids:=lapply(rbp_group, function(r){ hsa$gene_id[ match(r, hsa$gene_name) ] })]
rm(ame, fimo.all, fimo, fimo.control.all, fimo.control)
gc()


#  for each circRNA and RBP count the number of non-overlaping motifs it is found enriched in
#  N.B. RBP counts are shared across motifs now, i.e. RBP1 and RBP2 that had a common PWM from a non-overlapping motif will both get a +1 contribution 
#       from that motif in their corresponding counts. On the other hand, the counts for each circRNA and RBP are still from non-overlapping motifs!
rbp<-rbp.motifs[, c('circ_name', 'rbp_group')][, .(gene_name=unlist(rbp_group)), by=.(circ_name)][, .(ncount=.N), by=.(circ_name, gene_name)]
rbp<-rbp[ order(-ncount) ]


#  GOBP, MSigDB C2 enrichment analysis using as background all RBPs
#{{{

#  count number of circRNAs with enrichments in given RBP's motifs
B<-rbp[, sort(table(gene_name), decreasing=T)]
B<-setNames(as.numeric(B), names(B))


#  define signal group to be those with at least 5 circRNA counts
G<-setNames( sub('\\.[0-9]*$', '', RBP[ match(names(B[B>=5]), RBP$gene_name), 'gene_id' ]), names(B[B>=5]) )
U<-setNames( sub('\\.[0-9]*$', '', RBP$gene_id), RBP$gene_name )
rbp.goms<-enrichment_analysis(G=G, U=U, MS=setNames('/fast/groups/ag_schulte/work/reference/annotation/MSigDB/c2.all.v7.0.symbols.gmt', 'C2'))
rm(G,U,B)

#}}}


#  [date checked: 20201030] check NCBI for neuroblastoma-based publications associated with the RBPs 
#
#                           To individually check a PMID:
#
#                               https://www.ncbi.nlm.nih.gov/pubmed/?term=31436914%5Bpmid%5D
#
#                           XML entry of given PMID:
#
#                               https://www.ncbi.nlm.nih.gov/pubmed/?term=28560387%5Buid%5D&report=xml&format=text
#{{{

PREFIX<-'https://www.ncbi.nlm.nih.gov/pubmed?term=(neuroblastoma%5BTitle%2FAbstract%5D)%20AND%20'
SUFFIX<-'%5BTitle%2FAbstract%5D&report=uilist&format=text'
PMID<-list()
for(p in rbp[, unique(gene_name)]){
    NCBI<-paste0(PREFIX, p, SUFFIX, collapse='')
    cat('\nquerying', p, '...')
    PMID[[p]]<-setdiff(gsub('<.?pre>', '', system(paste0("wget -O - --quiet '", NCBI, "' | sed -n '/^<pre>/,/^<\\/pre>/p'"), intern=T)), '')
}
PMID<-PMID[ lengths(PMID)>0 ]
rm(PREFIX,SUFFIX)


#  check them all at once
browseURL(paste0('https://www.ncbi.nlm.nih.gov/pubmed/', paste0(unique(unlist(PMID)), collapse=',')))


#  check only those with at least 3 motifs
browseURL(paste0('https://www.ncbi.nlm.nih.gov/pubmed/', paste0(unique(unlist(PMID[ rbp[ncount>=3, unique(gene_name)] ])), collapse=',')))

#}}}


#  save 
save(RBP, RBP.tpm, rbp, rbp.motifs, rbp.goms, nb.circs, nb.meta, PMID, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_neuroblastoma-specific_RBPs.RData')

#}}}


#  load back
#  get the counts of circRNAs enriched in give RBP's motifs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_neuroblastoma-specific_RBPs.RData')
B<-rbp[, sort(table(gene_name), decreasing=T)]
B<-setNames(as.numeric(B), names(B))


#  percentage of circRNAs enriched in given RBP's motifs
B.p<-100*B/rbp[, length(unique(circ_name))]
browseURL(paste0('https://www.ncbi.nlm.nih.gov/pubmed/', paste0(unique(unlist(PMID[ names(PMID) %in% names(B.p[ B.p>=30 ])])), collapse=',')))


#  barplot of percentage of circRNAs enriched in given RBP's motifs
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(9.0, 7.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#YTICK<-pretty(c(0, max(B)), 4)
YTICK<-pretty(c(0, 100), 4)
bp<-barplot(B.p, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(YTICK[1], tail(YTICK, 1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B.p, border='white', col='darkgrey', axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
axis(2, at=YTICK, line=0, cex.axis=2.4)
#mtext('Count', side=2, line=5, padj=+0.3, las=0, cex=2.4)
mtext('Percentage', side=2, line=5, padj=+0.3, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=bp, las=2, adj=1, cex=1.8, col='black')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_circRNAs_neuroblastoma-specific_RBPs.svg', width=45, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(14.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  check RBP expression across risk-groups for all RBPs
for(p in names(B)){ cat('plotting:', p, '\n'); markers_boxplot(PICK=p, fileroot=''); readline() }


#  plot and save specific ones
PICK<-c('PTBP1', 'SFPQ', 'YBX1', 'SRSF3', 'HNRNPA1', 'KHSRP', 'HNRNPK', 'SNRPA', 'ZFP36', 'TIA1', 'RBM5', 'TARDBP','ELAVL4', 'AGO2', 'SRSF4', 'EIF4B', 'PABPC3', 'SNRNP70', 'RBM28', 'IGF2BP1', 'RBM24', 'AGO1', 'AKAP1', 'CELF4', 'CELF5', 'CELF6', 'G3BP2', 'RBM3', 'RBM6', 'NONO', 'PABPC1')
for(p in PICK){ 
    cat('plotting and saving:', p, '(', B[p], ')\n')
    markers_boxplot(PICK=p, fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_') 
    readline('')
}


#  hierachical clustering of tumors
#{{{

#  all RBPs
m<-nb.meta
x<-RBP.tpm[, m$bid]
colnames(x)<-sub('-11-R01', '', colnames(x))
d<-dist(t(x), method='euclidean')
hc<-hclust(d, method='ward.D2')
ex<-data.frame(Type=factor(m$risk_group, exclude=F), row.names=colnames(x))
cl<-setNames(list(setNames( unique(m$col), unique(m$risk_group) )), colnames(ex))
grid.newpage()
ph<-pheatmap(as.matrix(d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
grid.force()
#dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_tpm_RBPs_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)



#  above-picked RBPs
m<-nb.meta
x<-RBP.tpm[rownames(RBP.tpm) %in% RBP$gene_id[ match(PICK, RBP$gene_name) ], m$bid, drop=F]
colnames(x)<-sub('-11-R01', '', colnames(x))
d<-dist(t(x), method='euclidean')
hc<-hclust(d, method='ward.D2')
ex<-data.frame(Type=factor(m$risk_group, exclude=F), row.names=colnames(x))
cl<-setNames(list(setNames( unique(m$col), unique(m$risk_group) )), colnames(ex))
grid.newpage()
ph<-pheatmap(as.matrix(d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
grid.force()
rm(m,x,d,ex,hc,ph,cl)

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_neuroblastoma-specific_RBPs.RData



#  [SFPQ] circRNAs with motif enrichments
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(grid)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  run once
#{{{

#  load neuroblastoma-specific circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_neuroblastoma-specific.RData')
nb.specific<-CIRCS
rm(list=setdiff(ls(), 'nb.specific'))


#  load circRNA RBP motif enrichment results and isolate SFPQ-related ones
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
ci.lengths<-fimo[, c('circ_name', 'length')]
ci.all<-fimo.all[ sapply(rbp_group, function(r){ any(grepl('SFPQ', r)) }), ]
rm(fimo, fimo.all, fimo.control, fimo.control.all, ame)


#  collect all motifs per circRNA
ci.iso<-ci.all[, .(ncount=.N, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(circ_name)][ order(-ncount) ]


#  neuroblastoma-specific
ci.nb<-ci.iso[ sapply(circ_name, function(x){ any(x %in% nb.specific$circ_name) }) ]
ci.nb<-ci.nb[, c('g', 'rest'):=tstrsplit(circ_name, '_', fixed=T)][, c('gene_id', 'gene_name'):=tstrsplit(g, '|', fixed=T)][, c('g', 'rest'):=NULL]


#  identify corresponding gene_name and gene_id 
#  summarize circRNA isoforms by maximum and mean
ci<-ci.iso[, c('g', 'rest'):=tstrsplit(circ_name, '_', fixed=T)][, c('gene_id', 'gene_name'):=tstrsplit(g, '|', fixed=T)][, c('g', 'rest'):=NULL][, .(max=max(ncount), mean=mean(ncount), motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group))), circ_name=list(unique(unlist(circ_name)))), by=.(gene_id, gene_name)][ order(-max) ]


#  save
save(nb.specific, ci.nb, ci, ci.all, ci.iso, ci.lengths, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_enriched_in_SFPQ.RData')

#}}}


#  load the results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_enriched_in_SFPQ.RData')


#  recycle
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(9.0, 7.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  [maximum motifs found across the isoforms] barplot 
CUT<-ci[, quantile(max, 0.95)]
B<-setNames(ci[max>=CUT, max], ci[max>=CUT, gene_name])
B.cl<-rep('black', length(B))
B.cl[ names(B) %in% ci.nb[ncount>=CUT, gene_name] ]<-'steelblue3'
YTICK<-pretty(c(0, max(B)), 4)
bp<-barplot(B, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(YTICK[1], tail(YTICK, 1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col=B.cl, axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext('Counts', side=2, line=5, padj=+0.3, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=bp, las=2, adj=1, cex=1.8, col=B.cl)
#  WILL LOOK CUT ON DEVICE!
legend(x=par('usr')[2]*0.60, y=par('usr')[4]*1.05, legend=c('neuroblastoma-specific'), col=c('steelblue3'), bty='n', lty=1, lwd=25, pch=NA, cex=2.0, xpd=T, y.intersp=0.50, x.intersp=0.2, seg.len=0.2)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_circRNAs_enriched_in_SFPQ_motifs.svg', width=40, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_enriched_in_SFPQ.RData



#  [KHSRP] circRNAs with motif enrichments
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(grid)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  run once
#{{{

#  load neuroblastoma-specific circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_neuroblastoma-specific.RData')
nb.specific<-CIRCS
rm(list=setdiff(ls(), 'nb.specific'))


#  load circRNA RBP motif enrichment results and isolate KHSRP-related ones
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
ci.lengths<-fimo[, c('circ_name', 'length')]
ci.all<-fimo.all[ sapply(rbp_group, function(r){ any(grepl('KHSRP', r)) }), ]
rm(fimo, fimo.all, fimo.control, fimo.control.all, ame)


#  collect all motifs per circRNA
ci.iso<-ci.all[, .(ncount=.N, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(circ_name)][ order(-ncount) ]


#  neuroblastoma-specific
ci.nb<-ci.iso[ sapply(circ_name, function(x){ any(x %in% nb.specific$circ_name) }) ]
ci.nb<-ci.nb[, c('g', 'rest'):=tstrsplit(circ_name, '_', fixed=T)][, c('gene_id', 'gene_name'):=tstrsplit(g, '|', fixed=T)][, c('g', 'rest'):=NULL]


#  identify corresponding gene_name and gene_id 
#  summarize circRNA isoforms by maximum and mean
ci<-ci.iso[, c('g', 'rest'):=tstrsplit(circ_name, '_', fixed=T)][, c('gene_id', 'gene_name'):=tstrsplit(g, '|', fixed=T)][, c('g', 'rest'):=NULL][, .(max=max(ncount), mean=mean(ncount), motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group))), circ_name=list(unique(unlist(circ_name)))), by=.(gene_id, gene_name)][ order(-max) ]


#  save
save(nb.specific, ci.nb, ci, ci.all, ci.iso, ci.lengths, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_enriched_in_KHSRP.RData')

#}}}


#  load the results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_enriched_in_KHSRP.RData')


#  recycle
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(9.0, 7.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  [maximum motifs found across the isoforms] barplot 
CUT<-ci[, quantile(max, 0.95)]
B<-setNames(ci[max>=CUT, max], ci[max>=CUT, gene_name])
B.cl<-rep('black', length(B))
B.cl[ names(B) %in% ci.nb[ncount>=CUT, gene_name] ]<-'steelblue3'
YTICK<-pretty(c(0, max(B)), 4)
bp<-barplot(B, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(YTICK[1], tail(YTICK, 1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col=B.cl, axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext('Counts', side=2, line=5, padj=+0.3, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=bp, las=2, adj=1, cex=1.8, col=B.cl)
#  WILL LOOK CUT ON DEVICE!
legend(x=par('usr')[2]*0.60, y=par('usr')[4]*1.05, legend=c('neuroblastoma-specific'), col=c('steelblue3'), bty='n', lty=1, lwd=25, pch=NA, cex=2.0, xpd=T, y.intersp=0.50, x.intersp=0.2, seg.len=0.2)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_circRNAs_enriched_in_KHSRP_motifs.svg', width=30, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_enriched_in_KHSRP.RData



#  [ELAVL1-3] circRNAs with motif enrichments
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(grid)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  run once
#{{{

#  load neuroblastoma-specific circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_neuroblastoma-specific.RData')
nb.specific<-CIRCS
rm(list=setdiff(ls(), 'nb.specific'))


#  load circRNA RBP motif enrichment results and isolate ELAVL-related ones
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
ci.lengths<-fimo[, c('circ_name', 'length')]
ci.all<-fimo.all[ sapply(rbp_group, function(r){ any(grepl('ELAVL', r)) }), ]
rm(fimo, fimo.all, fimo.control, fimo.control.all, ame)


#  collect all motifs per circRNA
ci.iso<-ci.all[, .(ncount=.N, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(circ_name)][ order(-ncount) ]


#  neuroblastoma-specific
ci.nb<-ci.iso[ sapply(circ_name, function(x){ any(x %in% nb.specific$circ_name) }) ]
ci.nb<-ci.nb[, c('g', 'rest'):=tstrsplit(circ_name, '_', fixed=T)][, c('gene_id', 'gene_name'):=tstrsplit(g, '|', fixed=T)][, c('g', 'rest'):=NULL]


#  identify corresponding gene_name and gene_id 
#  summarize circRNA isoforms by maximum and mean
ci<-ci.iso[, c('g', 'rest'):=tstrsplit(circ_name, '_', fixed=T)][, c('gene_id', 'gene_name'):=tstrsplit(g, '|', fixed=T)][, c('g', 'rest'):=NULL][, .(max=max(ncount), mean=mean(ncount), motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group))), circ_name=list(unique(unlist(circ_name)))), by=.(gene_id, gene_name)][ order(-max) ]


#  save
save(nb.specific, ci.nb, ci, ci.all, ci.iso, ci.lengths, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_enriched_in_ELAVL.RData')

#}}}


#  load the results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_enriched_in_ELAVL.RData')


#  recycle
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(9.0, 7.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  [maximum motifs found across the isoforms] barplot 
CUT<-ci[, quantile(max, 0.95)]
B<-setNames(ci[max>=CUT, max], ci[max>=CUT, gene_name])
B.cl<-rep('black', length(B))
B.cl[ names(B) %in% ci.nb[ncount>=CUT, gene_name] ]<-'steelblue3'
YTICK<-pretty(c(0, max(B)), 4)
bp<-barplot(B, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(YTICK[1], tail(YTICK, 1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col=B.cl, axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext('Counts', side=2, line=5, padj=+0.3, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=bp, las=2, adj=1, cex=1.8, col=B.cl)
#  WILL LOOK CUT ON DEVICE!
legend(x=par('usr')[2]*0.60, y=par('usr')[4]*1.05, legend=c('neuroblastoma-specific'), col=c('steelblue3'), bty='n', lty=1, lwd=25, pch=NA, cex=2.0, xpd=T, y.intersp=0.50, x.intersp=0.2, seg.len=0.2)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_circRNAs_enriched_in_ELAVL_motifs.svg', width=40, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_enriched_in_ELAVL.RData



#  [MNA vs HR_nMNA] motif counts for the DE splice factors on the DE circRNA exons/introns
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(grid)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load DE results for genes and circRNAs and split to up-/downregulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_MNA_HR_nMNA.RData')
up.gns<-subset(RES, log2FoldChange>0)
down.gns<-subset(RES, log2FoldChange<0)
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs_MNA_HR_nMNA.RData')
up.circs<-subset(RES, log2FoldChange>0)
down.circs<-subset(RES, log2FoldChange<0)
rm(RES, CND, DDS, PCA, VE, VSC)


#  load union of splice-factors and identify the up-/downregulated 
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/SpliceAid-F+GO.RData')
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/MSigDB_c2_v7.0_splicing.RData')
db<-unique(rbind(data.frame(db[, c('gene_name', 'gene_id')]), msigdbSF[, c('gene_name', 'gene_id')]))
up.sf<-up.gns[ rownames(up.gns) %in% db$gene_id, , drop=F]
down.sf<-down.gns[ rownames(down.gns) %in% db$gene_id, , drop=F]
rm(db,msigdbSF)


#  load circRNA sequence RBP enrichments
#  keep only cases with 3 or more significant motifs 
#  summarize counts at the gene level by taking the mean and median
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
rbp.exons<-fimo[ ncount>=3 ][, c('g', 'rest'):=tstrsplit(circ_name, '_', fixed=T)][, c('gene_id', 'gene_name'):=tstrsplit(g, '|', fixed=T)][, c('g', 'rest'):=NULL][, .(circ_names=list(circ_name), rbp_ids=list(unique(unlist(gene_ids))), mean=mean(as.numeric(unlist(ncount))), median=median(as.numeric(unlist(ncount)))), by=.(gene_id, gene_name)]   #  as.numeric() because it complains about integer type
rm(ame, fimo.control, fimo.control.all, fimo.all, fimo)


#  load circRNA flanking introns RBP enrichments
#  keep only flanking intron cases with 3 or more significant motifs 
#  summarize counts of both flanking introns and isoforms at the gene level by taking the mean and median
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results.RData')
rbp.introns<-fimo[ ncount>=3 ][, c('g', 'rest'):=tstrsplit(circ_name, '_', fixed=T)][, c('gene_id', 'gene_name'):=tstrsplit(g, '|', fixed=T)][, c('g', 'rest'):=NULL][, .(circ_names=list(circ_name), rbp_ids=list(unique(unlist(gene_ids))), mean=mean(as.numeric(unlist(ncount))), median=median(as.numeric(unlist(ncount)))), by=.(gene_id, gene_name)]   #  as.numeric() because it complains about integer type
rm(ame, fimo.control, fimo.control.all, fimo.all, fimo)


#  [motif counts of significantly upregulated splice factor on significantly upregulated circRNA exons] Mann-Whitney U test
x.grp<-rownames(subset(up.circs, padj<0.05))                                            #  significantly upregulated circRNAs
y.grp<-rownames(subset(up.sf, padj<0.05))                                               #  significantly upregulated SF
y.rest<-c(rownames(down.sf), rownames(up.sf[ ! rownames(up.sf) %in% y.grp, , drop=F]))  #  not significantly upregulated SF
wilcox.test(x=rbp.exons[ gene_id %in% x.grp & sapply(rbp_ids, function(r){ any(r %in% y.grp)}), unlist(mean)], 
            y=rbp.exons[ gene_id %in% x.grp & sapply(rbp_ids, function(r){ any(r %in% y.rest)}), unlist(mean)], alternative='greater')$p.value
rm(x.grp, y.grp, y.rest)
#
#  => 0.503874809


#  [motif counts of significantly downregulated splice factor on significantly downregulated circRNA exons] Mann-Whitney U test
x.grp<-rownames(subset(down.circs, padj<0.05))                                            #  significantly downregulated circRNAs
y.grp<-rownames(subset(down.sf, padj<0.05))                                               #  significantly downregulated SF
y.rest<-c(rownames(up.sf), rownames(down.sf[ ! rownames(down.sf) %in% y.grp, , drop=F]))  #  not significantly downregulated SF
wilcox.test(x=rbp.exons[ gene_id %in% x.grp & sapply(rbp_ids, function(r){ any(r %in% y.grp)}), unlist(mean)], 
            y=rbp.exons[ gene_id %in% x.grp & sapply(rbp_ids, function(r){ any(r %in% y.rest)}), unlist(mean)], alternative='greater')$p.value
rm(x.grp, y.grp, y.rest)
#
#  => 0.0009768746504


#  [motif counts of significantly upregulated splice factor on significantly upregulated circRNA introns] Mann-Whitney U test
x.grp<-rownames(subset(up.circs, padj<0.05))                                            #  significantly upregulated circRNAs
y.grp<-rownames(subset(up.sf, padj<0.05))                                               #  significantly upregulated SF
y.rest<-c(rownames(down.sf), rownames(up.sf[ ! rownames(up.sf) %in% y.grp, , drop=F]))  #  not significantly upregulated SF
wilcox.test(x=rbp.introns[ gene_id %in% x.grp & sapply(rbp_ids, function(r){ any(r %in% y.grp)}), unlist(mean)], 
            y=rbp.introns[ gene_id %in% x.grp & sapply(rbp_ids, function(r){ any(r %in% y.rest)}), unlist(mean)], alternative='greater')$p.value
rm(x.grp, y.grp, y.rest)
#
#  => 0.5038733179


#  [motif counts of significantly downregulated splice factor on significantly downregulated circRNA introns] Mann-Whitney U test
x.grp<-rownames(subset(down.circs, padj<0.05))                                            #  significantly downregulated circRNAs
y.grp<-rownames(subset(down.sf, padj<0.05))                                               #  significantly downregulated SF
y.rest<-c(rownames(up.sf), rownames(down.sf[ ! rownames(down.sf) %in% y.grp, , drop=F]))  #  not significantly downregulated SF
wilcox.test(x=rbp.introns[ gene_id %in% x.grp & sapply(rbp_ids, function(r){ any(r %in% y.grp)}), unlist(mean)], 
            y=rbp.introns[ gene_id %in% x.grp & sapply(rbp_ids, function(r){ any(r %in% y.rest)}), unlist(mean)], alternative='greater')$p.value
rm(x.grp, y.grp, y.rest)
#
#  => 0.5000592517

#}}}



#  [MNA vs HR_nMNA] motif counts for the DE splice factors on the neuroblastoma-specific circRNA exons/introns
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(grid)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load the neuroblastoma-specific circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_neuroblastoma-specific.RData')
nb.circs<-CIRCS
rm(CIRCS, CIRCS.hg19, hb.circs, vt.circs, nb.genes, hb.genes, vt.genes, hb.meta, vt.meta)


#  load DE results for genes andsplit to up-/downregulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_MNA_HR_nMNA.RData')
up.gns<-subset(RES, log2FoldChange>0)
down.gns<-subset(RES, log2FoldChange<0)
rm(RES, CND, DDS, PCA, VE, VSC)


#  load the union of splice-factors and identify the up-/downregulated 
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/SpliceAid-F+GO.RData')
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/MSigDB_c2_v7.0_splicing.RData')
db<-unique(rbind(data.frame(db[, c('gene_name', 'gene_id')]), msigdbSF[, c('gene_name', 'gene_id')]))
up.sf<-up.gns[ rownames(up.gns) %in% db$gene_id, , drop=F]
down.sf<-down.gns[ rownames(down.gns) %in% db$gene_id, , drop=F]
rm(db, msigdbSF)


#  load circRNA sequence RBP enrichments
#  keep only cases with 3 or more significant motifs for the neuroblastoma-specific circRNAs
#  summarize counts at the gene level by taking the mean and median
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
rbp.exons<-fimo[ circ_name %in% nb.circs$circ_name & ncount>=3 ][, c('g', 'rest'):=tstrsplit(circ_name, '_', fixed=T)][, c('gene_id', 'gene_name'):=tstrsplit(g, '|', fixed=T)][, c('g', 'rest'):=NULL][, .(circ_names=list(circ_name), rbp_ids=list(unique(unlist(gene_ids))), mean=mean(as.numeric(unlist(ncount))), median=median(as.numeric(unlist(ncount)))), by=.(gene_id, gene_name)]   #  as.numeric() because it complains about integer type
rm(ame, fimo.control, fimo.control.all, fimo.all, fimo)


#  load circRNA flanking introns RBP enrichments
#  keep only flanking intron cases with 3 or more significant motifs for the neuroblastoma-specific circRNAs
#  summarize counts of both flanking introns and isoforms at the gene level by taking the mean and median
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results.RData')
rbp.introns<-fimo[ circ_name %in% nb.circs$circ_name & ncount>=3 ][, c('g', 'rest'):=tstrsplit(circ_name, '_', fixed=T)][, c('gene_id', 'gene_name'):=tstrsplit(g, '|', fixed=T)][, c('g', 'rest'):=NULL][, .(circ_names=list(circ_name), rbp_ids=list(unique(unlist(gene_ids))), mean=mean(as.numeric(unlist(ncount))), median=median(as.numeric(unlist(ncount)))), by=.(gene_id, gene_name)]   #  as.numeric() because it complains about integer type
rm(ame, fimo.control, fimo.control.all, fimo.all, fimo)


#  [motif counts of significantly upregulated splice factor on circRNA exons] Mann-Whitney U test
y.grp<-rownames(subset(up.sf, padj<0.05))                                               #  significantly upregulated SF
y.rest<-c(rownames(down.sf), rownames(up.sf[ ! rownames(up.sf) %in% y.grp, , drop=F]))  #  not significantly upregulated SF
wilcox.test(x=rbp.exons[ sapply(rbp_ids, function(r){ any(r %in% y.grp)}), unlist(mean)], 
            y=rbp.exons[ sapply(rbp_ids, function(r){ any(r %in% y.rest)}), unlist(mean)], alternative='greater')$p.value
rm(y.grp, y.rest)
#
#  => 0.5050236175


#  [motif counts of significantly downregulated splice factor on circRNA exons] Mann-Whitney U test
y.grp<-rownames(subset(down.sf, padj<0.05))                                               #  significantly downregulated SF
y.rest<-c(rownames(up.sf), rownames(down.sf[ ! rownames(down.sf) %in% y.grp, , drop=F]))  #  not significantly downregulated SF
wilcox.test(x=rbp.exons[ sapply(rbp_ids, function(r){ any(r %in% y.grp)}), unlist(mean)], 
            y=rbp.exons[ sapply(rbp_ids, function(r){ any(r %in% y.rest)}), unlist(mean)], alternative='greater')$p.value
rm(y.grp, y.rest)
#
#  => 0.2696635114

#}}}




###################################################
#
#
#  global intron/exon/transcript RBP motif analysis
#
#
###################################################




#  [run once] construct unique sequences of annotated exons and introns
#                   
#                   we remove chrM and chrY exons/introns
#                   we ask for at least 15nts long exons/introns
#                   we remove introns longer that 1e6 nts
#{{{
rm(list=ls())
library(data.table)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)


#  load reference
#  split to genes, exons, transcripts
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
mcols(hsa)<-mcols(hsa)[, c('type', 'gene_id', 'gene_name', 'gene_type', 'transcript_id', 'transcript_name', 'transcript_type', 'exon_id', 'exon_number')]
seqlevels(hsa, pruning.mode='coarse')<-setdiff(seqlevels(hsa), c('chrY', 'chrM'))
hsa.g<-hsa[ hsa$type %in% 'gene' ]
hsa.t<-hsa[ hsa$type %in% 'transcript' ]
hsa.e<-hsa[ hsa$type %in% 'exon' ]
hsa.e<-hsa.e[ width(hsa.e)>=15 & width(hsa.e)<=1e6 ]
rm(hsa)


#  extract all introns of all transcripts that have at least one
#  keep only those that are at least 15nts long
txdb<-loadDb('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/TxDb_GRCh38.gencode.v30.RData')
hsa.i<-intronsByTranscript(txdb, use.names=T)
seqlevels(hsa.i, pruning.mode='coarse')<-seqlevels(hsa.g)  #  if one hits the seqlevel, all will hit it
hsa.i<-unlist(hsa.i[lengths(hsa.i)>0])
hsa.i$transcript_id<-names(hsa.i)
names(hsa.i)<-NULL
hsa.i$gene_id<-hsa.t$gene_id[ match(hsa.i$transcript_id, hsa.t$transcript_id) ]
hsa.i$gene_name<-hsa.t$gene_name[ match(hsa.i$transcript_id, hsa.t$transcript_id) ]
hsa.i$gene_type<-hsa.t$gene_type[ match(hsa.i$transcript_id, hsa.t$transcript_id) ]
hsa.i$transcript_type<-hsa.t$transcript_type[ match(hsa.i$transcript_id, hsa.t$transcript_id) ]
hsa.i$transcript_name<-hsa.t$transcript_name[ match(hsa.i$transcript_id, hsa.t$transcript_id) ]
hsa.i<-hsa.i[ width(hsa.i)>=15 & width(hsa.i)<=1e6 ]


#  unique exons grouping the corresponding metadata together
x<-data.table(data.frame(hsa.e))[, c('width', 'pos'):=list(NULL, paste(seqnames, strand, start, end, sep='_'))][, .(gene_id=list(unique(gene_id)), gene_name=list(unique(gene_name)), gene_type=list(unique(gene_type)), transcript_id=list(unique(transcript_id)), transcript_name=list(unique(transcript_name)), transcript_type=list(unique(transcript_type)), exon_id=list(unique(exon_id))),  by=.(pos)][, c('seqnames', 'strand', 'start', 'end'):=tstrsplit(pos, '_', fixed=T)][, pos:=NULL]
EXONS<-GRanges(seqnames=x$seqnames, strand=x$strand, ranges=IRanges(start=as.integer(x$start), end=as.integer(x$end)), gene_id=List(x$gene_id), gene_name=List(x$gene_name), gene_type=List(x$gene_type), transcript_id=List(x$transcript_id), transcript_name=List(x$transcript_name), transcript_type=List(x$transcript_type), exon_id=List(x$exon_id))
rm(x)


#  unique introns grouping the corresponding metadata together
x<-data.table(data.frame(hsa.i))[, c('width', 'pos'):=list(NULL, paste(seqnames, strand, start, end, sep='_'))][, .(gene_id=list(unique(gene_id)), gene_name=list(unique(gene_name)), gene_type=list(unique(gene_type)), transcript_id=list(unique(transcript_id)), transcript_name=list(unique(transcript_name)), transcript_type=list(unique(transcript_type))),  by=.(pos)][, c('seqnames', 'strand', 'start', 'end'):=tstrsplit(pos, '_', fixed=T)][, pos:=NULL]
INTRONS<-GRanges(seqnames=x$seqnames, strand=x$strand, ranges=IRanges(start=as.integer(x$start), end=as.integer(x$end)), gene_id=List(x$gene_id), gene_name=List(x$gene_name), gene_type=List(x$gene_type), transcript_id=List(x$transcript_id), transcript_name=List(x$transcript_name), transcript_type=List(x$transcript_type))
rm(x, hsa.i, hsa.t, hsa.e, hsa.g, txdb)


#  create unique names based on genome positions
EXONS$exon_name<-paste0(as.character(seqnames(EXONS)),as.character(strand(EXONS)), start(EXONS), '-', end(EXONS))
INTRONS$intron_name<-paste0(as.character(seqnames(INTRONS)),as.character(strand(INTRONS)), start(INTRONS), '-', end(INTRONS))


#  get the sequences
EXONS.seqs<-setNames(getSeq(BSgenome.Hsapiens.UCSC.hg38, EXONS), EXONS$exon_name)
stopifnot( all.equal( width(EXONS), width(EXONS.seqs) ) )
INTRONS.seqs<-setNames(getSeq(BSgenome.Hsapiens.UCSC.hg38, INTRONS), INTRONS$intron_name)
stopifnot( all.equal( width(INTRONS), width(INTRONS.seqs) ) )


#  save
save(EXONS, EXONS.seqs, INTRONS, INTRONS.seqs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/unique_exons_introns.RData')
writeXStringSet(EXONS.seqs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/unique_exons.fa', format='fasta')
writeXStringSet(INTRONS.seqs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/unique_introns.fa', format='fasta')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/unique_exons_introns.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/unique_exons.fa
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/unique_introns.fa



#  [exons/introns] AME : ATtRACT human RBP motif enrichment for 6mers and above
#                 FIMO : count individual motif occurrences for 6mers and above 
#{{{
rm(list=ls())
library(data.table)
library(GenomicAlignments)
library(rtracklayer)
library(XML)


#  create the subdirectories
#  run in separate screen sessions the different analyses
#{{{

system2('mkdir', args='-p /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/ame/exons')
system2('mkdir', args='-p /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/ame/introns')
system2('mkdir', args='-p /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/fimo/exons')
system2('mkdir', args='-p /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/fimo/introns')
#
#  AME for exons:
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/ame/exons
#      fasta-get-markov -m 1 -rna -norc ../../../../unique_exons.fa background.model
#      ame --oc ./ --method fisher --scoring avg --evalue-report-threshold 0.05 --bfile background.model --control --shuffle-- --kmer 2 ../../../../unique_exons.fa /fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme
#
#  AME for introns:
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/ame/introns
#      fasta-get-markov -m 1 -rna -norc ../../../../unique_introns.fa background.model
#      ame --oc ./ --method fisher --scoring avg --evalue-report-threshold 0.05 --bfile background.model --control --shuffle-- --kmer 2 ../../../../unique_introns.fa /fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme
#
#  FIMO for exons (p-value threshold 1e-6 stricter than the circRNAs):
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/fimo/exons
#      fasta-get-markov -m 1 -rna -norc ../../../../unique_exons.fa background.model
#      fimo --oc ./ --norc --text --thresh 1e-6 --bfile background.model /data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme ../../../../unique_exons.fa > fimo.tsv
#
#  FIMO for introns (p-value threshold 1e-6 stricter than the circRNAs):
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/fimo/introns
#      fasta-get-markov -m 1 -rna -norc ../../../../unique_introns.fa background.model
#      fimo --oc ./ --norc --text --thresh 1e-6 --bfile background.model /data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme ../../../../unique_introns.fa > fimo.tsv

#}}}


#  load reference to add gene_id to the results
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]


#  [exons/introns] AME results
#                  convert p-values to numerics
#                  add FDR = FP/(FP+TP) column 
ame.exons<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/ame/exons/ame.tsv', header=T, sep='\t', strip.white=T, blank.lines.skip=F, drop=c('motif_DB', 'FASTA_max', 'PWM_min'))[, fdr:=round(FP/(FP+TP), 3)]  #  comments at the bottom of the file produce warning
colnames(ame.exons)<-c('rank', 'motif_id', 'motif_alt_id', 'consensus', 'pvalue', 'padj', 'evalue', 'tests', 'pos', 'neg', 'tp', '%tp', 'fp', '%fp', 'fdr')
ame.exons<-ame.exons[, c('pvalue', 'padj', 'evalue'):=list(as.numeric(pvalue), as.numeric(padj), as.numeric(evalue))]
ame.introns<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/ame/introns/ame.tsv', header=T, sep='\t', strip.white=T, blank.lines.skip=F, drop=c('motif_DB', 'FASTA_max', 'PWM_min'))[, fdr:=round(FP/(FP+TP), 3)]  #  comments at the bottom of the file produce warning
colnames(ame.introns)<-c('rank', 'motif_id', 'motif_alt_id', 'consensus', 'pvalue', 'padj', 'evalue', 'tests', 'pos', 'neg', 'tp', '%tp', 'fp', '%fp', 'fdr')
ame.introns<-ame.introns[, c('pvalue', 'padj', 'evalue'):=list(as.numeric(pvalue), as.numeric(padj), as.numeric(evalue))]


#  FIMO results
#{{{

#  FIMO cuts off motif_ids at 100 characters. We need to fix that.
load('/fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')
rm(pwm)

 
#  [exons/introns] load FIMO results (warning about last line skipping is normal, it contains comments about FIMO command run)
fimo.exons<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/fimo/exons/fimo.tsv', header=T, strip.white=T, blank.lines.skip=T, select=c('motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'score', 'p-value', 'q-value', 'matched_sequence'))
colnames(fimo.exons)<-c('motif_id', 'motif_alt_id', 'exon_name', 'start', 'end', 'score', 'pvalue', 'qvalue', 'seq')
fimo.introns<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/fimo/introns/fimo.tsv', header=T, strip.white=T, blank.lines.skip=T, select=c('motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'score', 'p-value', 'q-value', 'matched_sequence'))
colnames(fimo.introns)<-c('motif_id', 'motif_alt_id', 'intron_name', 'start', 'end', 'score', 'pvalue', 'qvalue', 'seq')


#  fix cut motif_ids
#{{{

#  identify the motif_ids 100 characters long, most likely they were cut, and the corresponding alternative motif_ids
#  construct the full-length motif_ids 
#  check that each cut motif_id matches one full length motif_id
#  
#  exons
m<-fimo.exons[nchar(motif_id)==100, unique(motif_id)]
a<-fimo.exons[ motif_id %in% m, unique(motif_alt_id)]
d<-db[ motif %in% a, paste(sapply(gene_name, paste, sep='', collapse=','), id, sep='_') ]
for(i in m){
    k<-grep(i, d, fixed=T)
    if (length(k)==0){
        stop(paste('cut motif:', i, 'does not have a home!'))
    } else if (length(k)>1){
        stop(paste('cut motif:', i, 'has many homes!'))
    }
    fimo.exons[ motif_id %in% i, motif_id:=d[k] ]
}
#
#  introns
#
m<-fimo.introns[nchar(motif_id)==100, unique(motif_id)]
a<-fimo.introns[ motif_id %in% m, unique(motif_alt_id)]
d<-db[ motif %in% a, paste(sapply(gene_name, paste, sep='', collapse=','), id, sep='_') ]
for(i in m){
    k<-grep(i, d, fixed=T)
    if (length(k)==0){
        stop(paste('cut motif:', i, 'does not have a home!'))
    } else if (length(k)>1){
        stop(paste('cut motif:', i, 'has many homes!'))
    }
    fimo.introns[ motif_id %in% i, motif_id:=d[k] ]
}

#}}}


#  identify the RBP groups for each motif
fimo.exons<-fimo.exons[, rbp_group:=lapply(strsplit(sub('_.*$', '', motif_id), ','), unique)]
fimo.introns<-fimo.introns[, rbp_group:=lapply(strsplit(sub('_.*$', '', motif_id), ','), unique)]


#  [~20mins] identify the overlapping motifs per sequence
fimo.exons[, o_group:={g<-seq_len(.N);
                 o<-mcols(reduce(IRanges(start=start, end=end), with.revmap=T))$revmap  #  God bless for reverse-map!!!!!
                 invisible(sapply(seq_along(o), function(n){ g[ o[[n]] ]<<-n }))        #  assign same group number to all overlapping motifs
                 g
                }, by=.(exon_name)]
fimo.introns[, o_group:={g<-seq_len(.N);
                 o<-mcols(reduce(IRanges(start=start, end=end), with.revmap=T))$revmap  #  God bless for reverse-map!!!!!
                 invisible(sapply(seq_along(o), function(n){ g[ o[[n]] ]<<-n }))        #  assign same group number to all overlapping motifs
                 g
                }, by=.(intron_name)]


#  group overlapping motifs per sequence
fimo.exons<-fimo.exons[, .(motif_id=list(unique(motif_id)), motif_alt_id=list(unique(motif_alt_id)), start=list(start), end=list(end), score=list(score), pvalue=list(pvalue), qvalue=list(qvalue), seq=list(seq), rbp_group=list(unique(unlist(rbp_group)))), by=.(exon_name, o_group)][, o_group:=NULL]
fimo.introns<-fimo.introns[, .(motif_id=list(unique(motif_id)), motif_alt_id=list(unique(motif_alt_id)), start=list(start), end=list(end), score=list(score), pvalue=list(pvalue), qvalue=list(qvalue), seq=list(seq), rbp_group=list(unique(unlist(rbp_group)))), by=.(intron_name, o_group)][, o_group:=NULL]


#  keep on the side all motifs grouped by overlaps
fimo.exons.all<-copy(fimo.exons)
fimo.introns.all<-copy(fimo.introns)


#  identify number of non-overlapping motifs per sequence collapsing motif_ids, motif_alt_ids, motif sequences and the RBPs involved
fimo.exons<-fimo.exons[, .(ncount=.N, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(exon_name)]
fimo.introns<-fimo.introns[, .(ncount=.N, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(intron_name)]


#  add gene_ids to the RBPs involved per sequence
fimo.exons[, gene_ids:=lapply(rbp_group, function(r){ hsa$gene_id[ match(r, hsa$gene_name) ] })]
fimo.introns[, gene_ids:=lapply(rbp_group, function(r){ hsa$gene_id[ match(r, hsa$gene_name) ] })]


#  include lengths of the exon/intron sequences 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/unique_exons_introns.RData')
l<-setNames(width(EXONS.seqs), names(EXONS.seqs))
fimo.exons[, length:=l[ exon_name ]]
fimo.exons.all[, length:=l[ exon_name ]]
l<-setNames(width(INTRONS.seqs), names(INTRONS.seqs))
fimo.introns[, length:=l[ intron_name ]]
fimo.introns.all[, length:=l[ intron_name ]]
rm(l)


#  compute count densities per 1K of sequence
fimo.exons[, density:=ncount/length*1e3]
fimo.introns[, density:=ncount/length*1e3]

#}}}


#  save
save(ame.exons, ame.introns, fimo.exons, fimo.introns, fimo.exons.all, fimo.introns.all, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/results_exons+introns.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/results_exons+introns.RData



#  [MNA vs HR_nMNA] are exons/introns of differentially spliced genes enriched in motifs of the differentially up-regulated RBPs
#                   as it was reported before (https://www.ncbi.nlm.nih.gov/pubmed/26683771)?
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(grid)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load ATtRACT human RBPs with 6mer and above motifs
load('/fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')
rbp<-db[ sapply(motif, function(m){ any(nchar(m)>=6) }), unique(unlist(gene_name))]
rm(db)


#  load RBP motifs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/results_exons+introns.RData')
fe<-fimo.exons
fi<-fimo.introns
rm(list=ls(pattern='fimo|ame'))


#  load DS results and DE (exons) results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/edgeR_exons_MNA_HR_nMNA.RData')
de.res<-DE.RES
ds.res<-DS.RES
rm(list=setdiff(ls(), c('rbp', 'fe', 'fi', 'ds.res', 'de.res')))


#  load DE results
#  isolate the significantly differentially up-/downregulated RBPs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_MNA_HR_nMNA.RData')
up.rbp<-subset(RES, log2FoldChange>0 & gene_name %in% rbp)
down.rbp<-subset(RES, log2FoldChange<0 & gene_name %in% rbp)
notde.rbp<-setdiff(rbp, c(up.rbp$gene_name, down.rbp$gene_name))  #  RBMY1A1 (not present altogether since it was not expressed)
rm(list=setdiff(ls(), c('fe', 'fi', 'ds.res', 'de.res', 'up.rbp', 'down.rbp')))


#  load exons/introns reference
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/unique_exons_introns.RData')
rm(INTRONS.seqs, EXONS.seqs)


#  isolate the exons/introns of differentially spliced genes and the rest
ds.exons<-EXONS[ sapply(EXONS$gene_id, function(i){ any( i %in% ds.res$gene_id )}) ]
notds.exons<-EXONS[ sapply(EXONS$gene_id, function(i){ !any( i %in% ds.res$gene_id )}) ]
ds.introns<-INTRONS[ sapply(INTRONS$gene_id, function(i){ any( i %in% ds.res$gene_id )}) ]
notds.introns<-INTRONS[ sapply(INTRONS$gene_id, function(i){ !any( i %in% ds.res$gene_id )}) ]


#  get the up-/downregulated RBP motif counts for exons/introns of differentially spliced genes and the rest
#  N.B. not upregulated in this case means downregulated and vice versa
ds.exons.up<-fe[ sapply(rbp_group, function(r){ any( r %in% up.rbp$gene_name )}) & exon_name %in% ds.exons$exon_name, ncount ]
notds.exons.up<-fe[ sapply(rbp_group, function(r){ any( r %in% up.rbp$gene_name )}) & exon_name %in% notds.exons$exon_name, ncount ]
ds.exons.down<-fe[ sapply(rbp_group, function(r){ any( r %in% down.rbp$gene_name )}) & exon_name %in% ds.exons$exon_name, ncount ]
notds.exons.down<-fe[ sapply(rbp_group, function(r){ any( r %in% down.rbp$gene_name )}) & exon_name %in% notds.exons$exon_name, ncount ]
ds.introns.up<-fi[ sapply(rbp_group, function(r){ any( r %in% up.rbp$gene_name )}) & intron_name %in% ds.introns$intron_name, ncount ]
notds.introns.up<-fi[ sapply(rbp_group, function(r){ any( r %in% up.rbp$gene_name )}) & intron_name %in% notds.introns$intron_name, ncount ]
ds.introns.down<-fi[ sapply(rbp_group, function(r){ any( r %in% down.rbp$gene_name )}) & intron_name %in% ds.introns$intron_name, ncount ]
notds.introns.down<-fi[ sapply(rbp_group, function(r){ any( r %in% down.rbp$gene_name )}) & intron_name %in% notds.introns$intron_name, ncount ]


#  one-sided Mann-Whitney test of RBP motif counts for exons/introns of differentially spliced genes and the rest
wilcox.test(x=ds.exons.up, y=notds.exons.up, alternative='greater')$p.value          #  8.871539593e-11
wilcox.test(x=ds.exons.down, y=notds.exons.down, alternative='greater')$p.value      #  8.53301384e-07
wilcox.test(x=ds.introns.up, y=notds.introns.up, alternative='greater')$p.value      #  1.52156208e-51
wilcox.test(x=ds.introns.down, y=notds.introns.down, alternative='greater')$p.value  #  5.824401163e-45

#}}}



#  [transcripts] AME : ATtRACT human RBP motif enrichment for 6mers and above
#               FIMO : count individual motif occurrences for 6mers and above 
#{{{
rm(list=ls())
library(data.table)
library(GenomicAlignments)
library(rtracklayer)
library(XML)


#  create the subdirectories
#  run in separate screen sessions the different analyses
#{{{

system2('mkdir', args='-p /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/ame/transcripts')
system2('mkdir', args='-p /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/fimo/transcripts')

#
#  AME for transcripts:
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/ame/transcripts
#      fasta-get-markov -m 1 -rna -norc /fast/groups/ag_schulte/work/reference/indices/GRCh38_STAR_2.7.1a/rsem.transcripts.fa background.model
#      ame --oc ./ --method fisher --scoring avg --evalue-report-threshold 0.05 --bfile background.model --control --shuffle-- --kmer 2 /fast/groups/ag_schulte/work/reference/indices/GRCh38_STAR_2.7.1a/rsem.transcripts.fa /fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme
#
#  FIMO for transcripts (p-value threshold 1e-5 one order of magnitude stricter than the circRNAs):
#
#      cd /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/fimo/transcripts
#      fasta-get-markov -m 1 -rna -norc /fast/groups/ag_schulte/work/reference/indices/GRCh38_STAR_2.7.1a/rsem.transcripts.fa background.model
#      fimo --oc ./ --norc --text --thresh 1e-5 --bfile background.model /data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme /fast/groups/ag_schulte/work/reference/indices/GRCh38_STAR_2.7.1a/rsem.transcripts.fa > fimo.tsv
#}}}


#  load reference to add gene_id to the RBPs and gene_id and gene_name to the transcripts 
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa.gns<-hsa[ hsa$type %in% 'gene' ]
hsa.txs<-hsa[ hsa$type %in% 'transcript' ]
mcols(hsa.gns)<-mcols(hsa.gns)[, c('gene_id', 'gene_type', 'gene_name')]
mcols(hsa.txs)<-mcols(hsa.txs)[, c('gene_id', 'gene_type', 'gene_name', 'transcript_id', 'transcript_type', 'transcript_name')]
rm(hsa)


#  AME results
#  convert p-values to numerics
#  add FDR = FP/(FP+TP) column 
ame.txs<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/ame/transcripts/ame.tsv', header=T, sep='\t', strip.white=T, blank.lines.skip=F, drop=c('motif_DB', 'FASTA_max', 'PWM_min'))[, fdr:=round(FP/(FP+TP), 3)]  #  comments at the bottom of the file produce warning
colnames(ame.txs)<-c('rank', 'motif_id', 'motif_alt_id', 'consensus', 'pvalue', 'padj', 'evalue', 'tests', 'pos', 'neg', 'tp', '%tp', 'fp', '%fp', 'fdr')
ame.txs<-ame.txs[, c('pvalue', 'padj', 'evalue'):=list(as.numeric(pvalue), as.numeric(padj), as.numeric(evalue))]


#  FIMO results
#{{{

#  FIMO cuts off motif_ids at 100 characters. We need to fix that.
load('/fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')
rm(pwm)

 
#  load FIMO results (warning about last line skipping is normal, it contains comments about FIMO command run)
fimo.txs<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/fimo/transcripts/fimo.tsv', header=T, strip.white=T, blank.lines.skip=T, select=c('motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'score', 'p-value', 'q-value', 'matched_sequence'))
colnames(fimo.txs)<-c('motif_id', 'motif_alt_id', 'transcript_id', 'start', 'end', 'score', 'pvalue', 'qvalue', 'seq')


#  fix cut motif_ids
#  identify the motif_ids 100 characters long, most likely they were cut, and the corresponding alternative motif_ids
#  construct the full-length motif_ids 
#  check that each cut motif_id matches one full length motif_id
m<-fimo.txs[nchar(motif_id)==100, unique(motif_id)]
a<-fimo.txs[ motif_id %in% m, unique(motif_alt_id)]
d<-db[ motif %in% a, paste(sapply(gene_name, paste, sep='', collapse=','), id, sep='_') ]
for(i in m){
    k<-grep(i, d, fixed=T)
    if (length(k)==0){
        stop(paste('cut motif:', i, 'does not have a home!'))
    } else if (length(k)>1){
        stop(paste('cut motif:', i, 'has many homes!'))
    }
    fimo.txs[ motif_id %in% i, motif_id:=d[k] ]
}


#  identify the RBP groups for each motif
fimo.txs<-fimo.txs[, rbp_group:=lapply(strsplit(sub('_.*$', '', motif_id), ','), unique)]


#  [~20mins] identify the overlapping motifs per sequence
fimo.txs[, o_group:={g<-seq_len(.N);
                 o<-mcols(reduce(IRanges(start=start, end=end), with.revmap=T))$revmap  #  God bless for reverse-map!!!!!
                 invisible(sapply(seq_along(o), function(n){ g[ o[[n]] ]<<-n }))        #  assign same group number to all overlapping motifs
                 g
                }, by=.(transcript_id)]


#  group overlapping motifs per sequence
fimo.txs<-fimo.txs[, .(motif_id=list(unique(motif_id)), motif_alt_id=list(unique(motif_alt_id)), start=list(start), end=list(end), score=list(score), pvalue=list(pvalue), qvalue=list(qvalue), seq=list(seq), rbp_group=list(unique(unlist(rbp_group)))), by=.(transcript_id, o_group)][, o_group:=NULL]


#  keep on the side all motifs grouped by overlaps
fimo.txs.all<-copy(fimo.txs)


#  identify number of non-overlapping motifs per sequence collapsing motif_ids, motif_alt_ids, motif sequences and the RBPs involved
fimo.txs<-fimo.txs[, .(ncount=.N, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(transcript_id)]


#  [~10mins] add gene_ids to the RBPs involved per sequence
fimo.txs[, rbp_ids:=lapply(rbp_group, function(r){ hsa.gns$gene_id[ match(r, hsa.gns$gene_name) ] })]


#  add gene_name, gene_id to the transcripts
fimo.txs.all[, c('gene_id', 'gene_name'):=list(hsa.txs$gene_id[ match(transcript_id, hsa.txs$transcript_id) ], hsa.txs$gene_name[ match(transcript_id, hsa.txs$transcript_id) ])]
fimo.txs[, c('gene_id', 'gene_name'):=list(hsa.txs$gene_id[ match(transcript_id, hsa.txs$transcript_id) ], hsa.txs$gene_name[ match(transcript_id, hsa.txs$transcript_id) ])]


#  include lengths of the transcript sequences 
txs.seqs<-readDNAStringSet('/fast/groups/ag_schulte/work/reference/indices/GRCh38_STAR_2.7.1a/rsem.transcripts.fa')
l<-setNames(width(txs.seqs), names(txs.seqs))
fimo.txs[, length:=l[ transcript_id ]]
fimo.txs.all[, length:=l[ transcript_id ]]
rm(l, txs.seqs)


#  compute count densities per 1K of sequence
fimo.txs[, density:=ncount/length*1e3]


#  reorder columns
colnames(fimo.txs)
setcolorder(fimo.txs, c('transcript_id', 'gene_id', 'gene_name', 'length', 'motif_id', 'motif_alt_id', 'seq', 'rbp_group', 'rbp_ids', 'ncount', 'density'))
colnames(fimo.txs.all)
setcolorder(fimo.txs.all, c('transcript_id', 'gene_id', 'gene_name', 'length', 'motif_id', 'motif_alt_id', 'start', 'end', 'score', 'pvalue', 'qvalue', 'seq', 'rbp_group'))

#}}}


#  save
save(ame.txs, fimo.txs, fimo.txs.all, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/results_transcripts.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/results_transcripts.RData



#  [KHSRP] distribution of mean motif counts across transcripts per gene
#          pick out interesting targets that are well-expressed in tumors
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(KEGGREST)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load gene members of different KEGG signaling pathways lifted over to GENCODE v30
load('/fast/groups/ag_schulte/work/Bock_CRISPR/202007/KEGG_signaling_pathways.RData')


#  load transcriptome counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/results_transcripts.RData')


#  load gene expression (TPMs) and keep only tumors
#  compute mean TPM across tumors
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')
meta<-meta[ !is.na(risk_group) ]
stopifnot( length(setdiff(meta[, bid], colnames(gns.tpm)))==0 )
gns.tpm<-gns.tpm[, meta[, bid]]
gns.tpm<-rowMeans(gns.tpm)
rm(list=setdiff(ls(), c(l, 'meta', 'hsa', 'gns.tpm')))


#  KHSRP
#{{{

PICK<-'KHSRP'
CUT<-5
STEP<-1
TPM_CUTOFF<-10


#  get the counts on all transcripts ONLY for KHSRP-related motifs even if they are shared among other RBPs (you need to use the fimo.txs.all object)
#  summarize them at the transcript level fist and then take the mean of counts and densities at the gene level
f<-fimo.txs.all[ sapply(rbp_group, function(r){ any(grepl(PICK, r)) }) ][, .(ncount=.N, density=.N/length*1e3, length=length, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(transcript_id, gene_id, gene_name)][, .(ncount=mean(ncount), density=mean(density), length=mean(length), rbp_group=list(unique(unlist(rbp_group))), motif_alt_id=list(unique(unlist(motif_alt_id)))), by=.(gene_id, gene_name)]


#  remove targets with mean TPM in tumors below the cutoff
#  order by mean counts
f<-f[ gene_name %in% names(gns.tpm[ gns.tpm>=TPM_CUTOFF ]) ][order(-ncount)]
f[1:10, c('gene_id', 'gene_name', 'ncount', 'density', 'length')]
#                gene_id  gene_name      ncount      density      length
#  1: ENSG00000181722.16     ZBTB20 8.000000000 0.2949308756 27125.00000
#  2: ENSG00000203930.11  LINC00632 7.222222222 0.3670722533 19249.22222
#  3: ENSG00000184226.15      PCDH9 7.000000000 0.2479895136 28227.00000
#  4:  ENSG00000259976.3 AC093010.3 7.000000000 0.4601025371 15214.00000
#  5: ENSG00000166233.15      ARIH1 6.250000000 0.3132591124 19477.62500
#  6: ENSG00000173611.17       SCAI 6.250000000 0.7360161295 10637.37500
#  7: ENSG00000139746.15      RBM26 6.000000000 1.6384489350  3662.00000
#  8: ENSG00000160131.13      VMA21 6.000000000 1.2719949120  4717.00000
#  9: ENSG00000170500.12     LONRF2 6.000000000 0.4312823462 13912.00000
# 10: ENSG00000061987.16       MON2 5.166666667 0.4584266480 11177.33333


#  which genes are members of the PI3K-AKT pathway and have at least 2 motifs?
f[ gene_name %in% KEGG[ pathways %in% 'pi3k_akt', unique(unlist(gene_names))] & ncount>1 ][order(-ncount)][, c('gene_id', 'gene_name', 'ncount', 'density', 'length', 'rbp_group')]
#                gene_id gene_name      ncount      density       length                                   rbp_group
#  1: ENSG00000135679.24      MDM2 4.000000000 0.3166561115 12632.000000 ELAVL1,ELAVL4,ZFP36,ZFP36L2,ELAVL2,TIA1,...
#  2: ENSG00000162409.11    PRKAA2 4.000000000 0.4293066588  9317.500000       ELAVL1,ELAVL4,ELAVL2,TIA1,TIAL1,KHSRP
#  3: ENSG00000105855.10     ITGB8 3.181818182 0.5773270307  5886.545455  ELAVL1,ELAVL4,ELAVL3,KHSRP,ELAVL2,TIA1,...
#  4: ENSG00000091409.15     ITGA6 3.000000000 0.5308047709  5652.000000  ELAVL1,ELAVL4,ELAVL2,TIA1,TIAL1,ELAVL3,...
#  5:  ENSG00000105810.9      CDK6 3.000000000 0.2583534275 11612.000000  ELAVL1,ELAVL4,ELAVL2,TIA1,TIAL1,ELAVL3,...
#  6: ENSG00000165699.14      TSC1 2.714285714 1.9714109252  1457.857143                                       KHSRP
#  7: ENSG00000115414.19       FN1 2.000000000 0.2180787264  9171.000000   ELAVL2,ELAVL4,TIA1,TIAL1,ELAVL3,KHSRP,...
#  8: ENSG00000170248.15   PDCD6IP 2.000000000 2.3752969121   842.000000       ELAVL2,ELAVL4,TIA1,TIAL1,ELAVL3,KHSRP
#  9: ENSG00000137801.10     THBS1 2.000000000 0.2572347267  7775.000000                         ELAVL3,ELAVL4,KHSRP
# 10:  ENSG00000135862.6     LAMC1 2.000000000 0.2522386177  7929.000000                                       KHSRP
# 11: ENSG00000187391.21     MAGI2 1.400000000 0.2904011644  4939.000000 ELAVL1,ELAVL4,ZFP36,ZFP36L2,ELAVL2,TIA1,...


#  histogram (full range)
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.5, 8.0, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
XLIM<-pretty(range(f[, ncount]), 4)
h<-hist(f[,ncount], breaks=seq(min(XLIM), max(XLIM), length.out=50), plot=F)
h$counts<-log10(1+h$counts)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
plot(h, col='blue4', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(XLIM), add=F)
mtext('Number of motifs', side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext(expression(log[10]('Frequency')), side=2, line=5, padj=+0.2, cex=2.4, las=0) 
mtext(PICK, side=3, line=0, padj=+0.2, cex=1.8, las=0) 
dev.off()


#  barplot (pooling counts above a certain value)
x11(width=20, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(8.5, 12.0, 2.0, 4.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=4.0, cex.axis=4.0)
x<-table(cut(f[, ncount], breaks=c(seq(0, CUT, STEP), f[, max(ncount)]), labels=c(seq(1, CUT, STEP), paste0('>', CUT)), right=T))
x<-log10(setNames(as.integer(x), names(x)))
YTICK<-pretty(c(0, max(x)), 5)
YMAX<-tail(YTICK, 1)
br<-barplot(x, axes=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(br)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
br<-barplot(x, col='grey39', border='white', ylim=c(0, YMAX), axisnames=F, xlab='', ylab='', yaxt='n', xaxt='n', axes=F, add=T)
#axis(1, at=br[ceiling(seq(1, length(x), length.out=8))], labels=names(x[ceiling(seq(1, length(x), length.out=8))]), line=0, tick=F, padj=+0.5, cex.axis=4.0, xpd=NA)
axis(1, at=br, labels=names(x), line=0, tick=F, padj=+0.5, cex.axis=4.0, xpd=NA)
axis(2, at=YTICK, labels=YTICK, cex.axis=4.0)
mtext(expression(log[10]('Frequency')), side=2, line=8, padj=+0.2, las=0, cex=4.0)
mtext('Number of motifs', side=1, line=6, padj=-0.1, las=0, cex=4.0)
mtext(PICK, side=3, line=0, padj=+0.2, cex=3.0, las=0) 
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_', PICK, '_mean_number_of_motifs_across_transcripts_per_gene.svg'), width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}

#}}}




######################
#
# 
# exploratory analysis
#
# 
######################




#  [tumors] Correlation of RBP (6mer motifs or longer) including SF expression with circRNA expression (based on CIRI2 counts).
#
#           N.B. Only ~2.3% of circRNAs originate from RBPs, there is no need to remove them from the correlation analysis.
#
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
library(grid)
library(ggplot2)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
options(scipen=0)



#  [run once] Spearman correlations between the union of RBPs and SFs and circRNA counts (based on CIRI2)
#             cluster-specific enrichment of RBP motifs in circRNAs and flanking introns 
#{{{

#  union of SF
l<-ls()
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/SpliceAid-F+GO.RData')
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/MSigDB_c2_v7.0_splicing.RData')
SF<-data.table(unique(rbind(data.frame(db[, c('gene_name', 'gene_id')]), msigdbSF[, c('gene_name', 'gene_id')])))
rm(list=setdiff(ls(), c(l, 'SF')))


#  load all RBPs 
#  inflate per motif the subgroups based on PWMs into one group of unique genes (motif-inflating needs to happen alone, different 1:many map than genes)
#  and keep those with 6mer or longer motifs
#  keep 6mers or longer
l<-ls()
load('/fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')
RBP<-unique(db[, .(gene_name=unique(unlist(gene_name)), gene_id=unique(unlist(gene_id)), motif=motif), by=.(id)][, .(motif=unlist(motif)), by=.(gene_name, gene_id)][ sapply(motif, nchar)>=6, ][, motif:=NULL])
rm(list=setdiff(ls(), c(l, 'RBP')))


#  join RBPs and SFs into one RBP list
length(setdiff(SF$gene_name, RBP$gene_name))  #  178
ALL<-unique(rbind(RBP, SF))
rm(RBP)


#  load prepackaged gene counts and metadata
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')
meta<-meta[ !(failed) & !is.na(risk_group) ]
ALL.cnt<-gns.cnt[ rownames(gns.cnt) %in% ALL$gene_id, colnames(gns.cnt) %in% meta$bid ]
#
ALL[ gene_id %in% setdiff(ALL$gene_id, rownames(ALL.cnt)) ]  #  one RBP not expressed
# 
#  RBMY1A1 ENSG00000234414.7
#
rownames(ALL.cnt)<-ALL[ match(rownames(ALL.cnt), gene_id), gene_name ]
rm(list=setdiff(ls(), c(l, 'meta', 'ALL.cnt')))


#  load CIRI2 counts and convert to DGE
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
CIR.cnt<-data.table(data.frame(mcols(nb.circs)[, c('circ_name', 'bid', 'jc_count')]))
CIR.cnt<-dcast(CIR.cnt, bid ~ circ_name, value.var='jc_count', fun.aggregate=mean, fill=0.0)
CIR.cnt<-t(as.matrix(data.frame(CIR.cnt[, -1], row.names=CIR.cnt[, bid], check.names=F)))
rm(list=setdiff(ls(), c(l, 'meta', 'CIR.cnt')))


#  fix the sample order 
stopifnot( length(setdiff( colnames(ALL.cnt), colnames(CIR.cnt) ))==0 )
CIR.cnt<-CIR.cnt[, colnames(ALL.cnt)]
stopifnot( all.equal( colnames(CIR.cnt), colnames(ALL.cnt) ) )


#  [RBPs+SFs] compute Spearman correlations and cluster based on Euclidean distance of 1 - Spearman correlations
CO.ALL<-cor(t(CIR.cnt), t(ALL.cnt), use='pairwise.complete.obs', method='spearman')
H.ALL.CIR<-hclust(dist(1-CO.ALL, method='euclidean'), method='ward.D2')
H.ALL<-hclust(dist(1-t(CO.ALL), method='euclidean'), method='ward.D2')


#  [RBPs-SFs] compute Spearman correlations and cluster based on Euclidean distance of 1 - Spearman correlations
CO.RBP<-cor(t(CIR.cnt), t(ALL.cnt[ !rownames(ALL.cnt) %in% SF$gene_name, ]), use='pairwise.complete.obs', method='spearman')
H.RBP.CIR<-hclust(dist(1-CO.RBP, method='euclidean'), method='ward.D2')
H.RBP<-hclust(dist(1-t(CO.RBP), method='euclidean'), method='ward.D2')


#  [SFs] compute Spearman correlations and cluster based on Euclidean distance of 1 - Spearman correlations
CO.SF<-cor(t(CIR.cnt), t(ALL.cnt[ SF$gene_name, ]), use='pairwise.complete.obs', method='spearman')
H.SF.CIR<-hclust(dist(1-CO.SF, method='euclidean'), method='ward.D2')
H.SF<-hclust(dist(1-t(CO.SF), method='euclidean'), method='ward.D2')


#  strong correlation cutoff
COR_CUTOFF<-0.40


#  [RBPs+SFs] identify strongly correlated pairs 
#             identify the self-correlations among them
#             identify all RBP-related correlations (including self)
#             cluster strong correlating pairs based on Euclidean distance of 1 - Spearman correlations
ALL.strong<-apply(CO.ALL, 2, function(x){ setNames(x[x>=COR_CUTOFF], names(x[x>=COR_CUTOFF])) })
ALL.strong<-ALL.strong[ lengths(ALL.strong)>0 ]
ALL.strong.self<-sapply(names(ALL.strong), function(n){ ALL.strong[[n]][grepl(n, names(ALL.strong[[n]]))] })
ALL.strong.self<-ALL.strong.self[ lengths(ALL.strong.self)>0 ]
a<-paste(ALL$gene_name, collapse='|')
ALL.strong.all<-sapply(ALL.strong, function(n){ n[grepl(a, names(n))] })
ALL.strong.all<-ALL.strong.all[lengths(ALL.strong.all)>0]
a<-unique(unname(unlist(lapply(ALL.strong, names))))
b<-names(ALL.strong)
CO.ALL.strong<-CO.ALL[ a, b]
H.ALL.strong.CIR<-hclust(dist(1-CO.ALL.strong, method='euclidean'), method='ward.D2')
H.ALL.strong<-hclust(dist(1-t(CO.ALL.strong), method='euclidean'), method='ward.D2')
rm(a,b)


#  [RBPs-SFs] identify strongly correlated pairs 
#             identify the self-correlations among them
#             identify all RBP-related correlations (including self)
#             cluster strong correlating pairs based on Euclidean distance of 1 - Spearman correlations
RBP.strong<-apply(CO.RBP, 2, function(x){ setNames(x[x>=COR_CUTOFF], names(x[x>=COR_CUTOFF])) })
RBP.strong<-RBP.strong[ lengths(RBP.strong)>0 ]
RBP.strong.self<-sapply(names(RBP.strong), function(n){ RBP.strong[[n]][grepl(n, names(RBP.strong[[n]]))] })
RBP.strong.self<-RBP.strong.self[ lengths(RBP.strong.self)>0 ]
a<-paste(ALL$gene_name, collapse='|')
RBP.strong.all<-sapply(RBP.strong, function(n){ n[grepl(a, names(n))] })
RBP.strong.all<-RBP.strong.all[lengths(RBP.strong.all)>0]
a<-unique(unname(unlist(lapply(RBP.strong, names))))
b<-names(RBP.strong)
CO.RBP.strong<-CO.RBP[ a, b]
H.RBP.strong.CIR<-hclust(dist(1-CO.RBP.strong, method='euclidean'), method='ward.D2')
H.RBP.strong<-hclust(dist(1-t(CO.RBP.strong), method='euclidean'), method='ward.D2')
rm(a,b)


#  [SFs] identify strongly correlated pairs 
#        identify the self-correlations among them
#        identify all SF-related correlations (including self)
#        cluster strong correlating pairs based on Euclidean distance of 1 - Spearman correlations
SF.strong<-apply(CO.SF, 2, function(x){ setNames(x[x>=COR_CUTOFF], names(x[x>=COR_CUTOFF])) })
SF.strong<-SF.strong[ lengths(SF.strong)>0 ]
SF.strong.self<-sapply(names(SF.strong), function(n){ SF.strong[[n]][grepl(n, names(SF.strong[[n]]))] })
SF.strong.self<-SF.strong.self[ lengths(SF.strong.self)>0 ]
a<-paste(ALL$gene_name, collapse='|')
SF.strong.all<-sapply(SF.strong, function(n){ n[grepl(a, names(n))] })
SF.strong.all<-SF.strong.all[lengths(SF.strong.all)>0]
a<-unique(unname(unlist(lapply(SF.strong, names))))
b<-names(SF.strong)
CO.SF.strong<-CO.SF[ a, b]
H.SF.strong.CIR<-hclust(dist(1-CO.SF.strong, method='euclidean'), method='ward.D2')
H.SF.strong<-hclust(dist(1-t(CO.SF.strong), method='euclidean'), method='ward.D2')
rm(a,b)


#  cut at this many clusters both RBPs and circRNAs
NCLUST<-2
w.x<-cutree(H.ALL.CIR, k=NCLUST)
w.y<-cutree(H.ALL, k=NCLUST)
table(w.x)
#
#    1    2 
# 3995 1208 
# 
table(w.y)
#
#   1   2 
# 113 202 



#  which clusters do circARID1A, KHSRP and DHX9 belong to?
w.x[ grep('ARID1A', names(w.x)) ]
#
# ENSG00000117713.20|ARID1A_chr1+26729651-26732792 
#                                                2 
w.y[ grep('KHSRP|DHX9', names(w.y)) ]
#
#  DHX9 KHSRP 
#     1     2 



#  summary of correlations for KHSRP and DHX9
summary(CO.ALL[, 'KHSRP'])
#
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.3686 -0.0031  0.0875  0.0932  0.1867  0.6519 
#
summary(CO.ALL[, 'DHX9'])
#
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  -0.216   0.138   0.216   0.216   0.295   0.584 
#
table(CO.ALL[, 'DHX9']<0)
# 
# FALSE  TRUE 
#  5016   187 



#  load circRNA RBP motif enrichment results
#  isolate the circRNAs and RBPs of interest
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
fimo.ex<-fimo.all[ circ_name %in% rownames(CO.ALL) & sapply(rbp_group, function(r){ any( r %in% colnames(CO.ALL) ) }) ][, c('motif_id', 'motif_alt_id', 'start', 'end', 'score', 'pvalue', 'qvalue'):=NULL] 
rm(list=setdiff(ls(), c(l, 'fimo.ex')))



#  load circRNA intron RBP motif enrichment results
#  isolate the circRNAs and RBPs of interest
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results.RData')
fimo.in<-fimo.all[ circ_name %in% rownames(CO.ALL) & sapply(rbp_group, function(r){ any( r %in% colnames(CO.ALL) ) }) ][, c('motif_id', 'motif_alt_id', 'start', 'end', 'score', 'pvalue', 'qvalue'):=NULL] 
rm(list=setdiff(ls(), c(l, 'fimo.in')))



#  count number of motifs per circRNAs in a cluster-specific way
#
#  N.B. The number of motifs is biased by the number of RBPs and circRNAs in each cluster. We need normalized numbers to compare.
#
ex.m<-in.m<-data.frame(circCluster=integer(0), RBPluster=integer(0), counts=integer(0), ncounts=double(0))
for (i in 1:NCLUST){      #  circRNA clusters
    N<-sum(w.x==i)
    for (j in 1:NCLUST){  #  RBP clusters
        M<-N*sum(w.y==j)
        x<-fimo.ex[ circ_name %in% names( w.x[ w.x==i ]) & sapply(rbp_group, function(r){ any( r %in% names( w.y[ w.y==j ] )) }), ][, .(counts=.N), by=.(circ_name)][, counts]
        ex.m<-rbind(ex.m, data.frame(circCluster=i, RBPCluster=j, counts=x, ncounts=x/M))
        x<-fimo.in[ circ_name %in% names( w.x[ w.x==i ]) & sapply(rbp_group, function(r){ any( r %in% names( w.y[ w.y==j ] )) }), ][, .(counts=.N), by=.(circ_name)][, counts]
        in.m<-rbind(in.m, data.frame(circCluster=i, RBPCluster=j, counts=x, ncounts=x/M))
    }
}
ex.m<-data.table(ex.m)
in.m<-data.table(in.m)
rm(i, N, j, M, x)


#  save
save(list=ls(), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_cor_circRNAs.RData')

#}}}



#  load back
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_cor_circRNAs.RData')



#  ggplot2 theme
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



#  functions
#{{{
wt<-function(circs=c(1, 2), rbps=c(1,2), m, alt='greater'){
    stopifnot( length(rbps)==2 | length(circs)==2 )
    if (length(rbps)==2 & length(circs)==1){
        wt<-wilcox.test(x=m[ circCluster %in% circs & RBPCluster %in% rbps[1], ncounts], y=m[ circCluster %in% circs & RBPCluster %in% rbps[2], ncounts], alternative=alt)$p.value
    } else if (length(circs)==2 & length(rbps)==1){
        wt<-wilcox.test(x=m[ circCluster %in% circs[1] & RBPCluster %in% rbps, ncounts], y=m[ circCluster %in% circs[2] & RBPCluster %in% rbps, ncounts], alternative=alt)$p.value
    } else if (length(circs)==length(rbps) & length(circs)==2){
        wt<-wilcox.test(x=m[ circCluster %in% circs[1] & RBPCluster %in% rbps[1], ncounts], y=m[ circCluster %in% circs[2] & RBPCluster %in% rbps[2], ncounts], alternative=alt)$p.value
    }
    return(wt)
}
#}}}



#  check: do the RBPs with which the circARID1A is highly correlated have significant motifs on it or on its introns?
#{{{

#  load circARID1A FIMO results 
#  N.B. RBPs have p-values but no filtering was done
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_RBPs_across_tissues.RData')
rbp.arid1a.ex<-rbp
rm(list=setdiff(ls(), c(l, 'rbp.arid1a.ex')))


#  load circARID1A intron FIMO results 
#  N.B. RBPs have p-values but no filtering was done
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_introns_RBPs_across_tissues.RData')
rbp.arid1a.in<-rbp
rm(list=setdiff(ls(), c(l, 'rbp.arid1a.in')))


#  highly correlated RBPs with circARID1A
x<-CO.ALL[ grepl('ARID1A', rownames(CO.ALL)), ]
x<-x[ x>=COR_CUTOFF ]


#  any of the RBPs of interest having significant motifs on circARID1A or its introns?
rbp.arid1a.ex[ gene_name %in% names(x) ]
#    gene_name            gene_id ncount     pv.vt     pv.hb    pvalue
# 1:     KHSRP ENSG00000088247.17      2 0.0005853 2.386e-14 0.0005853
# 2:     PCBP1  ENSG00000169564.6      5 0.0023819 1.386e-16 0.0023819
#
rbp.arid1a.in[ gene_name %in% names(x) ]
#     gene_name            gene_id ncount     pv.vt     pv.hb    pvalue
#  1:     NOVA2  ENSG00000104967.7     32 5.320e-46 1.419e-21 1.419e-21
#  2:     CELF5 ENSG00000161082.13     12 1.642e-39 4.897e-21 4.897e-21
#  3:     SRSF4 ENSG00000116350.17     17 1.081e-13 7.952e-24 1.081e-13
#  4:      RBM5 ENSG00000003756.16     76 3.151e-26 1.468e-11 1.468e-11
#  5:    TARDBP ENSG00000120948.17     32 1.781e-25 1.278e-10 1.278e-10
#  6:      PUM1 ENSG00000134644.15     53 1.928e-30 1.874e-10 1.874e-10
#  7:     CELF6 ENSG00000140488.16     18 6.449e-43 6.604e-08 6.604e-08
#  8:   ZFP36L2  ENSG00000152518.8      3 1.836e-07 8.254e-20 1.836e-07
#  9:      AGO1 ENSG00000092847.12     31 2.965e-15 2.888e-04 2.888e-04
# 10:     KHSRP ENSG00000088247.17     96 5.894e-04 2.483e-14 5.894e-04
# 11:     PCBP1  ENSG00000169564.6    142 2.423e-03 1.508e-16 2.423e-03
# 12:      FXR2 ENSG00000129245.11      9 1.607e-22 1.471e-02 1.471e-02
# 13:     CELF4 ENSG00000101489.20     12 1.586e-44 6.001e-02 6.001e-02
# 14:     ENOX1 ENSG00000120658.13      1 7.976e-47 1.111e-01 1.111e-01
# 15:     IFIH1  ENSG00000115267.7     34 8.836e-35 1.302e-01 1.302e-01
# 16:   SNRNP70 ENSG00000104852.15      1 6.448e-01 1.233e-02 6.571e-01
# 17:     RBM42 ENSG00000126254.12      2 4.453e-35 7.378e-01 7.378e-01
# 18:     SRP14 ENSG00000140319.10     26 8.783e-01 3.973e-07 8.783e-01

#}}}



#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
B.col<-list(RBPCluster=setNames(c('brown1', 'darkorange4'), 1:NCLUST), circCluster=setNames(c('darkgoldenrod', 'darkolivegreen'), 1:NCLUST))



#  histograms of correlations
#{{{

#  [ALL] histogram of all correlations
par(mar=c(5.0, 6.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#b<-pretty(c(min(CO.ALL), max(CO.ALL)), 10)
b<-pretty(c(-0.8, 1.0), 10)
h<-hist(CO.ALL, breaks=seq(min(b), max(b), length.out=50), plot=F)
h$counts<-log10(1+h$counts)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
h<-plot(h, col='darkgrey', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(b), xaxt='n', add=F)
axis(1, at=seq(min(b), max(b), by=0.4), line=0, cex.axis=2.4)
#abline(v=COR_CUTOFF, col='red4', lty=2, lwd=4, xpd=F)
mtext('Spearman correlation', side=1, line=4, padj=-0.5, las=0, cex=2.4)
mtext(expression(log[10](1+'Frequency')), side=2, line=2, padj=-0.4, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_cor_ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_circRNAs.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b,h,YMAX)



#  [SF] histogram of all correlations
par(mar=c(5.0, 6.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#b<-pretty(c(min(CO.SF), max(CO.SF)), 10)
b<-pretty(c(-0.8, 1.0), 10)
h<-hist(CO.SF, breaks=seq(min(b), max(b), length.out=50), plot=F)
h$counts<-log10(1+h$counts)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
h<-plot(h, col='darkgrey', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(b), xaxt='n', add=F)
axis(1, at=seq(min(b), max(b), by=0.4), line=0, cex.axis=2.4)
#abline(v=COR_CUTOFF, col='red4', lty=2, lwd=4, xpd=F)
mtext('Spearman correlation', side=1, line=4, padj=-0.5, las=0, cex=2.4)
mtext(expression(log[10](1+'Frequency')), side=2, line=2, padj=-0.4, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_cor_SpliceAid_SF+MSigDB_SF_circRNAs.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b,h,YMAX)

#}}}



#  heatmap 
#{{{

grid.newpage()
an.y<-data.frame('RBP clusters'=factor(cutree(H.ALL, k=NCLUST), levels=1:NCLUST), check.names=F) 
an.x<-data.frame('circRNA clusters'=factor(cutree(H.ALL.CIR, k=NCLUST), levels=1:NCLUST), check.names=F)
an.col<-list('RBP clusters'=setNames(c('brown1', 'darkorange4'), 1:NCLUST), 'circRNA clusters'=setNames(c('darkgoldenrod', 'darkolivegreen'), 1:NCLUST))
BREAKS<-c(-0.6, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8)
BREAKS.COL<-rev(colorRampPalette(brewer.pal(5,'RdBu'))((length(BREAKS)-1)*2))
BREAKS.COL<-BREAKS.COL[c(seq(1, by=2, length.out=length(BREAKS)-2), length(BREAKS.COL))]  #  maximum contrast for the last bin
lg<-packLegend(
    Legend(title='', title_gp=gpar(fontsize=18, fontface='bold'), labels_gp=gpar(fontsize=16), grid_height=unit(4,'mm'),
        col_fun=circlize::colorRamp2(BREAKS[-1], BREAKS.COL), at=c(min(BREAKS), 0.2, max(BREAKS)), ncol=1), 
    Legend(title=colnames(an.y), title_gp=gpar(fontsize=18, fontface='bold'), labels_gp=gpar(fontsize=16), grid_height=unit(4,'mm'),
        at=names(an.col[[colnames(an.y)]]), legend_gp=gpar(fill=an.col[[colnames(an.y)]]), ncol=1), 
    Legend(title=colnames(an.x), title_gp=gpar(fontsize=18, fontface='bold'), labels_gp=gpar(fontsize=16), grid_height=unit(4,'mm'),
        at=names(an.col[[colnames(an.x)]]), legend_gp=gpar(fill=an.col[[colnames(an.x)]]), ncol=1), 
    column_gap=unit(1, 'cm'), max_height=unit(1, 'npc'))
ph<-ComplexHeatmap::pheatmap(as.matrix(CO.ALL), color=BREAKS.COL, border_color=NA, scale='none', 
        breaks=BREAKS,
        use_raster=F, 
        cluster_rows=H.ALL.CIR,
        cluster_cols=H.ALL,
        cutree_row=NCLUST,
        cutree_col=NCLUST,
        legend=F,
        annotation_legend=F,
        annotation_names_row=F, annotation_names_col=F, 
        annotation_col=an.y, annotation_row=an.x, annotation_colors=an.col,
        show_rownames=F, show_colnames=F, 
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
ph@column_dend_param$show<-F
ph@row_dend_param$show<-F
ph<-draw(ph, annotation_legend_list=lg)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_cor_ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_circRNAs.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}



#  [exons] estimate significance in the differences observed
#
#  circRNA cluster 1
wt(circs=1, rbps=c(1, 2), m=ex.m, alt='greater')  #  7.458e-10
#
#  circRNA cluster 2
wt(circs=2, rbps=c(1, 2), m=ex.m, alt='greater')  #  9.42e-05
#
#  RBP cluster 1
wt(circs=c(1, 2), rbps=1, m=ex.m, alt='less')  #  2.986e-221
#
#  RBP cluster 2
wt(circs=c(1, 2), rbps=2, m=ex.m, alt='less')  #  1.588e-242
# 
#  circRNA cluster 1 + RBP cluster 2 vs circRNA cluster 2 + RBP cluster 1
wt(circs=c(1, 2), rbps=c(2, 1), m=ex.m, alt='less')  #  1.152e-267



#  [introns] estimate significance in the differences observed
#
#  circRNA cluster 1
wt(circs=1, rbps=c(1, 2), m=in.m, alt='greater')  #  1.736e-12
#
#  circRNA cluster 2
wt(circs=2, rbps=c(1, 2), m=in.m, alt='greater')  #  2.404e-05
#
#  RBP cluster 1
wt(circs=c(1, 2), rbps=1, m=in.m, alt='less')  #  3.036e-125
#
#  RBP cluster 2
wt(circs=c(1, 2), rbps=2, m=in.m, alt='less')  #  8.609e-121
# 
#  circRNA cluster 1 + RBP cluster 2 vs circRNA cluster 2 + RBP cluster 1
wt(circs=c(1, 2), rbps=c(2, 1), m=in.m, alt='less')  #  6.064e-167



#  boxplots of motif distribution densities
#{{{

B<-data.frame(rbind(cbind(ex.m, type='exons'), cbind(in.m, type='introns')))
B$type<-factor(B$type, levels=c('exons', 'introns'))
B$circCluster<-factor(as.character(B$circCluster), levels=as.character(seq_len(NCLUST)))
B$RBPCluster<-factor(as.character(B$RBPCluster), levels=as.character(seq_len(NCLUST)))
B$counts<-log10(B$counts)
B$ncounts<-log10(B$ncounts)
stopifnot( nrow(B)==nrow(ex.m)+nrow(in.m) )
YTICK<-pretty(range(B$ncounts, na.rm=T), 4)
#
#  circRNAs groups by RBPs
#
ggplot(B, aes(x=circCluster, ncounts, fill=RBPCluster)) + 
    stat_boxplot(geom='errorbar', width=0.2, position=position_dodge(width=0.8), linetype='solid', coef=Inf, size=0.4) +
    geom_boxplot(aes(ymin=..lower.., ymax=..upper..), position=position_dodge(width=0.8), coef=Inf, color='white') +
    facet_wrap(~type) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    scale_fill_manual(values=B.col[['RBPCluster']]) +
    labs(x='circRNA cluster', y=expression(log[10]('Normalized counts')), color='') +
    guides(fill=guide_legend('RBP cluster')) +
    coord_cartesian(clip='off') + 
    theme(axis.ticks.x=element_blank(),
          axis.text=element_text(size=45),
          axis.title.x=element_text(color='black', margin=margin(t=1, b=0.5, unit='cm')),  #  we need this or labels will be clipped
          axis.text.x=element_text(margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5, color=B.col[['circCluster']]),
          legend.title=element_text(margin=margin(l=0, b=0, t=0, r=0, unit='cm')),
          legend.justification=c(0, 1), 
          legend.position=c(0.05, 0.90)
    )
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_cor_ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_circRNAs_normalized_counts_by_RBPs.svg', width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  RBPs groups by circRNAs
#
ggplot(B, aes(x=RBPCluster, ncounts, fill=circCluster)) + 
    stat_boxplot(geom='errorbar', width=0.2, position=position_dodge(width=0.8), linetype='solid', coef=Inf, size=0.4) +
    geom_boxplot(aes(ymin=..lower.., ymax=..upper..), position=position_dodge(width=0.8), coef=Inf, color='white') +
    facet_wrap(~type) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    scale_fill_manual(values=B.col[['circCluster']]) +
    labs(x='RBP cluster', y=expression(log[10]('Normalized counts')), color='') +
    guides(fill=guide_legend('circRNA cluster')) +
    coord_cartesian(clip='off') + 
    theme(axis.ticks.x=element_blank(),
          axis.text=element_text(size=45),
          axis.title.x=element_text(color='black', margin=margin(t=1, b=0.5, unit='cm')),  #  we need this or labels will be clipped
          axis.text.x=element_text(margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5, color=B.col[['RBPCluster']]),
          legend.title=element_text(margin=margin(l=0, b=0, t=0, r=0, unit='cm')),
          legend.justification=c(0, 1), 
          legend.position=c(0.05, 0.90)
    )
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_cor_ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_circRNAs_normalized_counts_by_circRNAs.svg', width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  ecdfs of Spearman correlations per circRNA and RBP cluster
#{{{

par(mar=c(5.5, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
XLIM<-pretty(c(-0.5, 1.0), 10)
B<-list()
B[[1]]<-c(CO.ALL[ names(w.x[ w.x==1 ]), names(w.y[ w.y==1 ])])
B[[2]]<-c(CO.ALL[ names(w.x[ w.x==1 ]), names(w.y[ w.y==2 ])])
B[[3]]<-c(CO.ALL[ names(w.x[ w.x==2 ]), names(w.y[ w.y==1 ])])
B[[4]]<-c(CO.ALL[ names(w.x[ w.x==2 ]), names(w.y[ w.y==2 ])])
B.col<-colorRampPalette(c('lightsalmon3', 'lightgreen'))(4)
curve(ecdf(B[[1]])(x), from=XLIM[1], to=tail(XLIM, 1), n=100, ylab='', xlab='', pch=NA, col=B.col[1], lty=1, lwd=12, main='', xlim=range(XLIM), ylim=c(0, 1))
curve(ecdf(B[[2]])(x), from=XLIM[1], to=tail(XLIM, 1), n=100, pch=NA, col=B.col[2], lty=1, lwd=12, add=T)
curve(ecdf(B[[3]])(x), from=XLIM[1], to=tail(XLIM, 1), n=100, pch=NA, col=B.col[3], lty=1, lwd=12, add=T)
curve(ecdf(B[[4]])(x), from=XLIM[1], to=tail(XLIM, 1), n=100, pch=NA, col=B.col[4], lty=1, lwd=12, add=T)
mtext('Spearman correlation', side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.1, cex=2.4, las=0) 
legend('topleft', legend=c('circRNA 1 RBP 1', 'circRNA 1 RBP 2', 'circRNA 2 RBP 1', 'circRNA 2 RBP 2'), col=B.col, bty='n', pch=NA, lty=1, lwd=15, cex=1.6, y.intersp=0.60, x.intersp=0.4, seg.len=0.4)
text(x=XLIM[1], y=1.025, labels='Clusters', adj=-0.5, col='black', cex=1.6, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_cor_ATtRACT_RBP+SpliceAid_SF+MSigD_SF_circRNAs_by_cluster.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  check: differential expression and logFC distribution of the RBPs in the clusters
#{{{

#  load [MNA vs HR_nMNA] DE results
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_MNA_HR_nMNA.RData')
rm(list=setdiff(ls(), c(l, 'RES')))


#  isolate the RBPs
rbp<-subset(RES, gene_name %in% ALL$gene_name)[, c('baseMean', 'log2FoldChange', 'padj', 'gene_name')]
rbp<-data.table(gene_id=rownames(rbp), data.frame(rbp))[ order(-abs(log2FoldChange)) ][, cluster:=w.y[ gene_name ]]
rbp[ cluster==1, summary(log2FoldChange)]
#
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.3072 -0.2544 -0.0617 -0.0611  0.1963  0.7790 
#
rbp[ cluster==2, summary(log2FoldChange)]
#
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.4695  0.0504  0.2882  0.2853  0.4849  1.5664 
wilcox.test(x=rbp[ cluster==1, log2FoldChange], y=rbp[ cluster==2, log2FoldChange], alternative='less')$p.value        #  1.686e-14
#
#  slightly higher expressed RBPs in cluster 1
wilcox.test(x=rbp[ cluster==1, log10(baseMean)], y=rbp[ cluster==2, log10(baseMean)], alternative='greater')$p.value   #  0.02156
#
#  top 100 expressed similarly represented in the two clusters
rbp[ order(-baseMean) ][ 1:100, table(cluster) ]
#
#  1  2 
# 43 57 
#
#  top 100 absolute log2FCs of well-expressed RBPs are unequally represented in the two clusters
#  cluster 2 has the strongest absolute log2FCs but it is also bigger
rbp[ order(-abs(log2FoldChange)) ][ baseMean>=1000 ][1:100, table(cluster) ]
# cluster
#  1  2 
# 29 71 


#  [cluster 2, significantly up-DE genes] Fisher's exact test:
#
#                   |     up-DE        |      not up-DE       |
#  -----------------|------------------|----------------------|
#      cluster 2    |      aa          |          ab          |
#  -----------------|------------------|----------------------|
#      cluster 1    |      ba          |          bb          | 
#  -----------------|------------------|----------------------|
up<-subset(RES, log2FoldChange>0 & padj<0.05)$gene_name
aa<-length(intersect(names(w.y[ w.y==2 ]), up))
ba<-length(intersect(names(w.y[ w.y==1 ]), up))
notup<-setdiff(RES$gene_name, up)
ab<-length(intersect(names(w.y[ w.y==2 ]), notup))
bb<-length(intersect(names(w.y[ w.y==1 ]), notup))
fisher.test(data.frame('in'=c(aa, ba), 'out'=c(ab, bb)), alternative='greater')$p.value
#
#  => 5.423e-09 



#  [cluster 1, significantly down-DE genes] Fisher's exact test:
#
#                   |   down-DE        |    not down-DE       |
#  -----------------|------------------|----------------------|
#      cluster 1    |      aa          |          ab          |
#  -----------------|------------------|----------------------|
#      cluster 2    |      ba          |          bb          | 
#  -----------------|------------------|----------------------|
down<-subset(RES, log2FoldChange<0 & padj<0.05)$gene_name
aa<-length(intersect(names(w.y[ w.y==1 ]), down))
ba<-length(intersect(names(w.y[ w.y==2 ]), down))
notdown<-setdiff(RES$gene_name, down)
ab<-length(intersect(names(w.y[ w.y==1 ]), notdown))
bb<-length(intersect(names(w.y[ w.y==2 ]), notdown))
fisher.test(data.frame('in'=c(aa, ba), 'out'=c(ab, bb)), alternative='greater')$p.value
# 
#  => 1.707e-06 

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_cor_circRNAs.RData



#  [high risk tumors] Correlation of RBP (6mer motifs or longer) including SF expression with circRNA expression (based on CIRI2 counts).
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
library(grid)
library(ggplot2)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
options(scipen=0)



#  [run once] Spearman correlations between the union of RBPs and SFs and circRNA counts (based on CIRI2)
#             cluster-specific enrichment of RBP motifs in circRNAs and flanking introns 
#{{{

#  union of SF
l<-ls()
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/SpliceAid-F+GO.RData')
load('/fast/groups/ag_schulte/work/reference/annotation/human_splicing/MSigDB_c2_v7.0_splicing.RData')
SF<-data.table(unique(rbind(data.frame(db[, c('gene_name', 'gene_id')]), msigdbSF[, c('gene_name', 'gene_id')])))
rm(list=setdiff(ls(), c(l, 'SF')))


#  load all RBPs 
#  inflate per motif the subgroups based on PWMs into one group of unique genes (motif-inflating needs to happen alone, different 1:many map than genes)
#  and keep those with 6mer or longer motifs
#  keep 6mers or longer
l<-ls()
load('/fast/groups/ag_schulte/work/reference/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')
RBP<-unique(db[, .(gene_name=unique(unlist(gene_name)), gene_id=unique(unlist(gene_id)), motif=motif), by=.(id)][, .(motif=unlist(motif)), by=.(gene_name, gene_id)][ sapply(motif, nchar)>=6, ][, motif:=NULL])
rm(list=setdiff(ls(), c(l, 'RBP')))


#  join RBPs and SFs into one RBP list
length(setdiff(SF$gene_name, RBP$gene_name))  #  178
ALL<-unique(rbind(RBP, SF))
rm(RBP)


#  load prepackaged gene counts and metadata
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')
meta<-meta[ !(failed) & risk_group %in% c('HR_nMNA', 'MNA') ]
ALL.cnt<-gns.cnt[ rownames(gns.cnt) %in% ALL$gene_id, colnames(gns.cnt) %in% meta$bid ]
#
ALL[ gene_id %in% setdiff(ALL$gene_id, rownames(ALL.cnt)) ]  #  one RBP not expressed
# 
#  RBMY1A1 ENSG00000234414.7
#
rownames(ALL.cnt)<-ALL[ match(rownames(ALL.cnt), gene_id), gene_name ]
rm(list=setdiff(ls(), c(l, 'meta', 'ALL.cnt')))


#  load CIRI2 counts and convert to DGE
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
CIR.cnt<-data.table(data.frame(mcols(nb.circs)[, c('circ_name', 'bid', 'jc_count')]))
CIR.cnt<-dcast(CIR.cnt, bid ~ circ_name, value.var='jc_count', fun.aggregate=mean, fill=0.0)
CIR.cnt<-t(as.matrix(data.frame(CIR.cnt[, -1], row.names=CIR.cnt[, bid], check.names=F)))
rm(list=setdiff(ls(), c(l, 'meta', 'CIR.cnt')))


#  fix the sample order 
stopifnot( length(setdiff( colnames(ALL.cnt), colnames(CIR.cnt) ))==0 )
CIR.cnt<-CIR.cnt[, colnames(ALL.cnt)]
stopifnot( all.equal( colnames(CIR.cnt), colnames(ALL.cnt) ) )


#  [RBPs+SFs] compute Spearman correlations and cluster based on Euclidean distance of 1 - Spearman correlations
CO.ALL<-cor(t(CIR.cnt), t(ALL.cnt), use='pairwise.complete.obs', method='spearman')
H.ALL.CIR<-hclust(dist(1-CO.ALL, method='euclidean'), method='ward.D2')
H.ALL<-hclust(dist(1-t(CO.ALL), method='euclidean'), method='ward.D2')


#  [RBPs-SFs] compute Spearman correlations and cluster based on Euclidean distance of 1 - Spearman correlations
CO.RBP<-cor(t(CIR.cnt), t(ALL.cnt[ !rownames(ALL.cnt) %in% SF$gene_name, ]), use='pairwise.complete.obs', method='spearman')
H.RBP.CIR<-hclust(dist(1-CO.RBP, method='euclidean'), method='ward.D2')
H.RBP<-hclust(dist(1-t(CO.RBP), method='euclidean'), method='ward.D2')


#  [SFs] compute Spearman correlations and cluster based on Euclidean distance of 1 - Spearman correlations
CO.SF<-cor(t(CIR.cnt), t(ALL.cnt[ SF$gene_name, ]), use='pairwise.complete.obs', method='spearman')
H.SF.CIR<-hclust(dist(1-CO.SF, method='euclidean'), method='ward.D2')
H.SF<-hclust(dist(1-t(CO.SF), method='euclidean'), method='ward.D2')


#  strong correlation cutoff
COR_CUTOFF<-0.40


#  [RBPs+SFs] identify strongly correlated pairs 
#             identify the self-correlations among them
#             identify all RBP-related correlations (including self)
#             cluster strong correlating pairs based on Euclidean distance of 1 - Spearman correlations
ALL.strong<-apply(CO.ALL, 2, function(x){ setNames(x[x>=COR_CUTOFF], names(x[x>=COR_CUTOFF])) })
ALL.strong<-ALL.strong[ lengths(ALL.strong)>0 ]
ALL.strong.self<-sapply(names(ALL.strong), function(n){ ALL.strong[[n]][grepl(n, names(ALL.strong[[n]]))] })
ALL.strong.self<-ALL.strong.self[ lengths(ALL.strong.self)>0 ]
a<-paste(ALL$gene_name, collapse='|')
ALL.strong.all<-sapply(ALL.strong, function(n){ n[grepl(a, names(n))] })
ALL.strong.all<-ALL.strong.all[lengths(ALL.strong.all)>0]
a<-unique(unname(unlist(lapply(ALL.strong, names))))
b<-names(ALL.strong)
CO.ALL.strong<-CO.ALL[ a, b]
H.ALL.strong.CIR<-hclust(dist(1-CO.ALL.strong, method='euclidean'), method='ward.D2')
H.ALL.strong<-hclust(dist(1-t(CO.ALL.strong), method='euclidean'), method='ward.D2')
rm(a,b)


#  [RBPs-SFs] identify strongly correlated pairs 
#             identify the self-correlations among them
#             identify all RBP-related correlations (including self)
#             cluster strong correlating pairs based on Euclidean distance of 1 - Spearman correlations
RBP.strong<-apply(CO.RBP, 2, function(x){ setNames(x[x>=COR_CUTOFF], names(x[x>=COR_CUTOFF])) })
RBP.strong<-RBP.strong[ lengths(RBP.strong)>0 ]
RBP.strong.self<-sapply(names(RBP.strong), function(n){ RBP.strong[[n]][grepl(n, names(RBP.strong[[n]]))] })
RBP.strong.self<-RBP.strong.self[ lengths(RBP.strong.self)>0 ]
a<-paste(ALL$gene_name, collapse='|')
RBP.strong.all<-sapply(RBP.strong, function(n){ n[grepl(a, names(n))] })
RBP.strong.all<-RBP.strong.all[lengths(RBP.strong.all)>0]
a<-unique(unname(unlist(lapply(RBP.strong, names))))
b<-names(RBP.strong)
CO.RBP.strong<-CO.RBP[ a, b]
H.RBP.strong.CIR<-hclust(dist(1-CO.RBP.strong, method='euclidean'), method='ward.D2')
H.RBP.strong<-hclust(dist(1-t(CO.RBP.strong), method='euclidean'), method='ward.D2')
rm(a,b)


#  [SFs] identify strongly correlated pairs 
#        identify the self-correlations among them
#        identify all SF-related correlations (including self)
#        cluster strong correlating pairs based on Euclidean distance of 1 - Spearman correlations
SF.strong<-apply(CO.SF, 2, function(x){ setNames(x[x>=COR_CUTOFF], names(x[x>=COR_CUTOFF])) })
SF.strong<-SF.strong[ lengths(SF.strong)>0 ]
SF.strong.self<-sapply(names(SF.strong), function(n){ SF.strong[[n]][grepl(n, names(SF.strong[[n]]))] })
SF.strong.self<-SF.strong.self[ lengths(SF.strong.self)>0 ]
a<-paste(ALL$gene_name, collapse='|')
SF.strong.all<-sapply(SF.strong, function(n){ n[grepl(a, names(n))] })
SF.strong.all<-SF.strong.all[lengths(SF.strong.all)>0]
a<-unique(unname(unlist(lapply(SF.strong, names))))
b<-names(SF.strong)
CO.SF.strong<-CO.SF[ a, b]
H.SF.strong.CIR<-hclust(dist(1-CO.SF.strong, method='euclidean'), method='ward.D2')
H.SF.strong<-hclust(dist(1-t(CO.SF.strong), method='euclidean'), method='ward.D2')
rm(a,b)


#  cut at this many clusters both RBPs and circRNAs
NCLUST<-2
w.x<-cutree(H.ALL.CIR, k=NCLUST)
w.y<-cutree(H.ALL, k=NCLUST)
table(w.x)
#
#    1    2 
# 2260 2943 
# 
table(w.y)
#
#   1   2 
# 138 177 



#  which clusters do circARID1A, KHSRP and DHX9 belong to?
w.x[ grep('ARID1A', names(w.x)) ]
#
# ENSG00000117713.20|ARID1A_chr1+26729651-26732792 
#                                                1 
w.y[ grep('KHSRP|DHX9', names(w.y)) ]
#
#  DHX9 KHSRP 
#     2     2 



#  summary of correlations for KHSRP and DHX9
summary(CO.ALL[, 'KHSRP'])
#
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.5089 -0.0439  0.0662  0.0720  0.1847  0.6826 
#
summary(CO.ALL[, 'DHX9'])
#
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  -0.270   0.106   0.196   0.193   0.281   0.611 
#
table(CO.ALL[, 'DHX9']<0)
# 
# FALSE  TRUE 
#  4850   353 



#  load circRNA RBP motif enrichment results
#  isolate the circRNAs and RBPs of interest
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
fimo.ex<-fimo.all[ circ_name %in% rownames(CO.ALL) & sapply(rbp_group, function(r){ any( r %in% colnames(CO.ALL) ) }) ][, c('motif_id', 'motif_alt_id', 'start', 'end', 'score', 'pvalue', 'qvalue'):=NULL] 
rm(list=setdiff(ls(), c(l, 'fimo.ex')))



#  load circRNA intron RBP motif enrichment results
#  isolate the circRNAs and RBPs of interest
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/introns/results.RData')
fimo.in<-fimo.all[ circ_name %in% rownames(CO.ALL) & sapply(rbp_group, function(r){ any( r %in% colnames(CO.ALL) ) }) ][, c('motif_id', 'motif_alt_id', 'start', 'end', 'score', 'pvalue', 'qvalue'):=NULL] 
rm(list=setdiff(ls(), c(l, 'fimo.in')))



#  count number of motifs per circRNAs in a cluster-specific way
#
#  N.B. The number of motifs is biased by the number of RBPs and circRNAs in each cluster. We need normalized numbers to compare.
#
ex.m<-in.m<-data.frame(circCluster=integer(0), RBPluster=integer(0), counts=integer(0), ncounts=double(0))
for (i in 1:NCLUST){      #  circRNA clusters
    N<-sum(w.x==i)
    for (j in 1:NCLUST){  #  RBP clusters
        M<-N*sum(w.y==j)
        x<-fimo.ex[ circ_name %in% names( w.x[ w.x==i ]) & sapply(rbp_group, function(r){ any( r %in% names( w.y[ w.y==j ] )) }), ][, .(counts=.N), by=.(circ_name)][, counts]
        ex.m<-rbind(ex.m, data.frame(circCluster=i, RBPCluster=j, counts=x, ncounts=x/M))
        x<-fimo.in[ circ_name %in% names( w.x[ w.x==i ]) & sapply(rbp_group, function(r){ any( r %in% names( w.y[ w.y==j ] )) }), ][, .(counts=.N), by=.(circ_name)][, counts]
        in.m<-rbind(in.m, data.frame(circCluster=i, RBPCluster=j, counts=x, ncounts=x/M))
    }
}
ex.m<-data.table(ex.m)
in.m<-data.table(in.m)
rm(i, N, j, M, x)


#  save
save(list=ls(), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_cor_circRNAs_high-risk_tumors.RData')

#}}}



#  load back
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_cor_circRNAs_high-risk_tumors.RData')



#  ggplot2 theme
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



#  functions
#{{{
wt<-function(circs=c(1, 2), rbps=c(1,2), m, alt='greater'){
    stopifnot( length(rbps)==2 | length(circs)==2 )
    if (length(rbps)==2 & length(circs)==1){
        wt<-wilcox.test(x=m[ circCluster %in% circs & RBPCluster %in% rbps[1], ncounts], y=m[ circCluster %in% circs & RBPCluster %in% rbps[2], ncounts], alternative=alt)$p.value
    } else if (length(circs)==2 & length(rbps)==1){
        wt<-wilcox.test(x=m[ circCluster %in% circs[1] & RBPCluster %in% rbps, ncounts], y=m[ circCluster %in% circs[2] & RBPCluster %in% rbps, ncounts], alternative=alt)$p.value
    } else if (length(circs)==length(rbps) & length(circs)==2){
        wt<-wilcox.test(x=m[ circCluster %in% circs[1] & RBPCluster %in% rbps[1], ncounts], y=m[ circCluster %in% circs[2] & RBPCluster %in% rbps[2], ncounts], alternative=alt)$p.value
    }
    return(wt)
}
#}}}



#  check: do the RBPs with which the circARID1A is highly correlated have significant motifs on it or on its introns?
#{{{

#  load circARID1A FIMO results 
#  N.B. RBPs have p-values but no filtering was done
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_RBPs_across_tissues.RData')
rbp.arid1a.ex<-rbp
rm(list=setdiff(ls(), c(l, 'rbp.arid1a.ex')))


#  load circARID1A intron FIMO results 
#  N.B. RBPs have p-values but no filtering was done
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_introns_RBPs_across_tissues.RData')
rbp.arid1a.in<-rbp
rm(list=setdiff(ls(), c(l, 'rbp.arid1a.in')))


#  highly correlated RBPs with circARID1A
x<-CO.ALL[ grepl('ARID1A', rownames(CO.ALL)), ]
x<-x[ x>=COR_CUTOFF ]


#  any of the RBPs of interest having significant motifs on circARID1A or its introns?
rbp.arid1a.ex[ gene_name %in% names(x) ]
#    gene_name           gene_id ncount     pv.vt     pv.hb    pvalue
# 1:   HNRNPH2 ENSG00000126945.9      1 3.072e-24 1.238e-20 1.239e-20
# 2:     PCBP1 ENSG00000169564.6      5 2.382e-03 1.386e-16 2.382e-03
# 3:     ZFP36 ENSG00000128016.6      3 2.910e-02 2.597e-18 2.910e-02
rbp.arid1a.in[ gene_name %in% names(x) ]
#     gene_name            gene_id ncount     pv.vt     pv.hb    pvalue
#  1:     NOVA2  ENSG00000104967.7     32 5.320e-46 1.419e-21 1.419e-21
#  2:     CELF5 ENSG00000161082.13     12 1.642e-39 4.897e-21 4.897e-21
#  3:   HNRNPH2  ENSG00000126945.9    202 2.685e-24 1.203e-20 1.203e-20
#  4:      RBM4 ENSG00000173933.20     53 1.053e-20 8.575e-21 1.911e-20
#  5:      SFPQ ENSG00000116560.11    137 6.227e-29 3.635e-19 3.635e-19
#  6:    RBFOX1 ENSG00000078328.20     10 8.834e-47 1.733e-14 1.733e-14
#  7:     SRSF4 ENSG00000116350.17     17 1.081e-13 7.952e-24 1.081e-13
#  8:     RBM8A  ENSG00000265241.6      4 5.119e-12 2.023e-19 5.119e-12
#  9:    TARDBP ENSG00000120948.17     32 1.781e-25 1.278e-10 1.278e-10
# 10:      PUM1 ENSG00000134644.15     53 1.928e-30 1.874e-10 1.874e-10
# 11:     SRSF9  ENSG00000111786.9     45 8.588e-09 6.158e-18 8.588e-09
# 12:     CELF6 ENSG00000140488.16     18 6.449e-43 6.604e-08 6.604e-08
# 13:   ZFP36L2  ENSG00000152518.8      3 1.836e-07 8.254e-20 1.836e-07
# 14:     MBNL1 ENSG00000152601.17     71 1.027e-35 1.008e-06 1.008e-06
# 15:     CELF1 ENSG00000149187.18     33 5.166e-16 1.751e-05 1.751e-05
# 16:      AGO1 ENSG00000092847.12     31 2.965e-15 2.888e-04 2.888e-04
# 17:   IGHMBP2  ENSG00000132740.8     29 7.338e-04 7.108e-09 7.338e-04
# 18:     PCBP1  ENSG00000169564.6    142 2.423e-03 1.508e-16 2.423e-03
# 19:      FXR2 ENSG00000129245.11      9 1.607e-22 1.471e-02 1.471e-02
# 20:     ZFP36  ENSG00000128016.6     51 3.048e-02 2.815e-18 3.048e-02
# 21:     CSTF2 ENSG00000101811.13     36 5.503e-02 3.874e-05 5.506e-02
# 22:     CELF4 ENSG00000101489.20     12 1.586e-44 6.001e-02 6.001e-02
# 23:     IFIH1  ENSG00000115267.7     34 8.836e-35 1.302e-01 1.302e-01
# 24:    SRSF10 ENSG00000188529.14     47 5.322e-01 2.061e-11 5.322e-01
# 25:     DDX58 ENSG00000107201.10      2 5.166e-16 5.386e-01 5.386e-01
# 26:     RBM42 ENSG00000126254.12      2 4.453e-35 7.378e-01 7.378e-01
# 27:     SRP14 ENSG00000140319.10     26 8.783e-01 3.973e-07 8.783e-01
#     gene_name            gene_id ncount     pv.vt     pv.hb    pvalue
#}}}



#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
B.col<-list(RBPCluster=setNames(c('brown1', 'darkorange4'), 1:NCLUST), circCluster=setNames(c('darkgoldenrod', 'darkolivegreen'), 1:NCLUST))



#  histograms of correlations
#{{{

#  [ALL] histogram of all correlations
par(mar=c(5.0, 6.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#b<-pretty(c(min(CO.ALL), max(CO.ALL)), 10)
b<-pretty(c(-0.8, 1.0), 10)
h<-hist(CO.ALL, breaks=seq(min(b), max(b), length.out=50), plot=F)
h$counts<-log10(1+h$counts)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
h<-plot(h, col='darkgrey', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(b), xaxt='n', add=F)
axis(1, at=seq(min(b), max(b), by=0.4), line=0, cex.axis=2.4)
#abline(v=COR_CUTOFF, col='red4', lty=2, lwd=4, xpd=F)
mtext('Spearman correlation', side=1, line=4, padj=-0.5, las=0, cex=2.4)
mtext(expression(log[10](1+'Frequency')), side=2, line=2, padj=-0.4, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_cor_ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_circRNAs_high-risk_tumors.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b,h,YMAX)



#  [SF] histogram of all correlations
par(mar=c(5.0, 6.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#b<-pretty(c(min(CO.SF), max(CO.SF)), 10)
b<-pretty(c(-0.8, 1.0), 10)
h<-hist(CO.SF, breaks=seq(min(b), max(b), length.out=50), plot=F)
h$counts<-log10(1+h$counts)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
h<-plot(h, col='darkgrey', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(b), xaxt='n', add=F)
axis(1, at=seq(min(b), max(b), by=0.4), line=0, cex.axis=2.4)
#abline(v=COR_CUTOFF, col='red4', lty=2, lwd=4, xpd=F)
mtext('Spearman correlation', side=1, line=4, padj=-0.5, las=0, cex=2.4)
mtext(expression(log[10](1+'Frequency')), side=2, line=2, padj=-0.4, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_cor_SpliceAid_SF+MSigDB_SF_circRNAs_high-risk_tumors.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b,h,YMAX)

#}}}



#  heatmap 
#{{{

grid.newpage()
an.y<-data.frame('RBP clusters'=factor(cutree(H.ALL, k=NCLUST), levels=1:NCLUST), check.names=F) 
an.x<-data.frame('circRNA clusters'=factor(cutree(H.ALL.CIR, k=NCLUST), levels=1:NCLUST), check.names=F)
an.col<-list('RBP clusters'=setNames(c('brown1', 'darkorange4'), 1:NCLUST), 'circRNA clusters'=setNames(c('darkgoldenrod', 'darkolivegreen'), 1:NCLUST))
BREAKS<-c(-0.6, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8)
BREAKS.COL<-rev(colorRampPalette(brewer.pal(5,'RdBu'))((length(BREAKS)-1)*2))
BREAKS.COL<-BREAKS.COL[c(seq(1, by=2, length.out=length(BREAKS)-2), length(BREAKS.COL))]  #  maximum contrast for the last bin
lg<-packLegend(
    Legend(title='', title_gp=gpar(fontsize=18, fontface='bold'), labels_gp=gpar(fontsize=16), grid_height=unit(4,'mm'),
        col_fun=circlize::colorRamp2(BREAKS[-1], BREAKS.COL), at=c(min(BREAKS), 0.2, max(BREAKS)), ncol=1), 
    Legend(title=colnames(an.y), title_gp=gpar(fontsize=18, fontface='bold'), labels_gp=gpar(fontsize=16), grid_height=unit(4,'mm'),
        at=names(an.col[[colnames(an.y)]]), legend_gp=gpar(fill=an.col[[colnames(an.y)]]), ncol=1), 
    Legend(title=colnames(an.x), title_gp=gpar(fontsize=18, fontface='bold'), labels_gp=gpar(fontsize=16), grid_height=unit(4,'mm'),
        at=names(an.col[[colnames(an.x)]]), legend_gp=gpar(fill=an.col[[colnames(an.x)]]), ncol=1), 
    column_gap=unit(1, 'cm'), max_height=unit(1, 'npc'))
ph<-ComplexHeatmap::pheatmap(as.matrix(CO.ALL), color=BREAKS.COL, border_color=NA, scale='none', 
        breaks=BREAKS,
        use_raster=F, 
        cluster_rows=H.ALL.CIR,
        cluster_cols=H.ALL,
        cutree_row=NCLUST,
        cutree_col=NCLUST,
        legend=F,
        annotation_legend=F,
        annotation_names_row=F, annotation_names_col=F, 
        annotation_col=an.y, annotation_row=an.x, annotation_colors=an.col,
        show_rownames=F, show_colnames=F, 
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
ph@column_dend_param$show<-F
ph@row_dend_param$show<-F
ph<-draw(ph, annotation_legend_list=lg)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_cor_ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_circRNAs_high-risk_tumors.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}



#  [exons] estimate significance in the differences observed
#
#  circRNA cluster 1
wt(circs=1, rbps=c(1, 2), m=ex.m, alt='less')  #  1.608e-38
#
#  circRNA cluster 2
wt(circs=2, rbps=c(1, 2), m=ex.m, alt='less')  #  1.358e-47
#
#  RBP cluster 1
wt(circs=c(1, 2), rbps=1, m=ex.m, alt='greater')  #  1.842e-07
#
#  RBP cluster 2
wt(circs=c(1, 2), rbps=2, m=ex.m, alt='greater')  #  1.97e-09
#
#  circRNA cluster 1 + RBP cluster 1 vs circRNA cluster 2 + RBP cluster 2
wt(circs=c(1, 2), rbps=c(1, 2), m=ex.m, alt='less')  #  2.171e-14



#  [introns] estimate significance in the differences observed
#
#  circRNA cluster 1
wt(circs=1, rbps=c(1, 2), m=in.m, alt='less')  #  1.216e-06
#
#  circRNA cluster 2
wt(circs=2, rbps=c(1, 2), m=in.m, alt='less')  #  1.941e-09
#
#  RBP cluster 1
wt(circs=c(1, 2), rbps=1, m=in.m, alt='greater')  #  0.005729
#
#  RBP cluster 2
wt(circs=c(1, 2), rbps=2, m=in.m, alt='greater')  #  0.01416
# 
#  circRNA cluster 1 + RBP cluster 1 vs circRNA cluster 2 + RBP cluster 2
wt(circs=c(1, 2), rbps=c(1, 2), m=in.m, alt='less')  #  0.001944



#  boxplots of motif distribution densities
#{{{

B<-data.frame(rbind(cbind(ex.m, type='exons'), cbind(in.m, type='introns')))
B$type<-factor(B$type, levels=c('exons', 'introns'))
B$circCluster<-factor(as.character(B$circCluster), levels=as.character(seq_len(NCLUST)))
B$RBPCluster<-factor(as.character(B$RBPCluster), levels=as.character(seq_len(NCLUST)))
B$counts<-log10(B$counts)
B$ncounts<-log10(B$ncounts)
stopifnot( nrow(B)==nrow(ex.m)+nrow(in.m) )
YTICK<-pretty(range(B$ncounts, na.rm=T), 4)
#
#  circRNAs groups by RBPs
#
ggplot(B, aes(x=circCluster, ncounts, fill=RBPCluster)) + 
    stat_boxplot(geom='errorbar', width=0.2, position=position_dodge(width=0.8), linetype='solid', coef=Inf, size=0.4) +
    geom_boxplot(aes(ymin=..lower.., ymax=..upper..), position=position_dodge(width=0.8), coef=Inf, color='white') +
    facet_wrap(~type) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    scale_fill_manual(values=B.col[['RBPCluster']]) +
    labs(x='circRNA cluster', y=expression(log[10]('Normalized counts')), color='') +
    guides(fill=guide_legend('RBP cluster')) +
    coord_cartesian(clip='off') + 
    theme(axis.ticks.x=element_blank(),
          axis.text=element_text(size=45),
          axis.title.x=element_text(color='black', margin=margin(t=1, b=0.5, unit='cm')),  #  we need this or labels will be clipped
          axis.text.x=element_text(margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5, color=B.col[['circCluster']]),
          legend.title=element_text(margin=margin(l=0, b=0, t=0, r=0, unit='cm')),
          legend.justification=c(0, 1), 
          legend.position=c(0.05, 0.90)
    )
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_cor_ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_circRNAs_normalized_counts_by_RBPs_high-risk_tumors.svg', width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  RBPs groups by circRNAs
#
ggplot(B, aes(x=RBPCluster, ncounts, fill=circCluster)) + 
    stat_boxplot(geom='errorbar', width=0.2, position=position_dodge(width=0.8), linetype='solid', coef=Inf, size=0.4) +
    geom_boxplot(aes(ymin=..lower.., ymax=..upper..), position=position_dodge(width=0.8), coef=Inf, color='white') +
    facet_wrap(~type) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    scale_fill_manual(values=B.col[['circCluster']]) +
    labs(x='RBP cluster', y=expression(log[10]('Normalized counts')), color='') +
    guides(fill=guide_legend('circRNA cluster')) +
    coord_cartesian(clip='off') + 
    theme(axis.ticks.x=element_blank(),
          axis.text=element_text(size=45),
          axis.title.x=element_text(color='black', margin=margin(t=1, b=0.5, unit='cm')),  #  we need this or labels will be clipped
          axis.text.x=element_text(margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=1, hjust=0.5, color=B.col[['RBPCluster']]),
          legend.title=element_text(margin=margin(l=0, b=0, t=0, r=0, unit='cm')),
          legend.justification=c(0, 1), 
          legend.position=c(0.05, 0.90)
    )
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_cor_ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_circRNAs_normalized_counts_by_circRNAs_high-risk_tumors.svg', width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}



#  ecdfs of Spearman correlations per circRNA and RBP cluster
#{{{

par(mar=c(5.5, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
XLIM<-pretty(c(-0.5, 1.0), 10)
B<-list()
B[[1]]<-c(CO.ALL[ names(w.x[ w.x==1 ]), names(w.y[ w.y==1 ])])
B[[2]]<-c(CO.ALL[ names(w.x[ w.x==1 ]), names(w.y[ w.y==2 ])])
B[[3]]<-c(CO.ALL[ names(w.x[ w.x==2 ]), names(w.y[ w.y==1 ])])
B[[4]]<-c(CO.ALL[ names(w.x[ w.x==2 ]), names(w.y[ w.y==2 ])])
B.col<-colorRampPalette(c('lightsalmon3', 'lightgreen'))(4)
curve(ecdf(B[[1]])(x), from=XLIM[1], to=tail(XLIM, 1), n=100, ylab='', xlab='', pch=NA, col=B.col[1], lty=1, lwd=12, main='', xlim=range(XLIM), ylim=c(0, 1))
curve(ecdf(B[[2]])(x), from=XLIM[1], to=tail(XLIM, 1), n=100, pch=NA, col=B.col[2], lty=1, lwd=12, add=T)
curve(ecdf(B[[3]])(x), from=XLIM[1], to=tail(XLIM, 1), n=100, pch=NA, col=B.col[3], lty=1, lwd=12, add=T)
curve(ecdf(B[[4]])(x), from=XLIM[1], to=tail(XLIM, 1), n=100, pch=NA, col=B.col[4], lty=1, lwd=12, add=T)
mtext('Spearman correlation', side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.1, cex=2.4, las=0) 
legend('topleft', legend=c('circRNA 1 RBP 1', 'circRNA 1 RBP 2', 'circRNA 2 RBP 1', 'circRNA 2 RBP 2'), col=B.col, bty='n', pch=NA, lty=1, lwd=15, cex=1.6, y.intersp=0.60, x.intersp=0.4, seg.len=0.4)
text(x=XLIM[1], y=1.025, labels='Clusters', adj=-0.5, col='black', cex=1.6, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_cor_ATtRACT_RBP+SpliceAid_SF+MSigD_SF_circRNAs_by_cluster_high-risk_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ATtRACT_RBP+SpliceAid_SF+MSigDB_SF_cor_circRNAs_high-risk_tumors.RData




###################################
#
#
#  scratch/quick and dirty analyses
#
#
###################################




#  circARID1A backspliced junction counts correlation with KHSRP (ENSG00000088247.17) counts across all tumors and stratified by risk group
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/load_gene_expression.R')
source('~/bio/lib/draw_highlights.R')
source('~/bio/lib/enrichment_analysis.R')


#  load prepackaged gene counts and metadata
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')
meta<-meta[ !(failed) & !is.na(risk_group) ]
khsrp<-gns.cnt[ rownames(gns.cnt) %in% 'ENSG00000088247.17', meta$bid , drop=T]
rm(list=setdiff(ls(), c(l, 'meta', 'khsrp')))


#  load CIRI2 counts and isolate circARID1A
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
cir.cnt<-data.table(data.frame(mcols(nb.circs)[, c('circ_name', 'bid', 'jc_count')]))
cir.cnt<-dcast(cir.cnt, bid ~ circ_name, value.var='jc_count', fun.aggregate=mean, fill=0.0)
cir.cnt<-t(as.matrix(data.frame(cir.cnt[, -1], row.names=cir.cnt[, bid], check.names=F)))
arid1a<-cir.cnt[ rownames(cir.cnt) %in% 'ENSG00000117713.20|ARID1A_chr1+26729651-26732792', meta$bid, drop=T]
rm(list=setdiff(ls(), c(l, 'arid1a')))


#  correlations across all samples and per risk group
stopifnot( all.equal(names(khsrp), names(arid1a)) )
co<-cor(arid1a, khsrp, method='spearman')
co.risk<-mapply(function(x, y){ cor(x, y, method='spearman') }, split(arid1a, factor(meta$risk_group, levels=unique(meta$risk_group))), split(khsrp, factor(meta$risk_group, levels=unique(meta$risk_group))))
cl<-setNames(unique(meta$col), unique(meta$risk_group))
N<-meta[, .(count=.N), by=.(risk_group)]
N<-setNames(N$count, N$risk_group)


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel')


#  scatterplot across samples
par(mar=c(5.5, 7.5, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-log10(1+arid1a)
y<-log10(1+khsrp)
YTICK<-pretty(y, 4)
XTICK<-pretty(x, 4)
plot(x, y, type='p', pch=20, col=meta$col, main='', xlab='', ylab='', cex=2.4, xlim=range(XTICK), ylim=range(YTICK))
mtext(bquote(R == .(format(co, digits=2))), side=3, line=-1, padj=+0.5, cex=2.0, xpd=NA)
mtext(bquote(log[10](1+counts) ~ '[circARID1A]'), side=1, line=5, padj=-0.2, cex=2.4)
mtext(bquote(log[10](1+counts) ~ '[KHSRP]'), side=2, line=4, padj=-0.1, cex=2.4, las=0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/scatterplot_counts_circARID1A_vs_KHSRP.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  barplot of correlations per risk group
par(mar=c(10.0, 7.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(0, 1), 4)
bp<-barplot(co.risk, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(YTICK[1], tail(YTICK, 1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(co.risk, border='white', col=cl, axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext('Spearman correlation', side=2, line=5, padj=-0.1, las=0, cex=2.4)
#mtext(text=paste0(names(N), ' (', N, ')'), side=1, line=0, at=bp, las=2, adj=1, cex=2.2, col=cl)
mtext(text=names(N), side=1, line=1, at=bp, las=2, adj=1, cex=2.2, col=cl)
mtext('circARID1A vs KHSRP', side=3, line=-1, las=1, adj=0.5, cex=2.2, col='black')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_counts_correlation_per_risk-group_circARID1A_vs_KHSRP.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}



#  [explore transcript RBP results] check RBP motifs for select genes and RBPs
#{{{
rm(list=ls())
library(data.table)
library(GenomicAlignments)
library(rtracklayer)


#  load results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/results_transcripts.RData')


#  get the non-overlapping counts of all RBPs collected for each transcript
#  isolate the genes of interest
GENES<-c('ARID1A', 'SDHA')
g<-fimo.txs[ gene_name %in% GENES ][ order(gene_name) ]


#  SDHA
g[ gene_name %in% GENES[2], sort(unique(unlist(rbp_group)))]
#
#  AGO2, CELF1, ELAVL1, ELAVL2, ELAVL3, ELAVL4, FMR1, HNRNPA1, HNRNPA2B1, HNRNPD, HNRNPF, HNRNPH1, HNRNPH2, HNRNPH3, PCBP1, PTBP1, RBMX, SFPQ, SNRPA, 
#  SRSF1, SRSF2, TARDBP, TIA1, TIAL1, TRA2B, ZFP36, ZFP36L2


#  ARID1A
g[ gene_name %in% GENES[1]][ order(-ncount) ]
g[ gene_name %in% GENES[1], sort(unique(unlist(rbp_group)))]
#
#  AGO2, CELF1, CSTF2, ELAVL1, ELAVL2, ELAVL3, ELAVL4, HNRNPA1, HNRNPA2B1, HNRNPF, HNRNPH1, HNRNPH2, HNRNPH3, HNRNPK, HNRNPL, HNRNPLL, NELFE, PCBP1, 
#  PTBP1, QKI, RBM5, SFPQ, SRSF1, SRSF10, SRSF2, SRSF3, SRSF4, SRSF5, SRSF6, SRSF7, TARDBP, TIA1, TIAL1, TRA2A, TRA2B, YBX1, ZFP36


#  [KHSRP] get the counts on all transcripts ONLY for KHSRP-related motifs not shared among other RBPs (you need to use the fimo.txs.all object)
#          summarize them at the transcript level fist and then take the mean of counts and densities at the gene level
khsrp<-fimo.txs.all[ sapply(rbp_group, function(r){ all(grepl('KHSRP', r)) }) ][, .(ncount=.N, density=.N/length*1e3, length=length, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(transcript_id, gene_id, gene_name)][, .(ncount=mean(ncount), density=mean(density), length=mean(length), rbp_group=list(unique(unlist(rbp_group))), motif_alt_id=list(unique(unlist(motif_alt_id)))), by=.(gene_id, gene_name)][ order(-ncount) ]
khsrp[ ncount>2 ]


#  [ELAV] similar to KHSRP but ELAV-related 
elav<-fimo.txs.all[ sapply(rbp_group, function(r){ all(grepl('ELAV', r)) }) ][, .(ncount=.N, density=.N/length*1e3, length=length, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(transcript_id, gene_id, gene_name)][, .(ncount=mean(ncount), density=mean(density), length=mean(length), rbp_group=list(unique(unlist(rbp_group))), motif_alt_id=list(unique(unlist(motif_alt_id)))), by=.(gene_id, gene_name)][ order(-ncount) ]
elav[ ncount>2 ]

#}}}



#  [ELAV motifs on circARID1A] special analysis with no controls
#{{{
rm(list=ls())
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(data.table)


#  load FIMO results for the circRNA sequences and isolate the circARID1A ones
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/results.RData')
fimo<-fimo[ grep('ARID1A', circ_name) ]
fimo.all<-fimo.all[ grep('ARID1A', circ_name) ]
rm(list=setdiff(ls(), c(l, 'fimo', 'fimo.all')))


#  load ATtRACt db and isolate all ELAV motifs 6mers and above
l<-ls()
load('/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')
motifs<-DNAStringSet(db[ sapply(gene_name, function(g){ any(grepl('ELAV', g)) }), unique(unlist(motif)) ])
motifs<-motifs[ width(motifs)>=6 ]
rm(list=setdiff(ls(), c(l, 'motifs')))


#  load circARDI1A sequence
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')
arid1a<-CIRCS.exons.seqs[ grep('ARID1A', names(CIRCS.exons.seqs)) ][[1]]
rm(list=setdiff(ls(), c(l, 'arid1a')))


#  shuffle circARID1A sequence preserving dinucleotide frequencies
writeXStringSet(DNAStringSet(arid1a), file='~/Downloads/arid1a.fa')
#
#    shuffle -d -o ~/Downloads/arid1a_suffled.fa ~/Downloads/arid1a.fa 
#
arid1a.control<-readDNAStringSet('~/Downloads/arid1a_suffled.fa')[[1]]
oligonucleotideFrequency(arid1a, width=2, step=1)
oligonucleotideFrequency(arid1a.control, width=2, step=1)


#  look at all hits with up to one mismatch on circARID1A and control
hits<-hits.control<-List()
for(m in seq_along(motifs)){
    h<-matchPattern(motifs[[m]], arid1a, max.mismatch=1)
    if(length(h)>0){ hits[[toString(motifs[[m]])]]<-h }
    h<-matchPattern(motifs[[m]], arid1a.control, max.mismatch=1)
    if(length(h)>0){ hits.control[[toString(motifs[[m]])]]<-h }
}
hits<-unlist(hits)
hits.control<-unlist(hits.control)


#  are there any motifs not found in control?
length(d<-setdiff(hits@ranges@NAMES, hits.control@ranges@NAMES))    #  4
length(i<-intersect(hits@ranges@NAMES, hits.control@ranges@NAMES))  #  1
hits[ hits@ranges@NAMES %in% d ]
#      start end width
#  [1]   681 687     7 [TTCATTT]
#  [2]   681 687     7 [TTCATTT]
#  [3]   685 691     7 [TTTGGGT]
#  [4]   198 203     6 [CTCTTA]
#  [5]   200 205     6 [CTTATA]

#}}}



