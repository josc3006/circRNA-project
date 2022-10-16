###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




############################################################
#
#
#  collect results and do some basic statistics and plotting
#
#
############################################################




#  [featureCounts] collect results for genes and exons 
#                  unannotated but called exons are removed including exons on unplaced contigs
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)


#  functions
#{{{

parse_gene_counts<-function(B=''){
    require(data.table)


    #  load counts skipping first row which is a comment
    cn<-fread(B, sep='\t', skip=1, col.names=c('gene_id', 'chr', 'start', 'end', 'strand', 'length', 'counts'))[, c('chr', 'start', 'end', 'strand'):=list(NULL, NULL, NULL, NULL)]


    return(as.data.frame(cn))
}


parse_exon_counts<-function(B=''){
    require(data.table)


    #  load counts skipping first row which is a comment
    cn<-fread(B, sep='\t', skip=1, col.names=c('gene_id', 'chr', 'start', 'end', 'strand', 'length', 'counts'))


    #  remove identical exon entries
    cn<-cn[, .(counts=counts[1]), by=.(gene_id, chr, start, end, strand, length)]


    #  load junction counts
    jn<-fread(sub('$', '.jcounts', B), sep='\t', skip=0, header=T, col.names=c('gene_id_d', 'gene_id_a', 'chr_d', 'start_d', 'strand_d', 'chr_a', 'start_a', 'strand_a', 'count'))


    #  remove exons on unannotated features or exons on unplaced contigs
    jn<-jn[ !is.na(gene_id_d) ]


    #  convert the acceptors gene_id column to list 
    jn[, gene_id_a:=strsplit(gene_id_a, ',')]  


    return(list(exons=cn, junctions=jn))
}

#}}}


#  locate the gene counts
cnts<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/ -maxdepth 3 -type f -wholename \'*/counts/genes.tsv\' -print', stdout=T)


#  add ids
names(cnts)<-sub('^.*raw/(CB[^/]+)/counts/.*$', '\\1', cnts)


#  order them
cnts<-cnts[ order(names(cnts)) ]


#  collect the gene counts
totalrna<-List()
for (n in seq_along(cnts)){
    cat('\nprocessing: ', cnts[n], '\n')
    totalrna[[ names(cnts)[n] ]]<-parse_gene_counts(cnts[n])
}
rm(n)


#  save
save(totalrna, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')


#  collect the exon counts including the exon junctions
cnts<-sub('genes\\.tsv$', 'exons.tsv', cnts)
exns<-List()
for (n in seq_along(cnts)){
    cat('\nprocessing: ', cnts[n], '\n')
    exns[[ names(cnts)[n] ]]<-parse_exon_counts(cnts[n])
}
rm(n)


#  save
save(exns, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts_exons.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts_exons.RData



#  [kallisto] collect results 
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)


#  functions
#{{{

parse_results<-function(B=''){
    require(data.table)


    #  load counts 
    gn<-fread(B, sep='\t', header=T, col.names=c('transcript_id', 'length', 'effective_length', 'counts', 'tpm'))


    return(as.data.frame(gn))
}

#}}}


#  locate the counts 
cnts<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/ -maxdepth 3 -type f -wholename \'*/kallisto/abundance.tsv\' -print', stdout=T)


#  add ids
names(cnts)<-sub('^.*raw/(CB[^/]+)/kallisto/.*$', '\\1', cnts)


#  order them
cnts<-cnts[ order(names(cnts)) ]


#  append results
totalrna<-List()
for (n in seq_along(cnts)){
    cat('\nprocessing: ', cnts[n], '\n')
    totalrna[[ names(cnts)[n] ]]<-parse_results(cnts[n])
}
rm(n)


#  save
save(totalrna, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/kallisto.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/kallisto.RData



#  [CIRI2] collect results
#          transgenic circRNAs are split to a separate list
#          circRNAs on alternative contigs, chrM, and chrY are removed
#          circRNAs with at least 5 reads covering the junction in at least one sample are kept
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(openxlsx)


#  functions
#{{{

parse_results<-function(B=''){
    require(data.table)

    n<-tryCatch({
        cr<-fread(B, sep='\t', select=c(2:5, 7:11), col.names=c('seqnames', 'start', 'end', 'jc_count', 'non_jc_count', 'ratio', 'region', 'gene_id', 'strand'))

        #  remove last useless comma from gene_id
        cr$gene_id<-sub(',$', '', cr$gene_id)


        #  convert 'n/a' to empty string
        cr$gene_id<-sub('n/a', '', cr$gene_id)


        #  convert whole column to list
        cr[, gene_id:=strsplit(gene_id, ',', fixed=T)] 

        
        #  convert all character(0) that strsplit() returns when no comma is found to NA
        cr[ lengths(gene_id)==0, gene_id:=lapply(gene_id, function(g){ NA }) ]

        
        #  add gene_name column
        cr[, gene_name:=relist(hsa$gene_name[ match(unlist(cr$gene_id), hsa$gene_id) ], cr$gene_id)]

        
        #  convert to GRanges()
        cr<-GRanges(as.data.frame(cr))


        #  remove circRNAs on alternative contigs
        seqlevels(cr, pruning.mode='coarse')<-seqlevels(cr)[ grep('chr', seqlevels(cr)) ]


        #  remove circRNAs on chrM, chrY
        seqlevels(cr, pruning.mode='coarse')<-seqlevels(cr)[ grep('chrM|chrY', seqlevels(cr), invert=T) ]


        #  split transgenic circRNAs to a seprate list
        tr<-cr[ lengths(cr$gene_id)>1 ]
        cr<-cr[ lengths(cr$gene_id)==1 ]
        cr$gene_id<-unlist(cr$gene_id)
        cr$gene_name<-unlist(cr$gene_name)


        #  return GRanges object
        return(GRangesList(circs=cr, trans=tr))

    }, error=function(e){
        warning(e)
        return(GRangesList(circs=GRanges(), trans=GRanges()))
    }) 
}

#}}}


#  load the reference to add gene_name to gene_id
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
hsa<-mcols(hsa)[, c('gene_name', 'gene_id')]


#  locate the CIRI2 results from the non-failed samples
cir<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/ -maxdepth 3 -type f -wholename \'*/ciri/circRNAs.tsv\' -print', stdout=T)


#  add ids
names(cir)<-sub('^.*raw/(CB[^/]+)/ciri/.*$', '\\1', cir)


#  order them
cir<-cir[ order(names(cir)) ]


#  append results 
circ<-List()
trans<-List()
for (n in seq_along(cir)){
    cat('\nprocessing: ', cir[n], '\n')
    x<-parse_results(cir[n])
    circ[[ names(cir)[n] ]]<-x$circ
    trans[[ names(cir)[n] ]]<-x$trans
}
rm(n,x)


#  convert from List to GRanges and add bid to metadata
circ<-unlist(GRangesList(lapply(circ,c)))
circ$bid<-names(circ)
names(circ)<-NULL
trans<-unlist(GRangesList(lapply(trans,c)))
trans$bid<-names(trans)
names(trans)<-NULL


#  keep circRNAs with at least 5 reads covering the junction in at least one sample
circ<-data.table(as.data.frame(circ))
circ[, pass:=any(jc_count>=5), by=.(seqnames, start, end, strand)] 
circ<-circ[ pass %in% TRUE, ][, pass:=NULL]
circ<-GRanges(seqnames=circ$seqnames, strand=circ$strand, ranges=IRanges(start=circ$start, end=circ$end), data.frame(circ[, -c(1:5)]))
stopifnot( all(lengths(circ$gene_id)==1) )
stopifnot( all(lengths(circ$gene_name)==1) )
circ$gene_id<-unlist(circ$gene_id)
circ$gene_name<-unlist(circ$gene_name)
trans<-data.table(as.data.frame(trans))
trans[, pass:=any(jc_count>=5), by=.(seqnames, start, end, strand)] 
trans<-trans[ pass %in% TRUE, ][, pass:=NULL]
trans<-GRanges(seqnames=trans$seqnames, strand=trans$strand, ranges=IRanges(start=trans$start, end=trans$end), data.frame(trans[, -c(1:5)]))


#  save
save(circ, trans, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2.RData')
l<-circ
x<-cbind(data.frame(position=paste0(as.character(seqnames(l)), '(', as.character(strand(l)),'):', start(l), '-', end(l))), as.data.frame(l)[, -c(1:3,5)])
write.xlsx(x, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2.xlsx', col.names=T, row.names=F, sheetName='circRNAs', append=F)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2.{RData,xlsx}



#  [CIRI2] summarize common circRNA predictions across all samples 
#          plot number of circRNAs per gene 
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  functions
#{{{

granges2kable<-function(gr){
    circ_pos<-paste0('[', as.character(seqnames(gr)), '(', as.character(strand(gr)),'):', start(gr), '-', end(gr),
                     '](http://www.ensembl.org/Homo_sapiens/Location/View?r=', 
                     sub('^chr', '', as.character(seqnames(gr))), ':', start(gr), '-', end(gr), 
                     ')')
    gene_name<-paste0('[', gr$gene_name, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', gr$gene_name, '&keywords=', gr$gene_name, ')')
    gene_name[ grep('\\[NA\\]', gene_name) ]<-'not_annotated'
    gr.table<-data.table(  
        locus=circ_pos,
        #strand=paste0('\\', as.character(strand(gr))),
        width=width(gr),
        gene_name=gene_name,
        jc_count=gr$jc_count,
        non_jc_count=gr$non_jc_count,
        ratio=gr$ratio,
        region=gr$region,
        bid=gr$bid
    )

    #  unlist all lists and convert to comma-separated strings, make sure to add space so kable and line-break them
    gr.table[, ':='(locus=locus, 
        width=width, 
        gene_name=gene_name, 
        jc_count=sapply(jc_count, paste, sep='', collapse=', '), 
        non_jc_count=sapply(non_jc_count, paste, sep='', collapse=', '), 
        ratio=sapply(ratio, paste, sep='', collapse=', '), 
        region=region,
        bid=sapply(bid, paste, sep='', collapse=', '))]

    return(as.data.frame(gr.table))
}

#}}}


#  load collected results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2.RData')


#  remove failed samples and Pilot runs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
circ<-circ[ circ$bid %in% meta[!(failed) & !grepl('CBPilot', bid) , bid ] ]
rm(meta, meta.tum, meta.cel, meta.prefailed, trans)


#  summarize duplicate calls across samples
l<-circ
l<-data.table(as.data.frame(l))[, .(jc_count=list(jc_count), non_jc_count=list(non_jc_count), ratio=list(ratio), bid=list(bid)), by=.(seqnames, start, end, strand, gene_name, gene_id, region)]
l<-GRanges(as.data.frame(l))


#  create HTML table with proper links
l.table<-granges2kable(l)


#  save as R object
save(l, l.table, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2_summarized.RData')


#  save as Excel sheet after converting list columns to character vectors
x<-data.table(cbind(data.frame(position=paste0(as.character(seqnames(l)), '(', as.character(strand(l)),'):', start(l), '-', end(l))), as.data.frame(l)[, -c(1:3,5)]))
x[, ':='(jc_count=sapply(jc_count, paste, sep='', collapse=','), 
         non_jc_count=sapply(non_jc_count, paste, sep='', collapse=','), 
         ratio=sapply(ratio, paste, sep='', collapse=','), 
         bid=sapply(bid, paste, sep='', collapse=','))]
write.xlsx(as.data.frame(x), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2_summarized.xlsx', col.names=T, row.names=F, sheetName='summarized circRNAs', append=F)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2_summarized.{RData,xlsx}



#  [CIRI2] identify all unique circRNAs 
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(openxlsx)


#  functions
#{{{

granges2kable<-function(gr){
    circ_pos<-paste0('[', as.character(seqnames(gr)), '(', as.character(strand(gr)),'):', start(gr), '-', end(gr),
                     '](http://www.ensembl.org/Homo_sapiens/Location/View?r=', 
                     sub('^chr', '', as.character(seqnames(gr))), ':', start(gr), '-', end(gr), 
                     ')')
    gene_name<-paste0('[', gr$gene_name, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', gr$gene_name, '&keywords=', gr$gene_name, ')')
    gene_name[ grep('\\[NA\\]', gene_name) ]<-'not_annotated'
    gr.table<-data.table(  
        locus=circ_pos,
        #strand=paste0('\\', as.character(strand(gr))),
        width=width(gr),
        gene_name=gene_name,
        jc_count=gr$jc_count,
        non_jc_count=gr$non_jc_count,
        ratio=gr$ratio,
        region=gr$region,
        bid=gr$bid
    )

    #  unlist all lists and convert to comma-separated strings, make sure to add space so kable and line-break them
    gr.table[, ':='(locus=locus, 
        width=width, 
        gene_name=gene_name, 
        jc_count=sapply(jc_count, paste, sep='', collapse=', '), 
        non_jc_count=sapply(non_jc_count, paste, sep='', collapse=', '), 
        ratio=sapply(ratio, paste, sep='', collapse=', '), 
        region=region,
        bid=sapply(bid, paste, sep='', collapse=', '))]

    return(as.data.frame(gr.table))
}

#}}}


#  load collected results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2.RData')


#  remove failed samples and Pilot samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
circ<-circ[ circ$bid %in% meta[!(failed) & !grepl('CBPilot', bid) , bid ] ]
rm(meta, meta.tum, meta.cel, meta.prefailed, trans)


#  unique circRNAs
#{{{

#  sort by expression
circ.u<-unique(circ)
circ.u<-circ.u[ order(circ.u$jc_count, decreasing=T) ]  


#  create HTML table with proper links
l<-circ.u
l.table<-granges2kable(l)


#  remove bid since it does not make sense to keep
l$bid<-NULL
l.table$bid<-NULL


#  statistics
table(l$region)
# 
#              exon intergenic_region            intron 
#             33459              2138              4244 
length(l)   #  39841


#  save 
save(l, l.table, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2_unique.RData')
x<-cbind(data.frame(position=paste0(as.character(seqnames(l)), '(', as.character(strand(l)),'):', start(l), '-', end(l))), as.data.frame(l)[, -c(1:3,5)])
write.xlsx(x, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2_unique.xlsx', col.names=T, row.names=F, sheetName='unique circRNAs', append=F)
rm(x,l,l.table)

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2_unique.{RData,xlsx}



#  [HLA typing] collect results 
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)


#  functions
#{{{

parse_results<-function(B=''){
    require(data.table)

    #  load results
    hl<-as.data.frame(fread(B, sep='\t'))
    rownames(hl)<-hl[, 1]
    hl<-hl[, -1]

    return(hl)
}

#}}}


#  locate the gene estimated TPMs and counts, the transcripts will be looked for later on
hla<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/ -maxdepth 4 -type f -wholename \'*/optitype/*/*result.tsv\' -print', stdout=T)


#  add ids
names(hla)<-sub('^.*raw/(CB[^/]+)/optitype/.*$', '\\1', hla)


#  order them
hla<-hla[ order(names(hla)) ]


#  append results
hla.t<-List()
for (n in seq_along(hla)){
    cat('\nprocessing: ', hla[n], '\n')
    hla.t[[ names(hla)[n] ]]<-parse_results(hla[n])
}
hla<-hla.t
rm(n,hla.t)


#  save
save(hla, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/hla_typing.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/hla_typing.RData



#  [HLA typing] by eye informatics
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  functions
#{{{

scatterplots<-function(s1='', s2='', figs.dir=NULL){
    #  scatterplot based on counts and TPMs of commonly expressed genes between sample (s1) and (s2)


    #  isolate the samples of interest
    X<-totalrna[[ s1 ]]
    Y<-totalrna[[ s2 ]]


    #  remove mutually unexpressed genes
    keep<- X$counts!=0 & Y$counts!=0
    cat('\nTotal number of genes = ', nrow(X), ' of which ', sum(!keep), ' (', round(100*sum(!keep)/nrow(X), 1), '%) are mutually non-expressed and will be removed.\n\n', sep='') 
    X<-X[keep, ]
    Y<-Y[keep, ]
    rm(keep)


    #  add FPKM and TPM columns
    X$fpkm<-1e9*X$counts/X$length/sum(X$counts)
    Y$fpkm<-1e9*Y$counts/Y$length/sum(Y$counts)
    X$tpm<-1e6*X$fpkm/sum(X$fpkm)
    Y$tpm<-1e6*Y$fpkm/sum(Y$fpkm)


    #  isolate log10(1+counts)
    X.c<-setNames(X$counts, X$gene_id)
    Y.c<-setNames(Y$counts, Y$gene_id)
    X.c<-log10(1+X.c)
    Y.c<-log10(1+Y.c)


    #  isolate log10(1+TPMs)
    X.t<-setNames(X$tpm, X$gene_id)
    Y.t<-setNames(Y$tpm, Y$gene_id)
    X.t<-log10(1+X.t)
    Y.t<-log10(1+Y.t)


    #  linear regression
    l.c<-lm( Y.c ~ X.c )
    l.t<-lm( Y.t ~ X.t )


    #  [counts] scatterplot
    x11(width=11, height=11, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel')
    MAX<-ceiling(max(X.c, Y.c))
    par(mar=c(4.5,5.0,1.0,1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.8, cex.axis=1.8)
    den<-col2rgb(densCols(X.c, Y.c, colramp=colorRampPalette(c('black', 'white'))))[1, ]+1L
    cls<-colorRampPalette(c('#000099', '#00FEFF', '#45FE4F','#FCFF00', '#FF9400', '#FF3100'))(256)[den]
    plot(X.c[ order(den) ], Y.c[order(den)], type='p', pch=20, col=cls[order(den)], main='', xlab='', ylab='', cex=1.2, xlim=c(0, MAX), ylim=c(0, MAX))
    abline(l.c, lty=1, lwd=4, xpd=F)
    mtext(bquote(R^2 == .(format(summary(l.c)$r.squared, digits=2))), side=3, line=-1, padj=+0.5, cex=1.4, xpd=NA)
    mtext(bquote(log[10](1+counts) ~ group("[", .(s1), "]")), side=1, line=2, padj=+0.8, cex=1.8)
    mtext(bquote(log[10](1+counts) ~ group("[", .(s2), "]")), side=2, line=3, padj=+0.2, cex=1.8, las=0)
    if (!is.null(figs.dir)){
        #dev.print(device=pdf, file=paste0(figs.dir, '/scatterplot_counts_', s1, '_', s2, '.pdf'), width=16, height=16, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
        dev.print(device=svg, file=paste0(figs.dir, '/scatterplot_counts_', s1, '_', s2, '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    }

    #  [TPMs] scatterplot
    x11(width=11, height=11, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel')
    MAX<-ceiling(max(X.t, Y.t))
    par(mar=c(4.5,5.0,1.0,1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.8, cex.axis=1.8)
    den<-col2rgb(densCols(X.t, Y.t, colramp=colorRampPalette(c('black', 'white'))))[1, ]+1L
    cls<-colorRampPalette(c('#000099', '#00FEFF', '#45FE4F','#FCFF00', '#FF9400', '#FF3100'))(256)[den]
    plot(X.t[ order(den) ], Y.t[order(den)], type='p', pch=20, col=cls[order(den)], main='', xlab='', ylab='', cex=1.2, xlim=c(0, MAX), ylim=c(0, MAX))
    abline(l.t, lty=1, lwd=4, xpd=F)
    mtext(bquote(R^2 == .(format(summary(l.t)$r.squared, digits=2))), side=3, line=-1, padj=+0.5, cex=1.4, xpd=NA)
    mtext(bquote(log[10](1+TPM) ~ group("[", .(s1), "]")), side=1, line=2, padj=+0.8, cex=1.8)
    mtext(bquote(log[10](1+TPM) ~ group("[", .(s2), "]")), side=2, line=3, padj=+0.2, cex=1.8, las=0)
    if (!is.null(figs.dir)){
        #dev.print(device=pdf, file=paste0(figs.dir, '/scatterplot_tpms_', s1, '_', s2, '.pdf'), width=16, height=16, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
        dev.print(device=svg, file=paste0(figs.dir, '/scatterplot_tpms_', s1, '_', s2, '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    }
}

#}}}


#  load collected HLA typing results and keep best solution
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/hla_typing.RData')
h<-do.call(rbind, lapply(hla, function(x){ x[1, ] }))[, 1:6]
h$bid<-rownames(h)
h<-data.table(h)[ order(A1, A2, B1, B2, C1, C2) ]
h[, hla:=paste(A1, A2, B1, B2, C1, C2, sep='_')]


#  load sample metadata
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)


#  load count data
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')


#  indicate the identical HLA types
h[hla %in% names(table(hla)[ table(hla)>1]), ]


#  CB2019-11-R01, CB3037-11-R01 have identical HLA typing and are strongly correlated
scatterplots('CB2019-11-R01', 'CB3037-11-R01', figs.dir='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/')
meta[ bid %in% c('CB2019-11-R01', 'CB3037-11-R01') ]

#}}}



#  [sex determination] XIST, UTY, RPS4Y1 expressions 
#{{{
rm(list=ls())  
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(gplots)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  functions
#{{{

genes_barplot<-function(G=c(''), FE, names.trim='', figure.prefix='', PARAMS=NULL){
    source('~/bio/lib/trim_text.R')


    #  identify gene_id for given gene_name
    gid<-data.frame(mcols(hsa[ hsa$gene_name %in% G])[, c('gene_id', 'gene_name')])


    #  isolate genes of interest 
    #  add gene_name 
    #  sort by gene_name and bid
    #  compute log10(1+TPM)
    x<-FE[ gene_id %in% hsa$gene_id[ hsa$gene_name %in% G ] ]
    x$gene_name<-gid$gene_name[ match(x$gene_id, gid$gene_id) ]
    x<-x[, c('bid', 'gene_name', 'tpm')][ order(gene_name, bid) ]
    x$tpm<-log10(1+x$tpm)


    #  convert from:
    #
    #      BID , GENE_NAME , TPM  
    #  
    #  to:
    #
    #      BID , GENE_NAME_TPM      
    #
    #  for any number of genes provided!
    x<-dcast(x, bid ~ gene_name, value.var='tpm')
    y<-as.matrix(x[, -1])
    rownames(y)<-x$bid
    x<-y
    rm(y)

    
    #  define a color for each gene
    cl<-setNames(colorRampPalette(brewer.pal(8,'Dark2'))(length(G)), G)


    if(names.trim!=''){
        original.names<-rownames(x)
        rownames(x)<-sub(names.trim, '', rownames(x))
    }


    #  stacked barplot 
    YMAX<-ceiling(max(rowSums(x)))
    YTICK<-pretty(c(0, YMAX), 5)
    if(!is.null(PARAMS)){
        par(PARAMS$par)
    } else {
        par(mar=c(10.0, 8.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
    }
    bp<-barplot(t(x), beside=F, plot=F)
    plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(bp)+c(-1,2), axes=F, ann=F, xaxs='i', yaxs='i')
    bp<-barplot(t(x), border='white', col=cl[colnames(x)], axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=range(YTICK), add=T)
    axis(2, at=YTICK, line=0, cex.axis=2.4)
    mtext(text=rownames(x), side=1, line=0, at=bp, las=2, adj=1, cex=1.8)
    mtext(expression(log[10](1+TPM)), side=2, line=4, padj=+0.3, las=0, cex=2.4)
    legend(x=-par('usr')[2]*0.02, y=1.03*par('usr')[4], legend=names(cl) , col=cl, bty='n', lty=1, lwd=18, cex=2.0, y.intersp=0.60, x.intersp=0.1, seg.len=0.1)
    

    #  save only if directory is given
    if(nchar(figure.prefix)>0){
        dev.print(device=svg, file=paste0(figure.prefix, paste0(G, collapse='_'), '.svg'), width=60, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    }


    #  replace original names
    if(names.trim!=''){
        rownames(x)<-original.names
    }

    return(x)
}

#}}}


#  load annotation to identify gene_id from transcript_id
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_name', 'gene_type', 'level')]


#  load featureCount data 
#  unlist them to data.table
#  add TPMs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')
fe<-unlist(totalrna)
fe$bid<-sub('\\.[0-9]*$', '', rownames(fe))
rownames(fe)<-NULL
fe<-data.table(fe[, c('bid', 'gene_id', 'length', 'counts')])
fe<-fe[, tpm:=1e6*counts/length/sum(counts/length), by=.(bid)]
rm(totalrna)


#  load metadata and remove failed and Pilote samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta.tum<-meta.tum[ !(failed) ]
meta.cel<-meta.cel[ !(failed) & !grepl('CBPilot', bid) ]
fe.tum<-fe[ bid %in% meta.tum$bid ]
fe.cel<-fe[ bid %in% meta.cel$bid ]
rm(meta.cel, meta.prefailed, fe)


#  [tumors] stacked barplot
x11(width=22, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 
s<-genes_barplot(c('XIST', 'UTY', 'RPS4Y1'), FE=fe.tum, names.trim='-11-R01', figure.prefix='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_tumors_', PARAMS=list(par=list(mar=c(6.0, 8.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)))


#  save sex determination to add to metadata if necessary
s<-data.table(s, bid=rownames(s))
s[, sex:=if( RPS4Y1>XIST | UTY>XIST ){ 'M' }, by=.(bid)]
s[, sex:=if( RPS4Y1<XIST & UTY<XIST ){ 'F' }, by=.(bid)]
sex<-s
rm(s)
save(sex, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/sex_determination.RData')


#  does the determination agree with the clinical metadata?
x<-meta.tum[sex[, c('bid', 'sex'), with=F], on='bid'][!is.na(PAT_ID_BERLIN), c('bid', 'SEX', 'sex', 'risk_group')]
sex[ bid %in% x[ SEX != sex, bid ] ]  #  full agreement
table( x$risk_group, x$SEX )
#          
#            F  M
#   HR_nMNA  7 21
#   IMR      4  6
#   LR      11 19
#   MNA      8 14
#   ST4S     6  7

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/sex_determination.RData



#  [STAR] collect mapping statistics for all samples
#{{{
rm(list=ls())
library(data.table)


#  locate STAR Log.final.out files 
logs<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/ -maxdepth 3 -wholename \'*/star/Log.final.out\' -print', stdout=T)


#  populate mapping summary data.frame
map.s<-data.frame(logs=logs, dir='', bid='', raw_reads=0, unimapped=0, multimapped=0, unmapped=0, alignments=0)
for (l in seq_along(logs)){
    r<-readLines(logs[l]) 
    map.s$dir[l]<-dirname(dirname(logs[l]))
    map.s$bid[l]<-sub('^.*/(CB[^/]+)/.*$', '\\1',logs[l])
    map.s$raw_reads[l]<-as.integer(strsplit(r[ grep('Number of input reads', r) ], '\\t')[[1]][2])
    map.s$unimapped[l]<-as.integer(strsplit(r[ grep('Uniquely mapped', r) ], '\\t')[[1]][2])
    map.s$multimapped[l]<-as.integer(strsplit(r[ grep('Number of reads mapped to multiple loci', r) ], '\\t')[[1]][2])
    map.s$unmapped[l]<-map.s$raw_reads[l]-map.s$unimapped[l]-map.s$multimapped[l]
    map.s$alignments[l]<-as.integer(readLines(sub('Log.final.out', 'Aligned.sortedByCoord.out.bam.alignments', logs[l])))
}
rm(l,r,logs)


#  order them
setorder(map.s, bid)
map.s<-map.s[, c('dir', 'bid', 'raw_reads', 'unmapped', 'unimapped', 'multimapped', 'alignments')]
rownames(map.s)<-NULL


#  save
save(map.s, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/mapping_summary.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/mapping_summary.RData



#  [STAR] barplots including failed samples
#{{{
rm(list=ls())
library(data.table)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/trim_text.R')


#  load mapping summary and metadata
#  order libraries according to metadata order
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/mapping_summary.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
stopifnot( length(setdiff( map.s$bid , meta$bid ))==0 )  #  all sequenced libraries show up in the metadata
meta<-meta[ bid %in% map.s$bid ]
map.s<-map.s[ match(meta$bid, map.s$bid), ]


#  remove Pilot study samples
meta<-meta[ grep('CBPilote', bid, invert=T) ]
map.s<-map.s[ grep('CBPilote', map.s$bid, invert=T), ]


#  colors related to alignments and samples
cl.l<-data.frame(symbol=c('unimapped', 'multimapped','unmapped'), 
                  color=c('#228B22',  #  forestgreen
                          '#1874CD',  #  dodgerblue3
                          '#B22222')) #  firebrick


#  CLICK on it once to make sure it does not redraw, or options(scipen=-20) might be IGNORED
x11(width=40, height=18, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 


#  [tumors] stacked bars
par(mar=c(6.0, 12.0, 2.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
n<-map.s$bid %in% meta[ !is.na(risk_group), bid]
x<-map.s[n, c('unimapped', 'multimapped', 'unmapped')]
rownames(x)<-sub('-11-R01$', '', map.s[n, 'bid'])
cl.s<-meta[ match(rownames(x), PAT_ID_BERLIN), c('risk_group', 'col') ] 
cl.s<-setNames(cl.s$risk_group, cl.s$col)
YMAX<-max( rowSums(x) - rowSums(x)%%1e5 + 1e5 )
YTICK<-pretty(c(0, YMAX), 5)
YMAX<-tail(YTICK, 1)
bp<-barplot(t(x), plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(bp)+c(-1, +1), axes=F, ann=F, xaxs='i', yaxs='i') 
bp<-barplot(t(x), border='white', col=cl.l[match(colnames(x), cl.l$symbol), 2], axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, YMAX), add=T)
options(scipen=-20)
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext(text=rownames(x), side=1, line=0, at=bp, col=names(cl.s), las=2, adj=1, cex=1.8)
mtext('Number of raw reads', side=2, line=9, padj=+0.1, las=0, cex=2.4)
legend(x=par('usr')[2]*0.77, y=1.05*par('usr')[4], legend=cl.l$symbol[match(colnames(x), cl.l$symbol)] , col=cl.l$color[match(colnames(x), cl.l$symbol)], bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.40, x.intersp=0.2, seg.len=0.2)
legend(x=par('usr')[1]*0.5, y=1.05*par('usr')[4], legend=unique(cl.s) , col=unique(names(cl.s)),  bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.40, x.intersp=0.2, seg.len=0.2)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_mapping_statistics_alignments_tumors.svg', width=60, height=18, bg='white', antialias='subpixel', pointsize=20, family='Arial')
options(scipen=0)


#  [cell models] stacked bars
par(mar=c(22.0, 12.0, 2.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
n<-map.s$bid %in% meta[ !is.na(cell_model), bid]
x<-map.s[n, c('unimapped', 'multimapped', 'unmapped')]
rownames(x)<-map.s[n, 'bid']
cl.s<-meta[ match(rownames(x), bid), c('cell_model', 'treatment', 'col') ] 
cl.s<-setNames(paste(cl.s$cell_model, cl.s$treatment), cl.s$col)
YMAX<-max( rowSums(x) - rowSums(x)%%1e5 + 1e5 )
YTICK<-pretty(c(0, YMAX), 5)
YMAX<-tail(YTICK, 1)
bp<-barplot(t(x), plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(bp)+c(-1, +1), axes=F, ann=F, xaxs='i', yaxs='i') 
bp<-barplot(t(x), border='white', col=cl.l[match(colnames(x), cl.l$symbol), 2], axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, YMAX), add=T)
options(scipen=-20)
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext(text=sub('^CB-', '', sub('-R01$', '', rownames(x))), side=1, line=0, at=bp, col=names(cl.s), las=2, adj=1, cex=1.8)
mtext('Number of raw reads', side=2, line=9, padj=+0.1, las=0, cex=2.4)
legend(x=par('usr')[2]*0.77, y=1.05*par('usr')[4], legend=cl.l$symbol[match(colnames(x), cl.l$symbol)] , col=cl.l$color[match(colnames(x), cl.l$symbol)], bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.20, x.intersp=0.2, seg.len=0.2)
#legend(x=par('usr')[1]*0.5, y=1.05*par('usr')[4], legend=unique(cl.s) , col=unique(names(cl.s)),  bty='n', lty=1, lwd=15, cex=1.4, y.intersp=0.2, x.intersp=0.2, seg.len=0.2)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_mapping_statistics_alignments_cell_models.svg', width=60, height=18, bg='white', antialias='subpixel', pointsize=20, family='Arial')
options(scipen=0)


#  clean up
dev.off()

#}}}



#  boxplot per risk group of number of raw reads
#  boxplot per risk group of number of reads mapped/unmapped to features
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  keep only non-failed neuroblastoma tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta<-meta[ !(failed) & !is.na(risk_group) ]


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  boxplot per risk_group of number of raw reads
B<-split(meta$nreads/1e6, factor(meta$risk_group, levels=unique(meta$risk_group)))
B.cl<-setNames(unique(meta$col), unique(meta$risk_group))
par(mar=c(6.5, 7.5,0.1,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.8, cex.axis=1.8)
YTICK<-pretty(c(min(floor(sapply(B, min))), max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of raw reads (M)', side=2, line=5, padj=-0.4, las=0, cex=1.8)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, adj=1.0, cex=1.8, col=B.cl)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, adj=1.0, cex=1.8, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/boxplot_number_of_raw_reads_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  boxplot per risk_group of number of reads unmapped
B<-split(meta$unmapped/1e6, factor(meta$risk_group, levels=unique(meta$risk_group)))
B.cl<-setNames(unique(meta$col), unique(meta$risk_group))
par(mar=c(6.5, 7.5,0.1,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.8, cex.axis=1.8)
YTICK<-pretty(c(min(floor(sapply(B, min))), max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of unmapped reads (M)', side=2, line=5, padj=-0.4, las=0, cex=1.8)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, adj=1.0, cex=1.8, col=B.cl)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, adj=1.0, cex=1.8, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/boxplot_number_of_unmapped_reads_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  boxplot per risk_group of percentage of unmapped reads
B<-split(meta$p_unmapped, factor(meta$risk_group, levels=unique(meta$risk_group)))
B.cl<-setNames(unique(meta$col), unique(meta$risk_group))
par(mar=c(6.5, 6.5,0.1,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.8, cex.axis=1.8)
YTICK<-pretty(c(0, 100), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Percent of unmapped reads', side=2, line=4, padj=-0.4, las=0, cex=1.8)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, adj=1.0, cex=1.8, col=B.cl)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, adj=1.0, cex=1.8, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/boxplot_percent_of_unmapped_reads_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  boxplot per risk_group of percentage of mapped reads
B<-split(100-meta$p_unmapped, factor(meta$risk_group, levels=unique(meta$risk_group)))
B.cl<-setNames(unique(meta$col), unique(meta$risk_group))
par(mar=c(6.5, 6.5,0.1,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.8, cex.axis=1.8)
YTICK<-pretty(c(0, 100), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Percent of mapped reads', side=2, line=4, padj=-0.4, las=0, cex=1.8)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, adj=1.0, cex=1.8, col=B.cl)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, adj=1.0, cex=1.8, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/boxplot_percent_of_mapped_reads_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


dev.off()

#}}}



#  [featureCounts] barplots including failed samples
#{{{
rm(list=ls())
library(data.table)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/trim_text.R')


#  load sample metadata
#  remove Pilote study samples and prefailed samples that were never sequenced
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta<-meta[ !is.na(features) & ! grepl('CBPilote', bid) ]
rm(meta.tum, meta.cel, meta.prefailed)


#  colors related to alignments and samples
cl.l<-data.frame(symbol=c('failed','features','no_features','unmapped'), 
                  color=c('#000000',  #  black
                          '#228B22',  #  forestgreen
                          '#1874CD',  #  dodgerblue3
                          '#B22222')) #  firebrick


#  [tumors] stacked bars
x11(width=40, height=18, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 
par(mar=c(6.0, 8.0, 2.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-data.frame(meta[!is.na(risk_group), c('features', 'no_features', 'unmapped')], row.names=sub('-11-R01$', '', meta[!is.na(risk_group), bid]))
x<-data.frame(t(apply(x, 1, function(x){ x/1e6 })))  #  1M units
#
summary(x$features)
#
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 11.21471 26.39510 30.44305 34.70366 38.01140 84.58114 
#
cl.s<-meta[ match(rownames(x), PAT_ID_BERLIN), c('risk_group', 'col') ] 
cl.s<-setNames(cl.s$risk_group, cl.s$col)
MARK<-which(x$features<20)  #  less than 20M reads mapped on features
YMAX<-max( rowSums(x) )
YTICK<-pretty(c(0, YMAX), 5)
YMAX<-tail(YTICK, 1)
bp<-barplot(t(x), width=1, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(bp)+c(-1, +1), axes=F, ann=F, xaxs='i', yaxs='i') 
bp<-barplot(t(x), width=1, border='white', col=cl.l[match(colnames(x), cl.l$symbol), 2], axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, YMAX), add=T)
segments(x0=bp[MARK]-0.5, x1=bp[MARK]-0.5, y0=0, y1=rowSums(x)[MARK], lty=1, lwd=10, col='black')
segments(x0=bp[MARK]-0.5, x1=bp[MARK]+0.5, y0=rowSums(x)[MARK], y1=rowSums(x)[MARK], lty=1, lwd=10, col='black')
segments(x0=bp[MARK]+0.5, x1=bp[MARK]+0.5, y0=rowSums(x)[MARK], y1=0, lty=1, lwd=10, col='black')
segments(x0=bp[MARK]+0.5, x1=bp[MARK]-0.5, y0=0, y1=0, lty=1, lwd=10, col='black')
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext(text=rownames(x), side=1, line=0, at=bp, col=names(cl.s), las=2, adj=1, cex=1.8)
mtext('Number of reads (M)', side=2, line=5, padj=-0.5, las=0, cex=2.4)
legend(x=par('usr')[2]*0.77, y=1.05*par('usr')[4], legend=cl.l$symbol[match(colnames(x), cl.l$symbol)] , col=cl.l$color[match(colnames(x), cl.l$symbol)], bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.40, x.intersp=0.2, seg.len=0.2)
legend(x=par('usr')[1]*0.5, y=1.05*par('usr')[4], legend=unique(cl.s) , col=unique(names(cl.s)),  bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.40, x.intersp=0.2, seg.len=0.2)
#mtext(paste0('Number of samples = ', nrow(x), ' (failed = ', length(MARK),')'), side=3, line=0, padj=-0.4, las=0, cex=1.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_featureCounts_assignments_tumors.svg', width=50, height=18, bg='white', antialias='subpixel', pointsize=20, family='Arial')
options(scipen=0)


dev.off()

#}}}



#  [junk] percentage of junk genes across samples
#{{{
rm(list=ls())
library(data.table)
library(rtracklayer)
library(GenomicAlignments)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/trim_text.R')


#  define junk once and for all
JUNK<-'^RN7S|^RP[PLS]|^MT-|^RNVU[0-9]|^RNU[0-9]|^RF[0-9]+|^FP[0-9]+|^A[CFLP][0-9]+|^MTND[12]P|^MTCO[123]'


#  load sample metadata
#  remove Pilote samples, unsequenced or metafailed
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta<-meta[ !(failed) & !is.na(features) & ! grepl('CBPilote', bid) ]


#  load reference 
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene']
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]


#  load gene counts 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')
stopifnot(length(setdiff( meta$bid, names(totalrna) ))==0)
gns<-unlist(totalrna[ meta$bid ])
gns$bid<-sub('\\.[0-9]*$', '', rownames(gns))
rownames(gns)<-NULL
gns<-data.table(gns[, c('bid', 'gene_id', 'length', 'counts')])
gns<-gns[, .(gene_id=gene_id, counts=counts, tpm=1e6*counts/length/sum(counts/length), cpm=1e6*counts/sum(counts)), by=.(bid)]
gns[, gene_name:=hsa$gene_name[ match(gene_id, hsa$gene_id) ] ]


#  create the CPM matrix
#  remove unexpressed genes
#  order by top expression
gns.cpm<-dcast(gns, bid ~ gene_name, value.var='cpm', fun.aggregate=sum)
gns.cpm<-t(data.frame(gns.cpm[, -1], row.names=gns.cpm[, bid], check.names=F))
gns.cpm<-gns.cpm[ rowMeans(gns.cpm)>0, , drop=F]
gns.cpm<-gns.cpm[ order(rowMeans(gns.cpm), decreasing=T), ]


#  percentage of reads going to junk
p<-round(gns.cpm, digits=2)/1e6
j<-p[ grep(JUNK, rownames(p)), ]
j<-j[ order(rowMeans(j), decreasing=T), ,drop=F]
nrow(j)  #  21855


#  boxplot for cell lines and tumors separately
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 
par(mar=c(1.5, 8.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.8, cex.axis=2.8)
B<-colSums(j)
names(B)<-meta[match(colnames(j), bid), ifelse(is.na(cell_model), 'tumor', 'cell_line')]
B<-split(B, factor(names(B), levels=c('tumor', 'cell_line')))
B.cl<-c('cornflowerblue', 'coral')
YTICK<-pretty(c(0, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(text=c('tumors', 'cell lines'), side=1, line=0, at=seq_along(B), las=0, adj=0.5, cex=2.8, col='black')
mtext('Junk fraction', side=2, line=5, padj=-0.5, las=0, cex=2.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/boxplot_junk_fraction.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}




#########################################
#
#
#  distribution of MYCN target expression
#
#
#########################################




#  [tumors, cell lines] empirical distribution of TPMs for MYCN induced and repressed targets
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load MYCN targets
load('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/MYCN_targets.RData')


#  load gene expression
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')


#  remove failed samples and samples not yet integrated into the cohort
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta<-meta[ !(failed) ]
totalrna<-totalrna[ names(totalrna) %in% meta$bid ]
meta<-meta[ bid %in% names(totalrna) ]
rm(meta.prefailed)


#  compute TPMs
fe<-unlist(totalrna)
fe$bid<-sub('\\.[0-9]*$', '', rownames(fe))
rownames(fe)<-NULL
fe<-data.table(fe[, c('bid', 'gene_id', 'length', 'counts')])
fe<-fe[, tpm:=1e6*counts/length/sum(counts/length), by=.(bid)]
fe<-dcast(fe, bid ~ gene_id, value.var='tpm')
rm(totalrna)


#  isolate MYCN targets
fe<-fe[, c('bid', MYCN$gene_id), with=F]


#  separate tumors and cell lines
cel<-fe[ bid %in% meta[ !is.na(cell_model), bid ], ]
pat<-fe[ bid %in% meta[ !is.na(risk_group), bid ], ]
stopifnot( nrow(fe)==nrow(cel)+nrow(pat) )
rm(fe)


#  separate metadata as well and order samples according to metadata
pat.meta<-meta[ match(pat$bid, bid), ]
cel.meta<-meta[ match(cel$bid, bid), ]
stopifnot( nrow(meta)==nrow(cel.meta)+nrow(pat.meta) )
pat<-pat[ match(pat.meta$bid, bid),  ]
cel<-cel[ match(cel.meta$bid, bid), ]
rm(meta)


#  compute mean TPM within risk groups for each gene
#  invert and convert to matrix for easy use after that
pat$level<-pat.meta$risk_group
pat<-pat[, -1][, lapply(.SD, mean), by=.(level)]
x<-as.matrix(pat[, -1])
rownames(x)<-pat[, level]
pat<-t(x)[, c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA')]
rm(x)


#  compute mean TPM within cell line AND treatment for each gene
#  invert and convert to matrix for easy use after that
cel$level<-paste(cel.meta$cell_model, cel.meta$treatment, sep=' ')
cel<-cel[, -1][, lapply(.SD, mean), by=.(level)]
x<-as.matrix(cel[, -1])
rownames(x)<-cel[, level]
cel<-t(x)
rm(x)


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 
par(mar=c(5.5, 7.0,0.5,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, cex.lab=2.4, cex.axis=2.4, xpd=F, bty='n')


#  [induced targets] cumulative distributions of log10(1+mean TPMs) across risk groups/cell types
#{{{

#  induced targets
id<-MYCN[ type %in% 'induced', gene_id ]


#  tumors
X<-apply(log10(1+pat[id, ]), 2, function(tpm){ ecdf(tpm) })  #  tumors
h<-curve(X[[1]](x), from=-0.01, to=3.0, n=20, ylab='', xlab='', pch=NA, col=pat.meta[ match(names(X)[1], risk_group), col], lty=1, lwd=10, main='', ylim=c(0, 1), xlim=c(0, 3), xaxt='n')
for(n in seq_along(X)[-1]){
    curve(X[[n]](x), from=head(h$x,1), to=tail(h$x,1), n=length(h$x), pch=NA, col=pat.meta[ match(names(X)[n], risk_group), col], lty=1, lwd=10, main='', xaxt='n', add=T)
}
axis(1, at=pretty(c(0, max(h$x))), labels=pretty(c(0, max(h$x))))
mtext('Probability', side=2, line=4, padj=-0.5, cex=2.4, las=0)
mtext(expression(log[10](1+'mean TPM')), side=1, line=4, las=0, padj=+0.1, cex=2.4)
legend('topleft', legend=names(X), col=pat.meta[ match(names(X), risk_group), col], bty='n', lty=1, lwd=15, cex=1.6, y.intersp=0.6, x.intersp=0.5, seg.len=0.5, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/ecdf_MYCN_induced_targets_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(h,X)


#  cell lines
X<-apply(log10(1+cel[id, ]), 2, function(tpm){ ecdf(tpm) })  #  cell lines
cl<-cel.meta[, .(level=paste(cell_model, treatment), col)][, .(col=unique(col)), by=.(level)]  #  manually change colors to distinguish between treatments
cl<-setNames(colorRampPalette(unique(cl$col))(nrow(cl)), cl$level)
XMAX<-5.0
h<-curve(X[[1]](x), from=-0.01, to=XMAX, n=20, ylab='', xlab='', pch=NA, col=cl[ names(X)[1] ], lty=1, lwd=10, main='', ylim=c(0, 1), xlim=c(0, XMAX), xaxt='n')
for(n in seq_along(X)[-1]){
    curve(X[[n]](x), from=head(h$x,1), to=tail(h$x,1), n=length(h$x), pch=NA, col=cl[ names(X)[n] ], lty=1, lwd=10, main='', xaxt='n', add=T)
}
axis(1, at=pretty(c(0, max(h$x))), labels=pretty(c(0, max(h$x))))
mtext('Probability', side=2, line=4, padj=-0.5, cex=2.4, las=0)
mtext(expression(log[10](1+'mean TPM')), side=1, line=4, las=0, padj=+0.1, cex=2.4)
legend(x=par('usr')[2]*0.37, y=0.69*par('usr')[4], legend=names(X), col=cl[ names(X) ], bty='n', lty=1, lwd=15, cex=1.6, y.intersp=0.6, x.intersp=0.5, seg.len=0.5, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/ecdf_MYCN_induced_targets_cell_lines.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(h,X,cl)


#  Mann-Whitney between MYCN+Tet and MYCN-Tet 
wilcox.test(x=cel[id, 'MYCN +Tet 4h'], y=cel[id, 'MYCN -Tet (DMSO)'], alternative='greater', paired=F)$p.value  #  0.1036744562

#}}}


#  [repressed targets] cumulative distributions of log10(1+mean TPMs) across risk groups/cell types
#{{{

#  repressed targets
id<-MYCN[ type %in% 'repressed', gene_id ]


#  tumors
X<-apply(log10(1+pat[id, ]), 2, function(tpm){ ecdf(tpm) })  #  tumors
h<-curve(X[[1]](x), from=-0.01, to=3.0, n=20, ylab='', xlab='', pch=NA, col=pat.meta[ match(names(X)[1], risk_group), col], lty=1, lwd=10, main='', ylim=c(0, 1), xlim=c(0, 3), xaxt='n')
for(n in seq_along(X)[-1]){
    curve(X[[n]](x), from=head(h$x,1), to=tail(h$x,1), n=length(h$x), pch=NA, col=pat.meta[ match(names(X)[n], risk_group), col], lty=1, lwd=10, main='', xaxt='n', add=T)
}
axis(1, at=pretty(c(0, max(h$x))), labels=pretty(c(0, max(h$x))))
mtext('Probability', side=2, line=4, padj=-0.5, cex=2.4, las=0)
mtext(expression(log[10](1+'mean TPM')), side=1, line=4, las=0, padj=+0.1, cex=2.4)
legend('topleft', legend=names(X), col=pat.meta[ match(names(X), risk_group), col], bty='n', lty=1, lwd=15, cex=1.6, y.intersp=0.6, x.intersp=0.5, seg.len=0.5, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/ecdf_MYCN_repressed_targets_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(h,X)


#  cell lines
X<-apply(log10(1+cel[id, ]), 2, function(tpm){ ecdf(tpm) })  #  cell lines
cl<-cel.meta[, .(level=paste(cell_model, treatment), col)][, .(col=unique(col)), by=.(level)]  #  manually change colors to distinguish between treatments
cl<-setNames(colorRampPalette(unique(cl$col))(nrow(cl)), cl$level)
XMAX<-5.0
h<-curve(X[[1]](x), from=-0.01, to=XMAX, n=20, ylab='', xlab='', pch=NA, col=cl[ names(X)[1] ], lty=1, lwd=10, main='', ylim=c(0, 1), xlim=c(0, XMAX), xaxt='n')
for(n in seq_along(X)[-1]){
    curve(X[[n]](x), from=head(h$x,1), to=tail(h$x,1), n=length(h$x), pch=NA, col=cl[ names(X)[n] ], lty=1, lwd=10, main='', xaxt='n', add=T)
}
axis(1, at=pretty(c(0, max(h$x))), labels=pretty(c(0, max(h$x))))
mtext('Probability', side=2, line=4, padj=-0.5, cex=2.4, las=0)
mtext(expression(log[10](1+'mean TPM')), side=1, line=4, las=0, padj=+0.1, cex=2.4)
legend(x=par('usr')[2]*0.37, y=0.69*par('usr')[4], legend=names(X), col=cl[ names(X) ], bty='n', lty=1, lwd=15, cex=1.6, y.intersp=0.6, x.intersp=0.5, seg.len=0.5, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/ecdf_MYCN_repressed_targets_cell_lines.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(h,X,cl)


#  Mann-Whitney between MYCN+Tet and MYCN-Tet 
wilcox.test(x=cel[id, 'MYCN +Tet 4h'], y=cel[id, 'MYCN -Tet (DMSO)'], alternative='less', paired=F)$p.value  #  0.3070269006

#}}}

#}}}




#########################
#
#
#  stratification markers
#
#
#########################




#  [tumors, cell lines] stratify by:
#
#      percent of deduplicated reads
#      proliferative index based on TPMs which removes biases on gene length
#                       
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(ProliferativeIndex)
library(pheatmap)
library(cluster)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_ecdfs.R')


#  functions
#{{{
markers_boxplot<-function(PICKED, COLS, YMIN=0, fileroot=''){
    #  quick and dirty boxplot for the tumors

    B<-setNames(lapply(names(COLS), function(r){ c(PICKED[, colnames(PICKED) %in% r]) }), names(COLS))
    YTICK<-pretty(c(YMIN, sapply(B, max)), 5)
    plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
    bp<-boxplot(B, col=COLS, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=COLS, range=0, add=T)
    mtext('TERT', side=2, line=4, padj=-0.1, las=0, cex=2.4)
    mtext(text=paste0(names(B), ' (', lengths(B)/nrow(x), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=COLS)
    mtext(text=paste(rownames(PICKED), sep='', collapse=','), side=3, line=-1, padj=+0.2, cex=1.8)

    if (fileroot!=''){
        dev.print(device=svg, file=paste0(fileroot, paste(rownames(PICKED), sep='', collapse=','), '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    }
}

#}}}


#  load reference to convert gene_id to gene_name
#  remove chrM/chrY counts
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene']
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]
seqlevels(hsa, pruning.mode='coarse')<-seqlevels(hsa)[ grep('chrM|chrY', seqlevels(hsa), invert=T) ]


#  remove failed samples and Pilot samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta<-meta[ !(failed) & !grepl('CBPilote', bid) ]


#  load featureCounts 
#  remove chrM/chrY counts
#  compute TPMs
#  convert to matrix
#  log2-transform to regularize the range
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')
gns.tpm<-unlist(totalrna[ meta[, bid] ])
gns.tpm$bid<-sub('\\.[0-9]*$', '', rownames(gns.tpm))
rownames(gns.tpm)<-NULL
gns.tpm<-gns.tpm[ gns.tpm$gene_id %in% hsa$gene_id, ]
gns.tpm$gene_name<-hsa$gene_name[ match(gns.tpm$gene_id, hsa$gene_id) ]
gns.tpm<-data.table(gns.tpm)[, .(gene_name=gene_name, length=length, counts=counts, tpm=counts/length/sum(counts/length)*1e6), by=.(bid)]
gns.tpm<-dcast(gns.tpm, bid ~ gene_name, value.var='tpm', fun.aggregate=sum)
gns.tpm<-t(data.frame(gns.tpm[, -1], row.names=gns.tpm[, bid], check.names=F))
gns.tpm<-log2(1+gns.tpm)
rm(totalrna, hsa)


#  split back metadata to cell lines and tumors
tum<-meta[ !is.na(risk_group) ]
cel<-meta[ is.na(risk_group) ]
stopifnot( nrow(tum)+nrow(cel)==nrow(meta) )
rm(meta, meta.prefailed)


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  percent DUPLICATED reads, i.e. (100-dedup)
#{{{

#  [tumors] cluster samples
#           cut tree at three clusters
#           make sure to order them 
#           add discrete classification to the metadata
x<-as.matrix(setNames(100-tum$dedup, tum$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=3)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(3))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(3), labels=c('low', 'medium', 'high'))
tum$dup.class<-x.hc[ tum$bid ]
tum$dup.class.col<-factor(x.hc[ tum$bid ], levels=c('low', 'medium', 'high'), labels=c('#228B22', '#cccccc', '#b21f1f'))


#  [cell lines] cluster samples 
#               cut tree at three clusters
#               make sure to order them 
#               add discrete classification to the metadata
x<-as.matrix(setNames(100-cel$dedup, cel$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=3)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(3))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(3), labels=c('low', 'medium', 'high'))
cel$dup.class<-x.hc[ cel$bid ]
cel$dup.class.col<-factor(x.hc[ cel$bid ], levels=c('low', 'medium', 'high'), labels=c('#228B22', '#cccccc', '#b21f1f'))


#  tables
tum[, table(risk_group, dup.class)]
#           dup.class
# risk_group low medium high
#    HR_nMNA  12     15    2
#    IMR       8      2    0
#    LR       16      9    5
#    MNA       4     16    2
#    ST4S      3     10    0


#  save
save(tum, cel, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/dup_percentage.RData')


#  visualize the classification results
#{{{

#  load back the results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/dup_percentage.RData')


#  [tumors] boxplot of duplication percentage and class
par(mar=c(2.5, 7.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(100-tum$dedup, tum$dup.class)
B.cl<-setNames(levels(tum$dup.class.col), levels(tum$dup.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('% duplicated', side=2, line=4, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_dup_per_dup_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [tumors] heatmap of annotated samples clustered by Euclidean distance in duplication percentage
x<-t(as.matrix(setNames(100-tum$dedup, tum$bid), check.names=F))
x.ex<-data.frame('Risk group   '=tum$risk_group, Duplicated=tum$dup.class, row.names=tum$bid, check.names=F)
x.cl<-setNames(list(setNames( unique(tum$col), unique(x.ex$Risk)), setNames(levels(tum$dup.class.col), levels(tum$dup.class))), colnames(x.ex))
x.d<-dist(t(x), method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
ph<-pheatmap(as.matrix(x.d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=x.hc,
        cluster_cols=x.hc,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=x.ex, annotation_row=x.ex, annotation_colors=x.cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_dup_tumors.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=18, family='Arial')
rm(x, x.ex, x.cl, x.hc, ph)

#}}}

#}}}


#  PI
#{{{

#  split TPMs of marker genes to cell lines and tumors 
#  order according to metadata
tum.tpm<-gns.tpm[rownames(gns.tpm) %in% ProliferativeIndex:::metaPCNA2, tum$bid]
cel.tpm<-gns.tpm[rownames(gns.tpm) %in% ProliferativeIndex:::metaPCNA2, cel$bid]
stopifnot( ncol(tum.tpm)+ncol(cel.tpm)==ncol(gns.tpm) )


#  add proliferative index to metadata
tum$PI<-apply(tum.tpm, 2, median)
cel$PI<-apply(cel.tpm, 2, median)


#  [tumors] cluster samples by PI
#           cut tree at two cluster
#           add discrete PI classification to the metadata
x<-as.matrix(setNames(tum$PI, tum$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
tum$PI.class<-x.hc[ tum$bid ]
tum$PI.class.col<-factor(x.hc[ tum$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  [cell lines] cluster samples by PI
#               cut tree at two cluster
#               add discrete PI classification to the metadata
x<-as.matrix(setNames(cel$PI, cel$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
cel$PI.class<-x.hc[ cel$bid ]
cel$PI.class.col<-factor(x.hc[ cel$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  tables
tum[, table(risk_group, PI.class)]
#           PI.class
# risk_group low high
#    HR_nMNA   7   22
#    IMR       4    6
#    LR       16   14
#    MNA       0   22
#    ST4S      6    7


#  save
save(tum, tum.tpm, cel, cel.tpm, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/proliferative_index.RData')


#  visualize the classification results
#{{{

#  load back the results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/proliferative_index.RData')


#  [tumors] boxplot of PI values by PI class
par(mar=c(2.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(tum$PI, tum$PI.class)
B.cl<-setNames(levels(tum$PI.class.col), levels(tum$PI.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Proliferative index', side=2, line=4, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_PI_per_proliferative_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [tumors] boxplot of PI values per risk group
#           add horizontal line which shows class separation
par(mar=c(9.0, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(tum$PI, factor(tum$risk_group, levels=unique(tum$risk_group)))
B.cl<-setNames(unique(tum$col), unique(tum$risk_group))
#YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
YTICK<-pretty(range(sapply(B, range)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
abline(h=(tum[ PI.class %in% 'low', max(PI) ]+tum[ PI.class %in% 'high', min(PI) ])/2, lty=2, lwd=2, col='black', xpd=F)
mtext('Proliferative index', side=2, line=4, padj=-0.1, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_PI_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [tumors] boxplot of PI **high** values per risk group
par(mar=c(9.0, 7.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(tum[PI.class %in% 'high', PI], factor(tum[PI.class %in% 'high', risk_group], levels=unique(tum$risk_group)))
#
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='greater')$p.value  #  0.1787
wilcox.test(x=B[['MNA']], y=B[['LR']], alternative='greater')$p.value       #  0.009445
#
B.cl<-setNames(unique(tum$col), unique(tum$risk_group))
#YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
YTICK<-pretty(range(sapply(B, range)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Proliferative index', side=2, line=5, padj=-0.1, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_PI_high_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [tumors, MNA vs HR_nMNA] boxplot of PI  (separate plot for the report)
par(mar=c(1.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(tum[ risk_group %in% c('HR_nMNA', 'MNA'), PI], tum[ risk_group %in% c('HR_nMNA', 'MNA'), risk_group])
#
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='greater')$p.value  #  0.0129
#
B.cl<-setNames(tum[ risk_group %in% c('HR_nMNA', 'MNA'), unique(col)], tum[ risk_group %in% c('HR_nMNA', 'MNA'), unique(risk_group)])
YTICK<-pretty(range(sapply(B, range)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Proliferative index', side=2, line=4, padj=-0.1, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_PI_MNA_vs_HR_nMNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [tumors] heatmap of annotated samples clustered by Euclidean distance in log2-TPMs
x<-tum.tpm
stopifnot( all.equal( colnames(x), tum$bid ) )
x.ex<-data.frame('Risk group    '=tum$risk_group, PI=tum$PI.class, row.names=tum$bid, check.names=F)
x.cl<-setNames(list(setNames( unique(tum$col), unique(x.ex$Risk)), setNames(levels(tum$PI.class.col), levels(tum$PI.class))), colnames(x.ex))
x.d<-dist(t(x), method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
ph<-pheatmap(as.matrix(x.d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=x.hc,
        cluster_cols=x.hc,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=x.ex, annotation_row=x.ex, annotation_colors=x.cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_proliferative_index_tumors.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=18, family='Arial') 
rm(x, x.ex, x.cl, x.hc, ph)


#  [tumors] PCNA
par(mar=c(14.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-tum.tpm['PCNA', , drop=F]
colnames(x)<-tum$risk_group
x.cl<-setNames(tum[, unique(col)], tum[, unique(risk_group)])
markers_boxplot(PICKED=x, COLS=x.cl, YMIN=2, fileroot='')

#}}}

#}}}


#  MYCN expression
#{{{

#  split TPMs of MYCN to cell lines and tumors 
#  order according to metadata
tum.tpm<-gns.tpm[rownames(gns.tpm) %in% 'MYCN', tum$bid, drop=F]
cel.tpm<-gns.tpm[rownames(gns.tpm) %in% 'MYCN', cel$bid, drop=F]
stopifnot( ncol(tum.tpm)+ncol(cel.tpm)==ncol(gns.tpm) )


#  add gene expression to metadata
tum$mycn<-apply(tum.tpm, 2, median)
cel$mycn<-apply(cel.tpm, 2, median)


#  [tumors] cluster samples
#           cut tree at two clusters
#           make sure to order them 
#           add discrete classification to the metadata
x<-as.matrix(setNames(tum$mycn, tum$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
tum$mycn.class<-x.hc[ tum$bid ]
tum$mycn.class.col<-factor(x.hc[ tum$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  [cell lines] cluster samples 
#               cut tree at two clusters
#               make sure to order them 
#               add discrete classification to the metadata
x<-as.matrix(setNames(cel$mycn, cel$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
cel$mycn.class<-x.hc[ cel$bid ]
cel$mycn.class.col<-factor(x.hc[ cel$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  tables
tum[, table(risk_group, mycn.class)]
#           mycn.class
# risk_group low high
#    HR_nMNA  29    0
#    IMR      10    0
#    LR       30    0
#    MNA       0   22
#    ST4S     13    0


#  save
save(tum, tum.tpm, cel, cel.tpm, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYCN_expression.RData')


#  visualize the classification results
#{{{

#  load back the results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYCN_expression.RData')


#  [tumors] boxplot of MYCN expression by MYCN class
par(mar=c(2.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(tum$mycn, tum$mycn.class)
B.cl<-setNames(levels(tum$mycn.class.col), levels(tum$mycn.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('1+TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_MYCN_TPM_per_MYCN_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [tumors] heatmap of annotated samples clustered by Euclidean distance in log2-TPMs
x<-tum.tpm
colnames(x)<-tum$bid
x.ex<-data.frame('Risk group    '=tum$risk_group, MYCN=tum$mycn.class, row.names=tum$bid, check.names=F)
x.cl<-setNames(list(setNames( unique(tum$col), unique(x.ex$Risk)), setNames(levels(tum$mycn.class.col), levels(tum$mycn.class))), colnames(x.ex))
x.d<-dist(t(x), method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
ph<-pheatmap(as.matrix(x.d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=x.hc,
        cluster_cols=x.hc,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=x.ex, annotation_row=x.ex, annotation_colors=x.cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_MYCN_TPM_tumors.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=18, family='Arial')  #  workaround the cutting of legends?
rm(x, x.ex, x.cl, x.hc, ph)


#  [tumors] MYCN
par(mar=c(14.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-tum.tpm['MYCN', , drop=F]
colnames(x)<-tum$risk_group
x.cl<-setNames(tum[, unique(col)], tum[, unique(risk_group)])
markers_boxplot(PICKED=x, COLS=x.cl, YMIN=0, fileroot='')

#}}}

#}}}


#  MYC
#{{{

#  split TPMs of MYC to cell lines and tumors 
#  order according to metadata
tum.tpm<-gns.tpm[rownames(gns.tpm) %in% 'MYC', tum$bid, drop=F]
cel.tpm<-gns.tpm[rownames(gns.tpm) %in% 'MYC', cel$bid, drop=F]
stopifnot( ncol(tum.tpm)+ncol(cel.tpm)==ncol(gns.tpm) )
rm(gns.tpm)


#  add gene expression to metadata
tum$myc<-apply(tum.tpm, 2, median)
cel$myc<-apply(cel.tpm, 2, median)


#  [tumors] cluster samples
#           cut tree at two clusters
#           make sure to order them 
#           add discrete classification to the metadata
x<-as.matrix(setNames(tum$myc, tum$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
tum$myc.class<-x.hc[ tum$bid ]
tum$myc.class.col<-factor(x.hc[ tum$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  [cell lines] cluster samples 
#               cut tree at two clusters
#               make sure to order them 
#               add discrete classification to the metadata
x<-as.matrix(setNames(cel$myc, cel$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
cel$myc.class<-x.hc[ cel$bid ]
cel$myc.class.col<-factor(x.hc[ cel$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  tables
tum[, table(risk_group, myc.class) ]
#           myc.class
# risk_group low high
#    HR_nMNA   5   24
#    IMR       3    7
#    LR        4   26
#    MNA      16    6
#    ST4S      3   10


#  save
save(tum, tum.tpm, cel, cel.tpm, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYC_expression.RData')


#  visualize the classification results
#{{{

#  load back the results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYC_expression.RData')


#  [tumors] boxplot of MYC expression by MYC class
par(mar=c(2.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(tum$myc, tum$myc.class)
B.cl<-setNames(levels(tum$myc.class.col), levels(tum$myc.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('1+TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_MYC_TPM_per_MYC_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [tumors] heatmap of annotated samples clustered by Euclidean distance in log2-TPMs
x<-tum.tpm
colnames(x)<-tum$bid
x.ex<-data.frame('Risk group    '=tum$risk_group, MYC=tum$myc.class, row.names=tum$bid, check.names=F)
x.cl<-setNames(list(setNames( unique(tum$col), unique(x.ex$Risk)), setNames(levels(tum$myc.class.col), levels(tum$myc.class))), colnames(x.ex))
x.d<-dist(t(x), method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
ph<-pheatmap(as.matrix(x.d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=x.hc,
        cluster_cols=x.hc,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=x.ex, annotation_row=x.ex, annotation_colors=x.cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_MYC_TPM_tumors.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=18, family='Arial')  #  workaround the cutting of legends?
rm(x, x.ex, x.cl, x.hc, ph)


#  [tumors] MYC
par(mar=c(14.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-tum.tpm['MYC', , drop=F]
colnames(x)<-tum$risk_group
x.cl<-setNames(tum[, unique(col)], tum[, unique(risk_group)])
markers_boxplot(PICKED=x, COLS=x.cl, YMIN=0, fileroot='')

#}}}

#}}}


#  TERT 
#{{{

#  split TPMs of TERT to cell lines and tumors 
#  order according to metadata
tum.tpm<-gns.tpm[rownames(gns.tpm) %in% 'TERT', tum$bid, drop=F]
cel.tpm<-gns.tpm[rownames(gns.tpm) %in% 'TERT', cel$bid, drop=F]
stopifnot( ncol(tum.tpm)+ncol(cel.tpm)==ncol(gns.tpm) )
rm(gns.tpm)


#  add gene expression to metadata
tum$tert<-apply(tum.tpm, 2, median)
cel$tert<-apply(cel.tpm, 2, median)


#  [tumors] cluster samples
#           cut tree at two clusters
#           make sure to order them 
#           add discrete classification to the metadata
x<-as.matrix(setNames(tum$tert, tum$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
tum$tert.class<-x.hc[ tum$bid ]
tum$tert.class.col<-factor(x.hc[ tum$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  [cell lines] cluster samples 
#               cut tree at two clusters
#               make sure to order them 
#               add discrete classification to the metadata
x<-as.matrix(setNames(cel$tert, cel$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
cel$tert.class<-x.hc[ cel$bid ]
cel$tert.class.col<-factor(x.hc[ cel$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  tables
tum[, table(risk_group, tert.class)]
#           tert.class
# risk_group low high
#    HR_nMNA   5   24
#    IMR       3    7
#    LR        4   26
#    MNA      16    6
#    ST4S      3   10


#  save
save(tum, tum.tpm, cel, cel.tpm, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/TERT_expression.RData')


#  visualize the classification results
#{{{

#  load back the results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/TERT_expression.RData')


#  [tumors] boxplot of TERT expression by TERT class
par(mar=c(2.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(tum$tert, tum$tert.class)
B.cl<-setNames(levels(tum$tert.class.col), levels(tum$tert.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('1+TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=0, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_TERT_TPM_per_TERT_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [tumors] heatmap of annotated samples clustered by Euclidean distance in log2-TPMs
x<-tum.tpm
colnames(x)<-tum$bid
x.ex<-data.frame('Risk group    '=tum$risk_group, TERT=tum$tert.class, row.names=tum$bid, check.names=F)
x.cl<-setNames(list(setNames( unique(tum$col), unique(x.ex$Risk)), setNames(levels(tum$tert.class.col), levels(tum$tert.class))), colnames(x.ex))
x.d<-dist(t(x), method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
ph<-pheatmap(as.matrix(x.d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=x.hc,
        cluster_cols=x.hc,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=x.ex, annotation_row=x.ex, annotation_colors=x.cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_TERT_TPM_tumors.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=18, family='Arial')  #  workaround the cutting of legends?
rm(x, x.ex, x.cl, x.hc, ph)


#  [tumors] TERT
par(mar=c(14.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-tum.tpm['TERT', , drop=F]
colnames(x)<-tum$risk_group
x.cl<-setNames(tum[, unique(col)], tum[, unique(risk_group)])
markers_boxplot(PICKED=x, COLS=x.cl, YMIN=0, fileroot='')

#}}}

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/proliferative_index.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYCN_expresssion.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/dedup_percentage.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYC_expresssion.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/TERT_expresssion.RData



#  [tumors, cell lines ADRN/MES markers] expression among risk groups
#                                        calculate adrenergic/mesenchymal index based on TPMs
#                                        classify neuroblastomas by ADRN/MES score
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(cluster)
library(RColorBrewer)
library(pheatmap)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_ecdfs.R')
source('~/bio/lib/grouplist2boxplot.R')


#  functions
#{{{
markers_boxplot<-function(PICK, fileroot=''){
    #  internal function that depends on many global variables and plots a nice boxplot across risk groups

    x<-tum.tpm[rownames(tum.tpm) %in% PICK, , drop=F]
    colnames(x)<-tum[ match(colnames(x), bid), risk_group ]
    B<-setNames(lapply(c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'), function(r){ c(x[, colnames(x) %in% r]) }), c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))


    B.cl<-B.cl[ names(B) ]
    YTICK<-pretty(c(0.0, sapply(B, max)), 5)
    plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
    bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
    mtext('TPM', side=2, line=5, padj=-0.1, las=0, cex=2.4)
    mtext(text=paste0(names(B), ' (', lengths(B)/nrow(x), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
    mtext(text=paste(PICK, sep='', collapse=','), side=3, line=-1, padj=+0.2, cex=1.8)

    if (fileroot!=''){
    dev.print(device=svg, file=paste0(fileroot, paste(PICK, sep='', collapse=','), '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    }
}

#}}}


#  load GENCODE v30 markers
load('/fast/groups/ag_schulte/work/reference/annotation/MES_ADR_markers/Suppl Versteeg 2017 - MES ADR Genes_gencode_v30.RData')


#  load reference to convert gene_id to gene_name
#  remove chrM/chrY counts
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene']
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]
seqlevels(hsa, pruning.mode='coarse')<-seqlevels(hsa)[ grep('chrM|chrY', seqlevels(hsa), invert=T) ]


#  remove failed samples and Pilot samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta<-meta[ !(failed) & !grepl('CBPilote', bid) ]


#  load featureCounts 
#  remove chrM/chrY counts
#  compute TPMs
#  convert to matrix
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')
gns.tpm<-unlist(totalrna[ meta[, bid] ])
gns.tpm$bid<-sub('\\.[0-9]*$', '', rownames(gns.tpm))
rownames(gns.tpm)<-NULL
gns.tpm<-gns.tpm[ gns.tpm$gene_id %in% hsa$gene_id, ]
gns.tpm$gene_name<-hsa$gene_name[ match(gns.tpm$gene_id, hsa$gene_id) ]
gns.tpm<-data.table(gns.tpm)[, .(gene_name=gene_name, length=length, counts=counts, tpm=counts/length/sum(counts/length)*1e6), by=.(bid)]
gns.tpm<-dcast(gns.tpm, bid ~ gene_name, value.var='tpm', fun.aggregate=sum)
gns.tpm<-t(data.frame(gns.tpm[, -1], row.names=gns.tpm[, bid], check.names=F))
rm(totalrna, hsa)


#  split metadata to cell lines and tumors
tum<-meta[ !is.na(risk_group) ]
cel<-meta[ is.na(risk_group) ]
stopifnot( nrow(tum)+nrow(cel)==nrow(meta) )
rm(meta, meta.prefailed)


#  define cell types from bid
#  change colors by cell type
#{{{

#  all -TR terms are SKNAS-TR
cel[ grep('-TR', bid), bid ]   
#
#  [1] "CB-SKNAS-TR-MYCN--Early-Tet1" "CB-SKNAS-TR-MYCN--Early-Tet2" "CB-SKNAS-TR-MYCN--Early-Tet3" "CB-SKNAS-TR-MYCN--late-Tet1" 
#  [5] "CB-SKNAS-TR-MYCN--late-Tet2"  "CB-SKNAS-TR-MYCN--late-Tet3"  "CB-SKNAS-TR-MYCN-Early-ETOH1" "CB-SKNAS-TR-MYCN-Early-ETOH2"
#  [9] "CB-SKNAS-TR-MYCN-Early-ETOH3" "CB-SKNAS-TR-MYCN-late-ETOH1"  "CB-SKNAS-TR-MYCN-late-ETOH2"  "CB-SKNAS-TR-MYCN-late-ETOH3" 
# [13] "CB-MDM2-TR14-R01"            


#  all _TR terms are SHEP_TR
cel[ grep('_TR', bid), bid ]   
# [1] "CB-MYCN-SHEP_TR-DMSO1-R01" "CB-MYCN-SHEP_TR-DMSO2-R01" "CB-MYCN-SHEP_TR-DMSO3-R01" "CB-MYCN-SHEP_TR-Tet1-R01"  "CB-MYCN-SHEP_TR-Tet2-R01" 
# [6] "CB-MYCN-SHEP_TR-Tet3-R01" 


#  all Rh[0-9] terms are RMS
#
#  N.B. THESE WERE REMOVED 
#
#cel[ grep('Rh[0-9]', bid), bid ]   
# [1] "CB-RMS-Rh4-R01"  "CB-RMS-Rh30-R01"


#       TR should be renamed to SKNAS 
#  SHEP_TR should be renamed to SHEP
#  Rh[0-9] should be renamed to RMS
#
sapply(strsplit(cel$bid, '-', fixed=T), '[[', 3)
#
#  [1] "IMR5"    "IMR5"    "IMR5"    "IMR5"    "IMR5"    "IMR5"    "SHEP_TR" "SHEP_TR" "SHEP_TR" "SHEP_TR" "SHEP_TR" "SHEP_TR" "TR"      "TR"     
# [15] "TR"      "TR"      "TR"      "TR"      "TR"      "TR"      "TR"      "TR"      "TR"      "TR"      "LS"      "NB1691"  "NGP"     "TR14"   
# [29] "CLBGA"   "GIMEN"   "CHLA90"  "SKNFI"   "SKNBE2" 
#
cel[, cell_type:=sapply(strsplit(cel$bid, '-', fixed=T), '[[', 3)]
cel[grep('^TR$', cell_type), cell_type:='SKNAS']
cel[grep('SHEP_TR', cell_type), cell_type:='SHEP']
#cel[grep('Rh[0-9]', cell_type), cell_type:='RMS']


#  change colors by cell type
cl<-setNames( colorRampPalette(brewer.pal(8, 'Dark2'))(cel[, length(unique(cell_type))]), cel[, unique(cell_type)] )
cel[, col:=cl[ cell_type ]]
rm(cl)

#}}}


#  split TPMs of marker genes to cell lines and tumors 
#  order according to metadata
tum.tpm<-gns.tpm[rownames(gns.tpm) %in% mes.adr$gene_name, tum$bid]
cel.tpm<-gns.tpm[rownames(gns.tpm) %in% mes.adr$gene_name, cel$bid]
stopifnot( ncol(tum.tpm)+ncol(cel.tpm)==ncol(gns.tpm) )


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(14.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B.cl<-setNames(tum[, unique(col)], tum[, unique(risk_group)])  #  needed for markers_boxplot()


#  [ADRN] NOTCH1
mes.adr[ grep('NOTCH1', gene_name) ]
markers_boxplot(PICK=c('NOTCH1'), fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_')


#  [MES] NOTCH2-4
mes.adr[ grep('NOTCH[2-4]$', gene_name) ]
markers_boxplot(PICK=c('NOTCH2', 'NOTCH3', 'NOTCH4'), fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_')


#  [MES] PRRX1
mes.adr[ grep('PRRX1', gene_name) ]
markers_boxplot(PICK=c('PRRX1'), fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_')


#  [MES] SEMA3C
mes.adr[ grep('SEMA3C', gene_name) ]
markers_boxplot(PICK=c('SEMA3C'), fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_')


#  [MES] COL12A1
mes.adr[ grep('COL12A1', gene_name) ]
markers_boxplot(PICK=c('COL12A1'), fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_')


#  [MES] EGR1
mes.adr[ grep('EGR1', gene_name) ]
markers_boxplot(PICK=c('EGR1'), fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_')


#  [MES] MAML2
mes.adr[ grep('MAML2', gene_name) ]
markers_boxplot(PICK=c('MAML2'), fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_')


#  [ADRN] ALK
mes.adr[ grep('ALK', gene_name) ]
markers_boxplot(PICK=c('ALK'), fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_')


#  [ADRN] DLK1
mes.adr[ grep('DLK1', gene_name) ]
markers_boxplot(PICK=c('DLK1'), fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_')


#  [ADRN] GATA2-3 (the rest are lowly expressed)
mes.adr[ grep('GATA[2-3]', gene_name) ]
markers_boxplot(PICK=c('GATA2', 'GATA3'), fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_')


#  [ADRN] PHOX2A,2B
mes.adr[ grep('PHOX2[AB]', gene_name) ]
markers_boxplot(PICK=c('PHOX2A', 'PHOX2B'), fileroot='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_')


#  [tumors] all mesenchymal/adrenergic markers 
#{{{

#  split to mesenchymal/adrenergic 
x<-log2(1+tum.tpm)
colnames(x)<-tum[ match(colnames(x), bid), risk_group ]
m<-x[ rownames(x) %in% mes.adr[ cell_type %in% 'MES', gene_name ], ]
a<-x[ rownames(x) %in% mes.adr[ cell_type %in% 'ADRN', gene_name ], ]
stopifnot( nrow(a)+nrow(m)==nrow(x) )
B.m<-setNames(lapply(c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'), function(r){ c(m[, colnames(m) %in% r]) }), c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))
B.a<-setNames(lapply(c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'), function(r){ c(a[, colnames(a) %in% r]) }), c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))
rm(x)


#  all pairwise Mann-Whitney U tests
n<-matrix(names(B.m)[ combn(length(B.m), 2) ], nrow=2)
p.m<-p.a<-setNames(rep(NA, ncol(n)), apply(n, 2, paste0, collapse='-'))
for(i in seq_len(ncol(n))){
    p.m[i]<-wilcox.test(B.m[[ n[1, i] ]], B.m[[ n[2, i] ]], alternative='two.sided')$p.value
    p.a[i]<-wilcox.test(B.a[[ n[1, i] ]], B.a[[ n[2, i] ]], alternative='two.sided')$p.value
}
p.m[ setdiff(names(p.m), names(p.m[ p.m<0.05 ])) ]
#
# named numeric(0)
#
p.a[ setdiff(names(p.a), names(p.a[ p.a<0.05 ])) ]
#
#       LR-IMR 
# 0.4554079139 
#
rm(n, i)


#  mesenchymal ecdf
B<-B.m
N<-lengths(B)/nrow(m)
B.cl<-B.cl[ names(B) ]
my_ecdfs(B, B.cl, XLIM=c(0.0, 8), XLAB=expression(log[2](1+'TPM')), MAIN=paste(nrow(m), 'mesenchymal markers'), LTY=1, LWD=8, svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_across_risk_groups_TPMs_for_mesenchymal_all.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#
#  boxplot
#
#YTICK<-pretty(c(0.0, sapply(B, max)), 5)
YTICK<-pretty(c(0.0, 8), 4)
par(mar=c(14.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=1.5, add=T)
mtext(expression(log[2](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', N, ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=paste(nrow(m), 'mesenchymal markers'), side=3, line=-1, padj=+0.2, cex=1.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_mesenchymal_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  adrenergic boxplot
B<-B.a
N<-lengths(B)/nrow(a)
B.cl<-B.cl[ names(B) ]
my_ecdfs(B, B.cl, XLIM=c(0.0, 8), XLAB=expression(log[2](1+'TPM')), MAIN=paste(nrow(a), 'adrenergic markers'), LTY=1, LWD=8, svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_across_risk_groups_TPMs_for_adrenergic_all.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#
#  boxplot
#
#YTICK<-pretty(c(0.0, sapply(B, max)), 5)
YTICK<-pretty(c(0.0, 10.0), 4)
par(mar=c(14.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=1.5, add=T)
mtext(expression(log[2](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', N, ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=paste(nrow(a), 'adrenergic markers'), side=3, line=-1, padj=+0.2, cex=1.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_adrenergic_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] filtered mesenchymal/adrenergic markers 
#{{{

#  split to mesenchymal/adrenergic 
#  keep well-expressed markers with coefficient of variation of at least a given value
THRESHOLD<-3
CV<-0.2
x<-log2(1+tum.tpm)
colnames(x)<-tum[ match(colnames(x), bid), risk_group ]
m<-x[ rownames(x) %in% mes.adr[ cell_type %in% 'MES', gene_name ], ]
a<-x[ rownames(x) %in% mes.adr[ cell_type %in% 'ADRN', gene_name ], ]
stopifnot( nrow(a)+nrow(m)==nrow(x) )
cv.m<-sqrt(diag(var(t(m))))/rowMeans(m)
cv.a<-sqrt(diag(var(t(a))))/rowMeans(a)
m<-m[ rowMeans(m)>=THRESHOLD & cv.m>=CV, , drop=F]
a<-a[ rowMeans(a)>=THRESHOLD & cv.a>=CV, , drop=F]
cv.m<-cv.m[ rownames(m) ]
cv.a<-cv.a[ rownames(a) ]
B.m<-setNames(lapply(c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'), function(r){ c(m[, colnames(m) %in% r]) }), c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))
B.a<-setNames(lapply(c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'), function(r){ c(a[, colnames(a) %in% r]) }), c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))
rm(x, THRESHOLD, CV)


#  all pairwise Mann-Whitney U tests
n<-matrix(names(B.m)[ combn(length(B.m), 2) ], nrow=2)
p.m<-p.a<-setNames(rep(NA, ncol(n)), apply(n, 2, paste0, collapse='-'))
for(i in seq_len(ncol(n))){
    p.m[i]<-wilcox.test(B.m[[ n[1, i] ]], B.m[[ n[2, i] ]], alternative='two.sided')$p.value
    p.a[i]<-wilcox.test(B.a[[ n[1, i] ]], B.a[[ n[2, i] ]], alternative='two.sided')$p.value
}
p.m[ setdiff(names(p.m), names(p.m[ p.m<0.05 ])) ]
#
# named numeric(0)
#
p.a[ setdiff(names(p.a), names(p.a[ p.a<0.05 ])) ]
#
#     ST4S-IMR       LR-IMR   LR-HR_nMNA  HR_nMNA-MNA 
# 0.9226601578 0.1051103280 0.1367905728 0.1219189706 
rm(n, i)


#  mesenchymal ecdf
B<-B.m
N<-lengths(B)/nrow(m)
B.cl<-B.cl[ names(B) ]
my_ecdfs(B, B.cl, XLIM=c(0.0, 8), XLAB=expression(log[2](1+'TPM')), MAIN=paste(nrow(m), 'mesenchymal markers'), LTY=1, LWD=8, svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_across_risk_groups_TPMs_for_mesenchymal.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#
#  boxplot
#
#YTICK<-pretty(c(0.0, sapply(B, max)), 5)
YTICK<-pretty(c(0.0, 8.0), 4)  #  zoom-in and go back to default whisker range to see the significant difference by eye as well
par(mar=c(14.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=1.5, add=T)
mtext(expression(log[2](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', N, ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=paste(nrow(m), 'mesenchymal markers'), side=3, line=-1, padj=+0.2, cex=1.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_mesenchymal.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  adrenergic boxplot
B<-B.a
N<-lengths(B)/nrow(a)
B.cl<-B.cl[ names(B) ]
my_ecdfs(B, B.cl, XLIM=c(0.0, 8), XLAB=expression(log[2](1+'TPM')), MAIN=paste(nrow(a), 'adrenergic markers'), LTY=1, LWD=8, svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_across_risk_groups_TPMs_for_adrenergic.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#
#  boxplot
#
#YTICK<-pretty(c(0.0, sapply(B, max)), 5)
YTICK<-pretty(c(0.0, 10.0), 4)  #  zoom-in and go back to default whisker range to see the significant difference by eye as well
par(mar=c(14.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=1.5, add=T)
mtext(expression(log[2](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', N, ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=paste(nrow(a), 'adrenergic markers'), side=3, line=-1, padj=+0.2, cex=1.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_TPMs_for_adrenergic.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [cell lines] all mesenchymal/adrenergic markers 
#{{{

#  split to mesenchymal/adrenergic 
x<-log2(1+cel.tpm)
colnames(x)<-cel[ match(colnames(x), bid), cell_type ]
m<-x[ rownames(x) %in% mes.adr[ cell_type %in% 'MES', gene_name ], ]
a<-x[ rownames(x) %in% mes.adr[ cell_type %in% 'ADRN', gene_name ], ]
stopifnot( nrow(a)+nrow(m)==nrow(x) )
B.m<-setNames(lapply(cel[, unique(cell_type)], function(r){ c(m[, colnames(m) %in% r]) }), cel[, unique(cell_type)])
B.a<-setNames(lapply(cel[, unique(cell_type)], function(r){ c(a[, colnames(a) %in% r]) }), cel[, unique(cell_type)])
rm(x)


#  all pairwise Mann-Whitney U tests
n<-matrix(names(B.m)[ combn(length(B.m), 2) ], nrow=2)
p.m<-p.a<-setNames(rep(NA, ncol(n)), apply(n, 2, paste0, collapse='-'))
for(i in seq_len(ncol(n))){
    p.m[i]<-wilcox.test(B.m[[ n[1, i] ]], B.m[[ n[2, i] ]], alternative='two.sided')$p.value
    p.a[i]<-wilcox.test(B.a[[ n[1, i] ]], B.a[[ n[2, i] ]], alternative='two.sided')$p.value
}
p.m<-p.m[ p.m<0.05 ]
p.a<-p.a[ p.a<0.05 ]
rm(n, i)


#  mesenchymal ecdf
B<-B.m
N<-lengths(B)/nrow(m)
B.cl<-setNames(cel[, unique(col)], cel[, unique(cell_type)])
stopifnot( all.equal( names(B.cl), names(B) ) )
my_ecdfs(B, B.cl, XLIM=c(0.0, 8.0), XLAB=expression(log[2](1+'TPM')), MAIN=paste(nrow(m), 'mesenchymal markers'), LTY=1, LWD=8, LEGEND='bottomright', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_across_cell_types_TPMs_for_mesenchymal_all.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#
#  boxplot
#
#YTICK<-pretty(c(0.0, sapply(B, max)), 5)
YTICK<-pretty(c(0.0, 10.0), 4)  #  zoom-in and go back to default whisker range to see the significant difference by eye as well
par(mar=c(12.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=1.5, add=T)
mtext(expression(log[2](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', N, ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=paste(nrow(m), 'mesenchymal markers'), side=3, line=-1, padj=+0.2, cex=1.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_cell_types_TPMs_for_mesenchymal_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  adrenergic boxplot
B<-B.a
N<-lengths(B)/nrow(a)
B.cl<-setNames(cel[, unique(col)], cel[, unique(cell_type)])
stopifnot( all.equal( names(B.cl), names(B) ) )
my_ecdfs(B, B.cl, XLIM=c(0.0, 10), XLAB=expression(log[2](1+'TPM')), MAIN=paste(nrow(a), 'adrenergic markers'), LTY=1, LWD=8, LEGEND='bottomright', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_across_cell_types_TPMs_for_adrenergic_all.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#
#  boxplot
#
#YTICK<-pretty(c(0.0, sapply(B, max)), 5)
YTICK<-pretty(c(0.0, 10.0), 4)  #  zoom-in and go back to default whisker range to see the significant difference by eye as well
par(mar=c(12.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=1.5, add=T)
mtext(expression(log[2](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', N, ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=paste(nrow(a), 'adrenergic markers'), side=3, line=-1, padj=+0.2, cex=1.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_cell_types_TPMs_for_adrenergic_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [cell lines] DE mesencnymal/adrenergic markers in circARID1A KD
#{{{

#  load DE results
#  convert gene_ids to gene_names
env<-new.env()
local({
    load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_kable-ready_genes+circRNAs_circARID1A SI_circARID1A SCR.RData')
    up.mes.adr<-lapply(lapply(up.mes.adr, '[[', 1), function(g){ all.gns[g, 'gene_name'] })
    down.mes.adr<-lapply(lapply(down.mes.adr, '[[', 1), function(g){ all.gns[g, 'gene_name'] })
    }, envir=env)
up.mes.adr<-get('up.mes.adr', envir=env)
down.mes.adr<-get('down.mes.adr', envir=env)
rm(env)


#  split to mesenchymal/adrenergic 
x<-log2(1+cel.tpm)
colnames(x)<-cel[ match(colnames(x), bid), cell_type ]
m<-x[ rownames(x) %in% c(up.mes.adr$mes, down.mes.adr$mes), ]
a<-x[ rownames(x) %in% c(up.mes.adr$adr, down.mes.adr$adr), ]
B.m<-setNames(lapply(cel[, unique(cell_type)], function(r){ c(m[, colnames(m) %in% r]) }), cel[, unique(cell_type)])
B.a<-setNames(lapply(cel[, unique(cell_type)], function(r){ c(a[, colnames(a) %in% r]) }), cel[, unique(cell_type)])
rm(x)


#  mesenchymal ecdf
B<-B.m
N<-lengths(B)/nrow(m)
B.cl<-setNames(cel[, unique(col)], cel[, unique(cell_type)])
stopifnot( all.equal( names(B.cl), names(B) ) )
my_ecdfs(B, B.cl, XLIM=c(0.0, 8.0), XLAB=expression(log[2](1+'TPM')), MAIN=paste(nrow(m), 'mesenchymal markers'), LTY=1, LWD=8, LEGEND='bottomright', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_across_cell_types_TPMs_for_mesenchymal_circARID1A_KD.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#
#  boxplot
#
#YTICK<-pretty(c(0.0, sapply(B, max)), 5)
YTICK<-pretty(c(0.0, 10.0), 4)  #  zoom-in and go back to default whisker range to see the significant difference by eye as well
par(mar=c(12.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=1.5, add=T)
mtext(expression(log[2](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', N, ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=paste(nrow(m), 'mesenchymal markers'), side=3, line=-1, padj=+0.2, cex=1.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_cell_types_TPMs_for_mesenchymal_circARID1A_KD.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  adrenergic boxplot
B<-B.a
N<-lengths(B)/nrow(a)
B.cl<-setNames(cel[, unique(col)], cel[, unique(cell_type)])
stopifnot( all.equal( names(B.cl), names(B) ) )
my_ecdfs(B, B.cl, XLIM=c(0.0, 10), XLAB=expression(log[2](1+'TPM')), MAIN=paste(nrow(a), 'adrenergic markers'), LTY=1, LWD=8, LEGEND='bottomright', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_across_cell_types_TPMs_for_adrenergic_circARID1A_KD.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#
#  boxplot
#
#YTICK<-pretty(c(0.0, sapply(B, max)), 5)
YTICK<-pretty(c(0.0, 10.0), 4)  #  zoom-in and go back to default whisker range to see the significant difference by eye as well
par(mar=c(12.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=1.5, add=T)
mtext(expression(log[2](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', N, ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=paste(nrow(a), 'adrenergic markers'), side=3, line=-1, padj=+0.2, cex=1.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_cell_types_TPMs_for_adrenergic_circARID1A_KD.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  adrenergic/mesenchymal score similar to https://www.nature.com/articles/ng.3899#Sec2:
#
#      In each sample the genes with TPM>=1 were ranked by increasing TPM order. The median percentile of the ADRN and MES list of genes
#      corresponded to the ADRN and MES scores (highest score indicates higher overall expresssion for the corresonding group).
#{{{

#  compute ADRN/MES scores
AM.tum<-t(apply(tum.tpm, 2, function(tpms){ 
    keep<-tpms>=1
    gn<-rownames(tum.tpm)[ keep ]
    tpms<-tpms[ keep ]
    gn<-gn[ order(tpms, decreasing=F) ]
    c(adrn.score=median( which(gn %in% mes.adr[ cell_type %in% 'ADRN', gene_name])/length(gn) ),
      mes.score=median( which(gn %in% mes.adr[ cell_type %in% 'MES', gene_name])/length(gn) ) 
    )
    }))
AM.cel<-t(apply(cel.tpm, 2, function(tpms){ 
    keep<-tpms>=1
    gn<-rownames(cel.tpm)[ keep ]
    tpms<-tpms[ keep ]
    gn<-gn[ order(tpms, decreasing=F) ]
    c(adrn.score=median( which(gn %in% mes.adr[ cell_type %in% 'ADRN', gene_name])/length(gn) ),
      mes.score=median( which(gn %in% mes.adr[ cell_type %in% 'MES', gene_name])/length(gn) )
    )
    }))


#  add relevant metadata
AM.tum<-data.table(bid=rownames(AM.tum), AM.tum)
tum<-tum[, c('bid', 'risk_group', 'col'), with=F]
AM.tum<-tum[ AM.tum, on='bid']
AM.cel<-data.table(bid=rownames(AM.cel), AM.cel)
cel<-cel[, c('bid', 'cell_model', 'cell_type', 'col'), with=F]
AM.cel<-cel[ AM.cel, on='bid']
rm(tum, cel)


#  cluster tumors by MES score 
x<-as.matrix(setNames(AM.tum$mes.score, AM.tum$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
AM.tum$mes.class<-x.hc[ AM.tum$bid ]
AM.tum$mes.class.col<-factor(x.hc[ AM.tum$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  cluster tumors by ADRN score 
x<-as.matrix(setNames(AM.tum$adrn.score, AM.tum$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
AM.tum$adrn.class<-x.hc[ AM.tum$bid ]
AM.tum$adrn.class.col<-factor(x.hc[ AM.tum$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  cluster cell lines by MES score 
x<-as.matrix(setNames(AM.cel$mes.score, AM.cel$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
AM.cel$mes.class<-x.hc[ AM.cel$bid ]
AM.cel$mes.class.col<-factor(x.hc[ AM.cel$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  cluster cell lines by ADRN score 
x<-as.matrix(setNames(AM.cel$adrn.score, AM.cel$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-cutree(hclust(x.d, method='ward.D2'), k=2)
summary(silhouette(x.hc, dist=x.d))
o<-sapply(split(names(x.hc), x.hc), function(n){ mean(x[n, ,drop=F]) })
x.hc<-unclass(factor(x.hc, levels=names(sort(o)), labels=as.character(seq_len(2))))
attr(x.hc, 'levels')<-NULL
x.hc<-factor(x.hc, levels=seq_len(2), labels=c('low', 'high'))
AM.cel$adrn.class<-x.hc[ AM.cel$bid ]
AM.cel$adrn.class.col<-factor(x.hc[ AM.cel$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  tables
AM.tum[, table(risk_group, adrn.class)]
#           adrn.class
# risk_group low high
#    HR_nMNA   1   28
#    IMR       2    8
#    LR        7   23
#    MNA       1   21
#    ST4S      1   12
#
AM.tum[, table(risk_group, mes.class)]
#           mes.class
# risk_group low high
#    HR_nMNA  13   16
#    IMR       4    6
#    LR       11   19
#    MNA      18    4
#    ST4S      3   10


#  save
save(AM.tum, AM.cel, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/adrenergic-mesenchymal_score.RData')

#}}}


#  visualize the classification results
#{{{

#  load back the results 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/adrenergic-mesenchymal_score.RData')


#  tumors
#{{{

#  boxplot of scores stratified by risk groups
B<-data.frame(AM.tum[, c('risk_group', 'adrn.score', 'mes.score'), with=F])
colnames(B)<-c('risk_group', 'ADRN', 'MES')
B<-lapply(split(B[, c('ADRN', 'MES')], factor(B$risk_group, levels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))), t)
B.cl<-setNames(AM.tum[, unique(col)], AM.tum[, unique(risk_group)])[names(B)]
grouplist2boxplot(L=B, L.COL=B.cl, YLAB='Score', YLAB.CEX=2.4, XLAB.CEX=2.4, XLAS=1, XTEXT.ADJ=0.5, YTEXT.LINE=4, LEGEND='topright', mar=c(2.0, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_ADRN_MES_scores_per_risk_group.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B, B.cl)


#  boxplot of MES scores stratified by MES class
B<-data.frame(AM.tum[, c('mes.class', 'mes.score'), with=F])
colnames(B)<-c('class', 'MES')
B<-split(B[, 'MES'], B$class)
B.cl<-setNames(levels(AM.tum$mes.class.col), levels(AM.tum$mes.class))
YTICK<-pretty(c(min(sapply(B, min, na.rm=T)), max(sapply(B, max, na.rm=T))), 4)
par(mar=c(2.5, 8.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('MES score', side=2, line=5, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), '(', lengths(B), ')'), side=1, line=0, at=seq_along(B), las=1, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_MES_scores_per_MES_class.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B, B.cl, YTICK, bp)


#  boxplot of ADRN scores stratified by ADRN class
B<-data.frame(AM.tum[, c('adrn.class', 'adrn.score'), with=F])
colnames(B)<-c('class', 'ADRN')
B<-split(B[, 'ADRN'], B$class)
B.cl<-setNames(levels(AM.tum$adrn.class.col), levels(AM.tum$adrn.class))
YTICK<-pretty(c(min(sapply(B, min, na.rm=T)), max(sapply(B, max, na.rm=T))), 4)
par(mar=c(2.5, 8.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('ADRN score', side=2, line=5, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), '(', lengths(B), ')'), side=1, line=0, at=seq_along(B), las=1, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_ADRN_scores_per_ADRN_class.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B, B.cl, YTICK, bp)


#  [ADRN/MES scores together] heatmap of annotated samples clustered by Euclidean distance in ADRN/MES scores
x<-data.frame(AM.tum[, -1], row.names=AM.tum[, bid], check.names=F)
x.ex<-data.frame('Risk group    '=x$risk_group, ADRN=x$adrn.score, MES=x$mes.score, row.names=rownames(x), check.names=F)
x.cl<-factor(setNames(x$risk_group, x$col), levels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))
x.cl<-list(Risk=setNames(unique(names(x.cl)), levels(x.cl)),
    ADRN=colorRampPalette(brewer.pal(9,'YlOrRd'))(4),
    MES=colorRampPalette(brewer.pal(9,'YlOrRd'))(4))
x.d<-dist(x[, c('adrn.score', 'mes.score')], method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
ph<-pheatmap(as.matrix(x.d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=x.hc,
        cluster_cols=x.hc,
        annotation_names_row=T, annotation_names_col=T, annotation_legend=T,
        annotation_col=x.ex, annotation_row=x.ex, annotation_colors=x.cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_ADRN_MES_scores_across_risk_groups.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial') 
rm(x, x.ex, x.cl, x.hc, ph)

#}}}


#  cell lines
#{{{

#  barplot of ADRN/MES scores
B<-t(data.frame(AM.cel[, c('adrn.score', 'mes.score'), with=F], row.names=AM.cel$bid))  #  use bid for now because rownames needs to be unique
rownames(B)<-c('ADRN', 'MES')
colnames(B)<-AM.cel[, cell_type]  #  overwrite the colnames now
B.cl<-setNames( AM.cel[, col], colnames(B) )
o<-order(B['ADRN', ], decreasing=T)  #  order by ADRN score
B<-B[, o]
B.cl<-B.cl[o]
stopifnot(all.equal(names(B.cl), colnames(B)))
YTICK<-pretty(c(0, max(B)), 4)
par(mar=c(8.5, 6.5, 1.0, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
bp<-barplot(B, beside=T, border='white', col='white', axes=F, axisnames=F, xlab='', ylab='', las=1, xaxs='i', yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
for(i in 1:ncol(B)){ 
    b<-B
    b[, -i]<-NA    #  remove all values but the current column values
    if (i%%2){  #  odd-numbered columns should have this color scheme
        barplot(as.matrix(b), border='white', col=c('#b21f1f', '#cccccc'), axes=F, axisnames=F, beside=T, yaxt='n', add=T)
    } else {    #  even-numbered columns should have this color scheme
        barplot(as.matrix(b), border='white', col=c('#b21f1f', '#cccccc'), axes=F, axisnames=F, beside=T, yaxt='n', add=T)
    }
}
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext('Score', side=2, line=4, padj=-0.5, las=0, cex=2.4)
mtext(text=colnames(B), side=1, line=0, at=colMeans(bp), las=2, adj=1, cex=2.4, col=B.cl)
legend(x=par('usr')[2]*0.89, y=par('usr')[4]*1.00, legend=rownames(B), col=c('#b21f1f', '#cccccc'), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.65, x.intersp=0.2, seg.len=0.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_ADRN_MES_scores_cell_lines.svg', width=25, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B, B.cl, o, bp, i, b)


#  boxplot of scores stratified by cell type
B<-data.frame(AM.cel[, c('cell_type', 'adrn.score', 'mes.score'), with=F])
colnames(B)<-c('cell_type', 'ADRN', 'MES')
B<-lapply(split(B[, c('ADRN', 'MES')], factor(B$cell_type)), t)
B.cl<-setNames(AM.cel[, unique(col)], AM.cel[, unique(cell_type)])[ names(B) ]
grouplist2boxplot(L=B, L.COL=B.cl, YLAB='Score', YLAB.CEX=2.4, XLAB.CEX=2.4, XLAS=1, XTEXT.ADJ=0.5, YTEXT.LINE=4, LEGEND='topright', mar=c(2.0, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
rm(B, B.cl)


#  boxplot of MES scores stratified by MES class
B<-data.frame(AM.cel[, c('mes.class', 'mes.score'), with=F])
colnames(B)<-c('class', 'MES')
B<-split(B[, 'MES'], B$class)
B.cl<-setNames(levels(AM.cel$mes.class.col), levels(AM.cel$mes.class))
YTICK<-pretty(c(min(sapply(B, min, na.rm=T)), max(sapply(B, max, na.rm=T))), 4)
par(mar=c(2.0, 8.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('MES score', side=2, line=5, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), '(', lengths(B), ')'), side=1, line=0, at=seq_along(B), las=1, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_MES_scores_per_MES_class_cell_lines.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B, B.cl, YTICK, bp)


#  boxplot of ADRN scores stratified by ADRN class
B<-data.frame(AM.cel[, c('adrn.class', 'adrn.score'), with=F])
colnames(B)<-c('class', 'ADRN')
B<-split(B[, 'ADRN'], B$class)
B.cl<-setNames(levels(AM.cel$adrn.class.col), levels(AM.cel$adrn.class))
YTICK<-pretty(c(min(sapply(B, min, na.rm=T)), max(sapply(B, max, na.rm=T))), 4)
par(mar=c(2.0, 8.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('ADRN score', side=2, line=5, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), '(', lengths(B), ')'), side=1, line=0, at=seq_along(B), las=1, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_ADRN_scores_per_ADRN_class_cell_lines.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B, B.cl, YTICK, bp)


#  [ADRN/MES scores together] heatmap of annotated samples clustered by Euclidean distance in ADRN/MES scores
x<-data.frame(AM.cel[, -1], row.names=AM.cel[, bid], check.names=F)
x.ex<-data.frame('Cell type    '=x$cell_type, ADRN=x$adrn.score, MES=x$mes.score, row.names=rownames(x), check.names=F)
x.cl<-list('Cell type    '=setNames(unique(x$col), unique(x$cell_type)), 
    ADRN=colorRampPalette(brewer.pal(9,'Greys'))(4),
    MES=colorRampPalette(brewer.pal(9,'YlOrRd'))(4))
x.d<-dist(x[, c('adrn.score', 'mes.score')], method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
ph<-pheatmap(as.matrix(x.d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=x.hc,
        cluster_cols=x.hc,
        annotation_names_row=T, annotation_names_col=T, annotation_legend=T,
        annotation_col=x.ex, annotation_row=x.ex, annotation_colors=x.cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_ADRN_MES_scores_across_cell_types.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x, x.ex, x.cl, x.hc, ph)

#}}}

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/adrenergic-mesenchymal_score.RData



#  Compile a version of the metadata with all classifications included.
#  Clinical data explanations:
#
#      status : 
#
#          0 = Alive without progression.
#          1 = Died from disease.
#          2 = Alive after suffering a progression event.
#              For OS survival this should be in the status 0 bin.
#              For EFS survival this should be in the status 1 bin.
#
#      age:
#
#          young = <547 days
#            old = >546 days
#
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(ComplexHeatmap)
library(cluster)
library(RColorBrewer)
library(ggplot2)
library(survival)
library(survminer)
library(maxstat)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/cluster2groups.R')


#  load all the metrics and classifications already done
#{{{

#  load the nonfailed non-Pilote metadata
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta.tum<-meta.tum[ !(failed) & !grepl('Pilote', bid) ]
meta.cel<-meta.cel[ !(failed) & !grepl('Pilote', bid) ]
rm(meta.prefailed)


#  load duplication percentage class for the samples and annotate the data objects
l<-ls()
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/dup_percentage.RData')
m<-tum[, c('bid', 'dup.class', 'dup.class.col'), with=F]
meta.tum<-meta.tum[ m, on='bid' ]
m<-cel[, c('bid', 'dup.class', 'dup.class.col'), with=F]
meta.cel<-meta.cel[ m, on='bid' ]
rm(list=setdiff(ls(), l))


#  load PI class for the samples and annotate the data objects
l<-ls()
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/proliferative_index.RData')
m<-tum[, c('bid', 'PI', 'PI.class', 'PI.class.col'), with=F]
meta.tum<-meta.tum[ m, on='bid' ]
m<-cel[, c('bid', 'PI', 'PI.class', 'PI.class.col'), with=F]
meta.cel<-meta.cel[ m, on='bid' ]
rm(list=setdiff(ls(), l))


#  load ADRN/MES scores and MES classification for the tumor samples and annotate the data objects
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/adrenergic-mesenchymal_score.RData')
m<-AM.tum[, c('bid', 'adrn.score', 'mes.score', 'adrn.class', 'mes.class', 'adrn.class.col', 'mes.class.col'), with=F]
meta.tum<-meta.tum[ m, on='bid' ]
m<-AM.cel[, c('bid', 'cell_type', 'adrn.score', 'mes.score', 'adrn.class', 'mes.class', 'adrn.class.col', 'mes.class.col'), with=F]
meta.cel<-meta.cel[ m, on='bid' ]
rm(list=setdiff(ls(), l))


#  load MYCN class and annotate the data objects
l<-ls()
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYCN_expression.RData')
m<-tum[, c('bid', 'mycn', 'mycn.class', 'mycn.class.col'), with=F]
meta.tum<-meta.tum[ m, on='bid' ]
m<-cel[, c('bid', 'mycn', 'mycn.class', 'mycn.class.col'), with=F]
meta.cel<-meta.cel[ m, on='bid' ]
rm(list=setdiff(ls(), l))


#  load MYC class and annotate the data objects
l<-ls()
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYC_expression.RData')
m<-tum[, c('bid', 'myc', 'myc.class', 'myc.class.col'), with=F]
meta.tum<-meta.tum[ m, on='bid' ]
m<-cel[, c('bid', 'myc', 'myc.class', 'myc.class.col'), with=F]
meta.cel<-meta.cel[ m, on='bid' ]
rm(list=setdiff(ls(), l))


#  load TERT class and annotate the data objects
l<-ls()
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/TERT_expression.RData')
m<-tum[, c('bid', 'tert', 'tert.class', 'tert.class.col'), with=F]
meta.tum<-meta.tum[ m, on='bid' ]
m<-cel[, c('bid', 'tert', 'tert.class', 'tert.class.col'), with=F]
meta.cel<-meta.cel[ m, on='bid' ]
rm(list=setdiff(ls(), l))


#  [tumors] load sex determination to assign missing sex assignments
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/sex_determination.RData')
sex<-sex[, c('bid', 'sex')]
meta.tum<-meta.tum[ sex, on='bid']
meta.tum[ sex!=SEX | is.na(SEX) ]
#              bid sex PAT_ID_BERLIN risk_group        col INSS_STAGE COHORT TRIAL_PROTOCOL INRG_HR_risk_group MYCN_FC MYCN_INI_STATUS AGE_days
# 1: CB3009-11-R01   M        CB3009    HR_nMNA chocolate1          4   <NA>           NB43                 HR     1.5         non-amp     1419
#    EFS_days OS_days EFS_bin OS_bin status  SEX
# 1:      688    2661       1      0      2 <NA>
meta.tum[ sex!=SEX | is.na(SEX), SEX:=sex]
meta.tum[, sex:=NULL]
meta.tum[, table(SEX)]
#
# SEX
#  F  M 
# 36 68 
rm(sex)


#  load unified circRNAs
#  compute number of isoforms per sample
#  add them to the metadata
l<-ls()
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
iso<-data.table(data.frame(mcols(CIRCS.all)[, c('bid', 'circ_name')]))[, .(circRNAcounts=.N), by=.(bid)]
meta.tum<-meta.tum[ iso, on='bid' ]
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_cell_models.RData')
iso<-data.table(data.frame(mcols(CIRCS.all)[, c('bid', 'circ_name')]))[, .(circRNAcounts=.N), by=.(bid)][ bid %in% meta.cel$bid ]
meta.cel<-meta.cel[ iso, on='bid' ]
rm(list=setdiff(ls(), l))


#  sum of circRNA/mRNA CPMs per sample
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
tum<-unlist(lin.cir)[bid %in% meta.tum$bid][, c('c.cpm', 'l.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)][, .(c.cpm=sum(c.cpm, na.rm=T), l.cpm=sum(l.cpm, na.rm=T)), by=.(bid)]
meta.tum<-meta.tum[ tum, on='bid' ]
cel<-unlist(lin.cir)[bid %in% meta.cel$bid][, c('c.cpm', 'l.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)][, .(c.cpm=sum(c.cpm, na.rm=T), l.cpm=sum(l.cpm, na.rm=T)), by=.(bid)]
meta.cel<-meta.cel[ cel, on='bid' ]
rm(list=setdiff(ls(), l))


#  load linear vs circular junction quantification results
#  compute circARID1A and ARID1A CPMs based on the circular and linear junction counts, respectively, and the total (linear+circular) counts
#  compute circAKAP12 CPM based on the circular and linear junction counts, respectively, and the total (linear+circular) counts
#  add them to the metadata
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
N<-unlist(lin.cir)[bid %in% meta.tum$bid]
#
#  circARID1A
#
N<-unlist(lin.cir)[bid %in% meta.tum$bid]
N<-N[ circ_name %in% 'ENSG00000117713.20|ARID1A_chr1+26729651-26732792' ]
N<-N[, c('circARID1Acpm', 'ARID1Acpm'):=list(c.count/total*1e6, l.count.out/total*1e6)][, c('bid', 'circARID1Acpm', 'ARID1Acpm')]
meta.tum<-meta.tum[ N, on='bid' ]
meta.tum[ is.na(circARID1Acpm), circARID1Acpm:=0 ]
meta.tum[ is.na(ARID1Acpm), ARID1Acpm:=0 ]
#
#  circAKAP12
#
N<-unlist(lin.cir)[bid %in% meta.tum$bid]
N<-N[ circ_name %in% 'ENSG00000131016.17|AKAP12_chr6+151348711-151353752' ]
N<-N[, c('circAKAP12cpm', 'AKAP12cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)][, c('bid', 'circAKAP12cpm', 'AKAP12cpm')]
meta.tum<-meta.tum[ N, on='bid' ]
meta.tum[ is.na(circAKAP12cpm), circAKAP12cpm:=0 ]
meta.tum[ is.na(AKAP12cpm), AKAP12cpm:=0 ]
rm(list=setdiff(ls(), l))

#}}}


#  classification into low/high clusters
#{{{

#  [circARID1A, log10-transformed] classify tumors into clusters of low and high 
l<-ls()
d<-data.frame(meta.tum[, c('OS_days', 'OS_bin', 'circARID1Acpm'), with=F])
maxstat.test(Surv(OS_days, OS_bin) ~ circARID1Acpm, data=d, smethod='LogRank', pmethod='exactGauss', alpha=0.05, abseps=0.01)$estimate
# estimated cutpoint 
#              5.563 
#
#w<-setNames(ifelse(meta.tum$circARID1Acpm<5.6, 'low', 'high'), meta.tum$bid)
w<-cluster2groups(METRIC=log10(1+meta.tum$circARID1Acpm), NAMES=meta.tum$bid, META=meta.tum$risk_group, GROUPS=2, GROUP.NAMES=c('low', 'high'))
meta.tum$circARID1A.class<-factor(w[['km']][ meta.tum$bid ], levels=c('low', 'high'))
rm(list=setdiff(ls(), l))


#  [ARID1A, log10-transformed] classify tumors into clusters of low and high 
l<-ls()
d<-data.frame(meta.tum[, c('OS_days', 'OS_bin', 'ARID1Acpm'), with=F])
maxstat.test(Surv(OS_days, OS_bin) ~ ARID1Acpm, data=d, smethod='LogRank', pmethod='exactGauss', alpha=0.05, abseps=0.01)$estimate
# estimated cutpoint 
#              15.31 
#
#w<-setNames(ifelse(meta.tum$ARID1Acpm<15.3, 'low', 'high'), meta.tum$bid)
w<-cluster2groups(METRIC=log10(1+meta.tum$ARID1Acpm), NAMES=meta.tum$bid, META=meta.tum$risk_group, GROUPS=2, GROUP.NAMES=c('low', 'high'))
meta.tum$ARID1A.class<-factor(w[['km']][ meta.tum$bid ], levels=c('low', 'high'))
rm(list=setdiff(ls(), l))


#  [circAKAP12, log10-transformed] classify tumors into clusters of low and high 
l<-ls()
d<-data.frame(meta.tum[, c('OS_days', 'OS_bin', 'circAKAP12cpm'), with=F])
maxstat.test(Surv(OS_days, OS_bin) ~ circAKAP12cpm, data=d, smethod='LogRank', pmethod='exactGauss', alpha=0.05, abseps=0.01)$estimate
# estimated cutpoint 
#              1.169 
#
#w<-setNames(ifelse(meta.tum$circAKAP12cpm<1.2, 'low', 'high'), meta.tum$bid)
w<-cluster2groups(METRIC=log10(1+meta.tum$circAKAP12cpm), NAMES=meta.tum$bid, META=meta.tum$risk_group, GROUPS=2, GROUP.NAMES=c('low', 'high'))
meta.tum$circAKAP12.class<-factor(w[['km']][ meta.tum$bid ], levels=c('low', 'high'))
rm(list=setdiff(ls(), l))


#  [AKAP12, log10-transformed] classify tumors into clusters of low and high 
l<-ls()
d<-data.frame(meta.tum[, c('OS_days', 'OS_bin', 'AKAP12cpm'), with=F])
maxstat.test(Surv(OS_days, OS_bin) ~ AKAP12cpm, data=d, smethod='LogRank', pmethod='exactGauss', alpha=0.05, abseps=0.01)$estimate
# estimated cutpoint 
#               38.6 
#
#w<-setNames(ifelse(meta.tum$AKAP12cpm<38.6, 'low', 'high'), meta.tum$bid)
w<-cluster2groups(METRIC=log10(1+meta.tum$AKAP12cpm), NAMES=meta.tum$bid, META=meta.tum$risk_group, GROUPS=2, GROUP.NAMES=c('low', 'high'))
meta.tum$AKAP12.class<-factor(w[['km']][ meta.tum$bid ], levels=c('low', 'high'))
rm(list=setdiff(ls(), l))


#  [circRNA counts, log10-transformed] classify tumors into clusters of low and high 
#l<-ls()
#d<-data.frame(meta.tum[, c('OS_days', 'OS_bin', 'circRNAcounts'), with=F])
#maxstat.test(Surv(OS_days, OS_bin) ~ circRNAcounts, data=d, smethod='LogRank', pmethod='exactGauss', alpha=0.05, abseps=0.01)$estimate
# estimated cutpoint 
#               2698 
#
#w<-setNames(ifelse(meta.tum$circRNAcounts<2698, 'low', 'high'), meta.tum$bid)
#w<-cluster2groups(METRIC=log10(1+meta.tum$circRNAcounts), NAMES=meta.tum$bid, META=meta.tum$risk_group, GROUPS=2, GROUP.NAMES=c('low', 'high'))
#meta.tum$circRNAcounts.class<-factor(w[['km']][ meta.tum$bid ], levels=c('low', 'high'))
#rm(list=setdiff(ls(), l))


#  [circRNA sum of CPMs, log10-transformed] classify tumors into clusters of low and high
#l<-ls()
#d<-data.frame(meta.tum[, c('OS_days', 'OS_bin', 'c.cpm'), with=F])
#maxstat.test(Surv(OS_days, OS_bin) ~ c.cpm, data=d, smethod='LogRank', pmethod='exactGauss', alpha=0.05, abseps=0.01)$estimate
# estimated cutpoint 
#               4350 
#
#w<-setNames(ifelse(meta.tum$c.cpm<4350, 'low', 'high'), meta.tum$bid)
#w<-cluster2groups(METRIC=log10(1+meta.tum$c.cpm), NAMES=meta.tum$bid, META=meta.tum$risk_group, GROUPS=2, GROUP.NAMES=c('low', 'high'))
#meta.tum$c.cpm.class<-factor(w[['km']][ meta.tum$bid ], levels=c('low', 'high'))
#rm(list=setdiff(ls(), l))

#}}}


#  classification of circRNA counts and circRNA sum of CPMs into bins
#{{{

#  [circRNA counts] 
l<-ls()
h<-hist(meta.tum$circRNAcounts, breaks=c(190, 500, 800, 1200, 2000, seq(2500, 4300, 300), 5100), plot=F)$breaks
w<-setNames(cut(meta.tum$circRNAcounts, breaks=h, labels=tail(h, -1), right=T, ordered_result=F), meta.tum$bid)  #  remove lowest and right-close
meta.tum$circRNAcounts.class<-w
rm(list=setdiff(ls(), l))


#  [circRNA sum of CPMs, log10-transformed] classify tumors into clusters of low and high
l<-ls()
h<-hist(meta.tum$c.cpm, breaks=c(1200, seq(2200, 3800, 200), 4200, 4600, 6400, 12600), plot=F)$breaks
w<-setNames(cut(meta.tum$c.cpm, breaks=h, labels=tail(h, -1), right=T, ordered_result=F), meta.tum$bid)  #  remove lowest and right-close
meta.tum$c.cpm.class<-w
rm(list=setdiff(ls(), l))

#}}}


#  correlation circRNA counts with different metrics
meta.tum[, cor(nreads, circRNAcounts, method='spearman')]    #  0.4947
meta.tum[, cor(features, circRNAcounts, method='spearman')]  #  0.483


#  table of OS_bin and status
meta.tum[, table(OS_bin, status) ]
#       status
# OS_bin  0  1  2
#      0 63  0 15
#      1  0 26  0
#
meta.tum[, table(EFS_bin, status) ]
#        status
# EFS_bin  0  1  2
#       0 63  0  0
#       1  0 26 15


#  binarize age
meta.tum[, AGE_bin:=factor(ifelse(AGE_days<547, 'young', 'old'), levels=c('young', 'old'))]
meta.tum[, AGE_bin.col:=ifelse(AGE_bin=='young', '#cccccc', '#b21f1f')]


#  save
save(meta.tum, meta.cel, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata_with_classes.RData')

#}}}
#
#  => /fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata_with_classes.RData



#  [tumors] clinical data analysis
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(ComplexHeatmap)
library(cluster)
library(RColorBrewer)
library(ggplot2)
library(survival)
library(survminer)
library(maxstat)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load classifications
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata_with_classes.RData')


#  recycle
x11(width=20, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  stratification plots
#{{{

#  fix the ggplot2 theme
#{{{
theme_set(theme_bw(base_size=35) + theme(panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_blank(),
    plot.margin=margin(t=1.0, b=0.1, l=0.5, r=2.0, unit='cm'),
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
    legend.position=c(0.1, 1.10), 
    legend.margin=margin(),
    legend.key=element_blank(),
    legend.title=element_text(size=28, family='Arial', face='plain', margin=margin(b=0.5, unit='cm')),
    legend.text=element_text(size=28, family='Arial', face='plain', margin=margin(l=-0.2, unit='cm'))
))
#}}}


#  circRNA counts
#{{{

#  distributions of the low and high groups
B<-data.frame(meta.tum[, c('circRNAcounts', 'circRNAcounts.class')])
lh<-list(low=levels(B[, 2])[1:floor(length(levels(B[, 2]))/2)], high=levels(B[, 2])[ceiling(length(levels(B[, 2]))/2):length(levels(B[, 2]))]) 
XTICK<-pretty(c(0, max(as.numeric(levels(B[, 2])))), 5)
#
wilcox.test(x=B[ B[, 2] %in% lh[['low']], 1], y=B[ B[, 2] %in% lh[['high']], 1], alternative='less')$p.value   #  1.408e-18
#
ggplot(B, aes(circRNAcounts, fill=circRNAcounts.class)) +
    geom_histogram(breaks=as.numeric(levels(B[, 2])), color='white', position='identity', show.legend=F) +
    scale_y_continuous(limits=c(0, 18), n.breaks=7, expand=expansion(mult=c(0.01, 0))) +
    scale_fill_manual(values=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(B[, 2]))))) +
    scale_x_continuous(limits=range(XTICK), n.breaks=8, expand=expansion(mult=c(0.01, 0))) +
    coord_cartesian(clip='off') +
    labs(x='Number of circRNAs', y='Frequency')
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_per_stratum_circRNAcounts.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  distributions colored by status
B<-data.frame(meta.tum[, c('circRNAcounts', 'status')])
B$status<-factor(B$status, levels=c(0, 1, 2))
YTICK<-pretty(c(0, max(B[, 1], na.rm=T)), 5)
#
wilcox.test(x=B[ B[, 2] %in% 1, 1], y=B[ B[, 2] %in% 0, 1], alternative='less')$p.value     #  0.1823
wilcox.test(x=B[ B[, 2] %in% 1, 1], y=B[ B[, 2] %in% 2, 1], alternative='less')$p.value     #  0.5479
wilcox.test(x=B[ B[, 2] %in% 0, 1], y=B[ B[, 2] %in% 2, 1], alternative='greater')$p.value  #  0.1908
#
ggplot(B, aes(status, circRNAcounts, fill=status)) + 
    stat_boxplot(geom='errorbar', width=0.2, position=position_dodge(width=0.8), linetype='solid', coef=Inf, size=0.4) +
    geom_boxplot(aes(ymin=..lower.., ymax=..upper..), position=position_dodge(width=0.8), coef=Inf, color='white', show.legend=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    scale_fill_manual(values=c('#008b45', '#4682b4', '#ffa54f')) +
    labs(x='', y='Number of circRNAs', color='') +
    coord_cartesian(clip='off') + 
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text=element_text(size=45),
          #axis.title.x=element_text(color='black', margin=margin(t=5, b=5)),  #  we need this or labels will be clipped
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=0.5, hjust=0.5),
    )
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_per_status_circRNAcounts.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  circRNA sum of CPMs
#{{{

#  distributions of the low and high groups
B<-data.frame(meta.tum[, c('c.cpm', 'c.cpm.class')])
lh<-list(low=levels(B[, 2])[1:floor(length(levels(B[, 2]))/2)], high=levels(B[, 2])[ceiling(length(levels(B[, 2]))/2):length(levels(B[, 2]))]) 
XTICK<-pretty(c(0, max(as.numeric(levels(B[, 2])))), 5)
#
wilcox.test(x=B[ B[, 2] %in% lh[['low']], 1], y=B[ B[, 2] %in% lh[['high']], 1], alternative='less')$p.value   #  1.109e-18
#
ggplot(B, aes(c.cpm, fill=c.cpm.class)) +
    geom_histogram(breaks=as.numeric(levels(B[, 2])), color='white', position='identity', show.legend=F) +
    scale_y_continuous(limits=c(0, 14), n.breaks=7, expand=expansion(mult=c(0.01, 0))) +
    scale_fill_manual(values=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(B[, 2]))))) +
    scale_x_continuous(limits=range(XTICK), n.breaks=8, expand=expansion(mult=c(0.01, 0))) +
    coord_cartesian(clip='off') +
    labs(x='Expression of circRNAs (CPM)', y='Frequency', fill='')
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_per_stratum_circRNA_sum_CPM.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  distributions colored by status
B<-data.frame(meta.tum[, c('c.cpm', 'status')])
B$status<-factor(B$status, levels=c(0, 1, 2))
YTICK<-pretty(c(0, max(B[, 1], na.rm=T)), 5)
#
wilcox.test(x=B[ B[, 2] %in% 1, 1], y=B[ B[, 2] %in% 0, 1], alternative='less')$p.value     #  0.1307
wilcox.test(x=B[ B[, 2] %in% 1, 1], y=B[ B[, 2] %in% 2, 1], alternative='less')$p.value     #  0.3103
wilcox.test(x=B[ B[, 2] %in% 0, 1], y=B[ B[, 2] %in% 2, 1], alternative='greater')$p.value  #  0.4097
#
ggplot(B, aes(status, c.cpm, fill=status)) + 
    stat_boxplot(geom='errorbar', width=0.2, position=position_dodge(width=0.8), linetype='solid', coef=Inf, size=0.4) +
    geom_boxplot(aes(ymin=..lower.., ymax=..upper..), position=position_dodge(width=0.8), coef=Inf, color='white', show.legend=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    scale_fill_manual(values=c('#008b45', '#4682b4', '#ffa54f')) +
    labs(x='', y='Expression of circRNAs (CPM)', color='') +
    coord_cartesian(clip='off') + 
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text=element_text(size=45),
          #axis.title.x=element_text(color='black', margin=margin(t=5, b=5)),  #  we need this or labels will be clipped
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=0.5, hjust=0.5),
    )
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_per_status_circRNA_sum_CPM.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  circARID1A CPM
#{{{

#  distributions of the low and high groups
B<-data.frame(meta.tum[, c('circARID1Acpm', 'circARID1A.class')])
XTICK<-pretty(c(0, max(B[, 1], na.rm=T)), 5)
#
wilcox.test(x=B[ B[, 2] %in% 'low', 1], y=B[ B[, 2] %in% 'high', 1], alternative='less')$p.value   #  8.912e-18
#
ggplot(B, aes(circARID1Acpm, fill=circARID1A.class)) +
    geom_histogram(breaks=seq(min(XTICK), max(XTICK), length.out=36), color='white', position='identity') +
    scale_y_continuous(limits=c(0, 12), n.breaks=5, expand=expansion(mult=c(0.01, 0))) +
    scale_fill_manual(values=c('#cccccc', '#b21f1f')) +
    scale_x_continuous(limits=range(XTICK), n.breaks=8, expand=expansion(mult=c(0.01, 0))) +
    guides(fill=guide_legend(keywidth=unit(2, 'cm'), override.aes=list(size=14))) +
    labs(x='circARID1A CPM', y='Frequency', fill='') + 
    theme(legend.position=c(1.0, 1.10))
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_per_stratum_circARID1A.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  distributions colored by status
B<-data.frame(meta.tum[, c('circARID1Acpm', 'status')])
B$status<-factor(B$status, levels=c(0, 1, 2))
YTICK<-pretty(c(0, max(B[, 1], na.rm=T)), 5)
#
wilcox.test(x=B[ B[, 2] %in% 1, 1], y=B[ B[, 2] %in% 0, 1], alternative='less')$p.value   #  0.2713
wilcox.test(x=B[ B[, 2] %in% 1, 1], y=B[ B[, 2] %in% 2, 1], alternative='less')$p.value   #  0.569
wilcox.test(x=B[ B[, 2] %in% 0, 1], y=B[ B[, 2] %in% 2, 1], alternative='less')$p.value   #  0.7431
#
ggplot(B, aes(status, circARID1Acpm, fill=status)) + 
    stat_boxplot(geom='errorbar', width=0.2, position=position_dodge(width=0.8), linetype='solid', coef=Inf, size=0.4) +
    geom_boxplot(aes(ymin=..lower.., ymax=..upper..), position=position_dodge(width=0.8), coef=Inf, color='white', show.legend=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    scale_fill_manual(values=c('#008b45', '#4682b4', '#ffa54f')) +
    labs(x='', y='circARID1A CPM', color='') +
    coord_cartesian(clip='off') + 
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text=element_text(size=45),
          #axis.title.x=element_text(color='black', margin=margin(t=5, b=5)),  #  we need this or labels will be clipped
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=0.5, hjust=0.5),
    )
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_per_status_circARID1A.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  ARID1A CPM
#{{{

#  distributions of the low and high groups
B<-data.frame(meta.tum[, c('ARID1Acpm', 'ARID1A.class')])
XTICK<-pretty(c(0, max(B[, 1], na.rm=T)), 5)
#
wilcox.test(x=B[ B[, 2] %in% 'low', 1], y=B[ B[, 2] %in% 'high', 1], alternative='less')$p.value   #  1.109e-18
#
ggplot(B, aes(ARID1Acpm, fill=ARID1A.class)) +
    geom_histogram(breaks=seq(min(XTICK), max(XTICK), length.out=30), color='white', position='identity') +
    scale_y_continuous(limits=c(0, 20), n.breaks=5, expand=expansion(mult=c(0.01, 0))) +
    scale_fill_manual(values=c('#cccccc', '#b21f1f')) +
    scale_x_continuous(limits=range(XTICK), n.breaks=8, expand=expansion(mult=c(0.01, 0))) +
    guides(fill=guide_legend(keywidth=unit(2, 'cm'), override.aes=list(size=14))) +
    labs(x='ARID1A CPM', y='Frequency', fill='') + 
    theme(legend.position=c(1.0, 1.10))
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_per_stratum_ARID1A.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  distributions colored by status
B<-data.frame(meta.tum[, c('ARID1Acpm', 'status')])
B$status<-factor(B$status, levels=c(0, 1, 2))
YTICK<-pretty(c(0, max(B[, 1], na.rm=T)), 5)
#
wilcox.test(x=B[ B[, 2] %in% 1, 1], y=B[ B[, 2] %in% 0, 1], alternative='less')$p.value   #  0.4766
wilcox.test(x=B[ B[, 2] %in% 1, 1], y=B[ B[, 2] %in% 2, 1], alternative='less')$p.value   #  0.3691
wilcox.test(x=B[ B[, 2] %in% 0, 1], y=B[ B[, 2] %in% 2, 1], alternative='less')$p.value   #  0.3999
#
ggplot(B, aes(status, ARID1Acpm, fill=status)) + 
    stat_boxplot(geom='errorbar', width=0.2, position=position_dodge(width=0.8), linetype='solid', coef=Inf, size=0.4) +
    geom_boxplot(aes(ymin=..lower.., ymax=..upper..), position=position_dodge(width=0.8), coef=Inf, color='white', show.legend=F) +
    scale_y_continuous(limits=range(YTICK), breaks=YTICK, expand=c(0, 0)) +
    scale_fill_manual(values=c('#008b45', '#4682b4', '#ffa54f')) +
    labs(x='', y='ARID1A CPM', color='') +
    coord_cartesian(clip='off') + 
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text=element_text(size=45),
          #axis.title.x=element_text(color='black', margin=margin(t=5, b=5)),  #  we need this or labels will be clipped
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=0, margin=margin(l=0, b=0, t=0, r=0, unit='cm'), vjust=0.5, hjust=0.5),
    )
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_per_status_ARID1A.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}


#  survival analysis
#{{{

#  convert to data.frame and define some factors 
d<-data.frame(meta.tum)
d$SEX<-factor(d$SEX, levels=c('M', 'F'), labels=c('Male', 'Female'))
d$MYCN_INI_STATUS<-factor(d$MYCN_INI_STATUS, levels=c('non-amp', 'amp'))


#  [MYCN status] Kaplan-Meier survival curves
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ MYCN_INI_STATUS, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('non MNA', 'MNA'),  #  non-amp, amp level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_MYCN_INI_STATUS.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ MYCN_INI_STATUS, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('non MNA', 'MNA'),  #  non-amp, amp level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_MYCN_INI_STATUS.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [circRNA counts] Kaplan-Meier survival curves
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ circRNAcounts.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('circRNAs low', 'circRNAs high'),  #  low, high level order
    #palette=c('#E7B800', '#2E9FDF'),
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$circRNAcounts.class)))),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_circRNAcounts.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ circRNAcounts.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('circRNAs low', 'circRNAs high'),  #  low, high level order
    #palette=c('#E7B800', '#2E9FDF'),
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$circRNAcounts.class)))),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_circRNAcounts.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [circRNA sum of CPMs] Kaplan-Meier survival curves
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ c.cpm.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('CPMs low', 'CPMs high'),  #  low, high level order
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$c.cpm.class)))),
    #palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_circRNA_sum_CPMs.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ c.cpm.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('CPMs low', 'CPMs high'),  #  low, high level order
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$c.cpm.class)))),
    #palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_circRNA_sum_CPMs.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [circARID1A] Kaplan-Meier survival curve
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ circARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('circARID1A low', 'circARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_circARID1A.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ circARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('circARID1A low', 'circARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_circARID1A.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [ARID1A] Kaplan-Meier survival curve
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ ARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('ARID1A low', 'ARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_ARID1A.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ ARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('ARID1A low', 'ARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_ARID1A.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [circAKAP12] Kaplan-Meier survival curve
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ circAKAP12.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('circAKAP12 low', 'circAKAP12 high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_circAKAP12.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ circAKAP12.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('circAKAP12 low', 'circAKAP12 high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_circAKAP12.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [AKAP12] Kaplan-Meier survival curve
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ AKAP12.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('AKAP12 low', 'AKAP12 high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_AKAP12.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ AKAP12.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('AKAP12 low', 'AKAP12 high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_AKAP12.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}


#  [HR_nMNA+MNA] survival analysis
#{{{

#  prepare the data and define the factors
d<-data.frame(meta.tum[ risk_group %in% c('HR_nMNA', 'MNA')])
d$SEX<-factor(d$SEX, levels=c('M', 'F'), labels=c('Male', 'Female'))
d$MYCN_INI_STATUS<-factor(d$MYCN_INI_STATUS, levels=c('non-amp', 'amp'))


#  [MYCN status] Kaplan-Meier survival curves
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ MYCN_INI_STATUS, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('non MNA', 'MNA'),  #  non-amp, amp level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_MYCN_INI_STATUS_HR_nMNA+MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ MYCN_INI_STATUS, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('non MNA', 'MNA'),  #  non-amp, amp level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_MYCN_INI_STATUS_HR_nMNA+MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [circRNA counts] Kaplan-Meier survival curves
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ circRNAcounts.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('circRNAs low', 'circRNAs high'),  #  low, high level order
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$circRNAcounts.class)))),
    #palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_circRNAcounts_HR_nMNA+MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ circRNAcounts.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('circRNAs low', 'circRNAs high'),  #  low, high level order
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$circRNAcounts.class)))),
    #palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_circRNAcounts_HR_nMNA+MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [circRNA sum of CPMs] Kaplan-Meier survival curves
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ c.cpm.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('CPMs low', 'CPMs high'),  #  low, high level order
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$c.cpm.class)))),
    #palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_circRNA_sum_CPMs_HR_nMNA+MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ c.cpm.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('CPMs low', 'CPMs high'),  #  low, high level order
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$c.cpm.class)))),
    #palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_circRNA_sum_CPMs_HR_nMNA+MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [circARID1A] Kaplan-Meier survival curve
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ circARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('circARID1A low', 'circARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_circARID1A_HR_nMNA+MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ circARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('circARID1A low', 'circARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_circARID1A_HR_nMNA+MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [ARID1A] Kaplan-Meier survival curve
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ ARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('ARID1A low', 'ARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_ARID1A_HR_nMNA+MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ ARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('ARID1A low', 'ARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_ARID1A_HR_nMNA+MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}


#  [-MNA] survival analysis
#{{{

#  prepare the data and define the factors
d<-data.frame(meta.tum[ ! risk_group %in% c('MNA')])
d$SEX<-factor(d$SEX, levels=c('M', 'F'), labels=c('Male', 'Female'))
d$MYCN_INI_STATUS<-factor(d$MYCN_INI_STATUS, levels=c('non-amp', 'amp'))


#  [circRNA counts] Kaplan-Meier survival curves
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ circRNAcounts.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('circRNAs low', 'circRNAs high'),  #  low, high level order
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$circRNAcounts.class)))),
    #palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_circRNAcounts_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ circRNAcounts.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('circRNAs low', 'circRNAs high'),  #  low, high level order
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$circRNAcounts.class)))),
    #palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_circRNAcounts_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [circRNA sum of CPMs] Kaplan-Meier survival curves
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ c.cpm.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('CPMs low', 'CPMs high'),  #  low, high level order
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$c.cpm.class)))),
    #palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_circRNA_sum_CPMs_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ c.cpm.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=F,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    legend.labs=sub('^.*=', '', names(fit$strata)),
    #legend.labs=c('CPMs low', 'CPMs high'),  #  low, high level order
    palette=rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(d$c.cpm.class)))),
    #palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_circRNA_sum_CPMs_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [circARID1A] Kaplan-Meier survival curve
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ circARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('circARID1A low', 'circARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_circARID1A_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ circARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('circARID1A low', 'circARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_circARID1A_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [ARID1A] Kaplan-Meier survival curve
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ ARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('ARID1A low', 'ARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_ARID1A_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ ARID1A.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('ARID1A low', 'ARID1A high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_ARID1A_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [circAKAP12] Kaplan-Meier survival curve
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ circAKAP12.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('circAKAP12 low', 'circAKAP12 high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_circAKAP12_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ circAKAP12.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('circAKAP12 low', 'circAKAP12 high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_circAKAP12_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [AKAP12] Kaplan-Meier survival curve
#{{{

#  OS
fit<-survfit(Surv(OS_days, event=OS_bin) ~ AKAP12.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Overall survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('AKAP12 low', 'AKAP12 high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_OS_vs_AKAP12_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  EFS
fit<-survfit(Surv(EFS_days, event=EFS_bin) ~ AKAP12.class, data=d)
fit.summary<-surv_summary(fit)
p<-ggsurvplot(fit,
    pval=T,
    pval.size=10,
    pval.coord=c(median(seq(min(fit$time), max(fit$time), length.out=5)), 1.0),
    conf.int=F,
    linetype='solid', 
    size=4, 
    censor.shape='+',
    censor.size=20, 
    ylab='Probability',
    xlab='Event-free survival (days)',
    break.time.by=500,
    xlim=range(pretty(fit$time, 5)),
    ylim=c(0, 1),
    axes.offset=F,
    legend='top',
    legend.title='',
    #legend.labs=sub('^.*=', '', names(fit$strata)),
    legend.labs=c('AKAP12 low', 'AKAP12 high'),  #  low, high level order
    palette=c('#E7B800', '#2E9FDF'),
    font.x=40, 
    font.y=40,
    font.tickslab=25,
    ggtheme=theme_classic2(base_size=40) + 
            theme(plot.margin=margin(t=0.1, b=0.1, l=0.5, r=1.5, unit='cm'),
                  axis.line=element_line(color='black', size=0.5),
                  axis.ticks=element_line(color='black', size=0.5),
                  axis.text.x=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(t=0.2, b=0.8, unit='cm')), 
                  axis.text.y=element_text(color='black', size=40, family='Arial', face='plain', margin=margin(l=0.8, r=0.2, unit='cm')), 
                  axis.ticks.length=unit(0.1,  'cm'))
    )
gt<-ggplot_gtable(ggplot_build(p$plot))
gt$layout$clip<-'off'
grid.draw(gt)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/Kaplan-Meier_EFS_vs_AKAP12_-MNA.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}


dev.off()

#}}}



#  [tumors] visualize the classifications
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(ComplexHeatmap)
library(cluster)
library(RColorBrewer)
library(survival)
library(survminer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_ecdfs.R')


#  load the classifications
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata_with_classes.RData')


#  recycle
x11(width=20, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  heatmap of annotated samples clustered by Euclidean distance in log2-TPMs of MYCN expression
m<-copy(meta.tum)
x<-t(as.matrix(m$mycn))
colnames(x)<-m$bid
x.ex<-data.frame('Risk group    '=m$risk_group, MYCN=m$mycn.class, MYC=m$myc.class, TERT=m$tert.class, PI=m$PI.class, ADRN=m$adrn.class, row.names=m$bid, check.names=F)
x.cl<-setNames(list(
    setNames(unique(m$col), unique(m$risk_group)),
    setNames(levels(m$mycn.class.col), levels(m$mycn.class)),  #  extra 'moderate' level added because of cell lines
    setNames(levels(m$myc.class.col), levels(m$myc.class)),  #  extra 'moderate' level added because of cell lines
    setNames(levels(m$tert.class.col), levels(m$tert.class)),
    setNames(levels(m$PI.class.col), levels(m$PI.class)),  #  extra 'moderate' level added because of cell lines
    setNames(levels(m$adrn.class.col), levels(m$adrn.class))
    ), colnames(x.ex))
x.d<-dist(t(x), method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
ph<-pheatmap(as.matrix(x.d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=x.hc,
        cluster_cols=x.hc,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=x.ex, annotation_row=x.ex, annotation_colors=x.cl,
        show_rownames=F, show_colnames=F, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
print(ph)
rm(m, x, x.ex, x.cl, x.hc, ph)



#  sum of circRNA CPMs per risk group
#{{{

#  split into groups and compute log10(1+CPMs)
B<-split(log10(meta.tum$c.cpm), factor(meta.tum$risk_group, levels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA')))
B.cl<-unique(meta.tum[, c('risk_group', 'col')])
B.cl<-setNames(B.cl$col, B.cl$risk_group)
B.cl<-B.cl[ names(B) ]
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value  #  0.05641
wilcox.test(x=B[['MNA']], y=B[['LR']], alternative='less')$p.value       #  0.2876


#  boxplots
par(mar=c(9.0, 8.0, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(sapply(B, range, na.rm=T), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10]('circRNA expression')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, adj=1.00, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_sum_of_CPMs_per_sample_circular_junctions_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}



#  oncoPlot
#
#      risk group
#      INSS_STAGE
#      MYCN_INI_STATUS
#      MYCN expression
#      MYC expression
#      PI
#      OS_bin
#      EFS_bin
#      AGE_bin
#      circRNAcounts.class
#      c.cpm.class
#
#{{{

#  ordered by risk group primarily and circRNA numbers secondarily
#  show the sample sizes on the top annotation
#{{{

#  order by circRNA numbers
#  the risk group ordering is enforced when we ask for column splits but you have to define the factor level order
#  we also need to have cohort sizes added to the group annotation on top of the oncoPlot
r<-meta.tum[, table(risk_group)][c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA')]
r<-setNames(as.numeric(r), names(r))
m<-copy(meta.tum)[, risk_group:=factor(risk_group, levels=names(r), labels=paste0(names(r), ' (', r, ')'))][ order(circRNAcounts) ]
mat<-t(data.frame(m[, c('risk_group', 'INSS_STAGE', 'OS_bin', 'EFS_bin', 'AGE_bin', 'MYCN_INI_STATUS', 'PI.class', 'circRNAcounts.class', 'c.cpm.class'), with=F][, c('INSS_STAGE', 'OS_bin', 'EFS_bin', 'circRNAcounts.class', 'c.cpm.class'):=list(paste0('INSS', INSS_STAGE), paste0('OS', OS_bin), paste0('EFS', EFS_bin), paste0('counts', circRNAcounts.class), paste0('cpms', c.cpm.class))]))
rownames(mat)<-c('Risk group', 'INSS stage', 'OS', 'EFS', 'Age', 'MYCN status', 'Proliferative index', 'circRNA number', 'circRNA expression')
cl<-setNames( unique(m$col), unique(m$risk_group) )
cl<-append(cl, setNames( colorRampPalette(brewer.pal(9, 'Blues'))(12)[c(seq(6, 12, 2), 4)], c('INSS1', 'INSS2', 'INSS3', 'INSS4', 'INSS5') ) )
cl<-append(cl, setNames( c('#cccccc', '#b21f1f'), c('OS0', 'OS1') ) )
cl<-append(cl, setNames( c('#cccccc', '#b21f1f'), c('EFS0', 'EFS1') ) )
cl<-append(cl, setNames( unique(m$AGE_bin.col), unique(m$AGE_bin) ) )
cl<-append(cl, setNames( c('#cccccc', '#b21f1f'), c('non-amp', 'amp') ) )
cl<-append(cl, setNames( levels(m$PI.class.col), levels(m$PI.class) ) )
cl<-append(cl, setNames( rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(m$circRNAcounts.class)))), paste0('counts', levels(m$circRNAcounts.class)) ) )
cl<-append(cl, setNames( rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(m$c.cpm.class)))), paste0('cpms', levels(m$c.cpm.class)) ) )
N<-setNames(as.numeric(levels(m$circRNAcounts.class)), levels(m$circRNAcounts.class))
N<-N[c(1, floor(length(N)/2), length(N))]
#names(N)[seq(2, length(N), 2)]<-''
E<-setNames(as.numeric(levels(m$c.cpm.class)), levels(m$c.cpm.class))
E<-E[c(1, floor(length(E)/2), length(E))]
#names(E)[seq(2, length(E), 2)]<-''
lg<-packLegend(
    #Legend(title='Risk group', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'),
    #    at=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'), labels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'), 
    #    legend_gp=gpar(fill=cl[c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA')]), ncol=1),
    Legend(title='INSS stage', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'),
        at=c('INSS1', 'INSS2', 'INSS3', 'INSS4', 'INSS5'), labels=c('1', '2', '3', '4', '5'), 
        legend_gp=gpar(fill=cl[c('INSS1', 'INSS2', 'INSS3', 'INSS4', 'INSS5')]), ncol=1),
    Legend(title='OS/EFS event', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'), at=c('OS0', 'OS1'), 
        labels=c('0', '1'), legend_gp=gpar(fill=cl[c('OS0', 'OS1')]), ncol=1),
    Legend(title='Age', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'), at=c('young', 'old'),
        labels=c('<547 days', '>546 days'), legend_gp=gpar(fill=cl[c('young', 'old')]), ncol=1),
    Legend(title='MYCN status', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'),
        at=c('non-amp', 'amp'), labels=c('not amplified', 'amplified'), legend_gp=gpar(fill=cl[c('non-amp', 'amp')]), ncol=1),
    Legend(title='circRNA number', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), 
        grid_height=unit(2,'mm'), at=N, labels=names(N), legend_gp=gpar(fill=cl[paste0('counts', N)]), ncol=1),
    Legend(title='circRNA expression', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14),
        grid_height=unit(2,'mm'), at=E, labels=names(E), legend_gp=gpar(fill=cl[paste0('cpms', E)]), ncol=1),
    #Legend(title='circRNA number', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), 
    #    grid_height=unit(4,'mm'), at=as.numeric(levels(m$circRNAcounts.class)), labels=levels(m$circRNAcounts.class), 
    #    legend_gp=gpar(fill=cl[ paste0('counts', levels(m$circRNAcounts.class))]), ncol=2),
    #Legend(title='circRNA expression', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14),
    #    grid_height=unit(4,'mm'), at=as.numeric(levels(m$c.cpm.class)), labels=levels(m$c.cpm.class),
    #    legend_gp=gpar(fill=cl[ paste0('cpms', levels(m$c.cpm.class))]), ncol=2),
    Legend(title='Classifications', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'),
        at=c('low', 'high'), labels=c('low', 'high'), legend_gp=gpar(fill=cl[c('low', 'high')]), ncol=1), 
    column_gap=unit(0.5, 'cm'), max_height=unit(30, 'cm'))
op<-oncoPrint(mat, alter_fun=list(
        background=alter_graphic('rect', width=0.9, height=0.9, fill='#cccccc', col=NA), 
        'ST4S (13)'=alter_graphic('rect', width=0.9, height=0.9, fill=cl['ST4S (13)'], col=NA),
        'LR (30)'=alter_graphic('rect', width=0.9, height=0.9, fill=cl['LR (30)'], col=NA),
        'IMR (10)'=alter_graphic('rect', width=0.9, height=0.9, fill=cl['IMR (10)'], col=NA),
        'HR_nMNA (29)'=alter_graphic('rect', width=0.9, height=0.9, fill=cl['HR_nMNA (29)'], col=NA),
        'MNA (22)'=alter_graphic('rect', width=0.9, height=0.9, fill=cl['MNA (22)'], col=NA), 
        low=alter_graphic('rect', width=0.9, height=0.9, fill=cl['low'], col=NA), 
        high=alter_graphic('rect', width=0.9, height=0.9, fill=cl['high'], col=NA), 
        young=alter_graphic('rect', width=0.9, height=0.9, fill=cl['young'], col=NA), 
        old=alter_graphic('rect', width=0.9, height=0.9, fill=cl['old'], col=NA), 
        OS0=alter_graphic('rect', width=0.9, height=0.9, fill=cl['OS0'], col=NA), 
        OS1=alter_graphic('rect', width=0.9, height=0.9, fill=cl['OS1'], col=NA), 
        EFS0=alter_graphic('rect', width=0.9, height=0.9, fill=cl['EFS0'], col=NA), 
        EFS1=alter_graphic('rect', width=0.9, height=0.9, fill=cl['EFS1'], col=NA), 
        amp=alter_graphic('rect', width=0.9, height=0.9, fill=cl['amp'], col=NA), 
        'non-amp'=alter_graphic('rect', width=0.9, height=0.9, fill=cl['non-amp'], col=NA), 
        INSS1=alter_graphic('rect', width=0.9, height=0.9, fill=cl['INSS1'], col=NA), 
        INSS2=alter_graphic('rect', width=0.9, height=0.9, fill=cl['INSS2'], col=NA), 
        INSS3=alter_graphic('rect', width=0.9, height=0.9, fill=cl['INSS3'], col=NA), 
        INSS4=alter_graphic('rect', width=0.9, height=0.9, fill=cl['INSS4'], col=NA), 
        INSS5=alter_graphic('rect', width=0.9, height=0.9, fill=cl['INSS5'], col=NA),
        counts500=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts500'], col=NA),
        counts800=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts800'], col=NA),
        counts1200=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts1200'], col=NA),
        counts2000=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts2000'], col=NA),
        counts2500=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts2500'], col=NA),
        counts2800=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts2800'], col=NA),
        counts3100=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts3100'], col=NA),
        counts3400=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts3400'], col=NA),
        counts3700=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts3700'], col=NA),
        counts4000=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts4000'], col=NA),
        counts4300=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts4300'], col=NA),
        counts5100=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts5100'], col=NA),
        cpms2200=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms2200'], col=NA),
        cpms2400=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms2400'], col=NA),
        cpms2600=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms2600'], col=NA),
        cpms2800=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms2800'], col=NA),
        cpms3000=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms3000'], col=NA),
        cpms3200=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms3200'], col=NA),
        cpms3400=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms3400'], col=NA),
        cpms3600=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms3600'], col=NA),
        cpms3800=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms3800'], col=NA),
        cpms4200=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms4200'], col=NA),
        cpms4600=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms4600'], col=NA),
        cpms6400=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms6400'], col=NA),
        cpms12600=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms12600'], col=NA)
    ), alter_fun_is_vectorized=T, col=cl, show_pct=F, show_row_names=T, row_names_gp=gpar(fontsize=20), show_column_names=F,
    top_annotation=NULL, right_annotation=NULL,
    heatmap_legend_param=NULL,
    show_heatmap_legend=F,
#    width=unit(1.0, 'npc'),              #  size of heatmap body only
#    heatmap_width=unit(1.0, 'npc'),      #  size of full heatmap including components
    column_split=m[, risk_group],
    column_gap=unit(2, 'mm'),
    column_title_gp=gpar(fontsize=18)
)
op<-draw(op, annotation_legend_list=lg)
dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/oncoPlot_classes_tumors.svg', width=16, height=7.0, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [including MYC/MYCN mRNA expression] ordered by risk group primarily and circRNA numbers secondarily
#{{{

#  order by circRNA numbers, the risk group ordering is enforced when we ask for column splits but you have to define the factor level order
meta.tum<-meta.tum[, risk_group:=factor(risk_group, levels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))][ order(circRNAcounts) ]


mat<-t(data.frame(meta.tum[, c('risk_group', 'INSS_STAGE', 'OS_bin', 'EFS_bin', 'AGE_bin', 'MYCN_INI_STATUS', 'mycn.class', 'myc.class', 'PI.class', 'circRNAcounts.class', 'c.cpm.class'), with=F][, c('INSS_STAGE', 'OS_bin', 'EFS_bin', 'circRNAcounts.class', 'c.cpm.class'):=list(paste0('INSS', INSS_STAGE), paste0('OS', OS_bin), paste0('EFS', EFS_bin), paste0('counts', circRNAcounts.class), paste0('cpms', c.cpm.class))]))
rownames(mat)<-c('Risk group', 'INSS stage', 'OS', 'EFS', 'Age', 'MYCN status', 'MYCN mRNA', 'MYC mRNA', 'PI', 'circRNA number', 'circRNA expression')
cl<-setNames( unique(meta.tum$col), unique(meta.tum$risk_group) )
cl<-append(cl, setNames( colorRampPalette(brewer.pal(9, 'Blues'))(12)[c(seq(6, 12, 2), 4)], c('INSS1', 'INSS2', 'INSS3', 'INSS4', 'INSS5') ) )
cl<-append(cl, setNames( c('#cccccc', '#b21f1f'), c('OS0', 'OS1') ) )
cl<-append(cl, setNames( c('#cccccc', '#b21f1f'), c('EFS0', 'EFS1') ) )
cl<-append(cl, setNames( unique(meta.tum$AGE_bin.col), unique(meta.tum$AGE_bin) ) )
cl<-append(cl, setNames( c('#cccccc', '#b21f1f'), c('non-amp', 'amp') ) )
cl<-append(cl, setNames( levels(meta.tum$mycn.class.col), levels(meta.tum$mycn.class) ) )
cl<-append(cl, setNames( rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(meta.tum$circRNAcounts.class)))), paste0('counts', levels(meta.tum$circRNAcounts.class)) ) )
cl<-append(cl, setNames( rev(colorRampPalette(brewer.pal(9, 'RdYlGn'))(length(levels(meta.tum$c.cpm.class)))), paste0('cpms', levels(meta.tum$c.cpm.class)) ) )
N<-setNames(as.numeric(levels(meta.tum$circRNAcounts.class)), levels(meta.tum$circRNAcounts.class))
N<-N[c(1, floor(length(N)/2), length(N))]
#names(N)[seq(2, length(N), 2)]<-''
E<-setNames(as.numeric(levels(meta.tum$c.cpm.class)), levels(meta.tum$c.cpm.class))
E<-E[c(1, floor(length(E)/2), length(E))]
#names(E)[seq(2, length(E), 2)]<-''
lg<-packLegend(
    #Legend(title='Risk group', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'),
    #    at=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'), labels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'), 
    #    legend_gp=gpar(fill=cl[c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA')]), ncol=1),
    Legend(title='INSS stage', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'),
        at=c('INSS1', 'INSS2', 'INSS3', 'INSS4', 'INSS5'), labels=c('1', '2', '3', '4', '5'), 
        legend_gp=gpar(fill=cl[c('INSS1', 'INSS2', 'INSS3', 'INSS4', 'INSS5')]), ncol=1),
    Legend(title='OS/EFS event', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'), at=c('OS0', 'OS1'), 
        labels=c('0', '1'), legend_gp=gpar(fill=cl[c('OS0', 'OS1')]), ncol=1),
    Legend(title='Age', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'), at=c('young', 'old'),
        labels=c('<547 days', '>546 days'), legend_gp=gpar(fill=cl[c('young', 'old')]), ncol=1),
    Legend(title='MYCN status', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'),
        at=c('non-amp', 'amp'), labels=c('not amplified', 'amplified'), legend_gp=gpar(fill=cl[c('non-amp', 'amp')]), ncol=1),
    Legend(title='circRNA number', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), 
        grid_height=unit(2,'mm'), at=N, labels=names(N), legend_gp=gpar(fill=cl[paste0('counts', N)]), ncol=1),
    Legend(title='circRNA expression', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14),
        grid_height=unit(2,'mm'), at=E, labels=names(E), legend_gp=gpar(fill=cl[paste0('cpms', E)]), ncol=1),
    #Legend(title='circRNA number', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), 
    #    grid_height=unit(4,'mm'), at=as.numeric(levels(meta.tum$circRNAcounts.class)), labels=levels(meta.tum$circRNAcounts.class), 
    #    legend_gp=gpar(fill=cl[ paste0('counts', levels(meta.tum$circRNAcounts.class))]), ncol=2),
    #Legend(title='circRNA expression', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14),
    #    grid_height=unit(4,'mm'), at=as.numeric(levels(meta.tum$c.cpm.class)), labels=levels(meta.tum$c.cpm.class),
    #    legend_gp=gpar(fill=cl[ paste0('cpms', levels(meta.tum$c.cpm.class))]), ncol=2),
    Legend(title='Classifications', title_gp=gpar(fontsize=16, fontface='bold'), labels_gp=gpar(fontsize=14), grid_height=unit(4,'mm'),
        at=c('low', 'high'), labels=c('low', 'high'), legend_gp=gpar(fill=cl[c('low', 'high')]), ncol=1), 
    column_gap=unit(0.5, 'cm'), max_height=unit(30, 'cm'))
op<-oncoPrint(mat, alter_fun=list(
        background=alter_graphic('rect', width=0.9, height=0.9, fill='#cccccc', col=NA), 
        ST4S=alter_graphic('rect', width=0.9, height=0.9, fill=cl['ST4S'], col=NA), 
        LR=alter_graphic('rect', width=0.9, height=0.9, fill=cl['LR'], col=NA), 
        IMR=alter_graphic('rect', width=0.9, height=0.9, fill=cl['IMR'], col=NA), 
        HR_nMNA=alter_graphic('rect', width=0.9, height=0.9, fill=cl['HR_nMNA'], col=NA), 
        MNA=alter_graphic('rect', width=0.9, height=0.9, fill=cl['MNA'], col=NA), 
        low=alter_graphic('rect', width=0.9, height=0.9, fill=cl['low'], col=NA), 
        high=alter_graphic('rect', width=0.9, height=0.9, fill=cl['high'], col=NA), 
        young=alter_graphic('rect', width=0.9, height=0.9, fill=cl['young'], col=NA), 
        old=alter_graphic('rect', width=0.9, height=0.9, fill=cl['old'], col=NA), 
        OS0=alter_graphic('rect', width=0.9, height=0.9, fill=cl['OS0'], col=NA), 
        OS1=alter_graphic('rect', width=0.9, height=0.9, fill=cl['OS1'], col=NA), 
        EFS0=alter_graphic('rect', width=0.9, height=0.9, fill=cl['EFS0'], col=NA), 
        EFS1=alter_graphic('rect', width=0.9, height=0.9, fill=cl['EFS1'], col=NA), 
        amp=alter_graphic('rect', width=0.9, height=0.9, fill=cl['amp'], col=NA), 
        'non-amp'=alter_graphic('rect', width=0.9, height=0.9, fill=cl['non-amp'], col=NA), 
        INSS1=alter_graphic('rect', width=0.9, height=0.9, fill=cl['INSS1'], col=NA), 
        INSS2=alter_graphic('rect', width=0.9, height=0.9, fill=cl['INSS2'], col=NA), 
        INSS3=alter_graphic('rect', width=0.9, height=0.9, fill=cl['INSS3'], col=NA), 
        INSS4=alter_graphic('rect', width=0.9, height=0.9, fill=cl['INSS4'], col=NA), 
        INSS5=alter_graphic('rect', width=0.9, height=0.9, fill=cl['INSS5'], col=NA),
        counts500=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts500'], col=NA),
        counts800=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts800'], col=NA),
        counts1200=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts1200'], col=NA),
        counts2000=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts2000'], col=NA),
        counts2500=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts2500'], col=NA),
        counts2800=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts2800'], col=NA),
        counts3100=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts3100'], col=NA),
        counts3400=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts3400'], col=NA),
        counts3700=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts3700'], col=NA),
        counts4000=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts4000'], col=NA),
        counts4300=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts4300'], col=NA),
        counts5100=alter_graphic('rect', width=0.9, height=0.9, fill=cl['counts5100'], col=NA),
        cpms2200=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms2200'], col=NA),
        cpms2400=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms2400'], col=NA),
        cpms2600=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms2600'], col=NA),
        cpms2800=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms2800'], col=NA),
        cpms3000=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms3000'], col=NA),
        cpms3200=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms3200'], col=NA),
        cpms3400=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms3400'], col=NA),
        cpms3600=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms3600'], col=NA),
        cpms3800=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms3800'], col=NA),
        cpms4200=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms4200'], col=NA),
        cpms4600=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms4600'], col=NA),
        cpms6400=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms6400'], col=NA),
        cpms12600=alter_graphic('rect', width=0.9, height=0.9, fill=cl['cpms12600'], col=NA)
    ), alter_fun_is_vectorized=T, col=cl, show_pct=F, show_row_names=T, row_names_gp=gpar(fontsize=20), show_column_names=F,
    top_annotation=NULL, right_annotation=NULL,
    heatmap_legend_param=NULL,
    show_heatmap_legend=F,
#    width=unit(1.0, 'npc'),              #  size of heatmap body only
#    heatmap_width=unit(1.0, 'npc'),      #  size of full heatmap including components
    column_split=meta.tum[, risk_group],
    column_gap=unit(2, 'mm'),
    column_title_gp=gpar(fontsize=18)
)
op<-draw(op, annotation_legend_list=lg)
#dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/oncoPlot_classes_tumors.svg', width=16, height=7.0, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}

#}}}




############################
#
#
#  checks and quick analyses
#
#
############################




#  check count correlation between kallisto and featureCounts 
#{{{
rm(list=ls())  
library(GenomicAlignments)
library(rtracklayer)
library(data.table)


#  load annotation to identify gene_id from transcript_id
hsa<-import('/data/annotation/GRCh38/GRCh38.gencode.v27.gtf')


#  load featureCount data and massage the counts to a data.table
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')
fe<-unlist(totalrna)
fe$bid<-sub('\\.[0-9]*$', '', rownames(fe))
rownames(fe)<-NULL
fe<-data.table(fe[, c('bid', 'gene_id', 'counts')])


#  load kallisto data
#  add gene_id
#  take MEAN counts across isoforms to correspond to the gene count
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/kallisto.RData')
ka<-unlist(totalrna)
ka$bid<-sub('\\.[0-9]*$', '', rownames(ka))
rownames(ka)<-NULL
ka<-data.table(ka[, c('bid', 'transcript_id', 'counts')])
ka$gene_id<-hsa$gene_id[ match(ka$transcript_id, hsa$transcript_id) ]
ka<-ka[, .(counts=mean(counts)), by=.(bid, gene_id)]


#  remove failed samples throughout
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
fe<-fe[ bid %in% meta[!(failed), bid] ]
ka<-ka[ bid %in% meta[!(failed), bid] ]


#  compute Spearman correlations
x11(width=18, height=11, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
corr<-list()
for (s in unique(ka$bid)){
    x<-ka[ bid %in% s ]
    y<-fe[ bid %in% s ]
    stopifnot( length(setdiff(x$gene_id, y$gene_id))==0 )
    y<-y[ match(x$gene_id, gene_id), ]
    plot( log10(1+x$counts), log10(1+y$counts), pch=19, xlab='featureCounts', ylab='kallisto', cex=0.8)
    corr[[s]]<-cor( x$counts, y$counts, method='spearman' )
}
corr<-unlist(corr)


#  plot correlations
x11(width=18, height=11, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.5,6,0.5,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.8, cex.axis=1.6)
b<-barplot(corr, border=NA, col='grey39', axisnames=F, ylim=c(0, 1), xlab='', ylab='')
mtext(names(corr), side=1, las=2, at=b, adj=+1.05, cex=1.2)
mtext('Spearman correlation', side=2, line=3, padj=-0.8, cex=1.8, las=3)

#}}}



#  [HDAC11] check the expression of protein-coding isoforms in certain cell lines
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/draw_highlights.R')


#  load the reference
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'transcript' ]
mcols(hsa)<-mcols(hsa)[, c('gene_name', 'gene_id', 'gene_type', 'transcript_name', 'transcript_id', 'transcript_type')]


#  keep only specific non-failed cell lines
#  identify the cell line irrespective of the cell model
#  fix the colors to be unique per cell line
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta<-meta[ !(failed) & is.na(risk_group) & !grepl('Pilote', bid) & !grepl('Rhabdomyosarcoma|PAX3 FOXO1 fusion|Differentiation', cell_model) ]
meta[, cell_line:='']
meta[ grep('-IMR5-', bid), cell_line:='IMR5']
meta[ grep('-SHEP_TR-', bid), cell_line:='SHEP-TR']
meta[ grep('-SKNAS-TR-', bid), cell_line:='SKNAS-TR']
meta<-meta[ cell_line != '' ]
cl<-setNames( colorRampPalette(c('#00008B', '#006400', '#104E8B', '#2F4F4F', '#458B74', '#556B2F', '#8B1A1A', '#E9967A', '#FF8C00'))(meta[, length(unique(cell_line))]), meta[, unique(cell_line)] )
meta[, col:=cl[ cell_line ]]
rm(cl)


#  load kallisto TPMS
#  convert to matrix
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/kallisto.RData')
tra<-do.call(rbind, lapply(meta$bid, function(n){ data.table(totalrna[[ n ]])[, bid:=n] }))
tra<-dcast(tra, bid ~ transcript_id, value.var='tpm', fun.aggregate=sum)
tra<-t(data.frame(tra[, -1], row.names=tra[, bid], check.names=F))
rm(totalrna)


#  HDAC11 well-expressed isoforms 
#{{{

PICK<-hsa$transcript_id[ grepl('HDAC11', hsa$transcript_name) & hsa$transcript_type %in% 'protein_coding' ]
x<-tra[ PICK, , drop=F]
x<-x[ rowMeans(x)>0.01, ]
x<-x[ order(rowMeans(x), decreasing=T), ]
rownames(x)<-hsa$transcript_name[ match(rownames(x), hsa$transcript_id) ]
colnames(x)<-meta[ match(colnames(x), bid), cell_line ]
B<-list()
for(n in unique(colnames(x))){ 
    for(r in rownames(x)){
        B[[paste(n, r, sep='_')]]<-x[r, colnames(x) %in% n, drop=T]
    }
}
B.cl<-colorRampPalette(c('#00008B', '#006400', '#104E8B', '#2F4F4F', '#458B74', '#556B2F', '#8B1A1A', '#E9967A', '#FF8C00'))(nrow(x))


#  boxplot
x11(width=20, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(0.0, max(unlist(B))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('TPM', side=2, line=5, padj=-0.1, las=0, cex=2.4)
mtext(text=sub('^.*HDAC11-([0-9]+)$', '\\1', names(B)), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.0, col=B.cl)
mtext(text=unique(colnames(x)), side=1, at=seq(nrow(x)/2, length(B), nrow(x)), line=4, las=0, padj=-0.1, cex=2.4, col='black')
draw_highlights(L=length(B), STEP=nrow(x), YMAX=max(YTICK), YMIN=0.0)
legend('topright', legend=rownames(x), col=B.cl, bty='n', lty=1, lwd=15, pch=NA, cex=1.8, y.intersp=0.6, x.intersp=0.5, seg.len=0.5)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_cell_lines_TPMs_for_HDAC11_transcripts.svg'), width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}



#  [MYCN expression throughout the tumors] 
#{{{
rm(list=ls())  
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  fix the ggplot2 theme
#{{{
theme_set(theme_bw(base_size=35) + theme(panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_blank(),
    plot.margin=margin(t=1.0, b=0.1, l=0.5, r=1.0, unit='cm'),
    plot.title=element_text(hjust=0.5), 
    axis.line=element_line(color='black', size=0.5),
    axis.text=element_text(color='black', size=40, family='Arial', face='plain'), 
    axis.text.x=element_text(angle=90, margin=margin(t=0.2, b=0, unit='cm')),
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


#  load annotation to identify gene_id from transcript_id
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_name', 'gene_type', 'level')]


#  load metadata and keep tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-meta.tum[ !(failed) ]
rm(meta.tum, meta.cel, meta.prefailed)


#  load featureCount data 
#  unlist them to data.table
#  add TPMs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')
fe<-unlist(totalrna[ meta$bid ])
fe$bid<-sub('\\.[0-9]*$', '', rownames(fe))
rownames(fe)<-NULL
fe<-data.table(fe[, c('bid', 'gene_id', 'length', 'counts')])
fe<-fe[, tpm:=1e6*counts/length/sum(counts/length), by=.(bid)]


#  select MYCN
G<-'MYCN'
G<-setNames(hsa$gene_id[ hsa$gene_name %in% G ], G)


#  sort by TPM
#  beautify samples names
#  assign colors to samples by MNA status but single out CB2028, CB3016, CB3036
B<-data.frame(fe[ gene_id %in% G, ])[, c('bid', 'tpm')]
B<-B[ order(B$tpm, decreasing=T), , drop=F]
B$status<-meta[ match(B$bid, bid), MYCN_INI_STATUS]
B$bid<-sub('-[0-9]*.*$', '', B$bid)
B$bid<-factor(B$bid, levels=B$bid)
B$col<-'grey39'
B$col[ B$status %in% 'amp' ]<-'coral4'
B$col[ B$bid %in% c('CB2028', 'CB3016', 'CB3036') ]<-'darkgreen'
B$status[ B$bid %in% c('CB2028', 'CB3016', 'CB3036') ]<-'?'
cl<-setNames(unique(B$col), unique(B$status))


#  barplot
x11(width=22, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 
YTICK<-pretty(c(0, max(B[, 'tpm'])), 5)
YMAX<-tail(YTICK, 1)
ggplot(B, aes(bid, tpm, fill=status)) + 
    geom_bar(stat='identity', position='identity', aes(fill=status)) +
    scale_fill_manual(values=cl) +
    labs(title='MYCN expression') + 
    ylab('TPM') +
    scale_y_continuous(limits=c(0, YMAX), breaks=YTICK, expand=expansion(mult=c(0.01, 0))) +
    theme(axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text=element_text(size=48),
          axis.text.x=element_text(angle=90, vjust=0.5, hjust=1.0, color=B$col, size=40),
          axis.title.x=element_blank())
ggsave('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_MYCN.svg', device=svg, width=49, height=14, scale=1.0, bg='white', antialias='subpixel', family='Arial')

#}}}



#  Mass spectrometry plots
#{{{
rm(list=ls())
library(data.table)
library(openxlsx)
library(GenomicFeatures)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load the non-failed mass-spec results
ms<-read.xlsx('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/validations/rbp/mass_spec/Result_Pull_down_circRNA_SF_20200123.xlsx', cols=c(5, 24:28))
colnames(ms)<-c('gene_names', '-log10.p.value', 'q.value', 'diff', 'stat')
ms<-data.table(ms)[, c('gene_names', 'p.value'):=list(strsplit(gene_names, ';'), 10^(-`-log10.p.value`))][ order(p.value) ]


#  volcano plot
x11(width=12, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.5, 8.0, 0.1, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4)
#XLIM<-range(pretty(ms[, `diff`], 5))
XLIM<-c(-4, 4)
YLIM<-range(pretty(ms[, `-log10.p.value`], 5))
plot(ms[, `diff`], ms[, `-log10.p.value`], type='p', pch=21, xlim=XLIM, ylim=YLIM, main='', xlab='', ylab='', lwd=0, cex=1.2, col=adjustcolor('grey60', alpha.f=0.7), bg=adjustcolor('grey60', alpha.f=0.7))
points(ms[ p.value<0.05, `diff`], ms[ p.value<0.05, `-log10.p.value`], pch=21, lwd=1, cex=1.2, col='red3', bg='red3')
abline(v=0, lty=3, lwd=4, col='grey50')
mtext(expression(-log[10]('P-value')), side=2, line=5, padj=+0.2, cex=2.4, las=3)
mtext('Intensity fold-change', side=1, line=4, padj=-0.3, cex=2.4, las=1)
r<-ms[ p.value<0.05 ][, gene_name:=sapply(gene_names, '[[', 1)]
#text(r[, `diff`], r[, `-log10.p.value`], labels=r[, `gene_name`], adj=c(-0.15, -0.15), cex=1.2, col='red3')
identify(r[, `diff`], r[, `-log10.p.value`], labels=r[, `gene_name`], offset=0.35, cex=1.2, col='red3', xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/volcanoPlot_mass-spec_pull_down_circRNA_20200123.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}




###########
#
#
#  obsolete
#
#
###########




#  [DCC] collect results
#        keep circRNAs with at least 5 reads covering the junction in at least one sample
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(openxlsx)


#  functions
#{{{

parse_results<-function(B=''){
    require(data.table)

    #  load CircCoordinates
    n<-tryCatch({
        cr<-fread(B, sep='\t', col.names=c('seqnames', 'start', 'end', 'gene', 'junction_type', 'strand', 'region', 'whole_region'))
        if(nrow(cr)==0){
            return(GRanges())
        }


        #  load CircRNACount and LinearCount and add them
        x<-fread(sub('CircCoordinates', 'CircRNACount', B), sep='\t', col.names=c('seqnames', 'start', 'end', 'circ_count'))
        y<-fread(sub('CircCoordinates', 'LinearCount', B), sep='\t', col.names=c('seqnames', 'start', 'end', 'linear_count'))
        stopifnot( all.equal( cr[, 1:3] , x[, 1:3] ) )
        stopifnot( all.equal( cr[, 1:3] , y[, 1:3] ) )
        x$linear_count<-y$linear_count
        cr<-cbind(x, cr[, -c(1:3)])    #  add circ_counts, linear_counts right after coordinates


        #  load CircSkipJunctions and add them
        x<-fread(sub('CircCoordinates', 'CircSkipJunctions', B), sep='\t', skip=1, col.names=c('seqnames', 'start', 'end', 'strand', 'skipJ'))
        stopifnot( all.equal( cr[, 1:3] , x[, 1:3] ) )
        cr$skipJ<-x$skipJ

        
        #  break comma separated whole_region annotations, sort them, and join them again so that they are consistent throughout
        cr[, whole_region:=sapply(strsplit(whole_region, ','), function(x){ paste0(sort(x), collapse=',') }) ]

        
        #  return GRanges object
        return(GRanges(as.data.frame(cr)))

    }, error=function(e){
        warning(e)
        return(GRanges())
    }) 
}

#}}}


#  locate the DCC results from the non-failed samples
cir<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/ -maxdepth 3 -type f -wholename \'*/dcc/CircCoordinates\' -print', stdout=T)


#  add sequencing-sample-id (bid)
names(cir)<-sub('^.*(JS_SF_[0-9]+).*$', '\\1', cir)


#  order them
cir<-cir[ order(names(cir)) ]


#  append results
circ<-List()
for (n in seq_along(cir)){
    cat('\nprocessing: ', cir[n], '\n')
    circ[[ names(cir)[n] ]]<-parse_results(cir[n])
}
rm(n)


#  convert from List to GRanges and add bid to metadata
circ<-unlist(GRangesList(lapply(circ,c)))
circ$bid<-names(circ)
names(circ)<-NULL


#  keep circRNAs with at least 5 reads covering the junction in at least one sample
circ<-data.table(as.data.frame(circ))
circ[, pass:=any(circ_count>=5), by=.(seqnames, start, end, strand)] 
circ<-circ[ pass %in% TRUE, ][, pass:=NULL]
circ<-GRanges(seqnames=circ$seqnames, strand=circ$strand, ranges=IRanges(start=circ$start, end=circ$end), data.frame(circ[, -c(1:5)]))


#  save
save(circ, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_DCC.RData')
#l<-circ
#x<-cbind(data.frame(position=paste0(as.character(seqnames(l)), '(', as.character(strand(l)),'):', start(l), '-', end(l))), as.data.frame(l)[, -c(1:3,5)])
#write.xlsx(x, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_DCC.xlsx', col.names=T, row.names=F, sheetName='circRNAs', append=F)

#}}}



#  comparison of DCC and CIRI2
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  functions
#{{{
dcc_ciri_overlaps<-function(S=''){
    
    x<-dcc[ dcc$bid %in% S ]
    y<-ciri[ ciri$bid %in% S ]
    ov<-findOverlaps(x, y, type='equal', select='all', ignore.strand=F)

    dcc_ciri<-length(unique(queryHits(ov)))
    dcc_only<-length(x) - dcc_ciri
    ciri_only<-length(y) - dcc_ciri

    return(data.frame(dcc_ciri=dcc_ciri, dcc_only=dcc_only, ciri_only=ciri_only, row.names=S))

}

#}}}


#  load DCC predictions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_DCC.RData')
dcc<-circ
mcols(dcc)<-mcols(dcc)[, c('circ_count', 'linear_count', 'gene', 'junction_type', 'region', 'whole_region', 'bid')]         
seqlevels(dcc, pruning.mode='coarse')<-seqlevels(dcc)[ grep('chr', seqlevels(dcc)) ]
rm(circ)


#  load CIRI2 predictions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2.RData')
ciri<-circ
rm(circ)


#  load sample metadata
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)


#  identify common (identical) circRNA calls and either non-identical or not common 
#{{{

s<-sort(union(dcc$bid, ciri$bid))
ov<-dcc_ciri_overlaps(s[1])
for (n in 2:length(s)){
  ov<-rbind(ov, dcc_ciri_overlaps(s[n]))
}


#  CLICK on it once to make sure it does not redraw, or options(scipen=-20) might be IGNORED
x11(width=22, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 


#  colors related to the statistics
cl<-data.frame(symbol=c('common', 'DCC only', 'CIRI2 only'), 
               color=c('#666666',  #  grey40
                       '#A66753',  
                       '#1B9E77')) 


#  barplot
par(mar=c(5.5,4.5,1.5,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.4, cex.axis=1.4)
YMAX<-rowSums(ov) - rowSums(ov)%%1e3 + 1e3
YTICK<-pretty(c(0, YMAX), 5)
bp<-barplot(t(ov), border='white', col=cl$color, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=range(YTICK))
axis(2, at=YTICK, line=-1, cex.axis=1.2)
mtext(rownames(ov), side=1, line=0, at=bp, las=2, adj=1.01, cex=1.0)
mtext('circRNAs predicted', side=2, line=3, padj=-0.2, las=0, cex=1.4)
legend(x=par('usr')[2]*0.81, y=1.03*par('usr')[4], legend=cl$symbol , col=cl$color, bty='n', lty=1, lwd=8, cex=1.2)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_DCC_CIRI2_comparison.svg', width=22, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(cl,YMAX,YTICK,bp)

#}}}

#}}}



#  [tumors, cell lines] calculate the proliferative index (BASED ON VST COUNTS)
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(ProliferativeIndex)
library(pheatmap)
library(cluster)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_ecdfs.R')


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


#  split metadata to cell lines and tumors
tum<-meta[ !is.na(risk_group) ]
cel<-meta[ is.na(risk_group) ]
stopifnot( nrow(tum)+nrow(cel)==nrow(meta) )
rm(meta, meta.prefailed)


#  load counts of kept samples 
#  keep only the genes found in the trimmed reference
#  add gene_names
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')
totalrna<-unlist(totalrna[ union(tum$bid, cel$bid) ])
totalrna$bid<-sub('\\.[0-9]*$', '', rownames(totalrna))
rownames(totalrna)<-NULL
totalrna<-totalrna[ totalrna$gene_id %in% hsa$gene_id, , drop=F]
totalrna$gene_name<-hsa$gene_name[ match(totalrna$gene_id, hsa$gene_id ) ]


#  split to cell lines and tumors 
#  convert to matrices by summarizing gene_ids by gene_name 
#  order according to metadata
tum.cnt<-dcast(data.table(totalrna[ totalrna$bid %in% tum$bid, ]), bid ~ gene_name, value.var='counts', fun.aggregate=sum)
tum.cnt<-ceiling(t(data.frame(tum.cnt[, -1], row.names=tum.cnt[, bid], check.names=F)))
cel.cnt<-dcast(data.table(totalrna[ totalrna$bid %in% cel$bid, ]), bid ~ gene_name, value.var='counts', fun.aggregate=sum)
cel.cnt<-ceiling(t(data.frame(cel.cnt[, -1], row.names=cel.cnt[, bid], check.names=F)))
stopifnot( length( setdiff( totalrna$bid, union(colnames(tum.cnt), colnames(cel.cnt)) ) )==0 )
tum.cnt<-tum.cnt[, tum$bid]
cel.cnt<-cel.cnt[, cel$bid]
rm(totalrna)


#  apply variance-stabilizing transformations
tum.vst<-varianceStabilizingTransformation(tum.cnt, fitType='local', blind=T)
cel.vst<-varianceStabilizingTransformation(cel.cnt, fitType='local', blind=T)


#  keep only metaPCNA2 signature genes
tum.cnt<-tum.cnt[ rownames(tum.cnt) %in% ProliferativeIndex:::metaPCNA2, , drop=F]
tum.vst<-tum.vst[ rownames(tum.vst) %in% ProliferativeIndex:::metaPCNA2, , drop=F]
cel.cnt<-cel.cnt[ rownames(cel.cnt) %in% ProliferativeIndex:::metaPCNA2, , drop=F]
cel.vst<-cel.vst[ rownames(cel.vst) %in% ProliferativeIndex:::metaPCNA2, , drop=F]


#  add proliferative index to metadata
tum$PI<-calculatePI(list(vstData=tum.vst))
cel$PI<-calculatePI(list(vstData=cel.vst))


#  [tumors] cluster samples by PI
#           cut at 3 clusters 
#           make sure to order them by low, moderate, high
#           add discrete PI classification to the metadata
x<-as.matrix(setNames(tum$PI, tum$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
summary(silhouette(cutree(x.hc, k=2), dist=x.d))  #  two clusters have slightly higher mean silhouette width but hclust separates good at 3 as well
summary(silhouette(cutree(x.hc, k=3), dist=x.d))
w<-cutree(x.hc, k=3)
o<-sapply(split(names(w), w), function(n){ mean(x[n, ,drop=F]) })
for(n in names(o)){
    w[ w==n ]<-o[n]
}
w<-sort(w, decreasing=F)
x<-x[ names(w), ,drop=F ]
x.d<-dist(x, method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
w<-cutree(x.hc, k=3)
w[ w %in% '1' ]<-'low'
w[ w %in% '2' ]<-'moderate'
w[ w %in% '3' ]<-'high'
tum$PI.class<-factor(w[ tum$bid ], levels=c('low', 'moderate', 'high'))
tum$PI.class.col<-factor(w[ tum$bid ], levels=c('low', 'moderate', 'high'), labels=c('#cccccc', '#bfb610', '#b21f1f'))


#  [cell lines] cluster samples by PI
#               cut at 3 clusters 
#               make sure to order them by low, moderate, high
#               add discrete PI classification to the metadata
x<-as.matrix(setNames(cel$PI, cel$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
summary(silhouette(cutree(x.hc, k=2), dist=x.d))  #  two clusters have slightly higher mean silhouette width but hclust separates good at 3 as well
summary(silhouette(cutree(x.hc, k=3), dist=x.d))
w<-cutree(x.hc, k=3)
o<-sapply(split(names(w), w), function(n){ mean(x[n, ,drop=F]) })
for(n in names(o)){
    w[ w==n ]<-o[n]
}
w<-sort(w, decreasing=F)
x<-x[ names(w), ,drop=F ]
x.d<-dist(x, method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
w<-cutree(x.hc, k=3)
w[ w %in% '1' ]<-'low'
w[ w %in% '2' ]<-'moderate'
w[ w %in% '3' ]<-'high'
cel$PI.class<-factor(w[ cel$bid ], levels=c('low', 'moderate', 'high'))
cel$PI.class.col<-factor(w[ cel$bid ], levels=c('low', 'moderate', 'high'), labels=c('#cccccc', '#bfb610', '#b21f1f'))


#  save
save(tum, tum.cnt, tum.vst, cel, cel.cnt, cel.vst, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/proliferative_index.RData')


#  visualize the classification results
#{{{

#  load back the results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/proliferative_index.RData')


#  functions
#{{{
markers_boxplot<-function(PICKED, COLS, YMIN=0, fileroot=''){
    #  quick and dirty boxplot for the tumors

    B<-setNames(lapply(names(COLS), function(r){ c(PICKED[, colnames(PICKED) %in% r]) }), names(COLS))
    YTICK<-pretty(c(YMIN, sapply(B, max)), 5)
    plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
    bp<-boxplot(B, col=COLS, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=COLS, range=0, add=T)
    mtext('PI', side=2, line=4, padj=-0.1, las=0, cex=2.4)
    mtext(text=paste0(names(B), ' (', lengths(B)/nrow(x), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=COLS)
    mtext(text=paste(rownames(PICKED), sep='', collapse=','), side=3, line=-1, padj=+0.2, cex=1.8)

    if (fileroot!=''){
        dev.print(device=svg, file=paste0(fileroot, paste(rownames(PICKED), sep='', collapse=','), '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    }
}

#}}}


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  [tumors] boxplot of PI values by PI class
par(mar=c(12.0, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(tum$PI, tum$PI.class)
B.cl<-setNames(levels(tum$PI.class.col), levels(tum$PI.class))
YTICK<-pretty(c(7, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('PI', side=2, line=4, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_PI_per_proliferative_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [tumors] heatmap of annotated samples clustered by Euclidean distance in variance-stabilized counts of PI-related genes
x<-tum.vst
colnames(x)<-tum$bid
x.ex<-data.frame(Risk=tum$risk_group, PI=tum$PI.class, row.names=tum$bid)
x.cl<-setNames(list(setNames( unique(tum$col), unique(x.ex$Risk)), setNames(levels(tum$PI.class.col), levels(tum$PI.class))), colnames(x.ex))
x.d<-dist(t(x), method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
svg('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_proliferative_index_tumors.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=18, family='Arial')  #  workaround the cutting of legends?
ph<-pheatmap(as.matrix(x.d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=x.hc,
        cluster_cols=x.hc,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=x.ex, annotation_row=x.ex, annotation_colors=x.cl,
        drop_levels=F, show_rownames=F, show_colnames=F, 
        display_numbers=F, number_format='%.1f', number_color='grey39',
        fontsize=18, fontsize_row=18, fontsize_col=18, fontsize_number=10)
print(ph)
dev.off()
#dev.print(device=svg, file='/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_proliferative_index_tumors.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x, x.ex, x.cl, x.hc, ph)


#  [tumors] PCNA
par(mar=c(14.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-tum.vst['PCNA', , drop=F]
colnames(x)<-tum$risk_group
x.cl<-setNames(tum[, unique(col)], tum[, unique(risk_group)])
markers_boxplot(PICKED=x, COLS=x.cl, YMIN=7, fileroot='')

#}}}

#}}}



#  [CIRI2] fraction of genes producing circRNAs for tumors and cell lines
#          number of genes producing circRNAs stratified by risk_group and cell_model
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/trim_text.R')


#  load annotation
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')       
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_name', 'gene_type')]


#  load circRNAs 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2.RData')


#  remove failed samples and Pilot samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta<-meta[ !(failed) & !grepl('CBPilot', bid) ]
circ<-circ[ circ$bid %in% meta$bid ]
meta<-meta[ bid %in% unique(circ$bid) ]
rm(meta.tum, meta.cel, meta.prefailed, trans)


#  table of gene types with circRNAs
sort(table(hsa$gene_type[ hsa$gene_id %in% unique(circ$gene_id) ]), decreasing=T)
# 
#                     protein_coding                            lincRNA                          antisense transcribed_unprocessed_pseudogene 
#                               6733                                188                                120                                107 
#               processed_transcript             unprocessed_pseudogene               processed_pseudogene   transcribed_processed_pseudogene 
#                                 63                                 27                                 16                                 15 
#     transcribed_unitary_pseudogene                     sense_intronic      bidirectional_promoter_lncRNA                  sense_overlapping 
#                                 14                                  5                                  4                                  3 
#                           misc_RNA                             snoRNA                          IG_V_gene                       macro_lncRNA 
#                                  2                                  2                                  1                                  1 
#                                TEC 
#                                  1 


#  sum circRNA expression across isoforms per gene_id and per sample
#  convert to a matrix with (bid, gene_name_jc_count) structure by summing over identical gene_names but different gene_ids (e.g. CDR1)
circ<-data.table(data.frame(mcols(circ)[, c('gene_id', 'gene_name', 'jc_count', 'non_jc_count', 'bid')]))
circ<-circ[, .(jc_count=sum(jc_count), non_jc_count=sum(non_jc_count)), by=.(bid, gene_id, gene_name)]
circ<-dcast(circ, bid ~ gene_name, value.var='jc_count', fun.aggregate=sum)


#  fraction of genes producing circRNAs vs the fraction of cell lines/tumors
#{{{

#  isolate tumors and cell lines and remove circRNAs completely unexpressed within each group
cel<-circ[ bid %in% meta[ !is.na(cell_model), bid ] ]
pat<-circ[ bid %in% meta[ !is.na(risk_group), bid ] ]
stopifnot( nrow(circ)==nrow(cel)+nrow(pat) )
cel<-cel[ , c(1, 1+which(!sapply(cel[, -1], sum)==0)), with=F]  #  columns need to be referenced by number since we remove the bid column temporarily
pat<-pat[ , c(1, 1+which(!sapply(pat[, -1], sum)==0)), with=F]  #  columns need to be referenced by number since we remove the bid column temporarily
pat<-as.matrix(pat[, -1][ , lapply(.SD, function(s){ sum(s!=0)/length(s) })])[1, ]
cel<-as.matrix(cel[, -1][ , lapply(.SD, function(s){ sum(s!=0)/length(s) })])[1, ]


#  1 - cumulative distribution
x11(width=11, height=11, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 
par(mar=c(4.5,6.0,0.5,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, cex.lab=1.8, cex.axis=1.8, xpd=F, bty='n')
h<-curve(1-ecdf(pat)(x), from=-0.01, to=1.0, n=min(length(pat), 20), ylab='', xlab='', pch=NA, col='blue4', lty=1, lwd=8, main='', ylim=c(0, 1), xlim=c(0, 1), xaxt='n')
curve(1-ecdf(cel)(x), from=-0.01, to=1.0, n=min(nrow(cel), 10), pch=NA, col='grey39', lty=1, lwd=8, main='', xaxt='n', add=T)
axis(1, at=pretty(c(0, max(h$x))), labels=pretty(c(0, max(h$x))))
mtext('Fraction of genes', side=2, line=3, padj=-1.0, cex=1.8, las=0)
mtext('Fraction of samples', side=1, line=3, las=0, padj=+0.2, cex=1.8)
legend('topright', legend=c('tumors', 'cell lines'), col=c('blue4', 'grey39'), bty='n', lty=1, lwd=8, cex=1.6)
#lines(x=c(0.2, 0.2), y=c(-0.1, 1), lty=2, lwd=2, col='red3', xpd=F)
#dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/ecdf_circRNAs_CIRI2_frequencies_across_groups.pdf', width=14, height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/ecdf_circRNAs_CIRI2_frequencies_across_groups.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(h)
dev.off()

#}}}


#  number of genes producing circRNAs in cell lines/tumors
#{{{

#  isolate tumors and cell lines and remove circRNAs completely unexpressed within each group
cel<-circ[ bid %in% meta[ !is.na(cell_model), bid ] ]
pat<-circ[ bid %in% meta[ !is.na(risk_group), bid ] ]
stopifnot( nrow(circ)==nrow(cel)+nrow(pat) )
cel<-cel[ , c(1, 1+which(!sapply(cel[, -1], sum)==0)), with=F]  #  columns need to be referenced by number since we remove the bid column temporarily
pat<-pat[ , c(1, 1+which(!sapply(pat[, -1], sum)==0)), with=F]  #  columns need to be referenced by number since we remove the bid column temporarily
pat<-setNames(apply(as.data.frame(pat[, -1]), 1, function(x){ sum(x!=0) }), pat$bid)
cel<-setNames(apply(as.data.frame(cel[, -1]), 1, function(x){ sum(x!=0) }), cel$bid)


#  isolate metadata as well and reorder samples to follow metadata ordering
pat.meta<-meta[ !is.na(risk_group), ]
cel.meta<-meta[ !is.na(cell_model), ]
pat<-pat[ pat.meta$bid ]
cel<-cel[ cel.meta$bid ]


#  recycle
x11(width=25, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  [tumors] barplot of number of circRNA predictions
par(mar=c(6.5, 9.5, 4.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-pat
B.cl<-setNames(pat.meta$col, pat.meta$risk_group)
YTICK<-pretty(c(0, max(B)), 5)
bp<-barplot(B, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, tail(YTICK, 1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col='white', axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
for(i in seq_along(B)){ 
    b<-B
    b[-i]<-NA    #  remove all values but the current 
    barplot(b, border='white', col=B.cl[i], axes=F, axisnames=F, beside=F, yaxt='n', add=T)
}
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext(text=sub('-11-R01$', '', names(B)), side=1, line=0, at=bp, las=2, adj=1, cex=1.8, col=B.cl)
mtext('Number of genes producing circRNAs', side=2, line=7, padj=-0.1, las=0, cex=2.4)
legend(x=0.5*par('usr')[1], y=par('usr')[4]*1.01, legend=unique(names(B.cl)), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
#dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_circRNAs_CIRI2_number_of_predictions_tumors.pdf', width=40, height=12, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_circRNAs_CIRI2_number_of_predictions_tumors.svg', width=40, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [cell lines] barplot of number of circRNA predictions
par(mar=c(10.5, 9.5, 8.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-cel
B.cl<-setNames(cel.meta$col, cel.meta$cell_model)
YTICK<-pretty(c(0, max(B)), 5)
bp<-barplot(B, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, tail(YTICK, 1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col='white', axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
for(i in seq_along(B)){ 
    b<-B
    b[-i]<-NA    #  remove all values but the current
    barplot(b, border='white', col=B.cl[i], axes=F, axisnames=F, beside=F, yaxt='n', add=T)
}
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext(text=trim_text(names(B), 12), side=1, line=0, at=bp, las=2, adj=1, cex=1.8, col=B.cl)
mtext('Number of genes producing circRNAs', side=2, line=7, padj=-0.1, las=0, cex=2.4)
legend(x=0.5*par('usr')[1], y=par('usr')[4]*1.45, legend=unique(names(B.cl)), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
#dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_circRNAs_CIRI2_number_of_predictions_cell_lines.pdf', width=40, height=12, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/figures/barplot_circRNAs_CIRI2_number_of_predictions_cell_lines.svg', width=40, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}

#}}}




