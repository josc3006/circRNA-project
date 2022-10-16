###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




#######################################################
#
#
#  applied to original samples names (CH_JS_AM_[0-9]*):
#
#      create subfolder structures and symbolic links 
#
#{{{
rm(list=ls())


#  locate the first mates  
r1<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/FASTQ -maxdepth 1 -type f -wholename "*_R1_*.fastq.gz"', stdout=T)


#  identify the (provisional) subfolder names
names(r1)<-sub('(CH_JS_AM_[0-9]*)_.*$', '\\1', basename(r1))


#  create the subfolders and place symbolic links to the two mates in them
for(n in seq_along(r1)){
    d<-paste0(dirname(dirname(r1)[n]), '/', names(r1)[n])
    system2('mkdir', args=c('-p', d), stdout=T)
    system2('ln', args=c('-sf', r1[n], paste0(d, '/r1.fastq.gz')), stdout=T)
    system2('ln', args=c('-sf', sub('_R1_', '_R2_', r1[n]), paste0(d, '/r2.fastq.gz')), stdout=T)
}

#}}}
#
#
#######################################################




########################################################################################################################################################
#
#
#  applied to original samples names (CH_JS_AM_[0-9]*):
#
#      create the Makefile configuration
#      run the pipeline
#
#      manually copy the Makefile configuration from the parent directory and update all LIBS:= entries
#
#          sed '/#  raw read libraries/q' ../totalrna.conf > totalrna.conf
#          echo -ne 'LIBS:=' >> totalrna.conf
#          find /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/CH* -maxdepth 0 -type d -printf " %p" >> totalrna.conf
#          ln -s ../../lib/qsub.mf .
#
#      run the pipeline at the cluster: 
#      
#          make -f qsub.mf CONF=totalrna.conf qc optitype kallisto rrna star bwa
#          make -f qsub.mf CONF=totalrna.conf counts ciri lin_vs_circ
#
#      rename the KD samples according to their standardized bid by running parts of the code chunk below and then run MultiQC:
#
#          multiqc --interactive -o multiqc --ignore '*dcc*' --ignore '*ciri*' --ignore '*hg19*' --ignore '*lin_vs_circ*' -v -f -d -s ./
#
#
########################################################################################################################################################




###########################################################
#
#
#  define the metadata
#  rename samples according to their standardized bid names
#  add FASTQC and featureCount metrics to metadata
#
#{{{
rm(list=ls())
library(RColorBrewer)
library(data.table)


#  load and process the metadata for all libraries
meta<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/all_libraries.csv', header=T, sep=',', select=c(1, 5), col.names=c('id', 'cell_model'))[, treatment:=toupper(sapply(strsplit(cell_model, ' ', fixed=T), '[[', 2))][, cell_model:=sapply(strsplit(cell_model, ' ', fixed=T), '[[', 1)]
meta[treatment %in% 'SI', treatment:='circARID1A SI']
meta[treatment %in% 'SCR', treatment:='circARID1A SCR']


#  load the bids for the KD libraries and add them to the metadata
m<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/KD_libraries.csv', header=F, sep=',', select=1:2, col.names=c('bid', 'id'))
meta<-m[ meta, on='id']


#  add color palette for cell models
N<-meta[, unique(paste(cell_model, treatment, sep='_'))]
cl<-setNames(rev(colorRampPalette(brewer.pal(8, 'Dark2'))(length(N))), N)
for( l in seq_along(cl)){
    meta[ paste(cell_model, treatment, sep='_') %in% names(cl)[l], col:=cl[l] ] 
}
rm(N,cl)


#  replace empty bid by id
meta[ is.na(bid), bid:=id]


#  load FASTQC metrics from MultiQC report
#  identify bid
#  summarize %GC and total number of raw reads across read mates
#  add to metadata
fgc<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/multiqc/multiqc_data/multiqc_fastqc.txt')[, c('Sample', 'Filename', '%GC', 'Total Sequences')]
colnames(fgc)<-c('sample', 'filename', 'gc', 'nreads')
fgc[, bid:=sub('^.*(C[HB][^ ]+).*$', '\\1', sample) ]
fgc<-fgc[, .(gc=mean(gc), nreads=mean(nreads)), by=.(bid)]
meta<-fgc[meta, on='bid']
rm(fgc)


#  process featureCounts statistics and add to metadata
logs<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/ -maxdepth 3 -wholename \'*/counts/genes.tsv.summary\' -print', stdout=T)
names(logs)<-sub('^.*BATCH_201911/(C[HB][^/]+)/counts/.*$', '\\1', logs)
fc<-setNames(vector('list', length(logs)), names(logs))
for(l in logs){
    bid<-names(logs[ logs==l ])
    fc[[bid]]<-setNames(fread(l, header=T, sep='\t', data.table=F)[ c(1, 4:6), 2 ], c('features', 'no_features', 'unassigned_unmapped', 'unassigned_unmapped_mapq'))
}
fc<-as.data.frame(do.call(rbind, fc))
fc$bid<-rownames(fc)
rownames(fc)<-NULL
fc<-fc[, c(5, 1:4)]
fc$unmapped<-rowSums(fc[, 4:5])  #  add unmapped fields together and drop them (unmapped mates + alignments not passing MAPQ threshold asked) 
fc<-data.table(fc[, c(1:3, 6)])
meta<-fc[meta, on='bid']
rm(fc, l, logs, bid)


#  add percentage of the total number of alignments+reads that were dropped 
#  mark failed samples with at least 50% of dropouts in alignments+reads and number of reads covering features below the median
meta[, p_unmapped:=100*unmapped/(features+no_features+unmapped)]
meta[, failed:=ifelse(p_unmapped>=50 & features<median(features, na.rm=T), T, F)]
setcolorder(meta, c('bid', 'id', 'failed', 'gc', 'nreads', 'features', 'no_features', 'unmapped', 'p_unmapped', 'cell_model', 'treatment', 'col'))


#  save
save(meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/metadata.RData')


#  [run once] rename all samples according to their bid entries (if they have one)
m<-meta[ !is.na(bid) ]
for(s in seq_along(m$id)){
    system2('mv', args=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/', c(m[s, id], m[s, bid])))
    system2('sed', args=paste0('-i ', '"s/', m[s,id], '/', m[s, bid], '/g" ', '/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/totalrna.conf')) 
}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/metadata.RData
#
#
###########################################################




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
cnts<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/ -maxdepth 3 -type f -wholename \'*/counts/genes.tsv\' -print', stdout=T)


#  add ids
names(cnts)<-sub('^.*BATCH_201911/(C[BH][^/]+)/counts/.*$', '\\1', cnts)


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
save(totalrna, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/featureCounts.RData')


#  collect the exon counts including the exon junctions
cnts<-sub('genes\\.tsv$', 'exons.tsv', cnts)
exns<-List()
for (n in seq_along(cnts)){
    cat('\nprocessing: ', cnts[n], '\n')
    exns[[ names(cnts)[n] ]]<-parse_exon_counts(cnts[n])
}
rm(n)


#  save
save(exns, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/featureCounts_exons.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/featureCounts.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/featureCounts_exons.RData



#  [featureCounts, circARID1A KD samples] explore count results
#                                         save in Excel sheet TPMs and CPMs of expressed genes 
#{{{
rm(list=ls())
library(openxlsx)
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load reference 
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene']
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]


#  load metadata
#  keep circARID1A KD samples only
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/metadata.RData')
meta<-meta[ grep('circARID1A', bid) ]


#  load gene counts 
#  keep circARID1A KD samples only
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/featureCounts.RData')
gns<-unlist(totalrna[ meta$bid ])
gns$bid<-sub('\\.[0-9]*$', '', rownames(gns))
rownames(gns)<-NULL
gns<-data.table(gns[, c('bid', 'gene_id', 'length', 'counts')])
gns<-gns[, .(gene_id=gene_id, counts=counts, tpm=1e6*counts/length/sum(counts/length), cpm=1e6*counts/sum(counts)), by=.(bid)]
gns[, gene_name:=hsa$gene_name[ match(gene_id, hsa$gene_id) ] ]


#  create the TPM matrix
#  remove unexpressed genes
#  order by top expression
gns.tpm<-dcast(gns, bid ~ gene_name, value.var='tpm', fun.aggregate=sum)  #  TPM matrix
y<-as.data.frame(gns.tpm[, -1])
rownames(y)<-gns.tpm[, bid]
gns.tpm<-t(y)
gns.tpm<-gns.tpm[ rowMeans(gns.tpm)>0, , drop=F]
gns.tpm<-gns.tpm[ order(rowMeans(gns.tpm), decreasing=T), ]


#  create the CPM matrix
#  remove unexpressed genes
#  order by top expression
gns.cpm<-dcast(gns, bid ~ gene_name, value.var='cpm', fun.aggregate=sum)
y<-as.data.frame(gns.cpm[, -1])
rownames(y)<-gns.cpm[, bid]
gns.cpm<-t(y)
gns.cpm<-gns.cpm[ rowMeans(gns.cpm)>0, , drop=F]
gns.cpm<-gns.cpm[ order(rowMeans(gns.cpm), decreasing=T), ]


#  save as Excel sheet
wb<-createWorkbook()
addWorksheet(wb, 'TPMs of expressed genes')
addWorksheet(wb, 'CPMs of expressed genes')
writeDataTable(wb, 'TPMs of expressed genes', rowNames=T, as.data.frame(gns.tpm))
writeDataTable(wb, 'CPMs of expressed genes', rowNames=T, as.data.frame(gns.cpm))
setColWidths(wb, sheet=1, cols=seq_len(ncol(gns.tpm)), widths='auto')
setColWidths(wb, sheet=2, cols=seq_len(ncol(gns.cpm)), widths='auto')
saveWorkbook(wb, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/gene_TPMs_CPMs.xlsx', overwrite=T)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/gene_TPMs_CPMs.xlsx



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
cnts<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/ -maxdepth 3 -type f -wholename \'*/kallisto/abundance.tsv\' -print', stdout=T)


#  add ids
names(cnts)<-sub('^.*BATCH_201911/(C[HB][^/]+)/kallisto/.*$', '\\1', cnts)


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
save(totalrna, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/kallisto.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/kallisto.RData



#  [CIRI2] collect results
#          transgenic circRNAs are split to a separate list
#          circRNAs on alternative contigs, chrM, and chrY are removed
#          define circ_name
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
cir<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/ -maxdepth 3 -type f -wholename \'*/ciri/circRNAs.tsv\' -print', stdout=T)


#  add ids
names(cir)<-sub('^.*BATCH_201911/(C[HB][^/]+)/ciri/.*$', '\\1', cir)


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
#circ<-data.table(as.data.frame(circ))
#circ[, pass:=any(jc_count>=5), by=.(seqnames, start, end, strand)] 
#circ<-circ[ pass %in% TRUE, ][, pass:=NULL]
#circ<-GRanges(seqnames=circ$seqnames, strand=circ$strand, ranges=IRanges(start=circ$start, end=circ$end), data.frame(circ[, -c(1:5)]))
stopifnot( all(lengths(circ$gene_id)==1) )
stopifnot( all(lengths(circ$gene_name)==1) )
circ$gene_id<-unlist(circ$gene_id)
circ$gene_name<-unlist(circ$gene_name)
#trans<-data.table(as.data.frame(trans))
#trans[, pass:=any(jc_count>=5), by=.(seqnames, start, end, strand)] 
#trans<-trans[ pass %in% TRUE, ][, pass:=NULL]
#trans<-GRanges(seqnames=trans$seqnames, strand=trans$strand, ranges=IRanges(start=trans$start, end=trans$end), data.frame(trans[, -c(1:5)]))


#  name them (we do not remove anything since the unified cohort has already been defined)
circ$circ_name<-paste0(circ$gene_id, '|', circ$gene_name,'_', as.character(seqnames(circ)),as.character(strand(circ)), start(circ), '-', end(circ))


#  save
save(circ, trans, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/circRNAs_CIRI2.RData')
x<-data.frame(circ)[, -c(1:3,5)][, c(9, 1:8)]
x<-x[ order(x$jc_count, decreasing=T), ]
write.xlsx(x, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/circRNAs_CIRI2.xlsx', col.names=T, row.names=F, sheetName='circRNAs', append=F)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/circRNAs_CIRI2.{RData,xlsx}



#  [linear vs circular junction quantification] collect and process the quantification results
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)


#  functions
#{{{

parse_results<-function(B='', SSID=''){
    require(data.table)

    #  try to load the counts
    n<-tryCatch({
        cn<-fread(B, header=F, sep='\t', col.names=c('tx_name', 'count'))
        total<-cn[, sum(count)]


        #  separate circular from linear junctions
        ci<-cn[ grep('\\|jn_', tx_name, invert=T) ]
        colnames(ci)<-c('circ_name', 'c.count')
        li<-cn[ grep('\\|jn_', tx_name, invert=F) ]
        colnames(li)<-c('junction_name', 'l.count')
        rm(cn)


        #  add gene_id, gene_name, the linear junctions within the circRNA associated with the circular junction 
        #  and count==NA if the junction is not covered in this sample
        CI<-data.table(data.frame(mcols(CIRCS)[, c('circ_name', 'gene_id', 'gene_name', 'linear_junctions_within')]))
        ci<-ci[ CI, on='circ_name' ]


        #  add gene_id for the linear junction and collect all counts and junction names under the same gene_id
        li<-li[ LINEAR, on='junction_name' ][, .(l.counts=list(l.count), l.junct=list(junction_name)), by=.(gene_id)]


        #  add linear junction counts and junction names to the circular junctions
        ci<-ci[ li, on='gene_id']


        #  go over each circular junction and compute the mean/max counts across all linear junctions and across all linear junctions outside
        #  of the corresponding circRNA range
        ci<-ci[, .(c.count=c.count, gene_id=gene_id, gene_name=gene_name, l.count=mean(unlist(l.counts), na.rm=T), 
              l.count.max=as.numeric(max(unlist(l.counts), na.rm=T)), 
              l.count.out=mean(unlist(l.counts)[ unlist(l.junct) %in% setdiff(unlist(l.junct), unlist(linear_junctions_within)) ], na.rm=T),
              l.count.out.max=as.numeric(max(unlist(l.counts)[ unlist(l.junct) %in% setdiff(unlist(l.junct), unlist(linear_junctions_within)) ], na.rm=T))
              ), by=.(circ_name)]  #  warnings about "returning -Inf" are expected for unexpressed genes


        #  replace -Inf with NA for fully unexpressed genes or for genes with no outter linear junctions found expressed
        ci[ l.count.max %in% -Inf, l.count.max:=NA]
        ci[ l.count.out.max %in% -Inf, l.count.out.max:=NA]


        #  add total number of counts
        ci$total<-total


        #  add the bid for easy processing afterwards
        ci$bid<-SSID


        return(ci)

    }, error=function(e){
        warning(e)
        return(data.table())
    }) 
}

#}}}


#  locate the results
lin_cir<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/ -maxdepth 3 -type f -wholename \'*/lin_vs_circ/counts.tsv\' -print', stdout=T)


#  add sequencing-sample-id (bid)
names(lin_cir)<-sub('^.*/(C[HB][^/]+)/.*$', '\\1', lin_cir)


#  order them
lin_cir<-lin_cir[ order(names(lin_cir)) ]


#  load circular and linear junctions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular.RData')


#  collect results
lin.cir<-List()
for (n in seq_along(lin_cir)){
    cat('\nprocessing: ', lin_cir[n], '\n')
    lin.cir[[ names(lin_cir)[n] ]]<-parse_results(lin_cir[n], names(lin_cir)[n])  #  "returning -Inf" warnings expected for unexpressed junctions
}
rm(n,lin_cir)


#  save
save(lin.cir, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/circRNAs_linear_vs_circular_collected_results.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/circRNAs_linear_vs_circular_collected_results.RData



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
hla<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/ -maxdepth 4 -type f -wholename \'*/optitype/*/*result.tsv\' -print', stdout=T)


#  add ids
names(hla)<-sub('^.*BATCH_201911/(C[HB][^/]+)/optitype/.*$', '\\1', hla)


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
save(hla, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/hla_typing.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/hla_typing.RData



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
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/hla_typing.RData')
h<-do.call(rbind, lapply(hla, function(x){ x[1, ] }))[, 1:6]
h$bid<-rownames(h)
h<-data.table(h)[ order(A1, A2, B1, B2, C1, C2) ]
h[, hla:=paste(A1, A2, B1, B2, C1, C2, sep='_')]


#  load sample metadata
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/metadata.RData')


#  load count data
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/featureCounts.RData')


#  indicate the identical HLA types
h[hla %in% names(table(hla)[ table(hla)>1]), ]


#  scatterplot of a pair of samples
scatterplots('CB-IMR5-circARID1A_si6-R01-R1', 'CB-IMR5-circARID1A_si6-R01-R2', figs.dir='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/')

#}}}



#  [STAR] collect mapping statistics for all samples
#{{{
rm(list=ls())
library(data.table)


#  locate STAR Log.final.out files 
logs<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/ -maxdepth 3 -wholename \'*/star/Log.final.out\' -print', stdout=T)


#  populate mapping summary data.frame
map.s<-data.frame(logs=logs, dir='', bid='', raw_reads=0, unimapped=0, multimapped=0, unmapped=0, alignments=0)
for (l in seq_along(logs)){
    r<-readLines(logs[l]) 
    map.s$dir[l]<-dirname(dirname(logs[l]))
    map.s$bid[l]<-sub('^.*/(C[HB][^/]+)/.*$', '\\1',logs[l])
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
save(map.s, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/mapping_summary.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/mapping_summary.RData



#  [STAR] barplots 
#{{{
rm(list=ls())
library(data.table)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/trim_text.R')


#  load mapping summary and metadata 
#  order libraries according to metadata order
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/mapping_summary.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/metadata.RData')
stopifnot( length(setdiff( map.s$bid , meta$bid ))==0 )  #  all sequenced libraries show up in the metadata
meta<-meta[ bid %in% map.s$bid ]
map.s<-map.s[ match(meta$bid, map.s$bid), ]


#  colors related to alignments and samples
cl.l<-data.frame(symbol=c('unimapped', 'multimapped','unmapped'), 
                  color=c('#228B22',  #  forestgreen
                          '#1874CD',  #  dodgerblue3
                          '#B22222')) #  firebrick
cl.s<-setNames( unique(meta[, treatment]), unique(meta[, col]) )


#  CLICK on it once to make sure it does not redraw, or options(scipen=-20) might be IGNORED
x11(width=20, height=18, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 


#  stacked bars
par(mar=c(10.0, 12.0, 2.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-map.s[, c('unimapped', 'multimapped', 'unmapped')]
rownames(x)<-sub('-11-R01$', '', map.s$bid)
YMAX<-max( rowSums(x) - rowSums(x)%%1e5 + 1e5 )
YTICK<-pretty(c(0, YMAX), 5)
YMAX<-tail(YTICK, 1)
bp<-barplot(t(x), plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(bp)+c(-1, +1), axes=F, ann=F, xaxs='i', yaxs='i') 
bp<-barplot(t(x), border='white', col=cl.l[match(colnames(x), cl.l$symbol), 2], axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, YMAX), add=T)
options(scipen=-20)
axis(2, at=YTICK, line=-1, cex.axis=2.4)
mtext(text=sub('CH_JS_AM_|CB-IMR5-circARID1A_', '', rownames(x)), side=1, line=0, at=bp, col=meta$col, las=2, adj=1, cex=1.8)
mtext('Number of raw reads', side=2, line=9, padj=+0.1, las=0, cex=2.4)
legend(x=par('usr')[2]*0.72, y=1.05*par('usr')[4], legend=cl.l$symbol[match(colnames(x), cl.l$symbol)] , col=cl.l$color[match(colnames(x), cl.l$symbol)], bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.40, x.intersp=0.15, seg.len=0.5)
legend(x=par('usr')[1]*0.5, y=1.05*par('usr')[4], legend=cl.s , col=names(cl.s),  bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.40, x.intersp=0.15, seg.len=0.5)
#mtext(paste0('Number of samples = ', nrow(x)), side=3, line=-1, padj=-0.6, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/barplot_mapping_statistics_alignments.svg', width=20, height=18, bg='white', antialias='subpixel', pointsize=20, family='Arial')
options(scipen=0)


#  clean up
dev.off()

#}}}




###################################
#
#
#  differential expression analysis
#  enrichment analysis 
#
#
###################################




#  [run once] Prepare circRNAs and genes for the analysis. We compute:
#
#                 gene TPMs, raw counts and variance-stabilized log2-transformed counts
#                 circRNA raw counts and variance-stabilized log2-transformed counts based on size factors computed from gene counts
#
#                 PCA for genes and circRNAs is done THROUGHOUT THE SAMPLES using centered but not scaled variance-stabilized 
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


#  load reference 
#  discard chrM, chrY genes
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene']
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]
seqlevels(hsa, pruning.mode='coarse')<-seqlevels(hsa)[ grep('chrM|chrY', seqlevels(hsa), invert=T) ]


#  load metadata
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/metadata.RData')


#  load all of circRNAs
#  keep exon-exon ones and discard chrM, chrY genes
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/circRNAs_CIRI2.RData')
circs<-circ[ circ$region %in% 'exon' & circ$gene_id %in% hsa$gene_id ]
rm(circ)


#  load featureCounts 
#  discard chrM, chrY genes
#  compute TPMs, CPMs based on the filtered list of genes
load('/data/sequencing/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/featureCounts.RData')
gns<-unlist(totalrna[ meta$bid] )
gns$bid<-sub('\\.[0-9]*$', '', rownames(gns))
rownames(gns)<-NULL
gns<-gns[ gns$gene_id %in% hsa$gene_id, ]
gns<-data.table(gns[, c('bid', 'gene_id', 'length', 'counts')])
N<-gns[, .(N=sum(counts)), by=.(bid)]
N<-setNames(N$N, N$bid)
gns<-gns[, .(gene_id=gene_id, counts=counts, tpm=1e6*counts/length/sum(counts/length), cpm=1e6*counts/sum(counts)), by=.(bid)]
gns.cnt<-dcast(gns, bid ~ gene_id, value.var='counts', fun.aggregate=sum)   #  keep counts for all genes
y<-as.data.frame(gns.cnt[, -1])
rownames(y)<-gns.cnt[, bid]
gns.cnt<-t(y[names(N), ])
stopifnot( all.equal( names(N), colnames(gns.cnt) ) )
gns.tpm<-dcast(gns, bid ~ gene_id, value.var='tpm', fun.aggregate=sum)  #  TPM matrix
y<-as.data.frame(gns.tpm[, -1])
rownames(y)<-gns.tpm[, bid]
gns.tpm<-t(y[names(N), ])
stopifnot( all.equal( names(N), colnames(gns.tpm) ) )
circs$nreads<-N[ circs$bid ]
circs$cpm<-circs$jc_count/circs$nreads*1e6
rm(N,totalrna,y,gns)


#  remove unexpressed genes 
gns.cnt<-gns.cnt[rowSums(gns.cnt)!=0,  ]
gns.tpm<-gns.tpm[rowSums(gns.tpm)!=0,  ]


#  compute size factors across all samples (based on the filtered list of genes)
#  compute variance-stabilized gene counts
dds<-DESeqDataSetFromMatrix(countData=ceiling(gns.cnt), colData=data.frame(cell_model=factor(meta$cell_model), row.names=meta$bid), design=~1)
gns.sf<-sizeFactors(estimateSizeFactors(dds, type='poscounts'))
gns.vs<-assay(varianceStabilizingTransformation(dds, fitType='local'))
rm(dds)


#  PCA on variance-stabilized (and glog2-transformed) gene counts centered but not scaled
gns.pca<-prcomp(t(gns.vs), center=T, scale.=F)
gns.ve<-round(1000 * gns.pca$sdev^2/sum(gns.pca$sdev^2))/10 


#  summarize circRNA counts at the gene level 
#  force size factors to be those based from gene counts
#  variance-stabilize
#  PCA on variance-stabilized (and glog2-transformed) counts centered by not scaled
crs.cnt<-data.table(data.frame(mcols(circs)[, c('gene_id', 'jc_count', 'bid')]))[, .(jc_count=sum(jc_count)), by=.(bid, gene_id)]
crs.cnt<-dcast(crs.cnt, bid ~ gene_id, value.var='jc_count', fun.aggregate=sum)
crs.cnt<-t(as.matrix(data.frame(crs.cnt[, -1], row.names=crs.cnt$bid, check.names=F)))[, meta$bid]
dds<-DESeqDataSetFromMatrix(countData=crs.cnt, colData=data.frame(cell_model=factor(meta$cell_model), row.names=meta$bid), design=~1)
sizeFactors(dds)<-gns.sf
crs.vs<-assay(varianceStabilizingTransformation(dds, fitType='local'))
crs.pca<-prcomp(t(crs.vs), center=T, scale.=F)
crs.ve<-round(1000 * crs.pca$sdev^2/sum(crs.pca$sdev^2))/10 
rm(dds)


#  save all including the reference for convenience
save(list=c('hsa', ls(pattern='(gns|crs)\\.'), 'meta', 'circs'), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_circRNAs+genes_all_libraries.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_circRNAs+genes_all_libraries.RData



#  [circARID1A siRNA vs scrambled] clustering 
#                                  DESeq2 analysis
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
    m<-res[ intersect(annot$gene_id, rownames(res)), c('baseMean', 'log2FoldChange', 'gene_name')]
    m$col<-annot$col[ match( rownames(m), annot$gene_id ) ]
    m$cex<-annot$cex[ match( rownames(m), annot$gene_id ) ]
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
        points(r$baseMean , r$log2FoldChange, pch=21, lwd=6, col='black', bg=r$col, cex=r$cex)
        legend('topleft', legend=r$gene_name, bty='n', lty=0, lwd=0, pch=21, col='black', pt.bg=r$col, pt.cex=1.8, pt.lwd=4, cex=1.2, x.intersp=-0.4, y.intersp=0.4)
    }
    if(lfcThreshold>0){
        abline(h=c(-lfcThreshold, lfcThreshold), lty=1, lwd=4, col='cyan4')
    }
    #identify(RES$baseMean, RES$log2FoldChange, labels=RES$gene_name, cex=0.7, offset=0.2, xpd=T)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/plotMA_', Type, '_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    rm(YLIM,XLIM)
    dev.off()


    #  mean-sd plots to see if variance strongly depends on mean
    X11(width=12, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
    par(mar=c(5,4,0.1,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.4, cex.axis=1.4)
    meanSdPlot(VSC)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/meanSD_', Type, '_', sub(' vs ', '_', CND), '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
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
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/heatmap_', Type, '_vsc_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
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
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/heatmap_', Type, '_cor_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
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
    print(
        ggplot(x, aes(PC1, PC2, color=get(condition))) + 
        geom_point(size=10) + 
        geom_text_repel(aes(label=rownames(PCA$x)), size=10, box.padding=0.5, segment.size=1.0, min.segment.length=1.0) +
        scale_shape_manual(values=18) + 
        scale_fill_manual(name='sample', values=unique(x$col)) + 
        scale_color_manual(name='sample', values=unique(x$col)) + 
        theme(text = element_text(family='Arial'), axis.text.x=element_text(size=28), axis.title.x=element_text(size=28), 
            axis.title.y=element_text(size=28), axis.text.y=element_text(size=28), 
            legend.text=element_text(size=24), legend.title=element_text(size=24, face='plain'), aspect.ratio=1, panel.background = element_blank(),
            axis.line.x=element_line(size=0.2, colour = 'black', linetype='solid'), 
            axis.line.y=element_line(size=0.2, colour = 'black', linetype='solid'),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
        scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
        xlab(paste0('PC1: ', VE[1], '% variance')) + ylab(paste0('PC2: ', VE[2], '% variance')) 
    )
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/plotPCA_', Type, '_vsc_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    rm(XLIM,YLIM)
    dev.off()


    #  save DESeq results along with gene counts
    options(scipen=0)
    save(DDS,CND,RES,VSC,PCA,VE, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.RData'), compress=T)


    #  convert to table 
    x<-RES
    x$gene_id<-rownames(x)
    x<-x[, c('gene_id', 'gene_name', 'baseMean', 'log2FoldChange', 'pvalue', 'padj')]
    rownames(x)<-NULL
    x<-as.data.frame(x)
    write.table(x, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.tsv'), quote=F, sep='\t', row.names=F, col.names=T) 
    write.xlsx(x, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.xlsx'), row.names=F, col.names=T, sheetName=gsub(' ', '_', CND)) 
    rm(x)


   #readline("Hit ENTER to close all figures: ") 
   #for (n in setdiff(dev.list(), DEV_START)){ dev.off(n) }

   return(list(dds=DDS, res=RES, cnd=CND))
}

#}}}


#  annotate ARID1A at the MAplots
annot<-data.frame(gene_id='ENSG00000117713.20', gene_name='ARID1A', col='darkorange', cex=2.0)


#  load pre-prepared counts etc. for all samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_circRNAs+genes_all_libraries.RData')


#  [circRNAs + genes] clustering and PCA
#{{{

#  [circRNAs] heatmap based on Euclidean distance of variance stabilized (and glog2-transformed) counts 
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-crs.vs[, m$bid]
colnames(x)<-sub('R01-', '', colnames(x))
colnames(x)<-sub('CH_JS_AM_', '', colnames(x))
colnames(x)<-sub('CB-IMR5-', '', colnames(x))
d<-dist(t(x), method='euclidean')
hc<-hclust(d, method='ward.D2')
ex<-data.frame(Type=factor(m$treatment, exclude=F), row.names=colnames(x))
cl<-setNames(list(setNames( unique(m$col), unique(m$treatment) )), colnames(ex))
ph<-pheatmap(as.matrix(d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=14, fontsize_row=16, fontsize_col=16, fontsize_number=10.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/heatmap_circRNAs_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [genes] heatmap based on Euclidean distance of variance stabilized (and glog2-transformed) counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-gns.vs[, m$bid]
colnames(x)<-sub('R01-', '', colnames(x))
colnames(x)<-sub('CH_JS_AM_', '', colnames(x))
colnames(x)<-sub('CB-IMR5-', '', colnames(x))
d<-dist(t(x), method='euclidean')
hc<-hclust(d, method='ward.D2')
ex<-data.frame(Type=factor(m$treatment, exclude=F), row.names=colnames(x))
cl<-setNames(list(setNames( unique(m$col), unique(m$treatment) )), colnames(ex))
ph<-pheatmap(as.matrix(d), color=colorRampPalette(brewer.pal(9,'GnBu'))(20), border_color=NA, scale='none', 
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=14, fontsize_row=16, fontsize_col=16, fontsize_number=10.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/heatmap_genes_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [circRNAs] heatmaps based on Spearman correlations of raw counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-crs.cnt[, m$bid]
colnames(x)<-sub('R01-', '', colnames(x))
colnames(x)<-sub('CH_JS_AM_', '', colnames(x))
colnames(x)<-sub('CB-IMR5-', '', colnames(x))
d<-cor(x, method='spearman', use='pairwise.complete.obs')
hc<-hclust(as.dist(1-d), method='ward.D2')
ex<-data.frame(Type=factor(m$treatment, exclude=F), row.names=colnames(x))
cl<-setNames(list(setNames( unique(m$col), unique(m$treatment) )), colnames(ex))
ph<-pheatmap(as.matrix(d), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(20)), border_color=NA, scale='none', 
        breaks=seq(0, 1, length.out=21),
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=14, fontsize_row=16, fontsize_col=16, fontsize_number=10.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/heatmap_circRNAs_cor_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [genes] heatmaps based on Spearman correlations of raw counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-gns.cnt[, m$bid]
colnames(x)<-sub('R01-', '', colnames(x))
colnames(x)<-sub('CH_JS_AM_', '', colnames(x))
colnames(x)<-sub('CB-IMR5-', '', colnames(x))
d<-cor(x, method='spearman', use='pairwise.complete.obs')
hc<-hclust(as.dist(1-d), method='ward.D2')
ex<-data.frame(Type=factor(m$treatment, exclude=F), row.names=colnames(x))
cl<-setNames(list(setNames( unique(m$col), unique(m$treatment) )), colnames(ex))
ph<-pheatmap(as.matrix(d), color=rev(colorRampPalette(brewer.pal(9,'RdBu'))(20)), border_color=NA, scale='none', 
        breaks=seq(0.7, 1, length.out=21),  #  reduce range since samples are highly correlated
        cluster_rows=hc,
        cluster_cols=hc,
        #cutree_row=4, cutree_col=4,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=ex, annotation_row=ex, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=T, 
        #display_numbers=T, number_format='%.1f', number_color='grey39',
        fontsize=14, fontsize_row=16, fontsize_col=16, fontsize_number=10.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/heatmap_genes_cor_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x,m,d,ex,hc,ph,cl)
dev.off()


#  [circRNAs] PCA based on variance stabilized (and log2-transformed) counts centered but not scaled
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-crs.vs[, m$bid]
colnames(x)<-sub('R01-', '', colnames(x))
colnames(x)<-sub('CH_JS_AM_', '', colnames(x))
colnames(x)<-sub('CB-IMR5-', '', colnames(x))
p<-prcomp(t(x), center=T, scale.=F)
v<-round(1000 * p$sdev^2/sum(p$sdev^2))/10 
n<-cbind( data.frame( p$x[, c('PC1', 'PC2')] ), Type=m$treatment, bid=colnames(x))
cl<-list('Type'=setNames( unique(m$col), unique(m$treatment) ))
XLIM<-pretty(range(n[, 'PC1']), 2)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(n[, 'PC2']), 2)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(n[, c('PC1', 'PC2', 'Type', 'bid')], aes(PC1, PC2, color=Type)) + 
    geom_point(size=10) + 
    scale_shape_manual(values=18) + 
    scale_fill_manual(name='Type', values=cl$Type) + 
    scale_color_manual(name='Type', values=cl$Type) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=18), axis.title.x=element_text(size=22), 
        axis.title.y=element_text(size=22), axis.text.y=element_text(size=18), 
        legend.text=element_text(size=18), legend.title=element_text(size=18, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.2, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.2, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', v[1], '% variance')) + ylab(paste0('PC2: ', v[2], '% variance'))
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/plotPCA_circRNAs_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  [genes] PCA based on variance stabilized glog2-transformed counts centered but not scaled
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-gns.vs[, m$bid]
colnames(x)<-sub('R01-', '', colnames(x))
colnames(x)<-sub('CH_JS_AM_', '', colnames(x))
colnames(x)<-sub('CB-IMR5-', '', colnames(x))
p<-prcomp(t(x), center=T, scale.=F)
v<-round(1000 * p$sdev^2/sum(p$sdev^2))/10 
n<-cbind( data.frame( p$x[, c('PC1', 'PC2')] ), Type=m$treatment, bid=colnames(x))
cl<-list('Type'=setNames( unique(m$col), unique(m$treatment) ))
XLIM<-pretty(range(n[, 'PC1']), 2)
XLIM<-c(XLIM[1], XLIM[length(XLIM)])
YLIM<-pretty(range(n[, 'PC2']), 2)
YLIM<-c(YLIM[1], YLIM[length(YLIM)])
ggplot(n[, c('PC1', 'PC2', 'Type', 'bid')], aes(PC1, PC2, color=Type)) + 
    geom_point(size=10) + 
    scale_shape_manual(values=18) + 
    scale_fill_manual(name='Type', values=cl$Type) + 
    scale_color_manual(name='Type', values=cl$Type) + 
    theme(text=element_text(family='Arial'), axis.text.x=element_text(size=18), axis.title.x=element_text(size=22), 
        axis.title.y=element_text(size=22), axis.text.y=element_text(size=18), 
        legend.text=element_text(size=18), legend.title=element_text(size=18, face='plain'), 
        aspect.ratio=1, panel.background = element_blank(),
        axis.line.x=element_line(size=0.2, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.2, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(lim=XLIM, breaks=pretty(XLIM,5)) + 
    scale_y_continuous(lim=YLIM, breaks=pretty(YLIM,5)) + 
    xlab(paste0('PC1: ', v[1], '% variance')) + ylab(paste0('PC2: ', v[2], '% variance'))
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/plotPCA_genes_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,p,v,n,cl,XLIM,YLIM)
dev.off()

#}}}


#  [circRNAs] DESeq2 forcing library sizes from the gene counts 
m<-meta[grepl('CB-IMR5-circARID1A_', bid) ]
M<-crs.cnt[, m$bid ]
SF<-gns.sf[ colnames(M) ]
colnames(M)<-sub('R01-', '', colnames(M))
colnames(M)<-sub('CH_JS_AM_', '', colnames(M))
colnames(M)<-sub('CB-IMR5-', '', colnames(M))
names(SF)<-sub('R01-', '', names(SF))
names(SF)<-sub('CH_JS_AM_', '', names(SF))
names(SF)<-sub('CB-IMR5-', '', names(SF))
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=colnames(M))
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=ceiling(M), colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, 'circARID1A SCR')  #  define first level so that the comparison log2(level 2/level 1) is done
sizeFactors(DDS)<-SF[ colnames(DDS) ]
d<-do_deseq2(DDS=DDS, CT=CT, Type='circRNAs', condition='Treatment', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=T, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)


#  [genes] DESeq2 
m<-meta[grepl('CB-IMR5-circARID1A', bid) ]
M<-gns.cnt[, m$bid ]
colnames(M)<-sub('R01-', '', colnames(M))
colnames(M)<-sub('CH_JS_AM_', '', colnames(M))
colnames(M)<-sub('CB-IMR5-', '', colnames(M))
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=colnames(M))
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=ceiling(M), colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, 'circARID1A SCR')  #  define first level so that the comparison log2(level 2/level 1) is done
d<-do_deseq2(DDS=DDS, CT=CT, Type='genes', condition='Treatment', lfcThreshold=log2(1.0), lfcShrinkType='normal', include.annot=T, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_circRNAs_circARID1A SI_circARID1A SCR.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_genes_circARID1A SI_circARID1A SCR.RData



#  [circARID1A siRNA vs scrambled] process of the DESeq2 results into HTML-table-ready objects
#                                  do MSigDB C2, C3 enrichment analysis
#                                  do MSigDB GSEA analysis
#                                  do GO MF, BP analysis
#                                  do mesenchymal/adrenergic marker enrichment test
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(clusterProfiler)
library(topGO)
library(org.Hs.eg.db)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load circRNAs and genes used 
#  simplify the sample names like it is done in the DESeq2 objects
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_circRNAs+genes_all_libraries.RData')
circs<-circs[ grep('CB-IMR5-circARID1A', circs$bid) ]
circs$bid<-sub('CB-IMR5-', '', circs$bid)
circs$bid<-sub('R01-', '', circs$bid)


#  simplify circARID1A siRNA vs scrambled bids
meta$bid<-sub('CB-IMR5-', '', meta$bid)
meta$bid<-sub('R01-', '', meta$bid)


#  [circRNAs] load DE results
#             separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_circRNAs_circARID1A SI_circARID1A SCR.RData')
all.circ<-RES
up.circ<-subset(RES, log2FoldChange>0 & padj<0.1)    #  FDR cutoff: 0.1
down.circ<-subset(RES, log2FoldChange<0 & padj<0.1)  #  FDR cutoff: 0.1
bid.circ<-colnames(DDS)
meta.circ<-meta[ match( bid.circ, bid ), ]
rm(CND,DDS,PCA,RES,VE,VSC)


#  [genes] load DE results
#          separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_genes_circARID1A SI_circARID1A SCR.RData')
all.gns<-RES
up.gns<-subset(RES, log2FoldChange>0 & padj<0.05)
down.gns<-subset(RES, log2FoldChange<0 & padj<0.05)
bid.gns<-colnames(DDS)
meta.gns<-meta[ match( bid.gns, bid ), ]
rm(CND,DDS,PCA,RES,VE,VSC,meta)


#  C2,C3 MSigDB enrichment analysis (baseMean>=10)
#  C2,C3 MSigDB GSEA (baseMean>=10 and sum(log2FC) for genes with identical name)
#  GO MF/BP enrichment analysis (baseMean>=10)
#  mesenchymal/adrenergic gene marker enrichment test
#{{{

#  load the C2, C3 MSigDB gene sets
c2<-read.gmt('/fast/groups/ag_schulte/work/reference/annotation/MSigDB/c2.all.v7.0.symbols.gmt')
c3<-read.gmt('/fast/groups/ag_schulte/work/reference/annotation/MSigDB/c3.all.v7.0.symbols.gmt')


#  [genes] C2, C3 MSigDB enrichment analysis
#          universe is all genes with baseMean>=10 (if significantly DE genes have baseMean<10 let them drop out)
up.gns.c2<-enricher(up.gns$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=5, maxGSSize=200, universe=subset(all.gns, baseMean>=10)$gene_name, TERM2GENE=c2)  #  to access the main result: up.gns.c2@result
down.gns.c2<-enricher(down.gns$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=5, maxGSSize=200, universe=subset(all.gns, baseMean>=10)$gene_name, TERM2GENE=c2)
up.gns.c3<-enricher(up.gns$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=5, maxGSSize=200, universe=subset(all.gns, baseMean>=10)$gene_name, TERM2GENE=c3)  #  to access the main result: up.gns.c3@result
down.gns.c3<-enricher(down.gns$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=5, maxGSSize=200, universe=subset(all.gns, baseMean>=10)$gene_name, TERM2GENE=c3)


#  [circRNAs] C2, C3 MSigDB enrichment analysis
up.circ.c2<-enricher(up.circ$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=subset(all.circ, baseMean>=10)$gene_name, TERM2GENE=c2)  #  to access the main result: up.circ.c2@result
down.circ.c2<-enricher(down.circ$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=subset(all.circ, baseMean>=10)$gene_name, TERM2GENE=c2)
up.circ.c3<-enricher(up.circ$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=subset(all.circ, baseMean>=10)$gene_name, TERM2GENE=c3)  #  to access the main result: up.circ.c3@result
down.circ.c3<-enricher(down.circ$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=200, universe=subset(all.circ, baseMean>=10)$gene_name, TERM2GENE=c3)


#  [genes, baseMean>=10, add log2FC of genes with identical names] C2 MSigDB GSEA analysis
g<-data.table(data.frame(subset(all.gns, baseMean>=10)))[, .(log2FoldChange=sum(log2FoldChange)), by=.(gene_name)]
g<-sort( setNames( g$log2FoldChange , g$gene_name), decreasing=T ) 
gsea.gns.c2<-GSEA(g, pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=5, maxGSSize=200, TERM2GENE=c2)
gsea.gns.c3<-GSEA(g, pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=5, maxGSSize=200, TERM2GENE=c3)
rm(g)


#  [circRNAs, baseMean>=1]
g<-subset(all.circ, baseMean>=1)
g<-sort( setNames( g$log2FoldChange , g$gene_name), decreasing=T ) 
gsea.circ.c2<-GSEA(g, pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=5, maxGSSize=200, TERM2GENE=c2)
gsea.circ.c3<-GSEA(g, pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=5, maxGSSize=200, TERM2GENE=c3)
rm(g)


#  GO enrichment analysis (baseMean>=10)
#{{{

#  universe for BP and MF should be separate and keeping only genes with annotations so we can compute GeneRatios later on
uni<-sub('\\.[0-9]*$', '', unique(rownames(subset(all.gns, baseMean>=10))))  #  all genes to draw universes from
a<-factor(setNames( rep(0, length(uni)), uni), levels=c(0, 1))
uni.mf<-uni[ uni %in% unique(unlist(genesInTerm(new('topGOdata', ontology='MF', description='all genes', allGenes=a, annot=annFUN.org, mapping='org.Hs.eg.db', ID='Ensembl', nodeSize=1)))) ]
uni.bp<-uni[ uni %in% unique(unlist(genesInTerm(new('topGOdata', ontology='BP', description='all genes', allGenes=a, annot=annFUN.org, mapping='org.Hs.eg.db', ID='Ensembl', nodeSize=1)))) ]


#  upregulated
grp<-sub('\\.[0-9]*$', '', rownames(up.gns))
allG<-setNames(rep(0, length(uni.mf)), uni.mf)
allG[ names(allG) %in% grp  ]<-1
allG<-as.factor(allG)
up.mf.topgo<-new('topGOdata',description='',ontology='MF',allGenes=allG,annot=annFUN.org,mapping='org.Hs.eg.db',ID='Ensembl',nodeSize=1)
up.mf.test<-runTest(up.mf.topgo, algorithm='weight01', statistic='fisher')
up.mf<-data.table(GenTable(up.mf.topgo, p.value=up.mf.test, orderBy='p.value', topNodes=geneData(up.mf.test)['SigTerms'], numChar=120))
up.mf$p.value<-as.numeric(up.mf$p.value)
up.mf[, padj:=p.adjust(p.value, method='fdr')]
l<-lapply(genesInTerm(up.mf.topgo)[ up.mf[p.value<0.05, GO.ID] ], function(g){ intersect(g, grp) }) 
l<-data.table(GO.ID=names(l), gids=l)[, gene_names:=lapply(gids, function(i){ up.gns[ grep(paste(i, sep='', collapse='|'), rownames(up.gns)), 'gene_name' ]})][, gids:=NULL]
up.mf<-l[ up.mf, on='GO.ID']
allG<-setNames(rep(0, length(uni.bp)), uni.bp)
allG[ names(allG) %in% grp  ]<-1
allG<-as.factor(allG)
up.bp.topgo<-new('topGOdata',description='',ontology='BP',allGenes=allG,annot=annFUN.org,mapping='org.Hs.eg.db',ID='Ensembl',nodeSize=1)
up.bp.test<-runTest(up.bp.topgo, algorithm='weight01', statistic='fisher')
up.bp<-data.table(GenTable(up.bp.topgo, p.value=up.bp.test, orderBy='p.value', topNodes=geneData(up.bp.test)['SigTerms'], numChar=120))
up.bp$p.value<-as.numeric(up.bp$p.value)
up.bp[, padj:=p.adjust(p.value, method='fdr')]
l<-lapply(genesInTerm(up.bp.topgo)[ up.bp[p.value<0.05, GO.ID] ], function(g){ intersect(g, grp) }) 
l<-data.table(GO.ID=names(l), gids=l)[, gene_names:=lapply(gids, function(i){ up.gns[ grep(paste(i, sep='', collapse='|'), rownames(up.gns)), 'gene_name' ]})][, gids:=NULL]
up.bp<-l[ up.bp, on='GO.ID']
rm(l)
 

#  downregulated
grp<-sub('\\.[0-9]*$', '', rownames(down.gns))
allG<-setNames(rep(0, length(uni.mf)), uni.mf)
allG[ names(allG) %in% grp  ]<-1
allG<-as.factor(allG)
down.mf.topgo<-new('topGOdata',description='',ontology='MF',allGenes=allG,annot=annFUN.org,mapping='org.Hs.eg.db',ID='Ensembl',nodeSize=1)
down.mf.test<-runTest(down.mf.topgo, algorithm='weight01', statistic='fisher')
down.mf<-data.table(GenTable(down.mf.topgo, p.value=down.mf.test, orderBy='p.value', topNodes=geneData(down.mf.test)['SigTerms'], numChar=120))
down.mf$p.value<-as.numeric(down.mf$p.value)
down.mf[, padj:=p.adjust(p.value, method='fdr')]
l<-lapply(genesInTerm(down.mf.topgo)[ down.mf[p.value<0.05, GO.ID] ], function(g){ intersect(g, grp) }) 
l<-data.table(GO.ID=names(l), gids=l)[, gene_names:=lapply(gids, function(i){ down.gns[ grep(paste(i, sep='', collapse='|'), rownames(down.gns)), 'gene_name' ]})][, gids:=NULL]
down.mf<-l[ down.mf, on='GO.ID']
allG<-setNames(rep(0, length(uni.bp)), uni.bp)
allG[ names(allG) %in% grp  ]<-1
allG<-as.factor(allG)
down.bp.topgo<-new('topGOdata',description='',ontology='BP',allGenes=allG,annot=annFUN.org,mapping='org.Hs.eg.db',ID='Ensembl',nodeSize=1)
down.bp.test<-runTest(down.bp.topgo, algorithm='weight01', statistic='fisher')
down.bp<-data.table(GenTable(down.bp.topgo, p.value=down.bp.test, orderBy='p.value', topNodes=geneData(down.bp.test)['SigTerms'], numChar=120))
down.bp$p.value<-as.numeric(down.bp$p.value)
down.bp[, padj:=p.adjust(p.value, method='fdr')]
l<-lapply(genesInTerm(down.bp.topgo)[ down.bp[p.value<0.05, GO.ID] ], function(g){ intersect(g, grp) }) 
l<-data.table(GO.ID=names(l), gids=l)[, gene_names:=lapply(gids, function(i){ down.gns[ grep(paste(i, sep='', collapse='|'), rownames(down.gns)), 'gene_name' ]})][, gids:=NULL]
down.bp<-l[ down.bp, on='GO.ID']
rm(grp, uni, allG, l)

#}}}


#  mesenchymal/adrenergic gene marker enrichment tests
#{{{

#  load GENCODE v30 markers
load('/fast/groups/ag_schulte/work/reference/annotation/MES_ADR_markers/Suppl Versteeg 2017 - MES ADR Genes_gencode_v30.RData')


#  save the gene groups
up.mes.adr<-down.mes.adr<-list()


#  [upregulated, mesenchymal] Fisher's exact test:
#
#                   |      up          |        not up        |
#  -----------------|------------------|----------------------|
#     mesenchymal   |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-mesenchymal  |      x2          |          y2          | 
grp<-rownames( up.gns )                      
up.mes.adr[['mes']]<-list(genes=intersect(grp, mes.adr[ cell_type %in% 'MES', gene_id ]))
rst<-rownames( subset(all.gns, baseMean>=10 & rownames(all.gns) %in% setdiff( rownames(all.gns), grp )) )  #  only baseMean>=10
p<-fisher.test(data.frame('in'=c(length(intersect(grp, mes.adr[ cell_type %in% 'MES', gene_id ])), 
                              length(setdiff(grp, mes.adr[ cell_type %in% 'MES', gene_id ]))), 
                       'out'=c(length(intersect(rst, mes.adr[ cell_type %in% 'MES', gene_id ])), 
                               length(setdiff(rst, mes.adr[ cell_type %in% 'MES', gene_id ])))), alternative='greater')$p.value
up.mes.adr[['mes']]$pvalue<-unname(p)
rm(p)


#  [upregulated, adrenergic] Fisher's exact test:
#
#                   |      up          |        not up        |
#  -----------------|------------------|----------------------|
#     adrenergic    |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-adrenergic   |      x2          |          y2          | 
grp<-rownames( up.gns )                      
up.mes.adr[['adr']]<-list(genes=intersect(grp, mes.adr[ cell_type %in% 'ADRN', gene_id ]))
rst<-rownames( subset(all.gns, baseMean>=10 & rownames(all.gns) %in% setdiff( rownames(all.gns), grp )) )  #  only baseMean>=10 
p<-fisher.test(data.frame('in'=c(length(intersect(grp, mes.adr[ cell_type %in% 'ADRN', gene_id ])), 
                              length(setdiff(grp, mes.adr[ cell_type %in% 'ADRN', gene_id ]))), 
                       'out'=c(length(intersect(rst, mes.adr[ cell_type %in% 'ADRN', gene_id ])), 
                               length(setdiff(rst, mes.adr[ cell_type %in% 'ADRN', gene_id ])))), alternative='greater')$p.value
up.mes.adr[['adr']]$pvalue<-unname(p)
rm(p)


#  [downregulated, mesenchymal] Fisher's exact test:
#
#                   |     down         |        not down      |
#  -----------------|------------------|----------------------|
#     mesenchymal   |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-mesenchymal  |      x2          |          y2          | 
grp<-rownames( down.gns )                      
down.mes.adr[['mes']]<-list(genes=intersect(grp, mes.adr[ cell_type %in% 'MES', gene_id ]))
rst<-rownames( subset(all.gns, baseMean>=10 & rownames(all.gns) %in% setdiff( rownames(all.gns), grp )) )  #  only baseMean>=10
p<-fisher.test(data.frame('in'=c(length(intersect(grp, mes.adr[ cell_type %in% 'MES', gene_id ])), 
                              length(setdiff(grp, mes.adr[ cell_type %in% 'MES', gene_id ]))), 
                       'out'=c(length(intersect(rst, mes.adr[ cell_type %in% 'MES', gene_id ])), 
                               length(setdiff(rst, mes.adr[ cell_type %in% 'MES', gene_id ])))), alternative='greater')$p.value
down.mes.adr[['mes']]$pvalue<-unname(p)
rm(p)


#  [downregulated, adrenergic] Fisher's exact test:
#
#                   |     down         |        not down      |
#  -----------------|------------------|----------------------|
#     adrenergic    |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#  non-adrenergic   |      x2          |          y2          | 
grp<-rownames( down.gns )                      
down.mes.adr[['adr']]<-list(genes=intersect(grp, mes.adr[ cell_type %in% 'ADRN', gene_id ]))
rst<-rownames( subset(all.gns, baseMean>=10 & rownames(all.gns) %in% setdiff( rownames(all.gns), grp )) )  #  only baseMean>=10 
p<-fisher.test(data.frame('in'=c(length(intersect(grp, mes.adr[ cell_type %in% 'ADRN', gene_id ])), 
                              length(setdiff(grp, mes.adr[ cell_type %in% 'ADRN', gene_id ]))), 
                       'out'=c(length(intersect(rst, mes.adr[ cell_type %in% 'ADRN', gene_id ])), 
                               length(setdiff(rst, mes.adr[ cell_type %in% 'ADRN', gene_id ])))), alternative='greater')$p.value
down.mes.adr[['adr']]$pvalue<-unname(p)
rm(p)

#}}}

#}}}


#  [circRNAs] rounded mean circular junction coverage across samples for each circRNA isoform
#             build genomic locations as well (CIRI2 uses 1-based coordinate system)
circ<-data.table(data.frame(circs[ circs$bid %in% bid.circ ]))[, .(jc_count=ceiling(mean(jc_count)), gene_id=unique(gene_id), locus=paste0('[', sub('chr', '', seqnames), ':', start, '-', end, '](http://www.ensembl.org/Homo_sapiens/Location/View?r=', sub('chr', '', seqnames), ':', start, '-', end, ')')), by=.(gene_name, seqnames, start, end, strand)]


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
save(list=ls(pattern='up\\.|down\\.|all\\.|bid\\.|circ|crs|gns|gsea|uni'), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_kable-ready_genes+circRNAs_circARID1A SI_circARID1A SCR.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_kable-ready_genes+circRNAs_circARID1A SI_circARID1A SCR.RData



#  [circARID1A siRNA vs scrambled] plots 
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(clusterProfiler)
library(DESeq2)
library(topGO)
library(org.Hs.eg.db)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load post-processed data
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_kable-ready_genes+circRNAs_circARID1A SI_circARID1A SCR.RData')


#  recycle
x11(width=20, height=14, bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  [GO BP] upregulated genes (p-value cutoff<5e-3)
#{{{

x<-up.bp[p.value<5e-3, c('Term', 'Significant', 'p.value')][, GeneRatio:=Significant/length(sigGenes(up.bp.topgo))][ order(GeneRatio) ][, Term:=factor(Term, levels=Term)]
colnames(x)<-sub('Significant', 'Count', colnames(x))
print(ggplot(x, aes(GeneRatio, Term, size=Count, color=`p.value`)) + 
    geom_point() +
    scale_color_continuous(low='red', high='blue', name='p.value', breaks=pretty(x$p.value, n=5), 
        guide=guide_colorbar(reverse=T, draw.ulim=T, draw.llim=T, barheight=8)) +
    ylab(NULL) + 
    scale_size(range=c(3, 10)) +
    theme(text=element_text(family='Arial'), 
        axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/dotplot_genes_GOBP_upDE_circARID1A SI_circARID1A SCR.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [GO BP] downregulated genes (p-value cutoff<5e-3)
#{{{

x<-down.bp[p.value<5e-3, c('Term', 'Significant', 'p.value')][, GeneRatio:=Significant/length(sigGenes(down.bp.topgo))][ order(GeneRatio) ][, Term:=factor(Term, levels=Term)]
colnames(x)<-sub('Significant', 'Count', colnames(x))
print(ggplot(x, aes(GeneRatio, Term, size=Count, color=`p.value`)) + 
    geom_point() +
    scale_color_continuous(low='red', high='blue', name='p.value', breaks=pretty(x$p.value, n=5), 
        guide=guide_colorbar(reverse=T, draw.ulim=T, draw.llim=T, barheight=8)) +
    ylab(NULL) + 
    scale_size(range=c(3, 10)) +
    theme(text=element_text(family='Arial'), 
        axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/dotplot_genes_GOBP_downDE_circARID1A SI_circARID1A SCR.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [GO BP] upregulated genes (neuron/differentiation-specific p-value cutoff<0.05)
#{{{

x<-up.bp[p.value<0.05 & grepl('neuron|differentiation', Term), c('Term', 'Significant', 'p.value')][, GeneRatio:=Significant/length(sigGenes(up.bp.topgo))][ order(GeneRatio) ][, Term:=factor(Term, levels=Term)]
colnames(x)<-sub('Significant', 'Count', colnames(x))
print(ggplot(x, aes(GeneRatio, Term, size=Count, color=`p.value`)) + 
    geom_point() +
    scale_color_continuous(low='red', high='blue', name='p.value', breaks=pretty(x$p.value, n=5), 
        guide=guide_colorbar(reverse=T, draw.ulim=T, draw.llim=T, barheight=8)) +
    ylab(NULL) + 
    scale_size(range=c(3, 10)) +
    theme(text=element_text(family='Arial'), 
        axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/dotplot_genes_GOBP_neuron+differentiation_upDE_circARID1A SI_circARID1A SCR.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [GO BP] downregulated genes (neuron/differentiation-specific p-value cutoff<0.05)
#{{{

x<-down.bp[p.value<0.05 & grepl('neuron|differentiation', Term), c('Term', 'Significant', 'p.value')][, GeneRatio:=Significant/length(sigGenes(down.bp.topgo))][ order(GeneRatio) ][, Term:=factor(Term, levels=Term)]
colnames(x)<-sub('Significant', 'Count', colnames(x))
print(ggplot(x, aes(GeneRatio, Term, size=Count, color=`p.value`)) + 
    geom_point() +
    scale_color_continuous(low='red', high='blue', name='p.value', breaks=pretty(x$p.value, n=5), 
        guide=guide_colorbar(reverse=T, draw.ulim=T, draw.llim=T, barheight=8)) +
    ylab(NULL) + 
    scale_size(range=c(3, 10)) +
    theme(text=element_text(family='Arial'), 
        axis.text.x=element_text(size=28, color='black'), axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=28, color='black'), axis.title.y=element_text(size=30), 
        legend.text=element_text(size=25), legend.title=element_text(size=25, face='plain'), 
        panel.background = element_blank(),
        axis.line.x=element_line(size=0.5, colour = 'black', linetype='solid'), 
        axis.line.y=element_line(size=0.5, colour = 'black', linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/dotplot_genes_GOBP_neuron+differentiation_downDE_circARID1A SI_circARID1A SCR.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [GSEA C2]
#{{{

#  split tags, list and signal to separate columns and order by signal
g<-as.data.table(gsea.gns.c2@result)[, leading_edge:=lapply(strsplit(leading_edge, ', '), function(x){ as.numeric(sub('^.*=([0-9]+)%$', '\\1', x)) })]
g[, c('tags', 'list', 'signal'):=list(sapply(leading_edge, '[[', 1), sapply(leading_edge, '[[', 2), sapply(leading_edge, '[[', 3))][, c('leading_edge', 'Description'):=NULL]
g<-g[order(-signal, -setSize)]


#  histogram of signals
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.0, 7.0, 0.5, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
h<-hist(g[, signal], breaks=seq(0, 100, 5), col='grey39', border='white', xlab='', ylab='', main='', xlim=c(0, 100), ylim=c(0, 100), add=F)
mtext('Signal', side=1, line=3, padj=+0.3, las=0, cex=2.4)
mtext('Frequency', side=2, line=4, padj=-0.5, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/histogram_genes_C2_GSEA_signal.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  table to enrichments and depletions
g[, table(NES>0) ]
# 
# FALSE  TRUE 
#   238    54 


#  [enrichments] look at specific cases, save them if you want.
x<-g[NES>0, ]
n<-13; gseaplot(gsea.gns.c2, geneSetID=x[n, ID], title=x[n, ID])
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/gseaplot_genes_C2_', x[n, ID], '_circARID1A SI_circARID1A SCR.svg'), width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [depletions] look at specific cases, save them if you want.
x<-g[NES<0 & setSize>10, ]
n<-7; gseaplot(gsea.gns.c2, geneSetID=x[n, ID], title=x[n, ID])
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/gseaplot_genes_C2_', x[n, ID], '_circARID1A SI_circARID1A SCR.svg'), width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}



#  [differentiation markers] are they differentially expressed in KD?
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(clusterProfiler)
library(DESeq2)
library(topGO)
library(org.Hs.eg.db)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load post processed data
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_kable-ready_genes+circRNAs_circARID1A SI_circARID1A SCR.RData')


#  define differentiation markers
diff.markers<-list(positive=c('MART', 'DCX', 'CDH5', 'NTS', 'NEFL', 'GAP43'), 
                   negative=c('ASCL1', 'SOX2', 'SOX4', 'NTRK1'))


#  positive markers
all.gns[ all.gns$gene_name %in% diff.markers$positive, ]
#                      baseMean log2FoldChange     lfcSE      stat         pvalue         padj   gene_name
#                     <numeric>      <numeric> <numeric> <numeric>      <numeric>    <numeric> <character>
# ENSG00000077279.18 8610.96275     -0.2298611 0.0442677 -5.192456 0.000000207538 0.0000060186         DCX
# ENSG00000277586.2  6324.81304      0.2101744 0.0504469  4.166232 0.000030967574 0.0004454880        NEFL
# ENSG00000172020.12  136.30020     -0.1249035 0.1062115 -1.175857 0.239652033604 0.4279891720       GAP43
# ENSG00000133636.11    6.40412     -0.0538755 0.0643832 -0.838723 0.401624836398           NA         NTS
# ENSG00000179776.19    5.39708      0.0711557 0.0628159  1.137073 0.255507604700           NA        CDH5


#  negative markers
all.gns[ all.gns$gene_name %in% diff.markers$negative, ]
#                     baseMean log2FoldChange     lfcSE      stat         pvalue         padj   gene_name
#                    <numeric>      <numeric> <numeric> <numeric>      <numeric>    <numeric> <character>
# ENSG00000181449.4    237.199      -0.497327 0.0981693  -5.05632 0.000000427417 0.0000114822        SOX2
# ENSG00000139352.4   5124.783       0.224521 0.0478075   4.69629 0.000002649356 0.0000562568       ASCL1
# ENSG00000124766.6   6077.275       0.240387 0.0515582   4.66239 0.000003125540 0.0000652705        SOX4
# ENSG00000198400.11  1330.858       0.164977 0.0824743   2.00036 0.045461686339 0.1390881863       NTRK1

#}}}



#  [circARID1A RBPs] are they differentially expressed in KD?
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(clusterProfiler)
library(DESeq2)
library(topGO)
library(org.Hs.eg.db)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load post processed data
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_kable-ready_genes+circRNAs_circARID1A SI_circARID1A SCR.RData')


#  load RBPs with enriched motifs on circARID1A and flanking introns
env<-new.env()
local(load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_RBPs_across_tissues.RData'), envir=env)
rbp.circ<-get('rbp', envir=env)
local(load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_introns_RBPs_across_tissues.RData'), envir=env)
rbp.introns<-get('rbp', envir=env)
rm(env)


#  RBPs on circARID1A
subset(all.gns[ all.gns$gene_name %in% rbp.circ$gene_name, ], padj<0.05)
#                     baseMean log2FoldChange     lfcSE      stat      pvalue        padj   gene_name
#                    <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric> <character>
# ENSG00000096746.17 11713.022      -0.200705 0.0476808  -4.20931 2.56154e-05 0.000380843     HNRNPH3
# ENSG00000165119.21 30936.874      -0.149580 0.0446670  -3.34878 8.11683e-04 0.006445560      HNRNPK
# ENSG00000169813.16  8344.073       0.138285 0.0452632   3.05514 2.24958e-03 0.014485922      HNRNPF
# ENSG00000136450.13 18568.764      -0.116550 0.0395603  -2.94613 3.21779e-03 0.019241095       SRSF1
# ENSG00000065978.19 52101.303       0.138769 0.0474655   2.92356 3.46050e-03 0.020380555        YBX1
# ENSG00000135482.7    956.437      -0.198371 0.0694044  -2.85795 4.26385e-03 0.023968092      ZC3H10


#  RBPs on flanking introns
subset(all.gns[ all.gns$gene_name %in% rbp.introns$gene_name, ], padj<0.05)
#                     baseMean log2FoldChange     lfcSE      stat      pvalue        padj   gene_name
#                    <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric> <character>
# ENSG00000162374.17 11360.620      -0.269105 0.0515720  -5.21796 1.80906e-07 5.35992e-06      ELAVL4
# ENSG00000123908.12 10515.701      -0.197879 0.0406857  -4.86355 1.15296e-06 2.73478e-05        AGO2
# ENSG00000107105.15  5383.996      -0.251295 0.0522324  -4.81096 1.50207e-06 3.44344e-05      ELAVL2
# ENSG00000073792.15  1218.883      -0.289174 0.0624731  -4.62804 3.69134e-06 7.43432e-05     IGF2BP2
# ENSG00000096746.17 11713.022      -0.200705 0.0476808  -4.20931 2.56154e-05 3.80843e-04     HNRNPH3
# ENSG00000070756.16 40157.631       0.176266 0.0419816   4.19864 2.68518e-05 3.97457e-04      PABPC1
# ENSG00000138385.16  7996.656      -0.217522 0.0521856  -4.16817 3.07050e-05 4.42487e-04         SSB
# ENSG00000102081.14  3362.315      -0.228151 0.0563009  -4.05218 5.07429e-05 6.65783e-04        FMR1
# ENSG00000265241.6  18027.174       0.158158 0.0392331   4.03124 5.54838e-05 7.18227e-04       RBM8A
# ENSG00000151846.8    212.872       0.399355 0.0997386   4.00011 6.33139e-05 8.03117e-04      PABPC3
# ENSG00000152518.8    717.260      -0.266837 0.0761551  -3.50314 4.59816e-04 4.11528e-03     ZFP36L2
# ENSG00000165119.21 30936.874      -0.149580 0.0446670  -3.34878 8.11683e-04 6.44556e-03      HNRNPK
# ENSG00000122566.21 79247.080      -0.147216 0.0474154  -3.10480 1.90406e-03 1.27305e-02   HNRNPA2B1
# ENSG00000169813.16  8344.073       0.138285 0.0452632   3.05514 2.24958e-03 1.44859e-02      HNRNPF
# ENSG00000120948.17  8488.070      -0.133240 0.0442570  -3.01057 2.60754e-03 1.61938e-02      TARDBP
# ENSG00000136450.13 18568.764      -0.116550 0.0395603  -2.94613 3.21779e-03 1.92411e-02       SRSF1
# ENSG00000151923.17  5194.067      -0.157565 0.0536929  -2.93450 3.34084e-03 1.97966e-02       TIAL1
# ENSG00000065978.19 52101.303       0.138769 0.0474655   2.92356 3.46050e-03 2.03806e-02        YBX1
# ENSG00000164548.11  4952.126      -0.142257 0.0494440  -2.87710 4.01347e-03 2.29138e-02       TRA2A
# ENSG00000197451.12 11303.409      -0.115518 0.0402498  -2.87001 4.10460e-03 2.33124e-02     HNRNPAB
# ENSG00000135482.7    956.437      -0.198371 0.0694044  -2.85795 4.26385e-03 2.39681e-02      ZC3H10
# ENSG00000066044.15  5630.948      -0.127883 0.0463322  -2.76012 5.77796e-03 3.04558e-02      ELAVL1
# ENSG00000120658.13   804.028      -0.189112 0.0718049  -2.63336 8.45449e-03 4.05529e-02       ENOX1
# ENSG00000004534.15  7774.560      -0.124319 0.0473855  -2.62355 8.70177e-03 4.13284e-02        RBM6

#}}}



#  [expression plots] KHSRP: 
#                         are validated targets differentially expressed in KD?
#                         are putative expressed targets enriched in DE genes?
#                         ecdfs of TPMs of the putative targets different across tumors, MYCN Tet-induction, circARID1A KD
#
#                     (DHX9, FUS, QKI, HNRNPA1, CELF6) + (MYC, MYCN, MYCL):
#
#                         across circARID1A KD samples
#                         across tumors
#                         across MYCN Tet-inducible samples
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_ecdfs.R')
source('~/bio/lib/grouplist2boxplot.R')
source('~/bio/lib/my_stripchart.R')


#  KHSRP validated targets 
#{{{

#  load DE results and check
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_kable-ready_genes+circRNAs_circARID1A SI_circARID1A SCR.RData')
g<-c('KHSRP', 'HNRNPA1', 'HNRNPF', 'HNRNPAB', 'HNRNPC', 'H3F3A', 'GNAS', 'SORBS1', 'CTNNB1', 'PPP2CA')
subset(all.gns, gene_name %in% g)
#                     baseMean log2FoldChange     lfcSE       stat     pvalue      padj   gene_name
#                    <numeric>      <numeric> <numeric>  <numeric>  <numeric> <numeric> <character>
# ENSG00000169813.16  8344.073     0.13828516 0.0452632  3.0551366 0.00224958 0.0144859      HNRNPF
# ENSG00000197451.12 11303.409    -0.11551776 0.0402498 -2.8700092 0.00410460 0.0233124     HNRNPAB
# ENSG00000113575.10  5511.476    -0.14079611 0.0499770 -2.8171868 0.00484463 0.0264537      PPP2CA
# ENSG00000092199.17 34311.693    -0.11978921 0.0473631 -2.5291634 0.01143348 0.0505710      HNRNPC
# ENSG00000135486.17 53014.615    -0.08949755 0.0481939 -1.8570291 0.06330701 0.1742323     HNRNPA1
# ENSG00000095637.22   517.183     0.12446871 0.0863685  1.4410004 0.14958456 0.3145780      SORBS1
# ENSG00000163041.10  7167.830     0.06493923 0.0544162  1.1933844 0.23271886 0.4199019       H3F3A
# ENSG00000168036.18 21437.299    -0.03696141 0.0472266 -0.7826399 0.43383861 0.6150481      CTNNB1
# ENSG00000088247.17  7397.816     0.07005215 0.1054463  0.6643581 0.50646113 0.6777866       KHSRP
# ENSG00000087460.24 37589.133    -0.00365316 0.0402182 -0.0908335 0.92762493 0.9598459        GNAS
rm(list=setdiff(ls(), l))

#}}}


#  gene TPMs for the tumors and the MYCN Tet-inducible samples
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_circRNAs+genes_all_libraries.RData')
a.hsa<-hsa
a.meta<-meta[ !is.na(risk_group) | grepl( 'SKNAS-TR-MYCN', bid) ]
a.meta<-a.meta[ treatment %in% '+Tet 4h', col:='orange3']      #  manually change 
a.meta<-a.meta[ treatment %in% 'ETOH 4h', col:='palegreen4']   #  manually change 
a.meta<-a.meta[ treatment %in% '+Tet 48h', col:='orangered3']  #  manually change 
a.meta<-a.meta[ treatment %in% 'ETOH 48h', col:='seagreen4']   #  manually change
a.tpm<-gns.tpm[, a.meta$bid ]
stopifnot( all.equal( a.meta[, bid], colnames(a.tpm) ) )
rm(list=setdiff(ls(), c(l, ls(pattern='a\\.'))))


#  gene TPMs for the circARID1A KD system
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_circRNAs+genes_all_libraries.RData')
b.meta<-meta[ grep('circARID1A', treatment) ]
b.tpm<-gns.tpm[, b.meta$bid ]
stopifnot( all.equal( b.meta[, bid], colnames(b.tpm) ) )
rm(list=setdiff(ls(), c(l, ls(pattern='b\\.'))))


#  load FIMO on transcriptome results
#  identify all putative KHSRP target transcripts (in order to compute KHSRP-related counts only we need to use fimo.txs.all)
#  make also a stricter list that involve only KHSRP motifs
#  summarize at the gene level (mean counts and mean densities across isoforms) in both cases
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/meme/global/results_transcripts.RData')
targets<-fimo.txs.all[ sapply(rbp_group, function(r){ any(grepl('KHSRP', r)) }) ][, 
    .(ncount=.N, density=.N/length*1e3, length=length, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), 
    seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(transcript_id, gene_id, gene_name)][, 
    .(ncount=mean(ncount), density=mean(density), length=mean(length), rbp_group=list(unique(unlist(rbp_group))), 
    motif_alt_id=list(unique(unlist(motif_alt_id)))), by=.(gene_id, gene_name)]
targets<-setNames(targets$gene_name, targets$gene_id)
targets.strict<-fimo.txs.all[ sapply(rbp_group, function(r){ all(grepl('KHSRP', r)) }) ][,
    .(ncount=.N, density=.N/length*1e3, length=length, motif_id=list(unique(unlist(motif_id))), motif_alt_id=list(unique(unlist(motif_alt_id))), 
    seq=list(unique(unlist(seq))), rbp_group=list(unique(unlist(rbp_group)))), by=.(transcript_id, gene_id, gene_name)][,
    .(ncount=mean(ncount), density=mean(density), length=mean(length), rbp_group=list(unique(unlist(rbp_group))), 
    motif_alt_id=list(unique(unlist(motif_alt_id)))), by=.(gene_id, gene_name)]
targets.strict<-setNames(targets.strict$gene_name, targets.strict$gene_id)
length(targets)         #  7852
length(targets.strict)  #  3766
rm(list=setdiff(ls(), c(l, 'targets', 'targets.strict')))


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  KHSRP putative targets
#{{{

#  [all expressed targets] ecdfs in tumors
x<-a.tpm[ rownames(a.tpm) %in% targets, a.meta[ !is.na(risk_group), bid]]
nrow(x<-x[ rowMeans(x)>=1, ])  #  5664
colnames(x)<-a.meta[ !is.na(risk_group), risk_group]
B<-lapply(split(split(log10(1+x), col(x)), colnames(x)), unlist)[ unique(colnames(x)) ]
#
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value  #  4.111e-85
wilcox.test(x=B[['HR_nMNA']], y=B[['LR']], alternative='less')$p.value   #  3.788e-182
#
my_ecdfs(B=B, B.CL=a.meta[ !is.na(risk_group), unique(col)], XLAB=expression(log[10](1+'TPM')), XLIM=c(0, 3), YLINE=4, YPADJ=-0.5, svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_KHSRP_putative_targets_TPMs_tumors.svg')


#  [strict targets] ecdfs in tumors
x<-a.tpm[ rownames(a.tpm) %in% targets.strict, a.meta[ !is.na(risk_group), bid]]
nrow(x<-x[ rowMeans(x)>=1, ])  #  2740
colnames(x)<-a.meta[ !is.na(risk_group), risk_group]
B<-lapply(split(split(log10(1+x), col(x)), colnames(x)), unlist)[ unique(colnames(x)) ]
#
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value  #  2.271e-41
wilcox.test(x=B[['HR_nMNA']], y=B[['LR']], alternative='less')$p.value   #  9.191e-110
#
my_ecdfs(B=B, B.CL=a.meta[ !is.na(risk_group), unique(col)], XLAB=expression(log[10](1+'TPM')), XLIM=c(0, 3), YLINE=4, YPADJ=-0.5)


#  ecdfs in MYCN Tet-induction
x<-a.tpm[ rownames(a.tpm) %in% targets, a.meta[ is.na(risk_group), bid]]
nrow(x<-x[ rowMeans(x)>=1, ])  #  4055
colnames(x)<-a.meta[ is.na(risk_group), treatment]
B<-lapply(split(split(log10(1+x), col(x)), colnames(x)), unlist)[ unique(colnames(x)) ]
#
wilcox.test(x=B[['+Tet 4h']], y=B[['ETOH 4h']], alternative='less')$p.value    #  0.8167
wilcox.test(x=B[['+Tet 48h']], y=B[['ETOH 48h']], alternative='less')$p.value  #  1.093e-08
#
my_ecdfs(B=B, B.CL=a.meta[ is.na(risk_group), unique(col)], XLAB=expression(log[10](1+'TPM')), XLIM=c(0, 3), YLINE=4, YPADJ=-0.5, svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_KHSRP_putative_targets_TPMs_MYCN_Tet-inducible.svg')


#  ecdfs in circARID1A KD
x<-b.tpm[ rownames(b.tpm) %in% names(targets), ]
nrow(x<-x[ rowMeans(x)>=1, ])  #  4794
colnames(x)<-b.meta[, treatment]
B<-lapply(split(split(log10(1+x), col(x)), colnames(x)), unlist)[ unique(colnames(x)) ]
#
wilcox.test(x=B[['circARID1A SI']], y=B[['circARID1A SCR']], alternative='less')$p.value  #  3.237653513e-12
#
my_ecdfs(B=B, B.CL=b.meta[, unique(col)], XLAB=expression(log[10](1+'TPM')), XLIM=c(0, 3), YLINE=4, YPADJ=-0.5, svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/ecdf_KHSRP_putative_targets_TPMs.svg')
rm(x, B)


#  enrichment of putative KHSRP expressed targets in DE genes
#{{{

#  load DE results
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_kable-ready_genes+circRNAs_circARID1A SI_circARID1A SCR.RData')
b.gns<-subset(all.gns, baseMean>0)
rm(list=setdiff(ls(), c(l, 'b.gns')))


#  KHSRP expressed targets 
x<-b.tpm[ rownames(b.tpm) %in% names(targets), ]
nrow(x<-x[ rowMeans(x)>=1, ])  #  4794
KHSRP<-rownames(x)
rm(x)


#  [significantly up-DE genes] Fisher's exact test:
#
#                    |      DE          |        not DE        |
#  ------------------|------------------|----------------------|
#      KHSRP target  |      x1          |          y1          |
#  ------------------|------------------|----------------------|
#  non-KHSRP target  |      x2          |          y2          | 
grp<-rownames( subset(b.gns, log2FoldChange>0 & padj<0.05) )
rst<-setdiff( rownames(b.gns), grp )    #  not significantly up-DE which will include significantly down-DE for example
stopifnot( length(grp) + length(rst)==nrow(b.gns) )
fisher.test(data.frame('in'=c(length(intersect(grp, KHSRP)), 
                              length(setdiff(grp, KHSRP))), 
                       'out'=c(length(intersect(rst, KHSRP)), 
                               length(setdiff(rst, KHSRP)))), alternative='greater')$p.value
#
#  => 4.879674212e-18



#  [significantly down-DE genes] Fisher's exact test:
#
#                    |      DE          |        not DE        |
#  ------------------|------------------|----------------------|
#      KHSRP target  |      x1          |          y1          |
#  ------------------|------------------|----------------------|
#  non-KHSRP target  |      x2          |          y2          | 
grp<-rownames( subset(b.gns, log2FoldChange<0 & padj<0.05) )
rst<-setdiff( rownames(b.gns), grp )    #  not significantly down-DE which will include significantly up-DE for example
stopifnot( length(grp) + length(rst)==nrow(b.gns) )
fisher.test(data.frame('in'=c(length(intersect(grp, KHSRP)), 
                              length(setdiff(grp, KHSRP))), 
                       'out'=c(length(intersect(rst, KHSRP)), 
                               length(setdiff(rst, KHSRP)))), alternative='greater')$p.value
#
#  => 0

#}}}

#}}}


#  [circARID1A KD] KHSRP 
#{{{

#  pick the genes of interest
GENES<-c('KHSRP')
GENES<-setNames(GENES, a.hsa$gene_id[ match(GENES, a.hsa$gene_name) ])


#  [circARID1A KD] boxplots across circARID1A KD
x<-b.tpm[ rownames(b.tpm) %in% names(GENES), b.meta[, bid], drop=F]
rownames(x)<-GENES[ match(rownames(x), names(GENES)) ]
x<-x[ order(rowMeans(x), decreasing=T), , drop=F]
colnames(x)<-b.meta[, treatment]
x<-setNames(lapply(unique(colnames(x)), function(y){ log10(1+x[, colnames(x) %in% y, drop=F]) }), unique(colnames(x)))
grouplist2boxplot(L=x, L.COL=b.meta[, unique(col)], XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/boxplot_', paste0(GENES, collapse='+'), '_circARID1A_KD.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  [circARID1A KD] stripchart 
x<-log10(1+b.tpm[ rownames(b.tpm) %in% names(GENES), b.meta[, bid], drop=F])
rownames(x)<-GENES[ match(rownames(x), names(GENES)) ]
x<-x[ order(rowMeans(x), decreasing=T), , drop=F]
colnames(x)<-b.meta[, treatment]
B<-list()
for(n in rownames(x)){
    B[[n]]<-do.call(cbind, lapply(split(x[n, ], factor(colnames(x), levels=unique(colnames(x)))), function(x){ unname(x) }))
}
my_stripchart(B, JITTER=2.5, fig.prefix=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/stripchart_', paste0(GENES, collapse='+'), '_circARID1A_KD'),
    list(scipen=0, par=list(mar=c(6.0, 10.0, 0.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4), 
    gap=3, ylab=expression(log[10]('TPM')), yline=6, ycex=2.4, xline=0, xcex=1.8, pt.cex=2.0, cl=b.meta[, unique(col)], legend='bottomleft',
    width=16, height=14))

#}}}


#  [tumors, MYCN Tet-inducible, circARID1A KD] DHX9, FUS, QKI, HNRNPA1, CELF6
#{{{

#  pick the genes of interest
GENES<-c('DHX9', 'FUS', 'QKI', 'HNRNPA1', 'CELF6')
GENES<-setNames(GENES, a.hsa$gene_id[ match(GENES, a.hsa$gene_name) ])


#  boxplots across tumors
x<-a.tpm[ rownames(a.tpm) %in% GENES, a.meta[ !is.na(risk_group), bid]]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-a.meta[ !is.na(risk_group), risk_group]
x<-setNames(lapply(unique(colnames(x)), function(y){ log10(1+x[, colnames(x) %in% y]) }), unique(colnames(x)))
grouplist2boxplot(L=x, L.COL=a.meta[ !is.na(risk_group), unique(col)], XLAS=2, XTEXT.ADJ=1, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, mar=c(9.5, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_across_tumors.svg'), width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  boxplots across MYCN Tet-inducible samples
x<-a.tpm[ rownames(a.tpm) %in% GENES, a.meta[ is.na(risk_group), bid]]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-a.meta[ is.na(risk_group), treatment]
x<-setNames(lapply(unique(colnames(x)), function(y){ log10(1+x[, colnames(x) %in% y]) }), unique(colnames(x)))
grouplist2boxplot(L=x, L.COL=a.meta[ is.na(risk_group), unique(col)], XLAS=2, XTEXT.ADJ=1, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, mar=c(9.5, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_MYCN_Tet-inducible.svg'), width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  boxplots across circARID1A KD
x<-b.tpm[ rownames(b.tpm) %in% names(GENES), b.meta[, bid]]
rownames(x)<-GENES[ match(rownames(x), names(GENES)) ]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-b.meta[, treatment]
x<-setNames(lapply(unique(colnames(x)), function(y){ log10(1+x[, colnames(x) %in% y]) }), unique(colnames(x)))
grouplist2boxplot(L=x, L.COL=b.meta[, unique(col)], XLAS=2, XTEXT.ADJ=1, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, mar=c(9.5, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/boxplot_', paste0(GENES, collapse='+'), '_circARID1A_KD.svg'), width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  [circARID1A KD] stripchart 
x<-log10(1+b.tpm[ rownames(b.tpm) %in% names(GENES), b.meta[, bid]])
rownames(x)<-GENES[ match(rownames(x), names(GENES)) ]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-b.meta[, treatment]
B<-list()
for(n in rownames(x)){
    B[[n]]<-do.call(cbind, lapply(split(x[n, ], factor(colnames(x), levels=unique(colnames(x)))), function(x){ unname(x) }))
}
my_stripchart(B, JITTER=2.5, fig.prefix=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/stripchart_', paste0(GENES, collapse='+'), '_circARID1A_KD'),
    list(scipen=0, par=list(mar=c(8.0, 8.0, 0.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4), 
    gap=3, ylab=expression(log[10]('TPM')), yline=4, ycex=2.4, xline=0, xcex=1.8, pt.cex=2.0, cl=b.meta[, unique(col)], legend='bottomleft',
    width=16, height=14))

#}}}


#  [tumors, MYCN Tet-inducible, circARID1A KD] MYC, MYCN, MYCL
#{{{

#  pick the genes of interest
GENES<-c('MYC', 'MYCN', 'MYCL')
GENES<-setNames(GENES, a.hsa$gene_id[ match(GENES, a.hsa$gene_name) ])


#  boxplots across tumors
x<-a.tpm[ rownames(a.tpm) %in% GENES, a.meta[ !is.na(risk_group), bid]]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-a.meta[ !is.na(risk_group), risk_group]
x<-setNames(lapply(unique(colnames(x)), function(y){ log10(1+x[, colnames(x) %in% y]) }), unique(colnames(x)))
grouplist2boxplot(L=x, L.COL=a.meta[ !is.na(risk_group), unique(col)], YMAX=3, XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_across_tumors.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  boxplots across MYCN Tet-inducible samples
x<-a.tpm[ rownames(a.tpm) %in% GENES, a.meta[ is.na(risk_group), bid]]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-a.meta[ is.na(risk_group), treatment]
x<-setNames(lapply(unique(colnames(x)), function(y){ log10(1+x[, colnames(x) %in% y]) }), unique(colnames(x)))
grouplist2boxplot(L=x, L.COL=a.meta[ is.na(risk_group), unique(col)], XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_MYCN_Tet-inducible.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  boxplots across circARID1A KD
x<-b.tpm[ rownames(b.tpm) %in% names(GENES), b.meta[, bid]]
rownames(x)<-GENES[ match(rownames(x), names(GENES)) ]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-b.meta[, treatment]
x<-setNames(lapply(unique(colnames(x)), function(y){ log10(1+x[, colnames(x) %in% y]) }), unique(colnames(x)))
grouplist2boxplot(L=x, L.COL=b.meta[, unique(col)], XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/boxplot_', paste0(GENES, collapse='+'), '_circARID1A_KD.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  [circARID1A KD] stripchart 
x<-log10(1+b.tpm[ rownames(b.tpm) %in% names(GENES), b.meta[, bid]])
rownames(x)<-GENES[ match(rownames(x), names(GENES)) ]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-b.meta[, treatment]
B<-list()
for(n in rownames(x)){
    B[[n]]<-do.call(cbind, lapply(split(x[n, ], factor(colnames(x), levels=unique(colnames(x)))), function(x){ unname(x) }))
}
my_stripchart(B, JITTER=2.5, fig.prefix=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/stripchart_', paste0(GENES, collapse='+'), '_circARID1A_KD'),
    list(scipen=0, par=list(mar=c(8.0, 8.0, 0.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4), 
    gap=3, ylab=expression(log[10]('TPM')), yline=4, ycex=2.4, xline=0, xcex=1.8, pt.cex=2.0, cl=b.meta[, unique(col)], legend='bottomleft',
    width=16, height=14))

#}}}


#  [tumors, MYCN Tet-inducible, circARID1A KD] DHX9, ILF3
#{{{

#  pick the genes of interest
GENES<-c('DHX9', 'ILF3')
GENES<-setNames(GENES, a.hsa$gene_id[ match(GENES, a.hsa$gene_name) ])


#  boxplots across tumors
x<-a.tpm[ rownames(a.tpm) %in% GENES, a.meta[ !is.na(risk_group), bid]]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-a.meta[ !is.na(risk_group), risk_group]
x<-setNames(lapply(unique(colnames(x)), function(y){ log10(1+x[, colnames(x) %in% y]) }), unique(colnames(x)))
grouplist2boxplot(L=x, L.COL=a.meta[ !is.na(risk_group), unique(col)], YMAX=3, XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_across_tumors.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  boxplots across MYCN Tet-inducible samples
x<-a.tpm[ rownames(a.tpm) %in% GENES, a.meta[ is.na(risk_group), bid]]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-a.meta[ is.na(risk_group), treatment]
x<-setNames(lapply(unique(colnames(x)), function(y){ log10(1+x[, colnames(x) %in% y]) }), unique(colnames(x)))
grouplist2boxplot(L=x, L.COL=a.meta[ is.na(risk_group), unique(col)], XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=5, mar=c(2.0, 9.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', paste0(GENES, collapse='+'), '_MYCN_Tet-inducible.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  boxplots across circARID1A KD
x<-b.tpm[ rownames(b.tpm) %in% names(GENES), b.meta[, bid]]
rownames(x)<-GENES[ match(rownames(x), names(GENES)) ]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-b.meta[, treatment]
x<-setNames(lapply(unique(colnames(x)), function(y){ log10(1+x[, colnames(x) %in% y]) }), unique(colnames(x)))
grouplist2boxplot(L=x, L.COL=b.meta[, unique(col)], XLAS=0, XTEXT.ADJ=0.5, YLAB=expression(log[10](1+'TPM')), XLAB.CEX=2.4, YLAB.CEX=2.4, YTEXT.LINE=3, mar=c(2.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/boxplot_', paste0(GENES, collapse='+'), '_circARID1A_KD.svg'), width=14, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  [circARID1A KD] stripchart 
x<-log10(1+b.tpm[ rownames(b.tpm) %in% names(GENES), b.meta[, bid]])
rownames(x)<-GENES[ match(rownames(x), names(GENES)) ]
x<-x[ order(rowMeans(x), decreasing=T), ]
colnames(x)<-b.meta[, treatment]
B<-list()
for(n in rownames(x)){
    B[[n]]<-do.call(cbind, lapply(split(x[n, ], factor(colnames(x), levels=unique(colnames(x)))), function(x){ unname(x) }))
}
my_stripchart(B, JITTER=2.5, fig.prefix=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/stripchart_', paste0(GENES, collapse='+'), '_circARID1A_KD'),
    list(scipen=0, par=list(mar=c(8.0, 8.0, 0.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4), 
    gap=3, ylab=expression(log[10]('TPM')), yline=4, ycex=2.4, xline=0, xcex=1.8, pt.cex=2.0, cl=b.meta[, unique(col)], legend='bottomleft',
    width=16, height=14))

#}}}

#}}}



#  [miR-138-5p targets] check for enrichment in DE genes
#                       ecdf of their TPMs across conditions
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_ecdfs.R')
source('~/bio/lib/grouplist2boxplot.R')


#  gene TPMs for the circARID1A KD system
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_circRNAs+genes_all_libraries.RData')
c.meta<-meta[ grep('circARID1A', treatment) ]
c.tpm<-gns.tpm[, c.meta$bid ]
stopifnot( all.equal( c.meta[, bid], colnames(c.tpm) ) )
rm(list=setdiff(ls(), c(l, 'hsa', ls(pattern='c\\.'))))


#  we queried miRWalk on 2020-05-05
#  but we keep only the validated targets 
#  add gene_id to miR-138-5p targets (fix the two with old names...)
m138<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/miRWalk_miR-138-5p_target_sites.csv', header=T, sep=',')
m138<-m138[ grep('MIRT', validated) ]
m138<-m138[, .(counts=.N), by=.(genesymbol)][ order(-counts) ]
colnames(m138)<-c('gene_name', 'counts')
m138$gene_id<-hsa$gene_id[ match(m138$gene_name, hsa$gene_name) ]  
m138[ is.na(gene_id), gene_name]  #  FAM35A, FAM109A
m138[ gene_name %in% 'FAM35A', gene_name:='SHLD2']
m138[ gene_name %in% 'FAM109A', gene_name:='PHETA1']
m138$gene_id<-hsa$gene_id[ match(m138$gene_name, hsa$gene_name) ]  
m138[ is.na(gene_id), gene_name]  #  empty


#  load DE gene results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/DESeq2_genes_circARID1A SI_circARID1A SCR.RData')
up.gns<-subset(RES, padj<0.05 & log2FoldChange>0)
down.gns<-subset(RES, padj<0.05 & log2FoldChange<0)


#  [upregulated] Fisher's exact test:
#
#                   |      up          |        not up        |
#  -----------------|------------------|----------------------|
#       targets     |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#     not-targets   |      x2          |          y2          | 
grp<-up.gns$gene_name
rst<-subset(RES, ! gene_name %in% up.gns$gene_name)$gene_name
stopifnot( length(grp)+length(rst)==nrow(RES) )
fisher.test(data.frame('in'=c(length(intersect(grp, m138[, gene_name])), 
                              length(setdiff(grp, m138[, gene_name]))), 
                       'out'=c(length(intersect(rst, m138[, gene_name])), 
                               length(setdiff(rst, m138[, gene_name])))), alternative='greater')$p.value
#  => 0.03162582771


#  [downregulated] Fisher's exact test:
#
#                   |     down         |        not down      |
#  -----------------|------------------|----------------------|
#       targets     |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#     not-targets   |      x2          |          y2          | 
grp<-down.gns$gene_name
rst<-subset(RES, ! gene_name %in% down.gns$gene_name)$gene_name
stopifnot( length(grp)+length(rst)==nrow(RES) )
fisher.test(data.frame('in'=c(length(intersect(grp, m138[, gene_name])), 
                              length(setdiff(grp, m138[, gene_name]))), 
                       'out'=c(length(intersect(rst, m138[, gene_name])), 
                               length(setdiff(rst, m138[, gene_name])))), alternative='greater')$p.value
#  => 6.512017977e-08


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  ecdfs of expressed miR-138-5p targets
x<-c.tpm[ rownames(c.tpm) %in% m138[, gene_id], ]
x<-x[ rowMeans(x)>=1, ]
colnames(x)<-c.meta[, treatment]
B<-lapply(split(split(log10(1+x), col(x)), colnames(x)), unlist)[ unique(colnames(x)) ]
#
wilcox.test(x=B[['circARID1A SI']], y=B[['circARID1A SCR']], alterantive='less')$p.value  #  0.4096756754
#
par(mar=c(6.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
my_ecdfs(B=B, B.CL=c.meta[, unique(col)], XLAB=expression(log[10](1+'TPM')), YLINE=4, YPADJ=-0.5, svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/ecdf_miR-138-5p_validated_targets_TPMs.svg')
rm(x, B)

#}}}



#  [circARID1A siRNA vs scrambled] differential splicing analysis
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(clusterProfiler)
library(edgeR)
library(org.Hs.eg.db)
library(RColorBrewer)
library(pheatmap)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_plotSpliceDGE.R')


#  functions
#{{{

do_ds<-function(DGE, CT, lfcThreshold, plot.top.ds=F){
    #            DGE : DGEList()
    #             CT : metadata data.frame of only samples belonging to the two conditions compared
    #   lfcThreshold : log2FC threshold to use for the differential exon expression analysis
    #    plot.top.ds : FALSE or numeric indicating how many plots to do for the top differentially spliced genes
    require(clusterProfiler)
    require(edgeR)


    #  save open graphics devices by start
    DEV_START<-dev.list()


    #  factor level order defines the CND
    CND<-paste0(rev(levels(DGE$samples$group)), collapse=' vs ')


    #  keep exons with CPM>1 in at least 3 samples
    cat('\nkeeping exons with CPM>1 in at least 3 samples.\nnumber of exons that FAIL/PASS our filtering criterion:\n')
    print(summary(keep<-rowSums( cpm(DGE)>1 )>=3))
    DGE<-DGE[ keep, , keep.lib.sizes=F]
    rm(keep)


    #  calculate the normalization factors 
    #  define the model
    #  estimate the dispersions
    cat('\nestimating the dispersions...')
    DGE<-calcNormFactors(DGE)
    rg<-DGE$samples$group
    design<-model.matrix(~rg)
    colnames(design)<-sub('rg', '', colnames(design))
    DGE<-estimateDisp(DGE, design, robust=T)
    rm(rg)
    #x11()
    #plotBCV(DGE)


    #  fit the quasi-likelihood model
    cat('\nfitting the quasi-likelihood model...')
    FIT<-glmQLFit(DGE, design, robust=T)
    #x11()
    #plotQLDisp(FIT)


    #  test for differential exon expression above threshold
    cat('\ndoing the differential exon expression analysis...')
    DE<-glmTreat(FIT, coef=2, lfc=lfcThreshold)
    DE.RES<-as.data.frame(topTags(DE, n=Inf, p.value=0.05))


    #  test for alternative splicing using gene-level testing (more powerful for cases where several exons are differentially spliced)
    cat('\ndoing the differential splicing analysis...\n')
    DS<-diffSpliceDGE(FIT, coef=2, geneid='gene_id', exonid='exon_id')
    DS.RES<-topSpliceDGE(DS, test='gene', n=Inf, FDR=0.05)


    #  MSigDB C2 enrichment test of alternatively spliced genes
    if (nrow(DS.RES)>0){
        c2<-read.gmt('/data/annotation/MSigDB/c2.all.v7.0.symbols.gmt')  #  load C2 gene set of MSigDB
        DS.c2<-enricher(DS.RES$gene_name, pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=1, maxGSSize=200, universe=DS$gene.genes$gene_name, TERM2GENE=c2)
    } else {
        DS.c2<-NA
    }


    #  plot and save top differentially spliced genes 
    if( plot.top.ds>0 ){
        for(i in 1:min(c(plot.top.ds, nrow(DS.RES)))){
            my_plotSpliceDGE(DS, GENEID=DS.RES$gene_name[i], GENECOL='gene_name', FDR=0.05, CND=CND, FILEROOT='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/plotDS_')
        }
        rm(i)
    }


    #  MA-plot for the differentially expressed exons
    x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
    RES<-data.frame(AveLogCPM=DE$AveLogCPM, log2FoldChange=DE$table$logFC, sig=F, gene_id=DE$genes$gene_id, gene_name=DE$genes$gene_name, row.names=rownames(DE$table))
    RES$sig[ rownames(RES) %in% rownames(DE.RES) ]<-T
    YLIM<-range(pretty(range(RES$log2FoldChange, na.rm=T), 5))
    XLIM<-range(pretty(range(RES$AveLogCPM, na.rm=T), 5))
    par(mar=c(4.5, 7.5, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4)
    options(scipen=-3)
    geneplotter::plotMA(RES[, 1:3], xlab='Average log-CPM', ylab='', main='', ylim=YLIM, xlim=XLIM, cex=1.2, log='', colSig='red3')
    options(scipen=+10)
    mtext(CND, side=3, line=0, padj=+1.5, cex=2.4)
    mtext(expression(log[2]~'fold change'), side=2, line=4, padj=-0.2, cex=2.4, las=3)
    if(lfcThreshold>0){
        abline(h=c(-lfcThreshold, lfcThreshold), lty=1, lwd=4, col='cyan4')
    }
    #identify(RES$AveLogCPM, RES$log2FoldChange, labels=RES$gene_name, cex=0.9, offset=0.2, xpd=T)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/plotMA_exons_', sub(' vs ', '_', CND), '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    rm(YLIM,XLIM)


    #  save results
    cat('\nsaving the results...')
    FILE<-paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/edgeR_exons_', sub(' vs ', '_', CND), '.RData')
    save(CND, CT, DGE, FIT, DE, DE.RES, DS, DS.RES, DS.c2, file=FILE, compress=T)


    #  save the differential splicing results as TSV/XLSX 
    x<-as.data.frame(DS.RES)
    rownames(x)<-NULL
    x$gene.F<-NULL
    colnames(x)<-c('gene_id', 'gene_name', 'entrez', 'chr', 'strand', 'n_exons', 'pvalue', 'padj')
    write.table(x, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/edgeR_genes_DS_', sub(' vs ', '_', CND), '.tsv'), quote=F, sep='\t', row.names=F, col.names=T) 
    write.xlsx(x, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/edgeR_genes_DS_', sub(' vs ', '_', CND), '.xlsx'), row.names=F, col.names=T, sheetName=gsub(' ', '_', CND)) 
    rm(x)


    #  save as TSV the MSigDB C2 enrichment results as well
    if(nrow(DS.RES)>0){
        x<-DS.c2@result
        rownames(x)<-NULL
        x$Description<-NULL
        x$qvalue<-NULL
        colnames(x)<-c('msigdb_id', 'geneRatio', 'backgroundRatio', 'pvalue', 'padj', 'gene_names', 'count')
        write.table(x, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/edgeR_genes_DS_MSigDB_C2_', sub(' vs ', '_', CND), '.tsv'), quote=F, sep='\t', row.names=F, col.names=T) 
    }


   #  close all figures opened by the script
   for (n in setdiff(dev.list(), DEV_START)){ dev.off(n) }


   #  return the filename of the saved data
   cat('\n')
   invisible(FILE)
}

#}}}


#  load reference and add as many Entrez ids as possible
#{{{

#  import the reference 
#  keep genes only 
#  remove chrM, chrY
#  define gene_ids without subversion suffixes (there are _PAR_Y suffixes after the subversions as well)
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[hsa$type %in% 'gene']
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_name', 'gene_type')]
hsa$no_subversion<-sub('\\..*$', '', hsa$gene_id)
hsa<-hsa[ grep('chrY|chrM', seqnames(hsa), invert=T) ]


#  get the gene_ids
gid<-data.table(data.frame(mcols(hsa)[, c('gene_id', 'no_subversion', 'gene_name', 'gene_type')]))


#  identify by the gene_id and/or by the gene_name those genes with an Entrez id 
summary( exist<-gid$no_subversion %in% mappedRkeys(org.Hs.egENSEMBL2EG) | gid$gene_name %in% mappedRkeys(org.Hs.egSYMBOL2EG) )
#
#   FALSE    TRUE 
#   19872   38397 
#
gid[!exist, sort(table(gene_type), decreasing=T)]                                          #  510 protein_coding, 4961 lincRNA, ...
gid[!exist & gene_type %in% 'protein_coding' & !grepl('^[AFB][CLPDFX][0-9]', gene_name) ]  #  majority of protein-coding have provisional names 
gid<-gid[exist, ]
rm(exist)


#  merge the gene_id <-> entrez table with the entrez <-> gene_name table 
#  keep all Entrez ids without a gene_id but with a gene_name
g2e<-data.table(toTable(org.Hs.egENSEMBL2EG))
colnames(g2e)<-c('entrez', 'no_subversion')
e2s<-data.table(toTable(org.Hs.egSYMBOL2EG))
colnames(e2s)<-c('entrez', 'gene_name')
g2e<-g2e[ e2s, on='entrez']  
rm(e2s)


#  keep only the Entrez ids for those genes in our reference that can be identified by gene_id or gene_name
g2e<-g2e[ no_subversion %in% gid$no_subversion | gene_name %in% gid$gene_name ]


#  collect all Entrez ids under a given gene_id for those gene_ids found in our reference
#  update the Entrez ids for our reference 
x<-g2e[ no_subversion %in% gid$no_subversion ][, .(entrez=list(unique(entrez))), by=.(no_subversion)]
gid<-x[ gid, on='no_subversion']
rm(x)


#  collect all Entrez ids under a given gene_name for those gene_names found in our reference that did not have a gene_id in org.Hs.eg.db
#  keep only those gene_names without an Entrez id yet in our reference
y<-gid[ lengths(entrez)==0 ][, entrez:=NULL]
x<-g2e[ ! no_subversion %in% y$no_subversion ][, .(entrez=list(unique(entrez))), by=.(gene_name)][ gene_name %in% y$gene_name ]
y<-x[ y, on='gene_name' ]
stopifnot( all.equal( gid[ match(y$no_subversion, no_subversion),  c('no_subversion', 'gene_name')] , y[, c('no_subversion', 'gene_name')] ) )
gid[ match(y$no_subversion, no_subversion), entrez:=y$entrez ]
rm(x,y,g2e)


#  update the reference 
hsa$entrez<-vector('list', length(hsa))
hsa[ match(gid$no_subversion, hsa$no_subversion) ]$entrez<-gid$entrez
rm(gid)

#}}}


#  load metadata
#  keep circARID1A KD related samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/metadata.RData')
meta<-meta[ grep('ARID1A', bid) ]


#  load exon counts
#  keep samples of interest and order them according to metadata
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/featureCounts_exons.RData')
exns<-exns[ meta$bid ]
exns<-lapply(exns, '[[', 'exons')
gc()


#  isolate the fractional counts for each sample 
#  keep only exons/genes of the reference (removed chrY, chrM)
CNTS<-do.call(cbind, lapply(exns, function(x){ x[, 'counts'] }))
colnames(CNTS)<-names(exns)
CNTS<-cbind(exns[[1]][, 1:6], CNTS)  #  add the gene metadata
CNTS<-CNTS[ gene_id %in% hsa$gene_id ]
rm(exns)
gc()


#  identify genes with an Entrez id keeping the first entry among multiple entries
ez<-data.table(data.frame(mcols(hsa)[, c('gene_id', 'entrez')]))
ez<-ez[, .(entrez=entrez[[1]][1]), by=.(gene_id)]   #  incidentaly this removes gene_ids with no Entrez id
CNTS<-ez[CNTS, on='gene_id']                        #  but we want to keep them with NAs for Entrez id
rm(ez)


#  add gene_names and unique exon_ids
n<-data.table(data.frame(mcols(hsa)[, c('gene_id', 'gene_name')]))
CNTS<-n[CNTS, on='gene_id']
n<-CNTS[, 'gene_id', with=F][, exon_id:=paste0(gene_id, '_ex_', 1:.N), by=.(gene_id)] 
stopifnot( all.equal( n$gene_id, CNTS$gene_id ) )
CNTS<-cbind(n[, 'exon_id'], CNTS)  #  n[CNTS, on='gene_id'] fails because n$gene_id key is not unique anymore
rm(n)



#  [circARID1A siRNA vs scrambled]
CT<-meta
DGE<-DGEList(counts=ceiling(as.data.frame(CNTS[, CT$bid, with=F])), genes=CNTS[, 1:9], group=factor(CT$treatment, levels=c('circARID1A SCR', 'circARID1A SI')), samples=data.frame(CT[, c('gc', 'nreads', 'features', 'no_features', 'unmapped', 'p_unmapped'), with=F], row.names=CT$bid))
FILE<-do_ds(DGE, CT, lfcThreshold=log2(1.0), plot.top.ds=6)


#  further analysis and visualization of the results
#{{{

#  load the DS results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/edgeR_exons_circARID1A SI_circARID1A SCR.RData')


#  ARID1A differential exon usage
my_plotSpliceDGE(DS, GENEID='ARID1A', GENECOL='gene_name', FDR=0.05, CND=CND, FILEROOT='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/plotDS_')


#  EIF4G2 differential exon usage
my_plotSpliceDGE(DS, GENEID='EIF4G2', GENECOL='gene_name', FDR=0.05, CND=CND, FILEROOT='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/plotDS_')


#  HDAC9 differential exon usage
my_plotSpliceDGE(DS, GENEID='HDAC9', GENECOL='gene_name', FDR=0.05, CND=CND, FILEROOT='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/plotDS_')


#  DROSHA differential exon usage
my_plotSpliceDGE(DS, GENEID='DROSHA', GENECOL='gene_name', FDR=0.05, CND=CND, FILEROOT='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/plotDS_')


#  LINC00632 differential exon usage
my_plotSpliceDGE(DS, GENEID='LINC00632', GENECOL='gene_name', FDR=0.05, CND=CND, FILEROOT='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/plotDS_')


#  RPL37A differential exon usage
my_plotSpliceDGE(DS, GENEID='RPL37A', GENECOL='gene_name', FDR=0.05, CND=CND, FILEROOT='')

#}}}

#}}}




######################
#
#
#  proliferative index
#  ADRN/MES index
#
#
######################




#  [PI] calculate the proliferative index based on TPMs which removes biases on gene length
#
#       N.B. there is no difference across treatments! 
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


#  load reference to convert gene_id to gene_name
#  remove chrM/chrY counts
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene']
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]
seqlevels(hsa, pruning.mode='coarse')<-seqlevels(hsa)[ grep('chrM|chrY', seqlevels(hsa), invert=T) ]


#  load metadata and keep circARID1A KD related samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/metadata.RData')
meta<-meta[ grepl('circARID1A', treatment) ]


#  load featureCounts 
#  remove chrM/chrY counts
#  compute TPMs
#  convert to matrix
#  log2-transform to regularize the range
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/featureCounts.RData')
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


#  keep TPMs of marker genes
#  order according to metadata
gns.tpm<-gns.tpm[rownames(gns.tpm) %in% ProliferativeIndex:::metaPCNA2, ]


#  add proliferative index to metadata
meta$PI<-apply(gns.tpm, 2, median)


#  cluster samples by PI
#  cut at maximum silhouette width 
#  add discrete PI classification to the metadata
x<-as.matrix(setNames(meta$PI, meta$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
summary(silhouette(cutree(x.hc, k=2), dist=x.d))  #  best
summary(silhouette(cutree(x.hc, k=3), dist=x.d))
w<-cutree(x.hc, k=2)
o<-sapply(split(names(w), w), function(n){ mean(x[n, ,drop=F]) })
for(n in names(o)){
    w[ w==n ]<-o[n]
}
w<-sort(w, decreasing=F)
x<-x[ names(w), ,drop=F ]
x.d<-dist(x, method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
w<-cutree(x.hc, k=2)
w[ w %in% '1' ]<-'low'
w[ w %in% '2' ]<-'high'
meta$PI.class<-factor(w[ meta$bid ], levels=c('low', 'high'))
meta$PI.class.col<-factor(w[ meta$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  save
save(meta, gns.tpm, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/proliferative_index.RData')


#  visualize the classification results
#{{{

#  load back the results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/proliferative_index.RData')


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  heatmap of annotated samples clustered by Euclidean distance in log2-TPMs
x<-gns.tpm
stopifnot( all.equal( colnames(x), meta$bid ) )
x.ex<-data.frame(Treatment=meta$treatment, PI=meta$PI.class, row.names=meta$bid)
x.cl<-setNames(list(setNames( unique(meta$col), unique(x.ex$Treatment)), setNames(levels(meta$PI.class.col), levels(meta$PI.class))), colnames(x.ex))
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
rm(x, x.ex, x.cl, x.hc, ph)

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/proliferative_index.RData



#  [ADRN/MES] expression among conditions
#             calculate the adrenergic/mesenchymal index based on TPMs
#
#             N.B. there is no difference across treatments! 
#
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
    #  internal function that depends on many global variables and plots a nice boxplot 

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


#  load metadata and keep circARID1A KD related samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/metadata.RData')
meta<-meta[ grepl('circARID1A', treatment) ]


#  load featureCounts 
#  remove chrM/chrY counts
#  compute TPMs
#  convert to matrix
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/featureCounts.RData')
gns.tpm<-unlist(totalrna[ meta[, bid] ])
gns.tpm$bid<-sub('\\.[0-9]*$', '', rownames(gns.tpm))
rownames(gns.tpm)<-NULL
gns.tpm<-gns.tpm[ gns.tpm$gene_id %in% hsa$gene_id, ]
gns.tpm$gene_name<-hsa$gene_name[ match(gns.tpm$gene_id, hsa$gene_id) ]
gns.tpm<-data.table(gns.tpm)[, .(gene_name=gene_name, length=length, counts=counts, tpm=counts/length/sum(counts/length)*1e6), by=.(bid)]
gns.tpm<-dcast(gns.tpm, bid ~ gene_name, value.var='tpm', fun.aggregate=sum)
gns.tpm<-t(data.frame(gns.tpm[, -1], row.names=gns.tpm[, bid], check.names=F))
rm(totalrna, hsa)


#  keep TPMs of marker genes
#  order according to metadata
gns.tpm<-gns.tpm[rownames(gns.tpm) %in% mes.adr$gene_name, ]


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(14.0, 7.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B.cl<-setNames(meta[, unique(col)], meta[, unique(treatment)])


#  expression plots of mesenchymal/adrenergic markers 
#{{{

#  split to mesenchymal/adrenergic 
x<-log2(1+gns.tpm)
colnames(x)<-meta[ match(colnames(x), bid), treatment ]
m<-x[ rownames(x) %in% mes.adr[ cell_type %in% 'MES', gene_name ], ]
a<-x[ rownames(x) %in% mes.adr[ cell_type %in% 'ADRN', gene_name ], ]
stopifnot( nrow(a)+nrow(m)==nrow(x) )
B.m<-setNames(lapply(meta[, unique(treatment)], function(r){ c(m[, colnames(m) %in% r]) }), meta[, unique(treatment)])
B.a<-setNames(lapply(meta[, unique(treatment)], function(r){ c(a[, colnames(a) %in% r]) }), meta[, unique(treatment)])
rm(x)


#  all pairwise Mann-Whitney U tests
n<-matrix(names(B.m)[ combn(length(B.m), 2) ], nrow=2)
p.m<-p.a<-setNames(rep(NA, ncol(n)), apply(n, 2, paste0, collapse='-'))
for(i in seq_len(ncol(n))){
    p.m[i]<-wilcox.test(B.m[[ n[1, i] ]], B.m[[ n[2, i] ]], alternative='two.sided')$p.value
    p.a[i]<-wilcox.test(B.a[[ n[1, i] ]], B.a[[ n[2, i] ]], alternative='two.sided')$p.value
}
#p.m<-p.m[ p.m<0.05 ]
#p.a<-p.a[ p.a<0.05 ]
rm(n, i)


#  mesenchymal ecdf
B<-B.m
N<-lengths(B)/nrow(m)
B.cl<-B.cl[ names(B) ]
my_ecdfs(B, B.cl, XLIM=c(0.0, 8), XLAB=expression(log[2](1+'TPM')), MAIN=paste(nrow(m), 'mesenchymal markers'), LTY=1, LWD=8, svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/ecdf_TPMs_for_mesenchymal_all.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#
#  boxplot
#
#YTICK<-pretty(c(0.0, sapply(B, max)), 5)
YTICK<-pretty(c(0.0, 10), 4)
par(mar=c(4.5, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=1.5, add=T)
mtext(expression(log[2](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
#mtext(text=paste0(sub('circARID1A ', '', names(B)), ' (', N, ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=sub('circARID1A ', '', names(B)), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=paste(nrow(m), 'mesenchymal markers'), side=3, line=-1, padj=+0.2, cex=1.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/boxplot_TPMs_for_mesenchymal_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  adrenergic boxplot
B<-B.a
N<-lengths(B)/nrow(a)
B.cl<-B.cl[ names(B) ]
my_ecdfs(B, B.cl, XLIM=c(0.0, 8), XLAB=expression(log[2](1+'TPM')), MAIN=paste(nrow(a), 'adrenergic markers'), LTY=1, LWD=8, svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/ecdf_TPMs_for_adrenergic_all.svg', mar=c(6.0, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#
#  boxplot
#
#YTICK<-pretty(c(0.0, sapply(B, max)), 5)
YTICK<-pretty(c(0.0, 10), 4)
par(mar=c(4.5, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=1.5, add=T)
mtext(expression(log[2](1+'TPM')), side=2, line=3, padj=-0.1, las=0, cex=2.4)
#mtext(text=paste0(sub('circARID1A ', '', names(B)), ' (', N, ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=sub('circARID1A ', '', names(B)), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=paste(nrow(a), 'adrenerigc markers'), side=3, line=-1, padj=+0.2, cex=1.8)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/boxplot_TPMs_for_adrenergic_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  adrenergic/mesenchymal score similar to https://www.nature.com/articles/ng.3899#Sec2:
#
#      In each sample the genes with TPM>=1 were ranked by increasing TPM order. The median percentile of the ADRN and MES list of genes
#      corresponded to the ADRN and MES scores (highest score indicates higher overall expresssion for the corresonding group).
#{{{

#  compute ADRN/MES scores
AM<-t(apply(gns.tpm, 2, function(tpms){ 
    keep<-tpms>=1
    gn<-rownames(gns.tpm)[ keep ]
    tpms<-tpms[ keep ]
    gn<-gn[ order(tpms, decreasing=F) ]
    c(adrn.score=median( which(gn %in% mes.adr[ cell_type %in% 'ADRN', gene_name])/length(gn) ),
      mes.score=median( which(gn %in% mes.adr[ cell_type %in% 'MES', gene_name])/length(gn) ) 
    )
    }))


#  add relevant metadata
AM<-data.table(bid=rownames(AM), AM)
meta<-meta[, c('bid', 'treatment', 'col'), with=F]
AM<-meta[ AM, on='bid']
rm(meta)


#  cluster by MES score 
x<-as.matrix(setNames(AM$mes.score, AM$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
summary(silhouette(cutree(x.hc, k=2), dist=x.d))  #  best  
summary(silhouette(cutree(x.hc, k=3), dist=x.d))
w<-cutree(x.hc, k=2)
o<-sapply(split(names(w), w), function(n){ mean(x[n, ,drop=F]) })
for(n in names(o)){
    w[ w==n ]<-o[n]
}
w<-sort(w, decreasing=F)
x<-x[ names(w), ,drop=F ]
x.d<-dist(x, method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
w<-cutree(x.hc, k=2)
w[ w %in% '1' ]<-'low'
w[ w %in% '2' ]<-'high'
AM$mes.class<-factor(w[ AM$bid ], levels=c('low', 'high'))
AM$mes.class.col<-factor(w[ AM$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  cluster by ADRN score 
x<-as.matrix(setNames(AM$adrn.score, AM$bid), check.names=F)
x.d<-dist(x, method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
summary(silhouette(cutree(x.hc, k=2), dist=x.d))
summary(silhouette(cutree(x.hc, k=3), dist=x.d))  #  best but not convincing...
w<-cutree(x.hc, k=2)
o<-sapply(split(names(w), w), function(n){ mean(x[n, ,drop=F]) })
for(n in names(o)){
    w[ w==n ]<-o[n]
}
w<-sort(w, decreasing=F)
x<-x[ names(w), ,drop=F ]
x.d<-dist(x, method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
w<-cutree(x.hc, k=2)
w[ w %in% '1' ]<-'low'
w[ w %in% '2' ]<-'high'
AM$adrn.class<-factor(w[ AM$bid ], levels=c('low', 'high'))
AM$adrn.class.col<-factor(w[ AM$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))


#  save
save(AM, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/adrenergic-mesenchymal_score.RData')

#}}}


#  visualize the classification results
#{{{

#  load back the results 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/adrenergic-mesenchymal_score.RData')


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  boxplot of scores 
B<-data.frame(AM[, c('treatment', 'adrn.score', 'mes.score'), with=F])
colnames(B)<-c('treatment', 'ADRN', 'MES')
B<-lapply(split(B[, c('ADRN', 'MES')], factor(B$treatment, levels=AM[, unique(treatment)])), t)
B.cl<-setNames(AM[, unique(col)], AM[, unique(treatment)])[names(B)]
grouplist2boxplot(L=B, L.COL=B.cl, YLAB='Score', YLAB.CEX=2.4, XLAB.CEX=2.4, YTEXT.LINE=5, LEGEND='topright', mar=c(6.0, 7.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/figures/boxplot_ADRN_MES_scores.svg', width=16, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B, B.cl)


#  [ADRN/MES scores together] heatmap of annotated samples clustered by Euclidean distance in ADRN/MES scores
x<-data.frame(AM[, -1], row.names=AM[, bid], check.names=F)
x.ex<-data.frame(Treatment=x$treatment, ADRN=x$adrn.score, MES=x$mes.score, row.names=rownames(x))
x.cl<-factor(setNames(x$treatment, x$col), levels=AM[, unique(treatment)])
x.cl<-list(Treatment=setNames(unique(names(x.cl)), levels(x.cl)),
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
rm(x, x.ex, x.cl, x.hc, ph)

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/adrenergic-mesenchymal_score.RData




#############
#
#
#  scrap code
#
#
#############




#  ARID1A chr1:26729651-26732792 hsa_circ_0008494
arid1a<-'ENSG00000117713.20|ARID1A_chr1+26729651-26732792'



#  [lin_vs_circ] check mRNA and circRNA expression of ARID1A
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)


#  load all 
#  focus on circARID1A KD samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_201911/circRNAs_linear_vs_circular_collected_results.RData')
arid1a<-'ENSG00000117713.20|ARID1A_chr1+26729651-26732792'
arid1a<-do.call(rbind, lapply(lin.cir, function(x){ x[ circ_name %in% arid1a ] }))[ grep('CB-IMR5-circARID1A_', bid) ]
rm(lin.cir)


#  Mann-Whitney test
wilcox.test(x=arid1a[ grep('si6', bid), c.count ], y=arid1a[ grep('scr6', bid), c.count ], alternative='less')               #  p-value = 0.05
wilcox.test(x=arid1a[ grep('si6', bid), l.count.out ], y=arid1a[ grep('scr6', bid), l.count.out ], alternative='two.sided')  #  p-value = 0.7
wilcox.test(x=arid1a[ grep('si6', bid), l.count.max ], y=arid1a[ grep('scr6', bid), l.count.max ], alternative='two.sided')  #  p-value = 1.0
wilcox.test(x=arid1a[ grep('si6', bid), l.count ], y=arid1a[ grep('scr6', bid), l.count ], alternative='two.sided')          #  p-value = 0.7


#  t.test
t.test(x=arid1a[ grep('si6', bid), c.count ], y=arid1a[ grep('scr6', bid), c.count ], alternative='less')               #  p-value = 0.01
t.test(x=arid1a[ grep('si6', bid), l.count.out ], y=arid1a[ grep('scr6', bid), l.count.out ], alternative='two.sided')  #  p-value = 0.6
t.test(x=arid1a[ grep('si6', bid), l.count.max ], y=arid1a[ grep('scr6', bid), l.count.max ], alternative='two.sided')  #  p-value = 0.9
t.test(x=arid1a[ grep('si6', bid), l.count ], y=arid1a[ grep('scr6', bid), l.count ], alternative='two.sided')          #  p-value = 0.5

#}}}





