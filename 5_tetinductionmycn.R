###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################



#  BATCH_202004
#{{{


#################################################################################
#
#
#  define the metadata
#  create a preprocessed/ folder with cleaned up raw reads of at least 51nts long
#  create sample subfolders with symbolic links to raw FASTQ files
#  create Makefile and run the pipeline
#  add FASTQC and featureCount metrics to metadata
#
#{{{
rm(list=ls())
library(RColorBrewer)
library(data.table)


#  load metadata and process
meta<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/info/metadata.csv', header=F, sep=',', select=c(1:2, 4,6), col.names=c('bid', 'ids', 'library_prep', 'treatment'))[, cell_model:='MYCN']
meta[ grep('Tet', bid), treatment:='+Tet 120h']
meta[ grep('ETOH', bid), treatment:='ETOH 120h']


#  identify the sequencing ids from the Projektbericht PDF:
#
#      S2105Nr1 CD007
#      S2105Nr2 CD008
#      S2105Nr3 CD009
#      S2105Nr4 CD010
#      S2105Nr5 CD011
#      S2105Nr6 CD012
#  
#  add the sequencing ids to the metadata
stopifnot( all.equal( meta$ids, sprintf('CD0%02.0f', 7:12) ) )
meta[, sid:=paste0('S2105Nr', 1:6)]
setcolorder(meta, c('bid', 'sid', 'ids', 'library_prep', 'cell_model', 'treatment'))
meta[ treatment %in% '+Tet 120h', col:='orangered3']
meta[ treatment %in% 'ETOH 120h', col:='seagreen4']


#  provisional saving without FASTQC metrics and MultiQC metrics
#save(meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/metadata.RData')


#  locate the first mates
r1<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/fastq -maxdepth 1 -type f -wholename "*\\1\\.fastq.gz"', stdout=T)


#  identify the sequencing ids
names(r1)<-sub('(S2105Nr[0-9]*).*$', '\\1', basename(r1))


#  create the subfolders and place symbolic links to the two mates in them
root<-'/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/'
for(n in seq_along(r1)){
  d<-paste0(root, meta[ sid %in% names(r1)[n], bid])
  system2('mkdir', args=c('-p', d), stdout=T)
  system2('ln', args=c('-sf', r1[n], paste0(d, '/r1.fastq.gz')), stdout=T)
  system2('ln', args=c('-sf', sub('\\.1\\.fastq', '.2.fastq', r1[n]), paste0(d, '/r2.fastq.gz')), stdout=T)
}


#      copy the Makefile configuration from the parent directory and update all LIBS:= entries:
#
#          sed '/#  raw read libraries/q' ../totalrna.conf > totalrna.conf
#          echo -ne 'LIBS:=' >> totalrna.conf
#          find /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/CB* -maxdepth 0 -type d -printf " %p" >> totalrna.conf
#          ln -s ../../lib/qsub.mf .
#
#      run the pipeline at the cluster: 
#      
#          make -f qsub.mf CONF=totalrna.conf qc optitype kallisto rrna star bwa
#          make -f qsub.mf CONF=totalrna.conf counts ciri lin_vs_circ
#
#      run MultiQC:
#
#          multiqc --interactive -o multiqc --ignore '*dcc*' --ignore '*ciri*' --ignore '*hg19*' --ignore '*lin_vs_circ*' -v -f -d -s ./


#  load FASTQC metrics from MultiQC report
#  identify bid
#  summarize %GC and total number of raw reads across read mates
#  add to metadata
fgc<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/multiqc/multiqc_data/multiqc_fastqc.txt')[, c('Sample', 'Filename', '%GC', 'Total Sequences')]
colnames(fgc)<-c('sample', 'filename', 'gc', 'nreads')
fgc[, bid:=sub('^.*(C[HB][^ ]+).*$', '\\1', sample) ]
fgc<-fgc[, .(gc=mean(gc), nreads=mean(nreads)), by=.(bid)]
meta<-fgc[meta, on='bid']
rm(fgc)


#  process featureCounts statistics and add to metadata
logs<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/ -maxdepth 3 -wholename \'*/counts/genes.tsv.summary\' -print', stdout=T)
names(logs)<-sub('^.*BATCH_202004/(C[HB][^/]+)/counts/.*$', '\\1', logs)
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
setcolorder(meta, c('bid', 'sid', 'ids', 'library_prep', 'cell_model', 'treatment', 'failed', 'gc', 'nreads', 'features', 'no_features', 'unmapped', 'p_unmapped', 'col'))


#  save
save(meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/metadata.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/metadata.RData
#
#
#################################################################################




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
cnts<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/ -maxdepth 3 -type f -wholename \'*/counts/genes.tsv\' -print', stdout=T)


#  add ids
names(cnts)<-sub('^.*BATCH_202004/(C[BH][^/]+)/counts/.*$', '\\1', cnts)


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
save(totalrna, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/featureCounts.RData')


#  collect the exon counts including the exon junctions
cnts<-sub('genes\\.tsv$', 'exons.tsv', cnts)
exns<-List()
for (n in seq_along(cnts)){
  cat('\nprocessing: ', cnts[n], '\n')
  exns[[ names(cnts)[n] ]]<-parse_exon_counts(cnts[n])
}
rm(n)


#  save
save(exns, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/featureCounts_exons.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/featureCounts.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/featureCounts_exons.RData



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
cnts<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/ -maxdepth 3 -type f -wholename \'*/kallisto/abundance.tsv\' -print', stdout=T)


#  add ids
names(cnts)<-sub('^.*BATCH_202004/(C[HB][^/]+)/kallisto/.*$', '\\1', cnts)


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
save(totalrna, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/kallisto.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/kallisto.RData



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
cir<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/ -maxdepth 3 -type f -wholename \'*/ciri/circRNAs.tsv\' -print', stdout=T)


#  add ids
names(cir)<-sub('^.*BATCH_202004/(C[HB][^/]+)/ciri/.*$', '\\1', cir)


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
save(circ, trans, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/circRNAs_CIRI2.RData')
x<-data.frame(circ)[, -c(1:3,5)][, c(9, 1:8)]
x<-x[ order(x$jc_count, decreasing=T), ]
write.xlsx(x, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/circRNAs_CIRI2.xlsx', col.names=T, row.names=F, sheetName='circRNAs', append=F)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/circRNAs_CIRI2.{RData,xlsx}



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
lin_cir<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/ -maxdepth 3 -type f -wholename \'*/lin_vs_circ/counts.tsv\' -print', stdout=T)


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
save(lin.cir, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/circRNAs_linear_vs_circular_collected_results.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/circRNAs_linear_vs_circular_collected_results.RData



#  [STAR] collect mapping statistics for all samples
#{{{
rm(list=ls())
library(data.table)


#  locate STAR Log.final.out files 
logs<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/ -maxdepth 3 -wholename \'*/star/Log.final.out\' -print', stdout=T)


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
save(map.s, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/mapping_summary.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/mapping_summary.RData



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
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/mapping_summary.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/metadata.RData')
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
x11(width=18, height=16, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 


#  stacked bars
par(mar=c(9.0, 12.0, 2.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
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
mtext(text=sub('CB-SKNAS-TR-MYCN-', '', rownames(x)), side=1, line=0, at=bp, col=meta$col, las=2, adj=1, cex=1.8)
mtext('Number of raw reads', side=2, line=9, padj=+0.1, las=0, cex=2.4)
legend(x=par('usr')[2]*0.50, y=1.05*par('usr')[4], legend=cl.l$symbol[match(colnames(x), cl.l$symbol)] , col=cl.l$color[match(colnames(x), cl.l$symbol)], bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.40, x.intersp=0.15, seg.len=0.5)
legend(x=par('usr')[1]*0.50, y=1.05*par('usr')[4], legend=cl.s , col=names(cl.s),  bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.40, x.intersp=0.15, seg.len=0.5)
#mtext(paste0('Number of samples = ', nrow(x)), side=3, line=-1, padj=-0.6, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/barplot_mapping_statistics_alignments.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
options(scipen=0)


#  clean up
dev.off()

#}}}



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
hla<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/ -maxdepth 4 -type f -wholename \'*/optitype/*/*result.tsv\' -print', stdout=T)


#  add ids
names(hla)<-sub('^.*/(CB[^/]+)/optitype/.*$', '\\1', hla)


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
save(hla, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/hla_typing.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/hla_typing.RData



#  [HLA typing] by eye informatics comparing with previous SKNAS samples
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

scatterplots<-function(s1='', s2='', figs.dir=NULL, trim.name.root=''){
  #  scatterplot based on counts and TPMs of commonly expressed genes between sample (s1) and (s2)
  
  
  S1<-s1
  S2<-s2
  if(trim.name.root!=''){
    S1<-sub(trim.name.root, '', S1)
    S2<-sub(trim.name.root, '', S2)
  }
  
  
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
  mtext(bquote(log[10](1+counts) ~ group("[", .(S1), "]")), side=1, line=2, padj=+0.8, cex=1.8)
  mtext(bquote(log[10](1+counts) ~ group("[", .(S2), "]")), side=2, line=3, padj=+0.2, cex=1.8, las=0)
  if (!is.null(figs.dir)){
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
  mtext(bquote(log[10](1+TPM) ~ group("[", .(S1), "]")), side=1, line=2, padj=+0.8, cex=1.8)
  mtext(bquote(log[10](1+TPM) ~ group("[", .(S2), "]")), side=2, line=3, padj=+0.2, cex=1.8, las=0)
  if (!is.null(figs.dir)){
    dev.print(device=svg, file=paste0(figs.dir, '/scatterplot_tpms_', s1, '_', s2, '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
  }
}

#}}}


#  load collected HLA typing results and keep best solution
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/hla_typing.RData')
h<-do.call(rbind, lapply(hla, function(x){ x[1, ] }))[, 1:6]
h$bid<-rownames(h)
h<-data.table(h)[ order(A1, A2, B1, B2, C1, C2) ]
h[, hla:=paste(A1, A2, B1, B2, C1, C2, sep='_')]


#  load sample metadata
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/metadata.RData')


#  load count data
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/featureCounts.RData')


#  indicate the identical HLA types
h[hla %in% names(table(hla)[ table(hla)>1]), ]


#  scatterplot of identical HLA types
scatterplots('CB-SKNAS-TR-MYCN-120h-ETOH3', 'CB-SKNAS-TR-MYCN-120h-Tet1', figs.dir='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/', trim.name.root='CB-SKNAS-TR-MYCN-120h-')
meta[ bid %in% c('CB2019-11-R01', 'CB3037-11-R01') ]


#  check HLA types with the other SKNAS samples
h.this<-h
m.this<-meta
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/hla_typing.RData')
h<-do.call(rbind, lapply(hla, function(x){ x[1, ] }))[, 1:6]
h$bid<-rownames(h)
h<-data.table(h)[ order(A1, A2, B1, B2, C1, C2) ]
h[, hla:=paste(A1, A2, B1, B2, C1, C2, sep='_')]
h.prev<-h[ grep('CB-SKNAS-TR', bid) ]
rm(h, hla)

#}}}



#  [neuroblastoma-specific genes expression] MYCN, PHOX2B, ALK, BIRC5, CCND1, NTRK1, NRAS
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

genes_barplot<-function(G=c(''), names.trim='', plot.dir='', order.genes=T){
  source('~/bio/lib/trim_text.R')
  
  
  #  identify gene_id for given gene_name
  gid<-data.frame(mcols(hsa[ hsa$gene_name %in% G])[, c('gene_id', 'gene_name')])
  
  
  #  isolate genes of interest 
  #  add gene_name 
  #  sort by gene_name and bid
  #  compute log10(1+TPM)
  x<-fe[ gene_id %in% hsa$gene_id[ hsa$gene_name %in% G ] ]
  x$gene_name<-gid$gene_name[ match(x$gene_id, gid$gene_id) ]
  x<-x[, c('bid', 'gene_name', 'tpm')][ order(gene_name, bid) ]
  x$tpm<-log10(1+x$tpm)
  
  
  #  trim patient bids
  #x[, bid:=sub('-11-R01$', '', bid)]
  
  
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
  
  
  if(order.genes){
    x<-x[, order(colMeans(x), decreasing=T), drop=F]
  }
  
  
  #  define a color for each gene
  cl<-setNames(colorRampPalette(brewer.pal(8,'Dark2'))(ncol(x)), colnames(x))
  
  
  if(names.trim!=''){
    original.names<-rownames(x)
    rownames(x)<-sub(names.trim, '', rownames(x))
  }
  
  
  #  barplot 
  YTICK<-pretty(c(0, max(x)), 5)
  YMAX<-tail(YTICK, 1)
  par(mar=c(2.5, 8.0, 4.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
  bp<-barplot(t(x), beside=T, plot=F)
  plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(bp)+c(-1,2), axes=F, ann=F, xaxs='i', yaxs='i')
  bp<-barplot(t(x), beside=T, border='white', col=cl[colnames(x)], axisnames=F, xlab='', ylab='', las=1, yaxt='n', ylim=range(YTICK), add=T)
  axis(2, at=YTICK, line=0, cex.axis=2.4)
  mtext(text=trim_text(rownames(x), 12), side=1, line=1, at=colMeans(bp), las=0, adj=0.5, cex=1.8)
  mtext(expression(log[10](1+TPM)), side=2, line=5, padj=+0.1, las=0, cex=2.4)
  legend(x=par('usr')[2]*0.05, y=par('usr')[4]*1.15, legend=names(cl) , col=cl, bty='n', lty=1, lwd=18, cex=2.0, y.intersp=0.60, x.intersp=0.1, seg.len=0.1)
  
  
  #  save only if directory is given
  if(nchar(plot.dir)>0){
    dev.print(device=svg, file=paste0(plot.dir, '/barplot_', paste0(G, collapse='_'), '.svg'), width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
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
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/featureCounts.RData')
fe<-unlist(totalrna)
fe$bid<-sub('\\.[0-9]*$', '', rownames(fe))
rownames(fe)<-NULL
fe<-data.table(fe[, c('bid', 'gene_id', 'length', 'counts')])
fe<-fe[, tpm:=1e6*counts/length/sum(counts/length), by=.(bid)]


#  stacked barplot
x11(width=20, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 
GENES<-c('MYCN', 'PHOX2B', 'ALK', 'BIRC5', 'CCND1', 'NTRK1', 'NRAS')
s<-genes_barplot(GENES, names.trim='CB-SKNAS-TR-MYCN-120h-', plot.dir='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures')

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
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/metadata.RData')


#  from our cohort of circRNAs and identify their names
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
CIRCS.circ_name<-CIRCS$circ_name
rm(list=setdiff(ls(), c(l, 'CIRCS.circ_name')))


#  load all of circRNAs
#  keep only our cohort of circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/circRNAs_CIRI2.RData')
circs<-circ[ circ$region %in% 'exon' & circ$circ_name %in% CIRCS.circ_name ]
rm(circ, trans)


#  load featureCounts 
#  discard chrM, chrY genes
#  compute TPMs, CPMs based on the filtered list of genes
load('/data/sequencing/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/featureCounts.RData')
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
dds<-DESeqDataSetFromMatrix(countData=ceiling(gns.cnt), colData=data.frame(treatment=factor(meta$treatment), row.names=meta$bid), design=~1)
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
dds<-DESeqDataSetFromMatrix(countData=crs.cnt, colData=data.frame(treatment=factor(meta$treatment), row.names=meta$bid), design=~1)
sizeFactors(dds)<-gns.sf
crs.vs<-assay(varianceStabilizingTransformation(dds, fitType='local'))
crs.pca<-prcomp(t(crs.vs), center=T, scale.=F)
crs.ve<-round(1000 * crs.pca$sdev^2/sum(crs.pca$sdev^2))/10 
rm(dds)


#  save all including the reference for convenience
save(list=c('hsa', ls(pattern='(gns|crs)\\.'), 'meta', 'circs'), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_circRNAs+genes_all_libraries.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_circRNAs+genes_all_libraries.RData



#  clustering 
#  DESeq2 analysis
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
  plotMA(RES, xlab='', ylab='', main='', ylim=YLIM, xlim=XLIM, cex=1.2)
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
  identify(RES$baseMean, RES$log2FoldChange, labels=RES$gene_name, cex=0.7, offset=0.2, xpd=T)
  dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/plotMA_', Type, '_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
  rm(YLIM,XLIM)
  dev.off()
  
  
  #  mean-sd plots to see if variance strongly depends on mean
  X11(width=12, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
  par(mar=c(5,4,0.1,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.4, cex.axis=1.4)
  meanSdPlot(VSC)
  dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/meanSD_', Type, '_', sub(' vs ', '_', CND), '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
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
  dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/heatmap_', Type, '_vsc_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
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
  dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/heatmap_', Type, '_cor_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
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
  dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/plotPCA_', Type, '_vsc_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
  rm(XLIM,YLIM)
  dev.off()
  
  
  #  save DESeq results along with gene counts
  options(scipen=0)
  save(DDS,CND,RES,VSC,PCA,VE, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.RData'), compress=T)
  
  
  #  convert to table 
  x<-RES
  x$gene_id<-rownames(x)
  x<-x[, c('gene_id', 'gene_name', 'baseMean', 'log2FoldChange', 'pvalue', 'padj')]
  rownames(x)<-NULL
  x<-as.data.frame(x)
  write.table(x, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.tsv'), quote=F, sep='\t', row.names=F, col.names=T) 
  write.xlsx(x, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.xlsx'), row.names=F, col.names=T, sheetName=gsub(' ', '_', CND)) 
  rm(x)
  
  
  #readline("Hit ENTER to close all figures: ") 
  #for (n in setdiff(dev.list(), DEV_START)){ dev.off(n) }
  
  return(list(dds=DDS, res=RES, cnd=CND))
}

#}}}


#  load pre-prepared counts etc. for all samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_circRNAs+genes_all_libraries.RData')


#  [circRNAs + genes] clustering and PCA
#{{{

#  [circRNAs] heatmap based on Euclidean distance of variance stabilized (and glog2-transformed) counts 
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-crs.vs[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/heatmap_circRNAs_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [genes] heatmap based on Euclidean distance of variance stabilized (and glog2-transformed) counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-gns.vs[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/heatmap_genes_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [circRNAs] heatmaps based on Spearman correlations of raw counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-crs.cnt[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/heatmap_circRNAs_cor_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [genes] heatmaps based on Spearman correlations of raw counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-gns.cnt[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/heatmap_genes_cor_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x,m,d,ex,hc,ph,cl)
dev.off()


#  [circRNAs] PCA based on variance stabilized (and log2-transformed) counts centered but not scaled
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-crs.vs[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/plotPCA_circRNAs_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  [genes] PCA based on variance stabilized glog2-transformed counts centered but not scaled
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-gns.vs[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/plotPCA_genes_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,p,v,n,cl,XLIM,YLIM)
dev.off()

#}}}


#  [circRNAs] DESeq2 forcing library sizes from the gene counts 
m<-meta
M<-crs.cnt[, m$bid ]
SF<-gns.sf[ colnames(M) ]
colnames(M)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(M))
names(SF)<-sub('CB-SKNAS-TR-MYCN-', '', names(SF))
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=colnames(M))
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=ceiling(M), colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, 'ETOH 120h')  #  define first level so that the comparison log2(level 2/level 1) is done
sizeFactors(DDS)<-SF[ colnames(DDS) ]
d<-do_deseq2(DDS=DDS, CT=CT, Type='circRNAs', condition='Treatment', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)


#  [genes] DESeq2 
m<-meta
M<-gns.cnt[, m$bid ]
colnames(M)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(M))
names(SF)<-sub('CB-SKNAS-TR-MYCN-', '', names(SF))
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=colnames(M))
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=ceiling(M), colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, 'ETOH 120h')  #  define first level so that the comparison log2(level 2/level 1) is done
d<-do_deseq2(DDS=DDS, CT=CT, Type='genes', condition='Treatment', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_circRNAs_+Tet 120h_ETOH 120h.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_genes_+Tet 120h_ETOH 120h.RData



#  process of the DESeq2 results into HTML-table-ready objects
#  do MSigDB C2, C3 enrichment analysis
#  do MSigDB GSEA analysis
#  do GO MF, BP analysis
#  do mesenchymal/adrenergic marker enrichment test
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
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_circRNAs+genes_all_libraries.RData')
circs$bid<-sub('CB-SKNAS-TR-MYCN-', '', circs$bid)


#  simplify metadata bids
meta$bid<-sub('CB-SKNAS-TR-MYCN-', '', meta$bid)


#  [circRNAs] load DE results
#             separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_circRNAs_+Tet 120h_ETOH 120h.RData')
all.circ<-RES
up.circ<-subset(RES, log2FoldChange>0 & padj<0.1)    #  FDR cutoff: 0.1
down.circ<-subset(RES, log2FoldChange<0 & padj<0.1)  #  FDR cutoff: 0.1
bid.circ<-colnames(DDS)
meta.circ<-meta[ match( bid.circ, bid ), ]
rm(CND,DDS,PCA,RES,VE,VSC)


#  [genes] load DE results
#          separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_genes_+Tet 120h_ETOH 120h.RData')
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
save(list=ls(pattern='up\\.|down\\.|all\\.|bid\\.|circ|crs|gns|gsea|uni'), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_kable-ready_genes+circRNAs_+Tet 120h_ETOH 120h.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_kable-ready_genes+circRNAs_+Tet 120h_ETOH 120h.RData



#  gene set enrichment plots 
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
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/DESeq2_kable-ready_genes+circRNAs_+Tet 120h_ETOH 120h.RData')


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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/dotplot_genes_GOBP_upDE_+Tet 120h_ETOH 120h.svg', width=40, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/dotplot_genes_GOBP_downDE_+Tet 120h_ETOH 120h.svg', width=40, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/dotplot_genes_GOBP_neuron+differentiation_upDE_+Tet 120h_ETOH 120h.svg', width=40, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/dotplot_genes_GOBP_neuron+differentiation_downDE_+Tet 120h_ETOH 120h.svg', width=40, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}



#  cirular/(1+external linear) ratios
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(robustbase)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load circular and linear junction counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/circRNAs_linear_vs_circular_collected_results.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/metadata.RData')
m<-intersect(names(lin.cir), meta[ grepl('CB-SKNAS-TR-MYCN', bid) , bid ])
meta.tet<-meta[ bid %in% m ]
meta.tet<-meta.tet[, time.point:=factor(sub('^.*([0-9]+h)$', '\\1', treatment), levels=c('4h', '48h'))][ order(time.point), ][, time.point:=NULL]
lc.tet<-lin.cir[ meta.tet$bid ]
lc.tet<-unlist(lc.tet)
rm(m, lin.cir, meta)


#  convert all NA to zero
lc.tet[is.na(c.count), c.count:=0]
lc.tet[is.na(l.count.out), l.count.out:=0]
lc.tet[is.na(l.count.out.max), l.count.out.max:=0]


#  compute ratio of circular/(1+external linear) counts
#
#  N.B. ignore the warning about "invalid .internal.selfref detected"
#
lc.tet[, ratio:=c.count/(1+l.count.out)]


#  add risk_group metadata
#  add treatment metadata
lc.tet$treatment<-meta.tet$treatment[ match( lc.tet$bid , meta.tet$bid ) ]


#  load circular vs external linear junction correlation results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/cor_circular_vs_linear_junctions.RData')
length(R.high<-cl.cor[['all']][ cl.cor[['all']]>=0.6 ])  #  145
length(R.low<-cl.cor[['all']][ cl.cor[['all']]<0.6 ])    #  2053
lc.tet[, cor.class:='NA']
lc.tet[ gene_name %in% names(R.high), cor.class:='correlated']
lc.tet[ gene_name %in% names(R.low), cor.class:='uncorrelated']
lc.tet[, cor.class:=factor(cor.class, levels=c('correlated', 'uncorrelated', 'NA'))]
lc.tet[, table(cor.class)]
# cor.class
#   correlated uncorrelated           NA 
#         3186        25530           54 
#
lc.tet[, cor.class.col:='NA']
lc.tet[ cor.class %in% 'correlated', cor.class.col:='#b21f1f']
lc.tet[ cor.class %in% 'uncorrelated', cor.class.col:='#cccccc']
lc.tet[, cor.class.col:=factor(cor.class.col, levels=c('#b21f1f', '#cccccc', 'NA'))]
rm(cl.cor, R.high, R.low, ci, li, ci.cpm, li.cpm)


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  [MYCN Tet-inducible system]
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tet$ratio), factor(lc.tet$treatment, levels=unique(lc.tet$treatment)))
B.cl<-unique(meta.tet[, c('treatment', 'col')])
B.cl<-setNames(B.cl$col, B.cl$treatment)
B.cl['+Tet 120h']<-'orangered3'  #  manually change 
B.cl['ETOH 120h']<-'seagreen4'  #  manually change
wilcox.test(x=B[['+Tet 120h']], y=B[['ETOH 120h']], alternative='less')$p.value  #  1


#  boxplots
par(mar=c(13.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(meta.tet[, table(treatment)][names(B)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1.00, cex=2.4, col=B.cl)


#  [120h] ecdfs 
b<-B[c('ETOH 120h', '+Tet 120h')]
b.cl<-B.cl[c('ETOH 120h', '+Tet 120h')]
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(b[['ETOH 120h']])(x), from=min(YTICK), to=max(YTICK), n=min(length(b[['ETOH 120h']]), 20), ylab='', xlab='', pch=NA, col=b.cl['ETOH 120h'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(b[['+Tet 120h']])(x), from=min(YTICK), to=max(YTICK), n=min(length(b[['+Tet 120h']]), 20), pch=NA, col=b.cl['+Tet 120h'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
#legend('topleft', legend=paste0(names(b), ' (', as.numeric(meta.tet[, table(treatment)][names(b)]), ')'), col=b.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.90, x.intersp=0.2, seg.len=0.5)
legend('topleft', legend=names(b), col=b.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.90, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_+Tet_ETOH_120h.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b, b.cl)

#}}}


#  [MYCN Tet-inducible system] ratios of the correlated class
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tet$ratio), factor(paste(lc.tet$treatment, lc.tet$cor.class), levels=unique(paste(lc.tet$treatment, lc.tet$cor.class))))
B<-B[grep('NA$', names(B), invert=T)]  #  remove NA correlation class from circRNA/mRNA dropout pairs
B.cl<-unique(meta.tet[, c('treatment', 'col')])
B.cl<-setNames(B.cl$col, B.cl$treatment)
B.cl['+Tet 120h']<-'orangered3'  #  manually change 
B.cl['ETOH 120h']<-'seagreen4'  #  manually change
wilcox.test(x=B[['+Tet 120h correlated']], y=B[['ETOH 120h correlated']], alternative='greater')$p.value    #  0.01134584623


#  [120h] ecdfs 
b<-B[c('+Tet 120h correlated', 'ETOH 120h correlated')]
b.cl<-B.cl[c('+Tet 120h', 'ETOH 120h')]
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(b[[1]])(x), from=min(YTICK), to=max(YTICK), n=min(length(b[[1]]), 20), ylab='', xlab='', pch=NA, col=b.cl[1], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(b[[2]])(x), from=min(YTICK), to=max(YTICK), n=min(length(b[[2]]), 20), pch=NA, col=b.cl[2], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=names(b.cl), col=b.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.90, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202004/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_+Tet_ETOH_120h_correlated_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b, b.cl)

#}}}

#}}}

#}}}




##############################################################
#
#
#  BATCH_202006 below this point:
#
#
##############################################################




#  define the metadata
#  add FASTQC and featureCount metrics to metadata
#{{{
rm(list=ls())
library(RColorBrewer)
library(data.table)


#  load metadata and process
meta<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/P409_2a_Library_List_Illumina AG Schulte_27052020.tsv', header=F, sep='\t', col.names=c('bid', 'sid'))[, cell_model:='MYCN']
meta[ grep('Tet', bid), treatment:='+Tet 120h']
meta[ grep('ETOH', bid), treatment:='ETOH 120h']
meta[ treatment %in% '+Tet 120h', col:='orangered3']
meta[ treatment %in% 'ETOH 120h', col:='seagreen4']


#  provisional saving without FASTQC metrics and MultiQC metrics
#save(meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/metadata.RData')



#      run now MultiQC:
#
#          multiqc --interactive -o multiqc --ignore '*dcc*' --ignore '*ciri*' --ignore '*hg19*' --ignore '*lin_vs_circ*' -v -f -d -s ./



#  load FASTQC metrics from MultiQC report
#  identify bid
#  summarize %GC and total number of raw reads across read mates
#  add to metadata
fgc<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/multiqc/multiqc_data/multiqc_fastqc.txt')[, c('Sample', 'Filename', '%GC', 'Total Sequences')]
colnames(fgc)<-c('sample', 'filename', 'gc', 'nreads')
fgc[, bid:=sub('^.*(C[HB][^ ]+).*$', '\\1', sample) ]
fgc<-fgc[, .(gc=mean(gc), nreads=mean(nreads)), by=.(bid)]
meta<-fgc[meta, on='bid']
rm(fgc)


#  process featureCounts statistics and add to metadata
logs<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/ -maxdepth 3 -wholename \'*/counts/genes.tsv.summary\' -print', stdout=T)
names(logs)<-sub('^.*BATCH_202006/(C[HB][^/]+)/counts/.*$', '\\1', logs)
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
setcolorder(meta, c('bid', 'sid', 'cell_model', 'treatment', 'failed', 'gc', 'nreads', 'features', 'no_features', 'unmapped', 'p_unmapped', 'col'))


#  save
save(meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/metadata.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/metadata.RData




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
cnts<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/ -maxdepth 3 -type f -wholename \'*/counts/genes.tsv\' -print', stdout=T)


#  add ids
names(cnts)<-sub('^.*BATCH_202006/(C[BH][^/]+)/counts/.*$', '\\1', cnts)


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
save(totalrna, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/featureCounts.RData')


#  collect the exon counts including the exon junctions
cnts<-sub('genes\\.tsv$', 'exons.tsv', cnts)
exns<-List()
for (n in seq_along(cnts)){
  cat('\nprocessing: ', cnts[n], '\n')
  exns[[ names(cnts)[n] ]]<-parse_exon_counts(cnts[n])
}
rm(n)


#  save
save(exns, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/featureCounts_exons.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/featureCounts.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/featureCounts_exons.RData



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
cnts<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/ -maxdepth 3 -type f -wholename \'*/kallisto/abundance.tsv\' -print', stdout=T)


#  add ids
names(cnts)<-sub('^.*BATCH_202006/(C[HB][^/]+)/kallisto/.*$', '\\1', cnts)


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
save(totalrna, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/kallisto.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/kallisto.RData



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
cir<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/ -maxdepth 3 -type f -wholename \'*/ciri/circRNAs.tsv\' -print', stdout=T)


#  add ids
names(cir)<-sub('^.*BATCH_202006/(C[HB][^/]+)/ciri/.*$', '\\1', cir)


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
save(circ, trans, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/circRNAs_CIRI2.RData')
x<-data.frame(circ)[, -c(1:3,5)][, c(9, 1:8)]
x<-x[ order(x$jc_count, decreasing=T), ]
write.xlsx(x, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/circRNAs_CIRI2.xlsx', col.names=T, row.names=F, sheetName='circRNAs', append=F)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/circRNAs_CIRI2.{RData,xlsx}



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
lin_cir<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/ -maxdepth 3 -type f -wholename \'*/lin_vs_circ/counts.tsv\' -print', stdout=T)


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
save(lin.cir, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/circRNAs_linear_vs_circular_collected_results.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/circRNAs_linear_vs_circular_collected_results.RData



#  [STAR] collect mapping statistics for all samples
#{{{
rm(list=ls())
library(data.table)


#  locate STAR Log.final.out files 
logs<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/ -maxdepth 3 -wholename \'*/star/Log.final.out\' -print', stdout=T)


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
save(map.s, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/mapping_summary.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/mapping_summary.RData



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
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/mapping_summary.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/metadata.RData')
stopifnot( length(setdiff( map.s$bid , meta$bid ))==0 )  #  all sequenced libraries show up in the metadata
meta<-meta[ bid %in% map.s$bid ]
map.s<-map.s[ match(meta$bid, map.s$bid), ]


#  colors related to alignments and samples
cl.l<-data.frame(symbol=c('unimapped', 'multimapped','unmapped'), 
                 color=c('#1874CD',  #  dodgerblue3
                         '#EEB422',  #  goldenrod2
                         '#999999')) #  grey
cl.s<-setNames( unique(meta[, treatment]), unique(meta[, col]) )


#  CLICK on it once to make sure it does not redraw, or options(scipen=-20) might be IGNORED
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 


#  stacked bars
par(mar=c(12.0, 14.0, 2.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.8, cex.axis=2.8)
x<-map.s[, c('unimapped', 'multimapped', 'unmapped')]
rownames(x)<-sub('CB-SKNAS-TR-MYCN-', '', sub('-R01-R2$', '', map.s$bid))
YMAX<-max( rowSums(x) - rowSums(x)%%1e5 + 1e5 )
YTICK<-pretty(c(0, YMAX), 5)
YMAX<-tail(YTICK, 1)
bp<-barplot(t(x), plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(bp)+c(-1, +1), axes=F, ann=F, xaxs='i', yaxs='i') 
bp<-barplot(t(x), border='white', col=cl.l[match(colnames(x), cl.l$symbol), 2], axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, YMAX), add=T)
options(scipen=-20)
axis(2, at=YTICK, line=-1, cex.axis=2.8)
mtext(text=sub('CB-SKNAS-TR-MYCN-', '', rownames(x)), side=1, line=0, at=bp, col=meta$col, las=2, adj=1, cex=2.4)
mtext('Number of reads', side=2, line=10, padj=-0.2, las=0, cex=2.8)
legend(x=par('usr')[2]*0.50, y=1.14*par('usr')[4], legend=cl.l$symbol[match(colnames(x), cl.l$symbol)] , col=cl.l$color[match(colnames(x), cl.l$symbol)], bty='n', lty=1, lwd=18, cex=2.0, y.intersp=0.40, x.intersp=0.45, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/barplot_mapping_statistics_alignments.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
options(scipen=0)


#  clean up
dev.off()

#}}}



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
hla<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/ -maxdepth 4 -type f -wholename \'*/optitype/*/*result.tsv\' -print', stdout=T)


#  add bid
names(hla)<-sub('^.*/(CB[^/]+)/optitype/.*$', '\\1', hla)


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
save(hla, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/hla_typing.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/hla_typing.RData



#  [HLA typing] by eye informatics comparing with previous SKNAS samples
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

scatterplots<-function(s1='', s2='', figs.dir=NULL, trim.name.root=''){
  #  scatterplot based on counts and TPMs of commonly expressed genes between sample (s1) and (s2)
  
  
  S1<-s1
  S2<-s2
  if(trim.name.root!=''){
    S1<-sub(trim.name.root, '', S1)
    S2<-sub(trim.name.root, '', S2)
  }
  
  
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
  mtext(bquote(log[10](1+counts) ~ group("[", .(S1), "]")), side=1, line=2, padj=+0.8, cex=1.8)
  mtext(bquote(log[10](1+counts) ~ group("[", .(S2), "]")), side=2, line=3, padj=+0.2, cex=1.8, las=0)
  if (!is.null(figs.dir)){
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
  mtext(bquote(log[10](1+TPM) ~ group("[", .(S1), "]")), side=1, line=2, padj=+0.8, cex=1.8)
  mtext(bquote(log[10](1+TPM) ~ group("[", .(S2), "]")), side=2, line=3, padj=+0.2, cex=1.8, las=0)
  if (!is.null(figs.dir)){
    dev.print(device=svg, file=paste0(figs.dir, '/scatterplot_tpms_', s1, '_', s2, '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
  }
}

#}}}


#  load collected HLA typing results and keep best solution
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/hla_typing.RData')
h<-do.call(rbind, lapply(hla, function(x){ x[1, ] }))[, 1:6]
h$bid<-rownames(h)
h<-data.table(h)[ order(A1, A2, B1, B2, C1, C2) ]
h[, hla:=paste(A1, A2, B1, B2, C1, C2, sep='_')]


#  load sample metadata
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/metadata.RData')


#  load count data
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/featureCounts.RData')


#  indicate the identical HLA types
h[hla %in% names(table(hla)[ table(hla)>1]), ]


#  scatterplot of identical HLA types
#scatterplots('CB-SKNAS-TR-MYCN-120h-ETOH3', 'CB-SKNAS-TR-MYCN-120h-Tet1', figs.dir='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/', trim.name.root='CB-SKNAS-TR-MYCN-120h-')
#meta[ bid %in% c('CB2019-11-R01', 'CB3037-11-R01') ]

#}}}



#  [neuroblastoma-specific genes expression] MYCN, PHOX2B, ALK, BIRC5, CCND1, NTRK1, NRAS
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

genes_barplot<-function(G=c(''), names.trim='', plot.dir='', order.genes=T){
  source('~/bio/lib/trim_text.R')
  
  
  #  identify gene_id for given gene_name
  gid<-data.frame(mcols(hsa[ hsa$gene_name %in% G])[, c('gene_id', 'gene_name')])
  
  
  #  isolate genes of interest 
  #  add gene_name 
  #  sort by gene_name and bid
  #  compute log10(1+TPM)
  x<-fe[ gene_id %in% hsa$gene_id[ hsa$gene_name %in% G ] ]
  x$gene_name<-gid$gene_name[ match(x$gene_id, gid$gene_id) ]
  x<-x[, c('bid', 'gene_name', 'tpm')][ order(gene_name, bid) ]
  x$tpm<-log10(1+x$tpm)
  
  
  #  trim patient bids
  #x[, bid:=sub('-11-R01$', '', bid)]
  
  
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
  
  
  if(order.genes){
    x<-x[, order(colMeans(x), decreasing=T), drop=F]
  }
  
  
  #  define a color for each gene
  cl<-setNames(colorRampPalette(brewer.pal(8,'Dark2'))(ncol(x)), colnames(x))
  
  
  if(names.trim!=''){
    original.names<-rownames(x)
    rownames(x)<-sub(names.trim, '', rownames(x))
  }
  
  
  #  barplot 
  YTICK<-pretty(c(0, max(x)), 5)
  YMAX<-tail(YTICK, 1)
  par(mar=c(2.5, 8.0, 4.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
  bp<-barplot(t(x), beside=T, plot=F)
  plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(bp)+c(-1,2), axes=F, ann=F, xaxs='i', yaxs='i')
  bp<-barplot(t(x), beside=T, border='white', col=cl[colnames(x)], axisnames=F, xlab='', ylab='', las=1, yaxt='n', ylim=range(YTICK), add=T)
  axis(2, at=YTICK, line=0, cex.axis=2.4)
  mtext(text=trim_text(rownames(x), 12), side=1, line=1, at=colMeans(bp), las=0, adj=0.5, cex=1.8)
  mtext(expression(log[10](1+TPM)), side=2, line=5, padj=+0.1, las=0, cex=2.4)
  legend(x=par('usr')[2]*0.04, y=par('usr')[4]*1.15, legend=names(cl) , col=cl, bty='n', lty=1, lwd=18, cex=2.0, y.intersp=0.60, x.intersp=0.1, seg.len=0.1)
  
  
  #  save only if directory is given
  if(nchar(plot.dir)>0){
    dev.print(device=svg, file=paste0(plot.dir, '/barplot_', paste0(G, collapse='_'), '.svg'), width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
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
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/featureCounts.RData')
fe<-unlist(totalrna)
fe$bid<-sub('\\.[0-9]*$', '', rownames(fe))
rownames(fe)<-NULL
fe<-data.table(fe[, c('bid', 'gene_id', 'length', 'counts')])
fe<-fe[, tpm:=1e6*counts/length/sum(counts/length), by=.(bid)]


#  trim suffix from bids
fe[, bid:=sub('-R01-R2$', '', bid)]


#  stacked barplot
x11(width=20, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 
GENES<-c('MYCN', 'PHOX2B', 'ALK', 'BIRC5', 'CCND1', 'NTRK1', 'NRAS')
s<-genes_barplot(GENES, names.trim='CB-SKNAS-TR-MYCN-120h-', plot.dir='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures')

#}}}



#  [stratification markers] calculate the proliferative index based on TPMs which removes biases on gene length
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


#  load metadata
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/metadata.RData')
cel<-meta
rm(meta)


#  load featureCounts 
#  remove chrM/chrY counts
#  compute TPMs
#  convert to matrix
#  log2-transform to regularize the range
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/featureCounts.RData')
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


#  keep TPMs of the marker genes
#  order according to metadata
cel.tpm<-gns.tpm[rownames(gns.tpm) %in% ProliferativeIndex:::metaPCNA2, cel$bid]
rm(gns.tpm)


#  add proliferative index to metadata
cel$PI<-apply(cel.tpm, 2, median)


#  cluster samples by PI
#  cut at maximum silhouette width 
#  add discrete PI classification to the metadata
x<-as.matrix(setNames(cel$PI, cel$bid), check.names=F)
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
cel$PI.class<-factor(w[ cel$bid ], levels=c('low', 'high'))
cel$PI.class.col<-factor(w[ cel$bid ], levels=c('low', 'high'), labels=c('#cccccc', '#b21f1f'))



#  save
save(cel, cel.tpm, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/proliferative_index.RData')


#  visualize the classification results
#{{{

#  load back the results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/proliferative_index.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/proliferative_index.RData') #KH


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  boxplot of PI values by treatment
par(mar=c(7.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(cel$PI, cel$treatment)
#
wilcox.test(x=B[['+Tet 120h']], y=B[['ETOH 120h']], alternative='greater')$p.value  #  0.5
#
B.cl<-setNames(unique(cel[, col]), unique(cel[, treatment]))
YTICK<-pretty(c(4, max(sapply(B, max))), 4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('PI', side=2, line=4, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=0, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/boxplot_PI_per_treatment_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

cel$treatment_renamed = ifelse(
  cel$treatment == "ETOH 120h", 
  "Off", 
  ifelse(cel$treatment == "+Tet 120h", "On", NA))
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
boxplot(PI ~ treatment_renamed, data = cel, col = c("#22693c", "#86398b"), ylim = c(4,5.5), medcol='cyan', range = 0, ylab='', xlab='', frame.plot=F, xpd=F)#, show.names=F, frame.plot=F, boxwex=0.8, xpd=F, outline=F, color=B.cl, range=0, add=T)
vioplot(PI ~ treatment_renamed, data = cel, col = c("#22693c", "#86398b"), ylim = c(4,5.5))

stripchart(PI ~ as.factor(treatment_renamed), data = cel, vertical = T, ylim = c(4,5.5), col = c("#22693c", "#86398b"))

dev.off()
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
boxplot(PI ~ treatment_renamed, data = cel, col = NA, ylim = c(4,5.5), medcol=NA, range = 0, xlab='', frame.plot=F, xpd=F, border = NA)#, show.names=F, frame.plot=F, boxwex=0.8, xpd=F, outline=F, color=B.cl, range=0, add=T)
points(x = as.factor(cel$treatment_renamed), y = cel$PI, col = ifelse(cel$treatment_renamed == "Off", "#22693c", "#86398b"), pch = 19)
mtext('PI', side=2, line=4, padj=-0.1, las=0, cex=2.4)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/pointplot_PI_per_treatment_group_KH.pdf', width=8, height=8, bg='white', pointsize=20)

dev.off()
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
vioplot(PI ~ treatment_renamed, data = cel, col = c("#22693c", "#86398b"), ylim = c(4,5.5), medcol='cyan', range = 0, xlab='', frame.plot=F, xpd=F)#, show.names=F, frame.plot=F, boxwex=0.8, xpd=F, outline=F, color=B.cl, range=0, add=T)
dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/vioplot_PI_per_treatment_group_KH.pdf', width=8, height=8, bg='white', pointsize=20)



#  boxplot of PI values by PI class
par(mar=c(7.5, 6.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(cel$PI, cel$PI.class)
B.cl<-setNames(levels(cel$PI.class.col), levels(cel$PI.class))
YTICK<-pretty(range(sapply(B, max)), 4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('PI', side=2, line=4, padj=-0.1, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=0, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/boxplot_PI_per_proliferative_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  heatmap of annotated samples clustered by Euclidean distance in log2-TPMs
x<-cel.tpm
stopifnot( all.equal( colnames(x), cel$bid ) )
x.ex<-data.frame(Treatment=cel$treatment, PI=cel$PI.class, row.names=cel$bid)
x.cl<-setNames(list(setNames( unique(cel$col), unique(x.ex$Treatment)), setNames(levels(cel$PI.class.col), levels(cel$PI.class))), colnames(x.ex))
x.d<-dist(t(x), method='euclidean')
x.hc<-hclust(x.d, method='ward.D2')
dev.off()
svg('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/heatmap_proliferative_index.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=18, family='Arial')  #  workaround the cutting of legends?
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
rm(x, x.ex, x.cl, x.hc, ph)

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/proliferative_index.RData




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
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/metadata.RData')


#  from our cohort of circRNAs identify their names
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
CIRCS.circ_name<-CIRCS$circ_name
rm(list=setdiff(ls(), c(l, 'CIRCS.circ_name')))


#  load all of circRNAs
#  keep only our cohort of circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/circRNAs_CIRI2.RData')
circs<-circ[ circ$region %in% 'exon' & circ$circ_name %in% CIRCS.circ_name ]
rm(circ, trans)


#  load featureCounts 
#  discard chrM, chrY genes
#  compute TPMs, CPMs based on the filtered list of genes
load('/data/sequencing/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/featureCounts.RData')
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
dds<-DESeqDataSetFromMatrix(countData=ceiling(gns.cnt), colData=data.frame(treatment=factor(meta$treatment), row.names=meta$bid), design=~1)
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
dds<-DESeqDataSetFromMatrix(countData=crs.cnt, colData=data.frame(treatment=factor(meta$treatment), row.names=meta$bid), design=~1)
sizeFactors(dds)<-gns.sf
crs.vs<-assay(varianceStabilizingTransformation(dds, fitType='local'))
crs.pca<-prcomp(t(crs.vs), center=T, scale.=F)
crs.ve<-round(1000 * crs.pca$sdev^2/sum(crs.pca$sdev^2))/10 
rm(dds)


#  save all including the reference for convenience
save(list=c('hsa', ls(pattern='(gns|crs)\\.'), 'meta', 'circs'), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_circRNAs+genes_all_libraries.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_circRNAs+genes_all_libraries.RData



#  clustering 
#  DESeq2 analysis
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
  dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/plotMA_', Type, '_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
  rm(YLIM,XLIM)
  dev.off()
  
  
  #  mean-sd plots to see if variance strongly depends on mean
  X11(width=12, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
  par(mar=c(5,4,0.1,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=1.4, cex.axis=1.4)
  meanSdPlot(VSC)
  dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/meanSD_', Type, '_', sub(' vs ', '_', CND), '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
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
  dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/heatmap_', Type, '_vsc_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
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
  dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/heatmap_', Type, '_cor_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
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
  dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/plotPCA_', Type, '_vsc_', sub(' vs ', '_', CND), '.svg'), width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
  rm(XLIM,YLIM)
  dev.off()
  
  
  #  save DESeq results along with gene counts
  save(DDS,CND,RES,VSC,PCA,VE, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.RData'), compress=T)
  
  
  #  convert to table 
  x<-RES
  x$gene_id<-rownames(x)
  x<-x[, c('gene_id', 'gene_name', 'baseMean', 'log2FoldChange', 'pvalue', 'padj')]
  rownames(x)<-NULL
  x<-as.data.frame(x)
  write.table(x, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.tsv'), quote=F, sep='\t', row.names=F, col.names=T) 
  write.xlsx(x, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_', Type, '_', sub(' vs ', '_', CND), '.xlsx'), row.names=F, col.names=T, sheetName=gsub(' ', '_', CND)) 
  rm(x)
  
  
  #readline("Hit ENTER to close all figures: ") 
  #for (n in setdiff(dev.list(), DEV_START)){ dev.off(n) }
  
  return(list(dds=DDS, res=RES, cnd=CND))
}

#}}}


#  load pre-prepared counts etc. for all samples
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_circRNAs+genes_all_libraries.RData')


#  [circRNAs + genes] clustering and PCA
#{{{

#  [circRNAs] heatmap based on Euclidean distance of variance stabilized (and glog2-transformed) counts 
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-crs.vs[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/heatmap_circRNAs_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [genes] heatmap based on Euclidean distance of variance stabilized (and glog2-transformed) counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-gns.vs[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/heatmap_genes_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [circRNAs] heatmaps based on Spearman correlations of raw counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-crs.cnt[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/heatmap_circRNAs_cor_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,d,ex,hc,ph,cl)
dev.off()


#  [genes] heatmaps based on Spearman correlations of raw counts
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-gns.cnt[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/heatmap_genes_cor_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x,m,d,ex,hc,ph,cl)
dev.off()


#  [circRNAs] PCA based on variance stabilized (and log2-transformed) counts centered but not scaled
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-crs.vs[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/plotPCA_circRNAs_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  [genes] PCA based on variance stabilized glog2-transformed counts centered but not scaled
x11(width=14, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
m<-copy(meta)
x<-gns.vs[, m$bid]
colnames(x)<-sub('CB-SKNAS-TR-MYCN-', '', colnames(x))
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/plotPCA_genes_vsc_all.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(m,x,p,v,n,cl,XLIM,YLIM)
dev.off()

#}}}


#  [circRNAs] DESeq2 forcing library sizes from the gene counts 
m<-meta
M<-crs.cnt[, m$bid ]
SF<-gns.sf[ colnames(M) ]
colnames(M)<-sub('-R01-R2$', '', sub('CB-SKNAS-TR-MYCN-', '', colnames(M)))
names(SF)<-sub('-R01-R2$', '', sub('CB-SKNAS-TR-MYCN-', '', names(SF)))
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=colnames(M))
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=ceiling(M), colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, 'ETOH 120h')  #  define first level so that the comparison log2(level 2/level 1) is done
sizeFactors(DDS)<-SF[ colnames(DDS) ]
d<-do_deseq2(DDS=DDS, CT=CT, Type='circRNAs', condition='Treatment', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)


#  [genes] DESeq2 
m<-meta
M<-gns.cnt[, m$bid ]
colnames(M)<-sub('-R01-R2$', '', sub('CB-SKNAS-TR-MYCN-', '', colnames(M)))
names(SF)<-sub('-R01-R2$', '', sub('CB-SKNAS-TR-MYCN-', '', names(SF)))
CT<-data.frame(m[, c('cell_model', 'treatment', 'col'), with=F], row.names=colnames(M))
colnames(CT)<-c('cell_model', 'Treatment', 'col')
DDS<-DESeqDataSetFromMatrix(countData=ceiling(M), colData=data.frame(Treatment=factor(CT$Treatment), row.names=rownames(CT)), design=~Treatment)
DDS$Treatment<-relevel(DDS$Treatment, 'ETOH 120h')  #  define first level so that the comparison log2(level 2/level 1) is done
d<-do_deseq2(DDS=DDS, CT=CT, Type='genes', condition='Treatment', lfcThreshold=0.0, lfcShrinkType='normal', include.annot=F, names.trim='', list(par=list(mar=c(5.0, 8.0, 1.0, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.4, cex.axis=2.4, tcl=-0.2), ylab.line=4))
rm(m,M,CT,DDS)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_circRNAs_+Tet 120h_ETOH 120h.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_genes_+Tet 120h_ETOH 120h.RData



#  process of the DESeq2 results into HTML-table-ready objects
#  do MSigDB C2, C3 enrichment analysis
#  do MSigDB GSEA analysis
#  do GO MF, BP analysis
#  do mesenchymal/adrenergic marker enrichment test
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
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_circRNAs+genes_all_libraries.RData')
circs$bid<-sub('CB-SKNAS-TR-MYCN-', '', circs$bid)


#  simplify metadata bids
meta$bid<-sub('-R01-R2$', '', sub('CB-SKNAS-TR-MYCN-', '', meta$bid))


#  [circRNAs] load DE results
#             separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_circRNAs_+Tet 120h_ETOH 120h.RData')
all.circ<-RES
up.circ<-subset(RES, log2FoldChange>0 & padj<0.1)    #  FDR cutoff: 0.1
down.circ<-subset(RES, log2FoldChange<0 & padj<0.1)  #  FDR cutoff: 0.1
bid.circ<-colnames(DDS)
meta.circ<-meta[ match( bid.circ, bid ), ]
rm(CND,DDS,PCA,RES,VE,VSC)


#  [genes] load DE results
#          separate significantly up-regulated and down-regulated
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_genes_+Tet 120h_ETOH 120h.RData')
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
save(list=ls(pattern='up\\.|down\\.|all\\.|bid\\.|circ|crs|gns|gsea|uni'), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_kable-ready_genes+circRNAs_+Tet 120h_ETOH 120h.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_kable-ready_genes+circRNAs_+Tet 120h_ETOH 120h.RData



#  gene set enrichment plots 
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
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_kable-ready_genes+circRNAs_+Tet 120h_ETOH 120h.RData')


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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/dotplot_genes_GOBP_upDE_+Tet 120h_ETOH 120h.svg', width=40, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [GO BP] downregulated genes (p-value cutoff<5e-4)
#{{{

x<-down.bp[p.value<5e-4, c('Term', 'Significant', 'p.value')][, GeneRatio:=Significant/length(sigGenes(down.bp.topgo))][ order(GeneRatio) ][, Term:=factor(Term, levels=Term)]
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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/dotplot_genes_GOBP_downDE_+Tet 120h_ETOH 120h.svg', width=20, height=30, bg='white', antialias='subpixel', pointsize=20, family='Arial')

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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/dotplot_genes_GOBP_neuron+differentiation_upDE_+Tet 120h_ETOH 120h.svg', width=25, height=10, bg='white', antialias='subpixel', pointsize=20, family='Arial')

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
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/dotplot_genes_GOBP_neuron+differentiation_downDE_+Tet 120h_ETOH 120h.svg', width=20, height=10, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}



#  enrichment in MYCN targets
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(DESeq2)


#  load genes DE results
#  remove subversions
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_genes_+Tet 120h_ETOH 120h.RData')
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
#  => 2.095029649e-34


#  induced targets in the significantly down-DE genes
grp<-rownames( subset(gns, log2FoldChange<0 & padj<0.05) )
intersect(grp, IND)  #  ENSG00000166482, ENSG00000178035,ENSG00000172403


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
#  => 0.002880596799


#  repressed targets in the significantly up-DE genes
grp<-rownames( subset(gns, log2FoldChange>0 & padj<0.05) )
intersect(grp, REP)  #  ENSG00000078018, ENSG00000068366, ENSG00000132872

#}}}



#  cirular/(1+external linear) ratios
#  compare with MNA vs HR_nMNA results as well
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(robustbase)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/my_ecdfs.R')


#  [tumors] load ratios and keep MNA and HR_nMNA only
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_all_classifications_together.RData')
lc.hr<-lc.tumors[ risk_group %in% c('HR_nMNA', 'MNA') ]
rm(list=setdiff(ls(), c(l, 'lc.hr')))


#  load circular and linear junction counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/circRNAs_linear_vs_circular_collected_results.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/metadata.RData')
m<-intersect(names(lin.cir), meta[ grepl('CB-SKNAS-TR-MYCN', bid) , bid ])
meta.tet<-meta[ bid %in% m ]
meta.tet<-meta.tet[, time.point:=factor(sub('^.*([0-9]+h)$', '\\1', treatment), levels=c('4h', '48h'))][ order(time.point), ][, time.point:=NULL]
lc.tet<-lin.cir[ meta.tet$bid ]
lc.tet<-unlist(lc.tet)
rm(m, lin.cir, meta)


#  convert all NA to zero
lc.tet[is.na(c.count), c.count:=0]
lc.tet[is.na(l.count.out), l.count.out:=0]
lc.tet[is.na(l.count.out.max), l.count.out.max:=0]


#  compute ratio of circular/(1+external linear) counts
#
#  N.B. ignore the warning about "invalid .internal.selfref detected"
#
lc.tet[, ratio:=c.count/(1+l.count.out)]


#  add risk_group metadata
#  add treatment metadata
lc.tet$treatment<-meta.tet$treatment[ match( lc.tet$bid , meta.tet$bid ) ]


#  load circular vs external linear junction correlation results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/cor_circular_vs_linear_junctions.RData')
length(R.high<-cl.cor[['all']][ cl.cor[['all']]>=0.6 ])  #  139
length(R.low<-cl.cor[['all']][ cl.cor[['all']]<0.6 ])    #  2073
lc.tet[, cor.class:='NA']
lc.tet[ gene_name %in% names(R.high), cor.class:='correlated']
lc.tet[ gene_name %in% names(R.low), cor.class:='uncorrelated']
lc.tet[, cor.class:=factor(cor.class, levels=c('correlated', 'uncorrelated', 'NA'))]
lc.tet[, table(cor.class)]
# cor.class
#   correlated uncorrelated           NA 
#         2982        25914           54 
#
lc.tet[, cor.class.col:='NA']
lc.tet[ cor.class %in% 'correlated', cor.class.col:='#b21f1f']
lc.tet[ cor.class %in% 'uncorrelated', cor.class.col:='#cccccc']
lc.tet[, cor.class.col:=factor(cor.class.col, levels=c('#b21f1f', '#cccccc', 'NA'))]
rm(R.high, R.low, cl.cor, ci, li, ci.cpm, li.cpm)


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tet$ratio), factor(lc.tet$treatment, levels=unique(lc.tet$treatment)))
B.cl<-unique(meta.tet[, c('treatment', 'col')])
B.cl<-setNames(B.cl$col, B.cl$treatment)
wilcox.test(x=B[['+Tet 120h']], y=B[['ETOH 120h']], alternative='less')$p.value  #  1.171e-07


#  [120h] ecdfs 
b<-B[c('ETOH 120h', '+Tet 120h')]
b.cl<-B.cl[c('ETOH 120h', '+Tet 120h')]
my_ecdfs(b, b.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='topleft', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_+Tet_ETOH_120h.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
rm(b, b.cl)


#  ratios separately for the correlated/uncorrelated classes and together
#{{{

#  split groups by treatment and correlation class and compute log2(0.001 + ratio) and test for significant differences
B<-split(log2(1e-3+lc.tet$ratio), factor(paste(lc.tet$treatment, lc.tet$cor.class), levels=unique(paste(lc.tet$treatment, lc.tet$cor.class))))
B<-B[grep('NA$', names(B), invert=T)]  #  remove NA correlation class associated with circRNA/mRNA dropout pairs
B.cl<-unique(meta.tet[, c('treatment', 'col')])
B.cl<-setNames(B.cl$col, B.cl$treatment)
B.cl['+Tet 120h']<-'mediumslateblue'  #  manually change 
B.cl['ETOH 120h']<-'mediumpurple4'  #  manually change
B.cl[names(B)]<-colorRampPalette(brewer.pal(8, 'Paired'))(length(B))  #  add the uncorrelated/correlated as color pairs
wilcox.test(x=B[['+Tet 120h correlated']], y=B[['ETOH 120h correlated']], alternative='less')$p.value      #  0.006223402289
wilcox.test(x=B[['+Tet 120h uncorrelated']], y=B[['ETOH 120h uncorrelated']], alternative='less')$p.value  #  2.391706414e-06
wilcox.test(x=B[['+Tet 120h uncorrelated']], y=B[['+Tet 120h correlated']], alternative='less')$p.value    #  0.6994682988
wilcox.test(x=B[['ETOH 120h uncorrelated']], y=B[['ETOH 120h correlated']], alternative='less')$p.value    #  0.06432215358


#  [120h] ecdfs of correlated class only
b<-B[c('+Tet 120h correlated', 'ETOH 120h correlated')]
b.cl<-B.cl[names(b)]
my_ecdfs(b, b.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='bottomright', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_+Tet_ETOH_120h_correlated_group.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
rm(b, b.cl)


#  [120h] ecdfs of uncorrelated class only
b<-B[c('+Tet 120h uncorrelated', 'ETOH 120h uncorrelated')]
b.cl<-B.cl[names(b)]
my_ecdfs(b, b.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='bottomright', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_+Tet_ETOH_120h_uncorrelated_group.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
rm(b, b.cl)


#  [120h] ecdfs of both classes plotted together
b<-B[c('+Tet 120h correlated', 'ETOH 120h correlated', '+Tet 120h uncorrelated', 'ETOH 120h uncorrelated')]
b.cl<-B.cl[ names(b) ]
my_ecdfs(b, b.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='bottomright', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_+Tet_ETOH_120h_per_cor_group.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
rm(b, b.cl)

#}}}


# comparison with MNA vs HR_nMNA
#{{{


#  compute mean ratios across replicates/tumors and the difference between conditions
#
#  N.B. The Tet system has many more circRNA dropouts. This can be because it is a cell line and it also has much smaller cohort size (3+3)
#
hr<-lc.hr[, .(ratio=mean(ratio)), by=.(risk_group, circ_name)]
tet<-lc.tet[ treatment %in% c('+Tet 120h', 'ETOH 120h') ][, .(ratio=mean(ratio)), by=.(treatment, circ_name)]
tet<-tet[, .(ratio=list(ratio), diff=ratio[ treatment %in% '+Tet 120h']-ratio[ treatment %in% 'ETOH 120h']), by=.(circ_name)][ order(diff) ]
hr<-hr[, .(ratio=list(ratio), diff=ratio[ risk_group %in% 'MNA']-ratio[ risk_group %in% 'HR_nMNA']), by=.(circ_name)][order(diff)]


#  find commonly expressed circRNAs to be fair to the Tet system
length(common<-intersect(hr[ sapply(ratio, function(r){ any(r!=0) }), circ_name ], tet[ sapply(ratio, function(r){ any(r!=0) }), circ_name ]))  #  4253
hr<-hr[ circ_name %in% common  ]
tet<-tet[ circ_name %in% common ]
hr[, table(diff<0)]
# 
# FALSE  TRUE 
#   732  3521 
tet[, table(diff<0)]
# 
# FALSE  TRUE 
#  2061  2192 
#
wilcox.test(x=hr$diff, y=tet$diff, alternative='two.sided')$p.value   #  7.435e-281
#
x<-hr[ tet, on='circ_name' ]       #  i.diff belongs to tet
round(100*x[ diff<0 & i.diff<0, .N]/x[ diff<0, .N], 1)    #  51.5% of the downshifted in tumors were also downshifted in MYCN Tet-inducible system
round(100*x[ diff<0 & i.diff<0, .N]/x[ i.diff<0, .N], 1)  #  82.8% of the downshifted in MYCN Tet-inducible system were downshifted in the tumors


#  histogram
BREAKS<-pretty(range(c(hr$diff, tet$diff)), 200)
XLIM<-c(-1.0, 0.2)
h<-hist(hr$diff, breaks=BREAKS, border='white', col='seagreen', xlim=XLIM, freq=T, ylim=c(0, 2.5e3), xlab='', ylab='', main='')
hist(tet$diff, breaks=h$breaks, border='white', col=adjustcolor('lightsalmon', alpha.f=0.5), freq=T, add=T)


#  ecdfs of ratio difference
B<-list('MNA - HR_nMNA'=hr$diff, 'MYCN Tet - ETOH'=tet$diff)
B.cl<-setNames(c('darkseagreen', 'darksalmon'), names(B))
my_ecdfs(B, B.cl, XLIM=c(-0.3, 0.2), XLAB=expression(paste(Delta,'ratio')), XLINE=3, YLINE=4, MAIN='', LTY=1, LWD=12, LEGEND='topleft', 
         ABLINE=list(v=setNames(0.0, 'black'), h=setNames(c(ecdf(B[[1]])(0), ecdf(B[[2]])(0)), as.character(B.cl))), 
         svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_differences_circular_vs_linear_junctions_MNA_HR_nMNA+MYCN_Tet-inducible_+Tet_ETOH_120h.svg', 
         mar=c(4.5, 7.5, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  [downregulated in MYCN Tet enrichment in downregulated in MNA] Fisher's exact test:
#
#                    |     MYCN Tet down    |   MYCN Tet up   |
#  ------------------|----------------------|-----------------|
#         MNA down   |           x1         |       y1        |
#  ------------------|----------------------|-----------------|
#         MNA up     |           x2         |       y2        | 
#  ------------------|----------------------|-----------------|
x1<-length(intersect(tet[diff<0, circ_name], hr[ diff<0, circ_name]))
x2<-length(intersect(tet[diff<0, circ_name], hr[ diff>0, circ_name]))
y1<-length(intersect(tet[diff>0, circ_name], hr[ diff<0, circ_name]))
y2<-length(intersect(tet[diff>0, circ_name], hr[ diff>0, circ_name]))
fisher.test(data.frame('in'=c(x1, x2), 'out'=c(y1, y2)), alternative='greater')$p.value
# 
# => 0.5396 



#  [top N downregulated in MYCN Tet enrichment in top N downregulated in MNA] Fisher's exact test:
#
#                        |     MYCN Tet down top N   |   MYCN Tet not down top N  |
#  ----------------------|---------------------------|----------------------------|
#      MNA down top N    |              x1           |              y1            |
#  ----------------------|---------------------------|----------------------------|
#    MNA not down top N  |              x2           |              y2            | 
#  ----------------------|---------------------------|----------------------------|
N<-500
a1<-tet[diff<0][order(diff)][1:N, circ_name]  #  order again just in case
b1<-hr[diff<0][order(diff)][1:N, circ_name]   #  order again just in case
x1<-length(intersect(a1, b1))
x2<-length(intersect(a1, hr[ ! circ_name %in% b1, circ_name]))
y1<-length(intersect(tet[ ! circ_name %in% a1, circ_name], b1))
y2<-length(intersect(tet[ ! circ_name %in% a1, circ_name], hr[ ! circ_name %in% b1, circ_name]))
fisher.test(data.frame('in'=c(x1, x2), 'out'=c(y1, y2)), alternative='greater')$p.value
# 
# => 9.218e-58 

#}}}

#}}}




#################
#
#
#  extra analysis
#
#
#################




#  [MiSeq, circARID1A KD, Tet-inducible MYCN older runs] explore count results 
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


#  featureCounts with mate 2 forward strand as it should
#{{{

#  load gene counts with -s 2 orientation
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202005MISEQ/featureCounts_s2.RData')
gns.mi<-unlist(totalrna)
gns.mi$bid<-sub('\\.[0-9]*$', '', rownames(gns.mi))
rownames(gns.mi)<-NULL
gns.mi<-data.table(gns.mi[, c('bid', 'gene_id', 'length', 'counts')])
gns.mi<-gns.mi[, .(gene_id=gene_id, counts=counts, tpm=1e6*counts/length/sum(counts/length), cpm=1e6*counts/sum(counts)), by=.(bid)]
gns.mi[, gene_name:=hsa$gene_name[ match(gene_id, hsa$gene_id) ] ]


#  create the TPM matrix
#  remove unexpressed genes
#  order by top expression
gns.mi.tpm<-dcast(gns.mi, bid ~ gene_name, value.var='tpm', fun.aggregate=sum)  #  TPM matrix
y<-as.data.frame(gns.mi.tpm[, -1])
rownames(y)<-gns.mi.tpm[, bid]
gns.mi.tpm<-t(y)
gns.mi.tpm<-gns.mi.tpm[ rowMeans(gns.mi.tpm)>0, , drop=F]
gns.mi.tpm<-gns.mi.tpm[ order(rowMeans(gns.mi.tpm), decreasing=T), ]


#  create the CPM matrix
#  remove unexpressed genes
#  order by top expression
gns.mi.cpm<-dcast(gns.mi, bid ~ gene_name, value.var='cpm', fun.aggregate=sum)
y<-as.data.frame(gns.mi.cpm[, -1])
rownames(y)<-gns.mi.cpm[, bid]
gns.mi.cpm<-t(y)
gns.mi.cpm<-gns.mi.cpm[ rowMeans(gns.mi.cpm)>0, , drop=F]
gns.mi.cpm<-gns.mi.cpm[ order(rowMeans(gns.mi.cpm), decreasing=T), ]

#}}}


#  featureCounts Tet-inducible MYCN
#{{{

#  remove failed samples and Pilot samples
#  keep Tet-inducible MYCN system only
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-meta.cel[ !(failed) & grepl('CB-SKNAS-TR-MYCN', bid), ]
rm(meta.tum, meta.cel, meta.prefailed)


#  load gene counts 
#  keep Tet-inducible MYCN system only
load('/data/sequencing/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')
gns.te<-unlist(totalrna[ meta$bid ])
gns.te$bid<-sub('\\.[0-9]*$', '', rownames(gns.te))
rownames(gns.te)<-NULL
gns.te<-data.table(gns.te[, c('bid', 'gene_id', 'length', 'counts')])
gns.te<-gns.te[, .(gene_id=gene_id, counts=counts, tpm=1e6*counts/length/sum(counts/length), cpm=1e6*counts/sum(counts)), by=.(bid)]
gns.te[, gene_name:=hsa$gene_name[ match(gene_id, hsa$gene_id) ] ]


#  create the TPM matrix
#  remove unexpressed genes
#  order by top expression
gns.te.tpm<-dcast(gns.te, bid ~ gene_name, value.var='tpm', fun.aggregate=sum)  #  TPM matrix
y<-as.data.frame(gns.te.tpm[, -1])
rownames(y)<-gns.te.tpm[, bid]
gns.te.tpm<-t(y)
gns.te.tpm<-gns.te.tpm[ rowMeans(gns.te.tpm)>0, , drop=F]
gns.te.tpm<-gns.te.tpm[ order(rowMeans(gns.te.tpm), decreasing=T), ]


#  create the CPM matrix
#  remove unexpressed genes
#  order by top expression
gns.te.cpm<-dcast(gns.te, bid ~ gene_name, value.var='cpm', fun.aggregate=sum)
y<-as.data.frame(gns.te.cpm[, -1])
rownames(y)<-gns.te.cpm[, bid]
gns.te.cpm<-t(y)
gns.te.cpm<-gns.te.cpm[ rowMeans(gns.te.cpm)>0, , drop=F]
gns.te.cpm<-gns.te.cpm[ order(rowMeans(gns.te.cpm), decreasing=T), ]

#}}}


#  define junk
JUNK<-'^RN7S|^MT-|^AC[0-9]+|^AL[0-9]+|^MTND[12]P|^MTCO[123]'


#  percentage of reads going to junk
p.mi<-round(gns.mi.cpm, digits=2)/1e6
j.mi<-p.mi[ grep(JUNK, rownames(p.mi)), ]
j.mi<-j.mi[ order(rowMeans(j.mi), decreasing=T), ,drop=F]
colSums(j.mi)
#        Off         On 
# 0.22843751 0.23764826 


#  [Tet-inducible MYCN] percentage of reads going to junk
p.te<-round(gns.te.cpm, digits=2)/1e6
j.te<-p.te[ grep(JUNK, rownames(p.te)), ]
j.te<-j.te[ order(rowMeans(j.te), decreasing=T), ,drop=F]
colSums(j.te)
# CB-SKNAS-TR-MYCN--Early-Tet1 CB-SKNAS-TR-MYCN--Early-Tet2 CB-SKNAS-TR-MYCN--Early-Tet3  CB-SKNAS-TR-MYCN--late-Tet1  CB-SKNAS-TR-MYCN--late-Tet2 
#                   0.17964199                   0.19452938                   0.23270691                   0.21640543                   0.19550153 
#  CB-SKNAS-TR-MYCN--late-Tet3 CB-SKNAS-TR-MYCN-Early-ETOH1 CB-SKNAS-TR-MYCN-Early-ETOH2 CB-SKNAS-TR-MYCN-Early-ETOH3  CB-SKNAS-TR-MYCN-late-ETOH1 
#                   0.20839741                   0.17325671                   0.20910210                   0.20936928                   0.17184582 
#  CB-SKNAS-TR-MYCN-late-ETOH2  CB-SKNAS-TR-MYCN-late-ETOH3 
#                   0.19306515                   0.19822788 


#  what is capturing the non-junk reads?
nrow(nj.mi<-p.mi[ setdiff(rownames(p.mi), rownames(j.mi)), ])  #  21019
nrow(nj.te<-p.te[ setdiff(rownames(p.te), rownames(j.te)), ])  #  27219
nj.mi<-nj.mi[ order(rowMeans(nj.mi), decreasing=T), ,drop=F]
nj.te<-nj.te[ order(rowMeans(nj.te), decreasing=T), ,drop=F]


#  top 100 in each case
length(i<-intersect(rownames(nj.mi)[1:100], rownames(nj.te)[1:100]))   #  74
d<-setdiff(rownames(nj.mi)[1:100], rownames(nj.te)[1:100])
#
#  [1] "NEAT1"    "H19"      "EEF1A1P6" "UBC"      "UNC5C"    "UBR4"     "ENAH"     "RGS5"     "SFPQ"     "CENPF"    "HSP90B1"  "USP9X"    "SLC38A1" 
# [14] "CBX5"     "TPR"      "EIF3A"    "DDX21"    "NAP1L1"   "MACF1"    "DDX3X"    "CANX"     "EIF2S3"   "LARP1"    "DDX17"    "HNRNPH1"  "SIPA1L2" 


#   what about MYCN?
gns.te.cpm[ 'MYCN', ]
# CB-SKNAS-TR-MYCN--Early-Tet1 CB-SKNAS-TR-MYCN--Early-Tet2 CB-SKNAS-TR-MYCN--Early-Tet3  CB-SKNAS-TR-MYCN--late-Tet1  CB-SKNAS-TR-MYCN--late-Tet2 
#                 459.30980875                 497.09427673                 490.00476123                 444.40901617                 362.11996405 
#  CB-SKNAS-TR-MYCN--late-Tet3 CB-SKNAS-TR-MYCN-Early-ETOH1 CB-SKNAS-TR-MYCN-Early-ETOH2 CB-SKNAS-TR-MYCN-Early-ETOH3  CB-SKNAS-TR-MYCN-late-ETOH1 
#                 355.39196084                  64.23427419                  57.35490731                  68.56394248                  55.20539237 
#  CB-SKNAS-TR-MYCN-late-ETOH2  CB-SKNAS-TR-MYCN-late-ETOH3 
#                  54.10908353                  53.85903202 
#
#
gns.mi.cpm[ 'MYCN', ]
#          Off           On 
#  32.27824984 240.89197434 


#  barplot of MYCN and top 20 common+not-common non-junk genes averaged over Tet-induction groups
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial') 
i<-intersect(rownames(nj.mi)[1:20], rownames(nj.te)[1:20])
d<-setdiff(rownames(nj.mi)[1:20], rownames(nj.te)[1:20])
g<-gns.te.cpm[c('MYCN', i, d), , drop=F]
colnames(g)<-sub('CB-SKNAS-TR-MYCN--late-Tet[0-9]', 'On_48h', colnames(g))
colnames(g)<-sub('CB-SKNAS-TR-MYCN--Early-Tet[0-9]', 'On_4h', colnames(g))
colnames(g)<-sub('CB-SKNAS-TR-MYCN-late-ETOH[0-9]', 'Off_48h', colnames(g))
colnames(g)<-sub('CB-SKNAS-TR-MYCN-Early-ETOH[0-9]', 'Off_4h', colnames(g))
g<-t(do.call(rbind, setNames(lapply(unique(colnames(g)), function(x){ rowMeans(g[, colnames(g) %in% x]) }), unique(colnames(g)))))
x<-gns.mi.cpm[c('MYCN', i, d), , drop=F]
colnames(x)<-sub('Off', 'Off_120h', colnames(x))
colnames(x)<-sub('On', 'On_120h', colnames(x))
stopifnot(all.equal( rownames(g), rownames(x) ) )
g<-cbind(g, x)
g<-g[, c('Off_4h', 'Off_48h', 'Off_120h', 'On_4h', 'On_48h', 'On_120h')]
g<-g[ order(rowMeans(g), decreasing=T), ]
B<-log10(t(g))
B.cl<-setNames(colorRampPalette(brewer.pal(8,'Accent'))(nrow(B)), rownames(B))
YTICK<-pretty(c(0, max(B)), 5)
par(mar=c(12.0, 6.5, 2.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
bp<-barplot(B, beside=T, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, tail(YTICK, 1)), xlim=range(bp)+c(-1, +1), axes=F, ann=F, xaxs='i', yaxs='i') 
bp<-barplot(B, border='white', col=B.cl, axisnames=F, beside=T, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK, 1)), add=T)
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext(expression(log[10]('CPM')), side=2, line=2, padj=-0.3, las=0, cex=2.4)
mtext(text=colnames(B), side=1, line=1, at=colMeans(bp), col='black', las=2, adj=0.90, cex=2.4)
legend(x=par('usr')[2]*0.7, y=par('usr')[4]*1.10, legend=rownames(B), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.50, x.intersp=0.4, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202005MISEQ/figures/barplot_top_non-junk_CPMs.svg', width=20, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}




######################################
#
#
#  scrap code/quick and dirty analyses
#
#
######################################




#  list of circRNAs that could be seen as differentially expressed 
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(DESeq2)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load circular and linear junction counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/circRNAs_linear_vs_circular_collected_results.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/metadata.RData')
m<-intersect(names(lin.cir), meta[ grepl('CB-SKNAS-TR-MYCN', bid) , bid ])
meta.tet<-meta[ bid %in% m ]
meta.tet<-meta.tet[, time.point:=factor(sub('^.*([0-9]+h)$', '\\1', treatment), levels=c('4h', '48h'))][ order(time.point), ][, time.point:=NULL]
lc.tet<-lin.cir[ meta.tet$bid ]
lc.tet<-unlist(lc.tet)
rm(m, lin.cir, meta)


#  convert all NA to zero
lc.tet[is.na(c.count), c.count:=0]
lc.tet[is.na(l.count.out), l.count.out:=0]
lc.tet[is.na(l.count.out.max), l.count.out.max:=0]


#  compute ratio of circular/(1+external linear) counts
#
#  N.B. ignore the warning about "invalid .internal.selfref detected"
#
lc.tet[, ratio:=c.count/(1+l.count.out)]


#  add treatment metadata
#  compute mean ratio across all treatments per circRNA
lc.tet$treatment<-meta.tet$treatment[ match( lc.tet$bid , meta.tet$bid ) ]
ratios<-lc.tet[, .(ratio=mean(ratio)), by=.(circ_name)]
rm(lc.tet, meta.tet)


#  load the processed DE results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/BATCH_202006/DESeq2_kable-ready_genes+circRNAs_+Tet 120h_ETOH 120h.RData')


#  compute number of isoforms and mean CPMs across both conditions per circRNA isoform
#  group and rank isoforms by mean CPM
x<-data.table(data.frame(mcols(circs)[, c('circ_name', 'gene_name', 'bid', 'cpm')]))[, .(cpm=mean(cpm), gene_name=unique(gene_name)), by=.(circ_name)][, .(circs=list(circ_name[ order(cpm, decreasing=T) ]), cpms=list(sort(cpm, decreasing=T)), N=.N), by=.(gene_name) ]


#  declare first isoform among multiplets dominant if its expression is at least double of the second in rank
x[ N==1, dominant:=T]
x[ N>1, dominant:=sapply(cpms, function(x){ ifelse(x[1]/x[2]>=2.0, T, F) })]


#  keep only the list of dominant circRNAs 
dominant<-x[ (dominant) ][, .(circ_name=circs[[1]][1], cpm=cpms[[1]][1]), by=.(gene_name)]
rm(x)


#  order circRNA DE results by baseMean and abs(log2FC)
x<-data.table(data.frame(all.circ))
x[, gene_id:=rownames(all.circ)]
x[, c('lfcSE', 'stat', 'padj'):=NULL]
x<-x[ order(-baseMean, -abs(log2FoldChange)) ]


#  separate the validated (include MAN1A2 too)
#  add estimated cirRNA/mRNA ratio
validated<-x[ gene_name %in% c('ARID1A', 'EYA1', 'SETD3', 'HIPK3', 'SMARCA5', 'ZNF124', 'NFATC3', 'ASAP1', 'TET2', 'MAML3', 'AKAP12', 'LINC00632', 'MAN1A2')]
validated<-dominant[ validated, on='gene_name']
validated[, c('cpm', 'gene_id', 'pvalue'):=NULL]
validated<-ratios[ validated, on='circ_name' ]


#  rest of the well-expressed with strong fold-change circRNAs and a clearly dominant isoform
rest<-x[ ! gene_name %in% validated[, gene_name] & gene_name %in% dominant[, gene_name] & baseMean>=20 & abs(log2FoldChange)>=log2(1.2) ]
rest<-dominant[ rest, on='gene_name']
rest[, c('cpm', 'gene_id', 'pvalue'):=NULL]
rest<-ratios[ rest, on='circ_name' ]
rm(x)


#  create an Excel sheet with these two lists
wb<-createWorkbook()
addWorksheet(wb, 'validated')
addWorksheet(wb, 'rest')
writeDataTable(wb, 'validated', data.frame(validated))
writeDataTable(wb, 'rest', data.frame(rest))
setColWidths(wb, sheet=1, cols=seq_len(ncol(validated)), widths='auto')
setColWidths(wb, sheet=2, cols=seq_len(ncol(rest)), widths='auto')
saveWorkbook(wb, file='~/Downloads/MYCN_+Tet 120h_candidate_DE_circRNAs.xlsx', overwrite=T)

#}}}



