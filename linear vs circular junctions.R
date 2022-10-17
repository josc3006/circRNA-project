###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




#########################################################################
#
#
#  linear and circular junctions construction and expression estimation
#
#
#########################################################################




#  [run once] construct circularizing junctions and linear junctions reference:
#             
#                 circularizing junctions:
#             
#                     We used the unified cohort of circRNAs identified in the neuroblastoma tumors.
#                     We allowed for +1nt discrepancies in exon boundaries.
#                     In cases of ambiguity the exon with the longest part upstream of the 3' and downstream of the 5' sites was chosen.
#                     We spliced the corresponding exons so that the circular junctions are formed and kept those of at least 160nts.
#             
#                 linear junctions:
#                     
#                     For the genes associated with the kept circularized junctions we formed all unique annotated linear junctions and kept
#                     them at their natural lengths.
#                     This was done because the backspliced junctions contain variable sized chunks of exons sometimes larger than 80nts, 
#                     in which case the aligner chooses backspliced junctions instead of linear spliced ones!
#                     If a gene ended up with no linear junctions then all circularized junctions for that gene were also dropped
#
#{{{
rm(list=ls())
library(data.table)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)


#  load circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')


#  extract all exons per gene
#  exclude chrM exons
txdb<-loadDb('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/TxDb_GRCh38.gencode.v30.RData')
EXONS<-unlist(exonsBy(txdb, 'gene'))
EXONS$gene_id<-names(EXONS)
names(EXONS)<-NULL
seqlevels(EXONS, pruning.mode='coarse')<-seqlevels(CIRCS)


#  [~10min] create the circular junctions
#{{{

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


#  keep only circular junctions of at least 160nts long
keep<-P.UP$w+P.DOWN$w>=160
cat('Number of plus strand circular junctions of at least 160nts = ', sum(keep), ' (number of dropouts = ', sum(!keep), ')\n', sep='')
P<-P[ keep ]
P.UP<-P.UP[ keep ]
P.DOWN<-P.DOWN[ keep ]
keep<-M.UP$w+M.DOWN$w>=160
cat('Number of minus strand circular junctions of at least 160nts = ', sum(keep), ' (number of dropouts = ', sum(!keep), ')\n', sep='')
M<-M[ keep ]
M.UP<-M.UP[ keep ]
M.DOWN<-M.DOWN[ keep ]
rm(keep)


#  trim exons to appropriate lengths so that an 160nts circular junction is created
x<-data.table(up=P.UP$w, down=P.DOWN$w)
stopifnot( nrow(x[ up<80 & down<80])==0 )  #  if one side is short of nts, the other should always compensate
x[ up<80, c('up.w', 'down.w'):=list(as.integer(up), as.integer(160-up)) ]
x[ down<80, c('up.w', 'down.w'):=list(as.integer(160-down), as.integer(down)) ]
x[ up>=80 & down>=80, c('up.w', 'down.w'):=list(as.integer(80), as.integer(80)) ]
P.UP<-resize(P.UP, x$up.w, fix='end')
mcols(P.UP)<-mcols(P.UP)[, 'gene_id', drop=F]
P.DOWN<-resize(P.DOWN, x$down.w, fix='start')
mcols(P.DOWN)<-mcols(P.DOWN)[, 'gene_id', drop=F]
stopifnot( all.equal(end(P.UP), end(P)) )
stopifnot( all.equal(start(P.DOWN), start(P)) )
x<-data.table(up=M.UP$w, down=M.DOWN$w)
stopifnot( nrow(x[ up<80 & down<80])==0 )  #  if one side is short of nts, the other should always compensate
x[ up<80, c('up.w', 'down.w'):=list(as.integer(up), as.integer(160-up)) ]
x[ down<80, c('up.w', 'down.w'):=list(as.integer(160-down), as.integer(down)) ]
x[ up>=80 & down>=80, c('up.w', 'down.w'):=list(as.integer(80), as.integer(80)) ]
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

#}}}


#  [<1min] create the linear junctions
#{{{

#  get all annotated transcripts
#  N.B. transcrips on the minus strand have the exon order properly inverted
LINEAR<-exonsBy(txdb, 'tx', use.names=T)


#  keep only transcripts of genes in the circRNA list
gids<-select(txdb, names(LINEAR), columns=c('GENEID'), keytype='TXNAME')
colnames(gids)<-c('txid', 'gene_id')
LINEAR<-LINEAR[ unique(gids$txid[ gids$gene_id %in% CIRCS$gene_id ]) ]
gids<-gids[ match(names(LINEAR), gids$txid), ]


#  form all pairwise exon junctions per transcript
x<-unlist(LINEAR)
x$txid<-names(x)
names(x)<-NULL
x<-data.table(txid=x$txid, seqnames=as.character(seqnames(x)), start=start(x), end=end(x), strand=as.character(strand(x)))[, exon:=paste(seqnames, strand, start, end, sep='_')]
x<-x[, .(junctions=list(paste(head(exon, -1), tail(exon, -1), sep='|'))), by=.(txid)]


#  take unique exon junctions per gene across all annotated isoforms and inflate to one junction entry per row
x[, gene_id:=gids$gene_id[ match(txid, gids$txid) ] ]
x<-x[, .(junctions=unique(unlist(junctions))), by=.(gene_id)]


#  identify the genomic ranges of each exon involved in each junction
x[, c('j1', 'j2'):=list(sapply(strsplit(junctions, '|', fixed=T), '[[', 1), sapply(strsplit(junctions, '|', fixed=T), '[[', 2))]
x[, junctions:=NULL]
x[, c('seqnames1', 'strand1', 'start1', 'end1'):=list(sapply(strsplit(j1, '_', fixed=T), '[[', 1), sapply(strsplit(j1, '_', fixed=T), '[[', 2), sapply(strsplit(j1, '_', fixed=T), '[[', 3), sapply(strsplit(j1, '_', fixed=T), '[[', 4))]
x[, j1:=NULL]
x[, c('seqnames2', 'strand2', 'start2', 'end2'):=list(sapply(strsplit(j2, '_', fixed=T), '[[', 1), sapply(strsplit(j2, '_', fixed=T), '[[', 2), sapply(strsplit(j2, '_', fixed=T), '[[', 3), sapply(strsplit(j2, '_', fixed=T), '[[', 4))]
x[, j2:=NULL]
x[, c('start1', 'end1', 'start2', 'end2'):=list(as.integer(start1), as.integer(end1), as.integer(start2), as.integer(end2))]
stopifnot( all.equal( x$seqnames1, x$seqnames2 ) )
stopifnot( all.equal( x$strand1, x$strand2 ) )
LINEAR<-x
UP<-GRanges(seqnames=LINEAR$seqnames1, strand=LINEAR$strand1, ranges=IRanges(start=LINEAR$start1, end=LINEAR$end1), DataFrame(gene_id=LINEAR$gene_id))
DOWN<-GRanges(seqnames=LINEAR$seqnames2, strand=LINEAR$strand2, ranges=IRanges(start=LINEAR$start2, end=LINEAR$end2), DataFrame(gene_id=LINEAR$gene_id))
rm(x)


###################################################################################################################################################
#
#
#  N.B. If you trim linear junctions to 160nts then there are backspliced junctions with larger chunks of exons as compared to the linear junctions
#       in which case the aligner maps to the backspliced junctions reads that belong to the linear junction!
#       Therefore, we keep the linear junctions intact, even those less than 160nts long and give the aligner a chance to correctly decide
#       what belongs to where.
#
#  keep only junctions that can be at least 160nts
#LINEAR[, c('w1', 'w2'):=list(as.integer(end1-start1+1), as.integer(end2-start2+1))]
#keep<-LINEAR[, w1+w2>=160 ]
#cat('Number of linear junctions of at least 160nts = ', sum(keep), ' (number of dropouts = ', sum(!keep), ')\n', sep='')
#LINEAR<-LINEAR[ keep, ]
#UP<-GRanges(seqnames=LINEAR$seqnames1, strand=LINEAR$strand1, ranges=IRanges(start=LINEAR$start1, end=LINEAR$end1), DataFrame(gene_id=LINEAR$gene_id))
#DOWN<-GRanges(seqnames=LINEAR$seqnames2, strand=LINEAR$strand2, ranges=IRanges(start=LINEAR$start2, end=LINEAR$end2), DataFrame(gene_id=LINEAR$gene_id))
#rm(keep)


#  trim exons to appropriate lengths so that 160nts junctions are created
#x<-data.table(up=width(UP), down=width(DOWN))
#stopifnot( nrow(x[ up<80 & down<80])==0 )  #  if one side is short of nts, the other should always compensate
#x[ up<80, c('up.w', 'down.w'):=list(as.integer(up), as.integer(160-up)) ]
#x[ down<80, c('up.w', 'down.w'):=list(as.integer(160-down), as.integer(down)) ]
#x[ up>=80 & down>=80, c('up.w', 'down.w'):=list(as.integer(80), as.integer(80)) ]
#UP<-resize(UP, x$up.w, fix='end')
#DOWN<-resize(DOWN, x$down.w, fix='start')
#rm(x)


#  now that we trimmed junctions to 160nts we need to check again that unique linear junctions per gene are used
#x<-data.table(up=paste(as.character(seqnames(UP)), as.character(strand(UP)), start(UP), end(UP), sep='_'), down=paste(as.character(seqnames(UP)), as.character(strand(UP)), start(UP), end(UP), sep='_'), gene_id=UP$gene_id, i=seq_along(UP))[, junction:=paste(up, down, sep='|')]
#x<-x[, .(i=i[!duplicated(junction)], junction=junction[ !duplicated(junction) ]), by=.(gene_id) ]
#UP<-UP[ x$i ]
#DOWN<-DOWN[ x$i ]
#LINEAR<-LINEAR[ x$i, ]
#rm(x)
#
#
###################################################################################################################################################


#  construct the linear junctions 
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
LINEAR.junctions<-split(x, x$n)  #  we checked that the order of the exons is preserved
LINEAR<-LINEAR[, .(junction_name=paste0(gene_id, '|jn_', seq_len(.N))), by=.(gene_id)]
stopifnot( nrow(LINEAR)==length(LINEAR.junctions) )
names(LINEAR.junctions)<-LINEAR$junction_name
LINEAR.junctions.seqs<-extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, LINEAR.junctions)
rm(x,UP,DOWN)

#}}}


#  remove circular junctions of genes with no linear junctions at all
dropouts<-setdiff( CIRCS$gene_id, LINEAR$gene_id )
cat('dropout circular junctions with no linear junctions:\n')
CIRCS[ CIRCS$gene_id %in% dropouts ]   #  none, CDR1as is now LINC00632 and it has linear junctions!
dropouts<-CIRCS$circ_name[ CIRCS$gene_id %in% dropouts ]
CIRCS<-CIRCS[ ! CIRCS$circ_name %in% dropouts ]
CIRCS.junctions<-CIRCS.junctions[ ! names(CIRCS.junctions) %in% dropouts ]
CIRCS.junctions.seqs<-CIRCS.junctions.seqs[ ! names(CIRCS.junctions.seqs) %in% dropouts ]
rm(dropouts)


#  annotate the linear junctions that lie inside the circRNA bounds for later distinction between internal and external junctions
#
#  N.B. BOTH junction exons need to be found within, otherwise the junction should not count. 
#  N.B. all single-exon backspliced junctions should not have any linear junction within.
#
x<-unlist(LINEAR.junctions)
x$jc_name<-names(x)
names(x)<-NULL
o<-findOverlaps(x, CIRCS, select='all', type='within')
o<-data.table(data.frame(o[ x[queryHits(o)]$gene_id == CIRCS[ subjectHits(o) ]$gene_id ]))  #  keep only those with the same gene_id
o[, jn:=x$jc_name[ queryHits ]]                                                             #  add junction names to check if they appear twice per hit
o<-o[, .(n=.N), by=.(subjectHits, jn)][n==2][, .(junctions=list(jn)), by=.(subjectHits)]    #  keep only junctions with both exons within 
o[, circ:=CIRCS$circ_name[ subjectHits ]]                                                   #  name the backspliced junction hits
o<-o[ data.table(circ=CIRCS$circ_name), on='circ']                                          #  add back the backspliced junctions that dropped out
stopifnot( all.equal( o$circ, CIRCS$circ_name ) )
CIRCS$linear_junctions_within<-o[, junctions]
rm(x,o)


#  save
save(CIRCS, CIRCS.junctions, CIRCS.junctions.seqs, LINEAR, LINEAR.junctions, LINEAR.junctions.seqs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular.RData')
writeXStringSet(c(CIRCS.junctions.seqs, LINEAR.junctions.seqs), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular.fa', format='fasta')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular.{RData,fa}
#
#  index for bwa:
#
#      bwa index -a bwtsw -p circRNAs_linear_vs_circular.fa circRNAs_linear_vs_circular.fa



#  collect and process the quantification results
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
lin_cir<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/ -maxdepth 3 -type f -wholename \'*/lin_vs_circ/counts.tsv\' -print', stdout=T)


#  add sequencing-sample-id (bid)
names(lin_cir)<-sub('^.*/(CB[^/]+)/.*$', '\\1', lin_cir)


#  order them
lin_cir<-lin_cir[ order(names(lin_cir)) ]


#  remove failed samples or samples not yet integrated in the cohort
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
lin_cir<-lin_cir[ intersect( names(lin_cir), meta[ !(failed), bid] ) ]
rm(meta)


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
save(lin.cir, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData




#########################
#
#
#  correlation analyses
#
#
#########################




#  correlations
#
#      1) We compute the Spearman correlation between circular junction counts and CIRI2 circRNA counts.
#
#      2) The mean circular junction count per isoform and mean linear junction count over junctions outside the corresponding circular junctions
#         are used to summarize at the gene level the circRNA expression and mRNA expression and compute their Spearman correlation WITHIN SAMPLES
#         for all the commonly expressed pairs.
#
#      3) The mean circular junction count per isoform and mean linear junction count over junctions outside the corresponding circular junctions
#         are used to summarize at the gene level the circRNA expression and mRNA expression and compute their Spearman correlations ACROSS SAMPLES
#         for all the commonly expressed in at least 25% of the samples.
#
#      4) Similar to 3) but excluding MNA samples and reidentifying the circular+linear junctions commonly expressed in at least 25% 
#         of the samples.
#
#      5) Using only MNA samples we compute the correlations for the circRNA/mRNA pairs of 4).
#
#      6) Using only the HR_nMNA samples we compute the correlations similar to 4) and then for the same pairs we compute the correlations in the 
#         MNA samples. 
#
#      We find:
#
#          * In the MNA tumors there is no significant difference in the distribution of correlations for the correlated and uncorrelated groups 
#            defined across the rest of the tumors.
#
#          * In the MNA tumors there is a significant decrease in correlations both for the correlated and uncorrelated groups defined across 
#            the HR_nMNA tumors (~ 2 fold more correlated pairs than using all the non-MNA tumors).
#
#          * Irrespective of the way the group of correlated and anticorrelated pairs is determined, we always find the circRNA CPMs of the correlated
#            group to be HIGHER than the CPMs of the uncorrelated group and vice versa for the mRNA CPMs, i.e they are found LOWER in the correlated 
#            group compared to the uncorrelated group.
#
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(robustbase)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  [run once] compute correlations with CIRI2 counts
#             compute within sample correlations
#             compute across sample correlations
#{{{

#  load the results 
#  keep only the neuroblastoma tumors since the circular junctions are based on them
#  order them according to risk group
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
m<-intersect(names(lin.cir), meta[ !is.na(risk_group), bid ])
meta<-meta[ bid %in% m ]
lin.cir<-lin.cir[ meta$bid ]
rm(m, meta.prefailed)


#  load unified cohort of circRNAs
l<-ls()
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
CIRCS<-CIRCS.all[ CIRCS.all$bid %in% meta$bid ]
CIRCS<-data.table(data.frame(CIRCS))[, c('circ_name', 'gene_id', 'gene_name', 'jc_count', 'bid'), with=F]
rm(list=setdiff(ls(), c(l, 'CIRCS')))


#  compute the across isoforms mean circRNA and external linear junction counts 
#  remove circRNAs not found expressed anywhere
ci<-do.call(rbind, lin.cir)
ci<-dcast(ci, bid ~ gene_name, value.var='c.count', fun.aggregate=mean, na.rm=T)
ci<-t(data.frame(ci[, -1], row.names=ci[, bid], check.names=F))
ci<-ci[ apply(ci, 1, function(x){ any(!is.na(x)) }), names(lin.cir)]  #  keep the original sample order as well
li<-do.call(rbind, lin.cir)
li<-dcast(li, bid ~ gene_name, value.var='l.count.out', fun.aggregate=mean, na.rm=T)
li<-t(data.frame(li[, -1], row.names=li[, bid], check.names=F))
li<-li[ rownames(ci), colnames(ci) ]
stopifnot( all.equal( rownames(ci), rownames(li) ) ) 
stopifnot( all.equal( colnames(ci), colnames(li) ) ) 


#  compute the across isoforms mean circRNA and external linear junction CPMs 
x<-do.call(rbind, lin.cir)[, c('ci.cpm', 'li.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)]
ci.cpm<-dcast(x, bid ~ gene_name, value.var='ci.cpm', fun.aggregate=mean, na.rm=T)
ci.cpm<-t(data.frame(ci.cpm[, -1], row.names=ci.cpm[, bid], check.names=F))
li.cpm<-dcast(x, bid ~ gene_name, value.var='li.cpm', fun.aggregate=mean, na.rm=T)
li.cpm<-t(data.frame(li.cpm[, -1], row.names=li.cpm[, bid], check.names=F))
ci.cpm<-ci.cpm[rownames(ci), colnames(ci)]
li.cpm<-li.cpm[rownames(li), colnames(li)]
stopifnot( all.equal( rownames(ci.cpm), rownames(li.cpm) ) ) 
stopifnot( all.equal( colnames(ci.cpm), colnames(li.cpm) ) ) 
rm(x)


#  compute the across isoforms mean circRNA CIRI counts
#  keep only the genes whose circRNAs have been identified by the linear vs circular junction analysis
ci.ciri<-data.frame(dcast(CIRCS, bid ~ gene_name, value.var='jc_count', fun.aggregate=mean, na.rm=T), check.names=F)
rownames(ci.ciri)<-ci.ciri$bid
ci.ciri$bid<-NULL
ci.ciri<-t(ci.ciri)
ci.ciri<-ci.ciri[ rownames(ci), colnames(ci)]
stopifnot( all.equal( rownames(ci), rownames(ci.ciri) ) )
stopifnot( all.equal( colnames(ci), colnames(ci.ciri) ) )


#  [CIRI2] Spearman correlations WITHIN each sample and each risk group of the commonly expressed circular junction counts and CIRI2 counts 
#{{{

#R<-setNames(vector('numeric', length(lin.cir)), names(lin.cir))
#for( n in seq_along(lin.cir)){
#
#    #  circular junctions
#    x<-setNames(lin.cir[[n]]$c.count, lin.cir[[n]]$circ_name)
#
#    
#    #  CIRI2 circular junctions
#    y<-CIRCS[ bid %in% names(lin.cir)[n] & circ_name %in% names(x) ]
#    y<-setNames( y$jc_count, y$circ_name )
#
#
#    #  keep commonly expressed junctions
#    common<-intersect(names(x), names(y))
#    x<-x[ common ]
#    y<-y[ common ]
#    stopifnot(all.equal(names(x), names(y)))
#    rm(common)
#
#
#    #  compute robust R
#    #l<-lmrob(y~x, na.action=na.exclude, control=lmrob.control(max.it=100, maxit.scale=1000, k.max=1000, scale.tol=1e-6, solve.tol=1e-6, refine.tol=1e-6))
#    #R[n]<-sign(coef(l)[2])*sqrt(summary(l)$r.squared)
#
#
#    R[n]<-cor(x, y, use='complete.obs', method='spearman')
#}


#  compute the within-sample correlations
R<-diag(cor(ci, ci.ciri, use='pairwise.complete.obs', method='spearman'))


#  barplot of R
x11(width=25, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.5, 8.5, 4.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
stopifnot( all.equal(meta$bid, names(R) ) )
B<-R
B.cl<-setNames(meta$col, meta$risk_group)
YTICK<-pretty(c(0, 1.0), 4)
bp<-barplot(B, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, tail(YTICK,1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col='white', axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
for(i in seq_along(B)){ 
    b<-B
    b[-i]<-NA    #  remove all values but the current 
    barplot(b, border='white', col=B.cl[i], axes=F, axisnames=F, beside=F, yaxt='n', add=T)
}
axis(2, at=YTICK, line=0, cex.axis=2.8)
#mtext(text=sub('-11-R01', '', names(B)), side=1, line=0, at=bp, las=2, adj=1, cex=1.8, col=B.cl)
mtext('Correlation', side=2, line=5, padj=-0.2, las=0, cex=2.4)
legend(x=par('usr')[2]*0.001, y=par('usr')[4]*1.21, legend=unique(names(B.cl)), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
#dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_cor_circular_vs_CIRI2_per_sample.svg', width=40, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  boxplot per risk_group 
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(13.0, 8.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(R, factor(meta$risk_group, levels=unique(meta$risk_group)))
B.cl<-unique(meta[, c('risk_group', 'col')])
B.cl<-setNames(B.cl$col, B.cl$risk_group)
YTICK<-pretty(c(0, sapply(B, max)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Correlation', side=2, line=6, padj=+0.4, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_cor_circular_vs_CIRI2_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}


#  Spearman correlations WITHIN samples
#{{{

#R<-setNames(vector('numeric', length(lin.cir)), names(lin.cir))
#for( n in seq_along(lin.cir)){
#    #  commony expressed
#    l<-lin.cir[[n]][ !is.na(c.count) & !is.na(l.count.out) ]
#    x<-setNames(l$c.count, l$circ_name)
#    y<-setNames(l$l.count.out, l$circ_name)
#    #l<-lmrob(y~x, na.action=na.exclude, control=lmrob.control(max.it=100, maxit.scale=1000, k.max=1000, scale.tol=1e-6, solve.tol=1e-6, refine.tol=1e-6))
#    #R[n]<-sign(coef(l)[2])*sqrt(summary(l)$r.squared)
#    R[n]<-cor(x, y, use='complete.obs', method='spearman')
#}


#  compute the within-sample correlations
R<-diag(cor(ci, li, use='pairwise.complete.obs', method='spearman'))


#  barplot of R
x11(width=25, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(7.5, 8.0, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
stopifnot( all.equal(meta$bid, names(R) ) )
B<-R
B.cl<-setNames(meta$col, meta$risk_group)
YTICK<-pretty(range(B), 4)
bp<-barplot(B, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, tail(YTICK,1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col='white', axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=range(YTICK), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
for(i in seq_along(B)){ 
    b<-B
    b[-i]<-NA    #  remove all values but the current 
    barplot(b, border='white', col=B.cl[i], axes=F, axisnames=F, beside=F, yaxt='n', add=T)
}
axis(2, at=YTICK, line=0, cex.axis=2.4)
#mtext(text=sub('-11-R01', '', names(B)), side=1, line=0, at=bp, las=2, adj=1, cex=1.8, col=B.cl)
mtext('Correlation', side=2, line=5, padj=-0.6, las=0, cex=2.4)
legend(x=0.5*par('usr')[2], y=par('usr')[4]*1.001, legend=unique(names(B.cl)), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
#dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_cor_linear_vs_circular_junctions_per_sample.svg', width=40, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  boxplot per risk_group 
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(9.0, 8.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(R, factor(meta$risk_group, levels=unique(meta$risk_group)))
B.cl<-unique(meta[, c('risk_group', 'col')])
B.cl<-setNames(B.cl$col, B.cl$risk_group)
YTICK<-pretty(range(sapply(B, range)), 4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Correlation', side=2, line=7, padj=+0.4, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_cor_linear_vs_circular_junctions_circRNAs_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}


#  Spearman correlations ACROSS samples 
#{{{

#  correlation containers
cl.cor<-list()


#  [all samples] using only the commonly expressed circRNAs and mRNAs in at least 25% of the samples used for the analysis
#{{{

#  keep only circular and linear junctions found COMMONLY expressed in at least 25% of the samples
keep<-apply(t((apply(ci, 1, is.na) | apply(li, 1, is.na))), 1, function(x){ sum(x)/length(x)<0.75 })
table(keep)
# keep
# FALSE  TRUE 
#     8  2212 
ci.keep<-ci[ keep, ]
li.keep<-li[ keep, ]
stopifnot( all.equal( rownames(ci.keep), rownames(li.keep) ) )
rm(keep)


#  compute the correlations across the diagonal only ( much faster than diag(cor(...)) )
R<-setNames(rep(NA, nrow(ci.keep)), rownames(ci.keep))
for( n in seq_along(R)){
    R[n]<-cor(li.keep[n, ], ci.keep[n, ], use='pairwise.complete.obs', method='spearman')
}
table(strong_positive=R>=0.6, negative=R<0.0)  #  there is no strongly negative correlation
#                negative
# strong_positive FALSE TRUE
#           FALSE  2002   71
#           TRUE    139    0
cl.cor[['all']]<-R
rm(R)

#}}}


#  [no MNA samples] using only the commonly expressed circRNAs and mRNAs in at least 25% of the samples used for the analysis
#{{{

#  select the samples of interest
keep<-meta[ grep('MNA', risk_group, invert=T), bid ] 
ci.keep<-ci[, keep]
li.keep<-li[, keep]
rm(keep)


#  keep only circular and linear junctions found COMMONLY expressed in at least 25% of the samples
keep<-apply(t((apply(ci.keep, 1, is.na) | apply(li.keep, 1, is.na))), 1, function(x){ sum(x)/length(x)<0.75 })
table(keep)
# keep
# FALSE  TRUE 
#    30  2190 
ci.keep<-ci.keep[ keep, ]
li.keep<-li.keep[ keep, ]
stopifnot( all.equal( rownames(ci.keep), rownames(li.keep) ) )
rm(keep)


#  compute the correlations across the diagonal only ( much faster than diag(cor(...)) )
R<-setNames(rep(NA, nrow(ci.keep)), rownames(ci.keep))
for( n in seq_along(R)){
    R[n]<-cor(li.keep[n, ], ci.keep[n, ], use='pairwise.complete.obs', method='spearman')
}
table(strong_positive=R>=0.6, negative=R<0.0)  #  there is no strongly negative correlation
#                negative
# strong_positive FALSE TRUE
#           FALSE  1888  123
#           TRUE    179    0
cl.cor[['nMNA']]<-R
rm(R)

#}}}


#  [MNA samples] using the circRNAs and mRNAs of the "no MNA" sample list
#{{{

#  select the samples of interest
keep<-meta[ grep('MNA', risk_group), bid ] 
ci.keep<-ci[, keep]
li.keep<-li[, keep]
rm(keep)


#  keep only circular and linear junctions of the 'no MNA' list
keep<-names(cl.cor[['nMNA']])
ci.keep<-ci.keep[ keep, ]
li.keep<-li.keep[ keep, ]
stopifnot( all.equal( rownames(ci.keep), rownames(li.keep) ) )
rm(keep)


#  compute the correlations across the diagonal only ( much faster than diag(cor(...)) )
R<-setNames(rep(NA, nrow(ci.keep)), rownames(ci.keep))
for( n in seq_along(R)){
    R[n]<-cor(li.keep[n, ], ci.keep[n, ], use='pairwise.complete.obs', method='spearman')
}
table(strong_positive=R>=0.6, negative=R<0.0)  #  there is no strongly negative correlation
#                negative
# strong_positive FALSE TRUE
#           FALSE  1857  153
#           TRUE    180    0
cl.cor[['MNA']]<-R
rm(R)

#}}}


#  [HR_nMNA] using only the commonly expressed circRNAs and mRNAs in at least 25% of the samples used for the analysis
#{{{

#  select the samples of interest
keep<-meta[ grep('HR_nMNA', risk_group), bid ] 
ci.keep<-ci[, keep]
li.keep<-li[, keep]
rm(keep)


#  keep only circular and linear junctions found COMMONLY expressed in at least 25% of the samples
keep<-apply(t((apply(ci.keep, 1, is.na) | apply(li.keep, 1, is.na))), 1, function(x){ sum(x)/length(x)<0.75 })
table(keep)
# keep
# FALSE  TRUE 
#     9  2211 
ci.keep<-ci.keep[ keep, ]
li.keep<-li.keep[ keep, ]
stopifnot( all.equal( rownames(ci.keep), rownames(li.keep) ) )
rm(keep)


#  compute the correlations across the diagonal only ( much faster than diag(cor(...)) )
R<-setNames(rep(NA, nrow(ci.keep)), rownames(ci.keep))
for( n in seq_along(R)){
    R[n]<-cor(li.keep[n, ], ci.keep[n, ], use='pairwise.complete.obs', method='spearman')
}
table(strong_positive=R>=0.6, negative=R<0.0)  #  there is no strongly negative correlation
#                negative
# strong_positive FALSE TRUE
#           FALSE  1722  197
#           TRUE    292    0
cl.cor[['HR_nMNA']]<-R
rm(R)

#}}}


#  [MNA_HR_nMNA] using the circRNAs and mRNAs of the "HR_nMNA" sample list
#{{{

#  select the samples of interest
keep<-meta[ grep('MNA', risk_group), bid ] 
ci.keep<-ci[, keep]
li.keep<-li[, keep]
rm(keep)


#  keep only circular and linear junctions of the 'no MNA' list
keep<-names(cl.cor[['HR_nMNA']])
ci.keep<-ci.keep[ keep, ]
li.keep<-li.keep[ keep, ]
stopifnot( all.equal( rownames(ci.keep), rownames(li.keep) ) )
rm(keep)


#  compute the correlations across the diagonal only ( much faster than diag(cor(...)) )
R<-setNames(rep(NA, nrow(ci.keep)), rownames(ci.keep))
for( n in seq_along(R)){
    R[n]<-cor(li.keep[n, ], ci.keep[n, ], use='pairwise.complete.obs', method='spearman')
}
table(strong_positive=R>=0.6, negative=R<0.0)  #  there is no strongly negative correlation
#                negative
# strong_positive FALSE TRUE
#           FALSE  1873  156
#           TRUE    182    0
cl.cor[['MNA_HR_nMNA']]<-R
rm(R)

#}}}


#  save
save(ci, li, ci.cpm, li.cpm, meta, cl.cor, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/cor_circular_vs_linear_junctions.RData')

#}}}

#}}}


#  load back
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/cor_circular_vs_linear_junctions.RData')


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  identify the neuroblastoma-specific among the high-correlating
#{{{

#  load the neuroblastoma-specific circRNAs without polluting the session
e<-new.env()
local({load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_neuroblastoma-specific.RData')}, envir=e)
nb.circs<-get('CIRCS', envir=e)
rm(e)


#  which highly correlating are also neuroblastoma-specific?
R.high<-cl.cor[['all']][ cl.cor[['all']]>=0.6 ]
R.low<-cl.cor[['all']][ cl.cor[['all']]<0.6 ]
h<-sort( R.high[ names(R.high) %in% nb.circs$gene_name ], decreasing=T )
l<-sort( R.low[ names(R.low) %in% nb.circs$gene_name ], decreasing=T )


#  clean up
rm(nb.circs, l, h, R.high, R.low)

#}}}


#  [across sample means] scatterplot and correlations 
#{{{

#  take the mean across samples
x<-log10(rowMeans(ci, na.rm=T))
y<-log10(rowMeans(li, na.rm=T))
l<-lm( y ~ x )


#  scatterplot
x11(width=12, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
XLIM<-pretty(range(x), 5)
YLIM<-pretty(range(y), 5)
par(mar=c(6.5, 8.0, 1.5, 2.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.0, cex.axis=2.0)
den<-col2rgb(densCols(x, y, colramp=colorRampPalette(c('black', 'white'))))[1, ]+1L
cls<-colorRampPalette(c('#000099', '#00FEFF', '#45FE4F', '#FCFF00', '#FF9400', '#FF3100'))(256)[den]
plot(x[ order(den) ], y[order(den)], type='p', pch=19, cex=0.8, col=cls[order(den)], xlim=c(XLIM[1], tail(XLIM,1)), ylim=c(YLIM[1], tail(YLIM,1)), xlab='', ylab='')
abline(l, lty=1, lwd=4, col='black')
mtext(bquote(R^2 == .(format(summary(l)$r.squared, digits=2))), side=3, line=0, padj=+0.4, cex=2.0, xpd=NA)
mtext(expression(log[10]('mean counts')~'[linear]'), side=2, line=2, padj=-0.7, cex=2.0, las=3)
#mtext('mRNA expression', side=2, line=4, padj=+0.1, cex=1.8, las=3)  #  scaled down SVG for huge fonts
mtext(expression(log[10]('mean counts')~'[circular]'), side=1, line=4, padj=+0.2, cex=2.0, las=1)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/scatterplot_circular_vs_linear_junctions.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}


#  [all samples] histogram of across sample correlations 
#{{{

#  with highlight of the high correlated
par(mar=c(5.5, 7.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(min(cl.cor[['all']]), 1), 10)
h<-hist(cl.cor[['all']], breaks=seq(min(b), max(b), 0.1), plot=F)
R.high<-cl.cor[['all']][ cl.cor[['all']]>=0.6 ]
#YMAX<-max(pretty(c(1, max(h$counts)), 4))
YMAX<-500
h<-hist(cl.cor[['all']], breaks=h$breaks, col='grey39', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(b), add=F)
h.high<-hist(R.high, breaks=h$breaks, col='red3', border='white', add=T)
mtext('Correlation', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Number of circRNAs', side=2, line=5, padj=-0.2, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_cor_circular_vs_linear_junctions.svg', width=12, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  without highlight of the high correlated
par(mar=c(5.5, 7.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(min(cl.cor[['all']]), 1), 10)
h<-hist(cl.cor[['all']], breaks=seq(min(b), max(b), 0.1), plot=F)
#R.high<-cl.cor[['all']][ cl.cor[['all']]>=0.6 ]
#YMAX<-max(pretty(c(1, max(h$counts)), 4))
YMAX<-500
h<-hist(cl.cor[['all']], breaks=h$breaks, col='grey39', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(b), add=F)
#h.high<-hist(R.high, breaks=h$breaks, col='red3', border='white', add=T)
mtext('Correlation', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Number of circRNAs', side=2, line=5, padj=-0.2, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_cor_circular_vs_linear_junctions_no_highlight.svg', width=12, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [all samples] histogram of within sample correlations 
#{{{

#  compute the within-sample correlations
R<-diag(cor(ci, li, use='pairwise.complete.obs', method='spearman'))


#  histogram
par(mar=c(5.5, 6.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#b<-pretty(c(0, max(R)), 20)
b<-pretty(c(0, 0.4), 20)
h<-hist(R, breaks=b, plot=F)
#YMAX<-max(h$counts)
YMAX<-25
h<-hist(R, breaks=h$breaks, col='grey39', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(b), xaxp=c(0, max(b), 4), add=F)
mtext('Correlation', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Number of samples', side=2, line=4, padj=-0.4, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_cor-within_circular_vs_linear_junctions.svg', width=12, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [MNA vs the rest] are the correlated and uncorrelated pairs found by excluding the MNA samples differently correlated in the MNA tumors?
#{{{

#  correlated and uncorrelated pairs 
length(R.high<-cl.cor[['all']][ cl.cor[['all']]>=0.6 ])   #  139
length(R.low<-cl.cor[['all']][ cl.cor[['all']]<0.6 ])     #  2073


#  the corresponding MNA correlations
length(R.mna.high<-cl.cor[['MNA']][ names(cl.cor[['MNA']]) %in% names(R.high) ])   #  137
length(R.mna.low<-cl.cor[['MNA']][ names(cl.cor[['MNA']]) %in% names(R.low) ])     #  2053


#  Mann-Whitney U tests
wilcox.test(x=R.mna.high, y=R.high, alternative='two.sided')$p.value  #  0.532345572
wilcox.test(x=R.mna.low, y=R.low, alternative='two.sided')$p.value    #  0.07677253271


#  [correlated] histograms
par(mar=c(5.5, 7.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(min(R.high, R.mna.high), 1), 10)
h<-hist(c(R.high, R.mna.high), breaks=seq(min(b), max(b), 0.025), plot=F)
#YMAX<-max(pretty(c(1, max(h$counts)), 4))
YMAX<-60
XLIM<-range(b)
h.high<-hist(R.high, breaks=h$breaks, col='red3', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=XLIM, add=F, axes=F)
h.mna<-hist(R.mna.high[R.mna.high>min(XLIM)], breaks=h$breaks, col=adjustcolor(meta[ risk_group %in% 'MNA', col][1], alpha.f=0.95), border='white', add=T)  #  limit the x-axis
axis(1, at=seq(XLIM[1], XLIM[2], 0.1), labels=seq(XLIM[1], XLIM[2], 0.1))
axis(2, at=seq(0, YMAX, length.out=5), labels=seq(0, YMAX, length.out=5))
mtext('Correlation', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Frequency', side=2, line=5, padj=-0.1, las=0, cex=2.4)
legend('topleft', legend=c('nMNA', 'MNA'), col=c('red3', adjustcolor(meta[ risk_group %in% 'MNA', col][1], alpha.f=0.95)), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)


#  [correlated] ecdfs
par(mar=c(5.5, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(min(R.high, R.mna.high), 1), 10)
XLIM<-range(b)
h.high<-curve(ecdf(R.high)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(R.high), 100), ylab='', xlab='', pch=NA, col='red3', lty=1, lwd=12, main='', xlim=range(XLIM), ylim=c(0, 1))
h.mna<-curve(ecdf(R.mna.high[R.mna.high>min(XLIM)])(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(R.mna.high[R.mna.high>min(XLIM)]), 100), pch=NA, col=meta[ risk_group %in% 'MNA', col][1], lty=1, lwd=12, add=T)
mtext('Correlation', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Probability', side=2, line=5, padj=-0.1, las=0, cex=2.4)
legend('topleft', legend=c('nMNA', 'MNA'), col=c('red3', meta[ risk_group %in% 'MNA', col][1]), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_cor_circular_vs_linear_junctions_high_nMNA+MNA.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [uncorrelated] histograms
par(mar=c(5.5, 7.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(min(R.low, R.mna.low), 1), 10)
h<-hist(c(R.low, R.mna.low), breaks=seq(min(b), max(b), 0.025), plot=F)
#YMAX<-max(pretty(c(1, max(h$counts)), 4))
YMAX<-150
XLIM<-range(b)
h.low<-hist(R.low, breaks=h$breaks, col='#cccccc', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=XLIM, add=F, axes=F)
h.mna<-hist(R.mna.low[R.mna.low>min(XLIM)], breaks=h$breaks, col=adjustcolor(meta[ risk_group %in% 'MNA', col][1], alpha.f=0.95), border='white', add=T)  #  limit the x-axis
axis(1, at=seq(XLIM[1], XLIM[2], 0.2), labels=seq(XLIM[1], XLIM[2], 0.2))
axis(2, at=seq(0, YMAX, length.out=5), labels=seq(0, YMAX, length.out=5))
mtext('Correlation', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Frequency', side=2, line=5, padj=-0.1, las=0, cex=2.4)
legend('topleft', legend=c('nMNA', 'MNA'), col=c('#cccccc', adjustcolor(meta[ risk_group %in% 'MNA', col][1], alpha.f=0.95)), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)


#  [uncorrelated] ecdfs
par(mar=c(5.5, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(min(R.low, R.mna.low), 1), 10)
XLIM<-range(b)
h.low<-curve(ecdf(R.low)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(R.low), 100), ylab='', xlab='', pch=NA, col='#cccccc', lty=1, lwd=12, main='', xlim=range(XLIM), ylim=c(0, 1))
h.mna<-curve(ecdf(R.mna.low[R.mna.low>min(XLIM)])(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(R.mna.low[R.mna.low>min(XLIM)]), 100), pch=NA, col=meta[ risk_group %in% 'MNA', col][1], lty=1, lwd=12, add=T)
mtext('Correlation', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Probability', side=2, line=5, padj=-0.1, las=0, cex=2.4)
legend('topleft', legend=c('nMNA', 'MNA'), col=c('#cccccc', meta[ risk_group %in% 'MNA', col][1]), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_cor_circular_vs_linear_junctions_low_nMNA+MNA.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [MNA vs HR_nMNA] are the correlated and uncorrelated pairs found in HR_nMNA differently correlated in the MNA tumors?
#{{{

#  correlated and uncorrelated pairs 
length(R.high<-cl.cor[['HR_nMNA']][ cl.cor[['HR_nMNA']]>=0.6 ])   #  292
length(R.low<-cl.cor[['HR_nMNA']][ cl.cor[['HR_nMNA']]<0.6 ])     #  1919


#  the corresponding MNA correlations
length(R.mna.high<-cl.cor[['MNA_HR_nMNA']][ names(cl.cor[['MNA_HR_nMNA']]) %in% names(R.high) ])   #  292
length(R.mna.low<-cl.cor[['MNA_HR_nMNA']][ names(cl.cor[['MNA_HR_nMNA']]) %in% names(R.low) ])     #  1919


#  Mann-Whitney U tests
wilcox.test(x=R.mna.high, y=R.high, alternative='two.sided')$p.value  #  6.819851731e-24
wilcox.test(x=R.mna.low, y=R.low, alternative='two.sided')$p.value    #  0.003203832822


#  [correlated] histograms
par(mar=c(5.5, 7.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(min(R.high, R.mna.high), 1), 10)
h<-hist(c(R.high, R.mna.high), breaks=seq(min(b), max(b), 0.025), plot=F)
#YMAX<-max(pretty(c(1, max(h$counts)), 4))
YMAX<-80
XLIM<-range(b)
h.high<-hist(R.high, breaks=h$breaks, col='red3', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=XLIM, add=F, axes=F)
h.mna<-hist(R.mna.high[R.mna.high>min(XLIM)], breaks=h$breaks, col=adjustcolor(meta[ risk_group %in% 'MNA', col][1], alpha.f=0.95), border='white', add=T)  #  limit the x-axis
axis(1, at=seq(XLIM[1], XLIM[2], 0.1), labels=seq(XLIM[1], XLIM[2], 0.1))
axis(2, at=seq(0, YMAX, length.out=5), labels=seq(0, YMAX, length.out=5))
mtext('Correlation', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Frequency', side=2, line=5, padj=-0.1, las=0, cex=2.4)
legend('topleft', legend=c('HR_nMNA', 'MNA'), col=c('red3', adjustcolor(meta[ risk_group %in% 'MNA', col][1], alpha.f=0.95)), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)


#  [correlated] ecdfs
par(mar=c(5.5, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(min(R.high, R.mna.high), 1), 10)
XLIM<-range(b)
h.high<-curve(ecdf(R.high)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(R.high), 100), ylab='', xlab='', pch=NA, col='red3', lty=1, lwd=12, main='', xlim=range(XLIM), ylim=c(0, 1))
h.mna<-curve(ecdf(R.mna.high[R.mna.high>min(XLIM)])(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(R.mna.high[R.mna.high>min(XLIM)]), 100), pch=NA, col=meta[ risk_group %in% 'MNA', col][1], lty=1, lwd=12, add=T)
mtext('Correlation', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Probability', side=2, line=5, padj=-0.1, las=0, cex=2.4)
legend('topleft', legend=c('HR_nMNA', 'MNA'), col=c('red3', meta[ risk_group %in% 'MNA', col][1]), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_cor_circular_vs_linear_junctions_high_HR_nMNA+MNA.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [uncorrelated] histograms
par(mar=c(5.5, 7.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(min(R.low, R.mna.low), 1), 10)
h<-hist(c(R.low, R.mna.low), breaks=seq(min(b), max(b), 0.025), plot=F)
#YMAX<-max(pretty(c(1, max(h$counts)), 4))
YMAX<-120
XLIM<-range(b)
h.low<-hist(R.low, breaks=h$breaks, col='#cccccc', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=XLIM, add=F, axes=F)
h.mna<-hist(R.mna.low[R.mna.low>min(XLIM)], breaks=h$breaks, col=adjustcolor(meta[ risk_group %in% 'MNA', col][1], alpha.f=0.95), border='white', add=T)  #  limit the x-axis
axis(1, at=seq(XLIM[1], XLIM[2], 0.2), labels=seq(XLIM[1], XLIM[2], 0.2))
axis(2, at=seq(0, YMAX, length.out=5), labels=seq(0, YMAX, length.out=5))
mtext('Correlation', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Frequency', side=2, line=5, padj=-0.1, las=0, cex=2.4)
legend('topleft', legend=c('HR_nMNA', 'MNA'), col=c('#cccccc', adjustcolor(meta[ risk_group %in% 'MNA', col][1], alpha.f=0.95)), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)


#  [uncorrelated] ecdfs
par(mar=c(5.5, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(min(R.low, R.mna.low), 1), 10)
XLIM<-range(b)
h.low<-curve(ecdf(R.low)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(R.low), 100), ylab='', xlab='', pch=NA, col='#cccccc', lty=1, lwd=12, main='', xlim=range(XLIM), ylim=c(0, 1))
h.mna<-curve(ecdf(R.mna.low[R.mna.low>min(XLIM)])(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(R.mna.low[R.mna.low>min(XLIM)]), 100), pch=NA, col=meta[ risk_group %in% 'MNA', col][1], lty=1, lwd=12, add=T)
mtext('Correlation', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Probability', side=2, line=5, padj=-0.1, las=0, cex=2.4)
legend('topleft', legend=c('HR_nMNA', 'MNA'), col=c('#cccccc', meta[ risk_group %in% 'MNA', col][1]), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_cor_circular_vs_linear_junctions_low_HR_nMNA+MNA.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [all samples] distribution of circular/linear junction CPMs stratified by correlation group
#{{{

#  correlation groups
length(R.high<-names(cl.cor[['all']][ cl.cor[['all']]>=0.6 ]))  #  139
length(R.low<-names(cl.cor[['all']][ cl.cor[['all']]<0.6 ]))    #  2073
stopifnot( length(setdiff(R.high, rownames(ci.cpm)))==0 )
stopifnot( length(setdiff(R.low, rownames(ci.cpm)))==0 )


#  [circular junctions] boxplot of CPMs stratified by correlation group 
par(mar=c(2.5, 8.0, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-list(high=log10(1+c(ci.cpm[ R.high, ])), low=log10(1+c(ci.cpm[R.low, ])))
#
wilcox.test(x=B[['high']], y=B[['low']], alternative='greater')$p.value   #  0
#
B.cl<-setNames(c('red3', '#cccccc'), c('high', 'low'))
YTICK<-pretty(c(0, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'CPM')), side=2, line=6, padj=+0.4, las=0, cex=2.4)
mtext(names(B), side=1, line=0, at=seq(1, length(B), 1), las=0, adj=0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circular_junctions_cpms_cor_high-low.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  ecdf
#
par(mar=c(5.5, 8.0, 0.1, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
h<-curve(ecdf(B[['low']])(x), from=0, to=2, n=min(length(B[['low']]), 100), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(0, 2), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=0, to=2, n=min(length(B[['high']]), 100), ylab='', xlab='', pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=seq(0, 2, 0.4))
mtext(expression(log[10](1+'CPM')), side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=1.8, y.intersp=0.50, x.intersp=0.2, seg.len=0.5, xpd=NA)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_circular_junctions_cpms_cor_high-low.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [linear junctions] boxplot of CPMs stratified by correlation group 
par(mar=c(2.5, 6.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-list(high=log10(1+c(li.cpm[ R.high, ])), low=log10(1+c(li.cpm[R.low, ])))
#
wilcox.test(x=B[['high']], y=B[['low']], alternative='less')$p.value   #  2.338066784e-269
#
B.cl<-setNames(c('red3', '#cccccc'), c('high', 'low'))
YTICK<-pretty(c(0, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'CPM')), side=2, line=4, padj=+0.2, las=0, cex=2.4)
mtext(names(B), side=1, line=0, at=seq(1, length(B), 1), las=0, adj=0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_linear_junctions_cpms_cor_high-low.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  ecdf
#
par(mar=c(5.5, 8.0, 0.1, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
h<-curve(ecdf(B[['low']])(x), from=0, to=3, n=min(length(B[['low']]), 100), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(0, 3), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=0, to=3, n=min(length(B[['high']]), 100), ylab='', xlab='', pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=seq(0, 3, 0.5))
mtext(expression(log[10](1+'CPM')), side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=1.8, y.intersp=0.50, x.intersp=0.2, seg.len=0.5, xpd=NA)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_linear_junctions_cpms_cor_high-low.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [MNA vs the rest] distribution of circular/linear junction CPMs stratified by correlation group
#{{{

#  correlation groups
length(R.high<-names(cl.cor[['nMNA']][ cl.cor[['nMNA']]>=0.6 ]))  #  179
length(R.low<-names(cl.cor[['nMNA']][ cl.cor[['nMNA']]<0.6 ]))    #  2011
stopifnot( length(setdiff(R.high, rownames(ci.cpm)))==0 )
stopifnot( length(setdiff(R.low, rownames(ci.cpm)))==0 )


#  [circular junctions] boxplot of CPMs stratified by correlation group 
par(mar=c(2.5, 8.0, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-list(high=log10(1+c(ci.cpm[ R.high, ])), low=log10(1+c(ci.cpm[R.low, ])))
#
wilcox.test(x=B[['high']], y=B[['low']], alternative='greater')$p.value   #  0
#
B.cl<-setNames(c('red3', '#cccccc'), c('high', 'low'))
YTICK<-pretty(c(0, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'CPM')), side=2, line=6, padj=+0.4, las=0, cex=2.4)
mtext(names(B), side=1, line=0, at=seq(1, length(B), 1), las=0, adj=0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circular_junctions_cpms_cor_high-low_nMNA+MNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  ecdf
#
par(mar=c(5.5, 8.0, 0.1, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
h<-curve(ecdf(B[['low']])(x), from=0, to=2, n=min(length(B[['low']]), 100), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(0, 2), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=0, to=2, n=min(length(B[['high']]), 100), ylab='', xlab='', pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=seq(0, 2, 0.4))
mtext(expression(log[10](1+'CPM')), side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=1.8, y.intersp=0.50, x.intersp=0.2, seg.len=0.5, xpd=NA)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_circular_junctions_cpms_cor_high-low_nMNA+MNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [linear junctions] boxplot of CPMs stratified by correlation group 
par(mar=c(2.5, 6.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-list(high=log10(1+c(li.cpm[ R.high, ])), low=log10(1+c(li.cpm[R.low, ])))
#
wilcox.test(x=B[['high']], y=B[['low']], alternative='less')$p.value   #  8.071509396e-131
#
B.cl<-setNames(c('red3', '#cccccc'), c('high', 'low'))
YTICK<-pretty(c(0, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'CPM')), side=2, line=4, padj=+0.2, las=0, cex=2.4)
mtext(names(B), side=1, line=0, at=seq(1, length(B), 1), las=0, adj=0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_linear_junctions_cpms_cor_high-low_nMNA+MNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  ecdf
#
par(mar=c(5.5, 8.0, 0.1, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
h<-curve(ecdf(B[['low']])(x), from=0, to=3, n=min(length(B[['low']]), 100), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(0, 3), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=0, to=3, n=min(length(B[['high']]), 100), ylab='', xlab='', pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=seq(0, 3, 0.5))
mtext(expression(log[10](1+'CPM')), side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=1.8, y.intersp=0.50, x.intersp=0.2, seg.len=0.5, xpd=NA)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_linear_junctions_cpms_cor_high-low_nMNA+MNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [MNA vs HR_nMNA] distribution of circular/linear junction CPMs stratified by correlation group
#{{{

#  correlation groups
length(R.high<-names(cl.cor[['HR_nMNA']][ cl.cor[['HR_nMNA']]>=0.6 ]))  #  292
length(R.low<-names(cl.cor[['HR_nMNA']][ cl.cor[['HR_nMNA']]<0.6 ]))    #  1919
stopifnot( length(setdiff(R.high, rownames(ci.cpm)))==0 )
stopifnot( length(setdiff(R.low, rownames(ci.cpm)))==0 )


#  [circular junctions] boxplot of CPMs stratified by correlation group 
par(mar=c(2.5, 8.0, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-list(high=log10(1+c(ci.cpm[ R.high, ])), low=log10(1+c(ci.cpm[R.low, ])))
#
wilcox.test(x=B[['high']], y=B[['low']], alternative='greater')$p.value   #  0
#
B.cl<-setNames(c('red3', '#cccccc'), c('high', 'low'))
YTICK<-pretty(c(0, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'CPM')), side=2, line=6, padj=+0.4, las=0, cex=2.4)
mtext(names(B), side=1, line=0, at=seq(1, length(B), 1), las=0, adj=0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circular_junctions_cpms_cor_high-low_HR_nMNA+MNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  ecdf
#
par(mar=c(5.5, 8.0, 0.1, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
h<-curve(ecdf(B[['low']])(x), from=0, to=2, n=min(length(B[['low']]), 100), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(0, 2), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=0, to=2, n=min(length(B[['high']]), 100), ylab='', xlab='', pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=seq(0, 2, 0.4))
mtext(expression(log[10](1+'CPM')), side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=1.8, y.intersp=0.50, x.intersp=0.2, seg.len=0.5, xpd=NA)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_circular_junctions_cpms_cor_high-low_HR_nMNA+MNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [linear junctions] boxplot of CPMs stratified by correlation group 
par(mar=c(2.5, 6.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-list(high=log10(1+c(li.cpm[ R.high, ])), low=log10(1+c(li.cpm[R.low, ])))
#
wilcox.test(x=B[['high']], y=B[['low']], alternative='less')$p.value   #  2.371719907e-211
#
B.cl<-setNames(c('red3', '#cccccc'), c('high', 'low'))
YTICK<-pretty(c(0, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'CPM')), side=2, line=4, padj=+0.2, las=0, cex=2.4)
mtext(names(B), side=1, line=0, at=seq(1, length(B), 1), las=0, adj=0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_linear_junctions_cpms_cor_high-low_HR_nMNA+MNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  ecdf
#
par(mar=c(5.5, 8.0, 0.1, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
h<-curve(ecdf(B[['low']])(x), from=0, to=3, n=min(length(B[['low']]), 100), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(0, 3), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=0, to=3, n=min(length(B[['high']]), 100), ylab='', xlab='', pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=seq(0, 3, 0.5))
mtext(expression(log[10](1+'CPM')), side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=1.8, y.intersp=0.50, x.intersp=0.2, seg.len=0.5, xpd=NA)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_linear_junctions_cpms_cor_high-low_HR_nMNA+MNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/cor_circular_vs_linear_junctions.RData



#  correlation of circRNA CIRI2 counts summed at the gene level with featureCounts gene counts and 
#  correlation of circular junction counts with counts of linear junctions outside the circular junctions
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(robustbase)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load unified circRNAs and sum circRNA isoform counts per gene
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
ciri<-data.table(data.frame(mcols(CIRCS.all)[, c('jc_count', 'gene_name', 'gene_id', 'bid')])) 
ciri<-dcast(ciri, bid ~ gene_name, value.var='jc_count', fun.aggregate=sum)
rm(list=ls(pattern='GENE|gene|meta|circs'))


#  load featureCounts gene counts and keep only those associated with genes producing circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/featureCounts.RData')
fc<-unlist(totalrna[ unique(ciri$bid) ])
fc$bid<-sub('\\.[0-9]*$', '', rownames(fc))
rownames(fc)<-NULL
fc<-data.table(fc[, c('bid', 'gene_id', 'counts')])
gn<-setNames(colnames(ciri[, -1]), CIRCS$gene_id[ match(colnames(ciri[, -1]), CIRCS$gene_name) ])
fc<-fc[ gene_id %in% names(gn) ]
fc<-fc[, gene_name:=gn[ gene_id ]]
fc<-dcast(fc, bid ~ gene_name, value.var='counts', fun.aggregate=sum)
rm(totalrna, CIRCS, CIRCS.all, gn)


#  convert to matrices for easier use
x<-t(as.matrix(ciri[, -1]))
colnames(x)<-ciri$bid
ciri<-x
y<-t(as.matrix(fc[, -1]))
colnames(y)<-fc$bid
fc<-y
stopifnot( all.equal( colnames(ciri), colnames(fc) ) )
stopifnot( all.equal( rownames(ciri), rownames(fc) ) )
rm(x,y)


#  Spearman correlation across samples
R<-setNames(rep(NA, nrow(ciri)), rownames(ciri))
for( n in seq_along(R)){
    #l<-lmrob(fc[n, ]~ciri[n, ], na.action=na.exclude, control=lmrob.control(max.it=100, maxit.scale=1000, k.max=1000, scale.tol=1e-6, solve.tol=1e-6, refine.tol=1e-6))
    #R[n]<-sign(coef(l)[2])*sqrt(summary(l)$r.squared)
    R[n]<-cor(fc[n,], ciri[n,], use='complete.obs', method='spearman')
}


#  robust correlation
#R<-setNames(rep(NA, nrow(ciri)), rownames(ciri))
#for( n in seq_along(R)){
#    l<-lmrob(fc[n, ]~ciri[n, ], na.action=na.exclude, control=lmrob.control(max.it=100, maxit.scale=1000, k.max=1000, scale.tol=1e-6, solve.tol=1e-6, refine.tol=1e-6))
#    R[n]<-sign(coef(l)[2])*sqrt(summary(l)$r.squared)
#}


#  save 
save(ciri, fc, R, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/cor_CIRI2_vs_featureCounts.RData')


#  histograms and ecdfs together with the corresponding circular vs linear junction correlations summarized at the gene level as well
#{{{

#  load back
#  copy current correlations to avoid overwriting
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/cor_CIRI2_vs_featureCounts.RData')
R.<-R


#  load circular vs linear junction correlations summarized at the gene level
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/cor_circular_vs_linear_junctions.RData')
R<-cl.cor[['all']]
rm(li, ci, li.cpm, ci.cpm, cl.cor)


#  keep only the common genes between the two analyses 
common<-intersect(names(R.), names(R))
R<-R[ common ]
R.<-R.[ common ]
rm(common)


#  is the distribution of correlations based on mRNA expression significantly right-shifted?
wilcox.test(x=R., y=R, alternative='greater')$p.value  #  4.497343367e-127


#  histograms
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.5, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
h<-hist(R., freq=T, breaks=seq(-1, 1, 0.1), col=adjustcolor('burlywood4', alpha.f=0.8), border='white', xlab='', ylab='', main='', xlim=c(-0.3, 1), ylim=c(0, 600), xaxp=c(-0.2, 1, 6))
hist(R, freq=T, breaks=h$breaks, col=adjustcolor('darkorange', alpha.f=0.8), border='white', add=T) 
legend(x=par('usr')[1], y=par('usr')[4]*1.08, legend=c('circRNA vs gene', 'circular junction vs linear junction'), col=c('burlywood4', 'darkorange'), bty='n', lty=1, lwd=15, cex=1.8, y.intersp=0.50, x.intersp=0.1, seg.len=0.5, xpd=NA)
mtext('Correlation', side=1, line=4, las=0, padj=+0.2, cex=2.4, las=0)
mtext('Frequency', side=2, line=5, padj=-0.3, cex=2.4, las=0) 
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_cor_CIRI2_vs_featureCounts+circular_vs_linear.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  ecdf
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.5, 8.0, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
h<-curve(ecdf(R.)(x), from=-0.3, to=1.0, n=min(length(R.), 20), ylab='', xlab='', pch=NA, col='burlywood4', lty=1, lwd=12, main='', xlim=c(-0.3, 1), ylim=c(0, 1), xaxt='n')
curve(ecdf(R)(x), from=-0.3, to=1.0, n=min(length(R), 20), ylab='', xlab='', pch=NA, col='darkorange', lty=1, lwd=12, main='', xlim=c(-0.3, 1), ylim=c(0, 1), xaxt='n', add=T)
axis(1, at=seq(-0.2, 1, 0.2), labels=seq(-0.2, 1, 0.2))
mtext('Correlation', side=1, line=4, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
legend(x=par('usr')[1], y=par('usr')[4]*1.08, legend=c('circRNA vs gene', 'circular junction vs linear junction'), col=c('burlywood4', 'darkorange'), bty='n', lty=1, lwd=15, cex=1.8, y.intersp=0.50, x.intersp=0.2, seg.len=0.5, xpd=NA)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_cor_CIRI2_vs_featureCounts+circular_vs_linear.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/cor_CIRI2_vs_featureCounts.RData



#  [polyA] Spearman correlation of linear junction counts outside circular junctions with the corresponding polyA-seq counts
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(robustbase)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load the circular vs linear junction results and metadata
#  aggregate the linear junction counts under the gene_id, taking the mean would ignore as it should identical total linear junction counts 
#  that are present due to the different circRNA isoforms 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta.t<-meta.tum[ bid %in% names(lin.cir) ]
li.t<-unlist(lin.cir)
li.t<-li.t[, c('gene_id', 'l.count.out', 'total', 'bid'), with=F]
li.t<-dcast(li.t, bid ~ gene_id, value.var='l.count.out', fun.aggregate=mean)
rm(meta.tum,meta.cel,meta.prefailed,lin.cir)


#  load polyA-seq circular vs linear junction results and do the same
load('/fast/projects/peifer_wgs/work/work/2018-05-23_Fuchs_polyA/raw/circRNAs_linear_vs_circular_collected_results.RData')
load('/fast/projects/peifer_wgs/work/work/2018-05-23_Fuchs_polyA/raw/metadata.RData')
meta.p<-meta[ bid %in% names(lin.cir) & !is.na(risk_group) ]
li.p<-unlist(lin.cir)
li.p<-li.p[, c('gene_id', 'l.count.out', 'total', 'bid'), with=F]
li.p<-dcast(li.p, bid ~ gene_id, value.var='l.count.out', fun.aggregate=mean)
rm(lin.cir)


#  keep common totalRNA and polyA samples only
common<-intersect(meta.t$bid, meta.p$bid)
meta.t<-meta.t[match(common, bid), ]
meta.p<-meta.p[match(common, bid), ]
stopifnot( all.equal( meta.t$bid, meta.p$bid ) )
li.t<-li.t[ match(meta.t$bid, bid), ]
li.p<-li.p[ match(meta.p$bid, bid), ]
rm(common)


#  keep identical gene_ids and convert to matrices
common<-intersect(colnames(li.t[, -1]), colnames(li.p[, -1]))
x<-t(as.matrix(li.t[, -1]))
colnames(x)<-li.t$bid
li.t<-x[common, ]
x<-t(as.matrix(li.p[, -1]))
colnames(x)<-li.p$bid
li.p<-x[common, ]
stopifnot( all.equal( rownames(li.t), rownames(li.p) ) )
rm(x)


#  keep only junctions found COMMONLY expressed in at least 75% of samples (since we have only 10 tumors here)
stopifnot( all.equal( rownames(li.t), rownames(li.p) ) )
stopifnot( all.equal( colnames(li.t), colnames(li.p) ) )
keep<-apply(t((apply(li.t, 1, is.na) | apply(li.p, 1, is.na))), 1, function(x){ sum(x)/length(x)<0.25 })
li.t<-li.t[ keep, ]
li.p<-li.p[ keep, ]
stopifnot( all.equal( rownames(li.t), rownames(li.p) ) )
rm(keep)


#  Spearman correlation
R<-setNames(rep(NA, nrow(li.t)), rownames(li.t))
for( n in seq_along(R)){
    R[n]<-cor(li.t[n,], li.p[n,], use='complete.obs', method='spearman')
}
N<-names(sort(R[ R< -0.5 ] ))
summary(R)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.8666667  0.2121212  0.4424242  0.4058676  0.6363636  0.9757576 


#  histogram
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(-1, 1), 10)
h<-hist(R, breaks=seq(min(b), max(b), length.out=20), plot=F)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
h<-hist(R, breaks=h$breaks, col='grey39', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(b), add=F)
mtext('Correlation', side=1, line=4, padj=+0.1, las=0, cex=2.4)
mtext('Frequency', side=2, line=5, padj=-0.5, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_cor_linear_junctions_vs_polyA-seq_linear_junctions.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}




###################
#
#
#  analysis plots
#
#
###################




#  [tumors] distributions per risk group of circular and linear junction CPMs
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(robustbase)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load the circular vs linear junction results and metadata
#  aggregate the linear junction counts under the same gene_id by taking the mean which would ignore (as it should) identical 
#  total linear junction counts that are present due to the different circRNA isoforms 
#  aggregate the circular junction counts under the same gene_id by taking the maximum as representative isoform circRNA expression
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-meta.tum[ bid %in% names(lin.cir) ]
lin.cir<-unlist(lin.cir[ meta$bid ])
lin.cir<-lin.cir[, c('gene_id', 'c.count', 'l.count.out', 'total', 'bid'), with=F]
lin.cir<-lin.cir[, c('li_cpm', 'ci_cpm'):=list(l.count.out/total*1e6, c.count/total*1e6)]
li<-dcast(lin.cir, bid ~ gene_id, value.var='li_cpm', fun.aggregate=mean)
ci<-dcast(lin.cir, bid ~ gene_id, value.var='ci_cpm', fun.aggregate=max)
rm(lin.cir)


#  group samples in risk groups 
#  collapse CPMs into one pool per risk group and convert to a list
stopifnot( length(setdiff(li$bid, meta$bid))==0 )
stopifnot( length(setdiff(ci$bid, meta$bid))==0 )
li<-li[ match(meta$bid, bid),  ]
ci<-ci[ match(meta$bid, bid),  ]
stopifnot( all.equal(li$bid, meta$bid) )
stopifnot( all.equal(ci$bid, meta$bid) )
li[, risk_group:=meta$risk_group]
ci[, risk_group:=meta$risk_group]
li<-li[, bid:=NULL]
ci<-ci[, bid:=NULL]
li<-li[, .(pool=list(unlist(lapply(.SD, unlist)))), by=.(risk_group)]
ci<-ci[, .(pool=list(unlist(lapply(.SD, unlist)))), by=.(risk_group)]
li<-setNames(apply(li, 1, function(x){ x[['pool']] }), li$risk_group)
ci<-setNames(apply(ci, 1, function(x){ x[['pool']] }), ci$risk_group)


#  Mann-Whitney U tests
#
#  N.B. obviously when circRNA CPMs go down then mRNA CPMs go up for the MNA samples...
#
wilcox.test(x=ci[['HR_nMNA']], y=ci[['MNA']], alternative='greater')$p.value   #  1.393022218e-37
wilcox.test(x=li[['HR_nMNA']], y=li[['MNA']], alternative='less')$p.value      #  4.046180744e-08


#  recycle
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  collapse all risk groups (DFG grant version)
par(mar=c(2.0, 7.0, 0.8, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
options(scipen=0)
B<-list(ci=log10(1+unlist(ci, recursive=T)), li=log10(1+unlist(li, recursive=T)))
B.cl<-c('black', 'grey39')
YTICK<-pretty(c(0, max(ceiling(unlist(B)), na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.6), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10]('CPM')), side=2, line=4, padj=+0.2, las=0, cex=2.4)
mtext(text=c('circRNAs', 'mRNAs'), side=1, line=0, at=seq_along(B), las=1, adj=0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_linear_vs_circular_CPM.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  boxplots of circular and linear junctions groupped by risk group 
par(mar=c(8.5, 7.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
options(scipen=0)
B<-unlist(lapply(seq_along(li), function(n){ list(ci=log10(1+ci[[n]]), li=log10(1+li[[n]])) }), recursive=F, use.names=F)
B.cl<-rep(unique(meta$col), each=2)
YTICK<-pretty(c(0, max(ceiling(unlist(B)), na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 0.6), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10]('CPM')), side=2, line=4, padj=+0.2, las=0, cex=2.4)
mtext(text=rep(c('circRNAs', 'mRNAs'), length(B)), side=1, line=0, at=seq_along(B), las=2, adj=0.89, cex=2.4, col=B.cl)
legend('topleft', legend=unique(meta$risk_group), col=unique(meta$col), bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.1, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_linear_vs_circular_CPMs_per_risk_group.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}



#  [tumors] cirular/(1+external linear) ratios:
#
#               per risk group
#               per PI class
#               per ADRN/MES class
#               per MYCN expression class
#               per TERT expression class
#               per circRNA/mRNA correlation class
#
#  [MYCN Tet-induced system] cirular/(1+external linear) ratios across time points
#
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


#  [run once] collect all classifications together
#             save 
#{{{

#  load circular and linear junction counts
#  split the neuroblastoma tumors and MYCN Tet-inducible systems
#  order the tumors according to risk group
#  order the MYCN Tet-inducible system according to time points
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
m<-intersect(names(lin.cir), meta.tum[, bid ])
meta.tumors<-meta.tum[ bid %in% m ]
meta.tumors<-meta.tumors[, risk_group:=factor(risk_group, levels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))][ order(risk_group), ][, risk_group:=as.character(risk_group)]
lc.tumors<-unlist(lin.cir[ meta.tumors$bid ])
m<-intersect(names(lin.cir), meta.cel[ grepl('CB-SKNAS-TR-MYCN', bid) , bid ])
meta.tet<-meta.cel[ bid %in% m ]
meta.tet<-meta.tet[, time.point:=factor(sub('^.*([0-9]+h)$', '\\1', treatment), levels=c('4h', '48h'))][ order(time.point), ][, time.point:=NULL]
lc.tet<-unlist(lin.cir[ meta.tet$bid ])
rm(m, lin.cir, meta.tum, meta.cel, meta.prefailed)


#  convert all NA to zero
lc.tumors[is.na(c.count), c.count:=0]
lc.tumors[is.na(l.count.out), l.count.out:=0]
lc.tumors[is.na(l.count.out.max), l.count.out.max:=0]
lc.tet[is.na(c.count), c.count:=0]
lc.tet[is.na(l.count.out), l.count.out:=0]
lc.tet[is.na(l.count.out.max), l.count.out.max:=0]


#  compute ratio of circular/(1+external linear) counts
#
#  N.B. ignore the warning about "invalid .internal.selfref detected"
#
lc.tumors[, ratio:=c.count/(1+l.count.out)]  
lc.tet[, ratio:=c.count/(1+l.count.out)]


#  compute circular junction CPM and external linear junction CPM based on total counts of (linear+circular)
lc.tumors[, c('c.cpm', 'l.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)]
lc.tet[, c('c.cpm', 'l.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)]


#  add risk_group metadata
#  add treatment metadata
lc.tumors$risk_group<-meta.tumors$risk_group[ match( lc.tumors$bid , meta.tumors$bid ) ]
lc.tet$treatment<-meta.tet$treatment[ match( lc.tet$bid , meta.tet$bid ) ]


#  load PI class for the samples and annotate the data objects
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/proliferative_index.RData')
lc.tumors$pi.class<-tum[ match(lc.tumors$bid, bid), PI.class ]
lc.tumors$pi.class.col<-tum[ match(lc.tumors$bid, bid), PI.class.col ]
lc.tet$pi.class<-cel[ match(lc.tet$bid, bid), PI.class ]
lc.tet$pi.class.col<-cel[ match(lc.tet$bid, bid), PI.class.col ]
rm(cel, cel.tpm, tum.tpm, tum)


#  load ADRN/MES scores and MES classification for the tumor samples and annotate the data objects
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/adrenergic-mesenchymal_score.RData')
lc.tumors$adrn.score<-AM.tum[ match(lc.tumors$bid, bid), adrn.score ]
lc.tumors$mes.score<-AM.tum[ match(lc.tumors$bid, bid), mes.score ]
lc.tumors$adrn.class<-AM.tum[ match(lc.tumors$bid, bid), adrn.class ]
lc.tumors$adrn.class.col<-AM.tum[ match(lc.tumors$bid, bid), adrn.class.col ]
lc.tumors$mes.class<-AM.tum[ match(lc.tumors$bid, bid), mes.class ]
lc.tumors$mes.class.col<-AM.tum[ match(lc.tumors$bid, bid), mes.class.col ]
rm(AM.tum, AM.cel)


#  load MYCN class and annotate the data objects
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYCN_expression.RData')
lc.tumors$mycn.class<-tum[ match(lc.tumors$bid, bid), mycn.class ]
lc.tumors$mycn.class.col<-tum[ match(lc.tumors$bid, bid), mycn.class.col ]
lc.tet$mycn.class<-cel[ match(lc.tet$bid, bid), mycn.class ]
lc.tet$mycn.class.col<-cel[ match(lc.tet$bid, bid), mycn.class.col ]
rm(cel, cel.tpm, tum.tpm, tum)


#  load MYC class and annotate the data objects
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYC_expression.RData')
lc.tumors$myc.class<-tum[ match(lc.tumors$bid, bid), myc.class ]
lc.tumors$myc.class.col<-tum[ match(lc.tumors$bid, bid), myc.class.col ]
lc.tet$myc.class<-cel[ match(lc.tet$bid, bid), myc.class ]
lc.tet$myc.class.col<-cel[ match(lc.tet$bid, bid), myc.class.col ]
rm(cel, cel.tpm, tum.tpm, tum)


#  load TERT class and annotate the data objects
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/TERT_expression.RData')
lc.tumors$tert.class<-tum[ match(lc.tumors$bid, bid), tert.class ]
lc.tumors$tert.class.col<-tum[ match(lc.tumors$bid, bid), tert.class.col ]
lc.tet$tert.class<-cel[ match(lc.tet$bid, bid), tert.class ]
lc.tet$tert.class.col<-cel[ match(lc.tet$bid, bid), tert.class.col ]
rm(cel, cel.tpm, tum.tpm, tum)


#  load circular vs external linear junction correlation results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/cor_circular_vs_linear_junctions.RData')
length(R.high<-cl.cor[['all']][ cl.cor[['all']]>=0.6 ])  #  139
length(R.low<-cl.cor[['all']][ cl.cor[['all']]<0.6 ])    #  2073
lc.tumors[, cor.class:='NA']
lc.tumors[ gene_name %in% names(R.high), cor.class:='correlated']
lc.tumors[ gene_name %in% names(R.low), cor.class:='uncorrelated']
lc.tumors[, cor.class:=factor(cor.class, levels=c('correlated', 'uncorrelated', 'NA'))]
lc.tumors[, table(cor.class)]
# cor.class
#   correlated uncorrelated           NA 
#        51688       449176          936 
#
lc.tumors[, cor.class.col:='NA']
lc.tumors[ cor.class %in% 'correlated', cor.class.col:='#b21f1f']
lc.tumors[ cor.class %in% 'uncorrelated', cor.class.col:='#cccccc']
lc.tumors[, cor.class.col:=factor(cor.class.col, levels=c('#b21f1f', '#cccccc', 'NA'))]
lc.tet[, cor.class:='NA']
lc.tet[ gene_name %in% names(R.high), cor.class:='correlated']
lc.tet[ gene_name %in% names(R.low), cor.class:='uncorrelated']
lc.tet[, cor.class:=factor(cor.class, levels=c('correlated', 'uncorrelated', 'NA'))]
lc.tet[, table(cor.class)]
# cor.class
#   correlated uncorrelated           NA 
#         5964        51828          108 
#
lc.tet[, cor.class.col:='NA']
lc.tet[ cor.class %in% 'correlated', cor.class.col:='#b21f1f']
lc.tet[ cor.class %in% 'uncorrelated', cor.class.col:='#cccccc']
lc.tet[, cor.class.col:=factor(cor.class.col, levels=c('#b21f1f', '#cccccc', 'NA'))]
rm(R.high, R.low, cl.cor, ci, li, ci.cpm, li.cpm)


#  save
save(lc.tumors, lc.tet, meta.tumors, meta.tet, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_all_classifications_together.RData')

#}}}


#  load back
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_all_classifications_together.RData')


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  [tumors] ratios per risk group
#           ratios MNA vs HR_nMNA
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tumors$ratio), factor(lc.tumors$risk_group, levels=unique(lc.tumors$risk_group)))
B.cl<-unique(meta.tumors[, c('risk_group', 'col')])
B.cl<-setNames(B.cl$col, B.cl$risk_group)
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value  #  0
wilcox.test(x=B[['MNA']], y=B[['LR']], alternative='less')$p.value       #  1.002300741e-123


#  boxplots
par(mar=c(14.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(meta.tumors[, table(risk_group)][names(B)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1.00, cex=2.4, col=B.cl)


#  ecdfs 
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['ST4S']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['ST4S']]), 20), ylab='', xlab='', pch=NA, col=B.cl['ST4S'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['HR_nMNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['HR_nMNA'], lty=1, lwd=12, add=T)
curve(ecdf(B[['IMR']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['IMR'], lty=1, lwd=12, add=T)
curve(ecdf(B[['LR']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['LR'], lty=1, lwd=12, add=T)
curve(ecdf(B[['MNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['MNA'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
#legend('topleft', legend=paste0(names(B), ' (', as.numeric(meta.tumors[, table(risk_group)][names(B)]), ')'), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  ecdf for MNA vs HR_nMNA groups only
x<-lc.tumors[ lc.tumors$risk_group %in% c('MNA', 'HR_nMNA')]
B<-split(log2(1e-3+x$ratio), factor(x$risk_group, levels=unique(x$risk_group)))
B.cl<-unique(meta.tumors[, c('risk_group', 'col')])
B.cl<-setNames(B.cl$col, B.cl$risk_group)
B.cl<-B.cl[ names(B) ]
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['HR_nMNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['HR_nMNA']]), 20), ylab='', xlab='', pch=NA, col=B.cl['HR_nMNA'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['MNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['MNA'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_MNA_vs_HR_nMNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] mean ratios per risk group
#{{{

#  compute the mean ratio of each circRNA across the samples of the risk group
#  split into groups and compute log2(0.001 + mean ratio)
B<-lc.tumors[, .(ratio=mean(ratio)), by=.(risk_group, circ_name)]
B<-split(log2(1e-3+B$ratio), factor(B$risk_group, levels=unique(B$risk_group)))
B.cl<-unique(meta.tumors[, c('risk_group', 'col')])
B.cl<-setNames(B.cl$col, B.cl$risk_group)
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value  #  4.629e-66
wilcox.test(x=B[['MNA']], y=B[['LR']], alternative='less')$p.value       #  7.018e-30


#  boxplots
par(mar=c(14.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(meta.tumors[, table(risk_group)][names(B)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1.00, cex=2.4, col=B.cl)


#  ecdfs 
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['ST4S']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['ST4S']]), 20), ylab='', xlab='', pch=NA, col=B.cl['ST4S'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['HR_nMNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['HR_nMNA'], lty=1, lwd=12, add=T)
curve(ecdf(B[['IMR']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['IMR'], lty=1, lwd=12, add=T)
curve(ecdf(B[['LR']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['LR'], lty=1, lwd=12, add=T)
curve(ecdf(B[['MNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['MNA'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + Mean ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
#legend('topleft', legend=paste0(names(B), ' (', as.numeric(meta.tumors[, table(risk_group)][names(B)]), ')'), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_mean_ratio_circular_vs_linear_junctions_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  ecdf for MNA vs HR_nMNA groups only
B<-lc.tumors[ risk_group %in% c('HR_nMNA', 'MNA'), .(ratio=mean(ratio)), by=.(risk_group, circ_name)]
B<-split(log2(1e-3+B$ratio), factor(B$risk_group, levels=unique(B$risk_group)))
B.cl<-unique(meta.tumors[, c('risk_group', 'col')])
B.cl<-setNames(B.cl$col, B.cl$risk_group)
B.cl<-B.cl[ names(B) ]
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['HR_nMNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['HR_nMNA']]), 20), ylab='', xlab='', pch=NA, col=B.cl['HR_nMNA'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['MNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['MNA'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + Mean ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_mean_ratio_circular_vs_linear_junctions_MNA_vs_HR_nMNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] ratios PI high class per risk group
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tumors[ pi.class %in% 'high', ratio]), factor(lc.tumors[ pi.class %in% 'high', risk_group], levels=unique(lc.tumors$risk_group)))
B.cl<-unique(meta.tumors[, c('risk_group', 'col')])
B.cl<-setNames(B.cl$col, B.cl$risk_group)
#
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value  #  0
wilcox.test(x=B[['LR']], y=B[['MNA']], alternative='less')$p.value       #  2.262e-91


#  boxplots
par(mar=c(14.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(lc.tumors[ pi.class %in% 'high', .(risk_group=unique(risk_group)), by=.(bid)][, table(risk_group)][names(B)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1.00, cex=2.4, col=B.cl)


#  ecdfs 
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['ST4S']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['ST4S']]), 20), ylab='', xlab='', pch=NA, col=B.cl['ST4S'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['HR_nMNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['HR_nMNA'], lty=1, lwd=12, add=T)
curve(ecdf(B[['IMR']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['IMR'], lty=1, lwd=12, add=T)
curve(ecdf(B[['LR']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['LR'], lty=1, lwd=12, add=T)
curve(ecdf(B[['MNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['MNA'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=paste0(names(B), ' (', as.numeric(lc.tumors[ pi.class %in% 'high', .(risk_group=unique(risk_group)), by=.(bid)][, table(risk_group)][names(B)]), ')'), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
#legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_PI_high_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  ecdf for MNA vs HR_nMNA groups only
B<-split(log2(1e-3+lc.tumors[ pi.class %in% 'high' & risk_group %in% c('HR_nMNA', 'MNA'), ratio]), factor(lc.tumors[ pi.class %in% 'high' & risk_group %in% c('HR_nMNA', 'MNA'), risk_group], levels=c('HR_nMNA', 'MNA')))
B.cl<-unique(meta.tumors[, c('risk_group', 'col')])
B.cl<-setNames(B.cl$col, B.cl$risk_group)
B.cl<-B.cl[ names(B) ]
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['HR_nMNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['HR_nMNA']]), 20), ylab='', xlab='', pch=NA, col=B.cl['HR_nMNA'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['MNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['MNA'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_PI_high_MNA_vs_HR_nMNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] CPMs per risk group
#{{{

#  split into groups and compute log10(1+CPMs)
B<-split(log10(1+lc.tumors$c.cpm), factor(lc.tumors$risk_group, levels=unique(lc.tumors$risk_group)))
B.cl<-unique(meta.tumors[, c('risk_group', 'col')])
B.cl<-setNames(B.cl$col, B.cl$risk_group)
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value  #  0
wilcox.test(x=B[['MNA']], y=B[['LR']], alternative='less')$p.value       #  1.901e-94


#  boxplots
par(mar=c(10.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(sapply(B, range, na.rm=T), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1 +'CPM')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', as.numeric(meta.tumors[, table(risk_group)][names(B)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1.00, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1.00, cex=2.4, col=B.cl)


#  ecdfs 
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#YTICK<-pretty(sapply(B, range, na.rm=T), 5)
YTICK<-pretty(c(0, 1), 5)
h<-curve(ecdf(B[['ST4S']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['ST4S']]), 20), ylab='', xlab='', pch=NA, col=B.cl['ST4S'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['HR_nMNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['HR_nMNA'], lty=1, lwd=12, add=T)
curve(ecdf(B[['IMR']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['IMR'], lty=1, lwd=12, add=T)
curve(ecdf(B[['LR']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['LR'], lty=1, lwd=12, add=T)
curve(ecdf(B[['MNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['MNA'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[10](1+'CPM')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('bottomright', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_CPMs_circular_junctions_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  ecdf for MNA vs HR_nMNA groups only
x<-lc.tumors[ lc.tumors$risk_group %in% c('MNA', 'HR_nMNA')]
B<-split(log10(1+x$c.cpm), factor(x$risk_group, levels=unique(x$risk_group)))
B.cl<-unique(meta.tumors[, c('risk_group', 'col')])
B.cl<-setNames(B.cl$col, B.cl$risk_group)
B.cl<-B.cl[ names(B) ]
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(0, 1.0), 5)
h<-curve(ecdf(B[['HR_nMNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['HR_nMNA']]), 20), ylab='', xlab='', pch=NA, col=B.cl['HR_nMNA'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['MNA']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['MNA']]), 20), pch=NA, col=B.cl['MNA'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[10](1+'CPM')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('bottomright', legend=names(B), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_CPMs_circular_junctions_MNA_vs_HR_nMNA.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] ratios per PI class
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tumors$ratio), lc.tumors$pi.class)
B.cl<-setNames(levels(lc.tumors$pi.class.col), levels(lc.tumors$pi.class))
wilcox.test(x=B[['low']], y=B[['high']], alternative='greater')$p.value  #  0


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(pi.class)), by=.(bid)][, table(type)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs (only low/high)
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['low']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['low']]), 20), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['high']]), 20), pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(pi.class)), by=.(bid)][, table(type)]), ')'), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_per_proliferative_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] ratios per MES class
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tumors$ratio), lc.tumors$mes.class)
B.cl<-setNames(levels(lc.tumors$mes.class.col), levels(lc.tumors$mes.class))
wilcox.test(x=B[['low']], y=B[['high']], alternative='less')$p.value   #  1.35768085e-14


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(mes.class)), by=.(bid)][, table(type)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs (low/high)
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['low']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['low']]), 20), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['high']]), 20), pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(mes.class)), by=.(bid)][, table(type)]), ')'), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_per_MES_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] ratios per ADRN class
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tumors$ratio), lc.tumors$adrn.class)
B.cl<-setNames(levels(lc.tumors$adrn.class.col), levels(lc.tumors$adrn.class))
wilcox.test(x=B[['low']], y=B[['high']], alternative='less')$p.value   #  1.440847144e-173


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(adrn.class)), by=.(bid)][, table(type)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs (low/high)
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['low']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['low']]), 20), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['high']]), 20), pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(adrn.class)), by=.(bid)][, table(type)]), ')'), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_per_ADRN_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] ratios per MYCN class
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tumors$ratio), lc.tumors$mycn.class)
B.cl<-setNames(levels(lc.tumors$mycn.class.col), levels(lc.tumors$mycn.class))
wilcox.test(x=B[['low']], y=B[['high']], alternative='greater')$p.value  #  0


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(mycn.class)), by=.(bid)][, table(type)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs (low/high)
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['low']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['low']]), 20), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['high']]), 20), pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(mycn.class)), by=.(bid)][, table(type)]), ')'), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_per_MYCN_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] ratios PI high class per MYCN class
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tumors[ pi.class %in% 'high', ratio]), lc.tumors[ pi.class %in% 'high', mycn.class])
B.cl<-setNames(levels(lc.tumors$mycn.class.col), levels(lc.tumors$mycn.class))
wilcox.test(x=B[['low']], y=B[['high']], alternative='greater')$p.value  #  6.860157362e-191


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(lc.tumors[ pi.class %in% 'high', .(type=unique(mycn.class)), by=.(bid)][, table(type)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs (low/high)
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['low']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['low']]), 20), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['high']]), 20), pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=paste0(names(B), ' (', as.numeric(lc.tumors[ pi.class %in% 'high', .(type=unique(mycn.class)), by=.(bid)][, table(type)]), ')'), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_PI_high_per_MYCN_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] ratios per MYC class
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tumors$ratio), lc.tumors$myc.class)
B.cl<-setNames(levels(lc.tumors$myc.class.col), levels(lc.tumors$myc.class))
wilcox.test(x=B[['low']], y=B[['high']], alternative='less')$p.value  #  1.58683883e-150


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(myc.class)), by=.(bid)][, table(type)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs (low/high)
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['low']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['low']]), 20), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['high']]), 20), pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(myc.class)), by=.(bid)][, table(type)]), ')'), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_per_MYC_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] ratios MYCN low class per MYC class
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tumors[ mycn.class %in% 'low', ratio]), lc.tumors[ mycn.class %in% 'low', myc.class])
B.cl<-setNames(levels(lc.tumors$myc.class.col), levels(lc.tumors$myc.class))
wilcox.test(x=B[['high']], y=B[['low']], alternative='greater')$p.value  #  5.305356458e-34 (TESTING FOR MORE circRNAs in the MYC high class)


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(lc.tumors[ mycn.class %in% 'low', .(type=unique(myc.class)), by=.(bid)][, table(type)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs (low/high)
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['low']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['low']]), 20), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['high']]), 20), pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=paste0(names(B), ' (', as.numeric(lc.tumors[ mycn.class %in% 'low', .(type=unique(myc.class)), by=.(bid)][, table(type)]), ')'), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_low_per_MYC_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] ratios per TERT class
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tumors$ratio), lc.tumors$tert.class)
B.cl<-setNames(levels(lc.tumors$tert.class.col), levels(lc.tumors$tert.class))
wilcox.test(x=B[['low']], y=B[['high']], alternative='greater')$p.value  #  0


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(pi.class)), by=.(bid)][, table(type)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs (low/high)
par(mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, 0), 5)
h<-curve(ecdf(B[['low']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['low']]), 20), ylab='', xlab='', pch=NA, col=B.cl['low'], lty=1, lwd=12, main='', xlim=c(min(YTICK), max(YTICK)), ylim=c(0, 1), xaxt='n')
curve(ecdf(B[['high']])(x), from=min(YTICK), to=max(YTICK), n=min(length(B[['high']]), 20), pch=NA, col=B.cl['high'], lty=1, lwd=12, add=T)
axis(1, at=YTICK)
mtext('Probability', side=2, line=5, padj=-0.5, cex=2.4, las=0) 
mtext(expression(log[2]('0.001 + ratio')), side=1, line=4, padj=+0.4, las=0, cex=2.4)
legend('topleft', legend=paste0(names(B), ' (', as.numeric(lc.tumors[, .(type=unique(tert.class)), by=.(bid)][, table(type)]), ')'), col=B.cl, bty='n', lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_per_TERT_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}


#  [tumors] ratios per correlated class 
#{{{

#  split into groups and compute log2(0.001 + ratio)
#  remove NA class of circRNAs that dropped out from the correlation analysis
B<-split(log2(1e-3+lc.tumors$ratio), lc.tumors$cor.class)
B['NA']<-NULL
B.cl<-setNames(levels(lc.tumors$cor.class.col), levels(lc.tumors$cor.class))
B.cl<-B.cl[names(B)]
B.n<-lc.tumors[, .(cors=list(table(cor.class)[c('correlated', 'uncorrelated')])), by=.(bid)][1, unlist(cors)]  #  number of correlated/uncorrelated pairs
wilcox.test(x=B[['correlated']], y=B[['uncorrelated']], alternative='two.sided')$p.value  #  0


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs
my_ecdfs(B, B.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='topleft', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_per_cor_group.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)

#}}}


#  [tumors MNA vs HR_nMNA] ratios of the correlated and uncorrelated classes
#{{{


#  both classes
#{{{

#  split into groups and compute log2(0.001 + ratio)
x<-lc.tumors[ risk_group %in% c('HR_nMNA', 'MNA') ]
B<-split(x[, log2(1e-3+ratio)], x[, risk_group])
B.cl<-setNames(meta.tumors[risk_group %in% c('HR_nMNA', 'MNA'), unique(col)], meta.tumors[risk_group %in% c('HR_nMNA', 'MNA'), unique(risk_group)])
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='two.sided')$p.value  #  0


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs
my_ecdfs(B, B.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='topleft', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_MNA_vs_HR_nMNA.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)

#}}}


#  correlated
#{{{

#  split into groups and compute log2(0.001 + ratio)
x<-lc.tumors[ risk_group %in% c('HR_nMNA', 'MNA') & cor.class %in% 'correlated']
B<-split(x[, log2(1e-3+ratio)], x[, risk_group])
B.cl<-setNames(meta[risk_group %in% c('HR_nMNA', 'MNA'), unique(col)], meta[risk_group %in% c('HR_nMNA', 'MNA'), unique(risk_group)])
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='two.sided')$p.value  #  2.467212069e-120


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs
my_ecdfs(B, B.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='topleft', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_correlated_group_MNA_vs_HR_nMNA.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)

#}}}


#  uncorrelated
#{{{

#  split into groups and compute log2(0.001 + ratio)
x<-lc.tumors[ risk_group %in% c('HR_nMNA', 'MNA') & cor.class %in% 'uncorrelated']
B<-split(x[, log2(1e-3+ratio)], x[, risk_group])
B.cl<-setNames(meta[risk_group %in% c('HR_nMNA', 'MNA'), unique(col)], meta[risk_group %in% c('HR_nMNA', 'MNA'), unique(risk_group)])
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='two.sided')$p.value  #  2.467212069e-120


#  boxplots
par(mar=c(2.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=1, adj=0.5, cex=2.4, col=B.cl)


#  ecdfs
my_ecdfs(B, B.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='topleft', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_uncorrelated_group_MNA_vs_HR_nMNA.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)

#}}}

#}}}


#  [MYCN Tet-inducible system] ratios per condition
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tet$ratio), factor(lc.tet$treatment, levels=unique(lc.tet$treatment)))
B.cl<-unique(meta.tet[, c('treatment', 'col')])
B.cl<-setNames(B.cl$col, B.cl$treatment)
B.cl['+Tet 4h']<-'orange3'  #  manually change 
B.cl['ETOH 4h']<-'palegreen4'  #  manually change 
B.cl['+Tet 48h']<-'orangered3'  #  manually change 
B.cl['ETOH 48h']<-'seagreen4'  #  manually change
wilcox.test(x=B[['+Tet 4h']], y=B[['ETOH 4h']], alternative='less')$p.value    #  0.830640467
wilcox.test(x=B[['+Tet 48h']], y=B[['ETOH 48h']], alternative='less')$p.value  #  0.003901210755


#  boxplots
par(mar=c(13.0, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-10, sapply(B, max, na.rm=T)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', as.numeric(meta.tet[, table(treatment)][names(B)]), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1.00, cex=2.4, col=B.cl)


#  [altogether] ecdfs 
my_ecdfs(B, B.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='topleft', svg.file=NULL, mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)


#  [4h] ecdfs 
b<-B[c('ETOH 4h', '+Tet 4h')]
b.cl<-B.cl[c('ETOH 4h', '+Tet 4h')]
my_ecdfs(b, b.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='topleft', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_+Tet_ETOH_4h.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
rm(b, b.cl)


#  [48h] ecdfs 
b<-B[c('ETOH 48h', '+Tet 48h')]
b.cl<-B.cl[c('ETOH 48h', '+Tet 48h')]
my_ecdfs(b, b.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='topleft', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_+Tet_ETOH_48h.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
rm(b, b.cl)

#}}}


#  [MYCN Tet-inducible system] ratios of the correlated and uncorrelated classes
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tet$ratio), factor(paste(lc.tet$treatment, lc.tet$cor.class), levels=unique(paste(lc.tet$treatment, lc.tet$cor.class))))
B<-B[grep('NA$', names(B), invert=T)]  #  remove NA correlation class from circRNA/mRNA dropout pairs
B.cl<-unique(meta.tet[, c('treatment', 'col')])
B.cl<-setNames(B.cl$col, B.cl$treatment)
B.cl['+Tet 4h']<-'orange3'  #  manually change 
B.cl['ETOH 4h']<-'palegreen4'  #  manually change 
B.cl['+Tet 48h']<-'orangered3'  #  manually change 
B.cl['ETOH 48h']<-'seagreen4'  #  manually change
B.cl[names(B)]<-colorRampPalette(brewer.pal(8, 'Paired'))(length(B))  #  add the uncorrelated/correlated as color pairs
wilcox.test(x=B[['+Tet 4h correlated']], y=B[['ETOH 4h correlated']], alternative='less')$p.value        #  0.9672454597
wilcox.test(x=B[['+Tet 4h uncorrelated']], y=B[['ETOH 4h uncorrelated']], alternative='less')$p.value    #  0.6372411178
wilcox.test(x=B[['+Tet 48h correlated']], y=B[['ETOH 48h correlated']], alternative='less')$p.value      #  0.4296600023
wilcox.test(x=B[['+Tet 48h uncorrelated']], y=B[['ETOH 48h uncorrelated']], alternative='less')$p.value  #  0.002581112858


#  [ETOH all together] ecdfs
b<-B[c('ETOH 4h correlated', 'ETOH 4h uncorrelated', 'ETOH 48h correlated', 'ETOH 48h uncorrelated')]
b.cl<-B.cl[ names(b) ]
my_ecdfs(b, b.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='bottomright', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_ETOH_4h+48h_per_cor_group.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
rm(b, b.cl)


#  [+Tet all together] ecdfs
b<-B[c('+Tet 4h correlated', '+Tet 4h uncorrelated', '+Tet 48h correlated', '+Tet 48h uncorrelated')]
b.cl<-B.cl[ names(b) ]
my_ecdfs(b, b.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='bottomright', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_+Tet_4h+48h_per_cor_group.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
rm(b, b.cl)


#  [4h] ecdfs 
b<-B[c('+Tet 4h correlated', 'ETOH 4h correlated', '+Tet 4h uncorrelated', 'ETOH 4h uncorrelated')]
b.cl<-B.cl[ names(b) ]
my_ecdfs(b, b.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='bottomright', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_+Tet_ETOH_4h_per_cor_group.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
rm(b, b.cl)


#  [48h] ecdfs 
b<-B[c('+Tet 48h correlated', 'ETOH 48h correlated', '+Tet 48h uncorrelated', 'ETOH 48h uncorrelated')]
b.cl<-B.cl[ names(b) ]
my_ecdfs(b, b.cl, XLIM=c(-10, 5), XLAB=expression(log[2]('0.001 + ratio')), MAIN='', LTY=1, LWD=12, LEGEND='bottomright', svg.file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_ratio_circular_vs_linear_junctions_MYCN_Tet-inducible_+Tet_ETOH_48h_per_cor_group.svg', mar=c(6.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
rm(b, b.cl)

#}}}


#  [MYCN Tet-inducible system] look at circARID1A
#{{{

#  split into groups and compute log2(0.001 + ratio)
B<-split(log2(1e-3+lc.tet[ circ_name %in% 'ENSG00000117713.20|ARID1A_chr1+26729651-26732792', ratio]), factor(lc.tet[ circ_name %in% 'ENSG00000117713.20|ARID1A_chr1+26729651-26732792', treatment], levels=unique(lc.tet[ circ_name %in% 'ENSG00000117713.20|ARID1A_chr1+26729651-26732792', treatment])))
B.cl<-unique(meta.tet[, c('treatment', 'col')])
B.cl<-setNames(B.cl$col, B.cl$treatment)
B.cl['+Tet 4h']<-'orange3'  #  manually change 
B.cl['ETOH 4h']<-'palegreen4'  #  manually change 
B.cl['+Tet 48h']<-'orangered3'  #  manually change 
B.cl['ETOH 48h']<-'seagreen4'  #  manually change
wilcox.test(x=B[['+Tet 4h']], y=B[['ETOH 4h']], alternative='less')$p.value    #  0.65
wilcox.test(x=B[['+Tet 48h']], y=B[['ETOH 48h']], alternative='less')$p.value  #  0.35
B<-t(as.matrix(do.call(rbind, B)))


#  matplot
par(mar=c(9.5, 8.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(c(-3, +0.5), 4)
bp<-barplot(B, beside=T, plot=F)
plot(0:1, 0:1, type='n', ylim=range(YTICK), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col=rep(B.cl, each=nrow(B)), axes=F, axisnames=F, beside=T, xlab='', ylab='', las=1, yaxt='n', ylim=range(YTICK), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
axis(2, at=YTICK, line=0, cex.axis=2.8)
mtext(expression(log[2]('0.001 + ratio')), side=2, line=5, padj=+0.2, las=0, cex=2.4)
mtext(text=colnames(B), side=1, line=0, at=colMeans(bp), las=2, adj=1, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_ratio_circular_vs_linear_junctions_circARID1A_Tet-inducible_4h_48h.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}



#  distributions of coefficients of variation and Fano factors for circRNAs and mRNAs  
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  keep only non-failed neuroblastoma tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
nb.meta<-meta.tum[ !(failed) ]
rm(meta.tum, meta.cel, meta.prefailed)


#  compute mRNA CVs and Fano factors based on external linear junction counts 
#  compute circRNA CVs and Fano factors based on the circular junction counts
#
#  N.B. Sequencing depth is an issue but since this affects equally mRNAs and circRNAs the count-based comparisons still make sense. However, 
#       there are crazy cases like:
#
#           ENSG00000184007.21|PTP4A2_chr1-31915895-31919658
#           ENSG00000151779.13|NBAS_chr2-15461201-15504213
#
#       with CV>6 mostly due to CB2001-11-R01 (MNA).
#
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-unlist(lin.cir)[bid %in% nb.meta$bid]
nb<-lin.cir[, .(c.cv=sd(c.count, na.rm=T)/mean(c.count, na.rm=T), c.fano=var(c.count, na.rm=T)/mean(c.count, na.rm=T), l.cv=sd(l.count.out, na.rm=T)/mean(l.count.out, na.rm=T), l.fano=var(l.count.out, na.rm=T)/mean(l.count.out, na.rm=T)), by=.(circ_name, gene_name)]
rm(lin.cir)


#  summarize circRNAs, mRNAs at the gene level by taking the means
nb<-nb[, .(c.cv=mean(c.cv), c.fano=mean(c.fano), l.cv=mean(l.cv), l.fano=mean(l.fano)), by=.(gene_name)]


#  [CV] histogram
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.0, 8.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#
wilcox.test(x=nb[, c.cv], y=nb[, l.cv], alternative='greater')$p.value   #  0
#
XTICK<-pretty(c(0, max(nb[, c.cv], nb[, l.cv], na.rm=T)), 10)
hc<-hist(nb[, c.cv], breaks=seq(min(XTICK), max(XTICK), length.out=50), plot=F)
hl<-hist(nb[, l.cv], breaks=seq(min(XTICK), max(XTICK), length.out=50), plot=F)
YMAX<-max(pretty(c(0, max(hc$counts, hl$counts)), 4))
hc<-hist(nb[, c.cv], breaks=hc$breaks, col='lightblue4', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=c(0, 3), add=F)
hl<-hist(nb[, l.cv], breaks=hl$breaks, col='lightgoldenrod3', border='white', add=T)
mtext('CV', side=1, line=4, padj=-0.1, las=0, cex=2.4)
mtext('Frequency', side=2, line=6, padj=-0.4, las=0, cex=2.4)
legend('top', legend=c('circRNAs', 'mRNAs'), col=c('lightblue4', 'lightgoldenrod3'), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.70, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_CV_circular_vs_linear_junctions.svg', width=12, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  [Fano factors] barplot
x11(width=20, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
wilcox.test(x=nb[, c.fano], y=nb[, l.fano], alternative='less')$p.value   #  0
summary(nb[, c.fano])
#        Min.     1st Qu.      Median        Mean     3rd Qu.        Max.        NA's 
#    0.477311    3.032120    5.041318   14.257541    9.013561 9953.573702           3 
sd(nb[, c.fano], na.rm=T)  #  215.3762654
summary(nb[, l.fano])
#        Min.     1st Qu.      Median        Mean     3rd Qu.        Max. 
#     0.92508    14.04880    26.31533   111.01752    61.36200 46656.64860 
sd(nb[, l.fano], na.rm=T)  #  1298.380388
CUT<-50
STEP<-2
x<-table(cut(nb[, c.fano], breaks=c(seq(0, CUT+STEP, STEP), nb[, max(c.fano, l.fano, na.rm=T)]), labels=c(seq(0, CUT, STEP), paste0('>', CUT)), right=F))
y<-table(cut(nb[, l.fano], breaks=c(seq(0, CUT+STEP, STEP), nb[, max(c.fano, l.fano, na.rm=T)]), labels=c(seq(0, CUT, STEP), paste0('>', CUT)), right=F))
x<-setNames(as.integer(x), names(x))
y<-setNames(as.integer(y), names(y))
#
#  stacked barplot
#
n<-c(setdiff(seq(1, length(x), 5), length(x)-1), length(x))
#YTICK<-pretty(c(0, max(log10(1+x) + log10(1+y))), 5)
YTICK<-pretty(c(0, max(x + y)), 5)
YMAX<-tail(YTICK, 1)
#par(mar=c(5.5, 7.0, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
par(mar=c(5.5, 8.5, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#br<-barplot(log10(1+rbind(x, y)), beside=F, axes=F, plot=F)
br<-barplot(rbind(x, y), beside=F, axes=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(br)+c(-1, 5), axes=F, ann=F, xaxs='i', yaxs='i')
#br<-barplot(log10(1+rbind(x,y)), col=c('lightblue4', 'lightgoldenrod3'), border='white', ylim=c(0, YMAX), axisnames=F, xlab='', ylab='', yaxt='n', xaxt='n', axes=F, add=T)
br<-barplot(rbind(x,y), col=c('lightblue4', 'lightgoldenrod3'), border='white', ylim=c(0, YMAX), axisnames=F, xlab='', ylab='', yaxt='n', xaxt='n', axes=F, add=T)
axis(1, at=br[n], labels=names(x)[n], line=0, tick=T, lwd=0, lwd.ticks=1, tcl=-0.2, padj=0, cex.axis=2.4, xpd=NA)
axis(2, at=YTICK, cex.axis=2.4)
mtext('Fano factor', side=1, line=4, padj=+0.1, las=0, cex=2.4)
#mtext(expression(log[10]('Frequency')), side=2, line=6, padj=+0.8, las=0, cex=2.4)
mtext('Frequency', side=2, line=6, padj=-0.1, las=0, cex=2.4)
legend('top', legend=c('circRNAs', 'mRNAs'), col=c('lightblue4', 'lightgoldenrod3'), bty='n', lty=1, lwd=18, pch=NA, cex=2.0, xpd=T, y.intersp=0.80, x.intersp=0.2, seg.len=0.5)
#
#  overlayed barplots
#
n<-c(setdiff(seq(1, length(x), 5), length(x)-1), length(x))
#YTICK<-pretty(c(0, log10(1+max(x, y))), 5)
YTICK<-pretty(c(0, max(x, y)), 5)
YMAX<-tail(YTICK, 1)
par(mar=c(5.5, 8.5, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
#br<-barplot(log10(1+rbind(x, y)), beside=T, axes=F, plot=F)
br<-barplot(rbind(x, y), beside=T, axes=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(br)+c(-2, 2), axes=F, ann=F, xaxs='i', yaxs='i')
#br<-barplot(log10(1+rbind(x,y)), beside=T, col=c('lightblue4', 'lightgoldenrod3'), border='white', ylim=c(0, YMAX), axisnames=F, xlab='', ylab='', yaxt='n', xaxt='n', axes=F, add=T)
br<-barplot(rbind(x,y), beside=T, col=c('lightblue4', 'lightgoldenrod3'), border='white', ylim=c(0, YMAX), axisnames=F, xlab='', ylab='', yaxt='n', xaxt='n', axes=F, add=T)
axis(1, at=colMeans(br)[n], labels=names(x)[n], line=0, tick=T, lwd=0, lwd.ticks=1, tcl=-0.2, padj=0, cex.axis=2.4, xpd=NA)
axis(2, at=YTICK, cex.axis=2.4)
mtext('Fano factor', side=1, line=4, padj=+0.1, las=0, cex=2.4)
#mtext(expression(log[10]('Frequency')), side=2, line=6, padj=+0.5, las=0, cex=2.4)
mtext('Frequency', side=2, line=6, padj=-0.1, las=0, cex=2.4)
legend('top', legend=c('circRNAs', 'mRNAs'), col=c('lightblue4', 'lightgoldenrod3'), bty='n', lty=1, lwd=18, pch=NA, cex=2.0, xpd=T, y.intersp=0.60, x.intersp=0.4, seg.len=0.3)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_Fano_circular_vs_linear_junctions.svg', width=16, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}




####################################
#
#
#  neuroblastoma specific circRNAs
#
#
####################################




#  [neuroblastoma-specific circRNAs] mRNA expression estimated also from the linear junction counts normalized to CPMs
#                                    significance is estimated by two-sided Mann-Whitney tests
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  functions
#{{{

source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/add_sig_bars.R')

do_plots<-function(circ_name='', alternative='greater', do.plots=T){

    #  isolate the linear junction CPMs of the chosen circRNA
    nb.<-nb.lin[, c('bid', circ_name), with=F]
    colnames(nb.)<-c('bid', 'circ')
    hb.<-hb.lin[, c('bid', circ_name), with=F]
    colnames(hb.)<-c('bid', 'circ')
    vt.<-vt.lin[, c('bid', circ_name), with=F]
    colnames(vt.)<-c('bid', 'circ')
    nb.<-setNames(nb.$circ, nb.$bid)
    hb.<-setNames(hb.$circ, hb.$bid)
    vt.<-setNames(vt.$circ, vt.$bid)


    #  replace bid by tissue throughout
    names(nb.)<-rep('neuroblastoma', length(nb.))
    names(hb.)<-rep('brain tissue', length(hb.))
    names(vt.)<-rep('various tumors', length(vt.))


    #  do Mann-Whitney tests (does also human brain vs various tumors but this is ignored in the analysis)
    B<-list(nb=nb., hb=hb., vt=vt.)
    if(!all(is.na(nb.))){
        PV<-do.call(rbind, apply(combn(names(B), 2), 2, function(s){ data.frame(pv=wilcox.test(x=B[[s[1]]], y=B[[s[2]]], alternative=alternative)$p.value, row.names=paste(s[1], s[2], collapse=' ')) }))
    } else {
        PV<-data.frame(pv=c(NA, NA, NA), row.names=apply(combn(names(B), 2), 2, paste, collapse=' '))
    }
    PV$symbol<-'ns'
    PV$symbol[ PV$pv<1e-10 ]<-'***'
    PV$symbol[ PV$pv>=1e-10 & PV$pv<0.001 ]<-'**'
    PV$symbol[ PV$pv>=0.001 & PV$pv<0.05 ]<-'*'
    rm(B)

    
    if(do.plots){
        if(!all(is.na(nb.))){
            #  barplot of ordered CPMs 
            x11(width=40, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
            B<-c( sort(nb., decreasing=T, na.last=T), 
                  sort(hb., decreasing=T, na.last=T), 
                  sort(vt., decreasing=T, na.last=T) )
            B.cl<-c( setNames(rep('seagreen4', length(nb.)), names(nb.)), 
                     setNames(rep('cornflowerblue', length(hb.)), names(hb.)), 
                     setNames(rep('coral4', length(vt.)), names(vt.)) )
            par(mar=c(1.0, 5.5, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
            YTICK<-pretty(c(0, max(B, na.rm=T)), 5)
            bp<-barplot(B, beside=F, plot=F)
            plot(0:1, 0:1, type='n', ylim=c(0, tail(YTICK, 1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
            bp<-barplot(B, border='white', col=B.cl, axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
            axis(2, at=YTICK, line=0, cex.axis=2.4)
            legend('topright', legend=c('neuroblastoma', 'human brain tissue', 'various tumors'), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.70, x.intersp=0.2, seg.len=0.5)
            dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_across_tissues_linear_junction_CPMs_for_', circ_name, '.svg'), width=40, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
            dev.off()


            #  boxplot of CPMs
            options(scipen=+20)
            x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
            par(mar=c(5.0, 8.5, 2.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
            B<-list(nb=nb., hb=hb., vt=vt.)
            B.cl<-c('seagreen4', 'cornflowerblue', 'coral4')
            YTICK<-pretty(c(0, max(sapply(B, max, na.rm=T))), 4)
            plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
            bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
            YTICK<-pretty(c(0, 1.2*max(bp$stats[5, ])), 4)  #  recompute y-max to leave room for significance bars
            plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
            bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
            mtext('Counts per million', side=2, line=6, padj=-0.1, las=0, cex=2.4)
            mtext(text=c('neuroblastoma', 'brain\ntissue', 'various\ntumors'), side=1, line=-1, at=seq_len(length(B)), las=1, padj=1, cex=2.4, col=B.cl)
            add_sig_bars(bp, PV)
            dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_tissues_linear_junction_CPMs_for_', circ_name, '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
            options(scipen=0)
            dev.off()
        } else {
            cat('\n\n**** No mRNA expression for:', circ_name, '****\n\n')
        }
    }


    return(list(nb.=nb., hb.=hb., vt.=vt., pv=PV))
}

#}}}


#  load neuroblastoma-specific circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_neuroblastoma-specific.RData')


#  [neuroblastoma tumors] compute mRNA CPMs based on external linear junction counts and total counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-unlist(lin.cir)[bid %in% colnames(nb.circs)]
lin.cir<-lin.cir[ circ_name %in% rownames(nb.circs) ]
lin.cir<-lin.cir[, cpm:=l.count.out/total*1e6]
nb.lin<-dcast(lin.cir, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)
nb.lin<-nb.lin[ match(colnames(nb.circs), bid) ]
rm(lin.cir,nb.meta)


#  [human brain tissue] compute mRNA CPMs based on external linear junction counts and total counts
load('/fast/projects/Schulte_NB/work/downloads/paired-end_totalRNA_human_brain/neuroblastoma_circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-unlist(lin.cir)[ bid %in% colnames(hb.circs) ]  #  unnecessary since we keep all samples but good for consistency
lin.cir<-lin.cir[ circ_name %in% rownames(nb.circs) ]
lin.cir<-lin.cir[, cpm:=l.count.out/total*1e6]
hb.lin<-dcast(lin.cir, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)
hb.lin<-hb.lin[ match(colnames(hb.circs), bid) ]
rm(lin.cir,hb.meta)


#  [various tumors] compute mRNA CPMs based on external linear junction counts and total counts
#                   remove the two libraries with 52nt reads that are not included in the linear vs circular junction analysis
load('/fast/projects/Schulte_NB/work/downloads/paired-end_totalRNA_cancer/neuroblastoma_circRNAs_linear_vs_circular_collected_results.RData')
vt.circs<-vt.circs[, grep('SRR1617643|SRR1617644', colnames(vt.circs), invert=T)]
lin.cir<-unlist(lin.cir)[ bid %in% vt.meta[, unlist(srr)] ]
lin.cir<-lin.cir[ circ_name %in% rownames(nb.circs) ]
lin.cir<-lin.cir[, cpm:=l.count.out/total*1e6]
vt.lin<-dcast(lin.cir, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)
vt.lin<-vt.lin[ match(colnames(vt.circs), bid) ]
rm(lin.cir,vt.meta)


#  do all plots
#  NA values result in NA p-values
res<-List()
for (circ_name in colnames(nb.lin[, -1])){
    cat('\n', circ_name, ':\n', sep='')
    res[[circ_name]]<-do_plots(circ_name=circ_name, alternative='two.sided', do.plots=T)
    cat('\n')
}


#  save
save(res, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/linear_junction_CPMs_of_circRNAs_across_tissues_neuroblastoma-specific.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/linear_junction_CPMs_of_circRNAs_across_tissues_neuroblastoma-specific.RData



#  [neuroblastoma-specific circRNAs] plots including annotations for the mRNA expression computed above
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(grid)
library(viridis)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/add_sig_bars.R')
source('~/bio/lib/draw_highlights.R')


#  functions
#{{{
boxplot_per_risk_group<-function(X, META){
    colnames(X)<-META[ match(colnames(X), bid), risk_group ]
    B<-split(log10(1+X), colnames(X))[META[, unique(risk_group)]]
    B.cl<-META[ risk_group %in% names(B), unique(col) ] 
    YTICK<-pretty(c(0.0, sapply(B, max)), 5)
    plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
    bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
    mtext(expression(log[10](1+'CPM')), side=2, line=5, padj=+0.1, las=0, cex=2.4)
    #mtext(text=paste0(names(B), ' (', lengths(B)/nrow(X), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
    mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
}
#}}}


# load the neuroblastoma specific circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_neuroblastoma-specific.RData')
rm(list=setdiff(ls(), ls(pattern='circs|CIRCS|_|nb.meta')))  #  remove clutter


#  load the mRNA expression estimated from linear junction counts external to the corresponding circRNA and normalized to CPMs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/linear_junction_CPMs_of_circRNAs_across_tissues_neuroblastoma-specific.RData')


#  compare mean expressions to identify the consistently significantly up- or down- regulated mRNAs across both pairwise comparisons
up<-sapply(res, function(r){ 
    if(any(is.na(r$pv$pv))){ return(NA) }
    if(all(r$pv[c('nb hb', 'nb vt'), 'pv']<0.05)){
        if( mean(r$nb., na.rm=T)>mean(r$hb., na.rm=T) & mean(r$nb., na.rm=T)>mean(r$vt., na.rm=T) ){
            return(T)
        } else {
            return(F)
        }
    } else {
        return(F)
    }
})
down<-sapply(res, function(r){ 
    if(any(is.na(r$pv$pv))){ return(NA) }
    if(all(r$pv[c('nb hb', 'nb vt'), 'pv']<0.05)){
        if( mean(r$nb., na.rm=T)<mean(r$hb., na.rm=T) & mean(r$nb., na.rm=T)<mean(r$vt., na.rm=T) ){
            return(T)
        } else {
            return(F)
        }
    } else {
        return(F)
    }
})
down<-down[ rownames(nb.circs) ]
up_down<-up[ rownames(nb.circs) ]
up_down[down]<- -1
up_down[is.na(up_down)]<-0
rm(up,down)


#  order them by NB expression
n<-order(rowMeans(nb.circs), decreasing=T)
nb.circs<-nb.circs[n, ]
hb.circs<-hb.circs[n, ]
vt.circs<-vt.circs[n, ]
CIRCS<-CIRCS[n]
rm(n)


#  recycle
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  heatmap on the log10-transformed CPMs (best contrast, compared to z-scores of raw CPMs, or z-scores of log10-transformed CPMs)
X<-cbind(log10(1+nb.circs), log10(1+hb.circs), log10(1+vt.circs))
col.an<-data.frame(tissue=factor(rep(c('neuroblastoma', 'brain tissue', 'various tumors'), c(ncol(nb.circs), ncol(hb.circs), ncol(vt.circs))), exclude=F), row.names=colnames(X))
row.an<-data.frame(mRNA=factor(up_down, levels=c(-1, 1, 0), labels=c('down', 'up', 'neither')), row.names=rownames(X))
cl<-setNames(list(setNames( c('seagreen4', 'cornflowerblue', 'coral4'), c('neuroblastoma', 'brain tissue', 'various tumors') ),
                  setNames( c('lightsteelblue3', 'grey20', 'grey49'), c('down', 'up', 'neither') )), 
             c(colnames(col.an), colnames(row.an)))
ph<-pheatmap(X, color=viridis(10), border_color=NA, scale='none', 
        cluster_rows=F,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=col.an, annotation_row=row.an, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=F, 
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_circRNAs_across_tissues_NB-specific.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(X,col.an,row.an,cl,ph)


#  boxplots of log10-transformed CPMs
#  order them by NB expression
#  use circBase circ_id but replace hsa_circ with gene_name
options(scipen=0)
par(mar=c(15.0, 8.0, 0.5, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
stopifnot( all.equal(rownames(nb.circs), CIRCS$circ_name) )
stopifnot( all.equal(rownames(nb.circs), rownames(hb.circs)) )
stopifnot( all.equal(rownames(nb.circs), rownames(vt.circs)) )
NAMES<-sapply(apply(data.frame(mcols(CIRCS)[, c('gene_name', 'circ_id')]), 1, function(n){ sub('hsa_circ', paste0('circ', n[[1]]), n[[2]]) }) , paste, collapse=',\n')
B<-c(split(log10(1+nb.circs), row(nb.circs)), split(log10(1+hb.circs), row(hb.circs)), split(log10(1+vt.circs), row(vt.circs)))
n<-c(rbind(rbind(1:nrow(nb.circs), nrow(nb.circs)+(1:nrow(nb.circs)), 2*nrow(nb.circs)+(1:nrow(nb.circs)))))  #  order them in triplets
B<-B[n]
B.cl<-rep(c('seagreen4', 'cornflowerblue', 'coral4'), nrow(nb.circs))
m.cl<-up_down
m.cl[ m.cl==-1 ]<-'lightsteelblue3'
m.cl[ m.cl==1 ]<-'grey20'
m.cl[ m.cl==0 ]<-'grey49'
YTICK<-pretty(c(0, max(sapply(B, max, na.rm=T))), 4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'CPM')), side=2, line=4, padj=-0.2, las=0, cex=2.4)
mtext(text=sub('_.*$', '', NAMES), side=1, line=0, at= seq(2, length(B), 3), las=2, adj=1, cex=2.4, col=m.cl)  #  remove circBase id
draw_highlights(L=length(B), STEP=3, YMAX=max(YTICK))
legend('top', legend=c('neuroblastoma', 'brain tissue', 'various tumors'), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.40, x.intersp=0.3, seg.len=0.3)
legend('topright', legend=c('down', 'up', 'neither'), col=c('lightsteelblue3', 'grey20', 'grey49'), title='mRNA', bty='n', lty=1, lwd=15, pch=NA, cex=1.5, xpd=T, y.intersp=0.50, x.intersp=0.3, seg.len=0.3)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_across_tissues_NB-specific.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(YTICK,NAMES,B,n,B.cl,bp)


#  boxplots of log10-transformed CPMs per circRNA across risk groups
#{{{

#  plot all of them 
stopifnot( all.equal(rownames(nb.circs), CIRCS$circ_name) )
NAMES<-setNames(sapply(apply(data.frame(mcols(CIRCS)[, c('gene_name', 'circ_id')]), 1, function(n){ sub('hsa_circ', paste0('circ', n[[1]]), n[[2]]) }) , paste, collapse=',\n'), CIRCS$circ_name)
par(mar=c(10.0, 8.0, 1.0, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
for(p in rownames(nb.circs)){
    cat('plotting:', p, '\n')
    boxplot_per_risk_group(X=nb.circs[p, , drop=F], META=nb.meta)
    mtext(text=sub('_.*$', '', NAMES[p]), side=3, line=-1, las=1, padj=-0.1, cex=1.8)
    dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_', p, '_across_risk-groups.svg'), width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
    readline()
}

#}}}

#}}}



########################
#
#
#  scrap code
#  fast analyses/plots
#
#
########################




#  [MNA vs HR_nMNA DE] are the correlated/uncorrelated circRNA/mRNA pairs enriched in DE genes?
#                      
#                      N.B. I am not sure what this analysis tells us...
#
#{{{
library(DESeq2)

#  load DE results
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/DESeq2_genes_MNA_HR_nMNA.RData')


#  correlations across all samples 
#{{{

#  isolate correlated/uncorrelated groups
R.cor<-names(cl.cor[['all']][ cl.cor[['all']]>=0.6 ])
R.uncor<-names(cl.cor[['all']][ cl.cor[['all']]<0.6 ])


#  [DEPLETED in significantly up-DE genes with baseMean>100] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#     correlated    |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#   not correlated  |      x2          |          y2          | 
grp<-subset(RES, baseMean>100 & log2FoldChange>0 & padj<0.05)$gene_name
rst<-setdiff( subset(RES, baseMean>100)$gene_name, grp )
fisher.test(data.frame('in'=c(length(intersect(grp, R.cor)), 
                              length(setdiff(grp, R.cor))), 
                       'out'=c(length(intersect(rst, R.cor)), 
                               length(setdiff(rst, R.cor)))), alternative='less')$p.value  #  depleted
#
#  => 0.007311068114


#  [ENRICHED in significantly down-DE genes with baseMean>100] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#     correlated    |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#   not correlated  |      x2          |          y2          | 
grp<-subset(RES, baseMean>100 & log2FoldChange<0 & padj<0.05)$gene_name
rst<-setdiff( subset(RES, baseMean>100)$gene_name, grp )
fisher.test(data.frame('in'=c(length(intersect(grp, R.cor)), 
                              length(setdiff(grp, R.cor))), 
                       'out'=c(length(intersect(rst, R.cor)), 
                               length(setdiff(rst, R.cor)))), alternative='greater')$p.value  #  enriched
#
#  => 0.00005333347341


#  [DEPLETED in significantly up-DE genes with baseMean>100] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#   uncorrelated    |      x1          |          y1          |
#  -----------------|------------------|----------------------|
# not uncorrelated  |      x2          |          y2          | 
grp<-subset(RES, baseMean>100 & log2FoldChange>0 & padj<0.05)$gene_name
rst<-setdiff( subset(RES, baseMean>100)$gene_name, grp )
fisher.test(data.frame('in'=c(length(intersect(grp, R.uncor)), 
                              length(setdiff(grp, R.uncor))), 
                       'out'=c(length(intersect(rst, R.uncor)), 
                               length(setdiff(rst, R.uncor)))), alternative='less')$p.value  #  depleted
#
#  => 0.00001382068437


#  [ENRICHED in significantly down-DE genes with baseMean>100] Fisher's exact test:
#
#                   |      DE          |        not DE        |
#  -----------------|------------------|----------------------|
#   uncorrelated    |      x1          |          y1          |
#  -----------------|------------------|----------------------|
# not uncorrelated  |      x2          |          y2          | 
grp<-subset(RES, baseMean>100 & log2FoldChange<0 & padj<0.05)$gene_name
rst<-setdiff( subset(RES, baseMean>100)$gene_name, grp )
fisher.test(data.frame('in'=c(length(intersect(grp, R.uncor)), 
                              length(setdiff(grp, R.uncor))), 
                       'out'=c(length(intersect(rst, R.uncor)), 
                               length(setdiff(rst, R.uncor)))), alternative='greater')$p.value  #  enriched
#
#  => 0.00001884810333

#}}}

#}}}



#  [ARID1A mRNA, kallisto totalRNA-seq + polyA-seq] dominant isoform resolution
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/grouplist2boxplot.R')


#  keep only non-failed neuroblastoma tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
nb.meta<-meta.tum[ !(failed) ]
rm(meta.tum, meta.cel, meta.prefailed)


#  load the unified set of circRNAs and isolate circARID1A
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
int<-CIRCS[ grep('ARID1A', CIRCS$circ_name) ]
nb.circs<-nb.circs[ nb.circs$bid %in% nb.meta$bid & nb.circs$circ_name %in% int$circ_name ]
hb.circs<-hb.circs[ hb.circs$circ_name %in% int$circ_name ]
vt.circs<-vt.circs[ vt.circs$circ_name %in% int$circ_name ]
CIRCS<-int
rm(CIRCS.all, GENES, nb.genes, hb.genes, vt.genes, int)


#  remove the two libraries with 52nt reads that are not included in the linear vs circular junction analysis
vt.circs<-vt.circs[ grep('SRR1617643|SRR1617644', vt.circs$bid, invert=T) ]


#  load reference and identify all ARID1A isoforms
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'transcript' ]
mcols(hsa)<-mcols(hsa)[, c('gene_name', 'gene_id', 'gene_type', 'transcript_name', 'transcript_id', 'transcript_type')]
hsa<-hsa[ grep('ARID1A', hsa$gene_name) ]


#  load ARID1A isoform quantification by kallisto for both totalRNA-seq and polyA-seq
#  order isoforms by mean TPM expression
#{{{

#  totalRNA-seq
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/kallisto.RData')
tra<-lapply(totalrna[ nb.meta$bid ], function(d){ as.data.table(d)[ transcript_id %in% hsa$transcript_id ] } )
invisible(lapply(names(tra), function(n){ tra[[n]][, bid:=n] }))          #  add bid columns
tra<-do.call(rbind, tra)
tra<-dcast(tra, bid ~ transcript_id, value.var='tpm', fun.aggregate=sum)  #  data.table of TPMs
tra<-tra[ match(nb.meta$bid, bid), ]                                     #  order the samples again according to the metadata
stopifnot( all.equal( tra$bid, nb.meta$bid ) )
tra<-tra[, bid:=nb.meta$risk_group ]                                      #  replace bid column with risk_group
colnames(tra)<-sub('bid', 'risk_group', colnames(tra))
n<-sort(colMeans(tra[, -1]), decreasing=T)
tra<-tra[, c('risk_group', names(n)), with=F]
rm(n,totalrna)
gc()


#  tumor polyA-seq
load('/fast/projects/peifer_wgs/work/work/2018-05-23_Fuchs_polyA/raw/metadata.RData')
nb.meta.polya<-meta[ !is.na(risk_group) ]
load('/fast/projects/peifer_wgs/work/work/2018-05-23_Fuchs_polyA/raw/kallisto.RData')
tra.polya<-lapply(polya[ nb.meta.polya$bid ], function(d){ as.data.table(d)[ transcript_id %in% hsa$transcript_id ] } )
invisible(lapply(names(tra.polya), function(n){ tra.polya[[n]][, bid:=n] }))          #  add bid columns
tra.polya<-do.call(rbind, tra.polya)
tra.polya<-dcast(tra.polya, bid ~ transcript_id, value.var='tpm', fun.aggregate=sum)  #  data.table of TPMs
tra.polya<-tra.polya[ match(nb.meta.polya$bid, bid), ]                                #  order the samples again according to the metadata
stopifnot( all.equal( tra.polya$bid, nb.meta.polya$bid ) )
tra.polya<-tra.polya[, bid:=nb.meta.polya$risk_group ]                                #  replace bid column with risk_group
colnames(tra.polya)<-sub('bid', 'risk_group', colnames(tra.polya))
n<-sort(colMeans(tra.polya[, -1]), decreasing=T)
tra.polya<-tra.polya[, c('risk_group', names(n)), with=F]
rm(n,meta,polya)
gc()

#}}}


#  color-indicate protein-coding isoforms
C.cl<-setNames(hsa$transcript_type, hsa$transcript_name)
C.cl[ grep('protein_coding', C.cl) ]<-'firebrick1'
C.cl[ grep('firebrick1', C.cl, invert=T) ]<-'black'


#  boxplots of isoform expression across all risk groups
#{{{

#  totalRNA-seq
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
B<-as.data.frame(tra[, -1])
CL<-C.cl[ hsa$transcript_name[ match( colnames(tra)[-1], hsa$transcript_id ) ] ]
par(mar=c(10.5,5.5,0.5,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
YTICK<-pretty(c(0, ceiling(max(B))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col='grey39', ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol='grey39', range=0, add=T)
mtext('TPM', side=2, line=3, padj=-0.4, las=0, cex=2.0)
mtext(text=hsa$transcript_name[ match(names(B), hsa$transcript_id) ], side=1, line=0, at=seq_along(B), las=2, adj=0.99, cex=2.0, col=CL)
legend('topleft', legend='protein-coding', col=unique(CL[ grep('black', CL, invert=T)]), bty='n', lty=1, lwd=10, pch=NA, cex=1.8, y.intersp=0.8, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_ARID1A_isoforms_across_tumors_kallisto.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B,CL)
dev.off()


#  polyA-seq
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
B<-as.data.frame(tra.polya[, -1])
CL<-C.cl[ hsa$transcript_name[ match( colnames(tra.polya)[-1], hsa$transcript_id ) ] ]
par(mar=c(10.5,5.5,0.5,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
YTICK<-pretty(c(0, ceiling(max(B))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col='grey39', ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol='grey39', range=0, add=T)
mtext('TPM', side=2, line=3, padj=-0.4, las=0, cex=2.0)
mtext(text=hsa$transcript_name[ match(names(B), hsa$transcript_id) ], side=1, line=0, at=seq_along(B), las=2, adj=0.99, cex=2.0, col=CL)
legend('topleft', legend='protein-coding', col=unique(CL[ grep('black', CL, invert=T)]), bty='n', lty=1, lwd=10, pch=NA, cex=1.8, y.intersp=0.8, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2018-05-23_Fuchs_polyA/figures/boxplot_ARID1A_isoforms_across_tumors_kallisto.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B,CL)
dev.off()

#}}}


#  boxplots of isoform expression per risk group
#{{{

#  totalRNA-seq
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
B<-split(tra, tra$risk_group)
CL<-C.cl[ hsa$transcript_name[ match( colnames(tra)[-1], hsa$transcript_id ) ] ]
invisible(lapply(names(B), function(n){ B[[n]][, 'risk_group':=NULL,] } ))
B<-lapply(B, function(b){ b<-t(as.data.frame(b)); rownames(b)<-hsa$transcript_name[ match(rownames(b), hsa$transcript_id) ]; return(b); })
B<-B[ c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA') ]
B.cl<-setNames(nb.meta[, unique(col)], nb.meta[, unique(risk_group)])[na.omit(names(B))]
grouplist2boxplot(L=B, L.COL=B.cl, YLAB='TPM', XLAB.COL=CL, mar=c(11.0,5.5,0.5,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
legend('topleft', legend='protein-coding', col=unique(CL[ grep('black', CL, invert=T)]), bty='n', lty=1, lwd=10, pch=NA, cex=1.8, y.intersp=0.8, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_ARID1A_isoforms_across_tumors_risk_group-resolved_kallisto.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B,CL)
dev.off()


#  polyA-seq
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
B<-split(tra.polya, tra.polya$risk_group)
CL<-C.cl[ hsa$transcript_name[ match( colnames(tra.polya)[-1], hsa$transcript_id ) ] ]
invisible(lapply(names(B), function(n){ B[[n]][, 'risk_group':=NULL,] } ))
B<-lapply(B, function(b){ b<-t(as.data.frame(b)); rownames(b)<-hsa$transcript_name[ match(rownames(b), hsa$transcript_id) ]; return(b); })
B<-B[ c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA') ]
B.cl<-setNames(nb.meta.polya[, unique(col)], nb.meta.polya[, unique(risk_group)])[na.omit(names(B))]
grouplist2boxplot(L=B, L.COL=B.cl, YLAB='TPM', XLAB.COL=CL, mar=c(11.0,5.5,0.5,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
legend('topleft', legend='protein-coding', col=unique(CL[ grep('black', CL, invert=T)]), bty='n', lty=1, lwd=10, pch=NA, cex=1.8, y.intersp=0.8, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2018-05-23_Fuchs_polyA/figures/boxplot_ARID1A_isoforms_across_tumors_risk_group-resolved_kallisto.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(B,CL)
dev.off()

#}}}

#}}}



#  [ARID1A isoform resolution] dominant isoform resolution using linear junction counts only
#                              we compute the CPMs based on linear junction counts for all the ARID1A linear junctions that entered the analysis
#                              we remove the linear junctions shared with the circARID1A
#                              we resolve junction-membership of the ARID1A isoforms and average over the corresponding CPMs
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('~/bio/lib/grouplist2boxplot.R')


#  functions
#{{{

linear_counts_results<-function(B='', BID=''){
    require(data.table)

    #  try to load the counts
    n<-tryCatch({
        cn<-fread(B, header=F, sep='\t', col.names=c('tx_name', 'count')) 
        total<-cn[, sum(count)]


        #  keep the linear junction counts and compute CPMs
        cn<-cn[ grep('\\|jn_', tx_name, invert=F) ]
        colnames(cn)<-c('junction_name', 'count')
        cn[, cpm:=1e6*count/total]
        rm(total)


        #  add gene_id 
        cn<-cn[ LINEAR, on='junction_name' ]


        #  add the bid 
        cn$bid<-BID


        return(cn)

    }, error=function(e){
        warning(e)
        return(data.table())
    }) 
}

#}}}


#  locate the results
lin_cir<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/ -maxdepth 3 -type f -wholename \'*/lin_vs_circ/counts.tsv\' -print', stdout=T)


#  add sequencing-sample-id (bid)
names(lin_cir)<-sub('^.*/(CB[^/]+)/.*$', '\\1', lin_cir)


#  order them
lin_cir<-lin_cir[ order(names(lin_cir)) ]


#  load metadata to identify failed samples and exclude them from the list
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
nb.meta<-meta.tum[ !(failed) ]
lin_cir<-lin_cir[ names(lin_cir) %in% nb.meta$bid ]
rm(meta.tum, meta.cel, meta.prefailed)


#  load circular and linear junctions and isolate those from ARID1A
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular.RData')
CIRCS<-CIRCS[ CIRCS$gene_name %in% 'ARID1A' ]
CIRCS.junctions<-CIRCS.junctions[ grep(CIRCS$gene_id, names(CIRCS.junctions)) ]
LINEAR<-LINEAR[ gene_id %in% CIRCS$gene_id ]
LINEAR.junctions<-LINEAR.junctions[ LINEAR[, junction_name] ]
rm(CIRCS.junctions.seqs, LINEAR.junctions.seqs)


#  remove the linear junctions that are common with circARID1A
L<-unlist(LINEAR.junctions)
o<-findOverlaps(L, CIRCS, type='within', select='all')
o<-names(L[ queryHits(o) ])
o<-o[ duplicated(o) ]  #  valid junctions are those with both ends present
LINEAR.junctions<-LINEAR.junctions[ setdiff(names(LINEAR.junctions), o) ]
rm(o,L)


#  collect results
lin.cir<-List()
for (n in seq_along(lin_cir)){
    cat('\nprocessing: ', lin_cir[n], '\n')
    lin.cir[[ names(lin_cir)[n] ]]<-linear_counts_results(lin_cir[n], names(lin_cir)[n])  #  "returning -Inf" warnings expected for unexpressed junctions
}
lin.cir<-do.call(rbind, lin.cir)
rm(n,lin_cir,LINEAR)


#  load reference 
#  identify all ARID1A exons
#  group them by isoform
exn<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
exn<-exn[ exn$type %in% 'exon' ]
mcols(exn)<-mcols(exn)[, c('gene_name', 'gene_id', 'gene_type', 'transcript_name', 'transcript_id', 'transcript_type', 'exon_id', 'exon_number')]
exn<-exn[ grep('ARID1A', exn$gene_name) ]
exn<-split(exn, exn$transcript_name)


#  alocate exon-exon junctions to transcripts
#
#  N.B. I have checked this code only for + strand genes...
#
L<-unlist(LINEAR.junctions)
tra2j<-setNames(vector('list', length(exn)), names(exn))
for (n in seq_along(exn)){
    o<-findOverlaps(L, exn[[n]], type='within', select='all')
    o<-o[ start(L[queryHits(o)])==start(exn[[n]][subjectHits(o)]) | end(L[queryHits(o)])==end(exn[[n]][subjectHits(o)]) ]
    o<-names(L[ queryHits(o) ])
    tra2j[[n]]<-o[ duplicated(o) ]  #  count a junction only when both parts are present
}
rm(L,o,n)


#  query the junctions per sample and take the mean over their CPMs for each isoform
#
#  N.B. mean is more robust than sum and more forgiving for mistakes, such as when one shared junction drives the expression up in a transcript otherwise 
#       poorly expressed according to the rest of its junctions, or when annotated junctions differ by a single nucleotide and their CPMs are included 
#       separately in the counting
tra<-data.table(bid=nb.meta$bid)
for (n in seq_along(tra2j)){
    x<-lin.cir[ junction_name %in% tra2j[[n]] ][, .(cpm=mean(cpm, na.rm=T)), by=.(bid)][ !is.na(cpm) ]  #  returns NA if all values are NA so filter out
    colnames(x)<-c('bid', names(tra2j)[n])
    if(nrow(x)>0 ){
        tra<-x[tra, on='bid']  #  we are still going to get NA's here if samples are missing 
    }
}
rm(x,n)


#  order samples by metadata and replace bid by risk group 
tra<-tra[ match(nb.meta$bid, bid), ]
stopifnot( all.equal( tra$bid, nb.meta$bid ) )
tra<-tra[, bid:=nb.meta$risk_group ]
colnames(tra)<-sub('bid', 'risk_group', colnames(tra))


#  order isoforms by mean CPM expression
n<-sort(colMeans(tra[, -1], na.rm=T), decreasing=T)
tra<-tra[, c('risk_group', names(n)), with=F]
rm(n)


#  color-indicate protein-coding isoforms
C.cl<-sapply(exn, function(x){ x$transcript_type[1] })[ colnames(tra)[-1] ]
C.cl[ grep('protein_coding', C.cl) ]<-'firebrick1'
C.cl[ grep('firebrick1', C.cl, invert=T) ]<-'black'


#  boxplot of isoform expression across all risk groups
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
B<-as.data.frame(tra[, -1])
par(mar=c(10.5,6.0,0.5,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
YTICK<-pretty(c(0, ceiling(max(B, na.rm=T))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col='grey39', ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol='grey39', range=0, add=T)
mtext('CPM', side=2, line=3, padj=-0.9, las=0, cex=2.0)
mtext(text=names(B), side=1, line=0, at=seq_along(B), las=2, adj=0.99, cex=2.0, col=C.cl)
legend('topleft', legend='protein-coding', col=unique(C.cl[ grep('black', C.cl, invert=T)]), bty='n', lty=1, lwd=10, pch=NA, cex=1.8, y.intersp=0.8, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_ARID1A_isoforms_across_tumors_linear_vs_circular.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  boxplot of isoform expression per risk group
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
B<-split(tra, tra$risk_group)
invisible(lapply(names(B), function(n){ B[[n]][, 'risk_group':=NULL,] } ))
B<-lapply(B, function(b){ t(as.data.frame(b)) })
B<-B[ c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA') ]
B.cl<-setNames(nb.meta[, unique(col)], nb.meta[, unique(risk_group)])[names(B)]
grouplist2boxplot(L=B, L.COL=B.cl, YLAB='CPM', XLAB.COL=C.cl, mar=c(11.0,6.0,0.5,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.0, cex.axis=2.0)
legend('topleft', legend='protein-coding', col=unique(C.cl[ grep('black', C.cl, invert=T)]), bty='n', lty=1, lwd=10, pch=NA, cex=1.8, y.intersp=0.8, xpd=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_ARID1A_isoforms_across_tumors_risk_group-resolved_linear_vs_circular.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}



#  [circSLC45A4, circARID1A, across tissues] expression comparison of circRNAs of interest (+ their mRNA)
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  functions
#{{{

source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/add_sig_bars.R')

do_plots<-function(circ_name='', d=NULL, do.plots=T){
    #  replace bid by tissue throughout
    names(d$circ$nb)<-rep('neuroblastoma', length(d$circ$nb))
    names(d$circ$hb)<-rep('brain tissue', length(d$circ$hb))
    names(d$circ$vt)<-rep('various tumors', length(d$circ$vt))
    names(d$lin$nb)<-rep('neuroblastoma', length(d$lin$nb))
    names(d$lin$hb)<-rep('brain tissue', length(d$lin$hb))
    names(d$lin$vt)<-rep('various tumors', length(d$lin$vt))


    #  do the Mann-Whitney tests (also does human brain vs various tumors but this is ignored in the analysis)
    B<-list(nb=d$circ$nb, hb=d$circ$hb, vt=d$circ$vt)
    if(!all(is.na(d$circ$nb))){
        PV<-do.call(rbind, apply(combn(names(B), 2), 2, function(s){ data.frame(pv=wilcox.test(x=B[[s[1]]], y=B[[s[2]]], alternative=d$circ$alt)$p.value, row.names=paste(s[1], s[2], collapse=' ')) }))
    } else {
        PV<-data.frame(pv=c(NA, NA, NA), row.names=apply(combn(names(B), 2), 2, paste, collapse=' '))
    }
    PV$symbol<-'ns'
    PV$symbol[ PV$pv<1e-10 ]<-'***'
    PV$symbol[ PV$pv>=1e-10 & PV$pv<0.001 ]<-'**'
    PV$symbol[ PV$pv>=0.001 & PV$pv<0.05 ]<-'*'
    pv.circ<-PV
    B<-list(nb=d$lin$nb, hb=d$lin$hb, vt=d$lin$vt)
    if(!all(is.na(d$lin$nb))){
        PV<-do.call(rbind, apply(combn(names(B), 2), 2, function(s){ data.frame(pv=wilcox.test(x=B[[s[1]]], y=B[[s[2]]], alternative=d$lin$alt)$p.value, row.names=paste(s[1], s[2], collapse=' ')) }))
    } else {
        PV<-data.frame(pv=c(NA, NA, NA), row.names=apply(combn(names(B), 2), 2, paste, collapse=' '))
    }
    PV$symbol<-'ns'
    PV$symbol[ PV$pv<1e-10 ]<-'***'
    PV$symbol[ PV$pv>=1e-10 & PV$pv<0.001 ]<-'**'
    PV$symbol[ PV$pv>=0.001 & PV$pv<0.05 ]<-'*'
    pv.lin<-PV
    rm(B,PV)

    
    if(do.plots){
        #  circRNA boxplot
        if(!all(is.na(d$circ$nb))){
            options(scipen=+20)
            x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
            par(mar=c(5.0, 6.5, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
            B<-list(nb=d$circ$nb, hb=d$circ$hb, vt=d$circ$vt)
            B.cl<-c('seagreen4', 'cornflowerblue', 'coral4')
            YTICK<-pretty(c(0, max(sapply(B, max, na.rm=T))), 4) 
            plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
            bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
            YTICK<-pretty(c(0, 1.2*max(bp$stats[5, ])), 4)  #  recompute y-max to leave room for significance bars
            plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
            bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
            mtext('Counts per million', side=2, line=4, padj=-0.4, las=0, cex=2.4)
            mtext(text=c('neuroblastoma', 'brain\ntissue', 'various\ntumors'), side=1, line=-1, at=seq_len(length(B)), las=1, padj=1, cex=2.4, col=B.cl)
            add_sig_bars(bp, pv.circ)
            dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_tissues_circRNA_CPMs_for_', circ_name, '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
            options(scipen=0)
            dev.off()
        } else {
            cat('\n\n**** No neuroblastoma circRNA expression for:', circ_name, '****\n\n')
        }


        #  mRNA boxplot
        if(!all(is.na(d$lin$nb))){
            options(scipen=+20)
            x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
            par(mar=c(5.5, 8.5, 2.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
            B<-list(nb=d$lin$nb, hb=d$lin$hb, vt=d$lin$vt)
            B.cl<-c('seagreen4', 'cornflowerblue', 'coral4')
            YTICK<-pretty(c(0, max(sapply(B, max, na.rm=T))), 4) 
            plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
            bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
            YTICK<-pretty(c(0, 1.2*max(bp$stats[5, ])), 4)  #  recompute y-max to leave room for significance bars
            plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
            bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
            mtext('Counts per million', side=2, line=6, padj=-0.4, las=0, cex=2.4)
            mtext(text=c('neuroblastoma', 'brain\ntissue', 'various\ntumors'), side=1, line=-1, at=seq_len(length(B)), las=1, padj=1, cex=2.4, col=B.cl)
            add_sig_bars(bp, pv.lin)
            dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_tissues_linear_junction_CPMs_for_', circ_name, '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
            options(scipen=0)
            dev.off()
        } else {
            cat('\n\n**** No neuroblastoma mRNA expression for:', circ_name, '****\n\n')
        }
    }


    return(list(circ=list(nb=d$circ$nb, hb=d$circ$hb, vt=d$circ$vt, pv=pv.circ), lin=list(nb=d$lin$nb, hb=d$lin$hb, vt=d$lin$vt, pv=pv.lin)))
}

#}}}


#  keep only non-failed neuroblastoma tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
nb.meta<-meta.tum[ !(failed) ]
rm(meta.tum, meta.cel, meta.prefailed)


#  load the unified set of circRNAs and isolate the interesting ones and their CPMs based on CIRI2 across tissues
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
CIRCS<-CIRCS[ grep('SLC45A4|ARID1A', CIRCS$circ_name) ]
rm(CIRCS.all, GENES, nb.genes, hb.genes, vt.genes, nb.circs, hb.circs, vt.circs)


#  compute circRNA (mRNA) CPM based on circular (external linear junction counts) and the total number of junction counts
#{{{

#  [neuroblastoma tumors] make sure to remove failed/non-tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-unlist(lin.cir)[bid %in% nb.meta$bid]
nb<-lin.cir[ circ_name %in% CIRCS$circ_name ]
nb<-nb[, c('c.cpm', 'l.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)]
nb.cir<-data.frame(dcast(nb, bid ~ circ_name, value.var='c.cpm', fun.aggregate=sum), check.names=F)
nb.lin<-data.frame(dcast(nb, bid ~ circ_name, value.var='l.cpm', fun.aggregate=sum), check.names=F)
rm(lin.cir, nb)


#  [human brain tissue]
load('/fast/projects/Schulte_NB/work/downloads/paired-end_totalRNA_human_brain/neuroblastoma_circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-unlist(lin.cir)
hb<-lin.cir[ circ_name %in% CIRCS$circ_name ]
hb<-hb[, c('c.cpm', 'l.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)]
hb.cir<-data.frame(dcast(hb, bid ~ circ_name, value.var='c.cpm', fun.aggregate=sum), check.names=F)
hb.lin<-data.frame(dcast(hb, bid ~ circ_name, value.var='l.cpm', fun.aggregate=sum), check.names=F)
rm(lin.cir, hb)


#  [various tumors] make sure non-tumors are not present (they should not be)
load('/fast/projects/Schulte_NB/work/downloads/paired-end_totalRNA_cancer/neuroblastoma_circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-unlist(lin.cir)[ bid %in% vt.meta[, unlist(srr)] ]
vt<-lin.cir[ circ_name %in% CIRCS$circ_name ]
vt<-vt[, c('c.cpm', 'l.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)]
vt.cir<-data.frame(dcast(vt, bid ~ circ_name, value.var='c.cpm', fun.aggregate=sum), check.names=F)
vt.lin<-data.frame(dcast(vt, bid ~ circ_name, value.var='l.cpm', fun.aggregate=sum), check.names=F)
rm(lin.cir, vt)

#}}}


#  do all plots for circRNAs and mRNAs
#  NA values result in NA p-values
res<-List()
for (circ_name in CIRCS$circ_name){
    cat('\n', circ_name, '...', sep='')
    
    #  isolate circRNA + mRNA CPMs
    nb.c<-setNames(nb.cir[, circ_name, drop=T], nb.cir$bid)
    nb.l<-setNames(nb.lin[, circ_name, drop=T], nb.lin$bid)
    nb.l<-nb.l[ names(nb.c) ]
    stopifnot(all.equal(names(nb.c), names(nb.l)))
    hb.c<-setNames(hb.cir[, circ_name, drop=T], hb.cir$bid)
    hb.l<-setNames(hb.lin[, circ_name, drop=T], hb.lin$bid)
    hb.l<-hb.l[ names(hb.c) ]
    stopifnot(all.equal(names(hb.c), names(hb.l)))
    vt.c<-setNames(vt.cir[, circ_name, drop=T], vt.cir$bid)
    vt.l<-setNames(vt.lin[, circ_name, drop=T], vt.lin$bid)
    vt.l<-vt.l[ names(vt.c) ]
    stopifnot(all.equal(names(vt.c), names(vt.l)))

    
    #  call the function that will do both circRNA and mRNA plots and Mann-Whitney tests
    #  it returns a list of two separate lists of CPMs + pvalues for the circRNA and the mRNA, respectively
    #  
    #  N.B. we ask for neuroblastoma-specific significant differential expression here
    res[[circ_name]]<-do_plots(circ_name=circ_name, d=list(circ=list(nb=nb.c, hb=hb.c, vt=vt.c, alt='two.sided'), lin=list(nb=nb.l, hb=hb.l, vt=vt.l, alt='two.sided')), do.plots=T)
    cat('\n')
}

#  no need to save anything, execution time is very fast

#}}}



#  [circARID1A, ARID1A] expression among risk groups 
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  keep only non-failed neuroblastoma tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
nb.meta<-meta.tum[ !(failed) ]
rm(meta.tum, meta.cel, meta.prefailed)


#  load the unified set of circRNAs and isolate ARID1A 
#  add as zeros samples where ARID1A is not expressed
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
int<-CIRCS[ grep('ARID1A', CIRCS$circ_name) ]
nb.circs<-nb.circs[ nb.circs$bid %in% nb.meta$bid & nb.circs$circ_name %in% int$circ_name ]
CIRCS<-int
rm(CIRCS.all, GENES, nb.genes, hb.genes, vt.genes, int, hb.circs, vt.circs, hb.meta, vt.meta)


#  compute mRNA CPM based on external linear junction counts and the total (linear+circular) counts
#  compute circRNA CPM based on the circular junction counts and the total (linear+circular) counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-unlist(lin.cir)[bid %in% nb.meta$bid]
nb<-lin.cir[ circ_name %in% CIRCS$circ_name ]
nb<-nb[, c('c.cpm', 'l.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)]
nb.cir<-dcast(nb, bid ~ circ_name, value.var='c.cpm', fun.aggregate=sum)
nb.lin<-dcast(nb, bid ~ circ_name, value.var='l.cpm', fun.aggregate=sum)
nb.cir<-setNames(unname(unlist(nb.cir[, 2])), unlist(nb.cir[, 1]))
nb.lin<-setNames(unname(unlist(nb.lin[, 2])), unlist(nb.lin[, 1]))
rm(lin.cir)


#  correlation of circular and linear counts
nb[, summary(lm(l.count.out.max ~ c.count))]
# 
# Residuals:
#    Min     1Q Median     3Q    Max 
# -625.2 -175.5  -57.6   87.4 1140.3 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  330.906     39.342    8.41  2.9e-13 ***
# c.count        4.719      0.379   12.45  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 295 on 100 degrees of freedom
#   (2 observations deleted due to missingness)
# Multiple R-squared:  0.608,	Adjusted R-squared:  0.604 
# F-statistic:  155 on 1 and 100 DF,  p-value: <2e-16


#  boxplot of the circRNA across risk groups
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
options(scipen=0)
par(mar=c(10.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-nb.cir
names(x)<-nb.meta[ match(names(x), bid), risk_group ]
B<-split(x, names(x))[c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA')]
B.cl<-setNames(nb.meta[, unique(col)], nb.meta[, unique(risk_group)])
B.cl<-B.cl[ names(B) ]
YTICK<-pretty(c(0.0, sapply(B, max)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('CPM', side=2, line=3, padj=-0.9, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_circRNA_CPMs_for_', CIRCS$circ_name, '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x,B,B.cl,YTICK,bp)
dev.off()


#  boxplot of the mRNA across risk groups
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(10.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-nb.lin
names(x)<-nb.meta[ match(names(x), bid), risk_group ]
B<-split(x, names(x))[c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA')]
B.cl<-setNames(nb.meta[, unique(col)], nb.meta[, unique(risk_group)])
B.cl<-B.cl[ names(B) ]
YTICK<-pretty(c(0.0, sapply(B, max)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('CPM', side=2, line=5, padj=+0.4, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_linear_junction_CPMs_for_', CIRCS$circ_name, '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x,B,B.cl,YTICK,bp)
dev.off()

#}}}



#  [circSLC45A4, across risk groups] expression comparison of circRNA of interest
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  keep only non-failed neuroblastoma tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
nb.meta<-meta.tum[ !(failed) ]
rm(meta.tum, meta.cel, meta.prefailed)


#  load the unified set of circRNAs and isolate circSLC45A4
#  add as zeros samples where the circ is not expressed
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
int<-CIRCS[ grep('SLC45A4', CIRCS$circ_name) ]
nb.circs<-nb.circs[ nb.circs$bid %in% nb.meta$bid & nb.circs$circ_name %in% int$circ_name ]
CIRCS<-int
rm(CIRCS.all, GENES, nb.genes, hb.genes, vt.genes, int, hb.circs, vt.circs, hb.meta, vt.meta)


#  compute mRNA CPM based on external linear junction counts and the total (linear+circular) counts
#  compute circRNA CPM based on the circular junction counts and the total (linear+circular) counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-unlist(lin.cir)[bid %in% nb.meta$bid]
nb<-lin.cir[ circ_name %in% CIRCS$circ_name ]
nb<-nb[, c('c.cpm', 'l.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)]
nb.cir<-dcast(nb, bid ~ circ_name, value.var='c.cpm', fun.aggregate=sum)
nb.lin<-dcast(nb, bid ~ circ_name, value.var='l.cpm', fun.aggregate=sum)
nb.cir<-setNames(unname(unlist(nb.cir[, 2])), unlist(nb.cir[, 1]))
nb.lin<-setNames(unname(unlist(nb.lin[, 2])), unlist(nb.lin[, 1]))
rm(lin.cir, nb)


#  boxplot of the circRNA across risk groups
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
options(scipen=0)
par(mar=c(10.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-nb.cir
names(x)<-nb.meta[ match(names(x), bid), risk_group ]
B<-split(x, names(x))[c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA')]
B.cl<-setNames(nb.meta[, unique(col)], nb.meta[, unique(risk_group)])
B.cl<-B.cl[ names(B) ]
YTICK<-pretty(c(0.0, sapply(B, max)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('CPM', side=2, line=3, padj=-0.9, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_circRNA_CPMs_for_', CIRCS$circ_name, '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x,B,B.cl,YTICK,bp)
dev.off()


#  boxplot of the mRNA across risk groups
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(10.0, 6.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
x<-nb.lin
names(x)<-nb.meta[ match(names(x), bid), risk_group ]
B<-split(x, names(x))[c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA')]
B.cl<-setNames(nb.meta[, unique(col)], nb.meta[, unique(risk_group)])
B.cl<-B.cl[ names(B) ]
YTICK<-pretty(c(0.0, sapply(B, max)), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('CPM', side=2, line=5, padj=+0.4, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_risk_groups_linear_junction_CPMs_for_', CIRCS$circ_name, '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(x,B,B.cl,YTICK,bp)
dev.off()

#}}}



#  [HR_nMNA-specific circRNAs] mRNA expression estimated also from the linear junction counts normalized to CPMs
#                              significance is estimated by two-sided Mann-Whitney tests
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  functions
#{{{

source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/add_sig_bars.R')

do_plots<-function(circ_name='', alternative='greater', do.plots=T){

    #  isolate the linear junction CPMs of the chosen circRNA
    nb.<-nb.lin[, c('bid', circ_name), with=F]
    colnames(nb.)<-c('bid', 'circ')
    hb.<-hb.lin[, c('bid', circ_name), with=F]
    colnames(hb.)<-c('bid', 'circ')
    vt.<-vt.lin[, c('bid', circ_name), with=F]
    colnames(vt.)<-c('bid', 'circ')
    nb.<-setNames(nb.$circ, nb.$bid)
    hb.<-setNames(hb.$circ, hb.$bid)
    vt.<-setNames(vt.$circ, vt.$bid)


    #  replace bid by tissue throughout
    names(nb.)<-rep('neuroblastoma', length(nb.))
    names(hb.)<-rep('brain tissue', length(hb.))
    names(vt.)<-rep('various tumors', length(vt.))


    #  do Mann-Whitney tests (does also human brain vs various tumors but this is ignored in the analysis)
    B<-list(nb=nb., hb=hb., vt=vt.)
    if(!all(is.na(nb.))){
        PV<-do.call(rbind, apply(combn(names(B), 2), 2, function(s){ data.frame(pv=wilcox.test(x=B[[s[1]]], y=B[[s[2]]], alternative=alternative)$p.value, row.names=paste(s[1], s[2], collapse=' ')) }))
    } else {
        PV<-data.frame(pv=c(NA, NA, NA), row.names=apply(combn(names(B), 2), 2, paste, collapse=' '))
    }
    PV$symbol<-'ns'
    PV$symbol[ PV$pv<1e-10 ]<-'***'
    PV$symbol[ PV$pv>=1e-10 & PV$pv<0.001 ]<-'**'
    PV$symbol[ PV$pv>=0.001 & PV$pv<0.05 ]<-'*'
    rm(B)

    
    if(do.plots){
        if(!all(is.na(nb.))){
            #  barplot of ordered CPMs 
            x11(width=40, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
            B<-c( sort(nb., decreasing=T, na.last=T), 
                  sort(hb., decreasing=T, na.last=T), 
                  sort(vt., decreasing=T, na.last=T) )
            B.cl<-c( setNames(rep('seagreen4', length(nb.)), names(nb.)), 
                     setNames(rep('cornflowerblue', length(hb.)), names(hb.)), 
                     setNames(rep('coral4', length(vt.)), names(vt.)) )
            par(mar=c(1.0, 5.5, 1.0, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
            YTICK<-pretty(c(0, max(B, na.rm=T)), 5)
            bp<-barplot(B, beside=F, plot=F)
            plot(0:1, 0:1, type='n', ylim=c(0, tail(YTICK, 1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
            bp<-barplot(B, border='white', col=B.cl, axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
            axis(2, at=YTICK, line=0, cex.axis=2.4)
            legend('topright', legend=c('neuroblastoma', 'human brain tissue', 'various tumors'), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.70, x.intersp=0.2, seg.len=0.5)
            dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_across_tissues_linear_junction_CPMs_for_', circ_name, '.svg'), width=40, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
            dev.off()


            #  boxplot of CPMs
            options(scipen=+20)
            x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
            par(mar=c(5.0, 8.5, 2.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
            B<-list(nb=nb., hb=hb., vt=vt.)
            B.cl<-c('seagreen4', 'cornflowerblue', 'coral4')
            YTICK<-pretty(c(0, max(sapply(B, max, na.rm=T))), 4)
            plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
            bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
            YTICK<-pretty(c(0, 1.2*max(bp$stats[5, ])), 4)  #  recompute y-max to leave room for significance bars
            plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
            bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
            mtext('Counts per million', side=2, line=6, padj=-0.1, las=0, cex=2.4)
            mtext(text=c('neuroblastoma', 'brain\ntissue', 'various\ntumors'), side=1, line=-1, at=seq_len(length(B)), las=1, padj=1, cex=2.4, col=B.cl)
            add_sig_bars(bp, PV)
            dev.print(device=svg, file=paste0('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_across_tissues_linear_junction_CPMs_for_', circ_name, '.svg'), width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
            options(scipen=0)
            dev.off()
        } else {
            cat('\n\n**** No mRNA expression for:', circ_name, '****\n\n')
        }
    }


    return(list(nb.=nb., hb.=hb., vt.=vt., pv=PV))
}

#}}}


#  load HR_nMNA-specific circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_HR_nMNA-specific.RData')


#  [HR_nMNA tumors] compute mRNA CPMs based on external linear junction counts and total counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-unlist(lin.cir)[bid %in% colnames(nb.circs)]
lin.cir<-lin.cir[ circ_name %in% rownames(nb.circs) ]
lin.cir<-lin.cir[, cpm:=l.count.out/total*1e6]
nb.lin<-dcast(lin.cir, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)
nb.lin<-nb.lin[ match(colnames(nb.circs), bid) ]
rm(lin.cir,nb.meta)


#  [human brain tissue] compute mRNA CPMs based on external linear junction counts and total counts
load('/fast/projects/Schulte_NB/work/downloads/paired-end_totalRNA_human_brain/neuroblastoma_circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-unlist(lin.cir)[ bid %in% colnames(hb.circs) ]  #  unnecessary since we keep all samples but good for consistency
lin.cir<-lin.cir[ circ_name %in% rownames(nb.circs) ]
lin.cir<-lin.cir[, cpm:=l.count.out/total*1e6]
hb.lin<-dcast(lin.cir, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)
hb.lin<-hb.lin[ match(colnames(hb.circs), bid) ]
rm(lin.cir,hb.meta)


#  [various tumors] compute mRNA CPMs based on external linear junction counts and total counts
#                   remove the two libraries with 52nt reads that are not included in the linear vs circular junction analysis
load('/fast/projects/Schulte_NB/work/downloads/paired-end_totalRNA_cancer/neuroblastoma_circRNAs_linear_vs_circular_collected_results.RData')
vt.circs<-vt.circs[, grep('SRR1617643|SRR1617644', colnames(vt.circs), invert=T)]
lin.cir<-unlist(lin.cir)[ bid %in% vt.meta[, unlist(srr)] ]
lin.cir<-lin.cir[ circ_name %in% rownames(nb.circs) ]
lin.cir<-lin.cir[, cpm:=l.count.out/total*1e6]
vt.lin<-dcast(lin.cir, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)
vt.lin<-vt.lin[ match(colnames(vt.circs), bid) ]
rm(lin.cir,vt.meta)


#  do all plots
#  NA values result in NA p-values
res<-List()
for (circ_name in colnames(nb.lin[, -1])){
    cat('\n', circ_name, ':\n', sep='')
    res[[circ_name]]<-do_plots(circ_name=circ_name, alternative='two.sided', do.plots=F)
    cat('\n')
}


#  save
save(res, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/linear_junction_CPMs_of_circRNAs_across_tissues_HR_nMNA-specific.RData')

#}}}



#  [HR_nMNA-specific circRNAs] plots including annotations for the mRNA expression computed above
#{{{
rm(list=ls())
library(GenomicAlignments)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(grid)
library(viridis)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/add_sig_bars.R')
source('~/bio/lib/draw_highlights.R')


#  functions
#{{{
boxplot_per_risk_group<-function(X, META){
    colnames(X)<-META[ match(colnames(X), bid), risk_group ]
    B<-split(log10(1+X), colnames(X))[META[, unique(risk_group)]]
    B.cl<-META[ risk_group %in% names(B), unique(col) ] 
    YTICK<-pretty(c(0.0, sapply(B, max)), 5)
    plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
    bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
    mtext(expression(log[10](1+'CPM')), side=2, line=5, padj=+0.1, las=0, cex=2.4)
    mtext(text=paste0(names(B), ' (', lengths(B)/nrow(X), ')'), side=1, line=0, at=seq(1, length(B), 1), las=2, adj=1, cex=2.4, col=B.cl)
}
#}}}


# load the HR_nMNA specific circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_HR_nMNA-specific.RData')
rm(list=setdiff(ls(), ls(pattern='circs|CIRCS|_|nb.meta')))  #  remove clutter


#  load the mRNA expression estimated from linear junction counts external to the corresponding circRNA and normalized to CPMs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/linear_junction_CPMs_of_circRNAs_across_tissues_HR_nMNA-specific.RData')


#  compare mean expressions to identify the consistently significantly up- or down- regulated mRNAs across both pairwise comparisons
up<-sapply(res, function(r){ 
    if(any(is.na(r$pv$pv))){ return(NA) }
    if(all(r$pv[c('nb hb', 'nb vt'), 'pv']<0.05)){
        if( mean(r$nb., na.rm=T)>mean(r$hb., na.rm=T) & mean(r$nb., na.rm=T)>mean(r$vt., na.rm=T) ){
            return(T)
        } else {
            return(F)
        }
    } else {
        return(F)
    }
})
down<-sapply(res, function(r){ 
    if(any(is.na(r$pv$pv))){ return(NA) }
    if(all(r$pv[c('nb hb', 'nb vt'), 'pv']<0.05)){
        if( mean(r$nb., na.rm=T)<mean(r$hb., na.rm=T) & mean(r$nb., na.rm=T)<mean(r$vt., na.rm=T) ){
            return(T)
        } else {
            return(F)
        }
    } else {
        return(F)
    }
})
down<-down[ rownames(nb.circs) ]
up_down<-up[ rownames(nb.circs) ]
up_down[down]<- -1
up_down[is.na(up_down)]<-0
rm(up,down)


#  recycle
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  heatmap on the log10-transformed CPMs (best contrast, compared to z-scores of raw CPMs, or z-scores of log10-transformed CPMs)
X<-cbind(log10(1+nb.circs), log10(1+hb.circs), log10(1+vt.circs))
col.an<-data.frame(tissue=factor(rep(c('neuroblastoma', 'brain tissue', 'various tumors'), c(ncol(nb.circs), ncol(hb.circs), ncol(vt.circs))), exclude=F), row.names=colnames(X))
row.an<-data.frame(mRNA=factor(up_down, levels=c(-1, 1, 0), labels=c('down', 'up', 'neither')), row.names=rownames(X))
cl<-setNames(list(setNames( c('seagreen4', 'cornflowerblue', 'coral4'), c('neuroblastoma', 'brain tissue', 'various tumors') ),
                  setNames( c('lightsteelblue3', 'grey20', 'grey49'), c('down', 'up', 'neither') )), 
             c(colnames(col.an), colnames(row.an)))
ph<-pheatmap(X, color=viridis(10), border_color=NA, scale='none', 
        cluster_rows=F,
        cluster_cols=F,
        annotation_legend=T, annotation_names_row=T, annotation_names_col=T, 
        annotation_col=col.an, annotation_row=row.an, annotation_colors=cl,
        drop_levels=F, show_rownames=T, show_colnames=F, 
        fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/heatmap_circRNAs_across_tissues_HR_nMNA-specific.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(X,col.an,row.an,cl,ph)


#  boxplots of log10-transformed CPMs
#  use circBase circ_id but replace hsa_circ with gene_name
options(scipen=0)
par(mar=c(15.0, 8.0, 0.5, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
stopifnot( all.equal(rownames(nb.circs), CIRCS$circ_name) )
NAMES<-sapply(apply(data.frame(mcols(CIRCS)[, c('gene_name', 'circ_id')]), 1, function(n){ sub('hsa_circ', paste0('circ', n[[1]]), n[[2]]) }) , paste, collapse=',\n')
B<-c(split(log10(1+nb.circs), row(nb.circs)), split(log10(1+hb.circs), row(hb.circs)), split(log10(1+vt.circs), row(vt.circs)))
n<-c(rbind(rbind(1:nrow(nb.circs), nrow(nb.circs)+(1:nrow(nb.circs)), 2*nrow(nb.circs)+(1:nrow(nb.circs)))))  #  order them in triplets
B<-B[n]
B.cl<-rep(c('seagreen4', 'cornflowerblue', 'coral4'), nrow(nb.circs))
m.cl<-up_down
m.cl[ m.cl==-1 ]<-'lightsteelblue3'
m.cl[ m.cl==1 ]<-'grey20'
m.cl[ m.cl==0 ]<-'grey49'
YTICK<-pretty(c(0, max(floor(sapply(B, max, na.rm=T)))), 4)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'CPM')), side=2, line=4, padj=-0.2, las=0, cex=2.4)
#mtext(text=NAMES, side=1, line=-1, at= seq(2, length(B), 3), las=2, adj=1.09, cex=2.4, col=m.cl)
mtext(text=sub('_.*$', '', NAMES), side=1, line=0, at= seq(2, length(B), 3), las=2, adj=1, cex=2.4, col=m.cl)  #  remove circBase id
draw_highlights(L=length(B), STEP=3, YMIN=0.0, YMAX=max(YTICK))
legend('topleft', legend=c('neuroblastoma', 'brain tissue', 'various tumors'), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.40, x.intersp=0.3, seg.len=0.3)
legend('topright', legend=c('down', 'up', 'neither'), col=c('lightsteelblue3', 'grey20', 'grey49'), title='mRNA', bty='n', lty=1, lwd=15, pch=NA, cex=1.5, xpd=T, y.intersp=0.50, x.intersp=0.3, seg.len=0.3)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_across_tissues_HR_nMNA-specific.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(YTICK,NAMES,B,n,B.cl,bp)
dev.off()

#}}}





