###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




################################################################################
#
#
#  define the neuroblastoma specific circRNA cohort and everything related to it
#
#
################################################################################




#  [run once] define the circRNA cohort that will be used throughout the analysis:
#
#      remove failed samples and keep only neuroblastoma tumors
#      keep only exon-exon circRNAs (trans-gene + chrM + chrY circRNAs have been removed already)
#      allow +-1nt corresponding exon boundary discrepancies for the exons with the donor and acceptor sites provided the gene_id matches
#      keep only circRNAs supported by at least 25% of samples OR have backspliced junction coverage > 20 in at least three samples 
#
#  we query the circRNA and corresponding gene expression across neuroblastoma tumors, human brain tissue and various tumors
#
#  N.B. circRNA CPMs are based on the total counts on gene features 
#
#  For the different cell models we identify the same unified cohort of circRNAs 
#
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()
source('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/load_circ_gene_expression.R')


#  load circRNAs 
#  remove failed samples and keep only tumors (Pilot samples are removed as well)
#  keep only the exon-exon circRNAs
#  add circ_name
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-meta.tum[ !(failed) ]
CIRCS<-circ[ circ$bid %in% meta$bid & circ$region %in% 'exon' ]
META<-meta[ bid %in% CIRCS$bid ]
CIRCS$circ_name<-paste0(CIRCS$gene_id, '|', CIRCS$gene_name,'_', as.character(seqnames(CIRCS)),as.character(strand(CIRCS)), start(CIRCS), '-', end(CIRCS))
CIRCS.all<-CIRCS
CIRCS<-unique(CIRCS)
mcols(CIRCS)<-mcols(CIRCS)[, c('gene_name', 'gene_id', 'circ_name')]
rm(circ, meta.tum, meta.cel, meta.prefailed, trans)


#  load the annotated exons 
#  enforce the seqlevels of circRNAs (remove chrM)
txdb<-loadDb('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/TxDb_GRCh38.gencode.v30.RData')
EXONS<-exons(txdb, columns=c('EXONNAME', 'GENEID'))
colnames(mcols(EXONS))<-c('exon_name', 'gene_id')
seqlevels(EXONS, pruning.mode='coarse')<-seqlevels(CIRCS)
stopifnot( all(lengths(EXONS$gene_id)==1) )
EXONS$gene_id<-unlist(EXONS$gene_id)
rm(txdb)


#  allow for +-1nt discrepancies with the exon boundaries provided circRNA and donor/acceptor exon gene_ids match!
#{{{

#  allow +-1nt discrepancy at the donor site
ov<-findOverlaps(EXONS, resize(CIRCS, 2, fix='start'), type='any')
ex<-EXONS[ queryHits(ov) ]
ci<-CIRCS[ subjectHits(ov) ]
ov<-ov[ ex$gene_id==ci$gene_id & ( ( as.logical( strand(ci) %in% '+') & abs(start(ex)-start(ci))<=1 ) | 
                                   ( as.logical( strand(ci) %in% '-') & abs(end(ex)-end(ci))<=1 ) ) ]
f<-CIRCS[ setdiff(seq_along(CIRCS), subjectHits(ov)) ]  #  inspect those dropping out
CIRCS<-CIRCS[ unique(subjectHits(ov)) ]
rm(f,ov,ex,ci)


#  allow +-1nt discrepancy at the acceptor site
ov<-findOverlaps(EXONS, resize(CIRCS, 2, fix='end'), type='any')
ex<-EXONS[ queryHits(ov) ]
ci<-CIRCS[ subjectHits(ov) ]
ov<-ov[ ex$gene_id==ci$gene_id & ( ( as.logical( strand(ci) %in% '+') & abs(end(ex)-end(ci))<=1 ) | 
                                   ( as.logical( strand(ci) %in% '-') & abs(start(ex)-start(ci))<=1 ) ) ]
f<-CIRCS[ setdiff(seq_along(CIRCS), subjectHits(ov)) ]  #  inspect those dropping out
CIRCS<-CIRCS[ unique(subjectHits(ov)) ]
rm(f,ov,ex,ci,EXONS)

#}}}


#  keep only the circRNAs supported by at least 25% of samples OR have counts > 20 in at least three of the samples
x<-data.table(circ_name=CIRCS.all$circ_name, bid=CIRCS.all$bid, jc_count=CIRCS.all$jc_count)[ circ_name %in% CIRCS$circ_name ]
x<-dcast(x, bid ~ circ_name, value.var='jc_count')[, bid:=NULL]
keep<-unlist(c(as.data.frame(x[, lapply(.SD, function(x){ sum(is.na(x))/.N<0.75 | sum(x>20, na.rm=T)>=3 })])))
nf<-CIRCS.all[ CIRCS.all$circ_name %in% names(keep[ keep ]) ]
nf<-nf[ order(nf$jc_count, decreasing=T) ]
nf[ nf$gene_name %in% 'ARID1A' ]     #  chr1:26729651-26732792  
nf[ nf$gene_name %in% 'SETD3' ]      #  chr14:99458279-99465813 
nf[ nf$gene_name %in% 'EYA1' ]       #  chr8:71299047-71334174
nf[ nf$gene_name %in% 'TET2' ]       #  chr4:105233897-105237351
nf[ nf$gene_name %in% 'SMARCA5' ]    #  chr4:143543509-143543972
nf[ nf$gene_name %in% 'LINC00632' ]  #  CDR1as == LINC00632 chrX:140783175-140784659
nf[ nf$gene_name %in% 'HIPK3' ]      #  chr11:33286413-33287511
nf[ nf$gene_name %in% 'ATRX' ]       #  chrX:77652114-77656653 
nf[ nf$gene_name %in% 'HUWE1' ]      #  chrX:53614534-53615835, chrX:53645311-53654131
nf[ nf$gene_name %in% 'EZH2' ]       #  chr7:148846470-148847305
nf[ nf$gene_name %in% 'VRK1' ]       #  chr14:96833467-96860735
nf[ nf$gene_name %in% 'CHD7' ]       #  chr8:60794986-60801593, chr8:60741259-60743097
nf[ nf$gene_name %in% 'KDM1A' ]      #  chr1:23030469-23059167, chr1:23030469-23050520 
CIRCS.all<-nf
rm(nf,keep,x)


#  unique calls (removing bid, count information since it is irrelevant)
CIRCS<-unique(CIRCS.all)
mcols(CIRCS)<-mcols(CIRCS)[, c('gene_id', 'gene_name', 'circ_name')]


#  keep the unique genes also making sure to include Clara's housekeeping gene as well
GENES<-unique(mcols(CIRCS)[, c('gene_id', 'gene_name')])
if(nrow( GENES[ GENES$gene_id %in% 'ENSG00000073578.17' ,] )==0){
    GENES<-rbind(GENES, DataFrame(gene_id=c('ENSG00000073578.17'), gene_name=c('SDHA')))
}


#  import circRNA and gene expression in CPMs
l<-load_circ_gene_expression(CIRCS, GENES, META, nb.tumors.only=T, vt.tumors.only=T, gene.tpm=F)


#  inflate
nb.circs<-l$nb$circs
nb.genes<-l$nb$genes
nb.meta<-l$nb$meta
hb.circs<-l$hb$circs
hb.genes<-l$hb$genes
hb.meta<-l$hb$meta
vt.circs<-l$vt$circs
vt.genes<-l$vt$genes
vt.meta<-l$vt$meta
rm(l)


#  save 
save(CIRCS.all, CIRCS, GENES, hb.circs, hb.genes, hb.meta, nb.circs, nb.genes, nb.meta, vt.circs, vt.genes, vt.meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')


#  summarize circRNAs into an HTML table and Excel sheet
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


#  load circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')


#  summarize duplicate calls across samples
l<-CIRCS.all
l<-data.table(as.data.frame(l))[, .(jc_count=list(jc_count), non_jc_count=list(non_jc_count), ratio=list(ratio), bid=list(bid)), by=.(seqnames, start, end, strand, gene_name, gene_id, region, circ_name)]
l<-GRanges(as.data.frame(l))
l$jc_count<-as(l$jc_count, 'CompressedList')
l$non_jc_count<-as(l$non_jc_count, 'CompressedList')
l$ratio<-as(l$ratio, 'CompressedList')
l$bid<-as(l$bid, 'CompressedList')


#  create HTML table with proper links
l.table<-granges2kable(l)


#  save as R object
save(l, l.table, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_summarized.RData')


#  save as Excel sheet after converting list columns to character vectors
x<-data.table(cbind(data.frame(position=paste0(as.character(seqnames(l)), '(', as.character(strand(l)),'):', start(l), '-', end(l))), as.data.frame(l)[, -c(1:3,5)]))
x[, ':='(jc_count=sapply(jc_count, paste, sep='', collapse=','), 
         non_jc_count=sapply(non_jc_count, paste, sep='', collapse=','), 
         ratio=sapply(ratio, paste, sep='', collapse=','), 
         bid=sapply(bid, paste, sep='', collapse=','))]
write.xlsx(as.data.frame(x), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_summarized.xlsx', col.names=T, row.names=F, sheetName='summarized circRNAs', append=F)

#}}}


#  identify the same circRNA cohort in cell models
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  identify cell models
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-meta.cel[ ! (failed) & !is.na(cell_model) ]
rm(meta.tum, meta.cel, meta.prefailed)


#  load unified circRNAs based on the neuroblastoma tumors and keep their names
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
circ_name<-CIRCS$circ_name
rm(list=ls(pattern='\\.meta|\\.circs|\\.genes|GENES|CIRCS'))


#  load all CIRI2 circRNAs 
#  keep those defined by the tumor cohort found in the cell models
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/circRNAs_CIRI2.RData')
CIRCS.all<-circ[ circ$bid %in% meta$bid & circ$region %in% 'exon' ]
CIRCS.all$circ_name<-paste0(CIRCS.all$gene_id, '|', CIRCS.all$gene_name,'_', as.character(seqnames(CIRCS.all)),as.character(strand(CIRCS.all)), start(CIRCS.all), '-', end(CIRCS.all))
CIRCS.all<-CIRCS.all[ CIRCS.all$circ_name %in% circ_name ]
rm(circ,circ_name)


#  save
save(CIRCS.all, meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_cell_models.RData')

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_summarized.{RData,xlsx}
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_cell_models.RData



#  [run once] construct circRNA sequences and controls:
#             
#                 circRNA sequences:
#             
#                     we identify all exons within a given circRNA of identical gene_id as the circRNA and among overlapping exons we take the longest 
#                     we set the first/last exon boundaries to be identical to the circRNA boundaries, irrespective of possible +-1nt discrepancies
#             
#                 controls:
#                     
#                     identify all exons downstream of the last circRNA exon and upstream of the first circRNA exon 
#                     take the longest among overlapping ones in both lists 
#                     take the upstream or downstream list that produces the longest spliced sequence
#                     when possible trim the spliced sequence to identical length as the circRNA sequence by fixing the center (symmetrical trimming)
#                     discard empty lists when there is no upstream or downstream exon
#
#{{{
rm(list=ls())
library(data.table)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)


#  load circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
rm(list=c(ls(pattern='[nhv][tb]\\.*'), 'GENES', 'CIRCS.all'))


#  extract all annotated exons
#  take unique exon entries per gene
txdb<-loadDb('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/TxDb_GRCh38.gencode.v30.RData')
EXONS<-unlist(exonsBy(txdb, 'gene'))
EXONS$gene_id<-names(EXONS)
names(EXONS)<-NULL
x<-data.table( seqnames=as.character(seqnames(EXONS)), strand=as.character(strand(EXONS)), start=start(EXONS), end=end(EXONS), exon_name=EXONS$exon_name, gene_id=EXONS$gene_id)[, ex:=paste(seqnames, strand, start, end, sep='_')][, .(exon_name=exon_name[ !duplicated(ex) ], ex=ex[ !duplicated(ex) ]), by=.(gene_id)]
EXONS<-EXONS[ EXONS$exon_name %in% x$exon_name ]
mcols(EXONS)<-mcols(EXONS)[, 'gene_id', drop=F]
rm(x)


#  [~10min] circRNA sequences
#{{{

#  At this point the group of circRNAs has +-1nt agreement with the donor/acceptor exon boundaries,
#  so we identify now all exons in between making sure again that the gene_ids match
ov<-findOverlaps(EXONS, resize(CIRCS, width(CIRCS)+2, fix='center'), select='all', type='within')
ov<-ov[ EXONS$gene_id[ queryHits(ov) ]==CIRCS$gene_id[ subjectHits(ov) ] ]
print(CIRCS[ setdiff(seq_along(CIRCS), subjectHits(ov)) ])  #  there should not be any circRNAs rejected at this point


#  [<7min] group exons by circRNA
#          reverse the ordering for minus strand exons in order for extractTranscriptSeqs() to work properly
#          keep all non-overlapping and the longest among overlapping
EX<-split(EXONS[queryHits(ov)], subjectHits(ov))
CIRCS.exons<-GRangesList()
for(n in seq_along(EX)){
    X<-EX[[n]]
    R<-CIRCS[ as.integer(names(EX)[n]) ]
    M<-as.logical(strand(X)[1]=='-')  #  minus strand?

    #  reverse the order if on minus strand
    if(M){
        X<-X[ order(start(X), decreasing=T) ]
    }


    #  find out if there are overlapping exons and take the longest
    o<-findOverlaps(X, X, type='any', select='all')
    if(!all(queryHits(o)==subjectHits(o))){ 
        o<-data.frame(o[ queryHits(o)<subjectHits(o) ])
        w<-width(X)

        #  pairwise comparison of indices, if one exon overlaps many then sequentially the shortest is always discarded, so the same
        #  index might appear multiple times in the list to discard but this is fine
        X<-X[-apply(o, 1, function(x){ x[ which.min(w[x]) ] })]
    }


    #  if on plus (minus) strand make sure the start (end) of the first (last) exon and the end (start) of the last (first) exon are identical
    #  to those of the corresponding circRNA
    if(!M){
        start(X[1])<-start(R)
        end(X[length(X)])<-end(R)
    } else {
        end(X[1])<-end(R)
        start(X[length(X)])<-start(R)
    }

    CIRCS.exons[[n]]<-X
}
names(CIRCS.exons)<-names(EX)
stopifnot( length(CIRCS.exons)==length(CIRCS) )  #  split() has re-ordered them by ascending order of queryHits(), i.e. CIRCS order
rm(ov,n,EX)


#  check number of gene_ids per circRNA group of exons
g<-lapply(CIRCS.exons, function(x){ unique(x$gene_id) })
table(lengths(g))
stopifnot( all(lengths(g)==1) )
stopifnot( all.equal( CIRCS$gene_id, unname(unlist(g)) ) )
rm(g)


#  get the circRNA sequences
#
#  ----->!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!<-----
#  ----->!!!!!                                                                                     !!!!!<-----
#  ----->!!!!!  extractTranscriptSeqs() pastes together exons sequences IN THE ORDER THEY APPEAR   !!!!!<-----
#  ----->!!!!!  exonsBy(TxDb, 'tx') orders exons by rank which for minus-strand transcripts is     !!!!!<-----
#  ----->!!!!!  from right to left so the pasting would result in the correct transcript sequence  !!!!!<-----
#  ----->!!!!!                                                                                     !!!!!<-----
#  ----->!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!<-----
#
names(CIRCS.exons)<-CIRCS$circ_name
CIRCS.exons.seqs<-extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, CIRCS.exons)

#}}}


#  [~10min] control sequences
#{{{

CIRCS.controls<-endoapply(CIRCS.exons, function(X){
    ex<-EXONS[ EXONS$gene_id %in% X$gene_id[1] ]; 
    L<-length(ex);
    if (as.logical(strand(ex)[1]=='-')){
        ex<-ex[ order(start(ex), decreasing=T) ]
    };
    
    #  first and last circRNA exons are enough to identify the boundaries
    #  this will work for one-exon circRNAs as well, since first and last will be identical in this case
    o<-findOverlaps(X[c(1, length(X))], ex, select='all', type='any' )  
    up<-subjectHits(o)[1]-1
    down<-subjectHits(o)[2]+1  
    ex.up<-ex.down<-GRanges()
    if(up>0){
        ex.up<-ex[1:up]
        w<-width(ex.up);

        #  keep longest exons among overlapping exons
        o<-findOverlaps(ex.up, ex.up, type='any', select='all');
        if(!all(queryHits(o)==subjectHits(o))){
            o<-data.frame(o[ queryHits(o)<subjectHits(o) ]);
            ex.up<-ex.up[-apply(o, 1, function(x){ x[ which.min(w[x]) ] })]
        }
    };
    if(down<=L){
        ex.down<-ex[down:L]
        w<-width(ex.down);

        #  keep longest exons among overlapping exons
        o<-findOverlaps(ex.down, ex.down, type='any', select='all');
        if(!all(queryHits(o)==subjectHits(o))){
            o<-data.frame(o[ queryHits(o)<subjectHits(o) ]);
            ex.down<-ex.down[-apply(o, 1, function(x){ x[ which.min(w[x]) ] })]
        }
    };
    if(sum(width(ex.up))>=sum(width(ex.down))){
        return(ex.up)
    } else {
        return(ex.down)
    };
})


#  only ENSG00000281508.1|CDR1_chrX+140783176-140784660 drops out as expected there are no upstream/downsteam exons available
names(CIRCS.controls[ lengths(CIRCS.controls)==0 ])
CIRCS.controls<-CIRCS.controls[ lengths(CIRCS.controls)>0 ]


#  get the control sequences
#
#  ----->!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!<-----
#  ----->!!!!!                                                                                     !!!!!<-----
#  ----->!!!!!  extractTranscriptSeqs() pastes together exons sequences IN THE ORDER THEY APPEAR   !!!!!<-----
#  ----->!!!!!  exonsBy(TxDb, 'tx') orders exons by rank which for minus-strand transcripts is     !!!!!<-----
#  ----->!!!!!  from right to left so the pasting would result in the correct transcript sequence  !!!!!<-----
#  ----->!!!!!                                                                                     !!!!!<-----
#  ----->!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!<-----
#
CIRCS.controls.seqs<-extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, CIRCS.controls)

#}}}


#  statistics of circRNA sequence lengths vs control sequence lengths
x<-CIRCS.exons.seqs[ names(CIRCS.controls.seqs) ]
y<-CIRCS.controls.seqs
d<-width(x) - width(y)
summary(d)
#
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -75782.000  -7728.500  -4941.000  -5816.729  -2927.000   8654.000
hist( d, breaks=100, xlim=c(-4e4, 2e4))  #  most of control sequences are much larger than the circRNA sequences


#  resize control sequences to circRNA length if possible
i<-IRanges(start=1, end=width(CIRCS.controls.seqs))
ir<-resize(i, fix='center', width=width(CIRCS.exons.seqs[ names(CIRCS.controls.seqs) ]))
wrong<-start(ir)<=0 | end(ir)>width(i)          #  out of bound resizing either at the left or at the right
ir[ wrong ]<-i[ wrong ]                         #  should be corrected by replacing it with the original ranges
CIRCS.controls.seqs.resized<-unlist(extractAt(CIRCS.controls.seqs, as(ir, 'IRangesList')))
x<-CIRCS.exons.seqs[ names(CIRCS.controls.seqs.resized) ]
y<-CIRCS.controls.seqs.resized
d<-width(x) - width(y)
cat('Control sequences with identical length as circRNA sequences = ', sum(d==0), ' (', round(sum(d==0)/length(CIRCS.exons.seqs)*100, 1), '%)\n', sep='')
#
# Control sequences with identical length as circRNA sequences = 5104 (98.1%)


#  save everything and export circRNA sequences and resized control sequences
save(CIRCS, CIRCS.exons, CIRCS.exons.seqs, EXONS, CIRCS.controls, CIRCS.controls.seqs, CIRCS.controls.seqs.resized, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')
writeXStringSet(CIRCS.exons.seqs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.fa', format='fasta')
writeXStringSet(CIRCS.controls.seqs.resized, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_controls_sequences.fa', format='fasta')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.fa
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_controls_sequences.fa



#  [run once] construct upstream and downstream intronic sequences of circRNAs and controls
#
#      circRNA introns:
#
#          we pick the longest upstream and downstream introns (>=15nts) in each case when available
#
#      controls:
#
#          we identify all introns (>=15nts) not overlapping with the circRNA and its upstream and downstream introns
#          among overlapping introns in that list we take the longest
#          we keep as many introns as possible until their cumulative length is as long as the cumulative length of the circRNA introns
#
#{{{
rm(list=ls())
library(data.table)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)


#  load circRNA sequences
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')


#  load reference
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'transcript' ]
mcols(hsa)<-mcols(hsa)[, c('type', 'gene_id', 'gene_name', 'gene_type', 'transcript_id', 'transcript_name', 'transcript_type')]


#  extract all introns of all transcripts that have at least one
#  keep only those that are at least 15nts long
txdb<-loadDb('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/TxDb_GRCh38.gencode.v30.RData')
INTRONS<-intronsByTranscript(txdb, use.names=T)
seqlevels(INTRONS, pruning.mode='coarse')<-seqlevels(CIRCS)  #  if one hits the seqlevel, all will hit it
INTRONS<-unlist(INTRONS[lengths(INTRONS)>0])
INTRONS$transcript_id<-names(INTRONS)
names(INTRONS)<-NULL
INTRONS$gene_id<-hsa$gene_id[ match(INTRONS$transcript_id, hsa$transcript_id) ]
INTRONS<-INTRONS[ width(INTRONS)>=15 ]


#  keep only introns of genes that produce our circRNAs
INTRONS<-INTRONS[ INTRONS$gene_id %in% CIRCS$gene_id ]


#  identify the list of unique introns per gene
x<-data.table(data.frame(INTRONS))[, c('width', 'transcript_id', 'pos'):=list(NULL, NULL, paste(start, end, sep='-'))][, .(pos=unique(pos)), by=.(gene_id, seqnames,strand)][, c('start', 'end'):=tstrsplit(pos, '-', fixed=T)][, pos:=NULL]
INTRONS<-GRanges(seqnames=x$seqnames, strand=x$strand, ranges=IRanges(start=as.integer(x$start), end=as.integer(x$end)), gene_id=x$gene_id)
rm(x)


#  CDR1as is in GENCODE v30 LINC00632 which does have introns defined, so nothing drops out
CIRCS[ CIRCS$gene_id %in% setdiff(CIRCS$gene_id, INTRONS$gene_id) ]  #  empty


#  [~2min] identify upstream and downstream introns of circRNAs
#{{{

#  plus strand 5' introns (upstream):
P<-CIRCS[ strand(CIRCS) %in% '+' ]
IN<-INTRONS[ strand(INTRONS)=='+' ]
end(IN)<-end(IN)+2  #  move the intron end down by 2nts to make sure we hit the 5' donor site even if it is -1nt from the exon 5'
o<-findOverlaps(resize(IN, 1, fix='end'), resize(P, 2, fix='start'), type='any', select='all')
end(IN)<-end(IN)-2  #  revert back to original ranges
o<-o[ IN[ queryHits(o) ]$gene_id==P[ subjectHits(o) ]$gene_id ]
P.UP<-unlist(GRangesList(lapply(split(queryHits(o), subjectHits(o)), function(n){ IN[n][which.max(width(IN[n]))] })))
P.UP$circ_name<-P[ as.integer(names(P.UP)) ]$circ_name
#
#  plus strand 3' introns (downstream):
IN<-INTRONS[ strand(INTRONS)=='+' ]
start(IN)<-start(IN)-2  #  move the intron start up by 2nts to make sure we hit the 3' acceptor site even if it is -1nt from the exon 3'
o<-findOverlaps(resize(IN, 1, fix='start'), resize(P, 2, fix='end'), type='any', select='all')
start(IN)<-start(IN)+2  #  revert back to original ranges
o<-o[ IN[ queryHits(o) ]$gene_id==P[ subjectHits(o) ]$gene_id ]
P.DOWN<-unlist(GRangesList(lapply(split(queryHits(o), subjectHits(o)), function(n){ IN[n][which.max(width(IN[n]))] })))
P.DOWN$circ_name<-P[ as.integer(names(P.DOWN)) ]$circ_name
#
#  minus strand 5' introns (upstream):
M<-CIRCS[ strand(CIRCS) %in% '-' ]
IN<-INTRONS[ strand(INTRONS)=='-' ]
start(IN)<-start(IN)-2  #  move the intron end down by 2nts to make sure we hit the 5' donor site even if it is -1nt from the exon 5'
o<-findOverlaps(resize(IN, 1, fix='end'), resize(M, 2, fix='start'), type='any', select='all')
start(IN)<-start(IN)+2  #  revert back to original ranges
o<-o[ IN[ queryHits(o) ]$gene_id==M[ subjectHits(o) ]$gene_id ]
M.UP<-unlist(GRangesList(lapply(split(queryHits(o), subjectHits(o)), function(n){ IN[n][which.max(width(IN[n]))] })))
M.UP$circ_name<-M[ as.integer(names(M.UP)) ]$circ_name
#
#  minus strand 3' introns (downstream):
IN<-INTRONS[ strand(INTRONS)=='-' ]
end(IN)<-end(IN)+2  #  move the intron start up by 2nts to make sure we hit the 3' acceptor site even if it is -1nt from the exon 3'
o<-findOverlaps(resize(IN, 1, fix='start'), resize(M, 2, fix='end'), type='any', select='all')
end(IN)<-end(IN)-2  #  revert back to original ranges
o<-o[ IN[ queryHits(o) ]$gene_id==M[ subjectHits(o) ]$gene_id ]
M.DOWN<-unlist(GRangesList(lapply(split(queryHits(o), subjectHits(o)), function(n){ IN[n][which.max(width(IN[n]))] })))
M.DOWN$circ_name<-M[ as.integer(names(M.DOWN)) ]$circ_name


#  merge plus and minus strands back together
CIRCS.introns.up<-c(P.UP, M.UP)
CIRCS.introns.down<-c(P.DOWN, M.DOWN)
names(CIRCS.introns.up)<-NULL
names(CIRCS.introns.down)<-NULL
rm(P.UP, P.DOWN, M.UP, M.DOWN, IN, P, M, o)


#  announce to the world how many circRNAs with both up/down introns and how many with only up or only down
b<-length(intersect(intersect(CIRCS$circ_name, CIRCS.introns.up$circ_name), CIRCS.introns.down$circ_name))
u<-length(setdiff(CIRCS.introns.up$circ_name, CIRCS.introns.down$circ_name))
d<-length(setdiff(CIRCS.introns.down$circ_name, CIRCS.introns.up$circ_name))
n<-length(setdiff(setdiff(CIRCS$circ_name, CIRCS.introns.up$circ_name), CIRCS.introns.down$circ_name))
cat('Out of the', length(CIRCS), 'circRNAs:\n', 
'  ', b, '(', round(100*b/length(CIRCS), 2), '%) had both upstream and downstream introns\n',
'  ', u, '(', round(100*u/length(CIRCS), 2), '%) had upstream only introns\n',
'  ', d, '(', round(100*d/length(CIRCS), 2), '%) had downstream only introns\n',
'  ', n, '(', round(100*n/length(CIRCS), 2), '%) had neither upstream nor downstream introns.\n')
#
# Out of the 5203 circRNAs:
#     5188 ( 99.71 %) had both upstream and downstream introns
#     5 ( 0.1 %) had upstream only introns
#     10 ( 0.19 %) had downstream only introns
#     0 ( 0 %) had neither upstream nor downstream introns.
rm(b,u,d)


#  get the sequences
CIRCS.introns.up.seqs<-setNames(getSeq(BSgenome.Hsapiens.UCSC.hg38, CIRCS.introns.up), CIRCS.introns.up$circ_name)
CIRCS.introns.down.seqs<-setNames(getSeq(BSgenome.Hsapiens.UCSC.hg38, CIRCS.introns.down), CIRCS.introns.down$circ_name)

#}}}


#  [~15min] identify control introns not overlapping with the circRNA introns or the circRNA putative exons
#{{{

CIRCS.controls.introns<-GRangesList()
for (g in seq_along(CIRCS.controls)){
    n<-names(CIRCS.controls[g])
    i<-INTRONS[ INTRONS$gene_id==CIRCS.controls[[g]]$gene_id[1] ]


    #  combine upstream intron, backspliced junction and downstream intron into one range
    r<-reduce(c(CIRCS[ CIRCS$circ_name %in% n ], 
                CIRCS.introns.up[ CIRCS.introns.up$circ_name %in% n],
                CIRCS.introns.down[ CIRCS.introns.down$circ_name %in% n]))
    stopifnot(length(r)==1)


    #  find all the introns of this gene that do not overlap with the reduced range
    o<-findOverlaps(r, i, type='any', select='all')
    i<-i[-subjectHits(o)]

    
    #  find out if there are overlapping introns and take the longest
    o<-findOverlaps(i, i, type='any', select='all')
    if(!all(queryHits(o)==subjectHits(o))){ 
        o<-data.frame(o[ queryHits(o)<subjectHits(o) ])
        w<-width(i)

        #  pairwise comparison of indices, if one intron overlaps many then sequentially the shortest is always discarded, so the same
        #  index might appear multiple times in the list to discard but this is fine
        i<-i[-apply(o, 1, function(x){ x[ which.min(w[x]) ] })]
    }
    if (length(i)==0){
        CIRCS.controls.introns[[n]]<-i
        next
    }


    #  pick as many introns as possible of combined length at least as long as the combined length of the circRNA introns
    w<-sum(width(c(CIRCS.introns.up[ CIRCS.introns.up$circ_name %in% n], CIRCS.introns.down[ CIRCS.introns.down$circ_name %in% n])))
    i.w<-width(i)
    o<-order(i.w, decreasing=T)
    i<-i[o]
    i.w<-cumsum(i.w[o])
    o<-which(i.w>=w)[1]
    if(!is.na(o)){
        i<-i[1:o]
    }
    CIRCS.controls.introns[[n]]<-i
}
rm(g,n,i,r,o,w,i.w)


#  unlist back to GRanges since each intron should have its own sequence and no sequence splicing should happen when we construct the sequences
CIRCS.controls.introns<-unlist(CIRCS.controls.introns)
CIRCS.controls.introns$circ_name<-names(CIRCS.controls.introns)
names(CIRCS.controls.introns)<-NULL


#  define unique intron names by appending an index to the circ_name 
x<-data.table(data.frame(circ_name=CIRCS.controls.introns$circ_name))[, .(n=seq_len(.N)), by=.(circ_name)]
stopifnot(all.equal( CIRCS.controls.introns$circ_name, x$circ_name ) )
CIRCS.controls.introns$intron_name<-x[, paste(circ_name, n, sep='_')]


#  get the sequences
CIRCS.controls.introns.seqs<-setNames(getSeq(BSgenome.Hsapiens.UCSC.hg38, CIRCS.controls.introns), CIRCS.controls.introns$intron_name)

#}}}


#  save
save(CIRCS.introns.up, CIRCS.introns.down, CIRCS.introns.up.seqs, CIRCS.introns.down.seqs, CIRCS.controls.introns, CIRCS.controls.introns.seqs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_introns.RData')
writeXStringSet(CIRCS.introns.up.seqs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_introns_up.fa', format='fasta')
writeXStringSet(CIRCS.introns.down.seqs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_introns_down.fa', format='fasta')
writeXStringSet(CIRCS.controls.introns.seqs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_introns_controls.fa', format='fasta')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_introns.RData
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_introns_up.fa
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_introns_down.fa
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_introns_controls.fa




######################################
#
#
#  basic statistics and other analyses
#
#
######################################




#  [CIRCOS plot + analysis] circRNA per chromosome localizations
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(circlize)
library(ComplexHeatmap)
library(VariantAnnotation)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()



#  [run once] create hg38 cytoband
#             compute circRNA and gene number density per 1M of sequence
#             compute circRNA/mRNA ratio per gene
#             compute mean circular/(1+external junction) ratio per gene and per sample
#{{{

#  load annotation
#  exclude chrM and chrY
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')       
seqlevels(hsa, pruning.mode='coarse')<-setdiff(seqlevels(hsa), c('chrM', 'chrY'))
txs<-hsa[ hsa$type %in% 'transcript' ]
mcols(txs)<-mcols(txs)[, c('gene_id', 'gene_name', 'gene_type', 'transcript_id', 'transcript_name', 'transcript_type')]
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_name', 'gene_type')]


#  GRCh38 cytoBand
#  remove unplaced/decoy contigs, chrM, chrY
#  order them by number (besides chrX) and by start position
#  order chrX separately
cytoband<-fread('/data/genomes/GRCh38/GRCh38.cytoBand.tsv', header=T, sep='\t')
cytoband<-cytoband[ !grepl('_|chrM|chrY', seqnames) ][, index:=as.integer(sub('^chr', '', seqnames))][ order(index, start) ][, index:=NULL]  #  NAs introduced for chrX
cytoband[grep('chrX', seqnames)]<-cytoband[ grep('chrX', seqnames)][ order(seqnames, start) ]
cytoband<-cytoband[, c('start', 'end'):=list(as.numeric(start), as.numeric(end))]
cytoband.xlim<-cytoband[, .(start=min(start), end=max(end)), by=.(seqnames)]


#  load summarized circRNAs 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_summarized.RData')
circs<-l
rm(l, l.table)


#  compute density of circRNAs in 1Mb windows
circ.d<-cytoband[, .(length=tail(end, 1)), by=.(seqnames)]  #  chromosome lengths
circ.d<-circ.d[, .(start=list(seq(1, length, 1e6))), by=.(seqnames, length)]
circ.d<-circ.d[, .(start=unlist(start), end={i<-unlist(start)[-1]-1; append(i, length)}), by=.(seqnames)]
circ.d<-GRanges(seqnames=circ.d$seqnames, ranges=IRanges(start=circ.d$start, end=circ.d$end), n=rep(0, nrow(circ.d)))
ov<-findOverlaps(circs, circ.d, type='within', select='all')
m<-table(subjectHits(ov))
circ.d$n[ as.integer(names(m)) ]<-as.integer(m)
circ.d<-data.frame(circ.d)[, c('seqnames', 'start', 'end', 'n')]
rm(ov, m)


#  count number of genes per Mb per chromosome
gene.d<-cytoband[, .(length=tail(end, 1)), by=.(seqnames)]  #  chromosome lengths
gene.d<-gene.d[, .(start=list(seq(1, length, 1e6))), by=.(seqnames, length)]
gene.d<-gene.d[, .(start=unlist(start), end={i<-unlist(start)[-1]-1; append(i, length)}), by=.(seqnames)]
gene.d<-GRanges(seqnames=gene.d$seqnames, ranges=IRanges(start=gene.d$start, end=gene.d$end), n=rep(0, nrow(gene.d)))
ov<-findOverlaps(hsa, gene.d, type='within', select='all')
m<-table(subjectHits(ov))
gene.d$n[ as.integer(names(m)) ]<-as.integer(m)
gene.d<-data.frame(gene.d)[, c('seqnames', 'start', 'end', 'n')]
rm(ov, m)


#  compute number of circRNAs per gene and number of transcripts per gene
n<-data.table(data.frame(mcols(circs)[, c('gene_id', 'circ_name')]))[, .(n=.N), by=.(gene_id)]
circ.n<-hsa[ hsa$gene_id %in% n$gene_id ]
circ.n$number<-n[ match(circ.n$gene_id, gene_id), n]
n<-data.table(data.frame(mcols(txs)[, c('gene_id', 'transcript_id')]))[, .(n=.N), by=.(gene_id)]
gene.n<-hsa[ hsa$gene_id %in% n$gene_id ]
gene.n$number<-n[ match(gene.n$gene_id, gene_id), n]
rm(n)


#  compute ratio of circRNA number/transcript number for the genes that give rise the circRNA(s)
g<-gene.n[ match(circ.n$gene_id, gene.n$gene_id) ]
stopifnot( all.equal( g$gene_id, circ.n$gene_id) ) 
circ.n$ratio<-circ.n$number/g$number
rm(g)


#  annotate density peaks>10 for each chromosome and save
annot.d<-setNames(vector('list', length(unique(circ.d$seqnames))), unique(circ.d$seqnames))
for (chr in names(annot.d)){
    p<-subset(circ.d, seqnames %in% chr & n>10)
    if(nrow(p)==0){ next }
    p<-GRanges(seqnames=p$seqnames, ranges=IRanges(start=p$start, end=p$end), n=p$n)
    ov<-findOverlaps(circs, p, type='within', select='all')
    ov<-split(queryHits(ov), subjectHits(ov))
    ov<-lapply(ov, function(x){ unique(circs[ x ]$gene_name) })
    p$gene_names<-ov
    annot.d[[chr]]<-p
}
annot.d<-annot.d[ lengths(annot.d)>0 ]
annot.d<-unlist(GRangesList(annot.d))
names(annot.d)<-NULL
annot.d$gene_names<-as(annot.d$gene_names, 'CompressedCharacterList')
rm(chr,p,ov)


#  generate sids out of the bids that will match the genomics ids
circs$sid<-as(lapply(circs$bid, function(s){ sub('-11.*$', '', s) }), 'CompressedList')



#  identify the circRNAs residing on ratio peaks
circs$peaks<-F
circs$peaks[ circs$gene_name %in% circ.n[ circ.n$ratio>1.2]$gene_name ]<-T
table(circs$peaks)
# 
# FALSE  TRUE 
#  4837   366 



#  load circular and linear junction counts
#  convert all NA to zero
#  compute mean ratio of circular/(1+external linear) counts per gene and across samples
#
#  N.B. ignore the warning about "invalid .internal.selfref detected"
#  N.B. some genes with not external junctions available drop out and should become NAs (NOT ZEROS!) 
#
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-lin.cir[ unique(unlist(circs$bid)) ]
r<-unlist(lin.cir)[ circ_name %in% circs$circ_name ]
r[is.na(c.count), c.count:=0]
r[is.na(l.count.out), l.count.out:=0]
r[is.na(l.count.out.max), l.count.out.max:=0]
r<-r[, expression_ratio:=c.count/(1+l.count.out)][, .(expression_ratio=mean(expression_ratio)), by=.(bid, gene_name)][, .(expression_ratio=mean(expression_ratio)), by=.(gene_name)]
circ.n$expression_ratio<-r[ match(circ.n$gene_name, gene_name), expression_ratio]
rm(r, lin.cir)


#  anything interesting?
circ.n[ which(circ.n$expression_ratio>2) ]$gene_name  #  GUSBP1, LINC00632


#  save
save(annot.d, gene.d, circ.d, gene.n, circ.n, circs, cytoband, cytoband.xlim, hsa, txs, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/circos_circRNA_density.RData')

#}}}



#  load back
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/circos_circRNA_density.RData')



#  load all CNVs into a GRanges object
#  load the estimated tumor ploidy 
#  split CNVs per patient
#  complement the CNV gaps with the tumor ploidy as copy-number in order to form full genome copy-number calls
#{{{

#  identify all CNV calls
cn<-system2('find', args=c('/fast/projects/peifer_wgs/work/work/2017-09-25_BerlinWGS_hg38_jt/freec/ -maxdepth 2 -wholename \'*/*bam_CNVs\''), stdout=T)
names(cn)<-sub('^.*/freec/([^/]*)/.*$', '\\1', cn)


#  remove CB2009 which failed in totalRNA-seq
cn<-cn[ setdiff(names(cn), 'CB2009') ]
cn<-cn[ order(names(cn)) ]


#  collect the CNVs
p<-fread(sub('_CNVs$', '_info.txt', cn[1]), header=F, col.names=c('metric', 'value'))[ metric %in% 'Output_Ploidy', as.integer(value)]
cnv<-fread(cn[1], header=F, sep='\t', col.names=c('seqnames', 'start', 'end', 'copies', 'call'))[, call:=NULL][, seqnames:=paste0('chr', seqnames)][, sid:=names(cn[1])][, ploidy:=p]
for(n in 2:length(cn)){
    p<-fread(sub('_CNVs$', '_info.txt', cn[n]), header=F, col.names=c('metric', 'value'))[ metric %in% 'Output_Ploidy', as.integer(value)]
    cnv<-rbind(cnv, fread(cn[n], header=F, sep='\t', col.names=c('seqnames', 'start', 'end', 'copies', 'call'))[, call:=NULL][, seqnames:=paste0('chr', seqnames)][, sid:=names(cn[n])][, ploidy:=p])
}
cnv<-GRanges(seqnames=cnv$seqnames, strand='*', ranges=IRanges(start=cnv$start+1, end=cnv$end), data.frame(cnv[, c('copies', 'sid', 'ploidy'), with=F]))
seqlevels(cnv)<-c(paste0('chr', 1:22), 'chrX')
cnv<-sort(cnv)
rm(cn,n,p)


#  add chromosome lengths info to the GRanges and split by patient
seqlengths(cnv)<-cytoband.xlim[ match(names(seqlengths(cnv)), seqnames), end]
cnv<-GRangesList(split(cnv, cnv$sid))
cnv<-endoapply(cnv, function(x){ 
    g<-gaps(x)
    g<-g[ strand(g) %in% '*' ]
    copies<-rep(unique(x$ploidy), length(g))
    sid<-rep(unique(x$sid), length(g))
    mcols(g)<-DataFrame(copies=copies, sid=sid, ploidy=copies)
    sort(c(x, g))  #  checked it and it works!
})

#}}}



#  load circular and linear junction counts
#  compute log10(1+CPMs) for the external linear junctions to estimate host-gene expression
#
#  N.B. the list of genes with external junctions is a subset of the list of genes with circRNAs
#
#{{{

#  load the junction quantification results and keep only the tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-lin.cir[ unique(unlist(circs$bid)) ]


#  compute the across isoforms mean external linear junction CPMs 
x<-do.call(rbind, lin.cir)[, c('ci.cpm', 'li.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)]
li.cpm<-dcast(x, bid ~ gene_name, value.var='li.cpm', fun.aggregate=mean, na.rm=T)
li.cpm<-t(data.frame(li.cpm[, -1], row.names=li.cpm[, bid], check.names=F))
rm(x)


#  convert bids to sids
colnames(li.cpm)<-sub('-11.*$', '', colnames(li.cpm))
stopifnot( length(setdiff(unique(unlist(circs$sid)), colnames(li.cpm)))==0 )


#  transform to log10(1+CPM)
li.cpm<-log10(1+li.cpm)

#}}}



#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')



#  CIRCOS plot with densities per 1M
#{{{
circos.clear()
circos.par('start.degree'=90, 'gap.degree'=2, 'track.margin'=c(0.01,0.01), 'cell.padding'=c(0.0,0.05,0.00,0.05), 'unit.circle.segments'=500)
circos.initializeWithIdeogram(as.data.frame(cytoband), sort.chr=F)
circos.genomicTrackPlotRegion(circ.d, panel.fun=function(region, value, ...){
    circos.genomicLines(region, value, area=T, baseline=0, col='lightblue4', border='lightblue4', ...)
}, track.height=0.3, bg.border=NA)
circos.text(cytoband.xlim[ seqnames %in% 'chrX', end]*1.2, 28.0, 'circRNAs', sector.index='chrX', facing='bending', cex=1.4, col='lightblue4', niceFacing=T, font=2)
circos.genomicTrackPlotRegion(gene.d, panel.fun=function(region, value, ...){
    circos.genomicLines(region, value, area=T, baseline=0, col='lightgoldenrod3', border='lightgoldenrod3', ...)
}, track.height=0.3, bg.border=NA)
circos.text(cytoband.xlim[ seqnames %in% 'chrX', end]*1.2, 118.0, 'genes', sector.index='chrX', facing='bending', cex=1.4, col='lightgoldenrod3', font=2, niceFacing=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/circos_circRNA_density.svg', width=10, height=10, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  CIRCOS plot with number of circRNAs/transcripts per gene
#  distribution of ratio of circRNA number per gene to transcript number per gene
#{{{

circos.clear()
circos.par('start.degree'=90, 'gap.degree'=2, 'track.margin'=c(0.01,0.01), 'cell.padding'=c(0.0,0.05,0.00,0.05), 'unit.circle.segments'=500, 'canvas.xlim'=c(-0.5, 1.0), 'canvas.ylim'=c(-0.7, 1.0))
circos.initializeWithIdeogram(as.data.frame(cytoband), sort.chr=F)
circos.genomicTrackPlotRegion(data.frame(circ.n)[, c('seqnames', 'start', 'end', 'ratio')], panel.fun=function(region, value, ...){
    circos.genomicLines(region, value, area=T, baseline=0, col='red4', border='red4', ...)
}, track.height=0.2, bg.border=NA)
circos.text(cytoband.xlim[ seqnames %in% 'chrX', end]*1.2, 3.5, 'ratio', sector.index='chrX', facing='bending', cex=1.4, col='red4', niceFacing=T, font=2)
circos.genomicTrackPlotRegion(data.frame(circ.n)[, c('seqnames', 'start', 'end', 'number')], panel.fun=function(region, value, ...){
    circos.genomicLines(region, value, area=T, baseline=0, col='lightblue4', border='lightblue4', ...)
}, track.height=0.2, bg.border=NA)
circos.text(cytoband.xlim[ seqnames %in% 'chrX', end]*1.2, 24.0, 'circRNAs', sector.index='chrX', facing='bending', cex=1.4, col='lightblue4', niceFacing=T, font=2)
circos.genomicTrackPlotRegion(data.frame(gene.n[ gene.n$gene_id %in% circ.n$gene_id ])[, c('seqnames', 'start', 'end', 'number')], panel.fun=function(region, value, ...){
    circos.genomicLines(region, value, area=T, baseline=0, col='lightgoldenrod3', border='lightgoldenrod3', ...)
}, track.height=0.2, bg.border=NA)
circos.text(cytoband.xlim[ seqnames %in% 'chrX', end]*1.2, 118.0, 'mRNAs', sector.index='chrX', facing='bending', cex=1.4, col='lightgoldenrod3', font=2, niceFacing=T)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/circos_circRNA_number_per_gene.svg', width=10, height=10, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  distribution of the ratio of circRNA number per gene / transcript number per gene that produced circRNA(s)
summary(circ.n$ratio)
#
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.006   0.100   0.167   0.281   0.333   5.000 
table(circ.n$ratio>1.2)
# 
# FALSE  TRUE 
#  2257    45 
# 
round(sum(circ.n$ratio>1.2)/length(circ.n)*100, 1)  #  2%
#
par(mar=c(5.0, 8.5, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
b<-pretty(c(min(circ.n$ratio), max(circ.n$ratio)), 5)
h<-hist(circ.n$ratio, breaks=seq(min(b), max(b), length.out=50), plot=F)
h$counts<-log10(1+h$counts)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
h<-plot(h, col='darkgrey', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(b), add=F)
abline(v=1.2, col='red4', lty=2, lwd=4, xpd=F)
mtext('CircRNAs/mRNAs per gene', side=1, line=4, padj=-0.5, las=0, cex=2.4)
mtext(expression(log[10](1+'Frequency')), side=2, line=4, padj=-0.4, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_ratio_circRNA_over_transcript_number_per_gene.svg', width=18, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
rm(b,h,YMAX)

#}}}



#  [circRNA peaks] is there a significant number of involved genes found in COSMIC? 
#{{{


#  Fisher's exact test for circRNAs:
#
#                   |    on peaks      |       off peaks      |
#  -----------------|------------------|----------------------|
#      in COSMIC    |      x1          |          y1          |
#  -----------------|------------------|----------------------|
#    out of COSMIC  |                  |                      | 
#   protein-coding  |      x2          |          y2          | 
#                   |                  |                      | 
#  -----------------|------------------|----------------------|
co<-fread('/fast/groups/ag_schulte/work/reference/cosmic/Cosmic_Gene_Census_Complete_Apr2018.csv', header=T, sep=',')
co<-unique(unlist(co[, 1]))
not_co<-setdiff(hsa[hsa$gene_type %in% 'protein_coding']$gene_name , co)
on_peaks<-circ.n[ circ.n$ratio>1.2 ]$gene_name
off_peaks<-setdiff(circ.n$gene_name, on_peaks)
fisher.test(data.frame('on_peaks'=c(length(intersect(on_peaks, co)), 
                                    length(intersect(on_peaks, not_co))), 
                       'off_peaks'=c(length(intersect(off_peaks, co)), 
                                     length(intersect(off_peaks, not_co)))), alternative='greater')$p.value
#
#  => 0.2014
length(intersect(on_peaks, co))  #  5

#}}}



#  investigate structural variants, CNVs and other demons
#{{{
# 
#  All Berlin cohort analyses (64 CB20.. samples) reside in:
#
#      /fast/projects/peifer_wgs/work/work/2017-09-25_BerlinWGS_hg38_jt/
#
#  We do have whole-genome data of the 56 samples of the Peifer cohort but for those we have only poly-A RNA-Seq.
#
#  The single-nucleotide variants were called using Mutect2. The individual final result files are:
#
#      mutect2/*/somatic_final_ann.vcf.gz 
#
#  They can be imported using VariantAnnotation::readVcf.
#  Alternatively you can load rdata/vcfs.RData which is a list of data.tables with the variant calls.
#
#  The copy-number alteration files are: 
#
#      freec/*_CNVs
#
#  as the tool Control-FREEC was used to determine them. Alternatively you can load rdata/cnas_freec.RData which is a GRangesList with the calls.
#
#  The filtered structural variant calls is located in doc/BerlinWGS_NovobreakResults_allChr_hg38.txt.
#  The additional calls based on hg19 are located under barcelona/.



#  [all regions] compute frequencies of gains and losses on and off peaks
#
#  N.B. irrespective of the tumor ploidy a gain has copies>2 and a loss has copies<2
#  N.B. We merge multiple circRNA isoform hits at the gene level otherwise we are biasing the analysis by isoform counts
#
#{{{
res.circs<-res.linear<-list()
for(s in names(cnv)){

    #  isolate the circs expressed in the tumor
    x<-circs[ sapply(circs$sid, function(x){ any(x %in% s) }) ]  

    #  circRNAs on the peaks
    on.x<-x[ x$peaks ]
    on.hits<-findOverlaps(on.x, cnv[[s]], select='all', type='within')
    #  different isoforms of the same gene most likely will not hit multiple regions, but if they do we count them
    on.s<-unlist(lapply(split( subjectHits(on.hits), on.x[ queryHits(on.hits) ]$gene_name ), unique))
    on<-cnv[[s]][on.s]$copies
    on.neutral<-sum(on==2)/length(on)
    on.gain<-sum(on>2)/length(on)
    on.loss<-sum(on<2)/length(on)


    #  circRNAs off the peaks
    off.x<-x[ !x$peaks ]
    off.hits<-findOverlaps(off.x, cnv[[s]], select='all', type='within')
    #  different isoforms of the same gene most likely will not hit multiple regions, but if they do we count them
    off.s<-unlist(lapply(split( subjectHits(off.hits), off.x[ queryHits(off.hits) ]$gene_name ), unique))
    off<-cnv[[s]][off.s]$copies
    off.neutral<-sum(off==2)/length(off)
    off.gain<-sum(off>2)/length(off)
    off.loss<-sum(off<2)/length(off)


    #  collect the circRNA results
    res.circs[[s]]<-data.frame(on.neutral=on.neutral, off.neutral=off.neutral, on.gain=on.gain, off.gain=off.gain, on.loss=on.loss, off.loss=off.loss)


    #  mRNA CPM expressions
    res.linear[[s]]<-list(on.neutral=li.cpm[ rownames(li.cpm) %in% names(on.s[ on==2 ]), s], 
                          off.neutral=li.cpm[ rownames(li.cpm) %in% names(off.s[ off==2 ]), s],
                          on.gain=li.cpm[ rownames(li.cpm) %in% names(on.s[ on>2 ]), s], 
                          off.gain=li.cpm[ rownames(li.cpm) %in% names(off.s[ off>2 ]), s],
                          on.loss=li.cpm[ rownames(li.cpm) %in% names(on.s[ on<2 ]), s], 
                          off.loss=li.cpm[ rownames(li.cpm) %in% names(off.s[ off<2 ]), s])
}
res.circs<-do.call(rbind, res.circs)
on.gain<-na.omit(unlist(lapply(res.linear, '[[', 'on.gain')))
on.loss<-na.omit(unlist(lapply(res.linear, '[[', 'on.loss')))
on.neutral<-na.omit(unlist(lapply(res.linear, '[[', 'on.neutral')))
off.gain<-na.omit(unlist(lapply(res.linear, '[[', 'off.gain')))
off.loss<-na.omit(unlist(lapply(res.linear, '[[', 'off.loss')))
off.neutral<-na.omit(unlist(lapply(res.linear, '[[', 'off.neutral')))
res.linear<-list(on.neutral=on.neutral, off.neutral=off.neutral, on.gain=on.gain, off.gain=off.gain, on.loss=on.loss, off.loss=off.loss)
rm(on.gain, on.loss, on.neutral, off.gain, off.loss, off.neutral, x, s, on.x, on.hits, on.s, on, off.x, off.hits, off.s, off)


#  [circRNAs] one-sided Mann-Whitney U tests
wilcox.test(x=res.circs$on.neutral, y=res.circs$off.neutral, alternative='less')$p.value  #  0.5659
wilcox.test(x=res.circs$on.gain, y=res.circs$off.gain, alternative='greater')$p.value     #  0.5214
wilcox.test(x=res.circs$on.loss, y=res.circs$off.loss, alternative='less')$p.value        #  0.01069


#  [mRNA CPMs] one-sided Mann-Whitney U tests
wilcox.test(x=res.linear$on.neutral, y=res.linear$off.neutral, alternative='greater')$p.value  #  3.248e-246
wilcox.test(x=res.linear$on.gain, y=res.linear$off.gain, alternative='greater')$p.value        #  6.059e-123
wilcox.test(x=res.linear$on.loss, y=res.linear$off.loss, alternative='greater')$p.value        #  2.671e-39


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  [circRNAs] ecdfs of copy-number neutral, gains and losses
par(mar=c(5.5, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
XLIM<-pretty(range(res.circs), 4)
curve(ecdf(res.circs$on.gain)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), ylab='', xlab='', pch=NA, col='seagreen', lty=1, lwd=12, main='', xlim=range(XLIM), ylim=c(0, 1))
curve(ecdf(res.circs$off.gain)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), pch=NA, col='palegreen3', lty=1, lwd=12, add=T)
curve(ecdf(res.circs$on.loss)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), pch=NA, col='royalblue4', lty=1, lwd=12, add=T)
curve(ecdf(res.circs$off.loss)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), pch=NA, col='steelblue', lty=1, lwd=12, add=T)
curve(ecdf(res.circs$on.neutral)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), pch=NA, col='grey39', lty=1, lwd=12, add=T)
curve(ecdf(res.circs$off.neutral)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), pch=NA, col='darkgrey', lty=1, lwd=12, add=T)
mtext('Frequency', side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.1, cex=2.4, las=0) 
legend('bottomright', legend=c('on loss', 'off loss', 'on gain', 'off gain', 'on neutral', 'off neutral'), col=c('royalblue4', 'steelblue', 'seagreen', 'palegreen3', 'grey39', 'darkgrey'), bty='n', pch=NA, lty=1, lwd=15, cex=1.6, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_CNV_frequencies_on-off_circRNA_peaks.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [circRNAs] boxplots of copy-number neutral, gains and losses
par(mar=c(4.5, 7.0, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-res.circs
B.cl<-c('grey39', 'darkgrey', 'seagreen', 'palegreen3', 'royalblue4', 'steelblue')
YTICK<-pretty(range(B), 5)
plot(0:1, 0:1, xlim=c(0, ncol(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Frequency', side=2, line=4, padj=-0.5, las=0, cex=2.4)
mtext(text=rep(c('on', 'off'), 3), side=1, line=0, at=seq_len(ncol(B)), las=1, padj=+0.0, cex=2.2, col='black')
segments(x0=c(0.6, 2.6, 4.6), x1=c(2.4, 4.4, 6.4), y0=-0.1, lty=1, lwd=2, col='black')
mtext(text=c('neutral', 'gain', 'loss'), side=1, line=2, at=0.5+seq(1, ncol(B), 2), las=1, padj=+0.5, cex=2.4, col='black')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_CNV_frequencies_on-off_circRNA_peaks.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [mRNAs] boxplots of CPMs of gain and losses
par(mar=c(4.5, 7.0, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-res.linear
B.cl<-c('grey39', 'darkgrey', 'seagreen', 'palegreen3', 'royalblue4', 'steelblue')
YTICK<-pretty(range(B), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'CPM')), side=2, line=2, padj=-0.5, las=0, cex=2.4)
mtext(text=rep(c('on', 'off'), 3), side=1, line=0, at=seq_len(length(B)), las=1, padj=+0.5, cex=2.2, col='black')
segments(x0=c(0.6, 2.6, 4.6), x1=c(2.4, 4.4, 6.4), y0=-0.1, lty=1, lwd=2, col='black')
mtext(text=c('neutral', 'gain', 'loss'), side=1, line=2, at=0.5+seq(1, length(B), 2), las=1, padj=+0.5, cex=2.4, col='black')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_CNV_mRNA_CPMs_on-off_circRNA_peaks.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}



#  [excluding the MYCN locus] compute frequencies of gains and losses on and off peaks
#
#  N.B. irrespective of the tumor ploidy a gain has copies>2 and a loss has copies<2
#  N.B. We merge multiple circRNA isoform hits at the gene level otherwise we are biasing the analysis by isoform counts
#
#{{{

#  MYCN: chr2+:15940550-15947007
#
#  we define a generous MYCN locus +-800Kbps that encompasses NBAS, DDX1 etc
mycn<-GRanges(seqnames='chr2', strand='*', ranges=IRanges(start=15140550, end=16747007))


#  compute frequencies of gains and losses on and off peaks
#
#  N.B. irrespective of the tumor ploidy a gain has copies>2 and a loss has copies<2
#  N.B. We merge multiple circRNA isoform hits at the gene level otherwise we are biasing the analysis by isoform counts
#
res.circs<-res.linear<-list()
for(s in names(cnv)){

    #  isolate the circs expressed in the tumor
    x<-circs[ sapply(circs$sid, function(x){ any(x %in% s) }) ]  


    #  remove any overlaps with the MYCN locus
    o<-queryHits(findOverlaps(x, mycn, select='all', type='any'))  #  any overlap should be removed
    if( length(o)>0 ){
        x<-x[ -c(o) ]
    }

    
    #  remove the MYCN locus from CNV calls only if it is within the calls otherwise the large diploid regions around the locus 
    #  will be removed as well and that's too much.
    CNV<-cnv[[s]]
    o<-queryHits(findOverlaps(CNV, mycn, select='all', type='within'))
    if( length(o)>0 ){
        CNV<-CNV[ -c(o) ]
    }


    #  circRNAs on the peaks
    on.x<-x[ x$peaks ]
    on.hits<-findOverlaps(on.x, CNV, select='all', type='within')
    #  different isoforms of the same gene most likely will not hit multiple regions, but if they do we count them
    on.s<-unlist(lapply(split( subjectHits(on.hits), on.x[ queryHits(on.hits) ]$gene_name ), unique))
    on<-CNV[on.s]$copies
    on.neutral<-sum(on==2)/length(on)
    on.gain<-sum(on>2)/length(on)
    on.loss<-sum(on<2)/length(on)


    #  circRNAs off the peaks
    off.x<-x[ !x$peaks ]
    off.hits<-findOverlaps(off.x, CNV, select='all', type='within')
    #  different isoforms of the same gene most likely will not hit multiple regions, but if they do we count them
    off.s<-unlist(lapply(split( subjectHits(off.hits), off.x[ queryHits(off.hits) ]$gene_name ), unique))
    off<-CNV[off.s]$copies
    off.neutral<-sum(off==2)/length(off)
    off.gain<-sum(off>2)/length(off)
    off.loss<-sum(off<2)/length(off)


    #  collect the circRNA results
    res.circs[[s]]<-data.frame(on.neutral=on.neutral, off.neutral=off.neutral, on.gain=on.gain, off.gain=off.gain, on.loss=on.loss, off.loss=off.loss)


    #  mRNA CPM expressions
    res.linear[[s]]<-list(on.neutral=li.cpm[ rownames(li.cpm) %in% names(on.s[ on==2 ]), s], 
                          off.neutral=li.cpm[ rownames(li.cpm) %in% names(off.s[ off==2 ]), s],
                          on.gain=li.cpm[ rownames(li.cpm) %in% names(on.s[ on>2 ]), s], 
                          off.gain=li.cpm[ rownames(li.cpm) %in% names(off.s[ off>2 ]), s],
                          on.loss=li.cpm[ rownames(li.cpm) %in% names(on.s[ on<2 ]), s], 
                          off.loss=li.cpm[ rownames(li.cpm) %in% names(off.s[ off<2 ]), s])
}
res.circs<-do.call(rbind, res.circs)
on.gain<-na.omit(unlist(lapply(res.linear, '[[', 'on.gain')))
on.loss<-na.omit(unlist(lapply(res.linear, '[[', 'on.loss')))
on.neutral<-na.omit(unlist(lapply(res.linear, '[[', 'on.neutral')))
off.gain<-na.omit(unlist(lapply(res.linear, '[[', 'off.gain')))
off.loss<-na.omit(unlist(lapply(res.linear, '[[', 'off.loss')))
off.neutral<-na.omit(unlist(lapply(res.linear, '[[', 'off.neutral')))
res.linear<-list(on.neutral=on.neutral, off.neutral=off.neutral, on.gain=on.gain, off.gain=off.gain, on.loss=on.loss, off.loss=off.loss)
rm(on.gain, on.loss, on.neutral, off.gain, off.loss, off.neutral, x, s, on.x, on.hits, on.s, on, off.x, off.hits, off.s, off)


#  [circRNAs] one-sided Mann-Whitney U tests
wilcox.test(x=res.circs$on.gain, y=res.circs$off.gain, alternative='greater')$p.value     #  0.564
wilcox.test(x=res.circs$on.loss, y=res.circs$off.loss, alternative='less')$p.value        #  0.01311
wilcox.test(x=res.circs$on.neutral, y=res.circs$off.neutral, alternative='less')$p.value  #  0.5888


#  [mRNA CPMs] one-sided Mann-Whitney U tests
wilcox.test(x=res.linear$on.gain, y=res.linear$off.gain, alternative='greater')$p.value        #  8.301e-123
wilcox.test(x=res.linear$on.loss, y=res.linear$off.loss, alternative='greater')$p.value        #  2.085e-39
wilcox.test(x=res.linear$on.neutral, y=res.linear$off.neutral, alternative='greater')$p.value  #  1.325e-246


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  [circRNAs] ecdfs of copy-number neutral, gains and losses
par(mar=c(5.5, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
XLIM<-pretty(range(res.circs), 4)
curve(ecdf(res.circs$on.gain)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), ylab='', xlab='', pch=NA, col='seagreen', lty=1, lwd=12, main='', xlim=range(XLIM), ylim=c(0, 1))
curve(ecdf(res.circs$off.gain)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), pch=NA, col='palegreen3', lty=1, lwd=12, add=T)
curve(ecdf(res.circs$on.loss)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), pch=NA, col='royalblue4', lty=1, lwd=12, add=T)
curve(ecdf(res.circs$off.loss)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), pch=NA, col='steelblue', lty=1, lwd=12, add=T)
curve(ecdf(res.circs$on.neutral)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), pch=NA, col='grey39', lty=1, lwd=12, add=T)
curve(ecdf(res.circs$off.neutral)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), pch=NA, col='darkgrey', lty=1, lwd=12, add=T)
mtext('Frequency', side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.1, cex=2.4, las=0) 
legend('bottomright', legend=c('on loss', 'off loss', 'on gain', 'off gain', 'on neutral', 'off neutral'), col=c('royalblue4', 'steelblue', 'seagreen', 'palegreen3', 'grey39', 'darkgrey'), bty='n', pch=NA, lty=1, lwd=15, cex=1.6, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_CNV_frequencies_on-off_circRNA_peaks_no_MYCN-locus.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [circRNAs] boxplots of copy-number neutral, gains and losses
par(mar=c(4.5, 7.0, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-res.circs
B.cl<-c('grey39', 'darkgrey', 'seagreen', 'palegreen3', 'royalblue4', 'steelblue')
YTICK<-pretty(range(B), 5)
plot(0:1, 0:1, xlim=c(0, ncol(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Frequency', side=2, line=4, padj=-0.5, las=0, cex=2.4)
mtext(text=rep(c('on', 'off'), 3), side=1, line=0, at=seq_len(ncol(B)), las=1, padj=+0.0, cex=2.2, col='black')
segments(x0=c(0.6, 2.6, 4.6), x1=c(2.4, 4.4, 6.4), y0=-0.1, lty=1, lwd=2, col='black')
mtext(text=c('neutral', 'gain', 'loss'), side=1, line=2, at=0.5+seq(1, ncol(B), 2), las=1, padj=+0.5, cex=2.4, col='black')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_CNV_frequencies_on-off_circRNA_peaks_no_MYCN-locus.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  [mRNAs] boxplots of CPMs of gain and losses
par(mar=c(4.5, 7.0, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-res.linear
B.cl<-c('grey39', 'darkgrey', 'seagreen', 'palegreen3', 'royalblue4', 'steelblue')
YTICK<-pretty(range(B), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext(expression(log[10](1+'CPM')), side=2, line=2, padj=-0.5, las=0, cex=2.4)
mtext(text=rep(c('on', 'off'), 3), side=1, line=0, at=seq_len(length(B)), las=1, padj=+0.5, cex=2.2, col='black')
segments(x0=c(0.6, 2.6, 4.6), x1=c(2.4, 4.4, 6.4), y0=-0.1, lty=1, lwd=2, col='black')
mtext(text=c('neutral', 'gain', 'loss'), side=1, line=2, at=0.5+seq(1, length(B), 2), las=1, padj=+0.5, cex=2.4, col='black')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_CNV_mRNA_CPMs_on-off_circRNA_peaks_no_MYCN-locus.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}



#  proximity of circRNA/mRNA ratio > 1.2 peaks to breakpoints
#
#  N.B. the distance of the host gene to the breakpoint is collected to avoid biasing the statistics by the number of circRNA isoforms
#{{{

#  load the translocations only
brp<-fread('/fast/projects/peifer_wgs/work/work/2017-09-25_BerlinWGS_hg38_jt/doc/BerlinWGS_NovobreakResults_allChr_hg38.txt', sep='\t', select=c(1:5, 8))[ ALT %in% '<TRA>' ]


#  split start/end breakpoint positions and by sample
brp.s<-GRanges(seqnames=brp$CHR.A, ranges=IRanges(start=brp$POS.A, end=brp$POS.A), sid=brp$SAMPLE)
brp.s<-split(brp.s, brp.s$sid)
brp.e<-GRanges(seqnames=brp$CHR.B, ranges=IRanges(start=brp$POS.B, end=brp$POS.B), sid=brp$SAMPLE)
brp.e<-split(brp.e, brp.e$sid)
stopifnot(all.equal( lengths(brp.s), lengths(brp.e) ) )
stopifnot(all.equal( names(brp.s), names(brp.e) ) )


#  investigate distances from breakpoints sample by sample
res.circs<-list()
for(s in names(brp.s)){

    #  isolate the circs expressed in the tumor
    x<-circs[ sapply(circs$sid, function(x){ any(x %in% s) }) ]  


    #  isolate the genes on the peaks and compute distance in 1Mb
    on.x<-hsa[ hsa$gene_id %in% x[ x$peaks ]$gene_id ]
    on.dist.s<-setNames(lapply(levels(seqnames(on.x)), function(chr){
        c1<-on.x[ seqnames(on.x) %in% chr ]
        c2<-brp.s[[s]][ seqnames(brp.s[[s]]) %in% chr ]
        abs(apply(expand.grid(start(c1), start(c2)), 1, diff))/1e6
    }), levels(seqnames(on.x)))
    on.dist.e<-setNames(lapply(levels(seqnames(on.x)), function(chr){
        c1<-on.x[ seqnames(on.x) %in% chr ]
        c2<-brp.e[[s]][ seqnames(brp.e[[s]]) %in% chr ]
        abs(apply(expand.grid(start(c1), start(c2)), 1, diff))/1e6
    }), levels(seqnames(on.x)))
    on.dist<-c(unlist(on.dist.s), unlist(on.dist.e))
    

    #  isolate the genes off the peaks and compute distance in 1Mb
    off.x<-hsa[ hsa$gene_id %in% x[ !x$peaks ]$gene_id ]
    off.dist.s<-setNames(lapply(levels(seqnames(off.x)), function(chr){
        c1<-off.x[ seqnames(off.x) %in% chr ]
        c2<-brp.s[[s]][ seqnames(brp.s[[s]]) %in% chr ]
        abs(apply(expand.grid(start(c1), start(c2)), 1, diff))/1e6
    }), levels(seqnames(off.x)))
    off.dist.e<-setNames(lapply(levels(seqnames(off.x)), function(chr){
        c1<-off.x[ seqnames(off.x) %in% chr ]
        c2<-brp.e[[s]][ seqnames(brp.e[[s]]) %in% chr ]
        abs(apply(expand.grid(start(c1), start(c2)), 1, diff))/1e6
    }), levels(seqnames(off.x)))
    off.dist<-c(unlist(off.dist.s), unlist(off.dist.e))
    

    #  collect the distances
    res.circs[[s]]<-list(on=on.dist, off=off.dist)

}
rm(s, x, on.x, on.dist.s, on.dist.e, on.dist, off.x, off.dist.s, off.dist.e, off.dist)


#  collapse all sample distances in one pool we do not care for per sample statistics
on<-unlist(lapply(res.circs, '[[', 'on'))
off<-unlist(lapply(res.circs, '[[', 'off'))
res.circs<-list(on=on, off=off)
rm(on, off)


#  one-sided Mann-Whitney U test
summary(res.circs$on)
#
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    0.01   18.73   47.38   59.17   87.11  225.44 
#
summary(res.circs$off)
#
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     0.0    21.5    50.7    66.1    99.2   247.5 
wilcox.test(x=res.circs$on, y=res.circs$off, alternative='less')$p.value  #  8.094e-47


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  ecdfs
par(mar=c(5.5, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
XLIM<-pretty(range(res.circs), 4)
curve(ecdf(res.circs$on)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), ylab='', xlab='', pch=NA, col='seagreen', lty=1, lwd=12, main='', xlim=range(XLIM), ylim=c(0, 1))
curve(ecdf(res.circs$off)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(nrow(res.circs), 100), pch=NA, col='palegreen3', lty=1, lwd=12, add=T)
mtext('Distance to breakpoint (M)', side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.1, cex=2.4, las=0) 
legend('topleft', legend=c('on peak', 'off peak'), col=c('seagreen', 'palegreen3'), bty='n', pch=NA, lty=1, lwd=15, cex=1.6, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_distance_to_breakpoint_on-off_circRNA_peaks.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  boxplots of copy-number neutral, gains and losses
par(mar=c(2.5, 8.0, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-res.circs
B.cl<-c('seagreen', 'palegreen3')
YTICK<-pretty(range(B), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Distance to breakpoint (M)', side=2, line=5, padj=-0.5, las=0, cex=2.4)
mtext(text=c('on peak', 'off peak'), side=1, line=0, at=seq_len(length(B)), las=1, padj=+0.5, cex=2.4, col='black')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_distance_to_breakpoint_on-off_circRNA_peaks.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/circos_circRNA_density.RData



#  basic statistics 
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(plotrix)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  barplot of number of circRNAs per gene + including all genes without any (unified) circRNA
#{{{
rm(list=ls())


#  load summarized circRNAs
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_summarized.RData')


#  load annotation
#  exclude mitochondrial genes since CIRI2 ignores them
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')       
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_name', 'gene_type')]
hsa<-hsa[ ! seqnames(hsa) %in% 'chrM' ]


#  count all circRNA isoforms per gene and add zeros for all genes not producing circRNAs
#  bin them in specific bins 
x<-c(as.integer(table(l$gene_name, useNA='no')), rep(0, length(setdiff(hsa$gene_name, l$gene_name)))) 
x<-table(cut(x, breaks=c(seq(-1, 5, 1), 70), labels=c(seq(0, 5, 1), '>5'), right=T))
x<-setNames(as.integer(x), names(x))


#  barplot including zero
x11(width=18, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
#  below axis cut
par(fig=c(0.0, 1.0, 0.0, 0.8), mar=c(5.0, 9.0, 0.5, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.5, cex.axis=2.5, new=F)
YTICK<-pretty(c(0, max(x[-1])), 5)  #  compute ticks as if huge zero column did not exist
YMAX<-tail(YTICK, 1)
y<-x
y[1]<-YMAX  #  set for the bottom plot the huge zero column height to be the maximum height
br<-barplot(y, axes=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(br)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
br<-barplot(y, col='grey39', border='white', ylim=c(0, YMAX), axisnames=F, xlab='', ylab='', yaxt='n', xaxt='n', axes=F, add=T)
axis(1, at=br, labels=names(y), line=-1, tick=F, padj=+0.8, cex.axis=2.5, xpd=NA)
YTICK.LABELS<-YTICK
YTICK.LABELS[ length(YTICK.LABELS) ]<-''
axis(2, at=YTICK, labels=YTICK.LABELS, cex.axis=2.5)
axis.break(2, YMAX, style='slash', brw=0.02)
mtext('Number of genes', side=2, line=6, padj=-0.5, at=par('usr')[4]*0.7, las=0, cex=2.5)
mtext('Number of circRNAs per gene', side=1, line=4, padj=-0.4, las=0, cex=2.5)
#  above axis cut
par(fig=c(0.0, 1.0, 0.685, 1.0), mar=c(4.0, 9.0, 1.0, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.5, cex.axis=2.5, new=T)
y<-x
y[2:length(y)]<-NA
YTICK<-c(max(y, na.rm=T)-max(y, na.rm=T)%%100, max(y, na.rm=T)-max(y, na.rm=T)%%100) + c(-100, +100)
YMAX<-tail(YTICK, 1)
br<-barplot(y, axes=F, plot=F)
plot(0:1, 0:1, type='n', ylim=YTICK, xlim=range(br)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
br<-barplot(y-min(YTICK), col='grey39', border='white', ylim=YTICK, axisnames=F, xlab='', ylab='', yaxt='n', xaxt='n', axes=F, offset=min(YTICK), add=T)
YTICK.LABELS<-YTICK
YTICK.LABELS[1]<-''
axis(2, at=YTICK, labels=YTICK.LABELS, lwd.tick=0, cex.axis=2.5)
segments(x0=par('usr')[1], x1=par('usr')[1]-0.2, y0=par('usr')[4], lwd=1)  #  draw the fucking tickmark yourself
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_number_of_circRNAs_per_gene_with_zero.svg', width=14, height=10, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  barplot without zero
x11(width=18, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.0, 9.0, 1.0, 0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, new=F)
YTICK<-pretty(c(0, max(x[-1])), 5)
YMAX<-tail(YTICK, 1)
br<-barplot(x[-1], axes=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(br)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
br<-barplot(x[-1], col='grey39', border='white', ylim=c(0, YMAX), axisnames=F, xlab='', ylab='', yaxt='n', xaxt='n', axes=F, add=T)
axis(1, at=br, labels=names(x[-1]), line=-1, tick=F, padj=+0.8, cex.axis=2.4, xpd=NA)
axis(2, at=YTICK, cex.axis=2.4)
mtext('Number of genes', side=2, line=6, padj=-0.5, las=0, cex=2.4)
mtext('Number of circRNAs per gene', side=1, line=4, padj=-0.4, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_number_of_circRNAs_per_gene.svg', width=14, height=10, bg='white', antialias='subpixel', pointsize=20, family='Arial')  #  leave 14, 10 the axis-labels will be very nice like that
dev.off()

#}}}



#  [number of unified circRNA isoforms] per sample
#                                       per risk group
#                                       per duplicated read group
#                                       per sex group
#                                       per PI class
#                                       per ADRN/MES class
#                                       per MYCN expression class
#                                       per TERT expression class
#{{{
rm(list=ls())


#  load the classes
#{{{

#  load unified circRNAs
#  order by risk group
#  compute number of isoforms per sample
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')
iso<-data.table(data.frame(mcols(CIRCS.all)[, c('bid', 'circ_name')]))[, .(counts=.N), by=.(bid)]
nb.meta<-nb.meta[, risk_group:=factor(risk_group, levels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))][ order(risk_group), ][, risk_group:=as.character(risk_group)]
iso<-iso[ match(nb.meta$bid, bid), ]
stopifnot( all.equal( nb.meta$bid, iso$bid ) )
iso$risk_group<-nb.meta$risk_group



#  load sex determination
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/sex_determination.RData')
sex<-sex[ match(iso$bid, bid), c('bid', 'sex')]
stopifnot(all.equal( sex$bid, iso$bid ))
iso$sex<-sex$sex
rm(sex)



#  load PI class for the tumor samples and annotate the data objects
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/proliferative_index.RData')
iso$pi.class<-tum[ match(iso$bid, bid), PI.class ]
iso$pi.class.col<-tum[ match(iso$bid, bid), PI.class.col ]
rm(cel, cel.tpm, tum.tpm, tum)



#  load ADRN/MES scores and MES classification for the tumor samples and annotate the data objects
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/adrenergic-mesenchymal_score.RData')
iso$adrn.score<-AM.tum[ match(iso$bid, bid), adrn.score ]
iso$mes.score<-AM.tum[ match(iso$bid, bid), mes.score ]
iso$adrn.class<-AM.tum[ match(iso$bid, bid), adrn.class ]
iso$adrn.class.col<-AM.tum[ match(iso$bid, bid), adrn.class.col ]
iso$mes.class<-AM.tum[ match(iso$bid, bid), mes.class ]
iso$mes.class.col<-AM.tum[ match(iso$bid, bid), mes.class.col ]
rm(AM.tum, AM.cel)



#  load percent duplicated read class for the tumor samples and annotate the data objects
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/dup_percentage.RData')
iso$dup<-tum[ match(iso$bid, bid), 100-dedup ]
iso$dup.class<-tum[ match(iso$bid, bid), dup.class ]
iso$dup.class.col<-tum[ match(iso$bid, bid), dup.class.col ]
rm(cel, tum)



#  load MYCN class for the tumor samples and annotate the data objects
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYCN_expression.RData')
iso$mycn.class<-tum[ match(iso$bid, bid), mycn.class ]
iso$mycn.class.col<-tum[ match(iso$bid, bid), mycn.class.col ]
rm(cel, cel.tpm, tum.tpm, tum)



#  load MYC class for the tumor samples and annotate the data objects
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/MYC_expression.RData')
iso$myc.class<-tum[ match(iso$bid, bid), myc.class ]
iso$myc.class.col<-tum[ match(iso$bid, bid), myc.class.col ]
rm(cel, cel.tpm, tum.tpm, tum)



#  load TERT class for the tumor samples and annotate the data objects
load('/fast/work/projects/peifer_wgs/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/TERT_expression.RData')
iso$tert.class<-tum[ match(iso$bid, bid), tert.class ]
iso$tert.class.col<-tum[ match(iso$bid, bid), tert.class.col ]
rm(cel, cel.tpm, tum.tpm, tum)

#}}}



#  circRNA isoforms barplot per sample 
#{{{
x11(width=22, height=12, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(6.5, 9.5, 4.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-setNames(iso$counts, iso$bid)
B.cl<-setNames(nb.meta$col, nb.meta$risk_group)
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
mtext(text=sub('-11-R01', '', names(B)), side=1, line=0, at=bp, las=2, adj=1, cex=1.8, col=B.cl)
mtext('Number of circRNAs', side=2, line=6, padj=-0.7, las=0, cex=2.4)
legend(x=0.5*par('usr')[1], y=par('usr')[4]*1.20, legend=unique(names(B.cl)), col=unique(B.cl), bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.70, x.intersp=0.1, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_circRNAs_per_sample.svg', width=42, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()
#}}}



#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')



#  circRNA isoforms boxplot per risk_group 
#{{{
par(mar=c(9.0, 8.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso$counts, factor(iso$risk_group, levels=unique(iso$risk_group)))
#
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value   #  0.0006354
wilcox.test(x=B[['MNA']], y=B[['LR']], alternative='less')$p.value        #  0.07353
wilcox.test(x=B[['LR']], y=B[['HR_nMNA']], alternative='less')$p.value    #  0.0204
#
B.cl<-setNames(unique(nb.meta$col), unique(nb.meta$risk_group))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=7, padj=+0.4, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot for the PI high class per risk_group 
#{{{
par(mar=c(9.0, 8.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso[ pi.class %in% 'high', counts], factor(iso[ pi.class %in% 'high', risk_group], levels=unique(iso$risk_group)))
#
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value   #  0.004591
wilcox.test(x=B[['MNA']], y=B[['LR']], alternative='less')$p.value        #  0.5064
#
B.cl<-setNames(unique(nb.meta$col), unique(nb.meta$risk_group))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=7, padj=+0.4, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_PI_high_per_risk_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


#  tables
iso[ pi.class %in% 'high', table(risk_group, dup.class)]
#           dup.class
# risk_group low medium high
#    HR_nMNA   9     12    1
#    IMR       4      2    0
#    LR        7      4    3
#    MNA       4     16    2
#    ST4S      0      7    0

#}}}



#  circRNA isoforms boxplot per percent of duplicated reads
#{{{
par(mar=c(0.5, 8.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso$counts, factor(iso$dup.class, levels=levels(iso$dup.class)))
#
#  sequence duplication affects circRNA discovery?
#
wilcox.test(B[['low']], B[['medium']], alternative='greater')$p.value   #  0.04174
wilcox.test(B[['medium']], B[['high']], alternative='greater')$p.value  #  0.003524
#
B.cl<-setNames(levels(iso$dup.class.col), levels(iso$dup.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=7, padj=+0.4, las=0, cex=2.4)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=1, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_per_dup_class.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot for the low percent of duplicated reads class per risk group
#{{{
par(mar=c(9.5, 9.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso[ dup.class %in% 'low', counts], iso[ dup.class %in% 'low', risk_group])
B.cl<-setNames(unique(nb.meta$col), unique(nb.meta$risk_group))
B<-B[ names(B.cl) ]
#
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value  #  0.2231
#
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=7, padj=+0.4, las=0, cex=2.4)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_per_risk_group_dup_low.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot for the medium percent of duplicated reads class per risk group
#{{{
par(mar=c(9.0, 9.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso[ dup.class %in% 'medium', counts], iso[ dup.class %in% 'medium', risk_group])
B.cl<-setNames(unique(nb.meta$col), unique(nb.meta$risk_group))
B<-B[ names(B.cl) ]
#
wilcox.test(x=B[['MNA']], y=B[['HR_nMNA']], alternative='less')$p.value  #  0.002762
#
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=7, padj=+0.4, las=0, cex=2.4)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_per_risk_group_dup_medium.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot per sex group
#{{{
par(mar=c(2.0, 8.5, 0.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso$counts, factor(iso$sex, levels=unique(iso$sex)))
#
#  any sex-bias?
#
wilcox.test(B[[1]], B[[2]], alternative='two.sided')$p.value  #  0.3135
#
B.cl<-setNames(c('cornflowerblue', 'chocolate1'), c('M', 'F'))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=7, padj=+0.4, las=0, cex=2.4)
#mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=2, padj=+0.5, cex=2.4, col=B.cl)
mtext(text=names(B), side=1, line=-1, at=seq(1, length(B), 1), las=0, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_per_sex_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot per PI class
#{{{
par(mar=c(2.0, 9.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso$counts, iso$pi.class)
#
wilcox.test(x=B[['low']], y=B[['high']], alternative='greater')$p.value  #  0.0009645
#
B.cl<-setNames(levels(iso$pi.class.col), levels(iso$pi.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=6, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_per_proliferative_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot per MES class
#{{{
par(mar=c(2.0, 9.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso$counts, iso$mes.class)
#
wilcox.test(x=B[['low']], y=B[['high']], alternative='greater')$p.value  #  0.6252
#
B.cl<-setNames(levels(iso$mes.class.col), levels(iso$mes.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=6, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_per_MES_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot per ADRN class
#{{{
par(mar=c(2.0, 9.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso$counts, iso$adrn.class)
#
wilcox.test(x=B[['low']], y=B[['high']], alternative='less')$p.value  #  0.3859
#
B.cl<-setNames(levels(iso$adrn.class.col), levels(iso$adrn.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=6, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_per_ADRN_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot per MYCN class
#{{{
par(mar=c(2.0, 9.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso$counts, iso$mycn.class)
#
wilcox.test(x=B[['low']], y=B[['high']], alternative='greater')$p.value  #  0.005014
#
B.cl<-setNames(levels(iso$mycn.class.col), levels(iso$mycn.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=6, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_per_MYCN_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot for the high PI class per MYCN class
#{{{
par(mar=c(2.0, 9.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso[ pi.class %in% 'high', counts], iso[ pi.class %in% 'high', mycn.class])
#
wilcox.test(x=B[['low']], y=B[['high']], alternative='greater')$p.value  #  0.0361
#
B.cl<-setNames(levels(iso$mycn.class.col), levels(iso$mycn.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=6, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_PI_high_per_MYCN_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot per MYC class
#{{{
par(mar=c(2.0, 9.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso$counts, iso$myc.class)
#
wilcox.test(x=B[['low']], y=B[['high']], alternative='greater')$p.value  #  0.9192
#
B.cl<-setNames(levels(iso$myc.class.col), levels(iso$myc.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=6, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_per_MYC_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot for the MYCN low class per MYC class 
#{{{
par(mar=c(2.0, 9.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso[ mycn.class %in% 'low', counts], iso[ mycn.class %in% 'low', myc.class])
#
wilcox.test(x=B[['high']], y=B[['low']], alternative='greater')$p.value  #  0.2865 (TESTING FOR MORE circRNAs in the MYC high class)
#
B.cl<-setNames(levels(iso$myc.class.col), levels(iso$myc.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=6, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_MYCN_low_per_MYC_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}



#  circRNA isoforms boxplot per TERT class
#{{{
par(mar=c(2.0, 9.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-split(iso$counts, iso$tert.class)
#
wilcox.test(x=B[['low']], y=B[['high']], alternative='greater')$p.value  #  0.9192
#
B.cl<-setNames(levels(iso$tert.class.col), levels(iso$tert.class))
YTICK<-pretty(c(0, max(ceiling(sapply(B, max)))), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('Number of circRNAs', side=2, line=6, padj=-0.5, las=0, cex=2.4)
mtext(text=paste0(names(B), ' (', lengths(B),')'), side=1, line=-1, at=seq(1, length(B), 1), las=1, padj=+0.5, cex=2.4, col=B.cl)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circRNAs_per_TERT_group.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#}}}

#}}}



#  barplot of number of exons per circRNA
#{{{
rm(list=ls())


#  load circRNA sequences that contain the number of exons as well
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')


#  barplot 
x11(width=50, height=20, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')  #  scaled down SVG for huge fonts
par(mar=c(7.5, 11.0, 5.5, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=4.0, cex.axis=4.0)
x<-table(cut(lengths(CIRCS.exons), breaks=c(seq(0, 15, 1), 37), labels=c(seq(1, 15, 1), '>15'), right=T))
x<-setNames(as.integer(x), names(x))
YTICK<-pretty(c(0, max(x)), 5)
YMAX<-tail(YTICK, 1)
br<-barplot(x, axes=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(0, YMAX), xlim=range(br)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
br<-barplot(x, col='grey39', border='white', ylim=c(0, YMAX), axisnames=F, xlab='', ylab='', yaxt='n', xaxt='n', axes=F, add=T)
axis(1, at=br, labels=names(x), line=0, tick=F, padj=+0.5, cex.axis=3.0, xpd=NA)
axis(2, at=YTICK, labels=YTICK, cex.axis=3.0)
mtext('Number of circRNAs', side=2, line=6, padj=-0.5, las=0, cex=4.0)
mtext('Number of exons per circRNA', side=1, line=6, padj=-0.1, las=0, cex=4.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_number_of_exons_per_circRNA.svg', width=36, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}



#  ecdf of relative exon ranks of donor/acceptor sites
#{{{
rm(list=ls())


#  load circRNA sequences
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')


#  reduce to meta-exons all annotated exons per gene for the genes involved with circRNAs 
#  compute the ranks of the meta-exons
txdb<-loadDb('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/TxDb_GRCh38.gencode.v30.RData')
EXONS<-unlist(endoapply(exonsBy(txdb, 'gene')[ unique(CIRCS$gene_id) ], reduce))
EXONS$gene_id<-names(EXONS)
names(EXONS)<-NULL
EXONS<-data.table(data.frame(EXONS))[, c('exon_rank', 'max_rank'):=list(ifelse(strand=='+', seq_len(.N), rev(seq_len(.N))), rep(.N, .N)), by=.(gene_id)]
EXONS<-GRanges(data.frame(EXONS))
rm(txdb)


#  isolate the circRNA exons that have the donor/acceptor sites
donor<-unlist(endoapply(CIRCS.exons, function(x){ x[1] }))
donor$circ_name<-names(donor)
names(donor)<-NULL
acceptor<-unlist(endoapply(CIRCS.exons, function(x){ tail(x, 1) }))
acceptor$circ_name<-names(acceptor)
names(acceptor)<-NULL


#  find same gene overlaps of exonic parts with the donor/acceptor exons 
#
#  Remember: +-1nts descrepancies were allowed for circRNA exon boundaries 
#
donor.o<-findOverlaps(resize(donor, 2, fix='start'), EXONS, type='any', select='all')
donor.o<-donor.o[ donor$gene_id[ queryHits(donor.o) ] == EXONS$gene_id[ subjectHits(donor.o) ] ]
stopifnot( all.equal( seq_along(donor), queryHits(donor.o) ) )
acceptor.o<-findOverlaps(resize(acceptor, 2, fix='end'), EXONS, type='any', select='all')
acceptor.o<-acceptor.o[ acceptor$gene_id[ queryHits(acceptor.o) ] == EXONS$gene_id[ subjectHits(acceptor.o) ] ]
stopifnot( all.equal( seq_along(acceptor), queryHits(acceptor.o) ) )


#  add the absolute ranks
donor$exon_rank<-EXONS[subjectHits(donor.o)]$exon_rank
acceptor$exon_rank<-EXONS[subjectHits(acceptor.o)]$exon_rank


#  compute the relative ranks
donor.r<-setNames( EXONS[subjectHits(donor.o)]$exon_rank/EXONS[subjectHits(donor.o)]$max_rank, donor$circ_name)
acceptor.r<-setNames( EXONS[subjectHits(acceptor.o)]$exon_rank/EXONS[subjectHits(acceptor.o)]$max_rank, acceptor$circ_name)
stopifnot(all.equal(names(donor.r), names(acceptor.r)))


#  ecdf of relative ranks 
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.0, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
h.d<-curve(ecdf(donor.r)(x), from=0.0, to=1.0, n=min(length(donor.r), 100), ylab='', xlab='', pch=NA, col='darkolivegreen4', lty=1, lwd=12, main='', xlim=c(0.0, 1), ylim=c(0, 1), xaxt='n')
h.a<-curve(ecdf(acceptor.r)(x), from=0.0, to=1.0, n=min(length(acceptor.r), 100), ylab='', xlab='', pch=NA, col='brown2', lty=1, lwd=12, add=T)
curve(punif(x, min=0, max=1, log=F), from=0.0, to=1.0, n=min(length(donor.r), 100), ylab='', xlab='', pch=NA, col='black', lty=3, lwd=8, add=T)
axis(1, at=seq(0.0, 1, 0.2), labels=seq(0.0, 1, 0.2))
mtext('Relative exon rank', side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.1, cex=2.4, las=0) 
legend('topleft', legend=c('donors', 'acceptors', 'uniform'), col=c('darkolivegreen4', 'brown2', 'black'), bty='n', pch=NA, lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.5, seg.len=0.3)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_relative_exon_ranks.svg', width=12, height=12, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}



#  ecdf of circRNA sequence lengths and resized controls
#{{{
rm(list=ls())


#  load circRNA sequences 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')


#  lengths
i<-log10(width(CIRCS.exons.seqs))
j<-log10(width(CIRCS.controls.seqs.resized))


#  ecdfs
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.5, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
XLIM<-pretty(range(c(i,j)), 4)
h.j<-curve(ecdf(j)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(j), 100), ylab='', xlab='', pch=NA, col='darkorange', lty=1, lwd=12, main='', xlim=range(XLIM), ylim=c(0, 1))
h.i<-curve(ecdf(i)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(i), 100), pch=NA, col='blue4', lty=1, lwd=12, add=T)
mtext(expression(log[10]('Sequence length')), side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.1, cex=2.4, las=0) 
legend('topleft', legend=c('circRNA sequences', 'controls'), col=c('blue4', 'darkorange'), bty='n', pch=NA, lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_circRNAs_sequence_lengths.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}



#  ecdf of circRNA flanking intron lengths and unspliced and spliced controls 
#{{{
rm(list=ls())


#  load circRNA introns 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_introns.RData')


#  collapse upstream and downstream flanking intron lengths into a single vector
#  compute cumulative intron length for the controls
i<-log10(c(width(CIRCS.introns.up), width(CIRCS.introns.down)))
j<-log10(width(CIRCS.controls.introns))
k<-data.table(w=width(CIRCS.controls.introns), circ_name=CIRCS.controls.introns$circ_name)[, .(w=sum(w)), by=.(circ_name)][, log10(w)]


#  ecdfs
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
par(mar=c(5.5, 7.5, 1.5, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
XLIM<-pretty(range(c(i,j,k)), 4)
h.k<-curve(ecdf(k)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(k), 100), ylab='', xlab='', pch=NA, col='darkolivegreen', lty=1, lwd=12, main='', xlim=range(XLIM), ylim=c(0, 1))
h.i<-curve(ecdf(i)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(i), 100), pch=NA, col='blue4', lty=1, lwd=12, add=T)
h.j<-curve(ecdf(j)(x), from=XLIM[1], to=tail(XLIM, 1), n=min(length(j), 100), pch=NA, col='darkorange', lty=1, lwd=12, add=T)
mtext(expression(log[10]('Intron length')), side=1, line=3, las=0, padj=+0.5, cex=2.4, las=0)
mtext('Probability', side=2, line=5, padj=-0.1, cex=2.4, las=0) 
legend('topleft', legend=c('circRNA introns', 'controls unspliced', 'controls spliced'), col=c('blue4', 'darkorange', 'darkolivegreen'), bty='n', pch=NA, lty=1, lwd=15, cex=2.0, y.intersp=0.60, x.intersp=0.2, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_lengths_introns.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}



#  correlation of number of exons and number of circRNA isoforms
#{{{
rm(list=ls())


#  load unified circRNAs and order the metadata by risk group
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')


#  compute number of circRNA isoforms per gene and per sample (summing across cases with different gene_id but same gene_name)
circ<-data.table(data.frame(mcols(CIRCS.all)[, c('gene_id', 'gene_name', 'bid')]))[, .(N=.N), by=.(bid, gene_id, gene_name)]
circ<-dcast(circ, bid ~ gene_name, value.var='N', fun.aggregate=sum)


#  load the reference
#  count number of exons per transcript and combine them to lists per gene_name
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'exon' ]
mcols(hsa)<-mcols(hsa)[, c('type', 'gene_id', 'gene_name', 'gene_type', 'transcript_id', 'transcript_type', 'transcript_name', 'exon_number', 'exon_id')]
hsa<-data.table(data.frame(mcols(hsa)))
hsa<-hsa[, .(exons=.N), by=.(gene_name, transcript_id)][, .(transcript_id=list(transcript_id), exons=list(exons)), by=.(gene_name)]


#  maximum number of exons per gene_name ordered by the circRNA gene_names
hsa.max<-hsa[, .(N=max(unlist(exons))), by=.(gene_name)]
hsa.max<-setNames( hsa.max$N, hsa.max$gene_name )
hsa.max<-hsa.max[ colnames(circ[, -1]) ]
stopifnot( all.equal( colnames(circ[, -1]), names(hsa.max) ) )


#  square of Pearson's correlation coefficient per sample between number of circRNAs and maximum number of exons across all isoforms 
#  including only circRNAs found expressed in the sample
r2<-melt(circ, id=1, variable.name='gene_name', value.name='N')[, .(r={keep<-N!=0; cor(N[keep], hsa.max[keep])^2}), by=.(bid)]
r2<-r2[ match(nb.meta$bid, bid), ]
stopifnot( all.equal( r2$bid, nb.meta$bid ) )


#  recycle
x11(width=30, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  tumors
par(mar=c(6.5, 9.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-setNames(r2$r, r2$bid)
stopifnot(all.equal( names(B), nb.meta$bid ))
B.cl<-setNames(unique(nb.meta$col), unique(nb.meta$risk_group))
YTICK<-pretty(c(0.0, 0.1), 5)
bp<-barplot(B, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(YTICK[1], tail(YTICK, 1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col=nb.meta$col, axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext(text=sub('-11-R01', '', names(B)), side=1, line=0, at=bp, las=2, adj=1, cex=1.8, col=nb.meta$col)
mtext(expression(R^2), side=2, line=5, padj=-0.5, las=0, cex=2.4)
legend('topright', legend=names(B.cl), col=B.cl, bty='n', lty=1, lwd=15, pch=NA, cex=2.0, xpd=T, y.intersp=0.70, x.intersp=0.1, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_circRNAs_correlation_max_number_of_exons_per_gene.svg', width=40, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}



#  correlation of number of transcripts per gene and number of circRNAs per gene
#{{{
rm(list=ls())


#  load unified circRNAs and order the metadata by risk group
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')


#  compute number of circRNA isoforms per gene and per sample (summing across cases with different gene_id but same gene_name)
circ<-data.table(data.frame(mcols(CIRCS.all)[, c('gene_id', 'gene_name', 'bid')]))[, .(N=.N), by=.(bid, gene_id, gene_name)]
circ<-dcast(circ, bid ~ gene_name, value.var='N', fun.aggregate=sum)


#  load the reference
#  count number of transcripts per gene
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'transcript' ]
mcols(hsa)<-mcols(hsa)[, c('type', 'gene_id', 'gene_name', 'gene_type', 'transcript_id', 'transcript_type', 'transcript_name')]
hsa<-data.table(data.frame(mcols(hsa)))[, .(transcripts=.N, transcript_id=list(transcript_id)), by=.(gene_name)]
hsa<-hsa[ match(colnames(circ[, -1]), gene_name), ]
stopifnot( all.equal( hsa[, gene_name], colnames(circ[, -1]) ) )


#  number of transcripts per gene_name ordered by the circRNA gene_names
hsa.max<-setNames( hsa$transcripts, hsa$gene_name )


#  square of Pearson's correlation coefficient per sample between number of circRNAs and number of splice-variants based only on the 
#  circRNAs found expressed in the sample
r2<-melt(circ, id=1, variable.name='gene_name', value.name='N')[, .(r={keep<-N!=0; cor(N[keep], hsa.max[keep])^2}), by=.(bid)]
r2<-r2[ match(nb.meta$bid, bid), ]
stopifnot( all.equal( r2$bid, nb.meta$bid ) )


#  recycle
x11(width=30, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  tumors
par(mar=c(6.5, 9.5, 1.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-setNames(r2$r, r2$bid)
stopifnot(all.equal( names(B), nb.meta$bid ))
B.cl<-setNames(unique(nb.meta$col), unique(nb.meta$risk_group))
YTICK<-pretty(c(0.0, 0.1), 5)
bp<-barplot(B, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(YTICK[1], tail(YTICK, 1)), xlim=range(bp)+c(-1, 1), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col=nb.meta$col, axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4, add=T)
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext(text=sub('-11-R01', '', names(B)), side=1, line=0, at=bp, las=2, adj=1, cex=1.8, col=nb.meta$col)
mtext(expression(R^2), side=2, line=5, padj=-0.5, las=0, cex=2.4)
legend('topright', legend=names(B.cl), col=B.cl, bty='n', lty=1, lwd=10, pch=NA, cex=2.0, xpd=T, y.intersp=0.70, x.intersp=0.1, seg.len=0.5)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_circRNAs_correlation_number_of_splice-variants_per_gene.svg', width=40, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}



#  correlation of number of circRNAs per gene with gene expression
#{{{
rm(list=ls())


#  load number of circRNAs per gene estimated when doing the CIRCOS plots
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/circos_circRNA_density.RData')


#  load the junction quantification results and keep only the tumors
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
lin.cir<-lin.cir[ unique(unlist(circs$bid)) ]


#  compute the across isoforms mean external linear junction CPMs 
x<-do.call(rbind, lin.cir)[, c('ci.cpm', 'li.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)]
x[ is.na(li.cpm), li.cpm:=0]
li.cpm<-dcast(x, bid ~ gene_name, value.var='li.cpm', fun.aggregate=mean, na.rm=T)
li.cpm<-t(data.frame(li.cpm[, -1], row.names=li.cpm[, bid], check.names=F))
rm(x)


#  keep only the genes involved in this analysis
length(keep<-intersect( rownames(li.cpm), circs$gene_name ))   #  2221
li.cpm<-li.cpm[ keep, , drop=F]
stopifnot( all.equal( rownames(li.cpm), keep ) )
length(circs<-circs[ circs$gene_name %in% keep ])  #  5109
length(circ.n<-circ.n[ circ.n$gene_name %in% keep ])  #  2221
rm(keep)


#  mean CPM across tumors
li.cpm<-rowMeans(li.cpm)
stopifnot( length(setdiff(names(li.cpm), circ.n$gene_name))==0 )
li.cpm<-li.cpm[ circ.n$gene_name ]
stopifnot( all.equal( names(li.cpm), circ.n$gene_name ) )


#  Spearman correlation of circRNA numbers per gene and mean gene expression

#  circRNAs found expressed in the sample
cor(circ.n$number, li.cpm, method='spearman')  #  0.1467

#}}}



#  percentage of acceptor exons in 3'UTRs and CDSes
#{{{
rm(list=ls())


#  load circRNA sequences
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')


#  import gene_type and transcript_type annotations
#  remove chrM, chrY
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa.gns<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa.gns)<-mcols(hsa.gns)[, c('gene_id', 'gene_name', 'gene_type')]
seqlevels(hsa.gns, pruning.mode='coarse')<-setdiff(seqlevels(hsa.gns), c('chrY', 'chrM'))
hsa.txs<-hsa[ hsa$type %in% 'transcript' ]
mcols(hsa.txs)<-mcols(hsa.txs)[, c('gene_id', 'gene_name', 'gene_type', 'transcript_id', 'transcript_name', 'transcript_type')]
seqlevels(hsa.txs, pruning.mode='coarse')<-setdiff(seqlevels(hsa.txs), c('chrY', 'chrM'))
rm(hsa)


#  collect and annotate the 3'UTRs
#  remove chrM, chrY
#  annotate with transcript_type
txdb<-loadDb('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/TxDb_GRCh38.gencode.v30.RData')
txs<-unlist(transcriptsBy(txdb, 'gene'))
txs$gene_id<-names(txs)
names(txs)<-NULL
utr3<-unlist(threeUTRsByTranscript(txdb, use.names=T))
utr3$transcript_id<-names(utr3)
names(utr3)<-NULL
utr3$gene_id<-txs$gene_id[ match(utr3$transcript_id, txs$tx_name) ]
stopifnot( all(!is.na(utr3$gene_id)) )
mcols(utr3)<-mcols(utr3)[, c('transcript_id', 'gene_id')]
seqlevels(utr3, pruning.mode='coarse')<-setdiff(seqlevels(utr3), c('chrY', 'chrM'))
utr3$transcript_type<-hsa.txs$transcript_type[ match(utr3$transcript_id, hsa.txs$transcript_id) ]
sort(table(utr3$transcript_type), decreasing=T)
# 
# nonsense_mediated_decay          protein_coding  polymorphic_pseudogene               IG_C_gene               TR_C_gene               IG_D_gene 
#                   84275                   66723                      42                      24                       8                       2 
#               IG_V_gene 
#                       2 


#  collect and annotate the CDSes
#  remove chrM, chrY
cd<-unlist(cdsBy(txdb, by='tx', use.names=T))
cd$transcript_id<-names(cd)
names(cd)<-NULL
cd$gene_id<-txs$gene_id[ match(cd$transcript_id, txs$tx_name) ]
stopifnot( all(!is.na(cd$gene_id)) )
mcols(cd)<-mcols(cd)[, c('transcript_id', 'gene_id')]
seqlevels(cd, pruning.mode='coarse')<-setdiff(seqlevels(cd), c('chrY', 'chrM'))
cd$transcript_type<-hsa.txs$transcript_type[ match(cd$transcript_id, hsa.txs$transcript_id) ]
sort(table(cd$transcript_type), decreasing=T)
# 
#          protein_coding nonsense_mediated_decay          non_stop_decay  polymorphic_pseudogene               IG_V_gene               TR_V_gene 
#                  682408                   68630                     580                     326                     286                     210 
#               IG_C_gene               TR_J_gene               IG_D_gene               TR_C_gene               IG_J_gene               TR_D_gene 
#                      97                      79                      37                      21                      18                       4 
rm(txs)


#  keep only protein_coding 3'UTRs and CDSes
#  afterwards make sure to keep the unique list of each to avoid double-counting of hits
utr3<-unique(utr3[ utr3$transcript_type %in% 'protein_coding' ])
cd<-unique(cd[ cd$transcript_type %in% 'protein_coding' ])


#  sanity check: any overlaps need to be other transcripts/isoforms
o<-findOverlaps(utr3, cd, select='all', type='any')
stopifnot( all( utr3[ queryHits(o) ]$transcript_id != cd[ subjectHits(o) ]$transcript_id ) )
rm(o)


#  identify acceptor exons
acceptor<-unlist(endoapply(CIRCS.exons, function(x){ tail(x, 1) }))
acceptor$circ_name<-names(acceptor)
names(acceptor)<-NULL


#  find acceptor exon overlaps with 3'UTRs and CDSes 
#
#  Remember: +-1nts descrepancies were allowed for circRNA exon boundaries 
#
o.utr3<-findOverlaps(resize(acceptor, 2, fix='end'), utr3, type='any', select='all')
o.cd<-findOverlaps(resize(acceptor, 2, fix='end'), cd, type='any', select='all')


#  difference and intersect of the 3'UTR hits with the CDS hits
n.c<-length( setdiff(queryHits(o.cd), queryHits(o.utr3)) )                              #  4865
n.uc<-length( intersect(queryHits(o.utr3), queryHits(o.cd)) )                           #  103
n.u<-length( setdiff(queryHits(o.utr3), queryHits(o.cd)) )                              #  70
n.nuc<-length(setdiff(seq_along(acceptor), union(queryHits(o.cd), queryHits(o.utr3))))  #  165
cat('Percentage of acceptors on CDSes only: ', 100*round(n.c/length(acceptor), digits=3), '%\n', sep='')
cat('Percentage of acceptors on 3\'UTRs and CDSes: ', 100*round(n.uc/length(acceptor), digits=3), '%\n', sep='')
cat('Percentage of acceptors on 3\'UTRs only: ', 100*round(n.u/length(acceptor), digits=3), '%\n', sep='')
cat('Percentage of acceptors on neither 3\'UTRs nor CDSes: ', 100*round(n.nuc/length(acceptor), digits=3), '%\n', sep='')
#
#  =>  Percentage of acceptors on CDSes only:               93.5%
#  =>  Percentage of acceptors on 3'UTRs and CDSes:          2.0%
#  =>  Percentage of acceptors on 3'UTRs only:               1.3%
#  =>  Percentage of acceptors neither on 3'UTRs nor CDSes:  3.2%


#  what are the transcript annotations for those circRNAs with acceptors neither on 3'UTRs nor on CDSes?
x<-acceptor[ setdiff(seq_along(acceptor), union(queryHits(o.cd), queryHits(o.utr3))) ]
x$gene_type<-hsa.gns$gene_type[ match(x$gene_id, hsa.gns$gene_id) ]
sort(table(x$gene_type), decreasing=T)
# 
#                     protein_coding                            lincRNA transcribed_unprocessed_pseudogene                          antisense 
#                                 68                                 34                                 25                                 17 
#               processed_transcript             unprocessed_pseudogene      bidirectional_promoter_lncRNA     transcribed_unitary_pseudogene 
#                                 15                                  4                                  1                                  1 


#  protein-coding gene annotations MUST involve non-protein-coding transcripts or 5'UTRs!!!
x<-x[ x$gene_type %in% 'protein_coding' ]
y<-data.table(data.frame(mcols(hsa.txs[ hsa.txs$gene_id %in% x$gene_id ])))[, .(tx_types=list(unique(unlist(transcript_type)))), by=.(gene_id)]
x$tx_types<-as(y[ match(x$gene_id, gene_id), tx_types], 'CompressedList')
x[sapply(x$tx_types, function(s){ length(setdiff(s, 'protein_coding'))==0 }) ]  #  any acceptor with only protein-coding annotation?
#
#  ENSG00000138696.10|BMPR1B_chr4+94875831-94996134 => it is a spliced 5'UTR with circularizing exons!

#}}}

#}}}



#  [circARID1A, tumors] is circARID1A/ARID1A expression affected by copy number loss of the locus?
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(circlize)
library(ComplexHeatmap)
library(VariantAnnotation)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  GRCh38 cytoBand
#  remove unplaced/decoy contigs, chrM, chrY
#  order them by number (besides chrX) and by start position
#  order chrX separately
cytoband<-fread('/data/genomes/GRCh38/GRCh38.cytoBand.tsv', header=T, sep='\t')
cytoband<-cytoband[ !grepl('_|chrM|chrY', seqnames) ][, index:=as.integer(sub('^chr', '', seqnames))][ order(index, start) ][, index:=NULL]  #  NAs introduced for chrX
cytoband[grep('chrX', seqnames)]<-cytoband[ grep('chrX', seqnames)][ order(seqnames, start) ]
cytoband<-cytoband[, c('start', 'end'):=list(as.numeric(start), as.numeric(end))]
cytoband.xlim<-cytoband[, .(start=min(start), end=max(end)), by=.(seqnames)]


#  load all CNVs into a GRanges object
#  load the estimated tumor ploidy 
#  split CNVs per patient
#  complement the CNV gaps with the tumor ploidy as copy-number in order to form full genome copy-number calls
#{{{

#  identify all CNV calls
cn<-system2('find', args=c('/fast/projects/peifer_wgs/work/work/2017-09-25_BerlinWGS_hg38_jt/freec/ -maxdepth 2 -wholename \'*/*bam_CNVs\''), stdout=T)
names(cn)<-sub('^.*/freec/([^/]*)/.*$', '\\1', cn)


#  remove CB2009 which failed in totalRNA-seq
cn<-cn[ setdiff(names(cn), 'CB2009') ]
cn<-cn[ order(names(cn)) ]


#  collect the CNVs
p<-fread(sub('_CNVs$', '_info.txt', cn[1]), header=F, col.names=c('metric', 'value'))[ metric %in% 'Output_Ploidy', as.integer(value)]
cnv<-fread(cn[1], header=F, sep='\t', col.names=c('seqnames', 'start', 'end', 'copies', 'call'))[, call:=NULL][, seqnames:=paste0('chr', seqnames)][, sid:=names(cn[1])][, ploidy:=p]
for(n in 2:length(cn)){
    p<-fread(sub('_CNVs$', '_info.txt', cn[n]), header=F, col.names=c('metric', 'value'))[ metric %in% 'Output_Ploidy', as.integer(value)]
    cnv<-rbind(cnv, fread(cn[n], header=F, sep='\t', col.names=c('seqnames', 'start', 'end', 'copies', 'call'))[, call:=NULL][, seqnames:=paste0('chr', seqnames)][, sid:=names(cn[n])][, ploidy:=p])
}
cnv<-GRanges(seqnames=cnv$seqnames, strand='*', ranges=IRanges(start=cnv$start+1, end=cnv$end), data.frame(cnv[, c('copies', 'sid', 'ploidy'), with=F]))
seqlevels(cnv)<-c(paste0('chr', 1:22), 'chrX')
cnv<-sort(cnv)
rm(cn,n,p)


#  add chromosome lengths info to the GRanges and split by patient
seqlengths(cnv)<-cytoband.xlim[ match(names(seqlengths(cnv)), seqnames), end]
cnv<-GRangesList(split(cnv, cnv$sid))
cnv<-endoapply(cnv, function(x){ 
    g<-gaps(x)
    g<-g[ strand(g) %in% '*' ]
    copies<-rep(unique(x$ploidy), length(g))
    sid<-rep(unique(x$sid), length(g))
    mcols(g)<-DataFrame(copies=copies, sid=sid, ploidy=copies)
    sort(c(x, g))  #  checked it and it works!
})

#}}}


#  identify the tumors with 1p36 loss
g1p36<-cytoband[ seqnames %in% 'chr1' & grepl('^p36', name) ][, .(start=min(start)+1, end=max(end)), by=.(seqnames)]
g1p36<-GRanges(seqnames='chr1', strand='*', ranges=IRanges(start=g1p36$start, end=g1p36$end), name='1p36')
tum.1p36<-GRangesList()
for(n in seq_along(cnv)){
    s<-subjectHits(findOverlaps(g1p36, cnv[[n]], select='all', type='any'))
    tum.1p36[[names(cnv)[n]]]<-cnv[[n]][s]
}
tum.1p36<-data.table(data.frame(unlist(tum.1p36)))[, 'strand':=NULL]
tum.1p36[ copies<2 ]


#  identify the tumors with ARID1A loss
arid1a<-GRanges(seqnames='chr1', strand='*', ranges=IRanges(start=26693236, end=26782104), name='ARID1A')
tum.arid1a<-GRangesList()
for(n in seq_along(cnv)){
    s<-subjectHits(findOverlaps(arid1a, cnv[[n]], select='all', type='any'))
    tum.arid1a[[names(cnv)[n]]]<-cnv[[n]][s]
}
tum.arid1a<-data.table(data.frame(unlist(tum.arid1a)))[, 'strand':=NULL]
tum.arid1a[ copies<2 ]



#  load circular and linear junction counts
#  keep only the tumors with CNV calls
#  isolate circARID1A
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta<-meta[, sid:=sub('-11.*$', '', bid)][ sid %in% names(cnv) ]
meta<-meta[, risk_group:=factor(risk_group, levels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))][ order(risk_group), ][, risk_group:=as.character(risk_group)]
lc<-unlist(lin.cir[ meta$bid ])[ circ_name %in% 'ENSG00000117713.20|ARID1A_chr1+26729651-26732792' ]
rm(lin.cir, meta.prefailed)


#  convert all NA to zero
lc[is.na(c.count), c.count:=0]
lc[is.na(l.count.out), l.count.out:=0]
lc[is.na(l.count.out.max), l.count.out.max:=0]

#  compute ratio of circular/(1+external linear) counts
#
#  N.B. ignore the warning about "invalid .internal.selfref detected"
#
lc[, ratio:=c.count/(1+l.count.out)]  


#  compute circular junction CPM and external linear junction CPM
lc[, c('c.cpm', 'l.cpm'):=list(c.count/total*1e6, l.count.out/total*1e6)]


#  add risk_group and sid
lc$risk_group<-meta$risk_group[ match( lc$bid , meta$bid ) ]
lc$sid<-meta$sid[ match( lc$bid , meta$bid ) ]


#  add ARID1A locus region copy number
#  pool into diploid, gain, loss groups for higher statistical power
lc$copies<-tum.arid1a$copies[ match( lc$sid , tum.arid1a$sid ) ]
lc[ copies<2 , cnv:='loss']
lc[ copies==2 , cnv:='diploid']
lc[ copies>2 , cnv:='gain']


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  [circARID1A] boxplot of CPMs split by CNV type
B<-split(lc[, c.cpm], factor(lc[, cnv], levels=c('diploid', 'gain', 'loss')))
B.cl<-setNames(c('cadetblue4', 'darkcyan', 'darkolivegreen'), names(B))
#
wilcox.test(x=B[['diploid']], y=B[['loss']], alternative='greater')$p.value  #  0.0009984445355
wilcox.test(x=B[['diploid']], y=B[['gain']], alternative='less')$p.value     #  0.3702023464
#
par(mar=c(2.5, 6.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(range(B), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('CPM', side=2, line=4, padj=-0.1, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=seq_along(B), las=1, padj=+0.5, cex=2.4, col='black')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_circARID1A_CPM_by_CNV_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  [ARID1A] boxplot of CPMs split by CNV type
B<-split(lc[, l.cpm], factor(lc[, cnv], levels=c('diploid', 'gain', 'loss')))
B.cl<-setNames(c('cadetblue4', 'darkcyan', 'darkolivegreen'), names(B))
#
wilcox.test(x=B[['diploid']], y=B[['loss']], alternative='greater')$p.value  #  0.04829219285
wilcox.test(x=B[['diploid']], y=B[['gain']], alternative='less')$p.value     #  0.1028524791
#
par(mar=c(2.5, 6.5, 0.1, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
YTICK<-pretty(range(B), 5)
plot(0:1, 0:1, xlim=c(0, length(B))+c(0.5, 0.5), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
bp<-boxplot(B, col=B.cl, ylab='', xlab='', show.names=F, frame.plot=F, medcol='cyan', boxwex=0.8, xpd=F, outline=F, boxcol=B.cl, range=0, add=T)
mtext('CPM', side=2, line=4, padj=-0.1, las=0, cex=2.4)
mtext(text=names(B), side=1, line=0, at=seq_along(B), las=1, padj=+0.5, cex=2.4, col='black')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/boxplot_ARID1A_CPM_by_CNV_tumors.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')

#}}}




##################################
#
#
#  neuroblastoma specific circRNAs
#
#
##################################




#  identify neuroblastoma-specific circRNAs by comparing neuroblastoma to brain tissue and various tumors and looking for those that are significantly
#  differently expressed in both comparisons:
#
#      we identify the genes in the top 200 range of circRNA expression by summing CPMs across isoforms within 
#      samples and taking the mean CPM across samples
#
#      we do a one-sided Mann-Whitney test for each circRNA isoform comparing CPMs between NB and brain tissue and NB and various tumors asking
#      for the circRNA isoform to have CPM>=1 in at least 30% of the samples where it is expressed
#
#      after FDR correcting the p-values we keep cases which are significant in both comparisons (p-value<0.001) and we sum the two p-values to produce
#      a single p-value per circRNA isoform 
#
#      cirRNAs are ordered by increasing p-value, lifted over to hg19 and assigned corresponding circBase ID(s) if available
#
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load unified set of circRNAs
#  N.B. CPMs are computed based on the total number of feature counts per sample 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')


#  identify top 500 circRNA producing genes 
#  we sum CPMs across circRNA isoforms per sample and then take the mean across samples identifying the top 500 
g<-data.table(data.frame(mcols(nb.circs)[, c('gene_name', 'bid', 'cpm')]))[, .(cpm=sum(cpm)), by=.(bid, gene_name)][, .(cpm=mean(cpm)), by=.(gene_name)][ order(-cpm) ][1:500, gene_name]


#  all circRNA isoforms associated with these genes
CIRCS<-CIRCS[ CIRCS$gene_name %in% g ]
nb.circs<-nb.circs[ nb.circs$gene_name %in% g ]
nb.genes<-nb.genes[, g]
hb.circs<-hb.circs[ hb.circs$gene_name %in% g ]
hb.genes<-hb.genes[, g]
vt.circs<-vt.circs[ vt.circs$gene_name %in% g ]
vt.genes<-vt.genes[, g]
rm(g)


#  [a bit slow] go over each circRNA isoform and do Mann-Whitney tests with the alternative hypothesis being that NB is higher expressed
CIRCS$pv.hb<-CIRCS$pv.vt<-NA
NB<-unique(nb.circs$bid)
HB<-unique(hb.circs$bid)
VT<-unique(vt.circs$bid)
for (n in seq_along(CIRCS)){
    circ<-CIRCS[n]

    #  find the samples that express the circRNA isoform
    x<-mcols(nb.circs[ nb.circs$circ_name %in% circ$circ_name ])[, c('bid', 'cpm')]
    y<-mcols(hb.circs[ hb.circs$circ_name %in% circ$circ_name ])[, c('bid', 'cpm')]
    z<-mcols(vt.circs[ vt.circs$circ_name %in% circ$circ_name ])[, c('bid', 'cpm')]


    #  include the rest of the samples with zero CPM
    x<-rbind(x, DataFrame(bid=setdiff(NB, x$bid), cpm=rep(0.0, length(setdiff(NB, x$bid)))))
    y<-rbind(y, DataFrame(bid=setdiff(HB, y$bid), cpm=rep(0.0, length(setdiff(HB, y$bid)))))
    z<-rbind(z, DataFrame(bid=setdiff(VT, z$bid), cpm=rep(0.0, length(setdiff(VT, z$bid)))))


    #  we ask for at least 30% of the samples with the isoform expressed to have CPM>=1
    if(quantile(x$cpm[x$cpm!=0], 0.7)<1){ next }  


    #  Mann-Whitney test 
    p.xy<-wilcox.test(x=x$cpm, y=y$cpm, alternative='greater')$p.value
    p.xz<-wilcox.test(x=x$cpm, y=z$cpm, alternative='greater')$p.value
    CIRCS$pv.hb[n]<-p.xy 
    CIRCS$pv.vt[n]<-p.xz 
}
CIRCS<-CIRCS[ !is.na(CIRCS$pv.hb) ]  #  if one of them is not set then the other is not set either
rm(NB,HB,VT,n,x,y,z,circ,p.xy,p.xz)


#  FDR-adjust p-values
CIRCS$pv.hb<-p.adjust(CIRCS$pv.hb, method='BH')
CIRCS$pv.vt<-p.adjust(CIRCS$pv.vt, method='BH')


#  keep only circRNA isoforms that passed in both accounts
CIRCS<-CIRCS[ CIRCS$pv.hb<0.001 & CIRCS$pv.vt<0.001 ]


#  define a p-value which is the sum of p-values since if either test fails the circRNA is not included
CIRCS$pvalue<-CIRCS$pv.hb + CIRCS$pv.vt 


#  order them by p-value 
CIRCS<-CIRCS[ order(CIRCS$pvalue, decreasing=F) ]


#  keep only the circRNA isoforms that made it throughout
nb.circs<-nb.circs[ queryHits( findOverlaps( nb.circs, CIRCS, type='equal', select='all' ) ) ]
hb.circs<-hb.circs[ queryHits( findOverlaps( hb.circs, CIRCS, type='equal', select='all' ) ) ]
vt.circs<-vt.circs[ queryHits( findOverlaps( vt.circs, CIRCS, type='equal', select='all' ) ) ]


#  aggregate them into matrix
nb.circs<-data.table(data.frame(mcols(nb.circs)[, c('bid', 'circ_name', 'cpm')]))
hb.circs<-data.table(data.frame(mcols(hb.circs)[, c('bid', 'circ_name', 'cpm')]))
vt.circs<-data.table(data.frame(mcols(vt.circs)[, c('bid', 'circ_name', 'cpm')]))
nb.circs<-dcast(nb.circs, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)  #  sum-aggregate converts NAs to zeros
nb.circs<-t(as.matrix(data.frame(nb.circs[, -1], row.names=nb.circs$bid, check.names=F)))
hb.circs<-dcast(hb.circs, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)  #  sum-aggregate converts NAs to zeros
hb.circs<-t(as.matrix(data.frame(hb.circs[, -1], row.names=hb.circs$bid, check.names=F)))
vt.circs<-dcast(vt.circs, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)  #  sum-aggregate converts NAs to zeros
vt.circs<-t(as.matrix(data.frame(vt.circs[, -1], row.names=vt.circs$bid, check.names=F)))


#  add rows of zeros for circRNA isoforms not found in brain tissue or various tumors and order final matrices according to neuroblastoma order
n<-setdiff( rownames(nb.circs), rownames(hb.circs) )
hb.circs<-rbind(hb.circs, matrix(0, nrow=length(n), ncol=ncol(hb.circs), dimnames=list(n, colnames(hb.circs))))
hb.circs<-hb.circs[ rownames(nb.circs), ]
stopifnot( all.equal( rownames(nb.circs), rownames(hb.circs) ) )
n<-setdiff( rownames(nb.circs), rownames(vt.circs) )
vt.circs<-rbind(vt.circs, matrix(0, nrow=length(n), ncol=ncol(vt.circs), dimnames=list(n, colnames(vt.circs))))
vt.circs<-vt.circs[ rownames(nb.circs), ]
stopifnot( all.equal( rownames(nb.circs), rownames(vt.circs) ) )
stopifnot( length(setdiff(rownames(hb.circs), rownames(vt.circs)))==0 )  #  no circRNA isoform difference between brain tissue and various tumors
rm(n)


#  order them accordingly
stopifnot( length(setdiff(CIRCS$circ_name, rownames(nb.circs)))==0 )
nb.circs<-nb.circs[ CIRCS$circ_name, ]
hb.circs<-hb.circs[ CIRCS$circ_name, ]
vt.circs<-vt.circs[ CIRCS$circ_name, ]


#  manually lift the coordinates over to hg19 and query circBase
#
#  N.B. circBase uses wrongly 0-based format for the genome browser position
#
#  create BED format (0-based):
#
cat(paste(seqnames(CIRCS), start(CIRCS)-1, end(CIRCS), CIRCS$circ_name, 0, strand(CIRCS), sep='\t'), file='~/Downloads/circs.bed', sep='\n')
#
#  upload to:
#
#      https://genome.ucsc.edu/cgi-bin/hgLiftOver
#
#  and convert, then rename to 'circs_hg19.bed' and upload to circBase:
#
#      http://www.circbase.org
#
#  and do a list-search and download the CSV results


#  load the hg19 lifOvers conversions
CIRCS.hg19<-import('~/Downloads/circs_hg19.bed')
mcols(CIRCS.hg19)<-mcols(CIRCS.hg19)[, c('name')]
colnames(mcols(CIRCS.hg19))<-c('circ_name')


#  load the circBase results
#  combine circ_id of identical circRNAs into one 
#  convert to GRanges (MAKE SURE TO CONVERT THEM TO 1-BASED)
#
cb<-fread('~/Downloads/circBase_export.csv', sep=',', header=F, select=c(1:3, 11), col.names=c('position', 'strand', 'circ_id', 'gene_name'))
cb<-cb[, .(circ_id=list(circ_id)), by=.(position, strand, gene_name)]
cb<-cb[, c('seqnames', 'start', 'end'):=list(sub('^([^:]*):([^-]*)-([0-9]*)$', '\\1', cb$position), 
           1+as.integer(sub('^([^:]*):([^-]*)-([0-9]*)$', '\\2', cb$position)),  #  1-based
           as.integer(sub('^([^:]*):([^-]*)-([0-9]*)$', '\\3', cb$position)))]
cb<-GRanges(data.frame(cb)[, c('seqnames', 'start', 'end', 'strand')], gene_name=cb$gene_name, circ_id=as(cb$circ_id, 'CompressedList'))


#  MANUALLY FIX: end position of hsa_circ_0024339 (CADM1) to 115111140
end(cb[ cb$gene_name %in% 'CADM1' & sapply(cb$circ_id, function(i){ any(i %in% 'hsa_circ_0024339') }) ])<-115111140


#  check which circRNAs either have different name in circBase or are not present
print(g<-setdiff(sub('^[^|]+\\|([^_]+)_chr.*$', '\\1', CIRCS.hg19$circ_name), cb$gene_name))  #  Warning: isoforms are collapsed to one name
#
#  => AL121768.1, FIRRE, FMN1, LINC00632
#
ci<-CIRCS.hg19[ grep(paste(g, collapse='|'), CIRCS.hg19$circ_name ) ]
ov<-findOverlaps(ci, cb, type='equal', select='all')
cb[ subjectHits(ov) ]$gene_name<-sub('^[^|]+\\|([^_]+)_chr.*$', '\\1', ci$circ_name[ queryHits(ov) ])  #  this treats isoforms correctly 


#  MANUALLY ADD: FMN1 which is not present in circBase
ci[ setdiff(seq_along(ci), queryHits(ov)) ]$circ_name
# 
#  => ENSG00000248905.8|FMN1_chr15-32964107-32969477
#
cb<-append(cb, GRanges(data.frame(CIRCS.hg19[grep('FMN1', CIRCS.hg19$circ_name)])[, -6], circ_id='hsa_circ_NA', gene_name='FMN1'))
rm(ov, ci)


#  we are ready now to identify the circBase ids
ov<-findOverlaps(CIRCS.hg19, cb, type='equal', select='all')
stopifnot(all.equal(seq_along(CIRCS.hg19), queryHits(ov)))
mcols(CIRCS.hg19)<-cbind(mcols(CIRCS.hg19), mcols(cb[subjectHits(ov)]))
rm(ov, cb)


#  add circ_id to the hg38 circs for easy use
stopifnot( all.equal(CIRCS$circ_name, CIRCS.hg19$circ_name) )
CIRCS$circ_id<-CIRCS.hg19$circ_id


#  keep gene expression for only the relevant genes
nb.genes<-nb.genes[, colnames(nb.genes) %in% CIRCS$gene_name ]
hb.genes<-hb.genes[, colnames(hb.genes) %in% CIRCS$gene_name ]
vt.genes<-vt.genes[, colnames(vt.genes) %in% CIRCS$gene_name ]


#  save everything
save(CIRCS, CIRCS.hg19, hb.circs, hb.genes, hb.meta, nb.circs, nb.genes, nb.meta, vt.circs, vt.genes, vt.meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_neuroblastoma-specific.RData')


#  save also in Excel
x<-data.frame(mcols(CIRCS)[, c('circ_name', 'gene_name', 'gene_id', 'circ_id')])
x$circ_id<-sapply(x$circ_id, paste, collapse=', ')
write.xlsx(x, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_neuroblastoma-specific.xlsx', col.names=T, row.names=F, sheetName='NB-specific circRNAs', append=F)
rm(x)

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_neuroblastoma-specific.{RData,xlsx}




##############################
#
#
#  scrap/obsolete code section
#
#
##############################




#  [SHEP, SKNAS cell models] pick circRNA candidates for the MYCN induction/reversal qPCR panel
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load junction counts and metadata and keep SHEP cell model
load('/data/sequencing/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_linear_vs_circular_collected_results.RData')
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
meta<-rbind(meta.tum, meta.cel, fill=T)
meta.shep<-meta[ !(failed) & grepl('CB-MYCN-SHEP_TR', bid), ]
meta.sknas<-meta[ !(failed) & grepl('CB-SKNAS-TR-MYCN', bid), ]
shep<-unlist(lin.cir[ intersect(names(lin.cir), meta.shep$bid) ])
sknas<-unlist(lin.cir[ intersect(names(lin.cir), meta.sknas$bid) ])
rm(lin.cir, meta, meta.prefailed)


#  convert all NA to zero
shep[is.na(c.count), c.count:=0]
shep[is.na(l.count.out), l.count.out:=0]
shep[is.na(l.count.out.max), l.count.out.max:=0]
sknas[is.na(c.count), c.count:=0]
sknas[is.na(l.count.out), l.count.out:=0]
sknas[is.na(l.count.out.max), l.count.out.max:=0]


#  expression summaries of validated circRNAs by us
ours<-c('ENSG00000079785.15|DDX1_chr2+15603192-15607313',
        'ENSG00000198162.12|MAN1A2_chr1+117402186-117420649',
        'ENSG00000131016.17|AKAP12_chr6+151348711-151353752',
        'ENSG00000117713.20|ARID1A_chr1+26729651-26732792',
        'ENSG00000110422.12|HIPK3_chr11+33286413-33287511',
        'ENSG00000183576.13|SETD3_chr14-99458279-99465813',
        'ENSG00000168769.13|TET2_chr4+105233897-105237351',
        'ENSG00000104313.19|EYA1_chr8-71299047-71334174',
        'ENSG00000196782.12|MAML3_chr4-139889357-139890967',
        'ENSG00000196418.12|ZNF124_chr1-247159006-247159813',
        'ENSG00000072736.19|NFATC3_chr16+68121987-68126610',
        'ENSG00000153147.6|SMARCA5_chr4+143543509-143543972')


#  [SHEP] summarize circRNA expression across samples and define low, medium, and high expression ranges
shep.s<-shep[, .(c.count=mean(c.count), l.count.out=mean(l.count.out)), by=.(circ_name)]
shep.s[, quantile(100*c.count/max(c.count), probs=seq(0.1, 1.0, 0.1))]
#             10%             20%             30%             40%             50%             60%             70%             80%             90% 
#   0.00000000000   0.03993610224   0.11980830671   0.21964856230   0.35942492013   0.51916932907   0.77875399361   1.19808306709   2.37619808307 
#            100% 
# 100.00000000000 
#
shep.s[, r:=cut(100*c.count/max(c.count), breaks=c(0.0, 2.0, 8.0, 100), right=F, include.lowest=T, labels=c('low', 'medium', 'high'))]
shep.s[, table(r)]
#
#    low medium   high 
#   4228    471     96 
#
shep.s[ r %in% 'low', summary(c.count)]
#
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#  0.000000  0.500000  2.333333  3.656575  5.666667 16.666667 
#
shep.s[ r %in% 'medium', summary(c.count)]
#
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 16.83333 21.00000 26.33333 30.75088 38.75000 66.16667 
#
shep.s[ r %in% 'high', summary(c.count)]
#
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#  66.83333  81.37500  99.41667 147.05382 170.66667 834.66667 


#  [SKNAS] summarize circRNA expression across samples and define low, medium, and high expression ranges
sknas.s<-sknas[, .(c.count=mean(c.count), l.count.out=mean(l.count.out)), by=.(circ_name)]
sknas.s[, quantile(100*c.count/max(c.count), probs=seq(0.1, 1.0, 0.1))]
#
#             10%             20%             30%             40%             50%             60%             70%             80%             90% 
#   0.04198152813   0.16792611251   0.33585222502   0.50377833753   0.71368597817   1.00755667506   1.42737195634   2.26700251889   3.94626364400 
#            100% 
# 100.00000000000 
#
sknas.s[, r:=cut(100*c.count/max(c.count), breaks=c(0.0, 2.0, 8.0, 100), right=F, include.lowest=T, labels=c('low', 'medium', 'high'))]
sknas.s[, table(r)]
#
#    low medium   high 
#   3731    858    206 
#
sknas.s[ r %in% 'low', summary(c.count)]
#
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000000 0.3333333 0.9166667 1.2222148 1.8333333 3.9166667 
#
sknas.s[ r %in% 'medium', summary(c.count)]
#
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#  4.000000  5.000000  6.250000  7.258256  8.916667 15.833333 
#
sknas.s[ r %in% 'high', summary(c.count)]
#
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#  16.00000  20.02083  26.87500  37.71036  42.12500 198.50000 


#  check ranges of our circRNAs
shep.s[ circ_name %in% ours ]
sknas.s[ circ_name %in% ours ]


#  [SHEP+SKNAS] identify circRNAs that are commonly found in the low, medium and high ranges
low<-sknas.s[ r %in% 'low' ][ shep.s[ r %in% 'low' ], on='circ_name' ][ c.count>0 & i.c.count>0.0 ][ order(-c.count, -i.c.count, -l.count.out, -i.l.count.out) ]
medium<-sknas.s[ r %in% 'medium' ][ shep.s[ r %in% 'medium' ], on='circ_name' ][ c.count>0 & i.c.count>0.0 ][ order(-c.count, -i.c.count, -l.count.out, -i.l.count.out) ]
high<-sknas.s[ r %in% 'high' ][ shep.s[ r %in% 'high' ], on='circ_name' ][ c.count>0 & i.c.count>0.0 ][ order(-c.count, -i.c.count, -l.count.out, -i.l.count.out) ]


#  check our circRNAs
low[ circ_name %in% ours ]
#                                             circ_name     c.count l.count.out   r   i.c.count i.l.count.out i.r
# 1: ENSG00000131016.17|AKAP12_chr6+151348711-151353752 1.833333333 259.0902778 low 2.333333333   437.3055556 low
# 2:     ENSG00000079785.15|DDX1_chr2+15603192-15607313 1.000000000 415.1811328 low 1.833333333   310.5445419 low
#
medium[ circ_name %in% ours ]
# Empty data.table (0 rows and 7 cols): circ_name,c.count,l.count.out,r,i.c.count,i.l.count.out...
high[ circ_name %in% ours ]
#                                             circ_name      c.count   l.count.out    r    i.c.count i.l.count.out  i.r
# 1: ENSG00000198162.12|MAN1A2_chr1+117402186-117420649 159.16666667  489.95441919 high 834.66666667  310.84116162 high
# 2:   ENSG00000110422.12|HIPK3_chr11+33286413-33287511 147.16666667   95.30845588 high 380.66666667   65.66633987 high
# 3:  ENSG00000072736.19|NFATC3_chr16+68121987-68126610  69.08333333   31.01311469 high 371.50000000   65.67880755 high
# 4: ENSG00000153147.6|SMARCA5_chr4+143543509-143543972  55.33333333 1349.16319444 high  73.66666667  384.24728261 high
# 5:   ENSG00000183576.13|SETD3_chr14-99458279-99465813  43.33333333  241.15681818 high  68.66666667  157.54722222 high

#}}}



#  DOCK8 isoforms
#{{{
rm(list=ls())
library(rtracklayer)
hsa<-import('/data/annotation/GRCh38/GRCh38.gencode.v27.gtf')
hsa<-hsa[ hsa$type %in% 'transcript' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name', 'transcript_id', 'transcript_type', 'transcript_name')]
dock8<-hsa[ hsa$gene_name %in% 'DOCK8' ]
rm(hsa)
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/kallisto.RData')
kal<-unlist(totalrna[ nb.meta$bid ])
kal$bid<-sub('\\.[0-9]*$', '', rownames(kal))
rownames(kal)<-NULL
kal<-data.table(kal)
kal<-kal[ transcript_id %in% dock8$transcript_id ]
d<-dcast(kal, bid ~ transcript_id, value.var='tpm', fun.aggregate=sum)
sort(colMeans(d[, -1]), decreasing=T)
#}}}



#  HR_nMNA-specific circRNAs compared to brain tissue and various tumors 
#{{{
rm(list=ls())
library(Biostrings)
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load unified set of circRNAs
#  N.B. CPMs are computed based on the total number of feature counts per sample 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues.RData')


#  keep only HR_nMNA samples + associated circRNAs
nb.circs<-nb.circs[ nb.circs$risk_group %in% 'HR_nMNA' ]
CIRCS<-CIRCS[ CIRCS$circ_name %in% nb.circs$circ_name ]
hb.circs<-hb.circs[ hb.circs$circ_name %in% nb.circs$circ_name ]
vt.circs<-vt.circs[ vt.circs$circ_name %in% nb.circs$circ_name ]
GENES<-GENES[ GENES$gene_id %in% nb.circs$gene_id , ]
nb.genes<-nb.genes[ , colnames(nb.genes) %in% GENES$gene_name ]
hb.genes<-hb.genes[ , colnames(hb.genes) %in% GENES$gene_name ]
vt.genes<-vt.genes[ , colnames(vt.genes) %in% GENES$gene_name ]


#  identify top 500 circRNA producing genes 
#  we sum CPMs across circRNA isoforms per sample and then take the mean across samples identifying the top 500 
g<-data.table(data.frame(mcols(nb.circs)[, c('gene_name', 'bid', 'cpm')]))[, .(cpm=sum(cpm)), by=.(bid, gene_name)][, .(cpm=mean(cpm)), by=.(gene_name)][ order(-cpm) ][1:500, gene_name]


#  all circRNA isoforms associated with these genes
CIRCS<-CIRCS[ CIRCS$gene_name %in% g ]
nb.circs<-nb.circs[ nb.circs$gene_name %in% g ]
nb.genes<-nb.genes[, g]
hb.circs<-hb.circs[ hb.circs$gene_name %in% g ]
hb.genes<-hb.genes[, g]
vt.circs<-vt.circs[ vt.circs$gene_name %in% g ]
vt.genes<-vt.genes[, g]
rm(g)


#  [a bit slow] go over each circRNA isoform and do Mann-Whitney tests with the alternative hypothesis being that NB is higher expressed
CIRCS$pv.hb<-CIRCS$pv.vt<-NA
NB<-unique(nb.circs$bid)
HB<-unique(hb.circs$bid)
VT<-unique(vt.circs$bid)
for (n in seq_along(CIRCS)){
    circ<-CIRCS[n]

    #  find the samples that express the circRNA isoform
    x<-mcols(nb.circs[ nb.circs$circ_name %in% circ$circ_name ])[, c('bid', 'cpm')]
    y<-mcols(hb.circs[ hb.circs$circ_name %in% circ$circ_name ])[, c('bid', 'cpm')]
    z<-mcols(vt.circs[ vt.circs$circ_name %in% circ$circ_name ])[, c('bid', 'cpm')]


    #  include the rest of the samples with zero CPM
    x<-rbind(x, DataFrame(bid=setdiff(NB, x$bid), cpm=rep(0.0, length(setdiff(NB, x$bid)))))
    y<-rbind(y, DataFrame(bid=setdiff(HB, y$bid), cpm=rep(0.0, length(setdiff(HB, y$bid)))))
    z<-rbind(z, DataFrame(bid=setdiff(VT, z$bid), cpm=rep(0.0, length(setdiff(VT, z$bid)))))


    #  we ask for at least 30% of the samples with the isoform expressed to have CPM>=1
    if(quantile(x$cpm[x$cpm!=0], 0.7)<1){ next }  


    #  Mann-Whitney test 
    p.xy<-wilcox.test(x=x$cpm, y=y$cpm, alternative='greater')$p.value
    p.xz<-wilcox.test(x=x$cpm, y=z$cpm, alternative='greater')$p.value
    CIRCS$pv.hb[n]<-p.xy 
    CIRCS$pv.vt[n]<-p.xz 
}
CIRCS<-CIRCS[ !is.na(CIRCS$pv.hb) ]  #  if one of them is not set then the other is not set either
rm(NB,HB,VT,n,x,y,z,circ,p.xy,p.xz)


#  FDR-adjust p-values
CIRCS$pv.hb<-p.adjust(CIRCS$pv.hb, method='BH')
CIRCS$pv.vt<-p.adjust(CIRCS$pv.vt, method='BH')


#  keep only circRNA isoforms that passed in both accounts
CIRCS<-CIRCS[ CIRCS$pv.hb<0.001 & CIRCS$pv.vt<0.001 ]


#  define a p-value which is the sum of p-values since if either test fails the circRNA is not included
CIRCS$pvalue<-CIRCS$pv.hb + CIRCS$pv.vt 


#  order them by p-value 
CIRCS<-CIRCS[ order(CIRCS$pvalue, decreasing=F) ]


#  keep only the circRNA isoforms that made it throughout
nb.circs<-nb.circs[ queryHits( findOverlaps( nb.circs, CIRCS, type='equal', select='all' ) ) ]
hb.circs<-hb.circs[ queryHits( findOverlaps( hb.circs, CIRCS, type='equal', select='all' ) ) ]
vt.circs<-vt.circs[ queryHits( findOverlaps( vt.circs, CIRCS, type='equal', select='all' ) ) ]


#  aggregate them into matrix
nb.circs<-data.table(data.frame(mcols(nb.circs)[, c('bid', 'circ_name', 'cpm')]))
hb.circs<-data.table(data.frame(mcols(hb.circs)[, c('bid', 'circ_name', 'cpm')]))
vt.circs<-data.table(data.frame(mcols(vt.circs)[, c('bid', 'circ_name', 'cpm')]))
nb.circs<-dcast(nb.circs, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)  #  sum-aggregate converts NAs to zeros
nb.circs<-t(as.matrix(data.frame(nb.circs[, -1], row.names=nb.circs$bid, check.names=F)))
hb.circs<-dcast(hb.circs, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)  #  sum-aggregate converts NAs to zeros
hb.circs<-t(as.matrix(data.frame(hb.circs[, -1], row.names=hb.circs$bid, check.names=F)))
vt.circs<-dcast(vt.circs, bid ~ circ_name, value.var='cpm', fun.aggregate=sum)  #  sum-aggregate converts NAs to zeros
vt.circs<-t(as.matrix(data.frame(vt.circs[, -1], row.names=vt.circs$bid, check.names=F)))


#  add rows of zeros for circRNA isoforms not found in brain tissue or various tumors and order final matrices according to neuroblastoma order
n<-setdiff( rownames(nb.circs), rownames(hb.circs) )
hb.circs<-rbind(hb.circs, matrix(0, nrow=length(n), ncol=ncol(hb.circs), dimnames=list(n, colnames(hb.circs))))
hb.circs<-hb.circs[ rownames(nb.circs), ]
stopifnot( all.equal( rownames(nb.circs), rownames(hb.circs) ) )
n<-setdiff( rownames(nb.circs), rownames(vt.circs) )
vt.circs<-rbind(vt.circs, matrix(0, nrow=length(n), ncol=ncol(vt.circs), dimnames=list(n, colnames(vt.circs))))
vt.circs<-vt.circs[ rownames(nb.circs), ]
stopifnot( all.equal( rownames(nb.circs), rownames(vt.circs) ) )
stopifnot( length(setdiff(rownames(hb.circs), rownames(vt.circs)))==0 )  #  no circRNA isoform difference between brain tissue and various tumors
rm(n)


#  order them accordingly
stopifnot( length(setdiff(CIRCS$circ_name, rownames(nb.circs)))==0 )
nb.circs<-nb.circs[ CIRCS$circ_name, ]
hb.circs<-hb.circs[ CIRCS$circ_name, ]
vt.circs<-vt.circs[ CIRCS$circ_name, ]


#  manually lift the coordinates over to hg19 and query circBase
#
#  N.B. circBase uses wrongly 0-based format for the genome browser position
#
#  create BED format (0-based):
#
cat(paste(seqnames(CIRCS), start(CIRCS)-1, end(CIRCS), CIRCS$circ_name, 0, strand(CIRCS), sep='\t'), file='~/Downloads/circs.bed', sep='\n')
#
#  upload to:
#
#      https://genome.ucsc.edu/cgi-bin/hgLiftOver
#
#  and convert, then rename to 'circs_hg19.bed' and upload to circBase:
#
#      http://www.circbase.org
#
#  and do a list-search and download the CSV results


#  load the hg19 lifOvers conversions
CIRCS.hg19<-import('~/Downloads/circs_hg19.bed')
mcols(CIRCS.hg19)<-mcols(CIRCS.hg19)[, c('name')]
colnames(mcols(CIRCS.hg19))<-c('circ_name')


#  load the circBase results
#  combine circ_id of identical circRNAs into one 
#  convert to GRanges (MAKE SURE TO CONVERT THEM TO 1-BASED)
#
cb<-fread('~/Downloads/circBase_export.csv', sep=',', header=F, select=c(1:3, 11), col.names=c('position', 'strand', 'circ_id', 'gene_name'))
cb<-cb[, .(circ_id=list(circ_id)), by=.(position, strand, gene_name)]
cb<-cb[, c('seqnames', 'start', 'end'):=list(sub('^([^:]*):([^-]*)-([0-9]*)$', '\\1', cb$position), 
           1+as.integer(sub('^([^:]*):([^-]*)-([0-9]*)$', '\\2', cb$position)),  #  1-based
           as.integer(sub('^([^:]*):([^-]*)-([0-9]*)$', '\\3', cb$position)))]
cb<-GRanges(data.frame(cb)[, c('seqnames', 'start', 'end', 'strand')], gene_name=cb$gene_name, circ_id=as(cb$circ_id, 'CompressedList'))


#  MANUALLY FIX: end position of hsa_circ_0024339 (CADM1) to 115111140
end(cb[ cb$gene_name %in% 'CADM1' & sapply(cb$circ_id, function(i){ any(i %in% 'hsa_circ_0024339') }) ])<-115111140


#  check which circRNAs either have different name in circBase or are not present
print(g<-setdiff(sub('^[^|]+\\|([^_]+)_chr.*$', '\\1', CIRCS.hg19$circ_name), cb$gene_name))  #  Warning: isoforms are collapsed to one name
#
#  => AL121768.1, FIRRE
#
ci<-CIRCS.hg19[ grep(paste(g, collapse='|'), CIRCS.hg19$circ_name ) ]
ov<-findOverlaps(ci, cb, type='equal', select='all')
cb[ subjectHits(ov) ]$gene_name<-sub('^[^|]+\\|([^_]+)_chr.*$', '\\1', ci$circ_name[ queryHits(ov) ])  #  this treats isoforms correctly 


#  MANUALLY ADD: FMN1 which is not present in circBase
#ci[ setdiff(seq_along(ci), queryHits(ov)) ]$circ_name
# 
#  => ENSG00000248905.8|FMN1_chr15-32964107-32969477
#
#cb<-append(cb, GRanges(data.frame(CIRCS.hg19[grep('FMN1', CIRCS.hg19$circ_name)])[, -6], circ_id='hsa_circ_NA', gene_name='FMN1'))
rm(ov, ci)


#  we are ready now to identify the circBase ids
ov<-findOverlaps(CIRCS.hg19, cb, type='equal', select='all')
stopifnot(all.equal(seq_along(CIRCS.hg19), queryHits(ov)))
mcols(CIRCS.hg19)<-cbind(mcols(CIRCS.hg19), mcols(cb[subjectHits(ov)]))
rm(ov, cb)


#  add circ_id to the hg38 circs for easy use
stopifnot( all.equal(CIRCS$circ_name, CIRCS.hg19$circ_name) )
CIRCS$circ_id<-CIRCS.hg19$circ_id


#  keep gene expression for only the relevant genes
nb.genes<-nb.genes[, colnames(nb.genes) %in% CIRCS$gene_name ]
hb.genes<-hb.genes[, colnames(hb.genes) %in% CIRCS$gene_name ]
vt.genes<-vt.genes[, colnames(vt.genes) %in% CIRCS$gene_name ]


#  save everything
save(CIRCS, CIRCS.hg19, hb.circs, hb.genes, hb.meta, nb.circs, nb.genes, nb.meta, vt.circs, vt.genes, vt.meta, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_HR_nMNA-specific.RData')


#  save also in Excel
x<-data.frame(mcols(CIRCS)[, c('circ_name', 'gene_name', 'gene_id', 'circ_id')])
x$circ_id<-sapply(x$circ_id, paste, collapse=', ')
write.xlsx(x, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_across_tissues_HR_nMNA-specific.xlsx', col.names=T, row.names=F, sheetName='HR_nMNA-specific circRNAs', append=F)
rm(x)

#}}}





