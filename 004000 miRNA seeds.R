###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




#################################################################
#
#
#  miRNA seed counting on circRNAs with mutated seeds as controls
#
#
#################################################################




#  count 7mer miRNA seeds and their mutated controls on circRNAs
#  each time we count as one match any overlapping seeds of a given miRNA (same for the controls) 
#{{{
rm(list=ls())
library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)
library(data.table)


#  functions
#{{{

mutate<-function(M, S=3, E=4){
    #  mutate nucleotides in DNAStringSet object

    #  define the set of nucleotides
    SET<-c('A', 'C', 'G', 'T')


    #  extract all at once the nucleotides to be mutated 
    NTS<-do.call(rbind, strsplit(substr(M, S, E), ''))


    #  mutate them and combine them back to full string
    NTS.mut<-apply(NTS, 1, function(x){ paste0(c(sample(setdiff(SET, x[1]), 1), sample(setdiff(SET, x[2]), 1)), collapse='') } )


    #  replace original nucleotides with mutated nucleotides
    x<-sapply(M, toString)
    substr(x, S, E)<-NTS.mut
    M<-setNames(DNAStringSet(x), names(M))


    return(M)
}


countseeds<-function(SEEDS, SEQS){
    #  Wrapper function that runs countseeds.py.
    #  It temporarily saves SEEDS and SEQS as FASTA files.
    #  Calls countseeds.py and waits to finish so it can parse the counts.tsv outputs.
    #  When finished it deletes all temporary files and the output.
    require(Biostrings)
    require(data.table)

    #  create the temporary files
    targets<-tempfile('targets.fa.')
    seeds<-tempfile('seeds.fa.')
    output<-tempfile('counts.tsv.')

    
    # save FASTA
    writeXStringSet(SEQS, targets)
    writeXStringSet(SEEDS, seeds)


    #  call the python script
    r<-system2('countseeds.py', args=paste('--targets', targets, '--seeds', seeds, output), wait=T)

    
    #  parse the result into a data.table
    cnts<-fread(output, header=F, sep='\t', col.names=c('mir', 'circ_name', 'ncount'))

    
    #  delete the temporary files and the output
    system2('rm', args=paste(targets, seeds, output), wait=T)


    return(cnts)
}

#}}}


#  load the circRNA sequences 
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')
rm(list=ls(pattern='control'))


#  2-8nts 7mer seed-sites (reverse-complemented miRNA sequence)
#  2-7nts+A 7mer seed-sites (reverse-complemented miRNA sequence + A at the end)
MIRS<-readDNAStringSet('/data/annotation/GRCh38/hsa.mature.miRBase21.fa')
m.28<-reverseComplement(subseq(MIRS, start=2, end=8))
x<-reverseComplement(subseq(MIRS, start=2, end=7))
m.a27<-setNames(xscat(x,DNAString('A')), names(x))  #  concatenate two DNAStrings together
rm(x)


#  group all seeds together and create the mutated controls by mutating the 3rd and 4th nucleotides
mirs<-c(m.28, m.a27)
mirs.mut<-mutate(mirs, S=3, E=4)
rm(m.28,m.a27)


#  define miRNA families based on 2-7nts 6mers (reverse-complemented miRNA sequence)
#  add the corresponding miRNA seeds for each miRNA and the corresponding mutated miRNA seeds
fam<-reverseComplement(subseq(MIRS, start=2, end=7))
fam<-split(names(fam), sapply(fam, toString))
fam<-data.table(family.seed=names(fam), family=unname(fam))
fam<-fam[, .(mir=unlist(family), family=family), by=.(family.seed)]
s<-data.table(mir=names(mirs), seed=sapply(mirs, toString))[, .(seed=list(unique(seed))), by=.(mir) ]
fam<-s[ fam, on='mir']
s<-data.table(mir=names(mirs), seed.mut=sapply(mirs.mut, toString))[, .(seed.mut=list(unique(seed.mut))), by=.(mir) ]
fam<-s[ fam, on='mir']
setcolorder(fam, c('mir', 'family', 'family.seed', 'seed', 'seed.mut'))
rm(s)


#  count seed occurrences on circRNAs for both miRNA seeds and controls
CIRCS.sc<-countseeds(mirs, CIRCS.exons.seqs)
CIRCS.controls.sc<-countseeds(mirs.mut, CIRCS.exons.seqs)


#  add (seed, mutated seed, family, family.seed) to each miRNA interaction reported
CIRCS.sc<-fam[ CIRCS.sc, on='mir']
CIRCS.controls.sc<-fam[ CIRCS.controls.sc, on='mir']


#  summarize interactions to families by taking the maximum number of seed occurrences across family members 
#  keeping the unique list of corresponding seeds
CIRCS.sc.fam<-CIRCS.sc[, .(family=unique(family), seed=list(unique(unlist(seed))), seed.mut=list(unique(unlist(seed.mut))), ncount=max(ncount)), by=.(family.seed, circ_name)]
CIRCS.controls.sc.fam<-CIRCS.controls.sc[, .(family=unique(family), seed=list(unique(unlist(seed))), seed.mut=list(unique(unlist(seed.mut))), ncount=max(ncount)), by=.(family.seed, circ_name)]


#  save 
save(CIRCS.sc, CIRCS.sc.fam, CIRCS.controls.sc, CIRCS.controls.sc.fam, fam, mirs, mirs.mut, MIRS, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_hsa.mature.miRBase21_7mer_seed_counts.RData')


#  group transcripts by their miRNA interactions
#  are there any transcripts uniquely targeted?
x<-CIRCS.sc[, .(family=unique(family), circ_name=list(circ_name), seed=list(unique(seed)), seed.mut=list(unique(seed.mut)), ncount=list(unique(ncount))), by=.(family.seed, mir)]
y<-CIRCS.controls.sc[, .(family=unique(family), circ_name=list(circ_name), seed=list(unique(seed)), seed.mut=list(unique(seed.mut)), ncount=list(unique(ncount))), by=.(family.seed, mir)]
print(x[lengths(circ_name)==1]) 
print(y[lengths(circ_name)==1])

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_hsa.mature.miRBase21_7mer_seed_counts.RData



#  isolate the miRNA sponges 
#  recount seeds to non-overlapping seeds across all families involved
#  annotate them with AGO cluster coverage in 41nts windows centered around the seeds
#{{{
rm(list=ls())
library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
source('~/bio/lib/circ2genome.R')


#  load list of circRNAs sequences and spliced exons
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')
rm(list=ls(pattern='control'))


#  load the 7mer counts
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_hsa.mature.miRBase21_7mer_seed_counts.RData')
rm(CIRCS.sc, CIRCS.controls.sc, mirs, mirs.mut, MIRS)
gc()


#  load Huettelmaier miRNA expression
load('/fast/projects/Schulte_NB/work/Tumor_Neuroblastoma_Huettelmaier_Bell_Seq_97/miRNA/known_miRNAs_cpm.RData')
load('/fast/projects/Schulte_NB/work/Tumor_Neuroblastoma_Huettelmaier_Bell_Seq_97/metadata.RData')
hut.meta<-meta
hut.mirs<-data.table(mirs)
rm(meta, mirs)


#  3'UTRs
txdb<-loadDb('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/TxDb_GRCh38.gencode.v30.RData')
utr3<-unlist(threeUTRsByTranscript(txdb, use.names=T))
rm(txdb)


#  annotate 3'UTR-intersecting circRNAs
#  pass the 3'UTR annotation down to the 7mer counts both for signal and control (there can be cases where a circRNA does not have control hits)
ov<-unique(queryHits(findOverlaps(CIRCS, utr3, type='any', select='all')))
CIRCS$utr3<-F
CIRCS$utr3[ ov ]<-T
x<-data.table(as.data.frame(mcols(CIRCS)[, c('circ_name', 'utr3')]))
CIRCS.sc.fam<-x[ CIRCS.sc.fam, on='circ_name' ]
CIRCS.controls.sc.fam<-x[ CIRCS.controls.sc.fam, on='circ_name' ]
setcolorder(CIRCS.sc.fam, c('family.seed', 'family', 'circ_name', 'seed', 'seed.mut', 'ncount', 'utr3'))
setcolorder(CIRCS.controls.sc.fam, c('family.seed', 'family', 'circ_name', 'seed', 'seed.mut', 'ncount', 'utr3'))
rm(x, ov, utr3)


#  load AGO clusters
ago<-import('/fast/projects/Schulte_NB/work/downloads/ago-clip/bwa/sorted.bam.clusters.merged.bed')
seqlevels(ago, pruning.mode='coarse')<-seqlevels(CIRCS)


#  identify putative miRNA sponges:
#  
#      - putative sponges are defined as circRNAs with at least 10 miRNA binding sites for the same family 
#      - we remove miRNA-families with average CPM<5 across the Huettelmaier neuroblastoma tumors
#      - since before we enforced non-overlapping seed sites only within the same family, we collect again all interactions across families 
#        per circRNA isoform and re-identify the number of non-overlapping seed sites and their corresponding AGO coverage in a window of 41nts
#      - we do the same for the mutated seeds, although mutated seeds collected per family are many more compared to the original miRNA seeds per
#        family because in the former case each mutation was independently done for each miRNA member
#      - we summarize at the circRNA-gene level by collecting all circRNA isoforms and taking the maximum number of miRNA seed sites and the maximum
#        AGO coverage across them as representative
#
#  KNOWN BUG: AGO coverage of spliced seeds is not done correctly, fuck them, that would mean GRangesLists and lists of lists just for a single 
#             case in the controls. We drop AGO coverage for those suckers...anyway the controls are super conservative since each of the 
#             two possibilities of an 8mer seed is mutated separately!
#{{{

#  define sponges to have at least 10 miRNA binding sites for a given interaction
#  tabulate how many of them overlap with 3'UTRs
sponges<-CIRCS.sc.fam[ ncount>=10 ][ order(-ncount) ][, gene_name:=sub('^[^|]*\\|(.*)_chr.*$', '\\1', circ_name)]
sponges[, table(utr3, gene_name)]


#  introduce average expression of miRNA family-members across the Huettelmaier neuroblastoma tumors
#  keep interactions with average miRNA CPM>=5
#  order by CPM
sponges<-sponges[, cpm:=round(sapply(family, function(f){ mean(rowMeans(hut.mirs[ mir %in% f , -1])) }), 2)][ cpm>=5 ][ order(-cpm) ]


#  collect all interactions across families per circRNA isoform and sum over average family CPMs
sponges.circ<-sponges[, .(family.seed=list(unique(unlist(family.seed))), family=list(unique(unlist(family))), seed=list(unique(unlist(seed))), seed.mut=list(unique(unlist(seed.mut))), gene_name=unique(gene_name), cpm=sum(cpm)), by=.(circ_name)]


#  collect the circRNA isoform sequences 
sponges.seq<-CIRCS.exons.seqs[ sponges.circ$circ_name ]


#  add circRNA lengths 
stopifnot( all.equal( names(sponges.seq), sponges.circ$circ_name ) )
sponges.circ[, circ_width:=width(sponges.seq) ]


#  count non-overlapping seeds across all families for both miRNA and mutated seeds
#  annotate them with AGO coverage based on centered 41nts windows around the seed
sponges.circ[, c('ncount', 'seeds_hits', 'seeds_start', 'seeds_end', 'ago'):=list(0, list(), list(), list(), list())]
sponges.circ[, c('ncount.mut', 'seeds_hits.mut', 'seeds_start.mut', 'seeds_end.mut', 'ago.mut'):=list(0, list(), list(), list(), list())]
#sponges.circ[, c('ncount', 'seeds_hits', 'seeds_start', 'seeds_end', 'ago'):=list(NULL, NULL, NULL, NULL, NULL)]
#sponges.circ[, c('ncount.mut', 'seeds_hits.mut', 'seeds_start.mut', 'seeds_end.mut', 'ago.mut'):=list(NULL, NULL, NULL, NULL, NULL)]
for(n in 1:nrow(sponges.circ)){
    pd<-PDict(DNAStringSet(sponges.circ[n, seed][[1]]))
    m<-unlist(as(matchPDict(pd, sponges.seq[[ sponges.circ[n, circ_name] ]], min.mismatch=0, with.indels=F), 'CompressedIRangesList'))
    ov<-findOverlaps(m, m, type='any', select='all')
    if(!all(queryHits(ov)==subjectHits(ov))){ 
        ov<-ov[ queryHits(ov)<subjectHits(ov) ]
        m<-m[-subjectHits(ov)]
    }
    if(length(m)>0){
        coordinates<-circ2genome(CIR=m, EXN=CIRCS.exons[[ sponges.circ[n, circ_name] ]], SEQ=extractAt(sponges.seq[[n]], m))
        a<-sapply(coordinates$seqs, toString)
        b<-sapply(getSeq(BSgenome.Hsapiens.UCSC.hg38, coordinates), toString)
        i<-which(a!=b)
        a<-a[i]
        b<-b[i]
        if(length(a)>0){
            cat('signal seeds:', a, 'of circRNA:', sponges.circ[n, circ_name], 'are spliced and their AGO coverage is omitted!\n')
        }
        ov<-data.table(data.frame(findOverlaps(resize(coordinates, width=41, fix='center'), ago, type='any', select='all', ignore.strand=F)))
        colnames(ov)<-c('qh', 'sh')
        ov[, nsum:=ago$score[ sh ]]
        ov<-ov[, .(nsum=sum(nsum)), by=.(qh)]
        coordinates$ago<-0
        coordinates$ago[ ov$qh ]<-ov$nsum
        coordinates$ago[ i ]<-0  #  omit AGO coverage of spliced-seeds
        sponges.circ[n, c('ncount', 'seeds_hits', 'seeds_start', 'seeds_end', 'ago'):=list(length(m), list(sapply(coordinates$seqs, toString)), list(start(coordinates)), list(end(coordinates)), list(coordinates$ago))]
    } 


    #  do the same for the mutated seeds 
    pd<-PDict(DNAStringSet(sponges.circ[n, seed.mut][[1]]))
    m<-unlist(as(matchPDict(pd, sponges.seq[[ sponges.circ[n, circ_name] ]], min.mismatch=0, with.indels=F), 'CompressedIRangesList'))
    ov<-findOverlaps(m, m, type='any', select='all')
    if(!all(queryHits(ov)==subjectHits(ov))){ 
        ov<-ov[ queryHits(ov)<subjectHits(ov) ]
        m<-m[-subjectHits(ov)]
    }
    if(length(m)>0){
        coordinates<-circ2genome(CIR=m, EXN=CIRCS.exons[[ sponges.circ[n, circ_name] ]], SEQ=extractAt(sponges.seq[[n]], m))
        a<-sapply(coordinates$seqs, toString)
        b<-sapply(getSeq(BSgenome.Hsapiens.UCSC.hg38, coordinates), toString)
        i<-which(a!=b)
        a<-a[i]
        b<-b[i]
        if(length(a)>0){
            cat('control seeds:', a, 'of circRNA:', sponges.circ[n, circ_name], 'are spliced and their AGO coverage is omitted!\n')
        }
        ov<-data.table(data.frame(findOverlaps(resize(coordinates, width=41, fix='center'), ago, type='any', select='all', ignore.strand=F)))
        colnames(ov)<-c('qh', 'sh')
        ov[, nsum:=ago$score[ sh ]]
        ov<-ov[, .(nsum=sum(nsum)), by=.(qh)]
        coordinates$ago<-0
        coordinates$ago[ ov$qh ]<-ov$nsum
        coordinates$ago[ i ]<-0  #  omit AGO coverage of spliced-seeds
        sponges.circ[n, c('ncount.mut', 'seeds_hits.mut', 'seeds_start.mut', 'seeds_end.mut', 'ago.mut'):=list(length(m), list(sapply(coordinates$seqs, toString)), list(start(coordinates)), list(end(coordinates)), list(coordinates$ago))]
    } 
}


#  Fisher's exact test of significance of all families that made it into the sponges list against the rest of the families that did not make it
#  both for the signal and control subgroups
#
#  -------------------------------------------------------
#              |  these families  |   rest of families   |
#  -------------------------------------------------------
#  signal      |      x1          |         y1           |
#  -------------------------------------------------------
#  controls    |      x2          |         y2           |
#  -------------------------------------------------------
for(n in 1:nrow(sponges.circ)){
    x1<-sponges.circ[n, ncount]
    x2<-sponges.circ[n, ncount.mut]
    y1<-CIRCS.sc.fam[ circ_name %in% sponges.circ[n, circ_name] & ! family.seed %in% sponges.circ[n, unlist(family.seed)], sum(ncount)]
    y2<-CIRCS.controls.sc.fam[ circ_name %in% sponges.circ[n, circ_name] & ! family.seed %in% sponges.circ[n, unlist(family.seed)], sum(ncount)]
    if(all(!is.na(c(x1,x2,y1,y2)))){
        sponges.circ[n, pvalue:=fisher.test(data.frame('in'=c(x1, x2), 'out'=c(y1, y2)), alternative='greater')$p.value]
    }
}
sponges.circ[, padj:=p.adjust(pvalue, method='BH')]


#  summarize at the gene level by taking the max seed counts across circRNA isoforms and max AGO coverage
#  N.B. since circRNA isoforms overlap in sequence, summing the seed counts would not be appropriate, taking the mean would also bring down
#       the counts when short circRNA isoforms are present, the best would be to take the maximum as representative
sponges.gene<-sponges.circ[, .(circ_name=list(circ_name), 
                               family.seed=list(unique(unlist(family.seed))), 
                               family=list(unique(unlist(family))), 
                               seed=list(unique(unlist(seed))), 
                               seed.mut=list(unique(unlist(seed.mut))), 
                               cpm=sum(cpm), 
                               ncount=max(ncount), 
                               ncount.mut=max(ncount.mut), 
                               ago=max(unlist(ago))), 
                               by=.(gene_name)]


#  save
save(sponges, sponges.circ, sponges.seq, sponges.gene, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_miRNA_sponges.RData')

#}}}

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_miRNA_sponges.RData
#
#  =>      sponges : circRNA isofroms with at least 10 miRNA seed sites for a given miRNA family with miRNA family members mean expression CPM>=5
#  => sponges.circ : per circRNA isoform collection of miRNA family interactions with summing over corresponding mean CPMs, the number of miRNA seeds
#  =>                and mutated miRNA seeds has been recomputed to exclude overlaps of seeds within or across families
#  =>  sponges.seq : per circRNA isoform sequence
#  => sponges.gene : collection of all circRNA isoforms into a circRNA-gene summary with the number of miRNA seeds and mutated miRNA seeds taken 
#                    to be the maximum across the circRNA isoforms



#  miRNA sponges plots
#{{{
rm(list=ls())
library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)
library(data.table)
library(RColorBrewer)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load putative sponges
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_miRNA_sponges.RData')
#
#       sponges : circRNA isofroms with at least 10 miRNA seed sites for a given miRNA family with miRNA family members mean expression CPM>=5
#  sponges.circ : per circRNA isoform collection of miRNA family interactions with summing over corresponding mean CPMs, the number of miRNA seeds
#                 and mutated miRNA seeds has been recomputed to exclude overlaps of seeds across different families
#   sponges.seq : per circRNA isoform sequence
#  sponges.gene : collection of all circRNA isoforms into a circRNA-gene summary with the number of miRNA seeds and mutated miRNA seeds taken 
#                 to be the maximum across the circRNA isoforms


#  keep only the sponges found with significant number of miRNA seeds compared to controls
sponges.circ<-sponges.circ[ padj<0.05 ]
sponges.gene<-sponges.gene[ sapply(circ_name, function(n){ any(n %in% sponges.circ$circ_name)}) ] 
sponges.seq<-sponges.seq[ sponges.circ$circ_name ]
sponges<-sponges[ circ_name %in% sponges.circ$circ_name ]


#  load list of circRNA isoforms and sequences and keep only those found in sponges
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')
CIRCS<-CIRCS[ CIRCS$circ_name %in% unique(sponges$circ_name) ]
CIRCS.exons<-CIRCS.exons[ unique(sponges$circ_name) ]
CIRCS.exons.seqs<-CIRCS.exons.seqs[ unique(sponges$circ_name) ]
rm(EXONS)


#  recycle
x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  histogram of seed counts
par(mar=c(4.5,5.0,0.0,0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=1.8, cex.axis=1.6)
B<-setNames(sponges.gene$ncount, sponges.gene$gene_name)
XLIM<-pretty(c(10, max(B)), 4)
h<-hist(B, breaks=20, col='grey39', border='white', xlab='', ylab='', main='', xlim=range(XLIM), ylim=c(0, 8), add=F)
mtext('Number of miRNA seed sites', side=1, line=3, padj=+0.1, las=0, cex=1.8)
mtext('Number of circRNAs', side=2, line=3, padj=-0.2, las=0, cex=1.8)


#  barplot of seed counts per circRNA-gene showing representative miRNA and CPM as well
par(mar=c(11.5, 8.0, 3.0, 0.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
B<-as.data.frame(sponges.gene[, c('gene_name', 'ncount', 'ago'), with=F][order(-ncount)])
ago.cl<-ifelse(B$ago==0, 'black', 'grey39')
B<-setNames(B$ncount, B$gene_name)
YTICK<-pretty(c(0.0, max(B)), 5)
bp<-barplot(B, beside=F, plot=F)
plot(0:1, 0:1, type='n', ylim=c(YTICK[1], tail(YTICK, 1)), xlim=range(bp)+c(-1, 3), axes=F, ann=F, xaxs='i', yaxs='i')
bp<-barplot(B, border='white', col=ago.cl, axes=F, axisnames=F, beside=F, xlab='', ylab='', las=1, yaxt='n', ylim=c(0, tail(YTICK,1)), add=T)
axis(2, at=YTICK, line=0, cex.axis=2.4)
mtext(text=names(B), side=1, line=0, at=bp, las=2, adj=1, cex=2.4)
mtext('MiRNA seed sites', side=2, line=5, padj=-0.5, las=0, cex=2.4)
m<-sponges.gene[ match(names(B), gene_name), c('gene_name', 'family', 'cpm')][, .(family=sub('hsa-', '', sapply(family, '[[', 1)), cpm=round(cpm)), by=.(gene_name)]
m<-setNames(paste0(m$family, ' (', m$cpm, ')'), m$gene_name)
text(x=bp, y=B, labels=m, pos=4, offset=-0.1, srt=45, cex=1.0)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/barplot_circRNAs_sponges_miRNA_seed_sites.svg', width=18, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()


#  recycle
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')


#  log10(length) distribution of the sponges
par(mar=c(4.5, 7.0, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
L<-log10(width(CIRCS.exons.seqs))
XMAX<-max(L)
XMAX<-tail(pretty(c(1, XMAX), 5), 1)
#  ecdf
h<-curve(ecdf(L)(x), from=1, to=XMAX, n=min(L, 10), ylab='', xlab='', pch=NA, col='grey39', lty=1, lwd=8, main='', ylim=c(0, 1), xlim=c(1, XMAX), xaxt='n')
axis(1, at=pretty(c(1, max(h$x))), labels=pretty(c(1, max(h$x))))
mtext(expression(log[10]~'(Length)'), side=1, line=3, padj=+0.1, las=0, cex=2.4)
mtext('Frequency', side=2, line=4, padj=-0.4, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/ecdf_circRNAs_sponges_length.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#  histogram
XMAX<-pretty(range(L), 5)
b<-seq(XMAX[1], tail(XMAX,1)+0.2, 0.2)
h<-hist(L, breaks=b, plot=F)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
h<-hist(L, breaks=b, col='grey39', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), add=F)
mtext(expression(log[10]~'(Length)'), side=1, line=3, padj=+0.1, las=0, cex=2.4)
mtext('Frequency', side=2, line=3, padj=-0.7, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_circRNAs_sponges_length.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')



#  distribution of number of exons of the sponges
par(mar=c(4.5, 7.0, 0.1, 1.0), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=NA, bty='n', cex.lab=2.4, cex.axis=2.4)
L<-lengths(CIRCS.exons)
h<-hist(L, breaks=10, plot=F)
YMAX<-max(pretty(c(1, max(h$counts)), 4))
h<-hist(L, breaks=h$breaks, col='grey39', border='white', xlab='', ylab='', main='', ylim=c(0, YMAX), xlim=range(h$breaks), add=F, axes=F)
axis(1, at=h$breaks, cex.axis=2.4)
axis(2, at=pretty(c(0, YMAX), 5), cex.axis=2.4)
mtext('Number of exons', side=1, line=3, padj=+0.1, las=0, cex=2.4)
mtext('Frequency', side=2, line=3, padj=-0.7, las=0, cex=2.4)
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/histogram_circRNAs_sponges_number_of_exons.svg', width=14, height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')
dev.off()

#}}}



#  [ARID1A] non-overlapping counts of 7mer miRNA seeds
#{{{
rm(list=ls())
library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
source('~/bio/lib/circ2genome.R')


#  load list of circRNAs sequnces and keep circARID1A only
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.RData')
arid1a<-CIRCS[ grep('ARID1A', CIRCS$circ_name) ]
arid1a.exons<-CIRCS.exons[ grep('ARID1A', names(CIRCS.exons)) ]
arid1a.exons.seqs<-CIRCS.exons.seqs[ grep('ARID1A', names(CIRCS.exons.seqs)) ]
arid1a.exons.controls<-CIRCS.controls[ grep('ARID1A', names(CIRCS.controls)) ]
arid1a.controls.seqs<-CIRCS.controls.seqs[ grep('ARID1A', names(CIRCS.controls.seqs)) ]
arid1a.controls.seqs.resized<-CIRCS.controls.seqs.resized[ grep('ARID1A', names(CIRCS.controls.seqs.resized)) ]
rm(list=ls(pattern='CIRCS|EXONS'))


#  load the 7mer counts and keep only circARID1A related ones
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_hsa.mature.miRBase21_7mer_seed_counts.RData')
arid1a.fam<-CIRCS.sc.fam[ grep('ARID1A', circ_name) ]
arid1a.controls.fam<-CIRCS.controls.sc.fam[ grep('ARID1A', circ_name) ]
rm(list=ls(pattern='CIRCS|^fam$|mirs|MIRS'))
gc()


#  load Huettelmaier miRNA expression
load('/fast/projects/Schulte_NB/work/Tumor_Neuroblastoma_Huettelmaier_Bell_Seq_97/miRNA/known_miRNAs_cpm.RData')
load('/fast/projects/Schulte_NB/work/Tumor_Neuroblastoma_Huettelmaier_Bell_Seq_97/metadata.RData')
hut.meta<-meta
hut.mirs<-data.table(mirs)
rm(meta, mirs)


#  3'UTRs
txdb<-loadDb('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/TxDb_GRCh38.gencode.v30.RData')
utr3<-unlist(threeUTRsByTranscript(txdb, use.names=T))
rm(txdb)


#  annotate 3'UTR-intersecting circRNAs
#  pass the 3'UTR annotation down to the 7mer counts both for signal and control
ov<-unique(queryHits(findOverlaps(arid1a, utr3, type='any', select='all')))
arid1a$utr3<-F
arid1a$utr3[ ov ]<-T
x<-data.table(as.data.frame(mcols(arid1a)[, c('circ_name', 'utr3')]))
arid1a.fam<-x[ arid1a.fam, on='circ_name' ]
arid1a.controls.fam<-x[ arid1a.controls.fam, on='circ_name' ]
setcolorder(arid1a.fam, c('family.seed', 'family', 'circ_name', 'seed', 'seed.mut', 'ncount', 'utr3'))
setcolorder(arid1a.controls.fam, c('family.seed', 'family', 'circ_name', 'seed', 'seed.mut', 'ncount', 'utr3'))
rm(x, ov, utr3)


#  load AGO clusters
ago<-import('/fast/projects/Schulte_NB/work/downloads/ago-clip/bwa/sorted.bam.clusters.merged.bed')
seqlevels(ago, pruning.mode='coarse')<-seqlevels(arid1a)


#  introduce average expression of miRNA family-members across the Huettelmaier neuroblastoma tumors
#  keep interactions with average miRNA CPM>=5
#  order by CPM
arid1a.fam<-arid1a.fam[, cpm:=round(sapply(family, function(f){ mean(rowMeans(hut.mirs[ mir %in% f , -1])) }), 2)][ cpm>=5 ][ order(-cpm) ]
rm(list=ls(pattern='hut'))


#  collect all interactions across families and sum over average family CPMs
arid1a.fam.summarized<-arid1a.fam[, .(family.seed=list(unique(unlist(family.seed))), family=list(unique(unlist(family))), seed=list(unique(unlist(seed))), seed.mut=list(unique(unlist(seed.mut))), cpm=sum(cpm)), by=.(circ_name)]


#  collect the circRNA isoform sequences and controls
arid1a.seq<-arid1a.exons.seqs[ arid1a.fam.summarized$circ_name ]
arid1a.controls.seq<-arid1a.controls.seqs.resized[ intersect(arid1a.fam.summarized$circ_name, names(arid1a.controls.seqs.resized))  ]


#  add circRNA and control lengths 
stopifnot( all.equal( names(arid1a.seq), arid1a.fam.summarized$circ_name ) )
arid1a.fam.summarized[, circ_width:=width(arid1a.seq) ]
arid1a.fam.summarized[ match(names(arid1a.controls.seq), circ_name), controls_width:=width(arid1a.controls.seq)]


#  count non-overlapping seeds across all families for both miRNA and mutated seeds
#  annotate them with AGO coverage based on centered 41nts windows around the seed
arid1a.fam.summarized[, c('ncount', 'seeds_hits', 'seeds_start', 'seeds_end', 'ago'):=list(0, list(), list(), list(), list())]
arid1a.fam.summarized[, c('ncount.mut', 'seeds_hits.mut', 'seeds_start.mut', 'seeds_end.mut', 'ago.mut'):=list(0, list(), list(), list(), list())]
#arid1a.fam.summarized[, c('ncount', 'seeds_hits', 'seeds_start', 'seeds_end', 'ago'):=list(NULL, NULL, NULL, NULL, NULL)]
#arid1a.fam.summarized[, c('ncount.mut', 'seeds_hits.mut', 'seeds_start.mut', 'seeds_end.mut', 'ago.mut'):=list(NULL, NULL, NULL, NULL, NULL)]
for(n in 1:nrow(arid1a.fam.summarized)){
    pd<-PDict(DNAStringSet(arid1a.fam.summarized[n, seed][[1]]))
    m<-unlist(as(matchPDict(pd, arid1a.seq[[ arid1a.fam.summarized[n, circ_name] ]], min.mismatch=0, with.indels=F), 'CompressedIRangesList'))
    ov<-findOverlaps(m, m, type='any', select='all')
    if(!all(queryHits(ov)==subjectHits(ov))){ 
        ov<-ov[ queryHits(ov)<subjectHits(ov) ]
        m<-m[-subjectHits(ov)]
    }
    if(length(m)>0){
        coordinates<-circ2genome(CIR=m, EXN=arid1a.exons[[ arid1a.fam.summarized[n, circ_name] ]], SEQ=extractAt(arid1a.seq[[n]], m))
        a<-sapply(coordinates$seqs, toString)
        b<-sapply(getSeq(BSgenome.Hsapiens.UCSC.hg38, coordinates), toString)
        i<-which(a!=b)
        a<-a[i]
        b<-b[i]
        if(length(a)>0){
            cat('signal seeds:', a, 'of circRNA:', arid1a.fam.summarized[n, circ_name], 'are spliced and their AGO coverage is omitted!\n')
        }
        ov<-data.table(data.frame(findOverlaps(resize(coordinates, width=41, fix='center'), ago, type='any', select='all', ignore.strand=F)))
        colnames(ov)<-c('qh', 'sh')
        ov[, nsum:=ago$score[ sh ]]
        ov<-ov[, .(nsum=sum(nsum)), by=.(qh)]
        coordinates$ago<-0
        coordinates$ago[ ov$qh ]<-ov$nsum
        coordinates$ago[ i ]<-0  #  omit AGO coverage of spliced-seeds
        arid1a.fam.summarized[n, c('ncount', 'seeds_hits', 'seeds_start', 'seeds_end', 'ago'):=list(length(m), list(sapply(coordinates$seqs, toString)), list(start(coordinates)), list(end(coordinates)), list(coordinates$ago))]
    } 


    #  do the same for the mutated seeds 
    pd<-PDict(DNAStringSet(arid1a.fam.summarized[n, seed.mut][[1]]))
    m<-unlist(as(matchPDict(pd, arid1a.seq[[ arid1a.fam.summarized[n, circ_name] ]], min.mismatch=0, with.indels=F), 'CompressedIRangesList'))
    ov<-findOverlaps(m, m, type='any', select='all')
    if(!all(queryHits(ov)==subjectHits(ov))){ 
        ov<-ov[ queryHits(ov)<subjectHits(ov) ]
        m<-m[-subjectHits(ov)]
    }
    if(length(m)>0){
        coordinates<-circ2genome(CIR=m, EXN=arid1a.exons[[ arid1a.fam.summarized[n, circ_name] ]], SEQ=extractAt(arid1a.seq[[n]], m))
        a<-sapply(coordinates$seqs, toString)
        b<-sapply(getSeq(BSgenome.Hsapiens.UCSC.hg38, coordinates), toString)
        i<-which(a!=b)
        a<-a[i]
        b<-b[i]
        if(length(a)>0){
            cat('control seeds:', a, 'of circRNA:', arid1a.fam.summarized[n, circ_name], 'are spliced and their AGO coverage is omitted!\n')
        }
        ov<-data.table(data.frame(findOverlaps(resize(coordinates, width=41, fix='center'), ago, type='any', select='all', ignore.strand=F)))
        colnames(ov)<-c('qh', 'sh')
        ov[, nsum:=ago$score[ sh ]]
        ov<-ov[, .(nsum=sum(nsum)), by=.(qh)]
        coordinates$ago<-0
        coordinates$ago[ ov$qh ]<-ov$nsum
        coordinates$ago[ i ]<-0  #  omit AGO coverage of spliced-seeds
        arid1a.fam.summarized[n, c('ncount.mut', 'seeds_hits.mut', 'seeds_start.mut', 'seeds_end.mut', 'ago.mut'):=list(length(m), list(sapply(coordinates$seqs, toString)), list(start(coordinates)), list(end(coordinates)), list(coordinates$ago))]
    } 
}


#  save
save(arid1a, arid1a.exons.controls, arid1a.controls.fam, arid1a.controls.seq, arid1a.controls.seqs, arid1a.controls.seqs.resized, arid1a.exons, 
     arid1a.exons.seqs, arid1a.fam, arid1a.fam.summarized, arid1a.seq, 
     file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_hsa.mature.miRBase21_7mer_seed_counts.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/ARID1A_hsa.mature.miRBase21_7mer_seed_counts.RData




