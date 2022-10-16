###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




############################
#
#
#  aggregate sample metadata
#
#
############################




#  process cell lines and tumors metadata
#  group replicates under a common Berlin ID
#  separate pre-failed samples based on QC, or after mapping based on percentage of alignments + unmapped reads dropped from feature counting
#  define color palettes 
#
#  Some comments from Falk:
#
#       All tumors are 'initial' pre-biopsy tumors. The age at diagnosis is defined by the date of confirming the diagnosis, 
#       most of the time by histo-pahology or cytological evaluation (bone marrow) - in the first place, this means the tumor biopsy sample was 
#       taken days before the day of confirmed diagnosis. 
#       In other cases histology is somewhat behind and if the phenotype was confirmed by other imaging methods (e.g. ultrasound + mIBG scan and 
#       maybe bonemarrow infiltration) then maybe the hiostology on a tumor biopsy was taken later, to confirm - then the day of diagnosis may 
#       precede the tumor biopsy sampling. 
#
#{{{
rm(list=ls())
library(data.table)
library(RColorBrewer)
library(openxlsx)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  cell models
meta.cel<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/cell_models.csv', header=T, sep=',', select=1:6, col.names=c('bids', 'ids', 'library_prep', 'treatment', 'replicate', 'cell_model'))


#  tumors (remove mixups column, they have been converted to replicates)
meta.tum<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/tumor_samples.csv', sep=',', header=T, select=c(1:14), col.names=c('bids', 'ids', 'library_prep', 'pid', 'risk_group', 'group_v', 'stage', 'age', 'MYCN_FC', 'status', 'remission', 'survival', 'reason_of_death', 'LOH'))


#  change nMNA_HR to HR_nMNA 
#  change MNA-HET to MNA
#  patients without risk-group metadata should have risk-group NA for easy removal
#  remove question marks
meta.tum[, risk_group:=sub('nMNA_HR', 'HR_nMNA', risk_group)]
meta.tum[, risk_group:=sub('MNA-HET', 'MNA', risk_group)]
meta.tum[, risk_group:=sub('\\?', '', risk_group)]
meta.tum[ risk_group %in% '', risk_group:=NA]


#  [should be run at the BIH cluster] manually check symbolic links 
#{{{

#  identify the directories that have links
l<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/CB* -maxdepth 1 -type l -name "r1.fastq.gz" -printf "%h " -exec readlink -f {} \\;', stdout=T)
l<-data.frame(do.call(rbind, lapply(strsplit(l, ' ', fixed=T), function(x){ t(data.frame(x)) })))
rownames(l)<-NULL
colnames(l)<-c('dir', 'link')
l$dir<-basename(l$dir)
l$link<-basename(l$link)
l<-data.table(l)


#  identify the directories from the Berlin ID column as well
m<-rbind(meta.tum, meta.cel, fill=T)
m[, dir:=sub('-R01-R[0-9]$', '-R01', bids)]
m<-m[, .(bids=list(unique(bids)), ids=list(unique(ids))), by=.(dir)]
links<-m[ dir %in% l$dir ]
links<-l[ links, on='dir']
rest<-m[ ! dir %in% l$dir ]
stopifnot( nrow(links) + nrow(rest) == nrow(m) )
rm(m)


#  visually inspect one by one the symbolic links
for (n in 1:nrow(links)){
    cat(links[n, dir], ':', links[n, unlist(ids)], ':', links[n, link], '\n')
    readline()
}


#  visually inspect one by one the rest of the libraires (libraries not sequenced yet, or failed, or pooled)
for (n in 1:nrow(rest)){
    cat(rest[n, dir], ':', rest[n, unlist(ids)], '\n')
    readline()
}
rm(n, rest, links)

#}}}


#  group replicates under a single Berlin ID (failed but otherwise unrelated samples will group together)
#  keep distinctive information but collapse all identical information across replicates 
#  separate the metadata of failed samples and keep all information 
meta.tum[, bid:=sub('-R01-R[0-9]*$', '-R01', bids)]
meta.cel[, bid:=sub('-R01-R[0-9]*$', '-R01', bids)]
meta.prefailed<-rbind(meta.tum[ bid %in% 'failed' ], meta.cel[ bid %in% 'failed' ], fill=T)
meta.tum<-meta.tum[ !bid %in% 'failed' ]
meta.cel<-meta.cel[ !bid %in% 'failed' ]
meta.prefailed<-meta.prefailed[, lapply(.SD, list), by=.(bid)]
setcolorder(meta.prefailed, c('bid', 'bids', 'ids', 'library_prep', 'cell_model', 'treatment', 'risk_group', 'pid', 'group_v', 'stage', 'age', 'MYCN_FC', 'status', 'remission', 'survival', 'reason_of_death', 'LOH'))
meta.tum<-meta.tum[, .(bids=list(bids), ids=list(ids), library_prep=list(library_prep), risk_group=unique(risk_group), pid=unique(pid), group_v=unique(group_v), stage=unique(stage), age=unique(age), MYCN_FC=unique(MYCN_FC), status=unique(status), remission=unique(remission), survival=unique(survival), reason_of_death=unique(reason_of_death), LOH=unique(LOH)), by=.(bid)]
meta.cel<-meta.cel[, .(bids=list(bids), ids=list(ids), library_prep=list(library_prep), cell_model=unique(cell_model), treatment=unique(treatment)), by=.(bid)]


#  load FASTQC metrics from MultiQC report
#  identify bid
#  summarize %GC and total number of raw reads across read mates
#  add to metadata
fgc<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/multiqc/multiqc_data/multiqc_fastqc.txt')[, c('Sample', 'Filename', '%GC', 'Total Sequences', 'total_deduplicated_percentage')]
colnames(fgc)<-c('sample', 'filename', 'gc', 'nreads', 'dedup')
fgc[, bid:=sub('^.*(CB[^ ]+).*$', '\\1', sample) ]
fgc<-fgc[, .(gc=mean(gc), nreads=mean(nreads), dedup=mean(dedup)), by=.(bid)]
meta.tum<-fgc[meta.tum, on='bid']
meta.cel<-fgc[meta.cel, on='bid']
rm(fgc)


#  process featureCounts statistics and add to metadata
logs<-system2('find', args='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/ -maxdepth 3 -wholename \'*/counts/genes.tsv.summary\' -print', stdout=T)
names(logs)<-sub('^.*raw/(CB[^/]+)/counts/.*$', '\\1', logs)
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
meta.tum<-fc[meta.tum, on='bid']
meta.cel<-fc[meta.cel, on='bid']
rm(fc, l, logs, bid)


#  add percentage of the total number of alignments+reads that were dropped 
#  mark failed samples with at least 50% of dropouts in alignments+reads and number of reads covering features below the median
meta.tum[, p_unmapped:=100*unmapped/(features+no_features+unmapped)]
meta.cel[, p_unmapped:=100*unmapped/(features+no_features+unmapped)]
meta.tum[, failed:=ifelse(p_unmapped>=50 & features<median(features, na.rm=T), T, F)]
meta.cel[, failed:=ifelse(p_unmapped>=50 & features<median(features, na.rm=T), T, F)]
setcolorder(meta.tum, c('bid', 'failed', 'gc', 'nreads', 'dedup', 'features', 'no_features', 'unmapped', 'p_unmapped', 'bids', 'ids', 'library_prep', 'risk_group', 'pid', 'group_v', 'stage', 'age', 'MYCN_FC', 'status', 'remission', 'survival', 'reason_of_death', 'LOH'))
setcolorder(meta.cel, c('bid', 'failed', 'gc', 'nreads', 'dedup', 'features', 'no_features', 'unmapped', 'p_unmapped', 'bids', 'ids', 'library_prep', 'cell_model', 'treatment'))


#  manually designate CB2009-11-R01-R1 as failed
meta.tum[ bid %in% 'CB2009-11-R01', failed:=T ]


#  manually remove CB-p53-IMR5-DS3032b_1, CB-EWS-A4573-FCS-R01
meta.cel[ is.na(failed) ]
meta.cel<-meta.cel[ !is.na(failed) ]


#  inspect failed tumors and cell lines
meta.tum[ (failed) ]
meta.cel[ (failed) ]


#  integrate most recent clinical metadata
#
#  Falk writes:
#
#      CB2028 : het-MNA case which was only MNA in a relapse sample, while the initial diagnostic is and remains non-MNA.
#
#      CB3016 : het-MNA case with vast majority of MNA tumor cells and only minor fraction of tumor cells w/o MYCN-amp, keep it as MNA as we probably 
#               underestimate the detection of such cases among other MNA tumors.
#
#      CB3036 : het-MNA case considered as MNA when it comes to clinical/diagnostic risk evaluation but non-MNA when it comes to biological phenotype.
#               WE WILL REMOVE THIS SAMPLE FROM THE COHORT BY MARKING IT AS 'failed'.
#
#      CB2006 : correct age of diagnosis is 75
#
#      status : 0 = alive without progression, 1 = died from disease, 2 = alive after suffering a progression event
#
#{{{
m<-data.table(read.xlsx('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/CLINICAL_DATA.xlsx'))
m[ SUBGROUP %in% 'HR_MNA', SUBGROUP:='MNA']
m[ SUBGROUP %in% 'HR_non-MNA', SUBGROUP:='HR_nMNA']
m[ SUBGROUP %in% 'LR_4S', SUBGROUP:='ST4S']
m[, table(SUBGROUP)]
# SUBGROUP
# HR_nMNA     IMR      LR     MNA    ST4S 
#      29      10      30      23      13 


#  any clinical metadata missing?
setdiff(sub('-.*$', '', meta.tum[, bid]), m$PAT_ID_BERLIN)  #  CB3009


#  add PAT_ID_BERLIN key to tumors for easy integration
meta.tum[, PAT_ID_BERLIN:=sub('-.*$', '', bid) ]
meta.tum[, table(risk_group)]
# risk_group
# HR_nMNA     IMR      LR     MNA    ST4S 
#      30      10      29      24      13 


#  merge metadata by PAT_ID_BERLIN
meta.tum<-m[meta.tum, on='PAT_ID_BERLIN']  #  preserves CB3009


#  manually add missing clinical metadata for CB3009 from the previous version:
#
#      TRIAL_PROTOCOL:NB43, INSS_STAGE:=4, EFS_bin:1, EFS_days:688, OS_bin:0, OS_days:2661, STATUS:2, MYCN_ini_status:non-amp, subgroup:HR_nMNA
#
meta.tum[ PAT_ID_BERLIN %in% 'CB3009', c('TRIAL_PROTOCOL', 'INSS_STAGE', 'AGE_days', 'EFS_days', 'OS_days', 'EFS_bin', 'OS_bin', 'STATUS', 'MYCN_INI_STATUS', 'SUBGROUP', 'INRG_HR_risk_group', 'MYCN_FC'):=list('NB43', 4, 1419, 688, 2661, 1, 0, 2, 'non-amp', 'HR_nMNA', 'HR', 1.5)]


#  risk_group discrepancies!
meta.tum[ SUBGROUP!=risk_group ,c('PAT_ID_BERLIN', 'bid', 'INRG_HR_risk_group', 'SUBGROUP', 'MYCN_FC', 'failed', 'risk_group')]
#
#    PAT_ID_BERLIN           bid INRG_HR_risk_group SUBGROUP MYCN_FC failed risk_group
# 1:        CB3041 CB3041-11-R01             non-HR       LR       1  FALSE    HR_nMNA
# 2:        CB2047 CB2047-11-R01                 HR  HR_nMNA       1  FALSE        MNA


#  are there any discrepancies in the age of diagnosis?
meta.tum[ AGE_days!=age ,c('PAT_ID_BERLIN', 'bid', 'INRG_HR_risk_group', 'SUBGROUP', 'MYCN_FC', 'failed', 'risk_group', 'AGE_days', 'age')]
meta.tum[ AGE_days!=age , age-AGE_days ]
#
#  -6   10    1   -2  -17 -377    6    1    3    2    1    3   13


#  are there any staging discrepancies?
meta.tum[ INSS_STAGE!=stage ,c('PAT_ID_BERLIN', 'bid', 'INRG_HR_risk_group', 'SUBGROUP', 'MYCN_FC', 'failed', 'risk_group', 'INSS_STAGE', 'stage')]
#
#  => not really, there are ambiguous or missing entries only


#  status discrepancies!!
meta.tum[ STATUS!=status ,c('PAT_ID_BERLIN', 'bid', 'INRG_HR_risk_group', 'SUBGROUP', 'MYCN_FC', 'failed', 'risk_group', 'INSS_STAGE', 'stage', 'STATUS', 'status')]
#    PAT_ID_BERLIN           bid INRG_HR_risk_group SUBGROUP MYCN_FC failed risk_group INSS_STAGE stage STATUS status
# 1:        CB3001 CB3001-11-R01             non-HR     ST4S       1  FALSE       ST4S          5     5      0      2
# 2:        CB3004 CB3004-11-R01                 HR  HR_nMNA     1.5  FALSE    HR_nMNA          4     4      1      2
# 3:        CB3018 CB3018-11-R01                 HR      MNA      30  FALSE        MNA          4     4      2      1
# 4:        CB3022 CB3022-11-R01                 HR  HR_nMNA       1  FALSE    HR_nMNA          4     4      2      1
# 5:        CB3023 CB3023-11-R01                 HR  HR_nMNA       1  FALSE    HR_nMNA          4     4      2      1
# 6:        CB3028 CB3028-11-R01                 HR  HR_nMNA     1.5  FALSE    HR_nMNA          4     4      1      2
# 7:        CB3034 CB3034-11-R01                 HR  HR_nMNA       1  FALSE    HR_nMNA          4     4      1      2


#  with heterogeneous MNA state and classified as MNA (CB2028 is not part of them)
meta.tum[ MYCN_FC %in% 'het' | PAT_ID_BERLIN %in% 'CB2028', c('PAT_ID_BERLIN', 'bid', 'MYCN_INI_STATUS', 'MYCN_FC', 'failed', 'risk_group')]
#
#    PAT_ID_BERLIN           bid MYCN_INI_STATUS MYCN_FC failed risk_group
# 1:        CB2028 CB2028-11-R01         non-amp       1  FALSE    HR_nMNA
# 2:        CB3016 CB3016-11-R01             amp     het  FALSE        MNA
# 3:        CB3036 CB3036-11-R01             amp     het  FALSE        MNA


#  assume latest clinical metadata are the most up-to-date
#  fix CB2006 age back to 75
meta.tum[, risk_group:=SUBGROUP]
meta.tum[, status:=STATUS][, c('SUBGROUP', 'STATUS', 'CAUSE_OF_DEATH', 'COMMENT_THERA', 'COMMENT_MNA', 'pid', 'group_v', 'stage', 'age', 'i.MYCN_FC', 'remission', 'survival', 'reason_of_death', 'LOH'):=NULL]
setcolorder(meta.tum, c('PAT_ID_BERLIN', 'bid', 'bids', 'ids', 'library_prep', 'COHORT', 'TRIAL_PROTOCOL', 'failed', 'gc', 'nreads', 'dedup', 'features', 'no_features', 'unmapped', 'p_unmapped', 'SEX', 'status', 'INSS_STAGE', 'AGE_days', 'EFS_days', 'OS_days', 'EFS_bin', 'OS_bin', 'MYCN_INI_STATUS', 'INRG_HR_risk_group', 'MYCN_FC', 'risk_group'))
meta.tum[ PAT_ID_BERLIN %in% 'CB2006', AGE_days:=75]


#  mark CB3036 as 'failed' so that it is removed from the cohort
meta.tum[ PAT_ID_BERLIN %in% 'CB3036', failed:=T ]

#}}}


#  add set colors for tumors
meta.tum[, col:='black']
meta.tum[ risk_group %in% 'ST4S', col:='seagreen4' ] 
meta.tum[ risk_group %in% 'LR', col:='darkgreen' ] 
meta.tum[ risk_group %in% 'IMR', col:='cornflowerblue' ] 
meta.tum[ risk_group %in% 'HR_nMNA', col:='chocolate1' ] 
meta.tum[ risk_group %in% 'MNA', col:='coral4' ] 


#  add color palette for cell models
N<-meta.cel[ !is.na(cell_model), unique(paste(cell_model, treatment, sep='_'))]
cl<-setNames(colorRampPalette(brewer.pal(8, 'Dark2'))(length(N)), N)
for( l in seq_along(cl)){
    meta.cel[ paste(cell_model, treatment, sep='_') %in% names(cl)[l], col:=cl[l] ] 
}
rm(N,cl)


#  group tumors by risk group and cell lines by cell model and within each subgroup by bid
meta.tum<-meta.tum[, risk_group:=factor(risk_group, levels=c('ST4S', 'LR', 'IMR', 'HR_nMNA', 'MNA'))][ order(risk_group, bid), ][, risk_group:=as.character(risk_group)]
meta.cel<-meta.cel[, cell_model:=factor(cell_model, levels=as.character(na.omit(unique(cell_model))))][ order(cell_model, bid), ][, cell_model:=as.character(cell_model)]


#  save
save(meta.tum, meta.cel, meta.prefailed, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData
#
#  typical run of the pipeline modules and the final MultiQC summary:
#
#      make -f lib/qsub.mf CONF=raw/totalrna.conf qc optitype kallisto rrna star bwa
#      make -f lib/qsub.mf CONF=raw/totalrna.conf counts ciri lin_vs_circ
#
#      multiqc --interactive -o raw/multiqc --ignore '*dcc*' --ignore '*ciri*' --ignore '*hg19*' --ignore '*lin_vs_circ*' --ignore '*FAILED*' -v -f -d -s raw/




#########################
#
#
#  biomaRt GO annotations 
#
#
#########################




#  [run once] process the raw downloaded biomaRt GO annotations that exclude the
#             Inferred from Electronic Annotation (IEA) evidence code
#{{{
rm(list=ls())
library(data.table)


#  load raw biomaRt data
#  remove genes with no annotations
GO<-fread('/data/annotation/GRCh38/GRCh38.p12.biomaRt.GO.tsv', header=T, sep='\t', quote='', col.names=c('gene_id', 'gene_name', 'go_id', 'go_name', 'go_desc', 'go_ecode', 'go_domain'), fill=T)
GO<-GO[ go_domain!='' ]


#  split BP, MF, CC 
#  we remove gene_name since grouping by go_id might result in length(unique(gene_name)) != length(unique(gene_id)) 
#  we remove also the redundant go_domain column
BP<-GO[ go_domain %in% 'biological_process' ]
MF<-GO[ go_domain %in% 'molecular_function' ]
CC<-GO[ go_domain %in% 'cellular_component' ]
stopifnot( nrow(BP)+nrow(MF)+nrow(CC)==nrow(GO) )
BP<-BP[, c('gene_id', 'go_id', 'go_name', 'go_desc', 'go_ecode'), with=F]
MF<-MF[, c('gene_id', 'go_id', 'go_name', 'go_desc', 'go_ecode'), with=F]
CC<-CC[, c('gene_id', 'go_id', 'go_name', 'go_desc', 'go_ecode'), with=F]
rm(GO)


#  GO2GENE_ID grouping:
#
#      go_id:gene_id1, gene_id2, ...
#
BP.go2gid<-BP[, .(go_name=unique(go_name), go_ecode=list(unique(go_ecode)), gene_id=list(unique(gene_id))), by=.(go_id)]
MF.go2gid<-MF[, .(go_name=unique(go_name), go_ecode=list(unique(go_ecode)), gene_id=list(unique(gene_id))), by=.(go_id)]
CC.go2gid<-CC[, .(go_name=unique(go_name), go_ecode=list(unique(go_ecode)), gene_id=list(unique(gene_id))), by=.(go_id)]


#  GENE_ID2GO grouping:
#
#      gene_id:go_id1, go_id2, ...
#
BP.gid2go<-BP[, .(go_ecode=list(go_ecode), go_id=list(go_id)), by=.(gene_id)]
MF.gid2go<-MF[, .(go_ecode=list(go_ecode), go_id=list(go_id)), by=.(gene_id)]
CC.gid2go<-CC[, .(go_ecode=list(go_ecode), go_id=list(go_id)), by=.(gene_id)]


#  save them all
save(BP, MF, CC, BP.go2gid, MF.go2gid, CC.go2gid, BP.gid2go, MF.gid2go, CC.gid2go, file='/data/annotation/GRCh38/GRCh38.p12.biomaRt.GO.RData')

#}}}
#
#  => /data/annotation/GRCh38/GRCh38.p12.biomaRt.GO.RData




#######################
#
#
#  MYCN list of targets
#
#
#######################




#  [run once] process the MYCN list of targets from 
#             Valentjin et al. PNAS 2012 (https://www.pnas.org/content/109/47/19190)
#             and Westermann et al Genome Biology 2008 (https://doi.org/10.1186/gb-2008-9-10-r150)
#{{{
rm(list=ls())
library(data.table)
library(GenomicFeatures)
library(rtracklayer)


#  load gene annotation
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
mcols(hsa)<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]


#  Valentjin et al PNAS 2012 gene set
#  log2-FC values are associated with MYCN knockdown
#{{{

g<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/mycn_targets_pnas.1208215109.csv', header=T, sep=',')[, c(1, 3:4)]
colnames(g)<-c('gene_name', 'log2fc1', 'log2fc2')
stopifnot( length(g[ log2fc1*log2fc2<0, gene_name])==0 )  #  expect both measurements to agree how the gene level moves upon MYCN knockdown
g[, type:=ifelse( log2fc1<0, 'induced', 'repressed' )]    #  use first measurement to determine induction of repression activity of MYCN on the gene
g<-g[, c('log2fc1', 'log2fc2'):=list(NULL, NULL)]


#  remove star from gene names indicating the nerve/brain-specific genes found in another study
g[, gene_name:=sub('\\*', '', gene_name)]


#  manually fix old gene_names quering the HUGO database:
#
#      https://www.genenames.org/cgi-bin/symbol_checker
#
print(setdiff( g[, gene_name], hsa$gene_name ))
cat(setdiff( g[, gene_name], hsa$gene_name ), file='~/Downloads/genes')  #  upload to HUGO
h<-read.table('~/Downloads/hgnc-symbol-check.csv', header=T, sep=',', skip=1, col.names=c('input', 'type', 'gene_name', 'description', 'id', 'location'))
h<-h[ !h$type %in% 'Unmatched', ]
print(setdiff( h$gene_name, hsa$gene_name ))  #  FAM85A does not exist it is replaced with the opposite strand FAM85B
h[ h$input %in% 'FAM85A', 'gene_name']<-'FAM85B'
h<-h[ !duplicated(h$input), ]                 #  RAGE is somehow associated with MOK and AGER, both present in GENCODE v30
g[ match(h$input, gene_name), gene_name:=h$gene_name]


#  remove genes not found in HUGO or GENCODE v30
g<-g[ gene_name %in% intersect(gene_name, hsa$gene_name), ]

#}}}


#  initialize the target list
MYCN<-g


#  Westermann et al Genome Biology 2008
#{{{

g1<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/mycn_gb-2008-9-10-r150_class1,2.csv', header=T, sep=',')
g2<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/mycn_gb-2008-9-10-r150_class3.csv', header=T, sep=',')
g3<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/lib/mycn_gb-2008-9-10-r150_class4.csv', header=T, sep=',')
g<-unique(data.table(gene_name=c(g1$'Gene Name', g2$'Gene Name', g3$'Gene Name'), type='induced'))
rm(g1,g2,g3)


#  manually fix old gene_names quering the HUGO database:
#
#      https://www.genenames.org/cgi-bin/symbol_checker
#
print(setdiff( g[, gene_name], hsa$gene_name ))
cat(setdiff( g[, gene_name], hsa$gene_name ), file='~/Downloads/genes')  #  upload to HUGO
h<-read.table('~/Downloads/hgnc-symbol-check.csv', header=T, sep=',', skip=1, col.names=c('input', 'type', 'gene_name', 'description', 'id', 'location'))
h<-h[ !h$type %in% 'Unmatched', ]
print(setdiff( h$gene_name, hsa$gene_name ))  #  GAGE12I, TOMM40
h<-h[ !duplicated(h$input), ]                 #  NP is assiciated with ZNF384, PNP, CTF2P
g[ match(h$input, gene_name), gene_name:=h$gene_name]


#  remove genes not found in HUGO or GENCODE v30
g<-g[ gene_name %in% intersect(gene_name, hsa$gene_name), ]

#}}}


#  visually inspect common terms among the lists before removing
common<-intersect(MYCN$gene_name, g$gene_name)
MYCN[ gene_name %in% common ]
g[ gene_name %in% common ]
g<-g[ ! gene_name %in% common, ]


#  append to list
MYCN<-rbind(MYCN, g)
rm(g)


#  sort by type of regulation
MYCN<-MYCN[ order(type), ]


#  add gene_id
MYCN<-MYCN[, gene_id:=hsa$gene_id[ match( gene_name, hsa$gene_name ) ]]


#  reorder columns
MYCN<-MYCN[, c('gene_id', 'gene_name', 'type')]


#  save
save(MYCN, file='/fast/groups/ag_schulte/work/reference/annotation/GRCh38/MYCN_targets.RData')

#}}}
#
#  => /fast/groups/ag_schulte/work/reference/annotation/GRCh38/MYCN_targets.RData




###################################################
#
#
#  human RBPs: published data preparation/lift over
#
#
###################################################




#  human RBPs found upregulated in early brain development and downregulated in mature brain tissue (Gerstberger et al. 2014)
#  
#      https://www.nature.com/articles/nrg3813#f7
#
#  lift over to GENCODE v30
#  identify group I by repeating the procedure described in the paper
#{{{
rm(list=ls())
library(rtracklayer)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  do the lift over
#{{{

#  load annotation
hsa<-import('/fast/work/groups/ag_schulte/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
hsa<-data.table(data.frame(mcols(hsa)[, c('gene_name', 'gene_id', 'gene_type')]))


#  load the supplementary table
rbp<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/validations/rbp/brainspan_census_of_human_RBPs_nrg3813/nrg3813-s7.csv', sep=',', header=T)
colnames(rbp)<-c('target', 'old', 'gene_id', 'score', 'mean_RPKM', 'min', 'max', 'pcw_8', 'pcw_9', 'pcw_12', 'pcw_13', 'pcw_16', 'pcw_17', 'pcw_19', 'pcw_21', 'pcw_24', 'pcw_37', 'pcw_4month', 'pcw_1yr')
rbp[, gene_name:=old]


#  identify problematic genes, manually paste them to the HGNC database:
#
#      https://www.genenames.org/cgi-bin/symbol_checker
#
#  download the conversion as CSV
#
cat(setdiff( rbp$gene_name, hsa$gene_name ), file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/validations/rbp/brainspan_census_of_human_RBPs_nrg3813/genes.tsv', sep='\n')
g<-fread('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/validations/rbp/brainspan_census_of_human_RBPs_nrg3813/hgnc_genes.tsv', header=T, sep=',')
colnames(g)<-c('old', 'type', 'gene_name', 'description', 'id', 'location')
rbp[ match(g$old, old), gene_name:=g$gene_name ]
rm(g)


#  report those that were not found in HGNC so we can manually look them up
rbp[ gene_name %in% '', c('old', 'gene_name')]
rbp[ old %in% 'AC004381.6', gene_name:='REXO5' ]  #  google results point to ENST00000261377 transcript (REXO5-201)
rbp[ old %in% 'SF3B14', gene_name:='SF3B6' ]
rbp[ old %in% 'EIF2S3L', gene_name:='EIF2S3B' ]   #  EIF2S3L maps to EIF2S3B (ENSG00000180574) 
stopifnot( all(!is.na(rbp$gene_name)) )


#  update gene_id and remove old names
rbp[, gene_id:=hsa$gene_id[ match( gene_name , hsa$gene_name ) ]]
rbp[, c('old', 'min', 'max'):=list(NULL, NULL, NULL)]
setcolorder(rbp, c('target', 'gene_name', 'gene_id', 'score', 'mean_RPKM', 'pcw_8', 'pcw_9', 'pcw_12', 'pcw_13', 'pcw_16', 'pcw_17', 'pcw_19', 'pcw_21', 'pcw_24', 'pcw_37', 'pcw_4month', 'pcw_1yr'))
rm(hsa)

#}}}


#  follow the selection process as described in the supplement to identify groupI:
#
#      keep RBPs with sum of RPKMs across all Post Conception Week (PCW) timepoints >= 12
#      order by tissue specificity score and select top 200
#
#      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#      !!!!                                                                                  !!!!
#      !!!!  N.B. mean_RPKM column does not correspond to the mean across the time points    !!!!
#      !!!!                                                                                  !!!!
#      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#      compute mean across time points
#      divide RPKMs across time points by mean and take the log2(0.1+...) to compute log2FC 
#      cluster RBPs based on Pearson's correlations and the Ward.D2 method
#{{{

#  copy original data so that we can keep all RBPs
RBP<-copy(rbp)


#  select RBPs 
rbp<-rbp[ rowSums(as.data.frame(rbp[, c(6:ncol(rbp)), with=F]))>=12, ]  #  1402 RBPs should make it
rbp<-rbp[ order(-score) ][1:200, ]
rbp[, mean:=copy(rbp[, 6:17, with=F])[, n:=1:nrow(rbp)][, .(mean=rowMeans(.SD)), by=.(n)][, n:=NULL][, mean] ]
x<-log2( 0.1 + sweep(as.data.frame(rbp[, 6:17, with=F]), 1, rbp$mean, '/') )
rownames(x)<-rbp$gene_name
rm(rbp)


#  plot
x11(width=14, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
d<-cor(t(x), method='pearson', use='pairwise.complete.obs')
hc<-hclust(dist(d, method='euclidean'), method='ward.D2')
p<-pheatmap(as.matrix(d), color=colorRampPalette(c('#FFFFFF', '#FFFF33', '#FF0000'))(20), border_color=NA, scale='none', 
    breaks=seq(-1, 1, length.out=21),
    cluster_rows=hc,
    cluster_cols=hc,
    annotation_legend=F, annotation_names_row=T, annotation_names_col=T, 
    drop_levels=F, show_rownames=T, show_colnames=T, 
    fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
#dev.print(device=pdf, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/validations/rbp/brainspan_census_of_human_RBPs_nrg3813/heatmap_cor.pdf', width=16, height=16, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
dev.print(device=svg, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/validations/rbp/brainspan_census_of_human_RBPs_nrg3813/heatmap_cor.svg', width=16, height=16, bg='white', antialias='subpixel', pointsize=20, family='Arial')
#
#  looks good, very faithful to the published Figure 7B


#  identify groupI and check
ww<-cutree(hc, k=2)[ hc$order ]
ww<-names(ww[ ww==2 ])
xx<-x[ match(ww, rownames(x)), ]
dd<-cor(t(xx), method='pearson', use='pairwise.complete.obs')
p<-pheatmap(as.matrix(dd), color=colorRampPalette(c('#FFFFFF', '#FFFF33', '#FF0000'))(20), border_color=NA, scale='none', 
    breaks=seq(-1, 1, length.out=21),
    cluster_rows=F,
    cluster_cols=F,
    annotation_legend=F, annotation_names_row=T, annotation_names_col=T, 
    drop_levels=F, show_rownames=T, show_colnames=T, 
    fontsize = 12, fontsize_row=12, fontsize_col=12, fontsize_number=5.0)
rm(d, hc, dd, p, x, xx)


#  add groupI annotation to original data
RBP$groupI<-F
RBP[ gene_name %in% ww, groupI:=T ]
stopifnot( RBP[, sum(groupI)]==length(ww) )
rbp<-copy(RBP)
rm(ww,RBP)

#}}}


#  save
save(rbp, file='/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/validations/rbp/brainspan_census_of_human_RBPs_nrg3813/nrg3813-s7.RData')

#}}}
#
#  => /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/validations/rbp/brainspan_census_of_human_RBPs_nrg3813/nrg3813-s7.RData




########################################################################################################
#
#
#  ATtRACT analysis (preparation of the ATtRACT database is needed for the MEME motif discovery as well)
#
#
########################################################################################################




#  [run once] process
#             lift over to GENCODE v30
#             create MEME motif file based on the circRNA background nucleotide frequencies for all motifs
#             create MEME motif file based on the circRNA background nucleotide frequencies for 6mers or longer motifs
#{{{
rm(list=ls())
library(rtracklayer)
library(GenomicFeatures)
library(Biostrings)
library(data.table)


#  process the DB and the PWM into R-objects
#     
#      qscore : scores the binding affinity of given unambiguous motif as the product of probabilities for each nucleotide separately
#     
#  RPBs of identical motifs but of different PWMs coming from different experiments are kept
#  FUCKED up entries where motif length does not match PWM length are thrown out
#  dublicate PWMs with different consensus sequence motifs are thrown out and a unique PWM is kept
#
#{{{

#  load annotation
#  add column of gene_id without subversion 
hsa<-import('/fast/work/groups/ag_schulte/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
hsa<-mcols(hsa)[, c('gene_id', 'gene_type', 'gene_name')]
hsa$gid<-sub('\\.[0-9]*$', '', hsa$gene_id)


#  import the database
db<-fread('/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.tsv', sep='\t', header=T)[, c(1:2, 5, 12:13)]
colnames(db)<-c('gene_name', 'gene_id', 'motif', 'id', 'qscore')
db[, qscore:=as.numeric(sub('\\**$', '', qscore))]


#  convert U to T
db[, motif:=gsub('U', 'T', motif)]


#  add subversion to gene_id and update gene_names just in case they are old
stopifnot( length(setdiff( db$gene_id, hsa$gid ))==0 )
db[, gene_id:=hsa$gene_id[ match(gene_id, hsa$gid) ]]
db[, gene_name:=hsa$gene_name[ match(gene_id, hsa$gene_id) ]]
rm(hsa)


#  PWM is in FASTA format with PWM named by the matrix id
#  create a column which identifies all motifs belonging to the same secondary id
pwm<-fread('/data/annotation/ATtRACT/ATtRACT_db.pwm', sep='\t', header=F, fill=T, col.names=c('a_', 'c_', 'g_', 't_'))
n<-grep('>', pwm$a_) 
s<-sub('^>', '', rep(pwm[n, a_], c(diff(n), 1+nrow(pwm)-tail(n,1))))  #  repeat the matrix ids as many times as necessary to include given PWM entries
stopifnot( length(s)==nrow(pwm) )
stopifnot( length(s)==pwm[ grep('>', a_), sum(1+as.integer(c_))] )
pwm$id<-s
rm(s,n)


#  keep only the PWMs and not the FASTA-like name entries
#  split them to a list
#  keep only human
#  convert them to matrices holding numerical values and remove id column and column names
pwm<-pwm[ grep('>', a_, invert=T), ]
pwm<-split(pwm, by='id')
stopifnot( length(setdiff(db$id, names(pwm)))==0 )
pwm<-pwm[ db$id ]
pwm<-lapply(pwm, function(m){ unname(apply(as.matrix(m[,1:4]), 2, as.numeric)) })


#  SANITY CHECK: make sure all motifs and corresponding PWMs have the same length
#                ATtRACT database can have fucked up entries, e.g. 
#
#                    ELAVL2 AATTTATTTAA width=11 corresponds to the M329_0.6 PWM width of 9
#                    CNOT4 GACAGA width=6 corresponds to the M147_0.6 PWM with width of 7
TOREMOVE<-integer()
for(m in 1:nrow(db)){
    l<-db[m, nchar(motif)]
    if(nrow(pwm[[db[m, id]]])!=db[m, nchar(motif)]){
        cat('\n**** FUCKED UP motif:', db[m, gene_name], db[m, motif], db[m, id], 'removing... ****\n')
        TOREMOVE<-append(TOREMOVE, m)
    }
}
db<-db[-TOREMOVE, ]
pwm<-pwm[db[, id]]  #  keep only non-removed PWMs
rm(TOREMOVE)


#  SANITY CHECK: can the same motif_id have different PWMs?
p<-sapply(pwm, paste, sep='', collapse='') 
for(n in seq_along(p)){
    d<-which(p==p[n])
    if (!all(p[d]==p[d][1])){
        cat('FUCKED up motif_id:', names(p)[n], '\n')
    }
}  
#
#  this should check out and no motif_id should be associated with different PWMs!


#  collapse the PWMs to unique motif_id
pwm<-pwm[ match(unique(names(pwm)), names(pwm)) ]


#  identify the different motif_ids with identical PWMs
#  comma-separated merge of the motif_ids 
#  collapse the PWM to unique motif_ids again
#  replace motif_ids in the DB as well
p<-sapply(pwm, paste, sep='', collapse='') 
d<-list()
for(n in seq_along(p)){
    d[[names(p)[n]]]<-paste(names(p[ p==p[n] ]), sep='', collapse=',')
}
d<-unlist(d)
stopifnot( all.equal(names(p), names(d)) )
stopifnot( all.equal(names(p), names(pwm)) )
names(pwm)<-d
pwm<-pwm[ match(unique(names(pwm)), names(pwm)) ]
db[, id:=d[ id ]]


#  collect all RBPs under the same PWM
#  remove qscores
db<-db[, .(gene_name=list(unique(unlist((gene_name)))), gene_id=list(unique(unlist(gene_id))), motif=list(unique(unlist(motif)))), by=.(id)]


#  [obsolete and wrong]
#{{{
#  go over each motif separately and define subgroups based on the PWMs
#for(n in 1:nrow(db)){
#    p<-sapply(pwm[db[n, id][[1]]], paste, sep='', collapse='')                       #  quick identification of identical PWMs
#    i<-lapply(unname(split(names(p), p)), function(x){ match(x, db[n, id][[1]]) })   #  indices of the members of each group that have identical PWMs
#    db[n, c('gene_name', 'gene_id', 'id', 'qscore'):=list(gene_name=list(lapply(i, function(g){ unlist(gene_name)[g] })), 
#                                                          gene_id=list(lapply(i, function(g){ unlist(gene_id)[g] })),
#                                                          id=list(lapply(i, function(g){ unlist(id)[g] })),
#                                                          qscore=list(lapply(i, function(g){ unlist(qscore)[g] })))]
#}
#}}}


#  save
save(db, pwm, file='/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')

#}}}


#  load back
load('/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData')


#  manually create the MEME PWM motif file based on circRNA background nucleotide frequencies
#{{{

#  prepare the background model where we are going to read the nucleotide frequencies:
#
#      fasta-get-markov -m 0 -rna -norc /fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/unified/circRNAs_sequences.fa /data/annotation/ATtRACT/circRNA_sequences_background.model
#  
#  import them and write the MEME file header
bg<-read.table('/data/annotation/ATtRACT/circRNA_sequences_background.model', sep=' ', header=F)
bg<-setNames(bg[, 2], bg[, 1])
cat('MEME version 4\n\nALPHABET= ACGU\n\nBackground letter frequencies (from file \'/data/annotation/ATtRACT/circRNA_sequences_background.model\'):\n', paste(names(bg), bg, collapse=' '), '\n\n', file='/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.meme', sep='', append=F)
cat('MEME version 4\n\nALPHABET= ACGU\n\nBackground letter frequencies (from file \'/data/annotation/ATtRACT/circRNA_sequences_background.model\'):\n', paste(names(bg), bg, collapse=' '), '\n\n', file='/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme', sep='', append=F)
rm(bg)


#  for each unique PWM defined by the motif_id append the corresponding entry to the MEME motif file
for(n in 1:nrow(db)){
    #  grab the PWM
    p<-pwm[[db[n, id]]]
    
    
    #  define the motif names in the MEME file
    m.name<-sapply(db[n, gene_name], paste, sep='', collapse=',')
    

    #  define the alternative name from the consensus motifs
    a.name<-sapply(db[n, motif], paste, sep='', collapse=',')

    
    #  append the entries
    cat('MOTIF ', m.name, '_', db[n, id], ' ', a.name, '\n\n', sep='', file='/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.meme', append=T)
    cat('letter-probability matrix: ', 'alength= 4 w= ', nrow(p), ' nsites= 20 E=0\n', sep='', file='/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.meme', append=T)
    write.table(p, file='/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.meme', append=T, quote=F, row.names=F, col.names=F, eol='\n', sep=' ')
    cat('\n', sep='', file='/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.meme', append=T)
    
    
    #  write an entry to the strict DB if 6mer or longer
    if (nrow(p)>=6){
        cat('MOTIF ', m.name, '_', db[n, id], ' ', a.name, '\n\n', sep='', file='/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme', append=T)
        cat('letter-probability matrix: ', 'alength= 4 w= ', nrow(p), ' nsites= 20 E=0\n', sep='', file='/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme', append=T)
        write.table(p, file='/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme', append=T, quote=F, row.names=F, col.names=F, eol='\n', sep=' ')
        cat('\n', sep='', file='/data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme', append=T)
    }
}

#}}}

#}}}
#
#  => /data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.RData
#  => /data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db.meme
#  => /data/annotation/ATtRACT/Homo_sapiens_ATtRACT_db_strict.meme




###############################################################
#
#
#  mesenchymal/adrenergic gene set with GENCODE v30 annotations
#
#
###############################################################




#  load mesenchymal/adrenergic gene set
#  check gene names missing from GENCODE v30 and replace them with current gene names
#  add gene_id
#  save 
#{{{
rm(list=ls())
library(data.table)
library(rtracklayer)
library(GenomicAlignments)


#  import the reference
hsa<-import('/fast/groups/ag_schulte/work/reference/annotation/GRCh38/GRCh38.gencode.v30.gtf')
hsa<-hsa[ hsa$type %in% 'gene' ]
hsa<-data.table(data.frame(mcols(hsa)[, c('gene_name', 'gene_id')]))


#  load the original gene set
mes.adr<-fread('/fast/groups/ag_schulte/work/reference/annotation/MES_ADR_markers/Suppl Versteeg 2017 - MES ADR Genes.csv', header=F, sep=',', col.names=c('gene_name', 'cell_type'))


#  identify outdated names
#  search GeneCards (better than HUGO)
g<-setdiff( mes.adr$gene_name , hsa$gene_name )
#
#  C4orf32    =>  FAM241A
#  CYR61      =>  CCN1
#  FAM46A     =>  TENT5A
#  KIAA1462   =>  JCAD
#  KIRREL     =>  KIRREL1
#  LHFP       =>  LHFPL6
#  MRC2       =>  (for whatever reason MRC1 exists only, so MRC1 it is!)
#  NOTCH2NL   =>  NOTCH2NLA
#  PTRF       =>  CAVIN1
#  SPRY4-IT1  =>  (remove, it's basically SPRY4 which is part of the annotation)
#  ADRBK2     =>  GRK3
#  C7orf55    =>  FMC1
#  FAM60A     =>  SINHCAF
#  HMP19      =>  NSG2
#  HN1        =>  JPT1, NOTCH1 (add them both)
#  LOC100507194  =>  remove
#  LOC101928409  =>  remove
#
mes.adr[, gene_name:=lapply(gene_name, '[[', 1)]
mes.adr[ gene_name %in% 'C4orf32', gene_name:= 'FAM241A']
mes.adr[ gene_name %in% 'CYR61', gene_name:= 'CCN1']
mes.adr[ gene_name %in% 'FAM46A', gene_name:= 'TENT5A']
mes.adr[ gene_name %in% 'KIAA1462', gene_name:= 'JCAD']
mes.adr[ gene_name %in% 'KIRREL', gene_name:= 'KIRREL1']
mes.adr[ gene_name %in% 'LHFP', gene_name:= 'LHFPL6']
mes.adr[ gene_name %in% 'MRC2', gene_name:= 'MRC1']
mes.adr[ gene_name %in% 'NOTCH2NL', gene_name:= 'NOTCH2NLA']
mes.adr[ gene_name %in% 'PTRF', gene_name:= 'CAVIN1']
mes.adr[ gene_name %in% 'SPRY4-IT1', gene_name:= 'NA']
mes.adr[ gene_name %in% 'ADRBK2', gene_name:= 'GRK3']
mes.adr[ gene_name %in% 'C7orf55', gene_name:= 'FMC1']
mes.adr[ gene_name %in% 'FAM60A', gene_name:= 'SINHCAF']
mes.adr[ gene_name %in% 'HMP19', gene_name:= 'NSG2']
mes.adr[ gene_name %in% 'HN1', gene_name:= c('JPT1','NOTCH1')]
mes.adr[ gene_name %in% 'LOC100507194', gene_name:= 'NA']
mes.adr[ gene_name %in% 'LOC101928409', gene_name:= 'NA']
#
#mes.adr[, gene_name:=unlist(gene_name)]  #  this segmentfaults R!!!!
#
mes.adr<-mes.adr[, .(gene_name=unique(unlist(gene_name))), by=.(cell_type)][ ! gene_name %in% 'NA' ]


#  add gene_id
mes.adr<-hsa[ mes.adr, on='gene_name' ]  #  DIABLO has two gene_ids


#  save
save(mes.adr, file='/fast/groups/ag_schulte/work/reference/annotation/MES_ADR_markers/Suppl Versteeg 2017 - MES ADR Genes_gencode_v30.RData')

#}}}
#
#  => /fast/groups/ag_schulte/work/reference/annotation/MES_ADR_markers/Suppl Versteeg 2017 - MES ADR Genes_gencode_v30.RData




################################
#
#
#  scratch code and quick checks
#
#
################################




#  check MNA status of certain tumors
#{{{
rm(list=ls())
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(VariantAnnotation)
library(extrafont)  #  first time used need to run: font_import() , to load all for PDF device run: loadfonts(device='pdf')
loadfonts()


#  load metadata
load('/fast/projects/peifer_wgs/work/work/2017-11-08_Fuchs_totalRNAseq/raw/metadata.RData')
rm(meta.cel, meta.prefailed)


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
#  unlist back to GRanges for easy processing
#{{{

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
#  The filtered structural variant calls is located in doc/BerlinCohort_NovobreakResults_allChr_hg38.txt.
#  The additional calls based on hg19 are located under barcelona/.


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


#  unlist back
cnv<-unlist(cnv)

#}}}


#  MYCN: chr2+:15940550-15947007
#
#  we define a generous MYCN locus +-800Kbps that encompasses NBAS, DDX1 etc
mycn<-GRanges(seqnames='chr2', strand='*', ranges=IRanges(start=15140550, end=16747007))


#  check CNV on MYCN locus for particular samples
SAMPLES<-c('CB2028', 'CB3016', 'CB3036')
meta.tum[ PAT_ID_BERLIN %in% SAMPLES , c('PAT_ID_BERLIN', 'MYCN_INI_STATUS', 'MYCN_FC', 'risk_group')]
#
#    PAT_ID_BERLIN MYCN_INI_STATUS MYCN_FC risk_group
# 1:        CB2028         non-amp       1    HR_nMNA
# 2:        CB3016             amp     het        MNA
# 3:        CB3036             amp     het        MNA
#
x<-cnv[ cnv$sid %in% SAMPLES ]
unique(x$sid)  #  CB2028
x[queryHits(findOverlaps(x, mycn, type='any', select='all'))]
#
#          seqnames            ranges strand |    copies         sid    ploidy
#             <Rle>         <IRanges>  <Rle> | <integer> <character> <integer>
#   CB2028     chr2  3892246-15327585      * |         3      CB2028         2
#   CB2028     chr2 15327586-17497480      * |         4      CB2028         2
#   
mean(x[queryHits(findOverlaps(x, mycn, type='any', select='all'))]$copies)  #  3.5

#}}}






