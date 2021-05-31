library(edgeR, quietly=TRUE)
library(limma, quietly=TRUE)
library(tidyr, quietly=TRUE)
library(dplyr, quietly=TRUE)

mode = ifelse(length(commandArgs(trailingOnly=TRUE))==1, 'debug', 'batch')

if (mode=="debug") {
  #dir_work = '/Users/s229181/MSN/'
  #dir_ortho = paste0(dir_work, "OrthoFinder/Results_Feb09_2/Orthogroups")
  #dir_count = paste0(dir_work, "counts/")
  dir_work = '/Users/kef74yk/Dropbox_p/collaborators/Ken Naito/20210509_Vigna/gfe_data'
  dir_ortho = "/Users/kef74yk/Dropbox_p/collaborators/Ken Naito/20210509_Vigna/gfe_data/Orthogroups"
  dir_count = "/Users/kef74yk/Dropbox_p/collaborators/Ken Naito/20210509_Vigna/gfe_data/merge"
  setwd(dir_work)
} else if (mode=="batch") {
  args = commandArgs(trailingOnly=TRUE)
  dir_count = args[1]
  dir_ortho = args[2]
  dir_work = args[3]
}

get_spp_filled = function(dir_count, df_gc) {
  spp_filled = list.files(path=dir_count, pattern=".*count.\\.tsv")
  spp_filled = sub('_', '|', spp_filled)
  spp_filled = sub('_.*', '', spp_filled)
  spp_filled = sub('\\|', '_', spp_filled)
  is_missing_in_genecount = (!spp_filled %in% colnames(df_gc))
  if (sum(is_missing_in_genecount)) {
    for (sp in spp_filled[is_missing_in_genecount]) {
      warning(paste0('Species excluded. Not found in OrthoFinder\'s GeneCount table: ', sp))
    }
  }
  spp_filled = spp_filled[!is_missing_in_genecount]
  cat('Detected species:', spp_filled, '\n')
  return(spp_filled)
}

get_singlecopy_og = function(df_gc, spp_filled) {
  is_singlecopy = TRUE
  for (sp in spp_filled) {
    is_singlecopy = is_singlecopy & (df_gc[,sp]==1)
  }
  sc_og = df_gc[is_singlecopy,'Orthogroup']
  cat(length(sc_og), 'single-copy orthogroups were detected for the', length(spp_filled), 'species.\n')
  return(sc_og)
}

# set directory
if (!file.exists(dir_work)) {
  dir.create(dir_work)
}
setwd(dir_work)

# set output directory and create if not already there
dir_tmm = file.path(dir_work, 'cross_species_tmm_normalized_counts')
if (!file.exists(dir_tmm)) {
  dir.create(dir_tmm)
}

file_genecount = file.path(dir_ortho, 'Orthogroups.GeneCount.tsv')
file_orthogroup = file.path(dir_ortho, 'Orthogroups.tsv')

df_gc = read.table(file_genecount, header=TRUE, sep='\t', check.names=FALSE)
spp_filled = get_spp_filled(dir_count, df_gc)
single_orthogroups = get_singlecopy_og(df_gc, spp_filled)
df_og = read.table(file_orthogroup, header=TRUE, sep='\t', row.names=1, check.names=FALSE)
df_singleog = df_og[(rownames(df_og) %in% single_orthogroups), spp_filled, drop=FALSE]
df_og = NULL

uncorrected = list()
df_sog = df_singleog
for (sp in spp_filled) {
  infile = list.files(path = dir_count, pattern=paste0(sp,".*count.\\.tsv"))
  if (length(infile)> 1){
    stop(paste0("Multiple *count.tsv files found for: ", sp ,"\n"))
  } else if (length(infile)==0) {
    warning(paste0("Skipping. No *count.tsv files found for: ", sp ,"\n"))
    next
  }
  infile_path = file.path(dir_count, infile[1])
  if (file.exists(infile_path)) {
    cat('Input file found, reading:', infile[1], '\n')
    dat = read.delim(infile_path, header = T, row.names=1, sep='\t', check.names=FALSE)
    dat = dat[,(colnames(dat)!='length')]
    genus_cap = tolower(substring(sp,1,1))
    spec_name = strsplit(sp,'_')[[1]][2]
    colnames(dat) = paste(sp, colnames(dat), sep='_')
    uncorrected[[sp]] = dat
    df_sog = merge(df_sog, dat, by.x=sp, by.y="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)
    if (all(is.na(df_sog[,colnames(dat)]))) {
      warning('Gene IDs may not match between OrthoFinder outputs and amalgkit merge outputs. Please check.')
    }
  } else {
    cat('Input file not found:', infile, '\n')
  }
}
df_sog = df_sog[,-(1:length(spp_filled))]
rownames(df_sog) = rownames(df_singleog)

is_na_containing_row = apply(df_sog, 1, function(x){any(is.na(x))})
cat('Removing', sum(is_na_containing_row), 'out of', nrow(df_sog), 'orthogroups because missing values are observed in at least 1 species.\n')
is_no_count_col = apply(df_sog, 2, function(x){sum(x, na.rm=TRUE)==0})
cat('Removing', sum(is_no_count_col), 'out of', ncol(df_sog), 'RNA-seq samples because no read mapping values are all zero.\n')
df_nonzero = df_sog[!is_na_containing_row,!is_no_count_col]
cnf_in = edgeR::DGEList(counts = df_nonzero)

cnf_out1 = edgeR::calcNormFactors(cnf_in, method='TMM', refColumn=NULL)
x = cnf_out1[[2]][['norm.factors']]
cat('Round 1: Median normalization factor =', median(x), '\n')
median_value = sort(x)[ceiling(length(x)/2)]
median_index = (1:length(x))[x==median_value]

cnf_out2 = edgeR::calcNormFactors(cnf_in, method='TMM', refColumn=median_index)
cat('Round 2: Median normalization factor =', median(cnf_out2[[2]][['norm.factors']]), '\n')

df_nf = cnf_out2[[2]]
df_nf[['sample']] = rownames(df_nf)
df_nf = df_nf[,c('sample','group','lib.size','norm.factors')]
write.table(df_nf, 'normalization_factor.tsv', row.names=FALSE, sep='\t')

xlim = c(-2,2)
bins = seq(-2,2,0.1)
file_name='normalization_factor_histogram.pdf'
pdf(file_name, height=3.3, width=7.2) # full figure size = 9.7 x 7.2
x = log2(cnf_out2[[2]][['norm.factors']])
x[x>xlim[2]] = xlim[2]
x[x<xlim[1]] = xlim[1]
hist(x, xlab='log2(TMM normalization factor)', ylab='Count', main=NULL, col='black', xlim=xlim, breaks=bins, las=1)
graphics.off()

for (sp in names(uncorrected)) {
  cat('Applying TMM normalization factors:', sp, '\n')
  dat = uncorrected[[sp]]
  df_nf_sp = cnf_out2[[2]][startsWith(rownames(cnf_out2[[2]]),sp),]

  for (i in 1:length(df_nf_sp[,1])){
    # check if dat colnames start with the species name
    # and modify SRR variable to match the format of dat colnames
    # TODO: is this necessary? Is there any situation where the species name prefix is missing?
    if(all(startsWith(colnames(dat), sp))){
      SRR = as.character(row.names(df_nf_sp[i,]))
    }
    else{
      SRR = as.character(strsplit(row.names(df_nf_sp[i,]),'_')[[1]][3])
    }
    # manually apply normfactor
    tmm_normalization_factor = as.double(df_nf_sp[i,"norm.factors"])
    dat[,SRR] = dat[,SRR]/tmm_normalization_factor

  }
  dat_out = cbind(target_id=rownames(dat), dat)
  rownames(dat_out) = NULL
  colnames(dat_out) = sub(paste0(sp, '_'), '', colnames(dat_out))
  file_path = file.path(dir_tmm, paste0(sp, "_cstmm_counts.tsv"))
  write.table(dat_out, file_path, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)

}
cat('Done!\n')
