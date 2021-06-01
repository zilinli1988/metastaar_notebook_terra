rm(list=ls())
gc()

############################################################
#                     Load R-packages
############################################################
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(MetaSTAAR)
library(Matrix)
library(dplyr)
library(parallel)

############################################################
#                     User Input
############################################################
## study names for meta-analysis
study.names <- c("F8")
cov_maf_cutoff <- c(0.05)
## trait name
trait <- "LDL"
## sample size of each study
sample.sizes <- 61861
## directory of sparse weighted LD matrix
cov.dir <- c("gs://fc-282d7a75-f309-444a-9e1b-2679f8b67168/TOPMed_Summary_Stat/Cov/")
## directory of score statistics
sumstat.dir <- c("gs://fc-282d7a75-f309-444a-9e1b-2679f8b67168/TOPMed_Summary_Stat/Score/")


#### Sliding window parameter
## chr
chr <- 19
## start location of the chromosome
start_loc <- 60547
## end location of the chromosome
end_loc <- 58607590
## sliding window length
sliding_window_length <- 2000
## sliding window skip length
skip_length <- 1000

#### number of cores
user_cores <- 24

## google bucket storing output
cloud_dir_store <- "gs://fc-282d7a75-f309-444a-9e1b-2679f8b67168/TOPMed_F8_RVAT/"

## output file name
output_file_name <- "Sliding.Window.LDL.F8"

############################################################
#                    Preparation Step
############################################################
file.number <- c(498,485,397,381,363,342,319,291,277,268,271,267,229,214,204,181,167,161,118,129,94,102)

##### create folders
### Cov
if(!dir.exists("/home/jupyter-user/notebooks/Cov"))
{
	system("mkdir /home/jupyter-user/notebooks/Cov")
}

for(kk in 1:length(study.names))
{
	if(!dir.exists("/home/jupyter-user/notebooks/Cov"))
	{
		system(paste0("mkdir /home/jupyter-user/notebooks/Cov/",study.names[kk]))
	}
}

### Move files
for(kk in 1:length(study.names))
{
	for(segment in 1:file.number[chr])
	{
		if(!file.exists(paste0("/home/jupyter-user/notebooks/Cov/",study.names[kk],"/GTSinvG.rare.",trait,".",study.names[kk],".chr",chr,".segment",segment,".Rdata")))
		{
			system(paste0("gsutil cp ",cov.dir,study.names,"/GTSinvG.rare.",trait,".",study.names[kk],".chr",chr,".segment",segment,".Rdata /home/jupyter-user/notebooks/Cov/",study.names[kk],"/GTSinvG.rare.",trait,".",study.names[kk],".chr",chr,".segment",segment,".Rdata"))
		}
	}
}


### Score
if(!dir.exists("/home/jupyter-user/notebooks/Score"))
{
	system("mkdir /home/jupyter-user/notebooks/Score")
}

for(kk in 1:length(study.names))
{
	if(!dir.exists("/home/jupyter-user/notebooks/Score"))
	{
		system(paste0("mkdir /home/jupyter-user/notebooks/Score/",study.names[kk]))
	}
}

### Move files
for(kk in 1:length(study.names))
{
	for(segment in 1:file.number[chr])
	{
		if(!file.exists(paste0("/home/jupyter-user/notebooks/Score/",study.names[kk],"/summary.stat.",trait,".",study.names[kk],".chr",chr,".segment",segment,".Rdata")))
		{
			system(paste0("gsutil cp ",sumstat.dir,study.names,"/summary.stat.",trait,".",study.names[kk],".chr",chr,".segment",segment,".Rdata /home/jupyter-user/notebooks/Score/",study.names[kk],"/summary.stat.",trait,".",study.names[kk],".chr",chr,".segment",segment,".Rdata"))
		}
	}
}


#################################################################
#                       Main Step
################################################################# 
MetaSTAAR_RVAT_SW <- function(kk, chr, start_loc, end_loc, sliding_window_length, skip_length, study.names, sample.sizes, trait, cov_maf_cutoff)
{
	sliding_window_num <- floor((end_loc - start_loc)/skip_length) + 1

	subset.window.num <- 50
	subset.num <- floor(sliding_window_num/subset.window.num) + 1

	print(paste("kk =", kk))
    sub_start_loc <- start_loc + (kk-1)*skip_length*subset.window.num
    sub_end_loc <- sub_start_loc + (subset.window.num - 1)*skip_length + sliding_window_length
    # To avoid 'No such file or directory' error
    if (sub_end_loc > end_loc) {
        sub_end_loc <- end_loc
    }
	
	cov.dir <- "/home/jupyter-user/notebooks/Cov/"
	sumstat.dir <- "/home/jupyter-user/notebooks/Score/"
	
	sub_window_merge <- MetaSTAAR_merge(chr=chr,start.loc=sub_start_loc,end.loc=sub_end_loc,study.names=study.names,sample.sizes=sample.sizes,sumstat.dir=sumstat.dir,cov.dir=cov.dir,trait=trait,cov_maf_cutoff=cov_maf_cutoff)

	results <- c()
	if(!is.null(sub_window_merge))
    {
      if((length(sub_window_merge$U)>=2)&(sum(sub_window_merge$info$MAF<0.01)>=2))
      {
        for(k in 1:subset.window.num)
        {  
          region_start_loc <- sub_start_loc + (k-1)*skip_length
          region_end_loc <- region_start_loc + sliding_window_length - 1
          
          region.id <- (sub_window_merge$info$pos>=region_start_loc)&(sub_window_merge$info$pos<=region_end_loc)
          region_window_merge <- list(info = sub_window_merge$info[region.id,],
                                      U = sub_window_merge$U[region.id],
                                      cov = sub_window_merge$cov[region.id,region.id])
          
          results_temp <- c(chr,region_start_loc,region_end_loc)
          
          ## weights               
          pvalues <- 0
          try(pvalues <- MetaSTAAR(region_window_merge,annotation_phred=NULL))
          
          if(class(pvalues)=="list")
          {
            results_temp <- c(results_temp,pvalues$num_variant,pvalues$results_STAAR_S_1_25[1],pvalues$results_STAAR_S_1_1[1],
                              pvalues$results_STAAR_B_1_25[1],pvalues$results_STAAR_B_1_1[1],pvalues$results_STAAR_A_1_25[1],
                              pvalues$results_STAAR_A_1_1[1],pvalues$results_STAAR_O)
            
            results <- rbind(results,results_temp)
          }
		}
		if(!is.null(results))
		{
			colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
			colnames(results)[1:4] <- c("Chr","Start Loc","End Loc","#SNV")
			colnames(results)[11] <- c("STAAR-O")
		}
      }
      rm(sub_window_merge)
      gc()
    }
	return(results)
}

sliding_window_num <- floor((end_loc - start_loc)/skip_length) + 1
subset.window.num <- 50
subset.num <- floor(sliding_window_num/subset.window.num) + 1

a <- Sys.time()
out <- mclapply(1:subset.num, MetaSTAAR_RVAT_SW, chr=chr, start_loc=start_loc, end_loc=end_loc, 
				sliding_window_length=sliding_window_length, skip_length=skip_length, study.names=study.names, 
				sample.sizes=sample.sizes, trait=trait, cov_maf_cutoff=cov_maf_cutoff,mc.cores=user_cores)
b <- Sys.time()
## Computation time
b - a

results_sliding_window_chr <- do.call("rbind",out)

save(results_sliding_window_chr,file = paste0("/home/jupyter-user/notebooks/",output_file_name,".chr",chr,".Rdata"))

system(paste0("gsutil cp /home/jupyter-user/notebooks/",output_file_name,".chr",chr,".Rdata ",cloud_dir_store,output_file_name,".chr",chr,".Rdata"))




