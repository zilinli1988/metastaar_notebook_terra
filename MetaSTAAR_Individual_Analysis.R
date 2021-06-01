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
## trait name
trait <- "LDL"
## sample size of each study
sample.sizes <- 61861
## directory of score statistics
sumstat.dir <- c("gs://fc-282d7a75-f309-444a-9e1b-2679f8b67168/TOPMed_Summary_Stat/Score/")


#### Sliding window parameter
## chr
chr <- 19
## start location of the chromosome
start_loc <- 60547
## end location of the chromosome
end_loc <- 58607590
## MAF cut-off
common_maf_cutoff <- 1E-04

#### number of cores
user_cores <- 24

## google bucket storing output
cloud_dir_store <- "gs://fc-282d7a75-f309-444a-9e1b-2679f8b67168/TOPMed_F8_RVAT/"

## output file name
output_file_name <- "Individual.LDL.F8"

############################################################
#                    Preparation Step
############################################################
file.number <- c(498,485,397,381,363,342,319,291,277,268,271,267,229,214,204,181,167,161,118,129,94,102)

##### create folders
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
MetaSTAAR_merge_common <- function(chr,start.loc,end.loc,study.names,sample.sizes,sumstat.dir,common_maf_cutoff,trait,segment.size = 5e5){

  cov_maf_cutoff <- rep(0.5,length(study.names))		
  segment <- floor((start.loc - 1) / segment.size) + 1
  if (end.loc <= segment * segment.size) {
    ### summary statistics
    sumstat.files <- paste0(sumstat.dir,study.names,"/summary.stat.",trait,".",study.names,".chr",chr,".segment",segment,".Rdata")
    sumstat.list <- lapply(sumstat.files, function(x) {
      load(file = x)
      get(ls()[ls()!= "summary.stat"])
    })
    position.index <- mapply(function(x,y) {
      (x[(x$MAF<y)&(x$MAF>0),"pos"]>=start.loc)&(x[(x$MAF<y)&(x$MAF>0),"pos"]<=end.loc)
    }, x = sumstat.list, y = cov_maf_cutoff, SIMPLIFY = FALSE)
    sumstat.list <- lapply(sumstat.list, function(x) {
      x[(x$pos>=start.loc)&(x$pos<=end.loc),]
    })
    sumstat.varid.list <- mapply(function(x,y) {
      x[(x$MAF<y)&(x$MAF>0),1:4]
    }, x = sumstat.list, y = cov_maf_cutoff, SIMPLIFY = FALSE)
    sumstat.varid.merge <- do.call("rbind",sumstat.varid.list)
    sumstat.varid.nodup <- sumstat.varid.merge[!duplicated(sumstat.varid.merge),]
    if (is.null(sumstat.varid.nodup)) {
      return(NULL)
    }else if (dim(sumstat.varid.nodup)[1] == 0) {
      return(NULL)
    }
    sumstat.merge.list <- lapply(sumstat.list, function(x) {
      if (is.null(x)) {
        cbind(sumstat.varid.nodup,alt_AC=NA,MAC=NA,MAF=NA,N=NA,U=NA,V=NA,"1"=NA)
      }else {
        left_join(sumstat.varid.nodup,x,by=c("chr"="chr",
                                             "pos"="pos",
                                             "ref"="ref",
                                             "alt"="alt"))
      }
    })
    sumstat.merge.list <- mapply(function(x,y) {
      x[is.na(x[,"N"]),"N"] <- y
      x[is.na(x)] <- 0
      return(x)
    }, x = sumstat.merge.list, y = sample.sizes, SIMPLIFY = FALSE)
    sumstat.varid.nodup$index <- 1:dim(sumstat.varid.nodup)[1]
    sumstat.index.list <- lapply(sumstat.varid.list, function(x) {
      if (is.null(x)) {
        integer(0)
      }else{
        left_join(x,sumstat.varid.nodup,by=c("chr"="chr",
                                             "pos"="pos",
                                             "ref"="ref",
                                             "alt"="alt"))$index
      }
    })
    rm(sumstat.files,sumstat.list,sumstat.varid.list,sumstat.varid.merge)
    gc()

  }else if (end.loc <= (segment + 1) * segment.size) {
    ### summary statistics
    sumstat.files1 <- paste0(sumstat.dir,study.names,"/summary.stat.",trait,".",study.names,".chr",chr,".segment",segment,".Rdata")
    sumstat.list1 <- lapply(sumstat.files1, function(x) {
      load(file = x)
      get(ls()[ls()!= "summary.stat"])
    })
    position.index1 <- mapply(function(x,y) {
      (x[(x$MAF<y)&(x$MAF>0),"pos"]>=start.loc)&(x[(x$MAF<y)&(x$MAF>0),"pos"]<=end.loc)
    }, x = sumstat.list1, y = cov_maf_cutoff, SIMPLIFY = FALSE)
    sumstat.list1 <- lapply(sumstat.list1, function(x) {
      x[(x$pos>=start.loc)&(x$pos<=end.loc),]
    })
    sumstat.varid.list1 <- mapply(function(x,y) {
      x[(x$MAF<y)&(x$MAF>0),1:4]
    }, x = sumstat.list1, y = cov_maf_cutoff, SIMPLIFY = FALSE)

    sumstat.files2 <- paste0(sumstat.dir,study.names,"/summary.stat.",trait,".",study.names,".chr",chr,".segment",segment+1,".Rdata")
    sumstat.list2 <- lapply(sumstat.files2, function(x) {
      load(file = x)
      get(ls()[ls()!= "summary.stat"])
    })
    position.index2 <- mapply(function(x,y) {
      (x[(x$MAF<y)&(x$MAF>0),"pos"]>=start.loc)&(x[(x$MAF<y)&(x$MAF>0),"pos"]<=end.loc)
    }, x = sumstat.list2, y = cov_maf_cutoff, SIMPLIFY = FALSE)
    sumstat.list2 <- lapply(sumstat.list2, function(x) {
      x[(x$pos>=start.loc)&(x$pos<=end.loc),]
    })
    sumstat.varid.list2 <- mapply(function(x,y) {
      x[(x$MAF<y)&(x$MAF>0),1:4]
    }, x = sumstat.list2, y = cov_maf_cutoff, SIMPLIFY = FALSE)

    sumstat.list <- mapply(function(x,y) {
      rbind(x,y)
    }, x = sumstat.list1, y = sumstat.list2, SIMPLIFY = FALSE)
    sumstat.varid.list <- mapply(function(x,y) {
      rbind(x,y)
    }, x = sumstat.varid.list1, y = sumstat.varid.list2, SIMPLIFY = FALSE)
    sumstat.varid.merge <- do.call("rbind",sumstat.varid.list)
    sumstat.varid.nodup <- sumstat.varid.merge[!duplicated(sumstat.varid.merge),]
    if (is.null(sumstat.varid.nodup)) {
      return(NULL)
    }else if (dim(sumstat.varid.nodup)[1] == 0) {
      return(NULL)
    }
    sumstat.merge.list <- lapply(sumstat.list, function(x) {
      if (is.null(x)) {
        cbind(sumstat.varid.nodup,alt_AC=NA,MAC=NA,MAF=NA,N=NA,U=NA,V=NA,"1"=NA)
      }else {
        left_join(sumstat.varid.nodup,x,by=c("chr"="chr",
                                             "pos"="pos",
                                             "ref"="ref",
                                             "alt"="alt"))
      }
    })
    sumstat.merge.list <- mapply(function(x,y) {
      x[is.na(x[,"N"]),"N"] <- y
      x[is.na(x)] <- 0
      return(x)
    }, x = sumstat.merge.list, y = sample.sizes, SIMPLIFY = FALSE)
    sumstat.varid.nodup$index <- 1:dim(sumstat.varid.nodup)[1]
    sumstat.index.list1 <- lapply(sumstat.varid.list1, function(x) {
      if (is.null(x)) {
        integer(0)
      }else{
        left_join(x,sumstat.varid.nodup,by=c("chr"="chr",
                                             "pos"="pos",
                                             "ref"="ref",
                                             "alt"="alt"))$index
      }
    })
    sumstat.index.list2 <- lapply(sumstat.varid.list2, function(x) {
      if (is.null(x)) {
        integer(0)
      }else{
        left_join(x,sumstat.varid.nodup,by=c("chr"="chr",
                                             "pos"="pos",
                                             "ref"="ref",
                                             "alt"="alt"))$index
      }
    })
    rm(sumstat.files1,sumstat.files2,sumstat.list,sumstat.list1,sumstat.list2,
       sumstat.varid.list,sumstat.varid.list1,sumstat.varid.list2,
       sumstat.varid.merge)
    gc()
  }

  ### select rare variant based on the input cutoff
  alt_AC.merge <- as.integer(Reduce("+",lapply(sumstat.merge.list, function(x) {x$alt_AC})))
  N.merge.nonzero <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N * (x$MAC != 0)}))
  alt_AF.merge.nonzero <- alt_AC.merge / (2 * N.merge.nonzero)
  N.merge.zero <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N * (x$MAC == 0)}))
  alt_AC.merge <- alt_AC.merge + (alt_AF.merge.nonzero > 0.5) * (2 * N.merge.zero)
  N.merge <- Reduce("+",lapply(sumstat.merge.list, function(x) {x$N}))
  MAC.merge <- pmin(alt_AC.merge, 2 * N.merge - alt_AC.merge)
  MAF.merge <- MAC.merge / (2 * N.merge)
  cv.index <- (MAF.merge>common_maf_cutoff) & Reduce("*",mapply(function(x,y) {
    x$MAF<y
  }, x = sumstat.merge.list, y = cov_maf_cutoff, SIMPLIFY = FALSE))

  info <- cbind(sumstat.varid.nodup[,c("chr","pos","ref","alt")],
                MAC=MAC.merge,MAF=MAF.merge)[cv.index,]

  U.merge <- Reduce("+", lapply(sumstat.merge.list, function(x) {
    x <- x[cv.index,]
    return(x$U)
  }))
  
  V.merge <- Reduce("+", lapply(sumstat.merge.list, function(x) {
    x <- x[cv.index,]
    return(x$V)
  }))

  p.merge <- pchisq(U.merge^2/V.merge,df=1,lower.tail=FALSE)
  logp.merge <- -pchisq(U.merge^2/V.merge,df=1,lower.tail=FALSE,log.p=TRUE)	

  rm(list=setdiff(ls(), c("info","U.merge","V.merge","p.merge","logp.merge")))
  gc()

  return(data.frame(info=info,
              U=U.merge,
              V=V.merge,
			  p=p.merge,
			  logp=logp.merge))
}


MetaSTAAR_Individual <- function(kk, chr, start_loc, end_loc, study.names, sample.sizes, trait, common_maf_cutoff)
{
	subwindow_length <- 100000
	subwindow_num <- floor((end_loc - start_loc)/subwindow_length) + 1

	print(paste("kk =", kk))
    sub_start_loc <- start_loc + (kk-1)*subwindow_length
    sub_end_loc <- sub_start_loc + subwindow_length - 1
    # To avoid 'No such file or directory' error
    if (sub_end_loc > end_loc) {
        sub_end_loc <- end_loc
    }
	
	sumstat.dir <- "/home/jupyter-user/notebooks/Score/"
	
	common_merge <- MetaSTAAR_merge_common(chr=chr,start.loc=sub_start_loc,end.loc=sub_end_loc,study.names=study.names,sample.sizes=sample.sizes,sumstat.dir=sumstat.dir,trait=trait,common_maf_cutoff=common_maf_cutoff)

	return(common_merge)
}

subwindow_length <- 100000
subwindow_num <- floor((end_loc - start_loc)/subwindow_length) + 1

a <- Sys.time()
out <- mclapply(1:subwindow_num, MetaSTAAR_Individual, chr=chr, start_loc=start_loc, end_loc=end_loc, 
				study.names=study.names, sample.sizes=sample.sizes, trait=trait, common_maf_cutoff=common_maf_cutoff, mc.cores=user_cores)
b <- Sys.time()
## Computation time
b - a

results_individual_chr <- do.call("rbind",out)

save(results_sliding_window_chr,file = paste0("/home/jupyter-user/notebooks/",output_file_name,".chr",chr,".Rdata"))

system(paste0("gsutil cp /home/jupyter-user/notebooks/",output_file_name,".chr",chr,".Rdata ",cloud_dir_store,output_file_name,".chr",chr,".Rdata"))




