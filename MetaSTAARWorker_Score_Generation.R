rm(list=ls())
gc()

############################################################
#        Load R-packages
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
## chr
chr <- 19

## google bucket storing aGDS
aGDS_dir <- "gs://fc-secure-19d7e6e8-18bb-4649-9e66-b2ac84049a85/F8_aGDS"

## google bucket storing null model
null_model_dir <- "gs://fc-282d7a75-f309-444a-9e1b-2679f8b67168/Null_Model"

## Null model file name
null_model_name <- "obj.SMMAT.HDL.fulladj.sp.20190912.Rdata"

## google bucket storing output
cloud_dir_store <- "gs://fc-282d7a75-f309-444a-9e1b-2679f8b67168/TOPMed_Summary_Stat/Score/F8"

## output file name
output_file_name <- "summary.stat.HDL.F8"

## number of cores
user_cores <- 24


############################################################
#                    Preparation Step
############################################################

##### load Null model
if(!dir.exists("/home/jupyter-user/notebooks/Null_Model"))
{
	system("mkdir /home/jupyter-user/notebooks/Null_Model")
}

if(!file.exists(paste0("/home/jupyter-user/notebooks/Null_Model/",null_model_name)))
{
	system(paste0("gsutil cp ",null_model_dir,"/",null_model_name," /home/jupyter-user/notebooks/Null_Model/",null_model_name))
}
nullobj <- get(load(paste0("/home/jupyter-user/notebooks/Null_Model/",null_model_name)))

###### open aGDS
if(!dir.exists("/home/jupyter-user/notebooks/aGDS"))
{
	system("mkdir /home/jupyter-user/notebooks/aGDS")
}

if(!file.exists(paste0("/home/jupyter-user/notebooks/aGDS/freeze.8.chr",chr,".pass_and_fail.gtonly.minDP0.gds")))
{
	system(paste0("gsutil cp ",aGDS_dir,"/freeze.8.chr",chr,".pass_and_fail.gtonly.minDP0.gds /home/jupyter-user/notebooks/aGDS/freeze.8.chr",chr,".pass_and_fail.gtonly.minDP0.gds"))
}

### output folder
if(!dir.exists("/home/jupyter-user/notebooks/Summary_Stat_Score"))
{
	system("mkdir /home/jupyter-user/notebooks/Summary_Stat_Score")
}


######################################################
#                 Main Step
######################################################
### gds file
dir.geno <- "/home/jupyter-user/notebooks/aGDS"
gds.path <- paste(dir.geno,"/freeze.8.chr",chr,".pass_and_fail.gtonly.minDP0.gds",sep="")
genofile <- seqOpen(gds.path)

## get SNV id
filter <- seqGetData(genofile, "annotation/filter")
AVGDP <- seqGetData(genofile, "annotation/info/AVGDP")
SNVlist <- filter == "PASS" & AVGDP > 10 & isSNV(genofile)
rm(filter,AVGDP)
gc()

variant.id <- seqGetData(genofile, "variant.id")

## Position
position <- as.integer(seqGetData(genofile, "position"))
max_position <- max(position)

segment.size <- 5e5
segment.num <- ceiling(max_position/segment.size)

###  Generate Summary Stat Score
MetaSTAAR_worker_score_agds <- function(i, cloud_dir, chr, genofile, variant.id, position, nullobj, segment_size)
{
	print(paste0("Chromosome: ", chr, "; Segment: ", i))

	subsegment_num <- 25
	summary_stat <- NULL

    for(j in 1:subsegment_num)
	{
		### segment location 
		region_start_loc <- (i-1) * segment_size + (j-1) * (segment_size/subsegment_num) + 1
		region_end_loc <- (i-1) * segment_size + j * (segment_size/subsegment_num)
		
		### phenotype id
		phenotype.id <- as.character(nullobj$id_include)				
        
		is.in <- (SNVlist)&(position>=region_start_loc)&(position<=region_end_loc)
        seqSetFilter(genofile,variant.id=variant.id[is.in],sample.id=phenotype.id)
                
        pos <- as.integer(seqGetData(genofile, "position"))
        ref <- unlist(lapply(strsplit(seqGetData(genofile, "allele"),","),`[[`,1))
        alt <- unlist(lapply(strsplit(seqGetData(genofile, "allele"),","),`[[`,2))

		## genotype id
		id.genotype <- seqGetData(genofile,"sample.id")
        
		id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
		phenotype.id.merge <- data.frame(phenotype.id)
		phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
		id.genotype.match <- phenotype.id.merge$index

	    ##### Filtering all variants
        genotype <- seqGetData(genofile, "$dosage")
        genotype <- genotype[id.genotype.match,,drop=FALSE]
	
		if(!is.null(genotype))
		{
            variant_info <- data.frame(chr,pos,ref,alt)
                        
            results_temp <- NULL
            try(results_temp <- MetaSTAAR_worker_sumstat(genotype,nullobj,variant_info))
            summary_stat <- rbind(summary_stat,results_temp)
        }
	}
	
	## save results
	save(summary_stat,file = paste0("/home/jupyter-user/notebooks/Summary_Stat_Score/",output_file_name,".chr",chr,".segment",i,".Rdata"),
             compress = "xz")
	
	## move results to google bucket
	system(paste0("gsutil cp /home/jupyter-user/notebooks/Summary_Stat_Score/",output_file_name,".chr",chr,".segment",i,".Rdata ",cloud_dir,"/",output_file_name,".chr",chr,".segment",i,".Rdata"))

	## delete temp results to save storage
	system(paste0("rm /home/jupyter-user/notebooks/Summary_Stat_Score/",output_file_name,".chr",chr,".segment",i,".Rdata"))

    seqResetFilter(genofile)
    
	rm(genotype,results_temp,pos,ref,alt,variant_info,summary_stat)
    gc()
}

### Generate Summary Stat Score
a <- Sys.time()
mclapply(1:segment.num, MetaSTAAR_worker_score_agds, cloud_dir=cloud_dir_store, chr=chr, genofile=genofile, variant.id=variant.id, 
         position=position, nullobj=nullobj, segment_size=segment.size, mc.cores=user_cores)
b <- Sys.time()
## Computation time
b - a

seqClose(genofile)
	 
