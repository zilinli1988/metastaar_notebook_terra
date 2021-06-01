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
cloud_dir_store <- "gs://fc-282d7a75-f309-444a-9e1b-2679f8b67168/TOPMed_Summary_Stat/Cov/F8"

## output file name
output_file_name <- "GTSinvG.rare.HDL.F8"

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
if(!dir.exists("/home/jupyter-user/notebooks/Summary_Stat_Cov"))
{
	system("mkdir /home/jupyter-user/notebooks/Summary_Stat_Cov")
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

### Generate Summary Stat Cov
MetaSTAAR_worker_cov_agds <- function(i, cloud_dir, chr, genofile, variant.id, position, nullobj, segment_size, cov_maf_cutoff = 0.05)
{
	print(paste0("Chromosome: ", chr, "; Segment: ", i))
    
	### segment location 
	region_start_loc <- (i-1) * segment_size + 1
    region_midpos <- i * segment_size
    region_end_loc <- (i+1) * segment_size
	
	### phenotype id
	phenotype.id <- as.character(nullobj$id_include)	
	 
    is.in <- (SNVlist)&(position>=region_start_loc)&(position<=region_end_loc)
    variant_pos <- position[is.in]
    seqSetFilter(genofile,variant.id=variant.id[is.in],sample.id=phenotype.id)

    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
        
    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
	
	########################################################
	#                 Calculate MAF
	########################################################
		
	variant.id.sub <- seqGetData(genofile, "variant.id")
	### number of variants in each subsequence
	MAF_sub_snv_num <- 5000
	MAF_sub_seq_num <- ceiling(length(variant.id.sub)/MAF_sub_snv_num)

	genotype <- NULL
		    
	if(MAF_sub_seq_num > 0)
	{
		AF <- NULL
		      
		for(ii in 1:MAF_sub_seq_num)
		{
		    if(ii<MAF_sub_seq_num)
		    {
		       is.in.sub <- ((ii-1)*MAF_sub_snv_num+1):(ii*MAF_sub_snv_num)
		       invisible(capture.output(seqSetFilter(genofile,variant.id=variant.id.sub[is.in.sub],sample.id=phenotype.id)))
		    }
		        
		    if(ii==MAF_sub_seq_num)
		    {
		        is.in.sub <- ((ii-1)*MAF_sub_snv_num+1):length(variant.id.sub)
		        invisible(capture.output(seqSetFilter(genofile,variant.id=variant.id.sub[is.in.sub],sample.id=phenotype.id)))
		    }
		        
		    Geno <- seqGetData(genofile, "$dosage")
		    Geno <- Geno[id.genotype.match,,drop=FALSE]
		        
		    AF_sub <- apply(Geno,2,mean)/2
		    AF <- c(AF,AF_sub)
		        
		    rm(Geno)
		    gc()
		        
		    invisible(capture.output(seqResetFilter(genofile)))
		}
		MAF <- pmin(AF,1-AF)
		
		cov_maf_cutoff = 0.05
		### rare variant id
		RV_label <- (MAF<cov_maf_cutoff)&(MAF>0)
		variant.id.sub.rare <- variant.id.sub[RV_label]
		AF <- AF[RV_label]
		variant_pos <- variant_pos[RV_label]
		
		### Genotype
		RV_sub_num <- 5000
		RV_sub_seq_num <- ceiling(length(AF)/RV_sub_num)
		      
		for(jj in 1:RV_sub_seq_num)
		{
		    if(jj<RV_sub_seq_num)
		    {
				is.in.sub.rare <- ((jj-1)*RV_sub_num+1):(jj*RV_sub_num)
		        invisible(capture.output(seqSetFilter(genofile,variant.id=variant.id.sub.rare[is.in.sub.rare],sample.id=phenotype.id)))
		          
		        AF_sub <- AF[is.in.sub.rare]
		    }
		        
		    if(jj==RV_sub_seq_num)
		    {
		        is.in.sub.rare <- ((jj-1)*RV_sub_num+1):length(variant.id.sub.rare)
		        invisible(capture.output(seqSetFilter(genofile,variant.id=variant.id.sub.rare[is.in.sub.rare],sample.id=phenotype.id)))
		          
		        AF_sub <- AF[is.in.sub.rare]
		    }
		        
		    ## Genotype
		    Geno_sub <- seqGetData(genofile, "$dosage")
		    Geno_sub <- Geno_sub[id.genotype.match,,drop=FALSE]
		        
		    ## flip
		    if(sum(AF_sub>0.5)>0)
		    {
		        Geno_sub[,AF_sub>0.5] <- 2 - Geno_sub[,AF_sub>0.5]
		    }
		        
		    Geno_sub <- as(Geno_sub,"dgCMatrix")
		    genotype <- cbind(genotype,Geno_sub)
		        
		    rm(Geno_sub)
		    gc()
		        
		    invisible(capture.output(seqResetFilter(genofile)))
		}
	}
	
	GTSinvG_rare <- NULL
    try(GTSinvG_rare <- MetaSTAAR_worker_cov(genotype,obj_nullmodel = nullobj,cov_maf_cutoff = cov_maf_cutoff,variant_pos,
                                                 region_midpos,segment_size))

	## save results
	save(GTSinvG_rare,file = paste0("/home/jupyter-user/notebooks/Summary_Stat_Cov/",output_file_name,".chr",chr,".segment",i,".Rdata"),
		compress = "xz")

	## move results to google bucket
	system(paste0("gsutil cp /home/jupyter-user/notebooks/Summary_Stat_Cov/",output_file_name,".chr",chr,".segment",i,".Rdata ",cloud_dir,"/",output_file_name,".chr",chr,".segment",i,".Rdata"))

	## delete temp results to save storage
	system(paste0("rm /home/jupyter-user/notebooks/Summary_Stat_Cov/",output_file_name,".chr",chr,".segment",i,".Rdata"))

	seqResetFilter(genofile)
    
	rm(genotype,GTSinvG_rare)
    gc()
}

a <- Sys.time()
mclapply(1:segment.num, MetaSTAAR_worker_cov_agds, cloud_dir=cloud_dir_store, chr=chr, genofile=genofile, variant.id=variant.id, 
         position=position, nullobj=nullobj, segment_size=segment.size, cov_maf_cutoff = 0.05, mc.cores=user_cores)
b <- Sys.time()
## Computation time
b - a

seqClose(genofile)
	 
