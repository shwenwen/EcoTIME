library(dplyr)
library(Seurat)
library(Matrix)

ST_preprocess <- function(outDir, study_ID, file_suffix, genome) {
  # Preprocess ST data
  #
  # Args:
  #   outDir: Diretory for data saving. Choose from: '/fs/home/dongzhonghua/STARDUST/Data/' or '/fs/home/shixiaoying/Project/ST_LLM/Data/'.
  #   study_ID: GSE number. 
  #   file_suffix: Suffix of files. Example: '_gene_count.h5', '_filtered_feature_bc_matrix.h5', '.h5' and so on.
  #   genome: Genome version. Choose from: 'GRCh38', 'GRCh37', 'GRCm38', 'GRCm37'.
  # 
  # Return:
  #   QC-filtered .rds file.
  #   QC-filtered .h5 file.
  
  resDir <- file.path(outDir, study_ID, "suppl")
  allh5 <- list.files(resDir, pattern = '*\\.h5', recursive = T, full.names = T)
  suffix_list <- unique(c('_gene_count.h5', '_filtered_feature_bc_matrix.h5', '.h5', file_suffix))
  samples <- gsub(paste(suffix_list,  collapse= '|'), '', basename(allh5))
  print(samples)
  
  for (sample in samples) {
    print(sample)
    samplePath = allh5[grep(sample, allh5)]
    
    # Read in h5 data
    h5_file <- Read10X_h5(filename = samplePath)
    ST_object <- CreateSeuratObject(counts = h5_file, project = sample)
    
    # Filter out low quality spots
    ST_object[["percent.mt"]] <- PercentageFeatureSet(object = ST_object, pattern = "^MT-")
    ST_object <- subset(ST_object, subset = nFeature_RNA > 200 & percent.mt < 20)
    ST_object <- SCTransform(ST_object, assay = "RNA", verbose = FALSE)
    ST_object <- FindVariableFeatures(ST_object, selection.method = "vst", nfeatures = 2000)
    ST_object <- RunPCA(ST_object, features = VariableFeatures(object = ST_object))
    ST_object <- FindNeighbors(ST_object, dims = 1:30)
    ST_object <- FindClusters(ST_object, resolution = 0.5)
    ST_object <- RunUMAP(ST_object, dims = 1:30)

    # Save the .rds file for future analysis.
    saveRDS(ST_object, file = file.path(resDir, paste0(sample, "_filtered.rds")))

    # Save the QC-filtered .h5 file for STRIDE deconvolution.
    writeh5(h5path = file.path(resDir, paste0(sample, "_gene_count_QC.h5")), 
            count_data = ST_object@assays$RNA@data, 
            genome_assembly = genome)
    
  }
}

writeh5 <- function(h5path, count_data, genome_assembly){
  library(rhdf5)
  library(hdf5r)
  
  h5createFile(h5path)
  h5createGroup(h5path,"matrix")
  h5write(count_data@Dimnames[[2]] , h5path, "matrix/barcodes")
  h5write(as.integer(count_data@x), h5path, "matrix/data")
  h5createGroup(h5path,"matrix/features")
  h5write("genome", h5path, "matrix/features/_all_tag_keys")
  Genes <- rep("Gene Expression", length(count_data@Dimnames[[1]]))
  h5write(Genes,h5path, "matrix/features/feature_type")
  Genome <- rep(genome_assembly, length(count_data@Dimnames[[1]])) 
  h5write(Genome, h5path,"matrix/features/genome")
  h5write(count_data@Dimnames[[1]],h5path, "matrix/features/id")
  h5write(count_data@Dimnames[[1]],h5path, "matrix/features/name")
  h5write(count_data@i, h5path, "matrix/indices") 
  h5write(count_data@p, h5path, "matrix/indptr")
  h5write(count_data@Dim, h5path, "matrix/shape")
  h5closeAll()
}

# Usage
outDir <- "/fs/home/dongzhonghua/STARDUST/Data"
study_ID <- "GSE210041"
file_suffix <- "_filtered_feature_bc_matrix.h5"
genome <- "GRCh38"
ST_preprocess(outDir, study_ID, file_suffix, genome)
