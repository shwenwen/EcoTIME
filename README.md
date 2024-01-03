# Workflow to analysis ST data

## 0 Download ST data

This pipeline is mainly used to analyze sequencing-based spatial transcriptome data such as 10X Visium or Slide-seq data. Spot-gene expression matrix files and spot-location files are necessary. In some cases, spot-location information is included in spot-gene expression matrix as colnames or rownames of the matrix.

Spot meta information and HE images are also needed for some downstream analysis, you can download them from GEO supplementary file or original paper if available.

## 1 Preprocess ST data

### 1.1 Load ST data into R

ST expression data may have multiple file formats, such as HDF5(h5), plain text, and mtx format. This pipeline mainly uses h5 format. Expression files in other formats can be converted into h5 files through MAESTRO.

The following commands can be used for format conversion.&#x20;

    MAESTRO mtx-to-h5	      #Convert mtx format to HDF5 format
    MAESTRO count-to-h5       #Convert plain text format to HDF5 format

See <https://github.com/liulab-dfci/MAESTRO> for more information, or enter 'MAESTRO --help' on the command line to ask for help.

We use R package "Seurat" to preprocess ST data. First we need load ST expression data from HDF5 file and create seurat object.

    ST_expr <- Read10X_h5(filename = "Data/GSE210041/suppl/GSM6415705_GLMF1_filtered_feature_bc_matrix.h5")
    ST_object <- CreateSeuratObject(counts = ST_expr, project = "GSM6415705_GLMF1")

### 1.2 Quality control

We perform quality control to filter out low quality spots. Generally, we consider spots with less than 200 features or more than 20% mitochondrial gene percentage to be of low quality. You can adjust the quality control conditions according to the characteristics of the data.

    ST_object[["percent.mt"]] <- PercentageFeatureSet(object = ST_object, pattern = "^MT-")
    ST_object <- subset(ST_object, subset = nFeature_RNA > 200 & percent.mt < 20)

### 1.3 Normalization

    ST_object <- SCTransform(ST_object, assay = "RNA", verbose = FALSE)

### 1.4 Dimensionality reduction

    ST_object <- FindVariableFeatures(ST_object, selection.method = "vst", nfeatures = 2000)
    ST_object <- RunPCA(ST_object, features = VariableFeatures(object = ST_object))

### 1.5 Clustering based on expression

    ST_object <- FindNeighbors(ST_object, dims = 1:30)
    ST_object <- FindClusters(ST_object, resolution = 0.5)

### 1.6 Visualization

    ST_object <- RunUMAP(ST_object, dims = 1:30)

### 1.7 Load images slot

To import spatial image information, create a folder named 'sample\_spatial'. This folder should contain both images and spot location information. Acceptable image formats include 'tissue\_hires\_image.png', 'tissue\_lowres\_image.png', 'sample\_HE.tiff', 'sample\_HE.png' and 'sample\_HE.jpeg'. Acceptable location information formats include 'tissue\_positions\_list.csv' and 'sample\_spot\_location.txt'. Additionally, we can also read in 'scalefactors\_json.json' for additional image information if available.ST\_object

'Read\_image' function is used to load images slot into ST\_object, related code can be found in 'ST\_preprocess.R'

    image <- Read_image(sample = sample, spatial_dir = spatial_dir)
    DefaultAssay(object = image) <- 'RNA'
    image <- image[colnames(x = ST_object)]
    ST_object[[sample]] <- image

### 1.8 Save rds result and HDF5 file

Save preprocessed data for future analysis.&#x20;

    saveRDS(ST_object, file = file.path(resDir, paste0(sample, "_filtered.rds")))
    writeh5(h5path = file.path(resDir, paste0(sample, "_gene_count_QC.h5")), count_data = h5_object@assays$RNA@data, genome_assembly = genome)

**Note**: If you need to process multiple samples or multiple datasets, you can also perform batch processing through code in 'ST\_preprocess.R'.

## 2 Intergrate scRNA data and ST data

### 2.1 Download scRNA-seq data

Sequencing-based spatial transcriptome data has low resolution, and each spot contains several cells. In order to explore the cell type composition of each spot, we need to integrate ST data and single-cell data. In this pipeline, We use STRIDE to deconvolve and map ST data. The cell type composition of single-cell data will have a great impact on the deconvolution results, so it is best to use paired scRNA-seq data as reference. When paired scRNA-seq data is unavaliable, we can use scRNA data of the same cancer type or from the same tissue as reference.&#x20;

Single cell data can be found from the TISCH2 databases, or from GEO databases, published literatures, and other source.

### 2.2 Preprocess scRNA-seq data

We should perform basic analysis on scRNA data and prepare input file required by STRIDE. The scRNA data process pipeline can be found at <https://github.com/ytwang21/TISCH2/tree/master/code>. STRIDE deconvolution usually need three scRNA input files, including scRNA expression matrix file, celltype curated file, marker gene list file. Marker gene file is not necessary, if marker gene file is not provided, STRIDE deconvolution will find markers automatically but will cost longer time.&#x20;

### 2.3 Deconvolution by STRIDE

#### 2.3.1 Deconvolution

See <https://github.com/wanglabtongji/STRIDE> for more information, or enter 'STRIDE deconvolve --help' for help.

    STRIDE deconvolve --sc-count Data/GSE213699/suppl/scRNA/OV\_GSE154600\_count\_gene\_count.h5 \ 
    --sc-celltype Data/GSE213699/suppl/scRNA/OV\_GSE154600\_celltype\_curated.txt \ 
    --st-count Data/GSE213699/suppl/${sample}_gene_count_QC.h5 \ 
    --gene-use /Data/GSE213699/suppl/scRNA/OV_GSE154600_scRNA_top_marker_list.txt \ 
    --outdir Data/GSE213699/suppl/stride_result/${sample} 
    --outprefix ${sample} --normalize;


#### 2.3.2 Evaluate deconvolution results

We need to evaluate the results of deconvolution by calculating the correlation between the gene signature score and the deconvolution faction. Gene list for signature score is derived from single cell data. Generally, if the correlation is greater than 0.3, we believe that the deconvolution result meets the requirements. Otherwise, we need to re-run STRIDE deconvolution or change the single cell references. Gene signature score and correlation can be calculated according to the following code.

### 2.4 Mapping by STRIDE

In order to perform downstream spatial CCI and spatial GRN analysis, we need to map single cells to spatial locations.

#### 2.4.1 Mapping



#### 2.4.2 Single-cell transcriptome data with spatial location





## 3 Spatial domain



### 3.1 Spatial domains



### 3.2 Spatial variable genes (SVGs)



## 4 Spatial cell-cell interaction

## 5 Spatial
