# Workflow to analysis ST data

## 0 Download ST data

This pipeline is mainly used to analyze sequencing-based spatial transcriptome data such as 10X Visium or Slide-seq data. Spot-gene expression matrix file and spot-location file are necessary. In some cases, spot-location information is included in spot-gene expression matrix as colnames or rownames of the matrix.

Spot meta information and HE image are also needed for some downstream analysis, you can download them from GEO supplementary file or original paper if available.

## 1 Preprocess ST data

### 1.1 Load ST data into R

ST expression data may have multiple file formats, such as HDF5, plain text, and mtx format. This pipeline mainly uses h5 format. Expression files in other formats can be converted into h5 files through MAESTRO.

The following commands can be used for format conversion.

    MAESTRO mtx-to-h5	      #Convert mtx format to HDF5 format
    MAESTRO count-to-h5       #Convert plain text format to HDF5 format

We use R package "Seurat" to preprocess ST data. First we need load ST expression data from HDF5 file and create seurat object.

    ST_expr <- Read10X_h5(filename = "/fs/home/dongzhonghua/STARDUST/Data/GSE210041/suppl/GSM6415705_GLMF1_filtered_feature_bc_matrix.h5")
    ST_object <- CreateSeuratObject(counts = ST_expr, project = "GSM6415705_GLMF1")

### 1.2 Quality control

We perform quality control to filter out low quality spots. Generally, we consider spots with less than 200 features or more than 20% mitochondrial gene percentage to be of low quality. You can adjust the quality control conditions according to the characteristics of the data.

```
ST_object[["percent.mt"]] <- PercentageFeatureSet(object = ST_object, pattern = "^MT-")
ST_object <- subset(ST_object, subset = nFeature_RNA > 200 & percent.mt < 20)
    
```

### 1.3 Normalization

    ST_object <- SCTransform(ST_object, assay = "RNA", verbose = FALSE)

### 1.4 Dimensionality reduction

    ST_object <- FindVariableFeatures(ST_object, selection.method = "vst", nfeatures = 2000)
    ST_object <- RunPCA(ST_object, features = VariableFeatures(object = ST_object))

### 1.5 Clustering based on spatial expression

    ST_object <- FindNeighbors(ST_object, dims = 1:30)
    ST_object <- FindClusters(ST_object, resolution = 0.5)

### 1.6 Visualization

    ST_object <- RunUMAP(ST_object, dims = 1:30)

**Note**: If you need to process multiple samples or multiple datasets, you can also perform batch processing through the code in ST\_preprocess.R.

## &#x20;2 Intergrate scRNA data and ST data

### 2.1 Download scRNA-seq data

Sequencing-based spatial transcriptome data has low resolution, and each capture point contains multiple cells. In order to explore the cell type composition of each spot, we need to integrate ST data and single-cell data. In this pipeline, We use STRIDE to deconvolve spatial transcriptome data. The cell type composition in single-cell data will have a greater impact on the deconvolution results, so it is best to use single-cell data that is matched with idle data. However, this In rare cases, we can use single-cell data of the same cancer type or from the same tissue as a reference. In short, the closer the single-cell data and the spatial transcriptome data are, the better the deconvolution result will be.

Single cell data can be found from the TISCH2 website previously published by our laboratory, or from GEO , published literatures, and other public databases.

### 2.2 Preprocess scRNA-seq data

After downloading the single cell data, we need to do some basic processing and prepare the single cell data input file required by STRIDE. The scRNA data process pipeline can be found at <https://github.com/ytwang21/TISCH2/tree/master/code>
