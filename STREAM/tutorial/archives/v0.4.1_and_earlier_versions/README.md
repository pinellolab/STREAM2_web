[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/stream/README.html)

[![Build Status](https://travis-ci.org/pinellolab/STREAM.svg)](https://travis-ci.org/pinellolab/STREAM)

# STREAM (Latest version v0.4.1)

Latest News
-----------
> Jan 14, 2020

Version 0.4.1 is now available. We added support of feature `top_pcs` for *Mapping*

> Nov 26, 2019

Version 0.4.0 is now available. Numerous changes have been introduced. Please check [v0.4.0](https://github.com/pinellolab/STREAM/releases/tag/v0.4.0) for details.

Introduction
------------

STREAM (**S**ingle-cell **T**rajectories **R**econstruction, **E**xploration **A**nd **M**apping) is an interactive pipeline capable of disentangling and visualizing complex branching trajectories from both single-cell transcriptomic and epigenomic data.

STREAM is now published in *Nature Communications*! Please cite our paper [Chen H, et al. Single-cell trajectories reconstruction, exploration and mapping of omics data with STREAM.](https://www.nature.com/articles/s41467-019-09670-4)  *Nature Communications*, volume 10, Article number: 1903 (2019). if you find STREAM helpful for your research.

<img src="https://github.com/pinellolab/STREAM/blob/stream_python2/STREAM/static/images/Figure1.png">

STREAM is written using the class `anndata` [Wolf et al. Genome Biology (2018)](http://anndata.rtfd.io) and available as user-friendly open source software and can be used interactively as a web-application at [stream.pinellolab.org](http://stream.pinellolab.org/), as a bioconda package [https://bioconda.github.io/recipes/stream/README.html](https://bioconda.github.io/recipes/stream/README.html) and as a standalone command-line tool with Docker [https://github.com/pinellolab/STREAM](https://github.com/pinellolab/STREAM)

Installation with Bioconda (Recommended)
----------------------------------------

1)	If Anaconda (or miniconda) is already installed with **Python 3**, skip to 2) otherwise please download and install Python3 Anaconda from here: https://www.anaconda.com/download/

2)	Open a terminal and add the Bioconda channel with the following commands:

```sh
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```

3)	Create an environment named `myenv` , install **stream**, **jupyter**, and activate it with the following commands:

* *For single cell **RNA-seq** analysis*:
```sh
$ conda create -n myenv python stream jupyter
$ conda activate myenv
```
* *For single cell **ATAC-seq** analysis*:
```sh
$ conda create -n myenv python stream stream_atac jupyter
$ conda activate myenv
```

4)  To perform STREAM analyis in Jupyter Notebook as shown in **Tutorial**, type `jupyter notebook` within `myenv`:

```sh
$ jupyter notebook
```

You should see the notebook open in your browser.

Tutorial
--------

* Example for scRNA-seq: [1.1.STREAM_scRNA-seq (on 2D plane).ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/archives/v0.4.1_and_earlier_versions/1.1.STREAM_scRNA-seq%20%28on%202D%20plane%29.ipynb?flush_cache=true)

* Example for scRNA-seq: [1.2.STREAM_scRNA-seq (in 3D space).ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/archives/v0.4.1_and_earlier_versions/1.2.STREAM_scRNA-seq%20%28in%203D%20space%29.ipynb?flush_cache=true)

* Example for *mapping* feature: [2.STREAM_mapping.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/archives/v0.4.1_and_earlier_versions/2.STREAM_mapping.ipynb?flush_cache=true)

* Example for complex trajectories: [3.STREAM_complex_trajectories.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/archives/v0.4.1_and_earlier_versions/3.STREAM_complex_trajectories.ipynb?flush_cache=true)

* Example for scATAC-seq(using k-mers): [4.STREAM_scATAC-seq_k-mers.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/archives/v0.4.1_and_earlier_versions/4.STREAM_scATAC-seq_k-mers.ipynb?flush_cache=true)

* Example for scATAC-seq(using motifs): [5.STREAM_scATAC-seq_motifs.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/archives/v0.4.1_and_earlier_versions/5.STREAM_scATAC-seq_motifs.ipynb?flush_cache=true)

Installation with Docker
------------------------

With Docker no installation is required, the only dependence is Docker itself. Users will completely get rid of all the installation and configuration issues. Docker will do all the dirty work for you!

Docker can be downloaded freely from here: [https://store.docker.com/search?offering=community&type=edition](https://store.docker.com/search?offering=community&type=edition)

To get an image of STREAM, simply execute the following command:

```sh
$ docker pull pinellolab/stream
```

Basic usage of *docker run* 

```sh
$ docker run [OPTIONS] IMAGE [COMMAND] [ARG...]
```

OPTIONS:  
```
--publish , -p	Publish a container’s port(s) to the host  
--volume , -v	Bind mount a volume  
--workdir , -w	Working directory inside the container  
```

STREAM interactive website
--------------------------

In order to make STREAM user friendly and accessible to non-bioinformatician, we have created an interactive website: [http://stream.pinellolab.org](http://stream.pinellolab.org)

The website can also run on a local machine. More details can be found [https://github.com/pinellolab/STREAM_web](https://github.com/pinellolab/STREAM_web)


STREAM command line interface
-----------------------------

Please note that **STREAM command line is a streamlined script of 'stream' bioconda package. Some parameters and functions are not supported in STREAM command line.** To get a more flexible and advanced analysis, please check out our bioconda tutorial notebooks.

To run STREAM at the command-line interface:

* start a terminal session;

* enter ```docker run  -v ${PWD}:/data -w /data  pinellolab/stream --help [options]```

Users can specify the following options:
```
-m, --matrix  
input file name. Matrix is in .tsv or tsv.gz format in which each row represents a unique gene and each column is one cell. (default: None)
-l, --cell_labels  
file name of cell labels (default: None)
-c, --cell_labels_colors  
file name of cell label colors (default: None)
-s, --select_features  
LOESS, PCA, all: Select variable genes using LOESS or top principal components using PCA or keep all the gene (default: LOESS)
--TG  
detect transition genes automatically
--DE  
detect DE genes automatically
--LG  
etect leaf genes automatically
-g, --gene_list  
genes to visualize, it can either be filename which contains all the genes in one column or a set of gene names separated by comma (default: None)
-p, --use_precomputed  
use precomputed data files without re-computing structure learning part
--log2  
perform log2 transformation
--norm  
normalize data based on library size
--atac
indicate scATAC-seq data  
--n_processes  
Specify the number of processes to use. (default, all the available cores).
--loess_frac  
The fraction of the data used in LOESS regression (default: 0.1)
--pca_first_PC  
keep first PC
--pca_n_PC  
The number of selected PCs (default: 15)
--n_processes  
Specify the number of processes to use. The default uses all the cores available
--lle_neighbours  
LLE neighbour percent (default: 0.1)
--lle_components  
Number of components for LLE space (default: 3)
--clustering
Clustering method used for seeding the intial structure, choose from 'ap','kmeans','sc'.  
--damping  
Affinity Propagation: damping factor (default: 0.75)
--n_clusters  
Number of clusters for spectral clustering or kmeans
--EPG_n_nodes
Number of nodes for elastic principal graph (default: 50)
--EPG_lambda
lambda parameter used to compute the elastic energy (default: 0.02)
--EPG_mu
mu parameter used to compute the elastic energy (default: 0.1)
--EPG_trimmingradius
maximal distance of point from a node to affect its embedment (default: Inf)
--EPG_alpha  
positive numeric, alpha parameter of the penalized elastic energy (default: 0.02)
--disable_EPG_optimize  
disable optimizing branching  
--EPG_collapse  
Collapsing small branches
--EPG_collapse_mode  
the mode used to collapse branches. It can be 'PointNumber','PointNumber_Extrema', 'PointNumber_Leaves','EdgesNumber' or 'EdgesLength' (default:'PointNumber')
--EPG_collapse_par  
the control parameter used for collapsing small branches
--EPG_shift
shift branching point 
--EPG_shift_mode  
the mode to use to shift the branching points 'NodePoints' or 'NodeDensity' (default: NodeDensity)
--EPG_shift_DR  
positive numeric, the radius used when computing point density if EPG_shift_mode is 'NodeDensity' (default:0.05)
--EPG_shift_maxshift  
positive integer, the maximum distance (number of edges) to consider when exploring the branching point neighborhood (default:5)
--disable_EPG_ext  
disable extending leaves with additional nodes
--EPG_ext_mode  
the mode used to extend the graph. It can be 'QuantDists', 'QuantCentroid' or 'WeigthedCentroid'. (default: QuantDists)
--EPG_ext_par  
the control parameter used for contribution of the different data points when extending leaves with nodes (default: 0.5)
--DE_zscore_cutoff  
Differentially Expressed Genes z-score cutoff (default: 2)
--DE_logfc_cutoff  
Differentially Expressed Genes log fold change cutoff (default: 0.25)  
--TG_spearman_cutoff  
Transition Genes Spearman correlation cutoff (default: 0.4)
--TG_logfc_cutoff  
Transition Genes log fold change cutoff (default: 0.25)
--LG_zscore_cutoff  
Leaf Genes z-score cutoff (default: 1.5)
--LG_pvalue_cutoff  
Leaf Genes p value cutoff (default: 1e-2)
--umap  
Whether to use UMAP for visualization (default: No)  
-r
root node for subwaymap_plot and stream_plot (default:None)  
--stream_log_view
use log2 scale for y axis of stream_plot 
--for_web
Output files for website
-o, --output_folder  
Output folder (default: None)
--new  
file name of data to be mapped (default: None)
--new_l  
filename of new cell labels (default: None)
--new_c  
filename of new cell label colors (default: None)
```


Input file format
-----------------
#### **Transcriptomic data**

The main and required input file is a tab-separated gene expression matrix (raw counts or normalized values) in tsv file format. Each row represents a unique gene and each column is one cell.


For example, in python
```R
>import pandas as pd
>input_data = pd.read_csv('data_Guo.tsv.gz',sep='\t',header=0,index_col=0)
>input_data.iloc[0:5,0:5]
```

|        | HSPC_025      | HSC1.1    | HSPC_037    | LT-HSC_001    | HSPC_001   |
|--------|-----------|-----------|-----------|-----------|----------|
| Clec1b   | 0.000000  | 0.000000  | 0.000000  | 0.000000  | 0.000000 |
| Kdm3a | 4.891604 | 6.877725 | 0.000000 | 0.000000 | 0.000000 |
| Coro2b  | 1.426148  | 0.000000  | 6.913384  | 8.178374  | 9.475577 |
| 8430408G22Rik   | 0.000000 | 0.000000 | 0.000000  | 0.000000  | 0.000000 |
| Clec9a    | 0.000000  | 0.000000 | 0.000000  | 0.000000  | 0.000000 |

In addition, it is possible to provide these optional files in .tsv or .tsv.gz format: 

**cell_labels** file: .tsv or .tsv.gz format. Each item can be a putative cell type or sampling time point obtained from experiments. Cell labels are helpful for visually validating the inferred trajectory. The order of labels should be consistent with cell order in the gene expression matrix file. No header is necessary:

|        |
|--------|
| MPP    | 
| MPP    | 
| MPP    |
| HSC    |
| MPP    |

**cell_label_color** file: .tsv or .tsv.gz format. Customized colors to use for the different cell labels. The first column specifies cell labels and the second column specifies the color in the format of hex. No header is necessary:

|       |         | 
|-------|---------|
| HSC   | #40bdbd | 
| MPP   | #eea113 |
| CMP   | #d84f40 |
| GMP   | #10b460 |
| MEP   | #286ee1 |
| LMPP   | #7c5246 |

**gene_list** file: .tsv or .tsv.gz format. It contains genes that users may be interested in visualizing in subway map and stream plot. Genes are listed in one column. No header is necessary: 

|        |
|--------|
| Gata1 | 
| Mpo  | 
| Car2   |
| Prtn3   |
| Tmsb4x  |


#### **Epigenomic data**
To perform scATAC-seq trajectory inference analysis, the main input can be:   

1)a set of files including **count file**, **region file** and **sample file**. 

**count file**, .tsv or .tsv.gz format. A tab-delimited triplet file. It contains three columns. The first column specifies the rows indices (the regions, 1-based indexing) for non-zero entry. The second column specifies the columns indices (the samples, 1-based indexing) for non-zero entry. The last column contains the number of reads in a given region for a particular cell. No header is necessary:

|        |     |  |
|--------|-----|--|
| 3735   | 96  | 1|
| 432739 | 171 | 2|
| 133126 | 292 | 1|
| 219297 | 359 | 1|
| 284936 | 1222| 1|
| 442588 | 1580| 2|

**region file**, .bed or .bed.gz format. A tab-delimited .bed file with three columns. The first column specifies chromosome names. The second column specifies the start position of the region. The third column specifies the end position of the region. The order of regions should be consistent with the regions indices in the count file. No header is necessary:

|      |       |      |
|------|-------|------|
| chr1 | 10279 | 10779|
| chr1 | 13252 | 13752|
| chr1 | 16019 | 16519|
| chr1 | 29026 | 29526|
| chr1 | 96364 | 96864|

**sample file**, .tsv or .tsv.gz format. It has one column. Each row is a cell name.  The order of the cells should be consistent with the sample indices in count file. No header is necessary:

|                                    |
|------------------------------------|
| singles-BM0828-HSC-fresh-151027-1  | 
| singles-BM0828-HSC-fresh-151027-2  | 
| singles-BM0828-HSC-fresh-151027-3  |
| singles-BM0828-HSC-fresh-151027-4  |
| singles-BM0828-HSC-fresh-151027-5  |


2)a precomputed scaled z-score file by STREAM. Each row represents a k-mer DNA sequence. Each column represents one cell. Each entry is a scaled z-score of the accessibility of each k-mer across cells

|        | singles-BM0828-HSC-fresh-151027-1 | singles-BM0828-HSC-fresh-151027-2 | singles-BM0828-HSC-fresh-151027-3 |
|--------|-----------------------------------|-----------------------------------|-----------------------------------|
| AAAAAAA|-0.15973157637808505               | 0.18950966450007853               | 0.07713107176524692               | 
| AAAAAAG|-1.3630723054479532                | -0.04770034004421244              | 0.6387323857481045                |
| AAACACG|-0.2065161126378667                | -1.3375384076872765               | 0.2660278729402342                |
| AGCGTTA|-0.496859947462221                 | 0.7181918229050274                | 0.19603357892921522               |
| ATACTCA|-1.2127919166377426                | 0.7938414496478844                | -1.2665513250104594               |



Example
--------

All the datasets used in the following examples can also be downloaded using the following Dropbox link: 
[https://www.dropbox.com/sh/n8qq4m7w17i6b07/AAAro_qY_-q5VBDC1sZg-LE5a?dl=0](https://www.dropbox.com/sh/n8qq4m7w17i6b07/AAAro_qY_-q5VBDC1sZg-LE5a?dl=0)

Please note that for large dataset analysis it'll be necessary to increase the default allocated memory of container.

<img src="https://github.com/pinellolab/STREAM/blob/stream_python2/STREAM/static/images/docker.png" width="30%">

### **Transcriptomic data**

Here we we take a single cell RNA-seq dataset as an example,including data_Nestorowa.tsv.gz, cell_label.tsv.gz and cell_label_color.tsv.gz (Nestorowa, S. et al.,2016), and assuming that **they are in the current folder**, to perform trajectory inference analysis, users can simply run a single command:

*Using Bioconda:*
```sh
$ stream -m data_Nestorowa.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz
```
*Using Docker:*
```sh
$ docker run  -v ${PWD}:/data -w /data  pinellolab/stream -m data_Nestorowa.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz
```

If cell labels are not available or no customized cell label color file is available, **-l** or **-c** can also be omitted

*Using Bioconda:*
```sh
$ stream -m data_Nestorowa.tsv.gz
```
*Using Docker:*
```sh
$ docker run  -v ${PWD}:/data -w /data  pinellolab/stream -m data_Nestorowa.tsv.gz
```

To visualize genes of interest, user can provide a gene list file by adding **-g**, for example: gene_list.tsv.gz. Meanwhile, by adding the flag  **-p**, STREAM will use the precomputed file obtained from the first running (In this way, STREAM will import precomupted pkl file so the analysis will skip structure learning part and only execute the step of visualizing genes):

*Using Bioconda:*
```sh
$ stream -m data_Nestorowa.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz -g gene_list.tsv.gz -p
```
*Using Docker:*
```sh
$ docker run  -v ${PWD}:/data -w /data  pinellolab/stream -m data_Nestorowa.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz -g gene_list.tsv.gz -p
```

Users can also provide a set of gene names separated by comma or specify the root by adding **-r**:

*Using Bioconda:*
```sh
$ stream -m data_Nestorowa.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz -g Gata1,Mpo -r S1 -p
```
*Using Docker:*
```sh
$ docker run  -v ${PWD}:/data -w /data  pinellolab/stream -m data_Nestorowa.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz -g Gata1,Mpo -r S1 -p
```

To explore potential marker genes, it is possible to add the flags **--DE**, **--TG**, or **--LG** to detect DE (differentially expressed) genes, transition gens, and leaf genes respectively:

*Using Bioconda:*
```sh
$ stream -m data_Nestorowa.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz --DE --TG --LG -p
```
*Using Docker:*
```sh
$ docker run  -v ${PWD}:/data -w /data  pinellolab/stream -m data_Nestorowa.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz --DE --TG --LG -p
```

### **Mapping**

To explore the feature **mapping**, users need to provide two dataset, one is used for inferring trajectories. The other is the dataset that is going to be mapped to the inferred trajectories. Here we take data_Olsson.tsv.gz, data_perturbation.tsv (Olsson, A. et al.,2016) as an example. We assume that **all the datasets are in the current folder**.

Users first need to run the following command to get initial inferred trajetories from wild-type cells:

*Using Bioconda:*
```sh
$ stream -m data_Olsson.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz --lle_components 4 --EPG_shift 
```
*Using Docker:*
```sh
$ docker run  -v ${PWD}:/data -w /data  pinellolab/stream -m data_Olsson.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz --lle_components 4 --EPG_shift 
```

To map the genetically perturbed cells to the inferred trajectories, users can execute the following command:

*Using Bioconda:*
```sh
$ stream --new data_perturbation.tsv.gz --new_l cell_perturbation_label.tsv.gz --new_c cell_perturbation_label_color.tsv.gz 
```
*Using Docker:*
```sh
$ docker run  -v ${PWD}:/data -w /data  pinellolab/stream --new data_perturbation.tsv.gz --new_l cell_perturbation_label.tsv.gz --new_c cell_perturbation_label_color.tsv.gz 
```
After running this command,  a folder named **'mapping_result'** will be created under the current directory along with all the mapping analysis results.


### **scATAC-seq data**

To perform scATAC-seq trajectory inference analysis, three files are necessary, a .tsv file of counts in compressed sparse format, a sample file in .tsv format and a region file in .bed format. (Buenrostro, J.D. et al., 2018). We assume that **they are in the current folder**.

Using these three files, users can run `stream_atac` with the following command to preprocess sc-atac-seq data and get a z_score matrix file named **'zscore.tsv.gz'** (This step may take a couple of hours with a modest machine):

*Using Bioconda:*
```sh
$ stream_atac -c count_file.tsv.gz -s sample_file.tsv.gz -r region_file.bed.gz
```

Then, take z-score file as input to infer trajectories using `stream`:

*Using Bioconda:*
```sh
$ stream --atac -m zscore.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz --lle_components 4
```
*Using Docker:*
```sh
$ docker run  -v ${PWD}:/data -w /data pinellolab/stream --atac -m zscore.tsv.gz -l cell_label.tsv.gz -c cell_label_color.tsv.gz --lle_components 4
```

Output description
------------------

STREAM write all the results by default in the folder **stream_result**, unless a different directory is specified by the user with the flag **-o**. This folder mainly contains the following files and directories:

*   **std_vs_means.pdf**: selected most variable genes.
*   **dimension_reduction.pdf**: projected cells in the MLLE 3D space.
*   **seed_elastic_principal_graph_skeleton.pdf**: the initial structure skeleton with all the nodes and edges.
*   **seed_elastic_principal_graph.pdf**: the initial structure with cells.
*   **ElPiGraph_analysis.pdf**: the log of ElPiGraph strucuture learning.
*   **elastic_principal_graph_skeleton.pdf**: the elastic principal graph skeleton.
*   **elastic_principal_graph.pdf**: the elastic principal graph with cells.
*   **optimizing_elastic_principal_graph_skeleton.pdf**: the elastic principal graph skeleton after optimizing branching.
*   **optimizing_elastic_principal_graph.pdf**: the elastic principal graph with cells after optimizing branching.
*   **extending_elastic_principal_graph_skeleton.pdf**: the elastic principal graph with cells after extending leaf nodes.
*   **extending_elastic_principal_graph.pdf**: the elastic principal graph skeleton after extending leaf nodes.
*   **finalized_elastic_principal_graph_skeleton.pdf**: the finalized elastic principal graph skeleton.
*   **finalized_elastic_principal_graph.pdf**: the finalized elastic principal graph with cells.
*   **flat_tree.pdf**: flat tree plot.
*   **cell_info.tsv**: cell information file containing branch assignment id and pseudotime.
*	**stream_result.pkl**: stores anndata object from the analysis. It can be imported later to reproduce the whole analysis.
*   sub-folder **'transition_genes'** contains several files, one for each branch id, for example for (S0,S1):
    - **transition_genes_S0_S1.pdf**: Detected transition genes plot for branch S0_S1. Orange bars are genes whose expression values increase from state S0 to S1 and green bars are genes whose expression values decrease from S0 to S1
    - **transition_genes_S0_S1.tsv**: Table that stores information of detected transition genes for branch S1_S2.
*   sub-folder **'de_genes'** contains several files, one for each pair of branches, for example for (S0,S1) and (S0,S2):
    - **de_genes_S0_S1 and S0_S2.pdf**: Detected differentially expressed top 15 genes plot. Red bars are genes that have higher gene expression in branch S0_S1, blue bars are genes that have higher gene expression in branch S0_S2
    - **de_genes_greater_S0_S1 and S0_S2.tsv**: Table that stores information of DE genes that have higher expression in branch S0_S1.
    - **de_genes_less_S0_S1 and S0_S2.tsv**: Table that stores information of DE genes that have higher expression in branch S0_S2.
*   sub-folder **'leaf_genes'** contains several files:
    - **leaf_genes.tsv**: Table that stores information of leaf genes from all branches.
    - **leaf_genesS0_S1.tsv**: Table that stores information of leaf genes from branch S0_S1.
*	sub-folder **'S0'**: contains subway and stream plots for each of the cell states, for example, choosing S0 state as root state:   
    - **subway_map.pdf**: single-cell level cellular branches plot
    - **stream_plot.pdf**: density level cellular branches plot
    - **subway_map_gene.pdf**: gene expression pattern on subway map plot
    - **stream_plot_gene.pdf**: gene expression pattern on stream plot


Credits: H Chen, L Albergante, JY Hsu, CA Lareau, GL Bosco, J Guan, S Zhou, AN Gorban, DE Bauer, MJ Aryee, DM Langenau, A Zinovyev, JD Buenrostro, GC Yuan, L Pinello
