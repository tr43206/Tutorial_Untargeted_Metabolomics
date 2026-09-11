# Metabolomics_Tutorial
Output files of Compound Discoverer (v3.3 SP2), using Orbitrap LC-MS/MS data of Gut-Lung Axis project. 



## Step 0 : Select appropriate scaling method
(Use `01. scaling_selection_for_metabolomics_GutLungAxis.ipynb`)

Select between `none`, `auto`, `pareto`, `log`, `log_auto`, `log_pareto`.



## Step 1 : Preprocessing - Noise filtering and scaling
(Use `02. Preprocessing_GutLungAxis.ipynb`)

1) Filter noises based on subtraction of intensities between samples and blank.
2) Adjust filtering method based on `Step 0`.



## Step 2 : PCA
(Use `03. PCA_GutLungAxis.ipynb`)

Multivariate analysis using PCA method, and the statistical analysis is performed using Hotelling's T-squared test.



## Step 3 : Multiple testing and p-values adjustment
(Use `04. MeboAnalyst_OneFactor_format+bky_26.01.05.ipynb`)

Multiple testing with Mann-Whitney test (if only two groups are compared) or Kruskal-Wallis test (if more than three groups are compared).



## Step 4 : Pathway analysis
(1) Use `05. MetaboAnalyst_Pathways_format_rev_GutLungAxis.ipynb` in pre-processing step for MetaboAnanlyst 6.0 formats)
(2) Use `06. Pathway_Enrichment_GutLungAxis.R` for enrichment plot visualization after adjusting MetaboAnalyst 6.0 step)
(3) (Optional) If you want yo merge several groups in one enrichment plot, use `06. Pathway_Enrichment_merged_GutLungAxis.R`)

Mapping to KEGG pathway in MetaboAnalyst 6.0. This step is to exract significant features (FeatureID: `m/z`-`RT [min]`) based on `Step 3`.



## Step 5 : Adding Log 2 FoldChange criteria

(Use `07. Volcano_plot_GutLungAxis_26.04.26.R` to visualize differential metabolites after applying Log 2 FoldChange criteria)

Find out differentially abundant features by adding `lfc` criteria. Basically, it is set to 1, which means only features showing more than two times differences will be remained.



## Step 6 : Finding candidate biomarkers

(1) Use `08. Individual_metabolites_merged_v2_GutLungAxis.ipynb` to individually visualize differential metabolites to box plots)
(2) Or you can find candidate biomarkers using MFuzz package in R, which clusters metabolic features based on trending tendency.
Use `09. Mfuzz_GutLungAxis_26.04.26.R`)

Identify common features using venn-diagram, or Mfuzz clustering, or supervised methods-based multivariate analysis, such as PLS-DA, OPLS-DA, or sPLS-DA. Machine learning methods can be also used, such as Random Forest (RF) or Suppot Vector Machine (SVM).
e.g., Mfuzz clustering : Find out differential features which shows same trend (same cluster).
