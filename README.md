# Metabolomics_Tutorial

Output files of Compound Discoverer (v3.3 SP2), using Orbitrap LC-MS/MS data of Gut-Lung Axis project. 



# Step 0 : Select appropriate scaling method

Select between `none`, `auto`, `pareto`, `log`, `log_auto`, `log_pareto`.



# Step 1 : Preprocessing - Noise filtering and scaling

1) Filter noises based on subtraction of intensities between samples and blank.
2) Adjust filtering method based on `Step 0`.



# Step 2 : PCA

Multivariate analysis using PCA method, and the statistical analysis is performed using Hotelling's T-squared test.



# Step 3 : Multiple testing and p-values adjustment

Multiple testing with Mann-Whitney test (if only two groups are compared) or Kruskal-Wallis test (if more than three groups are compared).



# Step 4 : Pathway analysis

Mapping to KEGG pathway in MetaboAnalyst 6.0. This step is to exract significant features (FeatureID: `m/z`-`RT [min]`) based on `Step 3`.



# Step 5 : Adding Log 2 FoldChange criteria

Find out differentially abundant features by adding `lfc` criteria. Basically, it is set to 1, which means only features showing more than two times differences will be remained.



# Step 6 : Finding candidate biomarkers

Identify common features using venn-diagram, or Mfuzz clustering, or supervised methods-based multivariate analysis, such as PLS-DA, OPLS-DA, or sPLS-DA. Machine learning methods can be also used, such as Random Forest (RF) or Suppot Vector Machine (SVM).
e.g., Mfuzz clustering : Find out differential features which shows same trend (same cluster).
