# Metabolomics_Tutorial
+ Used Excel output file of Compound Discoverer (Thermo Fisher Scientific; v3.3 SP2), using Orbitrap LC-MS/MS data of `Gut-Lung Axis` project.



## Step 0 : Export Excel file of Compound Discoverer
1) `File` -> `Open Results...` -> Open `CDRESULTVIEW` or `CDRESULT` format file.
2) Check the `Area` box in `Field Chooser`. Its below the `Compounds` tag, on the left of `Tags` column.
3) Right-click the mouse on the metabolic table -> `Export` -> `As Excel...` -> change the `Path`, and click the `Export` button.
4) After exporting, change the `Area: ~` format columns into sample names.
   Separate different groups with same patterns is recommended (e.g., If your experimental groups are `Control`, `VNAM`, `Solvent`, and `MIA`, name them as `C01`, `V01`, `S01`, `M01`, .... Avoid using same starting alphabet.)


## Step 1 : Select appropriate scaling method
(Use `01. scaling_selection_for_metabolomics_GutLungAxis.ipynb`)

Select between `none`, `auto`, `pareto`, `log`, `log_auto`, `log_pareto`.
The more pooled QC samples get dense in PCA plot, it means scaling method is showing less bias.



## Step 2 : Preprocessing - Noise filtering and scaling
(Use `02. Preprocessing_GutLungAxis.ipynb`)

Filter noises based on subtraction of intensities between samples and blank.
If <`mean intensities of samples except pooled QC` - `mean intensities of blanks`> < 0, it means the metabolic feature could be interpreted as noise.
Adjust filtering method based on `Step 0`.



## Step 3 : PCA
(Use `03. PCA_GutLungAxis.ipynb`)

Observe the group separation using PCA (Principal Coordinate Analysis) method.
Statistical significance could also be performed using `Hotelling's T-squared test` (using PC values of cumulative sum proportion up to 80%).



## Step 4 : Multiple testing and p-values adjustment
(Use `04. MeboAnalyst_OneFactor_format+bky_26.01.05.ipynb`)

Multiple testing with Mann-Whitney test (if only two groups are compared) or Kruskal-Wallis test (if more than three groups are compared).
FDR (False Discovery Rate) adjustment could be performed using BKY method (Benjamini, Krieger, and Yekutieli, 2006).



## Step 5 : Pathway analysis
1) Use `05. MetaboAnalyst_Pathways_format_rev_GutLungAxis.ipynb` in pre-processing step for MetaboAnanlyst 6.0 formats.
2) Use `06. Pathway_Enrichment_GutLungAxis.R` for enrichment plot visualization after adjusting MetaboAnalyst 6.0 step.
3) (Optional) If you want to merge several groups in one enrichment plot, use `06. Pathway_Enrichment_merged_GutLungAxis.R`.

Mapping to KEGG pathway in MetaboAnalyst 6.0. This step is to exract significant features (FeatureID: `m/z`-`RT [min]`, without rounding) based on `Step 3`.



## Step 6 : Adding Log 2 FoldChange criteria

(Use `07. Volcano_plot_GutLungAxis_26.04.26.R` to visualize differential metabolites after applying Log 2 FoldChange criteria)

Find out differentially abundant features by adding `lfc` criteria. Basically, it is set to 1, which means only features showing more than two times differences will be remained.



## Step 7 : Finding candidate biomarkers

1) Use `08. Individual_metabolites_merged_v2_GutLungAxis.ipynb` to individually visualize differential metabolites to box plots.
2) Or you can find candidate biomarkers using MFuzz package in R, which clusters metabolic features based on trending tendency.
Use `09. Mfuzz_GutLungAxis_26.04.26.R`)

Identify common features using venn-diagram, or Mfuzz clustering, or supervised methods-based multivariate analysis, such as PLS-DA, OPLS-DA, or sPLS-DA. Machine learning methods can be also used, such as Random Forest (RF) or Suppot Vector Machine (SVM).
e.g., Mfuzz clustering : Find out differential features which shows same trend (same cluster).
