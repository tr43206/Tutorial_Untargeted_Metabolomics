# Metabolomics_Tutorial
+ Output Excel file of Compound Discoverer (Thermo Fisher Scientific; v3.3 SP2) was used for analysis.

+ `Gut-Lung Axis` project
+ Orbitrap LC-MS/MS
+ Phase A (Stationary phase) : Water with 1% FA (Formic acid)
+ Phase B (Mobile phase) : Methanol with 1% FA
+ Total of 20 minutes protocol was used.


## Step 0 : Export Excel file from Compound Discoverer
1) Initiate Compound Discoverer.
2) `File` -> `Open Results...` -> Open `CDRESULTVIEW` or `CDRESULT` format file.
3) Check the `Area` box in `Field Chooser` (its below the `Compounds` tag, next to the `Tags` column).
4) Right-click anywhere in the metabolic table -> `Export` -> `As Excel...` -> change the `Path`, and click the `Export` button.
5) After exporting, open exported excel file, and change `Area: ~` format columns to desiring sample names.
   I recommend to name them with different starting alphabet, to avoid errors in further steps (e.g., If your experimental groups are `Sucrose`, and `Sucralose`, name them as `A01`, `B01`, not like `So01`, `Sa01`.)


## Step 1 : Scaling method selection
+ Use `01. scaling_selection_for_metabolomics_GutLungAxis.ipynb`.

+ Select between `none`, `auto`, `pareto`, `log`, `log_auto`, `log_pareto`.
+ The more pooled QC samples become dense in PCA plot, it means its less biased.
+ If the raw data shows right-skewed distribution, Log 10 transformation is recommended.


## Step 2 : Noise filtering and scaling
+ Use `02. Preprocessing_GutLungAxis.ipynb`.

+ Filter noises based on subtraction between intensities of samples and blanks (only blanks used in extraction procedure is used in this step).
+ If <`mean intensities of samples except pooled QCs` - `mean intensities of blanks`> < 0, it means corresponding metabolic feature could be interpreted as noise.
+ Adjust filtering method from `Step 1`.


## Step 3 : PCA visualization
+ Use `03. PCA_GutLungAxis.ipynb`.

+ Statistical significance could also be identified using `Hotelling's T-squared test` (using PC values which cumulative sum proportion excess 80%).


## Step 4 : Statistical analysis with FDR adjusting
+ Use `04. MeboAnalyst_OneFactor_format+bky_26.01.05.ipynb`.

+ Statistical analysis with Mann-Whitney test (if only two groups are compared) or Kruskal-Wallis test (if more than two groups are compared).
+ FDR (False Discovery Rate) should be adjusted if multiple testing was performed (BKY method (Benjamini, Krieger, and Yekutieli, 2006) was used in the python script).


## Step 5 : Mapping to KEGG pathways
+ Use `05. MetaboAnalyst_Pathways_format_rev_GutLungAxis.ipynb` to make `MetaboAnanlyst 6.0` (https://www.metaboanalyst.ca/) format input files.
+ After processing above python script, Greek alphabets should be modified (e.g., α -> alpha, β -> beta, γ -> gamma).

1) Open up the MetaboAnalyst (link is above).
2) `Pathway Analysis` -> `Concentration table` ->
   Group Label : `Discrete`
   ID Type : `Compound Name`
   Data Format : `Samples in columns`
   Data File : (Output file from above python script)
   -> `Proceed` if sample count is correct
3) Sample normalization : `None`
   Data transformation : `None`
   Data scaling : `None`
4) Select a pathway library : Mammals : `Mus musculus (house mouse) (KEGG)` (if your samples are from human, select `Homo sapiens (KEGG)`. Other options are set to default.)
   -> Briefly check the result, then check `Submit`.
5) Download `pathway_result.csv`, then change their name to distinguish with other output files (e.g., `pathway_results_2w.csv`).

+ Use `06. Pathway_Enrichment_GutLungAxis.R` for visualization to enrichment plot.
+ (Optional) If you want to merge paired sample data in one enrichment plot, use `06. Pathway_Enrichment_merged_GutLungAxis.R`.


## Step 6 : Volcano plot visualization

+ Use `07. Volcano_plot_GutLungAxis_26.04.26.R` to visualize differential metabolites (default value is `q-value < 0.05` & `|Log2FoldChange| < 1`).
+ You can also show only specific metabolic features in interested pathways.


## Step 7 : Candidate biomarkers discovery

+ Use `08. Individual_metabolites_merged_v2_GutLungAxis.ipynb` to individually visualize differential metabolites into box plots.
+ (Optional) Or you can also find candidate biomarkers by using `MFuzz` package in `R`, which clusters metabolic features based on trending tendencies. In this case, use `09. Mfuzz_GutLungAxis_26.04.26.R`.

+ Identify candidate features using venn-diagram or Mfuzz clustering or supervised methods-based multivariate analysis, such as PLS-DA, OPLS-DA (to minimize risks of over-fitting), and sPLS-DA. Machine learning methods can be also used, such as Random Forest (RF; recommended for large dataset) or Suppot Vector Machine (SVM; recommended for small dataset).

