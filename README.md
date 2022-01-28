# Highly charged protein regions have secondary structure and conserved sequence properties

## Layout

### **`code`**  
 * **`analysis`**: Jupyter notebooks used to analyze and visualize raw data.
 * **`processing`**: Jupyter notebooks used to filter or extract subsets of raw data.
 * **`scripts`**: Python scripts containing functions used to perform analysis in `analysis` and `processing`.
 * **`tests`**: All test suites for code in `scripts`.

### **`data`** 
 * **`af_regions`**: Protein regions extracted from the AlphaFold dataset.
 * **`charged_regions`**: List of charged regions. `cr_raw.csv` contains the raw output from `extract_charged_regions.py`. `cr_filtered*.csv` contains the list of charged regions with a verified orf and a valid AYbRAH MSA. `cr_filtered_aflabel.csv` has an appended column characterizing the structure of the region based on AlphaFold prediction. `cr_filtered_lrlabel.csv` has an appended column characterizing the structure of the region based on prediction by the logistic regression model. `cr_filtered_unlabeled.csv` contains the list of regions whose structure cannot be defined based on the AlphaFold prediction.
 * **`misc`**: All other data files used.
 * **`ref_proteomes`**: Reference proteomes from SGD and 1000 Fungal Genomes Project.
 * **`sc_orfs`**: List of Sacharomyces cerevisae ORFs with different attributes.
 * **`uversky`**: Output files from Uversky analysis.

### **`figure`** 
Jupyter notebooks used to produce final figures for the manuscript.

### **`misc`** 
Notes and color schemes.