# LC-MS Automated Analysis Tool

This is a lightweight web application built with R Shiny and designed for laboratory researchers. By uploading a properly formatted CSV file, you can automatically remove outliers from LC-MS metabolite data, run statistical tests (`t-test` / `Wilcoxon test`), and generate publication-ready bar plots and scatter plots with a single click.

## Features

* **No-code workflow**: Use an intuitive web interface without modifying any R code.
* **Automated safeguards and outlier handling**: Built-in Z-score checking removes outliers automatically and safely handles samples that contain only `NA` values.
* **Dynamic statistical testing**: You can enter the name of the control group (`Standard`), and the system will automatically compare it with other genotypes and annotate significance levels with stars (`*`, `**`, `***`).
* **Flexible test methods**: Supports both Student's `t-test` and Mann-Whitney U test (`Wilcoxon`).
* **One-click bundled download**: After analysis, all metabolite plots can be downloaded together as a single `.zip` file.

## Data Format Requirements

The uploaded file must be in **`.csv`** format and include the following exact column names (case-sensitive):

1. `Genotype`: Genotype name (for example: `B73-W0`, `double-W2`)
2. `Treatment`: Treatment time or condition (for example: `KODA`, `KODA`)
3. `Rep`: Replicate sample ID (for example: `1`, `2`, `3`)
4. `Background`: Background information field
5. From the **6th column onward**: Quantitative values for each metabolite. These column names will automatically be used as plot titles.

> **Tip**: The system automatically sorts the X-axis based on the numbers in `Treatment` (such as `W0`, `W2`) and the control group you define (such as `B73`).

## How to Use

### Option 1: Online Use (Web App)
Use the following link to run the analysis directly:
  
👉 [Open the LC-MS Analysis Tool](https://ericbjchen.shinyapps.io/lcms_shiny/)

### Option 2: Run Locally
If you have an R environment set up, you can run the app locally on your computer:

1. Make sure the following R packages are installed:
   ```r
   install.packages(c("shiny", "data.table", "dplyr", "stringr", "ggplot2", "cowplot"))
   ```
2. Clone this repository or download the source code.
3. Open R or RStudio and set the working directory to the project folder.
4. Run the following command:
   ```r
   shiny::runApp()
   ```

## Output

After the analysis is complete, the application will generate:

* Publication-ready metabolite plots
* Statistical significance annotations
* A downloadable `.zip` archive containing all generated figures

## Notes

* Ensure the uploaded CSV file follows the required column naming and format exactly.
* If the control group name is entered incorrectly, statistical comparisons and plot ordering may not behave as expected.
