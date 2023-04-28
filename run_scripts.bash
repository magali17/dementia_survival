# this script runs all important scripts for the survival analysis
# to run this script in the terminal, enter: bash run_scripts.bash

Rscript 1_prep_data.R

Rscript -e 'rmarkdown::render("2_descriptives.Rmd", "html_document")'
Rscript -e 'rmarkdown::render("incidence.Rmd", "html_document")'

Rscript 3_run_models.R
Rscript -e 'rmarkdown::render("4_model_results.Rmd", "html_document")'



echo "done running scripts"
