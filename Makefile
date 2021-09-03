knit-clean-data:
	Rscript -e 'rmarkdown::render("Analysis/data_cleaning.Rmd", output_dir = "Output", encoding = "UTF-8")'

knit:
	Rscript -e 'rmarkdown::render("Analysis/Analysis.Rmd", output_dir = "Output", encoding = "UTF-8")'

knit-desc:
	Rscript -e 'rmarkdown::render("Analysis/descriptive.Rmd", output_dir = "Output", encoding = "UTF-8")'