

interfile_folder <- "D:/datasets/yunfu/interfile_folder/Normative_Model"
reference_interfile_folder <- "D:/datasets/yunfu/doublecheck_xy/results"

task_name <- "1-back"

deviationZ <- readRDS(file.path(interfile_folder, task_name, "back1Acc.deviations.rds"))
reference_deviationZ <- readRDS(file.path(reference_interfile_folder, "Oneback_Oneback_accdeviation.df.csv"))


