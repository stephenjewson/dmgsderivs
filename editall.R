#
# make derivative codes for fitdistcp, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
# List all R script files
files <- list.files(path = "./", pattern = "\\makederivs.R$", full.names = TRUE)

# Loop through each file, edit, and save
for (file in files) {
  # Read the script
  script <- readLines(file)

  # Modify the script (example: replace "old_function" with "new_function")
  script <- gsub("makederivatives", "dmgsderivs", script)

  # Write the modified script back to the file
  writeLines(script, file)
}
