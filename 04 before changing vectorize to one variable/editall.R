#
# edits all files (but doesn't work)
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/makederivatives/',sep=""))
#
# List all R script files
files <- list.files(path = "./", pattern = "\\makederivs.R$", full.names = TRUE)
#
# Loop through each file, edit, and save
for (file in files) {
  # Read the script
  script <- readLines(file)

  # Modify the script (example: replace "old_function" with "new_function")
  script <- gsub("\"x\"", "\\\"x\\\"", script)

  # Write the modified script back to the file
  writeLines(script, file)
  stop()
}
