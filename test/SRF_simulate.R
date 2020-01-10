#' 
#' 
#' 
#' 
#' 

file.choose()

fn <- "C:\\Users\\Shun Bi\\Documents\\MATLAB\\SRF.xlsx"

SRF_LIST <- list()
sheets <- excel_sheets(fn)
for(sheet in sheets){
  
  tmp <- read_excel(fn, sheet=sheet)
  SRF_LIST[[sheet]] <- tmp
  
  
  
}

