# Transfer docx data to txt
#
# @param
#
# @return
#
# @export

install.packages("docxtractr")

library(docxtractr)

importdocx <- function() {

  files <- list.files() #ディレクトリ内のファイル名をfilesに代入
  
  out <- matrix(0, ncol=2, nrow=length(files))


  for (i in 1:length(files)) {
    file.name <- files[i]

    if (regexpr('.docx$', file.name)  < 0) { # ファイル名の最後が '.txt'か？
      next                                 # そうでなければリストから削除してスキップ．
    }

    a <- read_docx(path=file.name)
    
    b <- docx_extract_tbl(docx = a, tbl_number = 1, header = FALSE)$V2[5]
    out[i,1] <- substr(b, 38, 41)
    
    c <- docx_extract_tbl(docx = a, tbl_number = 1, header = FALSE)$V2[3]
    out[i,2] <- substr(c, 37, 41)

  }
  return(out)
  }
  
  
result <- importdocx()

