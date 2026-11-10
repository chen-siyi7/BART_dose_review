args <- commandArgs(trailingOnly=TRUE)
stopifnot(!length(args))
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
root <- normalizePath(file.path(dirname(script),'..'))
required <- c('DESCRIPTION','NAMESPACE','README.md','R','man','inst','tests')
stopifnot(all(file.exists(file.path(root,required))))
files <- list.files(root,recursive=TRUE,all.files=TRUE,full.names=TRUE,no..=TRUE)
files <- files[!file.info(files)$isdir & !grepl('/[.]git/',files)]
r <- files[grepl('[.]R$',files,ignore.case=TRUE)]
for(p in r) {
  tokens <- getParseData(parse(file=p,keep.source=TRUE))
  stopifnot(!any(tokens$token=='COMMENT'))
}
stopifnot(!any(grepl('[.](pdf|png|jpg|zip|rds|rda|csv|xlsx|tex|py|sh|js)$',files,ignore.case=TRUE)))
description <- read.dcf(file.path(root,'DESCRIPTION'))
stopifnot(description[1L,'Package']=='bartDoseReview')
for(p in list.files(file.path(root,'man'),pattern='[.]Rd$',full.names=TRUE))
  tools::parse_Rd(p)
cat('R package layout and documentation verified;',length(r),'R files parse with no comments.\n')
