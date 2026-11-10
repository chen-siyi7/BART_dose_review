
argv<-commandArgs(trailingOnly=TRUE);stopifnot(length(argv)==1L,!file.exists(argv[1]))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE));root<-normalizePath(file.path(dirname(script),'..'))
dir.create(argv[1],recursive=TRUE);out<-normalizePath(argv[1])
for(f in c('reference.R','decisions.R','certificates.R','audits.R','feasibility.R','feasibility_audit.R'))source(file.path(root,'R',f))
result<-list(status='PASS',feasibility=audit_feasibility_study(file.path(root,'data'),out),
             weighted_report=audit_weighted_example(file.path(root,'data'),out))
save_json(result,file.path(out,'RESULTS.json'));cat('New R study and weighted report replay passed\n')
