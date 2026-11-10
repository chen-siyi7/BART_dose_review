

args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L,!file.exists(args[1]))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
root<-normalizePath(file.path(dirname(script),'..'))
out<-args[1];dir.create(out,recursive=TRUE);out<-normalizePath(out)
stopifnot(out!=root,!startsWith(out,paste0(root,'/')))
app<-file.path(root,'data/R_astodrimer_application')
stopifnot(dir.exists(app))
Rscript<-file.path(R.home('bin'),'Rscript')
rawout<-file.path(out,'raw_endpoints');analysisout<-file.path(out,'saved_chains')
rawlog<-file.path(out,'raw_endpoints.log');analysislog<-file.path(out,'saved_chains.log')
stopifnot(system2(Rscript,c(shQuote(file.path(app,'verification/audit_raw_endpoints.R')),
  shQuote(file.path(app,'raw')),shQuote(rawout)),stdout=rawlog,stderr=rawlog)==0L)
raw<-jsonlite::read_json(file.path(rawout,'RAW_ENDPOINT_AUDIT.json'),simplifyVector=TRUE)
stopifnot(raw$status=='PASS')
stopifnot(system2(Rscript,c(shQuote(file.path(app,'verification/audit_analysis.R')),
  shQuote(app),shQuote(analysisout)),stdout=analysislog,stderr=analysislog)==0L)
analysis<-jsonlite::read_json(file.path(analysisout,'ANALYSIS_AUDIT.json'),simplifyVector=TRUE)
stopifnot(analysis$status=='PASS',analysis$models==18L,analysis$efficacy_scenarios==17L)
contextout<-file.path(out,'posterior_context');contextlog<-file.path(out,'posterior_context.log')
stopifnot(system2(Rscript,c(shQuote(file.path(app,'descriptive_context/posthoc_context.R')),
  shQuote(app),shQuote(contextout)),stdout=contextlog,stderr=contextlog)==0L)
for(f in c('POSTERIOR_CONTEXT.json','CURE_CREDIBLE_INTERVALS.csv','PAIRED_DIFFERENCES.csv'))
  stopifnot(identical(digest::digest(file=file.path(contextout,f),algo='sha256'),
    digest::digest(file=file.path(app,'descriptive_context',f),algo='sha256')))
jsonlite::write_json(list(status='PASS',original_workbooks_audited=TRUE,original_protocol_supplied=TRUE,
  packed_application_audited=TRUE,models=18L,efficacy_scenarios=17L,
  posterior_context_audited=TRUE,posterior_context_separate_from_MC_family=TRUE,
  nominal_quantity_family=57L,nominal_family_error=.05,scalar_error=.05/57,
  raw_endpoint_audit=raw,saved_chain_audit=analysis,
  new_fits=FALSE,new_sampling=FALSE,new_timing=FALSE),file.path(out,'RESULTS.json'),
  auto_unbox=TRUE,pretty=TRUE,digits=17,na='null',null='null')
cat('Astodrimer original-data and saved-chain audits PASS; no model fitting\n')
