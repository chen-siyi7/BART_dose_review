
argv<-commandArgs(trailingOnly=TRUE);stopifnot(length(argv)==1L,!file.exists(argv[1]))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE));root<-normalizePath(dirname(script))
out<-argv[1];parent<-normalizePath(dirname(out));stopifnot(parent!=root,!startsWith(parent,paste0(root,'/')))
dir.create(out,recursive=TRUE);out<-normalizePath(out)
for(n in c('reference.R','decisions.R','certificates.R','sensitivity.R','calibration.R','audits.R','feasibility.R','feasibility_audit.R','reviewer_comparison_audit.R')) source(file.path(root,'R',n))
data<-file.path(root,'data');results<-list()
if(file.exists(file.path(root,'MANIFEST.json'))) results$manifest<-verify_manifest(root)
run<-function(name,fun) {
  cat('Checking',name,'\n');flush.console();value<-fun();results[[name]]<<-value
  save_json(value,file.path(out,paste0(name,'.json')));cat(name,'PASS\n');flush.console()
}
run('forest_accuracy',function()audit_forest_counts(data,out))
run('sampler',function()audit_sampler())
run('intercept',function()audit_intercept(file.path(data,'extended/exact_reference_benchmark_20261103'),out))
run('frontiers',function()audit_frontiers(data,out))
run('certificates',function()audit_certificates(data,out))
run('ordinary',function()replay_ordinary(file.path(data,'base/query/source/ordinary_review_benchmark_20261102'),out))
run('separate_calibration',function()replay_comparison(file.path(data,'extended/comparative_reference_20261103'),out))
run('prior',function()audit_prior(data,out))
run('prediction_points',function()audit_prediction_points(data))
run('original_workload',function()audit_workload(file.path(data,'worked_review')))
if(file.exists(file.path(data,'R_workload/RESULTS.json'))) run('R_workload',function()audit_workload(file.path(data,'R_workload')))
if(file.exists(file.path(data,'R_feasibility/COMPLETE.json'))) run('feasibility',function()audit_feasibility_study(data,out))
run('weighted_report',function()audit_weighted_example(data,out))
if(file.exists(file.path(data,'R_reviewer_comparisons/RESULTS.json'))) run('reviewer_comparisons',function()audit_reviewer_comparisons(root,out))
if(file.exists(file.path(data,'R_astodrimer_application/results/RESULTS.json'))) run('astodrimer',function() {
  clinical_out<-file.path(out,'astodrimer_replay');log<-file.path(out,'astodrimer.log')
  stopifnot(system2(file.path(R.home('bin'),'Rscript'),
    c(shQuote(file.path(root,'scripts/replay_astodrimer.R')),shQuote(clinical_out)),stdout=log,stderr=log)==0L)
  result<-read_json(file.path(clinical_out,'RESULTS.json'));stopifnot(result$status=='PASS');result
})
if(file.exists(file.path(root,'MANIFEST.json'))) results$final_manifest<-verify_manifest(root)
save_json(list(status='PASS',language='R',checks=results,session=capture.output(sessionInfo()),new_fits=FALSE,new_timing=FALSE),file.path(out,'RESULTS.json'))
cat('All R replay checks passed\n')
