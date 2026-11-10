


for(p in c('readxl','jsonlite','digest'))stopifnot(requireNamespace(p,quietly=TRUE))
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
verification <- normalizePath(dirname(script));argv <- commandArgs(trailingOnly=TRUE)
stopifnot(length(argv)<=2L)
out <- if(length(argv)==2L)argv[2]else verification
dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out)
raw <- if(length(argv)) normalizePath(argv[1]) else normalizePath(file.path(verification,'..','raw'))
sha <- function(p)digest::digest(file=p,algo='sha256')
expected <- c('pone.0232394.s002.xlsx'='11f77074b9008bc5fa67262aedfb128fb92dd2a1ed6adc18f17e97f0c9f26c4a',
 'pone.0232394.s003.xlsx'='6ae7e90a613ddaabb1859509b70cc9062588d37cb8efd05876a0b035ce5842f6',
 'pone.0232394.s005.xlsx'='eda1f419014529668fe3526131a6867b2cc78f540a95b3010ae366b5d72e5370')
for(f in names(expected))stopifnot(sha(file.path(raw,f))==expected[[f]])
p <- as.data.frame(readxl::read_excel(file.path(raw,names(expected)[1]),sheet='Populations'))
e <- as.data.frame(readxl::read_excel(file.path(raw,names(expected)[2]),sheet='Efficacy'))
a <- as.data.frame(readxl::read_excel(file.path(raw,names(expected)[3]),sheet='AEs'))

key <- function(x)vapply(strsplit(trimws(x),'/',fixed=TRUE),paste,collapse='',FUN.VALUE='')
pk <- key(p[['Site/ Subject']]);ek <- key(e[['Site/Subject']]);ak <- trimws(a$Patient)
stopifnot(!anyDuplicated(pk),!anyDuplicated(ek),all(ak%in%pk),all(ek%in%pk))
match_e <- match(pk,ek)
treatment <- p[['Treatment Group']]
safety <- p$Safety=='Yes';mitt <- p$mITT=='Yes';pp <- p$PP=='Yes'
stopifnot(sum(safety)==125L,sum(mitt)==111L,sum(pp)==107L,
 all(!is.na(match_e[safety])),all(is.na(match_e[!safety])),
 all(treatment[safety]==e$Treatment[match_e[safety]]),
 all(p$mITT[safety]==e$mITT[match_e[safety]]))
decode <- function(x) {
 stopifnot(all(is.na(x)|x%in%c('Yes','No')))
 ifelse(is.na(x),NA_integer_,as.integer(x=='Yes'))
}
toc <- decode(e[['CC TOC']])[match_e]
eot <- decode(e[['CC EOT']])[match_e]
rel <- c('POSSIBLY RELATED','PROBABLY RELATED','DEFINITELY RELATED')
related <- unique(ak[a$Relationship%in%rel])
any_ae <- unique(ak)
stopifnot(all(pk[pk%in%any_ae]%in%pk[safety]),length(any_ae)==84L,length(related)==30L,
 sum(is.na(toc[mitt]))==4L,sum(!is.na(toc[mitt]))==107L,
 sum(is.na(toc[safety]))==5L,sum(!is.na(toc[safety]))==120L)
miss <- which(mitt & is.na(toc))
stopifnot(all(!is.na(eot[miss])),identical(eot[miss],c(1L,1L,1L,0L)))
stopifnot(identical(pk[miss],c('10026','10029','14009','14038')),
 identical(match_e[miss]+1L,c(19L,22L,108L,116L)))
arms <- c('HEC placebo','0.5% SPL7013','1% SPL7013','3% SPL7013')
summarize <- function(mask)do.call(rbind,lapply(arms,function(arm) {
 z <- mask & treatment==arm;v <- toc[z]
 data.frame(treatment=arm,n=sum(z),cure_yes=sum(v==1L,na.rm=TRUE),
  cure_no=sum(v==0L,na.rm=TRUE),cure_missing=sum(is.na(v)),
  related_AE=sum(pk[z]%in%related),any_AE=sum(pk[z]%in%any_ae),stringsAsFactors=FALSE)
}))
safety_summary <- summarize(safety);mitt_summary <- summarize(mitt)
stopifnot(identical(safety_summary$related_AE,c(7L,8L,6L,9L)),
 identical(safety_summary$any_AE,c(23L,22L,19L,20L)),
 identical(mitt_summary$related_AE,c(5L,8L,4L,7L)),
 identical(mitt_summary$cure_yes,c(3L,7L,12L,7L)),
 identical(mitt_summary$cure_missing,c(2L,0L,1L,1L)))
locf <- toc;locf[miss] <- eot[miss]
locf_successes <- vapply(arms,function(arm)sum(locf[mitt & treatment==arm]),0)
stopifnot(identical(unname(locf_successes),c(5,7,12,8)))

patterns <- expand.grid(rep(list(0:1),length(miss)),KEEP.OUT.ATTRS=FALSE)
completions <- do.call(rbind,lapply(seq_len(nrow(patterns)),function(i) {
 z <- toc;z[miss] <- as.integer(patterns[i,])
 row <- data.frame(pattern=paste0(as.integer(patterns[i,]),collapse=''),
  stringsAsFactors=FALSE)
 for(j in seq_along(arms))row[[c('placebo','dose_0_5','dose_1','dose_3')[j]]] <- sum(z[mitt & treatment==arms[j]])
 row
}))
stopifnot(nrow(completions)==16L,
 all(completions$dose_0_5==7L),all(completions$dose_1%in%12:13),
 all(completions$dose_3%in%7:8),all(completions$placebo%in%3:5))
protocol_file <- file.path(verification,'pone.0232394.s006_protocol.pdf')
protocol_pages <- if(file.exists(file.path(verification,'protocol_text.txt')))
 strsplit(paste(readLines(file.path(verification,'protocol_text.txt'),warn=FALSE),collapse='\n'),'\f',fixed=TRUE)[[1L]] else NULL
if(!is.null(protocol_pages))stopifnot(
 any(grepl('last observation carried forward (LOCF)',protocol_pages,fixed=TRUE)),
 any(grepl('Assuming Clinical Cure rates of 50% and 15%',protocol_pages,fixed=TRUE)),
 any(grepl('abnormally high number',protocol_pages,fixed=TRUE)))
report <- list(status='PASS',date='2026-11-09',
 source_hashes=as.list(expected),
 original_workbooks=list(population_rows=nrow(p),efficacy_rows=nrow(e),AE_event_rows=nrow(a)),
 all_identifiers_match=TRUE,unmatched_AE_or_efficacy_identifiers=0L,
 public_primary_efficacy=list(population='mITT with as-observed primary clinical cure',
  eligible=111L,observed=107L,missing=4L,
  definition='No abnormal discharge and at most one of the other three Amsel findings at study days21-30',
  outcome_field='CC TOC',missing_is_not_coded_as_No=TRUE),
 all_treated_efficacy=list(population='Safety/treated population',eligible=125L,observed=120L,missing=5L,
  role='Different analysis population; not the publication primary mITT population'),
 original_safety=list(population=125L,not_treated=7L,related_levels=rel,
  endpoint='One participant indicator for at least one potentially treatment-related AE, at any recorded severity',
  not_a_DLT_endpoint=TRUE,no_recorded_AE_for_safety_participant_means_zero=TRUE),
 safety_population=safety_summary,mITT_population=mitt_summary,
 per_protocol_warning=list(original_PP_flag_n=107L,observed_primary_mITT_n=107L,
  overlap=104L,PP_only=3L,observed_primary_mITT_only=3L,
  same_number_does_not_mean_same_participants=TRUE),
 missing_primary_cells=data.frame(normalized_source_ID=pk[miss],published_efficacy_excel_row=match_e[miss]+1L,
  treatment=treatment[miss],age=as.numeric(sub('/.*','',p[['Age/ Race']][miss])),
  CC_EOT=eot[miss],original_PP_flag=pp[miss]),
 protocol=list(doi='https://doi.org/10.1371/journal.pone.0232394.s006',
  direct_url='https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0232394.s006&type=supplementary',
  sha256=if(file.exists(protocol_file))sha(protocol_file)else NULL,
  version='Final V3.0,15July2010',
  primary_definition=list(section='10.2.1 and10.2.4.1',printed_pages=c(40L,41L),physical_PDF_pages=c(51L,52L)),
  missing_data=list(section='10.8 and10.9',printed_pages=45L,physical_PDF_page=56L,
   plan='LOCF as appropriate, plus as-observed sensitivity; mITT and per-protocol analyses'),
  safety_review=list(section='8.10 and10.7',printed_pages=c(36L,37L,45L),physical_PDF_pages=c(47L,48L,56L),
   criteria='Clinical blinded review when potentially related events appear abnormally frequent or severe; no fixed numerical AE-rate/Bayesian posterior threshold located'),
  dose_selection=list(section='10.4 and10.5',printed_pages=c(43L,44L),physical_PDF_pages=c(54L,55L),
   note='50%active and15%placebo cure assumptions justify sample size, not minimum efficacy criteria; the interim review guides doses for a later protocol without an absolute efficacy-safety selection formula')),
 sensitivity=list(protocol_linked_LOCF=list(population=111L,completion_pattern='1110',
   ordered_IDs=pk[miss],successes=as.list(locf_successes),
   imputation='Use the observed CC EOT for the four CC TOC blanks; this is an explicit sensitivity, not an observed outcome'),
  full_binary_completions=list(scope='All16 hypothetical binary assignments for the four mITT CC TOC blanks',
   bit_order='Published efficacy workbook row order as shown in missing_primary_cells',aggregate_successes=completions,
   recommendation='Keep model, priors, common prediction population and decision rule fixed; compare possible-dose sets across all completions',
   coverage='A simultaneous finite-sample numerical guarantee over all completion analyses requires error allocation over that family',
   application_plan_interpretation='The current17-scenario application uses fitting chains: its57-quantity error allocation is nominal/asymptotic under the stated chain central-limit and variance-consistency assumptions, not a finite-sample guarantee',
   interpretation='This addresses sensitivity to the four unobserved efficacy outcomes. It does not establish a missing-at-random mechanism, impute verified patient outcomes or give a population/model robustness bound')),
 reference_verification=list(DOD_BART=list(authors=c('Yunqi Zhao','Rachael Liu','Jianchang Lin','Andy Chi','Simon Davies'),
  doi='10.1080/10543406.2024.2429463',primary_metadata='https://pubmed.ncbi.nlm.nih.gov/39611660',current_main_author_list_correct=TRUE),
  Zhou_Ji=list(authors=c('Tianjian Zhou','Yuan Ji'),year=2021L,volume=40L,issue=28L,pages='6421-6442',
   doi='10.1002/sim.9191',primary_metadata='https://onlinelibrary.wiley.com/doi/10.1002/sim.9191')))
jsonlite::write_json(report,file.path(out,'RAW_ENDPOINT_AUDIT.json'),pretty=TRUE,auto_unbox=TRUE,digits=17,na='null',null='null')
specification <- list(status='FIXED_BEFORE_NEW_MODEL_FITS',scope='Independent endpoint and missing-outcome specification; full model/run settings belong in the main application plan',
 raw_source_hashes=as.list(expected),missing_primary_IDs=pk[miss],
 missing_source_order='Ascending data-row order in published efficacy workbook',
 missing_excel_rows=match_e[miss]+1L,LOCF_values=eot[miss],LOCF_pattern='1110',
 observed_primary_mITT=list(eligible=111L,observed=107L),
 scenarios=list(primary='As observed mITT107',sensitivity='All16binary completions of the four mITT TOC blanks',total_efficacy_fits=17L),
 safety=list(participants=125L,endpoint='Any POSSIBLY,PROBABLY,DEFINITELY related recorded AE per treated participant'),
 prediction_population='Common125treated participant ages',
 dose_candidates=c(.5,1,3),
 agreed_fitting_plan=list(covariates=c('concentration','placebo indicator','age'),offset=0,
  chains=4L,burnin_per_chain=8000L,retained_per_chain=32000L),
 nominal_numerical_family=list(total_error=.05,toxicity_means=3L,toxicity_screening_probabilities=3L,
  efficacy_means=17L*3L,total_quantities=57L,scalar_error=.05/57,
  interpretation='Nominal/asymptotic MonteCarlo intervals under fitting-chain assumptions; completion sensitivity is separate from numerical interval coverage'),
 original_trial_thresholds='No numerical efficacy/safety acceptability or Bayesian selection thresholds identified in the published protocol; any new threshold must be identified as a retrospective analysis setting')
sp <- file.path(out,'MISSING_OUTCOME_SPECIFICATION.json')
if(file.exists(sp))stopifnot(identical(jsonlite::fromJSON(sp,simplifyVector=FALSE),
 jsonlite::fromJSON(jsonlite::toJSON(specification,auto_unbox=TRUE,digits=17),simplifyVector=FALSE)))else
 jsonlite::write_json(specification,sp,pretty=TRUE,auto_unbox=TRUE,digits=17)
print(mitt_summary,row.names=FALSE)
cat('PASS: source hashes,participant joins,population denominators,AE coding and missing-data sensitivity verified.\n')
