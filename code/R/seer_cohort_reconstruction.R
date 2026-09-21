suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(splines)
})

clean_chr <- function(x) {
  z <- trimws(as.character(x))
  z[z %in% c("", "Blank(s)", "Not applicable", "NA")] <- NA_character_
  z
}

find_col <- function(d, pattern) {
  hit <- grep(pattern, names(d), ignore.case = TRUE, value = TRUE)
  if (!length(hit)) stop("missing column pattern: ", pattern)
  hit[1]
}

norm_combined <- function(x, prefix) {
  z <- toupper(clean_chr(x))
  z <- sub("^[CP]", "", z)
  ifelse(is.na(z), NA_character_, paste0(prefix, z))
}

first_int <- function(x) {
  z <- clean_chr(x)
  out <- suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", z)))
  out[is.na(z) | !grepl("[0-9]", z)] <- NA_real_
  out
}

age_mid <- function(x) {
  z <- clean_chr(x)
  lo <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", z)))
  hi <- suppressWarnings(as.numeric(sub("^[0-9]+-([0-9]+).*$", "\\1", z)))
  hi[!grepl("^[0-9]+-[0-9]+", z)] <- NA_real_
  out <- ifelse(!is.na(lo) & !is.na(hi), (lo + hi) / 2, lo)
  out[grepl("90\\+", z)] <- 92
  out[grepl("<1", z)] <- 0.5
  out
}

grade_group <- function(old, new, year) {
  a <- clean_chr(old); b <- clean_chr(new)
  out <- rep("Unknown", length(year))
  old_low <- grepl("^(Well differentiated|Moderately differentiated)", a, ignore.case = TRUE)
  old_high <- grepl("^(Poorly differentiated|Undifferentiated)", a, ignore.case = TRUE)
  new_low <- grepl("category \\([12]\\)|Well differentiated|Moderately differentiated|^Low grade$", b, ignore.case = TRUE)
  new_high <- grepl("category \\([34]\\)|Poorly differentiated|Undifferentiated|^High grade$", b, ignore.case = TRUE)
  out[year <= 2017 & old_low] <- "Low (I-II)"
  out[year <= 2017 & old_high] <- "High (III-IV)"
  out[year >= 2018 & new_low] <- "Low (I-II)"
  out[year >= 2018 & new_high] <- "High (III-IV)"
  factor(out, levels = c("Low (I-II)", "High (III-IV)", "Unknown"))
}

as_factor_unknown <- function(x, levels = NULL) {
  z <- clean_chr(x); z[is.na(z)] <- "Unknown"
  if (is.null(levels)) factor(z) else factor(z, levels = levels)
}

weighted_var <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  x <- x[ok]; w <- w[ok]
  m <- sum(w*x)/sum(w)
  sum(w*(x-m)^2)/sum(w)
}

smd <- function(x, g, w) {
  i0 <- g == 0 & is.finite(x) & is.finite(w)
  i1 <- g == 1 & is.finite(x) & is.finite(w)
  m0 <- sum(w[i0]*x[i0])/sum(w[i0]); m1 <- sum(w[i1]*x[i1])/sum(w[i1])
  abs(m1-m0)/sqrt((weighted_var(x[i0],w[i0])+weighted_var(x[i1],w[i1]))/2)
}

reconstruct <- function(path, cancer) {
  cat("READ", cancer, path, "\n")
  d <- fread(path, na.strings = c("", "NA"), showProgress = FALSE)
  raw_n <- nrow(d)
  yr <- find_col(d, "^Year of diagnosis$")
  t6 <- find_col(d, "Derived AJCC T, 6th")
  n6 <- find_col(d, "Derived AJCC N, 6th")
  m6 <- find_col(d, "Derived AJCC M, 6th")
  tc <- find_col(d, "Derived SEER Combined T")
  nc <- find_col(d, "Derived SEER Combined N")
  mc <- find_col(d, "Derived SEER Combined M")
  te <- find_col(d, "Derived EOD 2018 T Recode")
  ne <- find_col(d, "Derived EOD 2018 N Recode")
  me <- find_col(d, "Derived EOD 2018 M Recode")
  surg <- find_col(d, "RX Summ--Surg Prim Site \\(1998-2022")
  d[, year := suppressWarnings(as.integer(clean_chr(get(yr))))]
  d[, T := fcase(year <= 2015, toupper(clean_chr(get(t6))),
                 year <= 2017, norm_combined(get(tc), "T"),
                 year >= 2018, toupper(clean_chr(get(te))), default = NA_character_)]
  d[, N := fcase(year <= 2015, toupper(clean_chr(get(n6))),
                 year <= 2017, norm_combined(get(nc), "N"),
                 year >= 2018, toupper(clean_chr(get(ne))), default = NA_character_)]
  d[, M := fcase(year <= 2015, toupper(clean_chr(get(m6))),
                 year <= 2017, norm_combined(get(mc), "M"),
                 year >= 2018, toupper(clean_chr(get(me))), default = NA_character_)]
  d <- d[grepl("^T1([A-Z0-9]*)$", T)]
  t1_n <- nrow(d)
  d[, surg_code := suppressWarnings(as.integer(sub("^[^0-9]*([0-9]+).*$", "\\1", clean_chr(get(surg)))))]
  d[is.na(clean_chr(get(surg))) | !grepl("[0-9]", clean_chr(get(surg))), surg_code := NA_integer_]

  seqc <- find_col(d, "Sequence number")
  histc <- find_col(d, "Histologic Type ICD-O-3")
  diagc <- find_col(d, "^Diagnostic Confirmation$")
  chemoc <- find_col(d, "Chemotherapy recode")
  radc <- find_col(d, "Radiation recode")
  agec <- find_col(d, "Age recode with <1 year olds and 90")
  survdc <- find_col(d, "^Survival Days$")
  survfc <- find_col(d, "Survival months flag")
  cssc <- find_col(d, "SEER cause-specific death classification")
  othc <- find_col(d, "SEER other cause of death classification")
  vitalc <- find_col(d, "Vital status recode")
  sexc <- find_col(d, "^Sex$")
  racec <- find_col(d, "Race recode")
  maritc <- find_col(d, "Marital status at diagnosis")
  g17c <- find_col(d, "Grade Recode \\(thru 2017")
  g18c <- find_col(d, "Derived Summary Grade 2018")
  cszc <- find_col(d, "CS tumor size")
  sszc <- find_col(d, "Tumor Size Summary")
  d[, first_primary := grepl("One primary only|1st of", get(seqc), ignore.case=TRUE)]
  d[, hist_code := suppressWarnings(as.integer(substr(gsub("[^0-9]", "", as.character(get(histc))), 1, 4)))]
  if (cancer == "CRC") {
    hist_ok <- c(8140L,8141L,8143L,8144L,8145L,8147L,8210L,8211L,8220L,8221L,8261L,8262L,8263L,8480L,8481L,8490L)
    sitec <- find_col(d, "Site recode ICD-O-3")
    d[, site_label := clean_chr(get(sitec))]
    d[, site_ok := !is.na(site_label) & !site_label %in% c("Appendix", "Large Intestine, NOS")]
    d[, surgery_group := fcase(surg_code >=20 & surg_code <=29, "Local procedure",
                               surg_code >=30 & surg_code <=80, "Oncologic resection",
                               default="Excluded")]
  } else {
    hist_ok <- c(8120L,8130L,8082L,8122L,8131L)
    sitec <- find_col(d, "^Primary Site$")
    d[, primary_site := suppressWarnings(as.integer(clean_chr(get(sitec))))]
    d[, site_ok := !is.na(primary_site) & primary_site >=670 & primary_site <=679]
    d[, surgery_group := fcase(surg_code >=10 & surg_code <=27, "Local procedure",
                               surg_code %in% c(30,50,60:64,70:74,80), "Cystectomy",
                               default="Excluded")]
  }
  d[, hist_ok := hist_code %in% hist_ok]
  d[, hist_confirm := clean_chr(get(diagc)) == "Positive histology"]
  d[, M0 := toupper(clean_chr(M)) == "M0"]
  d[, N0NX := toupper(clean_chr(N)) %in% c("N0","NX")]
  d[, chemo0 := clean_chr(get(chemoc)) == "No/Unknown"]
  d[, rad0 := clean_chr(get(radc)) == "None/Unknown"]
  d[, age := age_mid(get(agec))]
  d[, surv_days := suppressWarnings(as.numeric(clean_chr(get(survdc))))]
  d[, valid_surv := !is.na(surv_days) & !grepl("Death Certificate Only|Autopsy Only", get(survfc), ignore.case=TRUE)]
  steps <- list(c("T1", nrow(d)))
  apply_step <- function(label, expr) {
    keep <- eval(expr, envir = d, enclos = parent.frame())
    d <<- d[which(keep)]
    steps[[length(steps)+1]] <<- c(label,nrow(d))
  }
  apply_step("first", quote(first_primary))
  apply_step("hist", quote(hist_ok))
  apply_step("hist_confirm", quote(hist_confirm))
  apply_step("site", quote(site_ok))
  apply_step("M0", quote(M0))
  apply_step("N0NX", quote(N0NX))
  apply_step("no_therapy", quote(chemo0 & rad0))
  apply_step("adult_survival", quote(age >=18 & valid_surv))
  apply_step("surgery", quote(surgery_group != "Excluded"))
  d[, exposure01 := as.integer(surgery_group != "Local procedure")]
  d[, exposure := factor(surgery_group, levels=unique(c("Local procedure",if(cancer=="CRC") "Oncologic resection" else "Cystectomy")))]

  if (cancer == "CRC") {
    smap <- c("Cecum"="Proximal colon","Ascending Colon"="Proximal colon","Hepatic Flexure"="Proximal colon","Transverse Colon"="Proximal colon","Splenic Flexure"="Distal colon","Descending Colon"="Distal colon","Sigmoid Colon"="Distal colon","Rectosigmoid Junction"="Rectum/rectosigmoid","Rectum"="Rectum/rectosigmoid")
    d[, sitegrp := factor(unname(smap[site_label]), levels=c("Proximal colon","Distal colon","Rectum/rectosigmoid"))]
  } else {
    d[, sitegrp := factor(fcase(primary_site==670,"Trigone",primary_site %in% 671:674,"Main wall",default="Special/NOS"),levels=c("Trigone","Main wall","Special/NOS"))]
  }
  d[, sex := as_factor_unknown(get(sexc),c("Female","Male","Unknown"))]
  d[, race := as_factor_unknown(get(racec))]
  mr <- clean_chr(d[[maritc]]); mr[is.na(mr)] <- "Unknown"
  d[, marital := factor(fcase(grepl("Single|never married",mr,ignore.case=TRUE),"Never married",grepl("^Married",mr,ignore.case=TRUE),"Married",grepl("Divorced|Separated",mr,ignore.case=TRUE),"Divorced/separated",grepl("Widowed",mr,ignore.case=TRUE),"Widowed",default="Unknown"),levels=c("Married","Never married","Divorced/separated","Widowed","Unknown"))]
  cs <- first_int(d[[cszc]]); ss <- first_int(d[[sszc]])
  d[, size := ifelse(year <=2015,cs,ss)]; d[size<=0 | size>500,size:=NA_real_]
  d[, sizegrp := factor(fcase(is.na(size),"Unknown",size<=20,"<=20 mm",size<=50,"21-50 mm",default=">50 mm"),levels=c("<=20 mm","21-50 mm",">50 mm","Unknown"))]
  d[, gradegrp := grade_group(get(g17c),get(g18c),year)]
  if(cancer=="CRC") d[, histgrp:=factor(fcase(hist_code==8490,"Signet-ring cell",hist_code%in%c(8480,8481),"Mucinous",default="Conventional/other adenocarcinoma"))]
  else d[, histgrp:=factor(fcase(hist_code==8120,"UC NOS",hist_code==8130,"Papillary UC",default="Other urothelial"))]
  d[, os_event := as.integer(grepl("Dead",get(vitalc),ignore.case=TRUE))]
  d[, css_event := as.integer(clean_chr(get(cssc))=="Dead (attributable to this cancer dx)")]
  d[, other_event := as.integer(clean_chr(get(othc))=="Dead (attributable to causes other than this cancer dx)")]
  d[surv_days<=0,surv_days:=0.5]
  fit <- glm(exposure01 ~ ns(age,df=4)+sex+race+marital+sitegrp+sizegrp+gradegrp+histgrp+ns(year,df=4),data=d,family=binomial())
  d[, ps:=pmin(pmax(predict(fit,type="response"),1e-6),1-1e-6)]
  d[, ow:=ifelse(exposure01==1,1-ps,ps)]
  endpoint <- function(ev,name){
    e <- as.integer(d$surv_days<=5*365.25 & ev==1)
    f <- coxph(Surv(pmin(d$surv_days/365.25,5),e)~d$exposure01,weights=d$ow,robust=TRUE,ties="efron")
    b<-coef(f)[1]; se<-sqrt(vcov(f)[1,1]); ph<-tryCatch(cox.zph(f,transform="km")$table[1,"p"],error=function(e)NA_real_)
    data.table(endpoint=name,events=sum(e),hr=exp(b),lower=exp(b-1.96*se),upper=exp(b+1.96*se),ph=ph)
  }
  outcomes <- rbind(endpoint(d$css_event,"cancer"),endpoint(d$other_event,"other"))
  z<-d[,.(pid=.I,exposure01,ow,t5=pmin(surv_days/365.25,5),css=as.integer(surv_days<=5*365.25&css_event==1),oth=as.integer(surv_days<=5*365.25&other_event==1))]
  st<-rbind(z[,.(pid,exposure01,ow,t5,endpoint="c",event=css)],z[,.(pid,exposure01,ow,t5,endpoint="o",event=oth)])
  st[,endpoint:=factor(endpoint,levels=c("c","o"))];st[,xc:=exposure01*(endpoint=="c")];st[,xo:=exposure01*(endpoint=="o")]
  jf<-coxph(Surv(t5,event)~xc+xo+strata(endpoint)+cluster(pid),data=st,weights=ow,robust=TRUE,ties="efron")
  b<-coef(jf)[c("xc","xo")];v<-vcov(jf)[c("xc","xo"),c("xc","xo")];lr<-b[1]-b[2];ser<-sqrt(c(1,-1)%*%v%*%c(1,-1))
  ratio <- c(est=exp(lr),lo=exp(lr-1.96*ser),hi=exp(lr+1.96*ser),cov=v[1,2])
  ess <- d[,.(n=.N,sumw=sum(ow),ess=sum(ow)^2/sum(ow^2),psmin=min(ps),psq01=quantile(ps,.01),psmed=median(ps),psq99=quantile(ps,.99),psmax=max(ps)),by=exposure]
  unk <- d[,.(n=.N,NX=sum(N=="NX"),unknown_size=sum(is.na(size)),unknown_grade=sum(gradegrp=="Unknown"),unknown_race=sum(race=="Unknown"),unknown_marital=sum(marital=="Unknown")),by=exposure]
  overall <- d[,.(n=.N,NX=sum(N=="NX"),unknown_size=sum(is.na(size)),unknown_grade=sum(gradegrp=="Unknown"),unknown_race=sum(race=="Unknown"),unknown_marital=sum(marital=="Unknown"),min_year=min(year),max_year=max(year))]
  # Balance maximum across each represented column/level.
  bvals <- c(smd(d$age,d$exposure01,d$ow),smd(d$year,d$exposure01,d$ow))
  for(vn in c("sex","race","marital","sitegrp","sizegrp","gradegrp","histgrp")) for(lv in levels(droplevels(d[[vn]]))) bvals<-c(bvals,smd(as.numeric(d[[vn]]==lv),d$exposure01,d$ow))
  cat("RESULT",cancer,"raw",raw_n,"t1",t1_n,"\n")
  print(as.data.table(do.call(rbind,steps)))
  print(overall);print(unk);print(ess);print(outcomes);print(ratio);cat("max weighted SMD",max(bvals,na.rm=TRUE),"\n")
  rm(d);gc()
  invisible(list(outcomes=outcomes,ratio=ratio,overall=overall,ess=ess))
}
