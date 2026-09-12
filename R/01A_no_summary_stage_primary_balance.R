
rm(list=ls())
options(stringsAsFactors=FALSE, scipen=999)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
D <- file.path(ROOT,"02_data")
R <- file.path(ROOT,"03_results")

pick <- function(x) { z <- x[file.exists(x)]; if(length(z)) z[1] else NA_character_ }

PIPE <- pick(c(file.path(D,"Aline_frozen_pipeline_v2.RData"),
               file.path(ROOT,"Aline_frozen_pipeline_v2.RData")))
RAWF <- pick(c(file.path(D,"analytic_cohort_A.csv"),
               file.path(ROOT,"analytic_cohort_A.csv")))

if(is.na(PIPE) || is.na(RAWF)) stop("Missing frozen pipeline or analytic cohort.")
load(PIPE)
RAW <- read.csv(RAWF, stringsAsFactors=FALSE, check.names=FALSE)
if(nrow(RAW)!=6951) stop("Cohort fingerprint failed.")
if(!exists("pipeline",mode="function")) stop("pipeline() missing.")
if(!exists("PS_FORMULA")) stop("PS_FORMULA missing.")

## ---------- 1. Remove Summary Stage from PS formula ----------
PS_ORIG <- PS_FORMULA
tt <- attr(terms(PS_ORIG),"term.labels")
stage_terms <- tt[grepl("stage",tt,ignore.case=TRUE)]
if(!length(stage_terms)) stop("No stage term found in PS_FORMULA.")
PS_NOSTAGE <- reformulate(setdiff(tt,stage_terms),
                          response=all.vars(PS_ORIG)[1],
                          env=environment(PS_ORIG))

d_orig <- pipeline(RAW)

## pipeline() resolves PS_FORMULA in its enclosing environment.
penv <- environment(pipeline)
old_in_env <- if(exists("PS_FORMULA",envir=penv,inherits=FALSE))
  get("PS_FORMULA",envir=penv,inherits=FALSE) else NULL
had_in_env <- exists("PS_FORMULA",envir=penv,inherits=FALSE)

assign("PS_FORMULA",PS_NOSTAGE,envir=.GlobalEnv)
assign("PS_FORMULA",PS_NOSTAGE,envir=penv)

restore_formula <- function(){
  assign("PS_FORMULA",PS_ORIG,envir=.GlobalEnv)
  if(had_in_env) assign("PS_FORMULA",old_in_env,envir=penv)
}
on.exit(restore_formula(), add=TRUE)

d0 <- pipeline(RAW)
if(max(abs(d0$nodal_risk-d_orig$nodal_risk),na.rm=TRUE)>1e-10)
  stop("Nodal-risk predictions changed unexpectedly.")
if(max(abs(d0$ps-d_orig$ps),na.rm=TRUE)<1e-8)
  stop("PS did not change; inspect pipeline formula environment.")

## ---------- 2. Frozen primary estimators ----------
aj_w <- function(time,ev1,ev2,wt,start=-1,tau=60){
  ok <- !is.na(time)&!is.na(ev1)&!is.na(ev2)&!is.na(wt)&is.finite(wt)
  time<-time[ok]; ev1<-ev1[ok]; ev2<-ev2[ok]; wt<-wt[ok]
  m<-time>start; time<-time[m]; ev1<-ev1[m]; ev2<-ev2[m]; wt<-wt[m]
  if(!length(time)||sum(wt)<=0) return(NA_real_)
  dc<-rowsum(wt*ev1,time,reorder=FALSE)
  dr<-rowsum(wt*ev2,time,reorder=FALSE)
  dw<-rowsum(wt,time,reorder=FALSE)
  ut<-as.numeric(rownames(dc)); o<-order(ut)
  ut<-ut[o]; dc<-dc[o]; dr<-dr[o]; dw<-dw[o]
  S<-1; cif<-0; ar<-sum(wt)
  for(k in seq_along(ut)){
    if(ut[k]>tau) break
    if(ar>0){
      cif<-cif+S*dc[k]/ar
      S<-S*(1-(dc[k]+dr[k])/ar)
    }
    ar<-ar-dw[k]
  }
  as.numeric(cif)
}

mkdd <- function(d){
  z<-d[!is.na(d$css_event)&!is.na(d$os_months)&!is.na(d$ow)&
       !is.na(d$treat)&!is.na(d$nodal_risk),]
  z$cancer<-as.integer(z$css_event==1 & z$os_months<=60)
  z$otherd<-as.integer(z$other_death==1 & z$os_months<=60)
  z
}

cuts5 <- function(x){
  q<-as.numeric(quantile(x,0:5/5,na.rm=TRUE,names=FALSE))
  q[1]<-q[1]-1e-8; q[6]<-q[6]+1e-8; q
}

est7 <- function(d){
  cuts<-cuts5(d$nodal_risk)
  p10<-as.numeric(quantile(d$nodal_risk,.10,na.rm=TRUE))
  z<-mkdd(d)
  g<-cut(z$nodal_risk,cuts,labels=FALSE,include.lowest=TRUE)
  masks<-c(lapply(1:5,function(q) g==q & !is.na(g)),
           list(z$nodal_risk<=p10,z$nodal_risk<.05))
  sapply(masks,function(m){
    x<-z[m,]; a1<-x[x$treat==1,]; a0<-x[x$treat==0,]
    if(nrow(a1)<10||nrow(a0)<10) return(NA_real_)
    c1<-aj_w(a1$os_months,a1$cancer,a1$otherd,a1$ow)
    c0<-aj_w(a0$os_months,a0$cancer,a0$otherd,a0$ow)
    c1-c0
  })
}

overall_rd <- function(d){
  z<-mkdd(d); a1<-z[z$treat==1,]; a0<-z[z$treat==0,]
  aj_w(a1$os_months,a1$cancer,a1$otherd,a1$ow) -
    aj_w(a0$os_months,a0$cancer,a0$otherd,a0$ow)
}

labs <- c("Q1","Q2","Q3","Q4","Q5","bottom_decile","risk_lt_5pct")
p_old <- est7(d_orig)
p_new <- est7(d0)
o_old <- overall_rd(d_orig)
o_new <- overall_rd(d0)

## ---------- 3. Full-pipeline bootstrap B=1000 ----------
set.seed(20260831)
B<-1000
bb<-matrix(NA_real_,B,8)
colnames(bb)<-c(labs,"overall")
t0<-Sys.time()

for(b in 1:B){
  db<-RAW[sample(seq_len(nrow(RAW)),nrow(RAW),replace=TRUE),]
  bb[b,]<-tryCatch({
    dx<-pipeline(db)
    c(est7(dx),overall_rd(dx))
  },error=function(e) rep(NA_real_,8))
  if(b%%100==0) cat("bootstrap",b,"/",B,
    round(difftime(Sys.time(),t0,units="mins"),1),"min\n")
}
bv<-bb[complete.cases(bb),,drop=FALSE]
if(nrow(bv)<900) stop("Effective B <900.")

pt<-c(p_new,o_new)
res<-data.frame(
  analysis=c(labs,"overall"),
  estimate=pt,
  rd_pp=pt*100,
  lo=apply(bv,2,quantile,.025,na.rm=TRUE),
  hi=apply(bv,2,quantile,.975,na.rm=TRUE),
  effective_B=nrow(bv),
  stringsAsFactors=FALSE
)
res$lo_pp<-res$lo*100
res$hi_pp<-res$hi*100
res$upper_possible_reduction_pp<-pmax(0,-res$lo_pp)
res$criterion_2pp<-ifelse(res$analysis%in%c("Q1","bottom_decile","risk_lt_5pct"),
                          res$upper_possible_reduction_pp<2,NA)
res$criterion_3pp<-ifelse(res$analysis%in%c("Q1","bottom_decile","risk_lt_5pct"),
                          res$upper_possible_reduction_pp<3,NA)
res$criterion_5pp<-ifelse(res$analysis%in%c("Q1","bottom_decile","risk_lt_5pct"),
                          res$upper_possible_reduction_pp<5,NA)

## ---------- 4. Within-stratum balance/support audit ----------
ess<-function(w){w<-w[is.finite(w)&w>0]; sum(w)^2/sum(w^2)}
smd<-function(x,tr,w){
  ok<-is.finite(x)&is.finite(tr)&is.finite(w)
  x<-x[ok];tr<-tr[ok];w<-w[ok]
  m1<-weighted.mean(x[tr==1],w[tr==1]);m0<-weighted.mean(x[tr==0],w[tr==0])
  v1<-weighted.mean((x[tr==1]-m1)^2,w[tr==1])
  v0<-weighted.mean((x[tr==0]-m0)^2,w[tr==0])
  s<-sqrt((v1+v0)/2); if(!is.finite(s)||s==0) 0 else (m1-m0)/s
}
mm<-function(d){
  gr<-suppressWarnings(as.numeric(d$grade))
  data.frame(
    age=d$age,age2=d$age^2,female=d$female,year=d$year,
    logsize=d$logsize,size_miss=d$size_miss,
    T2=as.numeric(d$T=="T2"),T3=as.numeric(d$T=="T3"),T4=as.numeric(d$T=="T4"),
    grade1=as.numeric(gr==1),grade2=as.numeric(gr==2),
    grade3=as.numeric(gr==3),grade4=as.numeric(gr==4),
    nonmuc=as.numeric(d$hist_group=="nonmucinous"),
    muc=as.numeric(d$hist_group=="mucinous"),
    SRCC=as.numeric(d$hist_group=="SRCC"),
    race_black=as.numeric(d$race4=="Non-Hispanic Black"),
    race_white=as.numeric(d$race4=="Non-Hispanic White"),
    race_other=as.numeric(d$race4=="Other/Unknown"),
    income_hi=d$income_hi,metro=d$metro,married=d$married,
    nodal_risk=d$nodal_risk
  )
}

z<-mkdd(d0)
ct<-cuts5(d0$nodal_risk)
p10<-as.numeric(quantile(d0$nodal_risk,.10,na.rm=TRUE))
gg<-cut(z$nodal_risk,ct,labels=FALSE,include.lowest=TRUE)
masks<-c(setNames(lapply(1:5,function(q) gg==q & !is.na(gg)),paste0("Q",1:5)),
         list(bottom_decile=z$nodal_risk<=p10,
              risk_lt_5pct=z$nodal_risk<.05,
              overall=rep(TRUE,nrow(z))))

diag_list<-list()
bal_list<-list()

for(nm in names(masks)){
  x<-z[masks[[nm]],,drop=FALSE]
  bmat<-mm(x)
  sv<-vapply(names(bmat),function(v)smd(bmat[[v]],x$treat,x$ow),numeric(1))
  a0<-x[x$treat==0,];a1<-x[x$treat==1,]
  c0<-aj_w(a0$os_months,a0$cancer,a0$otherd,a0$ow)
  c1<-aj_w(a1$os_months,a1$cancer,a1$otherd,a1$ow)
  diag_list[[nm]]<-data.frame(
    analysis=nm,n=nrow(x),limited=nrow(a0),colectomy=nrow(a1),
    cancer_events_limited=sum(a0$cancer),cancer_events_colectomy=sum(a1$cancer),
    cif_limited_pct=c0*100,cif_colectomy_pct=c1*100,rd_pp=(c1-c0)*100,
    ess_limited=ess(a0$ow),ess_colectomy=ess(a1$ow),
    max_weight=max(x$ow,na.rm=TRUE),
    max_abs_weighted_smd=max(abs(sv),na.rm=TRUE),
    worst_variable=names(which.max(abs(sv))),
    balance_pass=max(abs(sv),na.rm=TRUE)<.10
  )
  bal_list[[nm]]<-data.frame(analysis=nm,variable=names(sv),weighted_smd=sv)
}
diag<-do.call(rbind,diag_list)
bal<-do.call(rbind,bal_list)

cmp<-data.frame(
  analysis=c(labs,"overall"),
  original_point_pp=c(p_old,o_old)*100,
  no_stage_point_pp=c(p_new,o_new)*100
)
cmp$change_pp<-cmp$no_stage_point_pp-cmp$original_point_pp

## ---------- 5. Save ----------
dir.create(R,recursive=TRUE,showWarnings=FALSE)
fres<-file.path(R,"REDTEAM_01A_NoSummaryStage_Primary_B1000.csv")
fboot<-file.path(R,"REDTEAM_01A_NoSummaryStage_Bootstrap_B1000.csv")
fdiag<-file.path(R,"REDTEAM_01A_NoSummaryStage_Stratum_Diagnostics.csv")
fbal<-file.path(R,"REDTEAM_01A_NoSummaryStage_Balance_Components.csv")
fcmp<-file.path(R,"REDTEAM_01A_NoSummaryStage_Comparison.csv")
frep<-file.path(R,"REDTEAM_01A_NoSummaryStage_Report.txt")

write.csv(res,fres,row.names=FALSE)
write.csv(data.frame(replicate=1:B,bb,check.names=FALSE),fboot,row.names=FALSE)
write.csv(diag,fdiag,row.names=FALSE)
write.csv(bal,fbal,row.names=FALSE)
write.csv(cmp,fcmp,row.names=FALSE)

out<-c(
"======================================================================",
"A6 RED-TEAM STEP 1A — REMOVE SUMMARY STAGE FROM PS",
"======================================================================",
paste0("Original PS formula: ",paste(deparse(PS_ORIG),collapse=" ")),
paste0("Removed term(s): ",paste(stage_terms,collapse=" | ")),
paste0("No-stage PS formula: ",paste(deparse(PS_NOSTAGE),collapse=" ")),
"",
"== ORIGINAL vs NO-STAGE POINT ESTIMATES ==",
apply(cmp,1,function(x)sprintf("%-15s original %+.2f | no-stage %+.2f | change %+.2f pp",
                               x[1],as.numeric(x[2]),as.numeric(x[3]),as.numeric(x[4]))),
"",
paste0("== NO-STAGE B=",nrow(bv)," =="),
apply(res,1,function(x)sprintf("%-15s RD %+.1f pp (95%% CI %+.1f to %+.1f), upper reduction %.1f",
                               x["analysis"],as.numeric(x["rd_pp"]),
                               as.numeric(x["lo_pp"]),as.numeric(x["hi_pp"]),
                               as.numeric(x["upper_possible_reduction_pp"]))),
"",
"== STRATUM DIAGNOSTICS ==",
apply(diag,1,function(x)sprintf("%-15s n=%s (%s/%s) events=%s/%s CIF=%.1f/%.1f RD=%+.1f ESS=%.1f/%.1f max|SMD|=%.3f [%s] %s",
                               x["analysis"],x["n"],x["limited"],x["colectomy"],
                               x["cancer_events_limited"],x["cancer_events_colectomy"],
                               as.numeric(x["cif_limited_pct"]),as.numeric(x["cif_colectomy_pct"]),
                               as.numeric(x["rd_pp"]),as.numeric(x["ess_limited"]),
                               as.numeric(x["ess_colectomy"]),as.numeric(x["max_abs_weighted_smd"]),
                               x["worst_variable"],ifelse(as.logical(x["balance_pass"]),"PASS","FAIL"))),
"",
"Primary low-risk anchor remains BOTTOM DECILE; <5% remains supportive.",
"Do not switch the primary anchor after seeing the results.",
"======================================================================"
)
writeLines(out,frep,useBytes=TRUE)
cat(paste(out,collapse="\n"),"\n")
cat("\nREPORT: ",frep,"\n",sep="")
