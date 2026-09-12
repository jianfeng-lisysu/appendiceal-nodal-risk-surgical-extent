## =============================================================================
## A6 RED-TEAM STEP 5 — ADEQUATE NODAL-ASSESSMENT SENSITIVITY (B=1000)
##
## Sensitivity definition:
##   - Limited-resection arm: retain all eligible patients.
##   - Oncologic-colectomy arm: retain only patients with >=12 examined nodes.
##
## Corrected framework:
##   - Summary Stage excluded from PS.
##   - Nodal-risk model re-estimated in the sensitivity cohort.
##   - Quintile / bottom-decile cut points re-estimated.
##   - PS/OW refitted separately WITHIN each predicted-risk stratum.
##   - Bottom predicted-risk decile remains the primary low-risk anchor.
##   - risk<5% remains supportive.
##   - B=1000 full-pipeline bootstrap.
##
## Expected sensitivity-cohort N from the frozen analysis: 6,130.
##
## OUTPUTS
##   03_results/REDTEAM_05_AdequateNodal_Report.txt
##   03_results/REDTEAM_05_AdequateNodal_B1000.csv
##   03_results/REDTEAM_05_AdequateNodal_Bootstrap_B1000.csv
##   03_results/REDTEAM_05_AdequateNodal_Checkpoint_B1000.rds
##
## NOTE
##   No new manuscript/table/figure version files are generated here.
##   Once this sensitivity is locked, the existing FINAL supplementary table
##   can be overwritten in place.
## =============================================================================

rm(list=ls())
options(stringsAsFactors=FALSE, scipen=999)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
DATA_DIR <- file.path(ROOT,"02_data")
RESULT_DIR <- file.path(ROOT,"03_results")
dir.create(RESULT_DIR,recursive=TRUE,showWarnings=FALSE)

B <- 1000L
SEED <- 20260905L

pick_existing <- function(paths) {
  z <- paths[file.exists(paths)]
  if (!length(z)) return(NA_character_)
  z[1]
}

PIPE_FILE <- pick_existing(c(
  file.path(DATA_DIR,"Aline_frozen_pipeline_v2.RData"),
  file.path(ROOT,"Aline_frozen_pipeline_v2.RData")
))
RAW_FILE <- pick_existing(c(
  file.path(DATA_DIR,"analytic_cohort_A.csv"),
  file.path(ROOT,"analytic_cohort_A.csv")
))

if (is.na(PIPE_FILE)) stop("Missing Aline_frozen_pipeline_v2.RData")
if (is.na(RAW_FILE)) stop("Missing analytic_cohort_A.csv")

load(PIPE_FILE)
if (!exists("pipeline",mode="function")) stop("pipeline() missing")
if (!exists("PS_FORMULA")) stop("PS_FORMULA missing")

RAW <- read.csv(RAW_FILE,stringsAsFactors=FALSE,check.names=FALSE)
if (nrow(RAW)!=6951) stop("Frozen cohort fingerprint failed: N=",nrow(RAW))

needed_raw <- c("exposure","adequate_ln")
miss <- setdiff(needed_raw,names(RAW))
if (length(miss)) stop("RAW missing: ",paste(miss,collapse=", "))

## -----------------------------------------------------------------------------
## 1. Sensitivity cohort
## -----------------------------------------------------------------------------

is_colectomy <- RAW$exposure == "oncologic_colectomy"
if (!any(is_colectomy,na.rm=TRUE)) {
  ## backward-compatible fallback if exposure labels differ
  is_colectomy <- RAW$exposure %in% c("colectomy","Oncologic colectomy",1,"1")
}

is_limited <- !is_colectomy & !is.na(RAW$exposure)

keep_s9 <- is_limited | (is_colectomy & RAW$adequate_ln==1)

RAW_S9 <- RAW[keep_s9 %in% TRUE,,drop=FALSE]

if (nrow(RAW_S9)!=6130) {
  stop(
    "Adequate-nodal sensitivity cohort fingerprint failed. Expected N=6130; current N=",
    nrow(RAW_S9),
    ". Audit exposure/adequate_ln coding before continuing."
  )
}

## -----------------------------------------------------------------------------
## 2. Remove Summary Stage from PS
## -----------------------------------------------------------------------------

PS_ORIG <- PS_FORMULA
tt <- attr(terms(PS_ORIG),"term.labels")
stage_terms <- tt[grepl("stage",tt,ignore.case=TRUE)]
if (!length(stage_terms)) stop("No stage-related PS term found.")

PS_NOSTAGE <- reformulate(
  setdiff(tt,stage_terms),
  response=all.vars(PS_ORIG)[1],
  env=environment(PS_ORIG)
)

penv <- environment(pipeline)
had_env <- exists("PS_FORMULA",envir=penv,inherits=FALSE)
old_env <- if (had_env) get("PS_FORMULA",envir=penv,inherits=FALSE) else NULL

on.exit({
  assign("PS_FORMULA",PS_ORIG,envir=.GlobalEnv)
  if (had_env) assign("PS_FORMULA",old_env,envir=penv)
},add=TRUE)

assign("PS_FORMULA",PS_NOSTAGE,envir=.GlobalEnv)
assign("PS_FORMULA",PS_NOSTAGE,envir=penv)

## -----------------------------------------------------------------------------
## 3. Common estimators
## -----------------------------------------------------------------------------

aj_w <- function(time,ev1,ev2,wt,start=-1,tau=60) {
  ok <- !is.na(time)&!is.na(ev1)&!is.na(ev2)&!is.na(wt)&is.finite(wt)
  time<-time[ok]; ev1<-ev1[ok]; ev2<-ev2[ok]; wt<-wt[ok]
  m <- time>start
  time<-time[m]; ev1<-ev1[m]; ev2<-ev2[m]; wt<-wt[m]
  if (!length(time)||sum(wt)<=0) return(NA_real_)
  dc<-rowsum(wt*ev1,time,reorder=FALSE)
  dr<-rowsum(wt*ev2,time,reorder=FALSE)
  dw<-rowsum(wt,time,reorder=FALSE)
  ut<-as.numeric(rownames(dc)); o<-order(ut)
  ut<-ut[o]; dc<-dc[o]; dr<-dr[o]; dw<-dw[o]
  S<-1; cif<-0; ar<-sum(wt)
  for (k in seq_along(ut)) {
    if (ut[k]>tau) break
    if (ar>0) {
      cif<-cif+S*dc[k]/ar
      S<-S*(1-(dc[k]+dr[k])/ar)
    }
    ar<-ar-dw[k]
  }
  as.numeric(cif)
}

make_outcome <- function(d) {
  z <- d[
    !is.na(d$css_event)&
    !is.na(d$os_months)&
    !is.na(d$os_event)&
    !is.na(d$treat)&
    !is.na(d$nodal_risk),
    ,drop=FALSE
  ]
  z$cancer <- as.integer(z$css_event==1 & z$os_months<=60)
  z$otherd <- as.integer(z$other_death==1 & z$os_months<=60)
  z
}

make_cuts <- function(x) {
  z <- as.numeric(quantile(x,0:5/5,na.rm=TRUE,names=FALSE))
  z[1]<-z[1]-1e-8
  z[6]<-z[6]+1e-8
  if (any(!is.finite(z))||any(diff(z)<=0)) stop("Invalid quintile cuts")
  z
}

make_masks <- function(d) {
  cuts <- make_cuts(d$nodal_risk)
  p10 <- as.numeric(quantile(d$nodal_risk,.10,na.rm=TRUE))
  grp <- cut(d$nodal_risk,breaks=cuts,labels=FALSE,include.lowest=TRUE)
  list(
    Q1=grp==1 & !is.na(grp),
    Q2=grp==2 & !is.na(grp),
    Q3=grp==3 & !is.na(grp),
    Q4=grp==4 & !is.na(grp),
    Q5=grp==5 & !is.na(grp),
    bottom_decile=d$nodal_risk<=p10,
    risk_lt_5pct=d$nodal_risk<.05
  )
}

fit_subgroup_ow <- function(g) {
  if (nrow(g)<40) stop("Subgroup too small")
  if (sum(g$treat==0)<15||sum(g$treat==1)<15)
    stop("Insufficient treatment support")
  fit <- suppressWarnings(glm(PS_NOSTAGE,data=g,family=binomial()))
  ps <- suppressWarnings(predict(fit,newdata=g,type="response"))
  if (length(ps)!=nrow(g)||any(!is.finite(ps))) stop("Non-finite subgroup PS")
  ps <- pmin(pmax(ps,1e-6),1-1e-6)
  g$ps_sub <- ps
  g$ow_sub <- ifelse(g$treat==1,1-ps,ps)
  if (any(!is.finite(g$ow_sub))||sum(g$ow_sub)<=0) stop("Invalid subgroup OW")
  g
}

estimate_one <- function(g) {
  gs <- fit_subgroup_ow(g)
  a1 <- gs[gs$treat==1,,drop=FALSE]
  a0 <- gs[gs$treat==0,,drop=FALSE]
  c1 <- aj_w(a1$os_months,a1$cancer,a1$otherd,a1$ow_sub)
  c0 <- aj_w(a0$os_months,a0$cancer,a0$otherd,a0$ow_sub)
  c(limited=c0,colectomy=c1,rd=c1-c0)
}

LABS <- c("Q1","Q2","Q3","Q4","Q5","bottom_decile","risk_lt_5pct")

estimate_all <- function(d) {
  z <- make_outcome(d)
  masks <- make_masks(z)
  out <- rep(NA_real_,length(LABS)); names(out)<-LABS
  for (nm in LABS) {
    g <- z[masks[[nm]],,drop=FALSE]
    out[nm] <- tryCatch(estimate_one(g)["rd"],error=function(e)NA_real_)
  }
  out
}

ess <- function(w) {
  w<-w[is.finite(w)&w>0]
  if (!length(w)) return(NA_real_)
  sum(w)^2/sum(w^2)
}

weighted_smd <- function(x,tr,w) {
  ok<-is.finite(x)&is.finite(tr)&is.finite(w)
  x<-x[ok]; tr<-tr[ok]; w<-w[ok]
  if (!sum(tr==0)||!sum(tr==1)) return(NA_real_)
  m0<-weighted.mean(x[tr==0],w[tr==0])
  m1<-weighted.mean(x[tr==1],w[tr==1])
  v0<-weighted.mean((x[tr==0]-m0)^2,w[tr==0])
  v1<-weighted.mean((x[tr==1]-m1)^2,w[tr==1])
  den<-sqrt((v0+v1)/2)
  if (!is.finite(den)||den<=0) return(0)
  (m1-m0)/den
}

balance_design <- function(g) {
  X <- model.matrix(delete.response(terms(PS_NOSTAGE)),
                    data=g,na.action=na.pass)
  if ("(Intercept)"%in%colnames(X))
    X<-X[,colnames(X)!="(Intercept)",drop=FALSE]
  cbind(X,nodal_risk=g$nodal_risk)
}

## -----------------------------------------------------------------------------
## 4. Point estimates + balance
## -----------------------------------------------------------------------------

d0 <- pipeline(RAW_S9)
z0 <- make_outcome(d0)
masks0 <- make_masks(z0)

point <- rep(NA_real_,length(LABS)); names(point)<-LABS
diag_rows <- list()

for (nm in LABS) {
  g <- z0[masks0[[nm]],,drop=FALSE]
  gs <- fit_subgroup_ow(g)
  est <- estimate_one(g)
  point[nm] <- est["rd"]
  a0<-gs[gs$treat==0,,drop=FALSE]
  a1<-gs[gs$treat==1,,drop=FALSE]
  X<-balance_design(gs)
  smd<-vapply(seq_len(ncol(X)),function(j)
    weighted_smd(X[,j],gs$treat,gs$ow_sub),numeric(1))
  names(smd)<-colnames(X)
  mx<-max(abs(smd),na.rm=TRUE)
  diag_rows[[nm]] <- data.frame(
    analysis=nm,
    n=nrow(gs),
    limited=nrow(a0),
    colectomy=nrow(a1),
    cancer_events_limited=sum(a0$cancer==1,na.rm=TRUE),
    cancer_events_colectomy=sum(a1$cancer==1,na.rm=TRUE),
    CIF_limited_pct=est["limited"]*100,
    CIF_colectomy_pct=est["colectomy"]*100,
    RD_pp=est["rd"]*100,
    ESS_limited=ess(a0$ow_sub),
    ESS_colectomy=ess(a1$ow_sub),
    max_weight=max(gs$ow_sub,na.rm=TRUE),
    max_abs_weighted_SMD=mx,
    worst_variable=names(which.max(abs(smd))),
    balance_pass=is.finite(mx)&&mx<.10,
    stringsAsFactors=FALSE
  )
}

DIAG <- do.call(rbind,diag_rows)

## -----------------------------------------------------------------------------
## 5. B=1000 full-pipeline bootstrap + checkpoint
## -----------------------------------------------------------------------------

CHECKPOINT <- file.path(
  RESULT_DIR,"REDTEAM_05_AdequateNodal_Checkpoint_B1000.rds"
)

boot <- matrix(
  NA_real_,nrow=B,ncol=length(LABS),
  dimnames=list(NULL,LABS)
)

start_b <- 1L

if (file.exists(CHECKPOINT)) {
  cp <- tryCatch(readRDS(CHECKPOINT),error=function(e)NULL)
  if (!is.null(cp)&&is.matrix(cp$boot)&&
      identical(dim(cp$boot),dim(boot))&&
      identical(colnames(cp$boot),colnames(boot))) {
    boot<-cp$boot
    done<-which(apply(boot,1,function(x)any(is.finite(x))))
    start_b<-if(length(done)) max(done)+1L else 1L
    if (!is.null(cp$rng_seed))
      assign(".Random.seed",cp$rng_seed,envir=.GlobalEnv)
    cat("Checkpoint loaded. Resume at ",start_b,"/",B,"\n",sep="")
  }
}

if (start_b==1L) set.seed(SEED)

if (start_b<=B) {
  t0<-Sys.time()
  for (b in seq.int(start_b,B)) {
    db <- RAW_S9[
      sample(seq_len(nrow(RAW_S9)),nrow(RAW_S9),replace=TRUE),
      ,drop=FALSE
    ]
    boot[b,] <- tryCatch({
      dx <- pipeline(db)
      estimate_all(dx)
    },error=function(e)rep(NA_real_,length(LABS)))
    if (b%%50==0||b==B) {
      saveRDS(
        list(boot=boot,completed_through=b,rng_seed=.Random.seed),
        CHECKPOINT
      )
      cat("Bootstrap ",b,"/",B," | elapsed ",
          round(difftime(Sys.time(),t0,units="mins"),1)," min\n",sep="")
    }
  }
}

## -----------------------------------------------------------------------------
## 6. Final results
## -----------------------------------------------------------------------------

res <- lapply(seq_along(LABS),function(i) {
  v<-boot[,i]
  valid<-is.finite(v)
  if (sum(valid)<900)
    warning(LABS[i],": effective B <900 (",sum(valid),")")
  lo<-as.numeric(quantile(v[valid],.025,na.rm=TRUE))
  hi<-as.numeric(quantile(v[valid],.975,na.rm=TRUE))
  data.frame(
    analysis=LABS[i],
    RD=point[i],
    RD_pp=point[i]*100,
    lo=lo,lo_pp=lo*100,
    hi=hi,hi_pp=hi*100,
    upper_possible_reduction_pp=max(0,-lo*100),
    effective_B=sum(valid),
    stringsAsFactors=FALSE
  )
})
RESULT <- do.call(rbind,res)

lowrisk <- c("Q1","bottom_decile","risk_lt_5pct")
RESULT$criterion_2pp <- ifelse(
  RESULT$analysis%in%lowrisk,
  RESULT$upper_possible_reduction_pp<2,
  NA
)
RESULT$criterion_3pp <- ifelse(
  RESULT$analysis%in%lowrisk,
  RESULT$upper_possible_reduction_pp<3,
  NA
)
RESULT$criterion_5pp <- ifelse(
  RESULT$analysis%in%lowrisk,
  RESULT$upper_possible_reduction_pp<5,
  NA
)

## full corrected primary reference (Step 1B)
ref_rd <- c(
  Q1=-0.170196,
  Q2=-3.604407,
  Q3=-0.204167,
  Q4=-4.432059,
  Q5=-5.066558,
  bottom_decile=2.053016,
  risk_lt_5pct=2.926840
)

RESULT$full_corrected_RD_pp <- ref_rd[RESULT$analysis]
RESULT$point_change_vs_full_pp <- RESULT$RD_pp-RESULT$full_corrected_RD_pp

## old adequate-node sensitivity reference (old/global framework)
old_s9 <- c(
  Q1=-0.6,
  Q2=-5.0,
  Q3=-1.5,
  Q4=-6.6,
  Q5=-7.0,
  bottom_decile=NA_real_,
  risk_lt_5pct=NA_real_
)
RESULT$old_S9_RD_pp <- old_s9[RESULT$analysis]
RESULT$point_change_vs_old_S9_pp <- RESULT$RD_pp-RESULT$old_S9_RD_pp

RESULT <- merge(
  RESULT,
  DIAG[,c(
    "analysis","n","limited","colectomy",
    "cancer_events_limited","cancer_events_colectomy",
    "CIF_limited_pct","CIF_colectomy_pct",
    "ESS_limited","ESS_colectomy","max_weight",
    "max_abs_weighted_SMD","worst_variable","balance_pass"
  )],
  by="analysis",all.x=TRUE,sort=FALSE
)

RESULT <- RESULT[match(LABS,RESULT$analysis),,drop=FALSE]

## -----------------------------------------------------------------------------
## 7. Score / decision checks
## -----------------------------------------------------------------------------

bottom <- RESULT[RESULT$analysis=="bottom_decile",,drop=FALSE]
risk5 <- RESULT[RESULT$analysis=="risk_lt_5pct",,drop=FALSE]

max_abs_shift <- max(abs(RESULT$point_change_vs_full_pp),na.rm=TRUE)

score_alert <- if (max_abs_shift>2) {
  paste0(
    "SCORE ALERT: at least one adequate-nodal sensitivity point estimate moved >2 pp ",
    "versus the full corrected cohort; maximum absolute change = ",
    sprintf("%.2f",max_abs_shift)," pp."
  )
} else {
  paste0(
    "No >2-pp point-estimate shift versus the full corrected cohort; ",
    "maximum absolute change = ",sprintf("%.2f",max_abs_shift)," pp."
  )
}

balance_alert <- if (all(RESULT$balance_pass %in% TRUE)) {
  paste0(
    "PASS: all seven risk strata have max|SMD| <0.10; worst = ",
    sprintf("%.3f",max(RESULT$max_abs_weighted_SMD,na.rm=TRUE)),"."
  )
} else {
  paste0(
    "SCORE ALERT: at least one risk stratum fails balance; worst max|SMD| = ",
    sprintf("%.3f",max(RESULT$max_abs_weighted_SMD,na.rm=TRUE)),"."
  )
}

## -----------------------------------------------------------------------------
## 8. Save
## -----------------------------------------------------------------------------

OUT_CSV <- file.path(
  RESULT_DIR,"REDTEAM_05_AdequateNodal_B1000.csv"
)
OUT_BOOT <- file.path(
  RESULT_DIR,"REDTEAM_05_AdequateNodal_Bootstrap_B1000.csv"
)
OUT_REPORT <- file.path(
  RESULT_DIR,"REDTEAM_05_AdequateNodal_Report.txt"
)

write.csv(RESULT,OUT_CSV,row.names=FALSE)
write.csv(data.frame(replicate=seq_len(B),boot,check.names=FALSE),
          OUT_BOOT,row.names=FALSE)

report <- c(
  "======================================================================",
  "A6 RED-TEAM STEP 5 — ADEQUATE NODAL-ASSESSMENT SENSITIVITY",
  "======================================================================",
  paste0("Frozen full cohort N = ",nrow(RAW)),
  paste0("Sensitivity cohort N = ",nrow(RAW_S9)," (expected 6130)"),
  paste0("Outcome-evaluable sensitivity N = ",nrow(z0)),
  paste0("Removed PS term(s): ",paste(stage_terms,collapse=" | ")),
  "Limited-resection arm retained without LN-count restriction.",
  "Oncologic-colectomy arm restricted to >=12 examined nodes.",
  "PS/OW refitted separately within each predicted-risk stratum.",
  "Bottom decile remains primary low-risk anchor; risk<5% remains supportive.",
  "",
  "== CORRECTED ADEQUATE-NODAL B=1000 =="
)

for (i in seq_len(nrow(RESULT))) {
  z<-RESULT[i,]
  crit <- ""
  if (z$analysis%in%lowrisk) {
    crit <- sprintf(
      " | upper reduction %.2f pp | 2pp:%s 3pp:%s 5pp:%s",
      z$upper_possible_reduction_pp,
      ifelse(z$criterion_2pp,"MET","NOT MET"),
      ifelse(z$criterion_3pp,"MET","NOT MET"),
      ifelse(z$criterion_5pp,"MET","NOT MET")
    )
  }
  report<-c(report,paste0(
    sprintf(
      "%-15s RD %+.2f pp (95%% CI %+.2f to %+.2f; B=%d)",
      z$analysis,z$RD_pp,z$lo_pp,z$hi_pp,z$effective_B
    ),
    crit,
    sprintf(
      " | change vs full %+.2f pp | max|SMD| %.3f [%s] %s",
      z$point_change_vs_full_pp,
      z$max_abs_weighted_SMD,
      z$worst_variable,
      ifelse(z$balance_pass,"PASS","FAIL")
    )
  ))
}

report<-c(
  report,
  "",
  "== SUPPORT / EVENTS ==",
  apply(RESULT,1,function(z) {
    sprintf(
      "%-15s n=%d (%d/%d) | cancer events=%d/%d | CIF %.2f/%.2f%% | ESS %.1f/%.1f",
      z[["analysis"]],
      as.integer(z[["n"]]),
      as.integer(z[["limited"]]),
      as.integer(z[["colectomy"]]),
      as.integer(z[["cancer_events_limited"]]),
      as.integer(z[["cancer_events_colectomy"]]),
      as.numeric(z[["CIF_limited_pct"]]),
      as.numeric(z[["CIF_colectomy_pct"]]),
      as.numeric(z[["ESS_limited"]]),
      as.numeric(z[["ESS_colectomy"]])
    )
  }),
  "",
  "== RED-TEAM / SCORE CHECK ==",
  balance_alert,
  score_alert,
  paste0(
    "Bottom-decile 3-pp criterion: ",
    ifelse(isTRUE(bottom$criterion_3pp),"MET","NOT MET"),
    " | upper compatible reduction ",
    sprintf("%.2f",bottom$upper_possible_reduction_pp)," pp."
  ),
  paste0(
    "risk<5% 3-pp criterion: ",
    ifelse(isTRUE(risk5$criterion_3pp),"MET","NOT MET"),
    " | upper compatible reduction ",
    sprintf("%.2f",risk5$upper_possible_reduction_pp)," pp."
  ),
  "",
  "Interpretation:",
  "This is a sensitivity analysis addressing adequacy of nodal assessment in the oncologic-colectomy arm.",
  "It must not be described as a cohort in which all patients had >=12 examined nodes.",
  "The key robustness checks are point-estimate stability, within-stratum balance, and whether the extreme low-risk interpretation materially changes.",
  "======================================================================"
)

writeLines(report,OUT_REPORT,useBytes=TRUE)
cat(paste(report,collapse="\n"),"\n")
cat("\nREPORT: ",OUT_REPORT,"\n",sep="")
