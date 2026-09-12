## =============================================================================
## A6 RED-TEAM STEP 4 — CORRECTED CONTINUOUS FIGURE 4 (B=1000)
## No Summary Stage in PS + quintile-localized PS/OW
## Secondary/descriptive continuous-risk analysis; bottom decile remains primary.
## =============================================================================

rm(list=ls())
options(stringsAsFactors=FALSE, scipen=999)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
DATA_DIR <- file.path(ROOT,"02_data")
RESULT_DIR <- file.path(ROOT,"03_results")
FIG_DIR <- file.path(ROOT,"04_figures")
dir.create(RESULT_DIR,recursive=TRUE,showWarnings=FALSE)
dir.create(FIG_DIR,recursive=TRUE,showWarnings=FALSE)

B <- 1000L
SEED <- 20260904L
GRID <- c(.03,.04,.05,.07,.10,.14,.20,.28,.40,.50)
OLD_GRID_RD_PP <- c(2.6,1.3,0.3,-1.0,-2.3,-3.3,-3.8,-3.8,-3.4,-3.0)
OLD_DELTA_PP <- -3.29
OLD_INTERACTION_P <- 0.2264

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

for (pkg in c("survival","splines","ggplot2")) {
  if (!requireNamespace(pkg,quietly=TRUE))
    install.packages(pkg,repos="https://cloud.r-project.org")
  if (!requireNamespace(pkg,quietly=TRUE))
    stop("Required package unavailable: ",pkg)
}

load(PIPE_FILE)
if (!exists("pipeline",mode="function")) stop("pipeline() missing")
if (!exists("PS_FORMULA")) stop("PS_FORMULA missing")

RAW <- read.csv(RAW_FILE,stringsAsFactors=FALSE,check.names=FALSE)
if (nrow(RAW)!=6951) stop("Frozen cohort fingerprint failed: N=",nrow(RAW))

## ---- remove Summary Stage from PS -------------------------------------------
PS_ORIG <- PS_FORMULA
tt <- attr(terms(PS_ORIG),"term.labels")
stage_terms <- tt[grepl("stage",tt,ignore.case=TRUE)]
if (!length(stage_terms)) stop("No stage term found in PS_FORMULA")
PS_NOSTAGE <- reformulate(setdiff(tt,stage_terms),
                          response=all.vars(PS_ORIG)[1],
                          env=environment(PS_ORIG))

penv <- environment(pipeline)
had_env <- exists("PS_FORMULA",envir=penv,inherits=FALSE)
old_env <- if (had_env) get("PS_FORMULA",envir=penv,inherits=FALSE) else NULL
on.exit({
  assign("PS_FORMULA",PS_ORIG,envir=.GlobalEnv)
  if (had_env) assign("PS_FORMULA",old_env,envir=penv)
},add=TRUE)
assign("PS_FORMULA",PS_NOSTAGE,envir=.GlobalEnv)
assign("PS_FORMULA",PS_NOSTAGE,envir=penv)

## ---- helpers ----------------------------------------------------------------
make_quintile <- function(risk) {
  cuts <- as.numeric(quantile(risk,0:5/5,na.rm=TRUE,names=FALSE))
  cuts[1] <- cuts[1]-1e-8
  cuts[6] <- cuts[6]+1e-8
  if (any(!is.finite(cuts)) || any(diff(cuts)<=0))
    stop("Invalid predicted-risk quintile cuts")
  cut(risk,breaks=cuts,labels=FALSE,include.lowest=TRUE)
}

fit_local_ow <- function(d) {
  d$risk_quintile_local <- make_quintile(d$nodal_risk)
  d$ps_local <- NA_real_
  d$ow_local <- NA_real_
  for (q in 1:5) {
    ii <- which(d$risk_quintile_local==q)
    g <- d[ii,,drop=FALSE]
    if (sum(g$treat==0,na.rm=TRUE)<20 || sum(g$treat==1,na.rm=TRUE)<20)
      stop("Insufficient treatment support in Q",q)
    fit <- suppressWarnings(glm(PS_NOSTAGE,data=g,family=binomial()))
    ps <- suppressWarnings(predict(fit,newdata=g,type="response"))
    if (length(ps)!=nrow(g) || any(!is.finite(ps)))
      stop("Non-finite localized PS in Q",q)
    ps <- pmin(pmax(ps,1e-6),1-1e-6)
    d$ps_local[ii] <- ps
    d$ow_local[ii] <- ifelse(g$treat==1,1-ps,ps)
  }
  if (any(!is.finite(d$ow_local)) || any(d$ow_local<=0))
    stop("Localized OW construction failed")
  d
}

weighted_smd <- function(x,tr,w) {
  ok <- is.finite(x)&is.finite(tr)&is.finite(w)
  x<-x[ok]; tr<-tr[ok]; w<-w[ok]
  if (!sum(tr==0)||!sum(tr==1)) return(NA_real_)
  m0<-weighted.mean(x[tr==0],w[tr==0]); m1<-weighted.mean(x[tr==1],w[tr==1])
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
    X <- X[,colnames(X)!="(Intercept)",drop=FALSE]
  cbind(X,nodal_risk=g$nodal_risk)
}

audit_balance <- function(d) {
  p10 <- as.numeric(quantile(d$nodal_risk,.10,na.rm=TRUE))
  masks <- list(
    Q1=d$risk_quintile_local==1,Q2=d$risk_quintile_local==2,
    Q3=d$risk_quintile_local==3,Q4=d$risk_quintile_local==4,
    Q5=d$risk_quintile_local==5,
    bottom_decile_under_quintile_weights=d$nodal_risk<=p10,
    risk_lt_5pct_under_quintile_weights=d$nodal_risk<.05
  )
  do.call(rbind,lapply(names(masks),function(nm) {
    g <- d[masks[[nm]] & !is.na(masks[[nm]]),,drop=FALSE]
    X <- balance_design(g)
    smd <- vapply(seq_len(ncol(X)),function(j)
      weighted_smd(X[,j],g$treat,g$ow_local),numeric(1))
    names(smd)<-colnames(X)
    mx<-max(abs(smd),na.rm=TRUE)
    data.frame(
      analysis=nm,n=nrow(g),
      limited=sum(g$treat==0,na.rm=TRUE),
      colectomy=sum(g$treat==1,na.rm=TRUE),
      max_abs_weighted_SMD=mx,
      worst_variable=names(which.max(abs(smd))),
      balance_pass=is.finite(mx)&&mx<.10,
      stringsAsFactors=FALSE
    )
  }))
}

make_curve_data <- function(d) {
  dd <- d[!is.na(d$css_event)&!is.na(d$os_months)&!is.na(d$os_event)&
          !is.na(d$ow_local)&!is.na(d$treat)&!is.na(d$nodal_risk),,drop=FALSE]
  Tt <- dd$os_months
  dead <- dd$os_event==1
  cancer <- as.numeric(dd$css_event==1 & Tt<=60)
  known <- (Tt>=60)|dead
  G <- survival::survfit(survival::Surv(Tt,1-dd$os_event)~1)
  Gfun <- stepfun(G$time,c(1,G$surv))
  Gval <- Gfun(pmin(Tt,60)-1e-6)
  ipcw <- ifelse(known,1/pmax(Gval,1e-3),0)
  list(dd=dd,cancer=cancer,w=ipcw*dd$ow_local)
}

make_ns <- function(risk,kn)
  splines::ns(log(risk),knots=kn[2:3],Boundary.knots=kn[c(1,4)])

rd_fit_moving <- function(d) {
  kn <- as.numeric(quantile(log(d$nodal_risk),c(.05,.35,.65,.95),na.rm=TRUE))
  p10 <- as.numeric(quantile(d$nodal_risk,.10,na.rm=TRUE))
  p90 <- as.numeric(quantile(d$nodal_risk,.90,na.rm=TRUE))
  z <- make_curve_data(d); dd<-z$dd
  bsx <- make_ns(dd$nodal_risk,kn); tr<-dd$treat
  fit <- suppressWarnings(lm(z$cancer~tr*bsx,weights=z$w))
  cf<-coef(fit); ic<-cf[grep("^tr:",names(cf))]
  if (!is.finite(cf["tr"])||length(ic)!=3||any(!is.finite(ic)))
    stop("Invalid continuous model coefficients")
  rd <- function(rs) as.numeric(cf["tr"]+predict(bsx,log(rs))%*%ic)
  c(rd(GRID),rd(p10),rd(p90))
}

interaction_fixed <- function(d,kn_fixed) {
  z<-make_curve_data(d); dd<-z$dd
  bsx<-make_ns(dd$nodal_risk,kn_fixed); tr<-dd$treat
  fit<-suppressWarnings(lm(z$cancer~tr*bsx,weights=z$w))
  ic<-coef(fit)[grep("^tr:",names(coef(fit)))]
  if (length(ic)!=3||any(!is.finite(ic))) stop("Invalid fixed-basis interaction")
  as.numeric(ic)
}

wald_from_boot_cov <- function(beta,V) {
  V<-(V+t(V))/2
  ee<-eigen(V,symmetric=TRUE)
  tol<-max(ee$values)*1e-10
  keep<-ee$values>tol
  if (!any(keep)) return(list(W=NA_real_,df=NA_integer_,p=NA_real_))
  Q<-ee$vectors[,keep,drop=FALSE]; lam<-ee$values[keep]
  z<-as.numeric(t(Q)%*%beta)
  W<-sum(z^2/lam); df<-length(lam)
  list(W=W,df=df,p=pchisq(W,df,lower.tail=FALSE))
}

build_corrected <- function(raw_data) fit_local_ow(pipeline(raw_data))

## ---- original sample ---------------------------------------------------------
cat("============================================================\n")
cat("A6 CORRECTED CONTINUOUS FIGURE 4\n")
cat("============================================================\n")
cat("Removed PS term(s):",paste(stage_terms,collapse=" | "),"\n")

d0 <- build_corrected(RAW)
BALANCE <- audit_balance(d0)
print(BALANCE)

KN_FIXED <- as.numeric(quantile(log(d0$nodal_risk),
                                c(.05,.35,.65,.95),na.rm=TRUE))
P10 <- as.numeric(quantile(d0$nodal_risk,.10,na.rm=TRUE))
P90 <- as.numeric(quantile(d0$nodal_risk,.90,na.rm=TRUE))
rdf0 <- rd_fit_moving(d0)
ic0 <- interaction_fixed(d0,KN_FIXED)

cat("Corrected points:",
    paste(sprintf("%.0f%%=%+.2fpp",GRID*100,rdf0[1:10]*100),
          collapse=" | "),"\n")
cat(sprintf("P10=%.2f%% | P90=%.2f%%\n\n",P10*100,P90*100))

## ---- bootstrap with checkpoint ----------------------------------------------
CHECKPOINT_FILE <- file.path(
  RESULT_DIR,"REDTEAM_04_Corrected_Continuous_Checkpoint_B1000.rds"
)
rbf <- matrix(NA_real_,B,15)
colnames(rbf) <- c(paste0("RD_",GRID*100,"pct"),"RD_P10","RD_P90",
                   "interaction_1","interaction_2","interaction_3")
start_b <- 1L

if (file.exists(CHECKPOINT_FILE)) {
  cp <- tryCatch(readRDS(CHECKPOINT_FILE),error=function(e)NULL)
  if (!is.null(cp)&&is.matrix(cp$rbf)&&identical(dim(cp$rbf),dim(rbf))&&
      identical(colnames(cp$rbf),colnames(rbf))) {
    rbf<-cp$rbf
    done<-which(apply(rbf,1,function(x)any(is.finite(x))))
    start_b<-if(length(done)) max(done)+1L else 1L
    if (!is.null(cp$rng_seed)) assign(".Random.seed",cp$rng_seed,envir=.GlobalEnv)
    cat("Checkpoint loaded. Resume at",start_b,"/",B,"\n")
  }
}
if (start_b==1L) set.seed(SEED)

if (start_b<=B) {
  t0<-Sys.time()
  for (b in seq.int(start_b,B)) {
    db<-RAW[sample(seq_len(nrow(RAW)),nrow(RAW),replace=TRUE),,drop=FALSE]
    rbf[b,]<-tryCatch({
      dx<-build_corrected(db)
      c(rd_fit_moving(dx),interaction_fixed(dx,KN_FIXED))
    },error=function(e)rep(NA_real_,15))
    if (b%%50==0||b==B) {
      saveRDS(list(rbf=rbf,completed_through=b,rng_seed=.Random.seed),
              CHECKPOINT_FILE)
      cat("Bootstrap",b,"/",B,"| elapsed",
          round(difftime(Sys.time(),t0,units="mins"),1),"min\n")
    }
  }
}

valid_curve <- apply(rbf[,1:12,drop=FALSE],1,function(x)all(is.finite(x)))
valid_int <- apply(rbf[,13:15,drop=FALSE],1,function(x)all(is.finite(x)))
eff_curve <- sum(valid_curve); eff_int <- sum(valid_int)
if (eff_curve<900) stop("Effective curve B <900: ",eff_curve)
if (eff_int<900) stop("Effective interaction B <900: ",eff_int)

rb_curve<-rbf[valid_curve,1:12,drop=FALSE]
rb_int<-rbf[valid_int,13:15,drop=FALSE]

lo<-apply(rb_curve[,1:10,drop=FALSE],2,quantile,.025,na.rm=TRUE)
hi<-apply(rb_curve[,1:10,drop=FALSE],2,quantile,.975,na.rm=TRUE)

CURVE <- data.frame(
  risk=GRID,risk_pct=GRID*100,
  rd=rdf0[1:10],rd_pp=rdf0[1:10]*100,
  lo=lo,lo_pp=lo*100,hi=hi,hi_pp=hi*100,
  old_rd_pp=OLD_GRID_RD_PP,
  point_change_vs_old_pp=rdf0[1:10]*100-OLD_GRID_RD_PP,
  stringsAsFactors=FALSE
)

db <- (rb_curve[,12]-rb_curve[,11])*100
delta_point <- (rdf0[12]-rdf0[11])*100
delta_lo <- as.numeric(quantile(db,.025,na.rm=TRUE))
delta_hi <- as.numeric(quantile(db,.975,na.rm=TRUE))
delta_frac <- mean(db<0,na.rm=TRUE)
wald <- wald_from_boot_cov(ic0,cov(rb_int,use="complete.obs"))

## ---- checks ------------------------------------------------------------------
worst_q <- max(BALANCE$max_abs_weighted_SMD[
  BALANCE$analysis%in%paste0("Q",1:5)],na.rm=TRUE)
worst_tail <- max(BALANCE$max_abs_weighted_SMD[
  BALANCE$analysis%in%c("bottom_decile_under_quintile_weights",
                        "risk_lt_5pct_under_quintile_weights")],na.rm=TRUE)

max_shift <- max(abs(CURVE$point_change_vs_old_pp),na.rm=TRUE)
low_shift <- max(abs(CURVE$point_change_vs_old_pp[CURVE$risk<=.05]),na.rm=TRUE)
new_pattern <- mean(CURVE$rd_pp[CURVE$risk>=.20])-
               mean(CURVE$rd_pp[CURVE$risk<=.05])

balance_alert <- if (worst_q < .10) {
  sprintf(
    "PASS: all localized quintiles balanced; worst max|SMD| %.3f.",
    worst_q
  )
} else {
  sprintf(
    "SCORE ALERT: localized quintile imbalance; worst max|SMD| %.3f.",
    worst_q
  )
}

tail_alert <- if (worst_tail < .10) {
  sprintf(
    "PASS: extreme low-risk tail also <0.10 under inherited quintile weights; worst %.3f.",
    worst_tail
  )
} else {
  sprintf(
    paste0(
      "CAUTION: extreme low-risk tail imbalance under inherited quintile weights; ",
      "worst max|SMD| %.3f. Keep Figure 4 strictly secondary."
    ),
    worst_tail
  )
}

shape_alert <- if (is.finite(new_pattern) && new_pattern < 0) {
  sprintf(
    "Directional pattern preserved: mean RD at 20%%-50%% is %.2f pp more negative than at 3%%-5%%.",
    abs(new_pattern)
  )
} else {
  "SCORE ALERT: corrected curve no longer preserves the low-to-higher-risk directional pattern."
}

shift_alert <- if (low_shift > 2) {
  sprintf(
    "SCORE ALERT: low-risk continuous point shifted >2 pp; maximum %.2f pp.",
    low_shift
  )
} else {
  sprintf(
    "No >2-pp low-risk continuous point shift; maximum %.2f pp.",
    low_shift
  )
}

## ---- save data ---------------------------------------------------------------
OUT_REPORT <- file.path(RESULT_DIR,"REDTEAM_04_Corrected_Continuous_Report_B1000.txt")
OUT_CURVE <- file.path(RESULT_DIR,"REDTEAM_04_Corrected_Continuous_Curve_B1000.csv")
OUT_BOOT <- file.path(RESULT_DIR,"REDTEAM_04_Corrected_Continuous_Bootstrap_B1000.csv")
OUT_BAL <- file.path(RESULT_DIR,"REDTEAM_04_Corrected_Continuous_Balance.csv")

write.csv(CURVE,OUT_CURVE,row.names=FALSE)
write.csv(data.frame(replicate=seq_len(B),rbf,check.names=FALSE),
          OUT_BOOT,row.names=FALSE)
write.csv(BALANCE,OUT_BAL,row.names=FALSE)

## ---- corrected Figure 4 ------------------------------------------------------
ymin <- floor(min(CURVE$lo_pp,na.rm=TRUE)/5)*5
ymax <- ceiling(max(CURVE$hi_pp,na.rm=TRUE)/5)*5
if (ymin==ymax) { ymin<-ymin-5; ymax<-ymax+5 }

p <- ggplot2::ggplot(CURVE,ggplot2::aes(risk_pct,rd_pp))+
  ggplot2::geom_ribbon(ggplot2::aes(ymin=lo_pp,ymax=hi_pp),
                       fill="grey85",linewidth=0)+
  ggplot2::geom_hline(yintercept=0,linewidth=.65)+
  ggplot2::geom_vline(xintercept=P10*100,linetype=3,linewidth=.5)+
  ggplot2::geom_vline(xintercept=P90*100,linetype=3,linewidth=.5)+
  ggplot2::geom_line(linewidth=1)+
  ggplot2::geom_point(size=2.5)+
  ggplot2::annotate("text",x=P10*100,y=ymax-.4,
                    label=sprintf("P10 = %.1f%%",P10*100),
                    hjust=1.08,vjust=1,size=3)+
  ggplot2::annotate("text",x=P90*100,y=ymax-.4,
                    label=sprintf("P90 = %.1f%%",P90*100),
                    hjust=-.08,vjust=1,size=3)+
  ggplot2::annotate("label",x=max(CURVE$risk_pct)-1,y=ymax-.4,
                    label=if(is.finite(wald$p))
                      sprintf("Global interaction P = %.3f",wald$p)
                    else "Global interaction P = NA",
                    hjust=1,vjust=1,size=3,fill="white",label.size=.25)+
  ggplot2::scale_x_continuous(breaks=seq(0,50,10),
                              limits=range(CURVE$risk_pct),
                              expand=ggplot2::expansion(mult=c(.01,.01)))+
  ggplot2::scale_y_continuous(breaks=seq(ymin,ymax,5),
                              limits=c(ymin,ymax),
                              expand=ggplot2::expansion(mult=c(.01,.01)))+
  ggplot2::labs(
    title="Continuous association between predicted nodal risk and 5-year mortality difference",
    x="Predicted nodal metastasis risk, %",
    y=paste0("Absolute 5-year cancer-specific mortality risk difference, pp\n",
             "(oncologic colectomy - limited resection)")
  )+
  ggplot2::theme_classic(base_size=11)+
  ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=12),
                 plot.margin=ggplot2::margin(8,8,8,8))

PNG_FILE <- file.path(FIG_DIR,"Figure4_Continuous_RD_Curve_CORRECTED.png")
TIFF_FILE <- file.path(FIG_DIR,"Figure4_Continuous_RD_Curve_CORRECTED.tiff")
ggplot2::ggsave(PNG_FILE,p,width=8.2,height=5.8,units="in",dpi=600,bg="white")
ggplot2::ggsave(TIFF_FILE,p,width=8.2,height=5.8,units="in",dpi=600,
                device="tiff",compression="lzw",bg="white")

## ---- report ------------------------------------------------------------------
report <- c(
  "======================================================================",
  "A6 RED-TEAM — CORRECTED CONTINUOUS FIGURE 4",
  "======================================================================",
  paste0("Frozen cohort N = ",nrow(RAW)),
  paste0("Removed PS term(s): ",paste(stage_terms,collapse=" | ")),
  "Corrected weighting: no-Summary-Stage PS refitted separately within predicted-risk quintiles.",
  "Spline/IPCW estimator otherwise preserves the frozen continuous analysis.",
  paste0("Effective curve bootstrap B = ",eff_curve," / ",B),
  paste0("Effective interaction bootstrap B = ",eff_int," / ",B),
  "",
  "== CORRECTED CONTINUOUS RD(r) =="
)

for (i in seq_len(nrow(CURVE))) {
  z<-CURVE[i,]
  report<-c(report,sprintf(
    "risk=%5.1f%% | RD %+.2f pp (95%% CI %+.2f to %+.2f) | old %+.2f | change %+.2f pp",
    z$risk_pct,z$rd_pp,z$lo_pp,z$hi_pp,z$old_rd_pp,z$point_change_vs_old_pp))
}

report<-c(
  report,"","== P10 / P90 CONTRAST ==",
  sprintf("P10 = %.2f%% | P90 = %.2f%%",P10*100,P90*100),
  sprintf("Delta = RD(P90)-RD(P10): %+.2f pp (95%% CI %+.2f to %+.2f) | fraction Delta<0 %.3f",
          delta_point,delta_lo,delta_hi,delta_frac),
  sprintf("Old delta = %+.2f pp | corrected-minus-old change = %+.2f pp",
          OLD_DELTA_PP,delta_point-OLD_DELTA_PP),
  "","== GLOBAL TREATMENT x CONTINUOUS-RISK INTERACTION ==",
  sprintf("W = %.3f | df = %d | P = %.4f",wald$W,wald$df,wald$p),
  sprintf("Old interaction P = %.4f",OLD_INTERACTION_P),
  "","== LOCALIZED BALANCE AUDIT =="
)

for (i in seq_len(nrow(BALANCE))) {
  z<-BALANCE[i,]
  report<-c(report,sprintf(
    "%-42s n=%d (%d/%d) | max|SMD| %.3f [%s] | %s",
    z$analysis,z$n,z$limited,z$colectomy,z$max_abs_weighted_SMD,
    z$worst_variable,ifelse(z$balance_pass,"PASS","FAIL")))
}

report<-c(
  report,"","== RED-TEAM / SCORE CHECK ==",
  balance_alert,tail_alert,shape_alert,shift_alert,
  sprintf("Maximum absolute grid-point change versus old Figure 4 = %.2f pp.",max_shift),
  "",
  "Interpretation rule:",
  "Figure 4 is secondary/descriptive and must not be used to select a new threshold.",
  "The prespecified bottom-decile analysis remains the primary low-risk decision anchor.",
  "Restore Figure 4 only if the corrected pattern is clinically interpretable and localized balance is acceptable.",
  "======================================================================"
)

writeLines(report,OUT_REPORT,useBytes=TRUE)
cat(paste(report,collapse="\n"),"\n")
cat("\nREPORT: ",OUT_REPORT,
    "\nFIGURE PNG: ",PNG_FILE,
    "\nFIGURE TIFF: ",TIFF_FILE,"\n",sep="")
