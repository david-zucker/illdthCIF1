### Cumulative Incidence Function Estimation in the Illness-Death Model Using All Disease Cases
### David Zucker and Malka Gorfine

#' Apply methods for estimating the cumulative incidence function for disease in an illness-death model.
#'
#' @param age_recr age at recruitment
#' @param age_diag age at diagnosis
#' @param age_death age at death
#' @param age_end_fu age at end of follow-up (death or censoring)
#' @param status_end status at end of follow-up
#' @param tgrd grid of ages at which the CIF will be computed
#' @param tgrd.cb grid of ages over which the simultaneous confidence band will be computed (subset of tgrd)
#' @param covpr desired confidence interval coverage probability
#' @param nresam number of bootstrap replications for confidence band
#' @param dtgrd_width grid width for age at diagnosis in conditional survival calculation
#' @param bwselect flag for whether to carry out bandwidth selection
#' @param bwvec vector of candidate bandwidths
#' @param nfld number of folds in crossvalidation procedure for selecting the bandwidth
#' @param redo_bw flag for redoing bandwidth selection in bootstrap
#' @return description
#' @export



### Code by David Zucker
### Version of 31 March 2026

### R code for for applying to a dataset methods for estimating the
### cumulative incidence function for disease in an illness-death model.

### Methods provided:
### 1. The method of Aalen and Johansen (1978, Scandinavian Journal of Statistics)
### 2. The method of Gorfine, Zucker, and Shoham (2025, Biometrics)
### 3. The method of Zucker and Gorfine (2026, arXiv)

### Method #3 includes all disease cases.

# The call to the main function is as follows:

#FUNCTION COMPUTE ALL ESTIMATORS WITH CONFIDENCE INTERVALS AND BANDS
# cifcmp.full(age_recr, age_diag, age_death, age_end_fu, status_end, tgrd, tgrd.cb,
#   covpr, nresam, dtgrd_width, bwselect, bwvec, nfld, redo_bw) {

#The arguments are as follows:

# age_recr = age at recruitment
# age_diag = age at diagnosis
# age_death = age at death
# age_end_fu = age at end of follow-up (death or censoring)
  #note that occurrence of disease does not end follow-up
# status_end = status at end of follow-up
  #0 = alive without disease (i.e. censored)
  #1 = died without disease
  #2 = alive with disease
  #3 = died with disease
# tgrd = grid of ages at which the CIF will be computed
# tgrd.cb = grid of ages over which the simultaneous confidence band will be computed (subset of tgrd)
# covpr = desired confidence interval coverage probability
# nresam = number of bootstrap replications for confidence band
# dtgrd_width = grid width for age at diagnosis in conditional survival calculation
# bwselect = flag for whether to carry out bandwidth selection
# bwvec = vector of candidate bandwidths
# nfld = number of folds in crossvalidation procedure for selecting the bandwidth
# redo_bw = flag for redoing bandwidth selection in bootstrap

# when no bandwidth selection is done, the bandwidth used is bwvec[1]

library(survival)
didl = 1e-8

### FUNCTIONS ##################################################################

#COMPLEMENTARY LOG-LOG TRANSFORM AND INVERSE
trans = function(u) {
  ans = log(-log(1-u)+didl)
  return(ans)
}
trans.inv = function(v) {
  ans = 1-exp(-(exp(v)-didl))
  return(ans)
}

#TRIWEIGHT KERNEL (including boundary adjustment)
krnl = function (u,cc) {
  kk = 1.09375*(abs(u) <= 1)*(((1-u^2))^3)
  mu0 = 1.09375*((16/35)+cc-(cc^3)+0.6*(cc^5)-(1/7)*(cc^7))
  mu1 = 1.09375*(0.125-0.5*(cc^2)+0.75*(cc^4)-0.5*(cc^6)+(1/8)*(cc^8))
  mu2 = 1.09375*((16/315)+(1/3)*(cc^3)-0.6*(cc^5)+(3/7)*(cc^7)-(1/9)*(cc^9))
  numer = mu2 - u*mu1
  denom = mu0*mu2 - (mu1^2)
  ans = (numer/denom)*kk
  return(ans)
}

#KAPLAN-MEIER FUNCTION
KM = function(trecr,tfu,del,wt) {

  #IDENTIFY UNIQUE EVENT TIMES
  evtim = tfu[which(del==1)]
  time = unique(evtim)
  time = sort(time)
  ndist = length(time)

  #CALCULATIONS
  atrskvec = rep(0,ndist)
  km_surv = rep(0,ndist)
  skmcur = 1
  for (m in 1:ndist) {
    cur.tim = time[m]
    i_ev = which((trecr < cur.tim) & (tfu==cur.tim) & (del==1))
    i_atrsk = which((trecr < cur.tim) & (tfu >= cur.tim))
    d = sum(wt[i_ev])
    atrsk = sum(wt[i_atrsk])
    p = d/atrsk
    skmcur = skmcur * (1-p)
    atrskvec[m] = atrsk
    km_surv[m] = skmcur
  }
  ans = list(time=time, surv=km_surv, atrsk=atrskvec)

  return(ans)

}
#END KAPLAN-MEIER FUNCTION

#BERAN-TYPE CONDITIONAL SURVIVAL FUNCTION ESTIMATE
cndsrv = function(tstart,tfu,del,xsam,xval,xvmin,tval,bw,btwts) {

  #IDENTIFY UNIQUE EVENT TIMES
  evtim = tfu[which(del==1)]
  time = unique(evtim)
  time = sort(time)
  time = time[which(time >= xval)]
  ndist = length(time)
  time = c(time,Inf)

  #CALCULATIONS
  atrskvec = rep(0,ndist)
  km_surv = rep(1,ndist)
  surv1 = rep(1,length(tval))
  chcur = 0
  skmcur = 1
  zz = xsam - xval
  cc = pmin((xval-xvmin)/bw,1)
  wt = btwts*krnl(zz/bw,cc)
  for (m in 1:ndist) {
    cur.tim = time[m]
    i_atrsk = which((tstart < cur.tim) & (tfu >= cur.tim))
    s.ky = sum(wt[i_atrsk])
    if (s.ky > 0) {
      i_ev = which((tstart < cur.tim) & (tfu==cur.tim) & (del==1))
      s.kdn = sum(wt[i_ev])
      p = s.kdn/s.ky
      chcur = chcur + p
      skmcur = exp(-chcur)
    }
    atrskvec[m] = s.ky
    km_surv[m] = skmcur
    ixt = which((tval >= time[m]) & (tval < time[m+1]))
    surv1[ixt] = skmcur
  }

  time = time[1:ndist]
  ans = list(time=time, surv0=km_surv, atrsk=atrskvec, surv1=surv1)

  return(ans)

}
#END CONDITIONAL SURVIVAL FUNCTION

#FUNCTION FOR CIF COMPUTATION FOR AALEN-JOHANSEN ESTIMATOR
cifcmp.aj = function(age_recr, age_diag, age_death, age_end_fu, status_end, tgrd, wts) {

  #status_end
  #0 = alive without disease (i.e. censored)
  #1 = died without disease
  #2 = alive with disease
  #3 = died with disease

  #tgrd = grid of timepoints over which the CIF will be computed
  ntgrd = length(tgrd)

  #sample size
  nsam = length(age_recr)

  #SETUPS
  noncens = (status_end != 0)
  ix.cen = which(status_end == 0)
  diseased = (status_end >= 2)
  prev = ((age_diag <= age_recr) & diseased)
  ix.prev = which(prev)
  case = (age_diag <= age_death) & (!prev) & (noncens)
  ix.case = which(case)
  ncase = length(ix.case)
  age_diag[ix.cen] = Inf
  age_death[ix.cen] = Inf
  cifest = rep(NA,ntgrd)
  age_first = pmin(age_diag, age_death, age_end_fu)

  #REMOVE PREVALENT CASES FOR COMPUTATION OF KM FOR FIRST TRANSITION
  nprev = length(ix.prev)
  n.aj = nsam - length(ix.prev)
  age_recr.aj = age_recr
  age_diag.aj = age_diag
  age_death.aj = age_death
  age_end_fu.aj = age_end_fu
  status_end.aj = status_end
  age_first.aj = age_first
  wts.aj = wts
  noncens.aj = noncens
  if (nprev > 0) {
    age_recr.aj = age_recr[-ix.prev]
    age_diag.aj = age_diag[-ix.prev]
    age_death.aj = age_death[-ix.prev]
    age_end_fu.aj = age_end_fu[-ix.prev]
    status_end.aj = status_end[-ix.prev]
    age_first.aj = age_first[-ix.prev]
    wts.aj = wts[-ix.prev]
    noncens.aj = noncens[-ix.prev]
  }

  #KAPLAN-MEIER ESTIMATE FOR TIME TO FIRST TRANSITION
  #OMITTING PREVALENT CASES
  km_fst = KM(age_recr.aj, age_first.aj, noncens.aj, wts.aj)
  km_fst_times = km_fst$time
  km_fst_surv = km_fst$surv
  km_fst_surv1 = c(1,km_fst_surv)
  km_fst_at_rsk = km_fst$atrsk
  ycal.circ = km_fst_at_rsk/nsam
  nfst = length(km_fst_times)
  afun = km_fst_surv1[1:nfst]/ycal.circ

  #PRELIMINARIES
  at_rsk_i = matrix(0,nsam,nfst)
  ixf_vec = rep(0,nsam)
  for (ix.f in 1:nfst) {
    icur = which((age_first==km_fst_times[ix.f]) & (!prev) & (noncens))
    ixf_vec[icur] = ix.f
    icur1 = which((age_recr < km_fst_times[ix.f]) &
                    (age_first >= km_fst_times[ix.f]))
    at_rsk_i[icur1,ix.f] = 1
  }

  #CIF COMPUTATION
  x.main = matrix(0,ntgrd,nsam)
  for (ix.g in 1:ntgrd) {
    x.main[ix.g,ix.case] = afun[ixf_vec[ix.case]] * (age_diag[ix.case] <= tgrd[ix.g])
    cifest[ix.g] = sum(wts[ix.case]*x.main[ix.g,ix.case])/nsam
    cifest[ix.g] = min(c(cifest[ix.g],1))
  }

  #RETURN RESULT
  ans = list(cifest=cifest)
  return(ans)

}
#END OF FUNCTION TO COMPUTE AJ ESTIMATOR

#FUNCTION FOR CIF COMPUTATION FOR GZS ESTIMATOR
cifcmp.gzs = function(age_recr, age_diag, age_death, age_end_fu,
  status_end, tgrd, wts) {

  #status_end
  #0 = alive without disease (i.e. censored)
  #1 = died without disease
  #2 = alive with disease
  #3 = died with disease

  #sample size
  nsam = length(age_recr)

  #tgrd = grid of timepoints over which the CIF will be computed
  ntgrd = length(tgrd)

  #SETUPS
  ix.cen = which(status_end == 0)
  ix.noncen = which(status_end > 0)
  n_noncen = length(ix.noncen)
  died = as.numeric((status_end == 1) | (status_end == 3))
  ix.died = which(died==1)
  alive = 1 - died
  ix.alive = which(alive==1)
  n_alive = length(ix.alive)
  #disease cases
  case = (status_end >= 2)
  ix.case = which(case)
  ncase = length(ix.case)
  #disease cases who died
  dcase = (status_end == 3)
  ix.dcase = which(dcase)
  ndcase = length(ix.dcase)
  #additional setups
  age_diag[ix.cen] = Inf
  age_death[ix.cen] = Inf

  #KAPLAN-MEIER ESTIMATE OF DEATH TIME DISTN
  km_dth = survfit(Surv(age_recr, age_end_fu, died) ~ 1, timefix=FALSE, weights=wts)
  km_dth_times = km_dth$time
  km_dth_surv = km_dth$surv
  ndth = length(km_dth_times)

  #SURVIVAL CALCULATIONS NEEDED FOR ESTIMATE BASED ON DECEASED CASES ONLY
  km_dth_surv1 = c(1, km_dth_surv)
  intwts.numer = -diff(km_dth_surv1)

  #K VECTOR FOR ANALYSIS BASED ON DECEASED CASES ONLY
  Kvec1 = rep(0,nsam)
  ixd_vec = rep(0,nsam)
  for (ix.d in 1:ndth) {
    icur = which(age_death == km_dth_times[ix.d])
    ixd_vec[icur] = ix.d
    ndth.cur.nrm = sum(wts[icur])/nsam
    Kvec1[icur] = intwts.numer[ix.d] / ndth.cur.nrm
  }

  #CIF COMPUTATION
  x.main = matrix(0,ntgrd,nsam)
  cifest = rep(NA,ntgrd)
  for (ix.g in 1:ntgrd) {
    x.main[ix.g,ix.dcase] = Kvec1[ix.dcase] * (age_diag[ix.dcase] <= tgrd[ix.g])
    cifest[ix.g] = sum(wts[ix.dcase]*x.main[ix.g,ix.dcase])/nsam
    cifest[ix.g] = min(c(cifest[ix.g],1))
  }

  #RETURN RESULT
  ans = list(cifest=cifest)
  return(ans)

}
# END OF FUNCTION TO COMPUTE GZS ESTIMATOR

#FUNCTION TO COMPUTE CROSSVALIDATION CRITERION FOR A GIVEN BANDWIDTH
cv_err = function(bwcur,tstart,tfu,del,xsam,dtgrd1,xvmin,age_grd,grp,nfld) {

  ndisgrd1 = length(dtgrd1)
  n_age_grd = length(age_grd)
  gof_fit = 0

  for (k in 1:nfld) {

    #preliminaries
    s2mat_cv = matrix(1,ndisgrd1,n_age_grd)
    ixdat = which(grp != k)
    ixtst = which(grp == k)

    #set up test data
    ntst = length(ixtst)
    xsam1 = xsam[ixtst]
    tstart1 = tstart[ixtst]
    tfu1 = tfu[ixtst]
    del1 = del[ixtst]
    ordr = order(xsam1)
    xsam1 = xsam1[ordr]
    tstart1 = tstart1[ordr]
    tfu1 = tfu1[ordr]
    del1 = del1[ordr]
    xsam1cut = cut(xsam1, breaks=3, labels=FALSE)

    #estimate conditional survival over a grid
    for (ixdt in 1:ndisgrd1) {
      v1icur = dtgrd1[ixdt]
      cscur = cndsrv(tstart[ixdat],tfu[ixdat],del[ixdat],xsam[ixdat],xval=v1icur,
        xvmin,tval=age_grd,bw=bwcur,btwts=rep(1,length(ixdat)))
      s2mat_cv[ixdt,] = pmax(cscur$surv1,didl)
    }

    #evaluate on test data
    prd = rep(0,ntst)
    ch2vec = rep(0,ntst)
    for (itst in 1:ntst) {
      xvalcur = xsam1[itst]
      tvalcur1 = tstart1[itst]
      tvalcur2 = tfu1[itst]
      xvgrd = findInterval(xvalcur,dtgrd1)
      tvgrd1 = findInterval(tvalcur1,age_grd)
      tvgrd2 = findInterval(tvalcur2,age_grd)
      srv1cur = s2mat_cv[xvgrd,tvgrd1]
      srv2cur = s2mat_cv[xvgrd,tvgrd2]
      ch1 = -log(srv1cur)
      ch2 = -log(srv2cur)
      prd[itst] = ch2 - ch1
      ch2vec[itst] = ch2
    }

    #wrap-up
    ixcv = which(prd>0)
    gof_fit = gof_fit + sum(((del1[ixcv]-prd[ixcv])^2)/prd[ixcv])

  }

  return(gof_fit)

}
#END OF CROSSVALIDATION FUNCTION

#FUNCTION FOR CIF COMPUTATION FOR NEW ESTIMATOR
cifcmp.new = function(age_recr, age_diag, age_death, age_end_fu,
  status_end, tgrd, bwselect, bwvec, nfld, dtgrd_width, wts) {

  #status_end
  #0 = alive without disease (i.e. censored)
  #1 = died without disease
  #2 = alive with disease
  #3 = died with disease

  #sample size
  nsam = length(age_recr)

  #tgrd = grid of timepoints over which the CIF will be computed
  ntgrd = length(tgrd)

  #SETUPS
  ix.cen = which(status_end == 0)
  ix.noncen = which(status_end > 0)
  n_noncen = length(ix.noncen)
  died = as.numeric((status_end == 1) | (status_end == 3))
  ix.died = which(died==1)
  alive = 1 - died
  ix.alive = which(alive==1)
  n_alive = length(ix.alive)
  #disease cases
  case = (status_end >= 2)
  ix.case = which(case)
  ncase = length(ix.case)
  #disease cases who died
  dcase = (status_end == 3)
  ix.dcase = which(dcase)
  ndcase = length(ix.dcase)
  #additional setups
  age_diag[ix.cen] = Inf
  age_death[ix.cen] = Inf
  dtgrd_width1 = 2*dtgrd_width
  rtgrd_width = dtgrd_width
  age_grd_w_cv = 2*dtgrd_width

  #KAPLAN-MEIER ESTIMATE OF DEATH TIME DISTN
  km_dth = survfit(Surv(age_recr, age_end_fu, died) ~ 1, timefix=FALSE, weights=wts)
  km_dth_times = km_dth$time
  km_dth_surv = km_dth$surv
  ndth = length(km_dth_times)

  #SURVIVAL CURVE FOR TIME FROM RECRUITMENT TO CENSORING
  fu_time = age_end_fu - age_recr
  km_cens = survfit(Surv(fu_time, alive) ~ 1, timefix=FALSE, weights=wts)

  #DATA SETUP FOR CONDITIONAL SURVIVAL CALCULATION
  tstart = pmax(age_recr[ix.case],age_diag[ix.case])
  tfu = age_end_fu[ix.case]
  del = died[ix.case]
  xsam = age_diag[ix.case]

  #SET UP GRID OF DIAGNOSIS AGE VALUES AND RECRUITMENT AGE VALUES
  dtlo = min(xsam)
  dthi = max(xsam)
  dtgrd = seq(dtlo,dthi,dtgrd_width)
  dtgrd1 = seq(dtlo,dthi,dtgrd_width1)
  ndisgrd = length(dtgrd)
  L_R = min(age_recr)
  U_R = max(age_recr)
  rtgrd = seq(L_R,U_R,rtgrd_width) - didl
  nrecgrd = length(rtgrd)

  #CONDITIONAL SURVIVAL CALCULATION OVER A GRID OF DIAGNOSIS AGE VALUES
  #AND RECRUITMENT AGE VARIABLES
  s2mat = matrix(1,ndisgrd,nrecgrd)

  #WITHOUT BANDWIDTH SELECTION
  if (!bwselect) {
    bwval = bwvec[1]
    for (ixdt in 1:ndisgrd) {
      v1icur = dtgrd[ixdt]
      if (v1icur > U_R) break
      incl = which(rtgrd > v1icur)
      cscur = cndsrv(tstart,tfu,del,xsam,xval=v1icur,xvmin=dtlo,
        tval=rtgrd[incl],bw=bwval,btwts=wts[ix.case])
      s2mat[ixdt,incl] = cscur$surv1
    }
  }

  #WITH BANDWIDTH SELECTION
  if (bwselect) {

    #GRID OF AGE VALUES
    age_lo = floor(min(c(age_diag,age_recr)))
    age_hi = floor(max(age_end_fu))
    age_grd = seq(age_lo,age_hi,age_grd_w_cv)
    n_age_grd = length(age_grd)

    #SPLIT SAMPLE INTO BLOCKS
    nblk = floor(ncase/nfld)
    rmdr = ncase - nfld*nblk
    grp = NULL
    for (k in 1:nfld) {
      if (k <= rmdr) {grp = c(grp,rep(k,nblk+1))}
      else {grp = c(grp,rep(k,nblk))}
    }
    prm = sample.int(ncase,ncase)
    grp = grp[prm]

    #IDENTIFY OPTIMAL BANDWIDTH
    nbw = length(bwvec)
    cv_err_vec = rep(NA,nbw)
    for (ii in 1:nbw) {
      bwcur = bwvec[ii]
      cv_err_vec[ii] = cv_err(bwcur,tstart,tfu,del,xsam,dtgrd1,dtlo,age_grd,grp,nfld)
    }
    cvmin = min(cv_err_vec, na.rm=TRUE)
    if (cvmin >= .Machine$double.xmax) {
      bw_opt = max(bwvec)
    }
    else {
      ixbw = which(cv_err_vec==cvmin)
      if (length(ixbw>1)) ixbw = ixbw[1]
      bw_opt = bwvec[ixbw]
    }

    #FINAL CONDITIONAL SURVIVAL CURVE CALCULATION
    for (ixdt in 1:ndisgrd) {
      v1icur = dtgrd[ixdt]
      if (v1icur > U_R) break
      incl = which(rtgrd > v1icur)
      cscur = cndsrv(tstart,tfu,del,xsam,xval=v1icur,xvmin=dtlo,
        tval=rtgrd[incl],bw=bw_opt,btwts=wts[ix.case])
      s2mat[ixdt,incl] = cscur$surv1
    }
    bwval = bw_opt

  }

  #SURVIVAL CURVE FOR CENSORING OVER A GRID OF DIAGNOSIS AGE VALUES
  #AND RECRUITMENT AGE VALUES
  cens_fac_mat = matrix(1,ndisgrd,nrecgrd)
  for (ixdt in 1:ndisgrd) {
    v1icur = dtgrd[ixdt]
    fu_time_cur = v1icur - rtgrd
    incl = which(fu_time_cur > 0)
    fu_time_cur.a = rev(fu_time_cur[incl]) - didl
    if (length(incl) > 0) {
      cf_obj = summary(km_cens, times=fu_time_cur.a, extend=TRUE)
      cens_fac_mat[ixdt,incl] = rev(cf_obj[[6]])
    }
  }

  #K VECTOR FOR ANALYSIS INCLUDING ALL DISEASE CASES
  Kvec2 = rep(0,nsam)
  st2r = rep(1,nsam)
  rtime = age_recr - didl
  oo = order(rtime)
  dth_obj_r = summary(km_dth, times=rtime[oo], extend=TRUE)
  st2r[oo] = dth_obj_r[[6]]
  ixrcr = findInterval(age_recr,rtgrd)
  ixdiag = findInterval(age_diag,dtgrd)
  ixrcr[ixrcr == 0] = 1
  ixrcr[ixrcr > nrecgrd] = nrecgrd

  for (ic in 1:ncase) {
    i = ix.case[ic]
    ixv1i = ixdiag[i]
    censfac_cur = cens_fac_mat[ixv1i,ixrcr]
    K2fac = sum(wts*s2mat[ixv1i,ixrcr]*censfac_cur/st2r)/nsam
    Kvec2[i] = 1/K2fac
  }

  #CIF COMPUTATION
  x.main = matrix(0,ntgrd,nsam)
  cifest = rep(NA,ntgrd)
  for (ix.g in 1:ntgrd) {
    x.main[ix.g,ix.case] = Kvec2[ix.case] * (age_diag[ix.case] <= tgrd[ix.g])
    cifest[ix.g] = sum(wts[ix.case]*x.main[ix.g,ix.case])/nsam
    cifest[ix.g] = min(c(cifest[ix.g],1))
  }

  #RETURN RESULT
  ans = list(cifest=cifest, bwval=bwval)
  return(ans)

}
# END OF FUNCTION TO COMPUTE NEW ESTIMATOR

#FUNCTION FOR BOOTSTRAP CI'S
boot.ci = function(nresam, orig.est, boot.est, covpr, ixcb) {

  #SETUPS
  ciqpr_lo = (1-covpr)/2
  ciqpr_hi = 1 - ciqpr_lo
  zcrit = qnorm(covpr)
  ci_adj_fac1 = zcrit*sqrt(ciqpr_lo*ciqpr_hi/nresam)
  ciqpr_lo_mod = ciqpr_lo - ci_adj_fac1
  ciqpr_lo_mod = floor(nresam*ciqpr_lo_mod)/nresam
  ciqpr_hi_mod = ciqpr_hi + ci_adj_fac1
  ciqpr_hi_mod = ceiling(nresam*ciqpr_hi_mod)/nresam
  ci_adj_fac2 = zcrit*sqrt(covpr*(1-covpr)/nresam)
  covpr_mod = covpr + ci_adj_fac2
  covpr_mod = ceiling(nresam*covpr_mod)/nresam
  orig.est.mat = matrix(rep(orig.est,nresam), nrow=nresam, byrow=TRUE)
  orig.est.tr = trans(orig.est)
  orig.est.tr.mat = matrix(rep(orig.est.tr,nresam), nrow=nresam, byrow=TRUE)
  boot.est.tr = trans(boot.est)
  boot.stat = boot.est - orig.est.mat
  boot.stat.tr = boot.est.tr - orig.est.tr.mat

  #POINTWISE CONFIDENCE INTERVALS
  qntlo = apply(boot.stat.tr, 2, quantile, probs=ciqpr_hi_mod, type=4, na.rm=TRUE)
  qnthi = apply(boot.stat.tr, 2, quantile, probs=ciqpr_lo_mod, type=4, na.rm=TRUE)
  ptwise.ci.lo = trans.inv(orig.est.tr - qntlo)
  ptwise.ci.hi = trans.inv(orig.est.tr - qnthi)
  ptwise.ci.width = ptwise.ci.hi - ptwise.ci.lo

  #SIMULTANEOUS CONFIDENCE BAND
  boot.stat.max = apply(abs(boot.stat[, ixcb, drop=FALSE]), 1, max, na.rm=TRUE)
  qmax = quantile(boot.stat.max, probs=covpr_mod, type=4, na.rm=TRUE)
  band.lo = orig.est[ixcb] - qmax
  band.hi = orig.est[ixcb] + qmax
  band.width = band.hi - band.lo

  #ESTIMATE WITH BOOTSTRAP BIAS CORRECTION
  boot.mean = apply(boot.est, 2, mean, na.rm=TRUE)
  bbc.est = 2*orig.est - boot.mean

  ans = list(
    bbc.est = bbc.est,
    ptwise.ci.lo = ptwise.ci.lo,
    ptwise.ci.hi = ptwise.ci.hi,
    ptwise.ci.width = ptwise.ci.width,
    band.lo = band.lo,
    band.hi = band.hi,
    band.width = band.width)
  return(ans)

}
#END OF FUNCTION FOR BOOTSTRAP CI'S

#FUNCTION COMPUTE ALL ESTIMATORS WITH CONFIDENCE INTERVALS AND BANDS
cifcmp.full = function(age_recr, age_diag, age_death, age_end_fu, status_end,
  tgrd, tgrd.cb, covpr, nresam, dtgrd_width, bwselect, bwvec, nfld, redo_bw) {

  #PRELIMINARIES
  n = length(age_recr)
  ntgrd = length(tgrd)
  ixcb = which(tgrd %in% tgrd.cb)
  verbose = TRUE

  #BOOTSTRAP SAMPLES
  b_ix = NULL
  for (ib in 1:nresam) {
    ixcur = sample.int(n,n,replace=TRUE)
    b_ix = rbind(b_ix,ixcur)
  }

  #SET UP ARRAYS FOR BOOTSTRAP
  aj.boot.est = matrix(NA,nresam,ntgrd)
  gzs.boot.est = matrix(NA,nresam,ntgrd)
  new.boot.est = matrix(NA,nresam,ntgrd)

  #AALEN-JOHANSEN ESTMATOR
  if (verbose) print(noquote('Computing AJ estimator ...'))
  cifest.aj = cifcmp.aj(age_recr, age_diag, age_death, age_end_fu,
    status_end, tgrd, rep(1,n))
  aj.est = cifest.aj$cifest
  for (ib in 1:nresam) {
    if (verbose) print(noquote(paste0('Bootstrap Replication ', ib)))
    ixc = b_ix[ib,]
    cifest.boot.aj = cifcmp.aj(age_recr[ixc], age_diag[ixc], age_death[ixc],
      age_end_fu[ixc], status_end[ixc], tgrd, rep(1,n))
    aj.boot.est[ib,] = cifest.boot.aj$cifest
  }
  aj.ci.rslts = boot.ci(nresam, aj.est, aj.boot.est, covpr,ixcb)
  aj.est.bbc = aj.ci.rslts$bbc.est
  aj.ptwise.ci.lo = aj.ci.rslts$ptwise.ci.lo
  aj.ptwise.ci.hi = aj.ci.rslts$ptwise.ci.hi
  aj.ptwise.ci.width = aj.ci.rslts$ptwise.ci.width
  aj.band.lo = aj.ci.rslts$band.lo
  aj.band.hi = aj.ci.rslts$band.hi
  aj.band.width = aj.ci.rslts$band.width

  #GZS ESTIMATOR
  if (verbose) print(noquote('Computing GZS estimator ...'))
  cifest.gzs = cifcmp.gzs(age_recr, age_diag, age_death, age_end_fu,
    status_end, tgrd, rep(1,n))
  gzs.est = cifest.gzs$cifest
  for (ib in 1:nresam) {
    if (verbose) print(noquote(paste0('Bootstrap Replication ', ib)))
    ixc = b_ix[ib,]
    cifest.boot.gzs = cifcmp.gzs(age_recr[ixc], age_diag[ixc], age_death[ixc],
      age_end_fu[ixc], status_end[ixc], tgrd, rep(1,n))
    gzs.boot.est[ib,] = cifest.boot.gzs$cifest
  }
  gzs.ci.rslts = boot.ci(nresam, gzs.est, gzs.boot.est, covpr, ixcb)
  gzs.est.bbc = gzs.ci.rslts$bbc.est
  gzs.ptwise.ci.lo = gzs.ci.rslts$ptwise.ci.lo
  gzs.ptwise.ci.hi = gzs.ci.rslts$ptwise.ci.hi
  gzs.ptwise.ci.width = gzs.ci.rslts$ptwise.ci.width
  gzs.band.lo = gzs.ci.rslts$band.lo
  gzs.band.hi = gzs.ci.rslts$band.hi
  gzs.band.width = gzs.ci.rslts$band.width

  #NEW ESTIMATOR
  if (verbose) print(noquote('Computing new estimator ...'))
  cifest.new = cifcmp.new(age_recr, age_diag, age_death, age_end_fu,
    status_end, tgrd, bwselect, bwvec, nfld, dtgrd_width, rep(1,n))
  new.est = cifest.new$cifest
  new.bw = cifest.new$bwval
  if (redo_bw) {
    bwvboot = bwvec
  }
  else {
    bwvboot = new.bw
  }
  for (ib in 1:nresam) {
    if (verbose) print(noquote(paste0('Bootstrap Replication ', ib)))
    ixc = b_ix[ib,]
    cifest.boot.new = cifcmp.new(age_recr[ixc], age_diag[ixc], age_death[ixc],
      age_end_fu[ixc], status_end[ixc], tgrd, redo_bw, bwvboot, nfld, dtgrd_width, rep(1,n))
    new.boot.est[ib,] = cifest.boot.new$cifest
  }
  new.ci.rslts = boot.ci(nresam, new.est, new.boot.est, covpr, ixcb)
  new.est.bbc = new.ci.rslts$bbc.est
  new.ptwise.ci.lo = new.ci.rslts$ptwise.ci.lo
  new.ptwise.ci.hi = new.ci.rslts$ptwise.ci.hi
  new.ptwise.ci.width = new.ci.rslts$ptwise.ci.width
  new.band.lo = new.ci.rslts$band.lo
  new.band.hi = new.ci.rslts$band.hi
  new.band.width = new.ci.rslts$band.width

  #BANDWIDTH FOR NEW ESTIMATOR
  bwval = cifest.new$bwval

  ans = list(
    aj.est = aj.est,
    aj.ptwise.ci.width = aj.ptwise.ci.width,
    aj.ptwise.ci.lo = aj.ptwise.ci.lo,
    aj.ptwise.ci.hi = aj.ptwise.ci.hi,
    gzs.est = gzs.est,
    gzs.ptwise.ci.width = gzs.ptwise.ci.width,
    gzs.ptwise.ci.lo = gzs.ptwise.ci.lo,
    gzs.ptwise.ci.hi = gzs.ptwise.ci.hi,
    new.est = new.est,
    new.ptwise.ci.width = new.ptwise.ci.width,
    new.ptwise.ci.lo = new.ptwise.ci.lo,
    new.ptwise.ci.hi = new.ptwise.ci.hi,
    aj.band.width = aj.band.width,
    aj.band.lo = aj.band.lo,
    aj.band.hi = aj.band.hi,
    gzs.band.width = gzs.band.width,
    gzs.band.lo = gzs.band.lo,
    gzs.band.hi = gzs.band.hi,
    new.band.width = new.band.width,
    new.band.lo = new.band.lo,
    new.band.hi = new.band.hi,
    aj.est.bbc = aj.est.bbc,
    gzs.est.bbc = gzs.est.bbc,
    new.est.bbc = new.est.bbc,
    bwval = bwval)
  #ALTOGETHER 25 OUTPUT ITEMS

  return(ans)

}
#END OF FUNCTION COMPUTE ALL ESTIMATORS WITH CONFIDENCE INTERVALS AND BANDS
