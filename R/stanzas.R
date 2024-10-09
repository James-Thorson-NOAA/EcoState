
#' @export
make_stanza_data <-
function( settings ){

  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  # Indexing
  s_s2 = match(names(settings$stanza_groups), settings$taxa)
  n_s2 = length(s_s2)
  g2_s2 = match( settings$stanza_groups, settings$unique_stanza_groups )
  if( n_s2>0 ){
    t2_s2 = sapply( seq_along(g2_s2), FUN=\(i)cumsum(g2_s2[i]==g2_s2)[i] )
  }else{
    t2_s2 = vector()
  }
  n_g2 = settings$n_g2

  # Checks for inputs for each taxa
  K_g2 = settings$K[settings$unique_stanza_groups]
  d_g2 = settings$d[settings$unique_stanza_groups]
  Wmat_g2 = settings$Wmat[settings$unique_stanza_groups]
  Amat_g2 = settings$Amat[settings$unique_stanza_groups]
  SpawnX_g2 = settings$SpawnX[settings$unique_stanza_groups]
  Amax_s2 = settings$Amax[settings$multigroup_taxa]
  Wmatslope_g2 = settings$Wmatslope[settings$unique_stanza_groups]
  #Leading_s2 = unlist(tapply( Amax_s2, FUN=\(v) v==max(v), INDEX=g2_s2))
  Leading_s2 = settings$Leading[settings$multigroup_taxa]
  plusage_g2 = tapply( Amax_s2, INDEX=g2_s2, FUN=max )

  # Variable-specific parameters
  stanzainfo_s2z = cbind( "s" = s_s2,
                          "s2" = seq_along(s_s2),
                          "g2" = g2_s2,
                          "lead" = Leading_s2,
                          "t2" = t2_s2,
                          "amax" = Amax_s2 )

  tmp = stanzainfo_s2z[which(stanzainfo_s2z[,'lead']==1),,drop=FALSE]
  stanzainfo_g2z = cbind( "Wmat" = Wmat_g2,
                          "Amat" = Amat_g2,
                          "Wmatslope" = Wmatslope_g2,
                          "SpawnX" = SpawnX_g2,
                          "plusage" = plusage_g2,
                          "K" = K_g2,
                          "d" = d_g2,
                          "lead_s" = tmp[match(seq_len(n_g2),tmp[,'g2']),'s'] )

  # EASIER TO LOOP THROUGH STANZAS
  #X_zz = NULL
  X_zz_g2 = NULL
  for( g2 in seq_len(n_g2) ){
    # Make fractional-age starting at age=0
    AGE = (seq_len(plusage_g2[g2] * settings$STEPS_PER_YEAR) - 1) / settings$STEPS_PER_YEAR
    stanzainfo_t2z = stanzainfo_s2z[which(stanzainfo_s2z[,'g2']==g2),,drop=FALSE]   # drop=FALSE in case only one stanza for single multigroup

    #
    t2_a = sapply( AGE, FUN=\(a){sum(a>=stanzainfo_t2z[,'amax'])} ) + 1
    t2_a = ifelse( t2_a > sum(stanzainfo_s2z[,'g2']==g2), sum(stanzainfo_s2z[,'g2']==g2), t2_a )

    # Stack
    Xg2_zz = cbind(
      g2 = g2,
      AGE = AGE,
      age_class = floor(AGE),
      t2 = t2_a,
      is_lead = stanzainfo_s2z[which(stanzainfo_s2z[,'g2']==g2),'lead'][t2_a]
    )
    #X_zz = rbind( X_zz, Xg2_zz )
    X_zz_g2[[g2]] = Xg2_zz
  }

  # Add and output
  stanza_data = list(
    n_s2 = n_s2,
    n_g2 = n_g2,
    stanzainfo_s2z = stanzainfo_s2z,
    stanzainfo_g2z = stanzainfo_g2z,
    #X_zz = X_zz
    X_zz_g2 = X_zz_g2
  )
  return( stanza_data )
}

#' @export
fecundity_by_weight <-
function( W,
          Wmat,
          Wmatslope ){

  # Globals
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  #pos = function(x){
  #  "c" <- ADoverload("c")
  #  "[<-" <- ADoverload("[<-")
  #  0.5*(abs(x)+x)
  #}

  if( is.infinite(Wmatslope) ){
    out = W - Wmat
    out = 0.5 * ( abs(out) + out )
  }else{
    out = plogis( (W - Wmat) * Wmatslope ) * W
  }
  return(out)
}

#' @export
add_stanza_params <-
function( p,
          stanza_data,
          settings ){

  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  n_s2 = stanza_data$n_s2
  n_g2 = stanza_data$n_g2
  Amat_g2 = stanza_data$stanzainfo_g2z[,'Amat']
  Wmatslope_g2 = stanza_data$stanzainfo_g2z[,'Wmatslope']
  d_g2 = plogis( p$logit_d_g2 )
  inv1minus_d_g2 = 1 / (1 - d_g2)
  #Wmat_g2 = exp(p$log_winf_z[1]) * stanza_data$stanzainfo_g2z[,'Wmat']
  plusage_g2 = stanza_data$stanzainfo_g2z[,'plusage']
  stanzainfo_s2z = stanza_data$stanzainfo_s2z
  stanzainfo_g2z = stanza_data$stanzainfo_g2z
  Z_s2 = exp(p$logPB_i)[stanzainfo_s2z[,'s']]
  #X_zz = stanza_data$X_zz
  X_zz_g2 = stanza_data$X_zz_g2

  # Replace Wmat with Amat if available
  which_replace = which(!is.na(Amat_g2))
  if( length(which_replace) > 0 ){
    # Given:
    #   dW = H * W^d - k * W
    # Then:
    #   W(A) = W_inf * (1 - exp(-k * (1-d) * A)) ^ (1/(1-d))
    #
    # Given:
    #   W(A) = a * L(A)^b where b=3
    # Then:
    #   L(A) = L_inf * (1 - exp(-K * (1-m) * A)) ^ (1/(1-m))
    # Where:
    #   k = b * K = 3K
    #   m = d*b + 1 - b = 3*d - 2
    #
    # Therefore:
    #   Wmat = (1 - exp(-3*K*(1-d) * Amat)) ^ (1/(1-d))
    k = 3 * exp(p$log_K_g2[which_replace])  # log_K_g2 is the log of K in length, where k = b*K and W = a * L^b, and we assume b=3
    p$Wmat_g2[which_replace] = (1 - exp(-k * (1 - d_g2[which_replace]) * Amat_g2[which_replace])) ^ inv1minus_d_g2[which_replace]
  }

  # Globals
  vbm_g2 = (1 - 3 * exp(p$log_K_g2) / settings$STEPS_PER_YEAR)
  ################## EXPERIMENT WITH RTMB
  # EASIER TO LOOP THROUGH STANZAS
  #Y_zz = matrix(nrow=0, ncol=4)
  Y_zz_g2 = NULL
  baseEggsStanza = baseSpawnBio = baseRzeroS = rep(0, n_g2)
  leading_ratio_s2 = Q_s2 = Consumption_s2 = rep(0, n_s2)
  for( g2 in seq_len(n_g2) ){
    #
    Xg2_zz = X_zz_g2[[g2]]

    # Make fractional-age starting at age=0
    AGE = Xg2_zz[,'AGE'] * Z_s2[1] / Z_s2[1]   # Extra stuff ensures that it is class-advector

    #
    stanzainfo_t2z = stanzainfo_s2z[which(stanzainfo_s2z[,'g2']==g2),,drop=FALSE]
    which_s2 = stanzainfo_t2z[,'s2']
    which_leading = which(stanzainfo_t2z[,'lead']==1)
    leading_s2 = stanzainfo_t2z[which_leading,'s2']

    #
    t2_a = Xg2_zz[,'t2']

    # WageS and QageS
    k = exp(p$log_K_g2[g2]) * 3     # 3 because W = L^3, i.e. converting VB length param K to weight param k
    WageS = (1 - exp(-k * (1 - d_g2[g2]) * AGE)) ^ inv1minus_d_g2[g2]
    #WageS = exp(p$log_winf_z[1]) * (1 - exp(-k * (1 - d) * (AGE))) ^ (1 / (1 - d))
    QageS = WageS ^ d_g2[g2]

    # SurvS
    Zrate = Z_s2[which_s2[t2_a]] / settings$STEPS_PER_YEAR
    SurvS = exp(-1 * c(0, cumsum(Zrate)[-length(Zrate)] ))

    # Correct plus-group
    SurvS[length(SurvS)] = SurvS[length(SurvS)] / (1 - exp(-Zrate[length(Zrate)]))

    # Solve for R0 and NageS
    BperR = sum(SurvS * WageS * Xg2_zz[,'is_lead'])
    baseRzeroS[g2] = exp(p$logB_i)[stanzainfo_g2z[g2,'lead_s']] / BperR
    NageS = SurvS * baseRzeroS[g2]

    #
    baseEggsStanza[g2] = sum(NageS * fecundity_by_weight(WageS, p$Wmat_g2[g2], Wmatslope_g2[g2]) )
    baseSpawnBio[g2] = sum(NageS * fecundity_by_weight(WageS, p$Wmat_g2[g2], Wmatslope_g2[g2]) )

    #
    s = stanzainfo_t2z[which_leading,'s']
    Consumption_s2[leading_s2] = exp(p$logQB_i[s]) * exp(p$logB_i[s])

    # Solve for expected consumption
    for( t2 in seq_len(nrow(stanzainfo_t2z)) ){
      s = stanzainfo_t2z[t2,'s']
      which_a = which( t2_a == stanzainfo_t2z[t2,'t2'] )
      p$logB_i[s] = log(sum(NageS[which_a] * WageS[which_a], na.rm=TRUE))
      Q_s2[stanzainfo_t2z[t2,'s2']] = sum(NageS[which_a] * QageS[which_a], na.rm=TRUE)
    }

    # Loop again for Cons and QB
    for( t2 in seq_len(nrow(stanzainfo_t2z)) ){
      s = stanzainfo_t2z[t2,'s']
      s2 = stanzainfo_t2z[t2,'s2']
      leading_ratio_s2[s2] = Consumption_s2[leading_s2] / Q_s2[leading_s2]
      Consumption_s2[s2] = Q_s2[s2] * leading_ratio_s2[s2]
      p$logQB_i[s] = log(Consumption_s2[s2]) - p$logB_i[s]
    }

    # SplitAlpha = REco.sim$stanzas$SplitAlpha
    SplitAlpha = rep(0, length=plusage_g2[g2] * settings$STEPS_PER_YEAR)
    for( t2 in seq_len(nrow(stanzainfo_t2z)) ){
      s = stanzainfo_t2z[t2,'s']
      s2 = stanzainfo_t2z[t2,'s2']
      which_a = which( t2_a == stanzainfo_t2z[t2,'t2'] )
      if(which_a[length(which_a)]==length(WageS)) which_a = which_a[-length(which_a)]
      # W(a+1) = vbm_g2 * W(a) + (Food/Pred) * SplitAlpha
      # SplitAlpha = (W(a+1) - vbm_g2 * W(a)) / (Food/Pred)
      SplitAlpha[which_a] = (WageS[which_a+1] - vbm_g2[g2] * WageS[which_a]) / leading_ratio_s2[s2]
      SplitAlpha[length(SplitAlpha)] = SplitAlpha[length(SplitAlpha)-1]
    }

    # Stack
    Yg2_zz = cbind( WageS=WageS, log_NageS=log(NageS), QageS=QageS, SplitAlpha=SplitAlpha)
    #Y_zz = rbind( Y_zz, Yg2_zz )  # deparse.level=0 avoids RTMB error
    Y_zz_g2[[g2]] = Yg2_zz
  }

  # Add and output
  p = c(p, list(
    #Y_zz = Y_zz,
    Y_zz_g2 = Y_zz_g2,
    baseEggs_g2 = baseEggsStanza,
    baseSB_g2 = baseSpawnBio,
    baseR0_g2 = baseRzeroS
  ))
  return(p)
}

#' @export
update_stanzas <-
function( p,
          stanza_data,
          FoodGain_s,
          LossPropToB_s,
          F_s,
          #increase_age = TRUE,
          STEPS_PER_YEAR = 1 ){

  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  # Globals
  vbm_g2 = (1 - 3 * stanza_data$stanzainfo_g2z[,'K'] / STEPS_PER_YEAR)       # 3*K = k, where k is from dW = H*W^d - k*W
  #SB_g2 = exp(p$logPB_i)[stanza_data$stanzainfo_s2z[,'s']] # Get it to be class-advector
  SB_g2 = Eggs_g2 = rep(0, stanza_data$n_g2 )
  Wmat_g2 = p$Wmat_g2
  Amat_g2 = stanza_data$stanzainfo_g2z[,'Amat']
  Wmatslope_g2 = stanza_data$stanzainfo_g2z[,'Wmatslope']
  #Wmat_g2 = exp(p$log_winf_z[1]) * stanza_data$stanzainfo_g2z[,'Wmat']
  #Y_zz = p$Y_zz
  Y_zz_g2 = p$Y_zz_g2
  #X_zz = stanza_data$X_zz
  X_zz_g2 = stanza_data$X_zz_g2
  #d_g2 = stanza_data$stanzainfo_g2z[,'d']
  d_g2 = plogis( p$logit_d_g2 )

  # Replace Wmat with Amat if available
  #which_replace = which(!is.na(Amat_g2))
  #if( length(which_replace) > 0 ){
  #  # Wmat = (1 - exp(-K * Amat)) ^ (1/(1-d))
  #  Wmat[which_replace] = (1 - exp(-Amat_g2[which_replace] * exp(p$log_K_g2[which_replace]))) ^ (1/(1-d_g2[which_replace]))
  #}

  # Increase age given plus group
  fmax = function(a,b){
    (a + b + abs(a-b) ) / 2
  }
  logspace_add <- function(logx, logy) {
    # https://github.com/kaskr/adcomp/issues/236
    # https://stackoverflow.com/questions/65233445/how-to-calculate-sums-in-log-space-without-underflow
    fmax(logx, logy) + log1p(exp(-abs(logx - logy)))
  }
  increase_vector = function(vec, plus_group="average"){
    "c" <- ADoverload("c")  # Necessary in packages
    "[<-" <- ADoverload("[<-")
    out = c( 0, vec[-length(vec)] )
    if(plus_group=="add"){
      out[length(out)] = out[length(out)] + vec[length(out)]
    }
    if(plus_group=="logspace_add"){
      out[length(out)] = logspace_add( out[length(out)], vec[length(out)] )
      #out[length(out)] = log( exp(out[length(out)]) + exp(vec[length(out)]) )
    }
    if(plus_group=="average"){
      out[length(out)] = (out[length(out)] + vec[length(out)]) / 2
    }
    return(out)
  }
  #pos = function(x){
  #  "c" <- ADoverload("c")
  #  "[<-" <- ADoverload("[<-")
  #  0.5*(abs(x)+x)
  #}
  # Loop through stanza-variables
  # Only does a single STEP of STEPS_PER_YEAR: Allows p to have updated B_t for each step
  Z_s2 = QB_s2 = rep( 0, nrow(stanza_data$stanzainfo_s2z) )
  for( s2 in seq_len(nrow(stanza_data$stanzainfo_s2z)) ){
    g2 = stanza_data$stanzainfo_s2z[s2,'g2']
    t2 = stanza_data$stanzainfo_s2z[s2,'t2']
    s = stanza_data$stanzainfo_s2z[s2,'s']
    Xg2_zz = X_zz_g2[[g2]]
    Yg2_zz = Y_zz_g2[[g2]]

    #
    stanzainfo_t2z = stanza_data$stanzainfo_s2z[which(stanza_data$stanzainfo_s2z[,'g2']==g2),,drop=FALSE]
    which_z = which( Xg2_zz[,'t2']==stanza_data$stanzainfo_s2z[s2,'t2'] )

    #
    stanzaPred = sum( exp(Yg2_zz[which_z,'log_NageS']) * Yg2_zz[which_z,'QageS'] )
    state_Biomass = sum( exp(Yg2_zz[which_z,'log_NageS']) * Yg2_zz[which_z,'WageS'] )

    # Copy ecosim.cpp#L890
    Z_s2[s2] = (LossPropToB_s[s] / state_Biomass) + F_s[s]
    #Su = exp(-Zrate / STEPS_PER_YEAR);
    log_Su = -Z_s2[s2] / STEPS_PER_YEAR
    QB_s2[s2] = FoodGain_s[s] / stanzaPred

    # Vectorized version
    #Yg2_zz[which_z,'NageS'] = Yg2_zz[which_z,'NageS'] * Su;
    Yg2_zz[which_z,'log_NageS'] = Yg2_zz[which_z,'log_NageS'] + log_Su;
    Yg2_zz[which_z,'WageS'] = vbm_g2[g2] * Yg2_zz[which_z,'WageS'] + QB_s2[s2] * Yg2_zz[which_z,'SplitAlpha']
    Y_zz_g2[[g2]] = Yg2_zz
    # W(a+1) = vbm_g2 * W(a) + (Food/Pred) * SplitAlpha(a)
    # SplitAlpha = (W(a+1) - vbm_g2 * W(a)) / (Food/Pred)
  }

  # Loop through multi-stanza groups
  # Could replace sum(exp(x)) using: https://stackoverflow.com/questions/65233445/how-to-calculate-sums-in-log-space-without-underflow
  for( g2 in seq_len(stanza_data$n_g2) ){
    # Record SpawnBiomass
    Yg2_zz = Y_zz_g2[[g2]]
    #SB_g2[g2] = sum(Yg2_zz[,'NageS'] * pos(Yg2_zz[,'WageS'] - Wmat_g2[g2]))
    SB_g2[g2] = sum( exp(Yg2_zz[,'log_NageS']) * fecundity_by_weight(Yg2_zz[,'WageS'], Wmat_g2[g2], Wmatslope_g2[g2]) )

    # Plus-group
    #Yg2_zz[,'NageS'] = increase_vector(Yg2_zz[,'NageS'], plus_group="add")
    Yg2_zz[,'log_NageS'] = increase_vector(Yg2_zz[,'log_NageS'], plus_group="logspace_add")
    Yg2_zz[,'WageS'] = increase_vector(Yg2_zz[,'WageS'], plus_group="average")

    # Eggs
    Eggs_g2[g2] = SB_g2[g2] * p$SpawnX_g2[g2] / (p$SpawnX_g2[g2] - 1.0 + (SB_g2[g2] / p$baseSB_g2[g2]))
    # Apply to first age
    #Yg2_zz[1,'NageS'] = p$baseR0_g2[g2] * EggsStanza / p$baseEggs_g2[g2] * exp(p$phi_g2[g2])
    Yg2_zz[1,'log_NageS'] = log(p$baseR0_g2[g2] * Eggs_g2[g2] / p$baseEggs_g2[g2]) + p$phi_g2[g2]
    Yg2_zz[1,'WageS'] = 0
    # Update QageS ... needed to calculate expected ration
    Yg2_zz[,'QageS'] = Yg2_zz[,'WageS'] ^ d_g2[g2]
    Y_zz_g2[[g2]] = Yg2_zz
  }


  # Update and return
  out = list(
    Y_zz_g2 = Y_zz_g2,
    SB_g2 = SB_g2,
    Z_s2 = Z_s2,
    QB_s2 = QB_s2,
    Eggs_g2 = Eggs_g2
  )
  return(out)
}

#' @export
get_stanza_total <-
function( stanza_data,
          Y_zz_g2,
          what = c("Biomass","Abundance") ){

  # Loop through stanza-variables
  what = match.arg(what)
  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  Y_s2 = rep(0, stanza_data$n_s2)
  #X_zz = stanza_data$X_zz
  for( s2 in seq_len(stanza_data$n_s2) ){
    g2 = stanza_data$stanzainfo_s2z[s2,'g2']
    t2 = stanza_data$stanzainfo_s2z[s2,'t2']
    s = stanza_data$stanzainfo_s2z[s2,'s']
    Xg2_zz = stanza_data$X_zz_g2[[g2]]
    Yg2_zz = Y_zz_g2[[g2]]

    #
    stanzainfo_t2z = stanza_data$stanzainfo_s2z[which(stanza_data$stanzainfo_s2z[,'g2']==g2),,drop=FALSE]
    which_z = which( Xg2_zz[,'t2'] == stanzainfo_t2z[t2,'t2'] )
    #if(what=="Biomass") Y_s2[s2] = sum(Yg2_zz[which_z,'NageS'] * Yg2_zz[which_z,'WageS'], na.rm=TRUE)
    if(what=="Biomass") Y_s2[s2] = sum( exp(Yg2_zz[which_z,'log_NageS']) * Yg2_zz[which_z,'WageS'], na.rm=TRUE)
    #if(what=="Abundance") Y_s2[s2] = sum(Yg2_zz[which_z,'NageS'], na.rm=TRUE)
    if(what=="Abundance") Y_s2[s2] = sum( exp(Yg2_zz[which_z,'log_NageS']), na.rm=TRUE)
  }
  return(Y_s2)
}

#' @export
project_stanzas <-
function( p,
          stanza_data,
          y,
          #xset,
          #increase_age = TRUE,
          correct_errors = FALSE,
          record_steps = FALSE,
          STEPS_PER_YEAR ){

  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  xset = seq( 1, nrow(y), length=STEPS_PER_YEAR+1)
  xset = round( rowMeans(cbind(xset[-length(xset)],xset[-1])) )
  if(record_steps) record = NULL
  #Y_zz = p$Y_zz
  Y_zz_g2 = p$Y_zz_g2
  TotalSB_g2 = TotalEggs_g2 = rep( 0, stanza_data$n_g2 )
  TotalZ_s2 = rep( 0, stanza_data$n_s2 )

  # Project
  for( STEP in seq_len(STEPS_PER_YEAR) ){
    # Load back in for update
    p$Y_zz_g2 = Y_zz_g2
    # Get food gain
    dBdt_step = dBdt( Time = 0,
              State = y[xset[STEP],],
              #State = out$B_g2
              Pars = p,
              what = "all")
    FoodGain = colSums(dBdt_step$Q_ij)
    # Update numbers
    updated_values = update_stanzas(
                  p = p,
                  stanza_data = stanza_data,
                  FoodGain_s = FoodGain,
                  LossPropToB_s = dBdt_step$M_i,
                  F_s = exp(p$logF_i),
                  STEPS_PER_YEAR = STEPS_PER_YEAR )
    Y_zz_g2 = updated_values$Y_zz_g2
    TotalEggs_g2 = TotalEggs_g2 + updated_values$Eggs_g2
    TotalSB_g2 = TotalSB_g2 + updated_values$SB_g2
    TotalZ_s2 = TotalZ_s2 + updated_values$Z_s2
    #if(record_steps){
    #  B_s2 = get_stanza_total( stanza_data = stanza_data,
    #                             Y_zz = Y_zz )
    #  record = rbind(record, B_s2)
    #}
  }

  if(correct_errors){
    # Calculate ending biomass
    B_s2 = get_stanza_total( stanza_data = stanza_data,
                               #Y_zz = Y_zz )
                               Y_zz_g2 = Y_zz_g2 )
    # Loop through multi-stanza groups
    for( g2 in seq_len(stanza_data$n_g2) ){
      stanzainfo_t2z = stanza_data$stanzainfo_s2z[which(stanza_data$stanzainfo_s2z[,'g2']==g2),,drop=FALSE]
      error_t2 = B_s2[stanzainfo_t2z[,'s2']] / y[nrow(y),stanzainfo_t2z[,'s']]
      #Y_zz_g2[[g2]][,'NageS'] = Y_zz_g2[[g2]][,'NageS'] / error_t2[stanza_data$X_zz_g2[[g2]][,'t2']]
      Y_zz_g2[[g2]][,'log_NageS'] = Y_zz_g2[[g2]][,'log_NageS'] - log(error_t2[stanza_data$X_zz_g2[[g2]][,'t2']])
    }
  }

  # Calculate ending biomass
  B_s2 = get_stanza_total( stanza_data = stanza_data,
                             #Y_zz = Y_zz )
                             Y_zz_g2 = Y_zz_g2 )

  # BUndle and return
  out = list( Y_zz_g2 = Y_zz_g2,
              B_s2 = B_s2,
              TotalEggs_g2 = TotalEggs_g2,
              TotalSB_g2 = TotalSB_g2,
              TotalZ_s2 = TotalZ_s2 )
  #if(record_steps) out$record = record
  return(out)
}

#' @title Detailed control for stanza structure
#'
#' @description Define a list of control parameters.
#'
#'
#' @export
stanza_settings <-
function( taxa,
          stanza_groups,
          K,
          d,
          Wmat,
          Amax,
          SpawnX,
          Leading,
          fit_K = c(),
          fit_d = c(),
          Amat = NULL,
          Wmatslope,
          STEPS_PER_YEAR = 1,
          comp_weight = c("multinom","dir","dirmult"),
          correct_errors = FALSE,
          min_agecomp_prob = 0 ){

  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  #
  comp_weight = match.arg(comp_weight)
  if(missing(stanza_groups)) stanza_groups = vector()
  unique_stanza_groups = setdiff(stanza_groups[taxa], NA)
  multigroup_taxa = names(stanza_groups)[which(stanza_groups %in% unique_stanza_groups)]

  # More defaults
  if(missing(SpawnX)) SpawnX = array(2, dim=length(unique_stanza_groups), dimnames=list(unique_stanza_groups))
  if(missing(d)) d = array(2/3, dim=length(unique_stanza_groups), dimnames=list(unique_stanza_groups))
  if(missing(Wmatslope)) Wmatslope = array(Inf, dim=length(unique_stanza_groups), dimnames=list(unique_stanza_groups))

  if( length(unique_stanza_groups)==0 ){
    K = d = Wmat = Amax = vector()
  }
  if(missing(Leading)){
    Leading = unlist(tapply( Amax, FUN=\(v) v==max(v), INDEX=match(stanza_groups, unique_stanza_groups) ))
  }
  names(Leading) = names(Amax)

  #
  if( is.null(Amat) ){
    Amat = rep(NA, length(unique_stanza_groups))
    names(Amat) = unique_stanza_groups
  }

  # Return
  structure( list(
    taxa = taxa,
    stanza_groups = stanza_groups,
    unique_stanza_groups = unique_stanza_groups,
    multigroup_taxa = multigroup_taxa,
    K = K,
    d = d,
    Wmat = Wmat,
    Amat = Amat,
    Wmatslope = Wmatslope,
    Amax = Amax,
    fit_K = fit_K,
    fit_d = fit_d,
    Leading = Leading,
    SpawnX = SpawnX,
    STEPS_PER_YEAR = STEPS_PER_YEAR,
    comp_weight = comp_weight,
    n_g2 = length(unique_stanza_groups),
    correct_errors = correct_errors,
    min_agecomp_prob = min_agecomp_prob
  ), class = "stanza_settings" )
}

