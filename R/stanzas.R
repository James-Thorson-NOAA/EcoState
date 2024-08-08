
#' @export
make_stanza_data <-
function( settings ){

  #
  #unique_stanza_groups = setdiff(settings$stanza_groups[settings$taxa], NA)
  #multigroup_taxa = names(settings$stanza_groups)[which(settings$stanza_groups %in% unique_stanza_groups)]

  # Indexing
  s_s2 = match(names(settings$stanza_groups), taxa)
  g2_s2 = match( settings$stanza_groups, settings$unique_stanza_groups )
  t2_s2 = sapply( seq_along(g2_s2), FUN=\(i)cumsum(g2_s2[i]==g2_s2)[i] )
  n_s2 = length(s_s2)
  n_g2 = length(settings$unique_stanza_groups)

  # Checks for inputs for each taxa
  K_g2 = settings$K[settings$unique_stanza_groups]
  d_g2 = settings$d[settings$unique_stanza_groups]
  Wmat_g2 = settings$Wmat[settings$unique_stanza_groups]
  SpawnX_g2 = settings$SpawnX[settings$unique_stanza_groups]
  #Z_s2 = exp(p$logPB_i)[multigroup_taxa]
  Amax_s2 = settings$Amax[settings$multigroup_taxa]
  Leading_s2 = unlist(tapply( Amax_s2, FUN=\(v) v==max(v), INDEX=g2_s2))
  plusage_g2 = tapply( Amax_s2, INDEX=g2_s2, FUN=max )

  # Variable-specific parameters
  stanzainfo_s2z = cbind( "s" = s_s2,
                          "s2" = seq_along(s_s2),
                          "g2" = g2_s2,
                          "lead" = Leading_s2,
                          "t2" = t2_s2,
                          "amax" = Amax_s2,
                          #"Z" = exp(p$logPB_i)[multigroup_taxa],
                          "Consumption" = NA,
                          "Q" = NA,
                          "leading_ratio" = NA,
                          "Gf" = NA )

  stanzainfo_g2z = cbind( "Wmat" = Wmat_g2,
                          "SpawnX" = SpawnX_g2,
                          "plusage" = plusage_g2,
                          "K" = K_g2,
                          "d" = d_g2,
                          "lead_s" = stanzainfo_s2z[which(stanzainfo_s2z[,'lead']==1),'s'] )

  #n_s2 = settings$n_s2
  #n_g2 = settings$n_g2
  #K_g2 = settings$K_g2
  #d_g2 = settings$d_g2
  #Wmat_g2 = settings$Wmat_g2
  #SpawnX_g2 = settings$SpawnX_g2
  #plusage_g2 = settings$plusage_g2
  #stanzainfo_s2z = settings$stanzainfo_s2z

  # Globals
  vbm_g2 = (1 - 3 * K_g2 / settings$STEPS_PER_YEAR)

  # EASIER TO LOOP THROUGH STANZAS
  X_zz = NULL
  for( g2 in seq_len(n_g2) ){
    # Make fractional-age starting at age=0
    AGE = (seq_len(plusage_g2[g2] * settings$STEPS_PER_YEAR) - 1) / settings$STEPS_PER_YEAR
    stanzainfo_t2z = stanzainfo_s2z[which(stanzainfo_s2z[,'g2']==g2),]

    #
    t2_a = sapply( AGE, FUN=\(a){sum(a>stanzainfo_t2z[,'amax'])} ) + 1
    t2_a = ifelse( t2_a > sum(stanzainfo_s2z[,'g2']==g2), sum(stanzainfo_s2z[,'g2']==g2), t2_a )

    # Stack
    Xg2_zz = cbind(
      g2 = g2,
      AGE = AGE,
      t2 = t2_a,
      is_lead = Leading_s2[t2_a]
    )
    X_zz = rbind( X_zz, Xg2_zz )   # deparse.level=0 avoids RTMB error
  }

  # Add and output
  stanza_data = list(
    n_s2 = n_s2,
    n_g2 = n_g2,
    #SpawnX_g2 = SpawnX_g2,
    #d_g2 = d_g2,
    stanzainfo_s2z = stanzainfo_s2z,
    stanzainfo_g2z = stanzainfo_g2z,
    #plusage_g2 = plusage_g2,
    #K_g2 = K_g2,
    X_zz = X_zz
    #Wmat_g2 = Wmat_g2
  )
  return(stanza_data)
}

#' @export
add_stanza_params <-
function( p,
          stanza_data,
          settings ){


  #s_s2 = settings$s_s2
  #g2_s2 = settings$g2_s2
  #t2_s2 = settings$t2_s2
  n_s2 = stanza_data$n_s2
  n_g2 = stanza_data$n_g2
  K_g2 = stanza_data$stanzainfo_g2z[,'K']
  d_g2 = stanza_data$stanzainfo_g2z[,'d']
  Wmat_g2 = stanza_data$stanzainfo_g2z[,'Wmat']
  SpawnX_g2 = stanza_data$stanzainfo_g2z[,'SpawnX']
  #Z_s2 = settings$Z_s2
  #Amax_s2 = settings$Amax_s2
  #Leading_s2 = settings$Leading_s2
  #Z_s2 = exp(p$logPB_i)[multigroup_taxa]
  plusage_g2 = stanza_data$stanzainfo_g2z[,'plusage']
  stanzainfo_s2z = stanza_data$stanzainfo_s2z
  stanzainfo_g2z = stanza_data$stanzainfo_g2z
  Z_s2 = exp(p$logPB_i)[stanzainfo_s2z[,'s']]
  X_zz = stanza_data$X_zz

  # Globals
  vbm_g2 = (1 - 3 * K_g2 / settings$STEPS_PER_YEAR)
  #pos = function(x) ifelse(x>0,x,0)
  pos = function(x) (x + sqrt(x^2))/2

  # EASIER TO LOOP THROUGH STANZAS
  #Y_zz = matrix(nrow=0, ncol=3)
  baseRzeroS = baseEggsStanza = baseSpawnBio = rep(0, n_g2)
  #plusage_g2 = tapply( stanzainfo_s2z[,'amax'], INDEX=stanzainfo_s2z[,'g2'], FUN=max )
  for( g2 in seq_len(n_g2) ){
    #
    Xg2_zz = X_zz[which(X_zz[,'g2']==g2),]

    # Make fractional-age starting at age=0
    AGE = Xg2_zz[,'AGE']

    #
    stanzainfo_t2z = stanzainfo_s2z[which(stanzainfo_s2z[,'g2']==g2),]
    which_s2 = stanzainfo_t2z[,'s2']
    which_leading = which(stanzainfo_t2z[,'lead']==1)
    leading_s2 = stanzainfo_t2z[which_leading,'s2']

    #
    t2_a = Xg2_zz[,'t2']

    # WageS and QageS
    k = K_g2[g2] * 3
    d = d_g2[g2]
    WageS = (1 - exp(-k * (1 - d) * (AGE))) ^ (1 / (1 - d))
    QageS = WageS ^ d

    # SurvS
    Zrate = Z_s2[which_s2[t2_a]] / settings$STEPS_PER_YEAR
    SurvS = exp(-1 * c(0, cumsum(Zrate)[-length(Zrate)] ))

    ################## EXPERIMENT WITH RTMB
    # Correct plus-group
    SurvS[length(SurvS)] = SurvS[length(SurvS)] / (1 - exp(-Zrate[length(Zrate)]))
    ################## EXPERIMENT WITH RTMB

    # Solve for R0 and NageS
    #which_a = which( Xg2_zz[,'is_lead']==1 )
    BperR = sum(SurvS * WageS * Xg2_zz[,'is_lead'])
    #which_s = stanzainfo_t2z[which(stanzainfo_t2z[,'lead']==1),'s']
    ################## EXPERIMENT WITH RTMB
    baseRzeroS[g2] = exp(p$logB_i)[stanzainfo_g2z[g2,'lead_s']] / BperR
    #baseRzeroS[g2] = 4000
    ################## EXPERIMENT WITH RTMB
    NageS = SurvS * baseRzeroS[g2]

    #
    baseEggsStanza[g2] = sum(pos(NageS * (WageS - Wmat_g2[g2])))
    baseSpawnBio[g2] = sum(pos(NageS * (WageS - Wmat_g2[g2])))

    #
    s = stanzainfo_t2z[which_leading,'s']
    stanzainfo_s2z[leading_s2,'Consumption'] = exp(p$logQB_i[s]) * exp(p$logB_i[s])

    # Solve for expected consumption
    for( t2 in seq_len(nrow(stanzainfo_t2z)) ){
      s = stanzainfo_t2z[t2,'s']
      which_a = which( t2_a == stanzainfo_t2z[t2,'t2'] )
      #if( is.na(p$logB_i[s]) ){
        p$logB_i[s] = log(sum(NageS[which_a] * WageS[which_a], na.rm=TRUE))
      #}
      stanzainfo_s2z[stanzainfo_t2z[t2,'s2'],'Q'] = sum(NageS[which_a] * QageS[which_a], na.rm=TRUE)
    }

    # Loop again for Cons and QB
    for( t2 in seq_len(nrow(stanzainfo_t2z)) ){
      s = stanzainfo_t2z[t2,'s']
      s2 = stanzainfo_t2z[t2,'s2']
      stanzainfo_s2z[s2,'leading_ratio'] = stanzainfo_s2z[leading_s2,'Consumption'] / stanzainfo_s2z[leading_s2,'Q']
      stanzainfo_s2z[s2,'Consumption'] = stanzainfo_s2z[s2,'Q'] * stanzainfo_s2z[s2,'leading_ratio']
      #if( is.na(p$logQB_i[s]) ){
        p$logQB_i[s] = log(stanzainfo_s2z[s2,'Consumption']) - p$logB_i[s]
      #}
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
      SplitAlpha[which_a] = (WageS[which_a+1] - vbm_g2[g2] * WageS[which_a]) / stanzainfo_s2z[s2,'leading_ratio']
      SplitAlpha[length(SplitAlpha)] = SplitAlpha[length(SplitAlpha)-1]
    }

    # Stack
    ################## EXPERIMENT WITH RTMB
    Y_zz = matrix( nrow=length(WageS), ncol=4 )
    colnames(Y_zz) = c("WageS", "NageS", "QageS", "SplitAlpha")
    Y_zz[,'WageS'] = WageS
    Y_zz[,'NageS'] = NageS
    Y_zz[,'QageS'] = QageS
    Y_zz[,'SplitAlpha'] = SplitAlpha
    #Y_zz = rbind( Y_zz, Yg2_zz )  # deparse.level=0 avoids RTMB error
    ################## EXPERIMENT WITH RTMB
  }

  # Add and output
  p = c(p, list(
    Y_zz = Y_zz,
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

  # Globals
  vbm_g2 = (1 - 3 * stanza_data$stanzainfo_g2z[,'K'] / STEPS_PER_YEAR)
  SB_g2 = rep(0, stanza_data$n_g2)
  Wmat_g2 = stanza_data$stanzainfo_g2z[,'Wmat']
  SpawnX_g2 = stanza_data$stanzainfo_g2z[,'SpawnX']
  Y_zz = p$Y_zz
  X_zz = stanza_data$X_zz
  d_g2 = stanza_data$stanzainfo_g2z[,'d']

  # Increase age given plus group
  increase_vector = function(vec, plus_group="average"){
    out = c( 0, vec[-length(vec)] )
    out[length(out)] = out[length(out)] + vec[length(out)]
    if(plus_group=="average") out[length(out)] = out[length(out)] / 2
    return(out)
  }
  #pos = function(x) ifelse(x>0,x,0)
  pos = function(x) (x + sqrt(x^2))/2

  # Loop through stanza-variables
  # Only does a single STEP of STEPS_PER_YEAR: Allows p to have updated B_t for each step
  for( s2 in seq_len(nrow(stanza_data$stanzainfo_s2z)) ){
    g2 = stanza_data$stanzainfo_s2z[s2,'g2']
    t2 = stanza_data$stanzainfo_s2z[s2,'t2']
    s = stanza_data$stanzainfo_s2z[s2,'s']

    #
    ################## EXPERIMENT WITH RTMB
    #stanzainfo_t2z = stanza_data$stanzainfo_s2z[which(stanza_data$stanzainfo_s2z[,'g2']==g2),]
    stanzainfo_t2z = stanza_data$stanzainfo_s2z
    ################## EXPERIMENT WITH RTMB
    which_z = which( (X_zz[,'g2']==stanza_data$stanzainfo_s2z[s2,'g2']) & (X_zz[,'t2']==stanza_data$stanzainfo_s2z[s2,'t2']) )
    #which_z = which( (X_zz[,'g2']==g2) & (X_zz[z_g2,'t2']==stanzainfo_t2z[t2,'t2']) )
    #age_vec = which( X_zz[z_g2,'t2'] == stanzainfo_t2z[t2,'t2'] )
    #which_z = z_g2[age_vec]

    #
    stanzaPred = sum( Y_zz[which_z,'NageS'] * Y_zz[which_z,'QageS'] )
    state_Biomass = sum( Y_zz[which_z,'NageS'] * Y_zz[which_z,'WageS'] )

    # Copy ecosim.cpp#L890
    Zrate = (LossPropToB_s[s] / state_Biomass) + F_s[s]
    Su = exp(-Zrate / STEPS_PER_YEAR);
    Gf = FoodGain_s[s] / stanzaPred;

    # Vectorized version
    Y_zz[which_z,'NageS'] = Y_zz[which_z,'NageS'] * Su;
    Y_zz[which_z,'WageS'] = vbm_g2[g2] * Y_zz[which_z,'WageS'] + Gf * Y_zz[which_z,'SplitAlpha']
    # W(a+1) = vbm_g2 * W(a) + (Food/Pred) * SplitAlpha(a)
    # SplitAlpha = (W(a+1) - vbm_g2 * W(a)) / (Food/Pred)
  }

  # Loop through multi-stanza groups
  for( g2 in seq_len(stanza_data$n_g2) ){
    # increase_age=FALSE is useful for testing precision of approximation ... actually, it trims off rising cohorts and breaks equilibrium
    #if( isTRUE(increase_age) ){
      # Record SpawnBiomass
      which_z = which(X_zz[,'g2']==g2)
      SB_g2[g2] = sum(pos(Y_zz[which_z,'NageS'] * pos(Y_zz[which_z,'WageS'] - Wmat_g2[g2])))
      #if(isFALSE(inherits(SB_g2[g2],"advector"))) stop("check SB_g2")

      # Plus-group
      Y_zz[which_z,'NageS'] = increase_vector(Y_zz[which_z,'NageS'], plus_group="add")
      Y_zz[which_z,'WageS'] = increase_vector(Y_zz[which_z,'WageS'], plus_group="average")

      # Eggs
      EggsStanza = SB_g2[g2] * SpawnX_g2[g2] / (SpawnX_g2[g2] - 1.0 + (SB_g2[g2] / p$baseSB_g2[g2]))
      # Apply to first age
      Y_zz[which_z[1],'NageS'] = p$baseR0_g2[g2] * EggsStanza / p$baseEggs_g2[g2]
      Y_zz[which_z[1],'WageS'] = 0
    #}
    # Update QageS ... needed to calculate expected ration
    Y_zz[which_z,'QageS'] = Y_zz[which_z,'WageS'] ^ d_g2[g2]
  }

  # Update and return
  return(Y_zz)
}

get_stanza_total <-
function( stanza_data,
          Y_zz,
          what = c("Biomass","Abundance") ){

  # Loop through stanza-variables
  what = match.arg(what)
  Y_s2 = rep(0, nrow(stanza_data$stanzainfo_s2z))
  X_zz = stanza_data$X_zz
  for( s2 in seq_len(nrow(stanza_data$stanzainfo_s2z)) ){
    g2 = stanza_data$stanzainfo_s2z[s2,'g2']
    t2 = stanza_data$stanzainfo_s2z[s2,'t2']
    s = stanza_data$stanzainfo_s2z[s2,'s']

    #
    stanzainfo_t2z = stanza_data$stanzainfo_s2z[which(stanza_data$stanzainfo_s2z[,'g2']==g2),]
    age_vec = which( X_zz[which(X_zz[,'g2']==g2),'t2'] == stanzainfo_t2z[t2,'t2'] )
    which_z = which(X_zz[,'g2']==g2)[age_vec]
    if(what=="Biomass") Y_s2[s2] = sum(Y_zz[which_z,'NageS'] * Y_zz[which_z,'WageS'], na.rm=TRUE)
    if(what=="Abundance") Y_s2[s2] = sum(Y_zz[which_z,'NageS'], na.rm=TRUE)
  }
  return(Y_s2)
}

project_stanzas <-
function( p,
          stanza_data,
          y,
          #xset,
          #increase_age = TRUE,
          correct_errors = FALSE,
          record_steps = FALSE,
          STEPS_PER_YEAR ){

  #if( missing(xset) ){
    xset = seq( 1, nrow(y), length=STEPS_PER_YEAR+1)
    xset = round( rowMeans(cbind(xset[-length(xset)],xset[-1])) )
  #}else{
  #  STEPS_PER_YEAR = length(xset)
  #}
  if(record_steps) record = NULL
  Y_zz = p$Y_zz

  # Change time-scale for F and epsilon
  # (F is already changed in update_stamzas)
  #p$logF_i = p$logF_i - log(STEPS_PER_YEAR)
  #p$epsilon_i = p$epsilon_i - log(STEPS_PER_YEAR)

  #
  #data2 = local({
  #                type_i = type_i
  #                n_species = n_species
  #                F_type = control$F_type
  #                environment()
  #})
  #environment(dBdt) <- data2

  # Project
  for( STEP in seq_len(STEPS_PER_YEAR) ){
    dBdt_step = dBdt( Time = 0,
              State = y[xset[STEP],],
              #State = out$B_g2
              Pars = p,
              what = "all")
    Y_zz = update_stanzas( p = p,
                  stanza_data = stanza_data,
                  FoodGain_s = colSums(dBdt_step$Q_ij),
                  LossPropToB_s = dBdt_step$G_i,
                  F_s = exp(p$logF_i),
                  #increase_age = increase_age,
                  STEPS_PER_YEAR = STEPS_PER_YEAR )
    if(record_steps){
      B_s2 = get_stanza_total( stanza_data = stanza_data,
                                 Y_zz = Y_zz )
      record = rbind(record, B_s2)
    }
  }

  if(correct_errors){
    # Calculate ending biomass
    B_s2 = get_stanza_total( stanza_data = stanza_data,
                               Y_zz = Y_zz )
    # Loop through multi-stanza groups
    for( g2 in seq_len(stanza_data$n_g2) ){
      which_z = which(stanza_data$X_zz[,'g2']==g2)
      stanzainfo_t2z = stanza_data$stanzainfo_s2z[which(stanza_data$stanzainfo_s2z[,'g2']==g2),]
      error_t2 = B_s2[stanzainfo_t2z[,'s2']] / y[nrow(y),stanzainfo_t2z[,'s']]
      Y_zz[which_z,'NageS'] = Y_zz[which_z,'NageS'] / error_t2[stanza_data$X_zz[which_z,'t2']]
    }
  }

  # Calculate ending biomass
  B_s2 = get_stanza_total( stanza_data = stanza_data,
                             Y_zz = Y_zz )

  # BUndle and return
  out = list( Y_zz = Y_zz,
              B_s2 = B_s2 )
  if(record_steps) out$record = record
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
          STEPS_PER_YEAR = 1 ){

  #
  if(missing(stanza_groups)) stanza_groups = array(NA, dim=length(taxa), dimnames=list(taxa))
  unique_stanza_groups = setdiff(stanza_groups[taxa], NA)
  multigroup_taxa = names(stanza_groups)[which(stanza_groups %in% unique_stanza_groups)]

  # More defaults
  if(missing(SpawnX)) SpawnX = array(2, dim=length(unique_stanza_groups), dimnames=list(unique_stanza_groups))
  if(missing(d)) d = array(2/3, dim=length(unique_stanza_groups), dimnames=list(unique_stanza_groups))

  # Return
  structure( list(
    taxa = taxa,
    stanza_groups = stanza_groups,
    unique_stanza_groups = unique_stanza_groups,
    multigroup_taxa = multigroup_taxa,
    K = K,
    d = d,
    Wmat = Wmat,
    Amax = Amax,
    SpawnX = SpawnX,
    STEPS_PER_YEAR = STEPS_PER_YEAR
    #s_s2 = s_s2,
    #g2_s2 = g2_s2,
    #t2_s2 = t2_s2,
    #n_s2 = n_s2,
    #n_g2 = n_g2,
    #K_g2 = K_g2,
    #d_g2 = d_g2,
    #Wmat_g2 = Wmat_g2,
    #SpawnX_g2 = SpawnX_g2,
    #Z_s2 = Z_s2,
    #Amax_s2 = Amax_s2,
    #Leading_s2 = Leading_s2,
    #plusage_g2 = plusage_g2,
    #stanzainfo_s2z = stanzainfo_s2z
  ), class = "stanza_settings" )
}

