
#' @title
#' Compute values for multi-stanza groups
#'
#' @inheritParams ecostate
#' @inheritParams ecostate_control
#'
#' @description
#' Compute various objects used for multi-stanza groups
#'
#' @param p list of parameters
#'
#' @details
#' todo
#'
#' @export
add_stanzas <-
function( p,
          #taxa,
          #stanza_groups,
          #K,
          #d,
          #Wmat,
          #Amax,
          #SpawnX,
          #STEPS_PER_YEAR = 1,
          settings ){

  # Indexing
  s_s2 = match(names(settings$stanza_groups), settings$taxa)
  g2_s2 = match( settings$stanza_groups, settings$unique_stanza_groups )
  t2_s2 = sapply( seq_along(g2_s2), FUN=\(i)cumsum(g2_s2[i]==g2_s2)[i] )
  n_s2 = length(s_s2)
  n_g2 = length(settings$unique_stanza_groups)

  # Checks for inputs for each taxa
  K_g2 = settings$K[settings$unique_stanza_groups]
  d_g2 = settings$d[settings$unique_stanza_groups]
  Wmat_g2 = settings$Wmat[settings$unique_stanza_groups]
  SpawnX_g2 = settings$SpawnX[settings$unique_stanza_groups]
  Z_s2 = exp(p$logPB_i)[settings$multigroup_taxa]
  Amax_s2 = settings$Amax[settings$multigroup_taxa]
  Leading_s2 = unlist(tapply( Amax_s2, FUN=\(v) v==max(v), INDEX=g2_s2))

  # Variable-specific parameters
  stanzainfo_s2z = cbind( "s" = s_s2,
                              "s2" = seq_along(s_s2),
                              "g2" = g2_s2,
                              "lead" = Leading_s2,
                              "t2" = t2_s2,
                              "amax" = Amax_s2,
                              "Z" = Z_s2,
                              "Consumption" = NA,
                              "Q" = NA,
                              "leading_ratio" = NA,
                              "Gf" = NA )

  # Globals
  vbm_g2 = (1 - 3 * K_g2 / settings$STEPS_PER_YEAR)
  pos = function(x) ifelse(x>0,x,0)

  # EASIER TO LOOP THROUGH STANZAS
  Y_az_g2 = vector("list", length=n_g2)
  baseRzeroS = baseEggsStanza = baseSpawnBio = rep(NA, n_g2)
  plusage_g2 = tapply( stanzainfo_s2z[,'amax'], INDEX=stanzainfo_s2z[,'g2'], FUN=max )
  for( g2 in seq_len(n_g2) ){
    #
    AGE = (seq_len(plusage_g2[g2] * settings$STEPS_PER_YEAR) - 1) / settings$STEPS_PER_YEAR

    #
    stanzainfo_t2z = stanzainfo_s2z[which(stanzainfo_s2z[,'g2']==g2),]
    which_s2 = stanzainfo_t2z[,'s2']
    which_leading = which(stanzainfo_t2z[,'lead']==1)
    leading_s2 = stanzainfo_t2z[which_leading,'s2']

    #
    t2_a = sapply( AGE, FUN=\(a){sum(a>stanzainfo_t2z[,'amax'])} ) + 1
    t2_a = ifelse( t2_a > sum(g2_s2==g2), sum(g2_s2==g2), t2_a )

    # WageS and QageS
    k = K_g2[g2] * 3
    d = d_g2[g2]
    WageS = (1 - exp(-k * (1 - d) * (AGE))) ^ (1 / (1 - d))
    QageS = WageS ^ d

    # SurvS
    Zrate = Z_s2[which_s2[t2_a]] / settings$STEPS_PER_YEAR
    SurvS = exp(-1 * c(0, cumsum(Zrate)[-length(Zrate)] ))
    # Correct plus-group
    SurvS[length(SurvS)] = SurvS[length(SurvS)] / (1 - exp(-Zrate[length(Zrate)]))

    # Solve for R0 and NageS
    which_a = which( t2_a == stanzainfo_t2z[which(stanzainfo_t2z[,'lead']==1),'t2'] )
    BperR = sum(SurvS[which_a] * WageS[which_a], na.rm=TRUE)
    baseRzeroS[g2] = exp(p$logB_i[stanzainfo_t2z[which(stanzainfo_t2z[,'lead']==1),'s']]) / BperR
    NageS = SurvS * baseRzeroS[g2]

    #
    baseEggsStanza[g2] = sum(pos(NageS * (WageS - Wmat_g2[g2])), na.rm=TRUE)
    baseSpawnBio[g2] = sum(pos(NageS * (WageS - Wmat_g2[g2])), na.rm=TRUE)

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
    SplitAlpha = rep(NA, length=plusage_g2[g2] * settings$STEPS_PER_YEAR)
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
    Y_az_g2[[g2]] = cbind(
      AGE = AGE,
      t2 = t2_a,
      WageS = WageS,
      NageS = NageS,
      #SurvS = SurvS,
      QageS = QageS,
      SplitAlpha = SplitAlpha
    )
  }

  # Add and output
  p = c(p, list(
    n_g2 = n_g2,
    Y_az_g2 = Y_az_g2,
    baseEggs_g2 = baseEggsStanza,
    baseSB_g2 = baseSpawnBio,
    baseR0_g2 = baseRzeroS,
    SpawnX_g2 = SpawnX_g2,
    d_g2 = d_g2,
    stanzainfo_s2z = stanzainfo_s2z,
    plusage_g2 = plusage_g2,
    K_g2 = K_g2,
    Wmat_g2 = Wmat_g2
  ))
  return(p)
}

#' @title
#' Update values for multi-stanza groups
#'
#' @inheritParams ecostate
#' @inheritParams ecostate_control
#'
#' @description
#' Update various objects used for multi-stanza groups
#'
#' @param p list of parameters
#'
#' @details
#' todo
#'
#' @export
update_stanzas <-
function( p,
          FoodGain_s,
          LossPropToB_s,
          F_s,
          #increase_age = TRUE,
          STEPS_PER_YEAR = 1 ){

  # Globals
  vbm_g2 = (1 - 3 * p$K_g2 / STEPS_PER_YEAR)
  SB_g2 = rep(0, p$n_g2)
  Y_az_g2 = p$Y_az_g2
  Wmat_g2 = p$Wmat_g2

  # Increase age given plus group
  increase_vector = function(vec, plus_group="average"){
    out = c( NA, vec[-length(vec)] )
    out[length(out)] = out[length(out)] + vec[length(out)]
    if(plus_group=="average") out[length(out)] = out[length(out)] / 2
    return(out)
  }
  pos = function(x) ifelse(x>0,x,0)

  # Loop through stanza-variables
  # Only does a single STEP of STEPS_PER_YEAR: Allows p to have updated B_t for each step
  for( s2 in seq_len(nrow(p$stanzainfo_s2z)) ){
    g2 = p$stanzainfo_s2z[s2,'g2']
    t2 = p$stanzainfo_s2z[s2,'t2']
    s = p$stanzainfo_s2z[s2,'s']

    #
    #AGE = seq_len(p$plusage_g2[g2] * STEPS_PER_YEAR) - 1
    stanzainfo_t2z = p$stanzainfo_s2z[which(p$stanzainfo_s2z[,'g2']==g2),]
    #t2_a = sapply( AGE, FUN=\(a){sum(a>stanzainfo_t2z$amax)} ) + 1
    #t2_a = ifelse( t2_a > sum(p$stanzainfo_s2z$g2==g2), sum(p$stanzainfo_s2z$g2==g2), t2_a )
    age_vec = which( Y_az_g2[[g2]][,'t2'] == stanzainfo_t2z[t2,'t2'] )

    #
    stanzaPred = sum( Y_az_g2[[g2]][age_vec,'NageS'] * Y_az_g2[[g2]][age_vec,'QageS'], na.rm=TRUE )
    state_Biomass = sum( Y_az_g2[[g2]][age_vec,'NageS'] * Y_az_g2[[g2]][age_vec,'WageS'], na.rm=TRUE )

    # Copy ecosim.cpp#L890
    Zrate = (LossPropToB_s[s] / state_Biomass) + F_s[s]
    Su = exp(-Zrate / STEPS_PER_YEAR);
    Gf = FoodGain_s[s] / stanzaPred;

    # Vectorized version
    Y_az_g2[[g2]][age_vec,'NageS'] = Y_az_g2[[g2]][age_vec,'NageS'] * Su;
    Y_az_g2[[g2]][age_vec,'WageS'] = vbm_g2[g2] * Y_az_g2[[g2]][age_vec,'WageS'] + Gf * Y_az_g2[[g2]][age_vec,'SplitAlpha']
    # W(a+1) = vbm_g2 * W(a) + (Food/Pred) * SplitAlpha(a)
    # SplitAlpha = (W(a+1) - vbm_g2 * W(a)) / (Food/Pred)
  }

  # Loop through multi-stanza groups
  for( g2 in seq_len(p$n_g2) ){
    # increase_age=FALSE is useful for testing precision of approximation ... actually, it trims off rising cohorts and breaks equilibrium
    #if( isTRUE(increase_age) ){
      # Record SpawnBiomass
      SB_g2[g2] = sum(pos(Y_az_g2[[g2]][,'NageS'] * (Y_az_g2[[g2]][,'WageS'] - Wmat_g2[g2])), na.rm=TRUE)

      # Plus-group
      Y_az_g2[[g2]][,'NageS'] = increase_vector(Y_az_g2[[g2]][,'NageS'], plus_group="add")
      Y_az_g2[[g2]][,'WageS'] = increase_vector(Y_az_g2[[g2]][,'WageS'], plus_group="average")

      # Eggs
      EggsStanza = SB_g2[g2] * p$SpawnX_g2[g2] / (p$SpawnX_g2[g2] - 1.0 + (SB_g2[g2] / p$baseSB_g2[g2]))
      # Apply to first age
      Y_az_g2[[g2]][1,'NageS'] = p$baseR0_g2[g2] * EggsStanza / p$baseEggs_g2[g2]
      Y_az_g2[[g2]][1,'WageS'] = 0
    #}
    # Update QageS ... needed to calculate expected ration
    Y_az_g2[[g2]][,'QageS'] = Y_az_g2[[g2]][,'WageS'] ^ p$d_g2[g2]
  }

  # Update and return
  return(Y_az_g2)
}

get_stanza_total <-
function( stanzainfo_s2z,
          Y_az_g2,
          what = c("Biomass","Abundance") ){
  what = match.arg(what)

  # Loop through stanza-variables
  Y_s2 = rep(NA, nrow(stanzainfo_s2z))
  for( s2 in seq_len(nrow(stanzainfo_s2z)) ){
    g2 = stanzainfo_s2z[s2,'g2']
    t2 = stanzainfo_s2z[s2,'t2']
    s = stanzainfo_s2z[s2,'s']

    #
    stanzainfo_t2z = stanzainfo_s2z[which(stanzainfo_s2z[,'g2']==g2),]
    age_vec = which( Y_az_g2[[g2]][,'t2'] == stanzainfo_t2z[t2,'t2'] )
    if(what=="Biomass") Y_s2[s2] = sum(Y_az_g2[[g2]][age_vec,'NageS'] * Y_az_g2[[g2]][age_vec,'WageS'], na.rm=TRUE)
    if(what=="Abundance") Y_s2[s2] = sum(Y_az_g2[[g2]][age_vec,'NageS'], na.rm=TRUE)
  }
  return(Y_s2)
}

project_stanzas <-
function( p,
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

  # Change time-scale for F and epsilon
  # (F is already changed in update_stamzas)
  #p$logF_i = p$logF_i - log(STEPS_PER_YEAR)
  #p$epsilon_i = p$epsilon_i - log(STEPS_PER_YEAR)

  # Project
  for( STEP in seq_len(STEPS_PER_YEAR) ){
    dBdt_step = dBdt( Time = 0,
              State = y[xset[STEP],],
              #State = out$B_g2
              Pars = p,
              what = "all")
    p$Y_az_g2 = update_stanzas( p = p,
                  FoodGain_s = colSums(dBdt_step$Q_ij),
                  LossPropToB_s = dBdt_step$G_i,
                  F_s = exp(p$logF_i),
                  #increase_age = increase_age,
                  STEPS_PER_YEAR = STEPS_PER_YEAR )
    if(record_steps){
      B_s2 = get_stanza_total( stanzainfo_s2z = p$stanzainfo_s2z,
                                 Y_az_g2 = p$Y_az_g2 )
      record = rbind(record, B_s2)
    }
  }

  if(correct_errors){
    # Calculate ending biomass
    B_s2 = get_stanza_total( stanzainfo_s2z = p$stanzainfo_s2z,
                               Y_az_g2 = p$Y_az_g2 )
    # Loop through multi-stanza groups
    for( g2 in seq_len(p$n_g2) ){
      stanzainfo_t2z = p$stanzainfo_s2z[which(p$stanzainfo_s2z[,'g2']==g2),]
      error_t2 = B_s2[stanzainfo_t2z[,'s2']] / y[nrow(y),stanzainfo_t2z[,'s']]
      p$Y_az_g2[[g2]][,'NageS'] = p$Y_az_g2[[g2]][,'NageS'] / error_t2[p$Y_az_g2[[g2]][,'t2']]
    }
  }

  # Calculate ending biomass
  B_s2 = get_stanza_total( stanzainfo_s2z = p$stanzainfo_s2z,
                             Y_az_g2 = p$Y_az_g2 )

  # BUndle and return
  out = list( Y_az_g2 = p$Y_az_g2,
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
  ), class = "stanza_settings" )
}

