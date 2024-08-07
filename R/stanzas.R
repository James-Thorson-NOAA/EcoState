
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
          taxa,
          stanzas,
          K,
          Z,
          d,
          Wmat,
          Amax,
          SpawnX,
          STEPS_PER_YEAR = 1 ){

  #
  if(missing(stanzas)) stanzas = array(NA, dim=length(taxa), dimnames=list(taxa))

  #
  unique_stanzas = setdiff(unique(stanzas),NA)
  multigroup_taxa = taxa[which(stanzas %in% unique_stanzas)]
  n_g2 = length(unique_stanzas)
  s_s2 = which(stanzas %in% unique_stanzas)
  n_s2 = length(s_s2)
  g2_s2 = match( stanzas[s_s2], unique_stanzas )
  t2_s2 = sapply( seq_along(g2_s2), FUN=\(i)cumsum(g2_s2[i]==g2_s2)[i] )

  # More defaults
  if(missing(SpawnX)) SpawnX = array(2, dim=length(unique_stanzas), dimnames=list(unique_stanzas))
  if(missing(d)) d = array(2/3, dim=length(unique_stanzas), dimnames=list(unique_stanzas))

  # Checks for inputs for each taxa
  K_g2 = K[unique_stanzas]
  d_g2 = d[unique_stanzas]
  Wmat_g2 = Wmat[unique_stanzas]
  SpawnX_g2 = SpawnX[unique_stanzas]
  Z_s2 = Z[multigroup_taxa]
  Amax_s2 = Amax[multigroup_taxa]
  Leading_s2 = unlist(tapply( Amax_s2, FUN=\(v)v==max(v), INDEX=g2_s2))

  # Variable-specific parameters
  stanzainfo_s2z = data.frame( "s" = s_s2,
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
  vbm_g2 = (1 - 3 * K_g2 / STEPS_PER_YEAR)
  pos = function(x) ifelse(x>0,x,0)

  # EASIER TO LOOP THROUGH STANZAS
  Y_az_g2 = vector("list", length=n_g2)
  baseRzeroS = baseEggsStanza = baseSpawnBio = rep(NA, n_g2)
  plusage_g2 = tapply( stanzainfo_s2z[,'amax'], INDEX=stanzainfo_s2z[,'g2'], FUN=max )
  for( g2 in seq_len(n_g2) ){
    #
    AGE = seq_len(plusage_g2[g2] * STEPS_PER_YEAR) - 1

    #
    stanzainfo_t2z = stanzainfo_s2z[which(stanzainfo_s2z$g2==g2),]
    which_s2 = stanzainfo_t2z$s2
    which_leading = which(stanzainfo_t2z[,'lead'])
    leading_s2 = stanzainfo_t2z[which_leading,'s2']

    #
    t2_a = sapply( AGE, FUN=\(a){sum(a>stanzainfo_t2z$amax)} ) + 1
    t2_a = ifelse( t2_a > sum(g2_s2==g2), sum(g2_s2==g2), t2_a )

    # WageS and QageS
    k = VBGF_Ksp[g2] * 3 / STEPS_PER_YEAR
    d = vBGFd[g2]
    WageS = (1 - exp(-k * (1 - d) * (AGE))) ^ (1 / (1 - d))
    QageS = WageS ^ d

    # SurvS
    Zrate = Z[which_s2[t2_a]] / STEPS_PER_YEAR
    SurvS = exp(-1 * c(0, cumsum(Zrate)[-length(Zrate)] ))

    # Solve for R0 and NageS
    which_a = which( t2_a == stanzainfo_t2z[which(stanzainfo_t2z$lead),'t2'] )
    BperR = sum(SurvS[which_a] * WageS[which_a], na.rm=TRUE)
    baseRzeroS[g2] = exp(p$logB_i[stanzainfo_t2z[which(stanzainfo_t2z$lead),'s']]) / BperR
    NageS = SurvS * baseRzeroS[g2]

    #
    baseEggsStanza[g2] = sum(pos(NageS * (WageS - Wmat[g2])), na.rm=TRUE)
    baseSpawnBio[g2] = sum(pos(NageS * (WageS - Wmat[g2])), na.rm=TRUE)

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
    SplitAlpha = rep(NA, length=plusage_g2[g2]*STEPS_PER_YEAR)
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
    SpawnX = SpawnX,
    d_g2 = d_g2,
    stanzainfo_s2z = stanzainfo_s2z,
    plusage_g2 = plusage_g2,
    K_g2 = K_g2
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
          STEPS_PER_YEAR = 1 ){

  # Globals
  vbm_g2 = (1 - 3 * p$K_g2 / STEPS_PER_YEAR)
  SB_g2 = rep(0, p$n_g2)
  Y_az_g2 = p$Y_az_g2

  # Increase age given plus group
  increase_age = function(vec, plus_group="average"){
    out = c( NA, vec[-length(vec)] )
    out[length(out)] = out[length(out)] + vec[length(out)]
    if(plus_group=="average") out[length(out)] = out[length(out)] / 2
    return(out)
  }
  pos = function(x) ifelse(x>0,x,0)

  # Loop through stanza-variables
  for( s2 in seq_len(nrow(p$stanzainfo_s2z)) ){
    g2 = p$stanzainfo_s2z[s2,'g2']
    t2 = p$stanzainfo_s2z[s2,'t2']
    s = p$stanzainfo_s2z[s2,'s']

    #
    #AGE = seq_len(p$plusage_g2[g2] * STEPS_PER_YEAR) - 1
    stanzainfo_t2z = p$stanzainfo_s2z[which(p$stanzainfo_s2z$g2==g2),]
    #t2_a = sapply( AGE, FUN=\(a){sum(a>stanzainfo_t2z$amax)} ) + 1
    #t2_a = ifelse( t2_a > sum(p$stanzainfo_s2z$g2==g2), sum(p$stanzainfo_s2z$g2==g2), t2_a )
    age_vec = which( Y_az_g2[[g2]][,'t2'] == stanzainfo_t2z[t2,'t2'] )

    #
    stanzaPred = sum( Y_az_g2[[g2]][age_vec,'NageS'] * Y_az_g2[[g2]][age_vec,'QageS'], na.rm=TRUE )
    state_Biomass = sum( Y_az_g2[[g2]][age_vec,'NageS'] * Y_az_g2[[g2]][age_vec,'WageS'], na.rm=TRUE )

    # Copy ecosim.cpp#L890
    Z = (LossPropToB_s[s] / state_Biomass) + F_s[s]
    Su = exp(-Z / STEPS_PER_YEAR);
    Gf = FoodGain_s[s] / stanzaPred;

    # Vectorized version
    Y_az_g2[[g2]][age_vec,'NageS'] = Y_az_g2[[g2]][age_vec,'NageS'] * Su;
    Y_az_g2[[g2]][age_vec,'WageS'] = vbm_g2[g2] * Y_az_g2[[g2]][age_vec,'WageS'] + Gf * Y_az_g2[[g2]][age_vec,'SplitAlpha']
    # W(a+1) = vbm_g2 * W(a) + (Food/Pred) * SplitAlpha(a)
    # SplitAlpha = (W(a+1) - vbm_g2 * W(a)) / (Food/Pred)
  }

  # Loop through multi-stanza groups
  for( g2 in seq_len(p$n_g2) ){
    # Record SpawnBiomass
    SB_g2[g2] = sum(pos(Y_az_g2[[g2]][,'NageS'] * (Y_az_g2[[g2]][,'WageS'] - Wmat[g2])), na.rm=TRUE)

    # Plus-group
    Y_az_g2[[g2]][,'NageS'] = increase_age(Y_az_g2[[g2]][,'NageS'], plus_group="add")
    Y_az_g2[[g2]][,'WageS'] = increase_age(Y_az_g2[[g2]][,'WageS'], plus_group="average")

    # Eggs
    EggsStanza = SB_g2[g2] * p$SpawnX[g2] / (p$SpawnX[g2] - 1.0 + (SB_g2[g2] / p$baseSB_g2[g2]))
    # Apply to first age
    #NageS_new[1,g2] = RscaleSplit[g2] * baseRzeroS[g2] * pow(double(EggsStanza[g2] / baseEggsStanza[g2]), double(RecPower[g2]))
    Y_az_g2[[g2]][1,'NageS'] = p$baseR0_g2[g2] * EggsStanza / p$baseEggs_g2[g2]
    Y_az_g2[[g2]][1,'WageS'] = 0
    # Update QageS ... needed to calculate expected ration
    Y_az_g2[[g2]][,'QageS'] = Y_az_g2[[g2]][,'WageS'] ^ p$d_g2[g2]
  }

  # Loop through stanza-variables
  B_s2 = rep(NA, nrow(p$stanzainfo_s2z))
  for( s2 in seq_len(nrow(p$stanzainfo_s2z)) ){
    g2 = p$stanzainfo_s2z[s2,'g2']
    t2 = p$stanzainfo_s2z[s2,'t2']
    s = p$stanzainfo_s2z[s2,'s']

    #
    #AGE = seq_len(p$plusage_g2[g2] * STEPS_PER_YEAR) - 1
    stanzainfo_t2z = p$stanzainfo_s2z[which(p$stanzainfo_s2z$g2==g2),]
    #t2_a = sapply( AGE, FUN=\(a){sum(a>stanzainfo_t2z$amax)} ) + 1
    #t2_a = ifelse( t2_a > sum(p$stanzainfo_s2z$g2==g2), sum(p$stanzainfo_s2z$g2==g2), t2_a )
    age_vec = which( Y_az_g2[[g2]][,'t2'] == stanzainfo_t2z[t2,'t2'] )

    B_s2[s2] = sum(Y_az_g2[[g2]][age_vec,'NageS'] * Y_az_g2[[g2]][age_vec,'WageS'], na.rm=TRUE)
  }

  # Update and return
  out = list(
    Y_az_g2 = Y_az_g2,
    B_s2 = B_s2
  )
  return(out)
}

