#' Fit an Extended Children-of-Twins model with indirect assortative mating
#' Written by Hans Fredrik Sunde (hfsu@fhi.no)
#'
#' Fits an Extended Children-of-Twins (COTS) structural equation model in OpenMx,
#' including spouses and two offspring per nuclear family. The model supports
#' monozygotic (MZ) and dizygotic (DZ) twin families, with optional inclusion of
#' full-sibling (FS), half-sibling (HS), and unknown-zygosity (UZ) groups.
#'
#' The model decomposes variance in the parental (focal) phenotype into additive
#' genetic (A), sibling-shared environment (C), twin-shared environment (T; optional),
#' and unique environment (E). Intergenerational resemblance is modeled via
#' phenotypic transmission (p), genetic transmission/confounding (a1p), and passive
#' environmental transmission/confounding (c1p). Offspring-specific variance is
#' decomposed into A (a2), C (c2), and E (e2), with an optional residual shared-
#' environmental correlation between cousins (rc).
#'
#' Assortative mating can be modeled as direct assortment on the focal phenotype
#' or as indirect assortment via a sorting factor. In the indirect assortment case,
#' the paths from A/C/T/E to the sorting factor can be freely estimated (subject to
#' identification constraints) and the model can optionally adjust for correlated
#' mate preferences (Z parameters), which otherwise can inflate co-in-law
#' correlations.
#'
#' Both continuous and dichotomous (ordinal) outcomes are supported via
#' mxExpectationNormal, using thresholds when \code{parent_ordinal} and/or
#' \code{child_ordinal} are \code{TRUE}. Model fitting can be ML or WLS.
#'
#' @param dta Data.frame containing family-level data. Must include a column named
#'   \code{fam.type} indicating group membership (strings containing "MZ", "DZ", and
#'   optionally "FS", "HS", "UZ") and the observed variables specified by
#'   \code{selVars}.
#' @param selVars Character vector of length 8 specifying observed variable names
#'   in the following order: parental pair, spouses, and offspring:
#'   \code{c(FAM1_Sib, FAM2_Sib, FAM1_Spouse, FAM2_Spouse, FAM1_CH1, FAM1_CH2, FAM2_CH1, FAM2_CH2)}.
#' @param name Character scalar, name of the top-level OpenMx model.
#' @param ... All arguments after ... must be named.
#'
#' @param a1.free,c1.free,t1.free Logical; free/fix parental A/C/T paths for the focal phenotype.
#' @param a1.val,c1.val,t1.val,e1.val Numeric starting values for parental A/C/T/E paths.
#' @param p.free,a1p.free,c1p.free Logical; free/fix intergenerational paths:
#'   phenotypic transmission (\code{p}), genetic transmission/confounding (\code{a1p}),
#'   and passive environmental transmission/confounding (\code{c1p}).
#' @param a2.free,c2.free Logical; free/fix offspring-specific A/C paths.
#' @param c2.cousin Logical; estimate residual cousin-shared environment correlation (\code{rc}).
#'
#' @param mu.free Logical; free/fix assortative mating strength parameter (\code{mu_param}).
#' @param mu.val Numeric starting value for assortative mating strength.
#' @param U.f Numeric in \{0,1\}; 0 indicates first generation of assortment (disequilibrium),
#'   1 indicates intergenerational equilibrium.
#' @param indirect.assortment Logical; if \code{FALSE}, sorting factor is equated to the focal
#'   phenotype; if \code{TRUE}, includes a separate sorting factor with its own ACTE paths.
#' @param VarS_constraint Logical; constrain variance of the sorting factor to equal the focal
#'   phenotype variance (\code{var1}).
#' @param a1s.free,c1s.free,t1s.free,e1s.free Logical; free/fix ACTE paths to the sorting factor
#'   (used when \code{indirect.assortment = TRUE}; otherwise constrained to equal corresponding
#'   focal-phenotype paths).
#' @param a1s.val,c1s.val,t1s.val,e1s.val Optional numeric starting values for sorting-factor paths.
#'   If \code{NULL}, defaults to corresponding focal-phenotype values/labels as implemented in the function.
#'
#' @param cor.preferences Logical; include sibling-correlated mate preferences (affecting co-in-law covariance).
#' @param preferences.mz.dz Logical; allow preference-correlation to differ by relatedness (MZ/DZ/UZ/HS).
#'
#' @param parental_rGE Logical or numeric; controls the gene-environment correlation parameter \code{w}.
#'   If \code{FALSE}, fixes \code{w=0}. If \code{TRUE}, estimates \code{w} with a constraint linking it to
#'   the implied offspring rGE. If numeric, fixes \code{w} to that value (must be between -1 and 1).
#' @param rGE_lbound Numeric lower bound for \code{w}.
#' @param rGE_twostep Logical; if \code{TRUE} and \code{parental_rGE=TRUE}, runs a two-step procedure:
#'   first with \code{w=0}, then re-runs estimating \code{w} under the rGE constraint.
#'
#' @param combine.FS_DZ Logical; combine full-sibling and DZ groups into one group (disables T components for parents).
#' @param use.fullsibs Logical; include FS group if present.
#' @param use.halfsibs Logical; include HS group if present.
#' @param halfsib.cor Numeric; C-correlation for half siblings (used when \code{use.halfsibs=TRUE}).
#' @param use.unknownzyg Logical; include UZ group if present.
#' @param MZprop Numeric; assumed proportion of MZ pairs in the UZ group (used to compute genetic correlation).
#'
#' @param parent_ordinal,child_ordinal Logical; treat parental and/or offspring phenotypes as ordinal.
#' @param factor.levels Vector of factor levels used for ordinal variables (default \code{c(0,1)}).
#'
#' @param fitFunction Character; "ML" or "WLS".
#' @param compute.intervals Logical; compute confidence intervals for selected derived parameters.
#' @param TryHard Logical; use mxTryHard / mxTryHardOrdinal.
#' @param RunModel Logical; if \code{FALSE}, returns the unfit OpenMx model object.
#'
#' @return If \code{RunModel=TRUE}, returns an OpenMx fit object (\code{MxModel}) from \code{mxRun}.
#'   Otherwise returns the constructed (unfit) \code{MxModel}. The model includes derived parameter
#'   algebras for variance decompositions, implied correlations, and parent-offspring decomposition.
#'
#' @details
#' Expected data structure is one row per extended family with columns given by \code{selVars} and
#' a grouping variable \code{fam.type}. The function performs basic checks for required groups and
#' will error if requested groups are absent.
#'
#' @section Dependencies:
#' Requires the \pkg{OpenMx} package.
#'
#' @seealso \code{\link[OpenMx]{mxModel}}, \code{\link[OpenMx]{mxRun}}, \code{\link[OpenMx]{mxTryHard}}
#'
#' @export


iam.cots.fun <- function(
    dta,
    selVars = c('FAM1_Sib','FAM2_Sib','FAM1_Spouse','FAM2_Spouse', 'FAM1_CH1', 'FAM1_CH2','FAM2_CH1', 'FAM2_CH2'),
    name = "iAM_COTS",
    ...,

    # Parental generation (Focal Phenotype)
    a1.free = TRUE,  # Additive genetic effects
    c1.free = TRUE, # Sibling-shared environment
    t1.free = FALSE,  # Twin-shared environment
    a1.val = .2,
    c1.val = .2,
    t1.val = 0,
    e1.val = .8,

    # Intergenerational paths
    p.free   = TRUE,  # If FALSE, phenoypic transmission constrained to 0
    a1p.free = TRUE,  # If FALSE, genetic transmission constrained to 0
    c1p.free = TRUE, # If FALSE, passive environmental transmission constrained to 0

    # Offspring generation
    a2.free = TRUE,
    c2.free = TRUE,
    c2.cousin = FALSE, # Estimate correlation between cousins' offspring-specific shared environmental factors?

    # Assortative mating options
    mu.free = TRUE,
    mu.val = .3,

    U.f = 1, # Intergenerational equilibrium (set to 0 if first generation of assortment, 1 if equilibirum)
    indirect.assortment = FALSE, # If FALSE, model will equate sorting factor and focal phenotype

    VarS_constraint = FALSE, # Will fix variance of the sorting factor to equal the focal phenotype. The model tends to do better if you instead constrain one of the paths (e.g., e1s)
    a1s.free = TRUE,
    c1s.free = TRUE,
    t1s.free = FALSE,
    e1s.free = FALSE,

    a1s.val = NULL, # If NULL, will use value of a1. If NULL & a1s.free = FALSE, then will use label of a1.
    c1s.val = NULL, # If NULL, will use value of c1. If NULL & c1s.free = FALSE, then will use label of c1.
    t1s.val = NULL, # If NULL, will use value of t1. If NULL & t1s.free = FALSE, then will use label of t1.
    e1s.val = NULL, # If NULL, will use value of e1. If NULL & e1s.free = FALSE, then will use label of e1.

    # Adjust for correlated mate-preferences? (can inflate co-in-law correlation)
    cor.preferences = F, # Freely estimated sibling-correlated mate preferences (inflates co-in-law correlations, can bias results unless accounted for)
    preferences.mz.dz = T, # Allow preference-correlation to differ across relatedness


    # Gene-environment correlation in parent-generation
    parental_rGE = FALSE, # If TRUE, constrains parental rGE to equal offspring rGE
    # Can also be specified manually (e.g., parental_rGE = 0.3)
    rGE_lbound=-1,
    rGE_twostep = TRUE, #If TRUE, runs model with w = 0, then re-runs to estimate w while constraining it to be equal the childhood generation


    # Family type options
    combine.FS_DZ = FALSE, # Combine FS and DZ group
    use.fullsibs = TRUE, # Use full siblings (FS) group

    use.halfsibs = FALSE, # Use half siblings (HS) group
    halfsib.cor = 0, # C-Correlation for half-sibs

    use.unknownzyg=FALSE, # Use twins of unknown zygosity (UZ)
    MZprop = 0.5, # Proportion MZs in UZ group (used to calculate genotypic correlation)

    # Dichotomous data
    parent_ordinal = FALSE, child_ordinal = FALSE, factor.levels = c(0,1),

    # Model options
    fitFunction = "ML",
    compute.intervals = FALSE,
    TryHard = FALSE, RunModel = TRUE
) {

  timestamp()
  requireNamespace("OpenMx", quietly = TRUE)


  famData <- as.data.frame(dta)


  selVars_pa   <- selVars[1:4] # Names of parental variables
  selVars_ch   <- selVars[5:8] # Names of offspring variables
  ntv          <- length(selVars)      # number of total variables



  # Converts to factor variables
  if (parent_ordinal) {famData[,selVars_pa] <- mxFactor(famData[,selVars_pa], levels = factor.levels); e1.free = FALSE} else {e1.free = TRUE}
  if (child_ordinal) {famData[,selVars_ch] <- mxFactor(famData[,selVars_ch], levels = factor.levels); e2.free = FALSE} else {e2.free = TRUE}




  # Data checks
  if (!any(grepl("MZ", famData$fam.type)))   {stop("No data for 'MZ' group")}
  if (use.halfsibs & !any(grepl("HS", famData$fam.type)))   {stop("No data for 'HS' group")}
  if (use.unknownzyg & !any(grepl("UZ", famData$fam.type))) {stop("No data for 'UZ' group")}

  if (combine.FS_DZ) {
    if (!any(grepl("FS", famData$fam.type)) & !any(grepl("DZ", famData$fam.type)))   {stop("No data for 'FS' and 'DZ' group")}
    t1.free <- t1s.free <- use.fullsibs <- FALSE
  } else {
    if (use.fullsibs & !any(grepl("FS", famData$fam.type)))   {stop("No data for 'FS' group")}
    if (!any(grepl("DZ", famData$fam.type)))   {stop("No data for 'DZ' group")}
  }

  if (use.halfsibs) {
    message("WARNING: the use of half siblings are not fully developed yet, especially parts concerning indirect assortative mating, assumptions about disequilibrium, and gene-environment correlations.")
    if (parental_rGE) {message("WARNING: rGE for half siblings only approximated by weighing 'w' by 'halfsib.cor', which is set to ", halfsib.cor, " \n")}
  }

  if (!cor.preferences) {preferences.mz.dz <- F} # Can only be free if correlated preferences are included


  #####################################################
  ####           Creates Free Parameters           ####
  #####################################################

  paths <- list(

    # Misc set parameters
    mxMatrix( type="Symm",  nrow=1, ncol=1, free=FALSE, values=halfsib.cor, label="HS_c_value",   name="HS_c" ),   # Half-sibling C-correlation
    mxMatrix( type="Symm",  nrow=1, ncol=1, free=FALSE, values=U.f,         labels="U1",          name="U" ) ,     # Closeness to intergenerational equilibirum (AM)
    mxMatrix( type="Symm",  nrow=1, ncol=1, free=FALSE, values=MZprop,      label="MZprop_value", name="MZprop" ), # Proportion MZs in unknown group


    # ACTE (parental generation, focal phenotype)
    mxMatrix( type="Symm", nrow=1, ncol=1, free=a1.free,  values=a1.val,lbound=0, label="est_a1",   name="a1" ), # Genetic influences
    mxMatrix( type="Symm", nrow=1, ncol=1, free=c1.free,  values=c1.val,lbound=0, label="est_c1",   name="c1" ), # (Sibling) shared environmental influences
    mxMatrix( type="Symm", nrow=1, ncol=1, free=t1.free,  values=t1.val,lbound=0, label="est_t1",   name="t1" ), # Twin shared environmental influences
    mxMatrix( type="Symm", nrow=1, ncol=1, free=e1.free,  values=e1.val,lbound=0, label="est_e1",   name="e1" ), # Unique environmental influences

    # Intergenerational paths
    mxMatrix( type="Symm", nrow=1, ncol=1, free=a1p.free,  values=0, label="est_a1p",  name="a1p"),
    mxMatrix( type="Symm", nrow=1, ncol=1, free=  p.free,  values=0, label="est_p",    name="p"  ),
    mxMatrix( type="Symm", nrow=1, ncol=1, free=c1p.free,  values=0, label="est_c1p",  name="c1p"),

    # ACE (child generation)
    mxMatrix( type="Symm", nrow=1, ncol=1, free=a2.free,  values=0, label="est_a2",   name="a2" ),
    mxMatrix( type="Symm", nrow=1, ncol=1, free=c2.free,  values=0, label="est_c2",   name="c2" ),
    mxMatrix( type="Symm", nrow=1, ncol=1, free=c2.cousin,values=0,lbound=-.0001,ubound=1, label="est_rc",   name="rc" ), # Correlation between cousins' residual shared environment
    mxMatrix( type="Symm", nrow=1, ncol=1, free=e2.free,  values=1, label="est_e2",   name="e2" )

  )


  #####################################################
  ####    Creates Assortative Mating Components    ####
  #####################################################

  # Does checks and fills NULL-objects
  if (a1s.free & c1s.free & t1s.free & e1s.free & !VarS_constraint) {stop("Model cannot be run if all paths to sorting factor are free without constraining the variance of the sorting factor.
                                                                         \nFix one of the paths (e.g., 'e1s.free = FALSE'), or select 'VarS_constraint = TRUE'")}
  if (indirect.assortment) {
    # The following code adds labels and starting values for the paths to the sorting factor.
    # If the path is not freely estimated, it either constrains it to the value given in the function arguments
    # If no argument is given, it is constrained to equal the corresponding path to the phenotype
    if(!a1s.free & is.null(a1s.val)) {lab_a1s <- "est_a1"; a1s.free <- a1.free; a1s.val = a1.val} else {lab_a1s <- "est_a1s"; if (is.null(a1s.val)) {a1s.val = a1.val} }
    if(!c1s.free & is.null(c1s.val)) {lab_c1s <- "est_c1"; c1s.free <- c1.free; c1s.val = c1.val} else {lab_c1s <- "est_c1s"; if (is.null(c1s.val)) {c1s.val = c1.val} }
    if(!t1s.free & is.null(t1s.val)) {lab_t1s <- "est_t1"; t1s.free <- t1.free; t1s.val = t1.val} else {lab_t1s <- "est_t1s"; if (is.null(t1s.val)) {t1s.val = t1.val} }
    if(!e1s.free & is.null(e1s.val)) {lab_e1s <- "est_e1"; e1s.free <- e1.free; e1s.val = e1.val} else {lab_e1s <- "est_e1s"; if (is.null(e1s.val)) {e1s.val = e1.val} }

  } else {
    # If no indirect assortment, then all paths so sorting factor must equal the corresponding path to the focal phentotype
    lab_a1s <- "est_a1"; a1s.free <- a1.free; a1s.val = a1.val
    lab_c1s <- "est_c1"; c1s.free <- c1.free; c1s.val = c1.val
    lab_t1s <- "est_t1"; t1s.free <- t1.free; t1s.val = t1.val
    lab_e1s <- "est_e1"; e1s.free <- e1.free; e1s.val = e1.val
    VarS_constraint == FALSE
  }


  paths <- append(paths, list(

    # ACTE (parental generation, sorting factor)
    mxMatrix( type="Symm", nrow=1, ncol=1, free=a1s.free,lbound=0,  values=a1s.val, label=lab_a1s,  name="a1s" ),
    mxMatrix( type="Symm", nrow=1, ncol=1, free=c1s.free,lbound=0,  values=c1s.val, label=lab_c1s,  name="c1s" ),
    mxMatrix( type="Symm", nrow=1, ncol=1, free=t1s.free,lbound=0,  values=t1s.val, label=lab_t1s,  name="t1s" ),
    mxMatrix( type="Symm", nrow=1, ncol=1, free=e1s.free,lbound=0,  values=e1s.val, label=lab_e1s,  name="e1s" ),

    # Variance of sorting factor, and covariance between focal phenotype and sorting factor
    mxAlgebra(expression = a1s^2 + c1s^2 + e1s^2 + 2*a1s*w*c1s + t1s^2,   name="VarS" ),
    mxAlgebra(expression = a1*a1s + c1*c1s + e1*e1s + c1*w*a1s + a1*w*c1s + t1*t1s,   name="covPS" ),
    mxAlgebra(expression = mu*(a1s + c1s*w)^2                       , name="partner_rG" ),

    # Estimates partner covariance on sorting factor, and rescales to copath
    mxMatrix( type="Lower", nrow=1, ncol=1, free=mu.free, values=mu.val, labels="mu1", name="mu_param" ),
    mxAlgebra( expression= mu_param/VarS^2 , name="mu" ),


    # Parameters for adjusting for sibling-correlated mate-preferences
    mxMatrix( type="Lower", nrow=1, ncol=1, free=cor.preferences,   values=0, labels="z1",  name="Z_dz" ), # Sibling-correlated mate-preferences (require partners of siblings-in-law, and tons of power)
    mxMatrix( type="Lower", nrow=1, ncol=1, free=preferences.mz.dz, values=0, labels="z2",  name="Z_mz" ), # MZ-correlated mate-preferences (require partners of siblings-in-law, and tons of power)
    mxMatrix( type="Lower", nrow=1, ncol=1, free=preferences.mz.dz, values=0, labels="z2",  name="Z_uz" ), # MZ-correlated mate-preferences (require partners of siblings-in-law, and tons of power)
    mxMatrix( type="Lower", nrow=1, ncol=1, free=preferences.mz.dz, values=0, labels="z2",  name="Z_hs" ), # MZ-correlated mate-preferences (require partners of siblings-in-law, and tons of power)
    mxAlgebra( expression=  Z_dz^2, name="Zdz" ),       # Sibling-correlated mate-preferences
    mxAlgebra( expression=  Zdz + Z_mz^2, name="Zmz" ), # MZ-correlated mate-preferences
    mxAlgebra( expression=  Zdz + Z_uz^2, name="Zuz" ), # UZ-correlated mate-preferences
    mxAlgebra( expression=  Zdz - Z_hs^2, name="Zhs" ) # HS-correlated mate-preferences

  )   )

  # Adds constraint if desired
  if (VarS_constraint) {VarS_cons <- mxConstraint(expression = VarS == var1, name = "VarS_constraint")} else {VarS_cons <- NULL}



  #####################################################
  ####        Creates Variance Components          ####
  #####################################################

  variances <- list(


    # Parental Generation (focal phenotype)
    mxAlgebra( expression = a1^2                    , name="VA1" ),
    mxAlgebra( expression = c1^2                    , name="VC1" ),
    mxAlgebra( expression = t1^2                    , name="VT1" ),
    mxAlgebra( expression = 2*c1*w*a1               , name="VrGE1" ),
    mxAlgebra( expression = e1^2                    , name="VE1" ),

    # Parental Generation (sorting factor)
    mxAlgebra( expression = a1s^2                    , name="VA1_s" ),
    mxAlgebra( expression = c1s^2                    , name="VC1_s" ),
    mxAlgebra( expression = t1s^2                    , name="VT1_s" ),
    mxAlgebra( expression = 2*c1s*w*a1s              , name="VrGE1_s" ),
    mxAlgebra( expression = e1s^2                    , name="VE1_s" ),

    # Offspring Generation
    mxAlgebra( expression =
                 2*(a1p*0.5*0.5*a1p) + # Genes correlated with parents
                 2*(a1p*0.5*partner_rG*0.5*a1p) + # AM-induced variance
                 a1p*(1 - FS_rG)*a1p, # Residual/Recombination Variance (0.5 - AM-induced variance)

               # All this is just a complicated way to ensure that the genetic variance is 1 under equilibirum, and that
               # the explained by VA1P (this variable) equals 1 * a1p^2
               # In other words, I could have written a1p^2 instead.
               # However, writing it this way allows the option of not assuming equilibrium,
               # which under assortative mating increases the variance of VA1P to be >1* a1p^2
               # If assuming base population instead (i.e., first generation of AM), then FS_rG should be 0.5,
               # and thus the recombination variance works out to be 0.5 as well.
               # This increases the variance to be more than one, resulting in VA1P > a1p^2.
               # Note that the interpretation of a1p is slightly different depending on what is assumed,
               # although offspring heritability should still be ((VA1P + VA2) / var1)

               name="VA1P" ),


    # Family-related environmental variance in offspring generation
    mxAlgebra( expression =
                 2*(p*var1*p) +
                 2*(c1p^2) +
                 4*(c1p*c1*p) +
                 4*(c1p*w*a1*p) +
                 2*((p*covPS + c1p*(c1s + w*a1s)) * mu * (covPS*p + (c1s + w*a1s)*c1p)),
               name="VF"  ),


    mxAlgebra( expression = a2^2               , name="VA2" ),
    mxAlgebra( expression = c2^2               , name="VC2" ),
    mxAlgebra( expression = e2^2               , name="VE2" ),

    # G-E covariance
    mxAlgebra( expression =
                 4*(a1p*0.5*(a1 + w*c1)*p) +
                 4*(a1p*0.5*w*c1p) +
                 4*(a1p*0.5*(a1s + w*c1s)*mu*(covPS*p + (c1s + a1s*w)*c1p)), name="VrGE2")
  )



  #######################################################
  ####        Creates Derived Parameters             ####
  #######################################################

  # Genetic correlation between siblings and half-siblings
  FS_rG <- mxAlgebra( expression = (1 + U*partner_rG)/2                                                , name="FS_rG" )
  HS_rG <- mxAlgebra( expression = 0.25 + U*((2*mu*(a1s + c1s*w)^2 + ((mu)^2)*((a1s + c1s*w)^2))/4)    , name="HS_rG" )     #TODO Double-check HS... I haven't worked it through thoroughly!!
  UZ_rG <- mxAlgebra( expression = (1-MZprop)*FS_rG + MZprop                                           , name="UZ_rG" )





  # TODO: Simplify by making an object for offspring rGE, and refer to that
  if (is.logical(parental_rGE)) {
    if (rGE_twostep) { rGE1 <- mxMatrix( type="Symm", nrow=1, ncol=1, free=FALSE, values=0, lbound=rGE_lbound, label="w_est",   name="w" )
    } else if (parental_rGE) {
      rGE1 <- list(
        mxMatrix( type="Symm", nrow=1, ncol=1, free=T, values=0,lbound=rGE_lbound, label="w_est",   name="w" ),
        mxConstraint( expression =  w ==  ((2*(a1p*0.5*(a1 + w*c1)*p)
                                            + 2*(a1p*0.5*w*c1p)
                                            + 2*((p*covPS + c1p*(c1s + w*a1s))* mu*(a1s+w*c1s)*0.5*a1p)) # G-E covariance inflation due to AM
                                           /(sqrt(VA1P+VA2)*sqrt(VF+VC2)))
                      , name="rGE_constraint" ))
    } else { rGE1 <- mxMatrix( type="Symm", nrow=1, ncol=1, free=FALSE, values=0, label="w_est",   name="w" ) }
  } else if (is.numeric(parental_rGE)) {
    if (parental_rGE < 1 & parental_rGE > -1) { rGE1 <- mxMatrix( type="Symm", nrow=1, ncol=1, free=FALSE, values=parental_rGE, label="w_est",   name="w" )
    } else { stop("parental_rGE must be either logical or numeric between -1 and 1") }
  } else { stop("parental_rGE must be either logical or numeric between -1 and 1")  }




  #######################################################
  ####    Creates Algebra for Expected Variances     ####
  #######################################################
  var1     <- mxAlgebra( expression= VA1 + VC1 + VrGE1 + VE1 + VT1, name="var1" )

  var2     <- mxAlgebra( expression=
                           VA2 + VC2 + VE2           # Child-specific variance
                         + VA1P + VF                 # Child-variance shared with parental phenotype
                         + VrGE2                     # G-E covariance
                         , name="var2" )


  ######################################################
  ####   Creates Algebra for Expected Covariances    ####
  #######################################################







  ####################################
  ####    Parental Generation     ####
  ####################################

  # Between twins (or siblings or however the related pair in the parental generation is related)
  covSibs_MZ     <- mxAlgebra(expression =     1*VA1 + VC1 + VrGE1 + VT1,     name = "covSibs_MZ")  # Covariation between siblings
  covSibs_DZ     <- mxAlgebra(expression = FS_rG*VA1 + VC1 + VrGE1 + VT1,     name = "covSibs_DZ")  # Covariation between siblings
  covSibs_UZ     <- mxAlgebra(expression = UZ_rG*VA1 + VC1 + VrGE1 + VT1,     name = "covSibs_UZ")  # Covariation between unknown twins
  covSibs_FS     <- mxAlgebra(expression = FS_rG*VA1 + VC1 + VrGE1,           name = "covSibs_FS")  # Covariation between siblings
  covSibs_HS     <- mxAlgebra(expression = HS_rG*VA1 + HS_c*VC1 + HS_c*VrGE1, name = "covSibs_HS") # Covariation between siblings # TODO: Doublecheck rGE for half sibs...

  covSibIn_MZ  <- mxAlgebra(expression = covPS*mu* (    1*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + t1*t1s      )  , name = "covSibIn_MZ") # Covariation between sibling-inlaw
  covSibIn_DZ  <- mxAlgebra(expression = covPS*mu* (FS_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + t1*t1s      )  , name = "covSibIn_DZ") # Covariation between sibling-inlaw
  covSibIn_UZ  <- mxAlgebra(expression = covPS*mu* (UZ_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + t1*t1s      )  , name = "covSibIn_UZ") # Covariation between sibling-inlaw
  covSibIn_FS  <- mxAlgebra(expression = covPS*mu* (FS_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s               )  , name = "covSibIn_FS") # Covariation between sibling-inlaw
  covSibIn_HS  <- mxAlgebra(expression = covPS*mu* (HS_rG*a1*a1s + c1*HS_c*c1s + HS_c*c1s*w*a1 + HS_c*c1*w*a1s)  , name = "covSibIn_HS") # Covariation between sibling-inlaw

  covInIn_MZ  <- mxAlgebra(expression = covPS^2 * ( mu^2 *  (    1*a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2  ) + Zmz ), name = "covInIn_MZ") # Covariation between co-inlaws
  covInIn_DZ  <- mxAlgebra(expression = covPS^2 * ( mu^2 *  (FS_rG*a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2  ) + Zdz ), name = "covInIn_DZ") # Covariation between co-inlaws
  covInIn_UZ  <- mxAlgebra(expression = covPS^2 * ( mu^2 *  (UZ_rG*a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2  ) + Zuz ), name = "covInIn_UZ") # Covariation between co-inlaws
  covInIn_FS  <- mxAlgebra(expression = covPS^2 * ( mu^2 *  (FS_rG*a1s^2 + c1s^2 + 2*c1s*w*a1s          ) + Zdz ), name = "covInIn_FS") # Covariation between co-inlaws
  covInIn_HS  <- mxAlgebra(expression = covPS^2 * ( mu^2 *  (HS_rG*a1s^2 + HS_c*c1s^2 + HS_c*2*c1s*w*a1s) + Zhs ), name = "covInIn_HS") # Covariation between co-inlaws

  covMATE <- mxAlgebra(expression = covPS*mu*covPS, name="covMate") # Covariance between partners



  ####################################
  ####       Nuclear Family       ####
  ####################################

  # Parent offspring covariance (w/o assortative mating)
  covPO   <- mxAlgebra( expression =           # Parent - Offspring Covariance
                          var1 * p             # Phenotypic Transmission
                        + a1 * 0.5 * a1p       # Genetic Confounding
                        + c1 * c1p             # Extended Family Confounding
                        + a1 * w * c1p         # rGE
                        + c1 * w * 0.5 * a1p   # rGE
                        ,name="covPO" )


  # Covariance between parental sorting factor and offspring phenotype
  SortPO <- mxAlgebra(expression = covPS*p + (a1s + c1s*w)*0.5*a1p + (c1s + a1s*w)*c1p, name = "SortPO")


  # Parent offspring covariance, including effects of assortative mating
  covPO_am <- mxAlgebra(expression = covPO + covPS*mu*SortPO, name = "covPO_am")





  # Covariance etween Siblings
  covSIB <- mxAlgebra(expression =
                        0.5 * VA2

                      + 2*(a1p*0.5*0.5*a1p) # Regular correlation via shared genes
                      + 2*(a1p*0.5*partner_rG*0.5*a1p) # AM-induced covariance via shared genes

                      + VC2 + VF
                      + VrGE2         # G-E  correlations
                      ,name="covSib" )


  ####################################
  ####          Avuncular         ####
  ####################################

  covAV_MZ <- mxAlgebra(expression =
                          a1 * 1 * 0.5 * a1p +  # Genetic relatedness
                          covSibs_MZ * p +      # Phenotypic via parent
                          c1 * c1p +            # Indirect environmental
                          c1 * w * 0.5 * a1p +   # rGE
                          a1 * w * c1p           # rGE
                        ,name="covAV_MZ" )

  covAV_DZ <- mxAlgebra(expression =
                          a1 * FS_rG * 0.5 * a1p +
                          covSibs_DZ * p +
                          c1 * c1p +
                          c1 * w * 0.5 * a1p +   # rGE
                          a1 * w * c1p          # rGE
                        ,name="covAV_DZ" )

  covAV_UZ <- mxAlgebra(expression =
                          a1 * UZ_rG * 0.5 * a1p +
                          covSibs_UZ * p +
                          c1 * c1p +
                          c1 * w * 0.5 * a1p +   # rGE
                          a1 * w * c1p,                 # rGE
                        name="covAV_UZ" )

  covAV_FS <- mxAlgebra(expression =
                          a1 * FS_rG * 0.5 * a1p +
                          covSibs_FS * p +
                          c1 * c1p +
                          c1 * w * 0.5 * a1p +   # rGE
                          a1 * w * c1p          # rGE
                        ,name="covAV_FS" )

  covAV_HS <- mxAlgebra(expression =
                          a1 * HS_rG * 0.5 * a1p +
                          covSibs_HS * p +
                          c1 * c1p * HS_c +

                          c1 * w * HS_c * 0.5 * a1p +   # rGE
                          a1 * w * HS_c * c1p           # rGE
                        ,name="covAV_HS" )




  # Accounts for assortative mating (uncle > parent > co-parent > child)
  covAVam_MZ <- mxAlgebra(expression = covAV_MZ + (    1*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + t1*t1s       ) * mu * SortPO, name="covAVam_MZ" )
  covAVam_DZ <- mxAlgebra(expression = covAV_DZ + (FS_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + t1*t1s       ) * mu * SortPO, name="covAVam_DZ" )
  covAVam_UZ <- mxAlgebra(expression = covAV_UZ + (UZ_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + t1*t1s       ) * mu * SortPO, name="covAVam_UZ" )
  covAVam_FS <- mxAlgebra(expression = covAV_FS + (FS_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s                ) * mu * SortPO, name="covAVam_FS" )
  covAVam_HS <- mxAlgebra(expression = covAV_HS + (HS_rG*a1*a1s + c1*HS_c*c1s + HS_c*c1s*w*a1 + HS_c*c1*w*a1s ) * mu * SortPO, name="covAVam_HS" )

  # Avuncular In-law Correlations
  covAVIn_MZ <- mxAlgebra(expression = covPS*mu*  ((a1s*    1 + c1s*w)*0.5* a1p + (c1s + a1s*w)*c1p + p*(    1*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + t1*t1s       ) + ((    1*a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2        ) * mu * SortPO)) + SortPO*Zmz*covPS , name="covAVIn_MZ" )
  covAVIn_DZ <- mxAlgebra(expression = covPS*mu*  ((a1s*FS_rG + c1s*w)*0.5* a1p + (c1s + a1s*w)*c1p + p*(FS_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + t1*t1s       ) + ((FS_rG*a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2        ) * mu * SortPO)) + SortPO*Zdz*covPS , name="covAVIn_DZ" )
  covAVIn_UZ <- mxAlgebra(expression = covPS*mu*  ((a1s*UZ_rG + c1s*w)*0.5* a1p + (c1s + a1s*w)*c1p + p*(UZ_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + t1*t1s       ) + ((UZ_rG*a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2        ) * mu * SortPO)) + SortPO*Zuz*covPS , name="covAVIn_UZ" )
  covAVIn_FS <- mxAlgebra(expression = covPS*mu*  ((a1s*FS_rG + c1s*w)*0.5* a1p + (c1s + a1s*w)*c1p + p*(FS_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s                ) + ((FS_rG*a1s^2 + c1s^2 + 2*c1s*w*a1s                ) * mu * SortPO)) + SortPO*Zdz*covPS , name="covAVIn_FS" )
  covAVIn_HS <- mxAlgebra(expression = covPS*mu*  ((a1s*HS_rG + c1s*w)*0.5* a1p + (c1s + a1s*w)*c1p + p*(HS_rG*a1*a1s + c1*HS_c*c1s + HS_c*c1s*w*a1 + HS_c*c1*w*a1s ) + ((HS_rG*a1s^2 + HS_c*c1s^2 + HS_c*2*c1s*w*a1s ) * mu * SortPO))      + SortPO*Zhs*covPS , name="covAVIn_HS" )












  ####################################
  ####     Cousin Covariances    ####
  ####################################



  covCous_MZ    <- mxAlgebra(expression =
                               a1p * 0.5 * 1 * 0.5 * a1p          # Genetic Covaraince
                             + 1 * 0.25 * VA2          # Genetic Covaraince

                             + c1p * c1p                  # Passive Environmental Covariance
                             + p * covSibs_MZ * p              # Phenotypic Transmission Covariance
                             + rc*VC2                         # Shared environmental factors unrelated to parental phenotype

                             + 2*(p * a1 * 1 * 0.5 * a1p)      # Phenotypic*Genetic Covariance
                             + 2*(p * c1 * c1p)                # Phenotypic*Passive Environ Covariance

                             + 2*(p * a1 * w * c1p)              #rGE
                             + 2*(p * c1 * w * 0.5 * a1p)         #rGE
                             + 2*(c1p * w * 0.5 * a1p)           #rGE


                             + 2*SortPO*mu*  ((a1s*1 + c1s*w)*0.5* a1p + (c1s + a1s*w)*c1p + p*(    1*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + t1*t1s        ))
                             + SortPO^2 * (mu^2 *  (    1*a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2         ) + Zmz )# Assortative Mating (via co-in-laws)

                             , name="covCous_MZ" )


  covCous_DZ    <- mxAlgebra(expression =
                               a1p * 0.5 * FS_rG * 0.5 * a1p        # Genetic Covaraince
                             + 0.5 * 0.25 * VA2        # Genetic Covaraince

                             + c1p * c1p                  # Passive Environmental Covariance
                             + p * covSibs_DZ * p              # Phenotypic Transmission Covariance
                             + rc*VC2                         # Shared environmental factors unrelated to parental phenotype

                             + 2*(p * a1 * FS_rG * 0.5 * a1p)    # Phenotypic*Genetic Covariance
                             + 2*(p * c1 * c1p)                # Phenotypic*Passive Environ Covariance

                             + 2*(p * a1 * w * c1p)            #rGE
                             + 2*(p * c1 * w * 0.5 * a1p)         #rGE
                             + 2*(c1p * w * 0.5 * a1p)           #rGE

                             + 2*SortPO*mu*  ((a1s*FS_rG + c1s*w)*0.5* a1p + (c1s + a1s*w)*c1p + p*(FS_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + t1*t1s               ))
                             + SortPO^2 * (mu^2 *  (FS_rG*a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2          ) + Zdz ) # Assortative Mating (via co-in-laws)

                             , name="covCous_DZ" )


  covCous_UZ    <- mxAlgebra(expression =
                               a1p * 0.5 * UZ_rG * 0.5 * a1p        # Genetic Covaraince
                             + ((1-MZprop)*0.5 + MZprop) * 0.25 * VA2                     # Genetic Covaraince

                             + c1p * c1p                  # Passive Environmental Covariance
                             + p * covSibs_UZ * p              # Phenotypic Transmission Covariance
                             + rc*VC2                         # Shared environmental factors unrelated to parental phenotype

                             + 2*(p * a1 * UZ_rG * 0.5 * a1p)    # Phenotypic*Genetic Covariance
                             + 2*(p * c1 * c1p)                # Phenotypic*Passive Environ Covariance

                             + 2*(p * a1 * w * c1p)            #rGE
                             + 2*(p * c1 * w * 0.5 * a1p)         #rGE
                             + 2*(c1p * w * 0.5 * a1p)           #rGE

                             +  2*SortPO*mu*  ((a1s*UZ_rG + c1s*w)*0.5* a1p + (c1s + a1s*w)*c1p + p*(UZ_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s  + t1*t1s           ))
                             + SortPO^2 *( mu^2 *  (UZ_rG*a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2          ) + Zuz ) # Assortative Mating (via co-in-laws)
                             , name="covCous_UZ" )


  covCous_FS    <- mxAlgebra(expression =
                               a1p * 0.5 * FS_rG * 0.5 * a1p        # Genetic Covaraince
                             + 0.5 * 0.25 * VA2        # Genetic Covaraince

                             + c1p * c1p                  # Passive Environmental Covariance
                             + p * covSibs_FS * p              # Phenotypic Transmission Covariance
                             + rc*VC2                         # Shared environmental factors unrelated to parental phenotype

                             + 2*(p * a1 * FS_rG * 0.5 * a1p)    # Phenotypic*Genetic Covariance
                             + 2*(p * c1 * c1p)                # Phenotypic*Passive Environ Covariance

                             + 2*(p * a1 * w * c1p)            #rGE
                             + 2*(p * c1 * w * 0.5 * a1p)         #rGE
                             + 2*(c1p * w * 0.5 * a1p)           #rGE

                             + 2*SortPO*mu*  ((a1s*FS_rG + c1s*w)*0.5* a1p + (c1s + a1s*w)*c1p + p*(FS_rG*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s                ))
                             + SortPO^2 * (mu^2 *  (FS_rG*a1s^2 + c1s^2 + 2*c1s*w*a1s          ) + Zdz ) # Assortative Mating (via co-in-laws)

                             , name="covCous_FS" )


  covCous_HS    <- mxAlgebra(expression =
                               a1p * 0.5 * HS_rG * 0.5 * a1p  # Genetic Covaraince
                             + 0.25 * 0.25 * VA2              # Genetic Covaraince

                             + c1p * c1p * HS_c               # Passive Environmental Covariance
                             + p * covSibs_HS * p              # Phenotypic Transmission Covariance
                             + rc*VC2                         # Shared environmental factors unrelated to parental phenotype

                             + 2*(p * a1 * HS_rG * 0.5 * a1p)   # Phenotypic*Genetic Covariance
                             + 2*(p * c1 * c1p * HS_c)                # Phenotypic*Passive Environ Covariance

                             + 2*(p * a1 * w * HS_c * c1p)            #rGE
                             + 2*(p * c1 * w * HS_c  * 0.5 * a1p)         #rGE
                             + 2*(c1p * w * HS_c * 0.5 * a1p)           #rGE

                             + 2*SortPO*mu*  ((a1s*HS_rG + c1s*w)*0.5* a1p + (c1s + a1s*w)*c1p + p*(HS_rG*a1*a1s + c1*HS_c*c1s + HS_c*c1s*w*a1 + HS_c*c1*w*a1s ))
                             + SortPO^2 * (mu^2 *  (HS_rG*a1s^2 + HS_c*c1s^2 + HS_c*2*c1s*w*a1s)  + Zhs )# Assortative Mating (via co-in-laws)

                             , name="covCous_HS" )












  #######################################################
  ####     Creates Expected Covariance Matrices      ####
  #######################################################

  # Create Expected Covariance Matrix             Sib1          Sib2     Spouse1     Spouse1       CH1_FAM1     CH2_FAM1     CH1_FAM2      CH2_FAM2
  expCov_MZ  <- mxAlgebra( expression=
                             rbind( cbind(       var1,   covSibs_MZ,    covMate,covSibIn_MZ,      covPO_am,      covPO_am,    covAVam_MZ,    covAVam_MZ), #FAM1
                                    cbind( covSibs_MZ,         var1,covSibIn_MZ,    covMate,    covAVam_MZ,    covAVam_MZ,      covPO_am,      covPO_am), #FAM2
                                    cbind(    covMate,  covSibIn_MZ,       var1, covInIn_MZ,      covPO_am,      covPO_am,    covAVIn_MZ,    covAVIn_MZ), #FAM1_spouse
                                    cbind(covSibIn_MZ,      covMate, covInIn_MZ,       var1,    covAVIn_MZ,    covAVIn_MZ,      covPO_am,      covPO_am), #FAM2_spouse

                                    cbind(   covPO_am,   covAVam_MZ,   covPO_am, covAVIn_MZ,          var2,        covSib,    covCous_MZ,    covCous_MZ),
                                    cbind(   covPO_am,   covAVam_MZ,   covPO_am, covAVIn_MZ,        covSib,          var2,    covCous_MZ,    covCous_MZ),
                                    cbind( covAVam_MZ,     covPO_am, covAVIn_MZ,   covPO_am,    covCous_MZ,    covCous_MZ,          var2,        covSib),
                                    cbind( covAVam_MZ,     covPO_am, covAVIn_MZ,   covPO_am,    covCous_MZ,    covCous_MZ,        covSib,         var2)),
                           name="expCov_MZ" )


  expCov_DZ  <- mxAlgebra( expression=
                             rbind( cbind(       var1,   covSibs_DZ,    covMate,covSibIn_DZ,      covPO_am,      covPO_am,    covAVam_DZ,    covAVam_DZ), #FAM1
                                    cbind( covSibs_DZ,         var1,covSibIn_DZ,    covMate,    covAVam_DZ,    covAVam_DZ,      covPO_am,      covPO_am), #FAM2
                                    cbind(    covMate,  covSibIn_DZ,       var1, covInIn_DZ,      covPO_am,      covPO_am,    covAVIn_DZ,    covAVIn_DZ), #FAM1_spouse
                                    cbind(covSibIn_DZ,      covMate, covInIn_DZ,       var1,    covAVIn_DZ,    covAVIn_DZ,      covPO_am,      covPO_am), #FAM2_spouse

                                    cbind(   covPO_am,   covAVam_DZ,   covPO_am, covAVIn_DZ,          var2,        covSib,    covCous_DZ,    covCous_DZ),
                                    cbind(   covPO_am,   covAVam_DZ,   covPO_am, covAVIn_DZ,        covSib,          var2,    covCous_DZ,    covCous_DZ),
                                    cbind( covAVam_DZ,     covPO_am, covAVIn_DZ,   covPO_am,    covCous_DZ,    covCous_DZ,          var2,        covSib),
                                    cbind( covAVam_DZ,     covPO_am, covAVIn_DZ,   covPO_am,    covCous_DZ,    covCous_DZ,        covSib,         var2)),
                           name="expCov_DZ" )



  if (use.unknownzyg) {
    expCov_UZ  <- mxAlgebra( expression=
                               rbind( cbind(       var1,   covSibs_UZ,    covMate,covSibIn_UZ,      covPO_am,      covPO_am,    covAVam_UZ,    covAVam_UZ), #FAM1
                                      cbind( covSibs_UZ,         var1,covSibIn_UZ,    covMate,    covAVam_UZ,    covAVam_UZ,      covPO_am,      covPO_am), #FAM2
                                      cbind(    covMate,  covSibIn_UZ,       var1, covInIn_UZ,      covPO_am,      covPO_am,    covAVIn_UZ,    covAVIn_UZ), #FAM1_spouse
                                      cbind(covSibIn_UZ,      covMate, covInIn_UZ,       var1,    covAVIn_UZ,    covAVIn_UZ,      covPO_am,      covPO_am), #FAM2_spouse

                                      cbind(   covPO_am,   covAVam_UZ,   covPO_am, covAVIn_UZ,          var2,        covSib,    covCous_UZ,    covCous_UZ),
                                      cbind(   covPO_am,   covAVam_UZ,   covPO_am, covAVIn_UZ,        covSib,          var2,    covCous_UZ,    covCous_UZ),
                                      cbind( covAVam_UZ,     covPO_am, covAVIn_UZ,   covPO_am,    covCous_UZ,    covCous_UZ,          var2,        covSib),
                                      cbind( covAVam_UZ,     covPO_am, covAVIn_UZ,   covPO_am,    covCous_UZ,    covCous_UZ,        covSib,         var2)),
                             name="expCov_UZ" )
  }

  if (use.fullsibs) {
    expCov_FS  <- mxAlgebra( expression=
                               rbind( cbind(       var1,   covSibs_FS,    covMate,covSibIn_FS,      covPO_am,      covPO_am,    covAVam_FS,    covAVam_FS), #FAM1
                                      cbind( covSibs_FS,         var1,covSibIn_FS,    covMate,    covAVam_FS,    covAVam_FS,      covPO_am,      covPO_am), #FAM2
                                      cbind(    covMate,  covSibIn_FS,       var1, covInIn_FS,      covPO_am,      covPO_am,    covAVIn_FS,    covAVIn_FS), #FAM1_spouse
                                      cbind(covSibIn_FS,      covMate, covInIn_FS,       var1,    covAVIn_FS,    covAVIn_FS,      covPO_am,      covPO_am), #FAM2_spouse

                                      cbind(   covPO_am,   covAVam_FS,   covPO_am, covAVIn_FS,          var2,        covSib,    covCous_FS,    covCous_FS),
                                      cbind(   covPO_am,   covAVam_FS,   covPO_am, covAVIn_FS,        covSib,          var2,    covCous_FS,    covCous_FS),
                                      cbind( covAVam_FS,     covPO_am, covAVIn_FS,   covPO_am,    covCous_FS,    covCous_FS,          var2,        covSib),
                                      cbind( covAVam_FS,     covPO_am, covAVIn_FS,   covPO_am,    covCous_FS,    covCous_FS,        covSib,         var2)),
                             name="expCov_FS" )

  }

  if (use.halfsibs) {
    expCov_HS  <- mxAlgebra( expression=
                               rbind( cbind(       var1,   covSibs_HS,    covMate,covSibIn_HS,      covPO_am,      covPO_am,    covAVam_HS,    covAVam_HS), #FAM1
                                      cbind( covSibs_HS,         var1,covSibIn_HS,    covMate,    covAVam_HS,    covAVam_HS,      covPO_am,      covPO_am), #FAM2
                                      cbind(    covMate,  covSibIn_HS,       var1, covInIn_HS,      covPO_am,      covPO_am,    covAVIn_HS,    covAVIn_HS), #FAM1_spouse
                                      cbind(covSibIn_HS,      covMate, covInIn_HS,       var1,    covAVIn_HS,    covAVIn_HS,      covPO_am,      covPO_am), #FAM2_spouse

                                      cbind(   covPO_am,   covAVam_HS,   covPO_am, covAVIn_HS,          var2,        covSib,    covCous_HS,    covCous_HS),
                                      cbind(   covPO_am,   covAVam_HS,   covPO_am, covAVIn_HS,        covSib,          var2,    covCous_HS,    covCous_HS),
                                      cbind( covAVam_HS,     covPO_am, covAVIn_HS,   covPO_am,    covCous_HS,    covCous_HS,          var2,        covSib),
                                      cbind( covAVam_HS,     covPO_am, covAVIn_HS,   covPO_am,    covCous_HS,    covCous_HS,        covSib,         var2)),
                             name="expCov_HS" )
  }





  if (child_ordinal & parent_ordinal) {

    # Create Algebra for expected Mean Matrix
    meanG     <- mxMatrix( type="Full", nrow=1, ncol=8, free=F, values=0,
                           labels=c(rep("parent_mean",4),rep("child_mean",4)), name="meanG" )
    # # Create Algebra for expected Threshold Matrix
    threG     <- mxMatrix( type="Full", nrow=1, ncol=8, free=TRUE, values=0, labels=c(rep("thres_pa",4),rep("thres_ch",4)), name="threG")

    thres_vars <- selVars
    superflous <- NULL
  } else if (!child_ordinal & parent_ordinal) {

    # Create Algebra for expected Mean Matrix
    meanG     <- mxMatrix( type="Full", nrow=1, ncol=8, free=c(rep(F,4),rep(T,4)), values=0,
                           labels=c(rep("parent_mean",4),rep("child_mean",4)), name="meanG" )

    # # Create Algebra for expected Threshold Matrix
    threG     <- mxMatrix( type="Full", nrow=1, ncol=4, free=TRUE, values=0, labels="thres_pa", name="threG")

    thres_vars <- selVars_pa
    superflous <- mxMatrix( type="Symm", nrow=1, ncol=1, free=FALSE, values=0, label="thres_ch",   name="superflous" )


  } else if (child_ordinal & !parent_ordinal) {

    # Create Algebra for expected Mean Matrix
    meanG     <- mxMatrix( type="Full", nrow=1, ncol=8, free=c(rep(T,4),rep(F,4)), values=0,
                           labels=c(rep("parent_mean",4),rep("child_mean",4)), name="meanG" )

    # # Create Algebra for expected Threshold Matrix
    threG     <- mxMatrix( type="Full", nrow=1, ncol=4, free=TRUE, values=0, labels="thres_ch", name="threG")

    thres_vars <- selVars_ch
    superflous <- mxMatrix( type="Symm", nrow=1, ncol=1, free=FALSE, values=0, label="thres_pa",   name="superflous" )

  }



  if (any(child_ordinal,parent_ordinal)) {
    exp_MZ     <- mxExpectationNormal( covariance="expCov_MZ", means="meanG", thresholds = "threG", dimnames=selVars, threshnames = thres_vars)
    exp_DZ     <- mxExpectationNormal( covariance="expCov_DZ", means="meanG", thresholds = "threG", dimnames=selVars, threshnames = thres_vars)
    if (use.unknownzyg) {exp_UZ     <- mxExpectationNormal( covariance="expCov_UZ", means="meanG", thresholds = "threG", dimnames=selVars, threshnames = thres_vars)}
    if (use.fullsibs)   {exp_FS     <- mxExpectationNormal( covariance="expCov_FS", means="meanG", thresholds = "threG", dimnames=selVars, threshnames = thres_vars)}
    if (use.halfsibs)   {exp_HS     <- mxExpectationNormal( covariance="expCov_HS", means="meanG", thresholds = "threG", dimnames=selVars, threshnames = thres_vars)}
  } else  {
    # Create Algebra for expected Mean Matrix
    meanG     <- mxMatrix( type="Full", nrow=1, ncol=8, free=c(rep(T,4),rep(T,4)), values=0,
                           labels=c(rep("parent_mean",4),rep("child_mean",4)), name="meanG" )
    threG     <- NULL
    superflous <- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, values=0, label=c("thres_pa","thres_ch"),   name="superflous" )

    exp_MZ     <- mxExpectationNormal(covariance="expCov_MZ", means="meanG", dimnames=selVars)
    exp_DZ     <- mxExpectationNormal(covariance="expCov_DZ", means="meanG", dimnames=selVars)
    if (use.unknownzyg) {exp_UZ     <- mxExpectationNormal(covariance="expCov_UZ", means="meanG", dimnames=selVars)}
    if (use.fullsibs)   {exp_FS     <- mxExpectationNormal(covariance="expCov_FS", means="meanG", dimnames=selVars)}
    if (use.halfsibs)   {exp_HS     <- mxExpectationNormal(covariance="expCov_HS", means="meanG", dimnames=selVars)}
  }











  #######################################################
  ####          Creates Model Objects                ####
  #######################################################





  # Create Data Objects
  if (combine.FS_DZ) {
    data_DZ <- mxData( observed = subset(famData, grepl("DZ",fam.type) | grepl("FS",fam.type), selVars), type="raw" )
  } else {
    data_DZ    <- mxData( observed = subset(famData, grepl("DZ",fam.type), selVars), type="raw" )
    if (use.fullsibs)   {data_FS    <- mxData( observed = subset(famData, grepl("FS",fam.type), selVars), type="raw" )}
  }

  data_MZ    <- mxData( observed = subset(famData, grepl("MZ",fam.type), selVars), type="raw" )
  if (use.unknownzyg) {data_UZ    <- mxData( observed = subset(famData, grepl("UZ",fam.type), selVars), type="raw" )}
  if (use.halfsibs)   {data_HS    <- mxData( observed = subset(famData, grepl("HS",fam.type), selVars), type="raw" )}



  # Compile objects into lists (it eases handling)
  common_objects <- list(meanG, threG, superflous, SortPO,
                         var1,var2,
                         covMATE,covPO,covPO_am,covSIB, FS_rG, HS_rG, UZ_rG, rGE1)

  algebras_MZ <- list(covSibs_MZ, covSibIn_MZ, covInIn_MZ,covAV_MZ, covAVam_MZ,covAVIn_MZ,covCous_MZ)
  algebras_DZ <- list(covSibs_DZ, covSibIn_DZ, covInIn_DZ,covAV_DZ, covAVam_DZ,covAVIn_DZ,covCous_DZ)
  if (use.unknownzyg) {algebras_UZ <- list(covSibs_UZ, covSibIn_UZ, covInIn_UZ,covAV_UZ, covAVam_UZ,covAVIn_UZ,covCous_UZ)}
  if (use.fullsibs)   {algebras_FS <- list(covSibs_FS, covSibIn_FS, covInIn_FS,covAV_FS, covAVam_FS,covAVIn_FS,covCous_FS)}
  if (use.halfsibs)   {algebras_HS <- list(covSibs_HS, covSibIn_HS, covInIn_HS,covAV_HS, covAVam_HS,covAVIn_HS,covCous_HS)}



  if (fitFunction == "WLS") {
    if (any(child_ordinal,parent_ordinal)) {
      fitFun     <- mxFitFunctionWLS()
    } else {
      fitFun     <- mxFitFunctionWLS(allContinuousMethod='marginals')
    }
  } else if (fitFunction == "ML") {
    fitFun     <- mxFitFunctionML()
  } else {
    stop("fitFunction argument must be either 'WLS' or 'ML'")
  }


  model_MZ <- mxModel(paths, variances, common_objects, algebras_MZ, expCov_MZ, exp_MZ, data_MZ, fitFun, name = "MZ")
  model_DZ <- mxModel(paths, variances, common_objects, algebras_DZ, expCov_DZ, exp_DZ, data_DZ, fitFun, name = "DZ")
  if (use.unknownzyg) {model_UZ <- mxModel(paths, variances, common_objects, algebras_UZ, expCov_UZ, exp_UZ, data_UZ, fitFun, name = "UZ")}
  if (use.fullsibs)   {model_FS <- mxModel(paths, variances, common_objects, algebras_FS, expCov_FS, exp_FS, data_FS, fitFun, name = "FS")}
  if (use.halfsibs)   {model_HS <- mxModel(paths, variances, common_objects, algebras_HS, expCov_HS, exp_HS, data_HS, fitFun, name = "HS")}


  exp_cov_list <- list(expCov_MZ,expCov_DZ)
  model_list   <- list(model_MZ, model_DZ  )
  algebras_all <- list(algebras_MZ, algebras_DZ)
  groups <- c("MZ","DZ")
  if (use.fullsibs) {
    exp_cov_list <- append(exp_cov_list, expCov_FS)
    model_list   <- append(model_list  , model_FS)
    groups       <- append(groups      , "FS")
    algebras_all <- append(algebras_all, algebras_FS)
  }
  if (use.halfsibs) {
    exp_cov_list <- append(exp_cov_list, expCov_HS)
    model_list   <- append(model_list  , model_HS)
    groups       <- append(groups      , "HS")
    algebras_all <- append(algebras_all, algebras_HS)
  }
  if (use.unknownzyg) {
    exp_cov_list <- append(exp_cov_list, expCov_UZ)
    model_list   <- append(model_list  , model_UZ)
    groups       <- append(groups      , "UZ")
    algebras_all <- append(algebras_all, algebras_UZ)
  }
  multi <- mxFitFunctionMultigroup( groups )



  #######################################################
  ####          Creates CI - Objects                ####
  #######################################################

  # Create Confidence Interval Objects
  est_var_Col     <<- c("var1","var2",

                        "VA1","VC1","VT1","VE1","VrGE1",
                        "VA1_s","VC1_s","VT1_s","VE1_s","VrGE1_s",
                        "VA1P","VP","VC1P","VA2","VC2","VE2","VF","VF_cov","VrGE2", "Offspring Heritability", "Offspring Shared Environment",
                        "Gene-Environment Correlation","Genetic Correlation (par_pheno and off_pheno)","Partner Correlation (Sorting Factor)","Partner Correlation (Focal Phenotype)", "Partner rG","Sibling rG")
  est_var     <- mxAlgebra(expression=cbind(var1,
                                            var2,

                                            VA1/var1,
                                            VC1/var1,
                                            VT1/var1,
                                            VE1/var1,
                                            VrGE1/var1,

                                            VA1_s/VarS,
                                            VC1_s/VarS,
                                            VT1_s/VarS,
                                            VE1_s/VarS,
                                            VrGE1_s/VarS,

                                            VA1P/var2,
                                            (2*((p^2)*var1) + 2*(p*covPS*mu*covPS*p))/var2,
                                            (2*(c1p^2) + 2*(c1p*c1s*mu*c1s*c1p))/var2,
                                            VA2/var2,
                                            VC2/var2,
                                            VE2/var2,
                                            VF/var2,
                                            (VF-(2*(c1p^2) + 2*(c1p*c1s*mu*c1s*c1p))-(2*((p^2)*var1) + 2*(p*covPS*mu*covPS*p)))/var2,
                                            VrGE2/var2,
                                            (VA2+VA1P)/var2,
                                            (VC2+VF)/var2,

                                            w,
                                            (a1*a1p)/sqrt(VA1*(VA1P + VA2)),
                                            mu * VarS,
                                            covMate/var1,
                                            partner_rG,
                                            FS_rG),
                           name="est_var", dimnames=list("est",est_var_Col) )


  est_path_Col     <<- c("mu_param","a1s","c1s","t1s","e1s","w","a1","c1","t1","e1","a1p","p","c1p","a2","c2","e2","Z_mz","Z_dz",
                         "parent_mean","child_mean", "thres_pa","thres_pa_Z", "thres_ch","thres_ch_Z")
  est_path     <- mxAlgebra(expression=cbind(mu_param,
                                             a1s,
                                             c1s,
                                             t1s,
                                             e1s,
                                             w,
                                             a1,
                                             c1,
                                             t1,
                                             e1,
                                             a1p,
                                             p,
                                             c1p,
                                             a2,
                                             c2,
                                             e2,
                                             Z_dz,Z_mz,
                                             parent_mean,
                                             child_mean,
                                             thres_pa,
                                             thres_pa/sqrt(var1),
                                             thres_ch,
                                             thres_ch/sqrt(var2)),
                            name="est_path", dimnames=list("est",est_path_Col) )



  est_PO_Col     <<- c("Parent Offspring Correlation",
                       "Phenotypic",
                       "Passive Genetic",
                       "Passive Environmental")


  est_PO_raw     <- mxAlgebra(expression=cbind(covPO,
                                               var1 * p,
                                               (a1 + c1*w) * 0.5 * a1p,
                                               (c1 + a1*w) * c1p
  ) / (sqrt(var1)*sqrt(var2)) ,
  name="est_PO_raw" )


  est_PO_AMinf     <- mxAlgebra(expression=cbind(covPS*mu*SortPO,
                                                 covPS*mu*covPS*p,
                                                 covPS*mu*(a1s + c1s*w)*0.5*a1p,
                                                 covPS*mu*(c1s + a1s*w)*c1p
  )/(sqrt(var1)*sqrt(var2)),
  name="est_PO_AMinf" )


  est_PO_total     <- mxAlgebra(expression=  est_PO_raw + est_PO_AMinf , name="est_PO_total")

  est_PO <- mxAlgebra(rbind(est_PO_raw,est_PO_AMinf,est_PO_total),
                      name="est_PO", dimnames=list(c("raw", "AMinf", "total"),est_PO_Col))





  ciCOTS     <- mxCI( c("est_var", "est_PO") )

  CICOLNAMES <<- c(est_var_Col,paste(est_PO_Col, "(raw)"), paste(est_PO_Col, "(AMinf)"), paste(est_PO_Col, "(total)"))


  # Create implied correlations object for model fit plot
  imp_cors     <- mxAlgebra(expression=cbind(
    covMate       /var1,
    covSibs_MZ    /var1,
    covSibs_DZ    /var1,
    covSibs_FS    /var1,
    covSibIn_MZ   /var1,
    covSibIn_DZ   /var1,
    covSibIn_FS   /var1,
    covInIn_MZ    /var1,
    covInIn_DZ    /var1,
    covInIn_FS    /var1,

    covPO_am      /sqrt(var1*var2),
    covAVam_MZ    /sqrt(var1*var2),
    covAVam_DZ    /sqrt(var1*var2),
    covAVam_FS    /sqrt(var1*var2),
    covAVIn_MZ    /sqrt(var1*var2),
    covAVIn_DZ    /sqrt(var1*var2),
    covAVIn_FS    /sqrt(var1*var2),

    covSib        /var2,
    covCous_MZ    /var2,
    covCous_DZ    /var2,
    covCous_FS    /var2
  ), name="implied_correlations",
  dimnames=list("est",
                c("Partners",
                  "Monozygotic Twins", "Dizygotic Twins", "Full Siblings",
                  "sib-in-law (MZ)","sib-in-law (DZ)","sib-in-law (FS)",
                  "co-in-law (MZ)","co-in-law (DZ)","co-in-law (FS)",
                  "Parent-Offspring",
                  "avuncular (MZ)","avuncular (DZ)","avuncular (FS)",
                  "avuncular-in-law (MZ)","avuncular-in-law (DZ)","avuncular-in-law (FS)",
                  "Offspring Siblings",
                  "cousins (MZ)", "cousins (DZ)", "cousins (FS)")) )


  results_objects <- list(
    est_path,
    est_var,
    est_PO,
    est_PO_raw,est_PO_AMinf,est_PO_total,
    imp_cors
  )


  #######################################################
  ####             Builds the model !                ####
  #######################################################


  model_iamCOTS <- mxModel(paths, multi, results_objects, ciCOTS, VarS_cons,
                           common_objects, variances, algebras_all, model_list,exp_cov_list,
                           name = name)





  if (rGE_twostep) {
    if (is.logical(parental_rGE)) {

      if (parental_rGE) {
        if (RunModel) {
          if (TryHard) {
            if(parent_ordinal | child_ordinal) {
              fit_iamCOTS    <- mxTryHardOrdinal(model_iamCOTS, intervals = compute.intervals)
            } else {
              fit_iamCOTS    <- mxTryHard(model_iamCOTS, intervals = compute.intervals)
            }
          } else {
            fit_iamCOTS    <- mxRun(model_iamCOTS, intervals = compute.intervals)
          }
        } else {fit_iamCOTS <- model_iamCOTS}
        model_iamCOTS <- omxSetParameters(fit_iamCOTS,labels="w_est", free=T)
        model_iamCOTS <- mxModel(model_iamCOTS,
                                 mxConstraint( expression =
                                                 w ==
                                                 ((2*(a1p*0.5*a1*p)  + 2*(a1p*0.5*w*(c1p + c1*p))

                                                   # G-E covariance inflation due to AM
                                                   + 2*(p*covPS*    mu*(a1s+w*c1s)*0.5*a1p)  # Up through P
                                                   + 2*(c1p*w*a1s*  mu*(a1s+w*c1s)*0.5*a1p)   # Up through P
                                                   + 2*(c1p*c1s*    mu*(a1s+w*c1s)*0.5*a1p))/  # Up through C
                                                    (sqrt(VA1P+VA2)*sqrt(VF+VC2)))
                                               , name="rGE_constraint" ))
      }
    }

  }


  if (!RunModel) {
    return(model_iamCOTS)
  } else {
    if (TryHard) {
      if(parent_ordinal | child_ordinal) {
        fit_iamCOTS    <- mxTryHardOrdinal(model_iamCOTS, intervals = compute.intervals)
      } else {
        fit_iamCOTS    <- mxTryHard(model_iamCOTS, intervals = compute.intervals)
      }
    } else {
      fit_iamCOTS    <- mxRun(model_iamCOTS, intervals = compute.intervals)
    }
    return(fit_iamCOTS)
  }




}




