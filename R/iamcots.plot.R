
#' Plot iAM-COTS variance and parent-offspring decomposition results
#'
#' Creates a multi-panel summary plot for a fitted iAM-COTS OpenMx model.
#' Panels include (1) variance decomposition for the parental phenotype,
#' (2) decomposition of the parent-offspring correlation into components,
#' (3) variance decomposition for the offspring phenotype, and optionally
#' (4) variance decomposition for the parental sorting factor.
#'
#' Uncertainty is shown using confidence intervals when available. If the fit
#' contains OpenMx confidence intervals, those are used. Otherwise Wald-type
#' intervals are computed from standard errors, or bootstrap quantile intervals
#' can be used when the model was fitted with MxComputeBootstrap.
#'
#' If simulated values are provided, they are overlaid as additional points
#' on the relevant panels.
#'
#' @param fit A fitted OpenMx model (typically returned by iam.cots.fun) that
#'   contains est_var and est_PO results.
#' @param sim_vars Optional simulated variance-component results to overlay.
#'   Must be compatible with get.varcomps(), which expects it to be coercible
#'   via t(sim_vars) and to align with est_var column names.
#' @param sim_PO Optional simulated parent-offspring decomposition results to
#'   overlay. Must be compatible with get.POcor(), i.e., contain columns
#'   parameter, type, and estimate.
#' @param incl.sort Logical. If TRUE, include the sorting-factor variance
#'   decomposition panel.
#' @param boot.CI Logical. If TRUE and fit$compute is an MxComputeBootstrap
#'   object, use bootstrap quantile intervals instead of Wald-type intervals.
#'
#' @return A grob object returned by gridExtra::grid.arrange().
#'
#' @details
#' This function depends on ggplot2 and gridExtra. It also assumes the following
#' functions exist in the package namespace: get.POcor(), get.varcomps(), and
#' custom_theme(). The plotted parameters are hard-coded to specific est_var
#' and est_PO labels produced by the iAM-COTS model.
#'
#'
#' @export



iamcots.plot <- function(fit, sim_vars = NULL, sim_PO = NULL, incl.sort = T, boot.CI=F) {
  requireNamespace("tidyverse")
  requireNamespace("gridExtra")

  results_PO <- get.POcor(fit, sim_PO, boot.CI = boot.CI)
  results_vars <- get.varcomps(fit, sim_vars, boot.CI = boot.CI)

  p2 <-  ggplot(filter(results_PO , parameter != "Parent Offspring Correlation", type == "total"),
                aes(x=estimate, xmin=lbound,xmax=ubound,y=parameter, label=paste0(round(estimate,2)," (",round(lbound,2) ,"???",round(ubound,2) ,")"))) +
    geom_col(data = filter(results_PO, parameter != "Parent Offspring Correlation", type != "total"),
             aes(x=estimate,y=parameter, alpha = type, fill = parameter),color="black", inherit.aes = F,
    ) +
    geom_point(size=4) +
    geom_errorbar( width=0, linewidth =1) +
    geom_vline(xintercept = 0, linewidth=1.5) +
    geom_text(size = 3, hjust=-.2,vjust=-.5,
              fontface="bold",
              color="black") +

    custom_theme() +
    scale_fill_viridis_d( end=.8, guide="none", option = "C") +
    scale_alpha_manual(values = c(raw = 1, AMinf=.5), labels = c(raw = "From index parent", AMinf= "Via co-parent"), breaks = c("raw", "AMinf")) +
    labs(title = "Decomposition of Parent-Offspring Correlation",
         alpha = "",y="", x="Parent-Offspring Correlation") +
    #coord_cartesian(ylim=c(0,.3)) +
    theme(legend.position = "top")








  p1 <- filter(results_vars,
               parameter %in% c("VA1","VC1","VT1","VrGE1","VE1","VE1b","VE1w"),
               abs(estimate) > 0.001) %>%
    mutate(
      generation = factor(ifelse( parameter %in% c("VA1","VC1","VT1","VrGE1","VE1","VE1b","VE1w"), "Parent Phenotype", "Offspring Phenotype"), levels=c("Parent Phenotype", "Offspring Phenotype")),
      parameter = factor(parameter, levels = c("VA1P","VA1","VA2","VC1","VT1","VF","VC2","VE1","VE2","VrGE1","VrGE2","VE1b","VE1w"))) %>%
    ggplot(aes(x = parameter,
               y = estimate, ymin =lbound, ymax=ubound,
               fill = parameter,
               label=paste0(round(estimate*100),"%\n(",round(lbound*100),"%???",round(ubound*100),"%)")))  +
    geom_col(color="black") +
    geom_text(y = 1, size = 3, hjust=.5,vjust=1.1,
              fontface="bold",
              color="black") +
    geom_point(size=4) +
    geom_errorbar( width=0, linewidth =1) +

    scale_fill_viridis_d(end = .85, guide="none", option = "C") +
    scale_y_continuous(labels=function(x) paste0(round(x*100),"%")) +
    custom_theme() +
    geom_hline(yintercept = c(0,1), linewidth =1.5) +
    #guides(fill = guide_legend(override.aes = list(shape = NA)) ) +
    labs(title = "Variance Decomposition", subtitle = "Parental Phenotype",
         y= "",
         x = "",
         fill = ""
    )

  p3 <- filter(results_vars,
               parameter %in% c("VA1P","VA2","VF","VC2","VrGE2", "VE2"),
               abs(estimate) > 0.001) %>%
    mutate(
      generation = factor(ifelse( parameter %in% c("VA1","VC1","VT1","VrGE1","VE1","VE1b","VE1w"), "Parent Phenotype", "Offspring Phenotype"), levels=c("Parent Phenotype", "Offspring Phenotype")),
      parameter = factor(parameter, levels = c("VA1P","VA1","VA2","VC1","VT1","VF","VC2","VE1","VE2","VrGE1","VrGE2","VE1b","VE1w"))) %>%
    ggplot(aes(x = parameter,
               y = estimate, ymin =lbound, ymax=ubound,
               fill = parameter,
               label=paste0(round(estimate*100),"%\n(",round(lbound*100),"%???",round(ubound*100),"%)")))  +
    geom_col(color="black") +
    geom_text(y = 1, size = 3, hjust=.5,vjust=1.1,
              fontface="bold",
              color="black") +
    geom_point(size=4) +
    geom_errorbar( width=0, linewidth =1) +

    scale_fill_viridis_d(end = .85, guide="none", option = "C") +
    scale_y_continuous(labels=function(x) paste0(round(x*100),"%")) +
    custom_theme() +
    geom_hline(yintercept = c(0,1), linewidth =1.5) +
    #guides(fill = guide_legend(override.aes = list(shape = NA)) ) +
    labs(title = "Variance Decomposition", subtitle = "Offspring Phenotype",
         y= "",
         x = "",
         fill = ""
    )




  p4 <- filter(results_vars,
               parameter %in% c("VA1_s","VC1_s","VT1_s","VrGE1_s", "VE1_s"),
               abs(estimate) > 0.001) %>%
    mutate(parameter = factor(parameter, levels = c("VA1P","VA1_s","VA2","VC1_s","VT1_s","VF","VC2","VE1_s","VE2","VrGE1_s","VrGE2","VE1b","VE1w"))) %>%
    ggplot(aes(x = parameter,
               y = estimate, ymin =lbound, ymax=ubound,
               fill = parameter,
               label=paste0(round(estimate*100),"%\n(",round(lbound*100),"%???",round(ubound*100),"%)")))  +
    geom_col(color="black") +
    geom_text(y = 1, size = 3, hjust=.5,vjust=1.1,
              fontface="bold",
              color="black") +
    geom_point(size=4) +
    geom_errorbar( width=0, linewidth =1) +

    scale_fill_viridis_d(end = .85, guide="none", option = "C") +
    scale_y_continuous(labels=function(x) paste0(round(x*100),"%")) +
    custom_theme() +
    geom_hline(yintercept = c(0,1), linewidth =1.5) +
    #guides(fill = guide_legend(override.aes = list(shape = NA)) ) +
    labs(title = "Variance Decomposition", subtitle = "Parental Sorting Factor",
         y= "",
         x = "",
         fill = ""
    )

  if (!is.null(sim_vars)) {
    p1 <- p1 + geom_point(aes(y = simulated), size = 3, shape=13, color = "red", stroke=1)
    p3 <- p3 + geom_point(aes(y = simulated), size = 3, shape=13, color = "red", stroke=1)
    p4 <- p4 + geom_point(aes(y = simulated), size = 3, shape=13, color = "red", stroke=1)
  }
  if (!is.null(sim_PO)) {
    p2 <- p2 + geom_point(aes(x = simulated), size = 3, shape=13, color = "red", stroke=1) +
      geom_point(data = filter(results_PO, parameter != "Parent Offspring Correlation", type == "raw"),
                 aes(x = simulated), size = 3, shape=13, color = "blue", stroke=1)


  }

  if (incl.sort) {
    return(grid.arrange(p1,p2,p3,p4, layout_matrix=rbind(c(4,2),c(1,3))))
  } else {
    return(grid.arrange(p1,p2,p3, layout_matrix=rbind(c(2,2),c(1,3))))
  }




}
