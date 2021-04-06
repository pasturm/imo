# run_imo ----------------------------------------------------------------
#' Runs ion mirror optimization.
#'
#' \code{run_imo} runs ion mirror optimization.
#'
#' The experiment needs to be configured in a imo_config file (TOML file format).
#' See \code{system.file("GLPMtuneR_config.toml", package = "imo")} for 
#' a template.
#' 
#' If you want to start from a previous best point, set \code{resume = TRUE}.
#' If there is a bestpoint_run.txt file in the imo_results directory, the starting 
#' values will be taken from there. Otherwise the latest bestpoint_run.txt from
#' all subfolders of imo_results will be taken.
#'
#' @param imo_config path and name of imo_config file (.toml)
#' @param type Ion mirror type. Currently supported are \code{"GLPM"}: gridless
#' planar mirrors, \code{"ZEIM"}: cylindrical Zhang-Enke three-element ion 
#' mirrors, \code{"PIM"}: planar ion mirrors.
#' @param resume if \code{FALSE} (default) starting values are taken from the
#' imo_config file, if \code{TRUE} starting values are taken from the last
#' bestpoint_run.txt (see 'Details').
#' @param write Write output files (\code{TRUE} (default) or \code{FALSE}).
#' @param plot Plot results (\code{TRUE} (default) or \code{FALSE}).
#' @param digits Controls the number of decimal places to print when
#' printing the best point values.
#' 
#' @return A data.frame containing the optimized best points from the last run. 
#'
#' @examples
#' \dontrun{
#' imo_config = system.file("ZEIMimo_config.toml", package = "imo")
#' bestpoint = run_imo(imo_config, type = "ZEIM")
#' }
#' 
#' @export
run_imo = function(imo_config, type = c("GLPM", "ZEIM", "PIM"),
                   resume = FALSE, write = TRUE, plot = TRUE, digits = 3) {
  
  # configuration --------------------------------------------------------------
  
  # load config file
  config = RcppTOML::parseTOML(imo_config)
  
  # number of processes
  np = config$np
  
  # Create np copies of R running in parallel
  cl = parallel::makeCluster(np)
  doParallel::registerDoParallel(cl)
  
  # Energy range for which time-of-flight variations are minimized
  range_E = config$range_E
  
  # result path
  imo_dir = file.path(dirname(imo_config), "imo_results")
  if (write) {
    if (!dir.exists(imo_dir)) { dir.create(imo_dir) }
  }
  
  # response variables for optimization
  responses = do.call(rbind.data.frame, c(config$responses, stringsAsFactors = FALSE))
  
  # delete previous results file if run was aborted
  if (file.exists(file.path(imo_dir, "results.txt"))) {
    file.remove(file.path(imo_dir, "results.txt"))
  }
  
  # factors
  factors = do.call(rbind.data.frame, c(config$factors, stringsAsFactors = FALSE))
  factors$Enabled = as.logical(factors$Enabled)
  nfact0 = length(factors$Name)  # number of factors
  nfact = length(factors$Name[factors$Enabled])  # number of enabled factors
  
  # controls with start values
  controls = do.call(rbind.data.frame, c(config$controls, stringsAsFactors = FALSE))
  ncont = length(controls$Name)  # number of factors
  
  # possibly overwrite start values with previous bestpoint_run.txt
  if (resume) {
    allfiles = list.files(imo_dir, pattern = "bestpoint_run.txt", full.names = TRUE, recursive = TRUE)
    if (length(allfiles)>0) {
      lastfile = allfiles[length(allfiles)]
      warning(paste("Start values are taken from", lastfile), immediate. = TRUE)
      bestpoint_csv = utils::read.csv(lastfile, stringsAsFactors = FALSE, sep = "|")
      for (i in 2:(length(bestpoint_csv)-1)) {
        controls$StartValue[controls$Name==names(bestpoint_csv)[i]] = as.numeric(bestpoint_csv[i])
      }
    } else if (exists("bestpoint")) {
      for (i in 2:(length(bestpoint)-1)) {
        controls$StartValue[controls$Name==names(bestpoint)[i]] = as.numeric(bestpoint[i])
      }
    }
  }
  
  # assign starting values of controls to a variable
  for (i in 1:ncont) {
    if (!is.na(controls$StartValue[i])) {
      assign(controls$Name[i], controls$StartValue[i])
    }
  }
  
  # evaluate factor values
  for (i in 1:nfact0) {
    factors$Value[i] = eval(parse(text = factors$Transformation[i]))
  }
  
  # formula names for rsm
  xnam = paste0("x", 1:nfact)
  
  # sequence loop --------------------------------------------------------------
  timestring = format(Sys.time(), "%Y-%m-%d-%Hh%Mm%Ss")
  for (k in 1:config$n_repeats) {
    
    # generate Box-Behnken design ------------------------------------------------
    
    if (k>1) {
      # load factors from bestpoint of last run
      for (i in 1:length(bestpoint)) {
        factors$Value[factors$Name==names(bestpoint)[i]] = as.numeric(bestpoint[i])
      }
    }
    
    factors_enabled = factors[factors$Enabled,]
    
    # coerce to limits
    for (i in (1:nfact)) {
      if (is.na(factors_enabled$LowLimit[i]) & !is.na(factors_enabled$HighLimit[i])) {
        factors_enabled$Value[i] = min(factors_enabled$Value[i],
                                       factors_enabled$HighLimit[i] - factors_enabled$Range[i])
      }
      else if (!is.na(factors_enabled$LowLimit[i]) & is.na(factors_enabled$HighLimit[i])) {
        factors_enabled$Value[i] = max(factors_enabled$LowLimit[i] + factors_enabled$Range[i],
                                       factors_enabled$Value[i])
      }
      else if (!is.na(factors_enabled$LowLimit[i]) & !is.na(factors_enabled$HighLimit[i])) {
        minValue = min(factors_enabled$LowLimit[i] + factors_enabled$Range[i],
                       factors_enabled$HighLimit[i])
        maxValue = max(factors_enabled$HighLimit[i] - factors_enabled$Range[i],
                       factors_enabled$LowLimit[i])
        factors_enabled$Value[i] = min(max(minValue, factors_enabled$Value[i]), maxValue)
      }
    }
    
    # bbd coding (see ?bbd)
    coding = lapply(1:nfact, function(i)
      stats::as.formula(paste(xnam[i], "~ (", factors_enabled$Name[i], 
                              if (factors_enabled$Value[i] > 0) {"-"} else {"+"},
                              abs(factors_enabled$Value[i]), ") /", factors_enabled$Range[i])))
    design = rsm::bbd(nfact, n0 = 1, randomize = FALSE, coding = coding, block = FALSE)
    
    # run optimization runs -----------------------------------------------------
    # make imo_results result directory
    resultdir = paste(timestring, formatC(k, width = 2, flag = "0"), sep="_")
    if (write) {
      dir.create(file.path(imo_dir, resultdir))
    }
    
    # bbd_data
    bbd_data = rsm::decode.data(design[,3:(3+nfact-1)])
    
    # assign factor values to variables
    for (i in 1:nfact0) {
      assign(factors$Name[i], factors$Value[i])
    }
    
    # overwrite enabled factors with bbd_data values
    for (i in 1:nfact) {
      assign(names(bbd_data)[i], bbd_data[i])
    }
    
    # make runs data.frame
    runs = data.frame(run_no = seq_len(length(bbd_data[,1])))
    for (i in 1:ncont) {
      runs[controls$Name[i]] = tryCatch(eval(parse(text = controls$Transformation[i])),
                                        error = function(e) controls$StartValue[i])
    }
    
    if (write) {
      utils::write.table(signif(runs, 12), file = file.path(imo_dir, resultdir, "runs.txt"),
                         sep = "|", row.names = FALSE, col.names = TRUE, eol = "|\n")
    }
    
    print(paste0("Repeat ", k, ", run ", 1, " to ", length(runs[,1]),  " running..."))
    
    # run experiments
    E = seq(1-range_E/200, 1+range_E/200, 0.005)  # energies (keV)

    type = match.arg(type)
    if (type == "GLPM") {
      # `%dopar%` <- foreach::`%dopar%`
      result = foreach::foreach(i = 1:length(runs$run_no), .combine = "rbind") %dopar% {
        L = c(runs$L1[i], runs$L2[i], runs$L3[i], runs$L4[i], runs$L5[i], runs$L6[i])
        V = c(runs$V1[i], runs$V2[i], runs$V3[i], runs$V4[i], runs$V5[i], runs$V6[i])
        x1 = glpm_find_x1(L, V)
        tmp = glpm_tofperiod(E = E, x1 = x1, L, V)
        return(c(i, 1/stats::sd(tmp), x1, tmp))
      }
      # for (i in 1:length(runs$run_no)) {
      #   L = c(runs$L1[i], runs$L2[i], runs$L3[i], runs$L4[i], runs$L5[i], runs$L6[i])
      #   V = c(runs$V1[i], runs$V2[i], runs$V3[i], runs$V4[i], runs$V5[i], runs$V6[i])
      #   x1 = glpm_find_x1(L, V)
      #   tmp = glpm_tofperiod(E = E, x1 = x1, L, V)
      #   if (i==1) {
      #     result = c(i, 1/stats::sd(tmp), x1, tmp)
      #   } else {
      #     result = rbind(result, c(i, 1/stats::sd(tmp), x1, tmp))
      #   }
      # }
    } else if (type == "ZEIM") {
      result = foreach::foreach(i = 1:length(runs$run_no), .combine = "rbind") %dopar% {
        x1 = zeim_find_x1(runs$Z1[i], runs$Z2[i], runs$L[i], runs$V1[i], runs$V2[i], runs$R[i])
        tmp = zeim_tofperiod(E = E, runs$Z1[i], runs$Z2[i], runs$L[i], runs$V1[i], 
                             runs$V2[i], runs$R[i]) + x1*1/sqrt(E)
        return(c(i, 1/stats::sd(tmp), x1, tmp))
        }
    } else if (type == "PIM") {
      result = foreach::foreach(i = 1:length(runs$run_no), .combine = "rbind") %dopar% {
        z = c(runs$Z1[i], runs$Z2[i], runs$Z3[i], runs$Z4[i], runs$L[i])
        V = c(0, runs$V1[i], runs$V2[i], runs$V3[i], runs$V4[i])
        x1 = pim_find_x1(z, V, runs$D5[i], runs$U5[i])
        tmp = pim_totaltof(E = E, z, V, x1, runs$D5[i], runs$U5[i])
        return(c(i, 1/stats::sd(tmp), x1, tmp))
      }
    }
    
    parallel::stopCluster(cl)
    result = as.data.frame(result, row.names = NA)
    names(result) = c("no", "res", "x1")
    
    # plot
    if (plot) {
      tmp = result[,!(names(result) %in% c("no", "res", "x1"))]
      ymin = min((tmp-tmp[,ceiling(length(E)/2)])/tmp[,ceiling(length(E)/2)]*1e6)
      ymax = max((tmp-tmp[,ceiling(length(E)/2)])/tmp[,ceiling(length(E)/2)]*1e6)
      graphics::plot((E-1)*100, (tmp[1,]-tmp[1,ceiling(length(E)/2)])/tmp[1,ceiling(length(E)/2)]*1e6,
                     type = "n", main = paste("Repeat", k), ylim = c(ymin, ymax),
                     ylab = "",
                     xlab = expression(paste(Delta,"E") / "E" %.% 100))
      graphics::title(ylab = expression(paste(Delta,"t") / "t" %.% 10^{6}), line = 2.5)
      graphics::grid()
      for (i in 1:length(tmp[,1])) {
        graphics::lines((E-1)*100, (tmp[i,]-tmp[i,ceiling(length(E)/2)])/tmp[i,ceiling(length(E)/2)]*1e6, 
                        col = grDevices::rgb(0,0,0,0.3))
      }
    }
    
    result = result[c("no", "res", "x1")]
    
    if (write) {
      utils::write.table(signif(result, 12), file = file.path(imo_dir, resultdir, "results.txt"), 
                         sep = "|", row.names = FALSE, col.names = TRUE, eol = "|\n")
    }
    
    # fit quadratic model and optimize -------------------------------------------
    
    result = result[order(result$no),]  # order
    
    design$res = result$res
 
    # fit second order model
    tuner_rsm_R = rsm::rsm(stats::as.formula(paste("res ~ + SO(", paste(xnam, collapse = ","), ")")), 
                           data = design)
    
    # optimize (using L-BFGS-B method)
    
    # Pred_min/Pred_max values (used to define desirability function)
    Pred_min_R = stats::optim(rep(0, nfact), fn = cost_function, rsm_output = tuner_rsm_R,
                              nfact = nfact, xnam = xnam, method = "L-BFGS-B",
                              lower = rep(-1, nfact), upper = rep(1, nfact),
                              control = list(fnscale = 1))$value
    Pred_max_R = stats::optim(rep(0, nfact), fn = cost_function, rsm_output = tuner_rsm_R,
                              nfact = nfact, xnam = xnam, method = "L-BFGS-B",
                              lower = rep(-1, nfact), upper = rep(1, nfact),
                              control = list(fnscale = -1))$value
    
    # run optimizer
    factors_optim = stats::optim(rep(0, nfact), fn = desirability_overall,
                                 Target = responses$Target, w = responses$Weight,
                                 nfact = nfact, xnam = xnam,
                                 tuner_rsm = list(tuner_rsm_R),
                                 Pred_min = Pred_min_R,
                                 Pred_max = Pred_max_R, method = "L-BFGS-B",
                                 lower = rep(-1, nfact), upper = rep(1, nfact))$par
    
    factors_optim_coded = data.frame(t(factors_optim))
    colnames(factors_optim_coded) = xnam[1:nfact]
    bestpoint = rsm::code2val(factors_optim_coded, rsm::codings(design))
    bestpoint_predicted = data.frame(res = stats::predict(tuner_rsm_R, factors_optim_coded))
    
    # verify bestpoint -----------------------------------------------------------
    # assign bestpoint values to variables
    for (i in 1:nfact) {
      assign(names(bestpoint)[i], bestpoint[i])
    }
    
    # make runs data.frame
    bestpoint_run = data.frame(run_no = 1)
    for (i in 1:ncont) {
      bestpoint_run[controls$Name[i]] = tryCatch(eval(parse(text = controls$Transformation[i])),
                                                 error = function(e) controls$StartValue[i])
    }
    
    if (write) {
      utils::write.table(signif(bestpoint_run, 12), 
                         file = file.path(imo_dir, resultdir, "bestpoint_run.txt"), 
                         sep = "|", row.names = FALSE, col.names = TRUE, eol = "|\n")
    }
    
    # run experiment
    bestpoint_result = data.frame(no = 1, res = NA, x1 = NA)
    
    # adjust to match choosen parameters
    if (type == "GLPM") {
      L = c(bestpoint_run$L1, bestpoint_run$L2, bestpoint_run$L3, bestpoint_run$L4, bestpoint_run$L5, bestpoint_run$L6)
      V = c(bestpoint_run$V1, bestpoint_run$V2, bestpoint_run$V3, bestpoint_run$V4, bestpoint_run$V5, bestpoint_run$V6)
      
      x1 = glpm_find_x1(L, V)
      tmp = glpm_tofperiod(E = E, x1 = x1, L, V)
      # bestpoint_result$res = 1/stats::sd(tmp)  # resolution based on response variable
      
      # calculate time-of-flight variation (resolution) for 10 % energy variation.
      E_res = seq(0.95, 1.05, 0.005)  # energies (keV)
      tmp_res = glpm_tofperiod(E = E_res, x1 = x1, L, V)
      resolution = 1/diff(range((2*(tmp_res-tmp_res[11]))/tmp_res[11]))
      bestpoint_result$res = resolution
    } else if (type == "ZEIM") {
      x1 = zeim_find_x1(bestpoint_run$Z1, bestpoint_run$Z2, bestpoint_run$L,
                        bestpoint_run$V1, bestpoint_run$V2, bestpoint_run$R)
      tmp = zeim_tofperiod(E = E, bestpoint_run$Z1, bestpoint_run$Z2,
                           bestpoint_run$L, bestpoint_run$V1, bestpoint_run$V2,
                           bestpoint_run$R) + x1*1/sqrt(E)
      bestpoint_result$res = 1/stats::sd(tmp)
    } else if (type == "PIM") {
      z = c(bestpoint_run$Z1, bestpoint_run$Z2, bestpoint_run$Z3, bestpoint_run$Z4, bestpoint_run$L)
      V = c(0, bestpoint_run$V1, bestpoint_run$V2, bestpoint_run$V3, bestpoint_run$V4)
      x1 = pim_find_x1(z, V, bestpoint_run$D5, bestpoint_run$U5)
      tmp = pim_totaltof(E = E, z, V, x1, bestpoint_run$D5, bestpoint_run$U5)
      bestpoint_result$res = 1/stats::sd(tmp)
    }
    
    bestpoint_result$x1 = x1
    
    result_string = bestpoint_run[2:length(bestpoint_run)]
    print(paste0("Best point: ", paste0(names(result_string), "=", 
                                        round(result_string, digits), collapse = "|"),
                 "|x1=", round(x1,2),"|res=",floor(bestpoint_result$res)))
    # plot
    if (plot) {
      graphics::lines((E-1)*100, (tmp-tmp[ceiling(length(E)/2)])/tmp[ceiling(length(E)/2)]*1e6, 
                      col = grDevices::rgb(1,0,0,1))
      graphics::mtext(paste0(names(result_string), "=", 
                   round(result_string,digits), collapse = " | "), side = 1, line = -1, cex = 0.8)
    }
    
    if (write) {
      utils::write.table(signif(bestpoint_result, 12), 
                         file = file.path(imo_dir, resultdir, "bestpoint_results.txt"), 
                         sep = "|", row.names = FALSE, col.names = TRUE, eol = "|\n")
    
      # save all data as .RData
      save(list = ls(all.names = TRUE), 
           file = file.path(imo_dir, resultdir, "results.RData"))
    }
    
    # end of loop
  }
  bestpoint_run$x1 = x1
  return(bestpoint_run)
  
}
