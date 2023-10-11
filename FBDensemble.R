#!/usr/local/bin/Rscript
#
# This R script is an edited version of code provided by Louca et al. for their paper "Fundamental identifiability limits in molecular epidemiology" (2021).
# Edits (by Kate Truman) were made in order to conduct a similar simulation for FBD trees.
# Edits are as follows:
# Only generate one tree type of larger size (175,000 - 200,000 tips)
# Include extant and extinct tips in generated trees
# Remove multiple grid sizes and selection via AIC due to technical issues
# Focus on piecewise linear models only (not also skyline models)
# Specify set psi as in Louca et al.'s' paper
# Allow labelling of sampled ancestors in tree output.
#
# Original comments from Louca et al.:
#     This R script is provided as a Supplemental code to the paper:
#     Louca, S., McLaughlin, A., MacPherson, A., Joy, J.B., Pennell, M.W. (in review as of 2021). Fundamental identifiability limits in molecular epidemiology.
#     If you want to run a smaller number of simulations, modify the parameter ENSEMBLE_HBD_FITTING_NSIMS below. 
#     You can also reduce the number of fitting trials per tree, at the cost of fitting accuracy, through the parameter FITTING_NTRIALS.
#     If you have a machine with many cores, you can utilize those by adjusting the parameter NUMBER_OF_PARALLEL_THREADS.
#
#     LICENSE AGREEMENT
#     - - - - - - - - -
#     THIS CODE IS PROVIDED BY THE AUTHOR (STILIANOS LOUCA) "AS IS" AND ANY 
#     EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
#     OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
#     IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, 
#     INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
#     PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
#     HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
#     (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS CODE, 
#     EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#     - - - - - - - - -
#
# Stilianos Louca
# March 25, 2021

###################################
# OPTIONS

REQUIRED_PACKAGES	= c("Rcpp", "nloptr", "ape", "castor") # list any required packages here 
OUTPUT_DIR			= "output"


NUMBER_OF_PARALLEL_THREADS 	= 50 	# number of parallel threads to use for fitting
MODEL_ADEQUACY_NBOOTSTRAPS 	= 1000 	# number of bootstraps to use for evaluating the adequacy of fitted models. 
MODEL_ADEQUACY_MAX_RUNTIME	= 10 	# maximum runtime (seconds) per bootstrap simulation when testing model adequacy

# fitting models to trees
FITTING_NTRIALS					= 100 # number of fitting trials per tree. A larger number leads to more accurate fits (by avoiding non-global likelihood maxima), but takes longer to compute.
FITTING_NSTART_ATTEMPTS			= 100
FITTING_NITERATIONS 			= 500
FITTING_NEVALUATIONS 			= 1000
FITTING_REL_TOLERANCE 			= 1e-12
FITTING_CONDITIONING			= "auto" 	# either "crown" or "stem" or "auto"
FITTING_STEP_MIN				= 0.001
FITTING_HOMOGENOUS_GRID			= FALSE
fitting_Ntips2max_model_runtime = function(Ntips) max(2,Ntips/1e4) # runtime in seconds to allocate for likelihood evaluations during fitting, as a function of tree size

ENSEMBLE_HBD_FITTING_NSIMS 			 	 		 = 15 # number of trees to simulate and fit models to, in each of the categories "exp" and "OU"
ENSEMBLE_HBD_FITTING_MIN_NTIPS		 	 		 = 1000 
ENSEMBLE_HBD_FITTING_MAX_NTIPS			 		 = 5000
ENSEMBLE_HBD_FITTING_REPEAT_FAILED_TREES 		 = TRUE
ENSEMBLE_HBD_FITTING_SKYLINE_FIX_PRESENT_DAY_PSI = TRUE
ENSEMBLE_HBD_FITTING_PLINEAR_FIX_PRESENT_DAY_PSI = TRUE
INCLUDE_EXS = TRUE

# plot styles
BLACK_CURVE_COLOR		= "#303030"
BLUE_CURVE_COLOR		= "#005a96"
GREY_CURVE_COLOR		= "#707070"
CI50_CURVE_COLOR		= "#707070"
CI50_SHADE_COLOR		= "#909090"
CI95_CURVE_COLOR		= "#BBBBBB"
CI95_SHADE_COLOR		= "#DDDDDD"
REFERENCE_CURVE_COLOR	= "#606060"

PLOT_COLOR_PALETTE		= c(BLACK_CURVE_COLOR, BLUE_CURVE_COLOR, "red", "brown", "darkgreen", "#D3762A", "#618ACC", "#B80DC4", "#565656")
PLOT_LINE_TYPE_PALETTE	= rep(c(1,1,2,2,3,3),ages=5)
PLOT_LINE_WIDTH_PALETTE	= rep(c(2,1),ages=5)

DEFAULT_PLOT_WIDTH	= 2.2 	# inches
DEFAULT_PLOT_HEIGHT	= 1.8  	# inches
DEFAULT_PLOT_TOP_MARGIN	= 1	# inches
DEFAULT_PLOT_RIGHT_MARGIN_WITH_LEGEND = 3 # inches
PLOT_DOWNSAMPLING_RESOLUTION = 1000


options(expressions=100000) # increase recursion depth
options(max.print=1000000) # increase print count
options(warn=1) # print warnings as they occur


###################################
# AUXILIARY FUNCTIONS

# KT: generate_tree_hbds function from castor edited to allow for past sample labelling

# generate a random phylogenetic tree according to a homogenous birth-death-sampling process
# the speciation, extinction and continuous (Poissonian) sampling rate can each be time-dependent, and there may be additional discrete sampling times included
# The simulation proceeds in forward time (starting from the root) until one of the stopping criteria are met, OR once all lineages are extinct.
generate_tree_hbds_man = function(max_sampled_tips		= NULL, 
                                  max_sampled_nodes		= NULL, 
                                  max_extant_tips			= NULL,
                                  max_extinct_tips		= NULL,
                                  max_tips				= NULL, 	# integer, max number of tips of any type (extant + extinct + sampled). The simulation is halted once this number is reached.
                                  max_time				= NULL,
                                  include_extant			= FALSE,	# logical, whether to include extant non-sampled tips in the final tree
                                  include_extinct			= FALSE,	# logical, whether to include extinct non-sampled tips in the final tree
                                  as_generations			= FALSE,	# if FALSE, then edge lengths correspond to time. If TRUE, then edge lengths correspond to generations (hence if include_extant==true and include_extinct==true, all edges will have unit length).
                                  time_grid				= NULL,		# numeric vector listing grid times in ascending order. The time grid should generally cover the maximum possible simulation time, otherwise everything is polynomially extrapolated (according to splines_degree).
                                  lambda					= NULL,		# numeric vector of the same length as time_grid[], listing per-capita birth rates (speciation rates) at each time_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
                                  mu						= NULL,		# numeric vector of the same length as time_grid[], listing per-capita death rates (extinction rates) at each time_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
                                  psi						= NULL,		# numeric vector of the same length as time_grid[], listing per-capita sampling rates (Poissonian detection rates) at each time_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
                                  kappa					= NULL,		# numeric vector of the same length as time_grid[], listing the retention probability (upon sampling) at each time_grid point, i.e. the probability that a sampled lineage remains in the pool. If 0, then every sampled lineage becomes a tip.
                                  splines_degree			= 1,		# polynomial degree of time-dependent model parameters (lambda, mu, psi, kappa) between time-grid points
                                  CSA_times				= NULL,		# optional numeric vector listing concentrated sampling times, in ascending order
                                  CSA_probs				= NULL,		# optional numeric vector listing concentrated sampling probabilities, corresponding to CSA_times[]
                                  CSA_kappas				= NULL,		# optional numeric vector listing retention probabilities during concentrated sampling attempts, corresponding to CSA_times[]
                                  no_full_extinction		= FALSE,	# if true, then extinction of the entire tree is prevented. This is done by temporarily disabling extinctions when the number of extant tips is 1.
                                  max_runtime				= NULL,		# maximum time (in seconds) to allow for the computation; if the computation roughly exceeds this threshold, it is aborted. Use this as protection against badly parameterized models. If NULL or <=0, this option is ignored.
                                  tip_basename			= "",		# basename for tips (e.g. "tip."). 
                                  node_basename			= NULL,		# basename for nodes (e.g. "node."). If NULL, then nodes will not have any labels.
                                  edge_basename			= NULL,		# basename for edge (e.g. "edge."). If NULL, then edges will not have any labels.
                                  include_birth_times		= TRUE,
                                  include_death_times		= TRUE){
  # basic input checking
  if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))
  if(is.null(max_tips) && is.null(max_sampled_tips) && is.null(max_extant_tips) && is.null(max_extinct_tips) && is.null(max_sampled_nodes) && is.null(max_time)) return(list(success=FALSE, error="ERROR: At least one of {max_tips, max_sampled_tips, max_sampled_nodes, max_extant_tips, max_extinct_tips, max_time} must be non-NULL"));
  if(is.null(time_grid) || (length(time_grid)<=1)){
    if((!is.null(lambda)) && (length(lambda)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of lambda values (%d); since no time grid was provided, you must either provide a single (constant) birth rate or NULL",length(lambda))))
    if((!is.null(mu)) && (length(mu)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of mu values (%d); since no time grid was provided, you must either provide a single (constant) death rate or none",length(mu))))
    if((!is.null(psi)) && (length(psi)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of psi values (%d); since no time grid was provided, you must either provide a single (constant) sampling rate or none",length(psi))))
    if((!is.null(kappa)) && (length(kappa)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of kappa values (%d); since no time grid was provided, you must either provide a single (constant) retention probability or none",length(kappa))))
    # create dummy time grid
    NG = 2;
    time_grid = seq(from=0,to=1,length.out=NG)
    if(!is.null(lambda)){
      lambda = rep(lambda,times=NG)
    }else{
      lambda = rep(0,times=NG)
    }
    if(!is.null(mu)){
      mu = rep(mu,times=NG)
    }else{
      mu = rep(0,times=NG)
    }
    if(!is.null(psi)){
      psi = rep(psi,times=NG)
    }else{
      psi = rep(0,times=NG)
    }
    if(!is.null(kappa)){
      kappa = rep(kappa,times=NG)
    }else{
      kappa = rep(0,times=NG)
    }
  }else{
    NG = length(time_grid);
    if(is.null(lambda)){
      lambda = rep(0,times=NG)
    }else if(length(lambda)==1){
      lambda = rep(lambda,times=NG)
    }else if(length(lambda)!=NG){
      return(list(success=FALSE, error=sprintf("Expected either a single birth-rate lambda or exactly %d birth-rates (=time_grid length), but instead got %d",NG,length(lambda))))
    }
    if(is.null(mu)){
      mu = rep(0,times=NG)
    }else if(length(mu)==1){
      mu = rep(mu,times=NG)
    }else if(length(mu)!=NG){
      return(list(success=FALSE, error=sprintf("Expected either a single death-rate mu or exactly %d death-rates (=time_grid length), but instead got %d",NG,length(mu))))
    }
    if(is.null(psi)){
      psi = rep(0,times=NG)
    }else if(length(psi)==1){
      psi = rep(psi,times=NG)
    }else if(length(psi)!=NG){
      return(list(success=FALSE, error=sprintf("Expected either a single sampling-rate psi or exactly %d sampling-rates (=time_grid length), but instead got %d",NG,length(psi))))
    }
    if(is.null(kappa)){
      kappa = rep(0,times=NG)
    }else if(length(kappa)==1){
      kappa = rep(kappa,times=NG)
    }else if(length(kappa)!=NG){
      return(list(success=FALSE, error=sprintf("Expected either a single retention probability kappa or exactly %d probabilities (=time_grid length), but instead got %d",NG,length(kappa))))
    }
    if(any(diff(time_grid)<=0)) return(list(success = FALSE, error = sprintf("Values in time_grid must be strictly increasing")))
  }
  NCSA = (if(is.null(CSA_times)) 0 else length(CSA_times))
  if((NCSA==0) && (!is.null(CSA_probs)) && (length(CSA_probs)>0)) return(list(success=FALSE, error="CSA_times is missing while CSA_probs was provided; either provide both or none"))
  if((NCSA==0) && (!is.null(CSA_kappas)) && (length(CSA_probs)>0)) return(list(success=FALSE, error="CSA_times is missing while CSA_kappas was provided; either provide both or none"))
  if((NCSA>0) && is.null(CSA_probs)) return(list(success=FALSE, error="CSA_probs is missing while CSA_times was provided; either provide both or none"))
  if((NCSA>0) && is.null(CSA_kappas)) return(list(success=FALSE, error="CSA_kappas is missing while CSA_times was provided; either provide both or none"))
  if((NCSA>0) && (!is.null(CSA_probs)) && (!is.null(CSA_kappas))){
    if(length(CSA_probs)==1) CSA_probs 	 = rep(CSA_probs, times=NCSA)
    if(length(CSA_kappas)==1) CSA_kappas = rep(CSA_kappas, times=NCSA)
    if(length(CSA_times)!=length(CSA_probs)) return(list(success=FALSE, error="Number of CSA_times (%d) differs from number of CSA_probs (%d)",length(CSA_times),length(CSA_probs)))
    if(length(CSA_times)!=length(CSA_kappas)) return(list(success=FALSE, error="Number of CSA_times (%d) differs from number of CSA_kappas (%d)",length(CSA_times),length(CSA_kappas)))
    if(any(diff(CSA_times)<=0)) return(list(success=FALSE, error="CSA_times must be in strictly increasing order"))
    if(any(CSA_probs<0) || any(CSA_probs>1)) return(list(success=FALSE, error="CSA_probs must be true probabilities, and thus between 0 and 1"))
    if(any(CSA_kappas<0) || any(CSA_kappas>1)) return(list(success=FALSE, error="CSA_kappas must be true probabilities, and thus between 0 and 1"))
  }
  if(is.null(max_runtime)) max_runtime = 0
  
  # generate tree
  results = generate_random_tree_HBDS_CPP(max_sampled_tips		= (if(is.null(max_sampled_tips)) -1 else max_sampled_tips),
                                                      max_sampled_nodes		= (if(is.null(max_sampled_nodes)) -1 else max_sampled_nodes),
                                                      max_extant_tips			= (if(is.null(max_extant_tips)) -1 else max_extant_tips),
                                                      max_extinct_tips		= (if(is.null(max_extinct_tips)) -1 else max_extinct_tips),
                                                      max_tips				= (if(is.null(max_tips)) -1 else max_tips),
                                                      max_time				= (if(is.null(max_time)) -1 else max_time),
                                                      time_grid				= time_grid,
                                                      birth_rates				= lambda,
                                                      death_rates				= mu,
                                                      sampling_rates			= psi,
                                                      retention_probs			= kappa,
                                                      splines_degree			= splines_degree,
                                                      CSA_times				= (if(is.null(CSA_times)) numeric() else CSA_times),
                                                      CSA_probs				= (if(is.null(CSA_probs)) numeric() else CSA_probs),
                                                      CSA_kappas				= (if(is.null(CSA_kappas)) numeric() else CSA_kappas),
                                                      as_generations			= as_generations,
                                                      no_full_extinction		= no_full_extinction,
                                                      runtime_out_seconds		= max_runtime,
                                                      include_extant			= include_extant,
                                                      include_extinct			= include_extinct,
                                                      include_birth_times		= include_birth_times,
                                                      include_death_times		= include_death_times)
  if(!results$success) return(list(success=FALSE, error=results$error)); # something went wrong
  Ntips	= results$Ntips
  Nnodes 	= results$Nnodes
  Nedges 	= results$Nedges
  # allow labelling of past sampling events
  allTips = 1:Ntips
  tipLabels = ifelse(allTips %in% results$sampled_clades & !(allTips %in% results$retained_clades), paste("removal.", allTips, sep=""), paste(tip_basename, allTips, sep=""))
  
  tree = list(Nnode 		= Nnodes,
              tip.label 	= tipLabels,
              node.label 	= (if(is.null(node_basename)) NULL else paste(node_basename, 1:Nnodes, sep="")),
              edge.label 	= (if(is.null(edge_basename)) NULL else paste(edge_basename, 1:Nedges, sep="")),
              edge 		= matrix(results$tree_edge,ncol=2,byrow=TRUE) + 1,
              edge.length = results$edge_length,
              root 		= results$root+1)
  class(tree) = "phylo";
  attr(tree,"order") = NULL
  
  return(list(success				= TRUE,
              tree				= tree,
              root_time			= results$root_time,
              final_time			= results$final_time,
              root_age			= results$final_time - results$root_time,
              Nbirths		 		= results$Nbirths,
              Ndeaths				= results$Ndeaths,
              Nsamplings			= results$Nsamplings,
              Nretentions			= results$Nretentions,
              sampled_clades		= results$sampled_clades+1,
              retained_clades		= results$retained_clades+1,
              extant_tips			= (if(include_extant) results$extant_tips+1 else integer()),
              extinct_tips		= (if(include_extinct) results$extinct_tips+1 else integer())));
  
}

# KT: move repeated calls to function
# Save parameters for congruent case 1
set_true_results = function(results_df, sim_true){
  str(sim_true)
  	results_df$true_slope_lambda[sim]				= get_linear_slope(x=sim_true$ages, y=sim_true$lambda, include_intercept=TRUE)
		results_df$true_slope_mu[sim]					= get_linear_slope(x=sim_true$ages, y=sim_true$mu, include_intercept=TRUE)
		results_df$true_slope_psi[sim]					= get_linear_slope(x=sim_true$ages, y=sim_true$psi, include_intercept=TRUE)
		results_df$true_slope_Reff[sim]				= get_linear_slope(x=sim_true$ages, y=sim_true$Reff, include_intercept=TRUE)
		results_df$true_slope_removal_rate[sim]		= get_linear_slope(x=sim_true$ages, y=sim_true$removal_rate, include_intercept=TRUE)
		results_df$true_slope_sampling_proportion[sim]	= get_linear_slope(x=sim_true$ages, y=sim_true$sampling_proportion, include_intercept=TRUE)
		results_df$true_slope_net_growth_rate[sim]		= get_linear_slope(x=sim_true$ages, y=sim_true$diversification_rate, include_intercept=TRUE)
		#results_df$true_slope_branching_density[sim]	= get_linear_slope(x=sim_true$ages, y=sim_true$branching_density, include_intercept=TRUE)
		#results_df$true_slope_sampling_density[sim]	= get_linear_slope(x=sim_true$ages, y=sim_true$sampling_density, include_intercept=TRUE)
		results_df$true_mean_lambda[sim]				= mean(sim_true$lambda, na.rm=TRUE)
		results_df$true_mean_mu[sim]					= mean(sim_true$mu, na.rm=TRUE)
		results_df$true_mean_psi[sim]					= mean(sim_true$psi, na.rm=TRUE)
		results_df$true_mean_Reff[sim]					= mean(sim_true$Reff, na.rm=TRUE)
		results_df$true_mean_removal_rate[sim]			= mean(sim_true$removal_rate, na.rm=TRUE)
		results_df$true_mean_sampling_proportion[sim]	= mean(sim_true$sampling_proportion, na.rm=TRUE)
		results_df$true_mean_net_growth_rate[sim]		= mean(sim_true$diversification_rate, na.rm=TRUE)
		# results_df$true_mean_branching_density[sim]	= mean(sim_true$branching_density, na.rm=TRUE)
		#results_df$true_mean_sampling_density[sim]		= mean(sim_true$sampling_density, na.rm=TRUE)
    return(results_df)
}

# KT: move repeated calls to function
# Save parameters for congruent case 2
set_plinear_results = function(results_df, sim_true, sim_fit){
      results_df$plinear_lambda_R2[sim] 					= get_R2(xtrue=sim_true$ages, ytrue=sim_true$lambda, xfit=sim_fit$ages, yfit=sim_fit$lambda)
      results_df$plinear_mu_R2[sim] 						= get_R2(xtrue=sim_true$ages, ytrue=sim_true$mu, xfit=sim_fit$ages, yfit=sim_fit$mu)
      results_df$plinear_psi_R2[sim] 					= get_R2(xtrue=sim_true$ages, ytrue=sim_true$psi, xfit=sim_fit$ages, yfit=sim_fit$psi)
      results_df$plinear_Reff_R2[sim] 					= get_R2(xtrue=sim_true$ages, ytrue=sim_true$Reff, xfit=sim_fit$ages, yfit=sim_fit$Reff)
      results_df$plinear_removal_rate_R2[sim] 			= get_R2(xtrue=sim_true$ages, ytrue=sim_true$removal_rate, xfit=sim_fit$ages, yfit=sim_fit$removal_rate)
      results_df$plinear_sampling_proportion_R2[sim] 	= get_R2(xtrue=sim_true$ages, ytrue=sim_true$sampling_proportion, xfit=sim_fit$ages, yfit=sim_fit$sampling_proportion)
      results_df$plinear_net_growth_rate_R2[sim] 		= get_R2(xtrue=sim_true$ages, ytrue=sim_true$diversification_rate, xfit=sim_fit$ages, yfit=sim_fit$diversification_rate)
      results_df$plinear_nLTT_R2[sim]					= get_R2(xtrue=sim_true$ages, ytrue=sim_true$nLTT, xfit=sim_fit$ages, yfit=sim_fit$nLTT)
      results_df$plinear_lambda_MMNE[sim] 				= get_MMNE(xtrue=sim_true$ages, ytrue=sim_true$lambda, xfit=sim_fit$ages, yfit=sim_fit$lambda)
      results_df$plinear_mu_MMNE[sim] 					= get_MMNE(xtrue=sim_true$ages, ytrue=sim_true$mu, xfit=sim_fit$ages, yfit=sim_fit$mu)
      results_df$plinear_psi_MMNE[sim] 					= get_MMNE(xtrue=sim_true$ages, ytrue=sim_true$psi, xfit=sim_fit$ages, yfit=sim_fit$psi)
      results_df$plinear_Reff_MMNE[sim] 					= get_MMNE(xtrue=sim_true$ages, ytrue=sim_true$Reff, xfit=sim_fit$ages, yfit=sim_fit$Reff)
      results_df$plinear_removal_rate_MMNE[sim] 			= get_MMNE(xtrue=sim_true$ages, ytrue=sim_true$removal_rate, xfit=sim_fit$ages, yfit=sim_fit$removal_rate)
      results_df$plinear_sampling_proportion_MMNE[sim] 	= get_MMNE(xtrue=sim_true$ages, ytrue=sim_true$sampling_proportion, xfit=sim_fit$ages, yfit=sim_fit$sampling_proportion)
      results_df$plinear_net_growth_rate_MMNE[sim] 		= get_MMNE(xtrue=sim_true$ages, ytrue=sim_true$diversification_rate, xfit=sim_fit$ages, yfit=sim_fit$diversification_rate)
      results_df$plinear_nLTT_MMNE[sim]					= get_MMNE(xtrue=sim_true$ages, ytrue=sim_true$nLTT, xfit=sim_fit$ages, yfit=sim_fit$nLTT)
      results_df$plinear_slope_lambda[sim]				= get_linear_slope(x=sim_fit$ages, y=sim_fit$lambda, include_intercept=TRUE)
      results_df$plinear_slope_mu[sim]					= get_linear_slope(x=sim_fit$ages, y=sim_fit$mu, include_intercept=TRUE)
      results_df$plinear_slope_psi[sim]					= get_linear_slope(x=sim_fit$ages, y=sim_fit$psi, include_intercept=TRUE)
      results_df$plinear_slope_Reff[sim]					= get_linear_slope(x=sim_fit$ages, y=sim_fit$Reff, include_intercept=TRUE)
      results_df$plinear_slope_removal_rate[sim]			= get_linear_slope(x=sim_fit$ages, y=sim_fit$removal_rate, include_intercept=TRUE)
      results_df$plinear_slope_sampling_proportion[sim]	= get_linear_slope(x=sim_fit$ages, y=sim_fit$sampling_proportion, include_intercept=TRUE)
      results_df$plinear_slope_net_growth_rate[sim]		= get_linear_slope(x=sim_fit$ages, y=sim_fit$diversification_rate, include_intercept=TRUE)
      results_df$plinear_mean_lambda[sim]				= mean(sim_fit$lambda, na.rm=TRUE)
      results_df$plinear_mean_mu[sim]					= mean(sim_fit$mu, na.rm=TRUE)
      results_df$plinear_mean_psi[sim]					= mean(sim_fit$psi, na.rm=TRUE)
      results_df$plinear_mean_Reff[sim]					= mean(sim_fit$Reff, na.rm=TRUE)
      results_df$plinear_mean_removal_rate[sim]			= mean(sim_fit$removal_rate, na.rm=TRUE)
      results_df$plinear_mean_sampling_proportion[sim]	= mean(sim_fit$sampling_proportion, na.rm=TRUE)
      results_df$plinear_mean_net_growth_rate[sim]		= mean(sim_fit$diversification_rate, na.rm=TRUE)
      return(results_df)
}

# KT: Simulate deterministic values given specified evolutionary rates.
sim_plinear = function(fit,sim_true,age0,tree_LTT0,kappa, results_df, fit_dir){
  # Calculate deterministic properties of fitted model.
   sim_fit = simulate_deterministic_hbds(	age_grid		= fit$age_grid, 
                                                   lambda			= fit$param_fitted$lambda,
                                                   mu				= fit$param_fitted$mu,
                                                   psi				= fit$param_fitted$psi,
                                                   kappa			= kappa,
                                                   splines_degree	= 1,
                                                   requested_ages	= sim_true$ages,
                                                   age0			= age0,
                                                   LTT0			= tree_LTT0)
    if(!sim_fit$success){
      cat2(sprintf("      WARNING: Simulation failed: %s\n",sim_fit$error))
      return(list(FALSE))
    }else{
      # Save results
      results_df = set_plinear_results(results_df, sim_true, sim_fit)
      cat2(sprintf("    Plotting fitted HBDS model..\n"))
      # Comparison plot
      plot_fitted_vs_true_model(	plot_dir			= fit_dir,
                                 case_tag			= "plinear comparison",
                                 subtitle			= NULL,
                                 true_model_name	 	= sprintf("%s.sim_%d%s",scenario$name,sim,formatC(kappa, digits = 1, format = "f")),
                                 fit_model_name		= "plinear",
                                 sim_true			= sim_true,
                                 sim_fit				= sim_fit,
                                 tree_LTT			= tree_LTT,
                                 root_age			= root_age,
                                 time_unit			= scenario$time_units,
                                 verbose				= TRUE,
                                 verbose_prefix		= "      ")
      return(list(TRUE,results_df, sim_fit))
    }
}

# KT: Fit and plot a piecewise linear model to obtain congruent case 2 given the existing tree from congruent case 1.
plinear_fit_and_plot = function(sim_true, tree, properties, correct_psi = FALSE, kappa, results_df){
  root_age = properties[[1]]
  stem_age = properties[[2]]
  end_age = properties[[3]]
  tree_LTT = properties[[4]]
  age0 = properties[[5]]
  tree_LTT0 = properties[[6]]

  grid_to_fit = seq(from = 0, to = properties[[1]], length.out = scenario$fitting_grid_size)
  # Fit parameters on one grid. As the new parameters correspond to the same tree, the second scenario will be congruent to the first.
  # Note that we do not specify a minimum lambda or mu value when fitting as we allow these to be greater or equal to zero.
  if (correct_psi == FALSE){
    fixed_psi=NULL
    if(ENSEMBLE_HBD_FITTING_SKYLINE_FIX_PRESENT_DAY_PSI){
      fixed_psi = rep(NA,scenario$fitting_grid_size)
      fixed_psi[grid_to_fit<=root_age/1000] = present_day_psi
    }
    str(tree)
    fit = fit_hbds_model_on_grid(tree				= tree,
                                          root_age 			= root_age,
                                          oldest_age			= root_age,
                                          age_grid = grid_to_fit,
                                          max_lambda			= 100*max(sim_true$lambda),
                                          max_mu				= 100*max(sim_true$mu),
                                          min_psi				= 0.01*min(sim_true$psi),
                                          max_psi				= 100*max(sim_true$psi),
                                          fixed_psi			= fixed_psi,
                                          fixed_kappa			= kappa,
                                          splines_degree		= 1,
                                          condition			= FITTING_CONDITIONING,
                                          Ntrials				= FITTING_NTRIALS,
                                          max_start_attempts	= FITTING_NSTART_ATTEMPTS,
                                          Nthreads			= NUMBER_OF_PARALLEL_THREADS,
                                          max_model_runtime	= fitting_Ntips2max_model_runtime(Ntips),
                                          fit_control			= list(eval.max=FITTING_NEVALUATIONS, iter.max=FITTING_NITERATIONS, rel.tol=FITTING_REL_TOLERANCE, step.min=FITTING_STEP_MIN),
                                          verbose				= TRUE,
                                          verbose_prefix		= "      ")
  }
  else {
    fit = fit_hbds_model_on_grid(tree				= tree, 
                                          root_age 			= root_age,
                                          oldest_age			= root_age,
                                          age_grid = grid_to_fit,
                                          max_lambda			= 100*max(sim_true$lambda),
                                          max_mu				= 100*max(sim_true$mu),
                                          fixed_psi			= sim_true$psi[c(1, seq(from=length(sim_true$psi)/(scenario$fitting_grid_size-1),to=length(sim_true$psi),length.out=scenario$fitting_grid_size-1))],
                                          fixed_kappa			= kappa,
                                          splines_degree		= 1,
                                          condition			= FITTING_CONDITIONING,
                                          Ntrials				= FITTING_NTRIALS,
                                          max_start_attempts	= FITTING_NSTART_ATTEMPTS,
                                          Nthreads			= NUMBER_OF_PARALLEL_THREADS,
                                          max_model_runtime	= fitting_Ntips2max_model_runtime(Ntips),
                                          fit_control			= list(eval.max=FITTING_NEVALUATIONS, iter.max=FITTING_NITERATIONS, rel.tol=FITTING_REL_TOLERANCE, step.min=FITTING_STEP_MIN),
                                          verbose				= TRUE,
                                          verbose_prefix		= "      ")
  }
  if(!fit$success){
    cat2(sprintf("      ERROR: Fitting failed: %s\n",fit$error));
    if(ENSEMBLE_HBD_FITTING_REPEAT_FAILED_TREES){
      cat2(sprintf("        Repeating entire simulation %d and fitting\n",sim))
      unlink(sim_dir, recursive=TRUE)
      return(list(FALSE, fit$error))
    }
    }else{
    
      fit_dir = sprintf("%s/fitted_plinear_psi_specified%s",sim_dir, correct_psi)
      # We no longer need to select the best fit
      # fit = fit$best_fit
      dir.create(fit_dir, showWarnings = FALSE, recursive=TRUE);
      sink(file=sprintf("%s/fit_results.txt",fit_dir)); print(fit); sink(); # save fit results to text file
      results_df$plinear_Ngrid[sim] 			= scenario$fitting_grid_size
      
      cat2(sprintf("    Simulating fitted HBDS plinear model..\n"));
      sim_result = sim_plinear(fit,sim_true,age0,tree_LTT0,kappa, results_df, fit_dir)
      if (!(sim_result[[1]])){
        return(list(FALSE))
      }
      results_df = sim_result[[2]]
      sim_fit = sim_result[[3]]
      cat2(sprintf("    Assessing adequacy of fitted HBDS plinear model..\n"))
      adequacy_age_grid = seq(from=end_age,to=stem_age,length.out=1000)
      adequacy = assess_model_adequacy(	tree 				= tree, 
                                        models				= list(list(ages=adequacy_age_grid, stem_age=stem_age, end_age=end_age, lambda=approx(x=sim_fit$ages,y=sim_fit$lambda,xout=adequacy_age_grid,rule=2)$y, mu=approx(x=sim_fit$ages,y=sim_fit$mu,xout=adequacy_age_grid,rule=2)$y, psi=approx(x=sim_fit$ages,y=sim_fit$psi,xout=adequacy_age_grid,rule=2)$y)),
                                        tree_name 			= sprintf("%s.sim_%d",scenario$name,sim),
                                        model_name			= "plinear_fit",
                                        Nbootstraps			= MODEL_ADEQUACY_NBOOTSTRAPS, 
                                        report_file			= sprintf("%s/fit_model_adequacy.txt",fit_dir), 
                                        max_extant_tips		= Ntips*10,
                                        Nthreads			= NUMBER_OF_PARALLEL_THREADS)
      if(!adequacy$success){
        cat2(sprintf("      WARNING: Adequacy test failed: %s\n",adequacy$error))
        return(list(FALSE))
      }else{
        cat2(sprintf("      --> PnodeKS = %g, PedgeKS = %g, PtipKS = %g\n",adequacy$PnodeKS,adequacy$PedgeKS,adequacy$PtipKS))
        results_df$plinear_PedgeKS[sim] = adequacy$PedgeKS
        results_df$plinear_PtipKS[sim]  = adequacy$PtipKS
        results_df$plinear_PnodeKS[sim] = adequacy$PnodeKS
      }
    }
  return(list(TRUE,results_df, fit, fit_dir))
}


check_output_file = function(file_path,force_replace,verbose,verbose_prefix){
  if(file.exists(file_path)){
    if(force_replace){
      cat(sprintf("%sNote: Replacing output file '%s'.\n",verbose_prefix,file_path))
      file.remove(file_path);
    }else{
      sprintf("Add to existing file %s.\n", file_path);
      return()
    }
  }
  dir.create(dirname(file_path), showWarnings = FALSE, recursive=TRUE)
}


prepare_output_file = function(file_path, force_replace, verbose, verbose_prefix){
  check_output_file(file_path,force_replace,verbose,verbose_prefix)
  fout = (if(endsWith(tolower(file_path),".gz")) gzfile(file_path, "wt") else file(file_path,"wt"))
  return(fout)
}


get_non_existent_dir = function(parent_path, child_basename, digits=3){
  existing = list.files(path=parent_path, full.names=FALSE, recursive=FALSE)
  counter = 1
  child_name = sprintf(sprintf("%%s%%0%dd",digits),child_basename,counter)
  while(any(sapply(existing, FUN=function(x) startsWith(x,child_name)))){
    counter = counter + 1;
    child_name = sprintf(sprintf("%%s%%0%dd",digits),child_basename,counter)
  }
  child_path = file.path(parent_path,child_name)
  return(child_path);
}


get_linear_slope = function(x,y,include_intercept=TRUE){
  if(include_intercept){
    fit = stats::lm(y~x, data=data.frame(x,y))
    slope = fit$coefficients[2]
  }else{
    fit = stats::lm(y~x-1, data=data.frame(x,y))	
    slope = fit$coefficients[1]
  }
  return(slope)
}


get_R2 = function(xtrue, ytrue, xfit, yfit, minx=NULL, maxx=NULL){
  if(is.null(minx)) minx = min(xtrue)
  if(is.null(maxx)) maxx = max(xtrue)
  yfit_on_xtrue 	= suppressWarnings(approx(x=xfit, y=yfit, xout=xtrue)$y)
  valids 			= which(is.finite(yfit_on_xtrue) & is.finite(ytrue) & (xtrue>=minx) & (xtrue<=maxx)) # only consider points where the fitted model and the data are defined
  if(length(valids)<=1) return(NaN);
  ytrue 			= ytrue[valids]
  xtrue 			= xtrue[valids]
  yfit_on_xtrue 	= yfit_on_xtrue[valids]
  SSR 			= sum((yfit_on_xtrue-ytrue)**2)
  return(1 - SSR/(length(ytrue)*var(ytrue)))
}



# calculate mean modulus of normalized error (normalized by the mean modulus of the true value)
get_MMNE = function(xtrue, ytrue, xfit, yfit, minx=NULL, maxx=NULL){
  if(is.null(minx)) minx = min(xtrue)
  if(is.null(maxx)) maxx = max(xtrue)
  yfit_on_xtrue 	= suppressWarnings(approx(x=xfit, y=yfit, xout=xtrue)$y)
  valids 			= which(is.finite(yfit_on_xtrue) & is.finite(ytrue) & (xtrue>=minx) & (xtrue<=maxx)) # only consider points where the fitted model and the data are defined
  if(length(valids)<=1) return(NaN);
  ytrue 			= ytrue[valids]
  xtrue 			= xtrue[valids]
  yfit_on_xtrue 	= yfit_on_xtrue[valids]
  return(mean(abs(yfit_on_xtrue-ytrue), na.rm=TRUE)/mean(abs(ytrue), na.rm=TRUE))
}


expand_range = function(xmin,xmax,factor,logarithmic=FALSE){
  if(logarithmic){
    return(exp(expand_range(log(xmin),log(xmax),factor,logarithmic=FALSE)));
  }else{
    L = xmax - xmin;
    M = 0.5*(xmin+xmax);
    return(c(M-factor*L/2, M+factor*L/2))
  }
}



# save & plot variables over age
plot_curves = function(	file_basepath, 		# e.g. 'output/SILVA_curves_last_1000years'
                        xtype,				# e.g. 'age'
                        ytype,				# e.g. 'birth_rates'
                        case_tag,			# e.g. 'model 2'
                        curves, 			# list of size Ncurves, each entry of which is a sub-list with two elements (ages, values) specifying a separate curve to be plotted
                        curve_names,		# 1D character vector of size >=Ncurves
                        curve_colors		= NULL,			# either NULL or a 1D vector of size >=Ncurves, specifying the color of a curve. If NULL, colors are picked automatically.
                        curve_line_types	= NULL,			# either NULL or a 1D vector of size >=Ncurves, specifying the line type of a curve. If NULL, line types are picked automatically.
                        curve_widths		= NULL,			# either NULL or a 1D vector of size >=Ncurves, specifying the line width of a curve. If NULL, line widths are picked automatically.
                        show_curves			= NULL,			# either a 1D vector of booleans of size Ncurves, indicating whether a curve should be plotted, or NULL (plot all curves)
                        shadings			= NULL,			# optional list of shading specifications. Each specification is a list of 3 entries (the 2 curves to shade between and the shading color). Can also be NULL. For example, list(list(1,2,'red'),list(5,6,'blue')) will shade the area between curves 1 & 2 in red and the area between curves 5 & 6 in blue.
                        reverse_x			= FALSE,		# boolean, specifying if the x-axis should be plotted reversed
                        plot_minx			= NULL,			# minimum x-value to plot. If NULL, this is automatically set by the data
                        plot_maxx			= NULL,			# maximum x-value to plot. If NULL, this is automatically set by the data
                        plot_miny			= NULL,			# optional lower bound for the plots. If NULL, this is automatically set by the data
                        plot_maxy			= NULL,			# optional upper bound for the plots. If NULL, this is automatically set by the data
                        resolution			= NULL,			# either NULL or an integer, specifying the temporal resolution for downsampling curves. May be needed in order to avoid excessively large data & plot files.
                        xlabel				= "x", 			# e.g. 'age (years)'
                        ylabel				= "y",			# e.g. 'number of lineages'
                        plot_log_values		= FALSE,
                        horizontal_reference= NULL,			# optional numeric, specifying a horizontal reference line to show (e.g., at value 0)
                        legend_pos			= "outside",	# (string) specification of legened position, e.g. "outside". If "none", no legend will be shown
                        plot_title			= "",
                        plot_width			= NULL,			# numeric, plot width in inches. If NULL, this is set to the document's default.
                        plot_height			= NULL,			# numeric, plot height in inches. If NULL, this is set to the document's default.
                        data_file_comments	= "",
                        save_data			= TRUE,			# (boolean)
                        scatterpoints		= NULL,			# optional 2D numeric matrix of size NP x 2 (points only, X & Y) or NP x 4 (points & vertical error bars, X & Y & minY & maxY) or NP x 6 (points & vertical & horizontal error bars, X & Y & minX & maxX & minY & maxY), specifying a list of scatterpoints to add to the figure
                        verbose 			= TRUE,
                        verbose_prefix 		= "  "){
  curves  = curves[sapply(curves,FUN = function(l) !is.null(l))] # remove NULL elements
  Ncurves = length(curves);
  if(is.null(curve_colors)){ curve_colors = PLOT_COLOR_PALETTE[1 + (seq_len(Ncurves)-1)%%length(PLOT_COLOR_PALETTE)]; }
  else if(length(curve_colors)==1){ curve_colors = rep(curve_colors,times=Ncurves); }
  if(is.null(curve_line_types)){ curve_line_types = PLOT_LINE_TYPE_PALETTE[1 + (seq_len(Ncurves)-1)%%length(PLOT_LINE_TYPE_PALETTE)] }
  else if(length(curve_line_types)==1){ curve_line_types = rep(curve_line_types,times=Ncurves); }
  if(is.null(curve_widths)){ curve_widths = PLOT_LINE_WIDTH_PALETTE[1 + (seq_len(Ncurves)-1)%%length(PLOT_LINE_WIDTH_PALETTE)] }
  else if(length(curve_widths)==1){ curve_widths = rep(curve_widths,times=Ncurves); }
  if(is.null(show_curves)) show_curves = rep(TRUE,Ncurves)
  
  if(!is.null(resolution)){
    if(verbose) cat(sprintf("Downsampling curves to %d time points..\n",resolution))
    for(n in seq_len(Ncurves)){
      if(sum(!is.na(curves[[n]][[2]]))<2) next;
      X = curves[[n]][[1]];
      if(length(X)>resolution){
        Xnew = seq(from=X[1], to=tail(X,1), length.out=resolution);
        curves[[n]][[1]] = Xnew;
        curves[[n]][[2]] = approx(x=X, y=curves[[n]][[2]], xout=Xnew, method="linear", yleft=NaN, yright=NaN, rule = 1, f = 0, ties = mean)$y;
      }
    }
  }
  
  # save data as text file
  if(save_data){
    if(verbose) cat(sprintf("%sSaving %s to table (%s)..\n",verbose_prefix,ytype,case_tag))
    cat("file_basepath is ", file_basepath, "\n")
    output_table=sprintf("%s.tsv",file_basepath)
    cat("Output table is ", output_table, "\n")
    check_output_file(output_table,force_replace=FALSE,TRUE,"  ")
    cat(sprintf("# %s over %s (%s)\n# %s\n#\n",ylabel,xlabel,case_tag,stringr::str_replace_all(data_file_comments,"\n","\n# ")), file=output_table, append=TRUE)
    for(n in seq_len(Ncurves)){
      cat(sprintf("# %s\n%s\t%s\n",curve_names[[n]],xtype,ytype), file=output_table, append=TRUE);
      write.table(x=cbind(curves[[n]][[1]],curves[[n]][[2]]), file=output_table, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE);
      cat(sprintf("\n\n"), file=output_table, append=TRUE);
    } 
  }
  
  # omit invalid values from curves
  for(n in seq_len(Ncurves)){
    if(plot_log_values){
      valids = which(is.finite(log(curves[[n]][[2]])));
    }else{
      valids = which(is.finite(curves[[n]][[2]]));
    }
    curves[[n]][[1]] = curves[[n]][[1]][valids]
    curves[[n]][[2]] = curves[[n]][[2]][valids]
  }
  
  # determine X & Y ranges
  data_miny = NA
  data_maxy = NA
  data_minx = NA
  data_maxx = NA
  for(n in seq_len(Ncurves)){
    X = curves[[n]][[1]]
    Y = curves[[n]][[2]]
    if(!is.null(plot_minx)) Y = Y[X>=plot_minx]
    if(!is.null(plot_maxx)) Y = Y[X<=plot_maxx]
    if(plot_log_values) Y = Y[Y>0]
    if(length(Y)>0){
      data_miny = (if(is.na(data_miny)) min(Y,na.rm=TRUE) else min(data_miny,min(Y,na.rm=TRUE)));
      data_maxy = (if(is.na(data_maxy)) max(Y,na.rm=TRUE) else max(data_maxy,max(Y,na.rm=TRUE)));
      data_minx = (if(is.na(data_minx)) min(X,na.rm=TRUE) else min(data_minx,min(X,na.rm=TRUE)));
      data_maxx = (if(is.na(data_maxx)) max(X,na.rm=TRUE) else max(data_maxx,max(X,na.rm=TRUE)));
    }
  }
  plot_minx	= (if((!is.null(plot_minx)) && is.finite(plot_minx)) plot_minx else data_minx)
  plot_maxx	= (if((!is.null(plot_maxx)) && is.finite(plot_maxx)) plot_maxx else data_maxx)
  plot_miny 	= (if((!is.null(plot_miny)) && is.finite(plot_miny)) plot_miny else data_miny)
  plot_maxy 	= (if((!is.null(plot_maxy)) && is.finite(plot_maxy)) plot_maxy else data_maxy + (data_maxy-plot_miny)*0.1)
  if(is.na(data_miny) || is.na(data_maxy) || is.na(data_minx) || is.na(data_maxx)){
    if(verbose) cat(sprintf("%sWARNING: No valid points for plotting\n",verbose_prefix))
    return();
  }
  
  # prepare for plotting
  if(verbose) cat(sprintf("%sPlotting %s over age (%s)..\n",verbose_prefix,ytype,case_tag))
  curves 				= curves[show_curves]
  curve_names			= curve_names[show_curves]
  curve_colors		= curve_colors[show_curves]
  curve_line_types	= curve_line_types[show_curves]
  curve_widths		= curve_widths[show_curves]
  plot_file 			= sprintf("%s%s.pdf",file_basepath,case_tag)
  Ncurves				= length(curves)
  check_output_file(plot_file,force_replace=FALSE,TRUE,"  ")
  plot_width			= (if(is.null(plot_width)) DEFAULT_PLOT_WIDTH else plot_width)
  plot_height			= (if(is.null(plot_height)) DEFAULT_PLOT_HEIGHT else plot_height)
  margins				= c((if(legend_pos=="none") 1 else max(1,Ncurves*0.23-plot_height)),1,DEFAULT_PLOT_TOP_MARGIN,(if(legend_pos=="outside") DEFAULT_PLOT_RIGHT_MARGIN_WITH_LEGEND else 0.5))
  pdf(file=plot_file, width=plot_width+margins[2]+margins[4], height=plot_height+margins[1]+margins[3])
  par(mai=margins)
  
  # initialize with empty plot
  suppressWarnings(plot(	x		= c(),
                         y		= c(), 
                         main 	= plot_title, 
                         xlab 	= xlabel, 
                         ylab 	= NA, 
                         log 	= (if(plot_log_values) "y" else ""),
                         yaxt	= (if(plot_log_values) "n" else NULL), 
                         cex		= 1.1,
                         las		= 1,
                         xaxs	= "i",
                         xlim	= (if(reverse_x) c(plot_maxx, plot_minx) else c(plot_minx,plot_maxx)), 
                         ylim	= c(plot_miny, plot_maxy)))
  title(ylab=ylabel, line=4)
  
  # plot shadings
  if(!is.null(shadings)){
    for(sh in seq_len(length(shadings))){
      curve1 		= shadings[[sh]][[1]]
      curve2 		= shadings[[sh]][[2]]
      shade_color = shadings[[sh]][[3]]
      order1		= order(curves[[curve1]][[1]])
      order2		= order(curves[[curve2]][[1]])
      polygon(c(curves[[curve1]][[1]][order1],rev(curves[[curve2]][[1]][order2])),
              c(curves[[curve2]][[2]][order2],rev(curves[[curve1]][[2]][order1])), 
              col=shade_color,
              border=NA)
    }
  }
  
  # plot reference curve
  if(!is.null(horizontal_reference)){
    lines(	x	= c(plot_minx,plot_maxx),
           y	= c(horizontal_reference,horizontal_reference),
           type= "l",
           col	= REFERENCE_CURVE_COLOR,
           lty	= 1,
           lwd	= 1)
  }
  
  # plot curves
  for(n in seq_len(Ncurves)){
    lines(	x	= curves[[n]][[1]],
           y	= curves[[n]][[2]],
           type= "l",
           col	= curve_colors[n],
           lty	= curve_line_types[n],
           lwd	= curve_widths[n]);				
  }
  
  # plot scatterpoints
  if((!is.null(scatterpoints)) && (nrow(scatterpoints)>0)){
    points(x=scatterpoints[,1], y=scatterpoints[,2], pch=21, col="#303030", bg="#909090")
    if(ncol(scatterpoints)==4){
      # include vertical error bars
      arrows(scatterpoints[,1], scatterpoints[,3], scatterpoints[,1], scatterpoints[,4], length=0.05, angle=90, code=3, col="#303030")
    }else if(ncol(scatterpoints)==6){
      # include horizontal & vertical error bars
      arrows(scatterpoints[,3], scatterpoints[,2], scatterpoints[,4], scatterpoints[,2], length=0.05, angle=90, code=3, col="#303030") # horizontal bars
      arrows(scatterpoints[,1], scatterpoints[,5], scatterpoints[,1], scatterpoints[,6], length=0.05, angle=90, code=3, col="#303030") # vertical bars
    }
  }
  
  # add legend
  if(legend_pos!="none"){
    if(legend_pos=="outside"){
      if(reverse_x){
        legendx = plot_minx - 0.2*(plot_maxx-plot_minx)
      }else{
        legendx = plot_maxx + 0.2*(plot_maxx-plot_minx)
      }
      legendy = plot_maxy
      legend(x=legendx, y=legendy, legend = curve_names, col=curve_colors, lty=curve_line_types, lwd=curve_widths, xpd=NA);
    }else{
      legend(legend_pos, legend = curve_names, col=curve_colors, lty=curve_line_types, lwd=curve_widths);
    }
  }
  if(plot_log_values){
    # improve appearance of y-axis ticks if logarithmic
    aty = axTicks(2)
    axis(2,at=aty,labels=sapply(aty, function(x) sprintf("%g",x)), las=1)
  }
  invisible(dev.off());
}


pairwise_scatterplots = function(	output_basepath,	# e.g. 'output/all_simulations'
                                  case_tag,			# character, will be included as subtitle in the plots
                                  scattervalues, 		# 2D numeric matrix of size NP x ND, storing NP scattered values for each of ND data types. Each scatterplot will show those values for a pair of data types, e.g. Y[,i] vs Y[,j]
                                  data_names,			# character vector of size ND. Can be NULL, in which case the column names of scattervalues are used as names.
                                  point_names,		# character vector of size NP. Can be NULL, in which case the row names of scattervalues are used as names.
                                  logarithmic,		# boolean vector of size ND, indicating if data types should be plotted on a log axis. Can also be a single boolean, applying to all data types
                                  include_diagonal,	# (boolean)
                                  comment_line,		# character, optional comment line to include in TSV file
                                  verbose,
                                  verbose_prefix){
  ND 	= ncol(scattervalues);
  NP	= nrow(scattervalues);
  if(is.null(data_names)) data_names = colnames(scattervalues)
  if(is.null(point_names)){
    point_names = rownames(scattervalues)
  }else{
    rownames(scattervalues) = point_names
  }
  if(length(logarithmic)==1) logarithmic  = rep(logarithmic,times=ND)
  
  # save to file
  if(verbose) cat(sprintf("%sSaving scattervalues for '%s'..\n",verbose_prefix,case_tag))
  output_table=(if(endsWith(output_basepath,"/")) sprintf("%s/all_data.tsv",output_basepath) else sprintf("%s.tsv",output_basepath))
  check_output_file(output_table,force_replace=FALSE,TRUE,"  ")
  cat(sprintf("# Scatterdata for %d variables, '%s'\n# %s\n#\n# \t%s\n",ND,case_tag,comment_line,paste(data_names,collapse="\t")), file=output_table, append=FALSE)
  write.table(x=scattervalues, file=output_table, append=TRUE, sep="\t", row.names=(!is.null(point_names)), col.names=FALSE, quote=FALSE);
  
  # plot to PDF
  if(ND<=1) return(); # nothing to plot
  for(d1 in 1:ND){
    for(d2 in 1:ND){
      if(d1<=d2) next;
      if(verbose) cat(sprintf("%sGenerating scatterplot '%s' vs '%s', %s..\n",verbose_prefix,data_names[d1],data_names[d2],case_tag))
      valids = which(is.finite(scattervalues[,d1]) & is.finite(scattervalues[,d2]))
      if(length(valids)==0){
        if(verbose) cat(sprintf("%sWARNING: No valid points to plot, skipping scatterplot\n",verbose_prefix))
        next;
      }
      if((ND==2) && (!endsWith(output_basepath,"/"))){
        plot_file = sprintf("%s%s.pdf",output_basepath, case_tag)
      }else{
        plot_file = sprintf("%s%s%s%s_vs_%s.pdf",output_basepath, case_tag,(if(endsWith(output_basepath,"/")) "" else "_"),data_names[d1],data_names[d2])
      }
      check_output_file(plot_file,force_replace=FALSE,TRUE,"  ")
      margins = c(1,1,DEFAULT_PLOT_TOP_MARGIN,0.5)
      pdf(file=plot_file, width=DEFAULT_PLOT_WIDTH+margins[2]+margins[4], height=DEFAULT_PLOT_HEIGHT+margins[1]+margins[3])
      par(mai=margins)
      X	= scattervalues[valids,d1]
      Y	= scattervalues[valids,d2]
      graphics::plot(	x = X,
                      y = Y,
                      type = "p",
                      main = sprintf("%s vs %s\n%s",data_names[d1],data_names[d2],case_tag),
                      xlab=sprintf("%s",data_names[d1]),
                      ylab=sprintf("%s",data_names[d2]),
                      col = "#303030",
                      pch=1,
                      xlim=(if(logarithmic[d1]) NULL else expand_range(min(X),max(X),1.1)),
                      ylim=(if(logarithmic[d2]) NULL else expand_range(min(Y),max(Y),1.1)),
                      log=paste0((if(logarithmic[d1]) "x" else ""), (if(logarithmic[d2]) "y" else "")),
                      xaxt=(if(logarithmic[d1]) "n" else NULL), 
                      yaxt=(if(logarithmic[d2]) "n" else NULL));
      # add diagonal
      if(include_diagonal){
        graphics::abline(a=0, b=1, col="#909090")
      }
      # improve axis tics if logarithmic
      if(logarithmic[d1]){
        atx = axTicks(1)
        axis(1,at=atx,labels=sapply(atx, function(x) sprintf("%g",x)))
      }
      if(logarithmic[d2]){
        aty = axTicks(2)
        axis(2,at=aty,labels=sapply(aty, function(y) sprintf("%g",y)), las=2)
      }
      invisible(dev.off());
    }
  }
}


compare_tree_to_model = function(	plot_dir,
                                  case_tag,		# e.g. 'BACTERIA_EMBL - fitted model, grid size 10'
                                  subtitle,
                                  tree_name,
                                  model_name,
                                  lambda = NULL,
                                  mu = NULL,
                                  psi = NULL,
                                  sim,		# deterministic simulation of the compared model
                                  tree_LTT,	# list containing times[] and lineages[], each of length NT
                                  root_age,
                                  time_unit,		# e.g. 'Myr'
                                  verbose,
                                  verbose_prefix){
  if(!is.null(subtitle)){ subtitle = sprintf("%s\n%s",case_tag,subtitle); }else{ subtitle = case_tag; }
  plot_curves(file_basepath		= sprintf("%s/LTT%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'LTTs',
              case_tag			= case_tag,
              curves				= list(list(tree_LTT$ages,tree_LTT$lineages), list(sim$ages,sim$LTT)),
              curve_names			= c(tree_name, model_name),
              curve_colors		= c(BLACK_CURVE_COLOR, BLUE_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= "lineages",
              plot_log_values		= FALSE,
              legend_pos			= "outside",
              plot_title			= sprintf("LTT of tree (%s)\nvs. model (%s)\n%s",tree_name,model_name,subtitle),
              data_file_comments	= sprintf("LTT of tree (%s) vs dLTT of model (%s)",tree_name,model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  tree_LTT_AUC = sum(0.5 * (tree_LTT$lineages[2:length(tree_LTT$ages)]+tree_LTT$lineages[1:(length(tree_LTT$ages)-1)]) * abs(diff(tree_LTT$ages)))
  tree_LTT$nLTT = tree_LTT$lineages/tree_LTT_AUC
  plot_curves(file_basepath		= sprintf("%s/nLTT%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'nLTTs',
              case_tag			= case_tag,
              curves				= list(list(tree_LTT$ages,tree_LTT$nLTT), list(sim$ages,sim$nLTT)),
              curve_names			= c(tree_name, model_name),
              curve_colors		= c(BLACK_CURVE_COLOR, BLUE_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= "density",
              plot_log_values		= FALSE,
              legend_pos			= "outside",
              plot_title			= sprintf("nLTT of tree (%s)\nvs. model (%s)\n%s",tree_name,model_name,subtitle),
              data_file_comments	= sprintf("Normalized LTT of tree (%s) vs dnLTT of model (%s)",tree_name,model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/lambda%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'lambda',
              case_tag			= case_tag,
              curves				= list(list(tree_LTT$ages, lambda), list(sim$ages,sim$lambda)),
              curve_names			= c("tree lambda", "model lambda"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLUE_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("lambda (1/%s)",time_unit),
              plot_log_values		= FALSE,
              legend_pos			= "none",
              plot_title			= sprintf("Lambda of (%s)\n vs model (%s)%s",tree_name,model_name,subtitle),
              data_file_comments	= sprintf("Speciation rate (lambda) of tree(%s) vs model (%s)",tree_name,model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/mu%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'mu',
              case_tag			= case_tag,
              curves				= list(list(tree_LTT$ages, mu),list(sim$ages,sim$mu)),
              curve_names			= c("tree mu", "model mu"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLUE_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("mu (1/%s)",time_unit),
              plot_log_values		= FALSE,
              legend_pos			= "none",
              plot_title			= sprintf("Tree mu (%s) vs model mu (%s)\n%s",tree_name,model_name,subtitle),
              data_file_comments	= sprintf("Extinction rate (mu) of tree (%s) vs. model (%s)",tree_name,model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/psi%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'psi',
              case_tag			= case_tag,
              curves				= list(list(tree_LTT$ages, psi),list(sim$ages,sim$psi)),
              curve_names			= c("tree_psi", "model psi"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLUE_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("psi (1/%s)",time_unit),
              plot_log_values		= FALSE,
              legend_pos			= "none",
              plot_title			= sprintf("Tree psi (%s) vs. Model psi (%s)\n%s",tree_name, model_name,subtitle),
              data_file_comments	= sprintf("Sampling rate (psi) of tree (%s) vs model (%s)",tree_name, model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/PSR%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'PSR',
              case_tag			= case_tag,
              curves				= list(list(sim$ages,sim$PSR)),
              curve_names			= c("model PSR"),
              curve_colors		= c(BLUE_CURVE_COLOR),
              curve_line_types	= c(1),
              curve_widths		= c(1.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("PSR (1/%s)",time_unit),
              plot_log_values		= FALSE,
              legend_pos			= "none",
              plot_title			= sprintf("Model PSR (%s)\n%s",model_name,subtitle),
              data_file_comments	= sprintf("Pulled speciation rate (PSR) of model (%s)",model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/PDR%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'PDR',
              case_tag			= case_tag,
              curves				= list(list(sim$ages,sim$PDR)),
              curve_names			= c("model PDR"),
              curve_colors		= c(BLUE_CURVE_COLOR),
              curve_line_types	= c(1),
              curve_widths		= c(1.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= (if(any(sim$PDR<0)) NULL else 0),
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("PDR (1/%s)",time_unit),
              plot_log_values		= FALSE,
              legend_pos			= "none",
              plot_title			= sprintf("Model PDR (%s)\n%s",model_name,subtitle),
              data_file_comments	= sprintf("Pulled diversification rate (PDR) of model (%s)",model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/IPDR%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'IPDR',
              case_tag			= case_tag,
              curves				= list(list(sim$ages,sim$IPDR)),
              curve_names			= c("model IPDR"),
              curve_colors		= c(BLUE_CURVE_COLOR),
              curve_line_types	= c(1),
              curve_widths		= c(1.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("IPDR"),
              plot_log_values		= FALSE,
              legend_pos			= "none",
              plot_title			= sprintf("Integrated PDR (%s)\n%s",model_name,subtitle),
              data_file_comments	= sprintf("Age-integrated pulled diversification rate (IPDR) of model (%s)",model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  tree_LTT$Reff = lambda/(psi+mu)
  plot_curves(file_basepath		= sprintf("%s/Reff%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'Reff',
              case_tag			= case_tag,
              curves				= list(list(tree_LTT$ages, tree_LTT$Reff), list(sim$ages,sim$Reff)),
              curve_names			= c("tree Ref", "model Reff"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLUE_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= "Reff",
              plot_log_values		= FALSE,
              legend_pos			= "none",
              plot_title			= sprintf("Tree Reff (%s) vs. model Reff (%s)\n%s",tree_name,model_name,subtitle),
              data_file_comments	= sprintf("Effective reproduction ratio (Reff) of tree (%s) vs. model (%s)",tree_name,model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  tree_LTT$sampling_proportion = psi/(mu + psi)
  plot_curves(file_basepath		= sprintf("%s/sampling_proportion%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'sampling_proportion',
              case_tag			= case_tag,
              curves				= list(list(tree_LTT$ages, tree_LTT$sampling_proportion),list(sim$ages,sim$sampling_proportion)),
              curve_names			= c("tree sampling proportion","model sampling_proportion"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLUE_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= "sampling_proportion",
              plot_log_values		= FALSE,
              legend_pos			= "none",
              plot_title			= sprintf("Tree sampling_proportion (%s) vs. model (%s)\n%s",tree_name,model_name,subtitle),
              data_file_comments	= sprintf("Sampling proportion (psi/(mu+psi)) of tree (%s) vs. model (%s)",tree_name,model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  tree_LTT$removal_rate = mu + psi
  plot_curves(file_basepath		= sprintf("%s/removal_rate%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'removal_rate',
              case_tag			= case_tag,
              curves				= list(list(tree_LTT$ages, tree_LTT$removal_rate),list(sim$ages,sim$removal_rate)),
              curve_names			= c("tree_removal_rate","model removal_rate"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLUE_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("removal_rate (1/%s)",time_unit),
              plot_log_values		= FALSE,
              legend_pos			= "none",
              plot_title			= sprintf("Tree removal rate (%s) vs. model (%s)\n%s",tree_name,model_name,subtitle),
              data_file_comments	= sprintf("Removal rate (aka. become-uninfectious rate) of tree (%s) vs. model (%s)",tree_name,model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  tree_LTT$lambda_psi = lambda*psi
  plot_curves(file_basepath		= sprintf("%s/lambda_psi%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'event_density',
              case_tag			= case_tag,
              curves				= list(list(tree_LTT$ages, tree_LTT$lambda_psi),list(sim$ages,sim$lambda_psi)),
              curve_names			= c("tree event_density","model event_density"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLUE_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("event_density (1/%s^2)",time_unit),
              plot_log_values		= FALSE,
              legend_pos			= "none",
              plot_title			= sprintf("Tree event_density lambda*psi (%s) vs. model (%s)\n%s",tree_name,model_name,subtitle),
              data_file_comments	= sprintf("Event density (lambda*psi) of tree (%s) vs. model (%s)",tree_name,model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
}



# Show true and fitted model parameters on plot
plot_fitted_vs_true_model = function(	plot_dir,
                                      case_tag,		# e.g. 'BACTERIA_EMBL - fitted model, grid size 10'
                                      subtitle,
                                      true_model_name,
                                      fit_model_name,
                                      sim_true,		# deterministic simulation of the true model
                                      sim_fit,		# deterministic simulation of the fitted model
                                      tree_LTT,	# list containing times[] and lineages[], each of length NT
                                      root_age,
                                      time_unit,		# e.g. 'Myr'
                                      verbose,
                                      verbose_prefix){
  if(!is.null(subtitle)){ subtitle = sprintf("%s\n%s",case_tag,subtitle); }else{ subtitle = case_tag; }
  plot_curves(file_basepath		= sprintf("%s/LTT%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'LTTs',
              case_tag			= case_tag,
              curves				= list(	list(sim_true$ages,sim_true$LTT),
                                list(tree_LTT$ages,tree_LTT$lineages),
                                list(sim_fit$ages,sim_fit$LTT)),
              curve_names			= c("true", "tree", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2,1),
              curve_widths		= 1.5,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= "lineages",
              legend_pos			= "outside",
              plot_title			= sprintf("LTT of true tree (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("LTT of tree generated by true model (%s) vs deterministic LTT of fit piecewise model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  tree_LTT_AUC = sum(0.5 * (tree_LTT$lineages[2:length(tree_LTT$ages)]+tree_LTT$lineages[1:(length(tree_LTT$ages)-1)]) * abs(diff(tree_LTT$ages)))
  tree_LTT$nLTT = tree_LTT$lineages/tree_LTT_AUC
  plot_curves(file_basepath		= sprintf("%s/nLTT_vs_true_vs_tree%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'nLTTs',
              case_tag			= case_tag,
              curves				= list(	list(sim_true$ages,sim_true$nLTT), 
                                list(tree_LTT$ages,tree_LTT$nLTT), 
                                list(sim_fit$ages,sim_fit$nLTT)),
              curve_names			= c("true", "tree", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, GREY_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,1,2),
              curve_widths		= c(1.5,2,1.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= "density",
              legend_pos			= "outside",
              plot_title			= sprintf("nLTT of true tree (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Normalized LTT of tree generated by true model (%s) vs deterministic nLTT of fit piecewise model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/nLTT%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'nLTTs',
              case_tag			= case_tag,
              curves				= list(	list(sim_true$ages,sim_true$nLTT), 
                                list(sim_fit$ages,sim_fit$nLTT)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= 1.5,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= "density",
              legend_pos			= "outside",
              plot_title			= sprintf("nLTT of true tree (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Normalized LTT of tree generated by true model (%s) vs deterministic nLTT of fit piecewise model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/lambda%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'lambda',
              case_tag			= case_tag,
              curves				= list(list(sim_true$ages,sim_true$lambda), list(sim_fit$ages,sim_fit$lambda)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("lambda (1/%s)",time_unit),
              legend_pos			= "outside",
              plot_title			= sprintf("True lambda (%s)\nvs fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Speciation rate (lambda) of true model (%s) vs lambda of fit model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/mu%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'mu',
              case_tag			= case_tag,
              curves				= list(list(sim_true$ages,sim_true$mu), list(sim_fit$ages,sim_fit$mu)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("mu (1/%s)",time_unit),
              legend_pos			= "outside",
              plot_title			= sprintf("True mu (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Extinction rate (mu) of true model (%s) vs mu of fit model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/psi%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'psi',
              case_tag			= case_tag,
              curves				= list(list(sim_true$ages,sim_true$psi), list(sim_fit$ages,sim_fit$psi)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("psi (1/%s)",time_unit),
              legend_pos			= "outside",
              plot_title			= sprintf("True psi (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Sampling rate (psi) of true model (%s) vs psi of fit model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/removal_rate%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'removal_rate',
              case_tag			= case_tag,
              curves				= list(list(sim_true$ages,sim_true$removal_rate), list(sim_fit$ages,sim_fit$removal_rate)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("removal_rate (1/%s)",time_unit),
              legend_pos			= "outside",
              plot_title			= sprintf("True removal_rate (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Removal rate (aka. become-uninfectious rate, mu+psi) of true model (%s) vs fit model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/PSR%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'PSR',
              case_tag			= case_tag,
              curves				= list(list(sim_true$ages,sim_true$PSR), list(sim_fit$ages,sim_fit$PSR)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("PSR (1/%s)",time_unit),
              legend_pos			= "outside",
              plot_title			= sprintf("True PSR (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Pulled speciation rate (PSR) of true model (%s) vs fit model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/PDR%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'PDR',
              case_tag			= case_tag,
              curves				= list(list(sim_true$ages,sim_true$PDR), list(sim_fit$ages,sim_fit$PDR)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= (if(any(sim_true$PDR<0) || any(sim_fit$PDR<0)) NULL else 0),
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("PDR (1/%s)",time_unit),
              legend_pos			= "outside",
              plot_title			= sprintf("True PDR (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Pulled diversification rate (PDR) of true model (%s) vs fit model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/IPDR%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'IPDR',
              case_tag			= case_tag,
              curves				= list(list(sim_true$ages,sim_true$IPDR), list(sim_fit$ages,sim_fit$IPDR)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("IPDR"),
              legend_pos			= "outside",
              plot_title			= sprintf("True IPDR (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Age-integrated pulled diversification rate (IPDR) of true model (%s) vs fit model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")
  plot_curves(file_basepath		= sprintf("%s/Reff%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'Reff',
              case_tag			= case_tag,
              curves				= list(list(sim_true$ages,sim_true$Reff), list(sim_fit$ages,sim_fit$Reff)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= "Reff",
              legend_pos			= "outside",
              plot_title			= sprintf("True Reff (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Effective reproduction ratio (Reff) of true model (%s) vs Reff of fit model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")	
  plot_curves(file_basepath		= sprintf("%s/sampling_proportion%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'sampling_proportion',
              case_tag			= case_tag,
              curves				= list(list(sim_true$ages,sim_true$sampling_proportion), list(sim_fit$ages,sim_fit$sampling_proportion)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= "sampling_proportion",
              legend_pos			= "outside",
              plot_title			= sprintf("True sampling_proportion (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Sampling proportion of true model (%s) vs fit model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")	
  plot_curves(file_basepath		= sprintf("%s/lambda_psi%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'lambda_psi',
              case_tag			= case_tag,
              curves				= list(list(sim_true$ages,sim_true$lambda_psi), list(sim_fit$ages,sim_fit$lambda_psi)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("lambda_psi (1/%s^2)",time_unit),
              legend_pos			= "outside",
              plot_title			= sprintf("True lambda*psi (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("lambda*psi of true model (%s) vs fit model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")	
  plot_curves(file_basepath		= sprintf("%s/Pmissing%s",plot_dir, case_tag),
              xtype				= 'age',
              ytype				= 'Pmissing',
              case_tag			= case_tag,
              curves				= list(list(sim_true$ages,sim_true$Pmissing), list(sim_fit$ages,sim_fit$Pmissing)),
              curve_names			= c("true", "fit"),
              curve_colors		= c(BLUE_CURVE_COLOR, BLACK_CURVE_COLOR),
              curve_line_types	= c(1,2),
              curve_widths		= c(1.5,2.5),
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= root_age,
              plot_miny			= 0,
              plot_maxy			= 1,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_unit),
              ylabel				= sprintf("Pmissing"),
              legend_pos			= "outside",
              plot_title			= sprintf("True Pmissing (%s)\nvs. fit model (%s)\n%s",true_model_name,fit_model_name,subtitle),
              data_file_comments	= sprintf("Pmissing of true model (%s) vs fit model (%s)",true_model_name,fit_model_name),
              verbose 			= TRUE,
              verbose_prefix 		= "  ")	
}



# Plot deterministic values for specified model
plot_model = function(	model_name,			# (string) e.g. 'exp_lambda_const_mu'
                       sim,				# deterministic simulation of the model
                       plot_maxx,			# (numeric) maximum age to plot, typically the span of the tree
                       plot_basepath,		# (string)
                       time_units,			# (string)
                       verbose,
                       verbose_prefix){	# (string)
  if(verbose) cat2(sprintf("%sPlotting LTT of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%sLTT",plot_basepath),
              xtype				= 'age',
              ytype				= 'LTT',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$LTT)),
              curve_names			= model_name,
              curve_colors		= BLUE_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= "lineages",
              legend_pos			= "none",
              plot_title			= sprintf("LTT\n%s",model_name),
              data_file_comments	= sprintf("LTT of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))
  if(verbose) cat2(sprintf("%sPlotting nLTT of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%snLTT",plot_basepath),
              xtype				= 'age',
              ytype				= 'nLTT',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$nLTT)),
              curve_names			= model_name,
              curve_colors		= BLUE_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= "density",
              legend_pos			= "none",
              plot_title			= sprintf("nLTT\n%s",model_name),
              data_file_comments	= sprintf("Normalized LTT of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))
  if(verbose) cat2(sprintf("%sPlotting PSR of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%sPSR",plot_basepath),
              xtype				= 'age',
              ytype				= 'PSR',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$PSR)),
              curve_names			= model_name,
              curve_colors		= BLUE_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= sprintf("PSR (1/%s)",time_units),
              legend_pos			= "none",
              plot_title			= sprintf("Pulled speciation rate\n%s",model_name),
              data_file_comments	= sprintf("PSR of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))
  if(verbose) cat2(sprintf("%sPlotting PDR of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%sPDR",plot_basepath),
              xtype				= 'age',
              ytype				= 'PDR',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$PDR)),
              curve_names			= model_name,
              curve_colors		= BLACK_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= (if(any(na.omit(sim$PDR)<0)) NULL else 0),
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= sprintf("PDR (1/%s)",time_units),
              legend_pos			= "none",
              plot_title			= sprintf("Pulled diversification rate\n%s",model_name),
              data_file_comments	= sprintf("PDR of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))
  if(verbose) cat2(sprintf("%sPlotting IPDR of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%sIPDR",plot_basepath),
              xtype				= 'age',
              ytype				= 'IPDR',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$IPDR)),
              curve_names			= model_name,
              curve_colors		= BLUE_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= sprintf("IPDR"),
              legend_pos			= "none",
              plot_title			= sprintf("Integrated pulled diversification rate\n%s",model_name),
              data_file_comments	= sprintf("Age-integrated pulled diversification rate (IPDR) of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))
  if(verbose) cat2(sprintf("%sPlotting event density of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%slambda_psi",plot_basepath),
              xtype				= 'age',
              ytype				= 'event_density',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$lambda_psi)),
              curve_names			= model_name,
              curve_colors		= BLUE_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= sprintf("event density (1/%s^2)",time_units),
              legend_pos			= "none",
              plot_title			= sprintf("Event density (lambda*psi)\n%s",model_name),
              data_file_comments	= sprintf("Event density (lambda*psi) of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))
  if(verbose) cat2(sprintf("%sPlotting Reff of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%sReff",plot_basepath),
              xtype				= 'age',
              ytype				= 'Reff',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$Reff)),
              curve_names			= model_name,
              curve_colors		= BLUE_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= "Reff",
              legend_pos			= "none",
              plot_title			= sprintf("Reff of model '%s'",model_name),
              data_file_comments	= sprintf("Reff of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))
  if(verbose) cat2(sprintf("%sPlotting removal rate of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%sremoval_rate",plot_basepath),
              xtype				= 'age',
              ytype				= 'removal_rate',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$removal_rate)),
              curve_names			= model_name,
              curve_colors		= BLUE_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= sprintf("removal rate (1/%s)",time_units),
              legend_pos			= "none",
              plot_title			= sprintf("Removal rate (mu+psi)\n%s",model_name),
              data_file_comments	= sprintf("Removal rate (mu+psi) of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))
  if(verbose) cat2(sprintf("%sPlotting sampling proportion of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%ssampling_proportion",plot_basepath),
              xtype				= 'age',
              ytype				= 'sampling_proportion',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$psi/(sim$mu+sim$psi))),
              curve_names			= model_name,
              curve_colors		= BLUE_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= "sampling proportion",
              legend_pos			= "none",
              plot_title			= sprintf("Sampling proportion\n%s",model_name),
              data_file_comments	= sprintf("Sampling proportion (psi/(mu+psi)) of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))
  if(verbose) cat2(sprintf("%sPlotting lambda of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%slambda",plot_basepath),
              xtype				= 'age',
              ytype				= 'lambda',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$lambda)),
              curve_names			= model_name,
              curve_colors		= BLUE_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= sprintf("lambda (1/%s)",time_units),
              legend_pos			= "none",
              plot_title			= sprintf("speciation rate (lambda)\n%s",model_name),
              data_file_comments	= sprintf("lambda of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))
  if(verbose) cat2(sprintf("%sPlotting psi of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%spsi",plot_basepath),
              xtype				= 'age',
              ytype				= 'psi',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$psi)),
              curve_names			= model_name,
              curve_colors		= BLUE_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= sprintf("psi (1/%s)",time_units),
              legend_pos			= "none",
              plot_title			= sprintf("sampling rate (psi)\n%s",model_name),
              data_file_comments	= sprintf("psi of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))	
  if(verbose) cat2(sprintf("%sPlotting mu of model '%s'..\n",verbose_prefix,model_name))
  plot_curves(file_basepath		= sprintf("%smu",plot_basepath),
              xtype				= 'age',
              ytype				= 'mu',
              case_tag			= sprintf("%s",model_name),
              curves				= list(list(sim$ages,sim$mu)),
              curve_names			= model_name,
              curve_colors		= BLUE_CURVE_COLOR,
              curve_line_types	= 1,
              curve_widths		= 2,
              reverse_x			= TRUE,
              plot_minx			= 0,
              plot_maxx			= plot_maxx,
              plot_miny			= 0,
              resolution			= PLOT_DOWNSAMPLING_RESOLUTION,
              xlabel				= sprintf("age (%s)",time_units),
              ylabel				= sprintf("mu (1/%s)",time_units),
              legend_pos			= "none",
              plot_title			= sprintf("extinction rate (mu)\n%s",model_name),
              data_file_comments	= sprintf("mu of model '%s'",model_name),
              save_data			= TRUE,
              verbose 			= TRUE,
              verbose_prefix 		= paste0(verbose_prefix,"  "))					
}



# assess the adequacy of one or more (presumably fitted) HBDS models to explaining a given (presumably real) tree
assess_model_adequacy = function(	tree, 
                                  models,
                                  splines_degree,
                                  tree_name,
                                  model_name,
                                  Nbootstraps,
                                  report_file = "",
                                  max_extant_tips = NULL,
                                  Nthreads = 1){
  root_age = get_tree_span(tree)$max_distance
  adequacy = model_adequacy_hbds(	tree, 
                                          models				= models,
                                          splines_degree		= 1,
                                          Nbootstraps			= Nbootstraps,
                                          Nthreads			= Nthreads,
                                          extrapolate			= TRUE,
                                          max_extant_tips		= max_extant_tips,
                                          max_model_runtime	= MODEL_ADEQUACY_MAX_RUNTIME)
  if(!adequacy$success) return(list(success=FALSE, error=sprintf("Failed to test model adequacy: %s",adequacy$error)))
  
  # save report to a text file
  if(report_file!=""){
    dir.create(dirname(report_file), showWarnings = FALSE, recursive=TRUE)
    sink(file=report_file); print(adequacy); sink();
  }				
  return(adequacy)
}


###############################								
# PREPARATIONS


cat(sprintf("Loading %d required packages (%s)..\n",length(REQUIRED_PACKAGES),paste(REQUIRED_PACKAGES,collapse=", ")));
for(p in 1:length(REQUIRED_PACKAGES)){
  if(dir.exists(REQUIRED_PACKAGES[p])){
    # load this package with devtools
    cat(sprintf("  Loading '%s' via devtools..\n",REQUIRED_PACKAGES[p]))
    library("devtools", warn.conflicts = FALSE);
    suppressWarnings(suppressPackageStartupMessages(devtools::load_all(REQUIRED_PACKAGES[p],quiet=TRUE)))
  }else if(!suppressMessages(suppressPackageStartupMessages(require(REQUIRED_PACKAGES[p], quietly=TRUE, character.only=TRUE)))){
    cat(sprintf("  Note: Installing missing package '%s'..\n",REQUIRED_PACKAGES[p]))
    install.packages(REQUIRED_PACKAGES[p], dependencies=TRUE, repos="http://cran.r-project.org/");
    suppressMessages(suppressPackageStartupMessages(require(REQUIRED_PACKAGES[p],quietly=TRUE,character.only=TRUE,warn.conflicts = FALSE)))
  }
}


# prepare output dir
output_dir = get_non_existent_dir("output", "run_", 3);
dir.create(output_dir, showWarnings = FALSE, recursive=TRUE);
cat(sprintf("All output will be written to '%s'..\n",output_dir))


# save backup of this script
initial.options = commandArgs(trailingOnly = FALSE);
script_path <- sub("--file=", "", initial.options[grep("--file=", initial.options)]);
OK = file.copy(script_path, to = file.path(output_dir,"workflow_backup.R"));

# prepare log file
display_date_time = sprintf("%s, %s (%s)",format(Sys.Date(), "%B %d, %Y"),format(Sys.time(), "%H:%M"),Sys.timezone())
logfile = sprintf("%s/log.txt",output_dir)
cat(sprintf("Log of HBDS workflow\nGenerated on: %s\n",display_date_time),file=logfile,append=TRUE)
cat2 = function(message, file=logfile, append=TRUE){
  cat(message);
  cat(message, file=file, append=append);
}



# KT: Initialise data frame to store model parameters (repeated code from Louca)
# AIC grid selection and density parameters are currently not used. The latter cannot currently be calculated for sampled ancestor trees.
blank_df = function(){
  return(data.frame(	Ntips								= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            Nnodes								= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            Nevents								= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            mean_event_density					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS), # mean events per time
                            root_age							= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_slope_lambda					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_slope_mu						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_slope_psi						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_slope_Reff						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_slope_removal_rate				= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_slope_sampling_proportion		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_slope_net_growth_rate			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # true_slope_branching_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # true_slope_deterministic_branching_density = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # true_slope_deterministic_sampling_density = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # true_slope_sampling_density			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_mean_lambda					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_mean_mu						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_mean_psi						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_mean_Reff						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_mean_removal_rate				= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_mean_sampling_proportion		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            true_mean_net_growth_rate			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # true_mean_branching_density			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # true_mean_deterministic_branching_density = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # true_mean_deterministic_sampling_density = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # true_mean_sampling_density			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_Ngrid						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_loglikelihood				= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_AIC							= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_Niterations					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_Nevaluations				= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_lambda_R2					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_mu_R2						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_psi_R2						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_Reff_R2						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_removal_rate_R2				= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_sampling_proportion_R2		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_net_growth_rate_R2			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_total_diversity_R2			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_nLTT_R2						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_branching_density_R2		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_sampling_density_R2			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_deterministic_branching_density_R2		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_deterministic_sampling_density_R2			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_lambda_MMNE					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_mu_MMNE						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_psi_MMNE					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_Reff_MMNE					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_removal_rate_MMNE			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_sampling_proportion_MMNE	= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_total_diversity_MMNE		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_net_growth_rate_MMNE		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_nLTT_MMNE					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_branching_density_MMNE		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_deterministic_branching_density_MMNE = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_deterministic_sampling_density_MMNE = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_sampling_density_MMNE		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_slope_lambda				= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_slope_mu					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_slope_psi					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_slope_Reff					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_slope_removal_rate			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_slope_sampling_proportion	= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_slope_net_growth_rate		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_slope_branching_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_slope_sampling_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_slope_deterministic_branching_density = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_slope_deterministic_sampling_density = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_mean_lambda					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_mean_mu						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_mean_psi					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_mean_Reff					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_mean_removal_rate			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_mean_sampling_proportion	= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_mean_net_growth_rate		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_mean_branching_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_mean_deterministic_branching_density = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_mean_deterministic_sampling_density = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # skyline_mean_sampling_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            skyline_PedgeKS						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS), # P-value of Kolmogorov-Smirnov test of edge lengths
                            skyline_PtipKS						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS), # P-value of Kolmogorov-Smirnov test of tip ages
                            skyline_PnodeKS						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS), # P-value of Kolmogorov-Smirnov test of node ages
                            # Parameters not needed for one grid size
                            plinear_Ngrid						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_loglikelihood				= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_AIC							= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_Niterations					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_Nevaluations				= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_lambda_R2					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_mu_R2						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_psi_R2						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_Reff_R2						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_removal_rate_R2				= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_sampling_proportion_R2		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_net_growth_rate_R2			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_nLTT_R2						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_branching_density_R2		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_deterministic_branching_density_R2 = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_deterministic_sampling_density_R2 = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_sampling_density_R2			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_total_diversity_R2			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_lambda_MMNE					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_mu_MMNE						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_psi_MMNE					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_Reff_MMNE					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_removal_rate_MMNE			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_sampling_proportion_MMNE	= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_net_growth_rate_MMNE		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_total_diversity_MMNE		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_nLTT_MMNE					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_branching_density_MMNE		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_deterministic_branching_density_MMNE = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_deterministic_sampling_density_MMNE = rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_sampling_density_MMNE		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_slope_lambda				= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_slope_mu					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_slope_psi					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_slope_Reff					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_slope_removal_rate			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_slope_sampling_proportion	= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_slope_net_growth_rate		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_slope_branching_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_slope_deterministic_branching_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_slope_deterministic_sampling_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_slope_sampling_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_mean_lambda					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_mean_mu						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_mean_psi					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_mean_Reff					= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_mean_removal_rate			= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_mean_sampling_proportion	= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_mean_net_growth_rate		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_mean_branching_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_mean_deterministic_branching_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_mean_deterministic_sampling_density		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            # plinear_mean_sampling_densisty		= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_PedgeKS						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_PtipKS						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS),
                            plinear_PnodeKS						= rep(NA, times=ENSEMBLE_HBD_FITTING_NSIMS)))
}

# KT: Save properties kept consistent for all models fitted (repeated code)
set_consistent_df = function(results_df){
  	results_df$Ntips[sim] 							= Ntips
		results_df$Nnodes[sim]							= Nnodes
		results_df$Nevents[sim]						= Ntips+Nnodes
		results_df$mean_event_density[sim]				= (Ntips+Nnodes)/root_age
		results_df$root_age[sim]						= root_age
    return(results_df)
}

###################################
# MAIN SCRIPT BODY
ENSEMBLE_HBD_SCENARIOS=list(
  list(	name					= "OU",
        type 					= "OU", # possible options are 'OU' and 'exp'
        include					= FALSE,
        time_units				= "year",
        lambda_relaxation_rate	= function(){ runif(n=1, min=0.05, max=0.2) }, 	# random number generator for the OU relaxation rate of lambda
        lambda_stationary_mean	= function(){ runif(n=1, min=1, max=10) }, 		# random number generator for the OU stationary expectation of lambda
        lambda_stationary_rstd	= function(){ return(0.5) }, 					# random number generator for the OU relative std of lambda
        mu_relaxation_rate		= function(){ runif(n=1, min=0.05, max=0.2) }, 	# random number generator for the OU relaxation rate of mu
        epsilon_stationary_mean	= function(){ runif(n=1, min=0.1, max=1) }, 	# random number generator for epsilon:=mu_stationary_mean/lambda_stationary_mean
        mu_stationary_rstd		= function(){ return(0.5) }, 					# random number generator for the OU relative std of mu
        psi_relaxation_rate		= function(){ runif(n=1, min=0.05, max=0.2) }, 	# random number generator for the OU relaxation rate of psi
        psi_stationary_mean		= function(){ exp(runif(n=1, min=log(0.01), max=log(1))) },	# random number generator for the OU stationary expectation of psi
        psi_stationary_rstd		= function(){ return(0.5) }, 					# random number generator for the OU relative std of psi
        max_time				= 10, # duration of a simulation, in years
        random_seed				= 1010,
        fitting_grid_size		= 21, # range of grid sizes to consider when fitting skyline or piecewise linear models. The optimal grid size will be determined via AIC.
        fit_skyline				= TRUE,  # whether to fit skyline (piecewise constant) models to the simulated trees
        fit_plinear				= TRUE), # , # whether to fit piecewise linear models to the simulated trees
  list(	name					= "exp",
        type					= "exp",
        include					= TRUE,
        time_units				= "year",
        lambda_relaxation_rate	= function(){ runif(n=1, min=0.1, max=0.5) }, # random number generator for the exponential rate of lambda
        lambda_start			= function(){ runif(n=1, min=1, max=10) },
        lambda_end				= function(){ runif(n=1, min=1, max=10) },
        mu_relaxation_rate		= function(){ runif(n=1, min=0.1, max=0.5) },
        epsilon_start			= function(){ runif(n=1, min=0.1, max=1) },
        epsilon_end				= function(){ runif(n=1, min=0.1, max=1) },
        psi_relaxation_rate		= function(){ runif(n=1, min=0.1, max=0.5) },
        psi_start				= function(){ exp(runif(n=1, min=log(0.01), max=log(1))) },
        psi_end					= function(){ exp(runif(n=1, min=log(0.01), max=log(1))) },
        max_time				= 10,
        random_seed				= 1010,
        fitting_grid_size		= 21, # (grid size - 1) should be factor of age_grid size. i.e. grid size of 11 for age_grid size of 1000
        fit_skyline				= FALSE,
        fit_plinear				= TRUE)
)


for(e in seq_len(length(ENSEMBLE_HBD_SCENARIOS))){
  scenario = ENSEMBLE_HBD_SCENARIOS[[e]]
  if(!scenario$include){
    cat2(sprintf("Note: Skipping ensemble simulations under scenario '%s', as requested\n",scenario$name))
    next
  }
  cat2(sprintf("Simulating trees under scenario '%s' (of type '%s')..\n",scenario$name,scenario$type))
  scenario_dir = sprintf("%s/ensemble_HBDS_fitting/%s",output_dir,scenario$name)
  
  # reset random seed for this scenario
  if(is.null(scenario$random_seed)) scenario$random_seed = sample.int(n=1000000,size=1)
  cat2(sprintf("  Note: Scenario random seed = %d\n",scenario$random_seed))
  set.seed(scenario$random_seed)
  
  # Prepare data frame
  fit_results = blank_df()
  for(sim in seq_len(ENSEMBLE_HBD_FITTING_NSIMS)){
    cat2(sprintf("  Simulation %d (%s)..\n",sim,scenario$name))
    sim_dir=sprintf("%s/individual_simulations/sim_%d",scenario_dir,sim)
    Nsim_attempts = 0
    
    while(TRUE){
      # Which retention probabilities to use, in addition to zero.
      kappas = c(1,0.5)
      Nsim_attempts = Nsim_attempts + 1
      # simulate random lambda & mu & psi time series
      if(scenario$type=="OU"){
        # determine parameters for the OU processes generating lambda, mu, psi
        lambda_relaxation_rate	= scenario$lambda_relaxation_rate()
        lambda_stationary_mean 	= scenario$lambda_stationary_mean()
        lambda_stationary_rstd 	= scenario$lambda_stationary_rstd()
        psi_relaxation_rate		= scenario$psi_relaxation_rate()
        psi_stationary_mean 	= scenario$psi_stationary_mean()
        psi_stationary_rstd 	= scenario$psi_stationary_rstd()
        mu_relaxation_rate		= scenario$mu_relaxation_rate()
        mu_stationary_mean 		= lambda_stationary_mean * scenario$epsilon_stationary_mean()
        mu_stationary_rstd 		= scenario$mu_stationary_rstd()
        series_times  			= seq(from=0,to=scenario$max_time,length.out=20)
        # generate profiles of lambda, mu, psi
        series_lambda = generate_OU_time_series(	times 			= series_times,
                                                          start_value		= lambda_stationary_mean,
                                                          stationary_mean = lambda_stationary_mean,
                                                          stationary_std	= lambda_stationary_rstd*lambda_stationary_mean,
                                                          decay_rate		= lambda_relaxation_rate,
                                                          constrain_min	= 0.1*lambda_stationary_mean)$values
        series_mu = generate_OU_time_series(	times 			= series_times,
                                                      start_value		= mu_stationary_mean,
                                                      stationary_mean = mu_stationary_mean,
                                                      stationary_std	= mu_stationary_rstd*mu_stationary_mean,
                                                      decay_rate		= mu_relaxation_rate,
                                                      constrain_min	= 0.1*mu_stationary_mean)$values
        series_psi = generate_OU_time_series(	times 			= series_times,
                                                       start_value		= psi_stationary_mean,
                                                       stationary_mean = psi_stationary_mean,
                                                       stationary_std	= psi_stationary_rstd*psi_stationary_mean,
                                                       decay_rate		= psi_relaxation_rate,
                                                       constrain_min	= 0.1*psi_stationary_mean)$values
      }else if(scenario$type=="exp"){
        # determine parameters for the functional forms of lambda, mu, psi, of the type y(t) = A + B*exp(C*t)
        time_start 		= 0
        time_end		= scenario$max_time
        lambda_start	= scenario$lambda_start()
        lambda_end		= scenario$lambda_end()
        lambdaC			= (1-2*rbinom(n=1,size=1,prob=0.5)) * scenario$lambda_relaxation_rate()
        lambdaB 		= (lambda_start-lambda_end)/(exp(lambdaC*time_start) - exp(lambdaC*time_end))
        lambdaA			= lambda_start - lambdaB*exp(lambdaC*time_start)
        psi_start		= scenario$psi_start()
        psi_end			= scenario$psi_end()
        psiC			= (1-2*rbinom(n=1,size=1,prob=0.5)) * scenario$psi_relaxation_rate()
        psiB 			= (psi_start-psi_end)/(exp(psiC*time_start) - exp(psiC*time_end))
        psiA			= psi_start - psiB*exp(psiC*time_start)
        mu_start		= lambda_start * scenario$epsilon_start()
        mu_end			= lambda_end * scenario$epsilon_end()
        muC				= (1-2*rbinom(n=1,size=1,prob=0.5)) * scenario$mu_relaxation_rate()
        muB 			= (mu_start-mu_end)/(exp(muC*time_start) - exp(muC*time_end))
        muA				= mu_start - muB*exp(muC*time_start)
        # generate profiles of lambda, mu, psi
        series_times	= seq(from=0,to=scenario$max_time,length.out=100)
        series_lambda	= lambdaA + lambdaB * exp(lambdaC * series_times)
        series_mu		= muA + muB * exp(muC * series_times)
        series_psi		= psiA + psiB * exp(psiC * series_times)
      }
       # generate random tree based on the specific lambda & mu & psi
      tree_gen = generate_tree_hbds_man(max_time= scenario$max_time,
                                        max_tips		= ENSEMBLE_HBD_FITTING_MAX_NTIPS,
                                        include_extant			= INCLUDE_EXS,
                                        include_extinct			= INCLUDE_EXS,
                                        time_grid				= series_times,
                                        lambda					= series_lambda,
                                        mu						= series_mu,
                                        psi						= series_psi,
                                        kappa					= 0,
                                        splines_degree			= 1,
                                        no_full_extinction		= FALSE,
                                        tip_basename			= "tip.")
			if(!tree_gen$success) next
			tree	 = tree_gen$tree
			root_age = get_tree_span(tree)$max_distance
			Ntips	 = length(tree$tip.label)
			if((Ntips<ENSEMBLE_HBD_FITTING_MIN_NTIPS) || (Ntips>ENSEMBLE_HBD_FITTING_MAX_NTIPS)) next
			
			# calculate some basic properties of this tree
			# Note that age is counted backward, with 0 being at the last sampled tip
			root_age 		= castor::get_tree_span(tree)$max_distance
			stem_age		= tree_gen$root_time + root_age # age of the stem, i.e. when the process actually started
			end_age			= root_age - (tree_gen$final_time-tree_gen$root_time) # age at which the HBDS simulation was halted. This might be slightly negative, e.g. if the process halted after the last sampled tip.
			tree_LTT 	 	= castor::count_lineages_through_time(tree, Ntimes=200, include_slopes=TRUE)
			tree_LTT$ages 	= root_age - tree_LTT$times
			age0 			= tree_LTT$ages[which.max(tree_LTT$lineages)]
			tree_LTT0 		= approx(x=tree_LTT$ages,y=tree_LTT$lineages,xout=age0)$y	

      properties = list(root_age, stem_age, end_age, tree_LTT, age0, tree_LTT0)

			# simulate deterministic model
			age_grid = seq(from=0, to=stem_age, length.out=1000)
			sim_true = simulate_deterministic_hbds(	age_grid		= age_grid,
													lambda			= rev(approx(x=series_times,y=series_lambda,xout=tree_gen$final_time+end_age-age_grid)$y),
													mu				= rev(approx(x=series_times,y=series_mu,xout=tree_gen$final_time+end_age-age_grid)$y),
													psi				= rev(approx(x=series_times,y=series_psi,xout=tree_gen$final_time+end_age-age_grid)$y),
													kappa			= 0,
													requested_ages	= seq(from=0,to=root_age,length.out=1000),
													age0			= age0,
													LTT0			= tree_LTT0,
													splines_degree	= 1)
			if(!sim_true$success) next
			# all seems OK with this simulation
			break
		}
		Nnodes = tree$Nnode
		cat2(sprintf("    Note: Tree has %d tips, %d nodes, spans %g, Nsim_attempts=%d\n",Ntips,Nnodes,root_age,Nsim_attempts))
		dir.create(sim_dir, showWarnings = FALSE, recursive=TRUE)
		castor::write_tree(tree, file=sprintf("%s/tree.tre",sim_dir))        
          

		# save some basic stats about this simulated tree
    fit_results = set_consistent_df(fit_results)
    fit_results = set_true_results(fit_results, sim_true)
    fit_results_fixed_psi = fit_results

    plot_model(	model_name		= sprintf("%s.sim_%d_kappa_0_c1",scenario$name,sim),
                sim				= sim_true,
                plot_maxx		= NULL,
                plot_basepath	= sprintf("%s/deterministic_simulation_plots/",sim_dir),
                time_units		= scenario$time_units,
                verbose			= TRUE,
                verbose_prefix	= "      ")			
          
    # save deterministic curves of this model
    fout = prepare_output_file(file_path = sprintf("%s/deterministic_simulation.tsv",sim_dir), FALSE, verbose=FALSE, verbose_prefix="  ")
    cat(sprintf("# Deterministic simulation %d of scenario '%s'\n# Generated on: %s\n# Random seed for this model: %d\n",sim,scenario$name,display_date_time,scenario$random_seed), file=fout, append=FALSE)
    param_names = c("LTT","nLTT","Pmissing","lambda","mu","psi","PDR","IPDR","PSR","Reff","removal_rate","sampling_proportion","diversification_rate", "lambda_psi") #"branching_density", "deterministic_branching_density",
    cat(sprintf("age\t%s\n",paste(param_names,collapse="\t")), file=fout, append=TRUE)
    write.table(cbind(sim_true$ages,as.data.frame(do.call(cbind, sim_true[param_names]))), file=fout, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    close(fout)

    present_day_psi = approx(x=sim_true$ages,y=sim_true$psi,xout=0)$y
      # Fit congruent models 
      results = plinear_fit_and_plot(sim_true,  tree, properties, correct_psi=FALSE, 0, fit_results)
      if(!(results[[1]])){
        cat2(sprintf("      ERROR: Fitting failed: %s\n",results[[2]]));
				if(ENSEMBLE_HBD_FITTING_REPEAT_FAILED_TREES){
					cat2(sprintf("        Repeating entire simulation %d and fitting\n",sim))
					unlink(sim_dir, recursive=TRUE)
					next
				}
      }
        fit_results = as.data.frame(results[[2]])
        fit = results[[3]]
    # }
      # Also fit congruent model with correct sampling rate to see if the birth and death rates are also correct.    
      results = plinear_fit_and_plot(sim_true,  tree, properties, correct_psi=TRUE, 0, fit_results_fixed_psi)
      if(!(results[[1]])){
        cat2(sprintf("      ERROR: Fitting failed: %s\n",results[[2]]));
				if(ENSEMBLE_HBD_FITTING_REPEAT_FAILED_TREES){
					cat2(sprintf("        Repeating entire simulation %d and fitting\n",sim))
					unlink(sim_dir, recursive=TRUE)
					next
				}
      }
        fit_results_fixed_psi = as.data.frame(results[[2]])
        fit = results[[3]]
    #}
    fit_results = as.data.frame(fit_results)
    fit_results_fixed_psi = as.data.frame(fit_results_fixed_psi)
    fit_dir = results[[4]]

    first_fit = simulate_deterministic_hbds(	age_grid		= fit$age_grid, 
                                                   lambda			= fit$param_fitted$lambda,
                                                   mu				= fit$param_fitted$mu,
                                                   psi				= fit$param_fitted$psi,
                                                   kappa			= 0,
                                                   splines_degree	= 1,
                                                   requested_ages	= sim_true$ages,
                                                   age0			= age0,
                                                   LTT0			= tree_LTT0)

    plot_model(	model_name		= sprintf("%s.sim_%d_kappa_0_c2",scenario$name,sim),
                sim				= sim_true,
                plot_maxx		= NULL,
                plot_basepath	= sprintf("%s/deterministic_simulation_plots/",sim_dir),
                time_units		= scenario$time_units,
                verbose			= TRUE,
                verbose_prefix	= "      ")			
          
    # save deterministic curves of this model
    fout = prepare_output_file(file_path = sprintf("%s/deterministic_simulation.tsv",sim_dir), FALSE, verbose=FALSE, verbose_prefix="  ")
    cat(sprintf("# Deterministic simulation %d of scenario '%s'\n# Generated on: %s\n# Random seed for this model: %d\n",sim,scenario$name,display_date_time,scenario$random_seed), file=fout, append=FALSE)
    param_names = c("LTT","nLTT","Pmissing","lambda","mu","psi","PDR","IPDR","PSR","Reff","removal_rate","sampling_proportion","diversification_rate", "lambda_psi") #"branching_density", "deterministic_branching_density",
    cat(sprintf("age\t%s\n",paste(param_names,collapse="\t")), file=fout, append=TRUE)
    write.table(cbind(first_fit$ages,as.data.frame(do.call(cbind, first_fit[param_names]))), file=fout, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    close(fout)

        
    cat2(sprintf("Saving results from all simulations of scenario '%s' kappa %s..\n",scenario$name, 0))
    fout = prepare_output_file(file_path = sprintf("%s/all_simulation_results.tsv",scenario_dir), FALSE, verbose=FALSE, verbose_prefix="  ")
    cat(sprintf("# Summary results from all simulations of scenario '%s' %s \n# Generated on: %s\n# Random seed for this scenario: %d\n",scenario$name, 0, display_date_time,scenario$random_seed), file=fout, append=FALSE)
    cat(sprintf("%s\t%s\n",paste(colnames(fit_results),collapse="\t"), paste("fixed_psi",colnames(fit_results),collapse="\t")), file=fout, append=TRUE)
    write.table(cbind(fit_results, fit_results_fixed_psi), file=fout, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    close(fout)

    # For each non-zero retention probability we want to consider, simulate models using the rates from the kappa = 0 case. Plot relevant curves for each model.
    for (k in 1:length(kappas)){
      while(TRUE){
        kappa = kappas[k]
        cat("kappa ", kappa, "\n")
        fit_results = blank_df()
        fit_results = set_consistent_df(fit_results)
        cat("simulating deterministic\n")
        # Simulate deterministic model using original rates but current retention probability.
        sim_true = simulate_deterministic_hbds(	age_grid		= age_grid,
										lambda			= rev(approx(x=series_times,y=series_lambda,xout=tree_gen$final_time+end_age-age_grid)$y),
										mu				= rev(approx(x=series_times,y=series_mu,xout=tree_gen$final_time+end_age-age_grid)$y),
										psi				= rev(approx(x=series_times,y=series_psi,xout=tree_gen$final_time+end_age-age_grid)$y),
										kappa			= kappa,
									  requested_ages	= seq(from=0,to=root_age,length.out=1000),
										age0			= age0,
										LTT0			= tree_LTT0,
										splines_degree	= 1)                    
        str(sim_true)
        if(!sim_true$success) next
        # all seems OK with this simulation
        
        fit_results = set_true_results(fit_results, sim_true)
        cat("call sim_plinear\n")
        # Simulate deterministic model using rates under second scenario from kappa = 0 case, but now use current retention probability.
        sim_result = sim_plinear(fit,sim_true,age0,tree_LTT0,kappa, fit_results, fit_dir)
        if(!(sim_result[[1]])){
          next
        }
        break
      }
      fit_results = sim_result[[2]]
      sim_fit = sim_result[[3]]
      
      plot_model(	model_name		= sprintf("%s.sim_%d_kappa_%s_c1",scenario$name,sim,formatC(kappa, digits = 1, format = "f")),
                sim				= sim_true,
                plot_maxx		= NULL,
                plot_basepath	= sprintf("%s/deterministic_simulation_plots/",sim_dir),
                time_units		= scenario$time_units,
                verbose			= TRUE,
                verbose_prefix	= "      ")			
      plot_model(	model_name		= sprintf("%s.sim_%d_kappa_%s_c2",scenario$name,sim,formatC(kappa, digits = 1, format = "f")),
                sim				= sim_fit,
                plot_maxx		= NULL,
                plot_basepath	= sprintf("%s/deterministic_simulation_plots/",sim_dir),
                time_units		= scenario$time_units,
                verbose			= TRUE,
                verbose_prefix	= "      ")

    cat2(sprintf("Saving results from all simulations of scenario '%s' kappa %s..\n",scenario$name, formatC(kappa, digits = 1, format = "f")))
    fout = prepare_output_file(file_path = sprintf("%s/all_simulation_results.tsv",scenario_dir), FALSE, verbose=FALSE, verbose_prefix="  ")
    cat(sprintf("# Summary results from all simulations of scenario '%s' %s \n# Generated on: %s\n# Random seed for this scenario: %d\n",scenario$name, 0, display_date_time,scenario$random_seed), file=fout, append=FALSE)
    cat(sprintf("%s\t%s\n",paste(colnames(fit_results),collapse="\t"), paste("fixed_psi",colnames(fit_results),collapse="\t")), file=fout, append=TRUE)
    write.table(fit_results, file=fout, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    close(fout)
			
    }
  }
  break
}
cat2(sprintf("Done. All outputs were written to '%s'\n",output_dir));
