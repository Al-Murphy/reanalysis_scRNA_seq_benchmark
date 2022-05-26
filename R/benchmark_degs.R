#-------------------------------------
# Note this script will rerun the benchmarking analysis using the updated 
# version of hierarchicell which returns ROCs so the methods' sensitivity can be
# compared at the same type 1 error rates.
# 
# The seed is set here to give the exact same results as those produced
# in the analysis.
#
# NOTE this is the benchmarking based on an unequal number of differentially and
# non-differentially expressed genes.
#
# WARNING: This is both extremely memory and time intensive. Each number
# of DEG combination is rerun 50 times giving a total of 200 runs. For the 
# mixed model/two-part hurdle models each of these runs can take a very long 
# time and use a lot of memory.
#-------------------------------------

devtools::install_github("neurogenomics/hierarchicell")

#necessary steps:
#not passing in data uses pancreatic alpha cells
gene_counts_filtered <- hierarchicell::filter_counts()
data_summ <-
  hierarchicell::compute_data_summaries(gene_counts_filtered, type = "Norm")

#methods:
method <- c("MAST","MAST_RE","MAST_Combat","GLM_tweedie","GLMM_tweedie",
            "GEE1","Pseudobulk_mean","Pseudobulk_sum","Monocle","ROTS")
#done for 5,000 non DEGs and differing DEG proportions 50 times
set.seed(101)
seeds <- sample.int(10000, 50, replace = FALSE)
genes <- 5000
#done for 100 cells in 20 individuals
n_ind <- 20
n_cells <- 100

#genes,prop_diff_genes
genes <- 5000
prop_DEGs <- c(0.05,0.1,0.2,0.3)
#need to loop by n_ind,n_cells,seeds
#get runs
l <- list("seeds"=seeds,"n_ind"=n_ind,"n_cells"=n_cells,"prop_DEGs"=prop_DEGs)
runs <- expand.grid(l)
runs$id <- paste0(runs$seeds,":",runs$n_ind,":",runs$n_cells,":",runs$prop_DEGs)
runs$save_id <- paste0(runs$n_ind,"_",runs$n_cells,"_",genes,"_",gsub(".","",runs$prop_DEGs,fixed=TRUE),"_",runs$seeds,".csv")
#in case run stalled, remove the runs from the list that have already completed
#5_100_8507.csv
#filter out those already run
completed_runs <- list.files("./results_imbal_DEGs/", pattern=NULL, all.files=FALSE,full.names=FALSE)
runs <- runs[!(runs$save_id %in% completed_runs),]
print(paste0(nrow(runs)," runs to complete."))

#do it in parallel
#do it in parallel
suppressMessages(library(parallel))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))
suppressMessages(library(data.table))
suppressMessages(library(hierarchicell))
#TODO: Add number of cores on your machine to run in parallel
#Note the largest runs need up to 120gb RAM each 
number_cores <- #add here

#make temp dir to store each runs results - note there are also automatic 
#checkpoints in hierarchicell which will also store in this folder
dir.create("./temp", showWarnings = FALSE)

doParallel::registerDoParallel(min(number_cores,nrow(runs)))
run_return <-
  foreach::foreach(run_i = seq_len(nrow(runs)),
                   .final = function(run_i) setNames(run_i, runs$id),
                   .packages = c("data.table","hierarchicell")) %dopar% {
                      seed_i <- runs[run_i,]$seeds
                      n_ind_i <- runs[run_i,]$n_ind
                      n_cells_case_i <- runs[run_i,]$n_cell_cases
                      n_cells_control_i <- runs[run_i,]$n_cell_controls
                      print("#############################")
                      print(paste0(seed_i,":",n_ind_i,":",n_cells_case_i,":",
                                   n_cells_control_i))
                      print("############################")
                      #run each DE method
                      res_tbl <-
                        data.table("seed"=rep(seed_i,length(method)),
                                   "n_ind"=rep(n_ind_i,length(method)),
                                   "n_cells"=rep(n_cells_i,length(method)),
                                   "method"=method,
                                   "type_1_err"=0.00,
                                   "type_2_err"=0.00,
                                   "MCC"=0.00,
                                   "TP"=0.00,
                                   "FP"=0.00,
                                   "FN"=0.00,
                                   "TN"=0.00
                        )
                      res_raw <- vector(mode="list",length=length(method))
                      names(res_raw)<-method
                      for (i in seq_len(length(method))){
                        print("-----------------------")
                        print(method[[i]])
                        print("-----------------------")
                        partial_runs <- 
                          list.files("./temp/", pattern=NULL, all.files=FALSE,
                                     full.names=FALSE)
                        method_i <- method[[i]]
                        partial_i <- 
                          paste0(n_ind_i,"_",n_cells_i,"_",genes,"_",
                                 gsub(".", "",prop_diff_genes,fixed=TRUE),"_",
                                 seed_i,"_",method_i,".rds")
                        #if method completed before don't rerun
                        if(partial_i %in% partial_runs){
                          res <- readRDS(paste0("./temp/",partial_i))
                        }
                        else{
                          #increase deg_fc_interval to 1 given low nums DEGs
                          res <-
                            error_hierarchicell(data_summ,method = method_i,
                                                n_genes = genes,
                                                prop_diff_genes = 
                                                  prop_diff_genes,
                                                n_per_group = n_ind_i,
                                                cells_per_case = n_cells_i,
                                                cells_per_control = n_cells_i,
                                                ncells_variation_type="Poisson",
                                                pval = 0.05,
                                                fc_range_min = 1.1,
                                                fc_range_max = 10,
                                                deg_fc_interval = 1,
                                                seed = seed_i,
                                                checkpoint = "./temp")
                          saveRDS(res,paste0("./temp/",partial_i))
                        }
                        res_tbl[method==method[[i]],MCC:= res[["MCC"]]]
                        res_tbl[method==method[[i]],
                                type_1_err:= res[["type_1_err"]]]
                        res_tbl[method==method[[i]],
                                type_2_err:= res[["type_2_err"]]]
                        res_tbl[method==method[[i]],
                                TP:= res[["TP"]]]
                        res_tbl[method==method[[i]],
                                TN:= res[["TN"]]]
                        res_tbl[method==method[[i]],
                                FP:= res[["FP"]]]
                        res_tbl[method==method[[i]],
                                FN:= res[["FN"]]]
                        res_raw[[method[i]]] <- data.table(pred=res$pred,
                                                           actual=res$actual)
                        
                      }
                      #write each out in case job fails
                      fwrite(res_tbl,paste0("./results_imbal_DEGs/",n_ind_i,"_",
                                            n_cells_i,"_",genes,"_",
                                            gsub(".", "",prop_diff_genes,
                                                 fixed=TRUE),"_",seed_i,".csv"))
                      saveRDS(res_raw,paste0("./results_imbal_DEGs/",n_ind_i,
                                             "_",n_cells_i,"_",genes,"_",
                                             gsub(".", "",prop_diff_genes,
                                                  fixed=TRUE),"_",seed_i,".rds")
                              )
                      return(list(res_tbl,res_raw))
                   }
#we will use res_raw data to produce the ROCs
#load and combine into single csv
method <- c("MAST","MAST_RE","MAST_Combat","GLM_tweedie","GLMM_tweedie",
            "GEE1","Pseudobulk_mean","Pseudobulk_sum","Monocle","ROTS")

#add descriptive method name
method_desc <- c("Two-part hurdle: Default","Two-part hurdle RE",
                 "Two-part hurdle: Corrected","Tweedie: GLM","Tweedie: GLMM",
                 "GEE1","Pseudobulk: Mean","Pseudobulk: Sum","Tobit",
                 "Modified t")
#filter to just rds files which contain raw data for ROCs
completed_runs <- list.files("./results_imbal_DEGs/", pattern="\\.rds$",
                             all.files=FALSE,full.names=TRUE)
prop_DEGs <- c(0.05,0.1,0.2,0.3)
#store data for ROCs
dat_rocs <- vector(mode="list",length=length(prop_DEGs))
names(dat_rocs) <- prop_DEGs
for(run_i in seq_along(completed_runs)){
  #get proportion of DEGs from file name
  prop_deg <- basename(completed_runs[[run_i]])
  #remove bits after
  prop_deg <- sub("^([^_]*_[^_]*_[^_]*_[^_]*).*", "\\1", prop_deg)
  #remove bits before
  prop_deg <- sub("(.*_\\s*(.*$))", "\\2", prop_deg)
  #add back in decimal place and make integer
  prop_deg <- as.numeric(sub("(.{1})(.*)", "\\1.\\2", prop_deg))
  run_i_dat <- readRDS(completed_runs[[run_i]])
  #make single data table
  run_i_dat <- data.table::rbindlist(run_i_dat,idcol=TRUE)
  setnames(run_i_dat,".id","method")
  #going to combine diff seed runs into one ROC
  if(is.null(dat_rocs[[as.character(prop_deg)]])){
    dat_rocs[[as.character(prop_deg)]] <- run_i_dat
  }else{
    dat_rocs[[as.character(prop_deg)]] <- 
      data.table::rbindlist(list(dat_rocs[[as.character(prop_deg)]],run_i_dat))
  }
  
}

#combine everything into one datatable
dat_rocs <- data.table::rbindlist(dat_rocs,idcol=TRUE)
setnames(dat_rocs,".id","prop_DEGs")
#since predictions are p-values need to go 1 - pred
dat_rocs[,pred:=1-pred]
#now save results
fwrite(dat_rocs,"./results/method_performances_imbal_DEGs.csv")
