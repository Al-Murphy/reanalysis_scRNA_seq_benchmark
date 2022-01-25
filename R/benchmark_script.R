#-------------------------------------
# Note this script will rerun the benchmarking analysis using the updated 
# version of hierarchicell which returns a balanced performance measure:
# Mathews Correlation Coefficient (MCC) as well as the type 1 error.
# The seed is set here to give the exact same results as those produced
# in the analysis.
#
# WARNING: This is both extremely memory and time intensive. Each number
# of cells and number of individuals combination is rerun 50 times giving
# a total of 20,000 runs. For the mixed model/two-part hurdle models each
# of these runs can take a very long time and use a lot of memory. For 
# example, 1 run of GLMM Tweedie with 40 individuals and 500 cells will 
# take up to 120 gb RAM and will run for over a week
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
#set the seed for all runs which all use 5000 genes 
#(5000 non and 5000 differentially expressed genes)
set.seed(101)
seeds <- sample.int(10000, 50, replace = FALSE)
genes <- 5000
#done for increasing numebrs of individuals and cells
n_ind <- c(5,10,20,30,40)
n_cells <- c(50,100,250,500)

#need to loop by n_ind,n_cells,seeds
#get runs
l <- list("seeds"=seeds,"n_ind"=n_ind,"n_cells"=n_cells)
runs <- expand.grid(l)
runs$id <- paste0(runs$seeds,":",runs$n_ind,":",runs$n_cells)
runs$save_id <- paste0(runs$n_ind,"_",runs$n_cells,"_",runs$seeds,".csv")
#filter out those runs that already completed in case stalled
completed_runs <- 
  list.files("./results/", pattern=NULL, all.files=FALSE,full.names=FALSE)
runs <- runs[!(runs$save_id %in% completed_runs),]
print(paste0(nrow(runs)," runs to complete."))

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
                     n_cells_i <- runs[run_i,]$n_cells
                     print("#############################")
                     print(paste0(seed_i,":",n_ind_i,":",n_cells_i))
                     print("############################")
                     #run each DE method
                     res_tbl <-
                       data.table("seed"=rep(seed_i,length(method)),
                                  "n_ind"=rep(n_ind_i,length(method)),
                                  "n_cells"=rep(n_cells_i,length(method)),
                                  "method"=method,
                                  "type_1_err"=0.00,
                                  "MCC"=0.00)
                     for (i in seq_len(length(method))){
                       print("-----------------------")
                       print(method[[i]])
                       print("-----------------------")
                       partial_runs <- 
                         list.files("./temp/", pattern=NULL, all.files=FALSE,
                                      full.names=FALSE)
                       method_i <- method[[i]]
                       partial_i <- paste0(n_ind_i,"_",n_cells_i,"_",
                                            seed_i,"_",method_i,".rds")
                       #if method completed before don't rerun
                       if(partial_i %in% partial_runs){
                         res <- readRDS(paste0("./temp/",partial_i))
                       }
                       else{
                         res <-
                           error_hierarchicell(data_summ,method = method_i,
                                               n_genes = genes,
                                               n_per_group = n_ind_i,
                                               cells_per_case = n_cells_i,
                                               cells_per_control = n_cells_i,
                                               ncells_variation_type="Poisson",
                                               pval = 0.05,
                                               fc_range_min = 1.1,
                                               fc_range_max = 10,
                                               deg_fc_interval = 0.5,
                                               seed = seed_i,
                                               #save checkpoints
                                               checkpoint = "./temp/")
                         saveRDS(res,paste0("./temp/",partial_i))
                       }
                       res_tbl[method==method[[i]],MCC:= res[["MCC"]]]
                       res_tbl[method==method[[i]],
                               type_1_err:= res[["type_1_err"]]]
                       
                     }
                     #write each out in case job fails
                     fwrite(res_tbl,paste0("./results/",n_ind_i,"_",n_cells_i,
                                            "_",seed_i,".csv"))
                     return(res_tbl)
                   }

#load res
completed_runs <- 
  list.files("./results/", pattern=NULL, all.files=FALSE,full.names=TRUE)
res<-rbindlist(lapply(completed_runs, function(x) fread(x)))
#save results
fwrite(res,"./results/method_performances.csv")