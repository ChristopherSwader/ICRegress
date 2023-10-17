


#FUNCTION DEFINITIONS####

#' Benchmarking
#'
#' Run a random set of generated (fake) data and comparing the results across multiple methods
#' @param random_seeds Benchmarking parameter. A vector of integers to serve as random seeds.
#' @param generations Benchmarking parameter. A vector of integers specifying the number of generations for initial ICR genetic training.
#' @param solution_thresh Benchmarking parameter. A vector of positive numeric values specifying the discrete solution target on each side of the solution, where solutions may range from 0 to 1.
#' @param n_strands Benchmarking parameter. A vector of integers indicating the number of population that survive each generation of the genetic training algorithm.
#' @param real_models Benchmarking parameter. A vector of integers specifying the 'ground truth' number of underlying models in synthetic data.
#' @param error_sd Benchmarking parameter. A vector of numeric values indicating the standard deviation of error of synthetic data. Synthetic data range from 0 to 1.
#' @param n_IVs Benchmarking parameter. A vector of integers specifying the 'ground truth' number of underlying relevant independent variables for the given dependent variable in synthetic data.
#' @param n_cases Integer value. The number of overall cases in the synthetic data.
#' @param test_set_model_closeness Do not change. Will be removed as a parameter.
#' @param parallelize T indicates that the benchmarking will happen via multiple cores (sockets), whereas F indicates sequential.
#' @param detailed_plots_on Produces bi-dimensional X-Y plots for various steps of the ICR algorithm. Only useful for low-dimensional analysis.
#' @param boxplots Produces boxplots comparing method outcomes with one another at the end of the benchmarking routine.
#' @param starting_population An integer representing the number of starting population.
#' @param hold_back_cores Prevents X cores from being involved in the parallel processing, thus reserving computing power.
#' @keywords Benchmarking Comparison
#' @examples
#'  test <- benchmarking(random_seeds=100001:100025, n_cases=10000)
#' @export
benchmarking <- function( random_seeds=c(100001:100050),#:100050, #parameter for testing
                          generations=c(40,200,1000), #parameter for testing
                          solution_thresh = c( .01 ), #parameter for testing
                          n_strands = c(50), #parameter for testing
                          real_models = c( 20), #parameter for testing
                          error_sd = c(.05), #parameter for testing
                          n_IVs=c(30), #parameter for testing
                          n_cases= 5000, #n rows in dataset
                          mutation_rate=.2,
                          n_children=2,
                          test_set_model_closeness= c("residuals"),
                          parallelize=F,
                          detailed_plots_on=F,
                          boxplots=T,
                          starting_population,
                          hold_back_cores=2,#how many of the available cores are not used for parallel procedure?
                          death_by_ageing =60,
                          nelitism=10){






  n_parameters <- 8
  n_results <- 12
  parameter_testing <- data.frame(matrix(0, nrow = 0, ncol=n_parameters+n_results))
  colnames(parameter_testing) <- c("seed", "generations", "solution_thresh", "n_dna_strands", "real_model_n",
                                   "N_IVs", "error_sd", "test_set_model_closeness",
                                   "simple_r2","init_training_r2", "cond_Training_r2", "stream_training_r2", "cond_test_r2", "stream_test_r2",
                                   "casebased_training_r2", "casebased_testing_r2", "case_assign_score", "model_similarity_error",
                                   "nmodel_condense_streamline", "nmodel_case_based")








 #Parallelize?####

  if (parallelize==T){

    print("starting parallel parameter testing")
    closeAllConnections()
    n_cores <- parallel::detectCores() - hold_back_cores

    my_cluster <- parallel::makeCluster(
      n_cores,
      type = "PSOCK")

    doParallel::registerDoParallel(my_cluster)



    parallel::clusterEvalQ(my_cluster,  library(ICRegress))


  }else { #if sequential
    #SEQUENTIAL?####

    closeAllConnections()
    foreach::registerDoSEQ()
  }
  #number of combinations
  p_number <- length(random_seeds)*
    length(generations)*
    length(solution_thresh)*
    length(n_strands)*
    length(real_models)*
    length(n_IVs)*
    length(error_sd)*
    length(test_set_model_closeness)

  combos <- expand.grid(random_seeds, generations, solution_thresh,
                        n_strands, real_models, n_IVs,
                        error_sd, test_set_model_closeness)

  set.seed(774411)
par_results <-  foreach(i = 1:p_number, .combine = rbind) %dorng% {





    ra <- combos$Var1[i]
    ge <- combos$Var2[i]
    so <- combos$Var3[i]
    ns <- combos$Var4[i]
    re <- combos$Var5[i]
    ni <- combos$Var6[i]
    er <- combos$Var7[i]
    tc <- combos$Var8[i]



    compare_three_methods <- data.frame(matrix(nrow=0, ncol=15))


    parameter_set <- c(ra, ge, so, ns, re, ni, er, tc)
    # saveRDS(parameter_testing, paste0(ra,ge,so,ns,re, ni, er, tc, "parameter testing.RDS"))
    #print is for for debugging
    print("these parameters")
    print(parameter_set)

    #collect output row into parameter_testing table




      set.seed(ra)

    ##1B. Load fake data###


    sample.data <- sample_data(c.number = n_cases, #vs. 10000                    # number of cases
                               m.number = re,                    # number of distinct models
                               m_sizes = c(rep(floor(1000/re), re-1),(1000-sum(rep(floor(1000/re), re-1)))),        # how many cases each model explains
                               bin.v.number = round(ni/2), #binary                # number of binary variables
                               con.v.number = ni-round(ni/2), #continuous                # number of continuous variables
                               error = TRUE,                    # adds error term to the models
                               error.sd = er,   #.1 or .05               # sd of the error term
                               data.form = 'normalized',               # form of the data
                               messy_subgroups = F)


    #note that the output includes a build_list, this can be directly compared with the cut_list in the last step.###
    newdf <- sample.data$data
    original_data <<- newdf #before dummy splits. newdf will be transformed to after dummy splits
    models <- sample.data$models

    dv <<- colnames(newdf)[1]
    iv <<- colnames(newdf)[2:length(colnames(newdf))]
    variables <- c(dv, iv)
    dummy_sets <<- rep(F, length(variables))



    original_variables <- c(dv, iv)
    original_variables <- original_variables[!original_variables=="model"]
    original_variables <<- original_variables


if (T %in% dummy_sets){
  save_row_names <- rownames(newdf)
    df_full_dummies <- dummy_cols(newdf, select_columns = original_variables[dummy_sets], remove_first_dummy = F, ignore_na = T, remove_selected_columns = T)
    df_full_dummies<<- df_full_dummies
rownames(newdf) <- save_row_names
}else{
  df_full_dummies<<- newdf
}
    #need to build "original_variables" and dummy_Sets

    if (detailed_plots_on==T){
      simple <-  ggplot(newdf, aes(x=X1, y=y)) +
        geom_point() +
        geom_smooth(method=lm, se=FALSE)+
        ggtitle("Simple Regression")
      print(simple)
    }

    model_plot(data_with_models=newdf, label = "True Models",
               models = sample.data$models,
               base_line_plot = F,
               turned_on = detailed_plots_on,
               X_variable = 1,
               local_regression_line = F)


    #remove the models from the sample data

    newdf <- newdf[,!colnames(newdf)=="model"]




    ##2. SPLIT DATA INTO TRAINING AND TEST SETS###
    #newdf <- sample.data$data
    training_proportion <- 0.6
    training_set_N <- round(nrow(newdf)*training_proportion)
    training_indices <<- sample(nrow(newdf), training_set_N, replace=F)
    training_set <<- as.data.frame(newdf[training_indices, ])
    testing_set <<- as.data.frame(newdf[-training_indices, ])


    ##3. SPAWN CHILDREN###
    children <- children_spawn(n_children= starting_population,#0 #500 was used in parameter testing,
                               dna_length=ncol(newdf), #this should be the correct number of columns. We need one coefficient for every IV column (including dummies) and one forthe intercept (the DV column)
                               col_names=c('intercept',colnames(newdf)[2:length(colnames(newdf))]), #the DV name column is replaced with 'intercept'
                               normal=T, #normally distributed with the following SD. If not, then this_range is used
                               this_sd=.5, #standard deviation of the coefficients.. this shoudl be tailored to the data's form! .25 is a good default
                               this_range=c(-10, 10))

    ##4.  RUN TRAINING SET###

    cat("Running training set \n")

    solution_threshold <<- so


#Test children call####
    result_training <- test_children(dna_pool=children,
                                     input_data=training_set,
                                     closeness_threshold=solution_threshold, #defines how close a solution is to the found solution
                                     generations=ge,
                                     difficulty=.8,#currently unused
                                     death_rating_thresh=0,
                                     death_number_thresh=ns,
                                     uniqueness_power=8,
                                     death_by_ageing = death_by_ageing,
                                     number_of_children = n_children,
                                     min_strand_solution_proportion=0, #valid solutions must cover at least this proportion of cases.
                                     mutations = mutation_rate,
                                     mating = "random",
                                     mating_power = 2,
                                     selection = "probability",
                                     nelitism = nelitism)

    # print("training results")

    if (dir.exists("testing")){

    }else{
      dir.create("testing")
    }

  #  saveRDS(result_training,paste0("./testing/",i, " training results.RDS"))

    ### IF no training. training is NA###
    if (!is.na(result_training[1])){ #if the dna did not die out, then continue with condensation ONLY



      assign_after_training <- model_closeness(data = training_set, result_training$dna_profiles,
                                               models = result_training$dna_pool ,
                                               method = tc)


      model_plot(data_with_models=assign_after_training, label = "Post training",
                 models = sample.data$models,
                 base_line_plot = F,
                 turned_on = detailed_plots_on,
                 X_variable = 1,
                 local_regression_line = T)





      ##5. CONDENSE Training models, then streamline, and then re-run training set and try test set###
      #this is all linked together into a loop for troubleshooting and optimizing.. to find the best condensation-streamlining combination

      #only condense if more than one model was found!
      ###IF more than one training model found###
      if (nrow(result_training$dna_pool)>1){
        condense_max <- 6
        condense_min <- 2
        parameter_check <- 1:(length(condense_min:condense_max)+1)
        condense_number <- c( condense_min:condense_max, "NULL")
        method_style <- c( rep("best", length(condense_min:condense_max)),"cluster" ) #condensation method
        best_check_ratio <- 0
        minimum_models <- 2
        best_check_list <- list()
        condense_possibilities <- vector("list", length =length(parameter_check) )

        for (i in parameter_check){


          this_condense_number <- condense_number[i]
          this_method_style <- method_style[i]




          this_condense <- condense_models(data = result_training,
                                           training_set = training_set,
                                           method = this_method_style, #or "cluster" / "best"
                                           auto_cluster = TRUE, # indicates if algorithm should find the ideal number of cluster by itself
                                           nmodel = this_condense_number, # only relevant for method = best or when auto_cluster is set to FALSE
                                           solution_table=result_training$solution_table, #relevant only for genetic method
                                           no_correlation_thresh =-.1)

          print("this condense")
          print(this_condense$models)

          model_plot(data_with_models=this_condense$data_with_models,
                     label = paste0("Condense training . Parameter: ",i),
                     turned_on = detailed_plots_on,
                     X_variable=1)

          condensed_training <- predict_new(models=this_condense$models,
                                            profiles=this_condense$dna_profiles,
                                            new_data=this_condense$data_with_models,
                                            closeness_threshold=solution_threshold, #set same as in training set!
                                            assign_models = F)

          ##6. STREAMLINE Take solution table and calculate true regressions for underlying subpopulations###
          #This makes preditions more accurate and more clearly separates subpopulations



          streamlined_data_and_models <- streamline(models=this_condense$models,
                                                    profiles=this_condense$dna_profiles,
                                                    data=training_set)





          model_plot(data_with_models=streamlined_data_and_models$model_assignment,
                     label = paste0("Streamlined Training. Parameter: ",i),
                     turned_on = detailed_plots_on,
                     X_variable = 1,
                     base_line_plot = F,
                     models = sample.data$models)


          #now test how well we can predict training cases with the model

          result_streamlined_training <- predict_new(models=streamlined_data_and_models$dna_pool,
                                                     profiles=streamlined_data_and_models$dna_profiles,
                                                     new_data=streamlined_data_and_models$model_assignment,
                                                     closeness_threshold=solution_threshold, #set same as in training set!
                                                     assign_models = F)







          ##7.  RUN TEST SET###


          result_testing <- predict_new(models=streamlined_data_and_models$dna_pool,
                                        profiles=streamlined_data_and_models$dna_profiles,
                                        new_data=testing_set,
                                        closeness_threshold=solution_threshold, #set same as in training set!
                                        assign_models = T,
                                        method=tc)

          assign_after_testing <- model_closeness(data = testing_set,
                                                  profiles=streamlined_data_and_models$dna_profiles,
                                                  models = streamlined_data_and_models$dna_pool ,
                                                  method = tc)


          model_plot(data_with_models=assign_after_testing,
                     label=paste0("Testset Streamlined Profile. Parameter: ",i),
                     turned_on = detailed_plots_on)





          result_testing_condensed <- predict_new(models=this_condense$models,
                                                  profiles=this_condense$dna_profiles,
                                                  new_data=testing_set,
                                                  closeness_threshold=solution_threshold, #set same as in training set!
                                                  assign_models = T,
                                                  method=tc)

          assign_condensed_after_testing <- model_closeness(data = testing_set,
                                                            profiles=this_condense$dna_profiles,
                                                            models = this_condense$models ,
                                                            method = tc)


          model_plot(data_with_models=assign_condensed_after_testing,
                     label=paste0("Testset Condensed Profile. Parameter: ",i),
                     turned_on = detailed_plots_on)



          streamline_testset_dna_profiles <- create_profile(data_after_dummies=df_full_dummies[-training_indices,],
                                                            model_assignment_data=assign_after_testing,
                                                            dummy_vars=original_variables[dummy_sets],
                                                            model_labels = NULL)

          streamline_compare <- compare_data(true_sample=sample.data$profile,
                                             found_sample=streamline_testset_dna_profiles,
                                             true_case_assignments=sample.data$data,
                                             found_case_assignments = assign_after_testing,
                                             dna_pool=sample.data$models,
                                             found_models=streamlined_data_and_models$dna_pool)

          condense_possibilities[[i]] <-     list(condensed=this_condense,condensed_training=condensed_training, streamlined=streamlined_data_and_models,streamlined_training=result_streamlined_training, testset_streamlined= result_testing, testset_condensed= result_testing_condensed, streamline_compare=streamline_compare )


          ##8. Compare Results to Choose Condensation###


          cat("\nCondensation to",nrow(this_condense$models), "models using style:", this_method_style )
          cat("\nInitial Training with full models. (Discrete  %):", result_training$discrete_percent_solved, "(R-squared):", result_training$r_squared)

          cat("\nCondensed Training  (Discrete  %):", condensed_training$discrete_percent_solved, "(R-squared):", condensed_training$r_squared)


          cat("\nStreamlined Training (Discrete  %):", result_streamlined_training$discrete_percent_solved, "(R-squared):",result_streamlined_training$r_squared)
          cat("\nTest Set Condensed (Discrete  %):", result_testing_condensed$discrete_percent_solved, "(R-squared):",result_testing_condensed$r_squared)


          cat("\nTest Set Streamlined (Discrete  %):", result_testing$discrete_percent_solved, "(R-squared):",result_testing$r_squared)


          #here is the simple regression one-model comparison for test set###

          this_check_ratio <- mean(x=c((as.numeric(result_testing$discrete_percent_solved)/100) , as.numeric( gsub("%", "",  result_testing$r_squared))))/(nrow(this_condense$models)^(1/3))

          cat("\nR-squared to Models ratio:", round(this_check_ratio, digits = 2), "\n\n")




          if (nrow(this_condense$models)>1){

            if (this_check_ratio>best_check_ratio){
              best_check_ratio <- this_check_ratio
              best_check_list <- list(parameter_check=as.integer(i), ratio=this_check_ratio)
            }

          }else{ #if no models




          }



        } #end parameter check



        select_this_condensation <- best_check_list$parameter_check
        print(best_check_list$parameter_check)
        #select_this_condensation <- 3
        cat("\nChoosing this condensation", select_this_condensation,"\n\n")

        #IF condensation ends in only one model###
        if (!is.null( select_this_condensation)){
          streamlined_data_and_models <- condense_possibilities[[select_this_condensation]]$streamlined


        }else{ #if condensation ends only in one model
          #ADD NAs for ICR results
          compare_three_methods <- rbind(compare_three_methods,  c("ICR streamlined",NA, NA, NA, NA, parameter_set,1,  "only one model from condensation" ))
          compare_three_methods <- rbind(compare_three_methods, c("ICR case-based",NA, NA, NA, NA, parameter_set,1, "only one model from condensation" ))

        }


      }else{ #if only one training row, then select this condensation becomes null and all ICR becomes NA
        select_this_condensation <- NULL

        cat("\nNo condensations available. Use single model or adjust condensation selection parameters.\n")


        #ADD NAs for ICR results
        compare_three_methods <- rbind(compare_three_methods,  c("ICR streamlined",NA, NA, NA, NA, parameter_set,1,  "only one model survived training" ))
        compare_three_methods <- rbind(compare_three_methods, c("ICR case-based",NA, NA, NA, NA, parameter_set,1, "only one model survived training" ))


      }  #end if related to one more more training rows






      ##9. Interpretation and conversion to case-based model assignments###



    }else{ #if training was NA, no models
      cat("\nNo condensation conducted because no models survived training","\n")

           select_this_condensation <- NULL

      #ADD NAs for ICR results
      compare_three_methods <- rbind(compare_three_methods,  c("ICR streamlined",NA, NA, NA, NA, parameter_set,0, "no models survived training" ))
      compare_three_methods <- rbind(compare_three_methods, c("ICR case-based",NA, NA, NA, NA, parameter_set,0,"no models survived training" ))


    } #end if, for situation with training with no surviving models

    #For Data with dummy sets
    if (!is.null(select_this_condensation)){ #Continue only if condensation was conducted
      if (T %in% dummy_sets){
        interpret_data <- transform_dummies(training_set, group_terms =original_variables[dummy_sets], reference_cat = reference_categories)

        for (i in 1:length(dummy_set_ivs)){
          this_dummy <- dummy_set_ivs[i]
          interpret_data[,this_dummy] <- as.factor(interpret_data[,this_dummy])
        }


      }else{
        interpret_data <- training_set
      }

    }else{ #if condensation was unsuccessful/ not conducted
      streamlined_data_and_models <- NA
      true_regression_case_based <- NA #there was no solution, children died out
    }



    ##10 a. Agglomerative Subgroup Characteristic Tree###

    if (!is.null(select_this_condensation)){ #if condensation was conducted

      population_tree_results <<- agglom_tree(data=streamlined_data_and_models,
                                             original_data = original_data[training_indices,],
                                             re_use_variables = F)
    }else{
      population_tree_results <<- NA
    }

    ##11. CONVERT FROM POPULATION-BASED TO CASE-BASED Model-Assignments###
    # apply population differences as case-based differences to data set


    #ADD dna profiles to that output!!
    if (!is.null(select_this_condensation)){

      #is cutlist NA?
      if (length(population_tree_results$cutlist)>1){


        case_based_model_assignment <- population_to_case_tree(cutlist=population_tree_results$cutlist,
                                                               population_data=population_tree_results$data)

      }else if (is.null( population_tree_results$data[1])){ #if only one model




        case_based_model_assignment$data <-  training_set
        case_based_model_assignment$data$model <- 1
        # case_based_model_assignment$data$data <- dplyr::rename(case_based_model_assignment$data ,y="intercept")
        case_based_model_assignment$profiles <- streamlined_data_and_models$dna_profiles


      }else{

        case_based_model_assignment <- list()
        case_based_model_assignment$data <- population_tree_results$data

      }
    }#end if related to whether condensation happened

    #update the models by calculating true regressions for those new model populations




    #so we remove columns that are not the model and not in the newdf dataset
    if (!is.null(select_this_condensation)){
      columns_to_keep <- c(colnames(newdf), "model")
      case_based_model_assignment$data <- case_based_model_assignment$data[,colnames(case_based_model_assignment$data) %in% columns_to_keep]



      true_regression_case_based <- streamline(models=NULL,
                                               profiles=case_based_model_assignment$profiles,
                                               data=case_based_model_assignment$data,
                                               assign_models = F)







      result_training_case_based <- predict_new(models=true_regression_case_based$dna_pool,
                                                profiles=true_regression_case_based$dna_profiles,
                                                new_data=true_regression_case_based$model_assignment,
                                                closeness_threshold=solution_threshold, #set same as in training set!
                                                assign_models = F)

      model_plot(data_with_models=case_based_model_assignment$data,
                 label = "Case-based Training",
                 turned_on = detailed_plots_on)

      ##10 b. simplified model descriptions to this###




      ## 12. Run PREDICTIONS for new cases###



      result_testing_case_based <- predict_new(models=true_regression_case_based$dna_pool,
                                               profiles=true_regression_case_based$dna_profiles,
                                               new_data=testing_set,
                                               closeness_threshold=solution_threshold, #set same as in training set!
                                               assign_models = T,
                                               method="case-based",
                                               cutlist=population_tree_results$cutlist
      )

      assign_test_case_based <- result_testing_case_based$model_closeness_data


      model_plot(data_with_models=assign_test_case_based,
                 label = "Case-based Testing",
                 turned_on = detailed_plots_on)





      ##13. Compare results across steps###



      ICR_dna_profiles <-    create_profile(data_after_dummies=df_full_dummies[-training_indices,], #remove model column here
                                            model_assignment_data=assign_test_case_based,
                                            dummy_vars=NULL,
                                            model_labels = NULL)



      compare_ICR <- compare_data(true_sample=sample.data$profile,
                                  found_sample=ICR_dna_profiles,
                                  true_case_assignments=sample.data$data,
                                  found_case_assignments = assign_test_case_based,
                                  dna_pool=sample.data$models,
                                  found_models=true_regression_case_based$dna_pool)



      print("case based compare results")
      print(compare_ICR)
      if (is.na(compare_ICR$model_similarity_error)){
        #
        #investigate NA's in model similarity error
      }




      compare_three_methods <- rbind(compare_three_methods, c("ICR case-based",result_testing_case_based$r_squared, result_testing_case_based$discrete_percent_solved, compare_ICR$case_assign_score, compare_ICR$model_similarity_error, parameter_set,nrow(ICR_dna_profiles), "" ))

      ##ICR streamlined####




      this_discrete <- condense_possibilities[[select_this_condensation]]$testset_streamlined$discrete_percent_solved
      this_r2 <- condense_possibilities[[select_this_condensation]]$testset_streamlined$r_squared
      this_case_score <- condense_possibilities[[select_this_condensation]]$streamline_compare$case_assign_score
      this_model_error <- condense_possibilities[[select_this_condensation]]$streamline_compare$model_similarity_error


      #how many models?
      compare_three_methods <- rbind(compare_three_methods,  c("ICR streamlined",this_r2, this_discrete, this_case_score, this_model_error, parameter_set,nrow(ICR_dna_profiles) ,"" ))
    }else{


    } #end if to contininue with condensation, streamlining, etc.


    ##Simple regression####

    #now the simple regression
    this_sample <- testing_set
    last_column <- ncol(this_sample)
    dv <- colnames(this_sample)[1]
    ivs <- colnames(this_sample)[2:last_column]

    f <- as.formula(
      paste(dv,
            paste(ivs, collapse = " + "),
            sep = " ~ "))

    simple_regression <- eval(bquote(   lm(.(f), data = this_sample)   ))
    #summary(simple_regression)
    simple_model <- t(data.frame(simple_regression$coefficients))




    simple_testset_dna_profiles <- create_profile(data_after_dummies=df_full_dummies[-training_indices,],
                                                  model_assignment_data=cbind(this_sample, model=1),
                                                  dummy_vars=original_variables[dummy_sets],
                                                  model_labels = NULL)



    result_simple_testset<- predict_new(models=simple_model,
                                        profiles=simple_testset_dna_profiles,
                                        new_data=this_sample,
                                        closeness_threshold=solution_threshold,#set same as in training set!
                                        assign_models = T,
                                        method="residuals")


    rownames(simple_model) <- "1"

    compare_simple <- compare_data(true_sample=sample.data$profile,
                                   found_sample=simple_testset_dna_profiles,
                                   true_case_assignments=sample.data$data,
                                   found_case_assignments = cbind(this_sample, model=1),
                                   dna_pool=sample.data$models,
                                   found_models=simple_model)


    compare_three_methods <- rbind(compare_three_methods,   c("simple regression", result_simple_testset$r_squared, result_simple_testset$discrete_percent_solved, compare_simple$case_assign_score, compare_simple$model_similarity_error, parameter_set,1,"" ))












    #ADD back in later ?###
    #Model depiction###
    # substantive_threshold <- .05 #below this is considered not relevant to present
    #
    # depicted_models <- true_regression_case_based$dna_pool
    #
    # passing_threshhold <- 1*(abs(depicted_models[, 2:ncol(depicted_models)])>=substantive_threshold)
    #
    # depicted_models[, 2:ncol(depicted_models)] <- depicted_models[, 2:ncol(depicted_models)] * passing_threshhold
    #
    # if (nrow(depicted_models)>1){
    #   model_freq <- data.frame(table(streamlined_data_and_models$model_assignment$model))
    #
    # }else{
    #   #this works?###
    #   #define model_freq as 100%###
    #
    #   model_freq <- data.frame(matrix(nrow=1, ncol=2))
    #   colnames(model_freq) <- c("var", "Freq")
    #   model_freq$var <- rownames(depicted_models)[1]
    #   model_freq$Freq <- nrow(streamlined_data_and_models$solution_table)
    #
    #   # streamlined_data_and_models$solution_table[,1] <- 1
    #   # model_freq <-  data.frame(table(streamlined_data_and_models$solution_table[,1]))
    # }
    #
    # depicted_models$cases <- model_freq$Freq
    #
    # #the first column, named after the DV, is the intercept!
    # cat("\n\nModel Descriptions")
    # #print(depicted_models)
    # ##Manually add reference categories for interpretation###
    # cat("\nReference Categories:\n")
    # #cat("\nCountry: Eastern Europe")
    # #cat("\ngndr:male")
    # cat("\n")
    #
    # for (i in 1:nrow(depicted_models)){
    #
    #
    #   this_model <- depicted_models[i,-ncol(depicted_models)]
    #   this_model <- this_model[,!(this_model==0 | is.na(this_model)), drop=FALSE]
    #   print( this_model, row.names = T)
    # }






    ##Mixture Regression####




    data <- training_set #form the model on the training set, test on testset
    if (!is.null(select_this_condensation)){
      #
      this_many_models <- nrow(true_regression_case_based$dna_pool)
    }else{
      this_many_models <- 5
    }

    outcome <- colnames(testing_set)[1]


    f <- as.formula(paste(outcome,paste0(colnames(testing_set)[-1], collapse = "+"), sep=" ~ " ))

    ex1 <- flexmix::flexmix(f,data = data, k = this_many_models,
                            control = list(verb = 5, iter = 1000))

    #update to reflect how many models flexmix actually found
    this_many_models <- ex1@k

    #print(table(clusters(ex1), data$class))
    #get coefficients here.
    #parameters(ex1, component = 2)
    #summary(ex1)
    new_data <- testing_set


    mix_model_list <- list()

    for (m in 1:this_many_models){
      x1 <- data.frame(t(flexmix::parameters(ex1, component = m)))
      mix_model_list[[m]] <- x1

    }



    mixture_models <- list_rbind(mix_model_list)


    mixture_models <- mixture_models[,-which(colnames(mixture_models)=="sigma")]



    rownames(mixture_models) <- 1:this_many_models


    #assign with my function or with clusters function
    new_data <- model_closeness(data = new_data,
                                profiles = NULL,
                                models = mixture_models,
                                method = "residuals",
                                cutlist = NA)





    mixture_dna_profiles <-    create_profile(data_after_dummies=df_full_dummies, #remove model column here
                                              model_assignment_data=new_data,
                                              dummy_vars=NULL,
                                              model_labels = NULL)


    result_mixture<- predict_new(models=mixture_models,
                                 profiles=mixture_dna_profiles,
                                 new_data=testing_set,
                                 closeness_threshold=solution_threshold,#set same as in training set!
                                 assign_models = T,
                                 method="residuals")


    mixture_compare_results <- compare_data(true_sample=sample.data$profile,
                                            found_sample=mixture_dna_profiles,
                                            true_case_assignments=sample.data$data,
                                            found_case_assignments = new_data,
                                            dna_pool=sample.data$models,
                                            found_models=mixture_models)


    compare_three_methods <- rbind(compare_three_methods,   c("mixture regression", result_mixture$r_squared, result_mixture$discrete_percent_solved, mixture_compare_results$case_assign_score, mixture_compare_results$model_similarity_error, parameter_set ,this_many_models, ""))


    ##Clusters (k means) plus regression####


    data <- training_set

    #

    k <- kmeans(data, this_many_models)

    #now rotate through each cluster and get a model


    cluster_model_list <- list()

    for (m in 1:this_many_models){

      this_sample <- data[k$cluster==m,]
      last_column <- ncol(this_sample)
      dv <- colnames(this_sample)[1]
      ivs <- colnames(this_sample)[2:last_column]

      f <- as.formula(
        paste(dv,
              paste(ivs, collapse = " + "),
              sep = " ~ "))

      cluster_regression <- eval(bquote(   lm(.(f), data = this_sample)   ))
      #summary(simple_regression)
      cluster_model <- t(data.frame(cluster_regression$coefficients))

      cluster_model_list[[m]] <-as.data.frame( cluster_model)

    }

    #continue with testing of clusters plus regression.


    cluster_models <- list_rbind(cluster_model_list)





    rownames(cluster_models) <- 1:this_many_models
    #replace NA's with zero
    cluster_models[is.na(cluster_models)] <- 0


    #assign with my function or with clusters function
    new_data <- model_closeness(data = data,
                                profiles = NULL,
                                models = cluster_models,
                                method = "residuals",
                                cutlist = NA)


    #assign cluster models!
    clust <- k$cluster
    data$model <-clust



    cluster_dna_profiles <-    create_profile(data_after_dummies=df_full_dummies, #remove model column here
                                              model_assignment_data=data,
                                              dummy_vars=NULL,
                                              model_labels = NULL)



    result_cluster<- predict_new(models=cluster_models,
                                 profiles=cluster_dna_profiles,
                                 new_data=testing_set,
                                 closeness_threshold=solution_threshold,#set same as in training set!
                                 assign_models = T,
                                 method="residuals")


    cluster_compare_results <- compare_data(true_sample=sample.data$profile,
                                            found_sample=cluster_dna_profiles,
                                            true_case_assignments=sample.data$data,
                                            found_case_assignments = data,
                                            dna_pool=sample.data$models,
                                            found_models=cluster_models)




    compare_three_methods <- rbind(compare_three_methods,    c("cluster regression", result_cluster$r_squared, result_cluster$discrete_percent_solved, cluster_compare_results$case_assign_score, cluster_compare_results$model_similarity_error, parameter_set,nrow(cluster_dna_profiles), "" ))

    cat("\n Tested:",ra, ge, so, ns, re, ni, er, tc, "\n")

    #for all compared sets this parameter set










    #1. compare ICR and mixture in terms of model complexity (higher IVs, higher error )











    colnames(compare_three_methods) <- c("method", "weighted r2", "discrete solution", "case assignment score", "model similarity error","seed", "generations", "solution_thresh", "n_dna_strands", "real_model_n",
                                         "N_IVs", "error_sd", "test_set_model_closeness","n models", "notes")

  #  saveRDS(compare_three_methods, paste0( "./testing/", paste0(parameter_set,collapse = "_"), ".RDS"))
    return(compare_three_methods)
  }#end foreach loop####

  if (parallelize==T){
   # parallel::clusterEvalQ(my_cluster,  sink())
    parallel::stopCluster(my_cluster)
  }else{

  }



  par_results$parameter_ID <- cumsum(!duplicated(par_results[,6:13]))
  View(par_results)

  #saveRDS(par_results, "./testing/fresh parameter testing parallelized.RDS")

  #saveRDS(compare_three_methods, paste0("./testing/",ra, "new three method comparison D.RDS"))


  # saveRDS(parameter_testing, "generation testing parameter testing.RDS")



   these_parameter_files <- par_results


  #look for missing tests

 if (boxplots==T){
  ##boxplots####

  these_parameter_files$method <- as.factor( these_parameter_files$method )
  these_parameter_files$generations <- as.factor( these_parameter_files$generations )

  #now identify all NAs and their corresponding parameter IDs
  remove_these_ids <- these_parameter_files$parameter_ID[  which(is.na(these_parameter_files$`weighted r2`))]

  limited_data <- these_parameter_files[!these_parameter_files$parameter_ID %in% remove_these_ids,]
  limited_data <- suppressWarnings( type.convert(limited_data))
  limited_data$method <- as.factor( limited_data$method )
  limited_data$generations <- as.factor( limited_data$generations )




  #limited_data <- these_parameter_files[these_parameter_files$real_model_n==16,]

  ##Remove redundant parameter checks###
  #we don't need multiple checks of simple, cluster + regression, or LCA + regression, because they are unimpacted by generations (Besides a different random seed)

  keep_this_gen <- generations[length(generations)]
  apply_to_these_methods <- c("cluster regression","mixture regression", "simple regression"  )
  #extraneous parameter tests are removed for the above cases
  limited_data <- limited_data[!((limited_data$method %in% apply_to_these_methods) & !(limited_data$generations %in% keep_this_gen)),]



  ##boxplot###
  without_streamlined <- limited_data[!limited_data$method=="ICR streamlined",]

  x1 <- ggplot(without_streamlined, aes(x=method, y=without_streamlined$`weighted r2`,  fill=generations)) +
#    x1 <- ggplot(without_streamlined, aes(x=method, y=without_streamlined$`weighted r2`,  fill=generations, color=factor(starting_pop))) +
       geom_boxplot(notch=F)
  png("./testing/R2.png")
  print(x1)
  dev.off()
  print(x1)

    x2<-  ggplot(without_streamlined, aes(x=method, y=without_streamlined$`discrete solution`,  fill=generations)) +
     # x2<-  ggplot(without_streamlined, aes(x=method, y=without_streamlined$`discrete solution`,  fill=generations,  color=factor(starting_pop))) +

        geom_boxplot(notch=F)
  png("./testing/discrete.png")
  print(x2)
  dev.off()
  print(x2)

  without_regression <- without_streamlined[!without_streamlined$method=="simple regression",]
  x3<-  ggplot(without_regression, aes(x=method, y=without_regression$`case assignment score`,  fill=generations)) +
#    x3<-  ggplot(without_regression, aes(x=method, y=without_regression$`case assignment score`,  fill=generations, color=factor(starting_pop))) +
    geom_boxplot(notch=F)
  png("./testing/case accuracy.png")
  print(x3)
  dev.off()
  print(x3)

  x4<-  ggplot(without_streamlined, aes(x=method, y=without_streamlined$`model similarity error`,  fill=generations)) +
    #x4<-  ggplot(without_streamlined, aes(x=method, y=without_streamlined$`model similarity error`,  fill=generations, color=factor(starting_pop))) +

       geom_boxplot(notch=F)
  png("./testing/model error.png")
  print(x4)
  dev.off()
  print(x4)


 }#end boxplots condition

    return(par_results)


}


#' ICR
#'
#' Run the ICR algo
#' @param generations An integer specifying the number of generations for initial ICR genetic training.
#' @param solution_thresh A positive numeric value specifying the discrete solution target on each side of the solution, where solutions may range from 0 to 1.
#' @param n_strands Survivors. An integer indicating the number of population that survive each generation of the genetic training algorithm.
#' @param detailed_plots_on Produces bi-dimentional X-Y plots for various steps of the ICR algorithm. Only useful for low-dimensional analysis.
#' @param input_data A pure data.frame containing the dependent and independent variables. If not specified otherwise, the first column will be used as the dependent variable. If this df contains a column called 'model', it will be removed. All dummy variables must be coded as numerical and 0-1, factors of more than 3 categories should be entered as char columns.
#' @param dv The dependent variable. If not specified, the first column is used.
#' @param ivs The independent variable. If not specified, all variables besides the dv are used.
#' @param starting_pop Starting sets of random coefficients (population, children) in first generation of genetic algorithm.
#' @keywords ICR
#' @examples
#'  icr_results <- icr(generations=400,
#' input_data=ess2016,
#' dv="happy",
#' iv=c("female", "age", "cntry", "ppltrst", "sclmeet",      "income", "health"),
#' starting_pop=10000, #1000
#' n_strands = 100,
#' mutation_rate = .02,
#' solution_thresh = .003,
#' n_children=3,
#' death_by_ageing =60,
#' nelitism=10)
#' @export
icr <- function(generations=c(1000),
                          solution_thresh = c( .01 ),
                          n_strands = c(50), #parameter for testing
                          detailed_plots_on=F,
                          input_data,
                          dv,
                          iv,
                          starting_pop=10000,
                          mutation_rate=.2,
                          n_children=2,
                          death_by_ageing =60,
                          nelitism=10){

  options(scipen=999)

  input_data <- data.frame(input_data)





tc <- "residuals"

  ###data prep start
  original_variables <<- c(dv, iv)

  df <- input_data[,original_variables]
  dv<<-dv
  #throw error if DV is not numerical

if(  is.numeric(df[,dv])){

}else{
  stop("Dependent variable must be numeric.")
}

  #now automatically find which are dummysets

which_maybe_dummysets <-   sapply(df, class) == "character" |  sapply(df, class) == "factor"

#now verify more than 2 categories
more_than_two_cats <- sapply(unique(df[, which_maybe_dummysets,drop=FALSE]), FUN = function(x) length( unique((x))))>2
which_maybe_dummysets[names(more_than_two_cats)] <- more_than_two_cats
dummy_sets <<- which_maybe_dummysets #made global because this never changes from original data


dummy_set_ivs <- names(dummy_sets)[dummy_sets]
dummy_set_ivs <- dummy_set_ivs[ !dummy_set_ivs==dv]
  which_not_dummies <- which(!dummy_sets)
  #!sapply(df, class) == 'character'

  df[,which_not_dummies] <- sapply(df[,which_not_dummies], rescale)
  df <- data.frame(df)

  #cut out any rows with NA values in any column

  df <- df[complete.cases(df),]

    newdf <- df
    original_data <<- df #before dummy splits!




  ###Form dummysets###
    if (T %in% dummy_sets){
  save_row_names <- rownames(newdf)
  newdf <- dummy_cols(newdf, select_columns = original_variables[dummy_sets], remove_first_dummy = T, ignore_na = T, remove_selected_columns = T)
  df_full_dummies <- dummy_cols(df, select_columns = original_variables[dummy_sets], remove_first_dummy = F, ignore_na = T, remove_selected_columns = T)
  df_full_dummies<<- df_full_dummies
   #the above function reset the row names!

  rownames(newdf) <- save_row_names
  reference_categories <- as.list(df [,original_variables[dummy_sets], drop=F])
  reference_categories <-   lapply(reference_categories, FUN = function(x) sort(unique(x) ))
  reference_categories <- sapply(reference_categories,FUN = function(x) unlist(x)[1])
  reference_categories <<- reference_categories
  #newdf is AFTER dummy splits

    }else{
      df_full_dummies<<- newdf
    }
  ###data prep end








  #Analysis routine




    ge <- generations
    so <- solution_thresh
    ns <- n_strands




    parameter_set <- c( ge, so, ns)
    # saveRDS(parameter_testing, paste0(ra,ge,so,ns,re, ni, er, tc, "parameter testing.RDS"))
    #print is for for debugging
    print("these parameters")
    print(parameter_set)




    #collect output row into parameter_testing table




    if (detailed_plots_on==T){
      simple <-  ggplot(newdf, aes(x=newdf[, iv[1]], y=newdf[,dv])) +
        geom_point() +
        geom_smooth(method=lm, se=FALSE)+
        ggtitle("Simple Regression")
      print(simple)
    }

    #!FIX TO WORK WITH NON SYNTHETIC DATA###
    model_plot(data_with_models=newdf, label = "True Models",
               models = sample.data$models,
               base_line_plot = F,
               turned_on = detailed_plots_on,
               X_variable = 1,
               local_regression_line = F)



    #remove the models from the sample data

    newdf <- newdf[,!colnames(newdf)=="model"]




    ##2. SPLIT DATA INTO TRAINING AND TEST SETS###

    training_proportion <- 0.6
    training_set_N <- round(nrow(newdf)*training_proportion)
    training_indices <<- sample(nrow(newdf), training_set_N, replace=F)
    training_set <<- as.data.frame(newdf[training_indices, ])
    testing_set <<- as.data.frame(newdf[-training_indices, ])

    ##3. SPAWN CHILDREN###
    children <- children_spawn(n_children= starting_pop,
                               dna_length=ncol(newdf), #this should be the correct number of columns. We need one coefficient for every IV column (including dummies) and one forthe intercept (the DV column)
                               col_names=c('intercept',colnames(newdf)[2:length(colnames(newdf))]), #the DV name column is replaced with 'intercept'
                               normal=T, #normally distributed with the following SD. If not, then this_range is used
                               this_sd=.5, #standard deviation of the coefficients.. this shoudl be tailored to the data's form.
                               this_range=c(-10, 10))

    ##4.  RUN TRAINING SET###

    cat("Running training set \n")


    solution_threshold <<- so

#Test children call####
    result_training <- test_children(dna_pool=children,
                                     input_data=training_set,
                                     closeness_threshold=solution_threshold, #defines how close a solution is to the found solution
                                     generations=ge,
                                     difficulty=.8,#currently unused
                                     death_rating_thresh=0,
                                     death_number_thresh=ns,
                                     uniqueness_power=8,
                                     death_by_ageing = death_by_ageing,
                                     number_of_children = n_children,
                                     min_strand_solution_proportion=0, #valid solutions must cover at least this proportion of cases.
                                     mutations = .2,
                                     mating = "random",
                                     mating_power = 2,
                                     selection = "probability",
                                     nelitism = nelitism)


    # saveRDS(result_training,paste0("./testing/",i, " training results.RDS"))
    #
    ### IF no training. training is NA###
    if (!is.na(result_training[1])){ #if the dna did not die out, then continue with condensation ONLY



      assign_after_training <- model_closeness(data = training_set, result_training$dna_profiles,
                                               models = result_training$dna_pool ,
                                               method = tc)


      model_plot(data_with_models=assign_after_training, label = "Post training",
                 models = sample.data$models,
                 base_line_plot = F,
                 turned_on = detailed_plots_on,
                 X_variable = 1,
                 local_regression_line = T)





      ##5. CONDENSE Training models, then streamline, and then re-run training set and try test set###
      #this is all linked together into a loop for troubleshooting and optimizing.. to find the best condensation-streamlining combination



      #only condense if more than one model was found!
      ###IF more than one training model found###
      if (nrow(result_training$dna_pool)>1){
        condense_max <- 6
        condense_min <- 2
        parameter_check <- 1:(length(condense_min:condense_max)+1)
        condense_number <- c( condense_min:condense_max, "NULL")
        method_style <- c( rep("best", length(condense_min:condense_max)),"cluster" ) #condensation method
        best_check_ratio <- 0
        minimum_models <- 2
        best_check_list <- list()
        condense_possibilities <- vector("list", length =length(parameter_check) )

        for (i in parameter_check){


          this_condense_number <- condense_number[i]
          this_method_style <- method_style[i]




          this_condense <- condense_models(data = result_training,
                                           training_set = training_set,
                                           method = this_method_style, #or "cluster" / "best"
                                           auto_cluster = TRUE, # indicates if algorithm should find the ideal number of cluster by itself
                                           nmodel = this_condense_number, # only relevant for method = best or when auto_cluster is set to FALSE
                                           solution_table=result_training$solution_table, #relevant only for genetic method
                                           no_correlation_thresh =-.1)

        #  print("this condense")
        #  print(this_condense$models)

          model_plot(data_with_models=this_condense$data_with_models,
                     label = paste0("Condense training . Parameter: ",i),
                     turned_on = detailed_plots_on,
                     X_variable=1)

          condensed_training <- predict_new(models=this_condense$models,
                                            profiles=this_condense$dna_profiles,
                                            new_data=this_condense$data_with_models,
                                            closeness_threshold=solution_threshold, #set same as in training set!
                                            assign_models = F)

          ##6. STREAMLINE Take solution table and calculate true regressions for underlying subpopulations###
          #This makes predictions more accurate and more clearly separates subpopulations


          streamlined_data_and_models <- streamline(models=this_condense$models,
                                                    profiles=this_condense$dna_profiles,
                                                    data=this_condense$data_with_models,
                                                    assign_models = F)





          model_plot(data_with_models=streamlined_data_and_models$model_assignment,
                     label = paste0("Streamlined Training. Parameter: ",i),
                     turned_on = detailed_plots_on,
                     X_variable = 1,
                     base_line_plot = F,
                     models = sample.data$models)


          #now test how well we can predict training cases with the model

          result_streamlined_training <- predict_new(models=streamlined_data_and_models$dna_pool,
                                                     profiles=streamlined_data_and_models$dna_profiles,
                                                     new_data=streamlined_data_and_models$model_assignment,
                                                     closeness_threshold=solution_threshold, #set same as in training set!
                                                     assign_models = F)







          ##7.  TEST SET###


          result_testing <- predict_new(models=streamlined_data_and_models$dna_pool,
                                        profiles=streamlined_data_and_models$dna_profiles,
                                        new_data=testing_set,
                                        closeness_threshold=solution_threshold, #set same as in training set!
                                        assign_models = T,
                                        method=tc)

          assign_after_testing <- model_closeness(data = testing_set,
                                                  profiles=streamlined_data_and_models$dna_profiles,
                                                  models = streamlined_data_and_models$dna_pool ,
                                                  method = tc)


          model_plot(data_with_models=assign_after_testing,
                     label=paste0("Testset Streamlined Profile. Parameter: ",i),
                     turned_on = detailed_plots_on)





          result_testing_condensed <- predict_new(models=this_condense$models,
                                                  profiles=this_condense$dna_profiles,
                                                  new_data=testing_set,
                                                  closeness_threshold=solution_threshold, #set same as in training set!
                                                  assign_models = T,
                                                  method=tc)

          assign_condensed_after_testing <- model_closeness(data = testing_set,
                                                            profiles=this_condense$dna_profiles,
                                                            models = this_condense$models ,
                                                            method = tc)


          model_plot(data_with_models=assign_condensed_after_testing,
                     label=paste0("Testset Condensed Profile. Parameter: ",i),
                     turned_on = detailed_plots_on)


          streamline_testset_dna_profiles <- create_profile(data_after_dummies=df_full_dummies[-training_indices,],
                                                            model_assignment_data=assign_after_testing,
                                                            dummy_vars=original_variables[dummy_sets],
                                                            model_labels = NULL)

          streamline_compare <- NA #only valid for benchmarking function
          condense_possibilities[[i]] <-     list(condensed=this_condense,condensed_training=condensed_training, streamlined=streamlined_data_and_models,streamlined_training=result_streamlined_training, testset_streamlined= result_testing, testset_condensed= result_testing_condensed, streamline_compare=streamline_compare )


          ##8. Compare Results to Choose Condensation###
          #
          # cat("\nCondensation to",nrow(this_condense$models), "models using style:", this_method_style )
          # cat("\nInitial Training with full models. (Discrete  %):", result_training$discrete_percent_solved, "(R-squared):", result_training$r_squared)
          #
          # cat("\nCondensed Training  (Discrete  %):", condensed_training$discrete_percent_solved, "(R-squared):", condensed_training$r_squared)
          #
          #
          # cat("\nStreamlined Training (Discrete  %):", result_streamlined_training$discrete_percent_solved, "(R-squared):",result_streamlined_training$r_squared)
          # cat("\nTest Set Condensed (Discrete  %):", result_testing_condensed$discrete_percent_solved, "(R-squared):",result_testing_condensed$r_squared)
          #
          #
          # cat("\nTest Set Streamlined (Discrete  %):", result_testing$discrete_percent_solved, "(R-squared):",result_testing$r_squared)
          #

          #here is the simple regression one-model comparison for test set###

          this_check_ratio <- mean(x=c((as.numeric(result_testing$discrete_percent_solved)/100) , as.numeric( gsub("%", "",  result_testing$r_squared))))/(nrow(this_condense$models)^(1/3))

          # cat("\nR-squared to Models ratio:", round(this_check_ratio, digits = 2), "\n\n")
          #



          if (nrow(this_condense$models)>1){

            if (this_check_ratio>best_check_ratio){
              best_check_ratio <- this_check_ratio
              best_check_list <- list(parameter_check=as.integer(i), ratio=this_check_ratio)
            }

          }else{ #if no models




          }



        } #end parameter check



        select_this_condensation <- best_check_list$parameter_check
        print(best_check_list$parameter_check)
        #select_this_condensation <- 3
        cat("\nChoosing this condensation", select_this_condensation,"\n\n")


        #####IF condensation ends in only one model###
        if (!is.null( select_this_condensation)){
          streamlined_data_and_models <- condense_possibilities[[select_this_condensation]]$streamlined


        }else{ #if condensation ends only in one model
          #ADD NAs for ICR results
          compare_three_methods <- rbind(compare_three_methods,  c("ICR streamlined",NA, NA, NA, NA, parameter_set,1,  "only one model from condensation" ))
          compare_three_methods <- rbind(compare_three_methods, c("ICR case-based",NA, NA, NA, NA, parameter_set,1, "only one model from condensation" ))

        }


      }else{ #if only one training row, then select this condensation becomes null and all ICR becomes NA
        select_this_condensation <- NULL

        cat("\nNo condensations available. Use single model or adjust condensation selection parameters.\n")


        #ADD NAs for ICR results
        compare_three_methods <- rbind(compare_three_methods,  c("ICR streamlined",NA, NA, NA, NA, parameter_set,1,  "only one model survived training" ))
        compare_three_methods <- rbind(compare_three_methods, c("ICR case-based",NA, NA, NA, NA, parameter_set,1, "only one model survived training" ))


      }  #end if related to one more more training rows


      ####! selection for the best condensation-streamlining mix in terms of interpretability###




      ###9. Interpretation and conversion to case-based model assignments###



    }else{ #if training was NA, no models because dna population died out
      cat("\nNo condensation conducted because no models survived training","\n")

         select_this_condensation <- NULL
      #streamlined_data_and_models <- result_training

      #ADD NAs for ICR results
      compare_three_methods <- rbind(compare_three_methods,  c("ICR streamlined",NA, NA, NA, NA, parameter_set,0, "no models survived training" ))
      compare_three_methods <- rbind(compare_three_methods, c("ICR case-based",NA, NA, NA, NA, parameter_set,0,"no models survived training" ))


    } #end if, for situation with training with no surviving models

    #Post condensation

    #For Data with dummy sets
    if (!is.null(select_this_condensation)){ #Continue only if condensation was conducted
      if (T %in% dummy_sets){

        interpret_data <- transform_dummies(training_set, group_terms =original_variables[dummy_sets], reference_cat = reference_categories)


        for (i in 1:length(dummy_set_ivs)){

        this_dummy <- dummy_set_ivs[i]
          interpret_data[,this_dummy] <- as.factor(interpret_data[,this_dummy])
        }

         }else{
        interpret_data <- training_set
      }

    }else{ #if condensation was unsuccessful/ not conducted
      streamlined_data_and_models <- NA
      true_regression_case_based <- NA #there was no solution, children died out
    }



    ###10 a. Agglomerative Subgroup Characteristic Tree###

    if (!is.null(select_this_condensation)){ #if condensation was conducted
#adapt model names here!

      #change model names in the model assignment and in the dnaprofiles
      models <- rownames( streamlined_data_and_models$dna_profiles)
      new_model_names <- 1:length(models)
      new_model_data_order <- match(streamlined_data_and_models$model_assignment$model, models)
      rownames( streamlined_data_and_models$dna_profiles) <- new_model_names
      streamlined_data_and_models$model_assignment$model <- new_model_data_order

      population_tree_results <<- agglom_tree(data=streamlined_data_and_models,
                                             original_data = original_data[training_indices,],
                                             re_use_variables = F)


    }else{
      population_tree_results <<- NA
    }

    ##11. CONVERT FROM POPULATION-BASED TO CASE-BASED Model-Assignments###
    # apply population differences as case-based differences to data set


    #ADD dna profiles to that output!!
    if (!is.null(select_this_condensation)){

      #is cutlist NA?
      if (length(population_tree_results$cutlist)>1){

     #this is the training data being converted into case-base assignments from population-based ones.
        case_based_model_assignment <- population_to_case_tree(cutlist=population_tree_results$cutlist,
                                                               population_data=population_tree_results$data)

      }else if (is.null( population_tree_results$data[1])){ #if only one model




        case_based_model_assignment$data <-  training_set
        case_based_model_assignment$data$model <- 1
        # case_based_model_assignment$data$data <- dplyr::rename(case_based_model_assignment$data ,y="intercept")
        case_based_model_assignment$profiles <- streamlined_data_and_models$dna_profiles


      }else{

        case_based_model_assignment <- list()
        case_based_model_assignment$data <- population_tree_results$data

      }
    }#end if related to whether condensation happened

    #update the models by calculating true regressions for those new model populations (streamlining again)
    #Some models naturally lack cases in the conversion to case-based modeling???###

    #now we no longer need those extra interaction columns from the dummy sets
    # because the cases are already assigned based on the cutlist

    #so we remove columns that are not the model and not in the newdf dataset
    if (!is.null(select_this_condensation)){
      columns_to_keep <- c(colnames(newdf), "model")
      case_based_model_assignment$data <- case_based_model_assignment$data[,colnames(case_based_model_assignment$data) %in% columns_to_keep]


#check number of rows!

      true_regression_case_based <- streamline(models=NULL,
                                               profiles=case_based_model_assignment$profiles,
                                               data=case_based_model_assignment$data,
                                               assign_models = F)







      result_training_case_based <- predict_new(models=true_regression_case_based$dna_pool,
                                                profiles=true_regression_case_based$dna_profiles,
                                                new_data=true_regression_case_based$model_assignment,
                                                closeness_threshold=solution_threshold, #set same as in training set!
                                                assign_models = F)

      model_plot(data_with_models=case_based_model_assignment$data,
                 label = "Case-based Training",
                 turned_on = detailed_plots_on)

      ##10 b. simplified model descriptions to this###




      ## 12. Run PREDICTIONS for new cases###

    #  testing_set$model <- unique(population_tree_results$data$model)
      result_testing_case_based <- predict_new(models=true_regression_case_based$dna_pool,
                                               profiles=true_regression_case_based$dna_profiles,
                                               new_data=testing_set,
                                               closeness_threshold=solution_threshold, #set same as in training set!
                                               assign_models = T,
                                               method="case-based",
                                               cutlist=population_tree_results$cutlist
      )



      assign_test_case_based <- result_testing_case_based$model_closeness_data

      model_plot(data_with_models=assign_test_case_based,
                 label = "Case-based Testing",
                 turned_on = detailed_plots_on)




      ##ICR case-based###
      #get testset profiles!


      ICR_dna_profiles <-    create_profile(data_after_dummies=df_full_dummies[-training_indices,], #remove model column here
                                            model_assignment_data=assign_test_case_based,
                                            dummy_vars=dummy_set_ivs,
                                            model_labels = NULL)








          }else{

      # compare_three_methods <- rbind(compare_three_methods,  c("ICR streamlined",NA, NA, NA, NA, parameter_set ))
      #compare_three_methods <- rbind(compare_three_methods, c("ICR case-based",NA, NA, NA, NA, parameter_set ))


    } #end if to contininue with condensation, streamlining, etc.





    #ADD back in later for package?###
    #Model depiction###
    # substantive_threshold <- .05 #below this is considered not relevant to present
    #
    # depicted_models <- true_regression_case_based$dna_pool
    #
    # passing_threshhold <- 1*(abs(depicted_models[, 2:ncol(depicted_models)])>=substantive_threshold)
    #
    # depicted_models[, 2:ncol(depicted_models)] <- depicted_models[, 2:ncol(depicted_models)] * passing_threshhold
    #
    # if (nrow(depicted_models)>1){
    #   model_freq <- data.frame(table(streamlined_data_and_models$model_assignment$model))
    #
    # }else{
    #   #this works?###
    #   #define model_freq as 100%###
    #
    #   model_freq <- data.frame(matrix(nrow=1, ncol=2))
    #   colnames(model_freq) <- c("var", "Freq")
    #   model_freq$var <- rownames(depicted_models)[1]
    #   model_freq$Freq <- nrow(streamlined_data_and_models$solution_table)
    #
    #   # streamlined_data_and_models$solution_table[,1] <- 1
    #   # model_freq <-  data.frame(table(streamlined_data_and_models$solution_table[,1]))
    # }
    #
    # depicted_models$cases <- model_freq$Freq
    #
    # #the first column, named after the DV, is the intercept!
    # cat("\n\nModel Descriptions")
    # #print(depicted_models)
    # ##Manually add reference categories for interpretation###
    # cat("\nReference Categories:\n")
    # #cat("\nCountry: Eastern Europe")
    # #cat("\ngndr:male")
    # cat("\n")
    #
    # for (i in 1:nrow(depicted_models)){
    #
    #
    #   this_model <- depicted_models[i,-ncol(depicted_models)]
    #   this_model <- this_model[,!(this_model==0 | is.na(this_model)), drop=FALSE]
    #   print( this_model, row.names = T)
    # }

    #OUTPUT
    #test set: assign_test_case_based
    #true_regression_case_based$model assignment
    #dna profiles of entire model set
    #training and testing statistics

    original_data$model <- NA
    original_data$model[training_indices] <- true_regression_case_based$model_assignment$model
    original_data$model[-training_indices] <- assign_test_case_based$model

    original_after_dummies_df <- rbind(training_set, testing_set)
    integer_rows <- as.integer(rownames(original_after_dummies_df))
    new_order <- match(sort(integer_rows), integer_rows)

    original_after_dummies_df <- original_after_dummies_df[new_order,]
    original_after_dummies_df$model <- original_data$model




    #profile for overall results on original data
    overall_ICR_profiles <-    create_profile(data_after_dummies=df_full_dummies, #remove model column here
                                          model_assignment_data=original_after_dummies_df,
                                          dummy_vars=NULL,
                                          model_labels = NULL)



#simple regression comparison###
 #now the simple regression
 #cut out model columns
 cut_models <- grepl("model", colnames(original_after_dummies_df))

 this_sample <- original_after_dummies_df[,!cut_models]
 last_column <- ncol(this_sample)
 dv <- colnames(this_sample)[1]
 ivs <- colnames(this_sample)[2:last_column]

 f <- as.formula(
   paste(dv,
         paste(ivs, collapse = " + "),
         sep = " ~ "))

 simple_regression <- eval(bquote(   lm(.(f), data = this_sample)   ))
 #summary(simple_regression)
 simple_model <- t(data.frame(simple_regression$coefficients))
 # use create profile function!!



 simple_testset_dna_profiles <- create_profile(data_after_dummies=df_full_dummies[-training_indices,],
                                               model_assignment_data=cbind(this_sample, model=1),
                                               dummy_vars=original_variables[dummy_sets],
                                               model_labels = NULL)



 result_simple_testset<- predict_new(models=simple_model,
                                     profiles=simple_testset_dna_profiles,
                                     new_data=this_sample,
                                     closeness_threshold=solution_threshold,#set same as in training set!
                                     assign_models = T,
                                     method="residuals")


  return(list(original_data_with_new_models=original_data,  cutlist=population_tree_results$cutlist,agglomeration_tree=population_tree_results$tree_blueprint,models= true_regression_case_based$dna_pool,  model_subgroup_profiles=overall_ICR_profiles, training_results=result_training_case_based, testing_results=result_testing_case_based, simple_regression_results=result_simple_testset, training_perc_solved=result_training$discrete_percent_solved))


}


#' Create Sample Data
#'
#' Produces synthetic data
#' @export
#? Adjust coefficients for sample data so that their values have a wider possible range?###
sample_data <- function(c.number = 300,  # number of cases
                        m.number = 5, # number of distinct models
                        m_sizes = c(100,50,50,50,50), #  how many cases each model explains
                        bin.v.number = 5, # number of binary variables
                        con.v.number = 5, # number of continuous variables
                        error = TRUE, # adds error term to the models
                        error.sd = 1, # sd of the error term
                        data.form = c("raw", "standardized", "normalized"),# form of the data
                        messy_subgroups =T, #T means that subgroups might be defined only one variable. F means all variables are subgroup specific
                        min_variable_split=.25 #subgroup means split by at least this much. Also to start each variable may be split only ONCE in these sample data, despite whether it is continuous or binary
                        ){




  v.number <- bin.v.number + con.v.number
  # Create randomized IVs for database
  df <- data.frame(matrix(NA, nrow = c.number, ncol = v.number+1))
  profiles <- data.frame(matrix(0, nrow = m.number, ncol = v.number))
df$model <- 0


#implement a tree!
build_list <- data.frame(matrix(nrow=0, ncol = 6))
colnames(build_list) <- c("cut", "leaf1", "leaf2", "leaf1_mean", "leaf2_mean", "parent_node")

model_list <- 1:m.number
available_nodes <- model_list
variables <- 1:v.number #IVs
variable_vector <- c("continuous", "binary")[c(rep(1, con.v.number), rep(2, bin.v.number))] #true means this variable has a special mean, SD. False means random
available_variables <-  variables
master_parent <-  paste(  as.character( available_nodes), collapse = "&")

build_complete <- m.number-1



#this is the list that builds the subgroups for the diff models
for (i in 1:build_complete){ #so long as the build list is not complete

two_nodes <- sample(available_nodes, 2, replace = F)

this_variable_to_split <- sample(available_variables, 1)

node1_mean <-  sample(0:100, 1)/100
#the second mean is separated from the first by a gap
node2_mean <- sample(c(0:max(((node1_mean*100)- (min_variable_split*100)),0),min(((node1_mean*100)+ (min_variable_split*100)),100):100)  , 1)/100
parent_node <-  paste( as.character( two_nodes), collapse = "&")
this_build_row <- c(this_variable_to_split, two_nodes[1], two_nodes[2], node1_mean, node2_mean, parent_node)
build_list <-rbind(build_list,this_build_row )

#adapt the available nodes and available variables
available_nodes <- available_nodes[!available_nodes %in% two_nodes]
available_nodes <- c(available_nodes,parent_node )


available_variables <- available_variables[!available_variables== this_variable_to_split]


if (length( available_variables)==0){
  break
}

} #end build list construction

colnames(build_list) <- c("cut", "leaf1", "leaf2", "leaf1_mean", "leaf2_mean", "parent_node")


original_build_list <- build_list
df$model <- rep(1:m.number,m_sizes)

while(nrow(build_list)>0){
which_biggest_parent_node <- which(nchar(build_list$parent_node)==max(nchar(build_list$parent_node)))[1]
this_row <- build_list[which_biggest_parent_node,]
this_parent_node <- build_list$parent_node[which_biggest_parent_node]
IV_column_index <-as.integer( this_row$cut)+1 #index of this IV
leaf1 <- this_row$leaf1
leaf2 <- this_row$leaf2
l1_members <- as.integer(unlist(stringr::str_split(leaf1, "&")))
l2_members <- as.integer(unlist(stringr::str_split(leaf2, "&")))
l1_n <- nrow(df[df$model%in% l1_members, IV_column_index, drop=F])
l2_n<- nrow(df[df$model%in% l2_members, IV_column_index, drop=F])

this_variable_type <- variable_vector[as.integer( this_row$cut)]

#check whether binary or continuous


if (this_variable_type=="continuous"){


df[df$model %in% l1_members, IV_column_index] <-  rnorm(n=l1_n, mean = as.numeric( this_row$leaf1_mean), sd = .1)
df[df$model%in% l2_members, IV_column_index] <-  rnorm(n=l2_n, mean = as.numeric( this_row$leaf2_mean), sd = .1)

#compress anything above 1 or under 0 that that range
df[, IV_column_index][df[,IV_column_index]<0] <- 0
df[, IV_column_index][df[,IV_column_index]>1] <- 1

}else if (this_variable_type=="binary"){

  df[df$model %in% l1_members, IV_column_index] <-  rbinom(n=l1_n, size=1, prob= as.numeric( this_row$leaf1_mean))
  df[df$model%in% l2_members, IV_column_index] <-  rbinom(n=l2_n, size=1, prob = as.numeric( this_row$leaf2_mean))


}

build_list <- build_list[-which_biggest_parent_node,]

}

#make a second build list sequence for a proportion of the variables NOT attached to group-specific models
#start Second build list####
#implement a tree!
#this is the non-model build list!

nm_build_list <- data.frame(matrix(nrow=0, ncol = 6))
colnames(nm_build_list) <- c("cut", "leaf1", "leaf2", "leaf1_mean", "leaf2_mean", "parent_node")

#how many independent variables do NOT appear in the initial build list.
#only some of those enter the second list
ind_variables_left <- variables[!(variables %in% as.integer(original_build_list$cut))]
#number of non-model groups
n_nm_groups <- sample(1:length(ind_variables_left), 1, replace = F)


available_variables <-  ind_variables_left
master_parent <-  paste(  as.character( available_nodes), collapse = "&")
available_nodes <- 1:n_nm_groups
build_complete <- n_nm_groups-1

if (build_complete>0){
#this is the list that builds the subgroups for the diff models
for (i in 1:build_complete){ #so long as the build list is not complete

  two_nodes <- sample(available_nodes, 2, replace = F)

  this_variable_to_split <- sample(available_variables, 1)
  node1_mean <-  sample(0:100, 1)/100
  #the second mean is separated from the first by a gap
  node2_mean <- sample(c(0:max(((node1_mean*100)- (min_variable_split*100)),0),min(((node1_mean*100)+ (min_variable_split*100)),100):100)  , 1)/100
  parent_node <-  paste( as.character( two_nodes), collapse = "&")
  this_build_row <- c(this_variable_to_split, two_nodes[1], two_nodes[2], node1_mean, node2_mean, parent_node)
  nm_build_list <-rbind(nm_build_list,this_build_row )

  #adapt the available nodes and available variables
  available_nodes <- available_nodes[!available_nodes %in% two_nodes]
  available_nodes <- c(available_nodes,parent_node )


  available_variables <- available_variables[!available_variables== this_variable_to_split]


  if (length( available_variables)==0){
    break
  }

} #end build list construction


colnames(nm_build_list) <- c("cut", "leaf1", "leaf2", "leaf1_mean", "leaf2_mean", "parent_node")


original_nm_build_list <- nm_build_list


#the non_model assignment covers EVERY case
y<- c(rep(1:max((n_nm_groups-1), 1),round(c.number/n_nm_groups)), rep(n_nm_groups, c.number-(n_nm_groups-1)*round(c.number/n_nm_groups)))

df$nm_model <- y[!y==0]

while(nrow(nm_build_list)>0){
  which_biggest_parent_node <- which(nchar(nm_build_list$parent_node)==max(nchar(nm_build_list$parent_node)))[1]
  this_row <- nm_build_list[which_biggest_parent_node,]
  this_parent_node <- nm_build_list$parent_node[which_biggest_parent_node]
  IV_column_index <-as.integer( this_row$cut)+1 #index of this IV
  leaf1 <- this_row$leaf1
  leaf2 <- this_row$leaf2
  l1_members <- as.integer(unlist(stringr::str_split(leaf1, "&")))
  l2_members <- as.integer(unlist(stringr::str_split(leaf2, "&")))
  l1_n <- nrow(df[df$nm_model%in% l1_members, IV_column_index, drop=F])
  l2_n<- nrow(df[df$nm_model%in% l2_members, IV_column_index, drop=F])

  this_variable_type <- variable_vector[as.integer( this_row$cut)]

  #check whether binary or continuous


  if (this_variable_type=="continuous"){


    df[df$nm_model %in% l1_members, IV_column_index] <-  rnorm(n=l1_n, mean = as.numeric( this_row$leaf1_mean), sd = .1)
    df[df$nm_model%in% l2_members, IV_column_index] <-  rnorm(n=l2_n, mean = as.numeric( this_row$leaf2_mean), sd = .1)

    #compress anything above 1 or under 0 that that range
    df[, IV_column_index][df[,IV_column_index]<0] <- 0
    df[, IV_column_index][df[,IV_column_index]>1] <- 1

  }else if (this_variable_type=="binary"){

    df[df$nm_model %in% l1_members, IV_column_index] <-  rbinom(n=l1_n, size=1, prob= as.numeric( this_row$leaf1_mean))
    df[df$nm_model%in% l2_members, IV_column_index] <-  rbinom(n=l2_n, size=1, prob = as.numeric( this_row$leaf2_mean))


  }

  nm_build_list <- nm_build_list[-which_biggest_parent_node,]

}
} #end if build_complete is at least 1

#end Second build list####


#remove variable identifying the non-model-significant subgroups
if ("nm_model" %in% colnames(df)){
df <- df[,-which(colnames(df)=="nm_model")]
}

#now check the NAs, the variables not yet defined by their subgroups, to see if they are accurate and fill them.
#cycle through each IV, give it a mean/prob near .5, and then randomly fill the NAs
IV_columns <- variables+1
for (i in 1:length(IV_columns)){
  this_iv_type <- variable_vector[i]

  if (this_iv_type=="continuous"){

    #first establish the mean of the NA cells, those not members of crucial groups. Should be near .5
   rest_mean <- ( rnorm(n=1, mean = .5, sd = .05))
   this_n <- length(     df[is.na(df[,IV_columns[i]]),IV_columns[i]])
   if (this_n>0){
    df[is.na(df[,IV_columns[i]]),IV_columns[i]] <- scales::rescale(to=c(0,1), rnorm(n=this_n, mean = rest_mean, sd = .1))
   }

  }else if (this_iv_type=="binary"){

    rest_prob <-   ( rnorm(n=1, mean = .5, sd = .05))
    this_n <- length(     df[is.na(df[,IV_columns[i]]),IV_columns[i]])
    if (this_n>0){
    df[is.na(df[,IV_columns[i]]),IV_columns[i]] <-  rbinom(n=this_n, size=1, prob= rest_prob)

    }

  }

}

 for (m in 1:m.number) {


    #now iterate through the build list in order to construct the separate models


    # Calculate the DV depending on the underlying randomly created models
    # first calculate the coefficients. They should not be massively different from one another? should tend toward zero
   # only X number of coefficients should matter

    n_substantial_IVs <- sample(variables, size=1, prob =1/(abs(variables-(mean(variables)))+1))
   which_substantial_IVs <- sample(variables, n_substantial_IVs)

   #so the substantial coefficients should have a bigger possible cofficient size.
   #the others shoudl be close to zero
    substantial_coef <- rnorm(n_substantial_IVs, mean=0, sd = .30)  # creating random coefficients for the models
    nonsubstantial_coef <- rnorm(length(variables)-n_substantial_IVs, mean=0, sd = .05)

    #now put the cofficients in the correct order
    beta_order <- rep(0, length(variables))
  beta_order[which_substantial_IVs] <- substantial_coef
  beta_order[-which_substantial_IVs]  <- nonsubstantial_coef

    value_coef_calculations <- data.frame(mapply(`*`,df[df$model==m,2:(v.number+1)],beta_order))
    intercept <- sample(0:100, size=1, prob =log(1/abs(sort(rnorm(101, mean=.5, sd=.20))-.50)))/100




    if(error == TRUE){
      error.values <- rnorm(m_sizes[m], mean = 0, sd = error.sd)
    }
    DV <- rowSums(value_coef_calculations) + error.values + intercept # calculating DVS
    df[df$model==m,1] <- DV

  }



  df <- df %>%
    rename("y" = "X1")%>%
    rename_at(vars(2:(ncol(.)-1)), ~ paste("X", 1:length(.), sep = ""))


  if(data.form == "normalized"){

    df[,1] <- scales::rescale(to=c(0,1), x= df[,1])
      }

    # Calculate the underlying models
  models <- data.frame(matrix(0, nrow = m.number, ncol = v.number + 1))
  for (n in 1:m.number){

    model <- lm(y ~ ., data = df[df$model==n,1:(v.number+1)])

    #if a model has NA, meaning it is a constant. replace NA with zero
    model$coefficients[is.na(model$coefficients)] <- 0

    models[n,] <- model$coefficients
  }

  models <- models %>%
    rename("Intercept" = "X1")%>%
    rename_at(vars(2:ncol(.)), ~ paste("b", 1:length(.), sep = ""))



  # Calculate the profiles of the different subgroups
  #Take care of this with the special function

dna_profiles <-    create_profile(data_after_dummies=df[,-ncol(df)], #remove model column here
                                             model_assignment_data=df,
                                             dummy_vars=NULL,
                                             model_labels = NULL)


  output <- list(df, models, dna_profiles, original_build_list)
  names(output) <- c("data", "models", "profiles", "build_list")



  return(output)
}











#' Create New Children
#'
#' Generate random population (coefficients) for the genetic algorithm ('DNA' or  'children')
#' @export
children_spawn <- function(n_children,
                           dna_length=ncol(newdf),
                           col_names=c('intercept',colnames(newdf)[2:length(colnames(newdf))]), #the DV name column is replaced with 'intercept'
                           normal=T, #normally distributed with the following SD. If not, then this_range is used
                           this_sd=.5, #standard deviation of the coefficients.. this shoudl be tailored to the data's form! .25 is a good default
                           this_range=c(-10, 10)){ #this might need to be adjusted. In practice most surviving coefficients are small. But i guess leaving this large gives a bigger chance for larger coefficients to have at least a chance.

  if (normal==T){
    #this makes the children normally distributed, thus extreme coefficients harder to get

    #the fake data use this  coef <- rnorm(v.number, mean=0, sd = .35)

     dna_table <- data.frame( matrix(  rnorm(dna_length*n_children,mean=0, sd=this_sd) , nrow=n_children, ncol=dna_length ))

    dna_table[dna_table>1] <- 1
    dna_table[dna_table<(-1)] <- -1

  }else{
    dna_table <- data.frame( matrix(  sample(this_range[1]:this_range[2], dna_length*n_children, replace=T), nrow=n_children, ncol=dna_length ))

  }

  #these columns are extra.. .just storing the characteristics of each genome
  dna_table$DNA_age <- 0
  dna_table$DNA_rating <- 0

  colnames(dna_table)[1:(ncol(dna_table)-2)] <- col_names

  op <- options(digits.secs = 6)
  options(op)
  IDs <-  as.integer(as.POSIXct(Sys.time()))
  IDs <- c(IDs:(IDs+(nrow(dna_table)-1)))


  rownames(dna_table) <- IDs

  #and now redfine intercept around the mean of the outcome.

  #intercept cannot possibly be other than between 0 and 1. Using the normal distribution around the grand mean before created impossible cases as well as restricted the opening pool of possibilities.

  dna_table[,1] <- sample(1:1000,n_children, replace = T)/1000

  return(dna_table)


}

#' Genetic Algorithm Training
#'
#' Run genetic algorithm and store average value characteristics for children###
#' This first checks the exact algo for a solution. Then loops through children to see if they can find it. Then it ranks the children according to the uniqueness of their solution.
#' @export
test_children <- function(dna_pool=children,
                          input_data= training_set,
                          closeness_threshold=solution_threshold, #defines how close a solution is to the found solution
                          generations=10,
                          difficulty=.8,
                          death_rating_thresh=.5,
                          death_number_thresh=NULL,
                          uniqueness_power=2,
                          number_of_children=n_children,
                          death_by_ageing=999999999,
                          min_strand_solution_proportion=.05,
                          mutations=mutation_rate, #proportion of genes that are mutated
                          mating = c("pam", "dam", "random"), # pam = positive assortative mating, dam = disassortative mating
                          mating_power = 1, # indicates the strength of the effect of closeness on mating
                          selection = c("threshold", "probability"), #note that threshold is always applied to the final generation results.
                          rem_doppleganger=T, #removes doppleganger dna, keeping the youngest
                          nelitism = 0 #defines the number of best models that survive every generation
                          ){

  ###


  ###
  input_data <- data.frame(input_data)


  #initial report about dna pool individual indicators
  dna_length <-  ncol(dna_pool)-2

  initial_column_means <- apply(dna_pool[,1:dna_length], 2, mean)




  intercept_coef <- rep(1, nrow(input_data))
  input_data <- add_column(input_data, intercept_coef, .after = 1)

  dna_evolution <- NULL # for visualization


  dependent_variable <- input_data[,1]
  input_data <- input_data[,-1]


  generation <- 0

  #need raw DNA strings without the ending two columns

  for (g in 1:generations){




    highest_solved <- 0
    generation <- generation+1

    #Main tables###
    #first is the dna_pool, that stores a strands history and score
    #that is cumulative over a strand's history across generations
    #second is the table that will hold the relation between cases(rows) and the dna strands that solve them (columns)
    #that is only for the current generation: solution_table
    #third is the table that holds the profile of the solution set for a particular strand
    #that is also only for the current generation
    #that holds the characteristics of the values of the cases in the solution set that this dna_strand solved




    if (g==1){ #if the first generation, run through iteratively, because the matrices will be too big. Otherwise, use matrix operations
     cat("\n Running first generation. This may take a while. Further generations will be faster.")
        solution_table <- matrix(0, nrow=nrow(input_data), ncol=nrow(dna_pool))
#
        #put aside dna-string traickig information for the dna pool
        dna_tracking <- dna_pool[,(ncol(dna_pool)-1):ncol(dna_pool)]
        dna_pool <- dna_pool[,1:(ncol(dna_pool)-2)]

      for (dna in 1:nrow(dna_pool)){
        cat("Processing child", dna, "of", nrow(dna_pool),"\n")

        this_dna_strand <- dna_pool[dna,]
        #now run that dna against every member of the input_data, collect residuals

        x <-sweep(input_data, MARGIN=2,unlist( this_dna_strand), `*`) #input_data without dv column

        #now take row sums of x, and take their difference from happiness
        these_predictions <- rowSums(x)
        gap_from_Y <-  (these_predictions- dependent_variable)
        #smallest gap
        smallest_gap <- min(abs( gap_from_Y))
        best <- which( abs( gap_from_Y)==    smallest_gap)
        solution_set <-   which( abs( gap_from_Y) <= closeness_threshold) #indices, not case IDs!
       #if the solution set is not consisting of at least two points, set it to integer(0)
        if (length(solution_set)<2){
          solution_set <- vector("integer",0)
        }
         attr(solution_set, 'names') <- NULL

        #store the solutions (cases as rows) for this particular dna strand (columns)
        solution_table[solution_set,dna] <- 1



      } #end DNA loop

    }else{ #for every generation but the first





    input_data <- data.matrix(input_data)

        dna_pool <- data.matrix(dna_pool)

    #the DNA has already been trimmed!!

    #put aside dna-string traickig information for the dna pool
   dna_tracking <- dna_pool[,(ncol(dna_pool)-1):ncol(dna_pool)]
     dna_pool <- dna_pool[,1:(ncol(dna_pool)-2)]

    #now create a matrix where each row of input data is repeated 50 times (nrow(dnapool) and each row of dna is repeated in sequence once for each row of data

    large_data <- input_data[rep(1:nrow(input_data), each=nrow(dna_pool)), ]
   large_dependent_variable <- dependent_variable[rep(1:nrow(input_data), each=nrow(dna_pool))]
     large_dna_pool <- dna_pool[rep(1:nrow(dna_pool), times=nrow(input_data)),]
    large_product <- large_data*large_dna_pool
    large_dv_predictions <- rowSums(large_product)
    gap_from_real_dv <- abs(large_dv_predictions-large_dependent_variable)
    solved_or_not <- gap_from_real_dv<= closeness_threshold

    #now fill in solution table and then get dna profiles.

    large_solution_set <-   which( solved_or_not==T ) #indices of the LARGE matrix, not case IDs!
    attr(large_solution_set, "names") <- NULL
    large_solution_set_plotting <- data.frame(index=large_solution_set, row=ceiling(large_solution_set/nrow(dna_pool)),col= large_solution_set-(ceiling(large_solution_set/nrow(dna_pool))-1)*nrow(dna_pool) )


    #calculate overall indices for solution table based on above plotting
    #the indices in the matrix function fill by column.. then next column, then next
    large_solution_set_plotting$new_matrix_index <- (large_solution_set_plotting$col-1)*nrow(input_data)+large_solution_set_plotting$row


    #store the solutions (cases as rows, dna strands as columns)
    solution_table <- matrix(0, nrow=nrow(input_data), ncol=nrow(dna_pool))
    small_solution_set <- large_solution_set_plotting$new_matrix_index

    if (length(small_solution_set)<2){
      small_solution_set <- vector("integer",0)
    }

    solution_table[small_solution_set] <- 1






    }



    #Adding the IDs from the dna pool

    colnames(solution_table) <- rownames(dna_pool)
    rownames(solution_table) <- rownames(input_data)

    dna_pool <- cbind(dna_pool, dna_tracking)
    #update DNA_rating and DNA_age for each dna strand###
    dna_pool <- data.frame(dna_pool)
    dna_pool$DNA_age <- dna_pool$DNA_age +1


    coverage <- colSums(solution_table) #Coverage of each Model
    min_solved <- nrow(input_data)*min_strand_solution_proportion
    #now set to zero any columns below the threshold
    which_to_set_zero <- which (coverage<min_solved)
    solution_table[,which_to_set_zero] <- 0

    uniqueness <-   rowSums(solution_table) #should be low. how often each case was solved
    nonzeroes <- which(!uniqueness==0)
    #the line below inverts uniqueness so that a higher value means more unique and puts it on a 0 to 1 scale
    rescale_nonzeroes <- 1-scales::rescale(uniqueness[nonzeroes],to=c(0,1),  from = c(1, max(uniqueness)))
    uniqueness[nonzeroes] <- rescale_nonzeroes



    #penalize those with low uniquess!
    uniqueness <- uniqueness^uniqueness_power

    #now lets try simple weightings
    #so each row should be multiplied by its case weights.
    #and then the columns can be summed to create DNA_rating


    new_solution_table <- uniqueness*solution_table




       rating <- colSums(   new_solution_table) #should be HIGH. how successful each strand was overall

   #leaving rating in its raw form.




    dna_pool$DNA_rating <- rating



    #take the indices of the saved DNA and save the dna_profiles that match them!

    #save only those that solved the threshold amount
    which_saved <-which(  dna_pool$DNA_rating>death_rating_thresh)



    if (length(which_saved)==0){
      cat("\nNot enough children. Natural selection over.\n ")
      dna_profiles <- NA

      return(NA)
    }

    dna_pool <- dna_pool[which_saved,]

    solution_table <- as.data.frame(  solution_table[,which_saved])



    matching <- order(dna_pool$DNA_rating, decreasing = T)

    dna_pool <- dna_pool[matching,] #dna profiles and pool adapted so that the order is the same, matching ids

    solution_table <-  as.data.frame( solution_table[,matching])


    colnames(solution_table) <- rownames(dna_pool)
    rownames(solution_table) <- rownames(input_data)


    # Elitism

    which_elit <- NULL
    solution_table_elit <- solution_table
    nelitism <- min(nelitism,  ncol(solution_table))


if (nrow(dna_pool)>1){
    for(e in 1:nelitism){
  #same elite should only be selected once!!

      this_elit <- attr(sort(colSums(solution_table_elit[,1:length(solution_table_elit)]),decreasing=TRUE)[1], "names")
      #remove all rows solved by this elite
      these_rows <-       which(colnames(solution_table_elit)== this_elit)


      solution_table_elit <- solution_table_elit[solution_table_elit[,these_rows] != 1,]

      #populate this with model names only!! not indices
       which_elit[e] <- this_elit


    }


    # DNA_pool without the ones chosen as elite
    dna_pool_noelit <- dna_pool[!rownames(dna_pool)%in% which_elit,]
}else{
  dna_pool_noelit <- dna_pool
}

    #save models based on threshold (only the top x number)

    if(selection == "threshold" | g==generations){
      if (!is.null(death_number_thresh)){
        dna_pool_noelit <- dna_pool_noelit[1:min(nrow(dna_pool_noelit), (death_number_thresh-nelitism)),]

        model.names <- rownames(dna_pool_noelit)

        #which_selected <- which(rownames(dna_pool) %in% model.names)
        which_selected <- c(which_elit, model.names)

      }
    }else if (selection== "probability") { #save models based on of their rating as a probability
#if not the last generation
      if (!g==generations & nrow(dna_pool_noelit)>0){
        if (!is.null(death_number_thresh)){
          which_selected_noelit <- sample(nrow(dna_pool_noelit), min(nrow(dna_pool_noelit), (death_number_thresh-nelitism)), prob = dna_pool_noelit$DNA_rating)
          model.names <- rownames(dna_pool_noelit[which_selected_noelit,])
          #indices selected
          these_rows <- which(rownames(dna_pool) %in% model.names)
          #elites added back in
          which_selected <- c(which_elit, model.names)


        }

      }else{ #if the last generation

        if (!is.null(death_number_thresh)){

          which_selected_noelit <- vector("integer", 0)
          #indices selected
          which_selected <- which_selected_noelit
          #elites added back in
          which_selected <- c(which_elit, which_selected)


        }

      }


    } #end probability selection





    # Update dna_pool, dna_profiles and solution_table
    dna_pool <- dna_pool[rownames(dna_pool) %in% which_selected,,drop=F]
    solution_table <- as.data.frame(solution_table[,colnames(solution_table) %in% which_selected, drop=F])

    #we create the dna profiles here.

# we start with the first column and then loop through the others

    solution_table$model <- solution_table[,1]

    building_dna_profiles <-    create_profile(data_after_dummies=df_full_dummies[training_indices,],
                                               model_assignment_data=solution_table,
                                               dummy_vars=original_variables[dummy_sets],
                                               model_labels = colnames(solution_table)[1])



#filling in remaining models
if (ncol(solution_table)>2){
for (i in 2:(ncol(solution_table)-1)){ #minus one because there is a model column for this solution table


   solution_table$model <- solution_table[,i]

  adding_dna_profiles <-    create_profile(data_after_dummies=df_full_dummies[training_indices,],
                                             model_assignment_data=solution_table,
                                             dummy_vars=original_variables[dummy_sets],
                                             model_labels = colnames(solution_table)[i])


  building_dna_profiles <- rbind(building_dna_profiles, adding_dna_profiles)

}



  #the cases have overlapping profiles at this point because multiple models solve multiple cases.
  #so the dnaprofiles indeed need to be built column by column

  #this leaves some cases without a model of course



  }else{ #if there is only one model, assign it to all cases
  solution_table$model <- 1


  building_dna_profiles <-    create_profile(data_after_dummies=df_full_dummies[training_indices,],
                                             model_assignment_data=solution_table,
                                             dummy_vars=original_variables[dummy_sets],
                                             model_labels = colnames(solution_table)[1])
}

dna_profiles <- building_dna_profiles


#and then remove the model column here because it does not make sense because models overlap
solution_table <- solution_table[,-which(colnames(solution_table)== "model"), drop=F]
    # Only mate if not the last generation
    if (!g==generations){
    # If positive assortative mating is used:
    # aim is to make it more likely for models to mate with other models that are close to them
    # in order to enable the simultaneous exploration of different local maxima
    # First: Create a closeness table based on the sum of the standard deviation differences
    #         for all coefficients between two models
    if (mating == "pam" || mating == "dam"){
      closeness_table <- matrix(0, nrow = nrow(dna_pool), ncol = nrow(dna_pool))
      for(row in 1:nrow(closeness_table)) {
        #! Vectorize this to speed up####
        sd <- sapply(dna_pool[, 1:(ncol(dna_pool)-2)], sd) + 1
        for(col in 1:ncol(closeness_table)) {
            closeness_table[row,col] <- sum(abs((dna_pool[row,1:(ncol(dna_pool)-2)] - dna_pool[col,1:(ncol(dna_pool)-2)]) / sd))
        }
      }

      if(mating == "pam"){
        closeness_table <- max(closeness_table)-closeness_table  # reverse so that higher values indicate more closeness
      }

      if(sum(closeness_table) == 0){
        closeness_table <- closeness_table + 1 # prevent dividing by zero and a probability vector of only zeros in the following
      }

      diag(closeness_table) <- 0    #set closeness between the model and itself to zero
      RS <- rowSums(closeness_table)
        closeness_table <- closeness_table / RS # transform so that sum is 1

       closeness_table <- closeness_table^mating_power
    }


    #save that generation###
    #    write.csv(dna_pool, row.names = F, paste0(getwd(), "/QUEENS/configurational regression/",format(Sys.time(), "%d-%b-%Y %H_%M_%OS4"),"DNA generation ", g,".csv"))
    #   write.csv(dna_profiles, row.names = F, paste0(getwd(),"/QUEENS/configurational regression/",format(Sys.time(), "%d-%b-%Y %H_%M_%OS4"),"DNA_profiles generation ", g,".csv"))
    #  write.csv(solution_table, row.names = F, paste0(getwd(),"/QUEENS/configurational regression/",format(Sys.time(), "%d-%b-%Y %H_%M_%OS4"),"solution_table generation ", g,".csv"))

    #now let them mutate and mate
    if (mating == "pam" || mating == "dam"){
      # New parents table based on closeness
      parents_table <- matrix(0, nrow = floor(nrow(dna_pool)/2), ncol = 2)
      nrow_dna <- 1:nrow(dna_pool)
      nrow_dna_rm <- nrow_dna
      # creating a parents_table with two columns for first and second
      # parent. First parent is randomly chosen among the remaining
      # models from the dna.pool. The second parent is chosen with a prob
      # parameter that takes into account the closeness to the first
      # parent. However this means in the second column the same parent
      # could potentially be chosen several times.

       for(parents in 1:floor(nrow(dna_pool)/2)){
        parents_table[parents,1] <- sample(nrow_dna_rm,1)
        parents_table[parents,2] <- sample(nrow_dna,1, prob = closeness_table[parents_table[parents,1],nrow_dna])
        nrow_dna_rm <- nrow_dna[-c(parents_table)]}
      first_parents <- parents_table[,1]
      second_parents <- parents_table[,2]
    }else{
      first_parents <- sample(1:nrow(dna_pool), floor(nrow(dna_pool)/2) )#indices of first parents
      second_parents <- (1:nrow(dna_pool))[-first_parents][1:length(first_parents)]
    }



    if (length(first_parents)>0 & !is.na(second_parents[1])){

      for (child in 1:number_of_children){
        dna_pool <-   parents_mate( dna_set=dna_pool,
                                    first_parent_indices=first_parents,
                                    second_parent_indices=second_parents,
                                    mutation_rate=mutations,
                                    random_mate=T)
      }

    }else{
      cat("\nNot enough children. Natural selection over.\n ")
      break
    }


    #sort this by rating AND by age. so that highest rating, highest age is on top.

       #death by ageing happens now... before doppleganger removal, so that the lone parent doesn't die by age after the doppleganger is removed

    which_saved <-which(  dna_pool$DNA_age < death_by_ageing)
    if (length(which_saved)==0){
      cat("\nNot enough children. Natural selection over.\n ")
      break
    }
    dna_pool <- dna_pool[which_saved,]
    } #ends if statement for implementing mating if not last round.

    if (rem_doppleganger==T){ #remove dopplegangers



      original_dna <- rownames(dna_pool)
      which_unique <- rownames(unique(dna_pool[,1:(ncol(dna_pool)-2)])  )
      which_saved <- which(original_dna %in% which_unique)

      dna_pool <- dna_pool[which_saved,]

      if (length(which_saved)<=dim(dna_profiles)[2]){
            dna_profiles <- dna_profiles[which_saved,]
      }



      #the solution table is wrong!!
      if (length(which_saved)<=dim(solution_table)[2]){
      solution_table <- solution_table[,which_saved, drop=F]
      }


    }

    #just use the trimmed down solution table!#


#   Ecological solutions below are the total solved by ANY model.
solved <- round(length( which (rowSums(  solution_table)>0))/nrow(solution_table)*100, digits = 2)

#below are the model_specific solutions
 #   these_dna_active <- dna_pool[ rownames(dna_pool) %in%  rownames(dna_profiles),]
  #  solved <- predict_new(models=these_dna_active,profiles = dna_profiles , new_data = training_set)

    cat("\n\n Generation", generation, "complete. Solved %:", solved, "of",nrow(input_data), "cases, with", ncol(solution_table), "DNA strands\n\n")

    # only for the purpose of visualization
    dna_pool_gen <-cbind(dna_pool, generation = rep(g, nrow(dna_pool)))
    dna_evolution <- rbind(dna_evolution, dna_pool_gen)
  }#end generation loop





#add R-squared solution

training_perc_explained <- predict_new(models=dna_pool,
                                       profiles=dna_profiles,
                                       new_data = input_data,
                                       assign_models = T,
                                       method = "residuals")



  # dna_evolution is only for visualization!! Deprecated.
  return(list(dna_pool=dna_pool, dna_profiles=dna_profiles, solution_table=solution_table,  indicator_report=compare_means, dna_evolution = dna_evolution, discrete_percent_solved=training_perc_explained$discrete_percent_solved, r_squared=training_perc_explained$r_squared))

}

#' Genetic Crossover
#'
#' Creates new children based on the current successful dnapool
#' @export
parents_mate <- function(dna_set=dna_pool,
                         first_parent_indices,
                         second_parent_indices,
                         mutation_rate,
                         random_mate=T,
                         mutation_granularity=1, #higher number means smaller shifts in original values
                         percent_explained=solved,
                         mode="continuous" #or "binary"
                         ){
  #now join each set of parents, mutate the results, and rbind these children to the old dna_pool
  #cut off unneeded columns
  #higher mutation granularity means mutations are closer to the original

  if (random_mate==T){
    second_parent_indices <- sample(second_parent_indices)
  }

  dna_length <- ncol(dna_set)-2
  dna_set <- dna_set[,1:(dna_length+2)]

  for (new_child in 1:length(first_parent_indices)){

    first_dna_division <- sample(1:dna_length, sample(1:dna_length,1), replace = F) #this assigns 0's to random parent genes, so that I can add the genes together, leaving only one parent relevant
    second_dna_division <- (1:dna_length)[-first_dna_division]

    first_parent_contribution <-  dna_set[first_parent_indices[new_child],1:dna_length]
    first_parent_contribution[,second_dna_division] <- 0
    second_parent_contribution <-  dna_set[second_parent_indices[new_child],1:dna_length]

    second_parent_contribution[,first_dna_division] <- 0

    this_child <- first_parent_contribution+second_parent_contribution


    mutations <- (sample(1:100,dna_length)/100) <= mutation_rate #smaller than or equal to mutation rate.

    mutations <- c(mutations, F, F)

    mutation_number <- length(mutations[mutations==T])

    if (mode=="continuous"){

      #This should be adjusted to adjust around the origial value based on a normal distribution.
    #take the standard deviations of all coefficients that will be mutated from the whole pool
    sd_list <-     sapply(as.list(data.frame(dna_set[,mutations==T])), sd)
    #take those original values
    original_coeff <-     this_child[mutations==T]

    mutation_chain <- mapply (coef=unname(unlist(original_coeff)), sds= sd_list, FUN=function(coef, sds)     rnorm(n=1, mean=coef, sd=sds/mutation_granularity))
    }else if (mode=="binary"){
      #now flip the bit of the mutations.

      mutation_chain <- 1-    this_child[mutations==T]
      #make zeroes or ones.
    }

    this_child[,mutations==T] <- mutation_chain
    this_child$DNA_age <- 0
    this_child$DNA_rating <- 0

    dna_set <- rbind(dna_set, this_child)
  }

  op <- options(digits.secs = 6)
  options(op)
  IDs <-  as.integer(as.POSIXct(Sys.time()))
  IDs <- c(IDs:(IDs+(nrow(dna_set)-1)))

  rownames(dna_set) <- IDs

  return(new_dna_pool=dna_set)

}

#' Create Profile
#'
#'Create a descriptive profile of the subpopulations for each model, including all categories (reference categories are not left out)
#' @export
create_profile <- function(data_after_dummies=df,
                           model_assignment_data=NULL,#must have rownames for each case and a column called model. If NULL, it assumes one model
                           dummy_vars="cntry",#a vector of the categorical variables (3+ categories) for which dummies should be made
                           model_labels=NULL #whether model labels are specified or automatically assigned
                           ){

  #input... a. data in original form after split into dummies with relevant variables, rownames must match cases
  #         b. data that shows model assignments and rownames
  #output. dna profiles for each model


  if (length(dummy_vars)>0){
  save_rownames <- rownames(data_after_dummies)

  working_df <- data_after_dummies

  for (d in 1:length(dummy_vars)){
this_dummy <- dummy_vars[d]


  if (this_dummy %in% colnames(working_df)){

#in every case, get FULL dummy columns for the profiles, don't leave ref categories out.
  working_df <- dummy_cols(working_df, select_columns = this_dummy, remove_first_dummy = F, ignore_na = T, remove_selected_columns = T)
  rownames(working_df) <- save_rownames
  }else{

  }

  }#end for loop for each dummy variable
  dummy_df <- working_df

}else{
  dummy_df <- data_after_dummies
}


  if (is.null( model_assignment_data)){
    case_model_vector <- 1
    models <- 1
  }else{
    #how to handle cases where only one model with one case?

   models <-  unique(as.integer( model_assignment_data$model))
   models <- models[!models==0] #zero model is always removed


   case_model_vector <- as.integer( model_assignment_data$model)
  }

  if (is.null(model_labels)){
    model_labels <- as.character(models)
  }

  if ("model" %in% model_labels){
    model_labels <- model_labels[-which(model_labels=="model")]
  }

  if ("model" %in% colnames(dummy_df)){
    dummy_df <- dummy_df[,-which(colnames(dummy_df)=="model")]
  }

  dna_profiles <- data.frame(matrix(0, nrow=length(models), ncol=ncol(dummy_df)))
  colnames(dna_profiles) <-paste( colnames(dummy_df), "mean/std")
  rownames(dna_profiles) <- as.character( models)




  for (m in 1:length(models)){

  #create profiles for the current solution set! for that dna strand

  dna_profile_mean <- round(colMeans(dummy_df[case_model_vector==models[m],, drop=F ]) , digits=3)



  if (nrow(dummy_df[case_model_vector==models[m],,drop=F])>0){
    dna_profile_sd <- round(matrixStats::colSds(as.matrix(dummy_df[case_model_vector==models[m], ,drop=F])), digits=3)
  }else{
    dna_profile_sd <- NA
  }
  new_profile <- mapply(x=dna_profile_mean, y=dna_profile_sd, function(x,y) paste(x,"/", y))
  dna_profiles[m,] <- new_profile
   }



if (!nrow(dna_profiles)==length(model_labels)){

}


rownames(dna_profiles) <- as.character(model_labels)


    return(dna_profiles)
}

#' Deprecated?
#'
#' out of date.. Potentially use the proportional explanation table as a useful statistic
#' @export
interpret_results <- function(dna=dna_pool, profiles=dna_profiles, solution_table=solution_table){

  #cut out 0 rows off dna
  dna <- dna[!dna$DNA_age==0,]
  #cut away unneeded columns
  dna_condensed <- dna[,-c((ncol(dna)-1):ncol(dna))]
  profiles_condensed <- profiles[,-1]

  x <- dna_condensed*profiles_condensed


  absolute_row_sums <-   rowSums(  abs(x))
  proportional_explanation_table <-  round( abs(x)/absolute_row_sums, digits = 2)
}







#' Weighted cluster means
#'
#' Weighted cluster means
#' @export
  cluster_means <- function(data){
    means <- as.data.frame(matrix(0, ncol = ncol(data)-3, nrow = max(data$cluster)))
    for(c in 1:max(data$cluster)){
      this_data <- data[data$cluster == c,]
      mean.cluster <- lapply(this_data[ , 1:(ncol(this_data)-3)], weighted.mean,  w = this_data$DNA_rating)
      means[c,] <- mean.cluster
    }
    return(means)
  }

#' Condenses models.
#'
#' @export
  condense_models <- function(data = NULL,
                              training_set = training_set,
                              method = c("best", "cluster", "genetic"), #genetic is a bit broken... the cross problem.
                              nmodel = 10, #crucial parameter for interpretation. probably best to try out a few
                              auto_cluster=TRUE,
                              closeness_threshold = solution_threshold,
                              solution_table, #needed for genetic method
                              no_correlation_thresh=-.1, #for genetic algo, a correlation between binary strings in the solution table
                              pearson_corr_threshold=.7 #this compares the model cofficients directly and keeps the best, removing others over the threshold.
                              ){




    #make sure that the model rownames are the same in the data and the models themselves

    dna_pool <- data$dna_pool
    training_set <- as.data.frame(training_set)

    #removed correlated_models

    if ("DNA_rating" %in% colnames(dna_pool)){
    dna_length <- ncol(dna_pool)-2
    }else{
      dna_length <- ncol(dna_pool)
    }

    #make sure dna are sorted.

    dna_pool <- dna_pool[order(dna_pool$DNA_rating, decreasing = T),]
  dnapool_correlations <-   cor(t(dna_pool[,1:dna_length]))
  illegal_correlations <- dnapool_correlations>pearson_corr_threshold

  #now we loop through columns keeping the column name (which is the best model), and removing all illegals

  while (ncol(illegal_correlations)>0){
    this_column <- illegal_correlations[,1, drop=F]
    this_model <- colnames(this_column)


    models_to_remove <-  rownames(this_column)[ this_column==T]
    #make sure we don't remove the original comparison model
models_to_remove <- models_to_remove[!models_to_remove==this_model]

if (length(models_to_remove)>0){
  #remove those models from the dna pool as well as the correlations
  #this loop stops when the correlation table becomes too small

  dna_pool <- dna_pool[!rownames(dna_pool) %in% models_to_remove,]
illegal_correlations <- illegal_correlations[!rownames(illegal_correlations)%in% models_to_remove,!colnames(illegal_correlations)%in% models_to_remove, drop=F ]


}
#also shave off first column and row, which was just checked


illegal_correlations <- illegal_correlations[-1,-1, drop=F]

}


  #cut solution table to match the dna_pool

  these_models_alive <- rownames(dna_pool)

  solution_table <- solution_table[,colnames(solution_table) %in% these_models_alive, drop=F]
 #cut out the models here too?
  #where does 'models' come from?

   if(method == "best"){


     nmodel <- min(as.integer(nmodel), nrow(dna_pool))
      which_best <- NULL

      for(e in 1:nmodel){

  if(nrow(solution_table)==0){
  break
}
        this_best <- attr(sort(colSums(solution_table),decreasing=TRUE)[1], "names")
        this_best_index <- which(colnames(solution_table)== this_best)
        solution_table <- solution_table[solution_table[,this_best_index] != 1,, drop=F]
        which_best[e] <- this_best_index
        models <- dna_pool[which_best,] #this does not save the models results
      }

         #we need to now assign all cases to these new clusters
      #and then take profiles from those cases

      elitist_model_assignments <-       model_closeness(data = training_set,
                                                                   profiles = NULL, #we don't have profiles yet because cases not yet assigned to these new models
                                                                   models = models,
                                                                   method = "residuals") #'profiles' or 'residuals'


       dna_profiles <-  create_profile(data_after_dummies=df_full_dummies[training_indices,],
                                     model_assignment_data=elitist_model_assignments,
                                     dummy_vars=original_variables[dummy_sets],
                                     model_labels = NULL)


model_assignments <- elitist_model_assignments

      } else if(method == "cluster"  ){





        nmodel <- nrow(dna_pool)
        if ( nrow(dna_pool)>2){

      # Find the variables that best divide the clusters

        if (nrow(dna_pool)>1){
      coef <- colSums(abs(dna_pool[,1:(ncol(dna_pool)-2), drop=F])) / matrixStats::colSds(as.matrix(dna_pool[,1:(ncol(dna_pool)-2)]))
      #remove any infinities
      coef <- coef[!coef==Inf]

      max.coef <- sort(coef, decreasing = T)[c(1,2)]

      # Find the best number of cluster
      if(auto_cluster == TRUE & nrow(dna_pool)>3){
        unique_data <- nrow(unique(select(dna_pool, names(max.coef))))

            gap <-cluster::clusGap(select(dna_pool, names(max.coef)), FUN=kmeans, K.max=max(min(10, unique_data-1), 2), B=500) # based on gap statistics
        gap <- as.data.frame(gap$Tab)




       #just take the highest gap statistics such that gap (This) is bigger than or equal to gap (next) minus gap (next) standard error
       ncluster <-  cluster::maxSE(gap$gap, gap$SE.sim, method = "Tibs2001SEmax")


      }else{ #no autocluster
        unique_data <- nrow(unique(select(dna_pool, names(max.coef))))


        ncluster = min(unique_data, nmodel-1)
      }


      #calculate the clusters for all models
      k <- kmeans(select(dna_pool, names(max.coef)), ncluster)
      #add them to the dna_pool
      dna_pool$cluster <- k$cluster

      #calculate the cluster means with the above defined function
     #this is the moment when the models are renamed
       models <- cluster_means(dna_pool)
      colnames(models) <- colnames(dna_pool[,1:dna_length])
     #cut off cluster column for pool, no longer needed
      dna_pool <- dna_pool[,-ncol(dna_pool)]
      #cut off all rows and replace with the new model row
      dna_pool <-  dna_pool[0,]

      dna_pool[1:nrow(models),1:dna_length] <- models
      dna_pool$DNA_age <- 0
      dna_pool$DNA_rating <- 0
        }

      #leave data intact... we need a real DV for the model closeness function


      # we have new models now... we could assign those models to the old cases via residuals
      # and then we could get their profile.

  cluster_condensed_model_assignments <-       model_closeness(data = training_set,
                                  profiles = NULL, #we don't have profiles yet because cases not yet assigned to these new models
                                  models = models,
                                  method = "residuals") #'profiles' or 'residuals'


      dna_profiles <- create_profile(data_after_dummies=df_full_dummies[training_indices,],
                                     model_assignment_data=cluster_condensed_model_assignments,
                                     dummy_vars=original_variables[dummy_sets],
                                     model_labels = NULL)

      model_assignments <- cluster_condensed_model_assignments
        }else{ #if too few models for cluster analysis, still need the model assignments and profiles

          models <- dna_pool

          cluster_condensed_model_assignments <-       model_closeness(data = training_set,
                                                                       profiles = NULL, #we don't have profiles yet because cases not yet assigned to these new models
                                                                       models = models,
                                                                       method = "residuals") #'profiles' or 'residuals'

          dna_profiles <- create_profile(data_after_dummies=df_full_dummies[training_indices,],
                                         model_assignment_data=cluster_condensed_model_assignments,
                                         dummy_vars=original_variables[dummy_sets],
                                         model_labels = NULL)

          model_assignments <- cluster_condensed_model_assignments


        }

      } else if (method=="genetic"){


        genetic_condensation_results <- genetic_condense(model_dna_pool = dna_pool,
                                                         solution_table = solution_table,
                                                         no_correlation_thresh=no_correlation_thresh)


          }

  if (!exists("model_assignments")){
    #if there are no model assignments, create them based on input data

    these_models <- rownames(dna_pool)
      model_assignments <-       model_closeness(data = training_set,
                                                               profiles = data$dna_profiles[rownames(data$dna_profiles) %in% these_models,], #we don't have profiles yet because cases not yet assigned to these new models
                                                               models = dna_pool,
                                                               method = "profiles") #'profiles' or 'residuals'
  models <- dna_pool
  }

#make sure models, profiles, assignments all have same model name

if (nrow(dna_pool)==1){
  if ((rownames(dna_pool)!=rownames(dna_profiles))){
  rownames(dna_profiles) <- rownames(dna_pool)
  model_assignments$model <- rownames(dna_pool)
}
}
    return(list(data_with_models=model_assignments, models = models, dna_profiles = dna_profiles))
  }

#' Gather R-squared
#'
#' Used to gather multiple R-squared values
#' @export
  regression_R_squared <- function(data_with_models){

    model_vector <- unique(data_with_models$model)
    dv <- colnames(data_with_models)[1]
    ivs <- colnames(data_with_models)[2:(ncol(data_with_models)-1)]
    ivs <- ivs[!ivs=="model"]

    results <- vector("list", length(model_vector))

    for (i in 1:length(model_vector)){

      this_sample <- data_with_models[data_with_models$model==model_vector[i],-ncol(data_with_models)]
      #are any of the IV's without any variation?
      sapply(lapply(data_with_models, unique), length)


      f <- as.formula(
        paste(dv,
              paste(ivs, collapse = " + "),
              sep = " ~ "))

      this_regression <- eval(bquote(   lm(.(f), data = this_sample)   ))
      x <- summary(this_regression)

      results[i] <- x$adj.r.squared


    }
    names(results) <- as.character(model_vector)
    #convert all negative numbers or NaN to zero
    convert_to_zero <- lapply(results, FUN=function(x) is.na(x)|x<0)
    results[unlist(convert_to_zero)] <- 0
    return(results)
  }




#' Predict new
#'
#' Use ICR results to predict never-seen case dependent variables
#' @export
predict_new <- function(models=result_training$dna_pool,
                        profiles=result_training$dna_profiles,
                        new_data=testing_set, #if assign_models is F, then there needs to be a model column here
                        closeness_threshold=solution_threshold, #set same as in training set!
                        assign_models=T, #if true, this assigns models based on the method below. #if we have assign on for only one model, it automatically assigns that model across the cases
                        method="case-based", #'residuals' when profiles are not known. 'profiles' when they are known, or 'case-based', when case rules are already created,
                        cutlist=NA #if case-based method, then a cutlist must be provided
                        ){

#In this function new data rows are tested against their assigned models
imported_data <- new_data
  #Assess the size of the input data

  if (ncol(models)>ncol(profiles)){ #resize models only if it does not match profiles
  models <- models[,1:(ncol(models)-2)]
  }
  #this_row <- new_data[1,]*models[1,]

  if (ncol(profiles)>ncol(models)){
  profiles <- profiles[,2:ncol(profiles)] # we dont want the DV or intercept profiles
  }

  if (assign_models==T){

    #use model_closeness function with an option for either residuals or profiles
   if (nrow(models)==1){ #if there is just one model, simply assign it
     data_with_models <- imported_data
     rownames(models) <- "1"
     data_with_models$model <- rownames(models)
   }else{

     data_with_models <-     model_closeness(data = imported_data,
                    profiles = profiles, #we don't have profiles yet because cases not yet assigned to these new models
                    models = models,
                    method = method,
                    cutlist=cutlist)

}

  }else{ #if models are already assigned
    data_with_models <- imported_data
   model_assignment <- imported_data$model
  #now convert that model assignment to indices to match the rownames of the models
  these_models <- unique(model_assignment)
  #matchthe model assignment names the appropriate row numbers of the model table
  these_model_indices <- match(these_models, rownames(models))


  #now recode the model assignment to those indices!

  model_assignment <- these_model_indices[ match(model_assignment,these_models) ]
  data_with_models$model <- model_assignment
  rownames(models) <- as.character(1:nrow(models))




  }


#predictions



which_model_column <- which(colnames(data_with_models)=="model")
data_with_models <- data_with_models[,c((1:ncol(data_with_models))[-which_model_column], which_model_column)]
#now apply the above dna models to the data

model_assignment <-  data_with_models$model
# for that, we can construct a matrix of the models in the correct order
# and then overlay that against the original data, multiplying them



model_assignment <- as.integer( gsub(" ", "", model_assignment))
#use indices so that they can repeat

relevant_model_matrix <- models[match(model_assignment,rownames(models)  ) ,]



dependent_variable <- imported_data[,1]

#remove dv from original data
new_data <- new_data[,-1]
#add 1 vector for intercept
new_data_matrix <- cbind(1, new_data)

if ("model" %in% colnames(new_data_matrix)){
  new_data_matrix <- new_data_matrix[,-(ncol(new_data_matrix))]
}

#now multiply


product_matrix <- new_data_matrix*relevant_model_matrix
actual_predictions <- rowSums(product_matrix)

predictions_compared_to_real <- abs(actual_predictions- dependent_variable)

solved_or_not <- predictions_compared_to_real<=closeness_threshold



solved_or_not <- data.frame(table(solved_or_not))





options(warn=-1) #OFF
r_squared <- regression_R_squared(data=data_with_models)
options(warn=0) #warnings on
case_n <- data.frame(table(data_with_models$model))
row_order <- match(names(r_squared), case_n$Var1)
case_n <- case_n[row_order,]
weighted_r_squared <- sum(case_n$Freq/sum(case_n$Freq)*unlist(r_squared))



discrete_solution <-paste0(  unlist(unname(round(solved_or_not[solved_or_not$solved_or_not==T,2]/(sum(solved_or_not$Freq))*100, digits = 2) )  ))
r_squared_solution <- round(weighted_r_squared, digits = 2)


if  ( length(discrete_solution)==0){

  discrete_solution <- 0
}

                                   #cat("\n",unlist(round(solved_or_not[2,2]/solved_or_not[1,2]*100, digits = 2)), "% of set solved")
return(list(discrete_percent_solved=discrete_solution, r_squared=r_squared_solution, model_closeness_data=data_with_models))




}





#' Model closeness
#'
#' Assigns models to cases. (e.g. for predictions for new cases)
#' @export
model_closeness <- function(data = training_set, #Get input_data with a column for the assigned modeled based on the dna_profiles or the residuals
                            profiles = dna_profiles,
                            models = models,
                            method = NULL, #'residuals' when profiles are not known. 'profiles' when they are known, or 'case-based', when case rules are already created,
                            cutlist=cutlist
                            ){

  #models need their row names.
#make case-based variant based on case-based rules.






  data_model <- data

  if ("DNA_age" %in% colnames(models)){
  models <- models[,1:(ncol(models)-2)]
  }

  if ("sigma" %in% colnames(models)){
    models <- models[,-which(colnames(models)=="sigma")]
  }

  #adapted to ignore models with no cases.
 if(!is.null(nrow(profiles))){
    if (nrow(models)>nrow(profiles)){
    keep_these <- rownames(profiles)
    models <- models[rownames(models) %in% keep_these,]
    }
 }

  data_for_sweep <- data_model



  if(method == "profiles"){

    #remove dependent variable
    data_for_sweep <- data_for_sweep[,-1, drop=F]
    #check for the DV!

    if (T %in% grepl("y mean",colnames(profiles))){
      profiles <- profiles[,-1, drop=F]

    }
        #remove model column if it exists
    if ("model" %in% colnames(data_for_sweep)){
  data_for_sweep <- data_for_sweep[,-ncol(data_for_sweep), drop=F]
    }

    if (T %in% grepl("model",colnames(profiles))){
      profiles <- profiles[,-ncol(profiles), drop=F]

    }


    #convert profile values into matrix


    #now compare it against the profiles by standardizing the difference in means in each variable
    x <- unlist( profiles)
    x <-  stringr::str_split_fixed(x, "/", 2)


    #now the first column is the means, put that into a matrix
    profile_matrix <- matrix(as.numeric(x[,1]), nrow=nrow(profiles), ncol=ncol(profiles))

    for(p in 1:nrow(profile_matrix)){





      eu_dist <-apply(X=data_for_sweep, MARGIN=1,FUN=function(x)  (sum(abs(x - ( profile_matrix[p,])))))

       data_model <- cbind(data_model, eu_dist)
      colnames(data_model)[colnames(data_model) == "eu_dist"] <- paste("Diff_Modell", p, sep = "")
    }



    data_model$model <- apply(data_model[,(ncol(data)+1):(ncol(data)+nrow(profiles)), drop=F], MARGIN = 1, FUN = which.min)
    data_model <- data_model[,-((ncol(data_model) -nrow(profiles)):(ncol(data_model) - 1))]
  } else if (method == "residuals"){


    dependent_variable <- data_for_sweep[,1]
    data_for_sweep[,1] <- 1  #make a constant for the intercept

    for(m in 1:nrow(models)){

      #ADD an intercept column!!
      this_model <- models[m,]

  if (is.null(this_model)){
    this_model <- t(data.frame(this_model))
  }
      #now run that model against every member of the input_data, collect residuals
      x <-sweep(data_for_sweep[,, drop=F], MARGIN=2,unlist(this_model), `*`) #input_data without dv column

      #now take row sums of x, and take their difference from happiness
      these_predictions <- rowSums(x)
     #this is made absolute. was raw before leading to the lowest values surfacing as the most important
       gap_from_Y <-  abs(these_predictions- dependent_variable)

      data_model <- cbind(data_model, gap_from_Y)
      colnames(data_model)[colnames(data_model) == "gap_from_Y"] <- paste("Residual_M", rownames(this_model), sep = "")
    }

    #remove the residual below
    if (nrow(models)>1){
model_indices <- apply(data_model[,(ncol(data)+1):(ncol(data)+nrow(models))], MARGIN = 1, FUN = which.min)

    data_model <- data.frame(data_model)
    data_model$model <- rownames(models)[unname(model_indices)]

    #move any columns holding residuals
    which_res <- grepl("Residual", colnames(data_model))
 # Temporarily turn off
    data_model <- data_model[,!which_res]

     }else{
      data_model$model <- 1
      data_model <- data_model[,-((ncol(data_model) -nrow(models)):(ncol(data_model) - 1))]

    }
  }else if (method == "case-based"){

#we take the universal model string from the population tree results
    # and apply that to the data that actually need to be assigned to models


    data_for_sweep$model <- unique(population_tree_results$data$model)



results <-     population_to_case_tree(population_data = data_for_sweep,
                                          cutlist = cutlist)
    data_model <- results$data


    #models need to be assigned here!

    }

  return(data_model)
}


#' Transform dummies
#'
#' Transforms dummysets back to standard variables
#' @export
transform_dummies <-  function(x, reference = NULL, group_terms = NULL, reference_cat = NULL){

   x <- as.data.frame(x)
  new_data <- x
  #group terms are the original column names applied to the column names of the categories

  for (g in 1:length(group_terms)){
    this_var <- as.data.frame(x[ , grepl( group_terms[g] , names(x))])
    if(ncol(this_var)>1){
      #Y is the new factor style variable representing the categorical values
      y <- colnames(this_var)[unlist(apply(this_var,1,which.max, simplify = F))]
       y[which(rowSums(this_var) == 0)] <- reference_cat[g]
       y <- gsub(pattern = paste0(group_terms[g],"_"),"", x= y)
    }else{

      y <- rep(colnames(x)[grepl( group_terms[g] , names(x))], nrow(this_var))
      y[which(this_var == 0)] <- reference_cat[g]
    }

    #maybe one cannot enter binary data as 'group terms'
    new_data <- new_data[,!grepl(group_terms[g],names(new_data))]
    new_data$y <- y
    colnames(new_data)[ncol(new_data)] <- group_terms[g]
  }

  return(new_data)
}


#' Used by case-based decision tree.
#'
#'Used by case-based decision tree.
#' @export
interpret_by_case <- function(data,
                              models,
                              profiles,
                              categorical_variables=NULL,
                              reference_categories=NULL,
                              min_bucket=0,
                              max_depth=NULL,
                              eliminate_multileaves=T,
                              assignments=NULL){
  ### with case-based DECISION TREE###

  #Assign cases to models
  if (is.null(assignments)){
  data_model <- model_closeness(data = data, profiles = profiles,
                                models = models, method = "residuals")
  }else{
    model_order <- assignments$model
    data$model <- model_order
    data_model <- data
  }

if (!is.null(categorical_variables)){
  data_cat <- transform_dummies(data_model, group_terms = categorical_variables, reference_cat = reference_categories)
}else{
  data_cat <- data_model
}




 if (!is.null(categorical_variables)){


   #! 'data_cat' has lost the original column names from the input of the function. How to shift them, make them match?####
  factorize_these <- c(categorical_variables, "model")

  for (f in 1:length(factorize_these)){

  # Factorize categorical variables
  data_cat[,grep(factorize_these[f], colnames(data_cat))] <- as.factor(  data_cat[,grep(factorize_these[f], colnames(data_cat))] )



  }
 }


data_cat$model <- as.factor(data_cat$model)

  # Calculate the tree
  #minsplit/minbucket/maxdepth in the control argument is a great way to limit the leaves.


#ignore_model_and_dv <- c(1, which(colnames(data_cat) == "model"))


  tree <- rpart(model ~ ., data = data_cat[,-1], method="class",
                control = c(minbucket=min_bucket, maxdepth=max_depth))

  #deeper tree for debugging only
  '
  tree <- rpart(model ~ ., data = data_cat[,-1], method="class",
                control = c(minbucket=1, maxdepth=10))
  '
  '
  par(mar = rep(1, 4))
  rpart.plot(tree, type = 1, extra = 1, cex = NULL, tweak = 1.0, clip.facs = T)
  '

  if (eliminate_multileaves==T){
  #when we do pruning, use tree$frame
  #check with end leaves are doubled, and limit model end leaves to 1.
  # list the parents
  #reorder ascending order of rownames
  new_tree <- tree$frame[order(as.integer(rownames(tree$frame)), decreasing = F),]
  full_tree <- new_tree

  #the parent referred to is the actual row index (not name) of the produced tree
  parent_logic <- c(NA, rep(1:nrow(new_tree), each=2,length.out=nrow(new_tree)-1))

  new_tree$parent <- parent_logic
  #parent_split <- c(NA,rep())

  leaves <- new_tree[ which(new_tree$var=="<leaf>"),]
  bare_tree <- new_tree[ -which(new_tree$var=="<leaf>"),]
  which_doubled_indices <- duplicated(leaves$yval) #leave indices truth vector. Only returns last one.
  which_doubled_values <- unique(leaves$yval[which_doubled_indices])


  new_trimmed_tree <- tree
  checked_these_doubles <- NULL

  while(length(which_doubled_values)>0){

    original_tree <- new_trimmed_tree
    #always just take the next doubled value on the list
    i <- 1

    #these are the doubled leaf indices
    double_indices <- unlist(lapply(which_doubled_values[i], FUN = function(x) which(leaves$yval==x)))
    #which of these has more final cases of the class?
    leaves <- as.data.frame(leaves)

    original_double_indices <- double_indices



    #compare how many cases in each of these competing leaves

    compare_cuts <- unlist(mapply(double_indices,which_doubled_values[i], FUN=function(x,y) leaves$yval2[x,1+y], SIMPLIFY = F))
    keep_this_index <-  which(compare_cuts==max(compare_cuts))
     compare_cuts[order(compare_cuts, decreasing = T)]
     #this is the leaf index that is kept
    keep_this_value <- unlist(double_indices)[unlist(keep_this_index)] #this is a leaf index
    keep_this_leaf <- rownames(leaves)[keep_this_value]

    for( l in 1:length(double_indices)){


    #these are leaf node names to remove
    remove_these_exact <- rownames(leaves)[ unlist(double_indices)[-unlist(keep_this_index)]]
    this_parents_children_to_remove <- lapply(remove_these_exact,FUN=function(x) leaves$parent[rownames(leaves) %in% x] )
    snip_below_this_node <- rownames(bare_tree)[ unlist(this_parents_children_to_remove)]


    #supply the rownames of the node to snip... (not the leaves)
    new_trimmed_tree <- snip.rpart(new_trimmed_tree, toss=as.integer(snip_below_this_node))

    #now check the new tree that it has all five leaves. If not, revert back to original tree and try other index

    new_tree_contains_all_models <- all(1:nrow(models) %in% unique(new_trimmed_tree$frame$yval))

    if (new_tree_contains_all_models==F){
      new_trimmed_tree <- original_tree
     next
    }else{
      #double_indices <- double_indices[0]
      break
    }

    } #end checking the double indices for this model value

    #plot for debug
    # par(mar = rep(1, 4))
    #rpart.plot(new_trimmed_tree, type = 1, extra = 1, cex = NULL, tweak = 1.5, clip.facs = T)

    #update leaves and bare_tree

    #We also need fresh parent lists!
    new_tree <- new_trimmed_tree$frame[order(as.integer(rownames(new_trimmed_tree$frame)), decreasing = F),]


    #the parent referred to is the actual row index (not name) of the produced tree
    parent_logic <- c(NA, rep(1:nrow(new_tree), each=2,length.out=nrow(new_tree)-1))

    new_tree$parent <- parent_logic

    leaves <- new_tree[ which(new_tree$var=="<leaf>"),]
    bare_tree <- new_tree[ -which(new_tree$var=="<leaf>"),]
    #lets update only the doubled values that have not just been covered

    #update doubled values

    checked_these_doubles <- c(which_doubled_values[i],checked_these_doubles )
    which_doubled_indices <- duplicated(leaves$yval) #leave indices truth vector. Only returns last one.

    which_doubled_values <- unique(leaves$yval[which_doubled_indices])[!unique(leaves$yval[which_doubled_indices]) %in% checked_these_doubles]


  }


  }else{
    new_trimmed_tree <- tree
  }

  return(new_trimmed_tree)

  }


#' Streamline
#'
#' Streamline forces underlying subpopulations to 'choose' one model.  It then plots a 'true regression' on those subpopulations. The purpose is to both increase prediction accuracy but more importantly to allow subpopulations to become distinct to aid in interpretation the updated model-subpopulation matchings are then tested to see how many cases are correctly 'solved' by their new sole models
#' @export
streamline <- function(models,
                       profiles,
                       data,
                       assign_models=T){  #this will *automatically* assign models through the algorithm. False means they are already assigned within the model data

if (!is.null(models)){
  if ("DNA_age" %in% colnames(models)){
    models <- models[,1:(ncol(models)-2)]
  }
}
  # we will sweep the condensed models against all cases and then assign the best model for each case:
if (assign_models==T){
data_with_models <- model_closeness(data=data,
                       profiles=profiles,
                       models=models,
                       method="residuals")
}else{
  data_with_models <- data
  data_with_models$model <- as.integer(data_with_models$model )
}
  #now calculate true regressions for each of these subpopulations!

  model_vector <-  unique( data_with_models$model)
  last_column <- ncol(data_with_models)
  dv <- dv
  ivs <- colnames(data_with_models)[!colnames(data_with_models)==dv]
  ivs <- ivs[!ivs=="model"]

  dna_pool <- data.frame(matrix(0, nrow=length(model_vector), ncol=dim(data_with_models)[2]-1))

  colnames(dna_pool) <- colnames(data_with_models)[-ncol(data_with_models)]


  if (!is.null(models)){
  #adapted to ignore models with no cases.
  if (nrow(models)>nrow(profiles)){
    keep_these <- rownames(profiles)
    models <- models[rownames(models) %in% keep_these,]
  }
  }

  for (i in 1:length(model_vector)){

    this_sample <- data_with_models[data_with_models$model==model_vector[i],-last_column]


    f <- as.formula(
      paste(dv,
            paste(ivs, collapse = " + "),
            sep = " ~ "))

    true_regression <- eval(bquote(   lm(.(f), data = this_sample)   ))

    dna_pool[i,] <- true_regression$coefficients


  }
#replace NA's with zero coefficients.####
  dna_pool[ is.na(dna_pool)] <- 0
  rownames(dna_pool) <-as.character( model_vector)

  #create profiles using function here!


 dna_profiles <- create_profile(data_after_dummies=df_full_dummies[training_indices,],
                                 model_assignment_data=data_with_models,
                                dummy_vars=original_variables[dummy_sets],
                                model_labels = NULL)




  return(list(dna_pool=dna_pool,
              model_assignment=data_with_models,
              dna_profiles=dna_profiles
              ))

}

#' Agglomerative decision tree
#'
#' A (novel?) agglomerative decision tree. It is built bottom-up, with one final leaf per model. Used for interpreting ICR results.
#' @param original_data Data before the dummyset splits
#' @export
agglom_tree <- function(data,
                        original_data=original_data[training_indices,],  #this one is without the dummy variables
                        re_use_variables=F){


 if (length( unique(data$model_assignment$model))<2){
   #if there is only one model
   populations <- data$model_assignment
   cut_list <- NA

 }else{

   dummy_df <- combinatorial_columns(original_data)




  models_by_case <- data$model_assignment$model
  dummy_df$model <- models_by_case
  populations <- dummy_df
  models <- unique(populations$model     )



  #one way anova
  all_ivs <- colnames(populations)[-c(1,ncol(populations))]
  dv <- "model"
  populations$model <- as.factor(populations$model)

  #principle of binary agglomeration... take the biggest gap (F-value from anova) of all variables,
  #then combine the two groups that are most similar by splitting between their mean.
  #then find widest gap in that pair= difference.
  #then split them from one another by difference and from rest by similiarity.

  #Bottom up population based tree####

  cut_list <- data.frame(matrix(nrow = 0, ncol=11))
  colnames(cut_list) <- c("cut", "leaf1", "leaf2", "leaf1_relation", "leaf2_relation", "parent_node", "leaf1_level", "leaf2_level", "parentnode_level", "node_location1",  "node_location2")


  i <- 0

  populations$model <- as.character(populations$model )
  original_models <- unique(populations$model )
  already_considered <- vector("logical", length(all_ivs))


  while(length( unique(populations$model))>1){

    i <- nrow(cut_list)+1
    models <- unique(populations$model)

    greatest_gap <- biggest_anova_gap(all_ivs[!already_considered], models, populations)

    if (is.na( greatest_gap$models_means)[1]){
      break
    }

    #we can only do the two-step agglomeration, split, if we have more than two models remaining.
    #otherwise, we simply split by greatest gap
    if (length( unique(populations$model))>2){

    #now we introduce the first split.. we take the greatest F vactor, and choose the narrowest neighbors to combine them to one group, then repeat.


    pairwise_comparison <- combn(greatest_gap$models_means$model, 2)
    abs_pairwise_gap <- apply(pairwise_comparison,2,FUN=function(x) abs(greatest_gap$models_means$x[greatest_gap$models_means$model== x[1]]-greatest_gap$models_means$x[ greatest_gap$models_means$model==x[2]]))
    smallest_pair_gap_biggest_F <-pairwise_comparison[, which(abs_pairwise_gap==min(abs_pairwise_gap))]
    #smallest pair_gap_biggest F... this correctly chooses the two most similar model in the greatest F factor overall.
    #these belong together in a node.

    #so e.g. Gender has biggest F value.
    # model 1 and 2 populations are the closest in gender
    # but models 1 and 2 differ the most in (Western-Northern Europe)

    if (re_use_variables==F){

      already_considered <- c(already_considered | (all_ivs %in% c(greatest_gap$variable)))

      #if this variable belongs to a dummyset, remove it!
      if (T %in% dummy_sets){
      if(grepl(paste0(dummy_set_ivs,  collapse = "|"), greatest_gap$variable)){
        #then remove all columns from the alreadyconsidered list   referring to that variable
        which_dummysets_used <- dummy_set_ivs[unlist(lapply(dummy_set_ivs, FUN=function(x) grepl(x,greatest_gap$variable )))]
        already_considered <- c(already_considered|(grepl(paste0(which_dummysets_used,  collapse = "|"), all_ivs)))


      }
      }


      greatest_gap_in_similar_models_otherwise <- biggest_anova_gap(all_ivs[!already_considered], #all ivs except the one already considered . we only split each iv once for parsimony
                                                                    models[models %in% smallest_pair_gap_biggest_F], #only the key pair
                                                                    populations)
    }else{

      #if this variable belongs to a dummyset, remove it!
      if (T %in% dummy_sets){
      if(grepl(paste0(dummy_set_ivs,  collapse = "|"), greatest_gap$variable)){
        #then remove all columns from the alreadyconsidered list   referring to that variable
        which_dummysets_used <- dummy_set_ivs[unlist(lapply(dummy_set_ivs, FUN=function(x) grepl(x,greatest_gap$variable )))]
        already_considered <- c(already_considered| (grepl(paste0(which_dummysets_used,  collapse = "|"), all_ivs)))


         }
      }
      greatest_gap_in_similar_models_otherwise <- biggest_anova_gap(all_ivs[!already_considered], #all ivs except the one already considered . we only split each iv once for parsimony
                                                                    models[models %in% smallest_pair_gap_biggest_F], #only the key pair
                                                                    populations)
    }

    #REPEAT?####
    #if this variable belongs to a dummyset, remove it!
    if (T %in% dummy_sets){
    if(grepl(paste0(dummy_set_ivs,  collapse = "|"), greatest_gap_in_similar_models_otherwise$variable)){
      #then add all columns from the alreadyconsidered list   referring to that variable

      which_dummysets_used <- dummy_set_ivs[unlist(lapply(dummy_set_ivs, FUN=function(x) grepl(x,greatest_gap_in_similar_models_otherwise$variable )))]
      already_considered <- c(already_considered| (grepl(paste0(which_dummysets_used,  collapse = "|"), all_ivs)))

    }
    }

   #now we remove from re-consideration the  greatest_gap_in_similar_models_otherwise variable
    if (re_use_variables==F){

      already_considered <- c(already_considered| (all_ivs %in% c(greatest_gap_in_similar_models_otherwise$variable)))
    }

    #if the greatest_gap_above is NA, then stop building the tree!
    if (is.na(greatest_gap_in_similar_models_otherwise$models_means)[1]){
      break
    }

    greatest_gap_in_similar_models_otherwise$models_means <-   greatest_gap_in_similar_models_otherwise$models_means[greatest_gap_in_similar_models_otherwise$models_means$model %in% smallest_pair_gap_biggest_F,]
    #initial split relative to the leafs
    initial_split <- combine_leafs(split_variable=greatest_gap_in_similar_models_otherwise)


    cut_list[i,] <- c(initial_split$cut1, initial_split$leaf1, initial_split$leaf2,initial_split$leaf1_relation, initial_split$leaf2_relation, initial_split$parent_node, initial_split$leaf1_level, initial_split$leaf2_level, initial_split$parent_node_level, 0, 0)

      #next relative to the leaf ends

    #combine the populations based on the above!
    populations$model[populations$model %in% c(initial_split$leaf1, initial_split$leaf2)] <- initial_split$parent_node
    } #end if for first step (if more than two models remain)

    #then find the greatest gap pair
    models <- unique(populations$model)

    greatest_gap <- biggest_anova_gap(greatest_gap$variable, models, populations)
    pairwise_comparison <- combn(greatest_gap$models_means$model, 2)
    abs_pairwise_gap <- apply(pairwise_comparison,2,FUN=function(x) abs(greatest_gap$models_means$x[greatest_gap$models_means$model== x[1]]-greatest_gap$models_means$x[ greatest_gap$models_means$model==x[2]]))
    biggest_pair_gap_biggest_F <-pairwise_comparison[, which(abs_pairwise_gap==max(abs_pairwise_gap))]

    greatest_gap$models_means <-   greatest_gap$models_means[greatest_gap$models_means$model %in% biggest_pair_gap_biggest_F,]

     next_split <- combine_leafs(split_variable=greatest_gap)
    i=nrow(cut_list)+1


    populations$model[populations$model %in% biggest_pair_gap_biggest_F] <-next_split$parent_node

    #merge the leaf populations together and repeat the above

     cut_list[i,] <- c(next_split$cut1, next_split$leaf1, next_split$leaf2,next_split$leaf1_relation, next_split$leaf2_relation, next_split$parent_node, next_split$leaf1_level, next_split$leaf2_level, next_split$parent_node_level, 0, 0)

  } #end while loop that keeps building the bottom-up tree


  #then iterate through and draw the cutlist


  #leaf locations in the tree_blueprint

  node_list <- as.vector(t(cut_list[,2:3]))

  #leaf list is only end leafs
  which_not_leafs <- lapply (gregexpr('&', node_list), FUN=function(x) (x>0)[1])

  leaf_list <- node_list[!unlist(which_not_leafs)]




  leaf_locations <- expand.grid(c(1,5), 0:ceiling((nrow(cut_list)*2+1)/4) )
  leaf_locations$location <- leaf_locations$Var2*6 + leaf_locations$Var1
  leaf_locations <- leaf_locations[1:length(original_models),]

  tree_blueprint <- vector("list",leaf_locations$location[nrow(leaf_locations)])  # a list of text lines that will later be parsed to build the tree
  #if a list item is Null, replace it with "\n"
  nullify <- unlist( lapply(tree_blueprint, is.null))
  #tree_blueprint[nullify] <- "\n"
  tree_blueprint[nullify] <- ""


  leaf_stem_length <- 5
  join_stem_length <- 20
  node_label_size <- 8
  cut_descriptor_label_size <- 18
  relation_label_length <- 11 #length of MORE or LESS label
  #"    )" is the standard for ending cat parentheses



#this this to prevent cuts to the cutlist
cut_list$original_cut <-  cut_list$cut
  for (this_cut in 1:length(cut_list$cut)){
    if (nchar(cut_list$cut[this_cut])==cut_descriptor_label_size){ #same size  as max label length
      #skip
    }else if (nchar(cut_list$cut[this_cut])<cut_descriptor_label_size){
      #add spaces to make max length
      extended <- str_pad(cut_list$cut[this_cut], width = cut_descriptor_label_size, 'both', ' ')
      cut_list$cut[this_cut] <- extended
    }else if (nchar(cut_list$cut[this_cut])>cut_descriptor_label_size){
      #cut off excess
      trimmed <-     substr(cut_list$cut[this_cut],1,cut_descriptor_label_size)
      cut_list$cut[this_cut] <- trimmed
    }

  }

  #leaf1
  for (this_node in 1:length(cut_list$leaf1)){
    node_name <- cut_list$leaf1[this_node]

    if (nchar(node_name)==node_label_size){ #same size  as max label length
      #skip
    }else if (nchar(node_name)<node_label_size){
      #add spaces to make max length
      extended <- str_pad(node_name, width = node_label_size, 'both', ' ')
      cut_list$leaf1[this_node] <- extended
    }else if (node_name>node_label_size){
      #cut off excess? No... not for leafs. will cause errors if trimmed
    }

  }

  #leaf2
  for (this_node in 1:length(cut_list$leaf2)){
    node_name <- cut_list$leaf2[this_node]

    if (nchar(node_name)==node_label_size){ #same size  as max label length
      #skip
    }else if (nchar(node_name)<node_label_size){
      #add spaces to make max length
      extended <- str_pad(node_name, width = node_label_size, 'both', ' ')
      cut_list$leaf2[this_node] <- extended
    }else if (node_name>node_label_size){
      #cut off excess? No... not for leafs. will cause errors if trimmed
    }

  }


  #parentnode
  for (this_node in 1:length(cut_list$parent_node)){
    node_name <- cut_list$parent_node[this_node]

    if (nchar(node_name)==node_label_size){ #same size  as max label length
      #skip
    }else if (nchar(node_name)<node_label_size){
      #add spaces to make max length
      extended <- str_pad(node_name, width = node_label_size, 'both', ' ')
      cut_list$parent_node[this_node] <- extended
    }else if (node_name>node_label_size){
      #cut off excess? No... not for leafs. will cause errors if trimmed
    }

  }


  for (r in 1:nrow(cut_list)){

    this_aggregation <- cut_list[r,]


    if (r==1){
      #the first two leaves are built identically in any tree
      # keep each list item in the blueprint independently printable


      #make join node the proper size
      join_node <- cut_list$parent_node[r]



      tree_blueprint[1] <- paste0('cat( cut_list[1,]$leaf1,paste0(rep("-",',leaf_stem_length,'), collapse = ""),cut_list[1,]$leaf1_relation, cut_list[1,]$cut, sep="")')
      tree_blueprint[2] <- paste0('cat("\n", paste0(rep(" ",node_label_size+leaf_stem_length), collapse=""),"|", sep="")')
      tree_blueprint[3] <-paste0('cat("\n",paste0(rep(" ",node_label_size+leaf_stem_length), collapse=""),"|" ,paste0(rep("-",',join_stem_length,'), collapse="") ,cut_list$parent_node[',r,'], sep="") ')
      tree_blueprint[4] <- paste0('cat("\n",paste0(rep(" ",node_label_size+leaf_stem_length), collapse=""),"|", sep="")')
      tree_blueprint[5] <- paste0('cat("\n",cut_list[1,]$leaf2,paste0(rep("-",',leaf_stem_length,'),collapse=""), cut_list[1,]$leaf2_relation, cut_list[1,]$cut, sep="")')


      cut_list$node_location1[r] <- 1
      cut_list$node_location2[r] <- 5
    }else{

      if (this_aggregation$leaf1_level==0){
        leaf1_location <- leaf_locations$location[ which(leaf_list==gsub(" " , "", cut_list[r,]$leaf1 )) ] #calling the leaf's location in the blueprint

      }else{

        this_rows_parent_node_index <-  which(gsub(" " , "",cut_list$parent_node)==  gsub(" " , "", cut_list[r,]$leaf1 ))

        if (as.integer(cut_list$node_location1[this_rows_parent_node_index])>0 & as.integer( cut_list$node_location2[this_rows_parent_node_index])>0){
          #find its location in the parent nodes... and then add two levels to the first leaf in that parent node
          this_node_new_location <- round(mean(c(as.integer(cut_list$node_location1[this_rows_parent_node_index]),as.integer( cut_list$node_location2[this_rows_parent_node_index]))))
          leaf1_location <- this_node_new_location
        }else{
          first_leaf_this_row <- gsub(" " , "",cut_list[this_rows_parent_node_index,'leaf1'])
          that_leaf_location <-  which(node_list==first_leaf_this_row)
          this_node_new_location <- leaf_locations$location[that_leaf_location]+2
          leaf1_location <- this_node_new_location
        }


      }

      if (this_aggregation$leaf2_level==0){
        leaf2_location <-  leaf_locations$location[ which(leaf_list== gsub(" " , "",cut_list[r,]$leaf2 )) ]

      }else{
        this_rows_parent_node_index <-  which(gsub(" " , "",cut_list$parent_node)== gsub(" " , "", cut_list[r,]$leaf2 ))

        if (as.integer(cut_list$node_location1[this_rows_parent_node_index])>0 & as.integer( cut_list$node_location2[this_rows_parent_node_index])>0){
          #find its location in the parent nodes... and then add two levels to the first leaf in that parent node
          this_node_new_location <- round(mean(c(as.integer(cut_list$node_location1[this_rows_parent_node_index]),as.integer( cut_list$node_location2[this_rows_parent_node_index]))))
          leaf2_location <- this_node_new_location

        }else{
          first_leaf_this_row <- gsub(" " , "",cut_list[this_rows_parent_node_index,'leaf2'])
          that_leaf_location <-  which(node_list==first_leaf_this_row)
          this_node_new_location <- leaf_locations$location[that_leaf_location]+2
          leaf2_location <- this_node_new_location

        }

      }

      cut_list$node_location1[r] <- leaf1_location
      cut_list$node_location2[r] <- leaf2_location

      #now appropriately add both leaves
      #add leaf 1
      if ( tree_blueprint[[ leaf1_location]] ==""){
        if (this_aggregation$leaf1_level==0){
          # add a brand new leaf
          tree_blueprint[leaf1_location] <-paste0(  'cat("\n", cut_list[',r, ',]$leaf1,paste0(rep("-",',leaf_stem_length,'), collapse = ""),cut_list[',r,',]$leaf1_relation, cut_list[',r,',]$cut, sep="")')
        }else{
          stem_length <- leaf_stem_length+1+join_stem_length+node_label_size+leaf_stem_length
          tree_blueprint[leaf1_location] <-paste0(  'cat( "\n",cut_list[',r, ',]$leaf1,paste0(rep("-",',stem_length,'), collapse = ""),cut_list[',r,',]$leaf1_relation, cut_list[',r,',]$cut, sep="")')

        }

      }else{
        #extend the other leaf outward

        tree_blueprint[leaf1_location] <- paste0( gsub(', sep="")' , "", tree_blueprint[[leaf1_location]]) ,
                                                  ',paste0(rep("-",',leaf_stem_length,'),collapse=""),
                                              cut_list[',r,',]$leaf1_relation,
                                              cut_list[',r,',]$cut, sep="")'  )
      }

      #add leaf 2
      if ( tree_blueprint[[ leaf2_location]] ==""){

        if (this_aggregation$leaf2_level==0){
          # add a brand new leaf
          tree_blueprint[leaf2_location] <-paste0(  'cat( "\n",cut_list[',r, ',]$leaf2,paste0(rep("-",',leaf_stem_length,'), collapse = ""),cut_list[',r,',]$leaf2_relation, cut_list[',r,',]$cut, sep="")')

        }else{
          stem_length <- leaf_stem_length+1+join_stem_length+node_label_size+leaf_stem_length
          tree_blueprint[leaf2_location] <-paste0(  'cat( "\n",cut_list[',r, ',]$leaf2,paste0(rep("-",',stem_length,'), collapse = ""),cut_list[',r,',]$leaf2_relation, cut_list[',r,',]$cut, sep="")')


        }

      }else{
        #extend the other leaf outward

        tree_blueprint[leaf2_location] <-paste0( gsub(', sep="")' , "", tree_blueprint[[leaf2_location]]) ,
                                                 ',paste0(rep("-",',leaf_stem_length,'), collapse=""),
                                              cut_list[',r,',]$leaf2_relation,
                                              cut_list[',r,',]$cut, sep="")'  )

      }


      #now link the two nodes together.
      join_location <- min(leaf1_location, leaf2_location)+ round(abs( leaf2_location - leaf1_location)/2)
      join_level <- as.integer( max(this_aggregation$leaf1_level, this_aggregation$leaf2_level) )# max of the two node levels


      #Do actual join location first... then get the length of the base...and calculate the others based on that?
      #     length of this base in terms of position
      invisible_base_stem_length <-  node_label_size  +leaf_stem_length  +(join_stem_length+ node_label_size+leaf_stem_length)*join_level

      this_join_base <- paste0('cat(paste0(rep(" ",', invisible_base_stem_length
                               ,'), collapse="")',
                               ',"|")')

      base_length <-      max(nchar(capture.output(eval(parse(text=gsub("\n", "", this_join_base))))))

      if (tree_blueprint[[ join_location]] ==""){ #add new
        #invisible_base_stem_length <-  node_label_size  +leaf_stem_length  +(relation_label_length  +cut_descriptor_label_size +join_stem_length-(cut_descriptor_label_size+relation_label_length)+node_label_size+leaf_stem_length )*join_level

        tree_blueprint[join_location] <- paste0('cat("\n",paste0(rep(" ",', base_length
                                                ,'), collapse="")',
                                                ',"|",paste0(rep("-",',join_stem_length,'), collapse="") ,cut_list$parent_node[',r,'], sep="")')


      }else{

        current_length <-  max(nchar(capture.output(eval(parse(text=gsub("\n", "", tree_blueprint[[join_location]]))))))
        extension_length <- max(0, base_length-current_length)

        tree_blueprint[join_location] <- paste0( gsub(', sep="")' , "",tree_blueprint[[join_location]] ),
                                                 ',paste0(rep(" ",',extension_length,
                                                 '), collapse=""),"|",paste0(rep("-",',join_stem_length,'), collapse="") ,cut_list$parent_node[',r,'], sep="")'  )



      }

      #Now loop through the other base locations
      other_base_locations <-    leaf1_location:leaf2_location
      other_base_locations <- other_base_locations[which(!other_base_locations %in% c(join_location, leaf1_location, leaf2_location))]

      for (o in other_base_locations){
        if ( tree_blueprint[[o]] ==""){ #add new
          current_length <-  0
          extension_length <- max(0,base_length-current_length)

          tree_blueprint[o] <-   paste0( 'cat("\n",paste0(c(rep(" ",',extension_length,')',',"|"), collapse=""), sep="")'     )
        }else { #extend

          current_length <-  max(nchar(capture.output(eval(parse(text=gsub("\n", "", tree_blueprint[[o]]))))))
          extension_length <- max(0,base_length-current_length)

          #the positioning will be equal to the null case invisible base stem length minus the existing parts in this stem
          #added_invisible_base_stem_length <-  node_label_size  +leaf_stem_length  +(relation_label_length   +cut_descriptor_label_size +join_stem_length-(cut_descriptor_label_size+relation_label_length)+node_label_size+leaf_stem_length)*join_level - node_label_size- leaf_stem_length-relation_label_length-cut_descriptor_label_size


          tree_blueprint[o] <- paste0( gsub(', sep="")' , "", tree_blueprint[[o]]),
                                       ',paste0(rep(" ",',extension_length,'), collapse="")',',"|", sep="")')

        }

      }

      #now if the two leaf levels are different, the smaller one needs to be extended out to the base
      if ( this_aggregation$leaf1_level> this_aggregation$leaf2_level){


        current_length <-  max(nchar(capture.output(eval(parse(text=gsub("\n", "", tree_blueprint[[leaf2_location]]))))))
        extension_length <- max(0,base_length-current_length)
        tree_blueprint[leaf2_location] <- paste0( gsub(', sep="")' , "",  tree_blueprint[[leaf2_location]]),
                                                  ',paste0(rep("-",',extension_length,'), collapse=""), sep="")')


      }else if (this_aggregation$leaf2_level>this_aggregation$leaf1_level){

        current_length <-  max(nchar(capture.output(eval(parse(text=gsub("\n", "", tree_blueprint[[leaf1_location]]))))))
        extension_length <- max(0,base_length-current_length)

        tree_blueprint[leaf1_location] <- paste0( gsub(', sep="")' , "",  tree_blueprint[[leaf1_location]]),
                                                  ',paste0(rep("-",',extension_length,'), collapse=""), sep="")')

      }







    }




  }

  #remove all nulls from the blueprint








  #this prints out the tree in text
  sink("agglomeration_tree.txt")
  eval(parse(text=tree_blueprint))
  sink()



 }

  return(list(cutlist=cut_list, data=populations, tree_blueprint=parse(text=tree_blueprint)))

}


#' Create combinatorial columns
#'
#' Creates combinatorial columns for use in the agglomerative decision tree
#' @param input_data Dataframe to be converted into combinatorial columns.
#' @export
combinatorial_columns <- function(input_data){
  #for dummy variable sets, will build their combinations, so that splits can be made properly between multiple categories (e.g. Eastern and Southern Europe vs. the rest)


  dummy_set_ivs <<- original_variables[dummy_sets]

  #now take all possible combinations of that dummy variable and keep all but the full combinations (So that there can be a comparison)

  initial_levels_list <- vector("list", length(dummy_set_ivs))
  combos_list <- vector("list", length(dummy_set_ivs))

   for (i in 1:length(dummy_set_ivs)){

    if (length(dummy_set_ivs)==0){
      break
    }

    initial_levels <-   unique(original_data[,dummy_set_ivs[i]])


    combos <-  as.matrix(gtools::combinations(length(initial_levels), length(initial_levels)-1, v=initial_levels, repeats.allowed = T))

    this_is_duplicated <-   t(unlist(apply(combos, 1, stri_duplicated, simplify = T)))

    combos[this_is_duplicated ] <- NA

    initial_levels_list[[i]] <- initial_levels
    combos_list[[i]] <- combos

  }

  #now create dummy variables for each of these below combos and remove the original dummy variables

  if (length(dummy_set_ivs)>0){
    save_rownames <- rownames(input_data)

   #dummy columns should be re-done, because they ones in the data have reference categories removed!
#so remove dummy columns and replace with the original column

    col_to_remove <- colnames(input_data)[grepl(paste0(dummy_set_ivs, collapse = "|"), colnames(input_data))]
    col_to_remove <- col_to_remove[!col_to_remove %in% dummy_set_ivs]

    #all categorical variables are gone now.
    data_with_no_cat_variables <- input_data[,!colnames(input_data) %in% col_to_remove]

        if(!identical(dummy_set_ivs, intersect(dummy_set_ivs, colnames(input_data)))){

          #to add the original data back, we first figure out if these data are from the training or the test set.

          if (nrow(input_data)==nrow(training_set)){
            this_set <-  original_data[training_indices,dummy_set_ivs]
          }else if (nrow(input_data)==nrow(testing_set)){
            this_set <-  original_data[-training_indices,dummy_set_ivs]
          }

          data_with_no_cat_variables[,dummy_set_ivs] <- this_set
        }
        working_df <- data_with_no_cat_variables

        for (d in 1:length(dummy_set_ivs) ){
        this_dummy <- dummy_set_ivs[d]
        working_df <- dummy_cols(working_df, select_columns = this_dummy, remove_first_dummy = F, ignore_na = T, remove_selected_columns = T)



        }

        dummy_df <- working_df
        rownames(dummy_df) <- save_rownames

    }else{
      dummy_df <- input_data
    }

  #match categories to column numbers

  if (length(dummy_set_ivs)>0){



    for (d2 in 1:length(dummy_set_ivs)){
    these_columns <-  paste0( dummy_set_ivs[d2], "_", initial_levels_list[[d2]])
    combos <- combos_list[[d2]]

    for (combo in 1:nrow(combos)){
      this_combo <-  combos[combo,][!is.na(combos[combo,])]
      combo_subset <- dummy_df[,these_columns[ grepl(paste0(this_combo, collapse="|"), these_columns)], drop=F]
      this_new_colname <- paste0(c( dummy_set_ivs[d2], "_", paste0(this_combo, collapse="_OR_")), collapse="")
      new_values <- apply(combo_subset,1, FUN = function(x) as.integer(1 %in% x))

      if (this_new_colname %in% colnames(dummy_df)){

      }else{

        dummy_df <- eval(parse(text = paste0("cbind(dummy_df,", this_new_colname,"= new_values)" )))
      }
    }
    }

  }

  return(dummy_df)
}

#' Population tree to case tree
#'
#' Applies population differences (from the agglomerative tree) as case-based differences to data set
#' @param population_data Population-based decision tree data from agglom_tree function.
#' @param cutlist This is the blueprint for building the agglomeration tree, exported from the agglom_tree() function.
#' @export
population_to_case_tree <- function(cutlist,
                                    population_data){

  #determine if combinatorial columns are present!
  #If even one variable is converted, then they all are...
  #identify number of columns associated with the categorical variable
  dummy_set_ivs_columns <- sum( grepl(dummy_set_ivs[1], colnames(population_data)))

  if(!is.na(dummy_set_ivs_columns)){

 original_length <- length(unique( original_data[,dummy_set_ivs[1]]))
  if (dummy_set_ivs_columns== original_length-1){
    #if comb columns not present not use combinatorial_columns here... otherwise the final model assignments are not accurate.

    processed_data <- combinatorial_columns(population_data)
    processed_data$model <- population_data$model
  }else{
    #if combinatorial columns already present, leave columns alone
    processed_data <- population_data


  }


  }else{ #if no combinatorial columns present
    processed_data <- population_data
  }





  #cycle through the cutlist in an inverse order, converting population comparisons to case-based rules
  #adding new model assignment to the end.
  #testing the predictive power of this (is it at least much better than standard regression?)
 #sometimes the processed data ends with no model column!



  #start with models from agglom tree, as the model structure is already established there.


  if (is.na(cutlist)[1]){
    #if cutlist is NA, there is only one model!

    processed_data$model <- 1
  }else{





        if (!"model" %in% colnames(processed_data)){
      #which is the master parent node containing all models?
  main_parent_node_row <-  which(  nchar(gsub(" ", "", cutlist$parent_node))==max( nchar(gsub(" ", "", cutlist$parent_node))))

      processed_data$model <- cutlist$parent_node[ main_parent_node_row]
        }else{
      #pass population_data models to these processed data, which include now all dummy categories

          processed_data$model <- population_data$model
    }

  for (i in nrow(cutlist):1){
    this_cut <- cutlist[i,]
    this_variable <- gsub(" ", "", this_cut$original_cut)
    leaf1_part <- gsub(") MORE", "",this_cut$leaf1_relation)
    leaf1_part <- gsub(") LESS", "",leaf1_part)
    leaf1_part <- gsub("\\(", "",leaf1_part)
    leaf2_part <- gsub(") MORE", "",this_cut$leaf2_relation)
    leaf2_part <- gsub(") LESS", "",leaf2_part)
    leaf2_part <- gsub("\\(", "",leaf2_part)
    #this is a partition between the two subgroup means of leaf1 and leaf2
    threshhold <- mean(c(as.numeric( leaf1_part), as.numeric( leaf2_part)))
    leaf_1_action <-str_sub(this_cut$leaf1_relation, -4, -1)
    which_var_column <- which(grepl(this_variable, colnames(processed_data)))
    if (length(which_var_column)>1){ #if more than one variable matches, then look for an exact match or a minimal match
      x <- which(colnames(processed_data)==this_variable)
      if (length(x)==1){
        which_var_column <- x
      }else{
        x <- which(nchar(colnames(processed_data)[which_var_column])==min(nchar(colnames(processed_data)[which_var_column])))
        which_var_column <- x
        }
    }


    parent_mask <- (grepl(gsub(" ", "", this_cut$parent_node), processed_data$model))
    if (leaf_1_action=="MORE"){
      threshhold_mask <-  processed_data[,which_var_column] >threshhold
      processed_data$model[parent_mask & threshhold_mask ] <- this_cut$leaf1
      processed_data$model[parent_mask & !threshhold_mask] <- this_cut$leaf2
      #but now the cuts do not apply to all data but just some
    }else{
      threshhold_mask <-  processed_data[,which_var_column] <= threshhold
      processed_data$model[parent_mask & threshhold_mask ] <- this_cut$leaf1
      processed_data$model[parent_mask & !threshhold_mask] <- this_cut$leaf2
    }
  }

  }

  #use create_profile function!






  dna_profiles <- create_profile(data_after_dummies=df_full_dummies[training_indices,],
                                 model_assignment_data=processed_data,
                                 dummy_vars=original_variables[dummy_sets],
                                 model_labels = NULL)




  return(list(data= processed_data, profiles=dna_profiles))

}



#' Combine leaves
#'
#' This function combines leafs based on a particular variable split. In an agglomerative manner this builds a population-based tree from the bottom up. Part of the agglomerative tree building proces
#' @export
combine_leafs <- function(split_variable){  #variable plus the model means for the variable

  #STart by splitting the two models that are most similar in the overall greatest gap by their pairwise greatest difference
  if(split_variable$models_means$x[1]>split_variable$models_means$x[2]){

    leaf1_relation <- paste0("(",round(split_variable$models_means$x[1], digits = 2), ") ","MORE")
    leaf2_relation <- paste0("(",round(split_variable$models_means$x[2], digits = 2), ") ","LESS")
  }else{

    leaf1_relation <- paste0("(",round(split_variable$models_means$x[1], digits = 2), ") ","LESS")
    leaf2_relation <- paste0("(",round(split_variable$models_means$x[2], digits = 2), ") ","MORE")

  }

  #now split them by this difference, and then split them fromthe rest of the models one node above according to their similarity.
  #somehow save these leaves and cuts so the whole thing can be built later
  leaf1 <-  split_variable$models_means$model[1]
  leaf2 <-  split_variable$models_means$model[2]
  cut1 <- split_variable$variable
  parent_node <- paste(  as.character( c(leaf1, leaf2)), collapse = "&")

  #calculate the line spacing for each leaf and put it in table.
  #count number of "&" signs in the leafs. 0= base level, 1 = 1 level higher, etc.
  # then identify highest leaf (number of "&" signs.. that determines the line spacing for the parent node)
  leaf1_level <- lengths(regmatches(leaf1, gregexpr("&", leaf1)))
  leaf2_level <- lengths(regmatches(leaf2, gregexpr("&", leaf2)))
  parent_node_level <- ifelse(leaf1_level>=leaf2_level, leaf1_level, leaf2_level)+1
  return(list(cut1=cut1, leaf1=leaf1, leaf2=leaf2, leaf1_relation=leaf1_relation, leaf2_relation=leaf2_relation, parent_node=parent_node, leaf1_level=leaf1_level, leaf2_level=leaf2_level, parent_node_level=parent_node_level) )
}

#' Biggest anova gap
#'
#' Identifies largest ANOVA gap between subgroups. Part of agglomerative tree building process.
#' @export
biggest_anova_gap <- function(ivs, #the ivs that will be compared
                              these_models, #the model populations across which they will be compared
                              populations, #the population dataframe, that has models as a column
                              F_threshold = 1.25
                              ){


  all_ivs <- ivs
#dv is globally defined
    models <- these_models

  populations <- populations[populations$model %in% models,]

  anova_results <- vector("list", length(all_ivs))
  best_F <- 0
  choose_me <- NULL #index if i
  for (i in 1:length(all_ivs)){
    this_iv <- all_ivs[i]
    dv=dv
    f <- as.formula(
      paste(this_iv,
            paste(dv, collapse = " + "),
            sep = " ~ "))


    one.way <- aov(f, data = populations)
    #take variable with highest or lowest F value
    #summary(one.way)

    model_comparison <- aggregate(populations[,this_iv], by=list(model= populations$model[populations$model %in% models]), mean)
    F_value <- summary(one.way)[[1]]$`F value`[1]
    if (is.nan( F_value)){
      F_value <- 0
    }

    if (F_value>best_F){
      best_F <- F_value
      choose_me <- i
    }

    anova_results[[i]] <- list(this_var=this_iv, F_value=F_value, mean_comparison=model_comparison)

  }

if (anova_results[[choose_me]]$F_value>=F_threshold){

  biggest_F <- anova_results[[choose_me]]$mean_comparison
}else{ #lower than threshold
  biggest_F <- NA

}

  return(list(variable=  anova_results[[choose_me]]$this_var, models_means=biggest_F))
}

#' Model plot
#'
#' Plot group-based bivariate regression lines
#' @export
model_plot <- function(data_with_models,
                       label,
                       base_line_plot=F, #this plots the base models themselves in addition to the assigned point lines
                       models, #model data if baseline plot is T, this is required
                       turned_on=F,
                       X_variable=1, #which X variable will be used, by index
                       point_labels=T,
                       local_regression_line=F  #prints the local regression for these two variables
                       ){

  #if it is activated, it prints, if not, the command is ignored

  if (turned_on==T){
    ggplot_order <- vector("list",4) ## a list of text lines that will later be parsed to build plot
    ggplot_order[1] <- paste0('ggplot(data_with_models, aes(x=',colnames(data_with_models)[X_variable+1],', y=y, color=as.factor(model))) +')
    ggplot_order[2] <- paste0('geom_point(size=1.5) +')
    if (local_regression_line==T){
    ggplot_order[3] <- paste0('geom_smooth(method = "lm", fill = NA, linewidth=1.5)+')
    }
    ggplot_order[4] <- paste0(' ggtitle(label)')


  #my_colors <- c("#1170AA", "#55AD89", "#EF6F6A", "yellow", "black")


 if ( base_line_plot==T){
   n_models <- nrow(models)
existing_ggplot_order <- length(ggplot_order)

#first append a + to the previous line
ggplot_order[existing_ggplot_order] <- paste0(ggplot_order[existing_ggplot_order], "+")

for (i in 1:n_models){

  ggplot_order[existing_ggplot_order+i] <- paste0(' geom_textabline(slope=models[',i,',',X_variable+1,'], intercept= models[',i,',1] , colour="red", linewidth=.5, linetype="longdash", label=',rownames(models)[i],') +')
}

#and remove the plus from the last line

ggplot_order[existing_ggplot_order+n_models] <- gsub("\\+", "", ggplot_order[existing_ggplot_order+n_models])

 }

 ggplot_order<-ggplot_order[!sapply(ggplot_order,is.null)]
 print(eval(parse(text=ggplot_order)))

  }
}





#' Compare data
#'
#' Compare original sample characteristics and original model characteristics with outcome####
#' @export
compare_data <- function(true_sample=sample.data$profile,
                         found_sample=true_regression_case_based$dna_profiles,
                         true_case_assignments=sample.data$data$model,
                         found_case_assignments = true_regression_case_based$model_assignment,
                         dna_pool=sample.data$models,
                         found_models=true_regression_case_based$dna_pool){


  # Case assignments.####
  #just compare the true with found and give a percentage accuracy.

  #sort those found by row number so that we can easily match the sample data
  found_case_rows <- rownames(found_case_assignments)
  found_case_models <- found_case_assignments$model

  #now try to match the names of the found, vs. not found models

  true_case_assignments <- true_case_assignments[match(found_case_rows, rownames(true_case_assignments)),]

  case_assign_comparison <- data.frame(matrix(nrow=0, ncol=3))
  colnames(case_assign_comparison) <- c("original model", "new model", "percentage")

  for (i in unique(true_case_assignments$model)){

    these_cases <-     which(true_case_assignments$model==i)

    models_found_cases <- found_case_models[these_cases]

    models_to_cases <- data.frame(table(models_found_cases))

    models_to_cases <- cbind(i, models_to_cases)

    case_assign_comparison <- rbind(case_assign_comparison, models_to_cases)


  }

  case_assign_comparison$models_found_cases <- gsub(" ", "",case_assign_comparison$models_found_cases)
  case_assign_comparison$models_found_cases <- as.character( case_assign_comparison$models_found_cases)
 # colnames(case_assign_comparison)[1] <- "true models"

  #now I can  score the case_assign_comparison
  #a perfect assignment would mean.
  #1. all original cases are grouped together in the same new model.
  #2. the number of new models is the same
  #so 1. can be a percentage in each original model that share a new dominant model
  #and 2. can be abs( new-models minus old)/old models.. and that penaltyshould range between 0 and 1 (max)
  case_scores <- NULL
  model_sim_scores <- NULL
  for (i in 1:length(unique(case_assign_comparison$i))){
    this_group <- case_assign_comparison[case_assign_comparison$i==i,]

    #case assignment score
    group_score <-   max(this_group$Freq)/sum(this_group$Freq)
    case_scores <- append(case_scores, group_score)

    #model similarity score
    #model comparison####

    #now we can count the error between the original models and the assigned models, using the case assignment % as a weight.

    these_models <- this_group[,-3]
    colnames(these_models) <- c("original", "found")

    original_model <- dna_pool[  these_models[1,1],]

    new_models <- found_models[match(these_models$found,rownames(found_models)  ),]

    if (is.null(dim(new_models))){
    new_models <- t(data.frame(new_models))
    }

        colnames(new_models) <- colnames(original_model)

    #now sweep for the differences
    x <-sweep(new_models, MARGIN=2,unlist( original_model), `-`) #input_data without dv column
    x <- x^2
    x <- rowSums(x)
    #now multiply the error by its weight based on case assignments

    group_case_total <- sum(this_group$Freq)
    case_percentages <- this_group$Freq/group_case_total

    x <- x*case_percentages
    #now sum up the error, and we have the model_based error
    model_similarity_error <- sum(x)
    model_sim_scores <- append(model_sim_scores, model_similarity_error)

  }
  case_assign_score <- mean(case_scores)
  model_similarity_error <-  mean(model_sim_scores)

  if (is.na( model_similarity_error)){
    #
  }


  #and now apply the n_model penalty, limited it to a 0 to .9 range
  x <- min(abs(  length(unique(case_assign_comparison$models_found_cases)) - length(unique(case_assign_comparison$i)))/ length(unique(case_assign_comparison$i)), .9)
  case_assign_score <- case_assign_score*(1-x)

  #profile comparison perhaps not needed because we already look at exact case matches...










  return(list(case_assign_score=case_assign_score, case_assign_comparison_data=case_assign_comparison, model_similarity_error=model_similarity_error))

}



#' Deprecated?
#'
#' Track how well the profiles are being correctly found or not
#' @export
check_profile_performance <- function(original_profiles, comparison_profiles){
  origin <- sample.data$profiles
  found <- true_regression_case_based$dna_profiles[,-1]

  x1 <- unlist( origin)
  x1 <-  stringr::str_split_fixed(x1, "/", 2)


  #now the first column is the means, put that into a matrix
  origin_matrix <- matrix(as.numeric(x1[,1]), nrow=nrow(origin), ncol=ncol(origin))

  x2 <- unlist( found)
  x2 <-  stringr::str_split_fixed(x2, "/", 2)

  found_matrix <- matrix(as.numeric(x2[,1]), nrow=nrow(found), ncol=ncol(found))

  found_matrix
  origin_matrix
}


