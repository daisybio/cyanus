evap <- function(function_call) {
  call <- substitute(function_call)
  reactiveVals$call_list <<- append(reactiveVals$call_list, call) # in the shiny app, we don't have to make this global call, I think
  return_value <- eval(call)
  return(return_value)
}

setup_logfile <- function(){
  ##### function calls for setup that always have to be written / executed
  reactiveVals$call_list[[1]] <- quote("# Log file for the CYANUS Shiny App. The following code is compatible with R 4.2.1, renv 0.16, and UNIX systems (Linux, MacOS)")
  reactiveVals$call_list[[2]] <- substitute(system('git clone https://github.com/biomedbigdata/cyanus_functions.git'))
  reactiveVals$call_list[[3]] <- substitute(system('mv cyanus_functions/* ./'))
  reactiveVals$call_list[[4]] <- substitute(system('rm -r cyanus_functions'))
  reactiveVals$call_list[[5]] <- '# Choose option 1: Restore the project from the lockfile. [...] Do you want to proceed? [Y/n]: Y'
  reactiveVals$call_list[[6]] <- substitute(renv::init())
  reactiveVals$call_list[[7]] <- quote(reactiveVals <- list())
  reactiveVals$call_list[[8]] <- substitute(lapply(list.files('functions', full.names = T), source))
  reactiveVals$call_list[[9]] <- quote(reactiveVals$data <- reactiveVals$data <- list(upload = list(fcs=NULL, panel=NULL, md=NULL), example = list(fcs=NULL, panel=NULL, md=NULL) ))
}  