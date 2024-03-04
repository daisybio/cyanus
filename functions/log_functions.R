evap <- function(function_call, params = NULL) {
  call <- do.call("substitute", list(function_call, params))
  reactiveVals$call_list <<- append(reactiveVals$call_list, call)
  return_value <- eval(call, envir = parent.frame())
  return(return_value)
}

setup_logfile <- function(){
  ##### function calls for setup that always have to be written / executed
  reactiveVals$call_list[[1]] <- quote("# Log file for the CYANUS Shiny App. The following code is compatible with R 4.2.1, renv 1.0.3, and UNIX systems (Linux, MacOS)")
  reactiveVals$call_list[[2]] <- substitute(system('git clone https://github.com/biomedbigdata/cyanus_functions.git'))
  reactiveVals$call_list[[3]] <- substitute(system('mv cyanus_functions/* ./'))
  reactiveVals$call_list[[4]] <- substitute(system('rm -r cyanus_functions'))
  reactiveVals$call_list[[5]] <- '# Choose option 1: Restore the project from the lockfile. [...] Do you want to proceed? [Y/n]: Y'
  reactiveVals$call_list[[6]] <- substitute(renv::init())
  reactiveVals$call_list[[7]] <- quote(reactiveVals <- list())
  reactiveVals$call_list[[8]] <- substitute(lapply(list.files('functions', full.names = T), source))
  reactiveVals$call_list[[9]] <- ' '
  reactiveVals$call_list[[10]] <- '# Please change paths for input files manually!'
}  