ct.format <- function(ct) {
    
    r = gsub(' ', '-', ct)
    r = gsub('\\.', '', r)

}

ct.format.alt <- function(ct) {
    
    r = gsub(' ', '_', ct)
    r = gsub('\\.', '', r)
    
}

create_parent_dir <- function(file_path) {
    
    dir_path <- dirname(file_path)
    
    if (!dir.exists(dir_path)) {
        
        dir.create(dir_path, recursive = TRUE)
        
    }
    
}

# Error handling function

handle_error <- function(e, call.stack, quit = FALSE) {

  error_log <- list(
    "Error message" = conditionMessage(e),
    "Condition call" = if (!is.null(conditionCall(e))) deparse(conditionCall(e)) else "No call available",
    "Call stack" = lapply(call.stack, deparse)
  )
  
  print(error_log)
  
  if (quit) {

    quit(status = 1)

  }

}


assert <- function(condition, message = "Assertion failed") {
  if (!condition) stop(message, call. = FALSE)
}

warn <- function(condition, message = 'Warning!') {
	if (!condition) warning(message, call. = FALSE)
}
