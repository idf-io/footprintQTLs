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
handle_error <- function(err) {
    cat("Error:", conditionMessage(err), "\n")
    quit(status = 1)
}


