#' Write extremely large datasets to JSON for ingestion by Stan. 
#'
#' @name writeTempJSON
#' @author Jack R. Leary
#' @description In order to bypass R's vector length limits when writing to disk, this function uses chunking to write very large datasets to a temporary JSON file that can then be accessed by \code{CmdStan}. 
#' @param data.list A named list containing the data to be modeled. Defaults to NULL. 
#' @param temp.json The output of \code{tempfile(fileext = ".json")}. Defaults to NULL. 
#' @param chunk.size A double specifying the length above which each element of \code{data.list} will be streamed to the temporary JSON file in chunks. Defaults to one million (\code{1e6}). 
#' @importFrom cli cli_abort
#' @importFrom jsonlite toJSON
#' @return \code{invisible(NULL)}.
#' 

writeTempJSON <- function(data.list = NULL, 
                          temp.json = NULL, 
                          chunk.size = 1e6) {
  if (is.null(data.list) || !inherits(data.list, "list") || is.null(names(data.list))) { cli::cli_abort("You must provide a named data list to writeTempJSON().") }
  if (is.null(temp.json)) { cli::cli_abort("You must provide a temporary JSON file to write to.") }
  con <- file(temp.json, open = "w")
  on.exit(close(con))
  cat("{\n", file = con)
  keys <- names(data.list)
  n_keys <- length(keys)
  for (i in seq_along(keys)) {
    key <- keys[i]
    val <- data.list[[key]]
    cat('  "', key, '": ', sep = "", file = con)
    # stream to temp json in chunks if vector is too long
    if (is.atomic(val) && is.vector(val) && length(val) > chunk.size) {
      cat("[", file = con)
      n <- length(val)
      starts <- seq(1, n, by = chunk.size)
      for (j in seq_along(starts)) {
        start_idx <- starts[j]
        end_idx <- min(start_idx + chunk.size - 1, n)
        chunk_str <- paste(val[start_idx:end_idx], collapse = ",")
        cat(chunk_str, file = con)
        if (end_idx < n) cat(",", file = con)
      }
      cat("]", file = con)
    } else {
      # append to temp json normally if vector is shorter 
      cat(jsonlite::toJSON(val, auto_unbox = TRUE, digits = NA), file = con)
    }
    if (i < n_keys) {
      cat(",\n", file = con) 
    } else { 
      cat("\n", file = con)
    }
  }
  cat("}\n", file = con)
  return(invisible(NULL))
}
