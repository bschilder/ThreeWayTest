#' Get remote data
#'
#' Download remotely stored data via \link[piggyback]{pb_download}. 
#' @keywords internal
#' @inheritParams piggyback::pb_download
#' @importFrom tools R_user_dir
get_data <- function(file,
                     repo = "bschilder/ThreeWayTest",
                     tag = "latest",
                     storage_dir = tools::R_user_dir(
                       package = "ThreeWayTest",
                       which = "cache"
                     ),
                     overwrite = TRUE,
                     .token = gh::gh_token(),
                     check = FALSE
                     ){
  requireNamespace("piggyback")
  requireNamespace("gh")
  
  tmp <- file.path(storage_dir, file)
  Sys.setenv("piggyback_cache_duration" = 10)
  dir.create(storage_dir, showWarnings = FALSE, recursive = TRUE) 
  if(as.character(.token)=="")  {
    message(
        paste0(
            "Personal access token is missing and it is required to download ",
            names(file), "\n"))
    message("Please, configure git with R and try again \n")
    stop("Personal access token missing missing")
    } else {
      piggyback::pb_download(
        file = file,
        dest = storage_dir,
        repo = repo,
        overwrite = overwrite,
        .token = .token
      )
    } 
  return(tmp)
}


