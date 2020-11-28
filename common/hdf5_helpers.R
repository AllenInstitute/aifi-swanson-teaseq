#' Decide how many bits to use to store integer values in an HDF5 Dataset object.
#'
#' @param x a numeric vector
#'
#' @return an integer value. One of 16L, 32L, 64L, or NULL if larger than 64-bit.
#' @export
#'
choose_integer_bits <- function(x) {
  max_val <- max(x, na.rm = TRUE)
  
  if(max_val < 2^16) {
    16L
  } else if(max_val < 2^32) {
    32L
  } else if(max_val < 2^64) {
    64L
  } else {
    NULL
  }
  
}

#' Select a reasonable chunk size for an HDF5 dataset object
#'
#' If x has length > 1e6, chunks will be 1e5.
#' If x has length > 1e2, chunks will be one log10 smaller than the length of x.
#' If x has length <= 1e2, chunk size is the length of x.
#'
#' @param x A vector to store in an HDF5 file
#'
#' @return a numeric value containing a suggested chunk size
#' @export
#'
choose_chunk_size <- function(x) {
  x_len <- length(x)
  x_logs <- floor(log10(x_len))
  
  if(x_logs > 6) {
    10^5
  } else if (x_logs > 2) {
    10^(x_logs - 1)
  } else {
    length(x)
  }
}

#' Convert all 1D Arrays in a list object to vectors recursively
#'
#' @param x a list object
#'
#' @return a list object with all 1D arrays converted.
#' @export
#'
strip_1d_array_recursive <- function(x) {
  assertthat::assert_that(class(x) == "list")
  
  if(length(x) > 0) {
    for(n in seq_along(x)) {
      if(class(x[[n]]) == "list") {
        x[[n]] <- strip_1d_array_recursive(x[[n]])
      } else if(class(x[[n]]) == "array" & length(dim(x[[n]] == 1))) {
        x[[n]] <- as.vector(x[[n]])
      }
    }
  }
  
  x
}


#' Convert all "NA" character values to actual NAs recursively
#'
#' @param x a list object
#'
#' @return a list object with "NA"s converted to NA
#' @export
#'
convert_char_na_recursive <- function(x) {
  assertthat::assert_that(class(x) == "list")
  
  if(length(x) > 0) {
    for(n in seq_along(x)) {
      if(class(x[[n]]) == "list") {
        x[[n]] <- convert_char_na_recursive(x[[n]])
      } else if(class(x[[n]]) == "character") {
        x[[n]] <- convert_char_na(x[[n]])
      }
    }
  }
  
  x
}

#' Generate a set of attributes based on 10x Genomics defaults
#'
#' @param library_ids Library IDs to store as attributes. if NULL (default), will generate a random ID using date and ids::proquint().
#'
#' @return a list of attributes
#' @export
#'
h5_attr_list <- function(library_ids = NULL) {
  
  attr_list <- list(chemistry_description = "Single Cell 3' v3",
                    filetype = "matrix",
                    library_ids = paste(Sys.Date(),
                                        ids::proquint(1),
                                        sep = "-"),
                    original_gem_groups = 1,
                    version = 2)
  if(!is.null(library_ids)) {
    assertthat::assert_that(class(library_ids) == "character")
    attr_list$library_ids <- library_ids
  }
  
  attr_list
}

#' Write an h5_list, as created by rhdf5::h5dump(), to an .h5 file
#'
#' @param h5_list a list object, e.g. a list created by rhdf5::h5dump()
#' @param h5_file a character object specifying the location of a .h5 file to write to.
#' @param overwrite a logical value specifying whether or not to overwrite an existing .h5 file. Default is FALSE.
#' @param h5_handle an existing h5_handle created by H5Fopen(). Used for recursion. The default (NULL) should usually be used.
#' @param h5_target a base location within the HDF5 file to write to. Mainly used for recursion. The default ("/") should usually be used.
#' @param h5_attributes a list of attributes to add to an .h5 file to try to imitate 10x Genomics outputs. If NULL (default), will be skipped. "tenx" uses in-built data from 'h5_attr_list()'.
#' @param library_ids a character vector of library ids to add to attributes. Only used if h5_attributes != NULL. Default is NULL.
#'
#' @return Writes a file; no return to R.
#' @export
#'
write_h5_list <- function(h5_list,
                          h5_file,
                          overwrite = FALSE,
                          addon = FALSE,
                          h5_handle = NULL,
                          h5_target = "/",
                          h5_attributes = NULL,
                          library_ids = NULL) {
  
  assertthat::assert_that(is.list(h5_list))
  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)
  
  if(!is.null(h5_attributes)) {
    if(h5_attributes  == "tenx") {
      h5_attributes <- h5_attr_list()
    }
  }
  
  if(!is.null(library_ids)) {
    h5_attributes$library_ids <- library_ids
  }
  
  # Make sure the HDF5 file connection is closed if the function
  # exits due to an error.
  on.exit(expr = {
    if(h5_target == "/") {
      rhdf5::h5closeAll()
    }
  })
  
  if(is.null(h5_handle)) {
    if(!overwrite & !addon) {
      # If we aren't overwriting or adding, check for file and halt if it exists.
      # If it doesn't exist, create it.
      
      if(file.exists(h5_file)) {
        stop(paste(h5_file, "already exists."))
      } else {
        rhdf5::H5Fcreate(h5_file)
      }
      
    } else if(overwrite) {
      # If we're overwriting, remove the old file and make a new one.
      # Partial changes are performed by addon
      
      # If overwrite and addon, overwrite takes precedence.
      
      if(file.exists(h5_file)) {
        file.remove(h5_file)
      }
      
      rhdf5::H5Fcreate(h5_file)
      
    } else if(addon) {
      # If addon, we don't need to create or remove the file.
      # Keeping this condition to help with reasoning even though
      # nothing happens here.
      # We'll open below.
    } else {
      rhdf5::H5Fcreate(h5_file)
    }
    
    h5_handle <- rhdf5::H5Fopen(h5_file)
  }
  
  # Add file attributes to match cellranger
  if(!is.null(h5_attributes) & !addon) {
    if(h5_target != "/" & class(h5_attributes) == "list") {
      base_obj <- rhdf5::H5Dopen(h5_handle, "/")
      
      for(i in 1:length(h5_attributes)) {
        rhdf5::h5writeAttribute(h5_attributes[[i]],
                                base_obj,
                                names(h5_attributes)[i])
      }
      
      rhdf5::H5Dclose(base_obj)
    }
    h5_attributes <- NULL
  }
  
  
  h5_names <- names(h5_list)
  if(length(h5_names) > 0) {
    for(h5_name in h5_names) {
      new_object <- paste0(h5_target, h5_name)
      
      # Correct 1d arrays to vectors for storage
      if(class(h5_list[[h5_name]]) == "array") {
        if(length(dim(h5_list[[h5_name]])) == 1) {
          h5_list[[h5_name]] <- as.vector(h5_list[[h5_name]])
        }
      }
      
      # Correct numeric arrays to integers if they're all whole numbers
      if(class(h5_list[[h5_name]]) == "numeric") {
        
        x <- as.integer(h5_list[[h5_name]])
        if(isTRUE(all.equal(h5_list[[h5_name]], x, check.attributes = FALSE))) {
          h5_list[[h5_name]] <- x
        }
      }
      
      if(class(h5_list[[h5_name]]) == "list") {
        rhdf5::h5createGroup(h5_handle,
                             group = new_object)
        # Recurse function to write children of list
        write_h5_list(h5_list[[h5_name]],
                                h5_file = h5_file,
                                h5_handle = h5_handle,
                                h5_target = paste0(new_object,"/"),
                                h5_attributes = h5_attributes,
                                addon = addon)
      } else if(class(h5_list[[h5_name]]) == "numeric") {
        rhdf5::h5createDataset(h5_handle,
                               dataset = new_object,
                               dims = list(length(h5_list[[h5_name]])),
                               chunk = choose_chunk_size(h5_list[[h5_name]]),
                               storage.mode = storage.mode(h5_list[[h5_name]]))
        
        rhdf5::h5write(obj = h5_list[[h5_name]],
                       file = h5_handle,
                       name = new_object)
        
      } else if(class(h5_list[[h5_name]]) == "integer") {
        bits <- choose_integer_bits(h5_list[[h5_name]])
        h5_type <- paste0("H5T_NATIVE_UINT",bits)
        
        rhdf5::h5createDataset(h5_handle,
                               dataset = new_object,
                               dims = list(length(h5_list[[h5_name]])),
                               chunk = choose_chunk_size(h5_list[[h5_name]]),
                               H5type = h5_type)
        
        rhdf5::h5write(obj = h5_list[[h5_name]],
                       file = h5_handle,
                       name = new_object)
        
      } else if(class(h5_list[[h5_name]]) == "character") {
        h5_list[[h5_name]][is.na(h5_list[[h5_name]])] <- "NA"
        
        rhdf5::h5createDataset(h5_handle,
                               dataset = new_object,
                               dims = list(length(h5_list[[h5_name]])),
                               chunk = choose_chunk_size(h5_list[[h5_name]]),
                               storage.mode = storage.mode(h5_list[[h5_name]]),
                               size = max(nchar(h5_list[[h5_name]])) + 1)
        
        rhdf5::h5write(obj = h5_list[[h5_name]],
                       file = h5_handle,
                       name = new_object)
      }
    }
  }
  
}

#' Read the /matrix from a .h5 file as a sparse matrix
#'
#' @param h5_file the path to an .h5 file in 10x Genomics format
#' @param target a character object specifying the target matrix within the file. Default is "matrix".
#' @param feature_names a character object specifying whether to use "id" or "name" for row.names. Default is "id".
#' @param sample_names a character object specifying which values to use for col.names. If "barcodes", will use /target/barcodes. Other values will be read from /target/observations/
#' @param index1 a logical object specifying whether index vectors should start with 0 (FALSE) or 1 (TRUE). Default is TRUE.
#'
#' @return a dgCMatrix of gene expression values.
#' @export
#'
read_h5_dgCMatrix <- function(h5_file,
                              target = "matrix",
                              feature_names = "id",
                              sample_names = "barcodes",
                              index1 = TRUE) {
  
  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)
  
  assertthat::assert_that(is.character(target))
  assertthat::assert_that(length(target) == 1)
  
  if(grepl("^/",target)) {
    target <- sub("^/","",target)
  }
  
  assertthat::assert_that(is.character(feature_names))
  assertthat::assert_that(length(feature_names) == 1)
  
  assertthat::assert_that(is.character(sample_names))
  assertthat::assert_that(length(sample_names) == 1)
  
  if(!file.exists(h5_file)) {
    stop(paste(h5_file, "does not exist."))
  }
  
  # Make sure the HDF5 file connection is closed if the function
  # exits due to an error.
  on.exit(expr = {
    rhdf5::H5Fclose(h5_handle)
  })
  
  feature_names <- match.arg(arg = feature_names,
                             choices = c("id","name"))
  
  h5_handle <- rhdf5::H5Fopen(h5_file)
  
  if(sample_names == "barcodes") {
    colname_target <- paste0("/", target, "/barcodes")
  } else {
    colname_target <- paste0("/", target, "/observations/", sample_names)
  }
  
  if(index1) {
    mat <- Matrix::sparseMatrix(x = rhdf5::h5read(h5_handle, paste0("/",target,"/data")),
                                i = rhdf5::h5read(h5_handle, paste0("/",target,"/indices")) + 1,
                                p = rhdf5::h5read(h5_handle, paste0("/",target,"/indptr")),
                                index1 = index1,
                                dims = rhdf5::h5read(h5_handle, paste0("/",target,"/shape")),
                                dimnames = list(as.vector(rhdf5::h5read(h5_handle, paste0("/",target,"/features/",feature_names))),
                                                as.vector(rhdf5::h5read(h5_handle, colname_target))
                                )
    )
  } else {
    mat <- Matrix::sparseMatrix(x = rhdf5::h5read(h5_handle, paste0("/",target,"/data")),
                                i = rhdf5::h5read(h5_handle, paste0("/",target,"/indices")),
                                p = rhdf5::h5read(h5_handle, paste0("/",target,"/indptr")),
                                index1 = index1,
                                dims = rhdf5::h5read(h5_handle, paste0("/",target,"/shape")),
                                dimnames = list(as.vector(rhdf5::h5read(h5_handle, paste0("/",target,"/features/",feature_names))),
                                                as.vector(rhdf5::h5read(h5_handle, colname_target))
                                )
    )
  }
  
  mat
}

#' Read .h5 Cell Metadata
#'
#' @param h5_file the path to an .h5 file in 10x Genomics format
#' @param target A matrix object in the .h5 file with a /barcodes object and/or a /target/observations/ sub-group. Default is "matrix".
#'
#' @return A data.frame containing all feature metadata found in /target/barcodes and /target/observations/
#' @export
#'
read_h5_cell_meta <- function(h5_file,
                              target = "matrix") {
  
  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)
  
  target <- ifelse(grepl("^/",target),
                   target,
                   paste0("/",target))
  
  h5_contents <- h5ls(h5_file)
  target_contents <- h5_contents[grepl(paste0("^",target), h5_contents$group),]
  
  h5_meta_targets <- character()
  
  target_bcs <- paste0(target, "/barcodes")
  
  if(target_bcs %in% target_contents$full_name) {
    h5_meta_targets <- c(h5_meta_targets,
                         target_bcs)
  }
  
  target_obs <- paste0(target, "/observations")
  
  if(target_obs %in% target_contents$full_name) {
    h5_meta_targets <- c(h5_meta_targets,
                         target_contents$full_name[target_contents$group == target_obs])
  }
  
  if(length(h5_meta_targets) > 0) {
    meta_list <- lapply(h5_meta_targets,
                        function(h5_meta_target) {
                          rhdf5::h5read(h5_file,
                                        h5_meta_target)
                        })
    rhdf5::h5closeAll()
    
    names(meta_list) <- sub(".+/","",h5_meta_targets)
    
    meta_list <- strip_1d_array_recursive(meta_list)
    meta_list <- convert_char_na_recursive(meta_list)
    
    df <- as.data.frame(meta_list,
                        stringsAsFactors = FALSE)
    
    df
  } else {
    stop("No cell metadata found in h5_file.")
  }
}

#' Read .h5 Feature Metadata
#'
#' @param h5_file the path to an .h5 file in 10x Genomics format
#' @param target A matrix object in the .h5 file with a /features/ sub-group. Default is "matrix".
#'
#' @return a data.frame containing all feature metadata found in /target/features/
#' @export
read_h5_feature_meta <- function(h5_file,
                                 target = "matrix") {
  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)
  
  target <- ifelse(grepl("^/",target),
                   target,
                   paste0("/",target))
  
  h5_contents <- h5ls(h5_file)
  target_contents <- h5_contents[grepl(paste0("^",target), h5_contents$group),]
  
  h5_meta_targets <- character()
  
  target_feat <- paste0(target, "/features")
  
  if(target_feat %in% target_contents$full_name) {
    h5_meta_targets <- c(h5_meta_targets,
                         target_contents$full_name[target_contents$group == target_feat])
  }
  
  h5_meta_names <- sub(".+/","",h5_meta_targets)
  h5_meta_targets <- h5_meta_targets[!grepl("^_",h5_meta_names)]
  h5_meta_names <- h5_meta_names[!grepl("^_",h5_meta_names)]
  
  if(length(h5_meta_targets) > 0) {
    meta_list <- lapply(h5_meta_targets,
                        function(h5_meta_target) {
                          rhdf5::h5read(h5_file,
                                        h5_meta_target)
                        })
    rhdf5::h5closeAll()
    
    names(meta_list) <- h5_meta_names
    
    meta_list <- strip_1d_array_recursive(meta_list)
    meta_list <- convert_char_na_recursive(meta_list)
    
    df <- as.data.frame(meta_list,
                        stringsAsFactors = FALSE)
    
    df
  } else {
    stop("No cell metadata found in h5_file.")
  }
}

#' Convert the matrix in an h5_list from 10x Genomics data to a sparse matrix
#'
#' This is very useful for subsetting the matrix based on column names.
#'
#' @param h5_list a list object generated by running rhdf5::h5dump() on a 10x HDF5 file.
#' @param target the group name for the sparse matrix. Default is "matrix", which is used by 10x.
#' @param feature_names the group name in target/features to use for rownames. Must be either "id" or "name". Default is "id".
#'
#' @return a list object
#' @export
#'
h5_list_convert_to_dgCMatrix <- function(h5_list,
                                         target = "matrix",
                                         feature_names = "id") {
  
  assertthat::assert_that(class(h5_list) == "list")
  assertthat::assert_that(target %in% names(h5_list))
  assertthat::assert_that(feature_names %in% c("id","name"))
  
  target_dgCMatrix <- paste0(target,"_dgCMatrix")
  
  h5_list[[target_dgCMatrix]] <- Matrix::sparseMatrix(x = h5_list[[target]]$data,
                                                      i = h5_list[[target]]$indices,
                                                      index1 = FALSE,
                                                      p = h5_list[[target]]$indptr,
                                                      dims = h5_list[[target]]$shape,
                                                      dimnames = list(h5_list[[target]]$features[[feature_names]],
                                                                      h5_list[[target]]$barcodes))
  h5_list[[target]]$data <- NULL
  h5_list[[target]]$indices <- NULL
  h5_list[[target]]$indptr <- NULL
  h5_list[[target]]$shape <- NULL
  h5_list[[target]]$features$id <- NULL
  h5_list[[target]]$barcodes <- NULL
  
  h5_list
}

#' Convert the matrix in an h5_list from a sparse matrix back to its original structure
#'
#' This is very useful for preparing to write back out to an .h5 file.
#'
#' @param h5_list a list object generated by running rhdf5::h5dump() on a 10x HDF5 file, with a matrix converted by h5_list_convert_to_dgCMatrix().
#' @param target the group name of the original sparse matrix. Default is "matrix", which is used by 10x.
#'
#' @return a list object
#' @export
#'
h5_list_convert_from_dgCMatrix <- function(h5_list,
                                           target = "matrix") {
  
  assertthat::assert_that(class(h5_list) == "list")
  
  target_dgCMatrix <- paste0(target,"_dgCMatrix")
  assertthat::assert_that(target_dgCMatrix %in% names(h5_list))
  
  h5_list[[target]]$data <- h5_list[[target_dgCMatrix]]@x
  h5_list[[target]]$indices <- h5_list[[target_dgCMatrix]]@i
  h5_list[[target]]$indptr <- h5_list[[target_dgCMatrix]]@p
  h5_list[[target]]$shape <- dim(h5_list[[target_dgCMatrix]])
  h5_list[[target]]$barcodes <- colnames(h5_list[[target_dgCMatrix]])
  h5_list[[target]]$features$id <- rownames(h5_list[[target_dgCMatrix]])
  
  h5_list[[target_dgCMatrix]] <- NULL
  
  h5_list
}
