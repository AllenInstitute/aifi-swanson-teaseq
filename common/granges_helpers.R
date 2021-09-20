filter_gr <- function(query_gr,
                      target_gr,
                      mode = "remove",
                      ignore_strand = TRUE) {
  
  overlapping_fragments <- unique(
    S4Vectors::queryHits(
      GenomicRanges::findOverlaps(query_gr,
                                  target_gr,
                                  ignore.strand = ignore_strand)))
  
  if(mode == "remove") {
    query_gr[-overlapping_fragments]
  } else if(mode == "keep") {
    query_gr[overlapping_fragments]
  }
}

convert_fragment_list <- function(fragment_list,
                                  to = c("cut_centers","boundaries","footprints")) {
  out_list <- list()
  if(to == "cut_centers") {
    for(i in 1:length(fragment_list)) {
      prime5 <- fragment_list[[i]]
      BiocGenerics::start(prime5) <- BiocGenerics::start(prime5) + 5
      BiocGenerics::end(prime5) <- BiocGenerics::start(prime5)
      prime3 <- fragment_list[[i]]
      BiocGenerics::end(prime3) <- BiocGenerics::end(prime3) - 4
      BiocGenerics::start(prime3) <- BiocGenerics::end(prime3)
      out_list[[i]] <- GenomicRanges::sort(c(prime5,prime3))
      names(out_list)[i] <- names(fragment_list)[i]
    }
  } else if(to == "boundaries") {
    for(i in 1:length(fragment_list)) {
      prime5 <- fragment_list[[i]]
      BiocGenerics::end(prime5) <- BiocGenerics::start(prime5) + 19
      BiocGenerics::start(prime5) <- BiocGenerics::end(prime5)
      prime3 <- fragment_list[[i]]
      BiocGenerics::start(prime3) <- BiocGenerics::end(prime3) - 19
      BiocGenerics::end(prime3) <- BiocGenerics::start(prime3)
      out_list[[i]] <- GenomicRanges::sort(c(prime5,prime3))
      names(out_list)[i] <- names(fragment_list)[i]
    }
  } else if(to == "footprints") {
    for(i in 1:length(fragment_list)) {
      prime5 <- fragment_list[[i]]
      end(prime5) <- BiocGenerics::start(prime5) + 19
      start(prime5) <- BiocGenerics::start(prime5) -10
      prime3 <- fragment_list[[i]]
      start(prime3) <- BiocGenerics::end(prime3) - 18
      end(prime3) <- BiocGenerics::end(prime3) + 10
      both <- c(prime5,prime3)
      start(both)[start(both) < 1] <- 1
      out_list[[i]] <- BiocGenerics::sort(both)
      names(out_list)[i] <- names(fragment_list)[i]
    }
  }
  out_list
}

pileup_fragments <- function(query_gr,
                            target_gr,
                            target_width = NULL) {
  
  if(is.null(target_width)) {
    target_width <- width(target_gr)[1]
  }
  target_gr <- GenomicRanges::resize(target_gr, 
                                     width = target_width,
                                     fix = "center")
  
  pos_target <- target_gr[strand(target_gr) == "+"]
  neg_target <- target_gr[strand(target_gr) == "-"]
  
  pos_ol <- findOverlaps(query_gr, pos_target, ignore.strand = TRUE)
  neg_ol <- findOverlaps(query_gr, neg_target, ignore.strand = TRUE)
  
  pos_rel <- data.frame(start = start(query_gr)[S4Vectors::queryHits(pos_ol)] - start(pos_target)[S4Vectors::subjectHits(pos_ol)],
                        end = end(query_gr)[S4Vectors::queryHits(pos_ol)] - start(pos_target)[S4Vectors::subjectHits(pos_ol)])
  pos_rel$start[pos_rel$start < 0] <- 0
  pos_rel$end[pos_rel$end > target_width] <- target_width
  
  pos_cov <- integer(target_width + 1)
  for(i in 1:nrow(pos_rel)) {
    pos_cov[(pos_rel$start[i] + 1):(pos_rel$end[i] + 1)] <- pos_cov[(pos_rel$start[i] + 1):(pos_rel$end[i] + 1)] + 1
  }
  
  neg_rel <- data.frame(start = end(neg_target)[S4Vectors::subjectHits(neg_ol)] - end(query_gr)[S4Vectors::queryHits(neg_ol)],
                        end = end(neg_target)[S4Vectors::subjectHits(neg_ol)] - start(query_gr)[S4Vectors::queryHits(neg_ol)])
  neg_rel$start[neg_rel$start < 0] <- 0
  neg_rel$end[neg_rel$end > target_width] <- target_width
  
  neg_cov <- integer(target_width + 1)
  for(i in 1:nrow(neg_rel)) {
    neg_cov[(neg_rel$start[i] + 1):(neg_rel$end[i] + 1)] <- neg_cov[(neg_rel$start[i] + 1):(neg_rel$end[i] + 1)] + 1
  }
  
  data.frame(pos = 0:(target_width),
             cov = pos_cov + neg_cov)
  
}

convert_fragments_gr <- function(fragments,
                                 n_threads = 1) {
  if(n_threads == 1) {
    lapply(fragments,
           function(x) {
             GenomicRanges::GRanges(seqnames = x[["chr"]],
                                    IRanges::IRanges(start = x[["start"]],
                                                     end = x[["end"]]))
           })
  } else {
    starting_order <- names(fragments)
    res <- parallel::mclapply(fragments,
                              function(x) {
                                GenomicRanges::GRanges(seqnames = x[["chr"]],
                                                       IRanges::IRanges(start = x[["start"]],
                                                                        end = x[["end"]]))
                              },
                              mc.cores = n_threads)
    res <- res[starting_order]
    return(res)
  }
}
