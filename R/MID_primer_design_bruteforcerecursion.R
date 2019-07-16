# MID
# not exported

MIDmaker <- function(length, homopolymercutoff){

  MID <- paste0(c(sample(x=c("C","G"), size=1), sample(x=c("C","G", "A", "T"), size=length-1, replace = T)), collapse="") # can't start with A or T

  while(max(unlist(parallel::mclapply(c("A", "C", "T", "G"), function(x){longestConsecutive(MID, x)}))) >= homopolymercutoff | # avoid homopolymers
        length(Biostrings::findPalindromes(Biostrings::BString(MID), min.armlength=2)) != 0 # avoid palindromes
  ){
    MID <- paste0(c(sample(x=c("C","G"), size=1), sample(x=c("C","G", "A", "T"), size=length-1, replace = T)), collapse="")
  }

  return(MID)

}




#' @title MID Primer Design
#'
#' @author Nick Brazeau
#' @description This will produces MIDs
#' @details This is a brute-force approach with recursion as the number of new MIDs proposed is based on "how close we are to solving the problem" with the current sit of MIDs
#'
#' @param design_fasta is expected to have the forward, reverse, and target gene in it. Must have `forward` and `reverse` in the contig names for the forward and reverse primers
#' @param MID2targetmismatchesAllowed is the number of base-pair mismatches allowed in the MID before it would bind to ANY place on the target gene
#' @importFrom magrittr %>%
#' @export

MIDPrimerFinder <- function(design_fasta = NULL,
                            MID2MIDmatchesAllowed = NULL,
                            MIDlength = NULL,
                            MIDhomopolymerallowance = NULL,
                            MID2targetmismatchesAllowed = NULL,
                            MIDnum = NULL
){

  design_fasta <- Biostrings::readDNAStringSet(filepath = design_fasta, format = "fasta", use.names = T)

  ## Housekeeping
  fw <- design_fasta[ grepl("forward", tolower(names(design_fasta))) ]
  rv <- Biostrings::reverseComplement(design_fasta[ grepl("reverse", tolower(names(design_fasta))) ])
  tg <- design_fasta[ ! tolower(names(design_fasta)) %in% tolower(c(names(fw), names(rv))) ]

  # error handling

  hits <- lapply(list(fw, rv), function(x) {
    Biostrings::vmatchPattern(as.character(x), tg) # coercion here is fine since dnastrinset of length 1 ... becomes a named character vector
  })


  if( any(sapply(hits, function(x){  S4Vectors::elementNROWS(x) }) == 0 ) ){
    stop("Your primers do not bind to your target")
  }

  if( any(sapply(hits, function(x){ S4Vectors::elementNROWS(x) }) > 1 ) ){
    stop("Your primers have more than one binding site on your target")
  }

  if( all(sapply(c(MID2MIDmatchesAllowed, MIDlength, MIDhomopolymerallowance, MID2targetmismatchesAllowed, MIDnum), class) != "numeric") ){
    stop("Your MID specifications all must be of class numeric")
  }

  if( !all( c(MID2MIDmatchesAllowed, MIDlength, MIDhomopolymerallowance, MID2targetmismatchesAllowed, MIDnum) %in% sapply(c(MID2MIDmatchesAllowed, MIDlength, MIDhomopolymerallowance, MID2targetmismatchesAllowed, MIDnum), floor)) ){
    stop(" Your MID specifications must be whole numbers (no decimals)")
  }

  # housekeeping II
  # this is a very strange S4 object. I assume this will not always be backwards compatible...
  primerinfo <-
    tibble::tibble( primer = c("Forward", "Reverse"),
            start = c(
              hits[[1]]@ends[[1]] - hits[[1]]@width0 + 1,
              hits[[2]]@ends[[1]] - hits[[2]]@width0 + 1
            ),
            ends = c(
              hits[[1]]@ends[[1]],
              hits[[2]]@ends[[1]]
            ),
            width = c(
              hits[[1]]@width0,
              hits[[2]]@width0
            )
    )

  # START
  finalMIDs <- NULL # init N
  failedMIDs <- NULL
  n <- MIDnum # init N
  int <- 1 # housekeeping


  while( length(finalMIDs) < n ){

    propMIDs <- replicate((n - length(finalMIDs)), MIDmaker(MIDlength, MIDhomopolymerallowance)) # init


    # housekeeping for prop checks
    rev_culprits <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(propMIDs))) %in% as.character(Biostrings::DNAStringSet(propMIDs))  # catch any reverse complements
    ot_cultprits <- c( sapply(propMIDs, function(x){Biostrings::vcountPattern(x, tg, max.mismatch = MID2targetmismatchesAllowed) }) > 0 )  # catch overlap of MIDs with target sequence
    dup_culprits <- duplicated(propMIDs) # catch duplicated MIDs
    failed_culprits <- propMIDs %in% failedMIDs


    # repropose culprits
    while( sum(rev_culprits) > 0 | sum(ot_cultprits) > 0 | sum(dup_culprits) > 0 | sum(failed_culprits) > 0){

      if(sum(rev_culprits) > 0){
        propMIDs[rev_culprits] <- replicate(sum(rev_culprits), MIDmaker(MIDlength, MIDhomopolymerallowance))
      }

      if(sum(ot_cultprits) > 0){
        propMIDs[ot_cultprits] <- replicate(sum(ot_cultprits), MIDmaker(MIDlength, MIDhomopolymerallowance))
      }

      if(sum(dup_culprits) > 0){
        propMIDs[dup_culprits] <- replicate(sum(dup_culprits), MIDmaker(MIDlength, MIDhomopolymerallowance))
      }

      if(sum(failed_culprits) > 0){
        propMIDs[failed_culprits] <- replicate(sum(failed_culprits), MIDmaker(MIDlength, MIDhomopolymerallowance))
      }

      # housekeeping for CANDIDATE checks
      rev_culprits <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(propMIDs))) %in% as.character(Biostrings::DNAStringSet(propMIDs))  # catch any reverse complements
      ot_cultprits <- c( sapply(propMIDs, function(x){Biostrings::vcountPattern(x, tg, max.mismatch = MID2targetmismatchesAllowed) }) > 0 )  # catch overlap of MIDs with target sequence
      dup_culprits <- duplicated(propMIDs) # catch duplicated MIDs
      failed_culprits <- propMIDs %in% failedMIDs
      # ^ note above is strictly for readability and not memory. I would call this not-best-practice, as it would be easier to evaluate the condition on the fly w/in the while loop...but for clarity



    }




    if(int%%10 == 0){cat(paste0("Iteration  ", int, "\n ------- \n"))}
    # get pairs and pair.list
    pairs <- t(combn(propMIDs, m=2)) # choose 2
    pairs_list <- split(pairs, seq(nrow(pairs)))
    pairs_list_ret <- parallel::mclapply(pairs_list, function(x){stringdist::stringdistmatrix(x[1], x[2], method = c("hamming"))})

    mid.dists.init <- tibble::tibble(primer = c(pairs[,1], pairs[,2]), Hdist = rep(unlist(pairs_list_ret),2)) # TODO sloppy

    # find primers that had distance metrics above MID2MIDmatches for all MIDs currently in set
    # and keep failed MIDs so we don't waste time looking there anymore in our state space
    failedMIDs <- c(failedMIDs, mid.dists.init$primer[ mid.dists.init$Hdist < MIDlength-MID2MIDmatchesAllowed ])
    failedMIDs <- failedMIDs[!duplicated(failedMIDs)]

     # find passed MIDs
    passedMIDs <- propMIDs[! propMIDs %in% failedMIDs]

    finalMIDs <- c(passedMIDs, finalMIDs)

  }



  ret <- list(finalMIDs = finalMIDs,
              primerinfo = primerinfo,
              fwprimer = fw,
              rvprimer = rv)
  return(ret)


}
