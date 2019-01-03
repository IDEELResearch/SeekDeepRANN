
#' @title SeekDeepDat2ExonAnnotation
#'
#' @description annotates exons
#'
#' @param input A tibble of \class{SeekDeepDat}
#' @param gff path to reference gff
#' @param geneid GeneID that is contained within gff
#' @param ampliconrefseqpath Path for amplicon reference sequence. Dependent on `SeekDeep` file architecture, `SeekDeep` will provide this after extraction under targetRefSeqs
#' @param forwardprimerpath Path for amplicon forward primer sequence. Dependent on `SeekDeep` file architecture, `SeekDeep` will provide this after extraction under targetRefSeqs
#' @param reverseprimerpath Path for amplicon reverse primer sequence. Dependent on `SeekDeep` file architecture, `SeekDeep` will provide this after extraction under targetRefSeqs
#' @param ncbigeneticcode NCBI amino acid coding scheme
#'
#' @details This function makes several assumptions: 1\) Assumes that after you have trimmed your primers, you are only in an exonic region; 2\) No upstream frameshift mutations; 3\) geneid is within gff
#'
#' @return
#'
#' @export

SeekDeepDat2ExonAnnotation <- function(input,
                                       gff, geneid,
                                       ampliconrefseqpath, forwardprimerpath, reverseprimerpath,
                                       ncbigeneticcode=1){


  # error handle
  if(!is.SeekDeepDat(input)){
    stop("Input must be of class SeekDeepDat See the SeekDeepOutput2SeekDeepDat function.")
  }


  #####  READ IN GFF GENE INFORMATION
  geneid.end <- grep('##FASTA',readLines(gff))
  geneid.gff <- read.delim(file = gff,
                           nrows = geneid.end-1, comment= "#", header=F)
  colnames(geneid.gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "info")


  geneid.gff <- subset(geneid.gff, geneid.gff$feature == "gene") # subset to only genes, don't want the other mRNA, etc data

  ############################################
  ##### EXTRACT GENEID from GFF INFO   #######
  ############################################
  geneid.gff$GeneID <- stringr::str_split_fixed(geneid.gff$info, ";", n=2)[,1] # give it to columns to parse on
  geneid.gff$GeneID <- gsub("ID=", "", geneid.gff$GeneID, fixed=T)
  # These are the gene identifiers that we care about

  ############################################
  #####    Subset Gene from GFF        #######
  ############################################
  geneid.gff <- data.frame(geneid.gff[geneid.gff$GeneID == geneid, ], stringsAsFactors = F)

  # read in gff fasta
  gffseq <- seqinr::read.fasta(file=gff, seqtype = c("DNA"), forceDNAtolower = F, strip.desc = T) # of note this is a hacky solution...review it and improve

  # subset
  gffseq <- gffseq[names(gffseq) %in% c(geneid)]



  ####################################
  #####     REFERNET FASTA     #######
  ####################################
  ampliconrefseq <- Biostrings::readDNAStringSet(filepath = ampliconrefseqpath, format="fasta")
  ampliconrefseq <- ampliconrefseq[grepl("pf3d7", tolower(names(ampliconrefseq)))]


  ############################################
  #### Now Trim Forward and Reverse Primer  ##
  ############################################
  Lprimer <- Biostrings::readDNAStringSet(filepath=forwardprimerpath, format = "fasta")
  Rprimer <- Biostrings::reverseComplement(Biostrings::readDNAStringSet(filepath=reverseprimerpath, format = "fasta"))

  Pf3D7haplotypeRef <- Biostrings::trimLRPatterns(Lpattern = Lprimer[[1]],
                                                  Rpattern = Rprimer[[1]],
                                                  subject = ampliconrefseq)  # trim off primers
  ### SANITY CHECK
  ampliconfull <- ampliconrefseq@ranges@width
  ampliconprimertrim <- ampliconfull - Rprimer@ranges@width - Lprimer@ranges@width
  if(ampliconprimertrim != Pf3D7haplotypeRef@ranges@width){
    stop("The amplicon and haplotype are of different lengths. There was an issue in primer triming. Check your primers.")
  }

  #### THIS IS NOW THE REFERENT HAPLOTYPE (primers trimmed that we can compare/align to for variants)

  ############################################
  ####    Make a Dataframe of Variants     ###
  ############################################

  skdpconsens_dnalist <- split(input, factor(input$h_popUID))

  compare_DNA <- function(input,dnastringobject){
    x <- as.integer(Biostrings::DNAString(input$c_Consensus[1]))
    y <- as.integer(Biostrings::DNAString(dnastringobject))

    SNPpos <- which(x != y) #https://www.biostars.org/p/16880/
    SNP <- seqinr::s2c(input$c_Consensus[1])[which(x != y)] # nucleotide bp that are associated with snps

    SNPdf <- data.frame(h_popUID=rep(input$h_popUID[1], length(SNPpos)),
                        SNPpos=SNPpos, SNP=SNP)
    return(SNPdf)   # a single consensus cluster may have several variants -- need to account for this --> will do left_join
  }

  # find snps
  skdpconsens_SNPlist <- lapply(X=skdpconsens_dnalist, FUN=compare_DNA, Pf3D7haplotypeRef[[1]])

  # call snps
  skdpconsens_SNP <- do.call("rbind", skdpconsens_SNPlist)

  if(nrow(skdpconsens_SNP) !=0){   # error handling if no snps in amplicon



    # Identify First base in Amplicon  from "global" Gene Fasta
    RefSeqGenePos <- Biostrings::matchPattern(pattern = Pf3D7haplotypeRef[[1]],
                                              subject = Biostrings::DNAString(seqinr::c2s(gffseq[[1]]))) # find pos of amplicon in gene


    # find starting gene pos
    skdpconsens_SNP$GenePos <- skdpconsens_SNP$SNPpos + RefSeqGenePos@ranges@start - 1 # minus one here so the first base doesn't get counted twice in our amp count and our vcf count
    # Referent Amplicon w/in Gene (from the GFF gene-specific fasta sequence, the strand orientation is already taken care of)
    AArefseq <- seqinr::translate(gffseq[[1]], sens="F", numcode=ncbigeneticcode)
    # Make Mutant Amplicon Haplotype
    haplist <- split(skdpconsens_SNP, f=factor(skdpconsens_SNP$h_popUID))

    mutatehap_possense <- function(haplistobj){
      mutseq <- gffseq[[1]]
      mutseq[haplistobj$GenePos] <- as.character(haplistobj$SNP) # make mutation seq by putting in all SNPs
      AAmutseq <- seqinr::translate(mutseq, sens="F", numcode=ncbigeneticcode) # translate

      haplistobj$CODON <- ceiling(haplistobj$GenePos/3) # ceiling will so that 1/3, 2/3, 3/3 all got to #1 CODON
      haplistobj$AAREF <- AArefseq[haplistobj$CODON]
      haplistobj$AAALT <- AAmutseq[haplistobj$CODON]
      haplistobj$MUT_Type <- ifelse(haplistobj$AAREF == haplistobj$AAALT, "Syn", "Nonsyn")
      return(haplistobj)
    }
    skdpconsens_SNP <- do.call("rbind", lapply(haplist, mutatehap_possense))


    # return
    skdpvcf <- dplyr::left_join(input, skdpconsens_SNP, by=c("h_popUID"))

  } else {
    skdpvcf <- dplyr::left_join(input, skdpconsens_SNP)
    skdpvcf$GenePos <- NA
    skdpvcf$CODON <- NA
    skdpvcf$AAREF <- NA
    skdpvcf$AAALT <- NA
    skdpvcf$MUT_Type <- NA
  } # will just return empty AA change table

  return(skdpvcf)
}
