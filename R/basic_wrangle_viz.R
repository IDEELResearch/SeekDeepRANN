#' @title skdp_filter_simplifier
#'
#' @description Drops clusters without a lot of read support and recalculates haplotype/cluster fractions
#'
#' @param skdpclustinfo_df From SeekDeep Process Clusters -- selectedClustersInfo.tab.txt
#' @param  readcountcutoff: Read Depth Cutoff
#' @return Filtered tibble with columns subsetted to {s_Sample}, {c_AveragedFrac_adj}, {h_popUID}, {c_Consensus}, {filtered_c_ReadCnt_denom}, {c_ReadCnt}
#'
#' @export

SeekDeepOutput2SeekDeepDat <- function(input, readcountcutoff = 0){




  # Filter
  input_simp <- input[input$c_ReadCnt > readcountcutoff, ] # filter at the cluster level (which is the within sample haplotype)

  # Error handling
  if(nrow(input_simp) > 0){
    filtered_c_ReadCnt_sum <- aggregate(input_simp$c_ReadCnt ~ input_simp$s_Sample, function(x){sum(x)}, data = input_simp) # get denominator
    colnames(filtered_c_ReadCnt_sum) <- c("s_Sample", "filtered_c_ReadCnt_denom")
    input_simp <-  dplyr::left_join(input_simp, filtered_c_ReadCnt_sum, by=c("s_Sample"))



    input_simp$c_AveragedFrac_adj <- input_simp$c_ReadCnt/input_simp$filtered_c_ReadCnt_denom # adjusted average fraction by cluster


    input_simp <- input_simp[, c("s_Sample", "c_AveragedFrac_adj", "h_popUID", "c_Consensus", "filtered_c_ReadCnt_denom", "c_ReadCnt")] # keep specific columns

  } else{
    stop("There was an error filtering the reads. Contact the developer")
  }

  input_simp <- as.tibble(input_simp)
  class(input_simp) <- append(  class(input_simp) , "SeekDeepDat" )
  return(input_simp)


}


#---------------------------------------------------------------------------------

#' @title SeekDeepDat2HapPlotter
#'
#' @description plots
#'
#' @param input from skdp_filter_simplifier
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw scale_fill_manual stat
#' @export

SeekDeepDat2HapPlotter <- function(input, target="Target"){

  # error handle
  if(!is.SeekDeepDat(input)){
    stop("Input must be of class SeekDeepDat See the SeekDeepOutput2SeekDeepDat function.")
  }


  # Color setup
  # stackoverflow, # https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
  n <- length(unique(input$h_popUID))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] # pul out qualitative paletes
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col=sample(col_vector, n)


  # Drops clusters without a lot of read support and recalculates haplotype/cluster fractions.
  #
  # Args:
  #   skdpclustinfo_df_simp: Plots output from above
  #
  # Returns:
  #  stacked ggplot



  input %>%
    ggplot() +
    geom_bar(aes(y = c_AveragedFrac_adj, x = s_Sample, fill = h_popUID), stat="identity") +
    scale_fill_manual(values = col) +
    theme_bw()

}

