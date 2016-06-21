## internal sequencing function

sgseq = function(readmat, transcripts, paired, outdir, extras, weightsOnly=FALSE){

  # new errors
  stopifnot(length(extras$fraglen) == ncol(readmat))
  stopifnot(length(extras$fragsd) == ncol(readmat))
  if (!is.null(extras$fragGCBias)) stopifnot(length(extras$fragGCBiasData) == ncol(readmat))
  # end new errors

  weightsOut <- numeric(ncol(readmat))
  
    for(i in 1:ncol(readmat)){

        tObj = rep(transcripts, times=readmat[,i])
        
        #get fragments

        ### new code: iterate over fraglen, fragse and fraggcbias ###
        tFrags = generate_fragments(tObj, extras$fraglen[i], extras$fragsd[i],
            extras$readlen, extras$distr, extras$custdens, extras$bias,
            extras$fragGCBias, extras$fragGCBiasData[[i]])
        # keep track of how many missing fragments there are
        # weights is the amount we need to increase library size
        # to get constant expected library size after fragGCBias
        if (weightsOnly) { weightsOut[i] <- length(tObj)/length(tFrags) }
        ### end new code ###

        if (!weightsOnly) {
          #reverse_complement some of those fragments
          rctFrags = reverse_complement(tFrags)
          
          #get reads from fragments
          reads = get_reads(rctFrags, extras$readlen, paired)

          #add sequencing error
          if(extras$error_model == 'uniform'){
            errReads = add_error(reads, extras$error_rate)
          }else if(extras$error_model == 'custom'){
             errReads = add_platform_error(reads, 'custom', paired, extras$path)
           }else{
              errReads = add_platform_error(reads, extras$error_model, paired)
            }

          #write read pairs
          write_reads(errReads, readlen=extras$readlen, 
                      fname=paste0(outdir, '/sample_', sprintf('%02d', i)), paired=paired)
        }
    }
  
  return(weightsOut)
}
