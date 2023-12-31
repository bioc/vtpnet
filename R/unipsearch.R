setupGWPM = function(cores=11) {
   require(foreach)
   require(doParallel)
   registerDoParallel(cores=14)
   require(MotifDb)
}

matrc = function(x) {
  x[4:1, ncol(x):1]
}

gwMatchPWM = function( mdb1, tag, bsg=Hsapiens, nchr=24, thresh="75%", matTx=force ) {
#
# assumes you only want to process sequences 1:nchr in bsg
# bsg is any object answering seqnames and [[ seqname[i] ]] wih DNAString
#
# use matTx=matrc to match the reverse complement model
#
    if (!is(mdb1, "MotifList")) stop("operates only on MotifList instances")
    if (length(mdb1)>1) stop("input must be of MotifList of length 1")
    if (missing(tag)) stop("a tag (character string) must be supplied to label metadata on hits")
    sn = seqnames(bsg)[1:nchr]
    ans = foreach(i=1:nchr) %dopar% {
#     cat(i)
     subj = bsg[[ sn[i] ]]
     hits = matchPWM( matTx(mdb1[[1]]), subj, thresh )
     scores = PWMscoreStartingAt( matTx(mdb1[[1]]), subject(hits), start(hits) )
     GRanges( seqnames=sn[i], ranges=as(hits, "IRanges"), tag=rep(tag,length(hits)),
         score=scores)
    }
    do.call(c, ans)
   }

getMtags = function(mdb, fixer=function(x) gsub("\\/", "_", x)) {
   if (!is(mdb, "MotifList")) stop("requires MotifList input")
   utags = make.unique(values(mdb)$geneSymbol)
   longtags = names(mdb)
   list(utags = fixer(utags), longtags = fixer(longtags))
}

MdbSearchToDisk = function(mdb, bsg, thresh="80%", cores=14, ...) {
   setupGWPM(cores=cores)
   tagSet = getMtags(mdb)
   obn = tagSet$utags
   fn = paste0(obn, ".rda")
  
   for (i in 1:length(obn)) {
#    cat(i)
    assign( obn[i], gwMatchPWM( mdb1=mdb[i], tag=obn[i], bsg=bsg, thresh=thresh, ... ) )
    save(list=obn[i], file=fn[i] )
    gc()
   }
}
