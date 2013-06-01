setupGWPM = function(cores=11) {
   require(foreach)
   require(doParallel)
   registerDoParallel(cores=14)
   require(BSgenome.Hsapiens.UCSC.hg19)
   require(MotifDb)
}

gwMatchPWM = function( mdb1, tag, bsg=Hsapiens, nchr=24, thresh="75%" ) {
#
# assumes you only want to process sequences 1:nchr in bsg
#
    if (!is(mdb1, "MotifList")) stop("operates only on MotifList instances")
    if (length(mdb1)>1) stop("input must be of MotifList of length 1")
    if (missing(tag)) stop("a tag (character string) must be supplied to label metadata on hits")
    sn = seqnames(bsg)[1:nchr]
    ans = foreach(i=1:nchr) %dopar% {
#     cat(i)
     subj = bsg[[ sn[i] ]]
     hits = matchPWM( mdb1[[1]], subj, thresh )
     scores = PWMscoreStartingAt( mdb1[[1]], subject(hits), start(hits) )
     GRanges( seqnames=sn[i], ranges=as(hits, "IRanges"), tag=rep(tag,length(hits)),
         score=scores)
    }
    do.call(c, ans)
   }

getMtags = function(mdb) {
   if (!is(mdb, "MotifList")) stop("requires MotifList input")
   utags = make.unique(values(mdb)$geneSymbol)
   longtags = names(mdb)
   list(utags = utags, longtags = longtags)
}

MdbSearchToDisk = function(mdb, thresh="80%", cores=14) {
   setupGWPM(cores=cores)
   tagSet = getMtags(mdb)
   obn = tagSet$utags
   fn = paste0(obn, ".rda")
  
   for (i in 1:length(obn)) {
#    cat(i)
    assign( obn[i], gwMatchPWM( mdb[i], obn[i], thresh=thresh ) )
    save(list=obn[i], file=fn[i] )
    gc()
   }
}
