rng2facHits = function(rngsetgen = function() makeCurrentGwascat(), origin=c("fimo", "gwMatchPWM"),
    factorRngs = dir(patt="rda$"), ncore=12) {
#
# factorRngs is collection of GRanges of factor binding regions on disk
# rngsetgen() will provide a set of variants in GRanges whose
#   coincidence with factorRngs is to be tabulated
#
 require(GenomicRanges)
 require(gwascat)
 curgw = as( rngsetgen(), "GRanges" )
 urngs = factorRngs
 tags = gsub(".rda", "", urngs)
 require(doParallel)
 require(foreach)
 registerDoParallel(cores=ncore)
 unipans3 = foreach (i = 1:length(urngs)) %dopar% {
   cat(i)
   tmp = get(load(urngs[i]))
   maxsco = max(tmp$score)
   #try(curgw[ overlapsAny( curgw, tmp ) ])
   fo = findOverlaps( curgw, tmp )
   ans = curgw[ queryHits(fo) ]
   ans$score = tmp$score[ subjectHits(fo) ]
   if (origin == "fimo") {
     ans$pvalue = as.numeric(tmp$pvalue)[subjectHits(fo)]
     ans$qvalue = as.numeric(tmp$qvalue)[subjectHits(fo)]
     }
   metadata(ans)$maxscore = maxsco
   ans
   } 
 names(unipans3) = tags
 unipans3
}
   
