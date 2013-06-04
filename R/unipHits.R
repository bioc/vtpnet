getOneHits = function(tfbsrng, gwr, origin=c("fimo", "gwMatchPWM")) {
   maxsco = max(tfbsrng$score)
   fo = findOverlaps( gwr, tfbsrng )
   ans = gwr[ queryHits(fo) ]
   ans$score = tfbsrng$score[ subjectHits(fo) ]
   ans$tfstart = start(tfbsrng)[ subjectHits(fo) ]
   ans$tfend = end(tfbsrng)[ subjectHits(fo) ]
   if (origin == "fimo") {
     ans$pvalue = as.numeric(tfbsrng$pvalue)[subjectHits(fo)]
     ans$qvalue = as.numeric(tfbsrng$qvalue)[subjectHits(fo)]
     }
   metadata(ans)$maxscore = maxsco
   ans
}


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
   getOneHits( tmp, curgw, origin )
   } 
 names(unipans3) = tags
 unipans3
}
   
