runfimo = function( mdb, parmvec = " --max-stored-scores 10000000 ", fasta="/home/stvjc/hg19.fa" ) {
  tags = getMtags(mdb)
  jnk = foreach (i = 1:length(tags[[1]])) %dopar% {
    fn = paste0(tags[[1]][i], ".meme")
    odir = paste0(tags[[1]][i], "_out")
    MotifDb::export( mdb[tags[[2]][i]], con=fn, "meme")
    comm = paste0("fimo ", parmvec, " --oc ", odir, " ", fn, " ", fasta)
    print(comm)
    system(comm)
  }
}
