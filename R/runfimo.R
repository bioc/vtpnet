
runfimo = function (mdb, parmvec = " --max-stored-scores 10000000 ", fasta = "/home/stvjc/hg19.fa", 
    setuponly = FALSE, commprefix="", continue=FALSE)
{
    tags = getMtags(mdb)
    if (continue) {
        ideal = paste0(tags[[1]], "_out/fimo.gff")
        allout = dir(patt="_out") 
        actual = paste0(allout, "/fimo.gff")
        OK = file.exists(actual)
        todo = which(!(ideal %in% actual[OK]))
        tags[[1]] = tags[[1]][todo]
        tags[[2]] = tags[[2]][todo]
        }
    jnk = foreach(i = 1:length(tags[[1]])) %dopar% {
        fn = paste0(tags[[1]][i], ".meme")
        odir = paste0(tags[[1]][i], "_out")
        MotifDb::export(mdb[tags[[2]][i]], con = fn, "meme")
        comm = paste0("fimo ", parmvec, " --oc ", odir, " ", 
            fn, " ", fasta)
        comm = paste0(commprefix, " ", comm)
        if (setuponly) {
            writeLines(comm, paste0(tags[[1]][i], ".comm"))
        }
        else {
            print(comm)
            system(comm)
        }
    }
}

