filterGWASbyMap = function(gwr, map, DTattr="Disease.Trait") {
  if (!all(names(map) %in% c("original", "nodename"))) stop("map must have column names 'original' and 'nodename'")
  gwr = as(gwr, "GRanges")
  DT = values(gwr)[[DTattr]]
  kp = which(DT %in% map[,"original"])
  if (length(kp)==0) stop("no trait values in map$original")
  gwr = gwr[kp]
  dclass = map[,"nodename"]
  names(dclass) = map[,"original"]
  dclass = dclass[ as.character(values(gwr)[[DTattr]]) ]
  values(gwr)$dclass = dclass
  gwr
}
