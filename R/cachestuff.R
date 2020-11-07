
setupCache = function(cachedir) {
 cachepath = paste0(cachedir, "/ggmcache.rda")

 if (!(file.exists(cachedir))) {
   message("creating cache folder...")
   dir.create(cachedir)
   }
 if (!file.exists(cachepath)) {
   ggmcache = new.env()
   }
 else {
   ggmcache = get(load(cachepath))
   }
 incache = names(ggmcache)
 if (length(incache)>500) warning("ggmcache has > 500 entries, consider removing")
 ggmcache
}

getGeneModel = function( sym, ggmcache, cachepath="~/.ggmcache/ggmcache.rda" ) {
 knownModels = ls(ggmcache)
 if (length(knownModels)>0 && sym %in% knownModels) 
   m1 = get(sym, ggmcache)
 else {
   message("computing model...")
   m1 = genemodel(sym)
   assign(sym, m1, ggmcache)
   message("caching...")
   save(ggmcache, file=cachepath)
   message("done.")
 }
 m1
}

