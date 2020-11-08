

ggmDefC = function() get(load("~/.ggmcache/ggmcache.rda"))  # need reference object?

browseTracks = function( tlist, ylocs=1:length(tlist), padpct=20,
   rectWidth=.2 ) {}

 fn = dir(system.file(
   "extdata", package="RNAseqData.HNRNPC.bam.chr14"), full=TRUE,
   patt="bam$")

 bfv = GenomicFiles(,fn)

browseWreads = function(sym="ADSS1", doAlu=TRUE, alugr=alu14, aluy=2, 
     bamf = files(bfv)[1], ready=3,
     cachedir="~/.ggmcache", y=1, y2=1.2, ylim=c(0,6), padpct=20) {
 base = renderGene(sym, cachedir=cachedir, y=y, y2=y2, 
           ylim=ylim, padpct=padpct) # populates cache with relevant info
 rng = rangeOfGene(sym, ggmDefC())
 wid = end(rng) - start(rng)
 pad = padpct*wid/100
 left = start(rng) - pad
 if (doAlu) 
   return( base %>% 
      rectsAtY(sym=sym, gr=alugr, y=aluy) %>% addText(sym="Alu", x=left, y=aluy) %>%
      readsAtY(sym=sym, bamf=bamf, y=ready) %>% addText(sym="reads", x=left, y=ready) 
      )
 base %>% readsAtY(sym=sym, bamf=bamf, y=ready) %>% addText(sym="reads", x=left, y=ready) 
}

renderGene = function(sym="ADSS1", 
     cachedir="~/.ggmcache", y=1, y2=1.2, ylim=c(0,6), padpct=20)
  {
#
# key interface for browser: represent exons as rectangles on a
# horizontal path, establish communications for other elements
#
  ggmcache = setupCache( cachedir )
  gvm = ggvForGeneModel( sym, ggmcache, padpct=padpct ) 
  rng = rangeOfGene(sym, ggmcache)
  wid = end(rng) - start(rng)
  pad = padpct*wid/100
  gvm %>% exonRects(y=y,y2=y2) %>% 
   addText(sym=sym, x=start(rng)-pad, y=y ) %>%
   addText(sym=" ", x=end(rng)+pad, y=y ) %>%
   axesForGgv(chr=chrOfGene(sym, ggmcache), ylim=ylim)
  }


#   aluRects(y=y+1,y2=y2+1) %>% 
#   addText(sym="Alu", x=.low, y=y+1 ) %>%

addText = function(..., sym, x, y) {
   tdf = data.frame(z=x, w=y, sym=sym)
   layer_text(..., ~z, ~w, fontSize:= 18, text:= ~sym, data=tdf)
}

ggvForGeneModel = function( sym, ggmcache, padpct,
  aluRng=alu14, rsRng=b1g, cachepath="~/.ggmcache/ggmcache.rda" ) {
#
# initialize ggvis stream with exons from gene model
#
 require( GenomicRanges )
 m1 = getGeneModel( sym, ggmcache )  # getting model is costly
 chr = as.character(seqnames(m1)[1])
 rm1 = range(m1)
 assign(paste0("ggmRange_", sym), rm1, ggmcache)  # can use range for trimming
 save(ggmcache, file=cachepath)  
 assign(paste0("ggmChr_", sym), chr, ggmcache)  # communicate chr
 save(ggmcache, file=cachepath)  # update cache
 as(m1, "data.frame") %>% ggvis(~start)    # x is start coord of exons
}

exonRects = function( ..., y=1, y2=2) {
 # assumes part of pipe following ggvForGeneModel
 layer_rects(..., y=y, x2=~end, y2=y2) 
}

rangeOfGene = function(sym, ggmcache=ggmDefC()) 
  get(paste0("ggmRange_", sym), ggmcache)
chrOfGene = function(sym, ggmcache=ggmDefC()) 
  get(paste0("ggmChr_", sym), ggmcache)

gr2confinedDF = function(input, sym, ignore.strand=TRUE,
    ggmcache=ggmDefC(), tx=force) {
 target = get(paste0("ggmRange_", sym), ggmcache)
 fo = findOverlaps(input, target, ignore.strand=ignore.strand)
 if (length(queryHits(fo))<1) return(NULL)
 tx(as(input[queryHits(fo)], "data.frame"))
 }

rectsAtY = function(..., sym, gr, ignore.strand=TRUE, y=2, y2=2.2, tx=force, ggmcache=ggmDefC()) {
 layer_rects(..., x=~start, y=y, x2=~end, y2=y2, 
          data=gr2confinedDF(gr, sym, ignore.strand=ignore.strand, tx=tx,
                 ggmcache=ggmcache) ) }

rl2df = function(rl) {
 ans = data.frame(pos=cumsum(runLength(rl)), y=runValue(rl), yz=0)
 ans[-nrow(ans),]
}

readsAtY = function(..., sym, bamf, y=3, 
   txy=function(x) pmin(x,7)/7+y, ggmcache=ggmDefC()) {
  bp = ScanBamParam(which=rangeOfGene(sym, ggmcache))
  chr = chrOfGene(sym, ggmcache)
  fcov = rl2df(coverage( bamf, param=bp )[[chr]])
  fcov[,2] = txy(fcov[,2])
  rng = rangeOfGene(sym, ggmcache)
  left = start(rng)
  layer_points(..., x=~pos, y=~y, size=.5, data=fcov ) 
}

axesForGgv = function(..., chr, nxticks=4, nyticks=0, xgrid=FALSE, ygrid=FALSE, 
   ylim=c(0,6))
  add_guide_axis(vis=..., type="x", scale="x", orient="top", title=chr,
     grid=xgrid, ticks=nxticks) %>% 
  set_dscale("y", type="numeric", domain=ylim) %>%
  add_guide_axis(vis=..., type="y", scale="y", grid=ygrid, orient="left", ticks=nyticks,
     title="") 

renderGeneOLD = function(sym="ADSS1", 
     cachedir="~/.ggmcache", y=1, y2=1.2, ylim=c(0,6), padpct=20)
  {
#
# key interface for browser: represent exons as rectangles on a
# horizontal path, establish communications for other elements
#
  ggmcache = setupCache( cachedir )
  gvm = ggvForGeneModel( sym, ggmcache, padpct=padpct ) 
  rng = rangeOfGene(sym, ggmcache)
  wid = end(rng) - start(rng)
  pad = padpct*wid/100
  gvm %>% exonRects(y=y,y2=y2) %>% 
   addText(sym=sym, x=start(rng)-pad, y=y ) %>%
   addText(sym=" ", x=end(rng)+pad, y=y ) %>%
  set_dscale("y", type="numeric", domain=ylim) %>%
  add_guide_axis(vis=..., type="y", scale="y", grid=ygrid, orient="left", ticks=nyticks,
     title="") 
#   axesForGgv(chr=chrOfGene(sym, ggmcache), ylim=ylim)
  }

renderGene = function(sym="ADSS1",
     cachedir="~/.ggmcache", y=1, y2=1.2, ylim=c(0,6), padpct=20)
  {
#
# key interface for browser: represent exons as rectangles on a
# horizontal path, establish communications for other elements
#
  ggmcache = setupCache( cachedir )
  gvm = ggvForGeneModel( sym, ggmcache, padpct=padpct )
  rng = rangeOfGene(sym, ggmcache)
  wid = end(rng) - start(rng)
  pad = padpct*wid/100
  gvm %>% exonRects(y=y,y2=y2) %>%
   addText(sym=sym, x=start(rng)-pad, y=y ) %>%
   addText(sym=" ", x=end(rng)+pad, y=y ) %>%
   scale_numeric("y", range=ylim) %>% # set_dscale("y", type="numeric", domain=ylim) %>%
  add_axis(type="y", scale="y", grid=FALSE, orient="left", ticks=0,
     title="") %>%
  add_axis(type="x", scale="x", orient="top", title=as.character(seqnames(rng)[1]),
     grid=FALSE, ticks=4)
  }

