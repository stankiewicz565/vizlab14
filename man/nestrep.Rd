\name{nestrep}
\alias{nestrep}
\docType{data}
\title{
data on nested repetitive elements, from UCSC/AnnotationHub
}
\description{
data on nested repetitive elements, from UCSC/AnnotationHub
}
\usage{data(nestrep)}
\format{
  The format is: \cr
Formal class 'GRanges' [package "GenomicRanges"] with 6 slots \cr
  ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 4 slots \cr
  .. .. ..@ values         : Factor w/ 93 levels "chr1","chr2",..: 1 2 3 4 5 6 7 9 10 8 ... \cr
  .. .. ..@ lengths        : int [1:90] 47078 44332 36905 35056 33488 30949 30252 27161 24187 32427 ... \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ metadata       : list() \cr
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots \cr
  .. .. ..@ start          : int [1:575600] 23804 26583 30694 31293 38059 38256 43243 54465 61696 71369 ... \cr
  .. .. ..@ width          : int [1:575600] 645 555 438 462 2568 2039 2511 1144 534 6867 ... \cr
  .. .. ..@ NAMES          : NULL \cr
  .. .. ..@ elementType    : chr "integer" \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ metadata       : list() \cr
  ..@ strand         :Formal class 'Rle' [package "IRanges"] with 4 slots \cr
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 1 2 1 2 1 3 1 2 1 2 ... \cr
  .. .. ..@ lengths        : int [1:278779] 1 1 6 1 1 2 6 1 1 1 ... \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ metadata       : list() \cr
  ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots \cr
  .. .. ..@ rownames       : NULL \cr
  .. .. ..@ nrows          : int 575600 \cr
  .. .. ..@ listData       :List of 5 \cr
  .. .. .. ..$ name   : chr [1:575600] "L2b" "L2c" "MLT1A" "MIRc" ... \cr
  .. .. .. ..$ score  : num [1:575600] 49 0 508 218 361 157 476 6 166 228 ... \cr
  .. .. .. ..$ itemRgb: chr [1:575600] NA NA NA NA ... \cr
  .. .. .. ..$ thick  :Formal class 'IRanges' [package "IRanges"] with 6 slots \cr
  .. .. .. .. .. ..@ start          : int [1:575600] 23804 26583 30694 31293 38059 38256 43243 54465 61696 71369 ... \cr
  .. .. .. .. .. ..@ width          : int [1:575600] 645 555 438 462 2568 2039 2511 1144 534 6867 ... \cr
  .. .. .. .. .. ..@ NAMES          : NULL \cr
  .. .. .. .. .. ..@ elementType    : chr "integer" \cr
  .. .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. .. ..@ metadata       : list() \cr
  .. .. .. ..$ blocks :Formal class 'CompressedIRangesList' [package "IRanges"] with 5 slots \cr
  .. .. .. .. .. ..@ elementType    : chr "IRanges" \cr
  .. .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. .. ..@ metadata       : list() \cr
  .. .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots \cr
  .. .. .. .. .. .. .. ..@ end            : int [1:575600] 2 4 6 8 10 13 15 17 19 25 ... \cr
  .. .. .. .. .. .. .. ..@ NAMES          : NULL \cr
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer" \cr
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. .. .. .. ..@ metadata       : list() \cr
  .. .. .. .. .. ..@ unlistData     :Formal class 'IRanges' [package "IRanges"] with 6 slots \cr
  .. .. .. .. .. .. .. ..@ start          : int [1:1483482] 1 452 1 472 1 260 1 442 1 2275 ... \cr
  .. .. .. .. .. .. .. ..@ width          : int [1:1483482] 235 194 208 84 155 179 143 21 133 294 ... \cr
  .. .. .. .. .. .. .. ..@ NAMES          : NULL \cr
  .. .. .. .. .. .. .. ..@ elementType    : chr "integer" \cr
  .. .. .. .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. .. .. .. ..@ metadata       : list() \cr
  .. .. ..@ elementType    : chr "ANY" \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ metadata       : list() \cr
  ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomicRanges"] with 4 slots \cr
  .. .. ..@ seqnames   : chr [1:93] "chr1" "chr2" "chr3" "chr4" ... \cr
  .. .. ..@ seqlengths : int [1:93] 249250621 243199373 198022430 191154276 180915260 171115067 159138663 155270560 146364022 141213431 ... \cr
  .. .. ..@ is_circular: logi [1:93] NA NA NA NA NA NA ... \cr
  .. .. ..@ genome     : chr [1:93] "hg19" "hg19" "hg19" "hg19" ... \cr
  ..@ metadata       :List of 2 \cr
  .. ..$ trackLine:Formal class 'BasicTrackLine' [package "rtracklayer"] with 13 slots \cr
  .. .. .. ..@ itemRgb      : logi(0)  \cr
  .. .. .. ..@ useScore     : logi(0)  \cr
  .. .. .. ..@ group        : chr(0)  \cr
  .. .. .. ..@ db           : chr(0)  \cr
  .. .. .. ..@ offset       : num(0)  \cr
  .. .. .. ..@ url          : chr "NA" \cr
  .. .. .. ..@ htmlUrl      : chr(0)  \cr
  .. .. .. ..@ colorByStrand: num[0 , 0 ]  \cr
  .. .. .. ..@ name         : chr "nestedRepeats" \cr
  .. .. .. ..@ description  : chr "nestedRepeats" \cr
  .. .. .. ..@ visibility   : chr "1" \cr
  .. .. .. ..@ color        : int(0)  \cr
  .. .. .. ..@ priority     : num(0)  \cr
  .. ..$          :List of 10 \cr
  .. .. ..$ BiocVersion  :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots \cr
  .. .. .. .. ..@ elementType    : chr "character" \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots \cr
  .. .. .. .. .. .. ..@ end            : int 4 \cr
  .. .. .. .. .. .. ..@ NAMES          : chr "5928" \cr
  .. .. .. .. .. .. ..@ elementType    : chr "integer" \cr
  .. .. .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. .. .. ..@ metadata       : list() \cr
  .. .. .. .. ..@ unlistData     : chr [1:4] "2.12" "2.13" "2.14" "3.0" \cr
  .. .. ..$ DataProvider : chr "hgdownload.cse.ucsc.edu" \cr
  .. .. ..$ Description  : chr "GRanges object from UCSC track 'Interrupted Rpts'" \cr
  .. .. ..$ Genome       : chr "hg19" \cr
  .. .. ..$ Tags         :Formal class 'CompressedCharacterList' [package "IRanges"] with 5 slots \cr
  .. .. .. .. ..@ elementType    : chr "character" \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots \cr
  .. .. .. .. .. .. ..@ end            : int 6 \cr
  .. .. .. .. .. .. ..@ NAMES          : chr "5928" \cr
  .. .. .. .. .. .. ..@ elementType    : chr "integer" \cr
  .. .. .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. .. .. ..@ metadata       : list() \cr
  .. .. .. .. ..@ unlistData     : chr [1:6] "nestedRepeats" "UCSC" "track" "Gene" ... \cr
  .. .. ..$ SourceUrl    : chr "rtracklayer://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/nestedRepeats" \cr
  .. .. ..$ SourceVersion: chr "hg19" \cr
  .. .. ..$ Species      : chr "Homo sapiens" \cr
  .. .. ..$ RDataPath    : chr "http://annotationhub.bioconductor.org/ah/goldenpath/hg19/database/nestedRepeats_0.0.1.RData" \cr
  .. .. ..$ RDataName    : chr "goldenpath.hg19.database.nestedRepeats_0.0.1.RData"
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
see the RDataName element
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(nestrep)
nestrep
}
\keyword{datasets}
