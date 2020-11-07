\name{emtab1147sdrf}
\alias{emtab1147sdrf}
\docType{data}
\title{
sample-level data for the RNAseqData package for HNRNP/U2AF65 competition paper
}
\description{
sample-level data for the RNAseqData package for HNRNP/U2AF65 competition paper
}
\usage{data(emtab1147sdrf)}
\format{
  The format is: \cr
Formal class 'DataFrame' [package "IRanges"] with 6 slots \cr
  ..@ rownames       : chr [1:8] "ERR127302" "ERR127303" "ERR127304" "ERR127305" ... \cr
  ..@ nrows          : int 8 \cr
  ..@ listData       :List of 33 \cr
  .. ..$ Source.Name                     : chr [1:8] "HNRNPC knockdown 1 (KD1) - Replicate 1" "HNRNPC knockdown 1 (KD1) - Replicate 2" "HNRNPC knockdown 2 (KD2) - Replicate 1" "HNRNPC knockdown 2 (KD2) - Replicate 2" ... \cr
  .. ..$ Comment.ENA_SAMPLE.             : chr [1:8] "ERS140919" "ERS140921" "ERS140924" "ERS140922" ... \cr
  .. ..$ Material.Type                   : chr [1:8] "cell" "cell" "cell" "cell" ... \cr
  .. ..$ Characteristics.Organism.       : chr [1:8] "Homo sapiens" "Homo sapiens" "Homo sapiens" "Homo sapiens" ... \cr
  .. ..$ Characteristics.Genotype.       : chr [1:8] "HNRNPC knockdown with siRNA1" "HNRNPC knockdown with siRNA1" "HNRNPC knockdown with siRNA2" "HNRNPC knockdown with siRNA2" ... \cr
  .. ..$ Characteristics.CellLine.       : chr [1:8] "HeLa" "HeLa" "HeLa" "HeLa" ... \cr
  .. ..$ Characteristics.Sex.            : chr [1:8] "female" "female" "female" "female" ... \cr
  .. ..$ Protocol.REF                    : chr [1:8] "P-MTAB-26782" "P-MTAB-26782" "P-MTAB-26782" "P-MTAB-26782" ... \cr
  .. ..$ Protocol.REF.1                  : chr [1:8] "P-MTAB-26783" "P-MTAB-26783" "P-MTAB-26783" "P-MTAB-26783" ... \cr
  .. ..$ Protocol.REF.2                  : chr [1:8] "P-MTAB-26784" "P-MTAB-26784" "P-MTAB-26784" "P-MTAB-26784" ... \cr
  .. ..$ Extract.Name                    : chr [1:8] "HNRNPC knockdown 1 (KD1) - Replicate 1" "HNRNPC knockdown 1 (KD1) - Replicate 2" "HNRNPC knockdown 2 (KD2) - Replicate 1" "HNRNPC knockdown 2 (KD2) - Replicate 2" ... \cr
  .. ..$ Comment.ORIENTATION.            : chr [1:8] "5-3-3-5" "5-3-3-5" "5-3-3-5" "5-3-3-5" ... \cr
  .. ..$ Comment.NOMINAL_LENGTH.         : int [1:8] 175 175 175 175 175 175 175 175 \cr
  .. ..$ Comment.NOMINAL_SDEV.           : int [1:8] 25 25 25 25 25 25 25 25 \cr
  .. ..$ Comment.LIBRARY_LAYOUT.         : chr [1:8] "PAIRED" "PAIRED" "PAIRED" "PAIRED" ... \cr
  .. ..$ Comment.LIBRARY_SOURCE.         : chr [1:8] "TRANSCRIPTOMIC" "TRANSCRIPTOMIC" "TRANSCRIPTOMIC" "TRANSCRIPTOMIC" ... \cr
  .. ..$ Comment.LIBRARY_STRATEGY.       : chr [1:8] "RNA-Seq" "RNA-Seq" "RNA-Seq" "RNA-Seq" ... \cr
  .. ..$ Comment.LIBRARY_SELECTION.      : chr [1:8] "RANDOM" "RANDOM" "RANDOM" "RANDOM" ... \cr
  .. ..$ Protocol.REF.3                  : chr [1:8] "P-MTAB-26785" "P-MTAB-26785" "P-MTAB-26785" "P-MTAB-26785" ... \cr
  .. ..$ Performer                       : chr [1:8] "CRUK Cambridge Research Institute" "CRUK Cambridge Research Institute" "CRUK Cambridge Research Institute" "CRUK Cambridge Research Institute" ... \cr
  .. ..$ Assay.Name                      : chr [1:8] "SLX-2646.CRIRUN_441.s_1" "SLX-2647.CRIRUN_441.s_2" "SLX-2648.CRIRUN_441.s_3" "SLX-2649.CRIRUN_441.s_4" ... \cr
  .. ..$ Technology.Type                 : chr [1:8] "high_throughput_sequencing" "high_throughput_sequencing" "high_throughput_sequencing" "high_throughput_sequencing" ... \cr
  .. ..$ Comment.ENA_EXPERIMENT.         : chr [1:8] "ERX103443" "ERX103446" "ERX103450" "ERX103447" ... \cr
  .. ..$ Scan.Name                       : chr [1:8] "SLX-2646.CRIRUN_441.s_1.1_sequence.txt.gz" "SLX-2647.CRIRUN_441.s_2.1_sequence.txt.gz" "SLX-2648.CRIRUN_441.s_3.1_sequence.txt.gz" "SLX-2649.CRIRUN_441.s_4.1_sequence.txt.gz" ... \cr
  .. ..$ Comment.SUBMITTED_FILE_NAME.    : chr [1:8] "SLX-2646.CRIRUN_441.s_1.1_sequence.txt.gz" "SLX-2647.CRIRUN_441.s_2.1_sequence.txt.gz" "SLX-2648.CRIRUN_441.s_3.1_sequence.txt.gz" "SLX-2649.CRIRUN_441.s_4.1_sequence.txt.gz" ... \cr
  .. ..$ Comment.ENA_RUN.                : chr [1:8] "ERR127302" "ERR127303" "ERR127304" "ERR127305" ... \cr
  .. ..$ Comment.FASTQ_URI.              : chr [1:8] "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR127/ERR127302/ERR127302_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR127/ERR127303/ERR127303_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR127/ERR127304/ERR127304_1.fastq.gz" "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR127/ERR127305/ERR127305_1.fastq.gz" ... \cr
  .. ..$ Comment.SPOT_LENGTH.            : int [1:8] 72 72 72 72 72 72 72 72 \cr
  .. ..$ Comment.quality_scoring_system. : chr [1:8] "phred" "phred" "phred" "phred" ... \cr
  .. ..$ Comment.quality_encoding.       : chr [1:8] "ascii" "ascii" "ascii" "ascii" ... \cr
  .. ..$ Comment.ascii_offset.           : chr [1:8] "@" "@" "@" "@" ... \cr
  .. ..$ Comment.READ_INDEX_1_BASE_CCORD.: int [1:8] 73 73 73 73 73 73 73 73 \cr
  .. ..$ Factor.Value.RNAI.              : chr [1:8] "KD1" "KD1" "KD2" "KD2" ... \cr
  ..@ elementType    : chr "ANY" \cr
  ..@ elementMetadata: NULL \cr
  ..@ metadata       : list() 
}
\details{
taken verbatim
}
\source{
EBI ArrayExpress
}
\references{
PMID 23374342
}
\examples{
data(emtab1147sdrf)
emtab1147sdrf[1:4,1:4]
}
\keyword{datasets}
