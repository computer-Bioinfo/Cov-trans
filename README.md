# Cov-trans
Discontinuous transcription allows coronaviruses and other RNA viruses to efficiently replicate and transmit within host
cells, enhancing their adaptability and survival. Assembling viral transcripts is crucial for virology research and the
development of antiviral strategies. However, traditional transcripts assembly methods primarily designed for variable
alternative splicing events in eukaryotes are not suitable for the viral transcripts assembly problem. Existing algorithms
dedicated to viral transcripts assembly suffer from low accuracy. There is an urgent need to develop a highly accurate
viral transcripts assembly algorithm. Here, we propose Cov-trans, a reference-based transcripts assembler specifically
tailored for the discontinuous transcription of coronaviruses. Cov-trans first identifies canonical transcripts based on
discontinuous transcription mechanisms, mRNA translation start and stop codons, as well as read alignment information.
Subsequently, it formulates the assembly of non-canonical transcripts as path extraction on a segment graph, employing
mixed integer linear programming to recover these non-canonical transcripts. We benchmarked Cov-trans against the viral
transcripts assembler Jumper, as well as two state-of-the-art eukaryotic transcriptome assemblers, Stringtie and Scallop.
Experimental results show that Cov-trans outperforms others in both accuracy and recall, with a notable strength in
accurately identifying the boundaries of canonical transcripts. 

# Pre-requistites for Cov-trans to run
Python 3.7	

Gurobi 9.5.1

pysam

numpy

pandas

# Usage
--bam(or -b) BAM  
--fasta(or -f)
                      [--min-base-qual MIN_BASE_QUAL]
                      [--min-mapping-qual MIN_MAPPING_QUAL] [-w WIDTH]
                      [--samplingFrequency SAMPLINGFREQUENCY]
                      [--sj_threshold SJ_THRESHOLD] [-n NEDGES]
                      [--phasing_threshold PHASING_THRESHOLD]
                      [--greedy GREEDY] [--outputCSV OUTPUTCSV]
                      [--outputPhasing OUTPUTPHASING] [--inputCSV INPUTCSV]
                      [--inputPhasing INPUTPHASING]
                      [--inputBreakpoints INPUTBREAKPOINTS]
                      [--inputEdges INPUTEDGES] [--outputGraph OUTPUTGRAPH]
                      [--outputDOT OUTPUTDOT]
                      [--outputTranscripts OUTPUTTRANSCRIPTS]
                      [--outputBreakpoints OUTPUTBREAKPOINTS]
                      [--outputEdges OUTPUTEDGES]
                      [--outputDecomposition OUTPUTDECOMPOSITION]
                      [--outputMatching OUTPUTMATCHING]
                      [--outputGTF OUTPUTGTF] [--report REPORT] [--noverbose]
                      [--threads THREADS] [--timelimit TIMELIMIT]
                      [--maxIter MAXITER]




