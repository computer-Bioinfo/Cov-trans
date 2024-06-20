# Cov-trans
Discontinuous transcription allows coronaviruses and other RNA viruses to efficiently replicate and transmit within host cells, enhancing their adaptability and survival. Assembling viral transcripts is crucial for virology research and the development of antiviral strategies. However, traditional transcripts assembly methods primarily designed for variable alternative splicing events in eukaryotes are not suitable for the viral transcripts assembly problem. 
The current algorithms designed for assembling viral transcripts often struggle with low accuracy in determining the transcript boundaries. There is an urgent need to develop a highly accurate viral transcripts assembly algorithm. Here, we propose Cov-trans, a reference-based transcripts assembler specifically tailored for the discontinuous transcription of coronaviruses. Cov-trans first identifies canonical transcripts based on discontinuous transcription mechanisms, start and stop codons, as well as reads alignment information. Subsequently, it formulates the assembly of non-canonical transcripts as path extraction on a segment graph, employing mixed integer linear programming to recover these non-canonical transcripts. Experimental results show that Cov-trans outperforms other assemblers in both accuracy and recall, with a notable strength in accurately identifying the boundaries of canonical transcripts.
# Pre-requistites for Cov-trans to run
Python 3.7	

Gurobi 9.5.1

samtools

pysam

numpy

pandas

# Usage
--bam(or -b) .bam  

--fasta(or -f) .fasta

--Graph OUTPUTTRANSCRIPTS
                     
--calTransGTF .gtf

--noncalTransGTF .gtf

--threads THREADS 

-h help

# Example

~~~
python cov-trans_main.py -f ../data/sars-cov-2.fasta -b ../data/sars2-sample.bam --calTransGTF sars-cov2_cal.gtf --noncalTransGTF sars-cov2_nonncal.gtf --Graph sars2_segmentgraph.out
 ~~~

# Feedback and bug reports
Your comments, bug reports, and suggestions are very welcomed. They will help us to further improve Cov-trans. If you have any troubles running Cov-trans, please contact us (Email:2021020692@qdu.edu.cn).




