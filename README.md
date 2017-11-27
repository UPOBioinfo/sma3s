# sma3s
Sma3s (Sequence massive annotator using 3 modules) - annotate complete proteomes &amp; transcriptomes

Sma3s has low computing requirements and can be used on virtually any computer. It is written in Perl language and you need its interpreter (http://www.perl.com), which is preinstalled in Linux and Mac OS X (in Windows it will not be necessary). Additionally, you need to install the Blast+ package for your operating system.

Finally, you will our Sma3s program:

>>> Sma3s_v2.pl <<< (see Creative Commons license)

To annotate your sequence dataset, you only need the following files:

Your query sequences in multi-FASTA format,
The reference database, which you can download from our server: http://www.bioinfocabd.upo.es/sma3s/db/
Usual command line for annotating proteomes:
./sma3s.pl -i query_dataset.fasta -d uniref90.fasta -goslim

Usual command line for annotating transcriptomes:
 ./sma3s.pl -i query_dataset.fasta -d uniref90.fasta -nucl -goslim

Run "sma3s_v2.pl --help" for help with these and other advanced parameters.

Alternatively, you can use Sma3s with the whole UniProt database, if you are interested in a more sensitive, though more slowly, annotation. To do that, you must download a .dat file from UniProt (ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/), and install the Blast Legacy package (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/).

