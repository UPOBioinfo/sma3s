# Sma3s (Sequence massive annotator using 3 modules) - Annotate complete proteomes &amp; transcriptomes

Sma3s has low computing requirements and can be used on virtually any computer. It is written in Perl language and you need its interpreter (http://www.perl.com), which is preinstalled in Linux and Mac OS X (in Windows it will not be necessary). Additionally, you need to install the Blast+ package for your operating system.

To annotate your sequence dataset, you only need the following files:
- Your query sequences in multi-FASTA format,
- The reference database, which you can download from our server: http://www.bioinfocabd.upo.es/sma3s/db/

## Linux
Usual command line for annotating proteomes:
*./sma3s_v2.pl -i query_dataset.fasta -d uniref90.fasta -goslim*

Usual command line for annotating transcriptomes:
 *./sma3s_v2.pl -i query_dataset.fasta -d uniref90.fasta -nucl -goslim*

Run "sma3s_v2.pl --help" for help with these and other advanced parameters.

## Windows
Install Blast+ for Windows from: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/
Save all the necessary files into a folder called 'annotation', and use the Windows binary file: Sma3s_v2.exe.
Execute cmd.exe, and write 'cd \annotation' in the console. 

## Mac OS X
Install Blast+ for Mac from: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/
Open Terminal from Applications/Utilities.


## Customized databases
Alternatively, you can use Sma3s with the whole UniProt database, if you are interested in a more sensitive, though more slowly, annotation. To do that, you must download a .dat file from UniProt (ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/), and install the Blast Legacy package (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/).

## Sma3s Website (including a video-tutorial)
http://www.bioinfocabd.upo.es/sma3s/
