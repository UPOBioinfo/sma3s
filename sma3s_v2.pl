#!/usr/bin/perl -w
use strict;

# AJPerez, 2013-09-01
# Sma3s: a three-step modular sequence annotator for large sequence datasets
# v2, 2015-06-12: Blast+ & remove Bioperl, and other external libraries (CSanchez); evidence code use;
# Interpro/Interaction removed; included many performance and accuracy improvements, as well as summary results

# For HyGe library (used for biological enrichment)
my %LChooseCache;
my $HUGE_VAL = 1e100;
my $EPS = log( 1 / $HUGE_VAL );

# Annotation Types and Keyword categories, and GO Slim
my @TYPES = qw(GENENAME DESCRIPTION ENZYME GO KEYWORD PATHWAY);
my @KW_Biological_process = ("Acetoin biosynthesis", "Acetoin catabolism", "Acute phase", "Alginate biosynthesis", "Alkaloid metabolism", "Alkylphosphonate uptake", "Amino-acid biosynthesis", "Angiogenesis", "Antibiotic biosynthesis", "Antibiotic resistance", "Antiviral defense", "Apoptosis", "Arginine metabolism", "Aromatic hydrocarbons catabolism", "Arsenical resistance", "Ascorbate biosynthesis", "ATP synthesis", "Autoinducer synthesis", "Autophagy", "Auxin biosynthesis", "B-cell activation", "Bacteriocin immunity", "Behavior", "Biological rhythms", "Biomineralization", "Biotin biosynthesis", "Branched-chain amino acid catabolism", "Cadmium resistance", "Calvin cycle", "cAMP biosynthesis", "Carbohydrate metabolism","Carbon dioxide fixation", "Carnitine biosynthesis", "Carotenoid biosynthesis", "Catecholamine biosynthesis", "Catecholamine metabolism", "Cell adhesion", "Cell cycle", "Cell shape", "Cellulose biosynthesis", "cGMP biosynthesis", "Chemotaxis", "Chlorophyll biosynthesis", "Chromate resistance", "Chromosome partition", "Citrate utilization", "Cobalamin biosynthesis", "Coenzyme A biosynthesis", "Coenzyme M biosynthesis", "Collagen degradation", "Competence", "Conjugation", "Cytadherence", "Cytochrome c-type biogenesis", "Cytokinin biosynthesis", "Cytolysis", "Cytosine metabolism", "Deoxyribonucleotide synthesis", "Detoxification", "Diaminopimelate biosynthesis", "Differentiation", "Digestion", "DNA condensation", "DNA damage", "DNA excision", "DNA integration", "DNA recombination", "DNA replication", "DNA synthesis", "Endocytosis", "Enterobactin biosynthesis", "Erythrocyte maturation", "Ethylene biosynthesis", "Exocytosis", "Exopolysaccharide synthesis", "Fertilization", "Flagellar rotation", "Flavonoid biosynthesis", "Flight", "Flowering", "Folate biosynthesis", "Fruit ripening", "Galactitol metabolism", "Gaseous exchange", "Gastrulation", "Germination", "Gluconate utilization", "Gluconeogenesis", "Glutathione biosynthesis", "Glycerol metabolism", "Glycogen biosynthesis", "Glycolate pathway", "Glycolysis", "Glyoxylate bypass", "GPI-anchor biosynthesis", "Growth regulation", "Stress response", "Heme biosynthesis", "Hemolymph clotting", "Hemostasis", "Herbicide resistance", "Histidine metabolism", "Hydrogen peroxide", "Hypusine biosynthesis", "Immunity", "Inflammatory response", "Inositol biosynthesis", "Intron homing", "Iron storage", "Isoprene biosynthesis", "Karyogamy", "Keratinization", "Lactation", "Lactose biosynthesis", "Lactose metabolism", "Leukotriene biosynthesis", "Lignin biosynthesis", "Lignin degradation", "Lipid metabolism", "Lipopolysaccharide biosynthesis", "Luminescence", "Maltose metabolism", "Mandelate pathway", "Mast cell degranulation", "Meiosis", "Melanin biosynthesis", "Melatonin biosynthesis", "Menaquinone biosynthesis", "Mercuric resistance", "Methanogenesis", "Methanol utilization", "Methotrexate resistance", "Mineral balance", "Molybdenum cofactor biosynthesis", "mRNA processing", "Myogenesis", "Neurogenesis", "Neurotransmitter biosynthesis", "Neurotransmitter degradation", "Nitrate assimilation", "Nitrogen fixation", "Nodulation", "Nucleotide biosynthesis", "Nucleotide metabolism", "Nylon degradation", "One-carbon metabolism", "Pantothenate biosynthesis", "Pentose shunt", "Peptidoglycan synthesis", "PHA biosynthesis", "Phagocytosis", "PHB biosynthesis", "Phenylalanine catabolism", "Phenylpropanoid metabolism", "Pheromone response", "Phosphotransferase system", "Photorespiration", "Photosynthesis", "Phytochrome signaling pathway", "Plant defense", "Plasmid copy control", "Plasmid partition", "Plasminogen activation", "Polyamine biosynthesis", "Porphyrin biosynthesis", "Pregnancy", "Proline metabolism", "Protein biosynthesis", "Purine biosynthesis", "Purine metabolism", "Purine salvage", "Putrescine biosynthesis", "Pyridine nucleotide biosynthesis", "Pyridoxine biosynthesis", "Pyrimidine biosynthesis", "Queuosine biosynthesis", "Quinate metabolism", "Quorum sensing", "Restriction system", "Rhamnose metabolism", "Riboflavin biosynthesis", "Ribosome biogenesis", "RNA repair", "Viral RNA replication", "rRNA processing", "Self-incompatibility", "Sensory transduction", "Serotonin biosynthesis", "Spermidine biosynthesis", "Sporulation", "Starch biosynthesis", "Steroidogenesis", "Sulfate respiration", "Teichoic acid biosynthesis", "Tellurium resistance", "Terminal addition", "Tetrahydrobiopterin biosynthesis", "Thiamine biosynthesis", "Thiamine catabolism", "Tissue remodeling", "Transcription", "Translation regulation", "Transport", "Transposition", "Tricarboxylic acid cycle", "tRNA processing", "Tryptophan catabolism", "Tyrosine catabolism", "Ubiquinone biosynthesis", "Ubl conjugation pathway", "Unfolded protein response", "Urea cycle", "Virulence", "Nonsense-mediated mRNA decay", "Taxol biosynthesis", "Wnt signaling pathway", "Chlorophyll catabolism", "PQQ biosynthesis", "Chondrogenesis", "Osteogenesis", "Thyroid hormones biosynthesis", "Two-component regulatory system", "Hibernation", "Notch signaling pathway", "Interferon antiviral system evasion", "Auxin signaling pathway", "Hypersensitive response elicitation", "Cytokinin signaling pathway", "Ethylene signaling pathway", "Abscisic acid biosynthesis", "Abscisic acid signaling pathway", "Gibberellin signaling pathway", "RNA-mediated gene silencing", "Host-virus interaction", "Cell wall biogenesis/degradation", "Peroxisome biogenesis", "Cilium biogenesis/degradation", "Capsule biogenesis/degradation", "Insecticide resistance", "Nickel insertion", "Bacterial flagellum biogenesis", "Hearing", "Fimbrium biogenesis", "Brassinosteroid signaling pathway", "Cap snatching", "Virus entry into host cell", "Syncytium formation induced by viral infection", "Jasmonic acid signaling pathway", "Virus exit from host cell", "Viral DNA replication", "Viral transcription", "Archaeal flagellum biogenesis", "Necrosis", "Viral genome excision", "Viral latency");
my @KW_Cellular_component = ("Amyloid", "Antenna complex", "Apoplast", "Centromere", "CF(0)", "CF(1)", "Chlorosome", "Chromosome", "Chylomicron", "Nematocyst", "Copulatory plug", "Cuticle", "DNA-directed RNA polymerase", "Dynein", "Endoplasmic reticulum", "Exosome", "Fimbrium", "Glycosome", "Glyoxysome", "Golgi apparatus", "HDL", "Hydrogenosome", "Intermediate filament", "Keratin", "LDL", "Lysosome", "Membrane", "Membrane attack complex", "MHC I", "MHC II", "Microtubule", "Mitochondrion", "Nucleus", "Nucleomorph", "Lipid droplet", "Periplasm", "Peroxisome", "Photosystem I", "Photosystem II", "Phycobilisome", "Primosome", "Proteasome", "Reaction center", "Sarcoplasmic reticulum", "Signal recognition particle", "Signalosome", "Spliceosome", "Thick filament", "Thylakoid", "Viral occlusion body", "VLDL", "Vacuole", "Plastid", "Virion", "Cytoplasm", "Secreted", "Cell junction", "Cell projection", "Endosome", "Cytoplasmic vesicle", "Archaeal flagellum", "Bacterial flagellum", "Kinetochore", "Mitosome", "Host cell junction", "Host cell projection", "Host cytoplasm", "Host cytoplasmic vesicle", "Host endoplasmic reticulum", "Host endosome", "Host Golgi apparatus", "Host lipid droplet", "Host lysosome", "Host mitochondrion", "Host nucleus", "Host periplasm", "Host thylakoid", "Target cell cytoplasm");
my @KW_Developmental_stage = ("Early protein", "Fruiting body", "Heterocyst", "Late protein", "Merozoite", "Sporozoite", "Bradyzoite", "Tachyzoite", "Trophozoite");
my @KW_Disease = ("AIDS", "Albinism", "Allergen", "Alport syndrome", "Alzheimer disease", "Ectodermal dysplasia", "Tumor suppressor", "Atherosclerosis", "Autoimmune encephalomyelitis", "Autoimmune uveitis", "Bernard Soulier syndrome", "Cardiomyopathy", "Chronic granulomatous disease", "Cockayne syndrome", "Cone-rod dystrophy", "Crown gall tumor", "Cystinuria", "Deafness", "Dental caries", "Diabetes insipidus", "Diabetes mellitus", "Disease mutation", "Down syndrome", "Dwarfism", "Ehlers-Danlos syndrome", "Epidermolysis bullosa", "Gaucher disease", "Glutaricaciduria", "Glycogen storage disease", "Gangliosidosis", "Gout", "Hemophilia", "Hereditary hemolytic anemia", "Hereditary multiple exostoses", "Hereditary nonpolyposis colorectal cancer", "Hirschsprung disease", "Holoprosencephaly", "Hyperlipidemia", "Leber hereditary optic neuropathy", "Leigh syndrome", "Li-Fraumeni syndrome", "Lissencephaly", "Long QT syndrome", "Malaria", "Maple syrup urine disease", "Mucopolysaccharidosis", "Neurodegeneration", "Obesity", "Oncogene", "Phenylketonuria", "Neuropathy", "Proto-oncogene", "Pseudohermaphroditism", "Retinitis pigmentosa", "Rhizomelic chondrodysplasia punctata", "SCID", "Stargardt disease", "Stickler syndrome", "Systemic lupus erythematosus", "Thrombophilia", "Trypanosomiasis", "von Willebrand disease", "Whooping cough", "Williams-Beuren syndrome", "Xeroderma pigmentosum", "MELAS syndrome", "Epilepsy", "Cataract", "Congenital disorder of glycosylation", "Leber congenital amaurosis", "Primary microcephaly", "Parkinson disease", "Parkinsonism", "Bartter syndrome", "Congenital muscular dystrophy", "Age-related macular degeneration", "Fanconi anemia", "Progressive external ophthalmoplegia", "Short QT syndrome", "Mucolipidosis", "Limb-girdle muscular dystrophy", "Aicardi-Goutieres syndrome", "Familial hemophagocytic lymphohistiocytosis", "Congenital adrenal hyperplasia", "Glaucoma", "Kallmann syndrome", "Peroxisome biogenesis disorder", "Atrial septal defect", "Ichthyosis", "Primary hypomagnesemia", "Congenital hypothyroidism", "Congenital erythrocytosis", "Amelogenesis imperfecta", "Osteopetrosis", "Intrahepatic cholestasis", "Craniosynostosis", "Mental retardation", "Brugada syndrome", "Aortic aneurysm", "Congenital myasthenic syndrome", "Palmoplantar keratoderma", "Amyloidosis", "Dyskeratosis congenita", "Microphthalmia", "Congenital stationary night blindness", "Hypogonadotropic hypogonadism", "Atrial fibrillation", "Pontocerebellar hypoplasia", "Congenital generalized lipodystrophy", "Dystonia", "Diamond-Blackfan anemia", "Leukodystrophy", "Niemann-Pick disease", "Heterotaxy", "Nemaline myopathy", "Asthma", "Peters anomaly", "Myofibrillar myopathy", "Cushing syndrome", "Hypotrichosis", "Osteogenesis imperfecta", "Premature ovarian failure", "Emery-Dreifuss muscular dystrophy", "Hemolytic uremic syndrome", "Ciliopathy", "Schizophrenia", "Corneal dystrophy", "Dystroglycanopathy", "Autism spectrum disorder");
my %SLIM = ("GO:0000228"=>"C:nuclear chromosome", "GO:0000229"=>"C:cytoplasmic chromosome", "GO:0005576"=>"C:extracellular region", "GO:0005578"=>"C:proteinaceous extracellular matrix", "GO:0005615"=>"C:extracellular space", "GO:0005618"=>"C:cell wall", "GO:0005622"=>"C:intracellular", "GO:0005623"=>"C:cell", "GO:0005634"=>"C:nucleus", "GO:0005635"=>"C:nuclear envelope", "GO:0005654"=>"C:nucleoplasm", "GO:0005694"=>"C:chromosome", "GO:0005730"=>"C:nucleolus", "GO:0005737"=>"C:cytoplasm", "GO:0005739"=>"C:mitochondrion", "GO:0005764"=>"C:lysosome", "GO:0005768"=>"C:endosome", "GO:0005773"=>"C:vacuole", "GO:0005777"=>"C:peroxisome", "GO:0005783"=>"C:endoplasmic reticulum", "GO:0005794"=>"C:Golgi apparatus", "GO:0005811"=>"C:lipid particle", "GO:0005815"=>"C:microtubule organizing center", "GO:0005829"=>"C:cytosol", "GO:0005840"=>"C:ribosome", "GO:0005856"=>"C:cytoskeleton", "GO:0005886"=>"C:plasma membrane", "GO:0005929"=>"C:cilium", "GO:0009536"=>"C:plastid", "GO:0009579"=>"C:thylakoid", "GO:0016023"=>"C:cytoplasmic, membrane-bounded vesicle", "GO:0030312"=>"C:external encapsulating structure", "GO:0043226"=>"C:organelle", "GO:0043234"=>"C:protein complex", "GO:0000988"=>"F:transcription factor activity, protein binding", "GO:0001071"=>"F:nucleic acid binding transcription factor activity", "GO:0003677"=>"F:DNA binding", "GO:0003723"=>"F:RNA binding", "GO:0003729"=>"F:mRNA binding", "GO:0003735"=>"F:structural constituent of ribosome", "GO:0003924"=>"F:GTPase activity", "GO:0004386"=>"F:helicase activity", "GO:0004518"=>"F:nuclease activity", "GO:0004871"=>"F:signal transducer activity", "GO:0005198"=>"F:structural molecule activity", "GO:0008092"=>"F:cytoskeletal protein binding", "GO:0008134"=>"F:transcription factor binding", "GO:0008135"=>"F:translation factor activity, RNA binding", "GO:0008168"=>"F:methyltransferase activity", "GO:0008233"=>"F:peptidase activity", "GO:0008289"=>"F:lipid binding", "GO:0008565"=>"F:protein transporter activity", "GO:0016301"=>"F:kinase activity", "GO:0016491"=>"F:oxidoreductase activity", "GO:0016746"=>"F:transferase activity, transferring acyl groups", "GO:0016757"=>"F:transferase activity, transferring glycosyl groups", "GO:0016765"=>"F:transferase activity, transferring alkyl or aryl (other than methyl) groups", "GO:0016779"=>"F:nucleotidyltransferase activity", "GO:0016791"=>"F:phosphatase activity", "GO:0016798"=>"F:hydrolase activity, acting on glycosyl bonds", "GO:0016810"=>"F:hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds", "GO:0016829"=>"F:lyase activity", "GO:0016853"=>"F:isomerase activity", "GO:0016874"=>"F:ligase activity", "GO:0016887"=>"F:ATPase activity", "GO:0019843"=>"F:rRNA binding", "GO:0019899"=>"F:enzyme binding", "GO:0022857"=>"F:transmembrane transporter activity", "GO:0030234"=>"F:enzyme regulator activity", "GO:0030533"=>"F:triplet codon-amino acid adaptor activity", "GO:0030555"=>"F:RNA modification guide activity", "GO:0030674"=>"F:protein binding, bridging", "GO:0032182"=>"F:ubiquitin-like protein binding", "GO:0042393"=>"F:histone binding", "GO:0043167"=>"F:ion binding", "GO:0051082"=>"F:unfolded protein binding", "GO:0000003"=>"P:reproduction", "GO:0000902"=>"P:cell morphogenesis", "GO:0002376"=>"P:immune system process", "GO:0003013"=>"P:circulatory system process", "GO:0005975"=>"P:carbohydrate metabolic process", "GO:0006091"=>"P:generation of precursor metabolites and energy", "GO:0006259"=>"P:DNA metabolic process", "GO:0006397"=>"P:mRNA processing", "GO:0006399"=>"P:tRNA metabolic process", "GO:0006412"=>"P:translation", "GO:0006457"=>"P:protein folding", "GO:0006461"=>"P:protein complex assembly", "GO:0006464"=>"P:cellular protein modification process", "GO:0006520"=>"P:cellular amino acid metabolic process", "GO:0006605"=>"P:protein targeting", "GO:0006629"=>"P:lipid metabolic process", "GO:0006790"=>"P:sulfur compound metabolic process", "GO:0006810"=>"P:transport", "GO:0006913"=>"P:nucleocytoplasmic transport", "GO:0006914"=>"P:autophagy", "GO:0006950"=>"P:response to stress", "GO:0007005"=>"P:mitochondrion organization", "GO:0007009"=>"P:plasma membrane organization", "GO:0007010"=>"P:cytoskeleton organization", "GO:0007034"=>"P:vacuolar transport", "GO:0007049"=>"P:cell cycle", "GO:0007059"=>"P:chromosome segregation", "GO:0007067"=>"P:mitotic nuclear division", "GO:0007155"=>"P:cell adhesion", "GO:0007165"=>"P:signal transduction", "GO:0007267"=>"P:cell-cell signaling", "GO:0007568"=>"P:aging", "GO:0008219"=>"P:cell death", "GO:0008283"=>"P:cell proliferation", "GO:0009056"=>"P:catabolic process", "GO:0009058"=>"P:biosynthetic process", "GO:0009790"=>"P:embryo development", "GO:0015979"=>"P:photosynthesis", "GO:0016192"=>"P:vesicle-mediated transport", "GO:0019748"=>"P:secondary metabolic process", "GO:0021700"=>"P:developmental maturation", "GO:0022607"=>"P:cellular component assembly", "GO:0022618"=>"P:ribonucleoprotein complex assembly", "GO:0030154"=>"P:cell differentiation", "GO:0030198"=>"P:extracellular matrix organization", "GO:0030705"=>"P:cytoskeleton-dependent intracellular transport", "GO:0032196"=>"P:transposition", "GO:0034330"=>"P:cell junction organization", "GO:0034641"=>"P:cellular nitrogen compound metabolic process", "GO:0034655"=>"P:nucleobase-containing compound catabolic process", "GO:0040007"=>"P:growth", "GO:0040011"=>"P:locomotion", "GO:0042254"=>"P:ribosome biogenesis", "GO:0042592"=>"P:homeostatic process", "GO:0043473"=>"P:pigmentation", "GO:0044281"=>"P:small molecule metabolic process", "GO:0044403"=>"P:symbiosis, encompassing mutualism through parasitism", "GO:0048646"=>"P:anatomical structure formation involved in morphogenesis", "GO:0048856"=>"P:anatomical structure development", "GO:0048870"=>"P:cell motility", "GO:0050877"=>"P:neurological system process", "GO:0051186"=>"P:cofactor metabolic process", "GO:0051276"=>"P:chromosome organization", "GO:0051301"=>"P:cell division", "GO:0051604"=>"P:protein maturation", "GO:0055085"=>"P:transmembrane transport", "GO:0061024"=>"P:membrane organization", "GO:0065003"=>"P:macromolecular complex assembly", "GO:0071554"=>"P:cell wall organization or biogenesis", "GO:0071941"=>"P:nitrogen cycle metabolic process");

# Summary
my %STATS = ("Annotations"=>0, "a1"=>0, "a2"=>0, "a3"=>0, "a123"=>0, "a12"=>0);
foreach my $type (@TYPES) {
  $STATS{$type} = 0;
}
my %SLIMS;
my %PATHS;
my %KWS;
my $N = 0; # number of query sequences

# Evidence codes
my @ECO = qw(0000501 0000203 0000209 0000348 0000350 0000347 0000331 0000213 0000246 0000254 0000313 0000259 0000256 0000258 0000210 0000211 0000248 0000265 0000249 0000251 0000261 0000263 0000332);
my $GOEV = ":IEA";

# Gather input parameters
my $ANNOTATOR = 123;     # Annotator program (1=UniProt sequences; 2=orthologues; 3=by clustering)
my $QUERY_FILE;          # input fasta file
my $DAT_FILE;            # target dat database (UniProt)
my $FASTA_FILE;          # target fasta database (it is created from dat_file if it does exist)
my $BLAST_FILE;          # blast report
my $EXTENSION;           # extension for the output file
my $ROST = 20;           # Rost curve + (ROST_params)
my $PV = 0.1;            # p-value threshold (default 0.1)
my $TRAINING = 0;        # when it is equal to 1, the target ID is not gathered if it is equal to query ID
my $NOEMPTY;             # when it is equal to 1, unannotated IDs are reported
my $BLAST1 = "blastp";   # blast program for the first stage
my $BLAST2 = "blastp";   # blast program for the second stage
my $ANNOT_FILE;          # annotation final file
my $ID_UNIPROT = 90;     # minimal ID for considering sequence in UniProt database
my $ID_ORTHOLOGUE = 75;  # minimal ID for considering a ortholog sequence
my $COV_UNIPROT = 90;    # minimal coverage for considering sequence in UniProt database
my $COV_ORTHOLOGUE = 80; # minimal coverage for considering a ortholog sequence
my $EXTEND = 0;          # Extended output (1=add GO term names)
my $QUALITY = 0;         # On/Off quality mode (1=unreviewed annotations are avoided by their evidence code)
my $NOPRED = 0;          # On/Off no_predicted mode (1=predicted proteins are removed from the database)
my $SOURCE = 0;          # On/Off source columns (1=used annotator and IDs are reported)
my $UNIREF = 1;          # Database type (1=UniRef/0=UniProt)
my $GOSLIM = 0;          # GO Slim terms in the final report (1=Include GO Slim terms)
my $CPUS = 1;            # Number of threads for Blast runnnings
my $MAX_SCORE = 14;      # Maximum GeneName_Description Score
my $LOW_COMP = "no";     # Low complexity filter for Blast (1=yes/0=no)

my $u; # u=1 when user change identity or coverage
for (my $x = 0; $x <= $#ARGV; $x++) {
  if ($ARGV[$x] eq "--help") {
    &error(1);
  } elsif ($ARGV[$x] eq "-a") {
    $ANNOTATOR = $ARGV[$x+1];
    $EXTENSION .= "_a$ANNOTATOR";
  } elsif ($ARGV[$x] eq "-i") {
    $QUERY_FILE = $ARGV[$x+1];
    $ANNOT_FILE = $QUERY_FILE;
    $ANNOT_FILE =~ s/\.[A-Za-z0-9_-]+$/_/;
  } elsif ($ARGV[$x] eq "-d") {
    $DAT_FILE = $ARGV[$x+1];
    $FASTA_FILE = $DAT_FILE;
    $FASTA_FILE =~ s/\.dat/\.fasta/;
  } elsif ($ARGV[$x] eq "-b") {
    $BLAST_FILE = $ARGV[$x+1];
  } elsif ($ARGV[$x] eq "-r") {
    $ROST = $ARGV[$x+1] || &error(1);
    $EXTENSION .= "_r$ROST";
  } elsif ($ARGV[$x] eq "-p") {
    $PV = $ARGV[$x+1] || &error(1);
    $EXTENSION .= "_p$PV";
  } elsif ($ARGV[$x] eq "-cov1") {
    $COV_UNIPROT = $ARGV[$x+1] || &error(1);
    $EXTENSION .= "_cov1$COV_UNIPROT";
    $u = 1;
  } elsif ($ARGV[$x] eq "-cov2") {
    $COV_ORTHOLOGUE = $ARGV[$x+1] || &error(1);
    $EXTENSION .= "_cov2$COV_ORTHOLOGUE";
    $u = 1;
  } elsif ($ARGV[$x] eq "-id1") {
    $ID_UNIPROT = $ARGV[$x+1] || &error(1);
    $EXTENSION .= "_id1$ID_UNIPROT";
    $u = 1;
  } elsif ($ARGV[$x] eq "-id2") {
    $ID_ORTHOLOGUE = $ARGV[$x+1] || &error(1);
    $EXTENSION .= "_id2$ID_ORTHOLOGUE";
    $u = 1;
  } elsif ($ARGV[$x] eq "-num_threads") {
    $CPUS = $ARGV[$x+1] || &error(1);
  } elsif ($ARGV[$x] eq "-go") {
    $EXTEND = 1;
    $EXTENSION .= "_go";
    next;
  } elsif ($ARGV[$x] eq "-filter") {
    $LOW_COMP = "yes";
    $EXTENSION .= "_filter";
  } elsif ($ARGV[$x] eq "-nucl") {
    $BLAST1 = "blastx";
    $BLAST2 = "tblastn";
    next;
  } elsif ($ARGV[$x] eq "-training") {
    $TRAINING = 1;
    $EXTENSION .= "_training";
    next;
  } elsif ($ARGV[$x] eq "-noempty") {
    $NOEMPTY = 1;
    $EXTENSION .= "_noempty";
    next;
  } elsif ($ARGV[$x] eq "-nopred") {
    $NOPRED = 1;
    $EXTENSION .= "_nopred";
    next;
  } elsif ($ARGV[$x] eq "-quality") {
    $QUALITY = 1;
    $EXTENSION .= "_quality";
    next;
  } elsif ($ARGV[$x] eq "-source") {
    $SOURCE = 1;
    $EXTENSION .= "_source";
    next;
  } elsif ($ARGV[$x] eq "-uniprot") {
    $UNIREF = 0;
    next;
  } elsif ($ARGV[$x] eq "-goslim") {
    $GOSLIM = 1;
    push @TYPES, "GOSLIM";
    $EXTENSION .= "_goslim";
    next;
  } else {
    undef $QUERY_FILE; # Force parameter error
    last; 
  }
  $x++;
}

# Check input parameters
if (!$u && $BLAST1 eq "blastx") { # identity and coverage >=51% when nucleotide sequences
  $ID_UNIPROT = 51 if $ID_UNIPROT > 51;
  $ID_ORTHOLOGUE = 51 if $ID_ORTHOLOGUE > 51;
  $COV_UNIPROT = 51 if $COV_UNIPROT > 51;
  $COV_ORTHOLOGUE = 51 if $COV_ORTHOLOGUE > 51;
  $u = "Nucleotide queries selected: identity and coverage thresholds have been reduced to 51%\n";
}
if (!$QUERY_FILE || !$DAT_FILE || $ANNOTATOR !~ /^1?2?3?$/) {
  &error(1);
} elsif ($ROST < 0 || $ROST > 100 || $PV < 1e-100 || $PV > 1) {
  &error ("Rost parameter must be between 0-100, p-value threshold between 1e-100-1, and low-complexity filter F or T");
} elsif ($ID_UNIPROT < 1 || $ID_UNIPROT > 100 || $ID_ORTHOLOGUE < 1 || $ID_ORTHOLOGUE > 100 || 
         $COV_UNIPROT < 1 || $COV_UNIPROT > 100 || $COV_ORTHOLOGUE < 1 || $COV_ORTHOLOGUE > 100) {
  &error ("Identity and Coverage thresholds should be >=1 and <=100");
} elsif (! -e $QUERY_FILE) {
  &error ("$QUERY_FILE file does not exist");
} elsif (! -e $DAT_FILE && ! -e $FASTA_FILE) {
  &error ("$DAT_FILE file does not exist");
}
# Check if programs exist
my $null;
$null = `$BLAST1 -help`; &error("Check if ncbi-blast+ package is installed") if $? == -1;
$null = `$BLAST2 -help`; &error("Check if ncbi-blast+ package is installed") if $? == -1;
$null = `blastdbcmd -help`; &error("Check if ncbi-blast+ package is installed") if $? == -1;
$null = `makeblastdb -help`; &error("Check if ncbi-blast+ package is installed") if $? == -1;
if (!$UNIREF) {
  $null = `blastclust --help`; &error("Check if blastclust from blast2 package is installed") if $? == -1; 
}
undef $null;

# Check for Blosum62 file
if (! -e "BLOSUM62") { &create_blosum62 }

###########
# PROGRAM #
###########

# Output file with a header and the annotations
print <<HEADER;
  __      __    
 (_  _  _  _) _ 
 __)|||(_|__)_) 
HEADER
print "\nStarting Sma3s.v2 (annotator $ANNOTATOR)\n";
print $u if ($u && $u ne 1); # threshold reduced to 51%
(my $DAT_BRIEF = $DAT_FILE) =~ s/\.dat$//;
$DAT_BRIEF =~ s/.*[\/\\]//; # remove folder symbols
if (!$UNIREF) { 
  $ANNOT_FILE .= join "", map { substr ucfirst $_, 0, 3 } split/_/, $DAT_BRIEF;
} else {
  $DAT_BRIEF =~ s/\..+//g;
  $ANNOT_FILE .= $DAT_BRIEF;
}
$BLAST_FILE = $ANNOT_FILE . ".blast" if (!$BLAST_FILE);
$ANNOT_FILE .= $EXTENSION if $EXTENSION;
$ANNOT_FILE .= ".tsv";
open FILE, ">$ANNOT_FILE";
print FILE "#ID";
foreach my $type (@TYPES) {
  if ($type eq "GO" && $EXTEND) {
    print FILE "\tGO(P)ID\tGO(P)NAME\tGO(F)ID\tGO(F)NAME\tGO(C)ID\tGO(C)NAME";
  } else {
    print FILE "\t$type";
  }
}
print FILE "\tANNOTATOR\tUSED_UNIPROT_SEQUENCES" if ($SOURCE);
print FILE "\n";

# Create filtered Query FASTA
my $nucl; # for checking if query are nucleotide sequences and 'nucl' was not selected
my @id_q; # for makeblastcmd
my %MAPID; # Mapping original ID versus sXX
my $QUERY_FILEo = $QUERY_FILE; # Original query FASTA
$QUERY_FILE .= 2;
open QUERYo, $QUERY_FILEo || &error ("Problem opening $QUERY_FILEo");
open QUERY, ">$QUERY_FILE" || &error ("Problem creating $QUERY_FILE");
while (<QUERYo>) {
  chomp;

  if (/^>(.*)/) {
    $N++;
    my $id = "s$N";
    print QUERY ">$id\n";
    $MAPID{$id} = $1;
    push @id_q, $id;
  } else {
    if ($N <= 5 && $BLAST1 ne "blastx") { # check if queries are nucleotide sequences and 'nucl' was not selected
      if ($N == 5 && $nucl =~ /^[ACGTUN]+$/i) {
        unlink glob "$QUERY_FILE.*";
        &error ("First 4 queries seems to be nucleotide sequences, but 'nucl' argument was not selected");
      }
      $nucl .= $_;
    }
    print QUERY "$_\n";
  }
}
undef $nucl;
close QUERY;
close QUERYo;

# Blast & Blastclust parameters
my $params    = "-db $FASTA_FILE -evalue 1e-6 -seg $LOW_COMP -max_target_seqs 250 -out $BLAST_FILE -num_threads $CPUS -outfmt \"6 qseqid sseqid pident qcovs length evalue\"";
my $params2   = "-db $QUERY_FILE -evalue 1e-6 -seg $LOW_COMP -max_target_seqs 1 -outfmt \"6 sseqid\"";
my $bc_params = "-p T -L .95 -b F -S 95";

# Creating FASTA & ANNOT file
(my $REFDB = $FASTA_FILE) =~ s/\.fasta/\.annot/; # ANNOT_FILE from .fasta
if (! -e $FASTA_FILE || ! -e $REFDB) {
  print "Creating Fasta and Annot from UniProt file\n";
  &createFastaAnnotFromUniprot($DAT_FILE, $FASTA_FILE, $REFDB);
}

# Run formatdb for databases
if (! -e "$FASTA_FILE.psq" && ! -e "$FASTA_FILE.00.psq") { # Database
  print "Indexing $FASTA_FILE file\n";
  &run_formatdb($FASTA_FILE, "prot");
}
if ( ($BLAST1 eq "blastp" && ! -e "$QUERY_FILE.psq" && ! -e "$QUERY_FILE.00.psq") || ($BLAST1 eq "blastx" && ! -e "$QUERY_FILE.nsq" && ! -e "$QUERY_FILE.00.nsq") ) { # Query
  print "Indexing $QUERY_FILEo file\n";
  if ($BLAST1 eq "blastx") {
    &run_formatdb($QUERY_FILE, "nucl");
  } else {
    &run_formatdb($QUERY_FILE, "prot");
  }
}

# Run FIRST BLAST (only if blast_file does not exist)
#  and gather IDs
my $report_obj;
my %id_flag;
if (! -e $BLAST_FILE) {
  print "Running Blast\n";
  my $null = `$BLAST1 -query $QUERY_FILE $params`;
}

# Gather Blast, including a flag to free RAM from annots 
my @id_s; # for makeblastcmd
my %br; # blast reports
my $id_old = ""; # for not repeating subjects
open BLAST, $BLAST_FILE;
while (<BLAST>) {
  chomp;
  
  my ($id, $id2, @line) = split /\t/;
  my $id3 = "$id$id2";
  
  next if (($TRAINING && $MAPID{$id} eq $id2) || $id3 eq $id_old); # itself is avoided for training & also repeated subject for the same query
  push @id_s, $id2 if (!$id_flag{$id2});
  $id_flag{$id2} = 1;
  my $line = "$id\t$id2\t" . join "\t", @line;
  push @{$br{$id}}, $line;
  $id_old = $id3;
}
close BLAST;

# Gather lengths for query and hit sequences
my %lenq;
my %lens;
print "Making initial calculations\n";
my $lengths_q = &get_seq($QUERY_FILE, \@id_q, "len");
my $lengths_s = &get_seq($FASTA_FILE, \@id_s, "len");
my (@lengths_q) = split/\n/, $lengths_q;
my (@lengths_s) = split/\n/, $lengths_s;
for (my $x = 0; $x <= $#id_q; $x++) {
  $lenq{$id_q[$x]} = $lengths_q[$x];
}
for (my $x = 0; $x <= $#id_s; $x++) {
  $lens{$id_s[$x]} = $lengths_s[$x];
}
undef @id_q; undef @lengths_q; undef $lengths_q;
undef @id_s; undef @lengths_s; undef $lengths_s;

# Biological enrichment for annotator 3
my %FREQ;
my %N_ANNOT;
  
# Gather reference annotations
my %ANN; # Annotations for all hits in Blast Report
open IN, "$REFDB"; # REF annotation file for enrichment
print "Reading annotations from UniProt reference file\n";
while (<IN>) {
  chomp $_;

  my ($id, $annot, @annot) = split /\t/; # annot=provisional annotation (begin with geneName score, 2nd column)
  for (my $i = 0; $i <= $#TYPES; $i++) { # Calculate reference annotations for enrichment and gather annotation for each ID
    my $type = $TYPES[$i];
    my @term_final; # terms finally kept for a ID
    
    if (!$annot[$i]) { # if no terms in this annotation type
      $annot .= "\t";
      next;
    }

    my (@term) = split /\;/, $annot[$i];
    TERM: foreach my $term (@term) {
      if ($QUALITY) { # Quality mode (only reviewed annotations)
        if ($type eq "GO" && $term =~ /\{(.+?)\}/) {
          next if ($1 =~ /$GOEV/);
        } elsif (($i != 3 && $i != 4) && $term =~ /\{(.+?)\}/) { # only for annotations with ECO numbers
          my (@eco) = split/\|\|/, $1;
          foreach my $eco (@eco) {
            $eco =~ s/ ?ECO://;
            next TERM if (grep(/$eco/,@ECO));
          }
        }
      }

      if ($EXTEND && $type eq "GO") { # remove evidence codes (except in Extended mode)
        $term =~ s/^GO://; $term =~ s/\{/:/; $term =~ s/:\w{2,3}\}//; # edit: removed \w{3}
      } else {
        $term =~ s/\{.+?\}//; # remove evidence codes
      }
      push @term_final, $term;
      $FREQ{$type}{$term}++;
      $N_ANNOT{$type}++;
    }
    $annot .= "\t" . join ";", @term_final;
  }
  $annot =~ s/\t+$//g; # remove tabs in the latter columns
  $ANN{$id} = $annot if ($id_flag{$id}); # final annotation
}
close IN;
undef %id_flag;

# ANNOTATION
my @queries;
print "Reading sequences to annotate\n";
open FASTA, $QUERY_FILE; # Reading sequences to annotate
while (<FASTA>) {
  chomp;

  next unless /^>/;
  my $id = $_;
  $id =~ s/>//;
  push @queries, $id;
}
close FASTA;

print "Running Annotator $ANNOTATOR\n";
QUERY: foreach my $q (@queries) { # Go through the FIRST BLAST
  my ($o1, $o2);              # The best output from the 3 annotators
  my @s1 = (0,0); my $s2 = 0; # The best scores from the 3 annotators
  if (!$br{$q}[0]) {       # hits not found in Blast report
    print FILE "$MAPID{$q}\n" if (!$NOEMPTY); 
    next;
  }
  
  my($name, $name2, $id_hsp, $qc_q, $len_hsp, $evalue)  = split /\t/, $br{$q}[0]; # Blast parameters
  
  # Only ANNOTATOR 1
  ##################
  if ($ANNOTATOR =~ /1/) {
    my $pv_hsp = &pvalue_from_evalue ($evalue);
    my $qc_s = &calculate_qc_subject ($lenq{$name}, $qc_q, $lens{$name2});
    
    if ($id_hsp >= $ID_UNIPROT && $qc_s >= $COV_UNIPROT && $pv_hsp <= $PV) { # alignment: 90% id + 90% qcoverage + pvalue
      $o1 = ""; # initialize annotation
      my ($s, $gn1, $de1, @remaining) = split/\t/, $ANN{$name2};
      (@s1) = split/,/, $s; # scores
      if ($gn1) { ($gn1) = split/;/, $gn1; $o1 .= $gn1; } $o1 .= "\t"; # only first GN
      if ($de1) { ($de1) = split/;/, $de1; $o1 .= $de1; } $o1 .= "\t"; # only first DE
      $o1 .= join "\t" , @remaining; # functional annotation

      if ($ANNOTATOR == 1) { # only A1 selected
        print FILE &create_annot($MAPID{$name}, $o1, "a1", $name2);
        next QUERY; 
      }
    } elsif ($ANNOTATOR == 1) {
      print FILE "$MAPID{$name}\n" if (!$NOEMPTY); 
      next QUERY; 
    }
  }
    
  # Only ANNOTATOR 2
  ##################  
  if ($ANNOTATOR =~ /2/) {
    foreach my $br (@{$br{$q}}) {
      next if $s1[0] == 1 || $s1[1] == $MAX_SCORE; # next if already swiss-prot or maximum gn_de score
      my($name, $name2, $id_hsp, $qc_q, $len_hsp, $evalue)  = split /\t/, $br; # Blast parameters
      
      # Score filter (only evaluating new candidate if higher score than previous
      my ($s, $gn1, $de1, @remaining) = split/\t/, $ANN{$name2};
      my (@s3) = split/,/, $s; # scores
      
      next unless $s3[0] == 1 || $s3[1] > $s1[1]; 
      
      my $pv_hsp = &pvalue_from_evalue ($evalue);
      my $qc_s = &calculate_qc_subject ($lenq{$name}, $qc_q, $lens{$name2});
    
      # Alignment: % id or rost + qcoverage + pvalue threshold
      my $threshold = $ID_ORTHOLOGUE;
      $threshold = &calculate_rost ($ROST,$len_hsp) if (!$u && $BLAST1 ne "blastx"); # use Rost only if id/cov have not been changed by the user
      if ($id_hsp >= $threshold && $qc_s >= $COV_ORTHOLOGUE && $pv_hsp <= $PV) {
        
        # Run SECOND BLAST
        my $fasta2 = &get_seq($FASTA_FILE, $name2, "seq");
        my $temp_f = "./.temp_$name2\_" . rand(1000000);
        open TEMP, ">$temp_f" || &error ("Problem creating files in this folder");
        print TEMP ">$name2\n$fasta2";
        close TEMP;
        my $name3 = `$BLAST2 -query $temp_f $params2`;
        chomp $name3;
        unlink $temp_f;
       
        # Check Reciprocal Blast
        if ($name3 eq $name) {
          $o1 = ""; # initialize annotation
          if ($gn1) { ($gn1) = split/;/, $gn1; $o1 .= $gn1; } $o1 .= "\t"; # only first GN
          if ($de1) { ($de1) = split/;/, $de1; $o1 .= $de1; } $o1 .= "\t"; # only first DE
          $o1 .= join "\t" , @remaining; # functional annotation
          
          @s1 = @s3;
          if ($ANNOTATOR == 2 || $ANNOTATOR == 12) { # only A2 selected
            print FILE &create_annot($MAPID{$name}, $o1, "a2", $name2);
            next QUERY; 
          }
        }
      } elsif ($ANNOTATOR == 2 || $ANNOTATOR == 12) {
        print FILE "$MAPID{$name}\n" if (!$NOEMPTY); 
        next QUERY; 
      }
    }
  }
      
  # ANNOTATOR 3: Rost equation threshold
  ######################################
  my @id_fastas = (); # Gather Fastas from Blast report
  if ($ANNOTATOR =~ /3/) {
    foreach (@{$br{$q}}) {
      my ($name, $name2, $id_hsp, $qc_q, $len_hsp, $evalue)  = split /\t/;
      my $threshold = &calculate_rost ($ROST,$len_hsp);
      next unless ($id_hsp > $threshold);

      # Collect fasta sequences with right parameters
      push @id_fastas, $name2;
    } 

    # Create clusters and gather all annotation for each one
    if ($#id_fastas > 0) { # only annotate when FASTA has more than one sequence
      my $fastas = &get_seq($FASTA_FILE, \@id_fastas, "fasta");
      my @clusters = @id_fastas;
      if ($UNIREF == 0) {
        my $clusters = &run_blastclust($fastas, $bc_params, $name);
        @clusters = @{$clusters};
      }
      my %freq = (); # final annotation frequency for each annotation
      my %n_annot;
      foreach my $cluster (@clusters) {
        my (@ids) = split(/\s/, $cluster);
        my %annot = ();
        foreach my $id (@ids) { # Ids in a cluster
          for (my $i = 0; $i <= $#TYPES; $i++) { # annotation columns
            my $type = $TYPES[$i];
            my ($terms) = (split /\t/, $ANN{$id})[$i+1];

            next unless ($terms); # next if term type is empty
            foreach my $term (split/;/, $terms) {
              my $term_m = &filter_cs ($term);
              next if (grep/^$term_m$/,@{$annot{$type}}); # remove redundancy
              push @{$annot{$type}}, $term;
            }
          }
        }
          
        # Count annotations from the different clusters
        foreach my $type (keys %annot) {
          foreach my $annot (@{$annot{$type}}) { # each particular annot
            $freq{$type}{$annot}++;
            $n_annot{$type}++;
          }
        }
      }

      # Recover only annotations with a minimal support among the different clusters
      foreach my $type (@TYPES) {  
        my $min_pv = 10;
        my @output;  # provisional output
        $o2 .= "\t"; # separator
        foreach my $annot (sort keys %{$freq{$type}}) {
          next if (!$annot); # next if empty annotation
          
          my $pvalue = &ComputeHyperPValue ($freq{$type}->{$annot}, $n_annot{$type}, $FREQ{$type}->{$annot}, $N_ANNOT{$type});
          
          next unless ($pvalue <= $PV);
          if ($type eq "GENENAME") { # One only Gene Name (the most informative one)
            my $s_gn  = &calculate_gnscore($annot);
            if ($s_gn > $s2) { # if current is more informative than the previous
              $s2 = $s_gn;
              $output[0] = $annot;
            }
          } elsif ($type eq "DESCRIPTION") { # One only DEscription
            next unless ($pvalue < $min_pv);
            $min_pv = $pvalue;
            $output[0] = $annot;
          } else { # Remaining fields
            push @output, $annot;
          }
        }
        $o2 .= join ";", @output;
      }
      $o2 =~ s/^\t//;
    }
  
    # If only A3 was selected  
    if ($ANNOTATOR == 3) {
      print FILE &create_annot($MAPID{$name}, $o2, "a3", (join ";", @id_fastas)) if ($o2);
      print FILE "$MAPID{$name}\n" if (!$o2 && !$NOEMPTY);
      next QUERY;
    } 
  }

  # Final OUTPUT (A123): $o1=A1|A2; $o2=A3
  if (!$o1 && !$o2) { # sequence not annotated
    print FILE "$MAPID{$name}\n" if (!$NOEMPTY);
    next QUERY;
  }
  
  my $o3; # A123
  if    (!$o2) { print FILE &create_annot($MAPID{$name}, $o1, "a12", $name2) } # A1|A2
  elsif (!$o1) { print FILE &create_annot($MAPID{$name}, $o2, "a3", (join ";", @id_fastas)) } # A3
  else { # A1|A2 & A3 (merging annotations):Real A123
    for (my $x = 0; $x <= $#TYPES; $x++) {
      my ($column1) = (split/\t/, $o1)[$x];
      my ($column2) = (split/\t/, $o2)[$x];

      if ($x == 0 || $x == 1) { # GN & DE: check GN score and select A1|A2 or A3
        $o3 .= "\t";
        if ($s1[1] >= $s2 && $column1) { $o3 .= $column1 } elsif ($column2) { $o3 .= $column2 } # only if it exists
        next;
      }

      my (@o1, @o2); # any of them could be empty
      @o1 = split/;/, $column1 if $column1;
      @o2 = split/;/, $column2 if $column2;
      my @oo = (@o1, @o2);                                                                                                                                   
      my %hash = map { $_, 1 } @oo; @oo = keys %hash; # remove redundancy
      $o3 .= "\t";
      $o3 .= join ";", @oo;
    }
    $o3 =~ s/^\t//; # remove first tab (previous to gene_name)
    print FILE &create_annot($MAPID{$name}, $o3, "a123", (join ";", @id_fastas));
  }
}
close FILE;

# Print summary
(my $STAT_FILE = $ANNOT_FILE) =~ s/\.tsv$/_summary.tsv/;
open STAT, ">$STAT_FILE";
print STAT "#Annotation summary\n";
print STAT "Number of query sequences:\t$N\n";
foreach my $stat ("Annotations", @TYPES) { # Annotations
  print STAT "With $stat\t$STATS{$stat}\n" unless $stat eq "GOSLIM";
}
foreach my $stat ("a1", "a2", "a3", "a12", "a123") { # Annotator
  print STAT "Annotator $stat\t$STATS{$stat}\n";
}
print STAT "\n";
if ($GOSLIM) {
  my %cat = ("P" => "Biological process", "C" => "Cellular component", "F" => "Molecular function");
  print STAT "#GO Slim\n";
  foreach my $g (keys %SLIMS) {
    print STAT "#Category \"$cat{$g}\"\n";
    foreach my $slim (sort { $SLIMS{$g}{$b} <=> $SLIMS{$g}{$a} } keys %{$SLIMS{$g}}) {
      my ($name) = (split/:/, $SLIM{$slim})[1];
      print STAT "$slim\t$name\t$SLIMS{$g}{$slim}\n";
    }
    print STAT "\n";
  }
}
print STAT "#UniProt Pathways\n" if %PATHS;
foreach my $path (sort { $PATHS{$b} <=> $PATHS{$a} } keys %PATHS) { # UniProt Pathways
  print STAT "$path\t$PATHS{$path}\n";
}
print STAT "\n" if %PATHS;
print STAT "#UniProt Keyword categories\n";
foreach my $cat (sort keys %KWS) { # UniProt Keyword categories
  print STAT "#Category \"$cat\"\n";
  foreach my $kw (sort { $KWS{$cat}{$b} <=> $KWS{$cat}{$a} } keys %{$KWS{$cat}}) {
    print STAT "$kw\t$KWS{$cat}{$kw}\n";
  }
  print STAT "\n";
}
close STAT;

# The end
print "\nAnnotation file created (open with a spreadsheet): $ANNOT_FILE\n";
print "Summary file created (open with a spreadsheet): $STAT_FILE\n\n";
print "Please, don't forget citing the paper: http://www.ncbi.nlm.nih.gov/pubmed/24501397\n\n";
exit;

##############
# SUBRUTINES #
##############

# pvalue_from_evalue
# params: evalue
################
sub pvalue_from_evalue () {
  my ($evalue) = @_;
  return 1 - exp(-$evalue);
}

# calculate_qc_subject
# params: len_query, querycoverage_query, length_subject
# description: calculate query coverage for subject from query
##############################################################
sub calculate_qc_subject () {
  my ($len_query, $qc_query, $len_subject) = @_;
  return $len_query * $qc_query / $len_subject;
}

# calculate_rost
# params: rost constant, alignment length
# description: calculate the %identity threshold for a given length
###################################################################
sub calculate_rost () {
  my ($ROST, $len) = @_;
  my $id = $ROST + (480 * ($len ** (-0.32 * (1+exp(-$len/1000))) ) );
  return $id;
}

# create_annot
# params: annotation_line, a123, used_IDs (it also needs @TYPES, EXTEND, EMPTY & $STATS)
# description: prepare the annotation line for an ID and calculate/add summary values
########################################################################################
sub create_annot () {
  my ($ANNOT, $o, $a, $used_ids) = @_; # ANNOT is initially the ID
 
  unless ($o =~ /[^\t]/) { # if annotation is empty
    return "$ANNOT\n" if (!$NOEMPTY);
    return;
  }
  
  my (@ann) = split/\t/, $o; # annotations
  $STATS{"Annotations"}++;
  $STATS{$a}++;
  
  for (my $i = 0; $i <= $#TYPES; $i++) {
    my $type = $TYPES[$i];
    if ($ann[$i]) {
      $STATS{$type}++ unless $type eq "GOSLIM";

      # Extended mode (separate GO)
      if ($type eq "GO" && $EXTEND) {
        my %pfc_id = ("P"=>[], "F"=>[], "C"=>[]); # classify GO by class
        my %pfc_name = ("P"=>[], "F"=>[], "C"=>[]); # classify GO by class
        my (@goes) = split/;/, $ann[$i];
        foreach my $go (@goes) {
          my ($goid, $g, $name) = split/:/, $go;
          push @{$pfc_id{$g}}, "GO:$goid"; 
          push @{$pfc_name{$g}}, "$name"; 
        }
        foreach my $group ("P", "F", "C") {
          $ANNOT .= "\t" . (join ";", @{$pfc_id{$group}}) . "\t" . (join ";", @{$pfc_name{$group}});
        }
      
      } else {
        $ANNOT .= "\t$ann[$i]";
      }

      # Count Pathways
      if ($type eq "PATHWAY") {
        my (@paths) = split/;/, $ann[$i];
        map { $PATHS{$_}++ } @paths;
      } elsif ($type eq "KEYWORD") {
        my (@kws) = split/;/, $ann[$i];
        foreach my $kw (@kws) {
          $KWS{"Biological_process"}{$kw}++ if (grep/^$kw$/, @KW_Biological_process);
          $KWS{"Cellular_component"}{$kw}++ if (grep/^$kw$/, @KW_Cellular_component);
          $KWS{"Developmental_stage"}{$kw}++ if (grep/^$kw$/, @KW_Developmental_stage);
          $KWS{"Disease"}{$kw}++ if (grep/^$kw$/, @KW_Disease);
        }
      } elsif ($type eq "GOSLIM") { # GOSlim groups
        my (@slims) = split/;/, $ann[$i];
        foreach my $goslim (@slims) {
          my ($g, $name) = split/:/, $SLIM{$goslim};
          $SLIMS{$g}{$goslim}++; 
        }
      }
    } else {
      $ANNOT .= "\t";
      $ANNOT .= "\t\t\t\t\t" if ($type eq "GO" && $EXTEND); # add more tabs if EXTENDing GO column
    }
  }

  $ANNOT .= "\t$a\t$used_ids" if ($SOURCE);

  return "$ANNOT\n";
}

# get_seq
# params: DB, ID list
# description; get sequence from blast+ index
#############################################
sub get_seq() {
  my ($db, $ids, $format) = @_;
  my %f = ("seq" => "%s", "len" => "%l", "fasta" => ">%a/%s");
  if ($format eq "seq") {
    return `blastdbcmd -db $db -entry $ids -outfmt \"$f{$format}\"`;
  } else {
    my $temp = "./.temp" . rand(1000000) . ".ids";
    open IDS, ">$temp";
    print IDS join "\n", @$ids;
    close IDS;
    
    my $out = `blastdbcmd -db $db -entry_batch $temp -outfmt \"$f{$format}\"`;
    if (!$out) { 
	    &error ("Problem with blastdbcmd. It could be due to an empty or erroneous Blast report file, indexed fasta, or low RAM memory"); 
	  }
    unlink $temp;
    $out =~ s/\//\n/g; # for windows compatibility
    return $out;
  }
}

# run_blastclust
# params: fasta_seqs params(eg: -p T -L .95 -b F -S 95)
# description: run Blastclust from Blast package (clustering
############################################################
sub run_blastclust () {
  my ($fasta, $params, $name) = @_;
  
  # Random temporal file
  my $temp = "./temp_$name.fas" . rand(1000000);

  # Create FASTA file
  open (FASTA, ">$temp") || &error ("Problem creating $temp");
  print FASTA $fasta;
  close FASTA;
  
  # Run Blastclust
  my (@out) = `blastclust $params -i $temp`;
  shift @out; # remove header
  
  # Remove FASTA file
  unlink $temp;
  
  return \@out;
}

# createFastaAnnotFromUniprot
# params: uniprot_file.dat
# description: create FASTA and annot files from a UniProt file
###############################################################
sub createFastaAnnotFromUniprot () {
  my ($dat, $fasta, $annot) = @_;
  my $e; # Entry
  my %n = ("Reviewed" => 1, "Unreviewed" => 0); # Swiss-Prot || TrEMBL
  my @STOP_DE = ("Uncharacterized protein", "Putative uncharacterized protein"); 
  my @STOP_KW = ("3D-structure", "Allosteric enzyme", "Complete proteome", "Genetically modified food", "Hybridoma", "Multifunctional enzyme", 
                 "Pharmaceutical", "ERV", "Direct protein sequencing", "Extinct organism protein", "Proteomics identification", 
                 "Reference proteome", "Alternative initiation", "Alternative splicing", "Polymorphism", "Alternative promoter usage");

  unlink glob "$fasta*" if (-e "$fasta.psq" || -e "$fasta.00.psq"); # Remove indexed fasta

  # Gather whole entry 
  open IN, $dat;
  open FASTA, ">$fasta";
  open ANNOT, ">$annot";
  while (<IN>) {
    $e .= $_;
    next unless /^\/\//;

    my $out;     # Output line
    my $gnS = 0; # GN/DE Score
    my $deS = 0; # DE Score
    
    $e =~ s/\|[^\}\,\|]+([\}\,])/$1/gs; # simplify ECO term
    $e =~ s/ \{ECO:/\{ECO:/g;
    $e =~ s/ECO:(\d+)\|\w+:\d+/ECO:$1\|\|/gs;
    $e =~ s/ECO:(\d+)[,;] ?/ECO:$1\|\|/gs;
    
    # ID (swiss-prot or trembl) & FASTA
    $e =~ /^ID   (.+?)\s+(\w+);/; # ID
    my $id = $1;
    my $db = $n{$2};
    my $fasta = ">$id\n";
    my @seq = $e =~ /\n     (.+)/g;
    my $seq = join "\n", @seq;
    $seq =~ s/[\n\s]//g;
    $fasta .= "$seq\n";

    # GeneName
    my @gn = $e =~ /\nGN   (.+)/g;
    my $gn = join "", @gn;
    $gn =~ s/\w+=//g;
    $gn =~ s/;$//;
    @gn = split/[;,] */, $gn;
    $out .= join ";", @gn;
    $out .= "\t";
    $gnS = &calculate_gnscore($gn[0]) if ($gn[0]); # use first GN for scoring

    # Description & EC numbers
    my $first_de;
    my @de = $e =~ /\nDE   (.+)/g;
    my $de = join "", @de;
    $de =~ s/;Flags:.+//;
    $de =~ s/\w+: //g;
    my @ec = $de =~ /EC=([0-9\.\-]+\{?E?C?O?:?[^\}\;]*\}?);?/g;
    $de =~ s/EC=[0-9\.]+\{?E?C?O?:?[^\}\;]*\}?;?//g;
    $de =~ s/\w+=//g;
    @de = split/; */, $de;
    foreach (@de) {
      my $clean = $_;
      $clean =~ s/\{.*\}//;
      my ($clean_m) = &filter_cs($clean);
      next if (grep/^$clean_m$/, @STOP_DE);

      $deS = 4;
      $out .= ";" unless !$first_de; 
      $first_de = 1;
      $out .= $_;
    }
    $out .= "\t";
    $out .= join ";", @ec;
    $out .= "\t";
    
    # GO terms
    my @go;
    foreach my $go ($e =~ /\nDR   GO; (GO:\d+; [PCF]:.+; \w+):/g) {
       $go =~ s/; /{/;
       $go =~ s/; /:/;
       $go .= "}";
       push @go, $go;
    }
    $out .= join ";", @go;
    $out .= "\t";

    # Keywords
    my $first_kw;
    if ($e =~ /\nKW   ([^\.]+)\./s) {
      my $kw = $1;
      if ($kw) {
        $kw =~ s/\n/ /g;
        $kw =~ s/KW   //g;
        my @kw = split /; ?/, $kw;
        foreach (@kw) { # remove STOP KW
          my $clean = $_;
          $clean =~ s/\{.*\}//;
          my ($clean_m) = &filter_cs($clean);
          next if (grep/^$clean_m/, @STOP_KW);
          $out .= ";" unless !$first_kw; 
          $first_kw = 1;
          $out .= $_;
        }
      }
      $out .= "\t";
    }
    
    # Pathways
    my @pw;
    foreach my $pw ($e =~ /\nCC   -\!- PATHWAY: (.+?)\nCC   -/sg) { # finish in CC   -----------
      $pw =~ s/\nCC      //g;
      my ($eco_pw);
      if ($pw =~ s/(\{ECO:.+\})//) {
        $eco_pw = $1;
      }
      ($pw) = split/;/, $pw;
      $pw =~ s/\.//g;
      $pw .= $eco_pw if $eco_pw;
      push @pw, $pw;
    }
    my %hash = map { $_, 1 } @pw; @pw = keys %hash; # remove redundancy
    $out .= join ";", @pw;
    $out .= "\t";
    
    # Collect Evidence for protein (4:Predicted could be removed if user select this argument)
    my $pe = 0;
    $e =~ /\nPE   (\d):/;
    $pe = $1 if $NOPRED == 1;
    
    # End
    $gnS += $deS; # Add GeneName_score + Description_score
    
    
    if ((@go || $first_kw) && ($pe != 4)) {
      print ANNOT "$id\t$db,$gnS\t$out\n";
      print FASTA $fasta;
    } 
    $e = "";
  }
  close FASTA;
  close ANNOT;
  close IN;

  return;
}

# run_formatdb
# params: database, -p[T|F]
# run formatdb tool from blast package
######################################
sub run_formatdb () {
  my ($db, $p) = @_;

  my $null = `makeblastdb -dbtype $p -in $db -parse_seqids`;

  return;
}

# filter_cs
# params: string
# filter special characters
###########################
sub filter_cs () {
  my ($s) = @_;

  $s =~ s/([\+\*\.\/\\\[\]\(\)\^\$\?])/\\$1/g;

  return $s;
}

# calculate_gnscore
# params: geneName
# calculate gene name score
###########################
sub calculate_gnscore () {
  my ($g) = @_;
 
  return 0 if (!$g);
  
  my $gnS = 0;
  $g =~ s/\{.*//;
  my $gn_l = length $g;    
  $gnS += 4 if ($g !~ /[\._-]/);
  $gnS += 2 if (substr($g,0,1) =~ /[a-z]/);
  if ($gn_l <= 4) { $gnS += 4 } 
    elsif ($gn_l <= 5) { $gnS += 3 } 
    elsif ($gn_l <= 6) { $gnS += 2 } 
    elsif ($gn_l <= 8) { $gnS += 1 } 

  return $gnS;
}

# error
# params: 1=command line, or error text
# die and show an error message
#######################################
sub error () {
  my ($e) = @_;

die "\nError! $e\n\n" unless $e eq "1";

print <<HELP;

Please run as: $0 -i fasta_input_file -d target_database.fasta
 (target_database.dat should be used when using a new database)
 Two files will be created: one with the annotations, and another with a summary which can be used to generate graphics
 Both files can be opened in a spreadsheet program

 optional arguments
  -a  annotators to annotate query sequences
    default = 123 (1=sequences already in database, 2=orthologues, 3=exhaustive)
    (123 is recommended, but you should discard annotator 2, and select 13, if sequences coming from different organisms)

  -b  blast file (it is automatically created when it is not indicated)

  -p  pvalue threshold
    default = 0.1 (both for Blast and enrichment)

  -id1  identity threshold for annotator 1
    default = 90

  -id2  identity threshold for annotator 2
    default = rost formulae (it is changed to 75% or to the user value, if any of id1, id2, cov1, cov2 is modified)

  -cov1  database sequence coverage for annotator 1
    default = 90

  -cov2  database sequence coverage for annotator 2
    default = 80

  -r  'N' parameter in Rost formulae
    default = 20 (the higher it is, the more restrictive it will be)
    it is used in both 2 and 3 annotators

  -num_threads  number of threads (CPUs) to use for Blast
    default = 1

 argument without values (they are activated if you give the argument)
  -nucl     nucleotide query sequences (identity and coverage threshold will be reduced to 51, unless user indicates something else) 
  -filter   activate low complexity filter for Blast
  -source   two additional columns for each query in the results: used annotator, and used database sequences
  -go       the 3 GO categories will be set apart in the results, together with the GO names
  -uniprot  if you use a uniprot database (with .dat extension) instead of the default uniref database
  -goslim   the final report will include GO Slim terms (it is useful for comparing proteomes/transcriptomes)
  -noempty  unannotated query sequences will not appear in the results
  -nopred   predicted proteins from the database will be discarded
  -quality  low quality annotations (usually without manual curation) will be discarded
HELP

  die "\n";
}

# Create BLOSUM62 file
######################
sub create_blosum62 () {
  open BLOSUM, ">BLOSUM62";

print BLOSUM <<BLOSUM;
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
BLOSUM

        close BLOSUM;
  return;
}

# HyGe libraries
################
sub lgamma {
    my $a = shift;
    return &gammaln($a);
}

sub lchoose {
    my ($x, $y) = @_;
    die ("**lchoose $x $y") if( $x < 0 || $y < 0 || $x > $y );
    unless (exists $LChooseCache{"$x:$y"}) {
  $LChooseCache{"$x:$y"} = lgamma($y+1) - lgamma($x+1) - lgamma($y-$x+1);
    }
    return $LChooseCache{"$x:$y"};
}

sub AddLog {
    my ($x, $y) = @_;
    return $y if ($x == -$HUGE_VAL);
    return $x if ($y == -$HUGE_VAL);
    if ($x >= $y ) {
  $y -= $x;
    } else {
  my $t = $x;
  $x = $y;
  $y = $t - $x;
    }
    $y = $EPS if ($y < $EPS);
    return ($x + log(1+exp($y)));
}

sub ComputeHyperPValue {
    my ($k, $n, $K, $N) = @_;
    return 0 if ($N == 0);
    my $PVal = -$HUGE_VAL;

    for( ; $k >= 0 && $k <= $n && $k <= $K && ($n-$k) <= ($N-$K); $k ++) 
    {
  my $x = lchoose($k,$K) + lchoose($n-$k, $N-$K) - lchoose($n,$N);
  $PVal = AddLog($PVal, $x);
    }
    $PVal = 0 if ($PVal > 0);
    return exp($PVal);
}

# Gamma.pm library extract
##########################

use constant PI => 4 * atan2 1, 1;

sub gammaln {
    my $x = @_ ? shift : $_;

    if ( $x <= 0 ) {
        return log(PI) - log(sin(PI*$x)) - gammaln(1-$x);
    }
    return &_gammaln ($x);

}


sub _gammaln ($) {
#----------------------------------------------------------------------
    my $x = $_[0];
    my ($res, $corr, $xnum, $xden, $xm1, $xm2, $xm4);
#----------------------------------------------------------------------
#  Mathematical constants
#----------------------------------------------------------------------
    my $pnt68  = 0.6796875e0;
    my $sqrtpi = 0.9189385332046727417803297e0;
#----------------------------------------------------------------------
#  Machine dependent parameters
#----------------------------------------------------------------------
    my $xbig   = 2.55e305;
    my $xinf   = 1.79e308;
    my $eps    = 2.22e-16;
    my $frtbig = 2.25e76;
#----------------------------------------------------------------------
#  Numerator and denominator coefficients for rational minimax
#     approximation over (0.5,1.5).
#----------------------------------------------------------------------
    my $D1 = -5.772156649015328605195174e-1;
    my @P1 = (4.945235359296727046734888e0, 2.018112620856775083915565e2,
              2.290838373831346393026739e3, 1.131967205903380828685045e4,
              2.855724635671635335736389e4, 3.848496228443793359990269e4,
              2.637748787624195437963534e4, 7.225813979700288197698961e3);
    my @Q1 = (6.748212550303777196073036e1, 1.113332393857199323513008e3,
              7.738757056935398733233834e3, 2.763987074403340708898585e4,
              5.499310206226157329794414e4, 6.161122180066002127833352e4,
              3.635127591501940507276287e4, 8.785536302431013170870835e3);
#----------------------------------------------------------------------
#  Numerator and denominator coefficients for rational minimax
#     Approximation over (1.5,4.0).
#----------------------------------------------------------------------
    my $D2 = 4.227843350984671393993777e-1;
    my @P2 = (4.974607845568932035012064e0, 5.424138599891070494101986e2,
              1.550693864978364947665077e4, 1.847932904445632425417223e5,
              1.088204769468828767498470e6, 3.338152967987029735917223e6,
              5.106661678927352456275255e6, 3.074109054850539556250927e6);
    my @Q2 = (1.830328399370592604055942e2, 7.765049321445005871323047e3,
              1.331903827966074194402448e5, 1.136705821321969608938755e6,
              5.267964117437946917577538e6, 1.346701454311101692290052e7,
              1.782736530353274213975932e7, 9.533095591844353613395747e6);
#----------------------------------------------------------------------
#  Numerator and denominator coefficients for rational minimax
#     Approximation over (4.0,12.0).
#----------------------------------------------------------------------
    my $D4 = 1.791759469228055000094023e0;
    my @P4 = (1.474502166059939948905062e4, 2.426813369486704502836312e6,
              1.214755574045093227939592e8, 2.663432449630976949898078e9,
              2.940378956634553899906876e10, 1.702665737765398868392998e11,
              4.926125793377430887588120e11, 5.606251856223951465078242e11);
    my @Q4 = (2.690530175870899333379843e3, 6.393885654300092398984238e5,
              4.135599930241388052042842e7, 1.120872109616147941376570e9,
              1.488613728678813811542398e10, 1.016803586272438228077304e11,
              3.417476345507377132798597e11, 4.463158187419713286462081e11);
#----------------------------------------------------------------------
#  Coefficients for minimax approximation over (12, INF).
#----------------------------------------------------------------------
    my @C = (-1.910444077728e-03, 8.4171387781295e-04,
             -5.952379913043012e-04, 7.93650793500350248e-04,
             -2.777777777777681622553e-03, 8.333333333333333331554247e-02,
             5.7083835261e-03);
#----------------------------------------------------------------------
    my $y = $x;
    if (($y > 0) && ($y <= $xbig)) {
        if ($y <= $eps) {
            $res = -log($y);
        } elsif ($y <= 1.5) {
#----------------------------------------------------------------------
#  EPS < X <= 1.5
#----------------------------------------------------------------------
            if ($y < $pnt68) {
                $corr = -log($y);
                $xm1 = $y
            } else {
                $corr = 0;
                $xm1 = ($y - 0.5) - 0.5;
            }
            if (($y <= 0.5) || ($y >= $pnt68)) {
                $xden = 1;
                $xnum = 0;
                foreach my $i (0 .. 7) {
                    $xnum = $xnum*$xm1 + $P1[$i];
                    $xden = $xden*$xm1 + $Q1[$i];
                }
                $res = $corr + ($xm1 * ($D1 + $xm1*($xnum/$xden)));
            } else {
                $xm2 = ($y - 0.5) - 0.5;
                $xden = 1;
                $xnum = 0;
                foreach my $i (0 .. 7) {
                    $xnum = $xnum*$xm2 + $P2[$i];
                    $xden = $xden*$xm2 + $Q2[$i];
                }
                $res = $corr + $xm2 * ($D2 + $xm2*($xnum/$xden));
            }
        } elsif ($y <= 4) {
#----------------------------------------------------------------------
#  1.5 < X <= 4.0
#----------------------------------------------------------------------
            $xm2 = $y - 2;
            $xden = 1;
            $xnum = 0;
            foreach my $i (0 .. 7) {
                $xnum = $xnum*$xm2 + $P2[$i];
                $xden = $xden*$xm2 + $Q2[$i];
            }
            $res = $xm2 * ($D2 + $xm2*($xnum/$xden));
        } elsif ($y <= 12) {
#----------------------------------------------------------------------
#  4.0 < X <= 12.0
#----------------------------------------------------------------------
            $xm4 = $y - 4;
            $xden = -1;
            $xnum = 0;
            foreach my $i (0 .. 7) {
                $xnum = $xnum*$xm4 + $P4[$i];
                $xden = $xden*$xm4 + $Q4[$i];
            }
            $res = $D4 + $xm4*($xnum/$xden);
        } else {
#----------------------------------------------------------------------
#  Evaluate for argument .GE. 12.0,
#----------------------------------------------------------------------
            $res = 0;
            if ($y <= $frtbig) {
                $res = $C[6];
                my $ysq = $y * $y;
                foreach my $i (0 .. 5) {
                    $res = $res / $ysq + $C[$i];
                }
            }
            $res = $res/$y;
            $corr = log($y);
            $res = $res + $sqrtpi - 0.5*$corr;
            $res = $res + $y*($corr-1);
        }
    } else {
#----------------------------------------------------------------------
#  Return for bad arguments
#----------------------------------------------------------------------
        $res = $xinf;
    }
#----------------------------------------------------------------------
#  Final adjustments and return
#----------------------------------------------------------------------
    return $res;
}

