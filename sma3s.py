#!/usr/bin/env python3
"""
sma3s.py

This script annotates biological sequences against UniProt/UniRef-like reference
files and writes two outputs: a tab-separated annotation table and a summary report.
It keeps the original Sma3s annotator logic:

  Annotator 1: direct high-identity/high-coverage reference hit.
  Annotator 2: orthologue assignment using reciprocal best hit.
  Annotator 3: multi-hit consensus/enrichment-based annotation.

Reference input rules:
  * -d reference.dat.gz -> decompress/reuse reference.dat, then create/reuse
                            reference.fasta and reference.annot.
  * -d reference.dat    -> create/reuse reference.fasta and reference.annot.
  * -d reference.fasta  -> reference.annot must already exist next to the FASTA.

Only Python standard-library modules are required for the core workflow. MMseqs2
must be installed and available in PATH, or supplied through --mmseqs-bin. The
optional rapidgzip Python module enables parallel decompression of .dat.gz files;
without it, the script automatically falls back to the standard gzip module.

Example:
    python3 sma3s_v3.py -i proteins.faa -d reference.fasta -num_threads 16
"""
from __future__ import annotations

import argparse
import contextlib
import concurrent.futures
import dataclasses
import gzip
import math
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

SMA3S_VERSION = "3.2"
MIN_PYTHON_VERSION = (3, 8)
TYPES_BASE = ["GENENAME", "DESCRIPTION", "ENZYME", "GO", "KEYWORD", "PATHWAY"]
HUGE_VAL = 1e100
EPS = math.log(1 / HUGE_VAL)
GOEV = ":IEA"
MAX_SCORE_DEFAULT = 14

KW_Biological_process = ['Acetoin biosynthesis',
 'Acetoin catabolism',
 'Acute phase',
 'Alginate biosynthesis',
 'Alkaloid metabolism',
 'Alkylphosphonate uptake',
 'Amino-acid biosynthesis',
 'Angiogenesis',
 'Antibiotic biosynthesis',
 'Antibiotic resistance',
 'Antiviral defense',
 'Apoptosis',
 'Arginine metabolism',
 'Aromatic hydrocarbons catabolism',
 'Arsenical resistance',
 'Ascorbate biosynthesis',
 'ATP synthesis',
 'Autoinducer synthesis',
 'Autophagy',
 'Auxin biosynthesis',
 'B-cell activation',
 'Bacteriocin immunity',
 'Behavior',
 'Biological rhythms',
 'Biomineralization',
 'Biotin biosynthesis',
 'Branched-chain amino acid catabolism',
 'Cadmium resistance',
 'Calvin cycle',
 'cAMP biosynthesis',
 'Carbohydrate metabolism',
 'Carbon dioxide fixation',
 'Carnitine biosynthesis',
 'Carotenoid biosynthesis',
 'Catecholamine biosynthesis',
 'Catecholamine metabolism',
 'Cell adhesion',
 'Cell cycle',
 'Cell shape',
 'Cellulose biosynthesis',
 'cGMP biosynthesis',
 'Chemotaxis',
 'Chlorophyll biosynthesis',
 'Chromate resistance',
 'Chromosome partition',
 'Citrate utilization',
 'Cobalamin biosynthesis',
 'Coenzyme A biosynthesis',
 'Coenzyme M biosynthesis',
 'Collagen degradation',
 'Competence',
 'Conjugation',
 'Cytadherence',
 'Cytochrome c-type biogenesis',
 'Cytokinin biosynthesis',
 'Cytolysis',
 'Cytosine metabolism',
 'Deoxyribonucleotide synthesis',
 'Detoxification',
 'Diaminopimelate biosynthesis',
 'Differentiation',
 'Digestion',
 'DNA condensation',
 'DNA damage',
 'DNA excision',
 'DNA integration',
 'DNA recombination',
 'DNA replication',
 'DNA synthesis',
 'Endocytosis',
 'Enterobactin biosynthesis',
 'Erythrocyte maturation',
 'Ethylene biosynthesis',
 'Exocytosis',
 'Exopolysaccharide synthesis',
 'Fertilization',
 'Flagellar rotation',
 'Flavonoid biosynthesis',
 'Flight',
 'Flowering',
 'Folate biosynthesis',
 'Fruit ripening',
 'Galactitol metabolism',
 'Gaseous exchange',
 'Gastrulation',
 'Germination',
 'Gluconate utilization',
 'Gluconeogenesis',
 'Glutathione biosynthesis',
 'Glycerol metabolism',
 'Glycogen biosynthesis',
 'Glycolate pathway',
 'Glycolysis',
 'Glyoxylate bypass',
 'GPI-anchor biosynthesis',
 'Growth regulation',
 'Stress response',
 'Heme biosynthesis',
 'Hemolymph clotting',
 'Hemostasis',
 'Herbicide resistance',
 'Histidine metabolism',
 'Hydrogen peroxide',
 'Hypusine biosynthesis',
 'Immunity',
 'Inflammatory response',
 'Inositol biosynthesis',
 'Intron homing',
 'Iron storage',
 'Isoprene biosynthesis',
 'Karyogamy',
 'Keratinization',
 'Lactation',
 'Lactose biosynthesis',
 'Lactose metabolism',
 'Leukotriene biosynthesis',
 'Lignin biosynthesis',
 'Lignin degradation',
 'Lipid metabolism',
 'Lipopolysaccharide biosynthesis',
 'Luminescence',
 'Maltose metabolism',
 'Mandelate pathway',
 'Mast cell degranulation',
 'Meiosis',
 'Melanin biosynthesis',
 'Melatonin biosynthesis',
 'Menaquinone biosynthesis',
 'Mercuric resistance',
 'Methanogenesis',
 'Methanol utilization',
 'Methotrexate resistance',
 'Mineral balance',
 'Molybdenum cofactor biosynthesis',
 'mRNA processing',
 'Myogenesis',
 'Neurogenesis',
 'Neurotransmitter biosynthesis',
 'Neurotransmitter degradation',
 'Nitrate assimilation',
 'Nitrogen fixation',
 'Nodulation',
 'Nucleotide biosynthesis',
 'Nucleotide metabolism',
 'Nylon degradation',
 'One-carbon metabolism',
 'Pantothenate biosynthesis',
 'Pentose shunt',
 'Peptidoglycan synthesis',
 'PHA biosynthesis',
 'Phagocytosis',
 'PHB biosynthesis',
 'Phenylalanine catabolism',
 'Phenylpropanoid metabolism',
 'Pheromone response',
 'Phosphotransferase system',
 'Photorespiration',
 'Photosynthesis',
 'Phytochrome signaling pathway',
 'Plant defense',
 'Plasmid copy control',
 'Plasmid partition',
 'Plasminogen activation',
 'Polyamine biosynthesis',
 'Porphyrin biosynthesis',
 'Pregnancy',
 'Proline metabolism',
 'Protein biosynthesis',
 'Purine biosynthesis',
 'Purine metabolism',
 'Purine salvage',
 'Putrescine biosynthesis',
 'Pyridine nucleotide biosynthesis',
 'Pyridoxine biosynthesis',
 'Pyrimidine biosynthesis',
 'Queuosine biosynthesis',
 'Quinate metabolism',
 'Quorum sensing',
 'Restriction system',
 'Rhamnose metabolism',
 'Riboflavin biosynthesis',
 'Ribosome biogenesis',
 'RNA repair',
 'Viral RNA replication',
 'rRNA processing',
 'Self-incompatibility',
 'Sensory transduction',
 'Serotonin biosynthesis',
 'Spermidine biosynthesis',
 'Sporulation',
 'Starch biosynthesis',
 'Steroidogenesis',
 'Sulfate respiration',
 'Teichoic acid biosynthesis',
 'Tellurium resistance',
 'Terminal addition',
 'Tetrahydrobiopterin biosynthesis',
 'Thiamine biosynthesis',
 'Thiamine catabolism',
 'Tissue remodeling',
 'Transcription',
 'Translation regulation',
 'Transport',
 'Transposition',
 'Tricarboxylic acid cycle',
 'tRNA processing',
 'Tryptophan catabolism',
 'Tyrosine catabolism',
 'Ubiquinone biosynthesis',
 'Ubl conjugation pathway',
 'Unfolded protein response',
 'Urea cycle',
 'Virulence',
 'Nonsense-mediated mRNA decay',
 'Taxol biosynthesis',
 'Wnt signaling pathway',
 'Chlorophyll catabolism',
 'PQQ biosynthesis',
 'Chondrogenesis',
 'Osteogenesis',
 'Thyroid hormones biosynthesis',
 'Two-component regulatory system',
 'Hibernation',
 'Notch signaling pathway',
 'Interferon antiviral system evasion',
 'Auxin signaling pathway',
 'Hypersensitive response elicitation',
 'Cytokinin signaling pathway',
 'Ethylene signaling pathway',
 'Abscisic acid biosynthesis',
 'Abscisic acid signaling pathway',
 'Gibberellin signaling pathway',
 'RNA-mediated gene silencing',
 'Host-virus interaction',
 'Cell wall biogenesis/degradation',
 'Peroxisome biogenesis',
 'Cilium biogenesis/degradation',
 'Capsule biogenesis/degradation',
 'Insecticide resistance',
 'Nickel insertion',
 'Bacterial flagellum biogenesis',
 'Hearing',
 'Fimbrium biogenesis',
 'Brassinosteroid signaling pathway',
 'Cap snatching',
 'Virus entry into host cell',
 'Syncytium formation induced by viral infection',
 'Jasmonic acid signaling pathway',
 'Virus exit from host cell',
 'Viral DNA replication',
 'Viral transcription',
 'Archaeal flagellum biogenesis',
 'Necrosis',
 'Viral genome excision',
 'Viral latency']
KW_Cellular_component = ['Amyloid',
 'Antenna complex',
 'Apoplast',
 'Centromere',
 'CF(0)',
 'CF(1)',
 'Chlorosome',
 'Chromosome',
 'Chylomicron',
 'Nematocyst',
 'Copulatory plug',
 'Cuticle',
 'DNA-directed RNA polymerase',
 'Dynein',
 'Endoplasmic reticulum',
 'Exosome',
 'Fimbrium',
 'Glycosome',
 'Glyoxysome',
 'Golgi apparatus',
 'HDL',
 'Hydrogenosome',
 'Intermediate filament',
 'Keratin',
 'LDL',
 'Lysosome',
 'Membrane',
 'Membrane attack complex',
 'MHC I',
 'MHC II',
 'Microtubule',
 'Mitochondrion',
 'Nucleus',
 'Nucleomorph',
 'Lipid droplet',
 'Periplasm',
 'Peroxisome',
 'Photosystem I',
 'Photosystem II',
 'Phycobilisome',
 'Primosome',
 'Proteasome',
 'Reaction center',
 'Sarcoplasmic reticulum',
 'Signal recognition particle',
 'Signalosome',
 'Spliceosome',
 'Thick filament',
 'Thylakoid',
 'Viral occlusion body',
 'VLDL',
 'Vacuole',
 'Plastid',
 'Virion',
 'Cytoplasm',
 'Secreted',
 'Cell junction',
 'Cell projection',
 'Endosome',
 'Cytoplasmic vesicle',
 'Archaeal flagellum',
 'Bacterial flagellum',
 'Kinetochore',
 'Mitosome',
 'Host cell junction',
 'Host cell projection',
 'Host cytoplasm',
 'Host cytoplasmic vesicle',
 'Host endoplasmic reticulum',
 'Host endosome',
 'Host Golgi apparatus',
 'Host lipid droplet',
 'Host lysosome',
 'Host mitochondrion',
 'Host nucleus',
 'Host periplasm',
 'Host thylakoid',
 'Target cell cytoplasm']
KW_Developmental_stage = ['Early protein',
 'Fruiting body',
 'Heterocyst',
 'Late protein',
 'Merozoite',
 'Sporozoite',
 'Bradyzoite',
 'Tachyzoite',
 'Trophozoite']
KW_Disease = ['AIDS',
 'Albinism',
 'Allergen',
 'Alport syndrome',
 'Alzheimer disease',
 'Ectodermal dysplasia',
 'Tumor suppressor',
 'Atherosclerosis',
 'Autoimmune encephalomyelitis',
 'Autoimmune uveitis',
 'Bernard Soulier syndrome',
 'Cardiomyopathy',
 'Chronic granulomatous disease',
 'Cockayne syndrome',
 'Cone-rod dystrophy',
 'Crown gall tumor',
 'Cystinuria',
 'Deafness',
 'Dental caries',
 'Diabetes insipidus',
 'Diabetes mellitus',
 'Disease mutation',
 'Down syndrome',
 'Dwarfism',
 'Ehlers-Danlos syndrome',
 'Epidermolysis bullosa',
 'Gaucher disease',
 'Glutaricaciduria',
 'Glycogen storage disease',
 'Gangliosidosis',
 'Gout',
 'Hemophilia',
 'Hereditary hemolytic anemia',
 'Hereditary multiple exostoses',
 'Hereditary nonpolyposis colorectal cancer',
 'Hirschsprung disease',
 'Holoprosencephaly',
 'Hyperlipidemia',
 'Leber hereditary optic neuropathy',
 'Leigh syndrome',
 'Li-Fraumeni syndrome',
 'Lissencephaly',
 'Long QT syndrome',
 'Malaria',
 'Maple syrup urine disease',
 'Mucopolysaccharidosis',
 'Neurodegeneration',
 'Obesity',
 'Oncogene',
 'Phenylketonuria',
 'Neuropathy',
 'Proto-oncogene',
 'Pseudohermaphroditism',
 'Retinitis pigmentosa',
 'Rhizomelic chondrodysplasia punctata',
 'SCID',
 'Stargardt disease',
 'Stickler syndrome',
 'Systemic lupus erythematosus',
 'Thrombophilia',
 'Trypanosomiasis',
 'von Willebrand disease',
 'Whooping cough',
 'Williams-Beuren syndrome',
 'Xeroderma pigmentosum',
 'MELAS syndrome',
 'Epilepsy',
 'Cataract',
 'Congenital disorder of glycosylation',
 'Leber congenital amaurosis',
 'Primary microcephaly',
 'Parkinson disease',
 'Parkinsonism',
 'Bartter syndrome',
 'Congenital muscular dystrophy',
 'Age-related macular degeneration',
 'Fanconi anemia',
 'Progressive external ophthalmoplegia',
 'Short QT syndrome',
 'Mucolipidosis',
 'Limb-girdle muscular dystrophy',
 'Aicardi-Goutieres syndrome',
 'Familial hemophagocytic lymphohistiocytosis',
 'Congenital adrenal hyperplasia',
 'Glaucoma',
 'Kallmann syndrome',
 'Peroxisome biogenesis disorder',
 'Atrial septal defect',
 'Ichthyosis',
 'Primary hypomagnesemia',
 'Congenital hypothyroidism',
 'Congenital erythrocytosis',
 'Amelogenesis imperfecta',
 'Osteopetrosis',
 'Intrahepatic cholestasis',
 'Craniosynostosis',
 'Mental retardation',
 'Brugada syndrome',
 'Aortic aneurysm',
 'Congenital myasthenic syndrome',
 'Palmoplantar keratoderma',
 'Amyloidosis',
 'Dyskeratosis congenita',
 'Microphthalmia',
 'Congenital stationary night blindness',
 'Hypogonadotropic hypogonadism',
 'Atrial fibrillation',
 'Pontocerebellar hypoplasia',
 'Congenital generalized lipodystrophy',
 'Dystonia',
 'Diamond-Blackfan anemia',
 'Leukodystrophy',
 'Niemann-Pick disease',
 'Heterotaxy',
 'Nemaline myopathy',
 'Asthma',
 'Peters anomaly',
 'Myofibrillar myopathy',
 'Cushing syndrome',
 'Hypotrichosis',
 'Osteogenesis imperfecta',
 'Premature ovarian failure',
 'Emery-Dreifuss muscular dystrophy',
 'Hemolytic uremic syndrome',
 'Ciliopathy',
 'Schizophrenia',
 'Corneal dystrophy',
 'Dystroglycanopathy',
 'Autism spectrum disorder']
SLIM = {'GO:0000003': 'P:reproduction',
 'GO:0000228': 'C:nuclear chromosome',
 'GO:0000229': 'C:cytoplasmic chromosome',
 'GO:0000902': 'P:cell morphogenesis',
 'GO:0000988': 'F:transcription factor activity, protein binding',
 'GO:0001071': 'F:nucleic acid binding transcription factor activity',
 'GO:0002376': 'P:immune system process',
 'GO:0003013': 'P:circulatory system process',
 'GO:0003677': 'F:DNA binding',
 'GO:0003723': 'F:RNA binding',
 'GO:0003729': 'F:mRNA binding',
 'GO:0003735': 'F:structural constituent of ribosome',
 'GO:0003924': 'F:GTPase activity',
 'GO:0004386': 'F:helicase activity',
 'GO:0004518': 'F:nuclease activity',
 'GO:0004871': 'F:signal transducer activity',
 'GO:0005198': 'F:structural molecule activity',
 'GO:0005576': 'C:extracellular region',
 'GO:0005578': 'C:proteinaceous extracellular matrix',
 'GO:0005615': 'C:extracellular space',
 'GO:0005618': 'C:cell wall',
 'GO:0005622': 'C:intracellular',
 'GO:0005623': 'C:cell',
 'GO:0005634': 'C:nucleus',
 'GO:0005635': 'C:nuclear envelope',
 'GO:0005654': 'C:nucleoplasm',
 'GO:0005694': 'C:chromosome',
 'GO:0005730': 'C:nucleolus',
 'GO:0005737': 'C:cytoplasm',
 'GO:0005739': 'C:mitochondrion',
 'GO:0005764': 'C:lysosome',
 'GO:0005768': 'C:endosome',
 'GO:0005773': 'C:vacuole',
 'GO:0005777': 'C:peroxisome',
 'GO:0005783': 'C:endoplasmic reticulum',
 'GO:0005794': 'C:Golgi apparatus',
 'GO:0005811': 'C:lipid particle',
 'GO:0005815': 'C:microtubule organizing center',
 'GO:0005829': 'C:cytosol',
 'GO:0005840': 'C:ribosome',
 'GO:0005856': 'C:cytoskeleton',
 'GO:0005886': 'C:plasma membrane',
 'GO:0005929': 'C:cilium',
 'GO:0005975': 'P:carbohydrate metabolic process',
 'GO:0006091': 'P:generation of precursor metabolites and energy',
 'GO:0006259': 'P:DNA metabolic process',
 'GO:0006397': 'P:mRNA processing',
 'GO:0006399': 'P:tRNA metabolic process',
 'GO:0006412': 'P:translation',
 'GO:0006457': 'P:protein folding',
 'GO:0006461': 'P:protein complex assembly',
 'GO:0006464': 'P:cellular protein modification process',
 'GO:0006520': 'P:cellular amino acid metabolic process',
 'GO:0006605': 'P:protein targeting',
 'GO:0006629': 'P:lipid metabolic process',
 'GO:0006790': 'P:sulfur compound metabolic process',
 'GO:0006810': 'P:transport',
 'GO:0006913': 'P:nucleocytoplasmic transport',
 'GO:0006914': 'P:autophagy',
 'GO:0006950': 'P:response to stress',
 'GO:0007005': 'P:mitochondrion organization',
 'GO:0007009': 'P:plasma membrane organization',
 'GO:0007010': 'P:cytoskeleton organization',
 'GO:0007034': 'P:vacuolar transport',
 'GO:0007049': 'P:cell cycle',
 'GO:0007059': 'P:chromosome segregation',
 'GO:0007067': 'P:mitotic nuclear division',
 'GO:0007155': 'P:cell adhesion',
 'GO:0007165': 'P:signal transduction',
 'GO:0007267': 'P:cell-cell signaling',
 'GO:0007568': 'P:aging',
 'GO:0008092': 'F:cytoskeletal protein binding',
 'GO:0008134': 'F:transcription factor binding',
 'GO:0008135': 'F:translation factor activity, RNA binding',
 'GO:0008168': 'F:methyltransferase activity',
 'GO:0008219': 'P:cell death',
 'GO:0008233': 'F:peptidase activity',
 'GO:0008283': 'P:cell proliferation',
 'GO:0008289': 'F:lipid binding',
 'GO:0008565': 'F:protein transporter activity',
 'GO:0009056': 'P:catabolic process',
 'GO:0009058': 'P:biosynthetic process',
 'GO:0009536': 'C:plastid',
 'GO:0009579': 'C:thylakoid',
 'GO:0009790': 'P:embryo development',
 'GO:0015979': 'P:photosynthesis',
 'GO:0016023': 'C:cytoplasmic, membrane-bounded vesicle',
 'GO:0016192': 'P:vesicle-mediated transport',
 'GO:0016301': 'F:kinase activity',
 'GO:0016491': 'F:oxidoreductase activity',
 'GO:0016746': 'F:transferase activity, transferring acyl groups',
 'GO:0016757': 'F:transferase activity, transferring glycosyl groups',
 'GO:0016765': 'F:transferase activity, transferring alkyl or aryl (other than methyl) groups',
 'GO:0016779': 'F:nucleotidyltransferase activity',
 'GO:0016791': 'F:phosphatase activity',
 'GO:0016798': 'F:hydrolase activity, acting on glycosyl bonds',
 'GO:0016810': 'F:hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds',
 'GO:0016829': 'F:lyase activity',
 'GO:0016853': 'F:isomerase activity',
 'GO:0016874': 'F:ligase activity',
 'GO:0016887': 'F:ATPase activity',
 'GO:0019748': 'P:secondary metabolic process',
 'GO:0019843': 'F:rRNA binding',
 'GO:0019899': 'F:enzyme binding',
 'GO:0021700': 'P:developmental maturation',
 'GO:0022607': 'P:cellular component assembly',
 'GO:0022618': 'P:ribonucleoprotein complex assembly',
 'GO:0022857': 'F:transmembrane transporter activity',
 'GO:0030154': 'P:cell differentiation',
 'GO:0030198': 'P:extracellular matrix organization',
 'GO:0030234': 'F:enzyme regulator activity',
 'GO:0030312': 'C:external encapsulating structure',
 'GO:0030533': 'F:triplet codon-amino acid adaptor activity',
 'GO:0030555': 'F:RNA modification guide activity',
 'GO:0030674': 'F:protein binding, bridging',
 'GO:0030705': 'P:cytoskeleton-dependent intracellular transport',
 'GO:0032182': 'F:ubiquitin-like protein binding',
 'GO:0032196': 'P:transposition',
 'GO:0034330': 'P:cell junction organization',
 'GO:0034641': 'P:cellular nitrogen compound metabolic process',
 'GO:0034655': 'P:nucleobase-containing compound catabolic process',
 'GO:0040007': 'P:growth',
 'GO:0040011': 'P:locomotion',
 'GO:0042254': 'P:ribosome biogenesis',
 'GO:0042393': 'F:histone binding',
 'GO:0042592': 'P:homeostatic process',
 'GO:0043167': 'F:ion binding',
 'GO:0043226': 'C:organelle',
 'GO:0043234': 'C:protein complex',
 'GO:0043473': 'P:pigmentation',
 'GO:0044281': 'P:small molecule metabolic process',
 'GO:0044403': 'P:symbiosis, encompassing mutualism through parasitism',
 'GO:0048646': 'P:anatomical structure formation involved in morphogenesis',
 'GO:0048856': 'P:anatomical structure development',
 'GO:0048870': 'P:cell motility',
 'GO:0050877': 'P:neurological system process',
 'GO:0051082': 'F:unfolded protein binding',
 'GO:0051186': 'P:cofactor metabolic process',
 'GO:0051276': 'P:chromosome organization',
 'GO:0051301': 'P:cell division',
 'GO:0051604': 'P:protein maturation',
 'GO:0055085': 'P:transmembrane transport',
 'GO:0061024': 'P:membrane organization',
 'GO:0065003': 'P:macromolecular complex assembly',
 'GO:0071554': 'P:cell wall organization or biogenesis',
 'GO:0071941': 'P:nitrogen cycle metabolic process'}
ECO = ['0000501',
 '0000203',
 '0000209',
 '0000348',
 '0000350',
 '0000347',
 '0000331',
 '0000213',
 '0000246',
 '0000254',
 '0000313',
 '0000259',
 '0000256',
 '0000258',
 '0000210',
 '0000211',
 '0000248',
 '0000265',
 '0000249',
 '0000251',
 '0000261',
 '0000263',
 '0000332']

KW_CATEGORIES = {
    "Biological_process": set(KW_Biological_process),
    "Cellular_component": set(KW_Cellular_component),
    "Developmental_stage": set(KW_Developmental_stage),
    "Disease": set(KW_Disease),
}
ECO_SET = set(ECO)


@dataclasses.dataclass
class Hit:
    query: str
    target: str
    pident: float
    qcov: float
    alnlen: int
    evalue: float
    tcov: float = 0.0
    qlen: int = 0
    tlen: int = 0
    rank: int = 0


@dataclasses.dataclass
class Config:
    annotator: str
    query_file: Path
    dat_file: Path
    fasta_file: Path
    compressed_dat_file: Optional[Path]
    mmseqs_file: Path
    annot_file: Path
    rost: float = 20.0
    pv: float = 0.1
    training: bool = False
    noempty: bool = False
    nucl: bool = False
    extend_go: bool = False
    quality: bool = False
    nopred: bool = False
    source: bool = False
    uniref: bool = True
    goslim: bool = False
    cpus: int = 1
    annotation_workers: int = 1
    cluster_threads: int = 1
    decompression_threads: int = 1
    id_uniprot: float = 90.0
    id_orthologue: float = 75.0
    cov_uniprot: float = 90.0
    cov_orthologue: float = 80.0
    user_changed_thresholds: bool = False
    max_score: int = MAX_SCORE_DEFAULT
    low_complexity_filter: bool = False
    max_seqs: int = 250
    sensitivity: float = 7.5
    split_memory_limit: Optional[str] = None
    db_load_mode: Optional[int] = None
    tmpdir: Path = Path(".")
    force_search: bool = False
    keep_tmp: bool = False
    mmseqs_bin: str = "mmseqs"
    genus_filter: Optional[str] = None
    exclude_genus_filter: Optional[str] = None
    family_filter: Optional[str] = None
    exclude_family_filter: Optional[str] = None
    order_filter: Optional[str] = None
    exclude_order_filter: Optional[str] = None
    report_hit_metrics: bool = False
    force_clean_outputs: bool = False
    force_clean_db: bool = False
    clean_only: bool = False
    reciprocal_fasta_mode: str = "auto"
    reciprocal_stream_threshold: int = 50000

    @property
    def types(self) -> List[str]:
        types = list(TYPES_BASE)
        if self.goslim:
            types.append("GOSLIM")
        return types


def die(msg: str, code: int = 2) -> None:
    raise SystemExit(f"\nError: {msg}\n")


def run_cmd(cmd: Sequence[str], quiet: bool = False) -> None:
    if not quiet:
        print("[cmd] " + " ".join(map(str, cmd)), flush=True)
    try:
        subprocess.run(list(map(str, cmd)), check=True)
    except FileNotFoundError:
        die(f"Executable not found: {cmd[0]}")
    except subprocess.CalledProcessError as exc:
        die(f"Command failed with exit code {exc.returncode}: {' '.join(map(str, cmd))}")


def mmseqs_version(mmseqs_bin: str) -> Tuple[Optional[str], Optional[str]]:
    """Return (resolved executable, version text), or (None, None) on failure."""
    resolved = shutil.which(mmseqs_bin)
    if resolved is None:
        return None, None
    try:
        proc = subprocess.run(
            [resolved, "version"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return resolved, None
    if proc.returncode != 0:
        return resolved, None
    version = (proc.stdout or "").strip().splitlines()
    return resolved, version[0] if version else "version command succeeded"


def check_mmseqs(mmseqs_bin: str) -> None:
    resolved, version = mmseqs_version(mmseqs_bin)
    if resolved is None:
        die("MMseqs2 is not installed or is not in PATH. Try: conda install -c conda-forge -c bioconda mmseqs2")
    if version is None:
        die(f"MMseqs2 was found at {resolved}, but '{mmseqs_bin} version' failed")


def check_installation(mmseqs_bin: str, tmpdir: Path) -> int:
    """Check the runtime dependencies and basic filesystem/SQLite functionality."""
    print(f"Sma3s v{SMA3S_VERSION} installation check")
    print("=" * 38)
    failures: List[str] = []
    warnings: List[str] = []

    py_ok = sys.version_info >= MIN_PYTHON_VERSION
    py_version = ".".join(map(str, sys.version_info[:3]))
    print(f"[{'OK' if py_ok else 'FAIL'}] Python {py_version} "
          f"(required >= {MIN_PYTHON_VERSION[0]}.{MIN_PYTHON_VERSION[1]})")
    if not py_ok:
        failures.append("unsupported Python version")

    try:
        tmpdir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=tmpdir) as td:
            test_dir = Path(td)
            test_file = test_dir / "write_test.txt"
            test_file.write_text("sma3s\n", encoding="utf-8")
            gzip_test = gzip.decompress(gzip.compress(b"sma3s-gzip-check"))
            if gzip_test != b"sma3s-gzip-check":
                raise RuntimeError("unexpected gzip round-trip result")
            db_file = test_dir / "sqlite_test.sqlite"
            con = sqlite3.connect(db_file)
            con.execute("CREATE TABLE test(value INTEGER)")
            con.execute("INSERT INTO test VALUES (1)")
            value = con.execute("SELECT value FROM test").fetchone()[0]
            con.close()
            if value != 1:
                raise RuntimeError("unexpected SQLite test result")
        print(f"[OK] Writable temporary directory and SQLite: {tmpdir.resolve()}")
    except Exception as exc:
        print(f"[FAIL] Temporary directory/SQLite test: {exc}")
        failures.append("temporary directory or SQLite is not usable")

    resolved, version = mmseqs_version(mmseqs_bin)
    if resolved is None:
        print(f"[FAIL] MMseqs2 executable not found: {mmseqs_bin}")
        failures.append("MMseqs2 is missing")
    elif version is None:
        print(f"[FAIL] MMseqs2 found at {resolved}, but its version command failed")
        failures.append("MMseqs2 cannot be executed")
    else:
        print(f"[OK] MMseqs2: {resolved} ({version})")

    try:
        import rapidgzip  # type: ignore
        with tempfile.TemporaryDirectory(dir=tmpdir) as td:
            gz_test_path = Path(td) / "rapidgzip_test.dat.gz"
            with gzip.open(gz_test_path, "wb") as out:
                out.write(b"ID   SMA3S_CHECK\n//\n")
            with rapidgzip.open(str(gz_test_path), mode="rb", parallelization=2) as inp:
                if inp.read() != b"ID   SMA3S_CHECK\n//\n":
                    raise RuntimeError("unexpected rapidgzip test result")
        rapidgzip_version = getattr(rapidgzip, "__version__", "installed")
        print(f"[OK] rapidgzip: {rapidgzip_version} (parallel .dat.gz decompression enabled)")
    except ImportError:
        msg = ("rapidgzip is not installed; .dat.gz files will still work through "
               "the standard gzip module, but decompression will be single-threaded")
        print(f"[WARN] {msg}")
        warnings.append(msg)
    except Exception as exc:
        msg = f"rapidgzip is installed but its decompression test failed: {exc}"
        print(f"[WARN] {msg}")
        warnings.append(msg)

    if failures:
        print("\nInstallation check: FAILED")
        for item in failures:
            print(f"  - {item}")
        return 1

    print("\nInstallation check: OK")
    if warnings:
        print("The core workflow is functional, with the warning shown above.")
    return 0


def decompressed_dat_path(compressed: Path) -> Path:
    """Return the .dat path corresponding to a .dat.gz input."""
    if compressed.suffix.lower() != ".gz":
        die(f"Compressed database must end in .gz: {compressed}")
    dat = compressed.with_suffix("")
    if dat.suffix.lower() != ".dat":
        die(f"Compressed UniProt input must end in .dat.gz: {compressed}")
    return dat


def _copy_decompressed_stream(source, destination: Path) -> None:
    """Copy a decompressed binary stream to an atomic temporary output."""
    with destination.open("wb") as out:
        shutil.copyfileobj(source, out, length=8 * 1024 * 1024)


def ensure_decompressed_dat(cfg: Config) -> Path:
    """Create or reuse the .dat file derived from a .dat.gz reference input.

    rapidgzip is used when available and performs parallel decompression with the
    requested number of threads. The standard-library gzip module is a compatible
    single-threaded fallback. The final .dat is installed atomically.
    """
    compressed = cfg.compressed_dat_file
    if compressed is None:
        return cfg.dat_file
    if not compressed.exists():
        if cfg.dat_file.exists():
            print(f"Compressed database is absent; reusing existing decompressed file: {cfg.dat_file}", flush=True)
            return cfg.dat_file
        die(f"Compressed database does not exist: {compressed}")

    if (cfg.dat_file.exists() and cfg.dat_file.stat().st_size > 0
            and cfg.dat_file.stat().st_mtime_ns >= compressed.stat().st_mtime_ns):
        print(f"Reusing decompressed UniProt database: {cfg.dat_file}", flush=True)
        return cfg.dat_file

    cfg.dat_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_dat = cfg.dat_file.with_name(f"{cfg.dat_file.name}.tmp.{os.getpid()}")
    with contextlib.suppress(FileNotFoundError):
        tmp_dat.unlink()

    method = "standard gzip (single thread)"
    try:
        rapidgzip_error: Optional[Exception] = None
        try:
            import rapidgzip  # type: ignore
        except ImportError:
            rapidgzip = None

        if rapidgzip is not None:
            threads = max(1, int(cfg.decompression_threads))
            method = f"rapidgzip ({threads} threads)"
            print(f"Decompressing {compressed} -> {cfg.dat_file} with {method}", flush=True)
            try:
                with rapidgzip.open(str(compressed), mode="rb", parallelization=threads) as inp:
                    _copy_decompressed_stream(inp, tmp_dat)
            except Exception as exc:
                rapidgzip_error = exc
                with contextlib.suppress(FileNotFoundError):
                    tmp_dat.unlink()

        if rapidgzip is None or rapidgzip_error is not None:
            if rapidgzip_error is not None:
                print(
                    f"rapidgzip failed ({rapidgzip_error}); retrying with the standard gzip module.",
                    flush=True,
                )
            else:
                print(
                    f"Decompressing {compressed} -> {cfg.dat_file} with {method}. "
                    "Install rapidgzip to enable parallel decompression.",
                    flush=True,
                )
            with gzip.open(compressed, mode="rb") as inp:
                _copy_decompressed_stream(inp, tmp_dat)

        if not tmp_dat.exists() or tmp_dat.stat().st_size == 0:
            raise RuntimeError("the decompressed .dat file is empty")
        with tmp_dat.open("rb") as check:
            probe = check.read(1024 * 1024)
        if probe.lstrip().startswith(b">"):
            raise RuntimeError("the decompressed input looks like FASTA, not a UniProt .dat file")
        if not (probe.startswith(b"ID   ") or b"\nID   " in probe):
            raise RuntimeError("no UniProt ID entry was found near the start of the decompressed file")
        os.replace(tmp_dat, cfg.dat_file)
        compressed_mtime = compressed.stat().st_mtime_ns
        os.utime(cfg.dat_file, ns=(time.time_ns(), compressed_mtime))
    except Exception as exc:
        with contextlib.suppress(FileNotFoundError):
            tmp_dat.unlink()
        die(f"Could not decompress {compressed}: {exc}")

    print(f"Decompressed UniProt database created: {cfg.dat_file}", flush=True)
    return cfg.dat_file


def pvalue_from_evalue(evalue: float) -> float:
    # Same conversion used by the Perl code.
    return 1.0 - math.exp(-evalue)


def calculate_qc_subject(len_query: int, qcov_query_percent: float, length_subject: int) -> float:
    if length_subject <= 0:
        return 0.0
    return len_query * qcov_query_percent / length_subject


def calculate_rost(rost: float, aln_len: int) -> float:
    if aln_len <= 0:
        return float("inf")
    return rost + (480.0 * (aln_len ** (-0.32 * (1.0 + math.exp(-aln_len / 1000.0)))))


def filter_cs(s: str) -> str:
    return re.escape(s)


def calculate_gnscore(g: str) -> int:
    if not g:
        return 0
    g = re.sub(r"\{.*", "", g)
    gn_l = len(g)
    score = 0
    if not re.search(r"[._-]", g):
        score += 4
    if g[:1].islower():
        score += 2
    if gn_l <= 4:
        score += 4
    elif gn_l <= 5:
        score += 3
    elif gn_l <= 6:
        score += 2
    elif gn_l <= 8:
        score += 1
    return score


@lru_cache(maxsize=2_000_000)
def lchoose(x: int, y: int) -> float:
    # Perl function name is lchoose(x, y) but mathematically computes log(C(y, x)).
    if x < 0 or y < 0 or x > y:
        raise ValueError(f"lchoose invalid args: {x}, {y}")
    return math.lgamma(y + 1) - math.lgamma(x + 1) - math.lgamma(y - x + 1)


def add_log(x: float, y: float) -> float:
    if x == -HUGE_VAL:
        return y
    if y == -HUGE_VAL:
        return x
    if x >= y:
        y -= x
    else:
        x, y = y, x - y
    if y < EPS:
        y = EPS
    return x + math.log1p(math.exp(y))


def compute_hyper_p_value(k: int, n: int, K: int, N: int) -> float:
    # Direct port of Perl ComputeHyperPValue.
    if N == 0:
        return 0.0
    p_val = -HUGE_VAL
    while k >= 0 and k <= n and k <= K and (n - k) <= (N - K):
        x = lchoose(k, K) + lchoose(n - k, N - K) - lchoose(n, N)
        p_val = add_log(p_val, x)
        k += 1
    if p_val > 0:
        p_val = 0.0
    return math.exp(p_val)


def fasta_iter(path: Path) -> Iterator[Tuple[str, str]]:
    seq_id = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(chunks)
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(re.sub(r"\s+", "", line))
        if seq_id is not None:
            yield seq_id, "".join(chunks)


def normalize_query_fasta(input_fasta: Path, output_fasta: Path, nucl: bool) -> Tuple[Dict[str, str], Dict[str, int], int]:
    mapping: Dict[str, str] = {}
    lengths: Dict[str, int] = {}
    n = 0
    nucl_probe = ""
    with input_fasta.open("r", encoding="utf-8", errors="replace") as inp, output_fasta.open("w", encoding="utf-8") as out:
        current = None
        for raw in inp:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                n += 1
                sid = f"s{n}"
                mapping[sid] = line[1:]
                lengths[sid] = 0
                current = sid
                out.write(f">{sid}\n")
            elif current:
                seq_line = re.sub(r"\s+", "", line)
                if n <= 5 and not nucl:
                    nucl_probe += seq_line
                    if n == 5 and re.fullmatch(r"[ACGTUNacgtun]+", nucl_probe or ""):
                        die("The first sequences look nucleotide-like, but -nucl was not used")
                lengths[current] += len(seq_line)
                out.write(seq_line + "\n")
    return mapping, lengths, n


def output_prefix(query_file: Path, dat_file: Path, uniref: bool) -> Path:
    annot_prefix = Path(re.sub(r"\.[A-Za-z0-9_-]+$", "_", str(query_file)))
    dat_brief = dat_file.name
    if dat_brief.endswith(".dat"):
        dat_brief = dat_brief[:-4]
    if not uniref:
        suffix = "".join(part.capitalize()[:3] for part in dat_brief.split("_"))
    else:
        suffix = re.sub(r"\..+", "", dat_brief)
    return Path(str(annot_prefix) + suffix)


def safe_name_suffix(value: str) -> str:
    """Return a filesystem/filename-safe suffix from a user-provided taxon/genus name."""
    suffix = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    suffix = suffix.strip("._-")
    return suffix or "taxon"


def tax_filter_suffix(
    genus: Optional[str] = None,
    exclude_genus: Optional[str] = None,
    family: Optional[str] = None,
    exclude_family: Optional[str] = None,
    order: Optional[str] = None,
    exclude_order: Optional[str] = None,
    *,
    legacy_single_genus: bool = True,
) -> str:
    """Build a deterministic suffix for derived filtered FASTA/.annot files."""
    if legacy_single_genus and genus and not any([exclude_genus, family, exclude_family, order, exclude_order]):
        return f".{safe_name_suffix(genus)}"
    if legacy_single_genus and exclude_genus and not any([genus, family, exclude_family, order, exclude_order]):
        return f".without_{safe_name_suffix(exclude_genus)}"
    parts: List[str] = []
    if order:
        parts.append(f"order_{safe_name_suffix(order)}")
    if family:
        parts.append(f"family_{safe_name_suffix(family)}")
    if genus:
        parts.append(f"genus_{safe_name_suffix(genus)}")
    if exclude_order:
        parts.append(f"without_order_{safe_name_suffix(exclude_order)}")
    if exclude_family:
        parts.append(f"without_family_{safe_name_suffix(exclude_family)}")
    if exclude_genus:
        parts.append(f"without_genus_{safe_name_suffix(exclude_genus)}")
    return "." + ".".join(parts) if parts else ""


def has_tax_filter(cfg_or_args) -> bool:
    return any(getattr(cfg_or_args, name, None) for name in (
        "genus_filter", "exclude_genus_filter", "family_filter", "exclude_family_filter", "order_filter", "exclude_order_filter"
    ))


def has_tax_exclusion(cfg_or_args) -> bool:
    return any(getattr(cfg_or_args, name, None) for name in (
        "exclude_genus_filter", "exclude_family_filter", "exclude_order_filter"
    ))


def uniprot_oc_terms(entry: str) -> List[str]:
    """Return normalized UniProt OC lineage nodes from one text entry."""
    oc_text = " ".join(re.findall(r"^OC\s+(.+)$", entry, flags=re.M))
    return [x.strip().lower() for x in re.split(r"[;.]", oc_text) if x.strip()]


def uniprot_entry_matches_taxon(entry: str, taxon: Optional[str], *, allow_genus_os_fallback: bool = False) -> bool:
    """Return True when a UniProt text entry contains an exact OC lineage node.

    For genus filters, UniProt OC normally contains the genus as an exact node, e.g.
    "OC   ...; Enterobacteriaceae; Salmonella.". As a safety fallback for genus only,
    the first word of OS is accepted when OC is absent.
    """
    if not taxon:
        return True
    wanted = taxon.strip().lower()
    if not wanted:
        return True

    oc_terms = uniprot_oc_terms(entry)
    if oc_terms:
        return wanted in oc_terms

    if allow_genus_os_fallback:
        os_text = " ".join(re.findall(r"^OS\s+(.+)$", entry, flags=re.M)).strip()
        if os_text:
            first = re.split(r"\s+", os_text, maxsplit=1)[0].strip(" .;:,()[]{}").lower()
            return first == wanted
    return False


def uniprot_entry_matches_genus(entry: str, genus: Optional[str]) -> bool:
    return uniprot_entry_matches_taxon(entry, genus, allow_genus_os_fallback=True)


def parse_uniprot_entry(
    entry: str,
    nopred: bool,
    genus_filter: Optional[str] = None,
    exclude_genus_filter: Optional[str] = None,
    family_filter: Optional[str] = None,
    exclude_family_filter: Optional[str] = None,
    order_filter: Optional[str] = None,
    exclude_order_filter: Optional[str] = None,
) -> Optional[Tuple[str, str]]:
    """Parse one UniProt text entry and return (fasta_record, annot_line).

    Returns None when the entry has no usable sequence/annotation or is filtered out.
    Taxonomic filters are exact UniProt OC lineage-node matches. Inclusion filters are
    combined with AND; exclusion filters remove an entry if any excluded node matches.
    This function is top-level so it can be used by ProcessPoolExecutor.
    """
    if genus_filter and not uniprot_entry_matches_taxon(entry, genus_filter, allow_genus_os_fallback=True):
        return None
    if family_filter and not uniprot_entry_matches_taxon(entry, family_filter):
        return None
    if order_filter and not uniprot_entry_matches_taxon(entry, order_filter):
        return None
    if exclude_genus_filter and uniprot_entry_matches_taxon(entry, exclude_genus_filter, allow_genus_os_fallback=True):
        return None
    if exclude_family_filter and uniprot_entry_matches_taxon(entry, exclude_family_filter):
        return None
    if exclude_order_filter and uniprot_entry_matches_taxon(entry, exclude_order_filter):
        return None

    n_reviewed = {"Reviewed": 1, "Unreviewed": 0}
    stop_de = {"Uncharacterized protein", "Putative uncharacterized protein"}
    stop_kw = {
        "3D-structure", "Allosteric enzyme", "Complete proteome", "Genetically modified food", "Hybridoma",
        "Multifunctional enzyme", "Pharmaceutical", "ERV", "Direct protein sequencing", "Extinct organism protein",
        "Proteomics identification", "Reference proteome", "Alternative initiation", "Alternative splicing",
        "Polymorphism", "Alternative promoter usage",
    }

    e = re.sub(r"\|[^\},\|]+([\},])", r"\1", entry, flags=re.S)
    e = e.replace(" {ECO:", "{ECO:")
    e = re.sub(r"ECO:(\d+)\|\w+:\d+", r"ECO:\1||", e, flags=re.S)
    e = re.sub(r"ECO:(\d+)[,;] ?", r"ECO:\1||", e, flags=re.S)

    m = re.search(r"^ID   (.+?)\s+(\w+);", e, flags=re.M)
    if not m:
        return None
    uid, reviewed = m.group(1), m.group(2)
    db_score = n_reviewed.get(reviewed, 0)

    seq_lines = re.findall(r"^     (.+)$", e, flags=re.M)
    seq = re.sub(r"[\s\n]", "", "\n".join(seq_lines))
    if not seq:
        return None

    out_fields: List[str] = []

    # GeneName
    gn_lines = re.findall(r"^GN   (.+)$", e, flags=re.M)
    gn = "".join(gn_lines)
    gn = re.sub(r"\w+=", "", gn).rstrip(";")
    gns = [x for x in re.split(r"[;,] *", gn) if x]
    out_fields.append(";".join(gns))
    gn_score = calculate_gnscore(gns[0]) if gns else 0

    # Description + EC
    de_lines = re.findall(r"^DE   (.+)$", e, flags=re.M)
    de = "".join(de_lines)
    de = re.sub(r";Flags:.+", "", de)
    de = re.sub(r"\w+: ", "", de)
    ecs = re.findall(r"EC=([0-9\.\-]+\{?E?C?O?:?[^\};]*\}?);?", de)
    de = re.sub(r"EC=[0-9\.]+\{?E?C?O?:?[^\};]*\}?;?", "", de)
    de = re.sub(r"\w+=", "", de)
    de_out: List[str] = []
    de_score = 0
    for item in re.split(r"; *", de):
        item = item.strip()
        if not item:
            continue
        clean = re.sub(r"\{.*\}", "", item)
        if clean in stop_de:
            continue
        de_score = 4
        de_out.append(item)
    out_fields.append(";".join(de_out))
    out_fields.append(";".join(ecs))

    # GO terms
    gos: List[str] = []
    for go in re.findall(r"^DR   GO; (GO:\d+; [PCF]:.+; \w+):", e, flags=re.M):
        go2 = go.replace("; ", "{", 1).replace("; ", ":", 1) + "}"
        gos.append(go2)
    out_fields.append(";".join(gos))

    # Keywords
    kw_out: List[str] = []
    mkw = re.search(r"^KW   ([^\.]+)\.", e, flags=re.M | re.S)
    if mkw:
        kw = mkw.group(1).replace("\n", " ").replace("KW   ", "")
        for item in re.split(r"; ?", kw):
            item = item.strip()
            if not item:
                continue
            clean = re.sub(r"\{.*\}", "", item)
            if any(clean.startswith(x) for x in stop_kw):
                continue
            kw_out.append(item)
    out_fields.append(";".join(kw_out))

    # Pathways
    pathways: List[str] = []
    for pw in re.findall(r"^CC   -!- PATHWAY: (.+?)\nCC   -", e, flags=re.M | re.S):
        pw = pw.replace("\nCC      ", "")
        eco_pw = ""
        meco = re.search(r"(\{ECO:.+\})", pw)
        if meco:
            eco_pw = meco.group(1)
            pw = pw.replace(eco_pw, "")
        pw = pw.split(";")[0].replace(".", "")
        pathways.append(pw + eco_pw)
    out_fields.append(";".join(sorted(set(pathways))))

    pe = 0
    mpe = re.search(r"^PE   (\d):", e, flags=re.M)
    if nopred and mpe:
        pe = int(mpe.group(1))

    total_score = gn_score + de_score
    if not ((gos or kw_out) and pe != 4):
        return None

    fasta_record = f">{uid}\n{seq}\n"
    annot_line = f"{uid}\t{db_score},{total_score}\t" + "\t".join(out_fields) + "\n"
    return fasta_record, annot_line


def iter_uniprot_entries(dat: Path) -> Iterator[str]:
    entry_lines: List[str] = []
    with dat.open("r", encoding="utf-8", errors="replace") as inp:
        for line in inp:
            entry_lines.append(line)
            if line.startswith("//"):
                yield "".join(entry_lines)
                entry_lines = []
        if entry_lines:
            # In case the last entry does not end with //.
            yield "".join(entry_lines)


def batched_iterator(iterator: Iterable[str], batch_size: int) -> Iterator[List[str]]:
    batch: List[str] = []
    for item in iterator:
        batch.append(item)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def create_fasta_annot_from_uniprot(
    dat: Path,
    fasta: Path,
    annot: Path,
    nopred: bool,
    workers: int = 1,
    genus_filter: Optional[str] = None,
    exclude_genus_filter: Optional[str] = None,
    family_filter: Optional[str] = None,
    exclude_family_filter: Optional[str] = None,
    order_filter: Optional[str] = None,
    exclude_order_filter: Optional[str] = None,
) -> None:
    with dat.open("r", encoding="utf-8", errors="replace") as check:
        for line in check:
            if line.strip():
                if line.startswith(">"):
                    die(f"{dat} looks like FASTA. Use a UniProt .dat file or a FASTA file with its companion .annot file")
                break

    workers = max(1, int(workers))
    filter_msgs = []
    if order_filter: filter_msgs.append(f"order={order_filter}")
    if family_filter: filter_msgs.append(f"family={family_filter}")
    if genus_filter: filter_msgs.append(f"genus={genus_filter}")
    if exclude_order_filter: filter_msgs.append(f"exclude order={exclude_order_filter}")
    if exclude_family_filter: filter_msgs.append(f"exclude family={exclude_family_filter}")
    if exclude_genus_filter: filter_msgs.append(f"exclude genus={exclude_genus_filter}")
    filter_msg = "; filters: " + ", ".join(filter_msgs) if filter_msgs else ""
    print(f"Creating FASTA and .annot from {dat} with {workers} worker(s){filter_msg}", flush=True)

    # Write in the main process to avoid write collisions and preserve .dat entry order.
    # Process bounded batches to avoid loading the entire .dat file into memory.
    batch_size = 512
    n_entries = 0
    n_kept = 0
    with fasta.open("w", encoding="utf-8") as fa, annot.open("w", encoding="utf-8") as an:
        if workers <= 1:
            for entry in iter_uniprot_entries(dat):
                n_entries += 1
                parsed = parse_uniprot_entry(entry, nopred, genus_filter, exclude_genus_filter, family_filter, exclude_family_filter, order_filter, exclude_order_filter)
                if parsed is None:
                    continue
                fasta_record, annot_line = parsed
                fa.write(fasta_record)
                an.write(annot_line)
                n_kept += 1
                if n_entries % 100000 == 0:
                    print(f"  UniProt entries processed: {n_entries:,}; kept: {n_kept:,}", flush=True)
        else:
            with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                for batch in batched_iterator(iter_uniprot_entries(dat), batch_size):
                    for parsed in executor.map(
                            parse_uniprot_entry,
                            batch,
                            [nopred] * len(batch),
                            [genus_filter] * len(batch),
                            [exclude_genus_filter] * len(batch),
                            [family_filter] * len(batch),
                            [exclude_family_filter] * len(batch),
                            [order_filter] * len(batch),
                            [exclude_order_filter] * len(batch),
                            chunksize=32,
                        ):
                        n_entries += 1
                        if parsed is None:
                            continue
                        fasta_record, annot_line = parsed
                        fa.write(fasta_record)
                        an.write(annot_line)
                        n_kept += 1
                    if n_entries % 100000 < batch_size:
                        print(f"  UniProt entries processed: {n_entries:,}; kept: {n_kept:,}", flush=True)

    print(f"FASTA/.annot created: {n_kept:,} proteins kept from {n_entries:,} UniProt entries", flush=True)


def fasta_index_db_path(fasta: Path) -> Path:
    return fasta.with_suffix(fasta.suffix + ".sma3s_fasta_index.sqlite")


def _sqlite_table_exists(con: sqlite3.Connection, table: str) -> bool:
    row = con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
        (table,),
    ).fetchone()
    return row is not None


def _fasta_index_is_current(db_path: Path, fasta: Path) -> bool:
    """Return True when the cached FASTA offset index matches the FASTA file.

    New indexes store FASTA size and nanosecond mtime in a metadata table.
    Older indexes did not have metadata; for those, fall back to the simpler
    check used previously: index mtime >= FASTA mtime.
    """
    if not db_path.exists():
        return False
    try:
        fasta_stat = fasta.stat()
        con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        try:
            if not _sqlite_table_exists(con, "seqs"):
                return False
            if _sqlite_table_exists(con, "metadata"):
                rows = dict(con.execute("SELECT key, value FROM metadata"))
                return (
                    rows.get("fasta_size") == str(fasta_stat.st_size)
                    and rows.get("fasta_mtime_ns") == str(fasta_stat.st_mtime_ns)
                )
        finally:
            con.close()
        # Backwards-compatible fallback for indexes made by older script versions.
        return db_path.stat().st_mtime >= fasta_stat.st_mtime
    except Exception:
        return False


def ensure_fasta_index(fasta: Path) -> Path:
    db_path = fasta_index_db_path(fasta)
    if _fasta_index_is_current(db_path, fasta):
        print(f"Reusing FASTA offset index: {db_path}", flush=True)
        return db_path

    print(f"Indexing FASTA by byte offsets: {fasta}", flush=True)
    if db_path.exists():
        print(f"  Index is missing, incomplete or older than the FASTA; it will be rebuilt: {db_path}", flush=True)

    tmp_db_path = db_path.with_name(f"{db_path.name}.tmp.{os.getpid()}")
    if tmp_db_path.exists():
        tmp_db_path.unlink()

    con = sqlite3.connect(tmp_db_path)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("CREATE TABLE seqs(id TEXT PRIMARY KEY, len INTEGER NOT NULL, offset INTEGER NOT NULL)")
    with fasta.open("rb") as fh:
        current_id = None
        seq_start = None
        seq_len = 0
        batch = []
        while True:
            line = fh.readline()
            if not line:
                break
            if line.startswith(b">"):
                if current_id is not None:
                    batch.append((current_id, seq_len, seq_start))
                    if len(batch) >= 10000:
                        con.executemany("INSERT OR REPLACE INTO seqs VALUES (?,?,?)", batch)
                        con.commit()
                        batch.clear()
                header = line[1:].decode("utf-8", "replace").strip()
                current_id = header.split()[0]
                seq_start = fh.tell()
                seq_len = 0
            else:
                seq_len += len(re.sub(rb"\s+", b"", line))
        if current_id is not None:
            batch.append((current_id, seq_len, seq_start))
        if batch:
            con.executemany("INSERT OR REPLACE INTO seqs VALUES (?,?,?)", batch)
            con.commit()
    con.execute("CREATE INDEX IF NOT EXISTS idx_len ON seqs(len)")

    fasta_stat = fasta.stat()
    con.execute("CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    con.executemany(
        "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
        [
            ("fasta_path", str(fasta.resolve())),
            ("fasta_size", str(fasta_stat.st_size)),
            ("fasta_mtime_ns", str(fasta_stat.st_mtime_ns)),
            ("created_unix", f"{time.time():.6f}"),
            ("script", "sma3s_mmseqs"),
        ],
    )
    con.commit()
    con.close()

    os.replace(tmp_db_path, db_path)
    print(f"FASTA index created: {db_path}", flush=True)
    return db_path


class FastaStore:
    def __init__(self, fasta: Path):
        self.fasta = fasta
        self.db_path = ensure_fasta_index(fasta)
        self.con = sqlite3.connect(self.db_path)
        self.fh = fasta.open("rb")

    def close(self) -> None:
        with contextlib.suppress(Exception):
            self.fh.close()
        with contextlib.suppress(Exception):
            self.con.close()

    def length(self, seq_id: str) -> int:
        row = self.con.execute("SELECT len FROM seqs WHERE id=?", (seq_id,)).fetchone()
        return int(row[0]) if row else 0

    def lengths(self, ids: Iterable[str]) -> Dict[str, int]:
        return {sid: self.length(sid) for sid in ids}

    def sequence(self, seq_id: str) -> str:
        row = self.con.execute("SELECT offset FROM seqs WHERE id=?", (seq_id,)).fetchone()
        if not row:
            return ""
        self.fh.seek(int(row[0]))
        chunks: List[bytes] = []
        while True:
            line = self.fh.readline()
            if not line or line.startswith(b">"):
                break
            chunks.append(re.sub(rb"\s+", b"", line))
        return b"".join(chunks).decode("ascii", "ignore")

    def write_subset_fasta(self, ids: Iterable[str], output: Path, mode: str = "index", stream_threshold: int = 50000) -> int:
        """Write a FASTA subset.

        mode='index' uses the SQLite offset index and random seeks; this is fast for small
        subsets. mode='stream' scans the FASTA once and writes matching IDs; this is much
        faster for large subsets on network filesystems because it avoids millions of
        random seeks. mode='auto' chooses stream when the target set is large.
        """
        id_list = list(dict.fromkeys(ids))
        if not id_list:
            output.write_text("", encoding="utf-8")
            return 0
        if mode not in {"auto", "index", "stream"}:
            mode = "auto"
        chosen = mode
        if mode == "auto":
            chosen = "stream" if len(id_list) >= stream_threshold else "index"

        if chosen == "stream":
            return self.write_subset_fasta_streaming(id_list, output)

        n = 0
        with output.open("w", encoding="utf-8") as out:
            for sid in id_list:
                seq = self.sequence(sid)
                if seq:
                    n += 1
                    out.write(f">{sid}\n{seq}\n")
        return n

    def write_subset_fasta_streaming(self, ids: Iterable[str], output: Path) -> int:
        wanted = set(ids)
        n = 0
        write_current = False
        with self.fasta.open("r", encoding="utf-8", errors="replace") as inp, output.open("w", encoding="utf-8") as out:
            for line in inp:
                if line.startswith(">"):
                    sid = line[1:].strip().split()[0]
                    write_current = sid in wanted
                    if write_current:
                        n += 1
                        out.write(f">{sid}\n")
                elif write_current:
                    out.write(line)
        return n


def clean_term_for_quality(term: str, type_name: str, type_index: int, quality: bool) -> Optional[str]:
    if quality:
        m = re.search(r"\{(.+?)\}", term)
        if type_name == "GO" and m:
            if GOEV in m.group(1):
                return None
        elif type_index not in (3, 4) and m:  # Perl skips GO and KEYWORD for this ECO filtering branch
            ecos = m.group(1).split("||")
            for eco in ecos:
                eco = eco.replace("ECO:", "").strip()
                if eco in ECO_SET:
                    return None
    return term


def clean_annotation_terms(raw_terms: str, type_name: str, type_index: int, quality: bool, extend_go: bool) -> List[str]:
    if not raw_terms:
        return []
    out: List[str] = []
    for term in raw_terms.split(";"):
        term = term.strip()
        if not term:
            continue
        term = clean_term_for_quality(term, type_name, type_index, quality)
        if term is None:
            continue
        if extend_go and type_name == "GO":
            term = re.sub(r"^GO:", "", term)
            term = term.replace("{", ":", 1)
            term = re.sub(r":\w{2,3}\}$", "", term)
        else:
            term = re.sub(r"\{.+?\}", "", term)
        if term:
            out.append(term)
    return out


def goslim_from_go_terms(go_terms: Iterable[str], extend_go: bool) -> List[str]:
    slims = []
    for go in go_terms:
        goid = ""
        if extend_go:
            first = go.split(":", 1)[0]
            if first.isdigit():
                goid = "GO:" + first
        else:
            if go.startswith("GO:"):
                goid = go.split(";", 1)[0]
        if goid in SLIM:
            slims.append(goid)
    return sorted(set(slims))


def annot_cache_path(annot: Path, cfg: Config) -> Path:
    flags = f"q{int(cfg.quality)}_go{int(cfg.extend_go)}_goslim{int(cfg.goslim)}"
    return annot.with_suffix(annot.suffix + f".{flags}.sma3s_annot.sqlite")


def flush_freqs(con: sqlite3.Connection, counter: Counter, n_annot: Counter) -> None:
    if counter:
        con.executemany(
            "INSERT INTO freqs(type, term, count) VALUES(?,?,?) "
            "ON CONFLICT(type, term) DO UPDATE SET count=count+excluded.count",
            [(t, term, c) for (t, term), c in counter.items()],
        )
        counter.clear()
    if n_annot:
        con.executemany(
            "INSERT INTO n_annot(type, count) VALUES(?,?) "
            "ON CONFLICT(type) DO UPDATE SET count=count+excluded.count",
            list(n_annot.items()),
        )
        n_annot.clear()
    con.commit()


def ensure_annotation_cache(annot: Path, cfg: Config) -> Path:
    db_path = annot_cache_path(annot, cfg)
    if db_path.exists() and db_path.stat().st_mtime >= annot.stat().st_mtime:
        return db_path
    print(f"Creating SQLite annotation cache: {db_path.name}", flush=True)
    if db_path.exists():
        db_path.unlink()
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("CREATE TABLE annotations(id TEXT PRIMARY KEY, annot TEXT NOT NULL)")
    con.execute("CREATE TABLE freqs(type TEXT NOT NULL, term TEXT NOT NULL, count INTEGER NOT NULL, PRIMARY KEY(type, term))")
    con.execute("CREATE TABLE n_annot(type TEXT PRIMARY KEY, count INTEGER NOT NULL)")
    types = cfg.types
    freq_counter: Counter = Counter()
    n_counter: Counter = Counter()
    batch = []
    with annot.open("r", encoding="utf-8", errors="replace") as fh:
        for line_no, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            uid, score = parts[0], parts[1]
            raw_fields = parts[2:2 + len(TYPES_BASE)]
            cleaned_fields: List[str] = []
            go_terms_for_slim: List[str] = []
            for i, type_name in enumerate(TYPES_BASE):
                raw = raw_fields[i] if i < len(raw_fields) else ""
                terms = clean_annotation_terms(raw, type_name, i, cfg.quality, cfg.extend_go)
                if type_name == "GO":
                    go_terms_for_slim = terms
                cleaned_fields.append(";".join(terms))
                for term in terms:
                    freq_counter[(type_name, term)] += 1
                    n_counter[type_name] += 1
            if cfg.goslim:
                slim_terms = goslim_from_go_terms(go_terms_for_slim, cfg.extend_go)
                cleaned_fields.append(";".join(slim_terms))
                for term in slim_terms:
                    freq_counter[("GOSLIM", term)] += 1
                    n_counter["GOSLIM"] += 1
            annot_line = score + "\t" + "\t".join(cleaned_fields)
            batch.append((uid, annot_line))
            if len(batch) >= 10000:
                con.executemany("INSERT OR REPLACE INTO annotations VALUES(?,?)", batch)
                batch.clear()
            if line_no % 100000 == 0:
                if batch:
                    con.executemany("INSERT OR REPLACE INTO annotations VALUES(?,?)", batch)
                    batch.clear()
                flush_freqs(con, freq_counter, n_counter)
    if batch:
        con.executemany("INSERT OR REPLACE INTO annotations VALUES(?,?)", batch)
    flush_freqs(con, freq_counter, n_counter)
    con.close()
    return db_path


class AnnotationStore:
    def __init__(self, annot: Path, cfg: Config):
        self.db_path = ensure_annotation_cache(annot, cfg)
        self.con = sqlite3.connect(self.db_path)
        self.freq_cache: Dict[Tuple[str, str], int] = {}
        self.n_cache: Dict[str, int] = {}

    def close(self) -> None:
        with contextlib.suppress(Exception):
            self.con.close()

    def get(self, seq_id: str) -> Optional[str]:
        row = self.con.execute("SELECT annot FROM annotations WHERE id=?", (seq_id,)).fetchone()
        return row[0] if row else None

    def freq(self, type_name: str, term: str) -> int:
        key = (type_name, term)
        if key not in self.freq_cache:
            row = self.con.execute("SELECT count FROM freqs WHERE type=? AND term=?", key).fetchone()
            self.freq_cache[key] = int(row[0]) if row else 0
        return self.freq_cache[key]

    def n_annot(self, type_name: str) -> int:
        if type_name not in self.n_cache:
            row = self.con.execute("SELECT count FROM n_annot WHERE type=?", (type_name,)).fetchone()
            self.n_cache[type_name] = int(row[0]) if row else 0
        return self.n_cache[type_name]


def parse_mmseqs_float(x: str, percent: bool = False) -> float:
    try:
        v = float(x)
    except ValueError:
        return 0.0
    # MMseqs coverage can be fraction-like depending on version/format. Normalize to percent.
    if percent and 0.0 <= v <= 1.0:
        v *= 100.0
    return v


def run_mmseqs_easy_search(query: Path, target: Path, out_file: Path, cfg: Config, max_seqs: int, force: bool = False) -> None:
    if out_file.exists() and out_file.stat().st_size > 0 and not force:
        print(f"Using existing MMseqs2 search output: {out_file}", flush=True)
        return
    tmp = cfg.tmpdir / (out_file.name + ".tmp")
    tmp.mkdir(parents=True, exist_ok=True)
    fmt = "query,target,pident,qcov,alnlen,evalue,tcov,qlen,tlen"
    cmd = [
        cfg.mmseqs_bin, "easy-search", str(query), str(target), str(out_file), str(tmp),
        "--threads", str(cfg.cpus),
        "-e", "1e-6",
        "--max-seqs", str(max_seqs),
        "-s", str(cfg.sensitivity),
        "--format-output", fmt,
        "--mask", "1" if cfg.low_complexity_filter else "0",
    ]
    if cfg.split_memory_limit:
        cmd += ["--split-memory-limit", str(cfg.split_memory_limit)]
    if cfg.db_load_mode is not None:
        cmd += ["--db-load-mode", str(cfg.db_load_mode)]
    if cfg.nucl:
        cmd += ["--search-type", "2"]
    run_cmd(cmd)
    if not cfg.keep_tmp:
        shutil.rmtree(tmp, ignore_errors=True)



def hits_db_path(mmseqs_file: Path) -> Path:
    return mmseqs_file.with_suffix(mmseqs_file.suffix + ".hits.sqlite")


def load_hits_to_sqlite(mmseqs_file: Path, cfg: Config, mapping: Dict[str, str]) -> Path:
    db_path = hits_db_path(mmseqs_file)
    if db_path.exists() and db_path.stat().st_mtime >= mmseqs_file.stat().st_mtime:
        return db_path
    print("Filtering MMseqs2 output and loading hits into SQLite", flush=True)
    if db_path.exists():
        db_path.unlink()
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute(
        "CREATE TABLE hits(q TEXT NOT NULL, t TEXT NOT NULL, pident REAL, qcov REAL, alnlen INTEGER, "
        "evalue REAL, tcov REAL, qlen INTEGER, tlen INTEGER, rank INTEGER, UNIQUE(q,t))"
    )
    con.execute("CREATE TABLE hit_ids(id TEXT PRIMARY KEY)")
    rank_by_q: Dict[str, int] = defaultdict(int)
    batch_hits = []
    batch_ids = []
    with mmseqs_file.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                continue
            q, t = cols[0], cols[1]
            if cfg.training and mapping.get(q) == t:
                continue
            pident = parse_mmseqs_float(cols[2], percent=True)
            qcov = parse_mmseqs_float(cols[3], percent=True)
            alnlen = int(float(cols[4])) if cols[4] else 0
            evalue = float(cols[5]) if cols[5] else 1.0
            tcov = parse_mmseqs_float(cols[6], percent=True) if len(cols) > 6 else 0.0
            qlen = int(float(cols[7])) if len(cols) > 7 and cols[7] else 0
            tlen = int(float(cols[8])) if len(cols) > 8 and cols[8] else 0
            rank_by_q[q] += 1
            batch_hits.append((q, t, pident, qcov, alnlen, evalue, tcov, qlen, tlen, rank_by_q[q]))
            batch_ids.append((t,))
            if len(batch_hits) >= 50000:
                con.executemany("INSERT OR IGNORE INTO hits VALUES(?,?,?,?,?,?,?,?,?,?)", batch_hits)
                con.executemany("INSERT OR IGNORE INTO hit_ids VALUES(?)", batch_ids)
                con.commit()
                batch_hits.clear(); batch_ids.clear()
    if batch_hits:
        con.executemany("INSERT OR IGNORE INTO hits VALUES(?,?,?,?,?,?,?,?,?,?)", batch_hits)
        con.executemany("INSERT OR IGNORE INTO hit_ids VALUES(?)", batch_ids)
    con.execute("CREATE INDEX idx_hits_q_rank ON hits(q, rank)")
    con.execute("CREATE INDEX idx_hits_t ON hits(t)")
    con.commit()
    con.close()
    return db_path


class HitsStore:
    def __init__(self, db_path: Path):
        self.con = sqlite3.connect(db_path)

    def close(self) -> None:
        with contextlib.suppress(Exception):
            self.con.close()

    def get_hits(self, query: str) -> List[Hit]:
        rows = self.con.execute(
            "SELECT q,t,pident,qcov,alnlen,evalue,tcov,qlen,tlen,rank FROM hits WHERE q=? ORDER BY rank", (query,)
        ).fetchall()
        return [Hit(*row) for row in rows]

    def distinct_targets(self) -> Iterator[str]:
        cur = self.con.execute("SELECT id FROM hit_ids")
        for (sid,) in cur:
            yield sid

    def reciprocal_candidate_targets(self, cfg: Config, lenq: Dict[str, int]) -> List[str]:
        """Return only targets that could pass the preliminary A2 filters.

        The old batched reciprocal search used all distinct targets from the first MMseqs2
        result. For large --max-seqs values this can create a huge reciprocal FASTA.
        A2 only needs reciprocal checks for hits that already pass identity/Rost, coverage
        and p-value thresholds, so filtering here avoids unnecessary FASTA extraction and
        reciprocal searching while preserving the annotation logic.
        """
        candidates: Set[str] = set()
        total_rows = 0
        total_targets = 0
        for (n,) in self.con.execute("SELECT COUNT(*) FROM hit_ids"):
            total_targets = int(n)
        cur = self.con.execute("SELECT q,t,pident,qcov,alnlen,evalue,tlen FROM hits")
        for q, t, pident, qcov, alnlen, evalue, tlen in cur:
            total_rows += 1
            try:
                evalue_f = float(evalue)
            except Exception:
                evalue_f = 1.0
            if pvalue_from_evalue(evalue_f) > cfg.pv:
                continue
            alnlen_i = int(alnlen or 0)
            threshold = cfg.id_orthologue
            if (not cfg.user_changed_thresholds) and (not cfg.nucl):
                threshold = calculate_rost(cfg.rost, alnlen_i)
            if float(pident or 0.0) < threshold:
                continue
            qlen = int(lenq.get(q, 0))
            tlen_i = int(tlen or 0)
            qc_s = calculate_qc_subject(qlen, float(qcov or 0.0), tlen_i)
            if qc_s < cfg.cov_orthologue:
                continue
            candidates.add(t)
        print(
            f"A2 reciprocal candidates after preliminary filters: {len(candidates):,} "
            f"de {total_targets:,} distinct targets ({total_rows:,} hits)",
            flush=True,
        )
        return sorted(candidates)



def reciprocal_db_path(mmseqs_file: Path) -> Path:
    return mmseqs_file.with_suffix(mmseqs_file.suffix + ".reciprocal.sqlite")


class ReciprocalBestStore:
    def __init__(self, db_path: Optional[Path]):
        self.db_path = db_path
        self.con: Optional[sqlite3.Connection] = None
        if db_path is not None and db_path.exists():
            self.con = sqlite3.connect(db_path)

    def close(self) -> None:
        if self.con is not None:
            with contextlib.suppress(Exception):
                self.con.close()

    def get(self, target_id: str) -> Optional[str]:
        if self.con is None:
            return None
        row = self.con.execute("SELECT q FROM reciprocal_best WHERE t=?", (target_id,)).fetchone()
        return row[0] if row else None


def build_reciprocal_best_db(cfg: Config, query_norm: Path, fasta_store: FastaStore, target_ids: Iterable[str]) -> Optional[Path]:
    """Build a disk-backed target->best-query map for annotator 2.

    A dict can become very large with huge UniProt/TrEMBL searches. SQLite keeps this map
    shared on disk so annotation workers do not duplicate it in RAM.
    """
    if "2" not in cfg.annotator:
        return None
    db_path = reciprocal_db_path(cfg.mmseqs_file)
    if db_path.exists() and db_path.stat().st_mtime >= cfg.mmseqs_file.stat().st_mtime and not cfg.force_search:
        return db_path

    print("Preparing batched MMseqs2 reciprocal search for annotator 2", flush=True)
    if db_path.exists():
        db_path.unlink()
    with tempfile.TemporaryDirectory(dir=cfg.tmpdir) as td:
        tdpath = Path(td)
        subset_fasta = tdpath / "reciprocal_targets.fasta"
        target_id_list = list(dict.fromkeys(target_ids))
        mode = cfg.reciprocal_fasta_mode
        chosen = mode
        if mode == "auto":
            chosen = "stream" if len(target_id_list) >= cfg.reciprocal_stream_threshold else "index"
        print(
            f"Creating reciprocal FASTA with {len(target_id_list):,} targets "
            f"(mode {chosen}; stream threshold={cfg.reciprocal_stream_threshold:,})",
            flush=True,
        )
        n = fasta_store.write_subset_fasta(
            target_id_list,
            subset_fasta,
            mode=mode,
            stream_threshold=cfg.reciprocal_stream_threshold,
        )
        print(f"Reciprocal FASTA created: {n:,} sequences", flush=True)
        con = sqlite3.connect(db_path)
        con.execute("PRAGMA journal_mode=OFF")
        con.execute("PRAGMA synchronous=OFF")
        con.execute("CREATE TABLE reciprocal_best(t TEXT PRIMARY KEY, q TEXT NOT NULL)")
        if n == 0:
            con.commit(); con.close()
            return db_path
        out_file = tdpath / "reciprocal.m8"
        run_mmseqs_easy_search(subset_fasta, query_norm, out_file, cfg, max_seqs=1, force=True)
        batch = []
        with out_file.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 2:
                    batch.append((cols[0], cols[1]))
                    if len(batch) >= 50000:
                        con.executemany("INSERT OR IGNORE INTO reciprocal_best VALUES(?,?)", batch)
                        con.commit()
                        batch.clear()
        if batch:
            con.executemany("INSERT OR IGNORE INTO reciprocal_best VALUES(?,?)", batch)
        con.execute("CREATE INDEX IF NOT EXISTS idx_recip_t ON reciprocal_best(t)")
        con.commit()
        con.close()
    return db_path

def run_mmseqs_cluster(ids: List[str], fasta_store: FastaStore, cfg: Config, name: str) -> List[str]:
    # Replacement for blastclust when -uniprot is used. For UniRef, the original script skips clustering.
    if len(ids) <= 1:
        return ids
    with tempfile.TemporaryDirectory(dir=cfg.tmpdir) as td:
        tdpath = Path(td)
        subset = tdpath / f"cluster_{name}.fasta"
        if fasta_store.write_subset_fasta(ids, subset, mode="index") <= 1:
            return ids
        prefix = tdpath / "clu"
        tmp = tdpath / "tmp"
        cmd = [
            cfg.mmseqs_bin, "easy-cluster", str(subset), str(prefix), str(tmp),
            "--min-seq-id", "0.95", "-c", "0.95", "--cov-mode", "0",
            "--threads", str(cfg.cluster_threads),
        ]
        try:
            run_cmd(cmd, quiet=True)
        except SystemExit:
            # Fallback: no clustering if easy-cluster behaves differently on this installation.
            return ids
        tsv = Path(str(prefix) + "_cluster.tsv")
        if not tsv.exists():
            return ids
        clusters: Dict[str, List[str]] = defaultdict(list)
        with tsv.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 2:
                    clusters[cols[0]].append(cols[1])
        if not clusters:
            return ids
        return [" ".join(members) for members in clusters.values()]


def format_hit_metrics(hits: Sequence[Hit]) -> List[str]:
    """Return semicolon-joined hit metrics for optional output columns.

    Columns are kept positional: target;target, identity;identity, qcov;qcov, etc.
    This lets a downstream parser match each metric by index.
    """
    if not hits:
        return ["", "", "", "", "", ""]
    targets: List[str] = []
    identities: List[str] = []
    qcovs: List[str] = []
    tcovs: List[str] = []
    alnlen: List[str] = []
    evalues: List[str] = []
    for h in hits:
        targets.append(h.target)
        identities.append(f"{h.pident:.3f}")
        qcovs.append(f"{h.qcov:.3f}")
        tcovs.append(f"{h.tcov:.3f}")
        alnlen.append(str(h.alnlen))
        evalues.append(f"{h.evalue:.3g}")
    return [
        ";".join(targets),
        ";".join(identities),
        ";".join(qcovs),
        ";".join(tcovs),
        ";".join(alnlen),
        ";".join(evalues),
    ]


def unique_hits_by_target(hits: Sequence[Hit], target_ids: Sequence[str]) -> List[Hit]:
    """Keep hit objects for target_ids in requested order, removing duplicates."""
    by_target: Dict[str, Hit] = {}
    for h in hits:
        by_target.setdefault(h.target, h)
    out: List[Hit] = []
    seen: Set[str] = set()
    for tid in target_ids:
        if tid in seen:
            continue
        h = by_target.get(tid)
        if h is not None:
            out.append(h)
            seen.add(tid)
    return out


class OutputWriter:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.stats = Counter({"Annotations": 0, "a1": 0, "a2": 0, "a3": 0, "a123": 0, "a12": 0})
        for t in cfg.types:
            self.stats[t] = 0
        self.slims: Dict[str, Counter] = defaultdict(Counter)
        self.paths = Counter()
        self.kws: Dict[str, Counter] = defaultdict(Counter)
        self.fh = cfg.annot_file.open("w", encoding="utf-8")
        self.write_header()

    def write_header(self) -> None:
        fields = ["#ID"]
        for t in self.cfg.types:
            if t == "GO" and self.cfg.extend_go:
                fields += ["GO(P)ID", "GO(P)NAME", "GO(F)ID", "GO(F)NAME", "GO(C)ID", "GO(C)NAME"]
            else:
                fields.append(t)
        if self.cfg.source:
            fields += ["ANNOTATOR", "USED_UNIPROT_SEQUENCES"]
        if self.cfg.report_hit_metrics:
            fields += [
                "ANNOTATION_HIT_SEQUENCE",
                "ANNOTATION_HIT_IDENTITY",
                "ANNOTATION_HIT_QCOV",
                "ANNOTATION_HIT_TCOV",
                "ANNOTATION_HIT_ALNLEN",
                "ANNOTATION_HIT_EVALUE",
            ]
        self.fh.write("\t".join(fields) + "\n")

    def close(self) -> None:
        self.fh.close()

    def write_empty(self, original_id: str) -> None:
        if not self.cfg.noempty:
            self.fh.write(original_id + "\n")

    def create_annot(self, original_id: str, o: str, annotator: str, used_ids: str, hit_metrics: Optional[List[str]] = None) -> Optional[str]:
        if not re.search(r"[^\t]", o or ""):
            return f"{original_id}\n" if not self.cfg.noempty else None
        ann = o.split("\t")
        self.stats["Annotations"] += 1
        self.stats[annotator] += 1
        row = [original_id]
        for i, type_name in enumerate(self.cfg.types):
            value = ann[i] if i < len(ann) else ""
            if value:
                if type_name != "GOSLIM":
                    self.stats[type_name] += 1
                if type_name == "GO" and self.cfg.extend_go:
                    pfc_id = {"P": [], "F": [], "C": []}
                    pfc_name = {"P": [], "F": [], "C": []}
                    for go in value.split(";"):
                        parts = go.split(":", 2)
                        if len(parts) == 3 and parts[1] in pfc_id:
                            goid, group, name = parts
                            pfc_id[group].append("GO:" + goid)
                            pfc_name[group].append(name)
                    for group in ("P", "F", "C"):
                        row.append(";".join(pfc_id[group]))
                        row.append(";".join(pfc_name[group]))
                else:
                    row.append(value)
                if type_name == "PATHWAY":
                    for path in value.split(";"):
                        if path:
                            self.paths[path] += 1
                elif type_name == "KEYWORD":
                    for kw in value.split(";"):
                        for cat, words in KW_CATEGORIES.items():
                            if kw in words:
                                self.kws[cat][kw] += 1
                elif type_name == "GOSLIM":
                    for goslim in value.split(";"):
                        if goslim in SLIM:
                            g, _name = SLIM[goslim].split(":", 1)
                            self.slims[g][goslim] += 1
            else:
                row.append("")
                if type_name == "GO" and self.cfg.extend_go:
                    row += ["", "", "", "", ""]
        if self.cfg.source:
            row += [annotator, used_ids]
        if self.cfg.report_hit_metrics:
            row += hit_metrics if hit_metrics is not None else ["", "", "", "", "", ""]
        return "\t".join(row) + "\n"

    def write_annotation(self, original_id: str, o: str, annotator: str, used_ids: str, hit_metrics: Optional[List[str]] = None) -> None:
        line = self.create_annot(original_id, o, annotator, used_ids, hit_metrics)
        if line:
            self.fh.write(line)

    def write_summary(self, n_queries: int) -> Path:
        stat_file = Path(re.sub(r"\.tsv$", "_summary.tsv", str(self.cfg.annot_file)))
        with stat_file.open("w", encoding="utf-8") as stat:
            stat.write("#Annotation summary\n")
            stat.write(f"Number of query sequences:\t{n_queries}\n")
            for name in ["Annotations"] + self.cfg.types:
                if name != "GOSLIM":
                    stat.write(f"With {name}\t{self.stats[name]}\n")
            for name in ["a1", "a2", "a3", "a12", "a123"]:
                stat.write(f"Annotator {name}\t{self.stats[name]}\n")
            stat.write("\n")
            if self.cfg.goslim:
                cat_names = {"P": "Biological process", "C": "Cellular component", "F": "Molecular function"}
                stat.write("#GO Slim\n")
                for g in sorted(self.slims):
                    stat.write(f"#Category \"{cat_names.get(g, g)}\"\n")
                    for slim, count in self.slims[g].most_common():
                        name = SLIM.get(slim, ":").split(":", 1)[1]
                        stat.write(f"{slim}\t{name}\t{count}\n")
                    stat.write("\n")
            if self.paths:
                stat.write("#UniProt Pathways\n")
                for path, count in self.paths.most_common():
                    stat.write(f"{path}\t{count}\n")
                stat.write("\n")
            stat.write("#UniProt Keyword categories\n")
            for cat in sorted(self.kws):
                stat.write(f"#Category \"{cat}\"\n")
                for kw, count in self.kws[cat].most_common():
                    stat.write(f"{kw}\t{count}\n")
                stat.write("\n")
        return stat_file



@dataclasses.dataclass
class AnnotationResult:
    line: Optional[str]
    stats: Dict[str, int]
    paths: Dict[str, int]
    kws: Dict[str, Dict[str, int]]
    slims: Dict[str, Dict[str, int]]


class ResultWriter:
    """Per-query writer used by annotation worker processes.

    It mirrors OutputWriter's accounting, but stores everything in memory and returns it
    to the parent process. Each query should emit at most one final line.
    """
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.stats = Counter()
        self.paths = Counter()
        self.kws: Dict[str, Counter] = defaultdict(Counter)
        self.slims: Dict[str, Counter] = defaultdict(Counter)
        self.line: Optional[str] = None

    def write_empty(self, original_id: str) -> None:
        if not self.cfg.noempty:
            self.line = original_id + "\n"

    def create_annot(self, original_id: str, o: str, annotator: str, used_ids: str, hit_metrics: Optional[List[str]] = None) -> Optional[str]:
        if not re.search(r"[^\t]", o or ""):
            return f"{original_id}\n" if not self.cfg.noempty else None
        ann = o.split("\t")
        self.stats["Annotations"] += 1
        self.stats[annotator] += 1
        row = [original_id]
        for i, type_name in enumerate(self.cfg.types):
            value = ann[i] if i < len(ann) else ""
            if value:
                if type_name != "GOSLIM":
                    self.stats[type_name] += 1
                if type_name == "GO" and self.cfg.extend_go:
                    pfc_id = {"P": [], "F": [], "C": []}
                    pfc_name = {"P": [], "F": [], "C": []}
                    for go in value.split(";"):
                        parts = go.split(":", 2)
                        if len(parts) == 3 and parts[1] in pfc_id:
                            goid, group, name = parts
                            pfc_id[group].append("GO:" + goid)
                            pfc_name[group].append(name)
                    for group in ("P", "F", "C"):
                        row.append(";".join(pfc_id[group]))
                        row.append(";".join(pfc_name[group]))
                else:
                    row.append(value)
                if type_name == "PATHWAY":
                    for path in value.split(";"):
                        if path:
                            self.paths[path] += 1
                elif type_name == "KEYWORD":
                    for kw in value.split(";"):
                        for cat, words in KW_CATEGORIES.items():
                            if kw in words:
                                self.kws[cat][kw] += 1
                elif type_name == "GOSLIM":
                    for goslim in value.split(";"):
                        if goslim in SLIM:
                            g, _name = SLIM[goslim].split(":", 1)
                            self.slims[g][goslim] += 1
            else:
                row.append("")
                if type_name == "GO" and self.cfg.extend_go:
                    row += ["", "", "", "", ""]
        if self.cfg.source:
            row += [annotator, used_ids]
        if self.cfg.report_hit_metrics:
            row += hit_metrics if hit_metrics is not None else ["", "", "", "", "", ""]
        return "\t".join(row) + "\n"

    def write_annotation(self, original_id: str, o: str, annotator: str, used_ids: str, hit_metrics: Optional[List[str]] = None) -> None:
        self.line = self.create_annot(original_id, o, annotator, used_ids, hit_metrics)

    def result(self) -> AnnotationResult:
        return AnnotationResult(
            line=self.line,
            stats=dict(self.stats),
            paths=dict(self.paths),
            kws={cat: dict(counter) for cat, counter in self.kws.items()},
            slims={cat: dict(counter) for cat, counter in self.slims.items()},
        )


def merge_annotation_result(writer: OutputWriter, result: AnnotationResult) -> None:
    if result.line:
        writer.fh.write(result.line)
    writer.stats.update(result.stats)
    writer.paths.update(result.paths)
    for cat, values in result.kws.items():
        writer.kws[cat].update(values)
    for cat, values in result.slims.items():
        writer.slims[cat].update(values)


_WORKER_CFG: Optional[Config] = None
_WORKER_MAPPING: Dict[str, str] = {}
_WORKER_LENQ: Dict[str, int] = {}
_WORKER_HITS: Optional[HitsStore] = None
_WORKER_FASTA: Optional[FastaStore] = None
_WORKER_ANN: Optional[AnnotationStore] = None
_WORKER_RECIP: Optional[ReciprocalBestStore] = None


def init_annotation_worker(
    cfg: Config,
    hits_db: str,
    fasta_file: str,
    ref_annot: str,
    reciprocal_db: Optional[str],
    mapping: Dict[str, str],
    lenq: Dict[str, int],
) -> None:
    global _WORKER_CFG, _WORKER_MAPPING, _WORKER_LENQ, _WORKER_HITS, _WORKER_FASTA, _WORKER_ANN, _WORKER_RECIP
    _WORKER_CFG = cfg
    _WORKER_MAPPING = mapping
    _WORKER_LENQ = lenq
    _WORKER_HITS = HitsStore(Path(hits_db))
    _WORKER_FASTA = FastaStore(Path(fasta_file))
    _WORKER_ANN = AnnotationStore(Path(ref_annot), cfg)
    _WORKER_RECIP = ReciprocalBestStore(Path(reciprocal_db)) if reciprocal_db else ReciprocalBestStore(None)


def annotate_query_worker(index: int) -> AnnotationResult:
    if not all([_WORKER_CFG, _WORKER_HITS, _WORKER_FASTA, _WORKER_ANN, _WORKER_RECIP]):
        raise RuntimeError("Annotation worker was not initialized")
    q = f"s{index}"
    writer = ResultWriter(_WORKER_CFG)
    hits = _WORKER_HITS.get_hits(q)
    annotate_query(q, hits, _WORKER_CFG, _WORKER_MAPPING, _WORKER_LENQ, _WORKER_ANN, _WORKER_FASTA, _WORKER_RECIP, writer)
    return writer.result()

def score_tuple(score: str) -> Tuple[int, int]:
    try:
        a, b = score.split(",", 1)
        return int(float(a or 0)), int(float(b or 0))
    except Exception:
        return 0, 0


def first_or_empty(value: str) -> str:
    return value.split(";", 1)[0] if value else ""


def build_o_from_annotation(ann_line: str, types: List[str]) -> Tuple[str, Tuple[int, int]]:
    parts = ann_line.split("\t")
    score = score_tuple(parts[0] if parts else "")
    fields = parts[1:]
    # GN and DE: only first term, remaining fields unchanged.
    out = []
    out.append(first_or_empty(fields[0]) if len(fields) > 0 else "")
    out.append(first_or_empty(fields[1]) if len(fields) > 1 else "")
    for idx in range(2, len(types)):
        out.append(fields[idx] if idx < len(fields) else "")
    return "\t".join(out), score


def annotate_query(
    q: str,
    hits: List[Hit],
    cfg: Config,
    mapping: Dict[str, str],
    lenq: Dict[str, int],
    ann_store: AnnotationStore,
    fasta_store: FastaStore,
    reciprocal_best: Dict[str, str],
    writer: OutputWriter,
) -> None:
    original = mapping[q]
    if not hits:
        writer.write_empty(original)
        return
    o1 = ""
    o2 = ""
    s1 = (0, 0)
    s2 = 0
    best_name2 = hits[0].target
    best_hit2 = hits[0]

    # Annotator 1: direct UniProt/UniRef hit.
    if "1" in cfg.annotator:
        h = hits[0]
        pv_hsp = pvalue_from_evalue(h.evalue)
        tlen = h.tlen or fasta_store.length(h.target)
        qc_s = h.tcov if h.tcov else calculate_qc_subject(lenq.get(h.query, h.qlen), h.qcov, tlen)
        if h.pident >= cfg.id_uniprot and qc_s >= cfg.cov_uniprot and pv_hsp <= cfg.pv:
            ann_line = ann_store.get(h.target)
            if ann_line:
                o1, s1 = build_o_from_annotation(ann_line, cfg.types)
            if cfg.annotator == "1":
                writer.write_annotation(original, o1, "a1", h.target, format_hit_metrics([h]))
                return
        elif cfg.annotator == "1":
            writer.write_empty(original)
            return

    # Annotator 2: orthologues / reciprocal best hit. Batched reciprocal map replaces thousands of subprocesses.
    if "2" in cfg.annotator:
        for h in hits:
            if s1[0] == 1 or s1[1] == cfg.max_score:
                break
            ann_line = ann_store.get(h.target)
            if not ann_line:
                continue
            parts = ann_line.split("\t")
            s3 = score_tuple(parts[0] if parts else "")
            if not (s3[0] == 1 or s3[1] > s1[1]):
                continue
            pv_hsp = pvalue_from_evalue(h.evalue)
            tlen = h.tlen or fasta_store.length(h.target)
            qc_s = h.tcov if h.tcov else calculate_qc_subject(lenq.get(h.query, h.qlen), h.qcov, tlen)
            threshold = cfg.id_orthologue
            if not cfg.user_changed_thresholds and not cfg.nucl:
                threshold = calculate_rost(cfg.rost, h.alnlen)
            if h.pident >= threshold and qc_s >= cfg.cov_orthologue and pv_hsp <= cfg.pv:
                if reciprocal_best.get(h.target) == q:
                    o1, s1 = build_o_from_annotation(ann_line, cfg.types)
                    best_name2 = h.target
                    best_hit2 = h
                    if cfg.annotator in ("2", "12"):
                        writer.write_annotation(original, o1, "a2", h.target, format_hit_metrics([h]))
                        return
            elif cfg.annotator in ("2", "12"):
                writer.write_empty(original)
                return

    # Annotator 3: enrichment among high-similarity hits.
    id_fastas: List[str] = []
    if "3" in cfg.annotator:
        for h in hits:
            threshold = calculate_rost(cfg.rost, h.alnlen)
            if h.pident > threshold:
                id_fastas.append(h.target)
        if len(id_fastas) > 1:
            clusters: List[str]
            if not cfg.uniref:
                clusters = run_mmseqs_cluster(id_fastas, fasta_store, cfg, q)
            else:
                clusters = id_fastas[:]
            freq: Dict[str, Counter] = defaultdict(Counter)
            n_annot: Counter = Counter()
            for cluster in clusters:
                ids = cluster.split()
                annot_seen: Dict[str, set] = defaultdict(set)
                annot_terms: Dict[str, List[str]] = defaultdict(list)
                for sid in ids:
                    ann_line = ann_store.get(sid)
                    if not ann_line:
                        continue
                    ann_fields = ann_line.split("\t")
                    fields = ann_fields[1:]
                    for i, type_name in enumerate(cfg.types):
                        terms = fields[i] if i < len(fields) else ""
                        if not terms:
                            continue
                        for term in terms.split(";"):
                            term_m = filter_cs(term)
                            if term_m in annot_seen[type_name]:
                                continue
                            annot_seen[type_name].add(term_m)
                            annot_terms[type_name].append(term)
                for type_name, terms in annot_terms.items():
                    for term in terms:
                        freq[type_name][term] += 1
                        n_annot[type_name] += 1

            out_fields: List[str] = []
            for type_name in cfg.types:
                min_pv = 10.0
                output: List[str] = []
                for annot in sorted(freq.get(type_name, {})):
                    if not annot:
                        continue
                    pval = compute_hyper_p_value(
                        freq[type_name][annot],
                        n_annot[type_name],
                        ann_store.freq(type_name, annot),
                        ann_store.n_annot(type_name),
                    )
                    if pval > cfg.pv:
                        continue
                    if type_name == "GENENAME":
                        s_gn = calculate_gnscore(annot)
                        if s_gn > s2:
                            s2 = s_gn
                            output = [annot]
                    elif type_name == "DESCRIPTION":
                        if pval < min_pv:
                            min_pv = pval
                            output = [annot]
                    else:
                        output.append(annot)
                out_fields.append(";".join(output))
            o2 = "\t".join(out_fields)

        if cfg.annotator == "3":
            if o2:
                writer.write_annotation(original, o2, "a3", ";".join(id_fastas), format_hit_metrics(unique_hits_by_target(hits, id_fastas)))
            else:
                writer.write_empty(original)
            return

    if not o1 and not o2:
        writer.write_empty(original)
        return
    if not o2:
        writer.write_annotation(original, o1, "a12", best_name2, format_hit_metrics([best_hit2]))
    elif not o1:
        writer.write_annotation(original, o2, "a3", ";".join(id_fastas), format_hit_metrics(unique_hits_by_target(hits, id_fastas)))
    else:
        cols1 = o1.split("\t")
        cols2 = o2.split("\t")
        merged: List[str] = []
        for i, type_name in enumerate(cfg.types):
            c1 = cols1[i] if i < len(cols1) else ""
            c2 = cols2[i] if i < len(cols2) else ""
            if i in (0, 1):
                if s1[1] >= s2 and c1:
                    merged.append(c1)
                elif c2:
                    merged.append(c2)
                else:
                    merged.append("")
            else:
                terms = []
                if c1:
                    terms += c1.split(";")
                if c2:
                    terms += c2.split(";")
                merged.append(";".join(sorted(set(t for t in terms if t))))
        writer.write_annotation(original, "\t".join(merged), "a123", ";".join(id_fastas), format_hit_metrics(unique_hits_by_target(hits, [best_name2] + id_fastas)))


def validate_and_prepare(args: argparse.Namespace) -> Config:
    if not args.query_file or not args.database:
        die("Minimum usage: sma3s.py -i input.fasta -d uniprot.dat|target.fasta")
    annotator = str(args.annotator)
    if not re.fullmatch(r"1?2?3?", annotator) or not annotator:
        die("-a must be 1, 2, 3, 12, 13, 23 or 123")
    query = Path(args.query_file)
    database_input = Path(args.database)
    compressed_dat: Optional[Path] = None
    if database_input.suffix.lower() == ".gz":
        compressed_dat = database_input
        dat = decompressed_dat_path(database_input)
    else:
        dat = database_input
    if not query.exists():
        die(f"Query file does not exist: {query}")
    genus_filter = args.genus.strip() if args.genus else None
    exclude_genus_filter = args.exclude_genus.strip() if args.exclude_genus else None
    family_filter = args.family.strip() if args.family else None
    exclude_family_filter = args.exclude_family.strip() if args.exclude_family else None
    order_filter = args.order.strip() if args.order else None
    exclude_order_filter = args.exclude_order.strip() if args.exclude_order else None

    for label, include_value, exclude_value in [
        ("genus", genus_filter, exclude_genus_filter),
        ("family", family_filter, exclude_family_filter),
        ("order", order_filter, exclude_order_filter),
    ]:
        if include_value and exclude_value and include_value.lower() == exclude_value.lower():
            die(f"Cannot select and exclude the same {label}: {include_value}")

    tax_filters_active = any([genus_filter, exclude_genus_filter, family_filter, exclude_family_filter, order_filter, exclude_order_filter])
    tax_exclusion_active = any([exclude_genus_filter, exclude_family_filter, exclude_order_filter])
    if tax_filters_active and dat.suffix != ".dat":
        die("Taxonomic filters (--genus/--family/--order and their --exclude-* counterparts) can only be applied when -d is a UniProt .dat file. If you provide a FASTA file, it must already be filtered and must have its companion .annot file")
    fasta = Path(str(dat))
    if fasta.suffix == ".dat":
        suffix = tax_filter_suffix(
            genus_filter, exclude_genus_filter,
            family_filter, exclude_family_filter,
            order_filter, exclude_order_filter,
        )
        if suffix:
            fasta = dat.with_name(f"{dat.stem}{suffix}.fasta")
        else:
            fasta = fasta.with_suffix(".fasta")
    annot_companion = fasta.with_suffix(".annot")
    compressed_exists = compressed_dat is not None and compressed_dat.exists()
    if not dat.exists() and not fasta.exists() and not compressed_exists:
        requested = compressed_dat if compressed_dat is not None else dat
        die(f"Neither the requested database nor the derived FASTA exists: {requested} / {fasta}")
    if dat.exists() and dat.suffix != ".dat" and not annot_companion.exists():
        # User may pass target.fasta in -d.
        fasta = dat
        annot_companion = fasta.with_suffix(".annot")
    if fasta.exists() and not annot_companion.exists() and (not dat.exists() or dat.suffix != ".dat"):
        die(f"FASTA input was provided but the companion .annot file is missing: {annot_companion}")

    user_changed = any(x is not None for x in [args.cov1, args.cov2, args.id1, args.id2])
    id1 = float(args.id1) if args.id1 is not None else 90.0
    id2 = float(args.id2) if args.id2 is not None else 75.0
    cov1 = float(args.cov1) if args.cov1 is not None else 90.0
    cov2 = float(args.cov2) if args.cov2 is not None else 80.0
    if args.nucl and not user_changed:
        id1 = min(id1, 51.0); id2 = min(id2, 51.0); cov1 = min(cov1, 51.0); cov2 = min(cov2, 51.0)
        print("Nucleotide queries selected: identity/coverage thresholds reduced to 51%", flush=True)
    for label, value in [("id1", id1), ("id2", id2), ("cov1", cov1), ("cov2", cov2)]:
        if value < 1 or value > 100:
            die(f"{label} must be between 1 and 100")
    if args.rost < 0 or args.rost > 100 or args.pvalue < 1e-100 or args.pvalue > 1:
        die("-r must be between 0 and 100 and -p between 1e-100 and 1")

    cpus = max(1, int(args.num_threads))
    annotation_workers = int(args.annotation_workers) if args.annotation_workers is not None else max(1, cpus // 4)
    annotation_workers = max(1, min(annotation_workers, cpus))
    cluster_threads = int(args.cluster_threads) if args.cluster_threads is not None else max(1, cpus // annotation_workers)
    cluster_threads = max(1, min(cluster_threads, cpus))
    decompression_threads = (int(args.decompression_threads) if args.decompression_threads is not None else cpus)
    decompression_threads = max(1, min(decompression_threads, cpus))

    prefix = output_prefix(query, dat, uniref=not args.uniprot)
    ext = ""
    if args.annotator != "123": ext += f"_a{args.annotator}"
    if args.rost != 20: ext += f"_r{args.rost}"
    if args.pvalue != 0.1: ext += f"_p{args.pvalue}"
    if args.cov1 is not None: ext += f"_cov1{args.cov1}"
    if args.cov2 is not None: ext += f"_cov2{args.cov2}"
    if args.id1 is not None: ext += f"_id1{args.id1}"
    if args.id2 is not None: ext += f"_id2{args.id2}"
    if args.go: ext += "_go"
    if args.filter: ext += "_filter"
    if args.training: ext += "_training"
    if args.noempty: ext += "_noempty"
    if args.nopred: ext += "_nopred"
    if args.quality: ext += "_quality"
    if args.source: ext += "_source"
    if args.goslim: ext += "_goslim"
    if order_filter: ext += f"_order{safe_name_suffix(order_filter)}"
    if family_filter: ext += f"_family{safe_name_suffix(family_filter)}"
    if genus_filter: ext += f"_genus{safe_name_suffix(genus_filter)}"
    if exclude_order_filter: ext += f"_excludeOrder{safe_name_suffix(exclude_order_filter)}"
    if exclude_family_filter: ext += f"_excludeFamily{safe_name_suffix(exclude_family_filter)}"
    if exclude_genus_filter: ext += f"_excludeGenus{safe_name_suffix(exclude_genus_filter)}"
    if args.report_hit_metrics and not tax_exclusion_active: ext += "_hitmetrics"
    annot_file = Path(str(prefix) + ext + ".tsv")
    mmseqs_file = Path(args.mmseqs_file) if args.mmseqs_file else Path(str(prefix) + ".mmseqs.m8")
    tmpdir = Path(args.tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    return Config(
        annotator=annotator,
        query_file=query,
        dat_file=dat,
        fasta_file=fasta,
        compressed_dat_file=compressed_dat,
        mmseqs_file=mmseqs_file,
        annot_file=annot_file,
        rost=float(args.rost),
        pv=float(args.pvalue),
        training=bool(args.training),
        noempty=bool(args.noempty),
        nucl=bool(args.nucl),
        extend_go=bool(args.go),
        quality=bool(args.quality),
        nopred=bool(args.nopred),
        source=bool(args.source),
        uniref=not bool(args.uniprot),
        goslim=bool(args.goslim),
        cpus=cpus,
        annotation_workers=annotation_workers,
        cluster_threads=cluster_threads,
        decompression_threads=decompression_threads,
        id_uniprot=id1,
        id_orthologue=id2,
        cov_uniprot=cov1,
        cov_orthologue=cov2,
        user_changed_thresholds=user_changed,
        low_complexity_filter=bool(args.filter),
        max_seqs=int(args.max_seqs),
        sensitivity=float(args.sensitivity),
        split_memory_limit=args.split_memory_limit,
        db_load_mode=args.db_load_mode,
        tmpdir=tmpdir,
        force_search=bool(args.force_search),
        keep_tmp=bool(args.keep_tmp),
        mmseqs_bin=args.mmseqs_bin,
        genus_filter=genus_filter,
        exclude_genus_filter=exclude_genus_filter,
        family_filter=family_filter,
        exclude_family_filter=exclude_family_filter,
        order_filter=order_filter,
        exclude_order_filter=exclude_order_filter,
        report_hit_metrics=bool(args.report_hit_metrics or tax_exclusion_active),
        force_clean_outputs=bool(args.force_clean_outputs or args.force_clean),
        force_clean_db=bool(args.force_clean_db or args.force_clean),
        clean_only=bool(args.clean_only),
        reciprocal_fasta_mode=args.reciprocal_fasta_mode,
        reciprocal_stream_threshold=max(1, int(args.reciprocal_stream_threshold)),
    )


def format_elapsed(seconds: float) -> str:
    total = int(round(seconds))
    hours, rem = divmod(total, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours:02d}:{minutes:02d}:{secs:02d}"
    return f"{minutes:02d}:{secs:02d}"




def _remove_path(path: Path, label: str = "") -> bool:
    """Remove a file/symlink/directory if present and report it."""
    try:
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
            print(f"  removed {label or 'directory'}: {path}", flush=True)
            return True
        if path.exists() or path.is_symlink():
            path.unlink()
            print(f"  removed {label or 'file'}: {path}", flush=True)
            return True
    except Exception as exc:
        print(f"  warning: could not remove {path}: {exc}", flush=True)
    return False


def _remove_glob(parent: Path, pattern: str, label: str = "") -> int:
    n = 0
    if not parent.exists():
        return 0
    for path in parent.glob(pattern):
        if _remove_path(path, label):
            n += 1
    return n


def summary_path_for(annot_file: Path) -> Path:
    return Path(re.sub(r"\.tsv$", "_summary.tsv", str(annot_file)))


def cleanup_outputs(cfg: Config) -> int:
    """Remove outputs derived from the current query/database/options.

    This does not touch the input query FASTA, input .dat, reference FASTA, or .annot.
    """
    print("Removing outputs from this run", flush=True)
    n = 0
    n += _remove_path(cfg.annot_file, "TSV")
    n += _remove_path(summary_path_for(cfg.annot_file), "summary")
    n += _remove_path(cfg.mmseqs_file, "MMseqs m8")
    n += _remove_path(hits_db_path(cfg.mmseqs_file), "SQLite hits")
    n += _remove_path(reciprocal_db_path(cfg.mmseqs_file), "SQLite reciprocal")
    n += _remove_path(Path(str(cfg.query_file) + ".sma3s_ids.fasta"), "normalized FASTA")
    n += _remove_path(cfg.tmpdir / (cfg.mmseqs_file.name + ".tmp"), "MMseqs tmp")

    print(f"Outputs removed: {n}", flush=True)
    return n


def cleanup_db_derived_files(cfg: Config) -> int:
    """Remove derived reference FASTA/.annot and their caches.

    For safety this only deletes files derived from a UniProt .dat/.dat.gz input. If
    -d is an existing FASTA, it is treated as user-provided input and is not deleted.
    For .dat.gz input, the decompressed .dat is derived and is removed; the archive is kept.
    """
    is_uniprot_dat = cfg.dat_file.suffix.lower() == ".dat" and (
        cfg.dat_file.exists() or cfg.compressed_dat_file is not None
    )
    if not is_uniprot_dat:
        print("DB-derived files were not removed: -d is not a UniProt .dat/.dat.gz file; FASTA/.annot appear to be user-provided input.", flush=True)
        return 0

    ref_annot = cfg.fasta_file.with_suffix(".annot")
    print("Removing DB-derived FASTA/.annot files and associated caches", flush=True)
    n = 0
    # Caches first, then the main files.
    n += _remove_path(fasta_index_db_path(cfg.fasta_file), "FASTA index")
    n += _remove_glob(ref_annot.parent, ref_annot.name + ".*.sma3s_annot.sqlite", ".annot cache")
    n += _remove_glob(ref_annot.parent, ref_annot.name + ".*.sma3s_annot.sqlite.tmp*", "tmp .annot cache")
    n += _remove_glob(cfg.fasta_file.parent, cfg.fasta_file.name + ".sma3s_fasta_index.sqlite.tmp*", "tmp FASTA index")
    n += _remove_path(ref_annot, "derived .annot")
    n += _remove_path(cfg.fasta_file, "derived FASTA")
    if cfg.compressed_dat_file is not None:
        n += _remove_path(cfg.dat_file, "decompressed .dat")

    print(f"Derived DB files removed: {n}", flush=True)
    return n


def run_requested_cleanup(cfg: Config) -> None:
    if cfg.force_clean_outputs:
        cleanup_outputs(cfg)
    if cfg.force_clean_db:
        cleanup_db_derived_files(cfg)

def build_parser() -> argparse.ArgumentParser:
    """Create the command-line interface with documented public parameters."""
    epilog = """
Default annotation thresholds
-----------------------------
Initial MMseqs2 search: e-value <= 1e-6.
Annotator 1: identity >= 90%, target/database coverage >= 90%, p-value <= 0.1.
Annotator 2: identity >= Rost(20, alignment length) unless any identity/coverage
             threshold is manually supplied; target/database coverage >= 80%;
             p-value <= 0.1; reciprocal best hit required.
Annotator 3: identity > Rost(20, alignment length); no explicit coverage filter;
             final terms are kept by hypergeometric enrichment p-value <= 0.1.

Reference input rules
---------------------
-d reference.dat.gz -> decompress/reuse reference.dat, then create/reuse FASTA/.annot.
-d reference.dat    -> create/reuse reference.fasta and reference.annot.
-d reference.fasta  -> reference.annot must already exist next to the FASTA.
Taxonomic filters require a UniProt .dat or .dat.gz file because taxonomy is read
from OC lines. rapidgzip enables parallel .gz decompression; gzip is the fallback.

Cache behavior
--------------
FASTA offset indexes (*.sma3s_fasta_index.sqlite), cleaned annotation caches
(*.sma3s_annot.sqlite), hit databases and reciprocal-best-hit maps are reused when
present and current. Use --force-clean-outputs to clean only run outputs. Use
--force-clean-db only when you want to rebuild derived reference FASTA/.annot files
and their caches.
"""
    p = argparse.ArgumentParser(
        description=(
            "Sma3s functional annotation workflow implemented in Python with MMseqs2. "
            "It annotates query sequences against UniProt/UniRef-style references using "
            "the three original Sma3s annotators."
        ),
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"Sma3s {SMA3S_VERSION}")

    # Installation check and core workflow options.
    p.add_argument("--check-install", action="store_true",
                   help="Check Python, MMseqs2, SQLite, temporary-directory access and optional rapidgzip support, then exit. -i and -d are not required for this check.")
    p.add_argument("-a", dest="annotator", default="123", choices=["1", "2", "3", "12", "13", "23", "123"],
                   help="Annotator combination to run. 1=direct high-confidence hit; 2=orthologue with reciprocal best hit; 3=multi-hit enrichment/consensus. Default: 123.")
    p.add_argument("-i", dest="query_file", default=None,
                   help="Input query FASTA file. Protein FASTA is assumed unless -nucl is supplied.")
    p.add_argument("-d", dest="database", default=None,
                   help="Reference database: UniProt .dat or .dat.gz file, or FASTA file with a companion .annot file using the same prefix.")
    p.add_argument("-b", dest="mmseqs_file", default=None,
                   help="Existing or target MMseqs2 tabular output file. If omitted, the script creates <prefix>.mmseqs.m8.")

    # Taxonomic database filters. These create/reuse derived FASTA/.annot files.
    p.add_argument("--genus", default=None,
                   help="Keep only UniProt .dat entries whose OC lineage contains this genus, e.g. --genus Salmonella.")
    p.add_argument("--exclude-genus", default=None,
                   help="Remove UniProt .dat entries whose OC lineage contains this genus, e.g. --exclude-genus Salmonella.")
    p.add_argument("--family", default=None,
                   help="Keep only UniProt .dat entries whose OC lineage contains this family, e.g. --family Enterobacteriaceae.")
    p.add_argument("--exclude-family", default=None,
                   help="Remove UniProt .dat entries whose OC lineage contains this family, e.g. --exclude-family Enterobacteriaceae.")
    p.add_argument("--order", default=None,
                   help="Keep only UniProt .dat entries whose OC lineage contains this order, e.g. --order Enterobacterales.")
    p.add_argument("--exclude-order", default=None,
                   help="Remove UniProt .dat entries whose OC lineage contains this order, e.g. --exclude-order Enterobacterales.")

    # Optional output columns.
    p.add_argument("--report-hit-metrics", action="store_true",
                   help="Add annotation-hit metrics to the final TSV: target ID, identity, query coverage, target coverage, alignment length and e-value. Automatically enabled when any taxonomic exclusion is used.")
    p.add_argument("-source", action="store_true",
                   help="Add original Sma3s source columns: final annotator label and reference sequence IDs used for the annotation.")
    p.add_argument("-go", action="store_true",
                   help="Split GO terms into Biological Process, Molecular Function and Cellular Component columns, and keep GO names when available.")
    p.add_argument("-goslim", action="store_true",
                   help="Add a GOSLIM column and a GO Slim section in the summary file.")

    # Biological thresholds from the original Sma3s logic.
    p.add_argument("-r", dest="rost", type=float, default=20.0,
                   help="N parameter in the Rost identity curve. Larger values make annotators 2 and 3 more restrictive. Default: 20.")
    p.add_argument("-p", dest="pvalue", type=float, default=0.1,
                   help="Maximum p-value used by annotators and the enrichment test. Initial MMseqs2 search still uses e-value <= 1e-6. Default: 0.1.")
    p.add_argument("-id1", dest="id1", type=float, default=None,
                   help="Minimum percent identity for annotator 1. Default: 90.")
    p.add_argument("-cov1", dest="cov1", type=float, default=None,
                   help="Minimum target/database coverage for annotator 1. Default: 90.")
    p.add_argument("-id2", dest="id2", type=float, default=None,
                   help="Minimum percent identity for annotator 2. Default behavior is the Rost curve; if any -id/-cov threshold is manually supplied, the fixed default becomes 75 unless -id2 is set.")
    p.add_argument("-cov2", dest="cov2", type=float, default=None,
                   help="Minimum target/database coverage for annotator 2. Default: 80.")

    # Query/database interpretation and annotation-quality filters.
    p.add_argument("-nucl", action="store_true",
                   help="Treat query sequences as nucleotide sequences and run translated MMseqs2 search. If thresholds are not manually set, id/cov thresholds are reduced to 51%%.")
    p.add_argument("-training", action="store_true",
                   help="Training mode: ignore hits where the original query ID is identical to the target ID.")
    p.add_argument("-noempty", action="store_true",
                   help="Do not print unannotated query IDs in the final TSV.")
    p.add_argument("-nopred", action="store_true",
                   help="When building FASTA/.annot from UniProt .dat, discard entries with PE=4 (predicted protein evidence).")
    p.add_argument("-quality", action="store_true",
                   help="Filter lower-confidence annotation terms using the same ECO/IEA logic as Sma3s. Useful when TrEMBL is included.")
    p.add_argument("-uniprot", action="store_true",
                   help="Reference is UniProt rather than UniRef. This enables MMseqs2 clustering inside annotator 3, mimicking the original blastclust branch.")

    # MMseqs2 search options.
    p.add_argument("-filter", action="store_true",
                   help="Enable MMseqs2 low-complexity masking (--mask 1). By default masking is disabled (--mask 0), matching the old default behavior.")
    p.add_argument("--max-seqs", type=int, default=250,
                   help="Maximum number of MMseqs2 hits retained per query in the initial search. Lower values reduce disk/RAM and speed up annotators 2/3. Default: 250.")
    p.add_argument("--sensitivity", type=float, default=7.5,
                   help="MMseqs2 sensitivity (-s). Higher is more sensitive and slower/heavier. Typical values: 5.7 fast, 7.5 balanced, 8.0-8.5 very sensitive. Default: 7.5.")
    p.add_argument("--split-memory-limit", default=None,
                   help="Memory limit for large MMseqs2 prefilter structures, e.g. 40G, 80G, 180G or 1T. Reduces OOM risk by splitting the database but can slow searches.")
    p.add_argument("--db-load-mode", type=int, default=None, choices=[0, 1, 2],
                   help="Optional MMseqs2 --db-load-mode. The script does not force it by default. Use only if you understand your filesystem/RAM behavior.")
    p.add_argument("--mmseqs-bin", default="mmseqs",
                   help="MMseqs2 executable name or full path. Default: mmseqs.")
    p.add_argument("--force-search", action="store_true",
                   help="Re-run MMseqs2 even if the tabular output file already exists and is non-empty.")

    # Parallelism and reciprocal-search construction.
    p.add_argument("-num_threads", dest="num_threads", type=int, default=1,
                   help="Number of MMseqs2 threads. Annotation workers default to one quarter of this value. Default: 1.")
    p.add_argument("--annotation-workers", type=int, default=None,
                   help="Number of Python worker processes for final per-query annotation. Default: max(1, -num_threads // 4).")
    p.add_argument("--cluster-threads", type=int, default=None,
                   help="Threads per MMseqs2 clustering job inside annotator 3. Default: max(1, -num_threads // annotation-workers).")
    p.add_argument("--decompression-threads", type=int, default=None,
                   help="Threads used by rapidgzip when -d is a .dat.gz file. Default: -num_threads. The standard gzip fallback is single-threaded.")
    p.add_argument("--reciprocal-fasta-mode", choices=["auto", "index", "stream"], default="auto",
                   help="How to build the target FASTA for annotator-2 reciprocal search: index=random offset seeks; stream=scan FASTA once; auto chooses stream for large candidate sets. Default: auto.")
    p.add_argument("--reciprocal-stream-threshold", type=int, default=50000,
                   help="Number of reciprocal target sequences above which auto mode switches from random offset reads to a sequential FASTA scan. Default: 50000.")

    # Temporary files and cleanup policy.
    p.add_argument("--tmpdir", default="./sma3s_mmseqs_tmp",
                   help="Temporary directory for MMseqs2 and intermediate files. Prefer node-local scratch/NVMe when possible. Default: ./sma3s_mmseqs_tmp.")
    p.add_argument("--keep-tmp", action="store_true",
                   help="Keep MMseqs2 temporary directories instead of deleting them after successful steps.")
    p.add_argument("--force-clean", action="store_true",
                   help="Before running, remove both current-run outputs and DB-derived files/caches. Never deletes the original .dat/.dat.gz or a FASTA supplied directly with -d.")
    p.add_argument("--force-clean-outputs", action="store_true",
                   help="Before running, remove only current-run outputs/caches: TSV, summary, .mmseqs.m8, hits SQLite, reciprocal SQLite and normalized query FASTA.")
    p.add_argument("--force-clean-db", action="store_true",
                   help="Before running, remove FASTA/.annot files derived from a UniProt .dat/.dat.gz and their SQLite caches. For .dat.gz input, also removes the decompressed .dat but never the original archive.")
    p.add_argument("--clean-only", action="store_true",
                   help="Perform the requested cleanup and exit without running MMseqs2 or annotation.")
    return p

def main(argv: Optional[Sequence[str]] = None) -> int:
    start_time = time.perf_counter()
    args = build_parser().parse_args(argv)
    if args.check_install:
        return check_installation(args.mmseqs_bin, Path(args.tmpdir))
    cfg = validate_and_prepare(args)

    if cfg.force_clean_outputs or cfg.force_clean_db:
        run_requested_cleanup(cfg)
        if cfg.clean_only:
            elapsed = time.perf_counter() - start_time
            print(f"Total runtime: {format_elapsed(elapsed)} ({elapsed:.2f} seconds)\n")
            return 0

    check_mmseqs(cfg.mmseqs_bin)
    ensure_decompressed_dat(cfg)

    print(f"\nStarting Sma3s Python/MMseqs2 (annotator {cfg.annotator})", flush=True)
    ref_annot = cfg.fasta_file.with_suffix(".annot")
    active_filters = []
    if cfg.order_filter: active_filters.append(f"order={cfg.order_filter}")
    if cfg.family_filter: active_filters.append(f"family={cfg.family_filter}")
    if cfg.genus_filter: active_filters.append(f"genus={cfg.genus_filter}")
    if cfg.exclude_order_filter: active_filters.append(f"exclude order={cfg.exclude_order_filter}")
    if cfg.exclude_family_filter: active_filters.append(f"exclude family={cfg.exclude_family_filter}")
    if cfg.exclude_genus_filter: active_filters.append(f"exclude genus={cfg.exclude_genus_filter}")
    if active_filters:
        print("Active taxonomic filters: " + ", ".join(active_filters), flush=True)
        print(f"Expected filtered database: {cfg.fasta_file} + {ref_annot}", flush=True)
    if has_tax_exclusion(cfg):
        print("The final TSV will include the hit used for annotation and its identity/coverage metrics", flush=True)
    elif cfg.report_hit_metrics:
        print("Hit-metric output is enabled: the final TSV will include identity/coverage for the annotation hit", flush=True)
    if (not cfg.fasta_file.exists() or not ref_annot.exists()) and cfg.dat_file.exists() and cfg.dat_file.suffix == ".dat":
        create_fasta_annot_from_uniprot(
            cfg.dat_file,
            cfg.fasta_file,
            ref_annot,
            cfg.nopred,
            cfg.annotation_workers,
            cfg.genus_filter,
            cfg.exclude_genus_filter,
            cfg.family_filter,
            cfg.exclude_family_filter,
            cfg.order_filter,
            cfg.exclude_order_filter,
        )
    if not cfg.fasta_file.exists() or not ref_annot.exists():
        die(f"Reference FASTA or .annot file is missing: {cfg.fasta_file}, {ref_annot}")

    query_norm = Path(str(cfg.query_file) + ".sma3s_ids.fasta")
    mapping, lenq, n_queries = normalize_query_fasta(cfg.query_file, query_norm, cfg.nucl)


    # First all-vs-reference search.
    run_mmseqs_easy_search(query_norm, cfg.fasta_file, cfg.mmseqs_file, cfg, cfg.max_seqs, cfg.force_search)
    hits_db = load_hits_to_sqlite(cfg.mmseqs_file, cfg, mapping)
    hits_store = HitsStore(hits_db)
    fasta_store = FastaStore(cfg.fasta_file)
    ann_store = AnnotationStore(ref_annot, cfg)

    if "2" in cfg.annotator:
        reciprocal_targets = hits_store.reciprocal_candidate_targets(cfg, lenq)
    else:
        reciprocal_targets = []
    reciprocal_db = build_reciprocal_best_db(cfg, query_norm, fasta_store, reciprocal_targets)
    reciprocal_best = ReciprocalBestStore(reciprocal_db)

    writer = OutputWriter(cfg)
    try:
        print(
            f"Anotando sequences con {cfg.annotation_workers} worker(s) Python "
            f"(MMseqs2 usa {cfg.cpus} thread(s); clustering A3 usa {cfg.cluster_threads} thread(s) por worker)",
            flush=True,
        )
        if cfg.annotation_workers <= 1:
            for i in range(1, n_queries + 1):
                q = f"s{i}"
                hits = hits_store.get_hits(q)
                annotate_query(q, hits, cfg, mapping, lenq, ann_store, fasta_store, reciprocal_best, writer)
                if i % 10000 == 0:
                    print(f"  processed {i}/{n_queries}", flush=True)
        else:
            # Close parent handles before forking/spawning worker processes. Each worker opens
            # independent read-only SQLite/FASTA handles. Results are merged in original FASTA order.
            hits_store.close()
            fasta_store.close()
            ann_store.close()
            reciprocal_best.close()
            chunksize = max(1, min(1000, n_queries // (cfg.annotation_workers * 8) if n_queries else 1))
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=cfg.annotation_workers,
                initializer=init_annotation_worker,
                initargs=(
                    cfg,
                    str(hits_db),
                    str(cfg.fasta_file),
                    str(ref_annot),
                    str(reciprocal_db) if reciprocal_db else None,
                    mapping,
                    lenq,
                ),
            ) as executor:
                for i, result in enumerate(executor.map(annotate_query_worker, range(1, n_queries + 1), chunksize=chunksize), 1):
                    merge_annotation_result(writer, result)
                    if i % 10000 == 0:
                        print(f"  processed {i}/{n_queries}", flush=True)
    finally:
        writer.close()
        hits_store.close()
        fasta_store.close()
        ann_store.close()
        reciprocal_best.close()
        if not cfg.keep_tmp:
            with contextlib.suppress(Exception):
                query_norm.unlink()

    stat_file = writer.write_summary(n_queries)
    elapsed = time.perf_counter() - start_time
    print(f"\nAnnotation file created: {cfg.annot_file}")
    print(f"Summary file created: {stat_file}")
    print(f"Total runtime: {format_elapsed(elapsed)} ({elapsed:.2f} seconds)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
