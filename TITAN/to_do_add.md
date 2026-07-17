# Audit d'ajouts d'outils pour TITAN

Date de l'audit web: 2026-07-16

Objectif: proposer des outils d'annotation structurale ou fonctionnelle a integrer dans TITAN, adaptes au graphe actuel. Les nouveaux outils de prediction structurale ci-dessous doivent d'abord etre publies comme annotations additionnelles independantes, sans etre injectes dans le merge AEGIS. Une option future pourra permettre leur inclusion dans AEGIS apres validation biologique.

## Sources principales consultees

- Mikado: https://mikado.readthedocs.io/
- Helixer: https://github.com/usadellab/Helixer
- eggNOG-mapper: https://github.com/eggnogdb/eggnog-mapper et https://github.com/eggnogdb/eggnog-mapper/blob/main/USAGE.md
- TransDecoder: https://github.com/TransDecoder/TransDecoder
- FLAIR: https://github.com/BrooksLabUCSC/flair
- SQANTI3: https://github.com/ConesaLab/SQANTI3

Les tags BioContainers ci-dessous ont ete verifies via l'API Quay le 2026-07-16. Avant implementation, pinner les images par digest comme le reste du projet.

## Decision d'architecture

### Regle pour les annotations additionnelles

- Nouveau dossier public: `${params.output_dir}/additional_annotations/<tool>/`.
- Connecter ces sorties a `subworkflows/aegis.nf` quand validé.
- Ajouter une provenance dediee dans `titan_provenance` ou un nouveau process `additional_annotations_provenance`.
- Chaque module doit emettre `versions.yml`.
- Chaque outil doit avoir un parametre `params.run_<tool> = false` par defaut, sauf eggNOG-mapper qui peut etre active par defaut si la base eggNOG locale est disponible.

## P0 - Annotation fonctionnelle apres AEGIS

### `modules/eggnog_mapper.nf`

- But: ajouter une annotation fonctionnelle orthology-based en complement de Diamond2GO sur les proteines finales AEGIS.
- Position dans le graphe: dans `subworkflows/aegis.nf`, apres `aegis_short_reads` ou `aegis_long_reads`, en parallele de `diamond2go`.
- Input principal: `aegis.out.aegis_proteins_main` et optionnellement `aegis.out.aegis_gff`.
- Output public: `${params.output_dir}/EggNOG_outputs/`.
- Image candidate: `quay.io/biocontainers/eggnog-mapper:2.1.15--pyhdfd78af_0`.
- Parametres a ajouter:
  - `params.container_eggnog_mapper`
  - `params.eggnog_data_dir`
  - `params.eggnog_mapper_cpus`
  - `params.eggnog_mapper_tax_scope = false`
  - `params.eggnog_mapper_sensmode = 'sensitive'`
  - `params.run_eggnog_mapper = true`
- Commande type:

```bash
emapper.py \
  -i final_annotation_proteins_main.fasta \
  --itype proteins \
  -m diamond \
  --cpu ${task.cpus} \
  --data_dir ${eggnog_data_dir} \
  --output titan_eggnog_main \
  --output_dir . \
  --excel \
  --report_orthologs \
  --sensmode ${params.eggnog_mapper_sensmode}
```

- Sorties attendues:
  - `titan_eggnog_main.emapper.annotations`
  - `titan_eggnog_main.emapper.seed_orthologs`
  - `titan_eggnog_main.emapper.orthologs`
- Notes d'implementation:
  - Le mode `--decorate_gff` existe dans eggNOG-mapper, mais il faut le tester avec le GFF3 AEGIS. Au premier ajout, produire seulement le TSV proteique et ajouter un module ulterieur de decoration GFF3 si les IDs protein/GFF sont compatibles.
  - La base eggNOG ne doit pas etre telechargee dans chaque workdir. Exiger `--eggnog_data_dir` ou ajouter un module separe `download_eggnog_data` explicitement cache.
  - Ajouter `withName: eggnog_mapper { cpus = { params.eggnog_mapper_cpus as int } }`.

## P1 - Annotations RNA-seq additionnelles sans merge AEGIS

### `modules/mikado_prepare.nf`, `modules/mikado_serialise.nf`, `modules/mikado_pick.nf`

- But: selectionner les meilleurs transcrits parmi StringTie, PsiCLASS et long-read GTFs.
- Pourquoi: Mikado consolide plusieurs assemblies transcriptomiques et selectionne des modeles par locus; c'est plus adapte qu'un simple merge de GTFs pour produire une annotation RNA-seq candidate.
- Position: apres les merges STAR/StringTie, STAR/PsiCLASS, long reads et Scallop2.
- Inputs:
  - genome FASTA: `params.new_assembly`
  - liste GTF/GFF transcriptomiques deja produits
  - proteines externes optionnelles depuis `protein_samplesheet`
- Output public: `${params.output_dir}/additional_annotations/mikado/`.
- Image candidate: `quay.io/biocontainers/mikado:2.3.4--py310h8ea774a_2` plutot qu'un tag `rc`.
- Parametres:
  - `params.run_mikado = false`
  - `params.container_mikado`
  - `params.mikado_mode = 'permissive'`
- Commandes types:

```bash
mikado configure \
  --list transcript_inputs.tsv \
  --reference ${genome} \
  --mode ${params.mikado_mode} \
  --scoring plants.yaml \
  configuration.yaml

mikado prepare \
  --json-conf configuration.yaml

mikado serialise \
  --json-conf configuration.yaml \
  --orfs mikado_prepared.fasta.transdecoder.bed \
  --blast_targets protein_hits.xml

mikado pick \
  --json-conf configuration.yaml \
  --subloci-out mikado.subloci.gff3 \
  --loci-out mikado.loci.gff3
```

- Sorties:
  - `mikado.loci.gff3`
  - `mikado.subloci.gff3`
  - `mikado_prepared.fasta`
- Notes:
  - Mikado est interessant seulement si on ajoute TransDecoder et/ou des hits proteiques. Sinon il reste moins informatif.
  - Prevoir un petit generateur `transcript_inputs.tsv` dans `scripts/` plutot qu'une closure Groovy complexe.

### `modules/transdecoder_longorfs.nf` et `modules/transdecoder_predict.nf`

- But: predire les regions codantes dans les transcrits RNA-seq consolides par Mikado ou PASA.
- Pourquoi: TransDecoder transforme des transcriptomes en CDS/proteines candidates; utile pour comparer les proteines AEGIS avec une evidence RNA-seq non mergee.
- Position: apres `mikado_prepare` ou apres une extraction FASTA depuis GTF.
- Input: FASTA de transcrits.
- Output public: `${params.output_dir}/additional_annotations/transdecoder/`.
- Image candidate: `quay.io/biocontainers/transdecoder:6.0.0--pl5321hdfd78af_0`.
- Parametres:
  - `params.run_transdecoder = false`
  - `params.container_transdecoder`
- Commandes types:

```bash
TransDecoder.LongOrfs \
  -t transcripts.fasta

TransDecoder.Predict \
  -t transcripts.fasta
```

- Sorties:
  - `transcripts.fasta.transdecoder.pep`
  - `transcripts.fasta.transdecoder.cds`
  - `transcripts.fasta.transdecoder.gff3`
- Notes:
  - Pour un resultat plus robuste, ajouter plus tard `--retain_blastp_hits` et `--retain_pfam_hits`, mais cela implique BLAST/Pfam/HMMER et des bases externes.

### `modules/flair_isoforms.nf`

- But: produire une annotation isoforme long-read alternative a StringTie long reads.
- Pourquoi: FLAIR est specialise pour correction, definition d'isoformes et analyses de splicing sur reads longs PacBio/ONT.
- Position: branche long reads, apres preparation FASTQ/FASTA long reads et avec genome FASTA.
- Inputs:
  - long-read FASTQ/FASTA par sample
  - genome FASTA
  - annotation de reference optionnelle: Liftoff ou AEGIS si execute en post-hoc
- Output public: `${params.output_dir}/additional_annotations/flair/`.
- Image candidate: `quay.io/biocontainers/flair:3.0.0--pyhdfd78af_0`.
- Parametres:
  - `params.run_flair = false`
  - `params.container_flair`
- Commande type:

```bash
flair align \
  -g ${genome} \
  -r ${reads} \
  -o ${sample_ID}.flair

flair correct \
  -q ${sample_ID}.flair.bed \
  -g ${genome} \
  -f ${reference_annotation} \
  -o ${sample_ID}.flair.corrected

flair collapse \
  -g ${genome} \
  -r ${reads} \
  -q ${sample_ID}.flair.corrected.bed \
  -o ${sample_ID}.flair
```

- Sorties:
  - `${sample_ID}.flair.isoforms.gtf`
  - `${sample_ID}.flair.isoforms.fa`
- Notes:
  - FLAIR est pertinent uniquement si `has_long_reads == true`.
  - Au premier ajout, utiliser Liftoff comme reference annotation pour ne pas creer de dependance circulaire avec AEGIS.

### `modules/sqanti3_qc.nf`

- But: QC et classification des transcriptomes long-read produits par FLAIR ou StringTie2.
- Pourquoi: SQANTI3 caracterise et cure les modeles de transcrits long-read, utile avant de faire confiance aux isoformes.
- Position: apres FLAIR et/ou merge long reads.
- Inputs:
  - isoforms GTF/FASTA
  - genome FASTA
  - reference annotation GFF3/GTF, idealement Liftoff ou AEGIS en mode post-hoc.
- Output public: `${params.output_dir}/additional_annotations/sqanti3/`.
- Image candidate: `quay.io/biocontainers/sqanti3:6.0.1--hdfd78af_0`.
- Parametres:
  - `params.run_sqanti3 = false`
  - `params.container_sqanti3`
- Commande type:

```bash
sqanti3_qc.py \
  ${isoforms_gtf} \
  ${reference_annotation} \
  ${genome} \
  --dir . \
  --output ${sample_ID}.sqanti3
```

- Sorties:
  - `${sample_ID}.sqanti3_classification.txt`
  - `${sample_ID}.sqanti3_corrected.gtf`
  - rapports QC SQANTI3
- Notes:
  - Integrer d'abord en QC/reporting, pas comme evidence AEGIS.

## P1 - Prediction ab initio complementaire sans merge AEGIS

### `modules/helixer_prediction.nf`

- But: produire une prediction ab initio/deep learning independante des evidences RNA-seq/proteines.
- Pourquoi: Helixer predit des modeles geniques directement depuis le genome; utile comme controle orthogonal a BRAKER3.
- Position: apres validation des inputs et avant AEGIS, en branche independante.
- Input: `params.new_assembly` ou genome EDTA masque selon test empirique.
- Output public: `${params.output_dir}/additional_annotations/helixer/`.
- Image candidate:
  - Pas de tag public BioContainers accessible via Quay API pendant l'audit.
  - Prevoir `dockerfiles/helixer/Dockerfile` ou utiliser une image officielle Helixer/Apptainer validee localement.
- Parametres:
  - `params.run_helixer = false`
  - `params.container_helixer`
  - `params.helixer_model = 'land_plant'`
  - `params.helixer_use_gpu = false`
- Commande type a valider selon l'image choisie:

```bash
Helixer.py \
  --fasta-path ${genome} \
  --lineage ${params.helixer_model} \
  --gff-output-path helixer.gff3 \
  --temporary-dir helixer_tmp
```

- Sorties:
  - `helixer.gff3`
  - logs et metadata modele
- Notes:
  - Pour plantes, tester explicitement le modele land plant sur les fixtures realistes.
  - GPU peut accelerer fortement, mais ne doit pas etre obligatoire pour le profil `test`.

## Changements de configuration a prevoir

Ajouter dans `nextflow.config` apres validation:

```groovy
params {
  run_eggnog_mapper = true
  eggnog_data_dir = false
  eggnog_mapper_cpus = 8
  container_eggnog_mapper = "quay.io/biocontainers/eggnog-mapper:2.1.15--pyhdfd78af_0"

  run_mikado = false
  container_mikado = "quay.io/biocontainers/mikado:2.3.4--py310h8ea774a_2"

  run_transdecoder = false
  container_transdecoder = "quay.io/biocontainers/transdecoder:6.0.0--pl5321hdfd78af_0"

  run_flair = false
  container_flair = "quay.io/biocontainers/flair:3.0.0--pyhdfd78af_0"

  run_sqanti3 = false
  container_sqanti3 = "quay.io/biocontainers/sqanti3:6.0.1--hdfd78af_0"

  run_helixer = false
  container_helixer = false
  helixer_model = "land_plant"
}
```

Ajouter dans `conf/base.config`:

```groovy
process {
  withName: eggnog_mapper {
    cpus = { params.eggnog_mapper_cpus as int }
    memory = 64.GB
    time = 24.h
  }
}
```

## Ordre d'implementation recommande

1. Ajouter `eggnog_mapper` apres AEGIS et Diamond2GO, avec test stub et validation de `--eggnog_data_dir`.
2. Ajouter `mikado_prepare/serialise/pick` et `transdecoder` pour produire une annotation RNA-seq candidate complete.
3. Ajouter `flair_isoforms` puis `sqanti3_qc` pour exploiter les long reads au-dela de StringTie.
4. Ajouter `helixer_prediction` comme prediction ab initio orthogonale.
