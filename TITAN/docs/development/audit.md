# Audit initial TITAN

Date: 2026-07-16
Racine TITAN: `/home/vmadmin/workflows/TITAN`
Depot Git: `/home/vmadmin/workflows`
Branche de travail: `codex/titan-hardening`

## Etat Git initial

Commandes executees:

```bash
git -C /home/vmadmin/workflows status --short --branch
git -C /home/vmadmin/workflows branch --show-current
git -C /home/vmadmin/workflows remote -v
git -C /home/vmadmin/workflows log --oneline -n 10
git status --short -- TITAN
```

Resultat synthetique: branche initiale `main`, remote `origin git@github-amandine:Grapedia/workflows.git`, aucune modification non committee detectee dans `TITAN/`. Branche creee: `codex/titan-hardening`.

## Architecture actuelle

Point d'entree principal: `main.nf`. Depuis P1-002, `main.nf` inclut et appelle le workflow `TITAN` defini dans `workflows/titan.nf`.

Sous-workflows:

* `subworkflows/generate_evidence_data.nf`: genere les evidences transcriptomiques, ab initio, Liftoff et EDTA.
* `subworkflows/aegis.nf`: integre les evidences via Aegis puis lance Diamond2GO.

Modules detectes: 27 process Nextflow dans `modules/*.nf`.

Scripts auxiliaires actifs: wrappers StringTie/EDTA, nettoyage FASTA proteique BRAKER3 et validation des fixtures minimales. Les anciens scripts AEGIS, scripts `retrieve_*` et scripts monoexons non references ont ete supprimes apres P1-006; les modules AEGIS utilisent le CLI upstream.

Donnees exemple: `data_example/` contient des placeholders tres legers, majoritairement de taille 0; ce repertoire ne constitue pas un jeu de test executable fiable.

## Commandes historiques

Documentees dans `README.md` et `launch_TITAN_example.sh`:

```bash
module load nextflow/24.04.3
nextflow run main.nf -with-dag dag_titan.png
```

Le script historique contient un chemin absolu specifique a une machine: `/home/avelt/data2/.../workflows/TITAN`.

## Baseline executee

| Test | Commande | Resultat | Duree approximative | Remarque |
| --- | --- | --- | --- | --- |
| Version Nextflow | `nextflow -version` | Succes | <1 s | Nextflow 26.04.3 disponible; TITAN declare 24.04.3. |
| Configuration | `nextflow config` | Succes | ~2 s | Les parametres par defaut resolvent. |
| Lancement minimal historique initial | `nextflow run main.nf --workflow aegis -ansi-log false` | Echec | ~3 s | Erreur de compilation DSL2: statements top-level melanges avec declarations. |
| Profils runtime | `python3 scripts/validate_profiles.py` | Succes | ~6 s | Verifie local, test, Apptainer et Slurm sans soumission Slurm. |
| Lancement TITAN stub apres contrat unique | `nextflow run main.nf -profile test -stub-run -ansi-log false` | Succes | ~5 s | Compile et lance generation d'evidences puis Aegis dans le meme graphe. |
| Profil test P0-003 | `nextflow config -profile test`; `python3 scripts/validate_minimal_test_data.py` | Succes | ~3 s | Profil local, sans Slurm ni Docker, pointe vers `test-data/minimal/valid`, ecrit dans `test-results/` et `test-work/`. |

Erreur exacte principale:

```text
Error main.nf:11:1: Statements cannot be mixed with script declarations -- move statements into a process, workflow, or function
```

Corrections P0 appliquees:

* validation des parametres et creation des channels de `generate_evidence_data` deplacees dans le bloc `workflow`;
* construction des channels RNA/proteines dans le workflow TITAN unique;
* detection des long reads depuis `library_layout=long` dans la samplesheet RNA-seq;
* suppression des flags biologiques `EDTA`, `run_edta`, `run_egapx` et `use_long_reads` du contrat runtime;
* suppression du routage public Aegis-only par fichiers publies;
* remplacement d'une boucle `for` non supportee par Nextflow 26 dans `subworkflows/aegis.nf`;
* correction du tag `diamond2go` qui referenceait une variable inexistante;
* erreur explicite lorsque les evidences Aegis obligatoires manquent.

## Fonctionnement biologique compris

TITAN annote un nouvel assemblage de genome en combinant:

* transfert d'annotations precedentes via Liftoff;
* annotation EGAPx obligatoire et emise sous forme de sorties nommees;
* conversion AGAT GFF3 vers CDS FASTA pour inference de brin;
* preparation RNA-seq court/long depuis FASTQ/FASTA/SRA;
* trimming fastp;
* alignements STAR, HISAT2 et Minimap2;
* assemblages transcriptomiques StringTie et PsiCLASS;
* fusion GFFCompare et StringTie;
* masquage EDTA;
* prediction BRAKER3 avec ou sans long reads;
* integration Aegis chromosome par chromosome;
* annotation fonctionnelle Diamond2GO.

## Matrice comparative atcg-rnaseq / TITAN

| Domaine | atcg-rnaseq | TITAN actuel | Ecart | Adaptation recommandee |
| --- | --- | --- | --- | --- |
| Architecture Nextflow | `main.nf` leger + `workflows/`, `subworkflows/`, `modules/local/` | `main.nf` est leger; `workflows/titan.nf` porte l'orchestration unique et connecte les evidences nommees a Aegis | Modules encore trop couples aux dossiers publies | Migrer progressivement vers des inputs `path` stages et des emits nommes. |
| DSL2 | DSL2 structure valide | DSL2 active mais compilation cassee sous Nextflow 26 | Bloquant | Deplacer statements top-level dans `workflow`. |
| Point d'entree | Help et validation parametres | Pas de `--help`; validation top-level | UX et compilation | Ajouter help non executant et validation interne. |
| Parametres | YAML utilisateur + schema | Parametres dans `nextflow.config` | Peu portable | Introduire config YAML/schema sans casser les anciens `--param`. |
| Validation entrees | Scripts Python testes | `scripts/validate_inputs.py` integre au workflow | Partiel | Etendre vers schema nf-validation si besoin. |
| Modules | Process atomiques labels | Modules nombreux mais scripts/volumes couples | Fragile | Clarifier contrats et labels module par module. |
| Sous-workflows | Fonctions coherentes | Deux sous-workflows biologiques | Base utile | Garder mais typer les emits. |
| Dependances | Versions documentees et images | Images runtime centralisees et digest-pinned | Partiel | Etendre les `versions.yml` module par module. |
| Conteneurs | Apptainer profile | Docker active globalement | Inadapte cluster | Ajouter profil Apptainer et garder Docker optionnel. |
| Profils | local/test/slurm/apptainer | Aucun profil explicite | Bloquant tests | Ajouter `conf/*.config`. |
| Ressources | Labels centralises | CPU dans modules + config globale | Peu maintenable | Centraliser par labels. |
| Donnees test | Synthetic fixtures versionnees | `data_example` placeholders | Non executable | Creer `test-data/minimal`. |
| Tests | pytest, nf-test, scripts | Aucun test detecte | Risque regression | Ajouter tests statiques puis nf-test. |
| CI | Scripts et workflows | Pas de CI TITAN detectee | Non verifie | CI rapide sans gros calcul. |
| Logs | Rapports structures | Echo dans scripts, logs implicites | Peu exploitable | Standardiser logs par process. |
| Rapports | MultiQC, trace, timeline, DAG | DAG possible manuellement | Incomplet | Documenter `-with-*`. |
| Documentation | README + docs | README utile mais incomplet | Manque exploitation | Ajouter docs dev/utilisateur. |
| Reproductibilite | versions.yml, provenance | Versions README + EGAPx/AEGIS/provenance | Partiel | Etendre versions aux modules restants. |
| Provenance | Manifestes et checksums | `provenance/evidence_manifest.json` | Partiel | Etendre le manifeste aux outputs intermediaires si necessaire. |
| Gestion erreurs | Statuts structures | Warnings et erreurs parfois non bloquants | Risque silencieux | Echec explicite sur invalides. |
| Reprise | Tests `-resume` | Compatibilite non verifiee | Inconnu | Tester apres baseline. |
| Outputs | Contrat documente | Sorties publiques declarees + `saveAs` | Partiel | Stabiliser les modules restants sous le meme contrat. |
| Versionnement | Changelog, version | Manifest version 1.0 | Minimal | Versionner pipeline et outils. |

Bonnes pratiques reutilisables: separation profils, fixtures minimales, validation precoce, scripts testables, provenance, CI rapide, contrats de sorties. Elements specifiques RNA-seq non reutilisables directement: ENA accessions, matrices featureCounts/TPM/z-score, logique HISAT2 expression.

## Inventaire des outils et conteneurs

Commande d'inspection statique utilisee pour P0-002:

```bash
rg -n "^(process|workflow|include|\\s*container|\\s*publishDir|\\s*input:|\\s*output:)" main.nf subworkflows modules nextflow.config conf
for f in modules/*.nf; do awk '/^[[:space:]]*input:/{flag=1} /^[[:space:]]*script:|^[[:space:]]*shell:|^[[:space:]]*stub:/{flag=0} flag{print}' "$f"; done
for f in modules/*.nf; do awk '/^[[:space:]]*output:/{flag=1} /^[[:space:]]*script:|^[[:space:]]*shell:|^[[:space:]]*stub:/{flag=0} flag{print}' "$f"; done
```

| Outil | Etape | Version actuelle | Version verrouillee | Module | Conteneur | Test | Version collectee |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Aegis | Integration finale annotation | v0.3.25 image label | Digest-pinned | `aegis_short_reads`, `aegis_long_reads` | `params.container_aegis` | Stub + test CLI Docker | `versions.yml` |
| AGAT | GFF3 vers CDS FASTA | 1.2.0 | Digest-pinned | `agat_convert_gff3_to_cds_fasta` | `params.container_agat` | Stub | Non |
| BRAKER3 | Prediction ab initio | v3.0.8 documente | Digest-pinned | `braker3_prediction*` | `params.container_braker3` | Stub | `versions.yml` |
| Diamond2GO | Annotation fonctionnelle | commit documente | Digest-pinned | `diamond2go` | `params.container_diamond2go` | Stub | Non |
| EDTA | Masquage TE | GitHub latest documente | Digest-pinned | `EDTA` | `params.container_edta` | Stub | Non |
| EGAPx | Annotation NCBI | 0.5.2 | Digest-pinned | `egapx` | `params.container_egapx` | Stub | `versions.yml` |
| fastp | Trimming | 0.23.2 | Digest-pinned | `trimming_fastq` | `params.container_fastp` | Stub | Non |
| GFFCompare | Fusion PsiCLASS | 0.12.6 | Digest-pinned | `gffcompare` | `params.container_gffcompare` | Stub | `versions.yml` |
| HISAT2 | Index/alignment | 2.2.1 documente | Digest-pinned | `hisat2_*` | `params.container_hisat2` | Stub | Non |
| Liftoff | Transfert annotation | 1.5.1 | Digest-pinned | `liftoff_annotations` | `params.container_liftoff` | Stub | Non |
| Minimap2 | Alignement long reads | 2.28 documente | Digest-pinned | `minimap2_*` | `params.container_minimap2_samtools` | Stub | Non |
| PsiCLASS | Assemblage transcriptome | 1.0.2 | Digest-pinned | `assembly_transcriptome_star_psiclass` | `params.container_psiclass_samtools` | Stub | `versions.yml` |
| Salmon | Inference strandedness | 1.10.3 | Digest-pinned | `salmon_*` | `params.container_salmon` | Stub | Non |
| ENA API downloader | SRA/ENA FASTQ staging | Python stdlib | Digest-pinned | `prepare_RNAseq_fastq_files_*`, `scripts/download_sra_fastq.py` | `params.container_python` | Unit + Stub | Non |
| STAR | Index/alignment | 2.7.11b | Digest-pinned | `star_*` | `params.container_star` | Stub | Non |
| StringTie | Assemblage/fusion | 2.2.3 documente | Digest-pinned | `assembly_transcriptome_*`, `Stringtie_merging_*` | `params.container_stringtie` | Stub | Non |
| gffread | Conversion GTF/GFF3 | Non requis par les modules AEGIS actuels | Non | n/a | n/a | n/a | n/a |
| samtools | BAM | 1.9 documente selon images | Partiel | alignements | images combinees | Non | Non |
| DIAMOND | Diamond2GO/BRAKER3 | 2.1.9/2.1.11 selon docs | Partiel | `diamond2go`, `braker3` | images combinees | Non | Non |

### Inventaire des entrees

| Entree | Parametre/channel | Format attendu | Colonnes ou structure | Etapes consommatrices | Remarques |
| --- | --- | --- | --- | --- | --- |
| Assemblage cible | `params.new_assembly` | FASTA | fichier existant | Liftoff, AGAT, STAR, HISAT2, Minimap2, EDTA, BRAKER3, Aegis | Chemin separe en dossier + nom dans les modules. |
| Assemblage precedent | `params.previous_assembly` | FASTA | fichier existant | Liftoff | Doit correspondre a `previous_annotations`. |
| Annotation precedente | `params.previous_annotations` | GFF3 | fichier existant | Liftoff, puis AGAT | Source du transfert et de la FASTA CDS pour Salmon. |
| Samplesheet RNA-seq | `params.RNAseq_samplesheet` | CSV avec header | `sampleID`, `SRA_or_FASTQ`, `library_layout` | preparation reads, fastp, Salmon, STAR, HISAT2, Minimap2 | `library_layout` alimente les branches `single`, `paired` et `long`; la presence d'une ligne `long` active automatiquement la branche longue. |
| Samplesheet proteines | `params.protein_samplesheet` | CSV avec header | `organism`, `filename` | BRAKER3 | Les modules BRAKER3 montent aussi `${projectDir}/data/protein_data`; ce couplage reste a normaliser. |
| Parametres EGAPx | `params.egapx_paramfile` | YAML | parametres EGAPx | module `egapx` | Obligatoire; sorties nommees publiees sous `${output_dir}/egapx`. |
| Options outil | `params.edta_cpus`, `params.egapx_cpus`, `params.PSICLASS_*`, `params.STAR_genomeSAindexNbases`, `params.STAR_sjdbGTFfile` | entiers CPU, floats, options STAR | valeurs scalaires | EDTA, EGAPx, PsiCLASS, STAR | Les anciens flags biologiques ne pilotent plus EDTA, EGAPx ou les long reads. |
| Evidences Aegis | emits nommes de `generate_evidence_data` | chemins stages par Nextflow | `masked_genome`, `liftoff_annotation`, `egapx_gff3`, `braker_augustus_gff`, `braker_genemark_gtf`, STAR/StringTie, STAR/PsiCLASS, long reads si detectes | sous-workflow `aegis` | Le workflow public ne relit plus ces evidences depuis `output_dir`. |

### Inventaire des sorties par module

| Module/process | Entrees principales | Sorties declarees | Publication actuelle | Consommateur principal |
| --- | --- | --- | --- | --- |
| `prepare_RNAseq_fastq_files_short` | tuples RNA-seq courts, script ENA | tuples `sample_ID`, `SRA_or_FASTQ`, `library_layout`, `prepared_1.fastq.gz`, `prepared_2.fastq.gz` | workdir uniquement | `trimming_fastq` |
| `prepare_RNAseq_fastq_files_long` | tuples RNA-seq longs, script ENA | tuples `sample_ID`, `SRA_or_FASTQ`, `library_layout`, `read_format`, `long_read_input.fastq.gz`, `long_read_input.fasta` | workdir uniquement | `minimap2_alignment` |
| `trimming_fastq` | FASTQ/SRA courts prepares | `${sample_ID}_1.trimmed.fastq.gz`, `${sample_ID}_2.trimmed.fastq.gz`, rapports fastp JSON/HTML | `${output_dir}/intermediate_files/evidence_data/RNAseq_data/trimmed_data` | Salmon, STAR, HISAT2 |
| `liftoff_annotations` | assemblage cible, assemblage precedent, GFF3 precedent | `liftoff_previous_annotations.gff3`, `unmapped_features.txt` | `${output_dir}` | AGAT, Aegis |
| `agat_convert_gff3_to_cds_fasta` | assemblage cible, Liftoff GFF3 | `${genome}.CDS.fasta.gz` | `${output_dir}/intermediate_files/liftoff/gff3_to_cds_fasta` | Salmon index |
| `salmon_index` | CDS FASTA | `salmon_index/` | workdir/cache uniquement | Salmon strand inference |
| `salmon_strand_inference` | reads trimmes, index Salmon | strand info en valeur, `${sample_ID}.strand_info.classified`, `${sample_ID}.log` | `${output_dir}/intermediate_files/salmon_strand` si intermediaires publies | STAR/HISAT2 channel enrichi avec strandedness |
| `star_genome_indices` | assemblage cible | `${genome}_index` | `${output_dir}/intermediate_files/evidence_data/star_databases` | STAR alignment |
| `star_alignment` | index STAR, reads, strandedness | `${sample_ID}_Aligned.sortedByCoord.out.bam` | `${output_dir}/intermediate_files/evidence_data/RNAseq_alignments/STAR/{stranded,unstranded}` | StringTie, PsiCLASS, BRAKER3 |
| `hisat2_genome_indices` | assemblage cible | `hisat2_index/` | `${output_dir}/intermediate_files/evidence_data/hisat2_databases` | HISAT2 alignment |
| `hisat2_alignment` | index HISAT2, reads, strandedness | `${sample_ID}_Aligned.sort.bam`, `${sample_ID}_Aligned.sort.bam.bai` | `${output_dir}/intermediate_files/evidence_data/RNAseq_alignments/HISAT2/{stranded,unstranded}` | StringTie HISAT2 |
| `minimap2_genome_indices` | assemblage cible | `${genome}.mmi` | workdir uniquement | Minimap2 alignment |
| `minimap2_alignment` | index Minimap2, reads longs, format FASTQ/FASTA | `${sample_ID}_Aligned.sorted.bam` | `${output_dir}/intermediate_files/evidence_data/RNAseq_alignments/minimap2` | StringTie long reads, BRAKER3 long reads |
| `assembly_transcriptome_star_stringtie` | BAM STAR | tuple `sample_ID`, default GTF, alt GTF, `strand_type` | `${output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/STAR/{stranded,unstranded}` | StringTie STAR merge |
| `assembly_transcriptome_hisat2_stringtie` | BAM HISAT2 | tuple `sample_ID`, default GTF, alt GTF, `strand_type` | `${output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/HISAT2/{stranded,unstranded}` | StringTie HISAT2 merge |
| `assembly_transcriptome_star_psiclass` | BAM STAR | `${sample_ID}_vote.gtf` | `${output_dir}/intermediate_files/evidence_data/transcriptomes/STAR_PsiCLASS/{stranded,unstranded}` | GFFCompare |
| `assembly_transcriptome_minimap2_stringtie` | BAM Minimap2 | tuple `sample_ID`, default GTF, alt GTF | `${output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/long_reads` | StringTie long reads merge |
| `Stringtie_merging_short_reads_STAR` | listes stagees de GTF STAR/StringTie par strandedness et mode default/alt | `merged_transcriptomes.STAR.short_reads.*.gtf`, `versions.yml` | `${output_dir}/tmp` et noms historiques sous `${output_dir}` via `publishDir saveAs` | Aegis |
| `Stringtie_merging_short_reads_hisat2` | listes stagees de GTF HISAT2/StringTie par strandedness et mode default/alt | `merged_transcriptomes.hisat2.short_reads.*.gtf`, `versions.yml` | `${output_dir}/tmp` et noms historiques sous `${output_dir}` via `publishDir saveAs` | Public/provenance; non utilise par Aegis |
| `Stringtie_merging_long_reads` | listes stagees de GTF Minimap2/StringTie default/alt | `merged_transcriptomes.minimap2.long_reads.*.gtf`, `versions.yml` | `${output_dir}/tmp` et noms historiques sous `${output_dir}` via `publishDir saveAs` | Aegis long reads |
| `gffcompare` | listes stagees de GTF PsiCLASS par strandedness | `stranded_merged_output.combined.gtf`, `unstranded_merged_output.combined.gtf` vide si absent, `versions.yml` | `${output_dir}/tmp` et noms historiques sous `${output_dir}` via `publishDir saveAs` | Aegis |
| `EDTA` | assemblage cible | `*TElib.fa`, `*TEanno.gff3`, `*MAKER.masked` | `${output_dir}/tmp` plus copies script vers `${output_dir}` | Aegis |
| `braker3_prediction` | genome, proteines, BAM courts | `augustus.hints.gff3`, `genemark.gtf`, `genemark_supported.gtf`, `braker.gff3`, `versions.yml`; logs et manifeste sous `intermediate_files/braker3` si `publish_intermediates` | `${output_dir}` | Aegis |
| `braker3_prediction_with_long_reads` | genome, proteines, BAM courts + longs | `augustus.hints.gff3`, `genemark.gtf`, `genemark_supported.gtf`, `braker.gff3`, `versions.yml`; logs et manifeste sous `intermediate_files/braker3` si `publish_intermediates` | `${output_dir}` | Aegis long reads |
| `aegis_merge` | genome masque, Liftoff, EGAPx GFF3, BRAKER3, STAR/StringTie, GFFCompare, StringTie long reads si mode long | `final_annotation.gff3`, `final_annotation_proteins_all.fasta`, `final_annotation_proteins_main.fasta`, `versions.yml`; `aegis_inputs.tsv`, logs et dossiers AEGIS sous `intermediate_files/aegis` si `publish_intermediates` | `${output_dir}/aegis_outputs` | Diamond2GO |
| `diamond2go` | proteines Aegis all/main | `*-diamond*` | `${output_dir}/Diamond2GO_outputs` | sortie finale fonctionnelle |
| `egapx` | YAML EGAPx | `egapx.complete.genomic.gff3`, `egapx.complete.genomic.gtf`, `egapx.complete.proteins.faa`, `egapx.complete.cds.fna`, `egapx.complete.transcripts.fna`, `egapx.annotated_genome.asn`, `egapx_out/`, `versions.yml` | `${output_dir}/egapx` | Evidence generation; le GFF3 est transmis a `aegis merge` |

Points de vigilance P0-002:

* plusieurs modules publient via `publishDir`; les copies manuelles `/outputdir` ont ete supprimees des merges StringTie, GFFCompare et EDTA;
* les images runtime sont centralisees dans `nextflow.config` et verrouillees par digest depuis P1-008;
* EGAPx, AEGIS et TITAN provenance collectent un `versions.yml`; les autres modules restent a couvrir progressivement;
* le schema d'entrees est valide avant construction des channels depuis P2-001: FASTA, GFF3, samplesheets, chemins FASTQ/proteines, YAML EGAPx minimal et enums;
* EGAPx est documente, obligatoire, expose des sorties nommees et son GFF3 est consomme par AEGIS.

## Constats critiques

Synthese operationnelle P0: `docs/development/p0-hardening.md`.
Audit architecture/refactor: `docs/development/architecture-audit.md`.

P0:

* `main.nf` ne compilait pas avec Nextflow 26.04.3 avant correction P0; depuis le durcissement P1, TITAN compile et lance le graphe complet en stub.
* La validation P0-005 echoue explicitement avant calcul lourd pour parametre obligatoire vide et fichier d'entree absent.
* Profil `test` local ajoute et valide: resolution de configuration et commande minimale TITAN sans Slurm, Docker ni donnees volumineuses.
* Jeu de donnees synthetique minimal ajoute sous `test-data/minimal`, avec fixtures RNA-seq, proteines, Liftoff/Aegis et cas invalides.
* Les chemins par defaut pointent vers `data/` absent du depot, alors que `data_example/` contient des placeholders.
* Profils `local`, `test`, `apptainer` et `slurm` explicites; Apptainer force EGAPx imbrique en `singularity`.
* Aucune suite de tests TITAN detectee.

P1:

* Les conteneurs runtime et bases Dockerfile sont digest-pinned depuis P1-008.
* les `containerOptions` Docker-centriques ont ete retires des modules actifs depuis P1-009; les chemins genome, RNA-seq et scripts sont stages par Nextflow.
* les scans internes de dossiers publies ont ete retires des merges StringTie, GFFCompare, BRAKER3, AGAT, Salmon et des alignements STAR/HISAT2/Minimap2.
* Telechargements SRA pendant le workflow via metadonnees FASTQ ENA, avec retry et verification MD5 optionnelle.
* Versions non collectees dans les sorties.

P2/P3:

* Schema d'entrees Python present; schema Nextflow/nf-validation non encore introduit.
* Pas de contrat de sorties.
* Documentation des limites et troubleshooting incomplete.
* Pas de CI rapide detectee.

## Jeux de donnees de test

Decision P0-004: creer un jeu synthetique minimal sous `test-data/minimal`, documente comme tel. Il sert a valider les contrats d'entree TITAN, le parsing FASTA/GFF3/GTF/FASTQ, les samplesheets et les futurs tests stub. Il ne pretend pas valider biologiquement EDTA, BRAKER3, Liftoff, STAR, HISAT2, Minimap2 ou Aegis en production.

Cas couverts:

* assemblies `reference.fa` et `target.fa` avec deux seqids;
* annotation precedente `reference.gff3` valide avec gene, mRNA, exon, CDS;
* RNA-seq FASTQ gzip: un single-end, une paire paired-end et un long-read;
* samplesheet RNA-seq contenant `single`, `paired` et `long`;
* deux sources proteiques et un samplesheet proteines multi-organismes;
* parametre `input_egapx.yaml`;
* evidences pre-calculees minimales: EDTA masked assembly, Liftoff GFF3, BRAKER/AUGUSTUS GFF3, GeneMark GTF, STAR/StringTie, STAR/PsiCLASS et Minimap2/StringTie;
* GFF3 invalides: seqid absent, coordonnees invalides, Parent invalide;
* FASTA vide.

Validation:

```bash
python3 scripts/validate_minimal_test_data.py
sha256sum -c test-data/minimal/checksums.sha256
```

Checksums actuels: `test-data/minimal/checksums.sha256`.

## Profil test local P0-003

Le profil `test` est defini dans `conf/test.config` et inclus depuis `nextflow.config`. Il conserve les parametres historiques du profil standard et surcharge uniquement le contexte de validation rapide:

* executeur local;
* Docker desactive;
* `workDir = "${projectDir}/test-work"`;
* `output_dir = "${projectDir}/test-results"`;
* entrees pointees vers `test-data/minimal/valid`;
* `RNAseq_data_dir = "${projectDir}/test-data/minimal/valid/rnaseq"` pour decoupler les FASTQ de `data/RNAseq_data`;
* aucun mode `workflow` public: le profil test lance le graphe TITAN complet en stub;
* `edta_cpus = 2`, `egapx_cpus = 2` et `diamond2go_cpus = 2` pour permettre les tests stub sur une machine de developpement;
* ressources plafonnees a 2 CPU et 4 GB pour les labels utilises par les futurs tests.

Commandes validees:

```bash
nextflow config -profile test
nextflow run main.nf -profile test -stub-run -ansi-log false
```

Decision P1-004: TITAN n'expose plus de modes partiels. Le workflow lance la generation d'evidences puis Aegis dans le meme graphe, en passant `evidence_data.masked_genome` directement au sous-workflow `aegis`.

Limite: ce profil valide la resolution Nextflow, les shapes de channels et le lien EDTA -> Aegis en mode stub. Il ne valide pas biologiquement les outils lourds ni les conteneurs.

## Validation des parametres P0-005

Decision initiale P0: conserver la validation dans un bloc `workflow`, pour eviter les effets de bord de code executable top-level avec Nextflow 26. Depuis P1-002, cette validation vit dans `workflows/titan.nf`. Les controles couvrent:

* parametres obligatoires vides ou passes comme flags sans valeur;
* fichiers d'entree absents;
* rejet explicite de l'ancien parametre `--workflow`;
* `Channel.fromPath(..., checkIfExists: true)` pour les samplesheets.

Commandes validees:

```bash
nextflow config -profile test
nextflow run main.nf -profile test -stub-run -ansi-log false
nextflow run main.nf -profile test --RNAseq_samplesheet '' -stub-run -ansi-log false
nextflow run main.nf -profile test --previous_annotations test-data/minimal/valid/missing.gff3 -stub-run -ansi-log false
nextflow run main.nf -profile test --workflow aegis -stub-run -ansi-log false
```

Messages attendus verifies:

* `Missing required parameter(s): --RNAseq_samplesheet`
* `Required input file(s) not found`
* `--workflow is no longer supported`

## Validation avancee des entrees P2-001

Le validateur `scripts/validate_inputs.py` est appele au demarrage de `workflows/titan.nf` et par `launch_TITAN_example.sh`. Il bloque avant calcul lourd si les entrees ne respectent pas le contrat TITAN.

Commandes de validation:

```bash
python3 scripts/validate_inputs.py --project-dir . --new-assembly test-data/minimal/valid/target.fa --previous-assembly test-data/minimal/valid/reference.fa --previous-annotations test-data/minimal/valid/reference.gff3 --rnaseq-samplesheet test-data/minimal/valid/rnaseq_samplesheet.csv --rnaseq-data-dir test-data/minimal/valid/rnaseq --protein-samplesheet test-data/minimal/valid/protein_samplesheet.csv --egapx-paramfile test-data/minimal/valid/input_egapx.yaml
python3 scripts/test_validate_inputs.py
nextflow run main.nf -profile test --RNAseq_samplesheet test-data/minimal/invalid/rnaseq_bad_layout.csv -stub-run -ansi-log false
```

Contrats controles:

* FASTA non vides et IDs uniques;
* GFF3 a 9 colonnes, seqids presents dans FASTA, coordonnees valides, IDs uniques et Parents resolus;
* samplesheet RNA-seq avec colonnes `sampleID,SRA_or_FASTQ,library_layout`, valeurs enumerees, sample IDs uniques et fichiers locaux presents;
* samplesheet proteines avec colonnes `organism,filename` et FASTA presents;
* YAML EGAPx contenant `genome`, `taxid`, `organism` et genome lisible;
* options `egapx_executor`, `PSICLASS_vd_option`, `PSICLASS_c_option`.

## Strategie de modularisation

Structure cible progressive:

```text
workflows/titan.nf
subworkflows/local/evidence_generation/main.nf
subworkflows/local/aegis_integration/main.nf
modules/local/<outil>/main.nf
modules/local/<outil>/meta.yml
conf/base.config
conf/local.config
conf/test.config
conf/apptainer.config
conf/slurm.config
```

Migration recommandee: ne pas deplacer massivement les modules existants. Commencer par validation, profils, fixtures et un module simple; conserver les noms publics d'outputs.

## Strategie CI

CI rapide proposee:

* installer Nextflow;
* `scripts/run-tests.sh`, qui lance les validations statiques, les tests Python, `nextflow config -profile test`, le run stub minimal et un cas negatif d'entree invalide.

Ne pas tester Slurm reel dans GitHub Actions. Apptainer peut etre ajoute plus tard ou remplace par validation de configuration si l'environnement CI est insuffisant.

## Limitations environnementales

* Nextflow local est 26.04.3; le manifest TITAN declare 24.04.3.
* Les profils Docker/Apptainer/Slurm resolvent localement; les runs Slurm/Apptainer reels restent a valider sur le cluster cible.
* Aucune execution scientifique complete n'a ete lancee.
* Les fixtures `data_example` existantes sont majoritairement vides.
