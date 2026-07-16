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

Point d'entree principal: `main.nf`.

Sous-workflows:

* `subworkflows/generate_evidence_data.nf`: genere les evidences transcriptomiques, ab initio, Liftoff et EDTA.
* `subworkflows/aegis.nf`: integre les evidences via Aegis puis lance Diamond2GO.

Modules detectes: 27 process Nextflow dans `modules/*.nf`.

Scripts auxiliaires: `scripts/Aegis1.py`, `scripts/Aegis2.sh`, `scripts/Aegis3.py`, wrappers StringTie/EDTA, scripts de recuperation de chemins et scripts Python d'analyse/filtrage.

Donnees exemple: `data_example/` contient des placeholders tres legers, majoritairement de taille 0; ce repertoire ne constitue pas un jeu de test executable fiable.

## Commandes historiques

Documentees dans `README.md` et `launch_TITAN_example.sh`:

```bash
module load nextflow/24.04.3
nextflow run main.nf -with-dag dag_evidence_data.png --workflow generate_evidence_data
nextflow run main.nf -with-dag dag_aegis.png --workflow aegis
```

Le script historique contient un chemin absolu specifique a une machine: `/home/avelt/data2/.../workflows/TITAN`.

## Baseline executee

| Test | Commande | Resultat | Duree approximative | Remarque |
| --- | --- | --- | --- | --- |
| Version Nextflow | `nextflow -version` | Succes | <1 s | Nextflow 26.04.3 disponible; TITAN declare 24.04.3. |
| Configuration | `nextflow config` | Succes | ~2 s | Les parametres par defaut resolvent. |
| Lancement minimal historique initial | `nextflow run main.nf --workflow aegis -ansi-log false` | Echec | ~3 s | Erreur de compilation DSL2: statements top-level melanges avec declarations. |
| Profils locaux | `nextflow config -profile local >/dev/null && nextflow config -profile test >/dev/null && nextflow config -profile slurm,apptainer,test >/dev/null` | Succes | ~5 s | Verifie la resolution des profils ajoutes, sans soumission Slurm. |
| Lancement minimal apres P0 | `nextflow run main.nf --workflow aegis --output_dir test-data/minimal/valid --new_assembly test-data/minimal/valid/target.fa --previous_assembly test-data/minimal/valid/reference.fa --previous_annotations test-data/minimal/valid/reference.gff3 --RNAseq_samplesheet test-data/minimal/valid/rnaseq_samplesheet.csv --protein_samplesheet test-data/minimal/valid/protein_samplesheet.csv --egapx_paramfile test-data/minimal/valid/input_egapx.yaml --EDTA no --use_long_reads false -ansi-log false` | Succes | ~4 s | Compile et lance la branche Aegis sans process lourd; les evidences absentes sont signalees. |
| Reprise minimale | `nextflow run main.nf -profile test --workflow aegis -ansi-log false -resume` | Succes | ~4 s | Aucun process scientifique execute; valide la commande minimale et la creation de `output_dir`. |

Erreur exacte principale:

```text
Error main.nf:11:1: Statements cannot be mixed with script declarations -- move statements into a process, workflow, or function
```

Corrections P0 appliquees:

* validation des parametres et creation des channels de `generate_evidence_data` deplacees dans le bloc `workflow`;
* construction des channels RNA/proteines limitee a la branche `generate_evidence_data`;
* interpretation explicite de `params.use_long_reads` pour eviter qu'une chaine `"false"` soit truthy;
* creation de `output_dir` avant ecriture des placeholders `dev_null*` de la branche `aegis`;
* remplacement d'une boucle `for` non supportee par Nextflow 26 dans `subworkflows/aegis.nf`;
* correction du tag `diamond2go` qui referenceait une variable inexistante;
* emits Aegis vides lorsque `EDTA=no`, pour eviter de referencer des process non invoques.

## Fonctionnement biologique compris

TITAN annote un nouvel assemblage de genome en combinant:

* transfert d'annotations precedentes via Liftoff;
* annotation EGAPx prevue mais actuellement commentee dans le sous-workflow;
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
| Architecture Nextflow | `main.nf` leger + `workflows/`, `subworkflows/`, `modules/local/` | `main.nf` orchestre directement channels et branches | Orchestration trop chargee | Migrer progressivement vers `workflows/titan.nf`. |
| DSL2 | DSL2 structure valide | DSL2 active mais compilation cassee sous Nextflow 26 | Bloquant | Deplacer statements top-level dans `workflow`. |
| Point d'entree | Help et validation parametres | Pas de `--help`; validation top-level | UX et compilation | Ajouter help non executant et validation interne. |
| Parametres | YAML utilisateur + schema | Parametres dans `nextflow.config` | Peu portable | Introduire config YAML/schema sans casser les anciens `--param`. |
| Validation entrees | Scripts Python testes | Validation minimale | Risque calcul lourd invalide | Ajouter validateur FASTA/GFF3/samplesheets. |
| Modules | Process atomiques labels | Modules nombreux mais scripts/volumes couples | Fragile | Clarifier contrats et labels module par module. |
| Sous-workflows | Fonctions coherentes | Deux sous-workflows biologiques | Base utile | Garder mais typer les emits. |
| Dependances | Versions documentees et images | Mix BioContainers + images `latest` | Reproductibilite partielle | Verrouiller progressivement. |
| Conteneurs | Apptainer profile | Docker active globalement | Inadapte cluster | Ajouter profil Apptainer et garder Docker optionnel. |
| Profils | local/test/slurm/apptainer | Aucun profil explicite | Bloquant tests | Ajouter `conf/*.config`. |
| Ressources | Labels centralises | CPU dans modules + config globale | Peu maintenable | Centraliser par labels. |
| Donnees test | Synthetic fixtures versionnees | `data_example` placeholders | Non executable | Creer `test-data/minimal`. |
| Tests | pytest, nf-test, scripts | Aucun test detecte | Risque regression | Ajouter tests statiques puis nf-test. |
| CI | Scripts et workflows | Pas de CI TITAN detectee | Non verifie | CI rapide sans gros calcul. |
| Logs | Rapports structures | Echo dans scripts, logs implicites | Peu exploitable | Standardiser logs par process. |
| Rapports | MultiQC, trace, timeline, DAG | DAG possible manuellement | Incomplet | Documenter `-with-*`. |
| Documentation | README + docs | README utile mais incomplet | Manque exploitation | Ajouter docs dev/utilisateur. |
| Reproductibilite | versions.yml, provenance | Versions dans README seulement | Non collecte | Ajouter collecte versions. |
| Provenance | Manifestes et checksums | Peu de manifestes | Faible tracabilite | Ajouter manifestes inputs/outputs. |
| Gestion erreurs | Statuts structures | Warnings et erreurs parfois non bloquants | Risque silencieux | Echec explicite sur invalides. |
| Reprise | Tests `-resume` | Compatibilite non verifiee | Inconnu | Tester apres baseline. |
| Outputs | Contrat documente | Sorties publiees et copies manuelles | Risque collision | Contrat outputs TITAN. |
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
| Aegis | Integration finale annotation | v2025_05_20 | Partielle | `aegis_short_reads`, `aegis_long_reads` | `avelt/aegis:v2025_05_20` | Non | Non |
| AGAT | GFF3 vers CDS FASTA | 1.2.0 | Oui | `agat_convert_gff3_to_cds_fasta` | `quay.io/biocontainers/agat:1.2.0--pl5321hdfd78af_0` | Non | Non |
| BRAKER3 | Prediction ab initio | v3.0.8 documente | Non, image latest | `braker3_prediction*` | `avelt/braker3:latest` | Non | Non |
| Diamond2GO | Annotation fonctionnelle | commit documente | Non, image latest | `diamond2go` | `avelt/diamond2go:latest` | Non | Non |
| EDTA | Masquage TE | GitHub latest documente | Non | `EDTA` | `avelt/edta:latest` | Non | Non |
| EGAPx | Annotation NCBI | 0.3.2-alpha | Oui | `egapx` | `avelt/ncbi_egapx:0.3.2-alpha` | Non | Non |
| fastp | Trimming | 0.23.2 | Oui | `trimming_fastq` | `quay.io/biocontainers/fastp:0.23.2--hb7a2d85_2` | Non | Non |
| GFFCompare | Fusion PsiCLASS | 0.12.6 | Non, image latest | `gffcompare` | `avelt/gffcompare:latest` | Non | Non |
| HISAT2 | Index/alignment | 2.2.1 documente | Non, image latest | `hisat2_*` | `avelt/hisat2:latest` | Non | Non |
| Liftoff | Transfert annotation | 1.5.1 | Oui | `liftoff_annotations` | `quay.io/biocontainers/liftoff:1.5.1--py_0` | Non | Non |
| Minimap2 | Alignement long reads | 2.28 documente | Non, image latest | `minimap2_*` | `avelt/minimap2_samtools:latest` | Non | Non |
| PsiCLASS | Assemblage transcriptome | 1.0.2 | Non, image latest | `assembly_transcriptome_star_psiclass` | `avelt/psiclass_samtools:latest` | Non | Non |
| Salmon | Inference strandedness | 1.10.3 | Oui | `salmon_*` | `quay.io/biocontainers/salmon:1.10.3--haf24da9_3` | Non | Non |
| sra-tools | SRA download | 3.1.1 | Oui | `prepare_RNAseq_fastq_files_*` | `quay.io/biocontainers/sra-tools:3.1.1--h4304569_0` | Non | Non |
| STAR | Index/alignment | 2.7.11b | Oui | `star_*` | `quay.io/biocontainers/star:2.7.11b--h43eeafb_2` | Non | Non |
| StringTie | Assemblage/fusion | 2.2.3 documente | Non, image latest | `assembly_transcriptome_*`, `Stringtie_merging_*` | `avelt/stringtie:latest` | Non | Non |
| gffread | Conversion GTF/GFF3 | Incluse Aegis | A confirmer | `aegis_*` | `avelt/aegis:v2025_05_20` | Non | Non |
| samtools | BAM | 1.9 documente selon images | Partiel | alignements | images combinees | Non | Non |
| DIAMOND | Aegis/Diamond2GO/BRAKER3 | 2.1.9/2.1.11 selon docs | Partiel | `aegis_*`, `diamond2go`, `braker3` | images combinees | Non | Non |

### Inventaire des entrees

| Entree | Parametre/channel | Format attendu | Colonnes ou structure | Etapes consommatrices | Remarques |
| --- | --- | --- | --- | --- | --- |
| Assemblage cible | `params.new_assembly` | FASTA | fichier existant | Liftoff, AGAT, STAR, HISAT2, Minimap2, EDTA, BRAKER3, Aegis | Chemin separe en dossier + nom dans les modules. |
| Assemblage precedent | `params.previous_assembly` | FASTA | fichier existant | Liftoff | Doit correspondre a `previous_annotations`. |
| Annotation precedente | `params.previous_annotations` | GFF3 | fichier existant | Liftoff, puis AGAT | Source du transfert et de la FASTA CDS pour Salmon. |
| Samplesheet RNA-seq | `params.RNAseq_samplesheet` | CSV avec header | `sample_ID`, `SRA_or_FASTQ`, `library_layout`, `read_type` | preparation reads, fastp, Salmon, STAR, HISAT2, Minimap2 | `read_type == short` alimente la branche courte; `read_type == long` alimente la branche longue. |
| Samplesheet proteines | `params.protein_samplesheet` | CSV avec header | Structure consommee par les scripts BRAKER3/Aegis | BRAKER3, Aegis | Les modules montent aussi `${projectDir}/data/protein_data`; ce couplage reste a normaliser. |
| Parametres EGAPx | `params.egapx_paramfile` | YAML | parametres EGAPx | module `egapx` | Module present mais include et appel commentes dans `generate_evidence_data`. |
| Options biologiques | `params.EDTA`, `params.use_long_reads`, `params.PSICLASS_*`, `params.STAR_memory_per_job` | chaines, booleens, floats, entier bytes | valeurs scalaires | EDTA, branches long reads, PsiCLASS, STAR | `use_long_reads` est interprete differemment selon modules; a harmoniser. |
| Entrees Aegis-only | channel `input_data` construit dans `main.nf` | liste de paires cle/fichier | cles `masked_genome.masked_genome`, `braker3_results.*`, `previous_annotations.*`, `merged_*`, `gffcompare_out.*` | sous-workflow `aegis` | En mode minimal `EDTA=no`, des fichiers `dev_null*` remplacent les sorties unstranded absentes. |

### Inventaire des sorties par module

| Module/process | Entrees principales | Sorties declarees | Publication actuelle | Consommateur principal |
| --- | --- | --- | --- | --- |
| `prepare_RNAseq_fastq_files_short` | tuples RNA-seq courts | tuples `sample_ID`, `SRA_or_FASTQ`, `library_layout` | workdir uniquement | `trimming_fastq` |
| `prepare_RNAseq_fastq_files_long` | tuples RNA-seq longs | tuples `sample_ID`, `SRA_or_FASTQ`, `library_layout` | workdir uniquement | `minimap2_alignment` |
| `trimming_fastq` | FASTQ/SRA courts prepares | `*.trimmed.fastq.gz` | `${output_dir}/intermediate_files/evidence_data/RNAseq_data/trimmed_data` | Salmon, STAR, HISAT2 |
| `liftoff_annotations` | assemblage cible, assemblage precedent, GFF3 precedent | `liftoff_previous_annotations.gff3`, `unmapped_features.txt` | `${output_dir}` | AGAT, Aegis |
| `agat_convert_gff3_to_cds_fasta` | assemblage cible, Liftoff GFF3 | `${genome}.CDS.fasta.gz` | `${output_dir}/intermediate_files/liftoff/gff3_to_cds_fasta` | Salmon index |
| `salmon_index` | CDS FASTA | `salmon_index/` | `${output_dir}/intermediate_files/salmon_index` | Salmon strand inference |
| `salmon_strand_inference` | reads trimmes, index Salmon | `${sample_ID}.strand_info.classified` | `${output_dir}/intermediate_files/salmon_strand` | STAR/HISAT2 channel enrichi avec strandedness |
| `star_genome_indices` | assemblage cible | `${genome}_index` | `${output_dir}/intermediate_files/evidence_data/star_databases` | STAR alignment |
| `star_alignment` | index STAR, reads, strandedness | `${sample_ID}_Aligned.sortedByCoord.out.bam` | `${output_dir}/intermediate_files/evidence_data/RNAseq_alignments/STAR/{stranded,unstranded}` | StringTie, PsiCLASS, BRAKER3 |
| `hisat2_genome_indices` | assemblage cible | `${genome}.*.ht2` | `${output_dir}/intermediate_files/evidence_data/hisat2_databases` | HISAT2 alignment |
| `hisat2_alignment` | index HISAT2, reads, strandedness | `${sample_ID}_Aligned.sort.bam` | `${output_dir}/intermediate_files/evidence_data/RNAseq_alignments/HISAT2/{stranded,unstranded}` | StringTie HISAT2 |
| `minimap2_genome_indices` | assemblage cible | `${genome}.mmi` | `${output_dir}/intermediate_files/evidence_data/minimap2_databases` | Minimap2 alignment |
| `minimap2_alignment` | index Minimap2, reads longs | `${sample_ID}_Aligned.sorted.bam` | `${output_dir}/intermediate_files/evidence_data/RNAseq_alignments/minimap2` | StringTie long reads, BRAKER3 long reads |
| `assembly_transcriptome_star_stringtie` | BAM STAR | `*_transcriptome.gtf`, `*_transcriptome.AltCommands.gtf` | `${output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/STAR/{stranded,unstranded}` | StringTie STAR merge |
| `assembly_transcriptome_hisat2_stringtie` | BAM HISAT2 | `*_transcriptome.gtf`, `*_transcriptome.AltCommands.gtf` | `${output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/HISAT2/{stranded,unstranded}` | StringTie HISAT2 merge |
| `assembly_transcriptome_star_psiclass` | BAM STAR | `${sample_ID}_vote.gtf` | `${output_dir}/intermediate_files/transcriptomes/STAR_PsiCLASS/{stranded,unstranded}` | GFFCompare |
| `assembly_transcriptome_minimap2_stringtie` | BAM Minimap2 | `*_transcriptome.gtf`, `*_transcriptome.AltCommands.gtf` | `${output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/long_reads` | StringTie long reads merge |
| `Stringtie_merging_short_reads_STAR` | liste GTF STAR/StringTie | `merged_transcriptomes.STAR.short_reads.*.gtf` | `${output_dir}/tmp` plus copies script vers `${output_dir}` | Aegis |
| `Stringtie_merging_short_reads_hisat2` | liste GTF HISAT2/StringTie | `*.gtf` | `${output_dir}` | Non utilise actuellement par Aegis |
| `Stringtie_merging_long_reads` | liste GTF Minimap2/StringTie | `merged_transcriptomes.minimap2.long_reads.*.gtf` | `${output_dir}/tmp` plus copies script vers `${output_dir}` | Aegis long reads |
| `gffcompare` | liste GTF PsiCLASS | `stranded_merged_output.combined.gtf`, optionnel `unstranded_merged_output.combined.gtf` | `${output_dir}/tmp` | Aegis |
| `EDTA` | assemblage cible | `*TElib.fa`, `*TEanno.gff3`, `*MAKER.masked` | `${output_dir}/tmp` plus copies script vers `${output_dir}` | Aegis |
| `braker3_prediction` | genome, proteines, BAM courts | `augustus.hints.gff3`, `genemark.gtf`, `genemark_supported.gtf`, `braker.gff3` | `${output_dir}` | Aegis |
| `braker3_prediction_with_long_reads` | genome, proteines, BAM courts + longs | `augustus.hints.gff3`, `genemark.gtf`, `genemark_supported.gtf`, `braker.gff3` | `${output_dir}` | Aegis long reads |
| `aegis_short_reads` | genome masque, BRAKER3, Liftoff, STAR/StringTie, GFFCompare | `final_annotation.gff3`, `final_annotation_proteins_all.fasta`, `final_annotation_proteins_main.fasta` | `${output_dir}/aegis_outputs` | Diamond2GO |
| `aegis_long_reads` | idem short + StringTie long reads | `final_annotation.gff3`, `final_annotation_proteins_all.fasta`, `final_annotation_proteins_main.fasta` | `${output_dir}/aegis_outputs` | Diamond2GO |
| `diamond2go` | proteines Aegis all/main | `*-diamond*` | `${output_dir}/Diamond2GO_outputs` | sortie finale fonctionnelle |
| `egapx` | YAML EGAPx | `*` | `${output_dir}/egapx` | Non branche actuellement |

Points de vigilance P0-002:

* plusieurs modules publient via `publishDir` et copient aussi manuellement dans `${output_dir}`;
* les images `latest` rendent l'inventaire de versions non reproductible sans validation externe;
* aucune version n'est collectee dans un artefact de sortie;
* les schemas reels de samplesheets ne sont pas valides avant construction des channels;
* EGAPx est documente et le module existe, mais la branche est commentee dans le sous-workflow.

## Constats critiques

P0:

* `main.nf` ne compilait pas avec Nextflow 26.04.3 avant correction P0; la commande minimale `aegis` compile et termine maintenant avec `EDTA=no`.
* Pas de profil `test` local, donc pas de validation rapide sans Slurm/Docker lourd.
* Les chemins par defaut pointent vers `data/` absent du depot, alors que `data_example/` contient des placeholders.
* Docker est active globalement; Apptainer n'est pas configure.
* Aucune suite de tests TITAN detectee.

P1:

* Nombreuses images `latest`.
* `containerOptions` montent `projectDir/work`, `projectDir/data` ou `/outputdir`, ce qui fragilise Apptainer et `-resume`.
* Certains scripts utilisent des chemins historiques, glob fragile, `ls`, backticks, `continue` hors boucle utile dans des scripts de process.
* Telechargements SRA pendant le workflow avec `--max-size 100G`.
* Versions non collectees dans les sorties.

P2/P3:

* Pas de schema de parametres.
* Pas de contrat de sorties.
* Documentation des limites et troubleshooting incomplete.
* Pas de CI rapide detectee.

## Jeux de donnees de test

Decision initiale: creer un jeu synthetique minimal sous `test-data/minimal`, documente comme tel. Il sert d'abord a valider parsing FASTA/GFF3 et futurs tests statiques. Il ne pretend pas valider biologiquement EDTA, BRAKER3 ou Aegis en production.

Cas prevus:

* FASTA valide avec deux seqids;
* GFF3 valide avec gene, mRNA, exon, CDS;
* proteines minimales;
* samplesheets minimalistes;
* GFF3 invalides: seqid absent, coordonnees invalides, Parent invalide;
* FASTA vide.

Checksums actuels: `test-data/minimal/checksums.sha256`.

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
* `nextflow config`;
* `nextflow config -profile test`;
* tests statiques Bash/Python;
* `nextflow run main.nf -profile test -stub-run` quand le profil test est executable.

Ne pas tester Slurm reel dans GitHub Actions. Apptainer peut etre ajoute plus tard ou remplace par validation de configuration si l'environnement CI est insuffisant.

## Limitations environnementales

* Nextflow local est 26.04.3; le manifest TITAN declare 24.04.3.
* Docker/Apptainer n'ont pas encore ete valides pour TITAN pendant cette passe.
* Aucune execution scientifique complete n'a ete lancee.
* Les fixtures `data_example` existantes sont majoritairement vides.
