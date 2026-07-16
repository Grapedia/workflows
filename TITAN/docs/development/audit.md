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
