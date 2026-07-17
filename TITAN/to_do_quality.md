# Audit qualite Nextflow de TITAN

Date de l'audit: 2026-07-16

Perimetre lu: `README.md`, `nextflow.config`, `main.nf`, `workflows/titan.nf`, `subworkflows/*.nf`, `modules/*.nf`, configuration de base et strategie de tests.

## Priorites transversales

### P0 - Robustesse du graphe et reproductibilite

- [x] Eviter de modifier `params` dans `main.nf`. Les valeurs par defaut doivent vivre dans `nextflow.config` ou dans un schema de parametres. Muter `params` au runtime rend les profils et les overrides moins lisibles.
- [x] Reduire l'usage direct de `params` dans les modules. Les modules DSL2 devraient recevoir leurs options metier en `val` ou `path`, et garder seulement `container`, `label` et eventuellement `publishDir` comme configuration externe.
- [x] Ajouter `set -euo pipefail` a tous les scripts shell non triviaux, surtout ceux avec pipes (`samtools`, `grep`, `awk`, `STAR`, `HISAT2`, `fastp`, `liftoff`, `EDTA`).
- [x] Remplacer les sorties glob trop larges (`file("*.gtf")`, `file("*.trimmed.fastq.gz")`, `path("${sample_ID}*")`) par des noms explicites ou des tuples normalises. Les glob larges peuvent capturer des fichiers temporaires et casser les contrats de canal.
- [x] Normaliser les canaux vides. Plusieurs merges utilisent `collect().ifEmpty([])` puis un `path(...)`; c'est fragile parce que `[]` n'est pas un fichier. Preferer produire un fichier sentinelle explicite, ou scinder les branches optionnelles avec des workflows conditionnels.
- [x] Ajouter des `versions.yml` a tous les processus publics, pas seulement EGAPx, AEGIS et provenance. C'est une pratique courante nf-core et utile pour tracer les outils.
- [x] Ajouter des tests de garde avant refactor: prepare FASTQ, trimming, alignements, merges, BRAKER, AEGIS, provenance. `nf-test` n'est pas disponible dans l'environnement courant; le verrouillage ajoute `scripts/validate_nextflow_quality.py` a la suite rapide et conserve le `stub-run` complet.

### P1 - Maintenabilite DSL2

- [x] Renommer les processus en `UPPER_SNAKE_CASE` ou adopter une convention stable. Le projet melange `EDTA`, `Stringtie_merging_short_reads_STAR`, `aegis_short_reads`, etc. Convention retenue: conserver les noms existants pour ne pas casser `withName`, traces et resume; les nouveaux modules doivent utiliser lower snake case aligne sur le fichier module.
- [x] Factoriser les duplications: STAR/HISAT2/StringTie ont des scripts tres proches; AEGIS short/long aussi; BRAKER short/long aussi. AEGIS short/long utilise maintenant `scripts/run_aegis_merge.sh`; les trois modules StringTie utilisent `scripts/run_stringtie_transcriptome.sh`.
- [x] Sortir les scripts shell complexes vers `scripts/` avec tests unitaires. Les modules Nextflow devraient orchestrer, pas contenir trop de logique bash. Les helpers partages sont testes par `scripts/test_shared_shell_scripts.py`.
- [x] Remplacer les lectures de fichier dans des closures de workflow (`file(strand_file).text`) par une etape process ou par une structure de sortie deja calculee. `salmon_strand_inference` emet maintenant le strand type via `env('STRAND_INFO')`.
- [x] Revoir les `publishDir` pour separer outputs publics et intermediaires. Publier beaucoup d'intermediaires avec `mode: copy` augmente le cout disque et rend les resultats moins lisibles. Les publications intermediaires sont maintenant controlees par `params.publish_intermediates`; les sorties publiques restent publiees.

### P2 - Hygiene et ergonomie

- [x] Harmoniser indentation, espaces autour de `emit:`, `input:` et `output:`, fins de ligne et commentaires. Les fichiers Nextflow ont ete normalises et le validateur qualite bloque CRLF, trailing whitespace et `emit :`.
- [x] Preferer `path(...)` a `file(...)` dans les sorties DSL2 pour rester coherent avec Nextflow moderne. Le validateur qualite bloque maintenant `file(...)` dans les blocs `output:`.
- [x] Ajouter des `tag` plus informatifs aux processus qui fusionnent beaucoup de fichiers, par exemple nombre de GTFs. Les merges StringTie et gffcompare decrivent maintenant le type d'evidence fusionnee.
- [x] Documenter les contrats de tuple de chaque module dans des commentaires courts ou dans des tests. Les contrats principaux RNA-seq/proteines/evidence sont documentes dans les subworkflows et verifies par la suite de tests.

## Audit fichier par fichier

### `README.md`

- Le contrat actuel est clair: graphe unique complet, modes partiels supprimes, validation en amont.
- [x] A ameliorer: ajouter une section "Developer quality contract" listant les conventions DSL2 attendues: modules sans logique globale, sorties nommees, `versions.yml`, tests `nf-test`, labels obligatoires.
- [x] A ameliorer: indiquer quels outputs sont publics et quels outputs sont seulement des intermediaires. Cela aidera a rationaliser les `publishDir`.

### `nextflow.config`

- Les containers sont pin avec digest, bon point pour la reproductibilite.
- [x] A ameliorer: les chemins par defaut pointent vers `data/...` alors que les exemples actuels sont sous `data_example/` et `test-data/`. Eviter des defaults de production non valides; utiliser `false` ou des fixtures de test seulement dans le profil `test`.
- [x] A ameliorer: `docker.enabled = true` au niveau global peut surprendre avec les profils HPC. Preferer activer Docker dans `conf/local.config`, et laisser le runtime au profil.
- [x] A ameliorer: harmoniser `container_egapx` et `egapx_container`, `container_aegis` et `aegis_container`. Garder un seul nom canonique avec compatibilite temporaire si necessaire. `container_egapx` et `container_aegis` sont canoniques; les anciens noms restent des alias.
- [x] A ameliorer: ajouter un schema `nextflow_schema.json` pour typer et documenter les parametres.

### `conf/base.config`

- Bon point: ressources centralisees par labels.
- A ameliorer: eviter `errorStrategy = 'finish'` global si certains processus ont besoin de retry controle. Definir les retries par label ou par process.
- A ameliorer: ajouter `maxRetries`, `withName` pour les processus reseau (`prepare_RNAseq_*`, `egapx`) et eventuellement `scratch`/`stageInMode` selon HPC.
- A ameliorer: revoir les labels inutilises (`process_medium`, `process_high`) ou les documenter.

### `main.nf`

- A ameliorer: retirer les assignations `params.* = false` du script principal. Elles ecrasent mentalement le contrat des profils et ne remplacent pas une validation de schema.
- A ameliorer: deplacer les defaults EGAPx/AEGIS vers `nextflow.config`.
- A ameliorer: garder `main.nf` minimal: activation DSL2, include du workflow, appel du workflow.

### `workflows/titan.nf`

- Bon point: validation d'entrees avant lancement des etapes lourdes.
- A ameliorer: remplacer `command.execute()` pour `scripts/validate_inputs.py` par un process Nextflow dedie ou une validation launcher. Executer un process local dans le workflow echappe au modele Nextflow, aux containers, a la trace et a la portabilite.
- A ameliorer: parser le CSV avec une seule source de canaux, puis brancher en `long`, `single`, `paired`; actuellement le samplesheet est relu trois fois.
- A ameliorer: `rnaseqLocalFiles()` retourne parfois liste, parfois fichier, parfois liste vide. Normaliser le tuple, par exemple `[meta, reads]` avec `reads` toujours liste.
- A ameliorer: `samplesheetHasLongReads()` parse le CSV manuellement avec `split(',')`; cela casse les champs quotes. Reutiliser le validateur Python ou `splitCsv`.
- A ameliorer: passer les fichiers principaux (`new_assembly`, `previous_assembly`, etc.) explicitement au subworkflow au lieu de relire `params` dans `generate_evidence_data`.
- A ameliorer: ajouter un process de validation final des parametres EGAPx/AEGIS si ces parametres restent runtime.

### `subworkflows/generate_evidence_data.nf`

- A ameliorer: le subworkflow depend fortement de `params` et `projectDir`; passer les fichiers, scripts et options en `take:`.
- A ameliorer: `file(strand_file).text.trim()` dans une closure lit un fichier de sortie cote orchestration. Preferer que `salmon_strand_inference` emette directement `val(strand_type)` et conserve le fichier en output secondaire si besoin.
- A ameliorer: les nombreux `collect()` bloquent jusqu'a la fin de la branche et peuvent masquer les canaux vides. Ajouter des tests pour les cas "aucun unstranded", "aucun long read", "un seul sample".
- A ameliorer: `merged_hisat2_stringtie` est execute mais ses outputs ne sont pas emis ni utilises par AEGIS. Decider si HISAT2 est un output public, une evidence AEGIS manquante ou un calcul inutile.
- A ameliorer: `Channel.value([])` pour les sorties long reads absentes est fragile avec les modules qui attendent `path`. Utiliser fichiers sentinelles crees par process ou structurer `aegis` pour ne pas prendre ces inputs en mode short.
- A ameliorer: `protein_list.map { file(filename) }.collect()` devrait valider et normaliser les chemins avant BRAKER, avec un tuple stable et une erreur claire si liste vide.
- A ameliorer: factoriser les blocs de filtrage stranded/unstranded repetes pour STAR, HISAT2, PsiCLASS.

### `subworkflows/aegis.nf`

- A ameliorer: la branche conditionnelle appelle soit `aegis_long_reads`, soit `aegis_short_reads`, mais l'`emit` reference les deux symboles dans une expression ternaire. Cela peut devenir fragile selon evaluation DSL2. Preferer definir des canaux intermediaires communs dans chaque branche.
- A ameliorer: `diamond2go.out` est emis sans nommage precis; le module `diamond2go` devrait avoir des outputs nommes.
- A ameliorer: factoriser les modules AEGIS short/long pour eviter deux scripts presque identiques.

### `modules/EDTA.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: eviter deux `publishDir` sur le meme process si possible; separer outputs publics et temporaires explicitement.
- A ameliorer: les globs `*TElib.fa`, `*TEanno.gff3`, `*MAKER.masked` peuvent capturer plusieurs fichiers. Renommer/copier vers des noms fixes avant `output:`.

### `modules/liftoff_annotations.nf`

- Bon point: `cache 'deep'` est pertinent pour un mapping deterministe sur fichiers d'entree.
- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: verifier que `unmapped_features.txt` existe meme si Liftoff ne le produit pas dans certains cas.
- A ameliorer: simplifier les inputs `val(genome)` et `val(previous_assembly)` si seuls les noms servent au `tag`.

### `modules/egapx.nf`

- Bon point: script defensif avec `set -euo pipefail` et `versions.yml`.
- A ameliorer: le module telecharge du code depuis GitHub a l'execution. Pour une reproductibilite stricte, preferer un runner prepackaged/pinned ou un cache gere hors process.
- A ameliorer: le module orchestre un autre workflow Nextflow dans un process. Documenter clairement les implications de nesting Nextflow, executors, workdir, cache et logs.
- A ameliorer: ajouter un `container` ou documenter pourquoi le process hote doit avoir `curl`, `tar`, `python3` et Docker/Apptainer disponibles.
- A ameliorer: rendre les noms attendus dans `egapx_out/` configurables ou valides avec messages d'erreur explicites avant les `cp`.

### `modules/agat_convert_gff3_to_cds_fasta.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: remplacer la logique `grep | awk | sed` par un script teste, car elle encode une correction PN40024 specifique dans un module generique.
- A ameliorer: emettre un output nomme (`emit: cds_fasta`).
- A ameliorer: gerer le cas `to_remove.txt` vide sans risque de comportement ambigu de `grep -f`.
- A ameliorer: ajouter `versions.yml`.

### `modules/salmon_index.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: publier l'index seulement si c'est un output utile; sinon garder en workdir/cache.

### `modules/salmon_strand_inference.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: emettre directement `val(strand_info)` pour eviter une lecture de fichier dans le subworkflow.
- A ameliorer: gerer explicitement les layouts inattendus avec `exit 1`; actuellement le script peut continuer sans log Salmon.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: publier logs Salmon ou les inclure comme output optionnel pour debug.

### `modules/prepare_RNAseq_fastq_files_short.nf`

- Bon point: retries applicatifs pour les downloads ENA et `errorStrategy 'retry'`.
- A ameliorer: retirer `debug true` par defaut ou le rendre activable par parametre.
- A ameliorer: la sortie `path("${sample_ID}*.fastq.gz", includeInputs: true)` est trop large. Emettre explicitement paired ou single, idealement avec un tuple `meta`.
- A ameliorer: passer `download_sra_fastq.py` comme `path` input pour rendre la dependance visible dans la trace/cache.
- A ameliorer: ajouter `maxRetries` controle dans la config pour ce process.
- A ameliorer: ajouter `versions.yml`.

### `modules/prepare_RNAseq_fastq_files_long.nf`

- A ameliorer: retirer `debug true` par defaut.
- A ameliorer: `path("${sample_ID}*", includeInputs: true)` peut capturer FASTA, FASTQ, fichiers temporaires ou logs. Emettre une sortie normalisee.
- A ameliorer: gerer separement FASTA et FASTQ dans le tuple de sortie pour que minimap2 sache le type sans inferer par nom.
- A ameliorer: passer `download_sra_fastq.py` comme input.
- A ameliorer: ajouter `versions.yml`.

### `modules/trimming_fastq.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: pour layout inconnu, echouer au lieu de seulement afficher un warning, sinon l'output attendu manquera avec une erreur moins claire.
- A ameliorer: remplacer `file("*.trimmed.fastq.gz")` par sorties explicites single/paired.
- A ameliorer: ajouter les rapports fastp JSON/HTML comme outputs optionnels publies.
- A ameliorer: ajouter `versions.yml`.

### `modules/star_genome_indices.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: eviter `chmod -R 755` si non necessaire; cela modifie tout l'index et peut couter cher.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: verifier les parametres STAR importants (`genomeSAindexNbases`, annotation GTF optionnelle) pour petits genomes/tests.

### `modules/star_alignment.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: gerer explicitement les layouts inconnus avec `exit 1`.
- A ameliorer: indexer le BAM ou documenter pourquoi l'index n'est pas requis par les etapes aval.
- A ameliorer: `params.STAR_memory_per_job` devrait etre derive de `task.memory` ou passe comme option explicite.
- A ameliorer: ajouter `versions.yml`.

### `modules/hisat2_genome_indices.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: emettre l'index sous forme de repertoire fixe plutot qu'un glob `path("${genome}.*.ht2")`; c'est plus facile a stage et a passer a l'alignement.
- A ameliorer: eviter `chmod -R 755 ${genome}*` trop large.
- A ameliorer: ajouter `versions.yml`.

### `modules/hisat2_alignment.nf`

- A ameliorer: ajouter `set -euo pipefail` avec `pipefail`, indispensable car `hisat2 | samtools sort` peut masquer un echec d'alignement.
- A ameliorer: gerer explicitement les layouts et strand types inconnus avec `exit 1`.
- A ameliorer: indexer le BAM si reutilise par d'autres outils.
- A ameliorer: pour single-end strandness, verifier que `FR/RF` est biologiquement correct pour HISAT2; les options single-end attendues sont souvent `F` ou `R`.
- A ameliorer: ajouter `versions.yml`.

### `modules/minimap2_genome_indices.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: publier l'index seulement si necessaire aux utilisateurs.

### `modules/minimap2_alignment.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: les trois branches FASTQ/SRA/FASTA executent la meme commande; factoriser.
- A ameliorer: remplacer `samtools view -b ... | samtools sort - > file.bam` par `samtools sort -o file.bam -` avec `pipefail`.
- A ameliorer: supprimer le SAM temporaire via pipe direct minimap2 vers samtools pour reduire I/O.
- A ameliorer: ajouter `versions.yml`.

### `modules/assembly_transcriptome_star_stringtie.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: les scripts `Stringtie.sh` et `Stringtie_AltCommands.sh` devraient etre versionnes comme inputs ou integres dans un module plus explicite.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: factoriser avec les modules HISAT2/minimap2 StringTie.

### `modules/assembly_transcriptome_hisat2_stringtie.nf`

- A ameliorer: memes points que le module STAR/StringTie: `set -euo pipefail`, `versions.yml`, factorisation.
- A ameliorer: harmoniser le tag et les chemins publies avec STAR pour faciliter la comparaison.

### `modules/assembly_transcriptome_minimap2_stringtie.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: verifier si les options StringTie long reads sont adaptees a FASTA et FASTQ SRA de maniere uniforme.

### `modules/assembly_transcriptome_star_psiclass.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: emettre avec `path(...)` plutot que `file(...)`.
- A ameliorer: passer `PSICLASS_vd_option` et `PSICLASS_c_option` comme `val` inputs ou les documenter dans un schema.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: chemin publie manque `evidence_data` contrairement aux autres transcriptomes; harmoniser.

### `modules/Stringtie_merging_short_reads_STAR.nf`

- Bon point: sorties nommees et `set -euo pipefail`.
- A ameliorer: `printf '%s\n' ${stranded_default_gtfs}` casse si un chemin contient espace et produit une ligne vide si liste vide. Construire les listes via une boucle bash sur `"$@"` est difficile dans Nextflow; mieux vaut generer les listes dans un petit script ou garantir les canaux non vides en amont.
- A ameliorer: les outputs unstranded sont declares `optional: true`, mais le script cree toujours des fichiers vides. Choisir un seul contrat: optionnel absent ou fichier sentinelle vide.
- A ameliorer: ajouter `versions.yml`.

### `modules/Stringtie_merging_short_reads_hisat2.nf`

- A ameliorer: output `file("*.gtf")` trop vague et sans `emit`. Nommer les quatre sorties comme pour STAR.
- A ameliorer: ce module est actuellement execute mais ses resultats ne sont pas emis par le subworkflow. Clarifier son utilite.
- A ameliorer: harmoniser les `publishDir` avec STAR et ajouter `saveAs` si outputs publics.
- A ameliorer: ajouter `versions.yml`.

### `modules/Stringtie_merging_long_reads.nf`

- Bon point: sorties nommees et `set -euo pipefail`.
- A ameliorer: robustifier la generation des listes de GTFs comme pour STAR.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: verifier que l'extension `emit: default_args_gff` correspond au contenu GTF; renommer en `*_gtf`.

### `modules/gffcompare.nf`

- Bon point: sorties nommees et publication publique avec `saveAs`.
- A ameliorer: le script cree toujours `unstranded_merged_output.combined.gtf`, donc `optional: true` n'est pas coherent.
- A ameliorer: proteger le cas `stranded_gtfs` vide avec une erreur explicite.
- A ameliorer: ajouter `versions.yml`.

### `modules/braker3_prediction.nf`

- Bon point: `set -euo pipefail`.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: les chemins internes `/BRAKER-3.0.8`, `/ProtHint-2.6.0`, etc. devraient etre documentes dans le Dockerfile ou exposes comme constantes testees.
- A ameliorer: la boucle `for file in ${protein_fastas}` est fragile avec espaces. Preferer fichiers sans espaces valides en amont ou script Python de preparation.
- A ameliorer: verifier que `bam_short` n'est jamais vide et produire une erreur claire.
- A ameliorer: publier seulement les outputs finaux publics, et garder le reste de BRAKER en workdir ou dans un dossier intermediaire dedie.

### `modules/braker3_prediction_with_long_reads.nf`

- A ameliorer: factoriser avec `braker3_prediction.nf`; seule la construction de la liste BAM change.
- A ameliorer: memes points: `versions.yml`, robustesse des listes, validation BAM/proteines, documentation des chemins internes.
- A ameliorer: nommer clairement `bam_short`, `bam_long`, `bam_all` dans logs/provenance.

### `modules/aegis_short_reads.nf`

- Bon point: `set -euo pipefail`, verification des FASTA proteines produites, `versions.yml`.
- A ameliorer: le script de merge est presque duplique avec `aegis_long_reads`; factoriser.
- A ameliorer: les checks `basename != dev_null*` suggerent un ancien contrat de sentinelle. Remplacer par un contrat explicite de fichiers optionnels/vides.
- A ameliorer: ajouter une validation des inputs obligatoires non vides avant `aegis merge`.
- A ameliorer: inclure le repertoire AEGIS complet comme output optionnel si utile au debug, ou documenter qu'il reste dans `work/`.

### `modules/aegis_long_reads.nf`

- A ameliorer: memes points que `aegis_short_reads.nf`.
- A ameliorer: valider que les deux GTF long reads sont non vides avant merge.
- A ameliorer: factoriser la liste `merge_inputs` pour eviter divergence avec le module short.

### `modules/diamond2go.nf`

- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: output `path("*-diamond*")` trop vague et sans `emit`. Nommer les resultats all/main.
- A ameliorer: ajouter `versions.yml`.
- A ameliorer: utiliser `${task.cpus}` si Diamond2GO/DIAMOND le supporte; sinon documenter pourquoi `diamond2go_cpus` ne change que la ressource reservee.
- A ameliorer: verifier que les deux commandes ne s'ecrasent pas si Diamond2GO produit des noms bases sur le repertoire plutot que sur la query.

### `modules/validate_final_annotation.nf`

- A ameliorer: ajouter `container params.container_python` ou un container dedie pour ne pas dependre du Python hote.
- A ameliorer: passer `scripts/validate_final_annotation.py` comme input `path` pour rendre la dependance visible dans le cache.
- A ameliorer: ajouter `set -euo pipefail`.
- A ameliorer: ajouter `versions.yml` ou integrer le hash/version du script de validation.

### `modules/titan_provenance.nf`

- Bon point: manifest JSON avec tailles et SHA-256.
- A ameliorer: ajouter `container params.container_python`.
- A ameliorer: passer la configuration a la provenance via inputs explicites plutot que interpoler beaucoup de `params`.
- A ameliorer: `records()` split sur espaces, ce qui casse les chemins contenant espaces et les listes Nextflow rendues en chaine. Preferer fournir un fichier liste ou serializer les chemins en JSON avant le process.
- A ameliorer: ajouter au manifest les versions de tous les modules une fois `versions.yml` generalise.
- A ameliorer: inclure `workflow.revision`, `workflow.commitId`, `workflow.nextflow.version`, profils actifs et ligne de commande si disponible.

## Plan de refactor recommande

1. Stabiliser les contrats sans changer la biologie: outputs nommes, `set -euo pipefail`, `versions.yml`, containers Python pour validation/provenance.
2. Ajouter des tests `nf-test` sur modules critiques avec fixtures minimales et stub mode.
3. Normaliser les tuples RNA-seq et proteins dans `workflows/titan.nf`, puis passer les inputs explicitement aux subworkflows.
4. Refactoriser les duplications AEGIS, BRAKER et StringTie.
5. Revoir publication et provenance: distinguer outputs publics, intermediaires utiles et artefacts de debug.
