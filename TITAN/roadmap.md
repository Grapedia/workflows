# Roadmap TITAN

Cette roadmap est specifique a TITAN. Elle transpose les bonnes pratiques observees dans `/home/vmadmin/atcg-rnaseq` vers un workflow d'annotation de genomes, sans transformer TITAN en pipeline RNA-seq.

## TITAN-P0-001 - Etablir une baseline fonctionnelle
Priorite : P0
Statut : Fait
Risque : Faible

### Objectif
Documenter l'etat initial, les commandes historiques et les erreurs bloquantes avant refactor.
### Constat
`nextflow config` fonctionne. Le lancement minimal echouait initialement a la compilation DSL2; une correction P0 a permis une commande minimale `aegis`. Depuis le durcissement P1, Aegis echoue volontairement si les evidences obligatoires manquent.
### Fichiers concernes
`main.nf`, `nextflow.config`, `docs/development/audit.md`.
### Etapes d'implementation
Capturer les commandes, corriger uniquement le blocage de compilation, verifier que la configuration reste resolvable.
### Tests
`nextflow config`; `nextflow run main.nf --workflow aegis -ansi-log false` ou commande test adaptee.
### Criteres d'acceptation
La baseline est documentee et une commande minimale atteint la phase workflow sans erreur de syntaxe.
### Risques et retour arriere
Risque faible si la correction se limite au placement de logique dans `workflow`. Retour arriere par revert du commit.

## TITAN-P0-002 - Inventorier entrees, sorties et outils
Priorite : P0
Statut : Fait
Risque : Faible

### Objectif
Obtenir un inventaire exploitable des outils, versions, conteneurs, entrees et sorties.
### Constat
Les versions sont partiellement documentees dans `README.md`; plusieurs images utilisent `latest`. L'inventaire P0-002 est maintenant consolide dans `docs/development/audit.md` et resume dans `README.md`.
### Fichiers concernes
`README.md`, `modules/*.nf`, `docs/development/audit.md`.
### Etapes d'implementation
Extraire les process, commandes, images, publishDir et dependances implicites.
### Tests
Inspection statique par `rg` et `awk`.
### Criteres d'acceptation
Tableau outil/etape/version/module/conteneur/test present dans l'audit.
### Risques et retour arriere
Aucun impact runtime.

## TITAN-P0-003 - Ajouter un profil test local minimal
Priorite : P0
Statut : Fait
Risque : Moyen

### Objectif
Permettre une execution de validation rapide sans Slurm ni donnees volumineuses.
### Constat
Le profil `test` local est disponible dans `conf/test.config`; il pointe vers `test-data/minimal/valid`, desactive Docker, isole `test-work`/`test-results` et permet une commande minimale sans Slurm.
### Fichiers concernes
`nextflow.config`, `conf/test.config`, `conf/local.config`, `test-data/`.
### Etapes d'implementation
Ajouter profils local/test sans changer les parametres historiques par defaut; pointer le profil test vers fixtures minimales.
### Tests
`nextflow config -profile test`; `nextflow run main.nf -profile test --workflow aegis -ansi-log false`; `nextflow run main.nf -profile test --workflow aegis -ansi-log false -resume`.
### Criteres d'acceptation
Le profil test resout et ne depend pas de Slurm.
### Risques et retour arriere
Risque de conflit de precedences de profils; conserver les defaults historiques hors profil.

## TITAN-P0-004 - Creer un jeu de donnees test minimal
Priorite : P0
Statut : Fait
Risque : Faible

### Objectif
Fournir un jeu de donnees synthetique, petit et versionne qui represente les familles d'entrees de TITAN: assemblies, annotation precedente pour Liftoff, RNA-seq single/paired/long, proteines, parametre EGAPx et evidences pre-calculees pour la branche Aegis.
### Constat
`data_example/` contient surtout des fichiers vides ou placeholders et ne couvre pas les contrats d'entree reels. Le dataset minimal doit etre isole de ces donnees historiques et servir a valider parsing, chemins, samplesheets et futures executions stub sans dependance Slurm/conteneur.
### Fichiers concernes
`test-data/minimal/`, `scripts/validate_minimal_test_data.py`, `docs/development/audit.md`.
### Etapes d'implementation
Structurer `test-data/minimal/valid` avec assemblies FASTA, GFF3 de reference, RNA-seq FASTQ gzip single/paired/long, deux FASTA proteiques, samplesheets TITAN, `input_egapx.yaml` et fixtures d'evidences EDTA/Liftoff/BRAKER/StringTie/PsiCLASS/Minimap2. Conserver des cas invalides courts pour les validations negatives. Ajouter un validateur statique et regenerer les checksums.
### Tests
`python3 scripts/validate_minimal_test_data.py`; `sha256sum -c test-data/minimal/checksums.sha256`; `nextflow config -profile test`; `nextflow run main.nf -profile test --workflow aegis -ansi-log false`.
### Criteres d'acceptation
Fixtures versionnees, petites, documentees et verifiees par checksums. Le samplesheet RNA-seq couvre `single`, `paired` et `long`; le samplesheet proteique couvre plusieurs sources; les evidences Aegis/Liftoff existent sous forme minimale pour tests futurs.
### Risques et retour arriere
Risque faible; les donnees sont synthetiques et isolees.

## TITAN-P0-005 - Controler les parametres obligatoires
Priorite : P0
Statut : Fait
Risque : Moyen

### Objectif
Valider explicitement les parametres avant de construire les channels.
### Constat
La validation est maintenant executee dans `workflow`, avant la construction des channels. Les parametres obligatoires, les fichiers d'entree absents et les noms de workflow invalides produisent des erreurs lisibles avant tout process lourd.
### Fichiers concernes
`main.nf`, `docs/development/audit.md`.
### Etapes d'implementation
Ajouter des helpers de validation dans `main.nf`, verifier les parametres obligatoires, verifier l'existence des fichiers d'entree et utiliser `Channel.fromPath(..., checkIfExists: true)` pour les samplesheets de la branche `generate_evidence_data`.
### Tests
`nextflow config -profile test`; `nextflow run main.nf -profile test --workflow aegis -ansi-log false`; cas negatifs avec parametre vide, fichier absent et workflow invalide.
### Criteres d'acceptation
Erreur lisible avant calcul lourd pour parametre manquant, fichier absent et workflow invalide.
### Risques et retour arriere
Peut modifier le moment d'echec; documenter.

## TITAN-P0-006 - Auditer l'architecture Nextflow cible
Priorite : P0
Statut : Fait
Risque : Faible

### Objectif
Identifier les defauts structurels qui rendent TITAN fragile avant d'ajouter de nouvelles fonctionnalites.
### Constat
L'audit confirme que le probleme principal est le contrat entre couches: l'orchestration est maintenant dans `workflows/titan.nf`, mais Aegis relit encore les fichiers publies dans `output_dir`, EDTA depend encore d'un nom publie pour Aegis-only, et EGAPx doit encore exposer des sorties nommees pour Aegis.
### Fichiers concernes
`docs/development/architecture-audit.md`, `roadmap.md`.
### Etapes d'implementation
Inspecter `main.nf`, `subworkflows/*.nf`, `modules/*.nf`, `conf/*.config` et formaliser les risques et l'architecture cible.
### Tests
Audit statique par `rg`, lecture des modules critiques, `nextflow config -profile test`.
### Criteres d'acceptation
Document d'audit present avec constats, impacts, recommandations et ordre de migration.
### Risques et retour arriere
Aucun impact runtime.

## TITAN-P1-001 - Clarifier les modes workflow et normaliser les booleens
Priorite : P1
Statut : Fait
Risque : Moyen

### Objectif
Rendre les modes `generate_evidence_data`, `aegis` et futur `all` explicites. Cette etape a d'abord normalise les booleens historiques; le contrat a ensuite evolue: EDTA et EGAPx sont obligatoires et les long reads sont detectes depuis la samplesheet.
### Constat
`main.nf` rejette maintenant `all` explicitement tant que ce mode n'est pas implemente. Les flags biologiques `EDTA`, `run_edta`, `run_egapx` et `use_long_reads` ne pilotent plus l'execution: EDTA/EGAPx sont obligatoires dans `generate_evidence_data`, et `library_layout=long` active automatiquement la branche long reads.
### Fichiers concernes
`main.nf`, `nextflow.config`, `conf/test.config`, `subworkflows/*.nf`, documentation.
### Etapes d'implementation
Ajouter un bloc de validation de parametres, refuser `all` tant qu'il n'est pas implemente, supprimer les booleens biologiques du contrat runtime et passer la detection long reads aux sous-workflows.
### Tests
`nextflow config -profile test`; `nextflow run main.nf -profile test --workflow aegis` doit echouer clairement si les evidences Aegis obligatoires manquent; `python3 scripts/validate_minimal_test_data.py`; test negatif `--workflow nope`.
### Criteres d'acceptation
Chaque mode a un comportement documente et aucun module ne depend d'une chaine booleenne ambigue.
### Risques et retour arriere
Peut changer le comportement historique de configs qui passaient des chaines; documenter les alias conserves.

## TITAN-P1-002 - Introduire une architecture Nextflow standard
Priorite : P1
Statut : Fait
Risque : Moyen

### Objectif
Faire de `main.nf` un point d'entree leger et deplacer l'orchestration vers une couche `workflows/` plus testable.
### Constat
`main.nf` est maintenant un point d'entree leger. L'orchestration historique a ete deplacee dans `workflows/titan.nf`; elle valide les parametres, construit les channels, scanne encore `output_dir` pour le mode Aegis-only et appelle les sous-workflows.
### Fichiers concernes
`main.nf`, `workflows/titan.nf`, `subworkflows/`, `docs/development/architecture-audit.md`.
### Etapes d'implementation
Creer `workflows/titan.nf`, y deplacer validation/orchestration, garder `main.nf` limite aux defaults publics, a l'include et a l'appel du workflow principal, conserver les noms publics de parametres.
### Tests
`nextflow config -profile test`; `python3 scripts/validate_minimal_test_data.py`; `nextflow run main.nf -profile test --workflow generate_evidence_data -stub-run -ansi-log false`; `nextflow run main.nf -profile test --workflow aegis -stub-run -ansi-log false`.
### Criteres d'acceptation
`main.nf` ne contient plus de logique biologique ni de scan de fichiers; les tests stub restent verts et les process sont names sous le workflow `TITAN`.
### Risques et retour arriere
Migration par etapes; commit dedie pour pouvoir revert facilement.

## TITAN-P1-003 - Definir un contrat d'evidences nomme
Priorite : P1
Statut : Fait
Risque : Eleve

### Objectif
Remplacer le channel mixte `[string_key, file]` et les listes Groovy par des emits nommes et documentes.
### Constat
`generate_evidence_data` expose maintenant des emits nommes pour les familles d'evidences. `aegis.nf` prend des evidences explicites au lieu d'une liste `[string_key, file]`. Le mode Aegis-only relit encore les fichiers publics depuis `output_dir`, mais les passe ensuite au sous-workflow via un contrat nomme.
### Fichiers concernes
`subworkflows/generate_evidence_data.nf`, `subworkflows/aegis.nf`, `docs/development/architecture-audit.md`, tests.
### Etapes d'implementation
Emettre des channels nommes: Liftoff, EGAPx obligatoire, EDTA masked genome, BRAKER/AUGUSTUS, GeneMark, STAR/StringTie, STAR/PsiCLASS, long reads detectes automatiquement. Documenter les evidences requises et optionnelles.
### Tests
`nextflow config -profile test`; `python3 scripts/validate_minimal_test_data.py`; `nextflow run main.nf -profile test --workflow generate_evidence_data -stub-run -ansi-log false`; `nextflow run main.nf -profile test --workflow aegis -stub-run -ansi-log false`; cas Aegis-only sans evidences publiees qui echoue avec la liste des fichiers requis.
### Criteres d'acceptation
Aegis consomme des evidences nommees et echoue clairement si une evidence requise manque.
### Risques et retour arriere
Risque eleve car c'est le coeur du pipeline; ne pas changer les commandes scientifiques dans le meme commit.

## TITAN-P1-004 - Refaire le lien EDTA -> Aegis
Priorite : P1
Statut : Fait
Risque : Eleve

### Objectif
Modeliser explicitement le genome hard-masked requis par Aegis.
### Constat
EDTA n'est plus optionnel dans `generate_evidence_data`. Le mode `all` execute maintenant evidence generation puis Aegis dans le meme graphe Nextflow: le `masked_genome` emis par EDTA passe par le contrat nomme de `generate_evidence_data` et devient l'input explicite `masked_genome` du sous-workflow Aegis. Le mode Aegis-only reste compatible avec les runs separes en relisant `assembly_masked.EDTA.fasta` depuis `output_dir` et echoue clairement si cette evidence manque.
### Fichiers concernes
`workflows/titan.nf`, `subworkflows/generate_evidence_data.nf`, `subworkflows/aegis.nf`, `modules/EDTA.nf`, `modules/aegis_*.nf`, documentation.
### Etapes d'implementation
Supprimer les flags `EDTA`/`run_edta` du contrat runtime, lancer EDTA systematiquement dans `generate_evidence_data`, ajouter un vrai mode `all`, transmettre `evidence_data.masked_genome` directement a `aegis`, et faire echouer Aegis-only si le genome hard-masked requis n'est pas disponible.
### Tests
`nextflow config -profile test`; `python3 scripts/validate_minimal_test_data.py`; `nextflow run main.nf -profile test --workflow all -stub-run -ansi-log false`; `nextflow run main.nf -profile test --workflow generate_evidence_data -stub-run -ansi-log false`; `nextflow run main.nf -profile test --workflow aegis -stub-run -ansi-log false`.
### Criteres d'acceptation
Aegis ne peut plus etre saute par `EDTA=no`; son input hard-masked est obligatoire. En mode `all`, Aegis consomme la sortie EDTA par channel Nextflow et ne depend pas d'un scan de `output_dir`.
### Risques et retour arriere
Risque scientifique et UX; conserver le mode Aegis-only par fichiers publies pendant une version, puis le remplacer par un manifeste d'evidences.

## TITAN-P1-005 - Integrer EGAPx proprement
Priorite : P1
Statut : En cours
Risque : Eleve

### Objectif
Brancher EGAPx comme source d'annotation obligatoire et testable, puis integrer ses sorties nommees au contrat d'evidences Aegis.
### Constat
`modules/egapx.nf` est maintenant appele par `generate_evidence_data` avec un `path egapx_paramfile`, publie `${output_dir}/egapx` et fournit un stub minimal. Ses sorties restent larges (`egapx.results`) tant que les noms exacts des fichiers EGAPx a consommer par Aegis ne sont pas stabilises.
### Fichiers concernes
`modules/egapx.nf`, `subworkflows/generate_evidence_data.nf`, `subworkflows/aegis.nf`, `modules/aegis_*.nf`, `test-data/minimal/valid/evidence/`.
### Etapes d'implementation
Supprimer `run_egapx`, corriger le module pour prendre `path egapx_paramfile`, lancer EGAPx systematiquement dans `generate_evidence_data`, ajouter un stub EGAPx, puis identifier et emettre les vrais outputs EGAPx (`egapx_gff3`, proteines, rapports) pour les connecter a Aegis.
### Tests
`nextflow config -profile test`; `nextflow run main.nf -profile test --workflow generate_evidence_data -stub-run -ansi-log false`; `nextflow run main.nf -profile test --workflow aegis -stub-run -ansi-log false`; test d'integration Aegis avec fixture EGAPx nommee.
### Criteres d'acceptation
EGAPx n'est plus activable/desactivable par parametre, produit des emits nommes et peut etre consomme explicitement par Aegis.
### Risques et retour arriere
EGAPx embarque son propre Nextflow et peut etre lourd; continuer avec stub/fixture et stabiliser les outputs nommes avant de les faire consommer par Aegis.

## TITAN-P1-006 - Supprimer les scans de dossiers internes
Priorite : P1
Statut : A faire
Risque : Eleve

### Objectif
Faire utiliser aux modules leurs inputs declares au lieu de scanner `output_dir`, `work` ou `data`.
### Constat
Les modules StringTie merge, GFFCompare, BRAKER3, read preparation et Aegis utilisent des mounts et scripts de recuperation de chemins qui contournent les channels Nextflow.
### Fichiers concernes
`modules/Stringtie_merging_*.nf`, `modules/gffcompare.nf`, `modules/braker3_prediction*.nf`, `modules/prepare_RNAseq_fastq_files_*.nf`, `modules/trimming_fastq.nf`, scripts `retrieve_*`.
### Etapes d'implementation
Migrer module par module vers des `path` inputs stages, remplacer les scripts de scan par des listes de fichiers passees en channels, conserver les chemins historiques uniquement en sortie publique.
### Tests
nf-test ou `-stub-run` par module, comparaison des noms de fichiers publies.
### Criteres d'acceptation
Les process ne lisent plus d'inputs internes depuis `${params.output_dir}` ou `${projectDir}/data` sauf parametre explicite.
### Risques et retour arriere
Risque eleve; proceder par familles de modules avec fixtures dediees.

## TITAN-P1-007 - Normaliser les outputs et la provenance
Priorite : P1
Statut : A faire
Risque : Moyen

### Objectif
Rendre les sorties previsibles et collecter versions/manifestes.
### Constat
Plusieurs modules copient manuellement vers `/outputdir` en plus de `publishDir`; aucun `versions.yml` ou manifeste d'evidences n'est produit.
### Fichiers concernes
`modules/*.nf`, `scripts/`, `docs/`.
### Etapes d'implementation
Publier uniquement des outputs declares, ajouter un `evidence_manifest.json`, ajouter des emits `versions`, conserver les noms historiques via `publishDir saveAs` si necessaire.
### Tests
Verification presence des fichiers historiques et nouveaux manifestes.
### Criteres d'acceptation
Pas de regression de noms publics; provenance lisible apres chaque run.
### Risques et retour arriere
Risque de rupture de scripts utilisateurs; maintenir une couche de compatibilite.

## TITAN-P1-008 - Verrouiller les conteneurs critiques
Priorite : P1
Statut : A faire
Risque : Moyen

### Objectif
Remplacer progressivement les tags `latest` par des tags versionnes compatibles.
### Constat
BRAKER3, EDTA, Diamond2GO, GFFCompare, HISAT2, Minimap2, PsiCLASS et StringTie utilisent `latest`.
### Fichiers concernes
`modules/*.nf`, `dockerfiles/`, documentation.
### Etapes d'implementation
Verifier versions, licences, commandes et resultats avant chaque remplacement; ne pas remplacer toutes les images dans un seul commit.
### Tests
Test module par module avec petites donnees ou stub.
### Criteres d'acceptation
Versions collectees et images documentees.
### Risques et retour arriere
Risque scientifique; remplacement progressif seulement.

## TITAN-P1-009 - Finaliser profils local, apptainer et slurm
Priorite : P1
Statut : En cours
Risque : Moyen

### Objectif
Isoler Docker, Apptainer, Slurm et ressources locales.
### Constat
Les profils existent et resolvent, mais de nombreux `containerOptions --volume` restent Docker-centriques et les ressources sont encore partiellement hardcodees dans les modules.
### Fichiers concernes
`nextflow.config`, `conf/*.config`, `modules/*.nf`.
### Etapes d'implementation
Valider Apptainer/Slurm avec un stub, retirer les mounts inutiles au profit du staging Nextflow, definir des labels par classe de process.
### Tests
`nextflow config -profile local`; `test`; `slurm,apptainer,test`; stub Apptainer si environnement disponible.
### Criteres d'acceptation
Les profils resolvent et les modules critiques ne dependent pas de mounts Docker evitables.
### Risques et retour arriere
Risque HPC; ne pas changer simultanement profils et commandes scientifiques.

## TITAN-P2-001 - Ajouter schema et validation avancee des entrees
Priorite : P2
Statut : A faire
Risque : Moyen

### Objectif
Valider FASTA, GFF3, samplesheets, chemins proteiques, options et coherence seqid avant calcul lourd.
### Constat
P0-005 valide presence des parametres et existence des fichiers, pas encore le schema complet ni la coherence biologique de base.
### Fichiers concernes
`scripts/validate_minimal_test_data.py`, futur `scripts/validate_inputs.py`, `test-data/minimal/invalid/`.
### Etapes d'implementation
Ajouter validation CSV stricte, presence FASTQ/proteines, seqids GFF3/FASTA, coordonnees, Parent, options enum.
### Tests
Fixtures invalides existantes plus nouveaux cas samplesheets.
### Criteres d'acceptation
Chaque entree invalide echoue avant calcul lourd avec message actionnable.
### Risques et retour arriere
Risque de faux positifs sur donnees historiques; commencer en mode warning si necessaire.

## TITAN-P2-002 - Mettre en place la strategie de tests
Priorite : P2
Statut : A faire
Risque : Moyen

### Objectif
Ajouter tests statiques, unitaires, integration et bout en bout minimal.
### Constat
Aucun repertoire `tests/` TITAN dedie.
### Fichiers concernes
`tests/`, `scripts/run-tests.sh`, `nf-test.config`.
### Etapes d'implementation
Commencer par checks config et validations Python; ajouter nf-test apres stabilisation des contrats; couvrir les sous-workflows evidence et Aegis en stub.
### Tests
Pytest ou Bash statique; nf-test module cible; `-stub-run`.
### Criteres d'acceptation
Une commande unique de tests rapides existe et tourne en local.
### Risques et retour arriere
Risque de tests fragiles si les contrats ne sont pas d'abord clarifies.

## TITAN-P2-003 - CI GitHub Actions rapide
Priorite : P2
Statut : A faire
Risque : Moyen

### Objectif
Verifier syntaxe, config et tests rapides sans gros calcul.
### Constat
Pas de CI TITAN dediee detectee.
### Fichiers concernes
`.github/workflows/titan-ci.yml`.
### Etapes d'implementation
Installer Nextflow, lancer checks statiques, validateur fixtures et profil test/stub.
### Tests
Execution locale des commandes CI.
### Criteres d'acceptation
CI sans secrets ni donnees lourdes.
### Risques et retour arriere
Apptainer peut etre limite dans GitHub Actions; documenter.

## TITAN-P3-001 - Refonte documentaire utilisateur
Priorite : P3
Statut : A faire
Risque : Faible

### Objectif
README complet: quick start, entrees, sorties, profils, exemples, troubleshooting.
### Constat
README utile mais encore insuffisant pour expliquer les modes production, les evidences Aegis-only et les limites actuelles.
### Fichiers concernes
`README.md`, `docs/`.
### Etapes d'implementation
Documenter modes workflow, contrats samplesheets, outputs publics, reprise Aegis via manifeste, EGAPx, EDTA/masked assembly et troubleshooting.
### Tests
Verifier chaque commande documentee.
### Criteres d'acceptation
Un nouvel utilisateur peut lancer un test local et comprendre comment preparer un run production.
### Risques et retour arriere
Faible.

## TITAN-P4-001 - Optimiser ressources HPC
Priorite : P4
Statut : A faire
Risque : Moyen

### Objectif
Definir CPU, memoire, temps, retries et concurrence par type de processus.
### Constat
Ressources hardcodees dans modules et config globale `100GB/20 CPU`.
### Fichiers concernes
`conf/base.config`, `modules/*.nf`.
### Etapes d'implementation
Introduire labels (`process_index`, `process_alignment`, `process_transcriptome`, `process_prediction`, `process_merge`, `process_aegis`, `process_low`) et les appliquer progressivement.
### Tests
`nextflow config -flat`; test local; validation Slurm sur petit run/stub.
### Criteres d'acceptation
Ressources adaptables sans modifier les modules.
### Risques et retour arriere
Risque de sous-dimensionnement; tester sur petits jeux puis production.

## TITAN-P5-001 - Validation scientifique des annotations
Priorite : P5
Statut : A faire
Risque : Eleve

### Objectif
Valider GFF3 final, FASTA/proteines et coherence biologique.
### Constat
Aucune validation finale automatisee detectee.
### Fichiers concernes
`scripts/`, `modules/local/validate_annotation`, `tests/`.
### Etapes d'implementation
Verifier GFF3, Parent, phases CDS, coordonnees, seqids, duplications, sequences absentes, statistiques d'annotation.
### Tests
Fixtures valides/invalides et comparaison aux resultats historiques.
### Criteres d'acceptation
Rapport de validation final present et bloque les erreurs critiques.
### Risques et retour arriere
Risque de faux positifs; commencer en mode rapport avant mode bloquant.
