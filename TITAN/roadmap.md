# Roadmap TITAN

Cette roadmap est specifique a TITAN. Elle transpose les bonnes pratiques observees dans `/home/vmadmin/atcg-rnaseq` vers un workflow d'annotation de genomes, sans transformer TITAN en pipeline RNA-seq.

## TITAN-P0-001 - Etablir une baseline fonctionnelle
Priorite : P0
Statut : Fait
Risque : Faible

### Objectif
Documenter l'etat initial, les commandes historiques et les erreurs bloquantes avant refactor.
### Constat
`nextflow config` fonctionne. Le lancement minimal echouait initialement a la compilation DSL2; une correction P0 permet maintenant une commande minimale `aegis` avec `EDTA=no`.
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
Statut : A faire
Risque : Moyen

### Objectif
Valider explicitement les parametres avant de construire les channels.
### Constat
La validation actuelle est top-level et casse la compilation avec Nextflow 26.
### Fichiers concernes
`main.nf`, futurs scripts de validation.
### Etapes d'implementation
Deplacer la validation dans `workflow`, utiliser `Channel.fromPath(..., checkIfExists: true)` pour les branches concernees.
### Tests
Cas nominal, parametre manquant, fichier absent.
### Criteres d'acceptation
Erreur lisible avant calcul lourd.
### Risques et retour arriere
Peut modifier le moment d'echec; documenter.

## TITAN-P1-001 - Consolider l'architecture DSL2
Priorite : P1
Statut : A faire
Risque : Moyen

### Objectif
Separer workflow principal, sous-workflows et modules locaux avec contrats d'entree/sortie.
### Constat
Modules presents mais channels et listes Groovy restent fragiles.
### Fichiers concernes
`main.nf`, `workflows/titan.nf`, `subworkflows/`, `modules/`.
### Etapes d'implementation
Introduire `workflows/titan.nf`, deplacer progressivement l'orchestration, documenter les tuples.
### Tests
Compilation, tests de sous-workflows, `-resume`.
### Criteres d'acceptation
`main.nf` devient un point d'entree leger.
### Risques et retour arriere
Migration par etapes seulement.

## TITAN-P1-002 - Verrouiller les conteneurs critiques
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
Verifier versions, licences, commandes et resultats avant chaque remplacement.
### Tests
Test module par module avec petites donnees ou stub.
### Criteres d'acceptation
Versions collectees et images documentees.
### Risques et retour arriere
Risque scientifique; aucun remplacement global.

## TITAN-P1-003 - Ajouter profils local, apptainer et slurm
Priorite : P1
Statut : En cours
Risque : Moyen

### Objectif
Isoler Docker, Apptainer, Slurm et ressources locales.
### Constat
Docker est active globalement; pas de configuration Slurm/Apptainer separee.
### Fichiers concernes
`nextflow.config`, `conf/*.config`.
### Etapes d'implementation
Ajouter includes par profil, ressources par labels, variables Slurm optionnelles.
### Tests
`nextflow config -profile local`; `test`; `slurm,apptainer,test`.
### Criteres d'acceptation
Les profils resolvent sans soumission cluster.
### Risques et retour arriere
Risque de changement runtime si Docker global est modifie; conserver comportement standard initial.

## TITAN-P1-004 - Normaliser les outputs et la provenance
Priorite : P1
Statut : A faire
Risque : Moyen

### Objectif
Rendre les sorties previsibles et collecter versions/manifestes.
### Constat
Plusieurs modules copient manuellement vers `/outputdir` en plus de `publishDir`.
### Fichiers concernes
`modules/*.nf`, `scripts/`, `docs/`.
### Etapes d'implementation
Centraliser publication via `publishDir`, emettre versions, conserver compatibilite des noms historiques.
### Tests
Verification presence des fichiers historiques et nouveaux manifestes.
### Criteres d'acceptation
Pas de regression de noms publics.
### Risques et retour arriere
Risque eleve si fait trop vite; proceder module par module.

## TITAN-P2-001 - Mettre en place la strategie de tests
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
Commencer par checks config et tests de validation FASTA/GFF3; evaluer nf-test apres stabilisation des modules.
### Tests
Pytest ou Bash statique; nf-test module cible; `-stub-run`.
### Criteres d'acceptation
Une commande unique de tests rapides existe.
### Risques et retour arriere
Risque de tests fragiles si les contrats ne sont pas d'abord clarifies.

## TITAN-P2-002 - Ajouter tests negatifs d'entrees
Priorite : P2
Statut : A faire
Risque : Faible

### Objectif
Verifier FASTA invalide, GFF3 invalide, seqid incompatible, coordonnees invalides et Parent invalide.
### Constat
Ces erreurs ne sont pas controlees explicitement.
### Fichiers concernes
`test-data/minimal/invalid/`, futur validateur.
### Etapes d'implementation
Creer fixtures invalides et assertions d'erreur lisible.
### Tests
Tests unitaires du validateur.
### Criteres d'acceptation
Chaque entree invalide echoue avant calcul lourd.
### Risques et retour arriere
Faible.

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
Installer Nextflow, lancer checks statiques et profil test/stub.
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
README utile mais incomplet sur tests, profils et limitations.
### Fichiers concernes
`README.md`, `docs/`.
### Etapes d'implementation
Documenter commandes historiques et nouvelles commandes test.
### Tests
Verifier chaque commande documentee.
### Criteres d'acceptation
Un nouvel utilisateur peut lancer un test local.
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
Introduire labels (`process_low`, `process_alignment`, `process_prediction`, `process_merge`, `process_aegis`) et les appliquer progressivement.
### Tests
`nextflow config -flat`; test local.
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
