# Roadmap TITAN

Cette roadmap est specifique a TITAN. Elle transpose les bonnes pratiques observees dans `/home/vmadmin/atcg-rnaseq` vers un workflow d'annotation de genomes, sans transformer TITAN en pipeline RNA-seq.

## TITAN-P0-001 - Etablir une baseline fonctionnelle
Priorite : P0
Statut : Fait
Risque : Faible

### Objectif
Documenter l'etat initial, les commandes historiques et les erreurs bloquantes avant refactor.
### Constat
`nextflow config` fonctionne. Le lancement minimal echouait initialement a la compilation DSL2; une correction P0 a permis une commande minimale. Depuis le durcissement P1, TITAN lance obligatoirement la generation d'evidences puis Aegis dans le meme graphe.
### Fichiers concernes
`main.nf`, `nextflow.config`, `docs/development/audit.md`.
### Etapes d'implementation
Capturer les commandes, corriger uniquement le blocage de compilation, verifier que la configuration reste resolvable.
### Tests
`nextflow config`; `nextflow run main.nf -ansi-log false` ou commande test adaptee.
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
Les versions sont documentees dans `README.md` et les images runtime sont digest-pinned depuis P1-008. L'inventaire P0-002 est maintenant consolide dans `docs/development/audit.md` et resume dans `README.md`.
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
`nextflow config -profile test`; `nextflow run main.nf -profile test -stub-run -ansi-log false`; `nextflow run main.nf -profile test -stub-run -ansi-log false -resume`.
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
`python3 scripts/validate_minimal_test_data.py`; `sha256sum -c test-data/minimal/checksums.sha256`; `nextflow config -profile test`; `nextflow run main.nf -profile test -stub-run -ansi-log false`.
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
Ajouter des helpers de validation dans `workflows/titan.nf`, verifier les parametres obligatoires, verifier l'existence des fichiers d'entree et utiliser `Channel.fromPath(..., checkIfExists: true)` pour les samplesheets.
### Tests
`nextflow config -profile test`; `nextflow run main.nf -profile test -stub-run -ansi-log false`; cas negatifs avec parametre vide et fichier absent.
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
L'audit confirme que le probleme principal est le contrat entre couches: l'orchestration est maintenant dans `workflows/titan.nf`, Aegis consomme les evidences nommees directement, et EGAPx expose des sorties nommees dont le GFF3 est transmis a AEGIS comme preuve supplementaire.
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

## TITAN-P1-001 - Clarifier le contrat workflow et normaliser les booleens
Priorite : P1
Statut : Fait
Risque : Moyen

### Objectif
Supprimer les modes partiels comme contrat utilisateur et normaliser les booleens historiques. TITAN lance obligatoirement la generation d'evidences puis Aegis; EDTA et EGAPx sont obligatoires et les long reads sont detectes depuis la samplesheet.
### Constat
`params.workflow` n'est plus un parametre public. Les flags biologiques `EDTA`, `run_edta`, `run_egapx` et `use_long_reads` ne pilotent plus l'execution: EDTA/EGAPx sont obligatoires dans TITAN, et `library_layout=long` active automatiquement la branche long reads.
### Fichiers concernes
`main.nf`, `nextflow.config`, `conf/test.config`, `subworkflows/*.nf`, documentation.
### Etapes d'implementation
Ajouter un bloc de validation de parametres, supprimer le routage par `params.workflow`, supprimer les booleens biologiques du contrat runtime et passer la detection long reads aux sous-workflows.
### Tests
`nextflow config -profile test`; `nextflow run main.nf -profile test -stub-run -ansi-log false`; `python3 scripts/validate_minimal_test_data.py`.
### Criteres d'acceptation
Le pipeline a un seul comportement documente et aucun module ne depend d'une chaine booleenne ambigue.
### Risques et retour arriere
Peut changer le comportement historique de configs qui passaient des chaines; documenter les alias conserves.

## TITAN-P1-002 - Introduire une architecture Nextflow standard
Priorite : P1
Statut : Fait
Risque : Moyen

### Objectif
Faire de `main.nf` un point d'entree leger et deplacer l'orchestration vers une couche `workflows/` plus testable.
### Constat
`main.nf` est maintenant un point d'entree leger. L'orchestration historique a ete deplacee dans `workflows/titan.nf`; elle valide les parametres, construit les channels et appelle obligatoirement les sous-workflows evidence generation puis Aegis.
### Fichiers concernes
`main.nf`, `workflows/titan.nf`, `subworkflows/`, `docs/development/architecture-audit.md`.
### Etapes d'implementation
Creer `workflows/titan.nf`, y deplacer validation/orchestration, garder `main.nf` limite aux defaults publics, a l'include et a l'appel du workflow principal, conserver les noms publics de parametres.
### Tests
`nextflow config -profile test`; `python3 scripts/validate_minimal_test_data.py`; `nextflow run main.nf -profile test -stub-run -ansi-log false`.
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
`generate_evidence_data` expose maintenant des emits nommes pour les familles d'evidences. `aegis.nf` prend des evidences explicites au lieu d'une liste `[string_key, file]`. TITAN connecte ces emits directement a Aegis dans un seul graphe.
### Fichiers concernes
`subworkflows/generate_evidence_data.nf`, `subworkflows/aegis.nf`, `docs/development/architecture-audit.md`, tests.
### Etapes d'implementation
Emettre des channels nommes: Liftoff, EGAPx obligatoire, EDTA masked genome, BRAKER/AUGUSTUS, GeneMark, STAR/StringTie, STAR/PsiCLASS, long reads detectes automatiquement. Documenter les evidences requises et optionnelles.
### Tests
`nextflow config -profile test`; `python3 scripts/validate_minimal_test_data.py`; `nextflow run main.nf -profile test -stub-run -ansi-log false`.
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
EDTA n'est plus optionnel. TITAN execute evidence generation puis Aegis dans le meme graphe Nextflow: le `masked_genome` emis par EDTA passe par le contrat nomme de `generate_evidence_data` et devient l'input explicite `masked_genome` du sous-workflow Aegis. Les modes partiels `generate_evidence_data`, `aegis` et `all` ont ete retires du contrat public.
### Fichiers concernes
`workflows/titan.nf`, `subworkflows/generate_evidence_data.nf`, `subworkflows/aegis.nf`, `modules/EDTA.nf`, `modules/aegis_*.nf`, documentation.
### Etapes d'implementation
Supprimer les flags `EDTA`/`run_edta` du contrat runtime, lancer EDTA systematiquement, supprimer le routage par `params.workflow`, transmettre `evidence_data.masked_genome` directement a `aegis`.
### Tests
`nextflow config -profile test`; `python3 scripts/validate_minimal_test_data.py`; `nextflow run main.nf -profile test -stub-run -ansi-log false`.
### Criteres d'acceptation
Aegis ne peut plus etre saute par `EDTA=no`; son input hard-masked est obligatoire. Aegis consomme la sortie EDTA par channel Nextflow et ne depend pas d'un scan de `output_dir`.
### Risques et retour arriere
Risque scientifique et UX; les reruns partiels devront etre reintroduits plus tard via un contrat explicite de reprise ou un manifeste d'evidences, pas par un mode Aegis-only implicite.

## TITAN-P1-005 - Integrer EGAPx proprement
Priorite : P1
Statut : Fait
Risque : Eleve

### Objectif
Brancher EGAPx comme source d'annotation obligatoire et testable, puis integrer ses sorties nommees au contrat d'evidences TITAN.
### Constat
`modules/egapx.nf` utilise maintenant le runner officiel EGAPx `v0.5.2`, pilote l'image Docker officielle `ncbi/egapx` verrouillee par digest, publie `${output_dir}/egapx` et emet des outputs nommes: GFF3, GTF, proteines, CDS, transcrits, ASN, repertoire complet et versions. Les sorties EGAPx sont typees dans `generate_evidence_data`; `egapx_gff3` est transmis a AEGIS et fusionne comme preuve d'annotation supplementaire.
### Fichiers concernes
`modules/egapx.nf`, `subworkflows/generate_evidence_data.nf`, `subworkflows/aegis.nf`, `modules/aegis_*.nf`, `test-data/minimal/valid/evidence/`.
### Etapes d'implementation
Supprimer `run_egapx`, corriger le module pour prendre `path egapx_paramfile`, lancer EGAPx systematiquement, remplacer l'image custom par l'image NCBI versionnee, telecharger/executer le runner officiel versionne, ajouter un stub EGAPx complet, puis emettre les vrais outputs EGAPx (`egapx_gff3`, `egapx_gtf`, proteines, CDS, transcrits, ASN, repertoire complet, versions).
### Tests
`docker manifest inspect --verbose ncbi/egapx@sha256:bc657b232d93364d5f3b75ad3bfaf14b6267e46173672b609f26078d48a04298`; `nextflow config -profile test`; `nextflow run main.nf -profile test -stub-run -ansi-log false`.
### Criteres d'acceptation
EGAPx n'est plus activable/desactivable par parametre, utilise une image officielle versionnee, produit des emits nommes et fournit un GFF3 consomme par AEGIS.
### Risques et retour arriere
EGAPx embarque son propre Nextflow et peut etre lourd; le module reste isole mais son GFF3 fait maintenant partie du contrat AEGIS.

## TITAN-P1-006 - Supprimer les scans de dossiers internes
Priorite : P1
Statut : Fait
Risque : Eleve

### Objectif
Faire utiliser aux modules leurs inputs declares au lieu de scanner `output_dir`, `work` ou `data`.
### Constat
Les modules StringTie merge, GFFCompare, BRAKER3, read preparation, trimming, AGAT, Salmon, Minimap2, STAR et HISAT2 utilisaient des mounts ou scripts de recuperation de chemins qui contournaient les channels Nextflow. Ils consomment maintenant les fichiers declares par `path` quand ils dependent d'outputs internes.
### Fichiers concernes
`modules/Stringtie_merging_*.nf`, `modules/gffcompare.nf`, `modules/braker3_prediction*.nf`, `modules/prepare_RNAseq_fastq_files_*.nf`, `modules/trimming_fastq.nf`; anciens scripts `retrieve_*` supprimes apres migration.
### Etapes d'implementation
Migrer module par module vers des `path` inputs stages, remplacer les scripts de scan par des listes de fichiers passees en channels, conserver les chemins historiques uniquement en sortie publique. Les preuves unstranded optionnelles sont representees par des fichiers vides et ignorees par AEGIS si absentes.
### Tests
`nextflow config -profile test`; `python3 scripts/validate_minimal_test_data.py`; `nextflow run main.nf -profile test -stub-run -ansi-log false`; verification des outputs publies AEGIS/EGAPx/StringTie/GFFCompare.
### Criteres d'acceptation
Les process critiques ne lisent plus d'inputs internes depuis `${params.output_dir}`, `${projectDir}/work` ou `${projectDir}/data` pour retrouver les outputs amont. Les fichiers locaux RNA-seq restent lus via le parametre explicite `RNAseq_data_dir`.
### Risques et retour arriere
Risque eleve; proceder par familles de modules avec fixtures dediees.

## TITAN-P1-007 - Normaliser les outputs et la provenance
Priorite : P1
Statut : Fait
Risque : Moyen

### Objectif
Rendre les sorties previsibles et collecter versions/manifestes.
### Constat
Les derniers outputs publics qui dependaient de copies manuelles vers `/outputdir` sont maintenant publies depuis des outputs declares. TITAN publie un manifeste d'evidences et un `versions.yml` de provenance sous `${output_dir}/provenance`.
### Fichiers concernes
`modules/*.nf`, `scripts/`, `docs/`.
### Etapes d'implementation
Publier uniquement des outputs declares, ajouter un `evidence_manifest.json`, ajouter des emits `versions`, conserver les noms historiques via `publishDir saveAs` si necessaire.
### Tests
`nextflow config -profile test`; `python3 scripts/validate_minimal_test_data.py`; `nextflow run main.nf -profile test -stub-run -ansi-log false`; validation JSON de `${output_dir}/provenance/evidence_manifest.json`; verification presence des fichiers historiques et nouveaux manifestes.
### Criteres d'acceptation
Pas de regression de noms publics; provenance lisible apres chaque run.
### Risques et retour arriere
Risque de rupture de scripts utilisateurs; maintenir une couche de compatibilite.

## TITAN-P1-008 - Verrouiller les conteneurs critiques
Priorite : P1
Statut : Fait
Risque : Moyen

### Objectif
Remplacer les images runtime flottantes par des references digest-pinned et rendre le contrat testable.
### Constat
Les modules actifs consomment maintenant des `params.container_*` centralises dans `nextflow.config`, tous verrouilles en `image@sha256`. Les `Dockerfile` locaux ont aussi des `FROM` verrouilles. La validation automatique refuse les `:latest` dans les fichiers runtime et les `FROM` non digestes.
### Fichiers concernes
`nextflow.config`, `modules/*.nf`, `dockerfiles/`, `scripts/validate_container_pins.py`, documentation.
### Etapes d'implementation
Centraliser les images, remplacer les `latest` par les digests de manifestes Docker/OCI, verrouiller les bases Dockerfile, documenter l'inventaire et ajouter un validateur.
### Tests
`python3 scripts/validate_container_pins.py`; `python3 scripts/validate_minimal_test_data.py`; `nextflow config -profile test`; `nextflow run main.nf -profile test -stub-run -ansi-log false`; validation registry des digests par `docker manifest inspect --verbose`.
### Criteres d'acceptation
Plus aucun `latest` dans les modules/config runtime; les Dockerfiles sont digest-pinned; l'inventaire est documente dans `docs/development/container-locks.md`.
### Risques et retour arriere
Les images `avelt/*` restent des images historiques sans tags semantiques publies; le pin par digest preserve le contenu actuel mais ne remplace pas une future reconstruction versionnee propre.

## TITAN-P1-009 - Finaliser profils local, apptainer et slurm
Priorite : P1
Statut : Fait
Risque : Moyen

### Objectif
Isoler Docker, Apptainer, Slurm et ressources locales.
### Constat
Les profils `local`, `test`, `apptainer` et `slurm` sont explicites. Les modules actifs ne declarent plus de `containerOptions`; les chemins genome, RNA-seq et scripts sont stages via `path`. Le profil `apptainer` force EGAPx imbrique en `singularity`, et `slurm` expose queue, compte, QOS, limite de soumission et taille de file.
### Fichiers concernes
`nextflow.config`, `conf/*.config`, `modules/*.nf`, `workflows/titan.nf`, `subworkflows/generate_evidence_data.nf`, `launch_TITAN_example.sh`, `scripts/validate_profiles.py`.
### Etapes d'implementation
Retirer les mounts Docker au profit du staging Nextflow, separer les profils runtime, adapter EGAPx a Apptainer/Singularity, durcir le lanceur production et ajouter une validation statique des profils.
### Tests
`python3 scripts/validate_profiles.py`; `nextflow config -profile local`; `nextflow config -profile apptainer`; `nextflow config -profile slurm,apptainer`; `nextflow config -profile test,slurm,apptainer`; `nextflow run main.nf -profile test -stub-run -ansi-log false`; variante stub sans long reads.
### Criteres d'acceptation
Les profils resolvent sans warnings critiques; les modules critiques ne dependent pas de mounts Docker evitables; le lanceur production prepare Apptainer et valide les contrats avant lancement.
### Risques et retour arriere
Risque HPC residuel: le vrai run Slurm/Apptainer doit etre valide sur le cluster cible avec cache Apptainer partage et images disponibles.

## TITAN-P2-001 - Ajouter schema et validation avancee des entrees
Priorite : P2
Statut : Fait
Risque : Moyen

### Objectif
Valider FASTA, GFF3, samplesheets, chemins proteiques, options et coherence seqid avant calcul lourd.
### Constat
TITAN lance maintenant `scripts/validate_inputs.py` au demarrage du workflow, avant creation des channels et avant calcul lourd. Le validateur controle FASTA, GFF3, samplesheets RNA-seq/proteines, chemins FASTQ/FASTA/proteines, YAML EGAPx minimal, enums et options numeriques.
### Fichiers concernes
`scripts/validate_inputs.py`, `scripts/test_validate_inputs.py`, `workflows/titan.nf`, `launch_TITAN_example.sh`, `test-data/minimal/invalid/`, documentation.
### Etapes d'implementation
Ajouter validation CSV stricte, presence FASTQ/proteines, seqids GFF3/FASTA, coordonnees, Parent, options enum et integration au workflow/lanceur.
### Tests
`python3 scripts/validate_inputs.py ...`; `python3 scripts/test_validate_inputs.py`; `nextflow run main.nf -profile test -stub-run -ansi-log false`; cas negatif Nextflow avec `test-data/minimal/invalid/rnaseq_bad_layout.csv`.
### Criteres d'acceptation
Chaque entree invalide echoue avant calcul lourd avec message actionnable.
### Risques et retour arriere
Risque de faux positifs sur donnees historiques; adapter explicitement le schema si un format historique valide doit etre conserve.

## TITAN-P2-002 - Mettre en place la strategie de tests
Priorite : P2
Statut : Fait
Risque : Moyen

### Objectif
Ajouter tests statiques, unitaires, integration et bout en bout minimal.
### Constat
Une commande unique `scripts/run-tests.sh` orchestre les tests rapides locaux: validations statiques, profils, fixtures, cas unitaires du validateur, run Nextflow stub et cas negatif d'entree invalide.
### Fichiers concernes
`tests/`, `scripts/run-tests.sh`, `nf-test.config`.
### Etapes d'implementation
Commencer par checks config et validations Python; ajouter nf-test apres stabilisation des contrats; couvrir les sous-workflows evidence et Aegis en stub. La premiere passe ajoute le runner rapide et prepare la configuration nf-test sans rendre nf-test obligatoire.
### Tests
`scripts/run-tests.sh`; `nextflow run main.nf -profile test -stub-run -ansi-log false`; cas negatif Nextflow avec `test-data/minimal/invalid/rnaseq_bad_layout.csv`.
### Criteres d'acceptation
Une commande unique de tests rapides existe et tourne en local.
### Risques et retour arriere
Risque de tests fragiles si les contrats ne sont pas d'abord clarifies.

## TITAN-P2-003 - CI GitHub Actions rapide
Priorite : P2
Statut : Fait
Risque : Moyen

### Objectif
Verifier syntaxe, config et tests rapides sans gros calcul.
### Constat
Une CI GitHub Actions dediee lance les tests rapides TITAN sur push et pull request quand `TITAN/**` ou le workflow CI changent.
### Fichiers concernes
`.github/workflows/titan-ci.yml`.
### Etapes d'implementation
Installer Java et Nextflow 24.04.3, lancer `scripts/run-tests.sh` depuis le repertoire `TITAN`.
### Tests
`scripts/run-tests.sh`; verification statique du workflow GitHub Actions.
### Criteres d'acceptation
CI sans secrets ni donnees lourdes.
### Risques et retour arriere
Apptainer peut etre limite dans GitHub Actions; la CI rapide n'execute pas Slurm, Apptainer ni les conteneurs scientifiques.

## TITAN-P3-001 - Refonte documentaire utilisateur
Priorite : P3
Statut : Fait
Risque : Faible

### Objectif
README complet: quick start, entrees, sorties, profils, exemples, troubleshooting.
### Constat
README restructure comme guide utilisateur: contrat workflow unique, quick start, entrees, samplesheets, profils, outputs, reprise, validation, CI, troubleshooting et limites.
### Fichiers concernes
`README.md`, `docs/`.
### Etapes d'implementation
Documenter contrat workflow unique, contrats samplesheets, outputs publics, reprise future via manifeste, EGAPx, EDTA/masked assembly et troubleshooting.
### Tests
`scripts/run-tests.sh`; verification des commandes documentees par inspection et controles statiques.
### Criteres d'acceptation
Un nouvel utilisateur peut lancer un test local et comprendre comment preparer un run production.
### Risques et retour arriere
Faible.

## TITAN-P4-001 - Optimiser ressources HPC
Priorite : P4
Statut : Fait
Risque : Moyen

### Objectif
Definir CPU, memoire, temps, retries et concurrence par type de processus.
### Constat
Les ressources des modules actifs sont pilotees par labels dans `conf/base.config`; les directives `cpus` locales ont ete retirees des modules et l'ancien defaut global `100GB/20 CPU` a ete supprime.
### Fichiers concernes
`conf/base.config`, `modules/*.nf`.
### Etapes d'implementation
Introduire labels (`process_index`, `process_alignment`, `process_transcriptome`, `process_prediction`, `process_merge`, `process_aegis`, `process_low`) et les appliquer progressivement. Les modules actifs sont maintenant labellises; EDTA, EGAPx et Diamond2GO gardent leurs parametres CPU via `withName`.
### Tests
`python3 scripts/validate_profiles.py`; `nextflow config -flat -profile test`; `nextflow config -flat -profile slurm,apptainer`; `nextflow run main.nf -profile test -stub-run -ansi-log false`; `scripts/run-tests.sh`.
### Criteres d'acceptation
Ressources adaptables sans modifier les modules.
### Risques et retour arriere
Risque de sous-dimensionnement; tester sur petits jeux puis production.

## TITAN-P5-001 - Validation scientifique des annotations
Priorite : P5
Statut : Fait
Risque : Eleve

### Objectif
Valider GFF3 final, FASTA/proteines et coherence biologique.
### Constat
Une validation finale structurelle bloque les erreurs critiques apres AEGIS et publie des rapports JSON/TXT sous `${output_dir}/validation`.
### Fichiers concernes
`scripts/validate_final_annotation.py`, `scripts/test_validate_final_annotation.py`, `modules/validate_final_annotation.nf`, `workflows/titan.nf`, documentation.
### Etapes d'implementation
Verifier GFF3, Parent, phases CDS, coordonnees, seqids, duplications, sequences absentes, statistiques d'annotation et integrite FASTA proteique. La validation biologique approfondie reste une revue scientifique separee.
### Tests
`python3 scripts/test_validate_final_annotation.py`; `nextflow run main.nf -profile test -stub-run -ansi-log false`; `scripts/run-tests.sh`.
### Criteres d'acceptation
Rapport de validation final present et bloque les erreurs critiques.
### Risques et retour arriere
Risque de faux positifs; commencer en mode rapport avant mode bloquant.
