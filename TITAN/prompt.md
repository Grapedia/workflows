# Prompt de developpement pour TITAN

Tu travailles exclusivement sur TITAN, workflow Nextflow DSL2 d'annotation de genomes situe dans ce depot sous `TITAN/`. Le depot de reference `/home/vmadmin/atcg-rnaseq` peut etre consulte uniquement comme reference methodologique; il ne doit jamais etre modifie, reformate, supprime ou committe.

## Identite et objectif du projet

TITAN signifie `The Intensive Transcript ANnotation pipeline`. C'est un pipeline d'annotation de genomes eucaryotes qui combine des evidences de transcriptome, proteines, annotations transferees et predictions ab initio afin de produire une annotation structurale finale GFF3 et des proteines annotees.

Utilisateurs vises: bioinformaticiens, ingenieurs de plateformes et biologistes computationnels executant des annotations de genomes sur VM Linux ou cluster Slurm.

Entrees principales:

* assemblage cible FASTA;
* assemblage precedent ou de reference FASTA;
* annotation precedente GFF3;
* samplesheet RNA-seq court et long;
* samplesheet de proteines;
* fichier de parametres EGAPx si cette branche est activee;
* parametres d'execution et options biologiques.

Sorties principales:

* evidences intermediaires: annotations Liftoff, predictions BRAKER3, transcriptomes StringTie/PsiCLASS, genome masque EDTA;
* annotation finale Aegis en GFF3;
* FASTA proteiques `main` et `all`;
* annotations fonctionnelles Diamond2GO;
* rapports, traces, versions et provenance.

Exigences: reproductibilite, tracabilite, execution locale sans Slurm pour les tests, execution Slurm pour la production, compatibilite Apptainer, absence de chemins absolus propres a une machine dans le workflow principal.

## Perimetre autorise

Fichiers modifiables dans TITAN uniquement:

* `main.nf`, `nextflow.config`, `nextflow_schema.json`;
* `conf/`, `modules/`, `subworkflows/`, `workflows/`;
* `scripts/`, `bin/`, `assets/`, `docs/`;
* `tests/`, `test-data/`, fixtures legeres;
* `.github/workflows/`;
* documentation et fichiers de configuration.

Ne jamais modifier:

* `/home/vmadmin/atcg-rnaseq`;
* fichiers utilisateur non lies a la tache;
* donnees de production;
* secrets, cles SSH, tokens;
* gros artefacts Nextflow, `work/`, `.nextflow/`, images SIF, resultats volumineux.

Fichiers generes a ne pas versionner: `work/`, `.nextflow/`, `.nextflow.log*`, `OUTDIR/`, `results/`, `test-results/`, `test-work/`, fichiers temporaires, logs volumineux, images Apptainer volumineuses.

## Principes de developpement

Avant toute modification, inspecter le depot, identifier le point d'entree, les parametres, les entrees, les sorties, les modules, les conteneurs, les profils et les tests existants. Etablir une baseline reproductible.

Realiser de petits changements atomiques, testables et reversibles. Limiter chaque passe a une fonctionnalite coherente. Eviter les refactorings non lies. Ne pas ajouter de dependance sans justification. Ne pas modifier une interface publique sans compatibilite ou deprecation documentee. Conserver les anciens parametres lorsque c'est possible.

Ne jamais masquer une erreur, ignorer silencieusement une entree invalide, utiliser `|| true` pour faire passer un test, supprimer un controle qualite ou pretendre qu'un test a reussi sans l'avoir execute.

Toute rupture necessaire doit etre documentee avec justification, impact, solution de repli et test de non-regression.

## Regles Nextflow

* DSL2 obligatoire.
* `main.nf` doit rester un orchestrateur leger.
* Les modules doivent etre autonomes et avoir une responsabilite principale.
* Les sous-workflows groupent des etapes biologiquement coherentes.
* Les parametres doivent etre explicites, documentes et valides autant que possible.
* Les channels doivent etre nommes clairement.
* Les tuples doivent etre documentes dans les modules ou meta-informations.
* Les directives de ressources doivent etre centralisees par labels.
* Les labels de processus doivent etre coherents.
* Les versions des outils doivent etre explicites et collectees.
* Les outputs doivent etre declares, deterministes et previsibles.
* Le workflow doit etre compatible avec `-resume`.
* Aucun chemin absolu specifique a une machine ne doit etre requis.
* Aucune dependance implicite a Slurm ne doit exister pour les tests.
* Le profil `test` doit fonctionner avec l'executeur local.

## Regles pour les modules

Chaque outil doit idealement etre encapsule dans un module distinct. Chaque module doit preciser son role, ses entrees, ses sorties, ses parametres, son conteneur, la version de l'outil, ses ressources, sa commande, un bloc `stub` lorsque pertinent, son mecanisme de version, ses tests et ses erreurs attendues.

Eviter les scripts monolithiques melangeant telechargement, pretraitement, analyse, validation, aggregation et publication. Les scripts Bash doivent utiliser `set -euo pipefail`, citer les chemins et verifier les fichiers attendus.

## Conteneurs Apptainer

Utiliser des images versionnees et immuables. Ne pas utiliser `latest` pour les analyses reproductibles. Preferer BioContainers officiel, image officielle de l'outil, image projet controlee, puis image specifique reproductible.

Les images doivent etre compatibles Apptainer/Singularity. Le cache doit etre configurable. Aucun telechargement non maitrise de dependance ne doit avoir lieu pendant une analyse. Les versions des outils doivent etre collectees dans les resultats.

## Tests

Maintenir au minimum:

* tests unitaires de scripts et modules;
* tests d'integration de sous-workflows;
* test minimal de bout en bout;
* tests d'entrees invalides;
* test de reprise `-resume`;
* test local sans Slurm;
* verification des sorties attendues;
* comparaison par contenu quand pertinent.

Aucun test ne doit utiliser de donnees de production volumineuses.

## Jeux de donnees de test

Les fixtures doivent etre petites, documentees, biologiquement coherentes, publiques ou synthetiques, sans donnees confidentielles. Documenter source, licence, methode de reduction, commandes de generation, checksums et sorties attendues. Ne jamais telecharger automatiquement de grands jeux de donnees.

## Profils d'execution

Maintenir au minimum:

* `standard`: comportement par defaut conservateur;
* `local`: executeur local, ressources limitees, aucun Slurm;
* `test`: donnees minimales, execution rapide, aucun Slurm;
* `apptainer`: activation Apptainer sans Docker daemon;
* `slurm`: directives cluster isolees.

Le profil `test` doit fonctionner sans Slurm. Le profil `slurm` ne doit pas rendre les autres profils dependants du cluster.

## Documentation

Toute nouvelle fonctionnalite doit mettre a jour README, parametres, exemples de commandes, entrees, sorties, dependances, profils, tests, limitations connues et changelog si present.

## Controles avant validation

Commandes minimales a adapter selon la tache:

```bash
nextflow config
nextflow config -profile test
nextflow run main.nf -profile test -stub-run
nextflow run main.nf -profile test -stub-run -resume
```

Si une commande ne peut pas etre executee, documenter l'erreur exacte et ne pas annoncer de succes.

## Git

Ne jamais committer de secrets, gros fichiers generes ou fichiers ignores. Respecter `.gitignore`. Faire des changements atomiques avec messages descriptifs. Ne jamais reecrire l'historique. Ne jamais pousser automatiquement. Ne jamais ecraser des changements utilisateur non lies.

## Criteres de fin

Une tache est terminee seulement si le workflow ou le composant concerne se lance, les tests cibles passent, les tests existants pertinents passent, les sorties attendues sont presentes, les modifications sont documentees, aucune fonctionnalite anterieure n'est cassee et aucun fichier non lie n'a ete modifie.
