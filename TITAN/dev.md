# Prompt d'amorçage pour Codex (développement isolé)

À coller tel quel comme premier prompt à Codex avant toute session de développement sur les phases de `to_do_add.md`.

```
Contexte : tu développes dans un git worktree isolé, PAS dans le dépôt de production.

Ton dossier de travail exclusif :
/data2/avelt/2026_annotations_PN40024_T2T/workflows/.claude-worktrees/codex-dev/TITAN
Branche : codex-dev

Règles strictes, à ne jamais enfreindre :

1. Ne jamais lire, écrire, ou référencer un chemin absolu commençant par
   /data2/avelt/2026_annotations_PN40024_T2T/workflows/TITAN/
   (sans le .claude-worktrees/codex-dev/ au milieu) — c'est le dépôt de
   production avec un pipeline Nextflow actif dessus. Reste toujours dans
   ton propre worktree, y compris pour tout `cd`, tout chemin passé en
   argument, tout `-work-dir`, tout `-c`.

2. Pour toute exécution Nextflow (test, validation, stub-run) :
   - Utilise systématiquement `-profile test`.
   - Ne passe JAMAIS `-work-dir` ni `-c` pointant vers un chemin en dehors
     de ton worktree (donc jamais vers .../TITAN/data/work ni vers
     .../TITAN/data/slurm_apptainer.config).
   - Laisse Nextflow utiliser ses répertoires par défaut (./work, ./.nextflow)
     à l'intérieur du worktree.
   - Référence pour les commandes de test : scripts/run-tests.sh dans ce
     worktree.

3. Git :
   - Reste sur la branche codex-dev. Ne fais jamais `git checkout main`,
     `git merge`, ni aucune opération touchant la branche main.
   - Commit et push sur codex-dev uniquement.
   - Le merge vers main sera fait manuellement par moi ensuite, depuis le
     dépôt de production — ce n'est pas ta responsabilité.

4. Feuille de route : suis to_do_add.md dans ce worktree (13 phases, une
   par outil, avec jalons et critères de validation précis pour chacune).
   Implémente phase par phase dans l'ordre recommandé en bas du document,
   ne saute pas les jalons de validation.
```

## Pourquoi cette isolation

`.nextflow/history` et `data/work` (le cache Nextflow) sont tous les deux dans `.gitignore` — un worktree neuf n'en a aucune trace. Le launcher de production (`launch_TITAN_serveur_colmar.sh`) relance toujours avec `-resume` sans ID de session explicite : il cible systématiquement la *dernière* entrée de `.nextflow/history`. Si des `-stub-run` de développement s'exécutaient dans le dépôt de production, ils deviendraient cette "dernière entrée" et un `-resume` ultérieur recalculerait tout le pipeline de production depuis zéro au lieu de reprendre où il en était — déjà arrivé deux fois pendant le développement de la fonctionnalité BUSCO/AGAT/MultiQC, corrigé chaque fois en nettoyant manuellement l'historique.

Tant que Codex reste dans le worktree avec `-profile test` et sans chemins explicites vers le dépôt de production, ce risque est structurellement écarté plutôt que dépendant de la discipline de chacun.
