# Plan : utiliser le masquage EDTA là où c'est utile

Contexte : parmi les 14 476 gènes codants sans expression, 36,9 % ont un signal de
transposon (contre 18,7 % des gènes exprimés) et 45 % des gènes BRAKER non exprimés en
ont un. Cause : BRAKER3 (donc AUGUSTUS et GeneMark, ensuite donnés bruts à Mikado) recevait
le génome **non masqué**, alors que le flag `--softmasking` était déjà passé (sans effet,
l'assemblage n'a aucune minuscule).

## Ce qui est fait dans cette modification

| Étape | Génome reçu | Pourquoi |
|---|---|---|
| **BRAKER3** (AUGUSTUS + GeneMark-ES) | **soft-masqué EDTA** (`assembly_softmasked.EDTA.fasta`) | `--softmasking` est déjà actif : AUGUSTUS/GeneMark baissent le poids des répétitions sans effacer la séquence. Source du problème. |
| **Mikado** prepare/serialise/pick | **soft-masqué EDTA** | Mikado met les sites d'épissage en majuscules avant de les tester, donc les minuscules sont sans effet. Avec des `N`, les introns au bord d'un TE perdent leur GT-AG canonique. En production, 81 360 des 1 599 787 transcrits préparés (5,1 %) contiennent des `N`. Le FASTA de transcrits est remis en majuscules après `mikado prepare`. |
| Helixer | hard-masqué EDTA (inchangé) | Helixer met tout en majuscules à l'import du FASTA (`geenuff` : `seq.upper()`). Soft-masqué = non masqué pour lui, il prédirait de nouveau des gènes dans les TE. Seul le `N` l'en empêche. |
| AEGIS, validation | hard-masqué EDTA (inchangé) | déjà le cas |
| EGAPx | **inchangé** (génome brut du YAML) | voir « EGAPx » ci-dessous |
| Liftoff | brut | transfère l'ancienne annotation ; des `N` feraient perdre des gènes |
| STAR / HISAT2 / StringTie / PsiCLASS / Salmon | brut | des lectures de TE seraient réalignées ailleurs ; l'évidence de transcription doit refléter le génome entier |
| tRNAscan-SE, Infernal/Rfam | brut | des ncRNA peuvent chevaucher des répétitions |
| extract_final_transcripts, SQANTI3, provenance | brut | extraction de séquences réelles |

Réglages :

* `--mask_genome_for_prediction false` : retour au comportement précédent (BRAKER3 sur génome brut).
* `--edta_precomputed_dir <dossier>` : n'exécute pas EDTA, réutilise
  `assembly_masked.EDTA.fasta`, `edta.TEanno.gff3`, `edta.TElib.fa`. Le dossier
  `data/edta_precomputed/` contient déjà des liens vers les sorties du run de production.

Le soft-masquage (`scripts/soft_mask_genome.py`) abaisse en minuscules, dans l'assemblage
original, les bases que EDTA a remplacées par `N`. Il apparie les séquences par identifiant
(EDTA met chr00 en premier), vérifie les longueurs, et laisse les `N` natifs en majuscules.
Contrôle sur T2T_ref : 106 510 239 bases masquées sur 504 621 863 (21,11 %), séquence
identique à l'original une fois remise en majuscules.

## Reprise avec EDTA déjà fait

* `-resume` **ne peut pas** réutiliser EDTA : le dossier `data/work` du run de production n'existe plus.
* À la place, `--edta_precomputed_dir /…/TITAN/data/edta_precomputed` saute EDTA.
* Le reste est relancé : BRAKER3 change d'entrée, donc tout ce qui est en aval (Mikado,
  AEGIS, fonctionnel...) aussi. Le cache imbriqué d'EGAPx (`.egapx_work`) devrait permettre de
  reprendre EGAPx, car son génome n'a pas changé.

```bash
./launch_TITAN_serveur_colmar.sh -- --edta_precomputed_dir $PWD/data/edta_precomputed
```

À lancer dans un nouveau `--output_dir`/`TITAN_RUN_NAME` pour garder le run actuel comme référence.

## Points à connaître

1. **EGAPx.** Le README d'EGAPx dit que le génome n'a pas besoin d'être masqué : il lance
   WindowMasker et produit son propre soft-mask pour Gnomon. Le paramètre `softmask` existe
   dans `nf/subworkflows/ncbi/main.nf` mais n'est qu'un commentaire dans la v0.5.2 (Gnomon
   reçoit `win_softmask`). Un FASTA en minuscules n'aurait donc aucun effet, et un FASTA
   hard-masqué empêcherait Gnomon de prédire sous les `N`, y compris de vrais gènes.
2. **Mikado** recevait un génome hard-masqué : 5,1 % des transcrits préparés contenaient des `N`
   (ORF coupés pour TransDecoder, sites d'épissage non canoniques). Corrigé. Effet à mesurer :
   les transcrits RNA-seq dans des TE ne sont plus mutilés, ils peuvent donc être mieux soutenus.
   Je n'ai pas pu tester TransDecoder sur du texte en minuscules (image sans exécutable utilisable
   hors Nextflow) ; c'est pourquoi le FASTA est remis en majuscules.
3. **Risque de perte de vrais gènes.** 4 314 gènes non exprimés ont ≥ 50 % de CDS dans un TE, mais
   4 200 gènes exprimés sont dans le même cas (certains gènes réels contiennent des TE).
   Le soft-masquage n'interdit pas les gènes dans les répétitions (les hints RNA-seq et
   protéines peuvent les soutenir), mais il faut le mesurer.
4. La documentation disait « soft-masked » pour Helixer ; c'est le génome hard-masqué. Corrigé.

## Validation proposée (à faire après le run)

1. **Contrôle direct** : proportion de gènes BRAKER (`augustus.hints.gff3`, `genemark.gtf`)
   dont ≥ 50 % de CDS est dans un TE EDTA, avant/après.
2. **Mêmes métriques que le rapport** : nombre de gènes, classes `no_signal`,
   `te_overlap_unsupported`, `te_protein` (`scripts/run_unsupported_genes_analysis.sh`).
3. **Gènes réels conservés** : BUSCO, et part des gènes `expressed` et
   `conserved_or_functional` du run précédent retrouvés (même critère que
   `docs/user/mikado_braker_input_test.md` : pas plus de ~2-3 % de perte de gènes exprimés).
4. **Critère de décision** : adopter si les gènes TE-protéine/`no_signal` baissent nettement sans
   perte notable de gènes exprimés ni de BUSCO ; sinon `--mask_genome_for_prediction false`.

## Pistes ensuite (non faites)

* Bras « EGAPx sur génome hard-masqué » comme expérience contrôlée, pas par défaut.
* Séparer rétrotransposons / transposons à ADN en croisant `edta.TEanno.gff3` avec les gènes.

## Pourquoi Liftoff et les alignements RNA-seq restent sur le génome brut

* **Liftoff** aligne (minimap2) les gènes de l'ancienne annotation sur le nouveau génome, puis
  les projette. Avec des `N` : le gène lifté tombe dans un trou, il est perdu ou tronqué. Même
  en soft-masqué, rien à gagner : minimap2 ignore la casse. De plus Liftoff n'est pas la source
  du bruit TE ; c'est une évidence indépendante (l'ancienne annotation), utile à garder entière.
* **Alignements RNA-seq (STAR, HISAT2) et assemblages (StringTie, PsiCLASS)** : une lecture qui
  vient d'un TE exprimé ou d'un gène contenant un TE (introns, UTR, exons) doit s'aligner à sa vraie
  place. Si on masque, elle s'aligne ailleurs (autre copie non masquée), ou est perdue, et les
  assemblages deviennent faux. STAR et HISAT2 n'utilisent pas le soft-masking. Ce sont aussi
  ces données qui mesurent « exprimé / non exprimé » : masquer ferait perdre l'évidence
  qui sert à juger les gènes TE.
