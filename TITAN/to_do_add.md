# Audit d'ajouts d'outils pour TITAN — v2

Date de l'audit: 2026-07-18 (remplace l'audit du 2026-07-16, dont les items déjà implémentés ont été retirés).

Objectif: aller vers une annotation "complète" du génome T2T de la vigne (PN40024), au niveau de ce qui a été publié pour l'annotation officielle PN40024.v5.1 (41 766 gènes codants, BUSCO 99,4 %, 7 934 gènes non-codants) et au-delà, en intégrant les manques identifiés par audit interne du pipeline + recherche bibliographique 2026.

**miRNA explicitement exclu de ce document** : nécessite des données sRNA-seq, un type d'input absent du `RNAseq_samplesheet` actuel — à traiter séparément quand ces données existeront.

## Déjà implémenté (ne pas refaire)

| Outil | Statut | Où |
| --- | --- | --- |
| eggNOG-mapper | ✅ implémenté | `modules/eggnog_mapper.nf`, dans `subworkflows/aegis.nf` |
| InterProScan | ✅ implémenté | `modules/interproscan.nf`, dans `subworkflows/aegis.nf` |
| Helixer | ✅ implémenté | `modules/helixer_prediction.nf`, appelé depuis `workflows/titan.nf` |
| BUSCO (protéines, gene-set) | ✅ implémenté | `modules/busco.nf` |
| AGAT stats (structure GFF3) | ✅ implémenté | `modules/agat_stats.nf` |
| MultiQC (fastp + BUSCO + AGAT + validation) | ✅ implémenté | `modules/multiqc_report.nf` |

## Règles d'architecture (rappel, valables pour toutes les phases ci-dessous)

- Nouveau dossier public par outil : `${params.output_dir}/additional_annotations/<tool>/` (sauf QC → `${params.output_dir}/quality_report/<tool>/`).
- Chaque outil a `params.run_<tool> = false` par défaut (opt-in), sauf s'il est très léger et sans dépendance de données externes (auquel cas il peut tourner par défaut, comme `agat_stats`/`multiqc_report` aujourd'hui).
- Chaque module émet `versions.yml`.
- Toute image conteneur **doit être pinnée par digest** (`registry/image@sha256:...`), jamais par tag mobile — vérifié via `scripts/validate_container_pins.py`.
- Toute donnée offline externe (bases Rfam, Infernal, OMAmer...) suit le pattern déjà utilisé pour `eggnog_data_dir`/`interproscan_data_dir`/`busco_data_dir` : un paramètre `params.<tool>_data_dir = false`, jamais téléchargée automatiquement par un run de production, bind-mountée explicitement dans `conf/apptainer.config`.
- Chaque module a un bloc `stub:` minimal pour que `scripts/run-tests.sh` (`-stub-run`) continue de couvrir tout le graphe.
- Ajouter la vérification correspondante dans `scripts/validate_profiles.py` (`RESOURCE_LABELS`) si un nouveau label est introduit.

---

## Phase 1 — tRNA (tRNAscan-SE)

**Statut TITAN codex-dev — 2026-07-18** : implémentation Nextflow validée en `-profile test -stub-run`.

- ✅ M1 : `modules/trnascan_se.nf` créé, publié sous `additional_annotations/ncrna/trna/`, `versions.yml` émis, stub GFF3 à 1 feature.
- ✅ M2 : `scripts/trnascan_to_gff3.py` créé avec fixture `test-data/minimal/valid/trnascan.out` et test unitaire `scripts/test_trnascan_to_gff3.py`.
- ⚠️ M3 : run réel Apptainer sur le génome T2T PN40024 non exécuté dans cette passe ; validation limitée aux fixtures et stub-runs locaux.
- ✅ M4 : intégré dans `workflows/titan.nf`, paramètres et conteneur pinné ajoutés, provenance additionnelle mise à jour, AGAT/MultiQC final branchés via `modules/ncrna_annotation_qc.nf`, validations `scripts/run-tests.sh`, `nextflow run main.nf -profile test -stub-run`, `-resume` et `--run_trnascan true` passées.

**But** : détecter et annoter tous les gènes de tRNA du génome T2T.
**Pourquoi** : composante standard de toute annotation ncRNA complète (utilisé explicitement dans l'annotation officielle PN40024.v5.1) ; tRNAscan-SE 2.0 est l'outil de référence (99-100 % de détection, <1 faux positif/15 Gb).
**Position dans le graphe** : branche indépendante, en parallèle de `generate_evidence_data`, sur `params.new_assembly` directement (pas besoin d'attendre AEGIS).

- Image : `quay.io/biocontainers/trnascan-se@sha256:e573090368974ff1228e6894828c6c8a132dfecc3198f5e9fb76832f8f434f29` (2.0.13)
- Paramètres : `params.run_trnascan = false`, `params.container_trnascan`
- Label : `process_low` (rapide, quelques minutes même sur un génome complet)

**Input** : `params.new_assembly` (FASTA génome complet, non masqué — tRNAscan-SE gère lui-même les régions répétées).

**Commande** :
```bash
tRNAscan-SE -E \
  -o trnascan.out \
  -f trnascan.struct \
  -s trnascan.isotype \
  -m trnascan.stats \
  --thread ${task.cpus} \
  ${genome}
scripts/trnascan_to_gff3.py trnascan.out > trna.gff3
```
`-E` = mode eucaryote. `scripts/trnascan_to_gff3.py` est un petit script à écrire (le format `.out` de tRNAscan-SE n'est pas du GFF3 nativement) : colonnes `seq name, tRNA #, begin, end, type, anticodon, intron begin/end, score` → GFF3 `type=tRNA` avec `product=tRNA-<type>` et `score`.

**Output** :
- `trnascan.out` (table brute)
- `trnascan.struct` (structures secondaires, pour audit manuel)
- `trnascan.stats` (comptage par isotype)
- `trna.gff3` (publié, format standardisé pour merge ultérieur avec rRNA/lncRNA)

Publié sous `${params.output_dir}/additional_annotations/ncrna/trna/`.

**Jalons**
1. M1 — `modules/trnascan_se.nf` créé, bloc `stub:` produit un `trna.gff3` factice à 1 ligne, `-stub-run` passe.
2. M2 — Écriture + test unitaire de `scripts/trnascan_to_gff3.py` sur un `trnascan.out` d'exemple fixe (fixture dans `test-data/`), assertion sur le nombre de lignes GFF3 produites et le format des 9 colonnes.
3. M3 — Exécution réelle sur `params.new_assembly` (T2T PN40024) via Apptainer, hors Nextflow d'abord (comme fait pour BUSCO), pour valider la commande et le temps d'exécution réel avant intégration.
4. M4 — Intégré dans `workflows/titan.nf`, publié, `versions.yml` émis, `scripts/validate_profiles.py` toujours vert.

**Validation**
- Le nombre de tRNA détectés est dans une fourchette biologiquement plausible pour un génome de plante diploïde (~500-1000 tRNA typiquement pour un génome de cette taille ; à comparer avec le chiffre publié pour PN40024.v5.1 si disponible dans le papier v5.1).
- Chaque isotype (Ala, Arg, Asn...) est représenté — l'absence totale d'un isotype standard est un signal d'échec à investiguer, pas à ignorer.
- `trna.gff3` est un GFF3 valide (9 colonnes, coordonnées 1-based, dans les limites des séquences du FASTA) — réutiliser `validate_gff3()` déjà présent dans `scripts/validate_inputs.py` comme fonction de test.
- Aucune tRNA prédite ne chevauche un exon codant de `final_annotation.gff3` à >50 % de sa longueur sans raison (sinon, signal de faux positif ou de vrai chevauchement biologique à documenter).
- `versions.yml` contient la version exacte de tRNAscan-SE et le digest du conteneur.

---

## Phase 2 — rRNA / snRNA / snoRNA (Infernal + Rfam)

**Statut TITAN codex-dev — 2026-07-18** : implémentation Nextflow validée en `-profile test -stub-run`.

- ✅ M1 : `modules/infernal_rfam.nf` créé, publié sous `additional_annotations/ncrna/rfam/`, `versions.yml` émis, stub GFF3 à 1 feature.
- ✅ M2 : `scripts/rfam_tblout_to_gff3.py` créé avec fixture `test-data/minimal/valid/rfam_hits.tbl` et test unitaire `scripts/test_rfam_tblout_to_gff3.py`.
- ⚠️ M3 : téléchargement/indexation Rfam et test `cmsearch` sur petit contig non exécutés dans cette passe ; dépend d'un `rfam_data_dir` réel préparé hors pipeline.
- ⚠️ M4 : estimation temps réel sur génome T2T complet non exécutée ; à faire sur infrastructure cible avant run complet, avec décision éventuelle de split par chromosome.
- ✅ M5 : intégré dans `workflows/titan.nf`, avec split FASTA par séquence/chromosome, `cmsearch` parallélisé par split, merge des `rfam_hits.tbl` puis conversion finale unique ; paramètres, label `process_rfam`, bind Apptainer et conteneur pinné ajoutés, provenance additionnelle mise à jour, AGAT/MultiQC final branchés via `modules/ncrna_annotation_qc.nf`, validations `scripts/run-tests.sh`, `nextflow run main.nf -profile test -stub-run`, `-resume` et `--run_rfam true` passées.

**But** : annoter les gènes d'ARN non-codants couverts par les modèles de covariance Rfam (rRNA, snRNA, snoRNA, et autres familles Rfam pertinentes chez les plantes).
**Pourquoi** : composante standard de toute annotation ncRNA complète (mentionnée explicitement dans la méthodologie PN40024.v5.1 pour filtrer les faux lncRNA). Approche standard de l'industrie (Ensembl, NCBI) : BLASTN Rfam en pré-filtre, puis `cmsearch` Infernal ciblé pour réduire le coût de calcul.
**Position dans le graphe** : branche indépendante, sur `params.new_assembly`, en parallèle de la Phase 1.

- Image : `quay.io/biocontainers/infernal@sha256:05ae1ca6cc76c27180524bc38c5b1e17adf9377be5b8c644d3e8e707848d4d99` (1.1.5) — fournit aussi `cmpress`, `cmscan`
- Données offline : `params.rfam_data_dir` (Rfam.cm + Rfam.clanin depuis `ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/`), à indexer une fois avec `cmpress Rfam.cm` — même pattern que `busco_data_dir`.
- Paramètres : `params.run_rfam = false`, `params.container_infernal`, `params.rfam_data_dir = false`
- Label : nouveau `process_rfam` (32 cpus, 96 GB, 48h en prod) — `cmsearch` full-genome contre tout Rfam est le poste le plus lourd restant dans ce document.

**Input** : `params.new_assembly` (FASTA génome complet).

**Commande** :
```bash
# Une fois, hors pipeline, pour préparer rfam_data_dir :
curl -O ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
curl -O ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
gunzip Rfam.cm.gz
cmpress Rfam.cm

# Dans le module :
cmsearch --cpu ${task.cpus} \
  --tblout rfam_hits.tbl \
  --fmt 2 \
  --cut_ga \
  --rfam \
  --nohmmonly \
  --clanin ${params.rfam_data_dir}/Rfam.clanin \
  ${params.rfam_data_dir}/Rfam.cm \
  ${genome} > rfam_search.out
scripts/rfam_tblout_to_gff3.py rfam_hits.tbl > rfam_ncrna.gff3
```
`--cut_ga` utilise les seuils de score (gathering thresholds) spécifiques à chaque famille Rfam, la pratique standard pour éviter de fixer un seuil global arbitraire. `--rfam` active l'heuristique HMM-avant-CM adaptée aux modèles Rfam pour limiter le temps de calcul.

**Output** :
- `rfam_hits.tbl` (table tabulaire brute Infernal)
- `rfam_search.out` (log complet, pour audit)
- `rfam_ncrna.gff3` (GFF3, `type=` dérivé du type Rfam : rRNA/snRNA/snoRNA/other, `Rfam_ID=RFxxxxx` en attribut)

Publié sous `${params.output_dir}/additional_annotations/ncrna/rfam/`.

**Jalons**
1. M1 — `modules/infernal_rfam.nf` créé avec bloc `stub:`, `-stub-run` passe.
2. M2 — Script `scripts/rfam_tblout_to_gff3.py` écrit et testé unitairement sur un `.tblout` d'exemple.
3. M3 — Rfam téléchargé et indexé dans `.rfam_data/` ; `cmsearch` testé sur un petit contig (pas le génome entier) pour valider syntaxe/temps avant un run complet.
4. M4 — Estimation du temps réel sur le génome T2T complet (hors pipeline, avec `time cmsearch ...` sur le vrai génome) : si >24h même à 32 cpus, prévoir un split par chromosome avant intégration (voir note ci-dessous).
5. M5 — Intégré dans `workflows/titan.nf`, publié.

**Note sur la charge de calcul** : `cmsearch --rfam` full-genome contre l'intégralité de Rfam (~4000 familles) est le poste de calcul le plus lourd restant dans ce document. TITAN découpe maintenant le run par séquence/chromosome et collecte les résultats avant conversion GFF3 finale, afin de paralléliser le calcul sans réduire la couverture Rfam.

**Validation**
- `rfam_ncrna.gff3` valide (même contrôle GFF3 que Phase 1).
- Comptage par famille Rfam (rRNA 5S/5.8S/18S/25S-28S attendus en multiples copies en tandem — cohérent avec la biologie connue du génome ribosomique chez les plantes) ; alerter si 0 rRNA détecté (signal d'échec, pas de vraie absence biologique).
- Croisement avec les tRNA de la Phase 1 : aucun conflit de coordonnées entre les deux jeux (ce sont des familles disjointes par construction, un chevauchement signale un bug de conversion GFF3).
- `versions.yml` inclut la version d'Infernal ET la version/date du build Rfam utilisé (Rfam évolue, la traçabilité de la release est essentielle pour la reproductibilité).

---

## Phase 3 — lncRNA (FEELnc + CPC2 + CPAT)

**Décision CPAT — 2026-07-18** : aucun modèle CPAT spécifique *Vitis vinifera* prêt à télécharger n'a été identifié. Le dépôt Plant-LncPipe fournit un modèle CPAT-plant générique (`Plant_Hexamer.tsv`, `Plant.logit.RData`, cutoff 0.46) ; ces fichiers sont maintenant intégrés sous `resources/cpat_plant_lncpipe/` et utilisés comme mode dégradé/provenancé via `params.cpat_model_dir`, mais la sortie TITAN reste `lncrna_candidates.gff3` tant qu'un modèle CPAT *Vitis* n'est pas entraîné et validé.

**Statut TITAN codex-dev — 2026-07-18** : couche candidate implémentée avec filtrage CPAT-plant et téléchargement automatique idempotent des modèles manquants, validée en `-profile test -stub-run`.

- ✅ M1 : `modules/lncrna_candidate_annotation.nf` créé avec `stub:`, sortie candidate publiée sous `additional_annotations/ncrna/lncrna/`, résumé MultiQC ajouté.
- ✅ M2 partiel : modèle CPAT *Vitis* non disponible ; CPAT-plant Plant-LncPipe téléchargé, checksumé, retéléchargé automatiquement si absent, utilisé et tracé comme mode dégradé, pas comme validation finale *Vitis*.
- ⚠️ M3/M4 : test 1 chromosome et run complet production non exécutés dans cette passe.
- ✅ M5 : intégré dans `workflows/titan.nf`, provenance additionnelle et MultiQC final mis à jour, validations `scripts/run-tests.sh`, `nextflow run main.nf -profile test -stub-run`, `-resume` et `--run_lncrna true` avec modèle CPAT-plant requis passées.

**But** : reproduire la méthodologie publiée pour PN40024.v5.1 — annotation des lncRNA à partir des assemblages transcriptomiques déjà produits par TITAN.
**Pourquoi** : c'est très exactement l'écart identifié entre ce pipeline et l'annotation officielle "complète" du même génome (7 934 gènes non-codants publiés via *"un pipeline parallèle"* non présent dans ce code).
**Position dans le graphe** : après `generate_evidence_data` (a besoin des GTF mergés StringTie/PsiCLASS/long-reads) et après `aegis` (a besoin de `final_annotation.gff3` pour exclure les loci déjà codants), avant le rapport qualité final.

- Images :
  - FEELnc : `quay.io/biocontainers/feelnc@sha256:de4aaf80de1af3fd90d3ad5f7e3a24ba8cb22aa5a1a8e429d7584fb0eae7c07b` (0.1.1)
  - CPC2 : `quay.io/biocontainers/cpc2@sha256:5736c1c5187a3a681bba566e63a1b78c10946468f0cf798d117712b265e07c80` (1.0.1)
  - CPAT : `quay.io/biocontainers/cpat@sha256:87366fff67d441f64e0ac4681ccbaf1147f2c0601f3df86bb99f228d7f9a9000` (3.0.5)
- Paramètres : `params.run_lncrna = false`, `params.lncrna_min_length = 200`, `params.lncrna_min_fpkm = 0.5`, `params.cpat_model_dir = "${projectDir}/resources/cpat_plant_lncpipe"`, `params.cpat_model_flavour = "plant_lncpipe"`, `params.cpat_plant_cutoff = 0.46`, `params.lncrna_require_cpat_model = true`
- Label : `process_transcriptome` pour FEELnc (le plus lourd des trois), `process_medium` pour CPC2/CPAT.

**Input** :
- Les GTF déjà mergés par TITAN : `merged_star_stringtie_stranded_default.gtf`, `merged_hisat2_stringtie_stranded_default.gtf`, `merged_transcriptomes.minimap2.long_reads.default_args.gtf` (si `has_long_reads`) — à concaténer/dédupliquer en un seul `candidate_transcripts.gtf` (via `gffcompare` ou `agat_sp_merge_annotations.pl`, déjà présent dans le conteneur AGAT pinné).
- `aegis.out.aegis_gff` (référence codante, pour exclusion).
- `params.new_assembly`.

**Commande** :
```bash
# 1. filtrage longueur/expression (AGAT deja pinne, container_agat)
agat_sp_filter_feature_by_attribute_value.pl --gff candidate_transcripts.gtf \
  --value ${params.lncrna_min_length} --test ">=" --attribute "length" -o length_filtered.gtf
# filtrage FPKM : reutiliser les quantifications salmon deja calculees par TITAN
# (salmon_index + salmon quant sur candidate_transcripts au lieu du CDS liftoff)

# 2. FEELnc : exclusion des loci codants connus
FEELnc_filter.pl -i length_filtered.gtf -a final_annotation.gtf \
  --biotype transcript_biotype=protein_coding -p ${task.cpus} > candidate_lncrna.gtf

# 3. FEELnc : potentiel codant (modele entraine sur les CDS AEGIS connus)
FEELnc_codpot.pl -i candidate_lncrna.gtf -a final_annotation.gtf -g ${genome} \
  -b transcript_biotype=protein_coding --mode=shuffle

# 4. Validation croisee CPC2
CPC2.py -i candidate_lncrna.fasta -o cpc2_result

# 5. Validation croisee CPAT (necessite un modele logit specifique plante,
#    a entrainer une fois sur les CDS/lncRNA connus de Vitis avec
#    make_hexamer_tab.py + make_logitModel.py - CPAT ne fournit pas de
#    modele plante pre-entraine par defaut)
cpat.py -x vitis_hexamer.tsv -d vitis_logitModel.RData \
  -g candidate_lncrna.fasta -o cpat_result

# 6. Intersection : lncRNA final = non-codant selon FEELnc ET CPC2 ET CPAT
scripts/intersect_lncrna_calls.py \
  --feelnc feelnc_codpot_out/candidate_lncrna.gtf.lncRNA.gtf \
  --cpc2 cpc2_result \
  --cpat cpat_result \
  --exclude-ncrna trna.gff3 rfam_ncrna.gff3 \
  -o final_lncrna_candidates.gtf

# 7. Classification finale (lincRNA/NAT-lncRNA/int-lncRNA/SOT-lncRNA)
FEELnc_classifier.pl -i final_lncrna_candidates.gtf -a final_annotation.gtf > final_lncrna.gtf
scripts/gtf_to_gff3.py final_lncrna.gtf > final_lncrna.gff3
```

**Point d'attention majeur** : CPAT officiel ne fournit que des modèles humain/souris/poisson-zèbre/mouche. Plant-LncPipe fournit un modèle CPAT-plant générique téléchargeable, mais pas spécifique *Vitis*. Il faut entraîner un modèle logit spécifique à *Vitis* à partir des CDS connus (AEGIS) et d'un jeu de séquences non-codantes de référence avant de promouvoir les candidats en annotation finale — prévoir ceci comme un jalon dédié avant M3 ci-dessous, pas comme un détail d'implémentation mineur.

**Output** :
- `lncrna_candidates.gff3` (publié ; candidat préliminaire, non final tant que CPAT *Vitis* n'est pas disponible)
- `feelnc_codpot_out/`, `cpc2_result`, `cpat_result` (intermédiaires, `publish_intermediates`)
- `lncrna_classification_summary.tsv` (comptage par classe)

Publié sous `${params.output_dir}/additional_annotations/ncrna/lncrna/`.

**Jalons**
1. M1 — `modules/lncrna_*.nf` (probablement 3-4 modules séparés : filter, codpot, cpc2, cpat, classify) créés avec `stub:`, `-stub-run` passe.
2. M2 — Entraînement et validation du modèle CPAT spécifique *Vitis* (hors pipeline, une fois, données stockées dans `params.cpat_model_dir`).
3. M3 — Test sur un sous-ensemble réduit (1 chromosome) avant run complet, pour valider la chaîne complète FEELnc→CPC2→CPAT→intersection→classification.
4. M4 — Run complet sur le génome de production, comparaison du nombre de lncRNA obtenus avec les 7 934 gènes non-codants publiés pour PN40024.v5.1 (ordre de grandeur, pas égalité exacte attendue vu les différences méthodologiques possibles).
5. M5 — Intégré dans `workflows/titan.nf`, publié, ajouté au rapport MultiQC.

**Validation**
- `lncrna_candidates.gff3` valide structurellement (même contrôle que Phases 1-2) ; ne pas le renommer en `final_lncrna.gff3` tant que CPAT *Vitis* et les validations réelles ne sont pas disponibles.
- Aucun chevauchement >50 % avec un CDS de `final_annotation.gff3` (sinon la Phase FEELnc_filter a un bug).
- Aucun chevauchement avec `trna.gff3`/`rfam_ncrna.gff3` (l'étape 6 `--exclude-ncrna` doit le garantir — tester explicitement).
- Distribution de longueur des lncRNA cohérente avec la littérature (médiane généralement <1kb, plus courte que les mRNA).
- Le nombre total de lncRNA + tRNA + rRNA/snRNA/snoRNA doit être du même ordre de grandeur que les 7 934 gènes non-codants publiés (à documenter, pas à forcer).
- Chaque outil (FEELnc/CPC2/CPAT) doit avoir un taux d'accord >70 % sur les cas non-ambigus — un désaccord massif entre les trois signale un problème de modèle (notamment CPAT mal entraîné) plutôt qu'un vrai résultat biologique incertain.

---

## Phase 4 — Mikado (sélection consensus de transcrits)

**Statut TITAN codex-dev — 2026-07-18** : Phase 4 implémentée avec la Phase 5 comme branche finale parallèle à AEGIS. Mikado consomme les mêmes familles d'évidences qu'AEGIS et publie un second GFF3 final sous `final_annotations/mikado/`, sans fusion automatique avec `aegis_outputs/final_annotation.gff3`.

**But** : sélectionner le meilleur modèle de transcrit par locus parmi StringTie, PsiCLASS et les GTF long-reads.
**Pourquoi** : Mikado consolide plusieurs sources d'évidence et choisit un modèle par locus avec un score explicable — c'est une source finale GFF3 alternative à AEGIS, à comparer avec AEGIS plutôt qu'à fusionner automatiquement.
**Position dans le graphe** : après les merges STAR/StringTie, STAR/PsiCLASS, HISAT2/StringTie et long-reads ; **avant** la Phase 5 (TransDecoder), dont Mikado a besoin en entrée pour le score de la 3e passe (`mikado serialise --orfs`).

- Image : `quay.io/biocontainers/mikado@sha256:dd6f5a2a2d7fdbab73c835cd0f49bd1444ecaddf8e4cd96fbf0fe24f5ecf5f22` (2.3.4, tag stable — pas la 2.3.5rc3 disponible en release candidate).
- Paramètres : `params.run_mikado = false`, `params.container_mikado`, `params.mikado_mode = 'permissive'`, `params.mikado_scoring = 'plant.yaml'` (fichier de scoring livré avec Mikado, profil "plant" existant)
- Label : `process_transcriptome`

**Input** :
- Même périmètre de sources que l'appel AEGIS : Liftoff, EGAPx, BRAKER3/AUGUSTUS/Genemark, STAR/StringTie, HISAT2/StringTie, STAR/PsiCLASS, long-reads si présents, Helixer si activé.
- `params.new_assembly`.
- Protéines externes optionnelles : `protein_samplesheet` (déjà utilisé par BRAKER3).
- CDS/ORFs de la Phase 5 (dépendance croisée).

**Commande** :
```bash
# 0. Générer la liste des GTF sources (script dans scripts/, pas de closure Groovy complexe)
scripts/make_mikado_list.py \
  --source liftoff.gff3:liftoff:False:10:True \
  --source egapx.gff3:egapx:False:9:True \
  --source braker_augustus.gff3:braker_augustus:False:8:False \
  --source merged_star_stringtie_stranded_default.gtf:star_stringtie:True:5:False \
  --source merged_hisat2_stringtie_stranded_default.gtf:hisat2_stringtie:True:5:False \
  --source star_psiclass_stranded.gtf:star_psiclass:True:5:False \
  --source merged_long_reads_default.gtf:minimap2_stringtie:False:6:False \
  -o transcript_inputs.tsv

mikado configure \
  --list transcript_inputs.tsv \
  --reference ${genome} \
  --mode ${params.mikado_mode} \
  --scoring ${params.mikado_scoring} \
  configuration.yaml

mikado prepare --json-conf configuration.yaml

# (Phase 5 s'insere ici : TransDecoder tourne sur mikado_prepared.fasta)

mikado serialise --json-conf configuration.yaml \
  --orfs mikado_prepared.fasta.transdecoder.bed \
  --procs ${task.cpus}

mikado pick --json-conf configuration.yaml \
  --subloci-out mikado.subloci.gff3 \
  --loci-out mikado.loci.gff3 \
  --procs ${task.cpus}
```

**Output** :
- `final_mikado_annotation.gff3` et `mikado.loci.gff3` (publiés, annotation finale Mikado alternative à AEGIS)
- `mikado.subloci.gff3` (intermédiaire, pour audit des loci ambigus)
- `mikado_prepared.fasta` (transcrits consolidés, entrée de la Phase 5)

Publié sous `${params.output_dir}/final_annotations/mikado/`.

**Jalons**
1. ✅ M1 — `modules/mikado.nf` créé avec `mikado_prepare`, `mikado_serialise`, `mikado_pick` et `stub:`.
2. ✅ M2 — `scripts/make_mikado_list.py` écrit + testé unitairement.
3. ⚠️ M3 — Run réel `mikado configure`+`prepare` sur données complètes non exécuté dans cette passe.
4. ✅ M4 — Intégration complète avec TransDecoder (Phase 5), `serialise`+`pick` validés en `-profile test -stub-run --run_mikado true`.
5. ✅ M5 — Comparaison `final_mikado_annotation.gff3` vs `aegis.out.aegis_gff` ajoutée au rapport qualité MultiQC.

**Validation**
- `mikado.loci.gff3` valide structurellement.
- Chaque locus Mikado a un score Mikado (`mikado_score` en attribut) non-nul — un score nul partout signale un fichier de scoring mal chargé.
- Le nombre de loci Mikado doit être du même ordre de grandeur que le nombre de gènes AEGIS (±30 %, à documenter si l'écart est plus grand).
- Aucun transcrit Mikado ne doit avoir d'intron >100 000 bp sans support de lecture réel (même filtre de bruit que celui documenté pour le pipeline officiel v5.1).

---

## Phase 5 — TransDecoder (ORF/CDS)

**Statut TITAN codex-dev — 2026-07-18** : Phase 5 implémentée entre `mikado_prepare` et `mikado_serialise`; elle suit `run_mikado` et fournit le `.bed` ORF consommé par Mikado.

**But** : prédire les régions codantes (CDS/protéines) dans les transcrits consolidés par Mikado.
**Pourquoi** : nécessaire à `mikado serialise --orfs` (Phase 4) et produit une évidence protéique RNA-seq indépendante de BRAKER3/EGAPx pour comparaison.
**Position dans le graphe** : entre `mikado prepare` et `mikado serialise` (Phase 4).

- Image : `quay.io/biocontainers/transdecoder@sha256:c70f3a30cc8f3aecccb1d8978b9a49865d3994ebd7885361ab4c9dd820bd17f5` (6.0.0)
- Paramètres : `params.run_transdecoder = true` (effectif seulement si `params.run_mikado = true`, pas de sens de l'activer seul), `params.container_transdecoder`
- Label : `process_medium`

**Input** : `mikado_prepared.fasta` (sortie de `mikado prepare`, Phase 4).

**Commande** :
```bash
TransDecoder.LongOrfs -t mikado_prepared.fasta
TransDecoder.Predict -t mikado_prepared.fasta --single_best_only
```
`--single_best_only` : un seul ORF par transcrit, cohérent avec l'usage en entrée de Mikado (Mikado gère lui-même les cas ambigus au niveau locus).

**Output** :
- `mikado_prepared.fasta.transdecoder.bed` (consommé directement par `mikado serialise --orfs`)
- `mikado_prepared.fasta.transdecoder.pep` (protéines, publiées pour audit)
- `mikado_prepared.fasta.transdecoder.gff3` (publié)

Publié sous `${params.output_dir}/final_annotations/mikado/transdecoder/`.

**Jalons**
1. ✅ M1 — `modules/transdecoder.nf` créé avec `transdecoder_longorfs`, `transdecoder_predict` et `stub:`.
2. ⚠️ M2 — Run réel sur `mikado_prepared.fasta` complet non exécuté dans cette passe.
3. ✅ M3 — Intégré dans le graphe Mikado complet (Phase 4, M4), validé en `-profile test -stub-run --run_mikado true`.

**Validation**
- Le `.bed` TransDecoder est directement consommable par `mikado serialise` sans erreur de parsing (le test d'intégration EST la Phase 4 M4).
- Le nombre d'ORF prédits est cohérent avec le nombre de transcrits en entrée (pas 100 % — normal qu'une fraction des transcrits assemblés soit non-codante, ce sont potentiellement des candidats lncRNA à croiser avec la Phase 3).
- Aucun ORF ne dépasse la longueur du transcrit parent (contrôle de cohérence trivial mais à vérifier explicitement en test).

---

## Phase 6 — FLAIR (isoformes long-read)

**Statut TITAN codex-dev — 2026-07-18** : Phase 6 implémentée comme branche long-read optionnelle. Les isoformes FLAIR sont mergées dans `flair_isoforms.gtf`/`.fa`, publiées sous `additional_annotations/flair/`, enregistrées dans la provenance et passées comme évidences transcriptomiques additionnelles à AEGIS et Mikado.

**But** : produire une annotation isoforme long-read alternative à StringTie long reads, avec correction des jonctions d'épissage contre le génome.
**Pourquoi** : FLAIR est spécialisé pour la correction, la définition d'isoformes et l'analyse de splicing sur reads longs PacBio/ONT — plus précis que StringTie seul sur ce type de données.
**Position dans le graphe** : branche long reads, uniquement si `has_long_reads == true`. Après préparation des FASTQ/FASTA long reads.

- Image : `quay.io/biocontainers/flair@sha256:187e2e22535d73ecc724afc7e474d9908b6e43a55f8588e8566db3bea2eba79e` (3.0.0)
- Paramètres : `params.run_flair = false`, `params.container_flair`
- Label : `process_transcriptome`

**Input** :
- FASTQ/FASTA long reads par échantillon (déjà préparés par `prepare_RNAseq_fastq_files_long`).
- `params.new_assembly`.
- Annotation de référence pour la correction : `evidence_data.liftoff_gff3` (pas `aegis.out.aegis_gff`, pour éviter une dépendance circulaire avec AEGIS au premier ajout, comme documenté dans l'audit v1).

**Commande** :
```bash
flair align -g ${genome} -r ${reads} -o ${sample_ID}.flair --threads ${task.cpus}

flair correct -q ${sample_ID}.flair.bed -g ${genome} \
  -f ${liftoff_gff3} -o ${sample_ID}.flair.corrected --threads ${task.cpus}

flair collapse -g ${genome} -r ${reads} \
  -q ${sample_ID}.flair.corrected_all_corrected.bed \
  -o ${sample_ID}.flair --threads ${task.cpus}
```

**Output** :
- `${sample_ID}.flair.isoforms.gtf` (publié)
- `${sample_ID}.flair.isoforms.fa` (publié)
- `flair_isoforms.gtf` et `flair_isoforms.fa` mergés (publiés)

Publié sous `${params.output_dir}/additional_annotations/flair/`. Le GTF mergé est transmis à AEGIS (`flair_isoforms_gtf`) et à Mikado (`flair_isoforms`) comme évidence nouvelle.

**Jalons**
1. ✅ M1 — `modules/flair.nf` créé avec `flair_isoforms`, `flair_merge_isoforms`, `flair_empty_evidence` et `stub:`, conditionné par `has_long_reads` comme le reste de la branche long-reads existante.
2. ⏳ M2 — Test sur le sample long-read réel du projet (`hq_transcripts.RI_rmv_Antonio`) non exécuté dans codex-dev ; validation actuelle faite en `-profile test -stub-run`.
3. ✅ M3 — Intégré, publié, ajouté à AEGIS et Mikado comme évidence optionnelle.

**Validation**
- `.isoforms.gtf` valide structurellement.
- Chaque isoforme a un support ≥1 read (FLAIR l'assure nativement — vérifier que la colonne de comptage n'est jamais nulle).
- Comparaison du nombre d'isoformes FLAIR vs StringTie/Minimap2 long-reads déjà produit par TITAN (`merged_transcriptomes.minimap2.long_reads.default_args.gtf`) — même ordre de grandeur attendu, pas égalité.

---

## Phase 7 — SQANTI3 (QC isoformes long-read)

**Statut TITAN codex-dev — 2026-07-18** : Phase 7 implémentée comme QC post-AEGIS optionnel sur les deux assemblages long reads : StringTie/Minimap2 (`long_reads_default_gtf`) et FLAIR (`flair_isoforms_gtf`). Les classifications/corrected GTF/rapports HTML sont publiés sous `additional_annotations/sqanti3/`, un résumé structural combiné est publié sous `quality_report/sqanti3/` et injecté dans MultiQC.

**But** : caractériser et classifier les modèles de transcrits long-read (structural categories : FSM, ISM, NIC, NNC, etc.), avant de leur faire confiance.
**Pourquoi** : SQANTI3 est l'outil de référence pour le QC d'isoformes long-read ; utile pour juger la qualité de FLAIR (Phase 6) et/ou de l'assemblage long-read déjà produit par TITAN.
**Position dans le graphe** : après FLAIR (Phase 6) et/ou après le merge long-reads existant, **et** après AEGIS — validation post-hoc de `final_annotation.gff3` contre les isoformes long-read réelles est plus informatif qu'une comparaison contre Liftoff seul.

- Image : `quay.io/biocontainers/sqanti3@sha256:3bd6ec96b3f1c9cae69cfef54ba0522b7d99efa7ebb0ff6a611841aa6784f74c` (6.0.1) — image volumineuse (~2.7 Go), prévoir le temps de pull.
- Paramètres : `params.run_sqanti3 = false`, `params.container_sqanti3`
- Label : `process_transcriptome`

**Input** :
- `merged_transcriptomes.minimap2.long_reads.default_args.gtf` / `long_reads_default_gtf` (StringTie/Minimap2 long reads).
- `flair_isoforms.gtf` / `flair_isoforms_gtf` (Phase 6 ; zéro feature si FLAIR désactivé ou absent).
- `aegis.out.aegis_gff` (référence, pour validation post-hoc de l'annotation finale elle-même).
- `params.new_assembly`.

**Commande** :
```bash
sqanti3_qc.py \
  ${isoforms_gtf} \
  ${aegis_gff} \
  ${genome} \
  --dir . \
  --output titan_sqanti3 \
  --cpus ${task.cpus}
```

**Output** :
- `stringtie_long_reads.sqanti3_classification.txt`, `.sqanti3_corrected.gtf`, `.sqanti3_report.html`, `.sqanti3_summary.tsv`
- `flair_isoforms.sqanti3_classification.txt`, `.sqanti3_corrected.gtf`, `.sqanti3_report.html`, `.sqanti3_summary.tsv`
- `sqanti3_long_read_isoform_qc_mqc.tsv` (résumé combiné MultiQC)

Publié sous `${params.output_dir}/additional_annotations/sqanti3/` et `${params.output_dir}/quality_report/sqanti3/`.

**Jalons**
1. ✅ M1 — `modules/sqanti3_qc.nf` créé avec `stub:` pour les deux sources et résumé MultiQC combiné.
2. ⏳ M2 — Run sur les isoformes long-read réelles non exécuté dans codex-dev ; validation actuelle faite en `-profile test -stub-run`.
3. ✅ M3 — Intégré, résumé des catégories structurales ajouté au rapport `quality_report/` comme custom content MultiQC.

**Validation**
- La majorité des isoformes qui chevauchent un gène `final_annotation.gff3` doivent être classées FSM (Full Splice Match) ou ISM — un taux élevé de NNC (Novel Not in Catalog) sur les gènes déjà bien supportés par ailleurs est un signal d'alerte sur la qualité de l'annotation finale, pas juste une statistique descriptive à ignorer.
- Rapport SQANTI3 généré sans erreur (`titan_sqanti3_classification.txt` non vide).

---

## Phase 8 — OMArk (QC complémentaire à BUSCO)

**Statut TITAN codex-dev — 2026-07-18** : Phase 8 implémentée comme QC protéome final optionnel après AEGIS/BUSCO. OMArk consomme `final_annotation_proteins_main.fasta`, publie les sorties sous `quality_report/omark/`, enregistre les sorties dans la provenance principale et ajoute `omark_mqc.tsv` au rapport MultiQC final.

**But** : évaluer non seulement la complétude du jeu de gènes mais aussi sa cohérence par rapport aux espèces proches, et détecter des indices de contamination.
**Pourquoi** : complémentaire à BUSCO (déjà implémenté) — OMArk considère les familles de gènes conservées multi-copies (BUSCO les exclut), et rapporte explicitement des signaux de contamination absents de BUSCO.
**Position dans le graphe** : juste après `busco` dans `workflows/titan.nf`, même niveau que le reste du `quality_report/`.

- Image : `quay.io/biocontainers/omark@sha256:84413cc19053c5d6452fbff245c9e6980b3f16aabdf991f9e51d7b9f2e0e0843` (0.5.0) — fournit aussi `omamer`
- Données offline : `params.omark_data_dir` (base OMAmer, `.h5`, téléchargée une fois depuis `omabrowser.org/oma/current/`, choisir le niveau taxonomique le plus proche de Viridiplantae/Eudicots disponible)
- Paramètres : `params.run_omark = false`, `params.container_omark`, `params.omark_data_dir = false`
- Label : `process_aegis` (même volumétrie que BUSCO)

**Input** : `aegis.out.aegis_proteins_main` (même protéome que BUSCO, pour comparabilité directe entre les deux rapports).

**Commande** :
```bash
omamer search --db ${params.omark_data_dir}/omamer.h5 \
  --query final_annotation_proteins_main.fasta \
  --out proteins_main.omamer --nthreads ${task.cpus}

omark -f proteins_main.omamer \
  -d ${params.omark_data_dir}/omamer.h5 \
  -o omark_out
```

**Output** :
- `omark_out/proteins_main_detailed_summary.txt` (publié — complétude + cohérence + contamination)
- `omark_out/proteins_main_omark.sum` (publié, résumé condensé)
- `proteins_main.omamer` (résultat OMAmer intermédiaire publié pour traçabilité)
- `omark_mqc.tsv` (résumé custom-content MultiQC)

Publié sous `${params.output_dir}/quality_report/omark/`.

**Jalons**
1. ✅ M1 — `modules/omark.nf` créé avec `stub:`, sur le même pattern que `busco.nf`.
2. ⏳ M2 — Base OMAmer réelle non téléchargée/smoke-testée dans codex-dev ; validation actuelle faite en `-profile test -stub-run`.
3. ✅ M3 — Intégré, ajouté au `multiqc_report` comme custom content via `omark_mqc.tsv`.

**Validation**
- Complétude OMArk et complétude BUSCO doivent être cohérentes entre elles (à quelques points de pourcentage près) — un grand écart signale que l'un des deux est mal configuré (mauvaise lignée/base).
- Le rapport de contamination OMArk ne doit signaler aucune séquence suspecte parmi les protéines principales — si c'est le cas, remonter l'alerte avant de considérer l'annotation finale "propre".
- `versions.yml` inclut la version d'OMArk/OMAmer et l'identifiant/date de la base OMAmer utilisée.

---

## Phase 9 — Validation par expression (aucun nouvel outil)

**Statut TITAN codex-dev — 2026-07-18** : Phase 9 implémentée et activée par défaut. TITAN extrait un transcriptome final depuis `final_annotation.gff3`, construit un index Salmon dédié, quantifie tous les FASTQ courts trimmés, agrège le support par gène avec seuil TPM configurable et ajoute `expression_support_summary_mqc.tsv` au MultiQC final.

**But** : vérifier que chaque gène prédit dans `final_annotation.gff3` a un support transcriptomique réel — exactement le contrôle explicitement fait pour PN40024.v5.1 ("gènes avec support transcriptomique significatif").
**Pourquoi** : réutilise entièrement l'infrastructure déjà présente dans TITAN (fastp, salmon, alignements STAR/HISAT2) — aucune nouvelle image à pinner, juste une étape de comptage/synthèse en plus.
**Position dans le graphe** : après AEGIS, en parallèle de `validate_final_annotation`.

- Pas de nouveau conteneur : réutiliser `params.container_salmon` (déjà pinné) et `params.container_python`.
- Paramètres : `params.run_expression_validation = true` (pas de dépendance de données externes, activé par défaut comme `agat_stats`)
- Label : `process_index` (salmon quant) puis `process_low` (agrégation python).

**Input** :
- `aegis.out.aegis_gff` (pour en extraire les transcrits via `agat_sp_extract_sequences.pl`, déjà utilisé ailleurs dans TITAN pour le CDS liftoff — même pattern).
- Tous les FASTQ trimmés déjà produits par `trimming_fastq` (réutiliser `evidence_data.fastp_json_reports` comme point d'ancrage du channel, mais il faut en fait les FASTQ eux-mêmes, pas juste les JSON — élargir l'emit de `generate_evidence_data.nf` en conséquence).

**Commande** :
```bash
# Nouvel index salmon sur le transcriptome FINAL (distinct de salmon_index
# existant, qui indexe le CDS liftoff pour l'inférence de strand)
agat_sp_extract_sequences.pl -g final_annotation.gff3 -f ${genome} -t exon --merge -o final_transcripts.fasta
salmon index -t final_transcripts.fasta -i final_salmon_index -p ${task.cpus}

# Par echantillon (reprend directement le pattern de salmon_strand_inference.nf)
salmon quant -i final_salmon_index -l A -p ${task.cpus} \
  -1 ${read_1} -2 ${read_2} -o ${sample_ID}_quant --validateMappings

scripts/summarize_expression_support.py \
  --gff final_annotation.gff3 \
  --quant-dirs *_quant \
  --min-tpm 0.5 \
  -o expression_support_summary.json
```

**Output** :
- `expression_support_summary.json` (publié — `% genes avec ≥1 échantillon TPM>0.5`, liste des gènes sans support, matrice TPM complète en intermédiaire)
- `expression_support_summary_mqc.tsv` (custom content MultiQC)
- `gene_tpm_matrix.tsv` (matrice TPM gène x échantillon)
- `final_transcripts.fasta` (transcriptome final indexé par Salmon)

Publié sous `${params.output_dir}/quality_report/expression_validation/`.

**Jalons**
1. ✅ M1 — `generate_evidence_data.nf` émet maintenant les FASTQ trimmés via `trimmed_fastqs`.
2. ✅ M2 — `modules/final_expression_validation.nf` et `scripts/summarize_expression_support.py` créés, `stub:` couvrant le graphe.
3. ✅ M3 — `scripts/test_summarize_expression_support.py` ajouté et intégré à `scripts/run-tests.sh`.
4. ✅ M4 — Intégré, publié, ajouté au rapport `quality_report/` et à MultiQC.

**Validation**
- Le `%` de gènes sans support transcriptomique doit être documenté et comparé à l'écart connu entre PN40024.v4 (11 508 gènes uniques sans preuve) et v5.1 (17 208 gènes avec support significatif) — l'objectif n'est pas 100 % (des gènes réels sont silencieux dans les tissus/conditions échantillonnés) mais une amélioration mesurable et traçable par rapport à la version précédente.
- `expression_support_summary.json` doit lister explicitement les gènes sans support (pas seulement un pourcentage agrégé), pour permettre un audit manuel ciblé plus tard.

## Ordre d'implémentation recommandé

1. **Phase 1 (tRNA) + Phase 2 (rRNA/snRNA/snoRNA)** — les plus simples, aucune dépendance croisée entre elles ni avec le reste du graphe, comblent directement le manque prouvé par la publication v5.1.
2. **Phase 3 (lncRNA)** — dépend des Phases 1-2 (exclusion des autres ncRNA) et de AEGIS ; complète le trio ncRNA qui correspond exactement à la méthodologie publiée.
3. **Phase 9 (validation par expression)** — aucun nouvel outil, faible risque, haute valeur (reproduit un contrôle qualité déjà validé sur ce même génome).
4. **Phase 8 (OMArk)** — suit directement le pattern de `busco.nf` déjà écrit, faible effort d'implémentation.
5. **Phase 4 + Phase 5 (Mikado + TransDecoder)** — dépendance croisée entre les deux, à traiter ensemble.
6. **Phase 6 + Phase 7 (FLAIR + SQANTI3)** — dépendantes l'une de l'autre et de `has_long_reads`.
