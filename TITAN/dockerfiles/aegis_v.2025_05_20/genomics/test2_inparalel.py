#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import os
from modules.tools import pickle_load, pickle_save
from modules.geneclasses import Annotation
from modules.geneclasses import Genome, Annotation_geneless

abinitio_evidences = pickle_load('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/ab_initio_evidences_corrected.pkl')
abinitio_evidences.detect_gene_overlaps_huge_file()
pickle_save('/media/tomslab2/Storage/Antonio/annotation_pipeline/Evidences_BLAST_hits/ab_initio_evidences_corrected_overlaps2.pkl', abinitio_evidences)
