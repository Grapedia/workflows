#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from modules.genome import Genome

genome = Genome('genome', '/home/tomslab2/Descargas/bdistachyon_dapmod_scffree_organellefree.fasta')

genome.extract_peak_sequences(output_file_name='test.fa', DAPseq_output_file='/home/tomslab2/Descargas/DAP_000041_raw_dap00044_750_final_peaks.csv', output_folder='/home/tomslab2/Descargas', top=600)