from modules.genome import Genome
from os import listdir

genome = Genome('VCost', '/media/tomslab2/Storage/Antonio/TOMSBioLab/Jorge_families_analysis/PME_paper/TF_matrices/VCost_genome.fasta')

# genome.extract_peak_sequences(output_file_name='MYB14_top600_seqs.fa', DAPseq_output_file='/media/tomslab2/Storage/Antonio/TOMSBioLab/Jorge_families_analysis/PME_paper/TF_matrices/MYB14_r1r2_STSs_from_VCost_with_markers_TM29.csv', output_folder='/media/tomslab2/Storage/Antonio/TOMSBioLab/Jorge_families_analysis/PME_paper/TF_matrices/output')

files = listdir('/media/tomslab2/Storage/Antonio/TOMSBioLab/Colabs/Zenoni_plot/NACs_peaks_DAPseq/csvs')

for file in files:

    print(file)
    output_name = file.replace('.csv', '_peaks.fasta')
    genome.extract_peak_sequences(output_file_name=output_name, DAPseq_output_file=f'/media/tomslab2/Storage/Antonio/TOMSBioLab/Colabs/Zenoni_plot/NACs_peaks_DAPseq/csvs/{file}', output_folder='/media/tomslab2/Storage/Antonio/TOMSBioLab/Colabs/Zenoni_plot/NACs_peaks_DAPseq')