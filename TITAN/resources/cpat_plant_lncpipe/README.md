Plant-LncPipe CPAT-plant model files bundled for TITAN lncRNA candidate filtering.

Source repository: https://github.com/xuechantian/Plant-LncRNA-pipline

Downloaded from:
- https://raw.githubusercontent.com/xuechantian/Plant-LncRNA-pipline/master/Model/Plant_Hexamer.tsv
- https://raw.githubusercontent.com/xuechantian/Plant-LncRNA-pipline/master/Model/Plant.logit.RData

The published Plant-LncPipe cutoff is 0.46: coding probability >= 0.46 is treated as coding, and coding probability < 0.46 is treated as non-coding.

These files are generic plant models, not Vitis-specific CPAT models. TITAN therefore keeps the output named `lncrna_candidates.gff3` until a Vitis-trained CPAT model is available and validated.
