"""
Snakefile for single-cell endoderm differentiation project


Author: Davis McCarthy
Affiliation: EMBL-EBI
Study: Single-cell endoderm differentiation project
Date: Sunday 12 June 2016
Run: snakemake --jobs 1000 --latency-wait 30 --cluster-config cluster.json --cluster 'bsub -q {cluster.queue} -J {cluster.name} -n {cluster.n} -R "rusage[mem={cluster.memory}]" -M {cluster.memory} -o {cluster.output} -e {cluster.error}'
add --prioritize flage to prioritize particular files
Latest modification:
  - todo

snakemake -fn --jobs 1000 --latency-wait 30 --cluster-config cluster.json --cluster 'bsub -q {cluster.queue} -J {cluster.name} -n {cluster.n} -R "rusage[mem={cluster.memory}]" -M {cluster.memory} -o {cluster.output} -e {cluster.error}' data_processed/scrnaseq/sceset.run_19776.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_18190.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_20287.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_20416.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_20450.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_20727.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_20759.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_20794.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_21241.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_21554.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_21672.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_21673.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_21843.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_21965.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_21999.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_22139.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_22194.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_22492.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_22606.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_22607.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_22710.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_22841.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_22944.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_23362.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_23794.salmon.preqc_gene.rds data_processed/scrnaseq/sceset.run_24086.salmon.preqc_gene.rds

STUDY ID: 3963
STUDY and PROJECT TITLE: Single Cell RNAseq at various stages of HiPSCs differentiating toward definitive endoderm and endoderm derived lineages.
Study name abbreviation: SC-RNAseq-definitive endoderm
PROJECT ID: 2010
HMDMC number (ethical approval): 13/042, 15_074
Project Cost Code: S0901
n
STUDY ID: 4262
STUDY and PROJECT TITLE: Single cell RNA-seq and bisulfite-seq at various stares of HiPSCs differentiating toward definitive endoderm derived lineages
PROJECT ID: 2218
HMDMC number (ethical approval): 13/042, 15_074
Project Cost Code: S0901

FACS data uploaded from Vallier group to Ian Streeter. Organised FACS data in:
/nfs/research2/hipsci/drop/hip-drop/tracked/endodiff/
and copied to
data_raw/facs

Raw data from Sanger first downloaded to:
data_raw/scrnaseq/run_{run}/cram

Raw data from Sanger sequencing core with some additions then uploaded to:
/hps/nobackup/hipsci/drop/hip-drop/incoming/stegle

From there it is moved into the HipSci tracked data area by Ian Streeter:
/hps/nobackup/hipsci/drop/hip-drop/tracked/

Transient files and analyses, and working files to share can be put in the scratch directory:
/hps/nobackup/hipsci/drop/hip-drop/scratch/
"""

import glob
import os
from subprocess import run
import pandas as pd
import re

configfile: "config.yaml"
shell.prefix("set -euo pipefail;") 

## REFERENCE FILES
fasta = os.path.join(config['references_dir'], config['human_tx_fasta'])
fasta_unzipped = config['human_gx_fasta']
fasta_dict = fasta_unzipped.replace('fa', 'dict')
fasta_idx = fasta_unzipped + '.fai'
kallisto_idx = expand('{basedir}{index}', basedir=config['references_dir'], index=config['kallisto_idx'])
fasta_GRCh38 = '/hps/nobackup/stegle/datasets/references/human/GRCh38/Homo_sapiens.GRCh38.rel85.cdna.all.ERCC92.fa.gz'
fasta_unzipped_GRCh38 = '/hps/nobackup/stegle/datasets/references/human/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.ERCC92.fa'
fasta_dict_GRCh38 = '/hps/nobackup/stegle/datasets/references/human/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.ERCC92.dict'
fasta_idx_GRCh38 = fasta_unzipped_GRCh38 + '.fai'
STAR_GENOME_DIR = '/hps/nobackup/stegle/datasets/references/human/STAR_GRCh37.75_ERCC'
star_genome_files = ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt', 'chrStart.txt', 'exonInfo.tab', 'Genome', 'genomeParameters.txt', 'SA', 'SAindex', 'sjdbInfo.txt', 'sjdbList.fromGTF.out.tab', 'sjdbList.out.tab', 'transcriptInfo.tab']
## variant files
dbSnpVcf = '/hps/nobackup/stegle/datasets/references/human/dbsnp_138.hg19.vcf.gz'
dbSnpVcfSmall = '/hps/nobackup/stegle/users/mjbonder/ref/GenotypingInformation/dbsnp_138.hg19.biallelicSNPs.HumanCoreExome12.Top1000ExpressedIpsGenes.Maf0.01.HWE0.0001.HipSci.vcf.gz'
reAlignmentIntervals = '/hps/nobackup/stegle/users/mjbonder/ref/GenotypingInformation/knownIndels.intervals'
knownIndelsMills = '/hps/nobackup/stegle/users/mjbonder/ref/GenotypingInformation/1000G_phase1.indels.hg19.sites.vcf.gz'
knownIndels100G = '/hps/nobackup/stegle/users/mjbonder/ref/GenotypingInformation/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz'
rRNAIntervals = '/hps/nobackup/stegle/datasets/references/human/rRNA.intervals'
## HipSci
HIPSCI_VCF = '/hps/nobackup/hipsci/scratch/genotypes/imputed/2017-03-27/INFO_0.4/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20170327.genotypes.allchr.vcf.gz'
hipsci_chr_vcf = []
hipsci_chr_vcf_idx = []
for i in range(1, 23):
    hipsci_chr_vcf.append('/hps/nobackup/hipsci/scratch/genotypes/imputed/2017-03-27/INFO_0.4/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20170327.genotypes.chr%s.vcf.gz' %str(i))
    hipsci_chr_vcf_idx.append('/hps/nobackup/hipsci/scratch/genotypes/imputed/2017-03-27/INFO_0.4/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20170327.genotypes.chr%s.vcf.gz.tbi' %str(i))


## define commands
kallisto_cmd = config['kallisto_cmd']
read_kallisto_to_scesets_cmd = 'src/preprocessing/read_kallisto_to_scesets.R'
kallisto_idx = os.path.join(config['references_dir'], config['kallisto_idx'])
kallisto_idx_GRCh38 = '/hps/nobackup/stegle/datasets/references/human/GRCh38/Homo_sapiens.GRCh38.rel85.cdna.all.ERCC92.kallisto_v0.43_idx'
python_cmd = '/nfs/software/stegle/users/davis/conda-envs/py3/bin/python'
star_cmd = '/nfs/software/stegle/bin/STAR'
#salmon_cmd = '/nfs/software/stegle/users/davis/Salmon-0.7.2_linux_x86_64/bin/salmon'
salmon_cmd = '/nfs/software/stegle/users/davis/Salmon-0.8.2_linux_x86_64/bin/salmon'
#salmon_idx = '/hps/nobackup/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.rel75.cdna.all.ERCC92.salmon_v0.7.1_idx'
salmon_idx = '/hps/nobackup/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.rel75.cdna.all.ERCC92.salmon_v0.8.2_idx'
#salmon_idx_GRCh38 = '/hps/nobackup/stegle/datasets/references/human/GRCh38/Homo_sapiens.GRCh38.rel85.cdna.all.ERCC92.salmon_v0.7.1_idx'
salmon_idx_GRCh38 = '/hps/nobackup/stegle/datasets/references/human/GRCh38/Homo_sapiens.GRCh38.rel85.cdna.all.ERCC92.salmon_v0.8.2_idx'
read_salmon_to_scesets_cmd = 'src/preprocessing/read_salmon_to_scesets.R'
#gatk_cmd = 'java -jar /nfs/software/stegle/users/davis/GenomeAnalysisTK.jar'
gatk_cmd = '/nfs/software/stegle/users/dseaton/java/jdk1.8.0_112/bin/java -Xmx48g -Xms24g -Djava.io.tmpdir=/hps/nobackup/hipsci/scratch/singlecell_endodiff/tmp -jar /hps/nobackup/stegle/users/mjbonder/tools/GATK/GenomeAnalysisTK.jar'
picard_cmd='/nfs/software/stegle/users/dseaton/java/jdk1.8.0_112/bin/java -Xmx48g -Xms24g -Djava.io.tmpdir=/hps/nobackup/hipsci/scratch/singlecell_endodiff/tmp -jar /nfs/software/stegle/users/dseaton/picard/picard.jar'
Rscript_cmd='/nfs/software/stegle/users/davis/conda-envs/py3/bin/Rscript'
R_cmd='/nfs/software/stegle/users/davis/conda-envs/py3/bin/R'

## parameter objects and samples
RUNS = ['run_19776', 'run_18190', 'run_20287', 'run_20416', 'run_20450', 'run_20727', 'run_20759', 'run_20794', 'run_21241', 'run_21554', 'run_21672', 'run_21673', 'run_21843', 'run_21965', 'run_21999', 'run_22139', 'run_22194', 'run_22492', 'run_22606', 'run_22607', 'run_22710', 'run_22841', 'run_22944', 'run_23362', 'run_23794', 'run_24086']
##RUNS = ['run_23794', 'run_24086']
RUNS_10x = ['cellranger201_count_22950_3_1000Genomes_hs37d5-ensembl_75_transcriptome', 'cellranger201_count_22950_5_1000Genomes_hs37d5-ensembl_75_transcriptome', 'cellranger201_count_22950_6_1000Genomes_hs37d5-ensembl_75_transcriptome', 'cellranger201_count_22950_7_1000Genomes_hs37d5-ensembl_75_transcriptome', 'cellranger201_count_22950_8_1000Genomes_hs37d5-ensembl_75_transcriptome', 'cellranger201_count_22951_3_1000Genomes_hs37d5-ensembl_75_transcriptome']
## avoid 'cellranger201_count_22950_4_1000Genomes_hs37d5-ensembl_75_transcriptome', for now, as v large number of dubious barcodes

#
LINES_DICT =  {'run_19776': ['eika_2;eipl_1;kuxp_1;pahc_4;podx_1;xugn_1'], 'run_18190': ['eika_2;eipl_1;kuxp_1;pahc_4;podx_1;xugn_1'], 'run_20287': ['joxm_1'], 'run_20416': ['eika_2;eipl_1;kuxp_1;pahc_4;podx_1;xugn_1'], 'run_20450': ['bima_1;lexy_1;oapg_5;pamv_3;rozh_4;fpdm_2'], 'run_20727': ['bezi_1;bubh_1;cuhk_1;eika_2;fpdm_2;oaaz_2'], 'run_20759': ['bezi_1;bubh_1;cuhk_1;eika_2;fpdm_2;oaaz_2'], 'run_20794': ['bezi_1;bubh_1;cuhk_1;eika_2;fpdm_2;oaaz_2'], 'run_21241': ['fafq_1;garx_2;hayt_1;sebz_1;sojd_3;wopl_1'], 'run_21554': ['dixh_2;fawm_2;koqx_1;naju_1;oebj_1;wigw_2;eoxi_6;fawm_2;iudw_4;oebj_1;oojs_1;pulk_1'], 'run_21672': ['fasu_2;kegd_2;zerv_8;zoio_2;xojn_3;fuai_1;eevy_7;oaqd_3;paab_4;sita_1;toss_3;zoio_2;heth_1;jogf_2;pelm_3;vass_1;wibj_2;zapk_3'], 'run_21673': ['fasu_2;kegd_2;zerv_8;zoio_2;xojn_3;fuai_1;eevy_7;iudw_4;kajh_3;tout_1;tavh_2'], 'run_21843': ['fafq_1;hiaf_2;iisa_3;joxm_1;lexy_1;wuye_2;oaqd_3;paab_4;sita_1;toss_3;zoio_2'], 'run_21965': ['babz_3;guyj_2;iisa_1;oikd_2;walu_1'], 'run_21999': ['eevy_7;fasu_2;iudw_4;kajh_3;tout_1;tavh_2'], 'run_22139': ['babz_3;guyj_2;iisa_1;oikd_2;walu_1'], 'run_22194': ['guyj_2;pulk_1;qayj_3;seru_1'], 'run_22492': ['guyj_2;pulk_1;qayj_3;seru_1;babz_3;iisa_1;oikd_2;walu_1'], 'run_22606': ['letw_1;olig_3;quls_2;rutc_2,sohd_3;vazt_1;iiyk_4;laey_4;miaj_6;poih_4;guyj_2;pulk_1;qayj_3;seru_1'], 'run_22607': ['aowh_2;keui_1;meue_4;naah_2;poih_4;vils_1'], 'run_22710': ['letw_1;olig_3;quls_2;rutc_2,sohd_3;vazt_1;iiyk_4;laey_4;miaj_6;poih_4'], 'run_22841': ['aowh_2;oicx_6;sehl_6;suop_5;wahn_1;wibj_2'], 'run_22944': ['aowh_2;oicx_6;sehl_6;suop_5;wahn_1;wibj_2;cicb_2;cuhk_2;hegp_3;lepk_1;ueah_1;veku_2'], 'run_23362': ['cicb_2;cuhk_2;hegp_3;lepk_1;ueah_1;veku_2'], 'run_23794': ['guss_1;lepk_1;mita_2;nocf_2;oibg_1;datg_2;feec_2;nudd_1;paab_4;qorq_2'], 'run_24086': ['guss_1;lepk_1;mita_2;nocf_2;oibg_1;datg_2;feec_2;nudd_1;paab_4;qorq_2']}
SAMPLES_DICT = {}
for run in RUNS:
    tmp = glob.glob("data_raw/scrnaseq/{0}/cram/*.cram".format(run))
    SAMPLES_DICT[run] = [os.path.basename(w).replace('.cram', '') for w in tmp]
SAMPLES = glob.glob("data_raw/scrnaseq/*/cram/*.cram")
SAMPLES = [os.path.basename(w).replace('.cram', '') for w in SAMPLES]

## targets
star_genome_output = expand('{genome_dir}/{genome_files}', genome_dir=STAR_GENOME_DIR, genome_files=star_genome_files)
kallisto_results = []
kallisto_results_GRCh38 = []
salmon_results = []
salmon_results_GRCh38 = []
kallisto_results_dict = {}
kallisto_results_GRCh38_dict = {}
salmon_results_dict = {}
salmon_results_GRCh38_dict = {}
fastqc_html_reports = []
scesets = []
scater_first_html_reports = []
star_bam_output = []
picard_readgroups_bam_output = []
picard_dedup_bam_output = []
gatk_split_bam_output = []
gatk_filtered_vcf_files = []
gatk_ase_output = []
variant_donor_id_files = []
variant_donor_id_files_dict = {}
donor_id_all_files = []
sample_table_files = []
multiqc_reports = []
seq_metadata_output = []
seq_metadata_files_dict = {}
feature_counts_uniq_output_dict = {}
feature_counts_uniq_output = []
feature_counts_nonuniq_output_dict = {}
feature_counts_nonuniq_output = []
picard_metrics = []
for run in RUNS:
    kallisto_results.append(expand('data_raw/scrnaseq/{run}/quant_kallisto/{sample}/abundance.tsv', run=run, sample=SAMPLES_DICT[run]))
    kallisto_results_dict[run] = expand('data_raw/scrnaseq/{run}/quant_kallisto/{sample}/abundance.tsv', run=run, sample=SAMPLES_DICT[run])
    kallisto_results_GRCh38.append(expand('data_raw/scrnaseq/{run}/quant_kallisto_GRCh38/{sample}/abundance.tsv', run=run, sample=SAMPLES_DICT[run]))
    kallisto_results_GRCh38_dict[run] = expand('data_raw/scrnaseq/{run}/quant_kallisto_GRCh38/{sample}/abundance.tsv', run=run, sample=SAMPLES_DICT[run])
    salmon_results.append(expand('data_raw/scrnaseq/{run}/quant_salmon/{sample}/quant.sf', run=run, sample=SAMPLES_DICT[run]))
    salmon_results_dict[run] = expand('data_raw/scrnaseq/{run}/quant_salmon/{sample}/quant.sf', run=run, sample=SAMPLES_DICT[run])
    salmon_results_GRCh38.append(expand('data_raw/scrnaseq/{run}/quant_salmon_GRCh38/{sample}/quant.sf', run=run, sample=SAMPLES_DICT[run]))
    salmon_results_GRCh38_dict[run] = expand('data_raw/scrnaseq/{run}/quant_salmon_GRCh38/{sample}/quant.sf', run=run, sample=SAMPLES_DICT[run])
    fastqc_html_reports.append(expand('reports/fastqc/scrnaseq/{run}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr_fastqc.html', run=run, sample=SAMPLES_DICT[run]))
    # scesets.append(expand('data_processed/scrnaseq/sceset.{run}.kallisto.preqc_tx.rds', run=run))
    # scesets.append(expand('data_processed/scrnaseq/sceset.{run}.kallisto.preqc_gene.rds', run=run))
    scesets.append(expand('data_processed/scrnaseq/sceset.{run}.salmon.preqc_tx.rds', run=run))
    scesets.append(expand('data_processed/scrnaseq/sceset.{run}.salmon.preqc_gene.rds', run=run))
    scater_first_html_reports.append(expand('reports/first_qc/{run}.{quant_tool}.first_qc.html', run=run, quant_tool = ['kallisto', 'salmon']))
    star_bam_output.append(expand('data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.out.bam', run=run, sample=SAMPLES_DICT[run]))
    picard_readgroups_bam_output.append(expand('data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.rgadded.bam', run=run, sample=SAMPLES_DICT[run]))
    picard_dedup_bam_output.append(expand('data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.dedup.bam', run=run, sample=SAMPLES_DICT[run]))
    gatk_split_bam_output.append(expand('data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.bam', run=run, sample=SAMPLES_DICT[run]))
    gatk_filtered_vcf_files.append(expand('data_raw/scrnaseq/{run}/vcf/{sample}.filtered.vcf', run=run, sample=SAMPLES_DICT[run]))
    gatk_ase_output.append(expand('data_raw/scrnaseq/{run}/ase/high_thresh/{sample}.ase.highthresh.tsv', run=run, sample=SAMPLES_DICT[run]))
    gatk_ase_output.append(expand('data_raw/scrnaseq/{run}/ase/low_thresh/{sample}.ase.lowthresh.tsv', run=run, sample=SAMPLES_DICT[run]))
    donor_id_tmp = expand('data_raw/scrnaseq/{run}/donor_id/{sample}.donor_id.csv', run=run, sample=SAMPLES_DICT[run])
    variant_donor_id_files_dict[run] = donor_id_tmp
    variant_donor_id_files.append(donor_id_tmp)
    seq_metadata_files_dict[run] = expand('data_raw/scrnaseq/{run}/meta/{sample}.meta', run=run, sample=SAMPLES_DICT[run])
    feature_counts_uniq_output_dict[run] = expand('data_raw/scrnaseq/{run}/featureCounts/unique_counts/{sample}.gene.counts.unique.tsv', run=run, sample=SAMPLES_DICT[run])
    feature_counts_uniq_output.append(expand('data_raw/scrnaseq/{run}/featureCounts/unique_counts/{sample}.gene.counts.unique.tsv', run=run, sample=SAMPLES_DICT[run]))
    feature_counts_nonuniq_output_dict[run] = expand('data_raw/scrnaseq/{run}/featureCounts/total_counts/{sample}.gene.counts.tsv', run=run, sample=SAMPLES_DICT[run])
    feature_counts_nonuniq_output.append(expand('data_raw/scrnaseq/{run}/featureCounts/total_counts/{sample}.gene.counts.tsv', run=run, sample=SAMPLES_DICT[run]))
    picard_metrics.append(expand('data_raw/scrnaseq/{run}/picard_metrics/{sample}.picard.rna.metrics', run=run, sample=SAMPLES_DICT[run]))

for i in range(len(RUNS)):
    donor_id_all_files.append("data_processed/donor_id/donor_id_all.{0}.csv".format(RUNS[i]))
    multiqc_reports.append("reports/multiqc/{0}/multiqc_report.{0}.html".format(RUNS[i]))
    sample_table_files.append("metadata/scrnaseq/samples.{0}.tab".format(RUNS[i]))
    seq_metadata_output.append("metadata/scrnaseq/seq_metadata.{0}.tsv".format(RUNS[i]))


filtered_vcf_10x = []
donor_id_10x = []
barcodes_10x = {}
barcode_bams_10x = {}
variant_donor_id_files_10x_dict = {}
for run in RUNS_10x:
    # filtered_vcf_10x.append(expand('data_raw/10x/{run}/vcf/pooled10x.filtered.vcf.gz', run=run))
    # filtered_vcf_10x.append(expand('data_raw/10x/{run}/vcf/pooled10x.filtered.vcf.gz.csi', run=run))
    donor_id_10x.append(expand('data_processed/donor_id_10x/donor_id_all.{run}.csv', run=run))
    fname = expand('data_raw/10x/{run}/filtered_gene_bc_matrices/1000Genomes_hs37d5/barcodes.tsv', run=run)
    with open(fname[0]) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    barcodes_10x[run] = content
    barcode_bams_10x[run] = expand('data_raw/10x/{run}/bam_barcoded/{barcode}.bam', run=run, barcode=content)
    variant_donor_id_files_10x_dict[run] = expand('data_raw/10x/{run}/donor_id/10x.{sample}.donor_id.csv', run=run, sample=content)
    filtered_vcf_10x.append(expand('data_raw/10x/{run}/vcf/10x.{barcode}.filtered.vcf.gz', run=run, barcode=content))
    filtered_vcf_10x.append(expand('data_raw/10x/{run}/vcf/10x.{barcode}.filtered.vcf.gz.csi', run=run, barcode=content))
filtered_vcf_10x = [filename for elem in filtered_vcf_10x for filename in elem]
donor_id_10x = [filename for elem in donor_id_10x for filename in elem]


## flatten these lists
kallisto_results = [filename for elem in kallisto_results for filename in elem]
kallisto_results_GRCh38 = [filename for elem in kallisto_results_GRCh38 for filename in elem]
salmon_results = [filename for elem in salmon_results for filename in elem]
salmon_results_GRCh38 = [filename for elem in salmon_results_GRCh38 for filename in elem]
fastqc_html_reports = [filename for elem in fastqc_html_reports for filename in elem]
#scesets = [filename for elem in scesets for filename in elem]
scater_first_html_reports = [filename for elem in scater_first_html_reports for filename in elem]
star_bam_output = [filename for elem in star_bam_output for filename in elem]
picard_readgroups_bam_output = [filename for elem in picard_readgroups_bam_output for filename in elem]
picard_dedup_bam_output = [filename for elem in picard_dedup_bam_output for filename in elem]
gatk_split_bam_output = [filename for elem in gatk_split_bam_output for filename in elem]
gatk_filtered_vcf_files = [filename for elem in gatk_filtered_vcf_files for filename in elem]
gatk_ase_output = [filename for elem in gatk_ase_output for filename in elem]
variant_donor_id_files = [filename for elem in variant_donor_id_files for filename in elem]
feature_counts_uniq_output = [filename for elem in feature_counts_uniq_output for filename in elem]
feature_counts_nonuniq_output = [filename for elem in feature_counts_nonuniq_output for filename in elem]
picard_metrics = [filename for elem in picard_metrics for filename in elem]
#donor_id_all_files = [filename for elem in donor_id_all_files for filename in elem]
## define star bam index files
star_bam_index_output = [x + '.bai' for x in star_bam_output]
## input files for summarising data at experiment level
expt_sceset_reqd_files = {'expt01': ['data_processed/donor_id/donor_id_all.run_18190.csv', 'metadata/scrnaseq/seq_metadata.run_18190.tsv', 'data_processed/scrnaseq/sceset.run_18190.salmon.preqc_gene.rds', 'data_processed/donor_id/donor_id_all.run_19776.csv', 'metadata/scrnaseq/seq_metadata.run_19776.tsv', 'data_processed/scrnaseq/sceset.run_19776.salmon.preqc_gene.rds', 'data_processed/donor_id/donor_id_all.run_20416.csv', 'metadata/scrnaseq/seq_metadata.run_20416.tsv', 'data_processed/scrnaseq/sceset.run_20416.salmon.preqc_gene.rds'], 'expt03': ['data_processed/donor_id/donor_id_all.run_20759.csv', 'metadata/scrnaseq/seq_metadata.run_20759.tsv', 'data_processed/scrnaseq/sceset.run_20759.salmon.preqc_gene.rds'], 'expt06': ['data_processed/donor_id/donor_id_all.run_20450.csv', 'metadata/scrnaseq/seq_metadata.run_20450.tsv', 'data_processed/scrnaseq/sceset.run_20450.salmon.preqc_gene.rds'], 'expt08': ['data_processed/donor_id/donor_id_all.run_20287.csv', 'metadata/scrnaseq/seq_metadata.run_20287.tsv', 'data_processed/scrnaseq/sceset.run_20287.salmon.preqc_gene.rds'], 'expt09': ['data_processed/donor_id/donor_id_all.run_21843.csv', 'metadata/scrnaseq/seq_metadata.run_21843.tsv', 'data_processed/scrnaseq/sceset.run_21843.salmon.preqc_gene.rds'], 'expt10': ['data_processed/donor_id/donor_id_all.run_21241.csv', 'metadata/scrnaseq/seq_metadata.run_21241.tsv', 'data_processed/scrnaseq/sceset.run_21241.salmon.preqc_gene.rds'], 'expt12': ['data_processed/donor_id/donor_id_all.run_21672.csv', 'metadata/scrnaseq/seq_metadata.run_21672.tsv', 'data_processed/scrnaseq/sceset.run_21672.salmon.preqc_gene.rds'], 'expt18': ['data_processed/donor_id/donor_id_all.run_21843.csv', 'metadata/scrnaseq/seq_metadata.run_21843.tsv', 'data_processed/scrnaseq/sceset.run_21843.salmon.preqc_gene.rds', 'data_processed/donor_id/donor_id_all.run_21672.csv', 'metadata/scrnaseq/seq_metadata.run_21672.tsv', 'data_processed/scrnaseq/sceset.run_21672.salmon.preqc_gene.rds'], 'expt19': ['data_processed/donor_id/donor_id_all.run_21673.csv', 'metadata/scrnaseq/seq_metadata.run_21673.tsv', 'data_processed/scrnaseq/sceset.run_21673.salmon.preqc_gene.rds', 'data_processed/donor_id/donor_id_all.run_21672.csv', 'metadata/scrnaseq/seq_metadata.run_21672.tsv', 'data_processed/scrnaseq/sceset.run_21672.salmon.preqc_gene.rds'], 'expt20': ['data_processed/donor_id/donor_id_all.run_21673.csv', 'metadata/scrnaseq/seq_metadata.run_21673.tsv', 'data_processed/scrnaseq/sceset.run_21673.salmon.preqc_gene.rds', 'data_processed/donor_id/donor_id_all.run_21999.csv', 'metadata/scrnaseq/seq_metadata.run_21999.tsv', 'data_processed/scrnaseq/sceset.run_21999.salmon.preqc_gene.rds'], 'expt21': ['data_processed/donor_id/donor_id_all.run_21554.csv', 'metadata/scrnaseq/seq_metadata.run_21554.tsv', 'data_processed/scrnaseq/sceset.run_21554.salmon.preqc_gene.rds'], 'expt22': ['data_processed/donor_id/donor_id_all.run_21554.csv', 'metadata/scrnaseq/seq_metadata.run_21554.tsv', 'data_processed/scrnaseq/sceset.run_21554.salmon.preqc_gene.rds'], 'expt22': ['data_processed/donor_id/donor_id_all.run_21554.csv', 'metadata/scrnaseq/seq_metadata.run_21554.tsv', 'data_processed/scrnaseq/sceset.run_21554.salmon.preqc_gene.rds'], 'expt23': ['data_processed/donor_id/donor_id_all.run_22606.csv', 'metadata/scrnaseq/seq_metadata.run_22606.tsv', 'data_processed/scrnaseq/sceset.run_22606.salmon.preqc_gene.rds', 'data_processed/donor_id/donor_id_all.run_22194.csv', 'metadata/scrnaseq/seq_metadata.run_22194.tsv', 'data_processed/scrnaseq/sceset.run_22194.salmon.preqc_gene.rds', 'data_processed/donor_id/donor_id_all.run_22492.csv', 'metadata/scrnaseq/seq_metadata.run_22492.tsv', 'data_processed/scrnaseq/sceset.run_22492.salmon.preqc_gene.rds'], 'expt24': ['data_processed/donor_id/donor_id_all.run_22139.csv', 'metadata/scrnaseq/seq_metadata.run_22139.tsv', 'data_processed/scrnaseq/sceset.run_22139.salmon.preqc_gene.rds', 'data_processed/donor_id/donor_id_all.run_21965.csv', 'metadata/scrnaseq/seq_metadata.run_21965.tsv', 'data_processed/scrnaseq/sceset.run_21965.salmon.preqc_gene.rds', 'data_processed/donor_id/donor_id_all.run_22492.csv', 'metadata/scrnaseq/seq_metadata.run_22492.tsv', 'data_processed/scrnaseq/sceset.run_22492.salmon.preqc_gene.rds'], 'expt27': ['data_processed/donor_id/donor_id_all.run_22606.csv', 'metadata/scrnaseq/seq_metadata.run_22606.tsv', 'data_processed/scrnaseq/sceset.run_22606.salmon.preqc_gene.rds', 'data_processed/donor_id/donor_id_all.run_22710.csv', 'metadata/scrnaseq/seq_metadata.run_22710.tsv', 'data_processed/scrnaseq/sceset.run_22710.salmon.preqc_gene.rds'], 'expt28': ['data_processed/donor_id/donor_id_all.run_22606.csv', 'metadata/scrnaseq/seq_metadata.run_22606.tsv', 'data_processed/scrnaseq/sceset.run_22606.salmon.preqc_gene.rds', 'data_processed/donor_id/donor_id_all.run_22710.csv', 'metadata/scrnaseq/seq_metadata.run_22710.tsv', 'data_processed/scrnaseq/sceset.run_22710.salmon.preqc_gene.rds'], 'expt29': ['data_processed/donor_id/donor_id_all.run_22607.csv', 'metadata/scrnaseq/seq_metadata.run_22607.tsv', 'data_processed/scrnaseq/sceset.run_22607.salmon.preqc_gene.rds'], 'expt33': ['data_processed/donor_id/donor_id_all.run_23794.csv', 'metadata/scrnaseq/seq_metadata.run_23794.tsv', 'data_processed/scrnaseq/sceset.run_23794.salmon.preqc_gene.rds']}
## experiment sceset files
expt_scesets = []
for expt in ['expt01', 'expt03', 'expt06', 'expt08', 'expt09', 'expt10', 'expt12', 'expt18', 'expt19', 'expt20', 'expt21', 'expt22', 'expt23', 'expt24', 'expt27', 'expt28', 'expt29', 'expt33']:
    expt_scesets.append('data_processed/{0}/sceset_{0}_salmon_allmeta_allcells.rds'.format(expt))


rule all:
    input:
        # 'data_raw/grm/hipsci_grm.rel.gz', 'data_processed/eqtl/REL-2016-09.Kpop_maf05.h5',
        donor_id_all_files,
        seq_metadata_output,
        scesets, #scater_first_html_reports, 
        # multiqc_reports,
        # expt_scesets,
        #feature_counts_uniq_output, feature_counts_nonuniq_output,
        # picard_metrics,
        #filtered_vcf_10x, donor_id_10x


### 10x rules ###


rule identify_donor_collect_10x:
    input:
        files=lambda wildcards: variant_donor_id_files_10x_dict[wildcards.run]
    output:
        'data_processed/donor_id_10x/donor_id_all.{run}.csv'
    params:
        dir='data_raw/10x/{run}/donor_id/',
    run:
        import pandas as pd
        import glob
        import os
        input_files = glob.glob(os.path.join(params.dir, '*.donor_id.csv'))
        df_list = {}
        for infile in input_files:
            df_tmp = pd.read_csv(infile)
            df_list[infile] = df_tmp
        df_out = pd.concat(df_list)
        df_out.to_csv(output[0], index=False)
 

rule identify_donor_runs_10x:
    input:
        sc_vcf='data_raw/10x/{run}/vcf/10x.{sample}.filtered.vcf.gz',
        hipsci_vcf='data_raw/10x/{run}/vcf/10x.{sample}.filtered.hipsci.overlap.vcf.gz',
        hipsci_vcf_idx='data_raw/10x/{run}/vcf/10x.{sample}.filtered.hipsci.overlap.vcf.gz.csi'
    output:
        'data_raw/10x/{run}/donor_id/10x.{sample}.donor_id.csv'
    params:
        prefix='data_raw/10x/{run}/donor_id/10x.{sample}.donor_id',
        lines='cicb_2;cuhk_2;hegp_3;lepk_1;ueah_1;veku_2'
    shell:
        '{Rscript_cmd} src/R/identify_donor_small_vcf.R --input_file "{input.sc_vcf}" '
        '--donor_lines "{params.lines}" '
        '--donor_vcf {input.hipsci_vcf} '
        '--output_prefix "{params.prefix}" '


rule index_hipsci_overlap_vcf_10x:
    input:
        'data_raw/10x/{run}/vcf/10x.{sample}.filtered.hipsci.overlap.vcf.gz'
    output:
        temp('data_raw/10x/{run}/vcf/10x.{sample}.filtered.hipsci.overlap.vcf.gz.csi')
    shell:
        'bcftools index {input}'


rule filter_hipsci_overlap_variants_10x:
    input:
        sc_vcf='data_raw/10x/{run}/vcf/10x.{sample}.filtered.vcf.gz',
        sc_vcf_idx='data_raw/10x/{run}/vcf/10x.{sample}.filtered.vcf.gz.csi',
        variant_list='data_raw/10x/{run}/vcf/10x.{sample}.variant_list.txt',
        hipsci_vcf=HIPSCI_VCF
    output:
        tmp=temp('data_raw/10x/{run}/vcf/10x.{sample}.tmp.vcf.gz'),
        vcf=temp('data_raw/10x/{run}/vcf/10x.{sample}.filtered.hipsci.overlap.vcf.gz')
    shell:
        """
        bcftools view -o {output.tmp} -O z -l 9 -R {input.variant_list} {input.hipsci_vcf}
        vcf-sort {output.tmp} | bgzip -c > {output.vcf}
        """


rule singlecell_variant_list_10x:
    input:
        'data_raw/10x/{run}/vcf/10x.{sample}.filtered.vcf.gz'
    output:
        temp('data_raw/10x/{run}/vcf/10x.{sample}.variant_list.txt')
    shell:
        """
        set +euo pipefail
        echo -e "1\\t1\\tA\\tC" > {output}
        bcftools view -O v {input} | grep -v ^# | awk \'{{sub(/chr/,""); print $1"\\t"$2"\\t"$4"\\t"$5}}\' >> {output}
        set -euo pipefail
        """


rule index_bgzip_vcf_10x:
    input:
        'data_raw/10x/{run}/vcf/10x.{sample}.filtered.vcf.gz'
    output:
        'data_raw/10x/{run}/vcf/10x.{sample}.filtered.vcf.gz.csi'
    shell:
        'bcftools index {input}'


rule bgzip_vcf_10x:
    input:
        'data_raw/10x/{run}/vcf/10x.{sample}.filtered.vcf'
    output:
        'data_raw/10x/{run}/vcf/10x.{sample}.filtered.vcf.gz'
    shell:
        'vcf-sort {input} | bgzip -c > {output}'


rule filter_variants_gatk_10x:
    input:
        vcf='data_raw/10x/{run}/vcf/10x.{sample}.unfiltered.vcf',
        fasta=fasta_unzipped,
        fai=fasta_idx,
        fa_dict=fasta_dict
    output:
        temp('data_raw/10x/{run}/vcf/10x.{sample}.filtered.vcf')
    shell:
        '{gatk_cmd} -T VariantFiltration -R {input.fasta} -V {input.vcf} '
        '-window 35 -cluster 3 -filterName FS -filter "FS > 30.0" '
        '-filterName QD -filter "QD < 2.0" -o {output}'


rule call_variants_gatk_10x_pooled:
    input:
        fasta=fasta_unzipped,
        fai=fasta_idx,
        fa_dict=fasta_dict,
        dbSnp=dbSnpVcf,
        dbSnpSmall=dbSnpVcfSmall,
        bam='data_raw/10x/{run}/star_possorted_genome_bam.split.realigned.bqsr.bam'
    output:
        'data_raw/10x/{run}/vcf/pooled10x.unfiltered.vcf'
    shell:
        '{gatk_cmd} -T HaplotypeCaller -R {input.fasta} -I {input.bam} '
        '-dontUseSoftClippedBases '
        '-D {input.dbSnp} -gt_mode GENOTYPE_GIVEN_ALLELES '
        '-alleles {input.dbSnpSmall} -L {input.dbSnpSmall} -o {output}' 


rule call_variants_gatk_10x_barcode:
    input:
        fasta=fasta_unzipped,
        fai=fasta_idx,
        fa_dict=fasta_dict,
        dbSnp=dbSnpVcf,
        dbSnpSmall= dbSnpVcfSmall,
        bam='data_raw/10x/{run}/bam_barcoded/{barcode}.rgadded.bam',
        bai='data_raw/10x/{run}/bam_barcoded/{barcode}.rgadded.bai',
    output:
        'data_raw/10x/{run}/vcf/10x.{barcode}.unfiltered.vcf'
    shell:
        '{gatk_cmd} -T HaplotypeCaller -R {input.fasta} -I {input.bam} '
        '-dontUseSoftClippedBases '
        '-D {input.dbSnp} -gt_mode GENOTYPE_GIVEN_ALLELES '
        '-alleles {input.dbSnpSmall} -L {input.dbSnpSmall} -o {output}' 


rule index_barcode_bams_10x:
    input:
        'data_raw/10x/{run}/bam_barcoded/{barcode}.rgadded.bam'
    output:
        'data_raw/10x/{run}/bam_barcoded/{barcode}.rgadded.bai'
    shell:
        '{picard_cmd} BuildBamIndex I={input} '


rule picard_read_groups_barcode_10x:
    input:
        bam='data_raw/10x/{run}/bam_barcoded/{barcode}.raw.bam'
    output:
        'data_raw/10x/{run}/bam_barcoded/{barcode}.rgadded.bam'
    shell:
        '{picard_cmd} AddOrReplaceReadGroups I={input.bam} O={output} SO=coordinate '
        'RGID={wildcards.barcode} RGLB={wildcards.barcode} '
        'RGPL=ILLUMINA RGPU=MACHINE1 RGSM={wildcards.barcode}'


# split bam into separate files by barcode
# samtools view  in.bam | cut -f 12 | tr "\t" "\n"  | grep  "^CB:Z:"  | cut -d ':' -f 3 | sort | uniq | while read S; do samtools view -h in.bam |  awk -v tag="CR:Z:$S" '($0 ~ /^@/ || index($0,tag)>0)' > ${S}.sam ; done

rule split_bam_by_barcode_10x:
    input:
        bam='data_raw/10x/{run}/star_possorted_genome_bam.split.realigned.bqsr.bam'
    output:
        sam=temp('data_raw/10x/{run}/bam_barcoded/{barcode}.raw.sam'),
        bam='data_raw/10x/{run}/bam_barcoded/{barcode}.raw.bam'
    shell:
        """
        samtools view -h {input.bam} | awk -v tag="CB:Z:{wildcards.barcode}" '($0 ~ /^@/ || index($0,tag)>0)' > {output.sam}
        samtools view -Sb {output.sam} > {output.bam}
        """


rule picard_qc_metrics_10x:
    input:
        bam='data_raw/10x/{run}/star_possorted_genome_bam.split.realigned.bqsr.bam',
        ref_flat='/hps/nobackup/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.75_ERCC.ref_flat.txt',
        ribo_intervals=rRNAIntervals
    output:
        metrics='data_raw/10x/{run}/picard.rna.metrics',
        #plot_pdf='figures/picard_cov_vs_pos/{run}/{sample}.picard.cov.pos.pdf'
    shell:
        '{picard_cmd} CollectRnaSeqMetrics I={input.bam} O={output.metrics} '
        'REF_FLAT={input.ref_flat} STRAND_SPECIFICITY=NONE '
        'RIBOSOMAL_INTERVALS={input.ribo_intervals} '
        #'CHART_OUTPUT={output.plot_pdf}'


rule recalibrated_writer_gatk_10x:
    input:
        fasta=fasta_unzipped,
        bam='data_raw/10x/{run}/star_possorted_genome_bam.split.realigned.bam',
        bqsr='data_raw/10x/{run}/star_possorted_genome_bam.split.realigned.bqsr'
    output:
        'data_raw/10x/{run}/star_possorted_genome_bam.split.realigned.bqsr.bam'
    shell:
        '{gatk_cmd} -T PrintReads -R {input.fasta} -I {input.bam} '
        '-BQSR {input.bqsr} -nct 2 '
        '-o {output}' 


rule base_recalibrator_gatk_10x:
    input:
        fasta=fasta_unzipped,
        dbSnp=dbSnpVcf,
        bam='data_raw/10x/{run}/star_possorted_genome_bam.split.realigned.bam',
        known1=knownIndelsMills,
        known2=knownIndels100G
    output:
        temp('data_raw/10x/{run}/star_possorted_genome_bam.split.realigned.bqsr')
    shell:
        '{gatk_cmd} -T BaseRecalibrator -R {input.fasta} -I {input.bam} '
        '-knownSites {input.known1} -knownSites {input.known2} -knownSites {input.dbSnp} '
        '-nct 2 '
        '-o {output}' 


rule indel_realignment_gatk_10x:
    input:
        fasta=fasta_unzipped,
        bam='data_raw/10x/{run}/star_possorted_genome_bam.split.bam',
        targetIntervals=reAlignmentIntervals,
        known1=knownIndelsMills,
        known2=knownIndels100G
    output:
        temp('data_raw/10x/{run}/star_possorted_genome_bam.split.realigned.bam')
    shell:
        '{gatk_cmd} -T IndelRealigner -R {input.fasta} -I {input.bam} '
        '-targetIntervals {input.targetIntervals} -known {input.known1} -known {input.known2} '
        '-U ALLOW_N_CIGAR_READS --consensusDeterminationModel KNOWNS_ONLY --LODThresholdForCleaning 0.4  '
        '-o {output}' 


rule split_n_trim_gatk_10x:
    input:
        fasta=fasta_unzipped,
        fai=fasta_idx,
        fa_dict=fasta_dict,
        bam='data_raw/10x/{run}/star_possorted_genome_bam.dedup.bam',
        bai='data_raw/10x/{run}/star_possorted_genome_bam.dedup.bai'
    output:
        temp('data_raw/10x/{run}/star_possorted_genome_bam.split.bam')
    shell:
        '{gatk_cmd} -T SplitNCigarReads -R {input.fasta} -I {input.bam} -o {output} '
        '-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 '
        '-U ALLOW_N_CIGAR_READS '


rule picard_mark_dups_10x:
    input:
        bam='data_raw/10x/{run}/star_possorted_genome_bam.tagsadded_possort.bam',
        bai='data_raw/10x/{run}/star_possorted_genome_bam.tagsadded_possort.bai'
    output:
        bam=temp('data_raw/10x/{run}/star_possorted_genome_bam.dedup.bam'),
        bai=temp('data_raw/10x/{run}/star_possorted_genome_bam.dedup.bai'),
        metrics='data_raw/10x/{run}/picard.dedup.output.metrics'
    shell:
        '{picard_cmd} MarkDuplicates I={input.bam} O={output.bam} CREATE_INDEX=true '
        'VALIDATION_STRINGENCY=SILENT M={output.metrics} '


rule picard_sort_coord_10x:
    input:
        bam='data_raw/10x/{run}/star/star.2pass.Aligned.sortedByCoord.tagsadded.bam'
    output:
        bam='data_raw/10x/{run}/star_possorted_genome_bam.tagsadded_possort.bam',
        bai='data_raw/10x/{run}/star_possorted_genome_bam.tagsadded_possort.bai'
    shell:
        """
        {picard_cmd} SortSam I={input.bam} O={output.bam} SO=coordinate CREATE_INDEX=true
        """


rule picard_merge_bam_alignment_10x:
    input:
        aligned_bam='data_raw/10x/{run}/star_possorted_genome_bam.qnamesort.bam',
        unmapped_bam='data_raw/10x/{run}/unmapped_bam.bam',
        fasta=fasta_unzipped
    output:
        temp('data_raw/10x/{run}/star/star.2pass.Aligned.sortedByCoord.tagsadded.bam')
    shell:
        '{picard_cmd} MergeBamAlignment REFERENCE_SEQUENCE={input.fasta} '
        'UNMAPPED_BAM={input.unmapped_bam} ALIGNED_BAM={input.aligned_bam} '
        'OUTPUT={output} INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false '
        'VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true '


rule picard_sort_qname_10x:
    input:
        bam1='data_raw/10x/{run}/star/star.2pass.Aligned.sortedByCoord.out.bam',
    output:
        bam1=temp('data_raw/10x/{run}/star_possorted_genome_bam.qnamesort.bam'),
    shell:
        """
        {picard_cmd} SortSam I={input.bam1} O={output.bam1} SO=queryname
        """


rule index_star_bams_10x:
    input:
        'data_raw/10x/{run}/star/star.2pass.Aligned.sortedByCoord.out.bam'
    output:
        temp('data_raw/10x/{run}/star/star.2pass.Aligned.sortedByCoord.out.bam.bai')
    shell:
        'samtools index {input} '


rule make_unmapped_bam_10x:
    input:
        'data_raw/10x/{run}/possorted_genome_bam.bam'    
    output:
        temp('data_raw/10x/{run}/unmapped_bam.bam')
    shell:
        """
        {picard_cmd} RevertSam \
        I={input} \
        O={output} \
        SANITIZE=false \
        MAX_DISCARD_FRACTION=0.005 \
        ATTRIBUTE_TO_CLEAR=XT \
        ATTRIBUTE_TO_CLEAR=XN \
        ATTRIBUTE_TO_CLEAR=AS \
        ATTRIBUTE_TO_CLEAR=OC \
        ATTRIBUTE_TO_CLEAR=OP \
        SORT_ORDER=queryname \
        RESTORE_ORIGINAL_QUALITIES=true \
        REMOVE_DUPLICATE_INFORMATION=true \
        REMOVE_ALIGNMENT_INFORMATION=true
    """

    
rule align_with_star_2pass_10x:
    input:
        star_genome_output,
        genome_dir=STAR_GENOME_DIR,
        fq='data_raw/10x/{run}/fastq/merged.fq.gz'
    output:
        'data_raw/10x/{run}/star/star.2pass.Aligned.sortedByCoord.out.bam'
    params: 
        prefix='data_raw/10x/{run}/star/star.2pass.'
    threads: 4
    shell:
        '{star_cmd} --genomeDir {input.genome_dir} '
        '--readFilesIn {input.fq} '
        '--outFileNamePrefix {params.prefix} '
        '--outSAMtype BAM SortedByCoordinate ' 
        '--alignSJoverhangMin 8 ' 
        '--alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory '
        '--alignIntronMin 20 --alignIntronMax 1000000 '
        '--alignMatesGapMax 1000000 --sjdbScore 2 '
        '--outFilterType BySJout '
        '--outFilterMultimapNmax 20 --outFilterMismatchNmax 999 '
        '--outFilterMismatchNoverLmax 0.04 '
        '--outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 '
        '--outSAMstrandField intronMotif '
        '--outFilterIntronMotifs RemoveNoncanonical '
        '--outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within '
        '--runThreadN {threads} --twopassMode Basic '
        '--readFilesCommand zcat ' 
        '--limitBAMsortRAM  75000000000 '

# @PG     ID:STAR PN:STAR VN:STAR_2.5.1b  CL:STAR   --runThreadN 4   
# --genomeDir /nfs/gs02/repository/transcriptomes/Homo_sapiens/ensembl_75_transcriptome/1000Genomes_hs37d5/10X/star   
# --readFilesIn /nfs/gs02/IL_seq_data/cellranger/analysis/cellranger201_count_22950_8_1000Genomes_hs37d5-ensembl_75_transcriptome/SC_RNA_COUNTER_CS/SC_RNA_COUNTER/EXTRACT_READS/fork0/chnk0/files/reads.fastq/1.fastq      
# --readNameSeparator space      --outStd SAM   --outSAMtype SAM      --outSAMunmapped Within      --outSAMorder PairedKeepInputOrder   
# --outSAMattrRGline ID:cellranger201_count_22950_8_1000Genomes_hs37d5-ensembl_75_transcriptome:MissingLibrary:1:HJ7JWBBXX:8   
# SM:cellranger201_count_22950_8_1000Genomes_hs37d5-ensembl_75_transcriptome   LB:MissingLibrary.1   PU:cellranger201_count_22950_8_1000Genomes_hs37d5-ensembl_75_transcriptome:MissingLibrary:1:HJ7JWBBXX:8   PL:ILLUMINA      --outSAMmultNmax 18446744073709551615


rule bam2fastq_10x:
    input:
        bam='data_raw/10x/{run}/possorted_genome_bam.bam',
        bai='data_raw/10x/{run}/possorted_genome_bam.bam.bai'
    output:
        fq='data_raw/10x/{run}/fastq/merged.fq.gz'
    shell:
        """
        samtools view -u {input} | \
        samtools collate -uOn 128 - tmp-prefix | \
        samtools fastq -F 0xB00 -t -0 {output.fq} --barcode-tag CB -T CB -
        """
        # 'samtools fastq -t -0 {output.fq} --barcode-tag CB -T CB {input.bam}'


#################


### Smartseq2 rules ###

rule run_to_expt_scesets:
    input:
        files=lambda wildcards: expt_sceset_reqd_files[wildcards.expt]
    output:
        rds='data_processed/{expt}/sceset_{expt}_salmon_allmeta_allcells.rds',
        html='reports/{expt}/{expt}.salmon.qc.html'
    params:
        rmd='reports/{expt}/{expt}.salmon.qc.Rmd'
    shell:
        """
        {Rscript_cmd} -e 'rmarkdown::render("{params.rmd}", clean = TRUE, output_format = "html_document")'
        """

rule multiqc_report:
    input:
        'data_processed/donor_id/donor_id_all.{run}.csv'
    output:
        "reports/multiqc/{run}/multiqc_report.{run}.html"
    shell:
        'multiqc --force --filename {output} ' 
        'data_raw/scrnaseq/{run} reports/fastqc/scrnaseq/{run}'


rule rough_qc:
    input:
        'data_processed/scrnaseq/sceset.{run}.{quant_tool}.preqc_gene.rds'
    output:
        'reports/first_qc/{run}.{quant_tool}.first_qc.html'
    shell:
        '{Rscript_cmd} src/R/compile_report.R -i {input} -o {output} '
        '--template src/Rmd/rough_qc_template.Rmd '


rule kallisto_to_sceset:
    input:
        files=lambda wildcards: kallisto_results_dict[wildcards.run]
    output:
        'data_processed/scrnaseq/sceset.{run}.kallisto.preqc_tx.rds',
        'data_processed/scrnaseq/sceset.{run}.kallisto.preqc_gene.rds',
        'data_processed/scrnaseq/sceset.{run}.kallisto.preqc.feather'
    params:
        input_dir='data_raw/scrnaseq/{run}/quant_kallisto',
        output_prefix='data_processed/scrnaseq/sceset.{run}.kallisto.preqc'
    shell:
        '{Rscript_cmd} {read_kallisto_to_scesets_cmd} '
        '--input_dir {params.input_dir} '
        '--output_prefix {params.output_prefix} '
        '--biomart feb2014.archive.ensembl.org'


rule salmon_to_sceset:
    input:
        files=lambda wildcards: salmon_results_dict[wildcards.run]
    output:
        'data_processed/scrnaseq/sceset.{run}.salmon.preqc_tx.rds',
        'data_processed/scrnaseq/sceset.{run}.salmon.preqc_gene.rds',
        'data_processed/scrnaseq/sceset.{run}.salmon.preqc.feather'
    params:
        input_dir='data_raw/scrnaseq/{run}/quant_salmon',
        output_prefix='data_processed/scrnaseq/sceset.{run}.salmon.preqc'
    shell:
        '{Rscript_cmd} {read_salmon_to_scesets_cmd} '
        '--input_dir {params.input_dir} '
        '--output_prefix {params.output_prefix} '
        '--biomart feb2014.archive.ensembl.org'


rule identify_donor_collect:
    input:
        files=lambda wildcards: variant_donor_id_files_dict[wildcards.run]
    output:
        'data_processed/donor_id/donor_id_all.{run}.csv'
    params:
        dir='data_raw/scrnaseq/{run}/donor_id/',
    run:
        import pandas as pd
        import glob
        import os
        input_files = glob.glob(os.path.join(params.dir, '*.donor_id.csv'))
        df_list = {}
        for infile in input_files:
            df_tmp = pd.read_csv(infile)
            df_list[infile] = df_tmp
        df_out = pd.concat(df_list)
        df_out.to_csv(output[0], index=False)
 

rule identify_donor_runs:
    input:
        sc_vcf='data_raw/scrnaseq/{run}/vcf/{sample}.filtered.vcf.gz',
        hipsci_vcf='data_raw/scrnaseq/{run}/vcf/{sample}.filtered.hipsci.overlap.vcf.gz',
        hipsci_vcf_idx='data_raw/scrnaseq/{run}/vcf/{sample}.filtered.hipsci.overlap.vcf.gz.csi'
    output:
        'data_raw/scrnaseq/{run}/donor_id/{sample}.donor_id.csv'
    params:
        prefix='data_raw/scrnaseq/{run}/donor_id/{sample}.donor_id',
        lines=lambda wildcards: LINES_DICT[wildcards.run][0]
    shell:
        '{Rscript_cmd} src/R/identify_donor_small_vcf.R --input_file "{input.sc_vcf}" '
        '--donor_lines "{params.lines}" '
        '--donor_vcf {input.hipsci_vcf} '
        '--output_prefix "{params.prefix}" '


rule index_hipsci_overlap_vcf:
    input:
        'data_raw/scrnaseq/{run}/vcf/{sample}.filtered.hipsci.overlap.vcf.gz'
    output:
        temp('data_raw/scrnaseq/{run}/vcf/{sample}.filtered.hipsci.overlap.vcf.gz.csi')
    shell:
        'bcftools index {input}'


rule filter_hipsci_overlap_variants:
    input:
        sc_vcf='data_raw/scrnaseq/{run}/vcf/{sample}.filtered.vcf.gz',
        sc_vcf_idx='data_raw/scrnaseq/{run}/vcf/{sample}.filtered.vcf.gz.csi',
        variant_list='data_raw/scrnaseq/{run}/vcf/{sample}.variant_list.txt',
        hipsci_vcf=HIPSCI_VCF
    output:
        tmp=temp('data_raw/scrnaseq/{run}/vcf/{sample}.tmp.vcf.gz'),
        vcf=temp('data_raw/scrnaseq/{run}/vcf/{sample}.filtered.hipsci.overlap.vcf.gz')
    shell:
        """
        bcftools view -o {output.tmp} -O z -l 9 -R {input.variant_list} {input.hipsci_vcf}
        vcf-sort {output.tmp} | bgzip -c > {output.vcf}
        """


rule singlecell_variant_list:
    input:
        'data_raw/scrnaseq/{run}/vcf/{sample}.filtered.vcf.gz'
    output:
        temp('data_raw/scrnaseq/{run}/vcf/{sample}.variant_list.txt')
    shell:
        """
        set +euo pipefail
        echo -e "1\\t1\\tA\\tC" > {output}
        bcftools view -O v {input} | grep -v ^# | awk \'{{sub(/chr/,""); print $1"\\t"$2"\\t"$4"\\t"$5}}\' >> {output}
        set -euo pipefail
        """


rule index_bgzip_vcf:
    input:
        'data_raw/scrnaseq/{run}/vcf/{sample}.filtered.vcf.gz'
    output:
        'data_raw/scrnaseq/{run}/vcf/{sample}.filtered.vcf.gz.csi'
    shell:
        'bcftools index {input}'


rule bgzip_vcf:
    input:
        'data_raw/scrnaseq/{run}/vcf/{sample}.filtered.vcf'
    output:
        'data_raw/scrnaseq/{run}/vcf/{sample}.filtered.vcf.gz'
    shell:
        'vcf-sort {input} | bgzip -c > {output}'


rule filter_variants_gatk:
    input:
        vcf='data_raw/scrnaseq/{run}/vcf/{sample}.unfiltered.vcf',
        fasta=fasta_unzipped,
        fai=fasta_idx,
        fa_dict=fasta_dict
    output:
        temp('data_raw/scrnaseq/{run}/vcf/{sample}.filtered.vcf')
    shell:
        '{gatk_cmd} -T VariantFiltration -R {input.fasta} -V {input.vcf} '
        '-window 35 -cluster 3 -filterName FS -filter "FS > 30.0" '
        '-filterName QD -filter "QD < 2.0" -o {output}'


rule call_variants_gatk:
    input:
        fasta=fasta_unzipped,
        fai=fasta_idx,
        fa_dict=fasta_dict,
        dbSnp=dbSnpVcf,
        dbSnpSmall= dbSnpVcfSmall,
        bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr.bam'
    output:
        temp('data_raw/scrnaseq/{run}/vcf/{sample}.unfiltered.vcf')
    shell:
        '{gatk_cmd} -T HaplotypeCaller -R {input.fasta} -I {input.bam} '
        '-dontUseSoftClippedBases '
        '-D {input.dbSnp} -gt_mode GENOTYPE_GIVEN_ALLELES '
        '-alleles {input.dbSnpSmall} -L {input.dbSnpSmall} -o {output}' 


# rule call_ase_gatk_low_thresh:
#     input:
#         fasta=fasta_unzipped,
#         fai=fasta_idx,
#         fa_dict=fasta_dict,
#         vcf='/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2014-11_SS/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.20151005.genotypes.gdid.mac1.recode.WGS.maf0.1.with_chr.vcf',
#         bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.bam'
#     output:
#         'data_raw/scrnaseq/{run}/ase/low_thresh/{sample}.ase.lowthresh.tsv'
#     shell:
#         '{gatk_cmd} -T ASEReadCounter -R {input.fasta} -I {input.bam} -o {output} '
#         '-sites {input.vcf} -U ALLOW_N_CIGAR_READS -minDepth 5 '
#         '--minMappingQuality 5 --minBaseQuality 0 -drf DuplicateRead '


# rule call_ase_gatk_high_thresh:
#     input:
#         fasta=fasta_unzipped,
#         fai=fasta_idx,
#         fa_dict=fasta_dict,
#         vcf='/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2014-11_SS/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.20151005.genotypes.gdid.mac1.recode.WGS.maf0.1.with_chr.vcf',
#         bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr.bam'
#     output:
#         'data_raw/scrnaseq/{run}/ase/high_thresh/{sample}.ase.highthresh.tsv'
#     shell:
#         '{gatk_cmd} -T ASEReadCounter -R {input.fasta} -I {input.bam} -o {output} '
#         '-sites {input.vcf} -U ALLOW_N_CIGAR_READS -minDepth 20 '
#         '--minMappingQuality 10 --minBaseQuality 2'


rule fastqc_reports:
    input:
        'data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr.bam'
    output:
        'reports/fastqc/scrnaseq/{run}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr_fastqc.html'
    params:
        output_dir="reports/fastqc/scrnaseq/{run}/"
    shell:
        '/nfs/software/stegle/FastQC/fastqc -o {params.output_dir} {input}'


rule count_genes_unique:
    input:
        bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr.bam',
        anno='/hps/nobackup/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.75_ERCC.gtf'
    output:
        'data_raw/scrnaseq/{run}/featureCounts/unique_counts/{sample}.gene.counts.unique.tsv'
    priority: 4
    shell:
        '/hps/nobackup/stegle/users/mjbonder/tools/subread-1.5.2-Linux-x86_64/bin/featureCounts --ignoreDup -B -C -p --primary '
        '-a {input.anno} -g gene_id '
        '-o  {output} {input.bam}'


rule count_genes:
    input:
        bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr.bam',
        anno='/hps/nobackup/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.75_ERCC.gtf'
    output:
        'data_raw/scrnaseq/{run}/featureCounts/total_counts/{sample}.gene.counts.tsv'
    priority: 4
    shell:
        '/hps/nobackup/stegle/users/mjbonder/tools/subread-1.5.2-Linux-x86_64/bin/featureCounts -B -C -p --primary '
        '-a {input.anno} -g gene_id '
        '-o  {output} {input.bam}'


rule picard_qc_metrics:
    input:
        bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr.bam',
        ref_flat='/hps/nobackup/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.75_ERCC.ref_flat.txt',
        ribo_intervals=rRNAIntervals
    output:
        metrics='data_raw/scrnaseq/{run}/picard_metrics/{sample}.picard.rna.metrics',
        #plot_pdf='figures/picard_cov_vs_pos/{run}/{sample}.picard.cov.pos.pdf'
    shell:
        '{picard_cmd} CollectRnaSeqMetrics I={input.bam} O={output.metrics} '
        'REF_FLAT={input.ref_flat} STRAND_SPECIFICITY=NONE '
        'RIBOSOMAL_INTERVALS={input.ribo_intervals} '
        #'CHART_OUTPUT={output.plot_pdf}'


rule recalibrated_writer_gatk:
    input:
        fasta=fasta_unzipped,
        bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bam',
        bqsr='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr'
    output:
        'data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr.bam'
    shell:
        '{gatk_cmd} -T PrintReads -R {input.fasta} -I {input.bam} '
        '-BQSR {input.bqsr} -nct 2 '
        '-o {output}' 


rule base_recalibrator_gatk:
    input:
        fasta=fasta_unzipped,
        dbSnp=dbSnpVcf,
        bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bam',
        known1=knownIndelsMills,
        known2=knownIndels100G
    output:
        temp('data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr')
    shell:
        '{gatk_cmd} -T BaseRecalibrator -R {input.fasta} -I {input.bam} '
        '-knownSites {input.known1} -knownSites {input.known2} -knownSites {input.dbSnp} '
        '-nct 2 '
        '-o {output}' 


rule indel_realignment_gatk:
    input:
        fasta=fasta_unzipped,
        bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.bam',
        targetIntervals=reAlignmentIntervals,
        known1=knownIndelsMills,
        known2=knownIndels100G
    output:
        temp('data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bam')
    shell:
        '{gatk_cmd} -T IndelRealigner -R {input.fasta} -I {input.bam} '
        '-targetIntervals {input.targetIntervals} -known {input.known1} -known {input.known2} '
        '-U ALLOW_N_CIGAR_READS --consensusDeterminationModel KNOWNS_ONLY --LODThresholdForCleaning 0.4  '
        '-o {output}' 


rule split_n_trim_gatk:
    input:
        fasta=fasta_unzipped,
        fai=fasta_idx,
        fa_dict=fasta_dict,
        bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.dedup.bam'
    output:
        temp('data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.bam')
    shell:
        '{gatk_cmd} -T SplitNCigarReads -R {input.fasta} -I {input.bam} -o {output} '
        '-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 '
        '-U ALLOW_N_CIGAR_READS '


rule picard_mark_dups:
    input:
        'data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.rgadded.bam'
    output:
        bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.dedup.bam',
        metrics='data_raw/scrnaseq/{run}/star/{sample}/{sample}.output.metrics'
    shell:
        '{picard_cmd} MarkDuplicates I={input} O={output.bam} CREATE_INDEX=true '
        'VALIDATION_STRINGENCY=SILENT M={output.metrics} '


rule picard_read_groups:
    input:
        bam='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.out.bam',
        cram='data_raw/scrnaseq/{run}/cram/{sample}.cram'
    output:
        temp('data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.rgadded.bam')
    shell:
        '{picard_cmd} AddOrReplaceReadGroups I={input.bam} O={output} SO=coordinate '
        'RGID={wildcards.sample} RGLB={wildcards.sample} '
        'RGPL=ILLUMINA RGPU=MACHINE1 RGSM={wildcards.sample}'


rule index_star_bams:
    input:
        'data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.out.bam'
    output:
        'data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.out.bam.bai'
    shell:
        'samtools index {input} '


rule align_with_star_2pass:
    input:
        star_genome_output,
        genome_dir=STAR_GENOME_DIR,
        fq1='data_raw/scrnaseq/{run}/fastq/{sample}_1_val_1.fq.gz',
        fq2='data_raw/scrnaseq/{run}/fastq/{sample}_2_val_2.fq.gz'
    output:
        'data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.out.bam'
    params: 
        prefix='data_raw/scrnaseq/{run}/star/{sample}/{sample}.2pass.'
    threads: 4
    shell:
        '{star_cmd} --genomeDir {input.genome_dir} '
        '--readFilesIn {input.fq1} {input.fq2} '
        '--outFileNamePrefix {params.prefix} '
        '--outSAMtype BAM SortedByCoordinate ' 
        '--alignSJoverhangMin 8 ' 
        '--alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory '
        '--alignIntronMin 20 --alignIntronMax 1000000 '
        '--alignMatesGapMax 1000000 --sjdbScore 2 '
        '--outFilterType BySJout '
        '--outFilterMultimapNmax 20 --outFilterMismatchNmax 999 '
        '--outFilterMismatchNoverLmax 0.04 '
        '--outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 '
        '--outSAMstrandField intronMotif '
        '--outFilterIntronMotifs RemoveNoncanonical '
        '--outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within '
        '--runThreadN {threads} --twopassMode Basic '
        '--readFilesCommand zcat ' 


rule build_star_genome_indexes:
    input:
        fasta=fasta_unzipped,
        annotation='/hps/nobackup/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.75_ERCC.gtf'
    output:
        star_genome_output
    threads: 8
    shell:
        '{star_cmd} --runMode genomeGenerate --genomeDir {STAR_GENOME_DIR} '
        '--genomeFastaFiles {input.fasta} --runThreadN {threads} '
        '--sjdbGTFfile {input.annotation} --sjdbOverhang 100'


rule kallisto_quant_GRCh38:
    input:
        kidx=kallisto_idx_GRCh38,
        fq1='data_raw/scrnaseq/{run}/fastq/{sample}_1_val_1.fq.gz',
        fq2='data_raw/scrnaseq/{run}/fastq/{sample}_2_val_2.fq.gz'
    output:
        'data_raw/scrnaseq/{run}/quant_kallisto_GRCh38/{sample}/abundance.tsv'
    params:
        folder='data_raw/scrnaseq/{run}/quant_kallisto_GRCh38/{sample}/'
    shell:
        '{kallisto_cmd} quant -i {input.kidx} '
        '-o {params.folder} '
        '--bias {input.fq1} {input.fq2}'


rule kallisto_quant:
    input:
        kidx=kallisto_idx,
        fq1='data_raw/scrnaseq/{run}/fastq/{sample}_1_val_1.fq.gz',
        fq2='data_raw/scrnaseq/{run}/fastq/{sample}_2_val_2.fq.gz'
    output:
        'data_raw/scrnaseq/{run}/quant_kallisto/{sample}/abundance.tsv'
    params:
        folder='data_raw/scrnaseq/{run}/quant_kallisto/{sample}/'
    shell:
        '{kallisto_cmd} quant -i {input.kidx} '
        '-o {params.folder} '
        '--bias {input.fq1} {input.fq2}'


rule salmon_quant_GRCh38:
    input:
        sidx=salmon_idx_GRCh38,
        fq1='data_raw/scrnaseq/{run}/fastq/{sample}_1_val_1.fq.gz',
        fq2='data_raw/scrnaseq/{run}/fastq/{sample}_2_val_2.fq.gz'
    output:
        'data_raw/scrnaseq/{run}/quant_salmon_GRCh38/{sample}/quant.sf'
    threads: 8
    params:
        folder='data_raw/scrnaseq/{run}/quant_salmon_GRCh38/{sample}/'
    shell:
        '{salmon_cmd} quant -i {input.sidx} -l IU '
        '-1 {input.fq1} -2 {input.fq2} '
        '--seqBias --gcBias --threads {threads} --useVBOpt ' 
        '-o {params.folder}'


rule salmon_quant:
    input:
        sidx=salmon_idx,
        fq1='data_raw/scrnaseq/{run}/fastq/{sample}_1_val_1.fq.gz',
        fq2='data_raw/scrnaseq/{run}/fastq/{sample}_2_val_2.fq.gz'
    output:
        'data_raw/scrnaseq/{run}/quant_salmon/{sample}/quant.sf'
    threads: 4
    params:
        folder='data_raw/scrnaseq/{run}/quant_salmon/{sample}/'
    shell:
        '{salmon_cmd} quant -i {input.sidx} -l IU '
        '-1 {input.fq1} -2 {input.fq2} '
        '--seqBias --gcBias --threads {threads} --useVBOpt ' 
        '-o {params.folder}'


rule trim_fastq:
    input:
        fq1="data_raw/scrnaseq/{run}/fastq/{sample}_1.fastq",
        fq2="data_raw/scrnaseq/{run}/fastq/{sample}_2.fastq"
    output:
        fq1='data_raw/scrnaseq/{run}/fastq/{sample}_1_val_1.fq.gz',
        fq2='data_raw/scrnaseq/{run}/fastq/{sample}_2_val_2.fq.gz'
    shell:
        '/nfs/software/stegle/trim_galore/trim_galore --gzip --nextera '
        '--output_dir data_raw/scrnaseq/{wildcards.run}/fastq '
        '--length 40 --paired '
        '{input.fq1} {input.fq2}'


rule cram2fastq:
    input:
        'data_raw/scrnaseq/{run}/cram/{sample}.cram'
    output:
        fq1=temp('data_raw/scrnaseq/{run}/fastq/{sample}_1.fastq'),
        fq2=temp('data_raw/scrnaseq/{run}/fastq/{sample}_2.fastq')
    shell:
        """
        samtools view -u {input} | \
        samtools collate -uOn 128 - tmp/tmp-prefix-{wildcards.sample} | \
        samtools fastq -F 0xB00 -1 {output.fq1} -2 {output.fq2} -
        """


rule seq_metadata_collect:
    input:
        files=lambda wildcards: seq_metadata_files_dict[wildcards.run]
    output:
        'metadata/scrnaseq/seq_metadata.{run}.tsv'
    params:
        dir='data_raw/scrnaseq/{run}/meta/',
    run:
        import os
        import glob
        import re
        import pandas as pd
        metafiles = glob.glob(os.path.join(params.dir, '*.meta'))
        seq_sample_names = [os.path.basename(x).replace(".meta", "") for x in metafiles]
        metadict = {}
        for fname in metafiles:
            tmpdict = {}
            with open(fname) as f:
                current_attribute = ''
                for line in f:
                    if re.search('^attribute', line):
                        current_attribute = line.replace('attribute: ', '').strip()
                    if re.search('^value', line):
                        tmpdict[current_attribute] = line.replace('value: ', '').strip()
                metadict[os.path.basename(fname).replace(".meta", "")] = tmpdict
        df = pd.DataFrame(metadict).transpose()
        df['cram_id'] = df.index
        df.to_csv(output[0], index=False, sep='\t')


rule build_kallisto_index_GRCh38:
    input:
        fasta_GRCh38
    output:
        kallisto_idx_GRCh38
    shell:
        '{kallisto_cmd} index -i {output} {input}'


rule build_kallisto_index:
    input:
        fasta
    output:
        kallisto_idx
    shell:
        '{kallisto_cmd} index -i {output} {input}'


rule build_salmon_index_GRCh38:
    input:
        fasta_GRCh38
    output:
        salmon_idx_GRCh38
    shell:
        '{salmon_cmd} index -t {input} -i {output} --type quasi -k 31 '
        '--perfectHash'


rule build_salmon_index:
    input:
        fasta
    output:
        salmon_idx
    shell:
        '{salmon_cmd} index -t {input} -i {output} --type quasi -k 31 '
        '--perfectHash'


rule create_fasta_index_GRCh38:
    input:
        fasta_unzipped_GRCh38
    output:
        fasta_idx_GRCh38
    shell:
        'samtools faidx {input}'


rule create_fasta_index:
    input:
        fasta_unzipped
    output:
        fasta_idx
    shell:
        'samtools faidx {input}'


rule create_fasta_dict_GRCh38:
    input:
        fasta_unzipped_GRCh38
    output:
        fasta_dict_GRCh38
    shell:
        '{picard_cmd} CreateSequenceDictionary R={input} O={output}'


rule create_fasta_dict:
    input:
        fasta_unzipped
    output:
        fasta_dict
    shell:
        '{picard_cmd} CreateSequenceDictionary R={input} O={output}'


rule process_genotypes_eqtl:
    input:
        HIPSCI_VCF
    output:
        'data_processed/eqtl/REL-2016-09.INFO_0.4_filtered.20160912.genotypes.with_monogenic_diabetes_lines.{chrom}.h5'
    run:
        import os
        import vcf
        import h5py
        import numpy as np
        import pandas as pd
        from datetime import datetime
        from cyvcf2 import VCF
        ## define files and samples
        vcf_file = input[0]
        h5_file = output[0]
        vcf_reader = vcf.Reader(filename=vcf_file, compressed=True)
        samples = np.array([a.encode('utf8') for a in vcf_reader.samples])
        chrom = wildcards.chrom.replace('chr', '')
        ## read in data from vcf
        refs = {}
        alts = {}
        chroms = {}
        poss = {}
        starts = {}
        ends = {}
        ids = {}
        gt_types = {}
        vcf = VCF(vcf_file)
        for variant in vcf(chrom + ':0-2000000000'):
            if variant.is_snp and variant.ID is not None:
                # numpy arrays of specific things we pull from the sample fields.
                # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
                gts = np.array(variant.gt_types)
                # To make things more understandable change the UNKNOWN value from 2 to nan
                # and the HOM_ALT value from 3 to 2
                gts = gts.astype(float)
                gts[gts == 2] = np.nan
                gts[gts == 3] = 2.0
                if (np.isnan(gts).mean() < 0.05): # retain variants with <5% missing genotypes
                    var_id = str(variant.ID, "utf-8")
                    refs[var_id] = variant.REF
                    alts[var_id] = variant.ALT[0] # e.g. REF='A', ALT=['C', 'T']
                    chroms[var_id] = variant.CHROM
                    poss[var_id] = variant.POS
                    starts[var_id] = variant.start
                    ends[var_id] = variant.end
                    ids[var_id] = var_id
                    gt_types[var_id] = gts
        ## put data into dataframes for convenience
        snp_info_df = pd.DataFrame({'ID': ids, 'CHROM': chroms, 'POS': poss, 'START': starts, 'END': ends, 'REF': refs, 'ALT': alts})
        gt_df = pd.DataFrame(gt_types)
        # # Write to HDF5 file
        fout = h5py.File(h5_file)
        fout.attrs['time_stamp'] = datetime.now().isoformat()
        fout.attrs['author'] = 'Davis McCarthy'
        fout.attrs['institution'] = 'EMBL-EBI'
        # Add sampleIDs to HDF5 file.
        fout.create_dataset('sampleID', data=samples)
        fout.create_group('sample_info')
        fout.create_dataset('sample_info/sampleID', data=samples)
        # Add genotype matrix.
        fout.create_group('genotype')
        fout.create_dataset('/genotype/matrix', data=gt_df.values)
        # Add SNP Info to HDF5 file.
        fout.create_group(name='snp_info')
        fout.create_dataset('/snp_info/alt', data=np.array([a.encode('utf-8') for a in snp_info_df['ALT'].values]))
        fout.create_dataset('/snp_info/ref', data=np.array([a.encode('utf-8') for a in snp_info_df['REF'].values]))
        fout.create_dataset('/snp_info/chrom', data=np.array([a.encode('utf-8') for a in snp_info_df['CHROM'].values]))
        fout.create_dataset('/snp_info/pos', data=snp_info_df['POS'].values)
        fout.create_dataset('/snp_info/start', data=snp_info_df['START'].values)
        fout.create_dataset('/snp_info/end', data=snp_info_df['END'].values)
        fout.create_dataset('/snp_info/gdid', data=np.array([a.encode('utf-8') for a in snp_info_df['ID'].values]))
        fout.close()


rule compute_Kpop_h5:
    input: 
        'data_raw/grm/hipsci_grm.rel.gz'
    output:
        'data_processed/eqtl/REL-2016-09.Kpop_maf05.h5'
    shell:
        '{Rscript_cmd} '
        'src/eqtl/grm_to_Kpop_h5.R -i {input} -o {output}'


rule compute_grm_hipsci:
    input: 
        HIPSCI_VCF
    output:
        'data_raw/grm/hipsci_grm.rel.gz'
    shell:
        '/nfs/software/stegle/plink --vcf {HIPSCI_VCF} --maf 0.05 --double-id '
        '--make-rel square gz --out data_raw/grm/hipsci_grm' 


rule create_hipsci_vcf:
    input:
        'Data/REL-2017-39_genotype_vcf_info_0.4_file_list.txt'
    output:
        vcf=HIPSCI_VCF,
        vcf_idx=HIPSCI_VCF + '.tbi'
    shell:
        """
        bcftools concat -f {input} -O z > {output.vcf}
        tabix -f -p vcf {output.vcf}
        """

rule make_refflat_gene_annotation:
    input:
        '/hps/nobackup/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.75_ERCC.gtf'
    output:
        keep='/hps/nobackup/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.75_ERCC.ref_flat.txt',
        tmp=temp('/hps/nobackup/stegle/datasets/references/human/GRCh37/tmp.ref_flat.txt')
    shell:
        """
        /nfs/software/stegle/users/davis/gtfToGenePred -genePredExt -geneNameAsName2 {input} {output.tmp}
        paste <(cut -f 12 {output.tmp}) <(cut -f 1-10 {output.tmp}) > {output.keep} 
        """


rule generate_new_tabix_idx:
    input:
        hipsci_chr_vcf
    output:
        hipsci_chr_vcf_idx
    run:
        for chrom in hipsci_chr_vcf:
            print(chrom)
            cmd = 'tabix -f -p vcf ' + chrom
            print(cmd)
            subprocess.run(cmd, shell = True)


rule create_vcf_file_list:
    output:
        'Data/REL-2017-39_genotype_vcf_info_0.4_file_list.txt'
    run:
        with open(output[0], "w+") as outfile:
                for i in range(1, 23):
                        outfile.write('/hps/nobackup/hipsci/scratch/genotypes/imputed/2017-03-27/INFO_0.4/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.INFO_0.4_filtered.20170327.genotypes.chr%s.vcf.gz\n' %str(i))
