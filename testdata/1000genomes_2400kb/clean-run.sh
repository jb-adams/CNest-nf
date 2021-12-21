# ##################################################
# CLEAN NEXTFLOW WORKING DIRECTORIES
# ##################################################

rm -rf .nextflow* 
rm -rf work

# ##################################################
# DOWNLOAD GRCh38 (IF NOT PRESENT)
# ##################################################

REF_FILE_PATH="./ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
REF_FILE_FTP_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
if [ -f ${REF_FILE_PATH} ]; then
    echo "GRCh38 reference genome already downloaded";
else
    echo "GRCh38 reference genome not downloaded...downloading";
    wget -P ./ref ${REF_FILE_FTP_URL}
fi;

REF_INDEX_PATH="${REF_FILE_PATH}.fai"
REF_INDEX_FTP_URL="${REF_FILE_FTP_URL}.fai"
if [ -f ${REF_INDEX_PATH} ]; then
    echo "GRCh38 index already downloaded";
else
    echo "GRCh38 index not downloaded...downloading";
    wget -P ./ref ${REF_INDEX_FTP_URL}
fi;

# ##################################################
# RUN CNEST ON DOWNSAMPLED 1000 GENOMES DATA
# ##################################################
nextflow run ../../main.nf \
    --alnformat cram \
    --ref ./ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --bed ./bed/2400kb.bed \
    --cor 0 \
    --skipem \
    --name_0 HG00443 --aln_0 ./cram/HG00443.final.2400kb.cram --idx_0 ./cram/HG00443.final.2400kb.cram.crai \
    --name_1 HG00445 --aln_1 ./cram/HG00445.final.2400kb.cram --idx_1 ./cram/HG00445.final.2400kb.cram.crai \
    --name_2 HG00446 --aln_2 ./cram/HG00446.final.2400kb.cram --idx_2 ./cram/HG00446.final.2400kb.cram.crai \
    --name_3 HG00448 --aln_3 ./cram/HG00448.final.2400kb.cram --idx_3 ./cram/HG00448.final.2400kb.cram.crai \
    --name_4 HG00449 --aln_4 ./cram/HG00449.final.2400kb.cram --idx_4 ./cram/HG00449.final.2400kb.cram.crai \
    --name_5 HG00451 --aln_5 ./cram/HG00451.final.2400kb.cram --idx_5 ./cram/HG00451.final.2400kb.cram.crai \
    --name_6 NA12718 --aln_6 ./cram/NA12718.final.2400kb.cram --idx_6 ./cram/NA12718.final.2400kb.cram.crai \
    --name_7 NA12748 --aln_7 ./cram/NA12748.final.2400kb.cram --idx_7 ./cram/NA12748.final.2400kb.cram.crai \
    --name_8 NA12775 --aln_8 ./cram/NA12775.final.2400kb.cram --idx_8 ./cram/NA12775.final.2400kb.cram.crai \
    --name_9 NA12778 --aln_9 ./cram/NA12778.final.2400kb.cram --idx_9 ./cram/NA12778.final.2400kb.cram.crai
