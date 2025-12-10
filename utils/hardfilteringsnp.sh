
#Execute
module purge
module load gatk/4.3.0.0-gcccore-11.3.0-java-11
module load samtools/1.18-gcc-12.3.0


gatk VariantFiltration \
    -V Chr00_17.vcf.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O snps_filterd_Chr00_17.vcf
