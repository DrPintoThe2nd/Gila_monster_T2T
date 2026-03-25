#!bin/bash

cp SRR2889291_Gallus-female.quant/quant.genes.sf SRR2889291_Gallus-female.quants.genes.tsv
cp SRR2889292_Gallus-female.quant/quant.genes.sf SRR2889292_Gallus-female.quants.genes.tsv
cp SRR2889293_Gallus-female.quant/quant.genes.sf SRR2889293_Gallus-female.quants.genes.tsv
cp SRR2889295_Gallus-male.quant/quant.genes.sf SRR2889295_Gallus-male.quants.genes.tsv
cp SRR2889296_Gallus-male.quant/quant.genes.sf SRR2889296_Gallus-male.quants.genes.tsv
cp SRR2889297_Gallus-male.quant/quant.genes.sf SRR2889297_Gallus-male.quants.genes.tsv
cp SRR6344890_Gila_35-F.quant/quant.genes.sf SRR6344890_Gila_35-F.quants.genes.tsv
cp SRR6344891_Gila_L-F.quant/quant.genes.sf  SRR6344891_Gila_L-F.quants.genes.tsv
cp SRR6344892_Gila_10-M.quant/quant.genes.sf SRR6344892_Gila_10-M.quants.genes.tsv
cp SRR6344893_Gila_16-M.quant/quant.genes.sf SRR6344893_Gila_16-M.quants.genes.tsv
cp SRR6344900_Gila_KI01-M.quant/quant.genes.sf SRR6344900_Gila_KI01-M.quants.genes.tsv
cp SRR6344901_Gila_30-F.quant/quant.genes.sf SRR6344901_Gila_30-F.quants.genes.tsv

###################################################getting Gila monster data ready for R
#extract and clean gff records
awk ' $3 == "gene" ' Heloderma_suspectum.final.reference.PAR-masked.fixed.gff3 | sed 's/ID=//' | sed 's/\;/\t/g' | sed 's/Note=//' | sed 's/gene_biotype=//' | sed 's/ /_/g' | sed 's/Similar_to_//' | cut -f1-12 > Gila_annotations.txt

#combine salmon quants w/ associated gff records
ls *Gila_*.quants.genes.tsv > samples.txt
sed -i 's/.quants.genes.tsv//' samples.txt
cat samples.txt | while read line
do
echo "Another one..."
Rscript join_expression.R Gila_annotations.txt $line\.quants.genes.tsv $line\.tsv
done

#############################################################################FIX THIS
#combine individual TPM counts with gff records (except effective length)
#males
cut -f1-13,15 rna_10-002_S1_L007.tsv > tmp
sed -i 's/TPM/rna_10/' tmp
cut -f15 rna_16-004_S2_L007.tsv > tmp2
sed -i 's/TPM/rna_16/' tmp2
cut -f15 rna_KI01-012_S6_L007.tsv > tmp3
sed -i 's/TPM/rna_KI01/' tmp3
paste tmp tmp2 tmp3 > tmp4
awk ' $1 != "16-S_mtDNA" ' tmp4 > rna_males.gff.tsv

#females
cut -f1-13,15 rna_30-005_S3_L007.tsv > tmp
sed -i 's/TPM/rna_30/' tmp
cut -f15 rna_35-006_S4_L007.tsv > tmp2
sed -i 's/TPM/rna_35/' tmp2
cut -f15 rna_L-007_S5_L007.tsv > tmp3
sed -i 's/TPM/rna_L007/' tmp3
paste tmp tmp2 tmp3 > tmp4
awk ' $1 != "16-S_mtDNA" ' tmp4 > rna_females.gff.tsv

######################################################getting Chicken data ready for R

sed 's/description=.*;gbkey/gbkey/g' GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff3 > tmp
awk ' $3 == "gene" ' tmp | sed 's/ID=//' | sed 's/\;/\t/g' |  sed 's/gene_biotype=//' | sed 's/ /_/g' | cut -f1-9,11,14 > tmp2
awk ' $3 == "gene" ' tmp | sed 's/ID=//' | sed 's/\;/\t/g' |  sed 's/gene_biotype=//' | sed 's/ /_/g' | cut -f10 > tmp3
paste tmp2 tmp3 > Gallus_annotations.txt

#combine salmon quants w/ associated gff records
ls *Gallus*.quants.genes.tsv > samples.txt
sed -i 's/.quants.genes.tsv//' samples.txt
cat samples.txt | while read line
do
echo "Another one..."
Rscript join_expression.R Gallus_annotations.txt $line\.quants.genes.tsv $line\.tsv
done

