#! bin/bash

#move and rename quant files
cp ../ERR216301_Echis_coloratus_skin.quant/quant.genes.sf                  ERR216301_Echis_coloratus_skin.genes.tsv
cp ../ERR216302_Echis_coloratus_skin.quant/quant.genes.sf                  ERR216302_Echis_coloratus_skin.genes.tsv
cp ../ERR216311_Echis_coloratus_venom_gland.quant/quant.genes.sf           ERR216311_Echis_coloratus_venom_gland.genes.tsv
cp ../ERR216312_Echis_coloratus_venom_gland.quant/quant.genes.sf           ERR216312_Echis_coloratus_venom_gland.genes.tsv
cp ../ERR216319_Echis_coloratus_scent_gland.quant/quant.genes.sf           ERR216319_Echis_coloratus_scent_gland.genes.tsv
cp ../ERR216326_Echis_coloratus_scent_gland.quant/quant.genes.sf           ERR216326_Echis_coloratus_scent_gland.genes.tsv
cp ../ERR216304_Eublepharis_macularius_skin.quant/quant.genes.sf           ERR216304_Eublepharis_macularius_skin.genes.tsv
cp ../ERR216306_Eublepharis_macularius_skin.quant/quant.genes.sf           ERR216306_Eublepharis_macularius_skin.genes.tsv
cp ../ERR216315_Eublepharis_macularius_salivary_gland.quant/quant.genes.sf ERR216315_Eublepharis_macularius_salivary_gland.genes.tsv
cp ../ERR216316_Eublepharis_macularius_salivary_gland.quant/quant.genes.sf ERR216316_Eublepharis_macularius_salivary_gland.genes.tsv
cp ../ERR216322_Eublepharis_macularius_scent_gland.quant/quant.genes.sf    ERR216322_Eublepharis_macularius_scent_gland.genes.tsv
cp ../ERR216325_Eublepharis_macularius_scent_gland.quant/quant.genes.sf    ERR216325_Eublepharis_macularius_scent_gland.genes.tsv
cp ../ERR216303_Opheodrys_aestivus_skin.quant/quant.genes.sf               ERR216303_Opheodrys_aestivus_skin.genes.tsv
cp ../ERR216305_Opheodrys_aestivus_skin.quant/quant.genes.sf               ERR216305_Opheodrys_aestivus_skin.genes.tsv
cp ../ERR216313_Opheodrys_aestivus_salivary_gland.quant/quant.genes.sf     ERR216313_Opheodrys_aestivus_salivary_gland.genes.tsv
cp ../ERR216314_Opheodrys_aestivus_salivary_gland.quant/quant.genes.sf     ERR216314_Opheodrys_aestivus_salivary_gland.genes.tsv
cp ../ERR216320_Opheodrys_aestivus_scent_gland.quant/quant.genes.sf        ERR216320_Opheodrys_aestivus_scent_gland.genes.tsv
cp ../ERR216321_Opheodrys_aestivus_scent_gland.quant/quant.genes.sf        ERR216321_Opheodrys_aestivus_scent_gland.genes.tsv
cp ../ERR216298_Pantherophis_guttatus_skin.quant/quant.genes.sf            ERR216298_Pantherophis_guttatus_skin.genes.tsv
cp ../ERR216307_Pantherophis_guttatus_salivary_gland.quant/quant.genes.sf  ERR216307_Pantherophis_guttatus_salivary_gland.genes.tsv
cp ../ERR216308_Pantherophis_guttatus_salivary_gland.quant/quant.genes.sf  ERR216308_Pantherophis_guttatus_salivary_gland.genes.tsv
cp ../ERR216317_Pantherophis_guttatus_scent_gland.quant/quant.genes.sf     ERR216317_Pantherophis_guttatus_scent_gland.genes.tsv
cp ../ERR216323_Pantherophis_guttatus_scent_gland.quant/quant.genes.sf     ERR216323_Pantherophis_guttatus_scent_gland.genes.tsv
cp ../ERR216299_Python_regius_skin.quant/quant.genes.sf                    ERR216299_Python_regius_skin.genes.tsv
cp ../ERR216300_Python_regius_skin.quant/quant.genes.sf                    ERR216300_Python_regius_skin.genes.tsv
cp ../ERR216309_Python_regius_salivary_gland.quant/quant.genes.sf          ERR216309_Python_regius_salivary_gland.genes.tsv
cp ../ERR216310_Python_regius_salivary_gland.quant/quant.genes.sf          ERR216310_Python_regius_salivary_gland.genes.tsv
cp ../ERR216318_Python_regius_scent_gland.quant/quant.genes.sf             ERR216318_Python_regius_scent_gland.genes.tsv
cp ../ERR216324_Python_regius_scent_gland.quant/quant.genes.sf             ERR216324_Python_regius_scent_gland.genes.tsv

#prepare NCBI gffs for tsv-based merging
sed 's/description=.*;gbkey/gbkey/g' ../Eublepharis_macularius_reference.genome.gff > tmp
awk ' $3 == "gene" ' tmp | sed 's/ID=//' | sed 's/\;/\t/g' |  sed 's/gene_biotype=//' | sed 's/ /_/g' | cut -f1-9,11,14 > tmp2
awk ' $3 == "gene" ' tmp | sed 's/ID=//' | sed 's/\;/\t/g' |  sed 's/gene_biotype=//' | sed 's/ /_/g' | cut -f10 > tmp3
paste tmp2 tmp3 > Eublepharis_macularius_annotations.txt

sed 's/description=.*;gbkey/gbkey/g' ../Pantherophis_guttatus_reference.genome.gff > tmp
awk ' $3 == "gene" ' tmp | sed 's/ID=//' | sed 's/\;/\t/g' |  sed 's/gene_biotype=//' | sed 's/ /_/g' | cut -f1-9,11,14 > tmp2
awk ' $3 == "gene" ' tmp | sed 's/ID=//' | sed 's/\;/\t/g' |  sed 's/gene_biotype=//' | sed 's/ /_/g' | cut -f10 > tmp3
paste tmp2 tmp3 > Pantherophis_guttatus_annotations.txt

sed 's/description=.*;gbkey/gbkey/g' ../Python_molurus_reference.genome.gff > tmp
awk ' $3 == "gene" ' tmp | sed 's/ID=//' | sed 's/\;/\t/g' |  sed 's/gene_biotype=//' | sed 's/ /_/g' | cut -f1-9,11,14 > tmp2
awk ' $3 == "gene" ' tmp | sed 's/ID=//' | sed 's/\;/\t/g' |  sed 's/gene_biotype=//' | sed 's/ /_/g' | cut -f10 > tmp3
paste tmp2 tmp3 > Python_molurus_annotations.txt

sed 's/description=.*;gbkey/gbkey/g' ../Vipera_berus_reference.genome.gff > tmp
awk ' $3 == "gene" ' tmp | sed 's/ID=//' | sed 's/\;/\t/g' |  sed 's/gene_biotype=//' | sed 's/ /_/g' | cut -f1-9,11,14 > tmp2
awk ' $3 == "gene" ' tmp | sed 's/ID=//' | sed 's/\;/\t/g' |  sed 's/gene_biotype=//' | sed 's/ /_/g' | cut -f10 > tmp3
paste tmp2 tmp3 > Vipera_berus_annotations.txt

#combine salmon quants w/ associated gff records for each species/genome combination
ls ERR*_Echis_coloratus_*.genes.tsv > samples.txt
sed -i 's/.genes.tsv//' samples.txt
cat samples.txt | while read line
do
echo "Another one..."
Rscript join_expression.R Vipera_berus_annotations.txt $line\.genes.tsv $line\.tsv
done

ls ERR*_Eublepharis_macularius_*.genes.tsv > samples.txt
sed -i 's/.genes.tsv//' samples.txt
cat samples.txt | while read line
do
echo "Another one..."
Rscript join_expression.R Eublepharis_macularius_annotations.txt $line\.genes.tsv $line\.tsv
done

ls ERR*_Opheodrys_aestivus_*.genes.tsv > samples.txt
sed -i 's/.genes.tsv//' samples.txt
cat samples.txt | while read line
do
echo "Another one..."
Rscript join_expression.R Pantherophis_guttatus_annotations.txt $line\.genes.tsv $line\.tsv
done

ls ERR*_Pantherophis_guttatus_*.genes.tsv > samples.txt
sed -i 's/.genes.tsv//' samples.txt
cat samples.txt | while read line
do
echo "Another one..."
Rscript join_expression.R Pantherophis_guttatus_annotations.txt $line\.genes.tsv $line\.tsv
done

ls ERR*_Python_regius_*.genes.tsv > samples.txt
sed -i 's/.genes.tsv//' samples.txt
cat samples.txt | while read line
do
echo "Another one..."
Rscript join_expression.R Python_molurus_annotations.txt $line\.genes.tsv $line\.tsv
done

##############combine per-tissue expression profiles for each species
#combine individual TPM counts with gff records (except effective length)
Rscript TPM_merge.R ERR216301_Echis_coloratus_skin.tsv ERR216302_Echis_coloratus_skin.tsv Echis_coloratus_skin.tsv
Rscript TPM_merge.R ERR216311_Echis_coloratus_venom_gland.tsv ERR216312_Echis_coloratus_venom_gland.tsv Echis_coloratus_venom_gland.tsv
Rscript TPM_merge.R ERR216319_Echis_coloratus_scent_gland.tsv ERR216326_Echis_coloratus_scent_gland.tsv Echis_coloratus_scent_gland.tsv
Rscript TPM_merge.R ERR216304_Eublepharis_macularius_skin.tsv ERR216306_Eublepharis_macularius_skin.tsv Eublepharis_macularius_skin.tsv
Rscript TPM_merge.R ERR216315_Eublepharis_macularius_salivary_gland.tsv ERR216316_Eublepharis_macularius_salivary_gland.tsv Eublepharis_macularius_salivary_gland.tsv
Rscript TPM_merge.R ERR216322_Eublepharis_macularius_scent_gland.tsv ERR216325_Eublepharis_macularius_scent_gland.tsv Eublepharis_macularius_scent_gland.tsv
Rscript TPM_merge.R ERR216303_Opheodrys_aestivus_skin.tsv ERR216305_Opheodrys_aestivus_skin.tsv Opheodrys_aestivus_skin.tsv
Rscript TPM_merge.R ERR216313_Opheodrys_aestivus_salivary_gland.tsv ERR216314_Opheodrys_aestivus_salivary_gland.tsv Opheodrys_aestivus_salivary_gland.tsv
Rscript TPM_merge.R ERR216320_Opheodrys_aestivus_scent_gland.tsv ERR216321_Opheodrys_aestivus_scent_gland.tsv Opheodrys_aestivus_scent_gland.tsv
Rscript TPM_merge.R ERR216307_Pantherophis_guttatus_salivary_gland.tsv ERR216308_Pantherophis_guttatus_salivary_gland.tsv Pantherophis_guttatus_salivary_gland.tsv
Rscript TPM_merge.R ERR216317_Pantherophis_guttatus_scent_gland.tsv ERR216323_Pantherophis_guttatus_scent_gland.tsv Pantherophis_guttatus_scent_gland.tsv
Rscript TPM_merge.R ERR216299_Python_regius_skin.tsv ERR216300_Python_regius_skin.tsv Python_regius_skin.tsv
Rscript TPM_merge.R ERR216309_Python_regius_salivary_gland.tsv ERR216310_Python_regius_salivary_gland.tsv Python_regius_salivary_gland.tsv
Rscript TPM_merge.R ERR216318_Python_regius_scent_gland.tsv ERR216324_Python_regius_scent_gland.tsv Python_regius_scent_gland.tsv

sed 's/TPM/ERR216298/' ERR216298_Pantherophis_guttatus_skin.tsv | cut -f1-13,15 > Pantherophis_guttatus_skin.tsv
