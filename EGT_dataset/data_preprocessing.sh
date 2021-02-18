mkdir MitoCOGs
cd MitoCOGS/
mkdir raw_data
cd raw_data
wget -r ftp://ftp.ncbi.nih.gov/pub/koonin/MitoCOGs
mv ftp.ncbi.nih.gov/pub/koonin/MitoCOGs/* .
rm -rf ftp.ncbi.nih.gov/
cd ..

# LIST OF ERRORS IN/BETWEEN "raw_data/MitoCOGs.txt" AND "raw_data/mitochondrial_proteomes.txt":
# Zancudomyces culisetae	-> TaxNbID:	4889 (OLD: "mitochondrial_proteomes.txt") -> 1213189 (NEW: "MitoCOGs.txt")	=> TODO: CORRECTION OF mitochondrial_proteomes.txt
# Ogataea angusta			-> TaxNbID:	4905 (OLD: "mitochondrial_proteomes.txt" & "MitoCOGs.txt") -> 870730 (NEW)	=> Both ID used in "MitoCOGs.txt" (Let's think they are 2 species!!!)
# 
# WRONG WRITTING IN "mitochondrial_proteomes.txt" (Add "," introducing error in column delimitation !!!!):
# Contracaecum rudolphii B Bullini et al., 1986,646521,cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Nematoda; Chromadorea; Ascaridida; Ascaridoidea; Anisakidae; Contracaecum; Contracaecum rudolphii group
#                                        *
# Leptocephalus sp. 'type II larva' (Smith, 1989),556254,cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Actinopterygii; Actinopteri; Neopterygii; Teleostei; Elopocephalai; Elopocephala; Elopomorpha; Anguilliformes; Nettastomatidae; Leptocephalus
#                                         * 

cut -d, -f2 raw_data/nuclear_proteomes.txt > listNP
cut -d, -f2 raw_data/nuclear_encoded_orthologs.txt|sort -n|uniq> listNEO
wc -l listNP listNEO
# There 57 species in "nuclear_proteomes.txt" BUT 52 in "nuclear_encoded_orthologs.txt" -->
diff listNEO listNP
# List of species in "nuclear_proteomes.txt" not in "nuclear_encoded_orthologs.txt"
# 	-> 184922, 284813, 294381, 412133 & 578460
rm listNEO listNP


# Sort MitoCOGs by species ID and by MitoCOGs ID
sort -t "," -k2n,2n -k4d,4d raw_data/MitoCOGs.txt > MitoCOGs_sorted_by_species_ID.txt
sort -t "," -k4d,4d -k2n,2n raw_data/MitoCOGs.txt > MitoCOGs_sorted_by_MitoCOGs.txt
cut -d"," -f2 MitoCOGs_sorted_by_species_ID.txt|uniq > listSpeciesMito
cut -d"," -f4 MitoCOGs_sorted_by_MitoCOGs.txt|uniq > listMitoCOGs
wc -l listMitoCOGs listSpeciesMito
# There are 140 MitoCOGs among all mitochondrial-encoded proteins from 2487 different species
rm listMitoCOGs listSpeciesMito MitoCOGs_sorted_by_MitoCOGs.txt

# !!! 2 different ID used for Physcomitrella patens between
# "raw_data/nuclear_proteomes.txt" (145481) & "raw_data/mitochondrial_proteomes.txt" (3218)
# Change ",145481," by ",3218," in "raw_data/nuclear_encoded_orthologs.txt"
sed 's/,145481,/,3218,/g' raw_data/nuclear_encoded_orthologs.txt > nuclear_encoded_orthologs_corrected.txt
# Sort nuclear_encoded_orthologs by species ID and by MitoCOGs ID
sort -t "," -k2n,2n -k4d,4d nuclear_encoded_orthologs_corrected.txt > nuclear_encoded_orthologs_corrected_sorted_by_species_ID.txt
sort -t "," -k4d,4d -k2n,2n nuclear_encoded_orthologs_corrected.txt > nuclear_encoded_orthologs_corrected_sorted_by_MitoCOGs.txt
cut -d"," -f2 nuclear_encoded_orthologs_corrected_sorted_by_species_ID.txt|uniq > listSpeciesNuc
cut -d"," -f4 nuclear_encoded_orthologs_corrected_sorted_by_MitoCOGs.txt|uniq > listNucCOGs
wc -l listNucCOGs listSpeciesNuc
# There are 45 MitoCOGs among all nuclear-encoded proteins from 52 different species 
rm listNucCOGs nuclear_encoded_orthologs_corrected.txt nuclear_encoded_orthologs_corrected_sorted_by_MitoCOGs.txt


# sort -t "," -k1d,1d raw_data/MitoCOGs_annotation.txt | grep "MitoCOG" > MitoCOGs_annotation_sorted_by_ID.txt
# sort -t "," -k2d,2d MitoCOGs_annotation_sorted_by_ID.txt > MitoCOGs_annotation_sorted_by_names.txt
# grep -v "hypothetical" MitoCOGs_annotation_sorted_by_names.txt > MitoCOGs_annotation_sorted_by_names-hypothetical_protein.txt


### Add command line to change 145481 in 3218 (P. patens)
sort -t "," -k2n,2n raw_data/nuclear_proteomes.txt > nuclear_proteomes_sorted_by_species_ID.txt
# !!! 2 different ID used for Physcomitrella patens between
# "raw_data/nuclear_proteomes.txt" (145481) & "raw_data/mitochondrial_proteomes.txt" (3218)
# Change ",145481," by ",3218," in "nuclear_proteomes_sorted_by_species_ID.txt"
sed -i 's/,145481,/,3218,/g' nuclear_proteomes_sorted_by_species_ID.txt
grep -f listSpeciesNuc nuclear_proteomes_sorted_by_species_ID.txt > nuclear_proteomes_sorted_by_species_ID_filtered.txt
rm listSpeciesNuc nuclear_proteomes_sorted_by_species_ID.txt


mkdir dataset
# grep Opisthokonta nuclear_proteomes_sorted_by_species_ID_filtered.txt > 28_Opisthokonta.txt 
# grep Fungi 28_Opisthokonta.txt > 10_Fungi.txt
# grep Metazoa 28_Opisthokonta.txt > 17_Metazoa.txt

grep Viridiplantae nuclear_proteomes_sorted_by_species_ID_filtered.txt > 9_Viridiplantae.txt

# grep -v Opisthokonta nuclear_proteomes_sorted_by_species_ID_filtered.txt \
# | grep -v Viridiplantae | grep -v "Cyanophora paradoxa" | grep -v "Cyanidioschyzon merolae" \
# > 13_Others_Eukaryota.txt


# Add Cyanophora paradoxa & Cyanidioschyzon merolae to the 9_Viridiplantae
mkdir 11_Plantae_dataset
mv 9_Viridiplantae.txt 11_Plantae_dataset/11_Plantae.txt
grep "Cyanophora paradoxa" nuclear_proteomes_sorted_by_species_ID_filtered.txt >> 11_Plantae_dataset/11_Plantae.txt
grep "Cyanidioschyzon merolae" nuclear_proteomes_sorted_by_species_ID_filtered.txt >> 11_Plantae_dataset/11_Plantae.txt
sort -t "," -k2n,2n -o 11_Plantae_dataset/11_Plantae_dataset.txt 11_Plantae_dataset/11_Plantae.txt
rm nuclear_proteomes_sorted_by_species_ID_filtered.txt


### Get mitochondrial-encoded MitoCOGS for the 11_Plantae dataset
for spe in $(cut -d',' -f2 11_Plantae_dataset/11_Plantae.txt); do
	grep ","$spe"," MitoCOGs_sorted_by_species_ID.txt
done > 11_Plantae_dataset/11_Plantae_MitoCOGs_sorted_by_species_ID.txt
rm MitoCOGs_sorted_by_species_ID.txt
sort -t "," -k4d,4d -k2n,2n 11_Plantae_dataset/11_Plantae_MitoCOGs_sorted_by_species_ID.txt \
> 11_Plantae_dataset/11_Plantae_MitoCOGs_sorted_by_MitoCOGs.txt
cut -d"," -f4 11_Plantae_dataset/11_Plantae_MitoCOGs_sorted_by_MitoCOGs.txt | uniq \
> 11_Plantae_dataset/list11Plantae_MitoCGOS_mitochondrial-encoded
wc -l 11_Plantae_dataset/list11Plantae_MitoCGOS_mitochondrial-encoded
# 68 MitoCOGs with only mitochondrial-encoded proteins in the 11 Plantae dataset 

### Get nuclear-encoded COGs for the 11_Plantae dataset
for spe in $(cut -d',' -f2 11_Plantae_dataset/11_Plantae.txt); do
	grep ","$spe"," nuclear_encoded_orthologs_corrected_sorted_by_species_ID.txt
done > 11_Plantae_dataset/11_Plantae_nuclear_encoded_orthologs_corrected_sorted_by_species_ID.txt
rm nuclear_encoded_orthologs_corrected_sorted_by_species_ID.txt
sort -t "," -k4d,4d -k2n,2n 11_Plantae_dataset/11_Plantae_nuclear_encoded_orthologs_corrected_sorted_by_species_ID.txt \
> 11_Plantae_dataset/11_Plantae_nuclear_encoded_orthologs_corrected_sorted_by_MitoCOGs.txt
cut -d"," -f4 11_Plantae_dataset/11_Plantae_nuclear_encoded_orthologs_corrected_sorted_by_MitoCOGs.txt | uniq \
> 11_Plantae_dataset/list11Plantae_MitoCGOS_nuclear-encoded
wc -l 11_Plantae_dataset/list11Plantae_MitoCGOS_nuclear-encoded
# 41 MitoCOGs with only nuclear-encoded proteins in the 11 Plantae dataset

### Get list of MitoCOGs containing nuclear- AND mitochondrial-encoded proteins
comm -12 11_Plantae_dataset/list11Plantae_MitoCGOS_mitochondrial-encoded \
11_Plantae_dataset/list11Plantae_MitoCGOS_nuclear-encoded \
> 11_Plantae_dataset/list11Plantae_MitoCOGs_nuc+mito
# 28 MitoCOGs with nuclear- AND mitochondrial-encoded proteins in the 11 Plantae dataset

### Get the 
grep -f 11_Plantae_dataset/list11Plantae_MitoCOGs_nuc+mito \
11_Plantae_dataset/11_Plantae_nuclear_encoded_orthologs_corrected_sorted_by_MitoCOGs.txt \
> 11_Plantae_dataset/11_Plantae_NUC+mito.txt
# !!! ADD MANUALLY ",nucleus" to the end of each line !!!
grep -f 11_Plantae_dataset/list11Plantae_MitoCOGs_nuc+mito \
11_Plantae_dataset/11_Plantae_MitoCOGs_sorted_by_MitoCOGs.txt \
> 11_Plantae_dataset/11_Plantae_nuc+MITO.txt
# !!! ADD MANUALLY ",mitochondrion" to the end of each line !!!
cat 11_Plantae_dataset/11_Plantae_nuc+MITO.txt 11_Plantae_dataset/11_Plantae_NUC+mito.txt \
| sort -t "," -k4d,4d -k2n,2n -k5d,5d > 11_Plantae_dataset/11_Plantae_MitoCOGs_nuc+mito_sorted_by_MitoCOGs.txt
rm 11_Plantae_dataset/11_Plantae_nuc+MITO.txt 11_Plantae_dataset/11_Plantae_NUC+mito.txt
# 326 proteins: 142 nuclear-encoded & 184 mitochondrial-encoded


#### Corrected "raw_data/nuclear_encoded_orthologs_sequences.fa"
sed 's/\[ubiquinone\]/\(ubiquinone\)/g' raw_data/nuclear_encoded_orthologs_sequences.fa \
| sed 's/Cyanidioschyzon merolae strain 10D/Cyanidioschyzon merolae/g' \
| sed 's/Physcomitrella patens subsp. patens/Physcomitrella patens/g' \
> 11_Plantae_dataset/nuclear_encoded_orthologs_sequences_corrected.fa

### Copy "MitoCOGs_sequences.fa" in 11_Plantae_dataset/
cp raw_data/MitoCOGs_sequences.fa 11_Plantae_dataset/


### Apply EndoRex
python3 ../pipeline_for_endorex.py \
11_Plantae_dataset/11_Plantae_MitoCOGs_nuc+mito_sorted_by_MitoCOGs.txt \
11_Plantae_dataset/11_Plantae.txt \
11_Plantae_dataset/11_Plantae_species_tree.nwk \
11_Plantae_dataset/MitoCOGs_sequences.fa \
11_Plantae_dataset/nuclear_encoded_orthologs_sequences_corrected.fa \
11_Plantae_dataset/RESULTS \
&> 11_Plantae_dataset/11_Plantae_EndoRex.log


# "11_Plantae_dataset/11_Plantae_species_tree.nwk" was created manually based on topology in Kannan et al., 2014