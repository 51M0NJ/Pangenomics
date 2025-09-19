#This file points to all the files stored by allthebacteria as a TSV
wget -O all_atb_files.tsv https://osf.io/download/r6gcp/ [osf.io]
#Step 1 get E. coli accessions
#Get the URL for the sylph species calls* from the all_atb_files.tsv
grep sylph all_atb_files.tsv  
wget -O sylph.tsv.gz https://osf.io/download/nu5a6/ [osf.io]
gunzip sylph.tsv.gz
#You will want to filter the sylph calls based on the sylph output
#Column 1 contains the Bioproject accession
#Column 4 contains taxonomic abundance
#Column 6 contains the average nucleotide identity
#Column 16 contains the species call
#Possibly you may only want to look at cases where there's a 1:1 map between Bioproject and species call, otherwise it might contaminate your pangenome?
#IMO you probably want to only keep ones where the taxonomic abundance and sequence abundance is high
#Optionally, get the list of high quality samples and use this to filter sylph (defined on this page:  ttps://allthebacteria.readthedocs.io/en/latest/sample_metadata.html [allthebacteria.readthedocs.io])
wget -O hq_set.sample_list.txt.gz https://osf.io/download/m26zn/ [osf.io]
gunzip hq_set.sample_list.txt.gz
gunzip sylph.tsv.gz

# Just the HQ sylph calls
#grep -f hq_set.sample_list.txt sylph.tsv > sylph_hq.tsv

###
awk -F'\t' '
    NR == 1 { print; next } 
    $16 == "Escherichia coli" && #species
    $4  >  0.8 && #Taxonomic_abundance
    $5  >  0.8 && #Sequence_abundance
    $6  >= 95  && #Adjusted_ANI
    $7  >= 20  && #Eff_cov
    $13 >  0.9    #containment
' sylph.tsv > sylph_filtered_Ecoli.tsv
# not significantly more stringent than just taking HQ
grep -f hq_set.sample_list.txt -w sylph_filtered_Ecoli.tsv > sylph_filtered_Ecoli_HQ.tsv
####

#Now just keep the accessions; we will use this to get the correct annotation files.
cut -f 1 sylph_filtered_Ecoli_HQ.tsv > Escherichia_coli_accessions.tsv

#Step 2: Figure out which files contain the E. coli bakta annotations
#From all_atb_files.tsv you can get the list of annotation files for all species
grep --ignore-case bakta all_atb_files.tsv > annotation_files.tsv

#The status files will point you to what to download:
wget -O atb.bakta.incr_release.202408.status.tsv.gz https://osf.io/download/2skzy/ [osf.io]
wget -O atb.bakta.r0.2.status.tsv.gz https://osf.io/download/rxfks [osf.io]
#NB: you need both the r0.2 release and the August 2024 incr_release, but not the r0.1 release (r0.1 is contained in r0.2, see ATB documentation).

gunzip atb.bakta.r0.2.status.tsv.gz
gunzip atb.bakta.incr_release.202408.status.tsv.gz

#Combine these into a single file
cat atb.bakta.r0.2.status.tsv atb.bakta.incr_release.202408.status.tsv > bakta.status.combined.tsv

#Now filter these for the accessions we identified as E. coli above
grep -f Escherichia_coli_accessions.tsv bakta.status.combined.tsv \
grep "PASS" \
> Escherichia_coli.bakta.status.tsv

#The E. coli baktas are batched into ~200 separate xf archives. Sort/unique to get which ones to download
cut -f 5 Escherichia_coli.bakta.status.tsv | sort -u > Escherichia_coli.bakta.batches.tsv

#Now back to the annotation_files.tsv for the URLs of the batches
grep -f Escherichia_coli.bakta.batches.tsv annotation_files.tsv > Escherichia_coli.bakta.batches.URLs.tsv

#At this point it's up to you whether you want to download these all in one go (takes a lot of space), or do one batch at a time in the steps below.
#Step 3: Extract all gene/aa/aa_hexdigest pairs from each strain for a 1x5000ish matrix
#Step 3a: Store aa/aa_hexdigest pairs into a dictionary file for each gene (i.e. only add to file if a new aa/aa_hexdigest pair is found)
#Step 4: Full join all strains from step 3 for a 400k x ?? matrix (?? = pan-genome size?)
#Step 5: Determine which aa_hexdigests correspond to defective proteins from the dictionary files (e.g. frameshift, big deletion, possibly even defects at functional residues?? Danna has done this for DNA replication fidelity genes and antibiotic resistance genes) --> classify as functional/non-functional instead of presence/absence
#Step 6: Association statistics (see Fiona)
