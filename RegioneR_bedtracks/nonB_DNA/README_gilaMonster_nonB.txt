This folder contains the non-B prediction for the Gila monster reference genome. The predictions include five different non-B types, inverted repeats (IRs, form cruciforms), mirror repeats (filtered to triplexes: TRI), A-phased repeats (APRs; form bent DNA), Z-DNA and G-quadruplexes (G4s). The first three were predicted using non-B-gfa (https://github.com/abcsFrederick/non-B_gfa), and the latter two were predicted using Z-DNA Hunter (https://bioinformatics.ibp.cz/#/analyse/zdna) and G4 Discovery (https://github.com/makovalab-psu/g4Discovery) respectively. 

The results are submitted as BED files. The first three columns are typical [chr, start (0-based), end] of the BED format. However, the 4th column contains different information based on the type of the non-B:
- IR: Repeat-Spacer-Permutations
- TRI: Repeat-Spacer-Permutations
- APR: Repeat-Spacer-Tracts-PurineFraction
- ZDNA: Z-DNA Score
- G4: pqsfinder score

Among all the non-B, only G4s are strand-specific. The 7th column in the G4 prediction BED file shows the G4 Hunter score of the particular G4 motif. 