nextclade run \
   --input-dataset data/sars-cov-2 \
   --output-all=/home/soniali/Desktop/02_china_recom_github/3_recom/nextclade_output \
/home/soniali/Desktop/02_china_recom_github/3_recom/FR.1.1_XCN_EG.5.1.1_sequences.fasta

iqtree -s FR.1.1_XCN_EG.5.1.1_masked_right.fasta -bb 1000 -alrt 1000 -o Reference

# XCN	Recombinant lineage of FR.1.1 and EG.5.1.1 (breakpoint between 22630-22663), China (Anhui/Shanghai), from sars-cov-2-variants/lineage-proposals#842
# XDL	Recombinant lineage of EG.5.1.1 and non-EG.5 (breakpoint between 25573-29624) with S:L455F and S:I1114T, China, from #2351