% Read omics data from Davidi 2016 and Map BIGG reaction IDs and gene names to KEGG reaction IDs

% --------------------------------------------
% Flux data 

flux_data = sbtab_table_load([ cmb_basedir '/resources/escherichia_coli/data-omics/davidi_2016/Davidi_2016_flux_mmolgCDW_s.tsv']);

conversion_bigg_kegg = load_any_table([cmb_basedir '/resources/annotations/bigg_reaction_to_kegg.tsv']);

conversion_bigg_kegg = [conversion_bigg_kegg; ...                   
                    {'PFK','R04779'};...
                    {'G6PDH2r','R00835'}; ...
                    {'PGK_reverse','R01512'}];

flux_data = sbtab_table_add_annotation(flux_data,'Reaction','Reaction:Identifiers:kegg.reaction',conversion_bigg_kegg);
 
sbtab_table_save(flux_data, [ cmb_basedir '/resources/escherichia_coli/data-omics/davidi_2016/Davidi_2016_flux_mmolgCDW_s_KEGG_IDs.tsv']);

% --------------------------------------------
% Protein data 

protein_data = sbtab_table_load([ cmb_basedir '/resources/escherichia_coli/data-omics/davidi_2016/Davidi_2016_protein_abundance_mmolgCDW.tsv']);

conversion_gene_kegg_reaction = load_any_table([ cmb_basedir '/resources/annotations/escherichia_coli_orf2kegg_reaction.tsv']);

protein_data = sbtab_table_add_annotation(protein_data,'Gene:LocusName','Reaction:Identifiers:kegg.reaction',conversion_gene_kegg_reaction);

sbtab_table_save(protein_data, [ cmb_basedir '/resources/escherichia_coli/data-omics/davidi_2016/Davidi_2016_protein_abundance_mmolgCDW_KEGG_IDs.tsv']);
