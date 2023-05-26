# TCR Repertoire Analysis
T-cell receptor (TCR) repertoire analysis helps study the diversity and characteristics of T-cell receptors expressed by T cells. TCRs play a crucial role in the adaptive immune response by recognizing and binding to specific antigens presented by antigen-presenting cells.

TCR repertoire analysis involves examining the diversity, clonality, and distribution of TCR sequences within a given population of T cells. This analysis provides insights into the immune response, such as the recognition of pathogens, autoimmune diseases, cancer immunotherapy, and vaccine development. The advancement of next generation sequencing (NGS) and single-cell technologies have enabled such in depth characterization of immune cell repertoires. 

TCR repertoire analysis helps in understanding immune responses, identifying disease-associated TCR signatures, and developing immunotherapies. The field continues to evolve, and new techniques and approaches are being developed to further enhance our understanding of TCR diversity and function.

Common techniques and approaches used in TCR repertoire analysis:

1. High-throughput sequencing: Next-generation sequencing (NGS) technologies enable the sequencing of millions of TCR DNA or RNA sequences from a sample. This approach provides a comprehensive view of the TCR repertoire diversity and allows for quantitative analysis.

2. TCR sequencing protocols: Specific protocols are used to amplify and sequence the TCR regions of interest, such as the complementarity-determining region 3 (CDR3). CDR3 is the most variable region of the TCR and is crucial for antigen recognition.

3. Data analysis: Bioinformatics tools and algorithms are used to process and analyze the massive amount of sequencing data generated. This includes preprocessing, quality control, error correction, sequence alignment, clonotype identification, diversity analysis, and clonality assessment.

4. Clonotype identification: Clonotypes are unique TCR sequences derived from the same T-cell lineage. Clonotype identification involves grouping similar TCR sequences into clonally related clusters based on sequence similarity and other criteria.

5. Diversity analysis: Measures such as clonality, richness, evenness, and diversity indices are calculated to assess the diversity and distribution of TCR sequences within a sample or across different samples.

6. Repertoire comparison: TCR repertoire analysis can compare TCR repertoires between different individuals, tissues, or time points to understand changes in immune responses or disease progression.

7. Functional analysis: TCR repertoire analysis can be integrated with functional assays to investigate the functionality of TCRs, including antigen specificity, cytokine production, and cytotoxicity.


# GLIPH2 Differential Motif Analysis
[GLIPH2 (Grouping of Lymphocyte Interactions by Paratope Hotspots 2)]([http://50.255.35.37:8080/](https://www.nature.com/articles/s41587-020-0505-4)) is a computational tool that aids in grouping TCR sequences based on shared patterns in the complementarity-determining region 3 (CDR3), facilitating understanding of TCR diversity and antigen recognition patterns in immune responses. In this repository, I provide a code template for an in-depth characterization of T-cell repertoires. I utilize the output from GLIPH2 to perform downstream analysis, distinguishing clonally expanded motifs. GLIPH clusters are initially filtered for Fisher's score < 0.05. Clonally expanded and enriched motifs are identified by comparing the summed contribution scores of samples from each cluster (summed template frequency) using the poisson.test function for p-value statistics. The poisson.test function performs an exact test of a simple null hypothesis regarding the ratio between two rate parameters. The enriched motifs' log2-fold change and p-values are then plotted in volcano plots to identify expanded TCRs.

With enriched motifs, it is possible to track their amino acid sequences within the repertoire. V gene and V-J pairing preferences across donors can be visualized as chord diagrams to identify similar preferences across donors. Additionally, shared motifs among different groups can be depicted in a network diagram or link plot, illustrating the frequency of contribution from each group. Here, I've provided a simple code template to make conventional chord plots as well as link plots.

