This Python script processes genetic data from a Thousand Genome Project VCF (Variant Call Format) file, performs population analysis, and visualizes the results using dimensionality reduction techniques (PCA, t-SNE, and UMAP) with interactive Bokeh plots. Below is a detailed explanation of the code, broken down into key sections:
1. Imports and Setup
The script begins by importing necessary libraries for data processing, visualization, and dimensionality reduction:

Standard Libraries: collections, gzip, itertools, os, time for handling data structures, file operations, and timing.
Numerical and Data Processing: numpy for array operations, Counter and defaultdict for counting and dictionary management.
Visualization: matplotlib.pyplot (though not used in the provided code), bokeh for interactive plotting, and bokeh.palettes for color schemes.
Dimensionality Reduction: sklearn.decomposition.PCA, sklearn.manifold.TSNE, and umap for reducing high-dimensional genetic data to 2D for visualization.
Interactivity: ipywidgets for potential interactive features (not used in the provided code).
File Paths: Defines paths to the VCF file (ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz), population description file (20131219.populations.tsv), and population file (affy_samples.20141118.panel).
2. SNP Class
The snp class is defined to parse and process individual lines from the VCF file, which contains genetic variant data.

Initialization (__init__)
Input: A line from the VCF file (as a string) and optional parameters (select for specific individuals, autosome_only to restrict to autosomes).
Processing:
Splits the line into fields (e.g., chromosome, position, reference allele, alternate allele, genotypes).
Performs validation checks:
Skips header lines (starting with #).
Ensures at least 5 fields are present.
If autosome_only=True, restricts to chromosomes 1–22.
Ensures reference and alternate alleles are valid (A, C, G, T).
Checks if the SNP passes quality control (filter field is PASS or .).
Extracts genotype strings (e.g., 0/0, 0/1, 1/1) for each individual.
Genotype Array:
Creates a NumPy array (genotype_array) to store genotypes numerically:
0/0 → 0 (homozygous reference).
0/1 or 1/0 → 1 (heterozygous).
1/1 → 2 (homozygous alternate).
./. → -1 (missing data).
Counts missing genotypes and handles unexpected formats by marking the SNP as failed.
Key Attributes:
chrom: Chromosome number.
position: SNP position on the chromosome.
rid: SNP identifier.
ref_allele, alt_allele: Reference and alternate alleles.
genotype_array: Numeric representation of genotypes.
failed, failure_cause: Flags and reasons for invalid SNPs.
3. Reading and Processing the VCF File
This section reads the VCF file and constructs a genotype matrix.

File Reading:
Opens the gzipped VCF file using gzip.open in text mode ('rt').
Iterates through the file, processing every number_of_lines_to_skip (10) lines to reduce computational load.
Skips header lines (starting with #) and failed SNPs.
For valid SNPs, extracts the genotype array, SNP ID (rid), and position (chrom, position).
Data Storage:
genotype_matrix: List of genotype arrays (rows = SNPs, columns = individuals).
genotype_names: List of SNP IDs.
genotype_positions: List of [chromosome, position] pairs.
Error Tracking:
Counts errors (error_count) and prints the first 10 failure reasons for debugging.
Output:
Prints the shape of the genotype matrix, number of SNPs processed, number of errors, and runtime.
Transpose:
Transposes the genotype matrix (transposed_genotype_matrix) so rows represent individuals and columns represent SNPs, preparing it for dimensionality reduction.
4. Population Data Processing
This section reads population metadata and organizes it for visualization.

Population File (affy_samples.20141118.panel):
Reads individual-to-population mappings.
Creates:
population_by_individual: Maps individual IDs to population codes.
individuals_by_population: Maps population codes to lists of individuals.
Population Description File (20131219.populations.tsv):
Reads population metadata (e.g., population name, code, continent).
Creates:
name_by_code: Maps population codes to full names.
pop_by_continent: Maps continents to lists of population codes.
continent_by_population: Maps population codes to continents.
Color Assignment:
Assigns colors to populations using Category20b and custom palettes (Purples, Greens, PuOr, RdGy) for visualization.
Ensures distinct colors for populations, with manual overrides for specific populations (e.g., CHS, STU, LWK).
5. Extracting Individuals from VCF
Reads the VCF file again to extract individual IDs from the header line (starting with #CHROM).
Builds lspop, a list of population codes for each individual, defaulting to "missing" if no population is assigned.
Creates indices_of_population_members, a dictionary mapping population codes to lists of individual indices in the VCF file.
Prints the total number of populations and the maximum/minimum number of individuals per population.
6. Dimensionality Reduction and Visualization
The script applies three dimensionality reduction techniques to the transposed genotype matrix and visualizes the results using Bokeh.

PCA (Principal Component Analysis):
Method: sklearn.decomposition.PCA with 2 components.
Process:
Fits and transforms the transposed genotype matrix to project it into 2D (PC1, PC2).
Creates a dictionary (pca_df) with PC1, PC2, population labels, individual IDs, and colors.
Visualization:
Uses Bokeh to create an interactive scatter plot.
Each point represents an individual, colored by population.
Includes a legend (clickable to hide/show populations).
Saves the plot as pca_genotypes.html.
t-SNE (t-Distributed Stochastic Neighbor Embedding):
Method: sklearn.manifold.TSNE with 2 components, perplexity=30, n_iter=1000, and random_state=42.
Process:
Fits and transforms the genotype matrix to 2D.
Creates a ColumnDataSource with t-SNE coordinates, colors, and population labels.
Visualization:
Creates a Bokeh scatter plot with interactive tools (pan, zoom, hover, save).
Saves the plot as tsne_plot.html.
UMAP (Uniform Manifold Approximation and Projection):
Method: umap.UMAP with n_neighbors=15, min_dist=0.1, metric='euclidean', and random_state=42.
Process:
Fits and transforms the genotype matrix to 2D.
Creates a ColumnDataSource with UMAP coordinates, colors, and population labels.
Visualization:
Creates a Bokeh scatter plot similar to t-SNE.
Saves the plot as umap_plot.html.
Key Outputs
Genotype Matrix: A matrix where rows are SNPs and columns are individuals, with genotypes encoded as 0, 1, 2, or -1.
Population Metadata: Organized mappings of individuals to populations and populations to continents.
Visualizations: Three interactive HTML plots (pca_genotypes.html, tsne_plot.html, umap_plot.html) showing individuals in 2D space, colored by population.
