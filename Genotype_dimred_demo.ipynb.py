import matplotlib.pyplot as plt
import collections
from collections import defaultdict
import gzip
import itertools
import numpy as np
import os
import time

from ipywidgets import interact
import bokeh
import bokeh.io
from bokeh.io import push_notebook

# Import colour palettes for later on
from bokeh.palettes import Category20b
from bokeh.palettes import Purples
from bokeh.palettes import Greens
from bokeh.palettes import YlOrBr
from bokeh.palettes import YlOrRd
from bokeh.palettes import PuOr
from bokeh.palettes import RdGy

from bokeh.plotting import figure, save, output_file
from bokeh.models import ColumnDataSource
from bokeh.palettes import Category20b, Purples, Greens, PuOr, RdGy

# Dimension reduction tools
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.impute import SimpleImputer
import umap 

data_dir = '/home/omer'
# These are the names of the files we use
vcf_name = 'ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz'
pop_desc_name = '20131219.populations.tsv'
pop_file_name = 'affy_samples.20141118.panel'

vcf_file = os.path.join(data_dir, vcf_name)
population_description_file = os.path.join(data_dir, pop_desc_name)
population_file = os.path.join(data_dir, pop_file_name)

from collections import Counter

class snp(object):

    def __init__(self, line, select=False, autosome_only =True):
        """The initialization method takes in a line from the vcf file, as a string, 
        and records the relevant information. 
        line: a string from a vcf file
        select: a list of positions of individuals to be analyzed, where positions run from 0 to 
        nInd-1, the number of individuals
        """ 
        
        split_line = line.split()  #  First break down the line into a list of each field
        
        self.failed = False  # A label that we will set to True if something goes wrong.
        
        if line.startswith('#'):
            self.failed = True
            self.failure_cause = "line was a header line, not a snp"
            return
        
        if len(split_line)<=5:
            self.failed = True
            self.failure_cause = "incorrectly formatted line, should have at least 5 fields " + line
            return
          
        self.chrom = split_line[0]
        print(f"Processing chromosome: {self.chrom}")  # Debug print
        if autosome_only:
            if self.chrom not in ["%d" % (i,) for i in range(1,23)]:
                self.failed = True
                self.failure_cause = "not recognized as an autosome while autosome_only set to True"
                return
        
        self.chrom = int(split_line[0]) # Chromosome (numbered)
        self.position = int(split_line[1])  # The coordinates of the snp
        self.rid = split_line[2] # Name/Record ID
        self.ref_allele = split_line[3]
        self.alt_allele = split_line[4] # The alterate allele according to the vcf; also a string 
        # Only accept snps in ACGT. 
        if self.ref_allele not in ["A","C","G","T"] or self.alt_allele not in ["A","C","G","T"]:
            self.failed = True
            self.failure_cause = "ref or alt not in ACGT"
            return
        self.filter = split_line[6]  # See vcf format specifications for the interpretation of 
                                    # the filter field
        if self.filter not in ['PASS', '.'] :  # PASS indicates a SNP that passed all QC filters.
            self.failed = True
            self.failure_cause = self.filter
            return
              
        self.genotype_strings = split_line[9:]

        # Prepare a list that will contain the transformed genotypes. 
        # Since we already know how long the list will be, it makes sense 
        # to create an array of zeros of the same length as self.gtypes, 
        
        self.genotype_array = np.zeros(len(self.genotype_strings), dtype = np.int8)             

        # Count the number of each genotype. 
        # There may be different strings giving the same genotype so we increment the 
        # counts found so far for the genotype by the number of times the  
        # For example, "0/0" and "0\0" give homref, and "0|1" and "1|0" give het
        
        n_missing = 0
        for index,genotype_string in enumerate(self.genotype_strings):
            if genotype_string == './.':
                n_missing +=1 
                self.genotype_array[index]=-1
                continue # missing data will be left as 0
            allele_0 = genotype_string[0] # Get the first allele (as a string)
            allele_1 = genotype_string[2]
            if (allele_0=='1' and allele_1=='1'): # Use rstrip because windows machines will occasionally have extra \n
                self.genotype_array[index]=2
            elif ((allele_0=='0' and allele_1=='1') or (allele_0=='1' and allele_1=='0')):
                self.genotype_array[index]=1   
            elif (allele_0=='0' and allele_1=='0'):
                # The array was initialized to zero, so nothing to do here!
                continue
            else:
                print(("unknown genotype", genotype_string))
                self.failed=True
                self.failedreason="unknown genotype"
                return
            
number_of_lines_to_skip = 10

start_time = time.time()

genotype_matrix = []  # Will contain our numerical genotype matrix. 
genotype_positions = []
genotype_names = []
x = 0
error_count = 0

with gzip.open(vcf_file,'rt') as f:
    count = 0
    for line in f:
        count+=1
        if count % number_of_lines_to_skip == 0:
            if line.startswith("#") or snp(line).failed:
                if snp(line).failure_cause != "line was a header line, not a snp":
                    error_count += 1
                    if x < 10:
                        print('Failed: ' + snp(line).failure_cause)
                        x+=1
                continue
            
            return_snp = snp(line, autosome_only=False)
            genotype_matrix.append(return_snp.genotype_array)
            genotype_names.append(return_snp.rid)
            genotype_positions.append([return_snp.chrom, return_snp.position])

end_time = time.time()

print(f"Genotype matrix shape: {np.array(genotype_matrix).shape}")
print(f"Number of SNPs processed: {len(genotype_matrix)}")
print(f"Number of errors: {error_count}")

print("Run time in seconds: " + str(end_time - start_time))

transposed_genotype_matrix = np.array(genotype_matrix).transpose()

population_by_individual = defaultdict(int)
individuals_by_population = defaultdict(list)  # A dictionary containing all the individuals in a given population

for line in open(population_file,'r'):
    split_line = line.split()
    if split_line[0] == 'sample':  # header line
        continue

    sample_name = split_line[0]
    population_name = split_line[1]
    population_by_individual[sample_name] = population_name
    individuals_by_population[population_name].append(sample_name) 

populations = list(individuals_by_population.keys())

name_by_code = {}  # A dictionary giving the full name of each population code
pop_by_continent = {}  # A dictionary giving the code of each population within a continent  
continent_by_population = {}  # A dictionary giving the continent for each population code
for line in open(population_description_file,'r'):
    split_line = line.split('\t')
    if split_line[0] in ['Population Description','Total','']:  # header or footer
        continue
    name_by_code[split_line[1]] = split_line[0]
    continent_by_population[split_line[1]] = split_line[2]
    try: 
        pop_by_continent[split_line[2]].append(split_line[1])
    except KeyError:
        pop_by_continent[split_line[2]] = [split_line[1]]

continents = list(pop_by_continent.keys()) 
    
    
# Populations listed by continent
pops=[]
for continent in continents:
    pops.extend(pop_by_continent[continent])

color_dict = {}
for i, cont in enumerate(continents): 
    for j, pop in enumerate(pop_by_continent[cont]):
        color_dict[pop] = Category20b[20][4*i+j%4]

# Colour palette above only really supports groups of 4 so we have to manually specify a few colours for the 5th/6th
# members of a group

color_dict['CHS'] = Purples[9][4]# purple
color_dict['STU'] = Greens[9][6] # green
color_dict['LWK'] = PuOr[11][-1] # brown
color_dict['MSL'] = PuOr[11][-2] # rusty brown
color_dict['YRI'] = PuOr[11][-3] # cappucino w/ extra milk (stirred)
color_dict['CEU'] = RdGy[11][-3]

for line in gzip.open(vcf_file,'rt'):
    if line.startswith("#"):
        if not line.startswith("##"):
            # Extract the individuals for the population, as a list of strings
            # Windows users may have trailing \n characters
            individuals = line.split()[9:]
            # Once we've extracted the individuals, we can exit the loops. 
            break

# Build a list of populations for each indiviudal in the vcf file
lspop = []
for ind in individuals:
    pop = population_by_individual[ind]
    if pop == 0:
        lspop.append("missing")
    else:
        lspop.append(pop)

        
indices_of_population_members = defaultdict(list)

for index,individual in enumerate(individuals):
    try:
        indices_of_population_members[population_by_individual[individual]].append(index)
    except KeyError: # We do not have population info for this individual
        continue

count = 0

for p in pop_by_continent:
    count+=len(pop_by_continent[p])
    
print(count)

imax = 0
imin = 200
for i in indices_of_population_members:
    imax = max(len(indices_of_population_members[i]),imax)
    imin = min(len(indices_of_population_members[i]),imin)
    
print(imax, imin)

pca = PCA(n_components=2)
pca_result = pca.fit_transform(transposed_genotype_matrix)

pca_df = {
    'PC1': pca_result[:, 0],
    'PC2': pca_result[:, 1],
    'population': lspop,
    'individual': individuals,
    'color': [color_dict.get(pop, "#000000") for pop in lspop]
}

# Bokeh ile interaktif görselleştirme
source = ColumnDataSource(data=pca_df)

p = figure(title="PCA of Genotype Matrix", width=800, height=600)
p.scatter(x='PC1', y='PC2', source=source, color='color', legend_field='population', size=8, alpha=0.6)

p.legend.location = "top_right"
p.legend.click_policy = "hide"

output_file("pca_genotypes.html")
save(p)

tsne = TSNE(n_components=2, perplexity=30, n_iter=1000, random_state=42)
tsne_results = tsne.fit_transform(transposed_genotype_matrix)

# Plot için veri hazırlığı
colors = [color_dict[pop] if pop in color_dict else "#999999" for pop in lspop]

source_tsne = ColumnDataSource(data=dict(
    x=tsne_results[:, 0],
    y=tsne_results[:, 1],
    color=colors,
    label=lspop
))

# Bokeh ile plot
tsne_plot = figure(title="T-SNE Plot of Genotype Data", tools="pan,wheel_zoom,reset,hover,save",
                   width=800, height=600)
tsne_plot.circle('x', 'y', size=8, color='color', legend_field='label', source=source_tsne, alpha=0.6)
tsne_plot.legend.click_policy = "hide"

output_file("tsne_plot.html")
save(tsne_plot)


# UMAP
reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric='euclidean', random_state=42)
umap_results = reducer.fit_transform(transposed_genotype_matrix)

# Plot için veri hazırlığı
source_umap = ColumnDataSource(data=dict(
    x=umap_results[:, 0],
    y=umap_results[:, 1],
    color=colors,
    label=lspop
))

# Bokeh ile plot
umap_plot = figure(title="UMAP Plot of Genotype Data", tools="pan,wheel_zoom,reset,hover,save",
                   width=800, height=600)
umap_plot.circle('x', 'y', size=8, color='color', legend_field='label', source=source_umap, alpha=0.6)
umap_plot.legend.click_policy = "hide"

output_file("umap_plot.html")
save(umap_plot)