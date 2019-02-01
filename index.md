### Overview
Genetic Counseling is a component of the umbrella [OpenMendel](https://openmendel.github.io) project. This analysis option computes risks to individuals in pedigrees segregating Mendelian diseases. As a conditional probability, a genetic risk involves two likelihoods, a numerator likelihood with the riskee having an disease genotype or phenotype at the disease locus and a denominator likelihood with the riskee having an unknown or non-specific phenotype at the disease locus. Depending on the problem, complicating features such as age of onset, mutation, linked markers, and biochemical tests come into play.

### Appropriate Problems and Data Sets
Genetic Counseling will handle fairly complicated problems involving mutation, reduced penetrance, and linked markers. Genetic Counseling does not provide empiric risks or theoretical risks under models for genetic heterogeneity or polygenic inheritance.

### Installation
*Note: The three OpenMendel packages (1) [SnpArrays](https://openmendel.github.io/SnpArrays.jl/latest/), (2) [MendelSearch](https://openmendel.github.io/MendelSearch.jl), and (3) [MendelBase](https://openmendel.github.io/MendelBase.jl) must be installed before any other OpenMendel package will run. It is easiest if these three packages are installed in the above order and before any other OpenMendel package.*

Within Julia, use the package manager to install MendelGeneticCounseling:

    pkg> add https://github.com/OpenMendel/MendelGeneticCounseling.jl.git

This package supports Julia v1.0+

### Input Files
The Genetic Counseling analysis package uses the following input files. Example input files can be found in the [data](https://github.com/OpenMendel/MendelGeneticCounseling.jl/tree/master/data) subfolder of the Genetic Counseling project. (An analysis won't always need every file type below.)

* [Control File](#control-file): Specifies the names of your data input and output files and any optional parameters (*keywords*) for the analysis. (For a list of common keywords, see [Keywords Table](https://openmendel.github.io/MendelBase.jl/#keywords-table)).
* [Locus File](https://openmendel.github.io/MendelBase.jl/#locus-file): Names and describes the genetic loci in your data.
* [Pedigree File](https://openmendel.github.io/MendelBase.jl/#pedigree-file): Gives information about your individuals, such as name, sex, family structure, and ancestry.
* [Phenotype File](https://openmendel.github.io/MendelBase.jl/#phenotype-file): Lists the available phenotypes.
* [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file): Defines your SNPs with information such as SNP name, chromosome, position, allele names, allele frequencies.
* [SNP Data File](https://openmendel.github.io/MendelBase.jl/#snp-data-file): Holds the genotypes for your data set. Must be a standard binary PLINK BED file in SNP major format. If you have a SNP data file you must have a SNP definition file.

<a id="control-file"></a>
### Control file
The Control file is a text file consisting of keywords and their assigned values. The format of the Control file is:

	Keyword = Keyword_Value(s)

Below is an example of a simple Control file to run Genetic Counseling:

	#
	# Input and Output files.
	#
	locus_file = genetic counseling 2 LocusFrame.txt
	pedigree_file = genetic counseling 2 PedigreeFrame.txt
	phenotype_file = genetic counseling 2 PhenotypeFrame.txt
	output_file = genetic counseling 2 Output.txt
	#
	# Analysis parameters for Genetic Counseling option.
	#
	eliminate_genotypes	= false
	lump_alleles		= false

In the example above, there are six keywords. The first four keywords specify input and output files: *genetic counseling 2 LocusFrame.txt*, *genetic counseling 2 PedigreeFrame.txt*, *genetic counseling 2 PhenotypeFrame.txt*, and *genetic counseling 2 Output.txt*. The last two keywords specify analysis parameters: *eliminate_genotypes* and *lump_alleles*. The text after the '=' are the keyword values. A list of OpenMendel keywords common to most analysis package can be found [here](https://openmendel.github.io/MendelBase.jl/#keywords-table). The names of keywords are *not* case sensitive. (The keyword values *may* be case sensitive.)

### Data Files
Genetic Counseling requires a [Control file](https://openmendel.github.io/MendelBase.jl/#control-file), and a [Pedigree file](https://openmendel.github.io/MendelBase.jl/#pedigree-file). Genotype data can be included in the Pedigree file, in which case a [Locus file](https://openmendel.github.io/MendelBase.jl/#locus-file) is required. Alternatively, genotype data can be provided in a [SNP data file](https://openmendel.github.io/MendelBase.jl/#snp-data-file), in which case a [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file) is required. OpenMendel will also accept [PLINK format](http://zzz.bwh.harvard.edu/plink) FAM and BIM files. Details on the format and contents of the Control and data files can be found on the [MendelBase](https://openmendel.github.io/MendelBase.jl) documentation page. There are example data files in the Genetic Counseling [data](https://github.com/OpenMendel/MendelGeneticCounseling.jl/tree/master/data) folder.

### Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:

     julia> using MendelGeneticCounseling

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")

Finally, to run the analysis using the parameters in your Control file, for example, Control_file.txt, use the command:

     julia> GeneticCounseling("Control_file.txt")

*Note: The package is called* MendelGeneticCounseling *but the analysis function is called simply* GeneticCounseling.

<!--- ### Interpreting the results
 ... --->

### Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ### Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

### Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
