# asd-development-coexpression-2013

This code was originally posted at:
https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/developingcortex/

Manuscript
Parikshak, N. N., Luo, R., Zhang, A., Won, H., Lowe, J. K., Chandran, V., Horvath S., Geschwind D.H. (2013). Integrative Functional Genomic Analyses Implicate Specific Molecular Pathways and Circuits in Autism. Cell, 155(5), 1008–1021. doi:10.1016/j.cell.2013.10.031 PMCID: PMC3107252

Short Description
Genetic studies have identified dozens of autism spectrum disorder (ASD) susceptibility genes, raising two critical questions:
1) do these genetic loci converge on specific biological processes, and
2) where does the phenotypic specificity of ASD arise, given its genetic overlap with intellectual disability (ID)?
To address this, we mapped ASD and ID risk genes onto co-expression networks representing developmental trajectories. This code will re-construct the developing cortex network in the manuscript, and run the gene set enrichments.

This code will run the developmental network analysis if used with the related .Rdata file. The code was written by Neelroop Parikshak (neelroop [-at-] ucla [-dot-] edu).
The code creates the relevant folders and runs the entire analysis up to creating tables of the network membership and gene set enrichments.
It is important to run this code on a system with enough RAM. The topological overlap matrix computation step with about 22k genes takes up to 16G of RAM. One can set the blockSize parameter to 10000 or 5000 and obtain similar output, but it is most accurate to run the entire network as one block (current parameter is set to 25000). This also allows easy visualization of the topological overlap dendrogram.

Input
Load the following data file into R (e.g. by using load("AllSupData-11-1-2013.Rdata") or by double clicking the file: AllSupData-11-1-2013.Rdata

Output
This code outputs three folders:
./figures
./tables
./data
Each containing analysis output as delineated in the code below.

 
