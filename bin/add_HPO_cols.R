# MIT License

# Copyright (c) 2019, Michaela Mueller, Vicente Yepez, Christian Mertes, Daniela Andrade, Julien Gagneur

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal 
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Function to add HPO terms to results table

add_HPO_cols <- function(RES, sample_id_col = 'sampleID',
                            gene_name_col = 'hgncSymbol', hpo_file = NULL){
    require(data.table)

    filename <- ifelse(is.null(hpo_file),
                    'https://www.cmm.in.tum.de/public/paper/drop_analysis/resource/hpo_genes.tsv.gz',
                    hpo_file)

    hpo_dt <- fread(filename)

    # change column names
    setnames(RES, old = c(sample_id_col, gene_name_col), new = c('sampleID', 'hgncSymbol'))

    # Get HPO terms of all genes
    f2 <- merge(unique(RES[, .(sampleID, hgncSymbol)]),
                hpo_dt[,.(hgncSymbol, HPO_id, HPO_label, Ngenes_per_pheno, Nphenos_gene)],
                by = 'hgncSymbol', allow.cartesian = TRUE)
    sa[, HPO_TERMS := gsub(', ', ',', HPO_TERMS)]

    if(nrow(f2) > 0){
        f3 <- merge(f2, sa[,.(RNA_ID, HPO_TERMS)], by.x = 'sampleID', by.y = 'RNA_ID')
        f3 <- f3[HPO_TERMS != ''] # remove cases with no HPO terms to speed up
        if(nrow(f3) > 0){
            f3[, HPO_match := any(grepl(HPO_id, HPO_TERMS)), by = 1:nrow(f3)]
            f3 <- f3[HPO_match == TRUE] #only take those that match HPOs
            if(nrow(f3) > 0){
                f4 <- f3[, .(HPO_label_overlap = paste(HPO_label, collapse = ', '),
                            HPO_id_overlap = paste(HPO_id, collapse = ', '),
                            Ngenes_per_pheno = paste(Ngenes_per_pheno, collapse = ', '),
                            Nphenos_gene = unique(Nphenos_gene)),
                        by = .(sampleID, hgncSymbol)]
                RES <- merge(RES, f4, by = c('sampleID', 'hgncSymbol'), all.x = TRUE)
            }
        }
    }

    setnames(RES, old = c('sampleID', 'hgncSymbol'), new = c(sample_id_col, gene_name_col))
    return(RES)
}
