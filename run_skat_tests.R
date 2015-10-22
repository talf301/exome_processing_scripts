if(!loaded)
{
    library(SKAT)
    # Assume dir contains the appropriate directory
    variants_path<-paste(dir, "/variants.txt",sep="")
    genos_path<-paste(dir, "/genotype_matrix.txt",sep="")
    weights_path<-paste(dir, "/weights.txt",sep="")
    age_path<-paste(dir, "/ages_per_var.txt",sep="")
    genes_folder<-paste(dir, "/genes/",sep="")

    # Load genotype matrix into a (num_patients, num_variants) matrix, clearly previous step each time
    geno_lines<-readLines(genos_path)
    gt<-lapply(X=geno_lines,FUN=function(x){strsplit(x,',')[[1]]})
    rm(geno_lines)
    gc()
    sink()
    gt2<-lapply(gt, FUN=function(x){replace(x,which(x=="."),"9")})
    rm(gt)
    gc()
    sink()
    genos<-lapply(gt2, as.numeric)
    rm(gt2)
    gc()
    sink()
    geno_mat<-matrix(unlist(genos), length(genos[[1]]), length(genos))
    rm(genos)
    gc()
    sink()

    # Load variant weights and subject ages
    weight_lines<-readLines(weights_path)
    # Add epsilon so weights are non-zero when transformed
    weights<-as.numeric(weight_lines) + 0.001

    age_lines<-readLines(age_path)
    ages<-as.numeric(age_lines)

    loaded<-TRUE
}
print("Loaded!")
# MODIFY THESE LINES TO CHANGE AGE RANGE
case<-which(ages<46)
control<-which(ages>=49)
include<-which(ages<46 | ages>=49)


# Set up SKAT model and prepare genotypes
gts<-geno_mat[include,]
gc(geno_mat)
case_or_control<-as.numeric(lapply(include, FUN=function(x){x %in% case}))
obj<-SKAT_Null_Model(case_or_control ~ 1, out_type="D")

# Load in gene files (ignoring output files), run tests
gene_file_list<-list.files(genes_folder)
gene_file_list<-grep(glob2rx('*.txt'), gene_file_list, value=TRUE)
genes<-list()
pvals<-list()
pvals_opt<-list()
qvals<-list()
for (file in gene_file_list)
{
    # Set up file paths
    gene_file<-paste(genes_folder, file, sep="")
    gene_out_file<-paste(gene_file, ".out", sep="")

    # Get modified info specifically for gene
    gene_variants<-as.numeric(readLines(gene_file))
    # Get remaining variants after filtering for rare
    filtered_variants<-Filter(function(x) {sum(gts[,x] %% 9) <= 0.05 * 2 * length(which(gts[,x] < 9))}, gene_variants)

    gene_gts<-as.matrix(gts[,filtered_variants])
    gene_weights<-weights[filtered_variants]
    genes<-c(genes, file)
    sink(gene_out_file, append=TRUE, split=FALSE)
    if (length(filtered_variants) > 0)
    {
        pval<-SKAT(gene_gts, obj, weights=gene_weights)$p.value
        pval_opt<-SKAT(gene_gts, obj, method="optimal.adj", weights=gene_weights)$p.value
        qval<-SKAT(gene_gts, obj, weights=gene_weights)$Q
    }
    else
    {
        pval<-NA
        pval_opt<-NA
        qval<-NA
    }
    sink()
    pvals<-c(pvals, pval)
    pvals_opt<-c(pvals_opt, pval_opt)
    qvals<-c(qvals, qval)

    if(length(genes) %% 100 == 0)
    {
        print(length(genes))
    }
}
results_file<-paste(dir, res_name, sep="/")
m<-matrix(c(unlist(genes), unlist(pvals), unlist(pvals_opt), unlist(qvals)), nrow=length(genes), ncol=4)
to_write<-data.frame(m)
write.table(to_write, sep="\t", row.names=FALSE, col.names=FALSE, file=results_file, quote=FALSE)
