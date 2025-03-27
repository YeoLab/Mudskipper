library(DirichletMultinomial)
library(lattice)
library(xtable)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# INPUTs
fl=args[1] 
annotation = args[2] 
ip_col = args[3]
in_col = args[4]
basedir= args[5]
out_stem = args[6]
dir.create(basedir, showWarnings = FALSE, recursive = TRUE)

sample_cols = c(ip_col, in_col)
print(sample_cols)

# Parameters
options(width=70, digits=2)
full <- TRUE 
.qualitative <- DirichletMultinomial:::.qualitative
dev.off <- function(...) invisible(grDevices::dev.off(...))
min_component = 1
max_component = 5
component_gap = 1
comp_attempt = seq(min_component, max_component, component_gap)
min_read = 10 # IgG had only 1 component when min_read = 0, CC was fine
fr_over_total = 6

# filter read count
count_df = read_tsv(fl)
print(head(count_df))


count_df = count_df[rowSums(count_df[sample_cols])>min_read, ]
# select by average fraction of reads instead of total reads. This might benefit lowly sequenced RBP
# average_fraction <- rowSums(t(t(count_df[sample_cols])/colSums(count_df[sample_cols])))/length(sample_cols)
# count_df = count_df[average_fraction>fr_over_total/nrow(count_df), ]

count <- as.matrix(count_df[sample_cols]) # don't need to t because we have different orientation
print('count matrix nrows=')
print(nrow(count))

print('count matrix ncol=')
print(ncol(count))

# fit data k=1 to max_component, using a subset of data
library(parallel)
cores = min(detectCores(), length(comp_attempt))
if (full) {
fit <- mclapply(comp_attempt, dmn, count=count, verbose=TRUE, mc.cores = cores)
save(fit, file=file.path(basedir, paste0(out_stem, ".fit.rda")))
} else load(file = file.path(basedir, paste0(out_stem, ".fit.rda")))

################ MODEL SELECTION ################



# plot Laplace against k
lplc <- sapply(fit, laplace)
aic <- sapply(fit, AIC)
bic <- sapply(fit, BIC)
 pdf(file.path(basedir, paste0(out_stem, ".goodness_of_fit.pdf")))
 plot(aic, type="b", xlab="Number of Dirichlet Components(k)",ylab="AIC")
 plot(bic, type="b", xlab="Number of Dirichlet Components(k)",ylab="BIC")
 plot(lplc, type="b", xlab="Number of Dirichlet Components(k)",ylab="Model Fit(Laplace)")
 dev.off()

  # find the best model: the DMM object
 (best <- fit[[which.min(bic)]])

################ CLUSTER SIZE ################
# reports the weight $\pi$ and $\theta$
# theta = \sum alphas, higher, more concentrated cluster
weights = mixturewt(best)
write_tsv(data.frame(weights) %>% rownames_to_column(), file.path(basedir, paste0(out_stem, '.weights.tsv')))

################ CLUSTER LABELLING ################
# contribution of each taxonomic group to the Dirichlet components fitted
fitted_df = data.frame(fitted(best, assign = FALSE))%>% rownames_to_column()
write_tsv(fitted_df, file.path(basedir, paste0(out_stem, '.alpha.tsv')))

# write null
fitted_df_null = data.frame(fitted(fit[[1]], assign = FALSE))%>% rownames_to_column()
names(fitted_df_null) <- c('rowname', 'single_component_weight')
write_tsv(fitted_df_null, file.path(basedir, paste0(out_stem, '.null.alpha.tsv')))

# how does the model differ from a single component DMM
p0 <- fitted(fit[[1]], scale=TRUE) # scale by theta
p_best <- fitted(best, scale=TRUE)
colnames(p_best) <- paste("m", 1:ncol(p_best), sep="")
(meandiff <- colSums(abs(p_best - as.vector(p0)))) # the difference of each component to 1 single component (possibly the null)

# export the component
mixture_df = mixture(best, assign = FALSE) # sample by component matrix
row.names(mixture_df) <- count_df$name

# join annotation
anno_df = read_tsv(annotation)
annotated_mixture=merge(mixture_df, anno_df, by.x='row.names', by.y='name')

write_tsv(data.frame(annotated_mixture) %>% rownames_to_column(), file.path(basedir, paste0(out_stem,'.mixture_weight.tsv')))

