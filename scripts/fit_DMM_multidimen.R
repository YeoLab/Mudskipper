library(DirichletMultinomial)
library(lattice)
library(xtable)
library(tidyverse)
library(parallel)

args = commandArgs(trailingOnly=TRUE)

# INPUTs
fl=args[1]
annotation = args[2]
basedir= args[3]
out_stem = args[4] 
seed = args[5]
dir.create(basedir, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(basedir, "/plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(basedir, "/intermediates"), showWarnings = FALSE, recursive = TRUE)

# Set the seed:
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)

options(width=70, digits=2)
full <- TRUE 
.qualitative <- DirichletMultinomial:::.qualitative
dev.off <- function(...) invisible(grDevices::dev.off(...))

min_component = 1
max_component = 40 
component_gap = 5
min_read = 10 

count_df = read_tsv(fl)
sample_cols = colnames(count_df)
count_df$name=row.names(count_df)

print(head(count_df))
count_df = count_df[rowSums(count_df[,sample_cols])>min_read, ]
count <- as.matrix(count_df[, sample_cols]) 
print('count matrix nrows=')
print(nrow(count))

print('count matrix ncol=')
print(ncol(count))

################ MODEL SELECTION ################

# fit data k=1 to max_component, using a subset of data
cores = detectCores()
if (full) {
fit <- mclapply(seq(min_component, max_component, component_gap), dmn, count=count, seed = seed, verbose=TRUE, mc.cores = cores, mc.set.seed = TRUE)
save(fit, file=file.path(paste0(basedir, "/intermediates"), paste0(out_stem, ".fit.rda")))
} else load(file = file.path(paste0(basedir, "/intermediates"), paste0(out_stem, ".fit.rda")))

# plot Laplace against k
lplc <- sapply(fit, laplace)
aic <- sapply(fit, AIC)
bic <- sapply(fit, BIC)
 pdf(file.path(paste0(basedir, "/plots"), paste0(out_stem, ".goodness_of_fit.pdf")))
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
write_tsv(data.frame(weights) %>% rownames_to_column(), file.path(paste(basedir, "/intermediates", sep = ""), paste0(out_stem, '.weights.tsv')))

################ CLUSTER LABELLING ################
# contribution of each taxonomic group to the Dirichlet components fitted 
fitted_df = data.frame(fitted(best, assign = FALSE))%>% rownames_to_column()
write_tsv(fitted_df, file.path(paste0(basedir, "/intermediates"), paste0(out_stem, '.alpha.tsv')))

# write null
fitted_df_null = data.frame(fitted(fit[[1]], assign = FALSE))%>% rownames_to_column()
names(fitted_df_null) <- c('rowname', 'single_component_weight')
write_tsv(fitted_df_null, file.path(paste0(basedir, "/intermediates"), paste0(out_stem, '.null.alpha.tsv')))

# how does the model differ from a single component DMM
p0 <- fitted(fit[[1]], scale=TRUE) # scale by theta
p_best <- fitted(best, scale=TRUE)
colnames(p_best) <- paste("m", 1:ncol(p_best), sep="")
(meandiff <- colSums(abs(p_best - as.vector(p0)))) # the difference of each component to 1 single component

# export the component
mixture_df = mixture(best, assign = FALSE) # sample by component matrix
row.names(mixture_df) <- count_df$name

# join annotation
anno_df = read_tsv(annotation)
annotated_mixture=merge(mixture_df, anno_df, by.x='row.names', by.y='name')

write_tsv(data.frame(annotated_mixture) %>% rownames_to_column(), file.path(paste0(basedir, "/intermediates"), paste0(out_stem,'.mixture_weight.tsv')))

