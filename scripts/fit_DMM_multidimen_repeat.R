library(DirichletMultinomial)
library(lattice)
library(xtable)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

fl=args[1]
basedir= args[2]
out_stem = args[3] 
dir.create(basedir, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(basedir, "/plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(basedir, "/intermediates"), showWarnings = FALSE, recursive = TRUE)

options(width=70, digits=2)
full <- TRUE
.qualitative <- DirichletMultinomial:::.qualitative
dev.off <- function(...) invisible(grDevices::dev.off(...))

min_component = 1
max_component = 21 
component_gap = 4
min_read = 10 

# load counts
count_df = read_tsv(fl)
sample_cols = colnames(count_df)[-1]
name_col = colnames(count_df)[1]
count_df['name']=count_df[name_col]

# filter read 
print(head(count_df))
count_df = count_df[rowSums(count_df[,sample_cols])>min_read, ]
count <- as.matrix(count_df[, sample_cols]) 
print('count matrix nrows=')
print(nrow(count))

print('count matrix ncol=')
print(ncol(count))

# fit
library(parallel)
cores = detectCores()
if (full) {
fit <- mclapply(seq(min_component, max_component, component_gap), dmn, count=count, verbose=TRUE, mc.cores = cores)
save(fit, file=file.path(paste0(basedir, "/intermediates"), paste0(out_stem, ".fit.rda")))
} else load(file = file.path(paste0(basedir, "/intermediates"), paste0(out_stem, ".fit.rda")))

# find best fit
lplc <- sapply(fit, laplace)
aic <- sapply(fit, AIC)
bic <- sapply(fit, BIC)
 pdf(file.path(paste0(basedir, "/plots"), paste0(out_stem, ".goodness_of_fit.pdf")))
 plot(aic, type="b", xlab="Number of Dirichlet Components(k)",ylab="AIC")
 plot(bic, type="b", xlab="Number of Dirichlet Components(k)",ylab="BIC")
 plot(lplc, type="b", xlab="Number of Dirichlet Components(k)",ylab="Model Fit(Laplace)")
 dev.off()
(best <- fit[[which.min(bic)]])

# mixture weights
weights = mixturewt(best)
write_tsv(data.frame(weights) %>% rownames_to_column(), file.path(paste0(basedir, "/intermediates"), paste0(out_stem, '.weights.tsv')))

# parameters
fitted_df = data.frame(fitted(best, assign = FALSE))%>% rownames_to_column()
write_tsv(fitted_df, file.path(paste0(basedir, "/intermediates"), paste0(out_stem, '.alpha.tsv')))

# single component parameters
fitted_df_null = data.frame(fitted(fit[[1]], assign = FALSE))%>% rownames_to_column()
names(fitted_df_null) <- c('rowname', 'single_component_weight')
write_tsv(fitted_df_null, file.path(paste0(basedir, "/intermediates"), paste0(out_stem, '.null.alpha.tsv')))

# which component is most different from the mean
p0 <- fitted(fit[[1]], scale=TRUE) # scale by theta
p_best <- fitted(best, scale=TRUE)
colnames(p_best) <- paste("m", 1:ncol(p_best), sep="")
(meandiff <- colSums(abs(p_best - as.vector(p0)))) # the difference of each component to 1 single component.

# posterior
mixture_df = mixture(best, assign = FALSE) # sample by component matrix
row.names(mixture_df) <- count_df$name
write_tsv(data.frame(mixture_df) %>% rownames_to_column(), file.path(paste0(basedir, "/intermediates"), paste0(out_stem,'.mixture_weight.tsv')))
