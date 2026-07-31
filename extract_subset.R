
load('/Users/hellyr/Documents/ONPRCliverSitePercMeth5x1090n94.RData')
sd_vals <- apply(perc5X, 1, sd, na.rm=TRUE)
top_sites <- names(sort(sd_vals, decreasing=TRUE)[1:200])
sub_mat <- perc5X[top_sites, ]
write.csv(sub_mat, 'liver_subset.csv')
