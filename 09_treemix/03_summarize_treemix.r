# the plotting funcs script comes with treemix
source("plotting_funcs.R")

# get proportion variance explained by base tree and trees with migration edges
x <- get_f("certhia_contact_10000.treemix")
x <- c(x, get_f("certhia_contact_10000_m1.treemix"))
x <- c(x, get_f("certhia_contact_10000_m2.treemix"))
x <- c(x, get_f("certhia_contact_10000_m3.treemix"))
plot(x)


plot_tree("certhia_contact_10000.treemix")
plot_tree("certhia_contact_10000_m1.treemix")

plot_tree("certhia_contact_10000_m2.treemix")

