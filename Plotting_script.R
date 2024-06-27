# Plotting script for the figures from All_species_mono_phylo_islands and Cory_mono_phylo_islands
# Run this script after you have run both All_species_mono_phylo_islands and Cory_mono_phylo_islands.rmd

#Setting the working directory
setwd("/home/au543206/Documents/Coryphoideae/Figures_in_r/Island_counts")

# Plotting coefficients table
pdf(file = file.path("./Figures","Coefficients_bars.pdf"),
    width = 7,
    height = 8)

plot(prow_3)

dev.off()

# Plotting the scatterplot of speciation for all species
prow1 <- cowplot::plot_grid(cory_scatter, all_vasc_scatter+theme(legend.position = "none"), labels = c("a. ", "b."), align = "h", rel_widths = c(1,1), ncol = 1, nrow = 2)
legend_side_1 <- get_legend(p2)

pdf(file = file.path("./Figures","scatter_plots_speciation_all_sp_facet.pdf"),
    width = 10,
    height = 7)

cowplot::plot_grid(prow1,legend_side_1, ncol=2,nrow = 1, rel_widths = c(1,.2))

dev.off()


# Plotting coefficients table for only volcanic islands
pdf(file = file.path("./Figures","Coefficients_bars_volcanic.pdf"),
    width = 7,
    height = 4)

plot(prow_volc)

dev.off()

# Plotting coefficients table
pdf(file = file.path("./Figures","Coefficients_bars_all_square.pdf"),
    width = 8,
    height = 6.6)

plot(prow_5)

dev.off()


# Plotting Anagenesis hump shape
pdf(file = file.path("./Figures","Anagenesis_hump.pdf"),
    width = 8,
    height = 6.6)

plot(plot_col_hump)

dev.off()


# 
# # Plotting idea 1 for all phylogenies
# pdf(file = file.path("./Figures","Phylogenies_draft.pdf"),
#     width = 14,
#     height = 8)
# 
# plot(phylo4)
# 
# dev.off()
# 
pdf(file = file.path("./Figures","Speciation_of_coryphoideae_on_islands.pdf"),
    width = 14,
    height = 8)

plot(cory_phylo_dat_plot)

dev.off()


