# Plotting script for the figures from All_species_mono_phylo_islands and Cory_mono_phylo_islands
# Run this script after you have run both All_species_mono_phylo_islands and Cory_mono_phylo_islands.rmd

# Plotting coefficients table
pdf(file = file.path(output_folder,"Coefficients_bars.pdf"),
    width = 7,
    height = 8)

plot(prow_3)

dev.off()

# Plotting the scatterplot of speciation for all species
prow1 <- cowplot::plot_grid(cory_scatter, all_vasc_scatter+theme(legend.position = "none"), labels = c("a. ", "b."), align = "h", rel_widths = c(1,1), ncol = 1, nrow = 2)
legend_side_1 <- get_legend(p2)


pdf(file = file.path(output_folder,"scatter_plots_speciation_all_sp_facet.pdf"),
    width = 10,
    height = 7)

cowplot::plot_grid(prow1,legend_side_1, ncol=2,nrow = 1, rel_widths = c(1,.2))

dev.off()


# Plotting coefficients table for only volcanic islands
pdf(file = file.path(output_folder,"Coefficients_bars_volcanic.pdf"),
    width = 7,
    height = 4)

plot(prow_volc)

dev.off()

# Plotting coefficients table
pdf(file = file.path(output_folder,"Coefficients_bars_all_square.pdf"),
    width = 8,
    height = 6.6)

plot(prow_5)

dev.off()


# Plotting Anagenesis hump shape
pdf(file = file.path(output_folder,"Anagenesis_hump.pdf"),
    width = 8,
    height = 6.6)

plot(plot_col_hump)

dev.off()


# Plotting the phylogeny of Coryphoideae
pdf(file = file.path(output_folder,"Speciation_of_coryphoideae_on_islands.pdf"),
    width = 14,
    height = 8)

plot(cory_phylo_dat_plot)

dev.off()




# combine cory boxplot and all vasc boxplot using cowplot
prow2 <- cowplot::plot_grid(cory_boxplot, all_vasc_boxplot, labels = c("a. ", "b."), align = "h", rel_widths = c(1,1), ncol = 1, nrow = 2)
prow2

# and then adding the boxplot of the physical variables
prow3 <- cowplot::plot_grid(prow2, boxplot_phys_vars_with_legend, labels = c("a. ", "c."), align = "h", rel_widths = c(1,1), ncol = 1, nrow = 2)

pdf(file = file.path(output_folder,"Boxplots_all_species.pdf"),
    width = 12,
    height = 12)

plot(prow3)

dev.off()
