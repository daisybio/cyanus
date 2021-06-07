# CYTOGLMM SIMULATION EVALUATION

source('DEComparison/benchmarking_plots.R')

stats_table <- preparePlotData('simulatedCytoGLMM')
plot_cells_vs_elapsed(stats_table)
plot_sens_vs_pre(stats_table)
plot_sens_vs_spec(stats_table)
plot_f1_vs_elapsed(stats_table)
