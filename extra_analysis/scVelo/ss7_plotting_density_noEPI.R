suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

meta = read.csv('saved_object/meta_ss5_latent_time_paga_dpt_noEPI.csv',row.names=1)

meta$cell_cluster = as.factor(meta$cell_cluster)

#######################

plot_1 = ggplot(meta) + 
geom_density(aes(x = velocity_pseudotime, y=..count.., color = cell_cluster)) + 
geom_rug(aes(x = velocity_pseudotime, color = cell_cluster))+ 
theme_bw()+
ggtitle('scvelo velocity pseudotime density and rug plot noEPI')+ 
theme(plot.title = element_text(size=15),axis.title = element_text(size = 15)) +
xlab("scvelo velocity pseudotime values")

plot_2 = ggplot(meta) + 
geom_density(aes(x = latent_time, y=..count.., color = cell_cluster)) + 
geom_rug(aes(x = latent_time, color = cell_cluster))+ 
theme_bw()+
ggtitle('scvelo latent time density and rug plot noEPI')+ 
theme(plot.title = element_text(size=15),axis.title = element_text(size = 15)) +
xlab("scvelo latent time values")

plot_3 = ggplot(meta) + 
geom_density(aes(x = dpt_pseudotime, y=..count.., color = cell_cluster)) + 
geom_rug(aes(x = dpt_pseudotime, color = cell_cluster))+ 
theme_bw()+
ggtitle('scvelo dpt pseudotime density and rug plot noEPI')+ 
theme(plot.title = element_text(size=15),axis.title = element_text(size = 15)) +
xlab("scvelo dpt pseudotime values")

######################

pdf('saved_plot/ss7_latent_dpt_pseudotime_density_scvelo_noEPI.pdf',width = 11, height = 15)
plot_grid(
  plot_1, plot_2, plot_3,
  labels = "AUTO", ncol = 1
)

dev.off()