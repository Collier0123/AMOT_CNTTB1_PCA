library(FactoMineR)
library(ggplot2)
library(factoextra)
library(RColorBrewer)
library(dplyr)
# Data preparation
Apical_Lateral_sum_median$Group <- as.factor(Apical_Lateral_sum_median$Group)
Apical_Lateral_sum_median$Domain <- as.factor(Apical_Lateral_sum_median$Domain)
Apical_Lateral_sum_median$Inj.Status <- as.factor(Apical_Lateral_sum_median$Inj.Status)
Apical_Lateral_sum_median$Condition <- interaction(Apical_Lateral_sum_median$Domain, Apical_Lateral_sum_median$Inj.Status)
Apical_Lateral_sum_median$Condition <- factor(Apical_Lateral_sum_median$Condition)
Apical_Lateral_sum_median$GroupCondition <- interaction(Apical_Lateral_sum_median$Group, Apical_Lateral_sum_median$Domain, Apical_Lateral_sum_median$Inj.Status)
Apical_Lateral_sum_median$GroupCondition <- factor(Apical_Lateral_sum_median$GroupCondition, 
                                                   levels = c("siNTC.Apical.Non.inj.", "siNTC.Apical.Inj.", 
                                                              "siNTC.Lateral.Non.inj.", "siNTC.Lateral.Inj.",
                                                              "siTead4.Apical.Non.inj.", "siTead4.Apical.Inj.",
                                                              "siTead4.Lateral.Non.inj.", "siTead4.Lateral.Inj."))
# Perform PCA
pca_data_domain <- Apical_Lateral_sum_median[, c("Sum.AMOT","Length", "Sum.B.Cat","Median.Ratio")]
pca_result_domain <- PCA(pca_data_domain, ncp = 5, graph = FALSE)


#pre-plot analysis
library(corrplot)
pca_result_domain$eig
fviz_eig(pca_result_domain, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_var(pca_result_domain, repel = TRUE)
fviz_contrib(pca_result_domain,"var")
fviz_contrib(pca_result_domain,"var", axes = 2)
fviz_contrib(pca_result_domain,"var", axes = 2)
fviz_contrib(pca_result_domain,"var", axes = 3)
fviz_contrib(pca_result_domain, "var", axes = 1:2)
fviz_contrib(pca_result_domain, "var", axes = 1:3)
corrplot(pca_result_domain$var$coord)

#Enhanced feature scree-plot
scree_plot <- fviz_eig(pca_result_domain, addlabels = FALSE, barfill = "black", barcolor = "black") +
  labs(
    title = "Scree Plot",         
    x = "Dimensions",             
    y = "Percentage of explained 
    varience (%)"  
  ) +
  theme_minimal() +                         
  theme(
    plot.title = element_text(size = 45,face = "bold", hjust = 0.5, margin = margin(r=60,b =50)),
    plot.subtitle = element_text(size=40, hjust =0.5, margin = margin(b=50, r=60)),
    axis.title = element_text(size = 40, hjust = 0.5),                           
    axis.text = element_text(size = 33), 
    panel.grid.major = element_line(color = "gray"), 
    panel.grid.minor = element_blank(),                              
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 20))
  )


print(scree_plot)
write.csv(eig_values, "eigenvalues.csv", row.names = FALSE)
#Values for enhanced feature graphs
#dim vaues
eig_values <- pca_result_domain$eig
dim1_var <- format(round(eig_values[1, 2], 1), nsmall = 1)
dim2_var <- format(round(eig_values[2, 2], 1), nsmall = 1)
dim3_var <- format(round(eig_values[3, 2], 1), nsmall = 1)
dim4_var <- format(round(eig_values[4, 2], 1), nsmall = 1)

#Row lables - dynamic
variable_name_mapping <- c(
  "Sum.AMOT" = "AMOT TCCF",
  "Length" = "Length 
  (µm)",
  "Sum.B.Cat" = "CTNNB1 TCCF",
  "Median.Ratio" = "Median Ratio 
  (AMOT/CTNNB1)"
)

raw_variable_names <- rownames(pca_result_domain$var$coord)
domain_labels_row <- variable_name_mapping[raw_variable_names]

#Column lables - dynamic
domain_correlation_matrix <- pca_result_domain$var$coord
domain_labels_column <- paste0("Dim", 1:ncol(domain_correlation_matrix))


#Matrix mapping and final features
rownames(domain_correlation_matrix) <- domain_labels_row
colnames(domain_correlation_matrix) <- domain_labels_column

#Enhanced corr plot (CHECK THAT TITLES MATCH BASIC!)
mwb_palette <- colorRampPalette(c("magenta","white", "black"))
corrplot(domain_correlation_matrix, 
         col = mwb_palette(200),          
         method = "circle",               
         tl.col = "black",                 
         tl.offset = 0.9,
         tl.cex = 1.5,                     
         tl.srt = 0,                   
         cl.cex = 1,
         cl.pos = "r")   

#Contribution plots
#Dim1
fviz_contrib_plot1 <- fviz_contrib(pca_result_domain, "var", axes = 1)
fviz_contrib_data1 <- fviz_contrib_plot1$data
fviz_contrib_data1$variable <- variable_name_mapping[rownames(fviz_contrib_data1)]
fviz_contrib_data1$variable <- reorder(fviz_contrib_data1$variable, -fviz_contrib_data1$contrib)
ggplot(fviz_contrib_data1, aes(x = variable, y = contrib)) +
  geom_bar(stat = "identity", color="black", fill = "black", width = 0.3) +
  geom_hline(yintercept = mean(fviz_contrib_data1$contrib), linetype = "dashed", color = "red") +
  labs(
    title = "Variable Contributions to Dimension 1",
    x = NULL,
    y = "Contribution (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 26,face = "bold", hjust = 0.5, margin = margin(r=60, b=20)),
    plot.subtitle = element_text(size=26, hjust =0.5, margin = margin(b=50, r=60)),
    axis.title = element_text(size = 22, hjust = 0.5),                           
    axis.text = element_text(size = 20), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20, color = "black"),
    panel.grid.major = element_line(color = "gray"), 
    panel.grid.minor = element_blank(),                              
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 20))
  )
fviz_contrib(pca_result_domain,"var")

#Dim2
fviz_contrib_plot2 <- fviz_contrib(pca_result_domain, "var", axes = 2)
fviz_contrib_data2 <- fviz_contrib_plot2$data
fviz_contrib_data2$variable <- variable_name_mapping[rownames(fviz_contrib_data2)]
fviz_contrib_data2$variable <- reorder(fviz_contrib_data2$variable, -fviz_contrib_data2$contrib)
ggplot(fviz_contrib_data2, aes(x = variable, y = contrib)) +
  geom_bar(stat = "identity", color="black", fill = "black", width = 0.3) +
  geom_hline(yintercept = mean(fviz_contrib_data2$contrib), linetype = "dashed", color = "red") +
  labs(
    title = "Variable Contributions to Dimension 2",
    x = NULL,
    y = "Contribution (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 26,face = "bold", hjust = 0.5, margin = margin(r=60, b=20)),
    plot.subtitle = element_text(size=26, hjust =0.5, margin = margin(b=50, r=60)),
    axis.title = element_text(size = 22, hjust = 0.5),                           
    axis.text = element_text(size = 20), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20, color = "black"),
    panel.grid.major = element_line(color = "gray"), 
    panel.grid.minor = element_blank(),                              
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 20))
  )
fviz_contrib(pca_result_domain,"var", axes = 2)

#Dim3
fviz_contrib_plot3 <- fviz_contrib(pca_result_domain, "var", axes = 3)
fviz_contrib_data3 <- fviz_contrib_plot3$data
fviz_contrib_data3$variable <- variable_name_mapping[rownames(fviz_contrib_data3)]
fviz_contrib_data3$variable <- reorder(fviz_contrib_data3$variable, -fviz_contrib_data3$contrib)
ggplot(fviz_contrib_data3, aes(x = variable, y = contrib)) +
  geom_bar(stat = "identity", color="black", fill = "black", width = 0.3) +
  geom_hline(yintercept = mean(fviz_contrib_data3$contrib), linetype = "dashed", color = "red") +
  labs(
    title = "Variable Contributions to Dimension 3",
    x = NULL,
    y = "Contribution (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 26,face = "bold", hjust = 0.5, margin = margin(r=60, b=20)),
    plot.subtitle = element_text(size=26, hjust =0.5, margin = margin(b=50, r=60)),
    axis.title = element_text(size = 22, hjust = 0.5),                           
    axis.text = element_text(size = 20), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20, color = "black"),
    panel.grid.major = element_line(color = "gray"), 
    panel.grid.minor = element_blank(),                              
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 20))
  )
fviz_contrib(pca_result_domain, "var", axes = 3)

#Dim1-2
fviz_contrib_plot4 <- fviz_contrib(pca_result_domain, "var", axes = 1:2)
fviz_contrib_data4 <- fviz_contrib_plot4$data
fviz_contrib_data4$variable <- variable_name_mapping[rownames(fviz_contrib_data4)]
fviz_contrib_data4$variable <- reorder(fviz_contrib_data4$variable, -fviz_contrib_data4$contrib)
ggplot(fviz_contrib_data4, aes(x = variable, y = contrib)) +
  geom_bar(stat = "identity", color="black", fill = "black", width = 0.3) +
  geom_hline(yintercept = mean(fviz_contrib_data4$contrib), linetype = "dashed", color = "red") +
  labs(
    title = "Variable Contributions to Dimensions 1-2",
    x = NULL,
    y = "Contribution (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 26,face = "bold", hjust = 0.5, margin = margin(r=60, b=20)),
    plot.subtitle = element_text(size=26, hjust =0.5, margin = margin(b=50, r=60)),
    axis.title = element_text(size = 22, hjust = 0.5),                           
    axis.text = element_text(size = 20), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20, color = "black"),
    panel.grid.major = element_line(color = "gray"), 
    panel.grid.minor = element_blank(),                              
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 20))
  )
fviz_contrib(pca_result_domain, "var", axes = 1:2)

#Dim1-3
fviz_contrib_plot5 <- fviz_contrib(pca_result_domain, "var", axes = 1:3)
fviz_contrib_data5 <- fviz_contrib_plot5$data
fviz_contrib_data5$variable <- variable_name_mapping[rownames(fviz_contrib_data5)]
fviz_contrib_data5$variable <- reorder(fviz_contrib_data5$variable, -fviz_contrib_data5$contrib)
ggplot(fviz_contrib_data5, aes(x = variable, y = contrib)) +
  geom_bar(stat = "identity", color="black", fill = "black", width = 0.3) +
  geom_hline(yintercept = mean(fviz_contrib_data5$contrib), linetype = "dashed", color = "red") +
  labs(
    title = "Variable Contributions to Dimensions 1-2-3",
    x = NULL,
    y = "Contribution (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 26,face = "bold", hjust = 0.5, margin = margin(r=60, b=20)),
    plot.subtitle = element_text(size=26, hjust =0.5, margin = margin(b=50, r=60)),
    axis.title = element_text(size = 22, hjust = 0.5),                           
    axis.text = element_text(size = 20), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20, color = "black"),
    panel.grid.major = element_line(color = "gray"), 
    panel.grid.minor = element_blank(),                              
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 20))
  )
fviz_contrib(pca_result_domain, "var", axes = 1:3)


#Create levels
ordered_levels_norm <- c(sort(unique(Apical_Lateral_sum_median$GroupCondition[grep("siNTC", Apical_Lateral_sum_median$GroupCondition)])),
                         sort(unique(Apical_Lateral_sum_median$GroupCondition[grep("siTead4", Apical_Lateral_sum_median$GroupCondition)])))
Apical_Lateral_sum_median$GroupCondition <- factor(Apical_Lateral_sum_median$GroupCondition, levels = ordered_levels_norm)

# Define the specific color scheme
condition_colours <- c(
  "siNTC.Apical.Non.inj." = "#B2BEB5",   
  "siNTC.Apical.Inj." = "black",       
  "siNTC.Lateral.Non.inj." = "#8fd9fb",   
  "siNTC.Lateral.Inj." = "#00AFFF",       
  "siTead4.Apical.Non.inj." = "#FF0000", 
  "siTead4.Apical.Inj." = "#880808",     
  "siTead4.Lateral.Non.inj." = "#ffb343", 
  "siTead4.Lateral.Inj." = "#fe8205"      
)

#Legend labels
group_condition_labels <- c(
  "siNTC.Apical.Non.inj." = "siNTC Non-inj.: Apical domain",
  "siNTC.Apical.Inj." = "siNTC Inj.: Apical domain",
  "siNTC.Lateral.Non.inj." = "siNTC Non-inj.: Lateral domain",
  "siNTC.Lateral.Inj." = "siNTC Inj.: Lateral domain",
  "siTead4.Apical.Non.inj." = "siTead4 Non-inj.: Apical domain",
  "siTead4.Apical.Inj." = "siTead4 Inj.: Apical domain",
  "siTead4.Lateral.Non.inj." = "siTead4 Non-inj.: Lateral domain",
  "siTead4.Lateral.Inj." = "siTead4 Inj.: Lateral domain"
)


#2D 4-axis plot values
library(ggplot2)
ind_coords_Apical_Lateral_sum_median <- as.data.frame(pca_result_domain$ind$coord)
plot_data_Apical_Lateral_sum_median <- data.frame(
  Dim1 = ind_coords_Apical_Lateral_sum_median[, 1],
  Dim2 = ind_coords_Apical_Lateral_sum_median[, 2],
  Dim3 = ind_coords_Apical_Lateral_sum_median[, 3],
  Dim4 = ind_coords_Apical_Lateral_sum_median[, 4],
  GroupCondition = Apical_Lateral_sum_median$GroupCondition
)

centroid_coords <- as.data.frame(pca_result_domain$ind$coord) %>%
  dplyr::mutate(GroupCondition = Apical_Lateral_sum_median$GroupCondition) %>%
  dplyr::group_by(GroupCondition) %>%
  dplyr::summarize(
    Dim.1 = mean(Dim.1),
    Dim.2 = mean(Dim.2)
  )

#Basic plot data - NEED THIS TO AQCUIRE CORRECT ELLIPSES, make sure to use correct habillage
fviz_pca_plot <- fviz_pca_ind(pca_result_domain,
                              habillage = Apical_Lateral_sum_median$GroupCondition,  
                              palette = condition_colours,
                              addEllipses = TRUE,
                              ellipse.level = 0.95,
                              ellipse.alpha = 0,
                              geom = "point",
                              repel = FALSE
) +
  scale_shape_manual(values = rep(19, length(unique(Apical_Lateral_sum_median$GroupCondition)))) +
  guides(shape = "none")

#2Dims
fviz_pca_plot +
  geom_hline(yintercept = 0, linetype = "solid", color = "#C0C0C0") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "#C0C0C0") +
  geom_point(data = centroid_coords, aes(x = Dim.1, y = Dim.2),color = "black", fill = condition_colours, 
             size = 3, shape = 24, stroke = 0.5) +
  scale_color_manual(values = condition_colours, name = "Cell Profile",
                     labels = group_condition_labels[levels(Apical_Lateral_sum_median$GroupCondition)]
  ) +
  labs(
    x = paste0("Dim1 (", dim1_var, "%)"),
    y = paste0("Dim2 (", dim2_var, "%)"),
    title = "PCA of Apical and Lateral domains"
  ) +
  guides(
    color = guide_legend(title = "Cell Profile"),
    size = guide_legend(title = paste0("Dim3 (", dim3_var, "%)")),
    fill = "none"  
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, margin = margin(b=20, l=100)),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    plot.margin = margin(0, 0, 0, 0))
