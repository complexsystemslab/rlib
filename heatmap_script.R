##########################################
# Function to generate heatmaps
##########################################


#save(m_Metagene_metaheatmap_all_signatures_workdir_ATNF_curated_onlyatnf, file = "heatmap_file_microarray_intersect_scseq.RData")
save(df_column_pheatmap2, file="df_column_pheatmap2.RData")
save(ann_colors, file = "ann_colors.RData")

##########################################
# Load data
##########################################
load(file = "heatmap_file_microarray_intersect_scseq.RData")
load(file = "df_column_pheatmap2.RData")
load(file = "ann_colors.RData")


##########################################
# Scale matrix
##########################################
s_m_Metagene_metaheatmap_all_signatures_workdir_ATNF_curated_onlyatnf <- t(scale(t(m_Metagene_metaheatmap_all_signatures_workdir_ATNF_curated_onlyatnf)))

#low  <- quantile(abs(s_counts_filtered_df_curated),0.05, na.rm = TRUE)
high <- quantile(abs(s_m_Metagene_metaheatmap_all_signatures_workdir_ATNF_curated_onlyatnf), 0.95, na.rm = TRUE)

ncolors <- 500
colors <- colorRampPalette(c("blue","white","red"))(ncolors)
breaks <- seq(-high, high, length.out = ncolors + 1)
# colors_function <- colorRamp2(
#   breaks = seq(-high, high, length.out = ncolors),
#   colors = colors)

##################################    
# Reordered rows function
##################################
fn_reorder_dendrogram <- function(m_matrix,
                                  str_method_distance_metric, 
                                  str_method_clustering,
                                  p_minkowski)
{
          ######################################################################
          # generic dendrogram reordering function
          #   m_matrix is NOT SCALED or UNSCALED matrix
          # 
          #  fn_reorder_dendrogram(m_matrix, "manhattan", "average", 1.2) 
          #  fn_reorder_dendrogram(m_matrix, "minkowski", "average", 1.2) 
          # 
          # returns the reordered hclust object
          #
          ######################################################################
          
          require(cba)
          
          if (str_method_distance_metric == "minkowski")
          {
            data_dist <- dist(m_matrix, method="minkowski", p=p_minkowski)
          }
          else
          {
            data_dist <- dist(m_matrix, method=str_method_distance_metric)
          }
          #data_dist <- dist(t(m_to_be_plotted_withbar_subset_active_CDnon_UC_control_curated_onlyatnf), method="manhattan")#, p=1.44)
          
          # data_hclust <- hclust(data_dist, method="average")
          data_hclust <- hclust(data_dist, method=str_method_clustering)
          
          # plot the dendrogram without leaf ordering
          den <- as.dendrogram(data_hclust)
          #par(mfrow=c(1,2))
          #par(mar=c(20,2,2,2))
          #plot(den, main = "without leaf ordering")
          
          # calculate the optimal leaf order
          co <- order.optimal(data_dist, data_hclust$merge)
          
          # overwrite the hclust object with the optimal leaf order
          data_hclust$merge <- co$merge
          data_hclust$order <- co$order
          
          # plot the dendrogram with leaf ordering
          den_opt <- as.dendrogram(data_hclust)
          
          #return(den_opt)
          return(data_hclust)
  
}





hclust_row_m_Metagene_metaheatmap_all_signatures_workdir_ATNF_curated_onlyatnf = 
  fn_reorder_dendrogram(m_matrix=(s_m_Metagene_metaheatmap_all_signatures_workdir_ATNF_curated_onlyatnf),
                        str_method_distance_metric = "minkowski",
                        str_method_clustering = "average",
                        p_minkowski = 1.2) 


hclust_col_m_Metagene_metaheatmap_all_signatures_workdir_ATNF_curated_onlyatnf = 
  fn_reorder_dendrogram(m_matrix=(t(m_Metagene_metaheatmap_all_signatures_workdir_ATNF_curated_onlyatnf)),
                        str_method_distance_metric = "minkowski",
                        str_method_clustering = "average",
                        p_minkowski = 1.2) 


pheatmap(s_m_Metagene_metaheatmap_all_signatures_workdir_ATNF_curated_onlyatnf,
         #color=colorRampPalette(rev(brewer.pal(n=11, name="RdYlBu")))(100),
         #color=colorRampPalette(c("blue","white","red"))(20) ,
         color=colors,
         breaks = breaks,
         scale = "none",
         annotation_col = df_column_pheatmap2,
         show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = as.hclust(hclust_col_m_Metagene_metaheatmap_all_signatures_workdir_ATNF_curated_onlyatnf),
         cluster_rows = as.hclust(hclust_row_m_Metagene_metaheatmap_all_signatures_workdir_ATNF_curated_onlyatnf),
         #annotation_row = df_row_pheatmap,
         annotation_colors = ann_colors#,
         #clustering_callback = callback
)

