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

