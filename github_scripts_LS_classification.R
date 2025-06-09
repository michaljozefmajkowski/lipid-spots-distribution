# to classify cells basing of lipid spots distribution
setwd("your working directory")

spot_centroids <- read_excel("your_exel_file.xlsx") [, c(2,3,6,7,21,22,23,24,26,27,28,29)] # centroids form Harmony software
spot_centroids <- spot_centroids %>% rename("cell_number" = `Spots Selected - Object No in Nuclei_SS`)


spot_centroids <-  rbind(spot_centroids, spot_centroids2)

object_charact <- read_excel("your_exel_file_objects.xlsx")[,c(2,3,6,7,21,22,23,24,25,26,27)] # here goes cells descrition from Harmony software
object_charact <-  object_charact %>% rename("cell_number" = `Object ID`) # obj. is a cell (Nuclei_SS)


compute_within_cluster_distances <- function(df) {
  df %>%
    # filter(cluster != "0") %>% 
    group_by(cluster) %>%
    summarise(
      mean_distance_within_clust = mean(
        as.vector(dist(cbind(`Spots Selected - Spot Centroid X in Image [µm]`, `Spots Selected - Spot Centroid Y in Image [µm]`)))),  # Mean pairwise distance
      min_distance_within_clust = min(
        as.vector(dist(cbind(`Spots Selected - Spot Centroid X in Image [µm]`,  `Spots Selected - Spot Centroid Y in Image [µm]`)))),  # Min distance
      max_distance_within_clust = max(
        as.vector(dist(cbind(`Spots Selected - Spot Centroid X in Image [µm]`, `Spots Selected - Spot Centroid Y in Image [µm]`))))   # Max distance
    )
}

row <- sort(unique(spot_centroids$Row))
column <- sort(unique(spot_centroids$Column))
field <- sort(unique(spot_centroids$Field))

clusters_df_list3 <- list()
sum_low_ls_list_cell_num <- list()

#removal of rows with NA
spot_centroids <- na.omit(spot_centroids)

complete.cases(spot_centroids) %>% length()

#------------------------------------------------------------ First - df with number of clusters per cell    --------------------------------

sum_low_ls_list_cell_num <- list()

for(r in 1 : 2){ # r means row
  
  clusters_df_list2 <- list()
  
  for(c in 1 : 12){ # c means column
    
    clusters_df_list1 <- list()
    
    for(f in 1 : 24){ # f means field
      
      clusters_df_list <- list()
      
      single_field <-  spot_centroids %>% 
        filter(Row == row[r],
               Column == column[c],
               Field == field[f]) %>%
        select(Row,
               Column,
               Field,
               cell_number,
               `Spots Selected - Spot Centroid X in Image [µm]`, 
               `Spots Selected - Spot Centroid Y in Image [µm]`,
               `Spots Selected - Spot Nearest Neighbor Distance [µm]`,
               `Spots Selected - Spot Contact Area with Neighbors [%]`)
      
      cell_num_in_field <-  single_field %>% 
        select(cell_number) %>% 
        unique() %>% 
        pull()
      
      if(length(cell_num_in_field) > 0) {
        
        for(cell in 1 : length(cell_num_in_field)) {
          
          min_spots <- single_field %>% 
            filter(cell_number == cell_num_in_field[cell]) %>% 
            nrow()    
          
          if(min_spots <= 4) {
            
            low_ls_list_cell_num <- list()
            
            low_ls_cell_number <- vector()
            
            low_ls_cell_number <- single_field %>% 
              group_by(cell_number) %>% 
              select(cell_number) %>% 
              summarise(n=n()) %>% 
              filter(n == min_spots) %>% 
              select(cell_number) %>% 
              pull()
            
            for(low in 1 : length(low_ls_cell_number)) {
              
              low_ls_df <-  spot_centroids %>% 
                filter(Row == row[r],
                       Column == column[c],
                       Field == field[f],
                       cell_number == low_ls_cell_number[low]) %>% 
                select(Row,
                       Column,
                       Field,
                       cell_number)
              
              low_ls_list_cell_num[[low]] <- low_ls_df
            }
            
            sum_low_ls_list_cell_num <- append(sum_low_ls_list_cell_num, low_ls_list_cell_num)
            
          } else {
            
            single_cell <-  single_field %>% 
              filter(cell_number == cell_num_in_field[cell]) %>%  
              select(`Spots Selected - Spot Centroid X in Image [µm]`, 
                     `Spots Selected - Spot Centroid Y in Image [µm]`,
                     `Spots Selected - Spot Nearest Neighbor Distance [µm]`,
                     `Spots Selected - Spot Contact Area with Neighbors [%]`)
            
            # calculation of distances
            sorted_dist <- sort(kNNdist(single_cell[ , c(1:2)], k = 4))
            
            # Compute first and second derivative approximations
            first_derivative <- diff(sorted_dist)
            second_derivative <- diff(first_derivative)
            
            # Find the index of the maximum curvature
            opt_eps_index <- which.max(second_derivative)  
            opt_eps <- sorted_dist[opt_eps_index]  # Best `eps`
            
            clusters <- dbscan(single_cell[,c(1:2)], eps = opt_eps, minPts = 4)
            
            # number of clusters per cell / noise points removed
            clust_number <- data.frame(table(clusters$cluster)) %>% 
              filter(Var1 != 0) %>% 
              nrow()
            
            # df with number of clustres per cell
            clusters_df <- data.frame(cell_number = cell_num_in_field[cell],
                                      Row = row[r],
                                      Column = column[c],
                                      Field = field[f],
                                      num_of_clusters = clust_number)
            
            # list with df with number of clustres per cell
            clusters_df_list[[cell]] <- clusters_df
          }
        }
        
      }
      clusters_df1 <- do.call(rbind, clusters_df_list)
      clusters_df_list1[[f]] <- clusters_df1
    }
    clusters_df2 <-  do.call(rbind, clusters_df_list1)
    clusters_df_list2[[c]] <- clusters_df2
  }
  clusters_df3 <-  do.call(rbind, clusters_df_list2)
  clusters_df_list3[[r]] <- clusters_df3
}
no_of_clusters_cell <- do.call(rbind, clusters_df_list3)
write.csv(no_of_clusters_cell, "no_of_clusters_cell_row_F.csv")
distinct(no_of_clusters_cell)

# df contains info about cells with less of 4 lipid spots per cell
df <- do.call(
  rbind,
  unique(sum_low_ls_list_cell_num))

df_first <-  df %>% group_by(Row, Column, Field, cell_number) %>% 
  select(Row) %>% 
  summarise(n=n())
write.csv(df_first, "df_first.csv")
# -----------------------------------------------------------------------------------------  Second; df with clusters characteristics
#                                                                                  mean_NN_distance; mean_NN_contact; percentage of ls in clusters
clust_charact_list3 <- list()
sum_low_ls_list_cell_num_second <- list()

for(r in 1 : 2){ # r means row
  
  clust_charact_list2 <- list()
  
  for(c in 1 : 12){ # c means column
    
    clust_charact_list1 <- list()
    
    for(f in 1 : 24){ # f means field
      
      clust_charact_list <- list()
      
      single_field <-  spot_centroids %>% 
        filter(Row == row[r],
               Column == column[c],
               Field == field[f]) %>%
        select(Row,
               Column,
               Field,
               cell_number,
               `Spots Selected - Spot Centroid X in Image [µm]`, 
               `Spots Selected - Spot Centroid Y in Image [µm]`,
               `Spots Selected - Spot Nearest Neighbor Distance [µm]`,
               `Spots Selected - Spot Contact Area with Neighbors [%]`)
      
      cell_num_in_field <-  single_field %>% 
        select(cell_number) %>% 
        unique() %>% 
        pull()
      
      if(length(cell_num_in_field) > 0) {
        
        for(cell in 1 : length(cell_num_in_field)) {
          
          min_spots <- single_field %>% 
            filter(cell_number == cell_num_in_field[cell]) %>% 
            nrow()    
          
          if(min_spots <= 4) {
            
            low_ls_list_cell_num_second <- list()
            
            low_ls_cell_number_second <- vector()
            
            low_ls_cell_number_second <- single_field %>% 
              group_by(cell_number) %>% 
              select(cell_number) %>% 
              summarise(n=n()) %>% 
              filter(n == min_spots) %>% 
              select(cell_number) %>% 
              pull()
            
            for(low in 1 : length(low_ls_cell_number_second)) {
              
              low_ls_df <-  spot_centroids %>% 
                filter(Row == row[r],
                       Column == column[c],
                       Field == field[f],
                       cell_number == low_ls_cell_number_second[low]) %>% 
                select(Row,
                       Column,
                       Field,
                       cell_number)
              
              low_ls_list_cell_num_second[[low]] <- low_ls_df
            }
            
            sum_low_ls_list_cell_num_second <- append(sum_low_ls_list_cell_num_second, low_ls_list_cell_num_second)
            
          } else {
            
            single_cell <-  single_field %>% 
              filter(cell_number == cell_num_in_field[cell]) %>%  
              select(`Spots Selected - Spot Centroid X in Image [µm]`, 
                     `Spots Selected - Spot Centroid Y in Image [µm]`,
                     `Spots Selected - Spot Nearest Neighbor Distance [µm]`,
                     `Spots Selected - Spot Contact Area with Neighbors [%]`)
            
            # calculation of distances
            sorted_dist <- sort(kNNdist(single_cell[ , c(1:2)], k = 4))
            
            # Compute first and second derivative approximations
            first_derivative <- diff(sorted_dist)
            second_derivative <- diff(first_derivative)
            
            # Find the index of the maximum curvature
            opt_eps_index <- which.max(second_derivative)  
            opt_eps <- sorted_dist[opt_eps_index]  # Best `eps`
            
            clusters <- dbscan(single_cell[,c(1:2)], eps = opt_eps, minPts = 4)
            
            # adding to a df with single cell column with numbers
            # zero are noise poitns
            # other numbers defines clusters
            single_cell$cluster <- as.factor(clusters$cluster)
            
            # sum of points in a cell
            # needed to calculate percentage of LS in clusters
            sum <- sum(data.frame(table(clusters$cluster))[ , 2])
            
            # df with number defining cluster, number of spots in cluster, percentage of spots in cluster
            # + 2 x mean 
            clust_charact <- single_cell %>% 
              group_by(cluster) %>% 
              summarise(n = n(),
                        mean_NN_dist = mean(`Spots Selected - Spot Nearest Neighbor Distance [µm]`),
                        mean_NN_contact = mean(`Spots Selected - Spot Contact Area with Neighbors [%]`),
                        percent = n / sum * 100) %>% 
              add_column(cell_number = cell_num_in_field[cell],
                         Row = row[r],
                         Column = column[c],
                         Field = field[f])
            
            clust_charact_list[[cell]] <- clust_charact
            
          }
        }
        
      }
      clust_charact_df <- do.call(rbind, clust_charact_list)
      clust_charact_list1[[f]] <- clust_charact_df
    }
    clust_charact_df2 <- do.call(rbind, clust_charact_list1)
    clust_charact_list2[[c]] <-  clust_charact_df2
  }
  clust_charact_df3 <- do.call(rbind, clust_charact_list2)
  clust_charact_list3[[r]] <-  clust_charact_df3
}

clust_charact_df_row_E <- do.call(rbind, clust_charact_list3)
write.csv(clust_charact_df_row_E, "clust_charact_df_row_F.csv")


df <- do.call(
  rbind,
  unique(sum_low_ls_list_cell_num_second))

df_second <- df %>% group_by(Row, Column, Field, cell_number) %>% 
  select(Row) %>% 
  summarise(n=n())

write.csv(df_second, "df_second.csv")

# ---------------------------------------------------------- Third df with distances between clusters centroids
# ---------------------------------------------------------- df with distances within clusters

dist_centroids_list3 <- list()
dist_clust_list3 <- list()
sum_low_ls_list_cell_num_third <- list()

for(r in 1 : 2){ # r means row
  
  dist_centroids_list2 <- list()
  dist_clust_list2 <- list()
  
  
  for(c in 1 : 12){ # c means column
    
    dist_centroids_list1 <- list()
    dist_clust_list1 <- list()
    
    for(f in 1 : 24){ # f means field
      
      dist_centroids_list <- list()
      dist_clust_list <- list()
      
      single_field <-  spot_centroids %>% 
        filter(Row == row[r],
               Column == column[c],
               Field == field[f]) %>%
        select(Row,
               Column,
               Field,
               cell_number,
               `Spots Selected - Spot Centroid X in Image [µm]`, 
               `Spots Selected - Spot Centroid Y in Image [µm]`,
               `Spots Selected - Spot Nearest Neighbor Distance [µm]`,
               `Spots Selected - Spot Contact Area with Neighbors [%]`)
      
      cell_num_in_field <-  single_field %>% 
        select(cell_number) %>% 
        unique() %>% 
        pull()
      
      if(length(cell_num_in_field) > 0) {
        
        for(cell in 1 : length(cell_num_in_field)) {
          
          min_spots <- single_field %>% 
            filter(cell_number == cell_num_in_field[cell]) %>% 
            nrow()    
          
          if(min_spots <= 4) {
            
            low_ls_list_cell_num_third <- list()
            
            low_ls_cell_number_third <- vector()
            
            low_ls_cell_number_third <- single_field %>% 
              group_by(cell_number) %>% 
              select(cell_number) %>% 
              summarise(n=n()) %>% 
              filter(n == min_spots) %>% 
              select(cell_number) %>% 
              pull()
            
            for(low in 1 : length(low_ls_cell_number_third)) {
              
              low_ls_df <-  spot_centroids %>% 
                filter(Row == row[r],
                       Column == column[c],
                       Field == field[f],
                       cell_number == low_ls_cell_number_third[low]) %>% 
                select(Row,
                       Column,
                       Field,
                       cell_number)
              
              low_ls_list_cell_num_third [[low]] <- low_ls_df
            }
            
            sum_low_ls_list_cell_num_third <- append(sum_low_ls_list_cell_num_third, low_ls_list_cell_num_third)
            
          } else {
            
            single_cell <-  single_field %>% 
              filter(cell_number == cell_num_in_field[cell]) %>%  
              select(`Spots Selected - Spot Centroid X in Image [µm]`, 
                     `Spots Selected - Spot Centroid Y in Image [µm]`,
                     `Spots Selected - Spot Nearest Neighbor Distance [µm]`,
                     `Spots Selected - Spot Contact Area with Neighbors [%]`)
            
            # calculation of distances
            sorted_dist <- sort(kNNdist(single_cell[ , c(1:2)], k = 4))
            
            # Compute first and second derivative approximations
            first_derivative <- diff(sorted_dist)
            second_derivative <- diff(first_derivative)
            
            # Find the index of the maximum curvature
            opt_eps_index <- which.max(second_derivative)  
            opt_eps <- sorted_dist[opt_eps_index]  # Best `eps`
            
            clusters <- dbscan(single_cell[,c(1:2)], eps = opt_eps, minPts = 4)
            
            # adding to a df with single cell column with numbers
            # zero are noise poitns
            # other numbers defines clusters
            single_cell$cluster <- as.factor(clusters$cluster) 
            
            # caluculation of clusters centroids
            centroids <- single_cell %>%
              filter(cluster != "0") %>% 
              group_by(cluster) %>%
              summarise(X = mean(`Spots Selected - Spot Centroid X in Image [µm]`), 
                        Y = mean(`Spots Selected - Spot Centroid Y in Image [µm]`))
            
            # df with distances between cluster centroides
            dist_centroids_list[[cell]] <- data.frame(max_dist_between_clust_centroids =  max(as.matrix(dist(centroids[, 2:3]))),
                                                      cell_number = cell_num_in_field[cell],
                                                      Row = row[r],
                                                      Column = column[c],
                                                      Field = field[f])
            
            # df with distances within clusters
            dist_clust_list[[cell]] <- compute_within_cluster_distances(single_cell) %>% 
              add_column(cell_number = cell_num_in_field[cell],
                         Row = row[r],
                         Column = column[c],
                         Field = field[f])
          }
        }
      }
      
      dist_clust_centroids_df <- do.call(rbind, dist_centroids_list)
      dist_centroids_list1[[f]] <-  dist_clust_centroids_df
      
      dist_clust_df <- do.call(bind_rows, dist_clust_list)
      dist_clust_list1[[f]] <-  dist_clust_df
      
    }
    
    dist_clust_centroids_df1 <- do.call(rbind, dist_centroids_list1)
    dist_centroids_list2[[c]] <-  dist_clust_centroids_df1
    
    dist_clust_df1 <- do.call(bind_rows,  dist_clust_list1)
    dist_clust_list2[[c]] <-  dist_clust_df1
    
  }
  
  dist_clust_centroids_df2 <- do.call(rbind, dist_centroids_list2)
  dist_centroids_list3[[r]] <-  dist_clust_centroids_df2
  
  dist_clust_df2 <- do.call(bind_rows,  dist_clust_list2)
  dist_clust_list3[[r]] <-  dist_clust_df2
  
}

dist_clust_centroids <- do.call(rbind, dist_centroids_list3)
write.csv(dist_clust_centroids, "dist_clust_centroids_row_F.csv")

dist_clust <- do.call(bind_rows, dist_clust_list3)
write.csv(dist_clust, "dist_clust_row_F.csv")



df <- do.call(
  rbind,
  unique(sum_low_ls_list_cell_num_third))

df_third <- df %>% group_by(Row, Column, Field, cell_number) %>% 
  select(Row) %>% 
  summarise(n=n())

write.csv(df_third, "df_third_row_F.csv")



# ---------------------------------------------------------- classification part

setwd("your working directory")
no_of_clusters_cell <- read.csv("no_of_clusters.csv") # number of clusters per cell
clust_charact_df <- read.csv("clust_charact_df.csv")
dist_clust_centroids <- read.csv("dist_clust_centroids.csv")
dist_clust <- read.csv("dist_clust.csv")
object_charact <- read_excel("object_character.xlsx")
object_charact <-  object_charact[,c(2,3,6,7,21,22,23)]
object_charact <-  object_charact %>% rename("cell_number" = `Object ID`)  # obj. is a cell 

length(complete.cases(object_charact)) 

spots <- merge(
  merge(no_of_clusters_cell, dist_clust_centroids),
  merge(clust_charact_df, dist_clust), 
  by = c("Row", "Column", "Field", "cell_number"), all.x = T, all.y = T) 
colnames(object_charact)

ggplot(object_charact %>% filter(Column %in% c(7:8)), aes(`nuclei_input - Cell Length [µm]`)) + 
  geom_histogram() + 
  facet_wrap(~Column) +
  theme_minimal()

object_charact %>% 
  filter(Column %in% c(7:8)) %>% 
  group_by(Column) %>%
  summarise(n=n(),
            area = mean(`nuclei_input - Cell Area [µm²]`),
            length = mean(`nuclei_input - Cell Length [µm]`))

objects_spots <- merge(object_charact, spots, by = c("Row", "Column", "Field", "cell_number"), all.x = T, all.y = T) 



objects_spots <- objects_spots %>% 
  select(-starts_with("X")) %>% 
  rename(number_of_spots_in_cluster = n)

objects_spots <-  objects_spots %>% 
  mutate(cell_cluster_length_ratio = `nuclei_input - Cell Length [µm]`/max_dist_between_clust_centroids,
         cell_max_cluster_size_ration = `nuclei_input - Cell Length [µm]`/max_distance_within_clust)

clust_groups <-  objects_spots[,c(1:4,8,10,11,14,18,19)]



# only single cluster per cell
single_cluster <- clust_groups %>% 
  filter(num_of_clusters == 1, cluster != 0)

clust_groups_2 <- clust_groups %>% # 2 at the end means at least two clusters per cell
  filter(num_of_clusters > 1, cluster != 0)

# clust_groups_2 <- clust_groups_2[,c(1:4, 8:13, 5:7)]

#-------------------------------------------------------------------------------------------------


ggplot(clust_groups_2, aes(cell_cluster_length_ratio)) + geom_histogram() 
ggplot(clust_groups_2, aes(cell_max_cluster_size_ration)) + geom_histogram() 
ggplot(clust_groups_2, aes(cell_max_cluster_size_ration, cell_cluster_length_ratio)) + geom_point(alpha = 0.3) + theme_minimal()

quantile(clust_groups_2$cell_cluster_length_ratio, probs =  seq(0,1, 0.1), na.rm = T)

# number of clusters per cell
row <- sort(unique(clust_groups_2$Row))
column <- sort(unique(clust_groups_2$Column))
field <- sort(unique(clust_groups_2$Field))

list_row <- list()

for(r in 1 : 4){ # r means row
  
  list_column <- list()
  
  for(c in 1 : 2){ # c means column
    
    list_field <- list()
    
    for(f in 1 : 32){ # f means field
      
      list_cell <- list()
      
      cell_in_field <- clust_groups_2 %>% 
        filter(Row == row[r],
               Column == column[c],
               Field == field[f]) %>% 
        select(cell_number) %>% 
        unique() %>% 
        pull()
      
      if(length(cell_in_field) > 0) {
        
        for(cell in 1 : length(cell_in_field)) {
          
          single_cell_df <-   clust_groups_2 %>% 
            filter(Row == row[r],
                   Column == column[c],
                   Field == field[f],
                   cell_number == cell_in_field[cell])
          
          df_cell <- single_cell_df %>% 
            filter(percent == max(single_cell_df$percent))
          
          if(nrow(df_cell) >= 2){
            df_cell <-  df_cell %>% filter(cell_max_cluster_size_ration == min(df_cell$cell_max_cluster_size_ration))
            list_cell[[cell]] <- df_cell
          } else {
            
            list_cell[[cell]] <- df_cell
          }
        }
      }
      
      df_cell <- do.call(rbind, list_cell)
      list_field [[f]] <- df_cell
    }
    df_field <- do.call(rbind, list_field)
    list_column[[c]] <- df_field
  }
  df_column <- do.call(rbind, list_column)
  list_row[[r]] <- df_column
}

#number of clusters per cell; at least two clusters per cell
num_clust <- do.call(rbind, list_row)


# loop to classify cells based on spots distribution: clustered, dispersed, intermediate
row <- sort(unique(clust_groups_2$Row))
column <- sort(unique(clust_groups_2$Column))
field <- sort(unique(clust_groups_2$Field))

final_list <- list()

for(r in 1 : 4) { # r means row
  
  pseudo_raw_list <- list()
  
  for(c in 1 : 2) { # c means column
    
    pseudo_col_list <- list()
    
    for(f in 1 : 31) { # f means field
      
      field_list <- list()    
      
      single_field <-  clust_groups_2 %>% 
        filter(Row == row[r],
               Column == column[c],
               Field == field[f])
      
      cells <- sort(unique(single_field$cell_number))
      
      if(length(cells) > 0) {
        
        for(cell in 1 : length(cells)) { # cell means cell_number; length of "cells"
          
          single_cell_df <- single_field %>% # single_cell_df contains info from single cell; 
            filter(cell_number == cells[cell])
          
          # dispersed - 2 clusters
          if(
            # 
            max(single_cell_df$num_of_clusters) == 2   &    
            #
            max(single_cell_df$percent) < 90 &
            # 
            sum(single_cell_df$percent) > 70 &                  
            #                  
            max(single_cell_df$cell_cluster_length_ratio) <= 1.5  
          )                  
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "dispersed")}
          # intermediate - 2 clusters
          else if(
            # 
            max(single_cell_df$num_of_clusters) == 2   &    
            #
            max(single_cell_df$percent) < 90 &
            # 
            sum(single_cell_df$percent) > 70 &                  
            #                  
            max(single_cell_df$cell_cluster_length_ratio) >= 1.9 &
            #
            single_cell_df %>% 
            select(cell_max_cluster_size_ration) %>% 
            pull() %>% 
            min() < 2
          )
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "intermediate")}
          # clustered - 2 clusters
          else if(
            # 
            max(single_cell_df$num_of_clusters) == 2   &    
            #
            max(single_cell_df$percent) < 90 &
            # 
            sum(single_cell_df$percent) > 70 &                  
            #                  
            max(single_cell_df$cell_cluster_length_ratio) >= 1.9 &
            #
            single_cell_df %>% 
            select(cell_max_cluster_size_ration) %>% 
            pull() %>% 
            min() > 2
          )
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "clustered")}
          
          
          # intermediate - 2 clusters
          else if(
            #
            max(single_cell_df$num_of_clusters) == 2 &    
            #
            max(single_cell_df$percent) < 90 &
            #
            sum(single_cell_df$percent)  > 70 &                    
            #
            max(single_cell_df$cell_cluster_length_ratio) > 1.5 &  
            max(single_cell_df$cell_cluster_length_ratio) < 1.9 &
            #
            single_cell_df %>% 
            select(cell_max_cluster_size_ration) %>% 
            pull() %>% 
            min() > 1.6
          ) 
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "intermediate")}
          #  min ratio between single cluster size and cell size  
          # dispersed because single cluster size determines clasiffication
          else if(
            #
            max(single_cell_df$num_of_clusters) == 2 &    
            #
            max(single_cell_df$percent) < 90 &
            #
            sum(single_cell_df$percent)  > 70 &                    
            #
            max(single_cell_df$cell_cluster_length_ratio) > 1.5 &  
            max(single_cell_df$cell_cluster_length_ratio) < 1.9 &
            #
            single_cell_df %>% 
            select(cell_max_cluster_size_ration) %>% 
            pull() %>% 
            min() <= 1.6
          ) 
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "dispersed")}
          
          # dispersed 3 clusters       
          else if(
            #
            max(single_cell_df$num_of_clusters) == 3 &       
            #
            max(single_cell_df$percent) < 90 &
            #
            sum(single_cell_df$percent)  > 70 &                    
            #
            max(single_cell_df$cell_cluster_length_ratio) <= 1.6  
          ) 
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "dispersed")
          }
          # clustered - 3 clusters
          else if(
            # 
            max(single_cell_df$num_of_clusters) == 3   &    
            #
            max(single_cell_df$percent) < 90 &
            # 
            sum(single_cell_df$percent) > 70 &                  
            #                  
            max(single_cell_df$cell_cluster_length_ratio) >= 2 
          )
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "clustered")}
          # intermediate - 3 clusters
          else if(
            #
            max(single_cell_df$num_of_clusters) == 3 &     
            #
            max(single_cell_df$percent) < 90 &
            #
            sum(single_cell_df$percent)  > 70 &                    
            #
            max(single_cell_df$cell_cluster_length_ratio) > 1.6 &  
            max(single_cell_df$cell_cluster_length_ratio) < 2 &
            #
            single_cell_df %>% 
            select(cell_max_cluster_size_ration) %>% 
            pull() %>% 
            min() > 1.6
          ) 
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "intermediate")}
          #  min ratio between single cluster size and cell size  
          # dispersed because single cluster size determines clasiffication   
          else if(
            #
            max(single_cell_df$num_of_clusters) == 3 &     
            #
            max(single_cell_df$percent) < 90 &
            #
            sum(single_cell_df$percent)  > 70 &                    
            #
            max(single_cell_df$cell_cluster_length_ratio) > 1.6 &  
            max(single_cell_df$cell_cluster_length_ratio) < 2 &
            #
            single_cell_df %>% 
            select(cell_max_cluster_size_ration) %>% 
            pull() %>% 
            min() <= 1.6
          ) 
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "dispersed")}
          # Dispersed; 4 clusters         
          else if(
            #
            max(single_cell_df$num_of_clusters) >= 4 &   
            #
            max(single_cell_df$percent) < 90 &
            #
            sum(single_cell_df$percent)  > 70 &                    
            #
            max(single_cell_df$cell_cluster_length_ratio) > 1.1
          ) 
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "dispersed")}
          #   single cluster >=90% spots       
          else if(
            #
            max(single_cell_df$num_of_clusters) >= 2 &
            #
            max(single_cell_df$percent) >= 90 &
            #
            sum(single_cell_df$percent)  > 70 &                    
            #
            single_cell_df %>% 
            filter(percent == max(single_cell_df$percent)) %>%
            select(cell_max_cluster_size_ration) %>% 
            pull() %>% 
            max() <= 2
          ) 
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "dispersed")}  
          #   single cluster >=90% spots       
          else if(
            #
            max(single_cell_df$num_of_clusters) >= 2 &
            #
            max(single_cell_df$percent) >= 90 &
            #
            sum(single_cell_df$percent)  > 70 &                    
            #
            single_cell_df %>% 
            filter(percent == max(single_cell_df$percent)) %>%
            select(cell_max_cluster_size_ration) %>% 
            pull() %>% 
            max() > 2
          ) 
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "intermediate")}  
          #   between 70 and 40 % of LS in all clusters
          else if(
            #
            max(single_cell_df$num_of_clusters) >= 2 &
            #
            max(single_cell_df$percent) < 90 &
            #
            sum(single_cell_df$percent)  < 70 &
            sum(single_cell_df$percent)  > 40
          ) 
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "intermediate")}   
          
          #   less that 40 % of LS in all clusters  
          else if(
            #
            max(single_cell_df$num_of_clusters) >= 2 &
            #
            max(single_cell_df$percent) < 90 &
            #
            sum(single_cell_df$percent)  <= 40
            
          ) 
          {df1 <- data.frame(Row = row[r],
                             Column = column[c],
                             field = field[f],
                             cell_number = cells[cell],
                             num_of_clusters = max(single_cell_df$num_of_clusters),
                             
                             cell_class = "dispersed")}   
          
          field_list[[cell]] <- df1
        }
        df2 <- do.call(rbind, field_list)
        pseudo_col_list[[f]] <- df2 
      }
    }
    df3 <-  do.call(rbind, pseudo_col_list)
    pseudo_raw_list[[c]] <- df3
  }
  df4 <- do.call(rbind, pseudo_raw_list)
  final_list[[r]] <- df4
}

# classes for cells with at least two groups
classes <- do.call(rbind, final_list)

# analysis of df with single cluster per cell ---------------------------------------------------------------------------------------

head(single_cluster)

row <- sort(unique(single_cluster$Row))
column <- sort(unique(single_cluster$Column))
field <- sort(unique(single_cluster$Field))

final_list_single_clu <- list()

for(r in 1 : 4) { # r means row
  
  row_list_single_clu <- list()
  
  for(c in 1 : 2) { # c means column
    
    col_list_single_clu <- list()
    
    for(f in 1 : 32) { # f means field
      
      field_list_single_clu <- list()    
      
      single_field_single_clu <-  single_cluster %>% 
        filter(Row == row[r],
               Column == column[c],
               Field == field[f])
      
      cells_s_c <- sort(unique(single_field_single_clu$cell_number)) # s_c mean single cluster
      
      if(length(cells_s_c) > 0) {
        
        for(cell in 1 : length(cells_s_c)) { # cell means cell_number; length of "cells"
          
          single_cell_df_s_c <- single_field_single_clu %>% # single_cell_df contains info from single cell; 
            filter(cell_number == cells_s_c[cell])
          
          
          if(single_cell_df_s_c$cell_max_cluster_size_ration <= 1.5)     
            
          {df1_s_c <- data.frame(Row = row[r],
                                 Column = column[c],
                                 field = field[f],
                                 cell_number = cells_s_c[cell],
                                 num_of_clusters = max(single_cell_df_s_c$num_of_clusters),
                                 
                                 cell_class = "dispersed")}
          
          else if(single_cell_df_s_c$cell_max_cluster_size_ration >= 2)
            
          {df1_s_c <- data.frame(Row = row[r],
                                 Column = column[c],
                                 field = field[f],
                                 cell_number = cells_s_c[cell],
                                 num_of_clusters = max(single_cell_df_s_c$num_of_clusters),
                                 
                                 cell_class = "clustered")}
          
          else if(single_cell_df_s_c$cell_max_cluster_size_ration < 2 &
                  single_cell_df_s_c$cell_max_cluster_size_ration > 1.5)
            
          {df1_s_c <- data.frame(Row = row[r],
                                 Column = column[c],
                                 field = field[f],
                                 cell_number = cells_s_c[cell],
                                 num_of_clusters = max(single_cell_df_s_c$num_of_clusters),
                                 
                                 cell_class = "intermediate")}
          
          field_list_single_clu [[cell]] <- df1_s_c
        }
        df2_s_c <- do.call(rbind, field_list_single_clu)
        col_list_single_clu [[f]] <- df2_s_c 
      }
    }
    df3_s_c <-  do.call(rbind, col_list_single_clu)
    row_list_single_clu[[c]] <- df3_s_c
  }
  df4_s_c <- do.call(rbind, row_list_single_clu)
  final_list_single_clu[[r]] <- df4_s_c
}
#s_c - means single cell.
classes_s_c <- do.call(rbind, final_list_single_clu)

classes_all <-  rbind(classes, classes_s_c)

unique(classes_all[,2])

write.csv(classes_all, "classes_all_A375.csv")

colnames(classes_all)
classes_all %>% filter(Column == 8) %>% nrow()

ggplot(classes_all, aes(cell_class, y =..prop.., group = Column, fill = Column)) + 
  geom_bar(position = "dodge") + 
  xlab("cell class") +
  ylab("proportion of cells") +
  facet_wrap(~interaction(Column, Row), ncol = 2) +
  theme_minimal() +
  theme(strip.text = element_text(size = 5),
        axis.title.x = element_text(size = 10),  # Bigger x-axis label
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 0, hjust = 1),  # x-axis tick labels size
        axis.text.y = element_text(size = 6),
        legend.position = "none")

ggsave("classes_all.png")

#---------------------------------------------------------------
# extra some data miniing
obj_spots <- objects_spots %>% filter(cell_cluster_length_ratio != Inf)

ggplot(obj_spots, aes(`Nuclei_SS - Intensity Cell Alexa 488 Mean`, cell_cluster_length_ratio)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~ FHOD1_version) 

ggplot(single_cluster, aes(`Nuclei_SS - Intensity Cell Alexa 488 Mean`, cell_max_cluster_size_ration)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~ FHOD1_version) +
  xlab("Mean GFP intensity") +
  ylab("cell length / max dist. within cluster") +
  guides(color = guide_legend(override.aes = list(size = 20))) +
  theme_minimal() +
  theme(strip.text = element_text(size = 10),
        axis.title.x = element_text(size = 14),  # Bigger x-axis label
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # x-axis tick labels size
        axis.text.y = element_text(size = 12),
        legend.position = "none")
ggsave("single_clust.png")

ggplot(clust_groups_2, aes(`Nuclei_SS - Intensity Cell Alexa 488 Mean`, cell_cluster_length_ratio)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~ FHOD1_version) +
  xlab("Mean GFP intensity") +
  ylab("cell length / max dist. between cluster") +
  theme_minimal() +
  theme(strip.text = element_text(size = 10),
        axis.title.x = element_text(size = 14),  # Bigger x-axis label
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # x-axis tick labels size
        axis.text.y = element_text(size = 12))
ggsave("between_clusters.png")

ggplot(single_cluster, aes(x = transfection_status, y = mean_NN_contact, colour = transfection_status)) +
  geom_jitter(width = 0.4, alpha = 0.2, size = 3) +
  scale_colour_manual(values = c("negative" = "black", "positive" = "green4")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
  facet_wrap(~ FHOD1_version) +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 3, colour = "red") +
  stat_summary(fun = median, geom = "point", shape = 16, size = 3, colour = "blue") +
  scale_colour_manual(values = c("negative" = "black", "positive" = "green4")) +
  scale_fill_manual(values = c("negative" = "black", "positive" = "green4")) +
  xlab("transfection status") +
  ylab("mean spots contact ") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_minimal() +
  theme(strip.text = element_text(size = 6),
        axis.title.x = element_text(size = 14),  # Bigger x-axis label
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5, hjust = 1),  # x-axis tick labels size
        axis.text.y = element_text(size = 12))
ggsave("mean_spots_contact_tow_more.png")

write.csv(single_cluster, "single_cluster_EF.csv")


























