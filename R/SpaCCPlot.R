#' Visualize the clustering results
#'
NULL

##' Run SpaCCPlot
##'
##' Windows users need to install the ggplot2 package
##'
##' @param "feature" is the name of the metadata column that needs to be displayed in color on the spatial graph. The string "/" is mandatory. Common values include "SpaCC" (the function will internally change it to "C1", "C2", etc.), and it can also be "seurat_clusters", "celltype", etc.
##' @param colors: Discrete color palette, character vector/default NULL. If NULL is passed, the built-in 16-color universal discrete color palette of the function is used. Custom vectors can be overridden.
##' @param pt_size: spot size, value/default NULL, if NULL is passed, set to 2.7; It can be adjusted according to the spot density or resolution by yourself.
##'
##' 1. Internal temporary variable
##' df_raw: Data frame of spatial coordinates + metadata
##' SpaCC = paste0("C", SpaCC) : Only takes effect when feature == "SpaCC". Change the number tags to the form of "C1" and "C2" for easier reading of the legend.
##' colors (second appearance) : Hard-coded 16-color discrete color palette within the function
##' pt_size (second appearance) : Final point size
##'
##' 2. Drawing details
##' aes(x = imagecol, y = imagerow) : Plot the spatial position using the pixel coordinates calculated by Seurat.
##' colour =.data[[feature]] : Color discretely by feature column.
##' scale_color_manual(values = colors, name = feature) : Manually specify the color and set the legend title.
##' scale_y_reverse() : Align the Y-axis with the direction of the tissue section (positive when the Y-axis in the image coordinate system is downward).
##' theme_bw() + subsequent theme(...) Remove the grid lines, coordinate axis scales and text, and keep the black border to make the picture clean.
##'
##' 3. Return value
##' The p: ggplot object can be directly printed (p) to be displayed in the RStudio plotting window, and it also supports ggsave() to export high-resolution images.
##'
SpaCCPlot <- function(seo,
                      feature = NULL,
                      colors = NULL,
                      pt_size = NULL) {
  if (!is.null(feature)) {
    df_raw <- SpatialDimPlot(seo, group.by = feature)[[1]][["data"]]
    
    if (feature == "SpaCC") {
      df_raw <- df_raw %>% mutate(SpaCC = paste0("C", SpaCC))
    } else {
      df_raw <- df_raw
    }
  } else {
    stop("Please provide the feature")
  }
  
  colors <- c(
    "#E41A1C",
    "#377EB8" ,
    "#4DAF4A",
    "#984EA3",
    "#FF7F00" ,
    "#FFFF33" ,
    "#A65628",
    "#F781BF" ,
    "#999999",
    "#90EE90",
    "#00CD00",
    "#008B8B",
    "#6495ED",
    "#FFC1C1",
    "#CD5C5C",
    "#00F5FF"
  )
  
  pt_size <- if (is.null(pt_size))
    2.7
  else
    pt_size
  
  p <- ggplot(df_raw, aes(x = imagecol, y = imagerow, colour = .data[[feature]])) +
    geom_point(size = pt_size) +
    scale_color_manual(values = colors, name = feature) +
    scale_y_reverse() +
    labs(x = "Spatial_1", y = "Spatial_2") +
    theme_bw() +
    theme(
      panel.border     = element_rect(colour = "black", fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(size = 14),
      text             = element_text(size = 14),
      legend.title     = element_text(size = 12),
      legend.text      = element_text(size = 10),
      axis.ticks       = element_blank(),
      axis.text        = element_blank()
    )
  
  return(p)
}