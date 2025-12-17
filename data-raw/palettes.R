palette_cluster <- paletteer::paletteer_d("ggsci::default_igv")
usethis::use_data(palette_cluster, 
                  overwrite = TRUE, 
                  compress = "xz")
palette_heatmap <- paletteer::paletteer_d("MetBrewer::Hiroshige", direction = -1)
usethis::use_data(palette_heatmap, 
                  overwrite = TRUE, 
                  compress = "xz")
