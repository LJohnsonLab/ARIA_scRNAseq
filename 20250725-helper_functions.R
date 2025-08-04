### ggplot decoration for the density plots
density_decor <- function(x){
    list(
        geom_density(),
        theme_minimal(),
        theme(axis.text.y = element_blank(),
              axis.text.x = element_text(size=12),
              axis.title.x = element_text(size=14))
    )
}    
################################################    