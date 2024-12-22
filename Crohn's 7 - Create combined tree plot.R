# Code for generating figure 6.

library(magick)
library(gridExtra)
graphics.off()

image_list = list()

image_list[[1]] <- image_read("radial_tree 2 CTRL vs TN-CD P.value.png")
#image_list[[1]] <- image_read("radial_tree 2 control vs treatment naive-CD Samples.png")
image_list[[2]] <- image_read("radial_tree 2 control vs treatment naive-CD DEG Count_empirical.png")
# image_list[[2]] <- image_read("radial_tree 2 control vs treatment naive-CD DEGs0.png")
image_list[[3]] <- image_read("radial_tree 2 CTRL vs CD P.value.png")
#image_list[[3]] <- image_read("radial_tree 2 control vs CD Samples.png")
image_list[[4]] <- image_read("radial_tree 2 control vs CD DEG Count_empirical.png")
# image_list[[4]] <- image_read("radial_tree 2 control vs CD DEGs0.png")
image_list[[5]] <- image_read("radial_tree 2 TN-CD vs CD P.value.png")
#image_list[[5]] <- image_read("radial_tree 2 treatment naive-CD vs CD Samples.png")
image_list[[6]] <- image_read("radial_tree 2 treatment naive-CD vs CD DEG Count_empirical.png")
# image_list[[6]] <- image_read("radial_tree 2 treatment naive-CD vs CD DEGs0.png")

location = "0+10"

image_list[[1]] <- image_annotate(image_list[[1]], "(a)          CTRL vs TN-CD", gravity = "northwest", location = location, size = 45, color = "black")
image_list[[2]] <- image_annotate(image_list[[2]], "(d)          CTRL vs TN-CD", gravity = "northwest", location = location, size = 45, color = "black")
image_list[[3]] <- image_annotate(image_list[[3]], "(b)               CTRL vs CD", gravity = "northwest", location = location, size = 45, color = "black")
image_list[[4]] <- image_annotate(image_list[[4]], "(e)               CTRL vs CD", gravity = "northwest", location = location, size = 45, color = "black")
image_list[[5]] <- image_annotate(image_list[[5]], "(c)              TN-CD vs CD", gravity = "northwest", location = location, size = 45, color = "black")
image_list[[6]] <- image_annotate(image_list[[6]], "(f)              TN-CD vs CD", gravity = "northwest", location = location, size = 45, color = "black")


for (i in seq(6)){
  image_list[[i]] = rasterGrob(image_list[[i]])
}

lineGrob <- rectGrob(gp = gpar(col = "grey", fill = "grey"), width = unit(0.01, "npc"))

#grid_layout <- grid.arrange(image_list[[1]], lineGrob, image_list[[2]], image_list[[3]], lineGrob, image_list[[4]], image_list[[5]], lineGrob, image_list[[6]],
#                            widths = c(10, 1, 10), ncol = 3, top = "     Abundance                                              Gene Expression")
grid_layout <- grid.arrange(image_list[[1]], lineGrob, image_list[[2]], image_list[[3]], lineGrob, image_list[[4]], image_list[[5]], lineGrob, image_list[[6]],
                            widths = c(10, 1, 10), ncol = 3, top = "             Abundance                                  Differential Gene Expression")

ggsave("combined_images.png", grid_layout, width=6, height=9.2, dpi = 500)

pp1 <- image_trim(image_read(paste0("combined_images.png"), strip = TRUE))
# pp1 <- image_scale(pp1, "75%")
image_write(pp1, paste0("combined_images.eps"), format = 'eps')
image_write(pp1, paste0("combined_images.pdf"), format = 'pdf')

