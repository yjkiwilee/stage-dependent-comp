# Define custom themes

theme_jun1 <- function() {
  (theme_linedraw() +
     theme(
       plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
       axis.title = element_text(face = "bold", size = rel(1)),
       legend.position = "bottom",
       legend.box = "vertical",
       legend.spacing.y = unit(0, unit = "pt"),
       legend.box.margin = margin(-9, 0, 0, 0, unit = "pt"),
       legend.margin = margin(0, 0, 0, 0, unit = "mm"),
       strip.background = element_rect(fill = "#eee"),
       strip.text = element_text(color = "black", face = "bold")
     ))
}

theme_jun2 <- function() {
  (theme_linedraw() +
     theme(
       plot.title = element_text(face = "bold", size = rel(1.8), hjust = 0.5),
       axis.title = element_text(size = rel(1.2)),
       legend.position = "right",
       legend.box = "vertical",
       legend.box.margin = margin(-9, 0, 0, 0, unit = "pt"),
       legend.margin = margin(0, 0, 0, 0, unit = "mm"),
       strip.text = element_text(size = rel(1.2)),
       strip.background = element_rect(fill = "#eee"),
       strip.text = element_text(color = "black", face = "bold")
     ))
}
