complex_summary <- data.frame(
  Complex = c("I", "II", "III", "IV"),
  Total = c(64, 8, 16, 51),
  Found = c(34, 1, 10, 12)
)

complex_summary$Percentage <- round(100 * complex_summary$Found / complex_summary$Total, 1)
library(ggplot2)

plot <- ggplot(complex_summary, aes(x = Complex, y = Percentage)) +
  geom_col(fill = "lightpink") +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5, size = 8, fontface = "bold") +
  labs(title = "Percentage of Mitochondrial Proteins Identified per Complex",
       y = "Percentage Identified", x = "Complex") +
  scale_y_continuous(limits = c(0, 100))+
  theme_classic()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    axis.text = element_text(size = 16, face = "bold", colour = "black"),
    axis.text.x =element_text(size = 16, face = "bold", colour = "black"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16, face = "bold")
  )
plot
ggsave("Complex_I_IV.png", plot = plot , dpi = 600, bg = "transparent", limitsize = F)
