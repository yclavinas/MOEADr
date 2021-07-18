rm(list = ls(all = TRUE))

suppressPackageStartupMessages(library(smoof))
suppressPackageStartupMessages(library(MOEADps))
suppressPackageStartupMessages(library(CAISEr))
suppressPackageStartupMessages(library(feather))
# suppressPackageStartupMessages(library(car))
suppressPackageStartupMessages(library(ggplot2))


load("~/Downloads/results.rds")

plot(my.results)




algopairs <- paste(my.results$data.summary$Alg1,
                   my.results$data.summary$Alg2,
                   sep = " - ")


par(mfrow = c(1, 1))
df <- cbind(Comparison = algopairs, my.results$data.summary)

mp <- ggplot(df, aes(x = Comparison, y = Phi, fill = Comparison))
mp +
  geom_violin(alpha = 0.6,
              show.legend = FALSE,
              scale = "width") +
  geom_boxplot(
    alpha = 0,
    show.legend = FALSE,
    outlier.shape = NA,
    width = .15
  ) +
  geom_point(
    shape = 16,
    col = "black",
    fill = "black",
    alpha = 0.6,
    position = position_jitter(width = .15)
  ) +
  geom_abline(
    slope = 0,
    intercept = 0,
    col = "red",
    lty = 2
  ) +
  ylab("Percent difference in IGD") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_minimal()

summary(my.results, test = "wilcoxon")

print(ggplot(df, aes(x = Instance, y = Phi, colour = Comparison,
                     ymin = Phi - SE, ymax = Phi + SE)) + 
        geom_pointrange(show.legend = FALSE) + 
        geom_abline(slope = 0, intercept = 0, col = 1, lty = 2) + 
        facet_grid(Comparison ~ .) + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
        xlab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
)
