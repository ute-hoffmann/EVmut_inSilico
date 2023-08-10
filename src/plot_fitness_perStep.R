library(tidyverse)

replication_data <- read.csv("results/2000replicates_500steps_3trust.csv")

# Fold into rows
#df = df %>%
#  select(-sequences) %>%
#  gather(Data, Value, -accept, -step, -Replicate, -Temperature, -Comparison) %>%
#  # Use only accepted data
#  filter(accept)

df_new <- as_tibble(replication_data)

# For each step, calculate mean value
df_mean_new = bind_rows(lapply(
  0:500, # do this for every step
  function(s){ # define function with argument s
    df_new %>%
      filter(Step <= s) %>%
      group_by(Replicate, Temperature) %>%
      top_n(1, Step) %>%
      group_by(Temperature) %>%
      summarise(Fitness=mean(Fitness), .groups="keep") %>%
      mutate(Step = s) %>%
      ungroup()
  }
))


# Combine data
df_new$Replicate <- as.character(df_new$Replicate)
df_new$Temperature <- as.character(df_new$Temperature)
df_mean_new$Temperature <- as.character(df_mean_new$Temperature)
df_plot = bind_rows(
  df_new %>% mutate(Type = "Replicate"),
  df_mean_new %>% mutate(Replicate = "-1", Type = "Mean")
) %>%
  # Add grouping variable
  mutate(Group = paste(Replicate, Type)) %>%
  # Add T to temperature
  mutate(
    Temperature = factor(
      paste("T =", Temperature),
      levels=paste(
        "T =",
        arrange(., as.numeric(Temperature)) %>%
          pull(Temperature) %>%
          unique()
      )
    )
  )

# Plot it
df_plot$Type <- as.factor(df_plot$Type)
df_plot$Type <- factor(df_plot$Type, c("Mean", "Replicate"))

gp = ggplot(
  df_plot,
  aes(
    group=Type,
    alpha=Type,
    x=Step, y=Fitness#, colour=Temperature
  )
)
gp = gp + geom_line()
gp = gp + facet_grid(1~Temperature, scales="free_y", switch="y")
gp = gp + scale_colour_manual(values=c("#b0636f", "#2b5168"))
gp = gp + scale_alpha_manual(values=c(1,0.2))
gp = gp + theme_bw()
gp = gp + theme(
  axis.ticks=element_line(colour="black"),
  axis.text=element_text(colour="black"),
  strip.background=element_blank(),
  legend.position="top",
  axis.title.y=element_blank(),
  strip.placement = "outside"
)

ggsave("results/2000replicates_500steps_trust3.png", gp, h=12/2.54, w=20/2.54, dpi=200)

