#' ---
#' title: Data Visualization
#' output: github_document
#' ---

#Class 05: Data Visualization
#Today we are going to use ggplot2

#Load package
library(ggplot2)

#Use inbuilt "cars" dataset first, confirm content with head function
head(cars)

#All ggplots have at least data + aes + geom layers
ggplot(data=cars) + 
  aes(x=speed, y=dist) + 
  geom_point() +
  geom_smooth(method="lm") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="Stopping Distance of Cars", x="Speed (MPH)", y="Stopping Distance (ft)")

#ggplot is not the only graphics system, "base" R has one built in too
plot(cars)

#Pull in RNA-seq dataset
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)

#Interrogate features of genes data frame
#Number of genes
nrow(genes)
#Names and amounts of columns
colnames(genes)
ncol(genes)
#Pull numbers of downregulated, unchanging, and upregulated genes
table(genes$State)
#Calculate % of each state
round((table(genes$State)/nrow(genes))*100, 2)

#Generate RNA-seq plot
RNA_plot <- ggplot(genes) +
  aes(x=Condition1, y=Condition2, col=State) +
  geom_point()

#Display plot
RNA_plot

#Adjust colors and add labels
RNA_plot <- RNA_plot + scale_color_manual(values=c("blue", "gray", "red"))
RNA_plot <- RNA_plot + theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="Gene Expression Changes Upon Drug Treatment", 
       x="Control (no drug)", y="Drug Treatment")

#Display new plot
RNA_plot

#Exploring gapminder dataset
library(gapminder)
head(gapminder)

#Plot year vs life expectancy
ggplot(gapminder) +
  aes(x=year, y=lifeExp, col=continent) +
  geom_jitter(width=0.3, alpha=0.4) +
  geom_violin(aes(group=year), alpha=0.2, draw_quantiles=0.5)

#Commenting out plotly for PDF render
#Load plotly
#library(plotly)
#ggplotly()

