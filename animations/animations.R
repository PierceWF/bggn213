#Read in gapminder dataset
url <- "https://raw.githubusercontent.com/jennybc/gapminder/master/inst/extdata/gapminder.tsv"
gapminder <- read.delim(url)

#Filter to single year of data
library(dplyr)
gapminder_2007 <- gapminder %>% filter(year==2007)
gapminder_1992 <- gapminder %>% filter(year==1992)
gapminder_1992_Africa <- gapminder_1992 %>% filter(continent=="Africa")
gapminder_1957 <- gapminder %>% filter(year==1957)

#Plot gapminder_2007
library(ggplot2)
ggplot(gapminder_2007) +
  aes(x=gdpPercap, y=lifeExp, color=continent, size=pop) +
  geom_point(alpha=0.5)

#Gradient for population
ggplot(gapminder_2007) +
  aes(x=gdpPercap, y=lifeExp, color=pop) +
  geom_point(alpha=0.8)

#Size for population
ggplot(gapminder_2007) +
  aes(x=gdpPercap, y=lifeExp, size=pop) +
  geom_point(alpha=0.8) +
  scale_size_area(max_size=10)

#Graph 1957
ggplot(gapminder_1957) +
  aes(x=gdpPercap, y=lifeExp, color=continent, size=pop) +
  geom_point(alpha=0.7) +
  scale_size_area(max_size=10)

#Graph two years
gapminder_1957_2007 <- gapminder %>% filter(year==1957 | year==2007)
ggplot(gapminder_1957_2007) +
  aes(x=gdpPercap, y=lifeExp, color=continent, size=pop) +
  geom_point(alpha=0.7) +
  scale_size_area(max_size=10) +
  facet_wrap(~year)


