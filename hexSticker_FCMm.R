library(hexSticker)
library(showtext)
## Loading Google fonts (http://www.google.com/fonts)
# font_add_google("Gochi Hand", "gochi")
font_add_google("Montserrat", "mont")
## Automatically use showtext to render text for future devices
showtext_auto()

library(FCMm)
library(magrittr)
library(ggplot2)

data("Bi_clusters")
colnames(Rrs_clusters.default) %<>% gsub("X", "", .)
Area <- trapz2(Rrs_clusters.default)

w = c(1, 3, 5, 7, 6)
dt <- {Rrs_clusters.default[w, ] / Area[w]} %>% 
  tibble::rownames_to_column(., var = "ID") %>%
  reshape2::melt(., id = "ID")
dt$variable %<>% as.vector() %>% as.numeric()
dt$ID %<>% as.factor()

# hcl.pals()

ggplot(data = dt, aes(x = variable, y = value, color = ID, group = ID)) + 
  geom_path(size = 2, lineend = "round") + 
  scale_color_manual(values = length(w) %>%
                       # FCMm::Spectral()
                       hcl.colors(., "Cork")
                     ) + 
  theme_void() + 
  theme(legend.position = "none") +
  theme_transparent() -> p

p

## use the ggplot2 example
sticker(p, 
        package="FCMm", 
        p_size=26, p_x = 1, p_y = 1.5, 
        p_family = "mont", p_fontface = "bold",
        s_x=1, s_y=.85, s_width=1.3, s_height=1,
        h_fill = "#1D1B1B", h_color = "gray75", p_color = "#FFFFFF",
        # h_fill = "#8AAAE5", h_color = "gray75", p_color = "#FFFFFF",
        spotlight = TRUE, l_x = 1, l_y = 0.5, l_alpha = 0.4,
        filename = "logo.png"
        )

sticker(p, 
        package="FCMm", 
        p_size=9, p_x = 1, p_y = 1.5, 
        p_family = "mont", p_fontface = "bold",
        s_x=1, s_y=.85, s_width=1.3, s_height=1,
        h_fill = "#1D1B1B", h_color = "gray75", p_color = "#FFFFFF",
        # h_fill = "#8AAAE5", h_color = "gray75", p_color = "#FFFFFF",
        spotlight = TRUE, l_x = 1, l_y = 0.5, l_alpha = 0.4,
        filename = "logo.svg"
)


