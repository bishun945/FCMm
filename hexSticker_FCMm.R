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

w = c(1, 3, 5, 6)
# w = 3
dt <- {Rrs_clusters.default[w, ] / Area[w]} %>% 
  tibble::rownames_to_column(., var = "ID") %>%
  reshape2::melt(., id = "ID")
dt$variable %<>% as.vector() %>% as.numeric()
dt$ID %<>% as.factor()

# hcl.pals()
cp.fun <- ggsci::pal_aaas()

ggplot(data = dt, aes(x = variable, y = value, color = ID, group = ID)) + 
  geom_path(size = 1.5, lineend = "round") + 
  scale_color_manual(values = length(w) %>%
                       # FCMm::Spectral()
                       # hcl.colors(., "Cork")
                       cp.fun
                     ) +
  theme_void() + 
  theme(legend.position = "none") +
  theme_transparent() -> p

# p <- ggplot(data = dt, aes(x = variable, y = value, group = ID)) + 
#   geom_path(size = 1.5, lineend = "round", color = "navy") + 
#   theme_void()

p + theme(legend.position = "right") + 
  scale_color_manual(values = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#008280FF", "#631879FF"))

p <- p + scale_color_manual(values = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#008280FF", "#631879FF"))

## use the ggplot2 example
sticker(p, 
        package="FCMm", 
        p_size=28, p_x = 1, p_y = 1.4, 
        p_family = "mont", p_fontface = "bold",
        s_x=1, s_y=.85, s_width=1.3, s_height=0.8,
        h_fill = "#FDFFE1", h_color = "black", p_color = "black", h_size = 1, 
        # h_fill = "#007DFF", h_color = "black", p_color = "#FF8200",
        # h_fill = "#8AAAE5", h_color = "gray75", p_color = "#FFFFFF",
        spotlight = TRUE, l_x = 1, l_y = 0.5, l_alpha = 0.4,
        filename = "./man/figures/logo.png"
        )

sticker(p, 
        package="FCMm", 
        p_size=28, p_x = 1, p_y = 1.4, 
        p_family = "mont", p_fontface = "bold",
        s_x=1, s_y=.85, s_width=1.3, s_height=0.8,
        h_fill = "#FDFFE1", h_color = "black", p_color = "black", h_size = 1, 
        # h_fill = "#007DFF", h_color = "black", p_color = "#FF8200",
        # h_fill = "#8AAAE5", h_color = "gray75", p_color = "#FFFFFF",
        spotlight = TRUE, l_x = 1, l_y = 0.5, l_alpha = 0.4,
        filename = "./man/figures/logo.svg"
)






