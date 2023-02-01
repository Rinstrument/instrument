
library(magick)
## Linking to ImageMagick 6.9.9.14
## Enabled features: cairo, freetype, fftw, ghostscript, lcms, pango, rsvg, webp
## Disabled features: fontconfig, x11
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
library(hexSticker)

img <- image_read("../../../Desktop/theta.png")
res = img %>% 
  image_convert("png") %>% 
  image_resize("1080 x 200")%>% 
  image_fill(color="#2766c5", point="+45")
res

# wrap in plot to preview ie plot(sticker(...))
final_res<-sticker(res, package="theta2", p_size=40,
             p_y = 1.4,
             s_x=1, s_y=0.65, s_width=1.1,
             s_height = 14,
        filename="theta2_icon.png",h_fill="#d8d8d8",h_color = "#0f3266",
        p_color = "#060e1a")

plot(final_res)
