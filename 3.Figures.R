library(tidyverse)
library(raster)
library(rgdal)
library(cowplot)
library(RColorBrewer)
library(reshape2)
library(broom)
library(reghelper)
library(modelr)
library(scales)
library(patchwork)
library(quantreg)
library(viridis)
library(ggplotify)
library(rgeos)
library(sf)
library(ggh4x)

#################### Fig. 1) Scaled predictions of SOC depth distributions ####################
soc_preds <- read_rds("2.Output/socconc_param.rds")

range01 <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}

preds <- soc_preds %>%
  unnest(data) %>% 
  mutate(pdepth = ifelse(pred_depth <= SOCz, pred_depth, NA)) %>% 
  filter(!is.na(pdepth)) %>% 
  nest(-id, -EP, -use) %>% 
  mutate(soc_scale = map(data, ~range01(pull(., pred_soc)))) %>% unnest(soc_scale, data) %>% 
  nest(-id, -EP, -use, -beta, -rsq, -pred_depth, -soc_scale)

p1 <- ggplot(preds %>% filter(rsq >= 0.8) %>% arrange(beta), 
             aes(x = soc_scale, y = pred_depth, color = beta))+
  geom_line(aes(group = id))+
  scale_y_reverse("Depth, cm", expand = c(0.01, 0.01))+
  scale_x_continuous("Predicted SOC concentrations (scaled)", breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = c(0, 0.25, 0.50, 0.75, 1.00), expand = c(0.01, 0))+
  scale_color_gradientn(na.value = "transparent",colors = rev(viridis(10)),
                        values = rescale(c(0, 0.02, 0.05, 0.07, 0.1, 0.2, 0.5, 1, 3, 7)),
                        breaks = c(0, 0.02, 0.05, 0.07, 0.1, 0.2, 0.5, 1, 3, 7),
                        
                        limits = c(0, 7),
                        aesthetics = "color",
                        guide = "colorsteps")+
  theme_bw()+
  theme(legend.direction = "horizontal", 
        legend.position = c(0.6, 0.1), text = element_text(color = "black"),
        axis.text = element_text(size = 11, color = "black"))+
  guides(color = guide_colorsteps(title = expression(paste("Estimated ", beta, " values")),
                                  title.position = "top", 
                                  title.hjust = 0.5, ticks = TRUE, 
                                  label.theme = element_text(size = 8,lineheight = 0.8),
                                  barwidth = unit(8, "cm"),
                                  label = TRUE, frame.colour = "black",
                                  frame.linewidth = 0.5, direction = "horizontal",
                                  barheight = unit(0.3, "cm"), show.limits = F))

stocks_preds <- read_rds("2.Output/socstocks_param.rds")

preds_st <- stocks_preds %>% unnest(data) %>% 
  mutate(pdepth = ifelse(pred_depth <= SOCz, pred_depth, NA)) %>% 
  filter(!is.na(pdepth)) %>% 
  nest(-id, -EP, -use) %>% 
  mutate(soc_scale = map(data, ~range01(pull(., pred_soc)))) %>% unnest(soc_scale, data) %>% 
  nest(-id, -EP, -use, -beta, -rsq, -pred_depth, -soc_scale)

p2 <- ggplot(preds_st %>% filter(rsq >= 0.8), aes(x = soc_scale, y = pred_depth, color = beta))+
  geom_path(aes(group = id))+
  scale_y_reverse("Depth, cm", expand = c(0.01, 0.01))+
  scale_x_continuous("Predicted SOC stocks (scaled)", breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = c(0, 0.25, 0.50, 0.75, 1.00), expand = c(0.01, 0))+
  scale_color_gradientn(na.value = "transparent",colors = rev(viridis(10)),
                        values = rescale(c(0, 0.02, 0.05, 0.07, 0.1, 0.2, 0.5, 1, 3, 7)),
                        breaks = c(0, 0.02, 0.05, 0.07, 0.1, 0.2, 0.5, 1, 3),
                        limits = c(0, 7),
                        aesthetics = c("color", "fill"),
                        guide = guide_colorsteps(title = expression(paste("Estimated ", beta, " values")),
                                                 title.position = "top", 
                                                 title.hjust = 0.5, ticks = TRUE, 
                                                 label.theme = element_text(size = 8,lineheight = 0.8),
                                                 barwidth = unit(8, "cm"),
                                                 label = TRUE, frame.colour = "black",
                                                 frame.linewidth = 0.5, direction = "horizontal",
                                                 barheight = unit(0.3, "cm"), show.limits = F))+
  theme_bw()+
  theme(legend.direction = "horizontal", 
        legend.position = c(0.6, 0.1), text = element_text(color = "black"),
        axis.text = element_text(size = 11, color = "black"))


source("https://raw.githubusercontent.com/adrfantini/plot_discrete_cbar/master/plot_discrete_cbar.R")
leg <- plot_discrete_cbar(c(0, 0.02, 0.05, 0.07, 0.1, 0.2, 0.5, 1, 3, Inf), 
                          colors = viridis(9), direction = -1, spacing = "constant",
                          legend_title = expression(paste("Estimated ", beta, " values")),
                          width = 0.8, spacing_scaling = 10, font_size = 3, border_color = "black")
leg

p2 <- (p2 + guides(color = "none")) + annotation_custom(as.grob(leg), xmin = 0.3, xmax = 0.95,
                                                        ymin = -350, ymax = -550)

p <- (p1 + guides(color = "none")) + (p2) +
  plot_annotation(tag_levels = list(c("(a)", "(b)"))) &
  theme(plot.tag = element_text(face = "bold", hjust = 1))

ggsave("3.Figures/Figure1.png", p, width = 12, height = 6)
ggsave("3.Figures/pdfs/Figure1.pdf", p, width = 12, height = 6)
rm(list = ls())

#################### Fig. 2) PEDS SOC concentration data ####################
midwest <- readOGR("1.Data/Shapefiles/USA/midwestst.shp")  
midwest_df <- fortify(midwest)

soilsmw <- read_rds("2.Output/optimsocconcdata.rds") %>% 
  mutate(mxsoc = map_dbl(data, ~pull(., soc)[1]),
         socdeep = map_dbl(data, ~pull(., soc)[length(pull(., soc))]),
         exc = ifelse(mxsoc >= socdeep, "ok", "remove")) %>% 
  filter(exc == "ok", beta < 1000)

soilsmw <- soilsmw %>% mutate(res = ifelse(rsq >= 0.8, "ok", "no"))

epmid <- raster("1.Data/Raster/ep_midwest.tif")
epmid_df <- as.data.frame(epmid, xy = TRUE)

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))


graticule <- shapefile("1.Data/Shapefiles/NaturalEarth/ne_10m_graticules_5.shp")
graticule <- as(gIntersection(graticule, midwest), "SpatialLinesDataFrame")
graticule_df <- fortify(graticule)

plota <- ggplot() + 
  geom_raster(data = epmid_df, aes(x, y, fill = ep_midwest)) +
  scale_fill_gradientn(colors = pull(filter(custompal, round(int, 2) <= 
                                              round(max(epmid_df$ep_midwest, na.rm = T), 2)), colors), 
                       na.value = "transparent",
                       limits = range(pull(filter(custompal, round(int, 2) <= 
                                                    ceiling(round(max(epmid_df$ep_midwest, na.rm = T), 2))), int)), 
                       name = "Effective precipitation",
                       breaks = seq(-4, 0, 0.5),
                       aesthetics = c("fill"),
                       guide = guide_colourbar(available_aes = c("fill")))+
  geom_path(data = graticule_df, aes(long, lat, group = group),
            linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_point(data = soilsmw, aes(long, lat, color = use, shape = res, size = res), stroke = 1.5)+
  scale_shape_manual(values = c("ok" = 19, "no" = 2),
                     labels = expression("R"^2*"< 0.8", "R"^2*">=0.8"), name = NULL)+
  scale_size_manual(values = c("ok" = 3, 'no' = 5), guide = NULL)+
  geom_polygon(data = midwest_df, aes(x = long, y = lat, group = group), 
               color = "black", alpha = 0, linewidth = 0.25)+
  scale_color_manual(values = c("N" = "darkgreen", "A" = "darkorange"),
                     labels = c("N" = "Native prairie", "A" = "Agriculture"), name = NULL)+
  theme_classic(base_size = 12) +
  scale_y_continuous(breaks = seq(35, 50, 5), labels = paste0(seq(35, 50, 5), "\u00B0N"))+
  scale_x_continuous(breaks = seq(-105, -80, 5), labels = paste0(seq(105, 80, -5), "\u00B0W"))+
  coord_fixed(ratio = 1.35)+
  # coord_sf(label_graticule = "SW", datum = "WGS84", 
  #          crs = crs(projection(midwest))) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(color = "black", size = 17),
        axis.text.y = element_text(color = "black", size = 17),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        plot.title = element_text(size = 15),
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        legend.text.align = 0,
        legend.box.margin = margin(-5, -5, -5, -5),
        legend.box.spacing = unit(0.2, "cm"),
        legend.box.just = "center",
        legend.position ="bottom", legend.direction = "horizontal", 
        legend.box = "horizontal",
        legend.key.width = unit(2.5, "cm"), 
        legend.spacing.y = unit(0.3, "cm"), legend.spacing.x = unit(0.05, "cm"),
        legend.title = element_text(size = 19),
        legend.text = element_text(size = 19, margin = margin(t = -0.2, l = -0.5, unit = "lines"))) +
  guides(fill = guide_colorbar(title.position = "top",
                               even.steps = FALSE,
                               show.limits = FALSE, ticks.colour = "black",
                               title.hjust = 0.5, ticks = TRUE,
                               label.theme = element_text(size = 17, lineheight = 0.8),
                               label = TRUE, frame.colour = "black", order = 3),
         color = guide_legend(ncol = 1, order = 2, reverse = TRUE, 
                              label.position = "right", override.aes = list(size = 6)),
         shape = guide_legend(ncol = 1, order = 1, override.aes = list(size = 6)))+
  annotate(geom = "text", x = -105, y = 50, label = "(a)", fontface = 2, size = 9)

molls <- soilsmw %>% unnest(data) %>% 
  nest(-id, -lat, -long, -EP, -use, -beta, -rsq, -beta_sd) %>% 
  mutate(k = log(beta),
         sd = log(beta + beta_sd) - log(beta))

model <- rq(k ~ EP*use, data = molls %>% filter(rsq >= 0.8, beta < 1000) %>% 
              transform(use = factor(use, levels = c("A", "N")))) %>% 
  summary.rq(iid = T, se = "iid")
model$coefficients %>% round(3)
cor(molls$k[molls$use == "A"], molls$EP[molls$use == "A"], method = "spearman") %>% round(3)

pb <- ggplot(molls %>% filter(use == "A", rsq >= 0.8, beta < 1000), 
             aes(x = EP, y = k, fill = use, color = use))+
  facet_grid(cols = vars(use), scales = "free_y",
             labeller = labeller(use = c("A" = "Agriculture",
                                         "N" = "Native")))+
  geom_point(size = 5, pch = 21, stroke = 0.5, alpha = 0.5)+
  geom_quantile(quantiles = 0.5, aes(group = use), lwd = 1.5, color = "black")+
  scale_fill_manual(values = c("N" = "darkgreen", "A" = "darkorange"),
                    labels = c("N" = "Native prairie", "A" = "Agriculture"), 
                    aesthetics = c("fill", "colour"))+
  scale_x_continuous("Effective precipitation", limits = c(-3.5, 0))+
  scale_y_continuous(expression(paste("ln(", beta ,"), cm"^-1*"")), 
                     breaks = seq(-12, 10, 2))+
  coord_cartesian(ylim = c(-6, 5))+
  theme_bw()+
  theme(text = element_text(size = 19, color = "black"), panel.grid = element_blank(),
        axis.text = element_text(size = 17, color = "black"), 
        plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm"), 
        strip.text.x = element_text(size = 19),
        legend.position = "none", plot.title = element_blank()) +
  guides(fill = "none", color = "none")+
  annotate(geom = "text", x = -3.5, y = 5, label = "(b)", fontface = 2, size = 9)+
  annotate(geom = "text", x = -2.5, y = 4, 
           label = expression(paste(rho, "= 0.029, ", italic("p"), " < 0.001")), size = 7)+
  annotate(geom = "text", x = -2.5, y = 3.3, 
           label = expression(paste(Delta, "ln(", beta, ") = 0.084 \u00B1 0.021 cm"^-1*"")),
           size = 7)+
  annotate(geom = "text", x = -2.5, y = 2.55, size = 7,
           label = expression(paste(italic("n"), " = 1923 profiles")))

model <- rq(k ~ EP*use, data = molls %>% filter(rsq >= 0.8, beta < 1000) %>% 
              transform(use = factor(use, levels = c("N", "A")))) %>% 
  summary.rq(iid = T, se = "iid")
model$coefficients %>% round(3)
cor(molls$k[molls$use == "N"], molls$EP[molls$use == "N"], method = "spearman") %>% round(3)

pc <- ggplot(molls %>% filter(use == "N", rsq >= 0.8, beta < 1000), 
             aes(x = EP, y = k, fill = use, color = use))+
  facet_grid(cols = vars(use), scales = "free_y",
             labeller = labeller(use = c("A" = "Agriculture",
                                         "N" = "Native prairie")))+
  geom_point(size = 5, pch = 21, stroke = 0.5, alpha = 0.5)+
  geom_quantile(quantiles = 0.5, aes(group = use), lwd = 1.5, color = "black")+
  scale_fill_manual(values = c("N" = "darkgreen", "A" = "darkorange"),
                    labels = c("N" = "Native prairie", "A" = "Agriculture"), 
                    aesthetics = c("fill", "colour"))+
  scale_x_continuous("Effective precipitation")+
  scale_y_continuous(expression(paste("ln(", beta ,"), cm"^-1*"")), breaks = seq(-12, 10, 2))+
  coord_cartesian(ylim = c(-6, 3))+
  theme_bw()+
  theme(text = element_text(size = 19, color = "black"), panel.grid = element_blank(),
        axis.text = element_text(size = 17, color = "black"), 
        plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm"), 
        strip.text.x = element_text(size = 19),
        legend.position = "none", plot.title = element_blank()) +
  guides(fill = "none", color = "none")+
  annotate(geom = "text", x = -3.5, y = 3, label = "(c)", fontface = 2, size = 9)+
  annotate(geom = "text", x = -2.5, y = 2, 
           label = expression(paste(rho, "= -0.244, ", italic("p"), " < 0.001")), size = 7)+
  annotate(geom = "text", x = -2.5, y = 1.5, size = 7,
           label = expression(paste(Delta, "ln(", beta, ") = -0.245 \u00B1 0.030 cm"^-1*"")))+
  annotate(geom = "text", x = -2.5, y = 0.9, size = 7,
           label = expression(paste(italic("n"), " = 683 profiles")))

fig2 <- plota / (pb + pc) + 
  plot_layout(ncol = 1, nrow = 2, widths = c(1), heights = c(2, 1.5))

ggsave("3.Figures/Figure2.png", fig2, dpi = 500, height = 15, width = 15)
ggsave("3.Figures/pdfs/Figure2.pdf", fig2, dpi = 500, height = 15, width = 15)
rm(list = ls())


#################### Fig. 3) PEDS SOC stocks data ####################
midwest <- readOGR("1.Data/Shapefiles/USA/midwestst.shp")  
midwest_df <- fortify(midwest)

soilsbd <- read_rds("2.Output/optimsocstockdata.rds") 

soilsbd <- soilsbd %>% unnest(data) %>% 
  nest(-id, -lat, -long, -EP, -use, -beta, -beta_sd, -rsq) %>% 
  mutate(res = ifelse(rsq >= 0.8, "ok", "no"))

epmid <- raster("1.Data/Raster/ep_midwest.tif")
epmid_df <- as.data.frame(epmid, xy = TRUE)

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))


graticule <- shapefile("1.Data/Shapefiles/NaturalEarth/ne_10m_graticules_5.shp")
graticule <- as(gIntersection(graticule, midwest), "SpatialLinesDataFrame")
graticule_df <- fortify(graticule)

plota <- ggplot() + 
  geom_raster(data = epmid_df, aes(x, y, fill = ep_midwest)) +
  scale_fill_gradientn(colors = pull(filter(custompal, round(int, 2) <= 
                                              round(max(epmid_df$ep_midwest, na.rm = T), 2)), colors), 
                       na.value = "transparent",
                       limits = range(pull(filter(custompal, round(int, 2) <= 
                                                    ceiling(round(max(epmid_df$ep_midwest, na.rm = T), 2))), int)), 
                       name = "Effective precipitation",
                       breaks = seq(-4, 0, 0.5),
                       aesthetics = c("fill"),
                       guide = guide_colourbar(available_aes = c("fill")))+
  geom_path(data = graticule_df, aes(long, lat, group = group),
            linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_point(data = soilsbd, aes(long, lat, color = use, shape = res, size = res), stroke = 1.5)+
  scale_shape_manual(values = c("ok" = 19, "no" = 2),
                     labels = expression("R"^2*"< 0.8", "R"^2*">=0.8",), name = NULL)+
  scale_size_manual(values = c("ok" = 3, 'no' = 5), guide = NULL)+
  geom_polygon(data = midwest_df, aes(x = long, y = lat, group = group), 
               color = "black", alpha = 0, linewidth = 0.25)+
  scale_color_manual(values = c("N" = "darkgreen", "A" = "darkorange"),
                     labels = c("N" = "Native prairie", "A" = "Agriculture"), name = NULL)+
  theme_classic(base_size = 12) +
  scale_y_continuous(breaks = seq(35, 50, 5), labels = paste0(seq(35, 50, 5), "\u00B0N"))+
  scale_x_continuous(breaks = seq(-105, -80, 5), labels = paste0(seq(105, 80, -5), "\u00B0W"))+
  coord_fixed(ratio = 1.35)+
  # coord_sf(label_graticule = "SW", datum = "WGS84", 
  #          crs = crs(projection(midwest))) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(color = "black", size = 17),
        axis.text.y = element_text(color = "black", size = 17),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        plot.title = element_text(size = 15),
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        legend.text.align = 0,
        legend.box.margin = margin(-5, -5, -5, -5),
        legend.box.spacing = unit(0.2, "cm"),
        legend.box.just = "center",
        legend.position ="bottom", legend.direction = "horizontal", legend.box = "horizontal",
        legend.key.width = unit(2.5, "cm"),
        legend.spacing.y = unit(0.3, "cm"), legend.spacing.x = unit(0.1, "cm"),
        legend.title = element_text(size = 19),
        legend.text = element_text(size = 19, margin = margin(t = -0.2, l = -0.5, unit = "lines"))) +
  guides(fill = guide_colorbar(title.position = "top",
                               even.steps = FALSE,
                               show.limits = FALSE, ticks.colour = "black",
                               title.hjust = 0.5, ticks = TRUE,
                               label.theme = element_text(size = 17, lineheight = 0.8),
                               label = TRUE, frame.colour = "black", order = 3),
         color = guide_legend(ncol = 1, order = 2, reverse = TRUE, 
                              label.position = "right", override.aes = list(size = 6)),
         shape = guide_legend(ncol = 1, order = 1, override.aes = list(size = 6)))+
  annotate(geom = "text", x = -105, y = 50, label = "(a)", fontface = 2, size = 9)

stockterms <- soilsbd %>% mutate(k = log(beta), sd = log(beta + beta_sd) - log(beta))
model <- rq(k ~ EP*use, data = stockterms %>% filter(rsq >= 0.8) %>% 
              transform(use = factor(use, levels = c("N", "A")))) %>% 
  summary.rq(iid = T, se = "iid")
model$coefficients %>% round(3)
cor(stockterms$k[stockterms$use == "A"], stockterms$EP[stockterms$use == "A"], method = "spearman") %>% round(3)

pb <- ggplot(stockterms %>% filter(use == "A", rsq >= 0.8), 
             aes(x = EP, y = k, fill = use, color = use))+
  facet_grid(cols = vars(use), scales = "free_y",
             labeller = labeller(use = c("A" = "Agriculture",
                                         "N" = "Native prairie")))+
  geom_point(size = 5, pch = 21, stroke = 0.5, alpha = 0.5)+
  scale_fill_manual(values = c("N" = "darkgreen", "A" = "darkorange"),
                    labels = c("N" = "Native prairie", "A" = "Agriculture"), 
                    aesthetics = c("fill", "colour"))+
  scale_x_continuous("Effective precipitation", limits = c(-3.5, 0))+
  scale_y_continuous(expression(paste("ln(", beta ,"), cm"^-1*"")), 
                     breaks = seq(-12, 10, 2))+
  coord_cartesian(ylim = c(-6, 3))+
  theme_bw()+
  theme(text = element_text(size = 19, color = "black"), panel.grid = element_blank(),
        axis.text = element_text(size = 17, color = "black"), 
        plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm"), 
        strip.text.x = element_text(size = 19),
        legend.position = "none", plot.title = element_blank()) +
  guides(fill = "none", color = "none")+
  annotate(geom = "text", x = -3.5, y = 3, label = "(b)", fontface = 2, size = 9)+
  annotate(geom = "text", x = -1.5, y = 2.5, 
           label = expression(paste(rho,  "= 0.073, ", italic("p"), " = 0.066")), size = 7)+
  annotate(geom = "text", x = -1.5, y = 2, 
           label = expression(paste(Delta, "ln(", beta, ") = 0.080 \u00B1 0.043 cm"^-1*"")),
           size = 7)+
  annotate(geom = "text", x = -1.5, y = 1.3, size = 7,
           label = expression(paste(italic("n"), " = 352 profiles")))

model <- rq(k ~ EP*use, data = stockterms %>% filter(rsq >= 0.8) %>% 
              transform(use = factor(use, levels = c("N", "A")))) %>% 
  summary.rq(iid = T, se = "iid")
model$coefficients %>% round(3)
cor(stockterms$k[stockterms$use == "N"], stockterms$EP[stockterms$use == "N"], method = "spearman") %>% round(3)

pc <- ggplot(stockterms %>% filter(use == "N", rsq >= 0.8), 
             aes(x = EP, y = k, fill = use, color = use))+
  facet_grid(cols = vars(use), scales = "free_y",
             labeller = labeller(use = c("A" = "Agriculture",
                                         "N" = "Native prairie")))+
  geom_point(size = 5, pch = 21, stroke = 0.5, alpha = 0.5)+
  geom_quantile(quantiles = 0.5, aes(group = use), lwd = 1.5, color = "black")+
  scale_fill_manual(values = c("N" = "darkgreen", "A" = "darkorange"),
                    labels = c("N" = "Native prairie", "A" = "Agriculture"), 
                    aesthetics = c("fill", "colour"))+
  scale_x_continuous("Effective precipitation")+
  scale_y_continuous(expression(paste("ln(", beta ,"), cm"^-1*"")), breaks = seq(-12, 0, 2))+
  coord_cartesian(ylim = c(-6, 0))+
  theme_bw()+
  theme(text = element_text(size = 19, color = "black"), panel.grid = element_blank(),
        axis.text = element_text(size = 17, color = "black"), 
        plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm"), 
        strip.text.x = element_text(size = 19),
        legend.position = "none", plot.title = element_blank()) +
  guides(fill = "none", color = "none")+
  annotate(geom = "text", x = -3.5, y = 0, label = "(c)", fontface = 2, size = 9)+
  annotate(geom = "text", x = -2, y = -0.2, size = 7,
           label = expression(paste(rho, "= -0.234, ", italic("p"), " < 0.001")))+
  annotate(geom = "text", x = -2, y = -0.6, size = 7,
           label = expression(paste(Delta, "ln(", beta, ") = -0.257 \u00B1 0.074 cm"^-1*"")))+
  annotate(geom = "text", x = -2, y = -1, size = 7,
           label = expression(paste(italic("n"), " = 103 profiles")))

fig3 <- plota / (pb + pc) + 
  plot_layout(ncol = 1, nrow = 2, widths = c(2), heights = c(2, 1.5))

ggsave("3.Figures/Figure3.png", fig3, dpi = 500, height = 15, width = 15)
ggsave("3.Figures/pdfs/Figure3.pdf", fig3, dpi = 500, height = 15, width = 15)
rm(list = ls())


#################### Fig. 4) Comparison of PEDS betas from SOC concentration and stocks ####################
concterms <- read_rds("2.Output/optimsocconcdata.rds") %>% unnest(data) %>% 
  nest(-id, -use, -EP, -rsq, -beta, -beta_sd, -beta_se) %>% 
  mutate(k = log(beta), kse = log(beta + beta_sd) - log(beta))
stockterms <- read_rds("2.Output/optimsocstockdata.rds") %>% unnest(data) %>% 
  nest(-id, -use, -EP, -rsq, -beta, -beta_sd, -beta_se) %>% 
  mutate(k = log(beta), kse = log(beta + beta_sd) - log(beta))

fits <- left_join(stockterms, concterms, 
                  by = c("id", "EP", "use"), suffix = c("stock", "conc")) %>% 
  filter(rsqconc >= 0.8, rsqstock >= 0.8)

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))

cor.test(fits$kconc, fits$kstock, method = "spearman", exact = FALSE)

f4a <- ggplot(fits %>% filter(rsqconc >= 0.8, rsqstock >= 0.8), 
              aes(x = kconc, y = kstock, fill = EP, shape = use))+
  geom_hline(aes(yintercept = 0), color = "gray90")+
  geom_vline(aes(xintercept = 0,), color = "gray90")+
  geom_abline(slope = 1, intercept = 0, color = "gray25") +
  geom_point(size = 7)+
  scale_fill_gradientn(colors = custompal$colors, 
                       na.value = "transparent",
                       limits = range(custompal$int), 
                       name = "Effective precipitation",
                       breaks = seq(-4, 0, 0.5),
                       aesthetics = c("fill"),
                       guide = guide_colourbar(available_aes = c("fill")))+
  scale_shape_manual(values = c("N" = 21, "A" = 24),
                     labels = c("N" = "Native prairie", "A" = "Agriculture"), name = NULL)+
  scale_x_continuous(expression(paste("ln(", beta ,") from SOC concentration, cm"^-1*"")), 
                     breaks = seq(-20, 0, 2))+
  scale_y_continuous(expression(paste("ln(", beta ,") from SOC stocks, cm"^-1*"")), 
                     breaks = seq(-20, 4, 2))+
  coord_cartesian(ylim = c(-6, 3), xlim = c(-5.5, 2))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.line = element_blank(),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 15),
        plot.title = element_text(size = 15),
        plot.margin = unit(c(1, 1, 1, 1), "lines"), 
        legend.box.margin = margin(-5, -5, -5, -5), 
        legend.box.spacing = unit(0.2, "cm"), 
        legend.position = c(0.82, 0.15), legend.box = "vertical", legend.box.just = "right",
        legend.direction = "horizontal",
        legend.key.width = unit(1.5, "cm"),
        legend.spacing.y = unit(0.08, "cm"), legend.spacing.x = unit(0.5, "cm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11, margin = margin(t = -0.2, r = -2.5, unit = "lines"))) +
  guides(fill = guide_colorbar(title.position = "top",
                               even.steps = FALSE,
                               show.limits = FALSE, ticks.colour = "black",
                               title.hjust = 0.5, ticks = TRUE,
                               label.theme = element_text(size = 11, lineheight = 0.8),
                               label = TRUE, frame.colour = "black", order = 2),
         shape = guide_legend(ncol = 2, order = 1, label.position = "left", reverse = T, override.aes = list(size = 7)))+
  annotate(geom = "text", x = -5.7, y = 3, label = "(a)", fontface = 2, size = 6)+
  annotate(geom = "text", x = -5, y = 1.5, size = 5,
           label = expression(paste(rho, " = 0.926")))+
  annotate(geom = "text", x = -5, y = 1, size = 5,
           label = expression(paste(italic("p"), " < 0.001")))+
  annotate(geom = "text", x = -5, y = 0.5, size = 5,
           label = expression(paste(italic("n"), " = 445 profiles")))

profiles <- readRDS("1.Data/soils_midwest.rds")

native <- profiles %>% filter(id %in% c("93P0138")) %>% 
  unnest(data)
ag <- profiles %>% filter(id %in% c("40A1947")) %>% unnest(data)
profs <- bind_rows(native, ag) %>% mutate(socstock = soc*bd)

socpred <- read_rds( "2.Output/socconc_param.rds")
conc <- socpred %>% filter(id %in% c("93P0138", "40A1947"))

bdsocpred <- read_rds("2.Output/socstocks_param.rds")
bd <- bdsocpred %>% filter(id %in% c("93P0138", "40A1947"))

f4b <- ggplot(profs %>% filter(use == "A") %>% melt(id.vars = c("id", "lat", "long", "EP", "use", "depth"), 
                                                    measure.vars = c("soc", "socstock"), 
                                                    variable.name = "soccal", value.name = "soc"), 
              aes(x = soc, y = depth, shape = soccal))+
  geom_path(data = conc %>% filter(use == "A") %>% unnest(data) %>% mutate(soccal = "soc"),
            aes(x = pred_soc, y = pred_depth), color = "darkorange", show.legend = T)+
  geom_path(data =  bd %>% filter(use == "A") %>% unnest(data) %>% mutate(soccal = "socstock"),
            aes(x = pred_soc, y = pred_depth), color = "darkorange", lty = 2, show.legend = T)+
  scale_shape_manual("test", values = c("soc" = 22, "socstock" = 23),
                     labels = c("soc" = "SOC concentration, %", 
                                "socstock" = expression(paste("SOC stock, gC cm"^-3*""["soil"]*""))))+
  geom_point(size = 6, fill = "darkorange", color = "black")+
  scale_y_reverse("Depth, cm", breaks = seq(0, 600, 50), limits = c(240, 0))+
  scale_x_continuous("Soil organic carbon, %", limits = c(0, 4), position = "top",
                     sec.axis = sec_axis(trans = ~ . * 2, name = expression("Soil organic carbon, gC cm"^-3*""["soil"])))+
  coord_cartesian(xlim = c(0, NA))+
  
  annotate(geom = "text", x = 1.5, y = 60, label = "Soil profile: 40A1947",
           hjust = 0, fontface = 2,  size = 5,)+
  annotate(geom = "text", x = 1.5, y = 75,  size = 5,
           label = paste("Effective precipitation = -0.34"), hjust = 0)+
  
  annotate(geom = "text", x = 1.5, y = 90, hjust = 0, size = 5,
           label = expression(paste("R"^2*""["SOC, %"]*" = 0.99")))+
  annotate(geom = "text", x = 1.5, y = 105,  size = 5,
           label = expression(paste("ln(", beta, ")"["SOC, %"]*"= ", "-2.82 cm"^-1*"")), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 120,  size = 5,
           label = expression(paste("SOC"["i"]*" = 1.87%")), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 135,  size = 5,
           label = expression(paste("SOC"["f"]*" = 0.20%")), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 150,  size = 5,
           label = expression(paste("z"["SOC, %"]*" = 206 cm")), hjust = 0)+
  
  annotate(geom = "text", x = 1.5, y = 165, hjust = 0, size = 5,
           label = expression(paste("R"^2*""["SOC, gC cm"^-3*""["soil"]*""]*" = 0.99")))+
  annotate(geom = "text", x = 1.5, y = 180,  size = 5,
           label = expression(paste("ln(", beta, ")"["SOC, gC cm"^-3*""["soil"]*""]*"= ", "-2.80 cm"^-1*"")), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 195,  size = 5,
           label = expression(paste("SOC"["i"]*" = 3.02 gC cm"^-3*""["soil"])), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 210,  size = 5,
           label = expression(paste("SOC"["f"]*" = 0.37 gC cm"^-3*""["soil"])), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 225,  size = 5,
           label = expression(paste("z"["SOC, gC cm"^-3*""["soil"]*""]*" = 211 cm")), hjust = 0)+
  ggtitle("Agriculture")+
  annotate(geom = "text", x = 0, y = 0, label = "(b)", fontface = 2, size = 6)+
  theme_bw() +
  theme(strip.text.x = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 14), 
        plot.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 14))+
  scale_linetype_manual("test", values = c(1, 2),
                        labels = c("SOC concentration, %", 
                                   expression(paste("SOC stock, gC cm"^-3*""["soil"]*""))))+
  guides(shape = guide_legend(ncol = 2, order = 1, label.position = "right", 
                              override.aes = list(fill = "black", color = "black", 
                                                  linetype = c(1, 2)),
                              keywidth = unit(1.8, "cm")))

f4c <- ggplot(profs %>% filter(use == "N"), aes(x = soc, y = depth))+
  geom_path(data = conc %>% filter(use == "N") %>% unnest(data), 
            aes(x = pred_soc, y = pred_depth), color = "darkgreen")+
  geom_path(data = bd %>% filter(use == "N") %>% unnest(data), 
            aes(x = pred_soc, y = pred_depth), color = "darkgreen", lty = 2)+
  geom_point(size = 6, fill = "darkgreen", color = "black", pch = 22)+
  geom_point(data = profs %>% filter(use == "N"), aes(x = socstock, y = depth),
             size = 6, fill = "darkgreen", color = "black", pch = 23)+
  scale_y_reverse("Depth, cm", breaks = seq(0, 600, 50), limits = c(240, 0))+
  scale_x_continuous("Soil organic carbon, %", limits = c(0, 4), position = "top",
                     sec.axis = sec_axis(trans = ~ . * 2, name = expression("Soil organic carbon, gC cm"^-3*""["soil"])))+
  coord_cartesian(xlim = c(0, NA))+
  
  annotate(geom = "text", x = 1.5, y = 60, size = 5,
           label = "Soil profile: 93P0138", hjust = 0, fontface = 2)+
  annotate(geom = "text", x = 1.5, y = 75,  size = 5,
           label = paste("Effective precipitation = -0.46"), hjust = 0)+
  
  annotate(geom = "text", x = 1.5, y = 90, size = 5, hjust = 0,
           label = expression(paste("R"^2*""["SOC, %"]*" = 0.99")))+
  annotate(geom = "text", x = 1.5, y = 105, size = 5,
           label = expression(paste("ln(", beta, ")"["SOC, %"]*"= ", "-2.66 cm"^-1*"")), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 120, size = 5,
           label = expression(paste("SOC"["i"]*" = 5.86%")), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 135, size = 5,
           label = expression(paste("SOC"["f"]*" = 0.07%")), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 150, size = 5,
           label = expression(paste("z"["SOC, %"]*" = 197 cm")), hjust = 0)+
  
  annotate(geom = "text", x = 1.5, y = 165, size = 5, hjust = 0,
           label = expression(paste("R"^2*""["SOC, gC cm"^-3*""["soil"]*""]*" = 0.99")))+
  annotate(geom = "text", x = 1.5, y = 180, size = 5,
           label = expression(paste("ln(", beta, ")"["SOC, gC cm"^-3*""["soil"]*""]*"= ", "-2.80 cm"^-1*"")), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 195, size = 5,
           label = expression(paste("SOC"["i"]*" = 6.76 gC cm"^-3*""["soil"])), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 210, size = 5,
           label = expression(paste("SOC"["f"]*" = 0.11 gC cm"^-3*""["soil"])), hjust = 0)+
  annotate(geom = "text", x = 1.5, y = 225, size = 5,
           label = expression(paste("z"["SOC, gC cm"^-3*""["soil"]*""]*" = 239 cm")), hjust = 0)+
  ggtitle("Native prairie")+
  annotate(geom = "text", x = 0, y = 0, label = "(c)", fontface = 2, size = 6)+
  theme_bw() +
  theme(strip.text.x = element_blank(), panel.grid = element_blank(), 
        legend.position = "none", axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 14), plot.title = element_text(face = "bold", size = 16))

fig4 <- (f4a + plot_layout(guides = "keep")) / ((f4b + f4c) + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + 
  plot_layout(ncol = 1, nrow = 2, widths = c(0.6), heights = c(0.8, 0.8))

ggsave("3.Figures/Figure4.png", fig4, width = 12, height = 12, dpi = 450)  
ggsave("3.Figures/pdfs/Figure4.pdf", fig4, width = 12, height = 12, dpi = 450)
rm(list = ls())

#################### Fig. 5) Zsoc from PEDS ####################
zconc <- read_rds("2.Output/socconc_param.rds") %>% dplyr::select(-data) %>% 
  mutate(class = "conc")
zbd <- read_rds("2.Output/socstocks_param.rds") %>% dplyr::select(-data) %>% 
  mutate(class = "stock")

zsoc <- bind_rows(zconc, zbd) %>% 
  transform(class = factor(class, levels = c("conc", "stock"),
                           labels = c(expression(paste("SOC concentrations, %")),
                                      expression(paste("SOC stocks, gC cm"^-3*""["soil"]*"")))))

avgs <- read_csv("2.Output/resamplingttest_peds.csv") %>% 
  transform(class = factor(class, levels = c("conc", "stock"),
                           labels = c(expression(paste("SOC concentrations, %")),
                                      expression(paste("SOC stocks, gC cm"^-3*""["soil"]*"")))))


ggplot(zsoc %>% filter(rsq >= 0.8, SOCz < 500), aes(x = SOCz, color = use))+
  facet_wrap(vars(class),
             labeller = label_parsed, scales = "free_y", shrink = T)+
  geom_histogram(aes(fill = use), lwd = 1, show.legend = T, binwidth = 10)+
  geom_vline(data = avgs, aes(xintercept = avg, color = use), lty = 2)+
  geom_text(data = avgs %>% 
              mutate(label = paste0(round(avg, 1), " cm"),
                     y = c(55, 55, 13.5, 13.5)),
            aes(x = avg-12, y = y, label = label, color = use), angle = 90, fontface = 2)+
  scale_color_manual(values = c("N" = "darkgreen", "A" = "darkorange"),
                     labels = c("N" = "Native", "A" = "Agriculture"),
                     aesthetics = c("color", "fill"))+
  theme_bw()+
  theme(legend.position = c(0.95, 0.94), legend.title = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10), 
        axis.ticks.x.bottom = element_line(colour = "black"),
        text = element_text(size = 14, color = "black"), 
        panel.grid = element_blank(), legend.key.size = unit(0.5, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        axis.text = element_text(size = 11, color = "black"), 
        legend.background = element_rect(fill = "transparent"),
        axis.line = element_line(), axis.title = element_text(size = 13))+
  guides(colour = guide_legend(label.position = "left", size = 1, override.aes = list(size = 0.5)))+
  scale_y_continuous("Number of profiles", limits = c(0, NA), 
                     expand = c(0, 0, 0.05, 0.05))+
  scale_x_continuous(expression(paste(italic("z"), ""["SOC"]*", cm")), 
                     expand = c(0.025, 0.025))+
  geom_text(data = zsoc %>% nest(-class) %>%
              mutate(label = paste0("(", letters[1:2], ")" )),
            aes(label = label, x = 2, y = c(60, 15)), 
            color = "black", size = 5, fontface = 2)+
  geom_text(data = zsoc %>% nest(-class) %>% 
              mutate(x = c(50, 50), y = c(50, 12.5), use = "N"), size = 3.5,
            aes(x = x, y = y), label = expression(paste(italic("p"), " < 0.001")), parse = T,
            color = "black")

ggsave("3.Figures/Figure5.png", last_plot(), dpi = 500, height = 5, width = 10)
ggsave("3.Figures/pdfs/Figure5.pdf", last_plot(), dpi = 500, height = 5, width = 10)
rm(list = ls())


#################### Fig. 6) SOC across Kansas ####################
ep <- raster("1.Data/Raster/ep.tif")
ep_df <- as.data.frame(ep, xy = TRUE)

sites <- read_csv("1.Data/samplingsites.csv")
sites <- sites %>% mutate(Use = str_extract(plot, "^^\\w{1}")) %>% 
  dplyr::select(Site, Use, lat, long, soil_type) %>% group_by(Site, Use) %>% 
  summarise(lat = mean(lat), long = mean(long)) %>% 
  arrange(long) %>% nest(-Site) %>% 
  add_column(long_i = c(-101.5, -101.1, -101, -99.9, -99.5, -99.1, 
                        -98, -97.2, -96.4, -96, -96.3, -95.5), 
             lat_i = c(39.3, 38.2, 39.7, 39.5, 38.3, 39.7,
                       39.5, 38.3, 39.5, 39, 38.2, 39.8)) %>% 
  mutate(Site_lab = (ifelse(Site == "EKS" | Site == "KNZ" | Site == "HAY" | Site == "TRB",
                            paste0("bold(underline(paste(", Site, ")))"), 
                            paste0("paste(", Site, ")"))))

ks <- readOGR("1.Data/Shapefiles/Kansas/Kansas.shp")
ks_df <- fortify(ks)

graticule <- shapefile("1.Data/Shapefiles/NaturalEarth/ne_10m_graticules_1.shp")
graticule <- crop(graticule, ks)
graticule_df <- fortify(graticule)

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))

fig6a <- ggplot() + 
  geom_raster(data = ep_df, aes(x, y, fill = ep)) +
  scale_fill_gradientn(colors = pull(filter(custompal, round(int, 2) <= 
                                              round(max(ep_df$ep, na.rm = T), 2)), colors), 
                       na.value = "transparent",
                       limits = range(pull(filter(custompal, round(int, 2) <= 
                                                    round(max(ep_df$ep, na.rm = T), 2)), int)), 
                       name = "Effective precipitation",
                       breaks = seq(-4, 0, 0.5),
                       aesthetics = c("fill"),
                       guide = guide_colourbar(available_aes = c("fill")))+
  geom_path(data = graticule_df, aes(long, lat, group = group),
            linetype = "dashed", color = "black", size = 0.1) +
  geom_polygon(data = ks_df, aes(x = long, y = lat, group = group), 
               color = "black", alpha = 0, size = 0.25)+
  geom_point(data = sites %>% filter(!(Site %in% c("TRB", "SVR", "TRG", "HAY", "HNC", "KNZ"))) %>% unnest(data), 
             aes(long, lat, color = Use), pch = 19, size = 5, stroke = 0.5)+
  geom_point(data = sites %>% filter(Site %in% c("TRB", "SVR", "TRG", "HAY", "HNC", "KNZ")) %>% unnest(data), 
             aes(long, lat, color = Use), pch = 19, size = 5, stroke = 0.5, 
             position = position_dodge2(width = 0.2))+
  scale_color_manual(values = c("N" = "darkgreen", "A" = "darkorange"),
                     labels = c("N" = "Native prairie", "A" = "Agriculture"), name = NULL)+
  geom_text(data = sites %>% filter(!(Site %in% c("TRB", "HAY", "KNZ", "EKS"))), size = 3,
            aes(x = long_i, y = lat_i, label = Site_lab), parse = TRUE,
            check_overlap = TRUE, fontface = 2)+
  geom_text(data = sites %>% filter(Site %in% c("TRB", "HAY", "KNZ", "EKS")), size = 3.5,
            aes(x = long_i, y = lat_i, label = Site_lab), parse = TRUE,
            check_overlap = TRUE, fontface = 2)+
  geom_segment(data = sites %>% unnest(data) %>% 
                 mutate(lat_i = ifelse(lat_i > lat, lat_i - 0.1, lat_i + 0.1)) %>% 
                 group_by(Site) %>% mutate(mlat = mean(lat_i)),
               aes(x = ifelse(Site == "EKS", -96.1, 
                              ifelse(Site == "JEF", -95.8, long_i)), 
                   y = mlat, 
                   xend = long, yend = lat), 
               arrow = arrow(length = unit(0.25, "cm")))+
  theme_classic(base_size = 12) +
  scale_y_continuous(breaks = seq(37, 40, 1), labels = paste0(seq(37, 40, 1), "\u00B0N"))+
  scale_x_continuous(breaks = seq(-102, -96, 2), labels = paste0(seq(102, 96, -2), "\u00B0W"))+
  coord_fixed(ratio = 1.3)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(color = "black", size = 17),
        axis.text.y = element_text(color = "black", size = 17),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        plot.title = element_text(size = 15),
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines"), 
        legend.text.align = 0, legend.box.just = "center",
        legend.box.margin = margin(-5, -5, -5, -5),
        legend.box.spacing = unit(0.2, "cm"),
        legend.position ="bottom", legend.direction = "horizontal", 
        legend.box = "horizontal",
        legend.key.width = unit(2.5, "cm"), 
        legend.spacing.y = unit(0.3, "cm"), legend.spacing.x = unit(0.05, "cm"),
        legend.title = element_text(size = 19),
        legend.text = element_text(size = 19, margin = margin(t = -0.2, l = -0.5, unit = "lines"))) +
  guides(fill = guide_colorbar(title.position = "top",
                               even.steps = FALSE,
                               show.limits = FALSE, ticks.colour = "black",
                               title.hjust = 0.5, ticks = TRUE,
                               label.theme = element_text(size = 17, lineheight = 0.8),
                               label = TRUE, frame.colour = "black", order = 3),
         color = guide_legend(ncol = 1, order = 2, reverse = FALSE, 
                              label.position = "right", override.aes = list(size = 6)))+
  annotate(geom = "text", x = -102.3, y = 40, label = "(a)", fontface = 2, size = 5)

betas <- read_rds("2.Output/optimsoc_ks.rds") %>% 
  mutate(kconc = log(beta), kste = log(beta + beta_se) - log(beta)) %>% 
  filter(rsq >= 0.8)

model <- lm(kconc ~ ep*Use, data = betas %>% filter(rsq >= 0.8) %>% 
              transform(Use = factor(Use, levels = c("A", "N")))) 
model %>% residuals %>% shapiro.test()
modelres <- model %>% summary()
modelres
modelres$coefficients %>% round(3)

fig6b <- ggplot(betas %>% filter(rsq >= 0.8), aes(x = ep, y = kconc, fill = Use, color = Use))+
  geom_smooth(data = betas %>% filter(rsq >= 0.8),
              method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 6, pch = 21, color = "black", stroke = 0.5)+
  scale_fill_manual(values = c("N" = "darkgreen", "A" = "darkorange"),
                    labels = c("N" = "Native prairie", "A" = "Agriculture"), aesthetics = c("fill", "colour"))+
  theme_bw()+
  scale_y_continuous(expression(paste("ln(", beta, "), cm"^-1*"")))+
  scale_x_continuous("Effective precipitation", limits =  c(min(betas$ep), max(betas$ep)), breaks = seq(-3, 0, 0.5))+
  theme(legend.position = c(0.15, 0.1), 
        legend.title = element_blank(), legend.text = element_text(size = 12),
        text = element_text(size = 14, color = "black"), panel.grid = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        legend.background = element_rect(fill = "transparent"))+
  guides(fill = "none",
         colour = "none")+
  annotate(geom = "text", x = -3, y = 0.5, label = "(b)", fontface = 2, size = 5)+
  annotate(geom = "text", x = -1.4, y = 0.4, size = 4.5, 
           label = expression(paste(italic("p"), ""["ag"]*" = 0.002")))+
  annotate(geom = "text", x = -1.4, y = 0.2, size = 4.5, 
           label = expression(paste(italic("p"), ""["np"]*" < 0.001")))+
  annotate(geom = "text", x = -1.4, y = -0.1, size = 4.5, 
           label = expression(paste(Delta, "ln(", beta, ")"["ag"]*" = -0.937 \u00B1 0.252 cm"^-1*"")))+
  annotate(geom = "text", x = -1.4, y = -0.4, size = 4.5,
           label = expression(paste(Delta, "ln(", beta, ")"["np"]*" = -1.017 \u00B1 0.184 cm"^-1*"")))

fig6c <- ggplot(betas %>% filter(rsq >= 0.8), 
                aes(x = ep, y = SOCz, fill = Use, color = Use))+
  geom_hline(yintercept = 200, color = "gray")+
  geom_smooth(method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 6, pch = 21, color = "black", stroke = 0.5,
             position = position_jitter(width = 0, height = 15))+
  scale_fill_manual(values = c("N" = "darkgreen", "A" = "darkorange"),
                    labels = c("N" = "Native prairie", "A" = "Agriculture"), aesthetics = c("fill", "colour"))+
  theme_bw()+
  scale_y_reverse(expression(paste("z"["SOC"]*", cm")), breaks = seq(0, 500, 100))+
  coord_cartesian(ylim = c(500, -50))+
  scale_x_continuous("Effective precipitation", limits =  c(min(betas$ep), max(betas$ep)), breaks = seq(-3, 0, 0.5))+
  theme(legend.position = c(0.85, 0.93), legend.title = element_blank(), legend.text = element_text(size = 12),
        text = element_text(size = 14, color = "black"), panel.grid = element_blank(),
        axis.text = element_text(size = 13, color = "black"), legend.background = element_rect(fill = "transparent"))+
  guides(fill = "none",
         colour = "none")+
  annotate(geom = "text", x = -3.02, y = -47, label = "(c)", fontface = 2, size = 5)+
  annotate(geom = "text", x = -0.8, y = 0, size = 4.5, 
           label = expression(paste(italic("p"), ""["ag"]*" = 0.014")))+
  annotate(geom = "text", x = -0.8, y = 25, size = 4.5, 
           label = expression(paste(italic("p"), ""["np"]*" = 0.013")))+
  annotate(geom = "text", x = -0.8, y = 50, size = 4.5, 
           label = expression(paste(italic("p"), ""["ag-np"]*" = 0.573")))

fig6 <- fig6a / (fig6b + fig6c) + 
  plot_layout(ncol = 1, nrow = 2, widths = c(0.8, 1), heights = c(2, 2))

ggsave("3.Figures/Figure6.png", fig6, dpi = 500, height = 10, width = 11)
ggsave("3.Figures/pdfs/Figure6.pdf", fig6, dpi = 500, height = 10, width = 11)
rm(list = ls())


#################### Fig. 7&8) Rooting distributions ####################
rootings <- read_csv("1.Data/rootings.csv")

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], 
                                         col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))

## Normalized rooting abundance depth-distribution
normal_rootab <- rootings %>% group_by(Site, Use) %>% 
  mutate(normalizedFine = SumFine/sum(TotalRoot, na.rm = T)*100,
         normalizedCoarse = SumCoarse/sum(TotalRoot, na.rm = T)*100,
         normalizedTotal = TotalRoot/sum(TotalRoot, na.rm = T)*100) %>% 
  melt(id.vars = c("Site", "Use", "Depth", "ppt", "ari", "pet", "ep", "pptpet"),
       measure.vars = c("FracFine", "normalizedFine", "FracCoarse", "normalizedCoarse"),
       variable.name = "roots", value.name = "normalroots") %>% 
  filter(roots != "normalizedTotal") %>% 
  transform(roots = factor(roots, levels = c("FracFine", "normalizedFine", "FracCoarse", "normalizedCoarse"),
                           labels = c("Absolute\nfine root abundance",
                                      "Normalized\nfine root abundance",
                                      "Absolute\ncoarse root abundance",
                                      "Normalized\ncoarse root abundance")),
            Use = factor(Use, levels = c("A", "N"),
                         labels = c("Agriculture", "Native prairie")))

## Plots
for(rt in c("fine", "coarse")){
  if(rt == "fine"){
    lab1 <- expression(paste("Fraction of fine roots cm"^-1*""['soil']*""))
  } else{
    lab1 <- expression(paste("Fraction of coarse roots cm"^-1*""['soil']*""))
  }
  r1 <- ggplot(normal_rootab %>% filter(str_detect(roots, "Absolute"), 
                                        str_detect(roots, rt)), 
               aes(x = normalroots,  y = Depth, color = ep))+
    facet_grid(cols = vars(Use), rows = vars(roots),
               labeller = label_wrap_gen(width = 40), scales = 'free_x', drop = T)+
    geom_point(aes(fill = ep), color = "black", size = 3, pch = 21, stroke = 0.5) +
    scale_x_continuous(lab1, 
                       limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
    scale_y_reverse("Depth, cm")+
    scale_fill_gradientn(colors = pull(filter(custompal, 
                                              round(int, 2) >= round(min(normal_rootab$ep), 2), 
                                              round(int, 2) <= round(max(normal_rootab$ep), 2)), colors),
                         limits = c(min(normal_rootab$ep), max(normal_rootab$ep)),
                         breaks = seq(round(min(normal_rootab$ep)), round(max(normal_rootab$ep)), 0.5),
                         aesthetics = c("fill", "colour"),
                         guide = guide_colourbar(even.steps = FALSE,
                                                 show.limits = FALSE,
                                                 title.hjust = 0.5, ticks = TRUE,
                                                 label.theme = element_text(size = 8,lineheight = 0.8),
                                                 label = TRUE, frame.colour = "black", available_aes = c("fill")))+
    theme_bw()+
    theme(legend.position = "bottom", legend.title = element_text(size = 11),
          panel.grid = element_blank(),
          strip.text = element_text(size = 12),
          text = element_text(size = 12), legend.key.width = unit(2, "cm"))+
    guides(fill = guide_colourbar(title = "Effective precipitation", title.position = "top",
                                  title.hjust = 0.5, frame.colour = "black"), 
           colour = "none")+
    geom_text(data = normal_rootab %>% filter(str_detect(roots, "Absolute"), 
                                              str_detect(roots, rt)) %>% nest(-Use, - roots) %>% arrange(Use) %>% 
                mutate(label = paste0("(", LETTERS[1:2], ")" )),
              aes(label = label, x = 1.09, y = 190), color = "black", size = 4, fontface = 2)
  
  if(rt == "fine"){
    lab2 <- expression(paste("Fine roots cm"['soil']*" ", "Total roots"^-1*""["whole profile"]*", %"))
    lm <- 2
  } else{
    lab2 <- expression(paste("Coarse roots cm"['soil']*" ", "Total roots"^-1*""["whole profile"]*", %"))
    lm <- 1.3
  }
  
  
  r2 <- ggplot(normal_rootab %>% filter(str_detect(roots, "Normalized"), 
                                        str_detect(roots, rt)), 
               aes(x = normalroots,  y = Depth, color = ep))+
    facet_grid(cols = vars(Use), rows = vars(roots),
               labeller = label_wrap_gen(width = 40), scales = 'free_x', drop = T)+
    geom_point(aes(fill = ep), color = "black", size = 3, pch = 21, stroke = 0.5) +
    scale_x_continuous(lab2,
                       limits = c(0, lm), breaks = seq(0, 2, 0.4)) +
    scale_y_reverse("Depth, cm")+
    scale_fill_gradientn(colors = pull(filter(custompal, 
                                              round(int, 2) >= round(min(normal_rootab$ep), 2), 
                                              round(int, 2) <= round(max(normal_rootab$ep), 2)), colors),
                         limits = c(min(normal_rootab$ep), max(normal_rootab$ep)),
                         breaks = seq(round(min(normal_rootab$ep)), round(max(normal_rootab$ep)), 0.5),
                         aesthetics = c("fill", "colour"),
                         guide = guide_colourbar(even.steps = FALSE,
                                                 show.limits = FALSE,
                                                 title.hjust = 0.5, ticks = TRUE,
                                                 label.theme = element_text(size = 8,lineheight = 0.8),
                                                 label = TRUE, frame.colour = "black", available_aes = c("fill")))+
    theme_bw()+
    theme(legend.position = "bottom", legend.title = element_text(size = 11),
          strip.text = element_text(size = 12),
          panel.grid = element_blank(), strip.text.y = element_text(hjust = 0.5, vjust = 0.5),
          text = element_text(size = 12), legend.key.width = unit(2, "cm"))+
    guides(fill = guide_colourbar(title = "Effective precipitation", title.position = "top",
                                  title.hjust = 0.5, frame.colour = "black", barheight = 0.5), 
           colour = "none")+
    geom_text(data = normal_rootab %>%  filter(str_detect(roots, "Normalized"), 
                                               str_detect(roots, rt)) %>% nest(-Use, - roots) %>% arrange(Use) %>% 
                mutate(label = paste0("(", LETTERS[3:4], ")" )),
              aes(label = label, x = lm - 0.1, y = 190), color = "black", size = 4, fontface = 2)
  
  tmp <- (r1 + theme(legend.position = "none"))/(r2 + theme(strip.text.x = element_blank()))
  
  ggsave(paste0("3.Figures/Figure", 6 + which(rt == c("fine", "coarse")), ".png"), tmp, width = 7, height = 8)
  ggsave(paste0("3.Figures/pdfs/Figure", 6 + which(rt == c("fine", "coarse")), ".pdf"), tmp, width = 7, height = 8)
}
rm(list = ls())


#################### Fig. 9) EOC per unit of total root abundance ####################
eoc <- read_csv("1.Data/eoc.csv") %>% filter(!is.na(Horizon))

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))

plota <- ggplot(eoc %>% filter(Use == "A"), aes(x = (eoctotal), y = MidDepth, color = ep, fill = ep))+
  facet_grid(cols = vars(Use), labeller = labeller(Use = c("N" = "Native prairie",
                                                           "A" = "Agriculture")),
             scales = "free")+
  geom_path(aes(group = desc(ep)), lwd = 1.5)+
  geom_point(size = 6, stroke = 0.5, color = "black", pch = 21, show.legend = F)+
  scale_y_reverse("Depth, cm", limits = c(200, -3))+
  scale_x_continuous(expression("Extractable organic carbon (EOC) per unit of total root abundance"),
                     breaks = c(0, 0.5, 1, 2),
                     minor_breaks =c(0.025, 0.050, 0.075, 0.1), guide = "axis_minor")+
  coord_cartesian(xlim = c(0, NA))+
  scale_fill_gradientn(colors = pull(filter(custompal, 
                                            round(int, 2) >= round(min(eoc$ep), 2), 
                                            round(int, 2) <= round(max(eoc$ep), 2)), colors),
                       limits = c(min(eoc$ep), max(eoc$ep)),
                       breaks = seq(round(min(eoc$ep)), round(max(eoc$ep)), 0.5),
                       aesthetics = c("fill", "colour"),
                       guide = guide_colourbar(even.steps = FALSE,
                                               show.limits = FALSE,
                                               title.hjust = 0.5, ticks = TRUE,
                                               label.theme = element_text(size = 8,lineheight = 0.8),
                                               label = TRUE, frame.colour = "black", available_aes = c("fill")))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_text(size = 11),
        axis.text = element_text(color = "black", size = 12),
        text = element_text(size = 12), legend.key.width = unit(2, "cm"), 
        ggh4x.axis.ticks.length.minor = rel(-1),
        legend.key.height = unit(0.4, "cm"))+
  guides(colour = guide_colourbar(title = "Effective precipitation", title.position = "top",
                                  title.hjust = 0.5, frame.colour = "black"))+ 
  annotate(geom = "text", x = 0, y = -2, label = "(a)", fontface = 2, size = 5, vjust = 0, hjust = 0.2)+
  annotation_custom(ggplotGrob(ggplot(eoc %>% filter(Use == "A"), aes(x = (eoctotal), y = MidDepth, color = ep, fill = ep))+
                                 geom_path(aes(group = desc(ep)), lwd = 1.5)+
                                 geom_point(size = 6, stroke = 0.5, color = "black", pch = 21, show.legend = F)+
                                 scale_y_reverse("Depth, cm")+
                                 scale_x_continuous(expression("Extractable organic carbon (EOC) per unit of total root abundance"),
                                                    minor_breaks = c(0.025, 0.050, 0.075),
                                                    breaks = seq(0, 0.05, 0.01),
                                                    guide = "axis_minor")+
                                 coord_cartesian(xlim = c(0, 0.055))+
                                 scale_fill_gradientn(colors = pull(filter(custompal, 
                                                                           round(int, 2) >= round(min(eoc$ep), 2), 
                                                                           round(int, 2) <= round(max(eoc$ep), 2)), colors),
                                                      limits = c(min(eoc$ep), max(eoc$ep)),
                                                      breaks = seq(round(min(eoc$ep)), round(max(eoc$ep)), 0.5),
                                                      aesthetics = c("fill", "colour"),
                                                      guide = guide_colourbar(even.steps = FALSE,
                                                                              show.limits = FALSE,
                                                                              title.hjust = 0.5, ticks = TRUE,
                                                                              label.theme = element_text(size = 8,lineheight = 0.8),
                                                                              label = TRUE, frame.colour = "black", available_aes = c("fill")))+
                                 theme_bw()+
                                 theme(legend.position = "bottom", legend.title = element_text(size = 11),
                                       axis.text = element_text(color = "black", size = 10),
                                       text = element_text(size = 10), legend.key.width = unit(2, "cm"), 
                                       legend.key.height = unit(0.4, "cm"),
                                       axis.title = element_blank(),
                                       ggh4x.axis.ticks.length.minor = rel(-1))+
                                 guides(colour = "none")),
                    xmin = 1, xmax = 2.7, ymin = 0, ymax = -150)

plotb <- ggplot(eoc %>% filter(Use == "N"), aes(x = (eoctotal), y = MidDepth, color = ep, fill = ep))+
  facet_grid(cols = vars(Use), labeller = labeller(Use = c("N" = "Native prairie",
                                                           "A" = "Agriculture")),
             scales = "free")+
  geom_path(aes(group = desc(ep)), lwd = 1.5)+
  geom_point(size = 6, stroke = 0.5, color = "black", pch = 21, show.legend = F)+
  scale_y_reverse("Depth, cm", position = "right", limits = c(200, -3))+
  scale_x_continuous(expression("Extractable organic carbon (EOC) per unit of total root abundance"))+
  coord_cartesian(xlim = c(0, NA))+
  scale_fill_gradientn(colors = pull(filter(custompal, 
                                            round(int, 2) >= round(min(eoc$ep), 2), 
                                            round(int, 2) <= round(max(eoc$ep), 2)), colors),
                       limits = c(min(eoc$ep), max(eoc$ep)),
                       breaks = seq(round(min(eoc$ep)), round(max(eoc$ep)), 0.5),
                       aesthetics = c("fill", "colour"),
                       guide = guide_colourbar(even.steps = FALSE,
                                               show.limits = FALSE,
                                               title.hjust = 0.5, ticks = TRUE,
                                               label.theme = element_text(size = 8,lineheight = 0.8),
                                               label = TRUE, frame.colour = "black", available_aes = c("fill")))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_text(size = 11),
        axis.text = element_text(color = "black", size = 12),
        text = element_text(size = 12), legend.key.width = unit(2, "cm"), 
        legend.key.height = unit(0.4, "cm"))+
  guides(colour = guide_colourbar(title = "Effective precipitation", title.position = "top",
                                  title.hjust = 0.5, frame.colour = "black"))+ 
  annotate(geom = "text", x = 0.0, y = -2, label = "(b)", fontface = 2, size = 5, vjust = 0, hjust = 0.2)
  # annotation_custom(ggplotGrob(ggplot(eoc %>% filter(Use == "N"), aes(x = (eoctotal), y = MidDepth, color = ep, fill = ep))+
  #                                geom_path(aes(group = desc(ep)), lwd = 1.5)+
  #                                geom_point(size = 6, stroke = 0.5, color = "black", pch = 21, show.legend = F)+
  #                                scale_y_reverse("Depth, cm")+
  #                                scale_x_continuous(expression("Extractable organic carbon (EOC) per unit of total root abundance"))+
  #                                coord_cartesian(xlim = c(0, 0.05))+
  #                                scale_fill_gradientn(colors = pull(filter(custompal, 
  #                                                                          round(int, 2) >= round(min(eoc$ep), 2), 
  #                                                                          round(int, 2) <= round(max(eoc$ep), 2)), colors),
  #                                                     limits = c(min(eoc$ep), max(eoc$ep)),
  #                                                     breaks = seq(round(min(eoc$ep)), round(max(eoc$ep)), 0.5),
  #                                                     aesthetics = c("fill", "colour"),
  #                                                     guide = guide_colourbar(even.steps = FALSE,
  #                                                                             show.limits = FALSE,
  #                                                                             title.hjust = 0.5, ticks = TRUE,
  #                                                                             label.theme = element_text(size = 8,lineheight = 0.8),
  #                                                                             label = TRUE, frame.colour = "black", available_aes = c("fill")))+
  #                                theme_bw()+
  #                                theme(legend.position = "bottom", legend.title = element_text(size = 11),
  #                                      axis.text = element_text(color = "black", size = 10),
  #                                      text = element_text(size = 10), legend.key.width = unit(2, "cm"), 
  #                                      legend.key.height = unit(0.4, "cm"),
  #                                      axis.title = element_blank())+
  #                                guides(colour = "none")),
  #                   xmin = 0.050, xmax = 0.09, ymin = -50, ymax = -200)

xlb <- ggplot()+
  annotate(geom = "text", x = 1, y = 1, size = 5,
           label = "Extractable organic carbon (EOC) per unit of total root abundance")+
  coord_cartesian(clip = "off")+
  theme_void()

(plota + theme(axis.title.x = element_blank()) | plotb + theme(axis.title.x = element_blank()))/xlb +
  plot_layout(guides = "collect", heights = c(1, .01)) &
  theme(legend.position = "bottom", strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        # axis.title = element_text(size = 13), axis.text = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))

ggsave("3.Figures/Figure9.png", plot = last_plot(), width = 12, height = 7)
ggsave("3.Figures/pdfs/Figure9.pdf", plot = last_plot(), width = 15, height = 8)
rm(list = ls())


#################### Fig. 10) Depth-distribution of modern C ####################
radiocarbon <- read_csv("1.Data/radiocarbon.csv")

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))

ggplot(radiocarbon, aes(x = pMC/100,  y = Depth_med, color = ep, shape = Use))+
  scale_x_continuous("Fraction modern C", limits = c(0, 1))+
  geom_path(aes(group = ep))+
  geom_point(aes(fill = ep), color = "black", size = 6, stroke = 0.5, guide = NULL)+
  scale_shape_manual(values = c("N" = 21, "A" = 24),
                     labels = c("N" = "Native prairie", "A" = "Agriculture"))+
  scale_y_reverse("Depth, cm")+
  scale_fill_gradientn(colors = pull(filter(custompal, 
                                            round(int, 2) >= round(min(radiocarbon$ep), 2), 
                                            round(int, 2) <= round(max(radiocarbon$ep), 2)), colors),
                       limits = c(min(radiocarbon$ep), max(radiocarbon$ep)),
                       breaks = seq(round(min(radiocarbon$ep)), round(max(radiocarbon$ep)), 0.5),
                       aesthetics = c("fill", "colour"),
                       guide = guide_colourbar(even.steps = FALSE,
                                               show.limits = FALSE,
                                               title.hjust = 0.5, ticks = TRUE,
                                               label.theme = element_text(size = 8,lineheight = 0.8),
                                               label = TRUE, frame.colour = "black", available_aes = c("fill")))+
  theme_bw()+
  theme(legend.position = c(0.75, 0.125), legend.title = element_text(size = 10),
        text = element_text(size = 12), legend.key.width = unit(1.7, "cm"),
        axis.text = element_text(size = 13), legend.direction = "horizontal",
        panel.grid = element_blank(), legend.box.just = "right", legend.spacing.y = unit(0.11, "cm"),
        legend.box = "vertical", legend.spacing.x = unit(0.5, "cm"),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 10, margin = margin(r = -3, unit = "lines")))+
  guides(fill = guide_colourbar(title = "Effective precipitation", title.position = "top",
                                title.hjust = 0.5, frame.colour = "black"),
         shape = guide_legend(title = NULL, order = 1, reverse = T,  label.position = "left"),
         colour = "none")
ggsave("3.Figures/Figure10.png", last_plot(), width = 8, height = 6)
ggsave("3.Figures/pdfs/Figure10.pdf", last_plot(), width = 8, height = 6)
rm(list = ls())


#################### Fig. 11) Depth-distribution of aggregates ####################
clim <- read_csv("1.Data/sampsites_climvars.csv")
avgaggs <- read_csv("1.Data/avgaggs.csv")
avgaggs <- left_join(avgaggs, clim, by = c("Site", "Use"))

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))

## Figure 10a
fig11a <- ggplot(avgaggs, aes(x = avgMAD,  y = MidDepth, color = ep, shape = Use))+
  facet_grid(cols = vars(Use), labeller = labeller(Use = c("A" = "Agriculture",
                                                           "N" = "Native prairie")),
             scales = "free_x")+
  scale_x_continuous("Mean aggregate diameter, mm", breaks = seq(0, 1.6, 0.4),
                     limits = c(0, NA))+
  scale_shape_manual(values = c("N" = 21, "A" = 24),
                     labels = c("N" = "Native prairie", "A" = "Agriculture"))+
  geom_path(aes(group = ep))+
  geom_point(aes(fill = ep), color = "black", size = 5, stroke = 0.5)+
  geom_errorbar(aes(xmin = avgMAD - sdMAD, xmax = avgMAD + sdMAD), color = "black")+
  scale_y_reverse("Depth, cm")+
  scale_fill_gradientn(colors = pull(filter(custompal, 
                                            round(int, 2) >= round(min(avgaggs$ep), 2), 
                                            round(int, 2) <= round(max(avgaggs$ep), 2)), colors),
                       limits = c(min(avgaggs$ep), max(avgaggs$ep)),
                       breaks = seq(round(min(avgaggs$ep)), round(max(avgaggs$ep)), 0.5),
                       aesthetics = c("fill", "colour"),
                       guide = guide_colourbar(even.steps = FALSE,
                                               show.limits = FALSE,
                                               title.hjust = 0.5, ticks = TRUE,
                                               label.theme = element_text(size = 8,lineheight = 0.8),
                                               label = TRUE, frame.colour = "black", available_aes = c("fill")))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_text(size = 10),
        text = element_text(size = 12), legend.key.width = unit(1.7, "cm"),
        axis.text = element_text(size = 13), legend.direction = "horizontal",
        panel.grid = element_blank(), legend.box.just = "center", legend.spacing.y = unit(0.11, "cm"),
        legend.box = "horizontal", legend.spacing.x = unit(-0.5, "cm"),
        axis.line = element_line(color = "black"))+
  guides(fill = guide_colourbar(title = "Effective precipitation", title.position = "top",
                                title.hjust = 0.5, frame.colour = "black"),
         shape = 'none',
         colour = "none")+
  geom_text(data = avgaggs %>% nest(-Use) %>% arrange(Use) %>% 
              mutate(label = paste0("(", letters[1:2], ")" )),
            aes(label = label, x = c(1.4, 1.85), y = c(195, 195)), color = "black", size = 4, fontface = 2)

## Figure 10b
fig11b <- ggplot(avgaggs, aes(x = intmacroOC,  y = MidDepth, color = ep, shape = Use))+
  facet_grid(cols = vars(Use), labeller = labeller(Use = c("A" = "Agriculture",
                                                           "N" = "Native prairie")),
             scales = "free_x")+
  scale_x_continuous(expression("g"["soil in aggregates"]*" g"^-1*""["organic carbon"]*""),
                     breaks = seq(0, 900, 200))+
  scale_shape_manual(values = c("N" = 21, "A" = 24),
                     labels = c("N" = "Native prairie", "A" = "Agriculture"))+
  geom_path(aes(group = ep))+
  geom_point(aes(fill = ep), color = "black", size = 5, stroke = 0.5)+
  geom_errorbar(aes(xmin = intmacroOC - sdintmacroOC, xmax = intmacroOC + sdintmacroOC), color = "black")+
  scale_y_reverse("Depth, cm")+
  scale_fill_gradientn(colors = pull(filter(custompal, 
                                            round(int, 2) >= round(min(avgaggs$ep), 2), 
                                            round(int, 2) <= round(max(avgaggs$ep), 2)), colors),
                       limits = c(min(avgaggs$ep), max(avgaggs$ep)),
                       breaks = seq(round(min(avgaggs$ep)), round(max(avgaggs$ep)), 0.5),
                       aesthetics = c("fill", "colour"),
                       guide = guide_colourbar(even.steps = FALSE,
                                               show.limits = FALSE,
                                               title.hjust = 0.5, ticks = TRUE,
                                               label.theme = element_text(size = 8,lineheight = 0.8),
                                               label = TRUE, frame.colour = "black", available_aes = c("fill")))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_text(size = 10),
        text = element_text(size = 12), legend.key.width = unit(1.7, "cm"),
        axis.text = element_text(size = 13), legend.direction = "horizontal",
        axis.title.x = element_text(size = 14),
        panel.grid = element_blank(), legend.box.just = "center", legend.spacing.y = unit(0.11, "cm"),
        legend.box = "horizontal", legend.spacing.x = unit(-0.5, "cm"),
        axis.line = element_line(color = "black"))+
  guides(fill = guide_colourbar(title = "Effective precipitation", title.position = "top",
                                title.hjust = 0.5, frame.colour = "black"),
         shape = "none", 
         colour = "none")+
  geom_text(data = avgaggs %>% nest(-Use) %>% arrange(Use) %>% 
              mutate(label = paste0("(", letters[3:4], ")" )),
            aes(label = label, x = c(820, 950), y = c(195, 195)), color = "black", size = 4, fontface = 2)

prow <- plot_grid(fig11a + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), legend.position = "none", plot.title = element_blank()),
                  fig11b + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), legend.position = "none", plot.title = element_blank()),
                  align = "h",
                  byrow = T,
                  # hjust = -3, 
                  nrow = 2,
                  ncol = 1)

lg <- get_legend(fig11a)
prow2 <- plot_grid(prow, lg, ncol = 1, rel_heights = c(1, .1))
save_plot("3.Figures/Figure11.png", prow2, dpi = 450, base_width = 7, base_height = 4, nrow = 2, ncol = 1)
save_plot("3.Figures/pdfs/Figure11.pdf", prow2, dpi = 450, base_width = 7, base_height = 4, nrow = 2, ncol = 1)
rm(list = ls())
