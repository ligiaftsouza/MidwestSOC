library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(broom)
library(scales)
library(reghelper)

#################### Supplementary Fig. 1 ####################
profiles <- readRDS("1.Data/soils_midwest.rds")
profiles <- profiles %>% filter(id %in% c("95P0518", "40A2396", "64IL059002")) %>% 
  unnest(data)

socpred <- read_rds( "2.Output/socconc_param.rds")
predprof <- socpred %>% filter(id %in% c("95P0518", "40A2396", "64IL059002"))

for(p in c("95P0518", "40A2396", "64IL059002")){
  mx <- max(profiles$soc[profiles$id == p], na.rm = T)
  ep <- unique(round(profiles$EP[profiles$id == p], 2))
  if(p == "95P0518"){
    bt <-  annotate(geom = "text", x = 1.5, y = 330, label = expression(paste("ln(", beta, ") = ", "-3.24 cm"^-1*"")), hjust = 0)
    isoc <- annotate(geom = "text", x = 1.5, y = 360, label = expression(paste("SOC"["i"]*" = 2.30%")), hjust = 0)
    fsoc <- annotate(geom = "text", x = 1.5, y = 390, label = expression(paste("SOC"["f"]*" = 0.35%")), hjust = 0)
    zsoc <- annotate(geom = "text", x = 1.5, y = 420, label = expression(paste("z"["SOC"]*" = 346 cm")), hjust = 0)
  } else if(p == "40A2396"){
    bt <-  annotate(geom = "text", x = 1.5, y = 330, label = expression(paste("ln(", beta, ") = ", "-3.39 cm"^-1*"")), hjust = 0)
    isoc <- annotate(geom = "text", x = 1.5, y = 360, label = expression(paste("SOC"["i"]*" = 3.48%")), hjust = 0)
    fsoc <- annotate(geom = "text", x = 1.5, y = 390, label = expression(paste("SOC"["f"]*" = 0.02%")), hjust = 0)
    zsoc <- annotate(geom = "text", x = 1.5, y = 420, label = expression(paste("z"["SOC"]*" = 412 cm")), hjust = 0)
  } else{
    bt <- annotate(geom = "text", x = 1.5, y = 330, label = expression(paste("ln(", beta, ") = ", "-3.90 cm"^-1*"")), hjust = 0)
    isoc <- annotate(geom = "text", x = 1.5, y = 360, label = expression(paste("SOC"["i"]*" = 1.40%")), hjust = 0)
    fsoc <- annotate(geom = "text", x = 1.5, y = 390, label = expression(paste("SOC"["f"]*" = 0.09%")), hjust = 0)
    zsoc <- annotate(geom = "text", x = 1.5, y = 420, label = expression(paste("z"["SOC"]*" = 500 cm")), hjust = 0)
  }
  pp <- ggplot(profiles %>% filter(id == p), 
               aes(x = soc, y = depth))+
    facet_grid(cols = vars(id), scales = "free_x")+
    geom_path(data = predprof %>% filter(id == p) %>% unnest(data), 
              aes(x = pred_soc, y = pred_depth), color = "cornsilk4")+
    geom_point(size = 5, fill = "cornsilk4", color = "black", pch = 21)+
    scale_y_reverse("Depth, cm", breaks = seq(0, 500, 50), limits = c(500, 0))+
    scale_x_continuous("Soil organic carbon, %",  breaks = seq(0, 3, 1))+
    coord_cartesian(xlim = c(0, 3.5))+
    annotate(geom = "text", x = 1.5, y = 300, label = paste("Effective precipitation = ", ep), hjust = 0)+
    bt +
    isoc +
    fsoc +
    zsoc +
    ggtitle(paste0("Soil profile ID: ", p))+
    theme_bw() +
    theme(strip.text.x = element_blank(), panel.grid = element_blank(), 
          legend.position = "none", axis.text = element_text(colour = "black", size = 11),
          axis.title = element_text(size = 12))
  assign(paste0("p", which(c("95P0518", "40A2396", "64IL059002") == p)), pp)
}

soc <- read_rds("2.Output/optimsoc_ks.rds") %>% 
  filter(rsq >= 0.8, Use == "N") %>% unnest(data) %>% unnest(data) %>% 
  nest(-Site, -Use, -soc, -depth)

kssocpred <- read_rds("2.Output/optimsoc_ks.rds")
kssocpred <- kssocpred %>% filter(Use == "N")

for(ks in c("LGN", "RKS", "EKS")){
  mx <- max(soc$soc[soc$Site == ks], na.rm = T)
  ep <- soc %>% filter(Site == ks) %>% unnest(data) %>% pull(ep) %>% round(2) %>% unique 
  if(ks == "LGN"){
    bt <-  annotate(geom = "text", x = 2, y = 330, label = expression(paste("ln(", beta, ") = ", "-0.29 cm"^-1*"")), hjust = 0)
    isoc <- annotate(geom = "text", x = 2, y = 360, label = expression(paste("SOC"["i"]*" = 2.24%")), hjust = 0)
    fsoc <- annotate(geom = "text", x = 2, y = 390, label = expression(paste("SOC"["f"]*" = 1.68%")), hjust = 0)
    zsoc <- annotate(geom = "text", x = 2, y = 420, label = expression(paste("z"["SOC"]*" = 16 cm")), hjust = 0)
  } else if(ks == "RKS"){
    bt <-  annotate(geom = "text", x = 2, y = 330, label = expression(paste("ln(", beta, ") = ", "-1.89 cm"^-1*"")), hjust = 0)
    isoc <- annotate(geom = "text", x = 2, y = 360, label = expression(paste("SOC"["i"]*" = 3.09%")), hjust = 0)
    fsoc <- annotate(geom = "text", x = 2, y = 390, label = expression(paste("SOC"["f"]*" = 0.75%")), hjust = 0)
    zsoc <- annotate(geom = "text", x = 2, y = 420, label = expression(paste("z"["SOC"]*" = 83 cm")), hjust = 0)
  } else{
    bt <- annotate(geom = "text", x = 2, y = 330, label = expression(paste("ln(", beta, ") = ", "-2.99 cm"^-1*"")), hjust = 0)
    isoc <- annotate(geom = "text", x = 2, y = 360, label = expression(paste("SOC"["i"]*" = 3.00%")), hjust = 0)
    fsoc <- annotate(geom = "text", x = 2, y = 390, label = expression(paste("SOC"["f"]*" = 0.19%")), hjust = 0)
    zsoc <- annotate(geom = "text", x = 2, y = 420, label = expression(paste("z"["SOC"]*" = 269 cm")), hjust = 0)
  }
  pp <- ggplot(soc %>% filter(Site == ks), 
               aes(x = soc, y = depth))+
    facet_grid(cols = vars(Site), scales = "free_x")+
    geom_path(data = kssocpred %>% filter(Site == ks) %>% unnest(data), 
              aes(x = soc_pred, y = pred_depth), color = "cornsilk4")+
    geom_point(size = 5, fill = "cornsilk4", color = "black", pch = 21)+
    scale_y_reverse("Depth, cm", breaks = seq(0, 500, 50), limits = c(510, 0))+
    scale_x_continuous("Soil organic carbon, %", limits = c(0, 4))+
    coord_cartesian(xlim = c(0, NA))+
    annotate(geom = "text", x = 2, y = 300, 
             label = paste("Effective precipitation = ", ep), hjust = 0)+
    bt +
    isoc +
    fsoc +
    zsoc +
    ggtitle(paste0("Soil profile in ", c("Logan County, KS", "Rooks County, KS", "Eastern Kansas, KS")[which(c("LGN", "RKS", "EKS") == ks)]))+
    theme_bw() +
    theme(strip.text.x = element_blank(), panel.grid = element_blank(), 
          legend.position = "none", axis.text = element_text(colour = "black", size = 11),
          axis.title = element_text(size = 12))
  assign(paste0("p", 3 + which(c("LGN", "RKS", "EKS") == ks)), pp)
}

plot <- plot_grid(p1 + theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm")),
          p2 + theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm")),
          p3 + theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm")),
          p4 + theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm")),
          p5 + theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm")),
          p6 + theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm")),
          align = 'h',
          labels = c("(a)", "(d)", "(b)", "(e)", "(c)", "(f)"),
          byrow = F,
          hjust = -1, 
          nrow = 3,
          ncol = 2)
save_plot("3.Figures/FigureS1.png", plot, dpi = 450, 
          base_width = 6, base_height = 4.5, nrow = 3, ncol = 2)
save_plot("3.Figures/pdfs/FigureS1.pdf", plot, dpi = 450, 
          base_width = 6, base_height = 4.5, nrow = 3, ncol = 2)
rm(list = ls())


#################### Supplementary Fig. 2 ####################
clim <- read_csv("1.Data/sampsites_climvars.csv")
avgaggs <- read_csv("1.Data/avgaggs.csv")
avgaggs <- left_join(avgaggs, clim, by = c("Site", "Use"))

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))


maggs <- left_join(avgaggs %>% 
                     melt(id.vars = c("Site", "Use", "Horizon", "MidDepth", "ep"),
                          measure.vars = c("micro", "intmacro", "macro"),
                          variable.name = "agg", value.name = "avgagg"), 
                   avgaggs %>% 
                     melt(id.vars = c("Site", "Use", "Horizon", "MidDepth", "ep"),
                          measure.vars = c("sdmicro", "sdintmacro", "sdmacro"),
                          variable.name = "agg", value.name = "sdagg") %>% 
                     mutate(agg = str_extract(agg, "(?<=sd)\\w{5,8}")),
                   id = c("Site", "Use", "Horizon", "MidDepth", "ep", "agg")) %>% 
  transform(agg = factor(agg, levels = c("micro", "intmacro", "macro"), 
                         labels = c("Microaggregates (< 0.21 mm)",
                                    "Intermediate\nmacroaggregates (0.21 - 4.76 mm)",
                                    "Macroaggregates (> 4.76 mm)")),
            Use = factor(Use, levels = c("A", "N"),
                         labels = c("Agriculture", 
                                    expression(paste("Native prairie")))))

ggplot(maggs, 
       aes(x = avgagg,  y = MidDepth, color = ep, shape = Use))+
  facet_grid(cols = vars(Use), rows = vars(agg),
             labeller = labeller(Use = label_parsed,
                                 agg = label_wrap_gen(width = 18)),
             scales = "free_x")+
  scale_x_continuous(expression(paste("Aggregate size distribution, g"["aggregate"]*" g"^-1*""["soil"]*"")), 
                     breaks = seq(0, 1, 0.2),
                     limits = c(-0.1, 1))+
  # scale_shape_manual(values = c("N" = 21, "A" = 24),
  #                    labels = c("N" = "Native prairie", "A" = "Agriculture"))+
  scale_shape_manual(values = c(24, 21))+
  geom_path(aes(group = ep))+
  geom_point(aes(fill = ep), color = "black", size = 5, stroke = 0.5)+
  geom_errorbar(aes(xmin = avgagg - sdagg, xmax = avgagg + sdagg))+
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
        legend.key.height = unit(0.5, "cm"),
        axis.text = element_text(size = 13), legend.direction = "horizontal",
        panel.grid = element_blank(), legend.box.just = "center", 
        legend.spacing.y = unit(0.11, "cm"),
        legend.box = "vertical", legend.spacing.x = unit(-0.5, "cm"),
        axis.line = element_line(color = "black"),
        strip.text.y = element_text(hjust = 0.5, vjust = 0.5))+
  guides(fill = guide_colourbar(title = "Effective precipitation", title.position = "top",
                                title.hjust = 0.5, frame.colour = "black"),
         shape = "none",
         colour = "none")+
  geom_text(data = maggs %>% nest(-Use, -agg) %>% arrange(Use) %>% 
              mutate(label = paste0("(", letters[1:6], ")" )),
            aes(label = label, x = -0.08, y = 5), color = "black", size = 4, fontface = 2)

ggsave("3.Figures/FigureS2.png", last_plot(), width = 7, height = 8, dpi = 400)
ggsave("3.Figures/pdfs/FigureS2.pdf", last_plot(), width = 7, height = 8, dpi = 400)
rm(list = ls())

#################### Supplementary Fig. 3 ####################
eoc <- read_csv("1.Data/eoc.csv")

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))

ggplot(eoc, aes(x = eoc, y = MidDepth, color = ep, fill = ep))+
  facet_grid(cols = vars(Use), labeller = labeller(Use = c("N" = "Native prairie",
                                                           "A" = "Agriculture")))+
  geom_path(aes(group = desc(ep)), lwd = 1.2)+
  geom_point(size = 4, stroke = 0.5, color = "black", pch = 21, show.legend = F)+
  scale_y_reverse("Depth, cm")+
  scale_x_continuous(expression("Extractable organic carbon (EOC), mg g"^-1*""["dry soil"]), 
                     breaks = c(0.0, 0.05, 0.10, 0.15), limits = c(0, 0.15))+
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
        panel.grid = element_blank(), axis.text = element_text(color = "black", size = 10),
        text = element_text(size = 11), legend.key.width = unit(2, "cm"), 
        legend.key.height = unit(0.4, "cm"))+
  guides(colour = guide_colourbar(title = "Effective precipitation", title.position = "top",
                                  title.hjust = 0.5, frame.colour = "black"))
ggsave("3.Figures/FigureS3.png", plot = last_plot(), width = 7, height = 5)
ggsave("3.Figures/pdfs/FigureS3.pdf", plot = last_plot(), width = 7, height = 5)
rm(list = ls())



#################### Supplementary Fig. 4 ####################
clim <- read_csv("1.Data/sampsites_climvars.csv")
nrcs <- read_csv("1.Data/nrcs.csv")
nrcs <- left_join(nrcs, clim, by = c("Site", "Use"))

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))

ggplot(nrcs, aes(x = COLE, y = MidDepth, color = ep, fill = ep))+
  facet_grid(cols = vars(Use), labeller = labeller(Use = c("N" = "Native prairie",
                                                           "A" = "Agriculture")))+
  geom_path(aes(group = desc(ep)), lwd = 1.5)+
  geom_point(size = 6, stroke = 0.5, color = "black", pch = 21, show.legend = F)+
  scale_y_reverse("Depth, cm", limits = c(210, 0))+
  scale_x_continuous("Coefficient of linear extensibility (COLE)",
                     breaks = c(0, 0.03, 0.06, 0.09, 0.12), limits = c(0, 0.13))+
  scale_fill_gradientn(colors = pull(filter(custompal, 
                                            round(int, 2) >= round(min(nrcs$ep), 2), 
                                            round(int, 2) <= round(max(nrcs$ep), 2)), colors),
                       limits = c(min(nrcs$ep), max(nrcs$ep)),
                       breaks = seq(round(min(nrcs$ep)), round(max(nrcs$ep)), 0.5),
                       aesthetics = c("fill", "colour"),
                       guide = guide_colourbar(even.steps = FALSE,
                                               show.limits = FALSE,
                                               title.hjust = 0.5, ticks = TRUE,
                                               label.theme = element_text(size = 8,lineheight = 0.8),
                                               label = TRUE, frame.colour = "black", available_aes = c("fill")))+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_text(size = 11),
        legend.key.height = unit(0.4, "cm"),
        panel.grid = element_blank(), axis.text = element_text(color = "black", size = 12),
        text = element_text(size = 13), legend.key.width = unit(2, "cm"))+
  geom_segment(aes(x = 0, xend = 0.03, y = 205, yend = 205), color = "black",
               arrow = arrow(ends = "both", angle = 90, length = unit(0.1, "cm")))+
  annotate(geom = "text", x = 0.015, y = 208, label = "low")+
  geom_segment(aes(x = 0.031, xend = 0.06, y = 205, yend = 205), color = "black",
               arrow = arrow(ends = "both",  angle = 90, length = unit(0.1, "cm")))+
  annotate(geom = "text", x = 0.045, y = 208, label = "moderate", hjust = 0.5)+
  geom_segment(aes(x = 0.061, xend = 0.09, y = 205, yend = 205), color = "black",
               arrow = arrow(ends = "both",  angle = 90, length = unit(0.1, "cm")))+
  annotate(geom = "text", x = 0.075, y = 208, label = "high", hjust = 0.5)+
  geom_segment(aes(x = 0.091, xend = 0.115, y = 205, yend = 205),color = "black",
               arrow = arrow(ends = "both",  angle = 90, length = unit(0.1, "cm")))+
  annotate(geom = "text", x = 0.105, y = 208, label = "very high",  hjust = 0.6)+
  guides(colour = guide_colourbar(title = "Effective precipitation", title.position = "top",
                                  title.hjust = 0.5, frame.colour = "black"))
ggsave("3.Figures/FigureS4.png", plot = last_plot(), width = 9, height = 6)
ggsave("3.Figures/pdfs/FigureS4.pdf", plot = last_plot(), width = 9, height = 6)
rm(list = ls())


#################### Supplementary Fig. 5 ####################
nrcs <- read_csv("1.Data/nrcs.csv")
clim <- read_csv("1.Data/sampsites_climvars.csv")
nrcs <- left_join(nrcs, clim, by = c("Site", "Use")) %>% filter(!str_detect(soc, "tr")) %>% 
  transform(soc = as.numeric(soc))

custompal <- data.frame(cbind(int = as.numeric(c(-4, -3.83, -3.07, -2.5, -1.87, -1.5, 
                                                 -0.754, -0.458, -0.275, 0.5)), 
                              colors = c(brewer.pal(11, "RdYlBu")[c(1:7, 10:11)], col2hcl("navy")))) %>% 
  transform(int = as.numeric(int))


coc <- nrcs %>% dplyr::select(c(-COLE, -WaterCont, -totalporosity, -effporosity)) %>% 
  mutate(coc = as.numeric(ifelse((soc < clay*0.163) == TRUE, soc, clay*0.163)),
         ncoc = as.numeric(ifelse((coc - soc) == TRUE, coc - soc, 0)),
         cc = as.numeric(ifelse((soc/0.163 < clay) == TRUE, soc/0.163, clay)),
         ncc = as.numeric(ifelse((clay - cc > 0) == TRUE, clay - cc, 0)),
         cmax = as.numeric(clay*0.163),
         paoc = ifelse((cmax-coc > 0) == TRUE, (cmax-coc), 0),
         ### COC indices above are in %
         ### COC indices in areal basis, multiplying by BD, layer height, and 10 to convert from g/cm2 to kg/m2
         coc_area= (coc/100)*BD*(horizon_bottom - horizon_top)*10,
         cmax_area = (cmax/100)*BD*(horizon_bottom - horizon_top)*10,
         paoc_area = (paoc/100)*BD*(horizon_bottom - horizon_top)*10,
         pits = paste0(Site, Use))

ggplot(coc, aes(x = (cc + coc)/(soc + clay)*100,  y = MidDepth, color = ep, shape = Use, fill = ep))+
  scale_shape_manual(values = c("N" = 21, "A" = 24),
                     labels = c("N" = "Native prairie", "A" = "Agriculture"))+
  geom_path(aes(group = ep))+
  geom_point(aes(fill = ep), color = "black", size = 6, stroke = 0.5) +
  scale_x_continuous(expression("Complexed clay and organic carbon, %")) +
  scale_y_reverse("Depth, cm")+
  scale_fill_gradientn(colors = pull(filter(custompal, 
                                            round(int, 2) >= round(min(coc$ep), 2), 
                                            round(int, 2) <= round(max(coc$ep), 2)), colors),
                       limits = c(min(coc$ep), max(coc$ep)),
                       breaks = seq(round(min(coc$ep)), round(max(coc$ep)), 0.5),
                       aesthetics = c("fill", "colour"),
                       guide = guide_colourbar(even.steps = FALSE,
                                               show.limits = FALSE,
                                               title.hjust = 0.5, ticks = TRUE,
                                               label.theme = element_text(size = 8,lineheight = 0.8),
                                               label = TRUE, frame.colour = "black", available_aes = c("fill")))+
  theme_bw()+
  theme(legend.position = c(0.7, 0.15), legend.title = element_text(size = 10),
        text = element_text(size = 12), legend.key.width = unit(1.7, "cm"),
        legend.key.height = unit(0.5, "cm"),
        axis.text = element_text(size = 13, color = "black"), legend.direction = "horizontal",
        panel.grid = element_blank(), legend.box.just = "right", legend.spacing.y = unit(0.11, "cm"),
        legend.box = "vertical", legend.spacing.x = unit(0.5, "cm"),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 10, margin = margin(r = -3, unit = "lines")))+
  guides(fill = guide_colourbar(title = "Effective precipitation", title.position = "top",
                                title.hjust = 0.5, frame.colour = "black"),
         shape = guide_legend(title = NULL, order = 1, reverse = T,  label.position = "left"),
         colour = "none")
ggsave("3.Figures/FigureS5.png", last_plot(), width = 7, height = 6)
ggsave("3.Figures/pdfs/FigureS5.pdf", last_plot(), width = 7, height = 6)
rm(list = ls())
