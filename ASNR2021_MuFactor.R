#' ---
#' title: "R code for the analysis conducted in 'Geodesic path analysis of neural networks in the Alzheimer\\'s disease connectome project'"
#' author: "Nagesh Adluru"
#' output:
#'    html_document:
#'      toc: true
#'      toc_float: true
#' ---
#'
#' # Initialization
# Loading the libraries =========
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(purrr)
library(forcats)
library(stringr)
library(broom)
library(lemon)
library(latex2exp)
library(tidytext)
library(circlize)
library(grDevices)
library(scales)
library(igraph)
library(ComplexHeatmap)
library(progress)

# Initializing variables ======
rm(list = ls(all = T))
csvroot = 'CSVs/'
figroot = 'Figures/'

# IIT Desikan node info=======
nodes = read.csv(paste0(csvroot, 'LUT_GM_Desikan_COG_Mass.csv'), as.is = T, header = T)
hemispheres = c('Left', 'Right')
lobes = nodes$Lobe %>% unique
LobeLobeList = map(6:1, ~map(.x:1, ~paste0(lobes[[.y]], "_", lobes[[.x]]), .y = .x)) %>% unlist
LobeLobeMirroredList = map(6:1, ~map(.x:1, ~paste0(lobes[[.x]], "_", lobes[[.y]]), .y = .x)) %>% unlist

# ggplot theme =======
dodge = position_dodge(width = 0.9)
txtSize = 12
gtheme = theme(
  legend.text = element_text(size = txtSize),
  legend.title = element_blank(),
  legend.key = element_blank(),
  legend.background = element_blank(),
  legend.position = "top", 
  strip.text.x = element_text(size = txtSize),
  strip.text.y = element_text(size = txtSize),
  strip.background = element_blank(),
  axis.text = element_text(colour = "black", size = txtSize),
  axis.title.x = element_text(size = txtSize),
  axis.title.y = element_text(size = txtSize),
  plot.title = element_text(size = txtSize),
  axis.title = element_text(size = txtSize),
  axis.line = element_line(),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "gray"),
  panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "gray"),
  # ticks
  axis.ticks.length = unit(0.25, "cm")
)

#' # Basic demographics
basicdemo = read.csv(paste0(csvroot, 'BasicDemo_ASNR2021.csv')) %>% mutate(ConsensusDiagnosis = fct_relevel(ConsensusDiagnosis, 'AD', after = Inf)) %>% mutate(ConsensusDiagnosis = fct_relevel(ConsensusDiagnosis, 'AD', after = Inf)) %>% 
  filter(!(SubjectID %in% c('ADCP_2066'))) # brain masking left eye tissue
basicdemosummary = basicdemo %>% group_by(Sex, ConsensusDiagnosis) %>% 
  summarise(ma = mean(visit_age), 
            n = n(), 
            sa = sd(visit_age),
            .groups = 'drop')

#+ fig.width=5.65, fig.height=6.20, warning=F
p = basicdemo %>% 
  filter(visit_age < 82) %>% # anonymity filter 
  ggplot(aes(x = visit_age, y = fct_reorder(SubjectID, visit_age))) + 
  geom_line(aes(group = SubjectID), alpha = 0.2) + 
  geom_point(size = 1, shape = 21, fill = NA) + 
  geom_vline(data = basicdemosummary, aes(xintercept = ma), linetype = 'longdash') +
  geom_vline(data = basicdemosummary, 
             aes(xintercept = ma - sa), 
             alpha = 0.75) + 
  geom_vline(data = basicdemosummary, 
             aes(xintercept = ma + sa), 
             alpha = 0.75) + 
  geom_rect(data = basicdemosummary, 
            aes(xmin = ma - sa, 
                xmax = ma + sa, 
                ymin = -Inf, 
                ymax = Inf), 
            alpha = 0.1, 
            fill = 'gray64', inherit.aes = F) +
  geom_text(data = basicdemosummary, 
            aes(x = Inf, y = -Inf, 
                label = paste0('n = ', n)), 
            size = txtSize/1.5, hjust = 1, vjust = -0.5) +
  facet_rep_grid(ConsensusDiagnosis ~ Sex, scales = 'free_y') + 
  labs(x = 'Age y', y = 'Subjects ordered by their age') + 
  gtheme + 
  scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.x = element_line(size = 0.2, 
                                          linetype = 'solid', 
                                          color = 'gray'), 
        panel.spacing.y = unit(-1.0, 'lines'),
        axis.ticks.length.y = unit(0.125, 'cm'))
p
ggsave(paste0(figroot, 'BasicDemographics_ASNR2021', '.pdf'),
       width = 5.65,
       height = 6.20,
       p)

#' # Geodesic distance matrix
pb = progress_bar$new(total = basicdemo %>% nrow)
distancedf = pmap(list(basicdemo$ConnectomePaths,
                       basicdemo$Subject.ID.Number), ~{
                         pb$tick()
                         wtc = read.csv(..1,
                                        sep = ",",
                                        header = F) %>% as.matrix
                         wtc = (wtc + t(wtc)) * (paste0(csvroot, 'ADCP_', ..2, '_mu.txt') %>% scan)
                         net = graph_from_adjacency_matrix(wtc^(-1),
                                                           weighted = T,
                                                           diag = F,
                                                           mode = "undirected")
                         distancedf = graph_from_adjacency_matrix(
                           distances(net) %>% as.matrix,
                           weighted = T,
                           diag = F,
                           mode = "undirected"
                         ) %>%
                           igraph::as_data_frame("edges") %>%
                           left_join(nodes %>% select(Vertex, Name, Hemisphere, Lobe),
                                     by = c("from" = "Vertex")) %>%
                           left_join(nodes %>% select(Vertex, Name, Hemisphere, Lobe),
                                     by = c("to" = "Vertex")) %>%
                           mutate(
                             HemHem = plyr::mapvalues(
                               paste0(Hemisphere.x, "_", Hemisphere.y),
                               from = "Right_Left",
                               to = "Left_Right"
                             ),
                             LobeLobe = plyr::mapvalues(paste0(Lobe.x, "_", Lobe.y),
                                                        from = LobeLobeList,
                                                        to = LobeLobeMirroredList)
                           ) %>%
                           pivot_longer(c(HemHem, LobeLobe),
                                        names_to = 'RegionalClass',
                                        values_to = 'RegionalName') %>%
                           group_by(RegionalName) %>%
                           summarise(MeanDistance = mean(weight, na.rm = T)) %>%
                           mutate(ID = ..2)
                       }) %>%
  bind_rows %>%
  left_join(basicdemo, by = c('ID' = 'Subject.ID.Number')) %>%
  rename(Age = visit_age,
         Consensus.Diagnosis = ConsensusDiagnosis)

# Removing cerebellar regions and normalizing the mean distance by the mean distance of the CU group
distancedf %>%
  left_join(distancedf %>%
              filter(Consensus.Diagnosis %>% str_detect('CU') &
                       is.finite(MeanDistance)) %>%
              group_by(RegionalName) %>%
              summarise(MD = mean(MeanDistance, na.rm = T))) %>%
  mutate(MeanDistance = MeanDistance / MD) %>%
  filter(!str_detect(RegionalName, "Cerebell.*") &
           !is.infinite(MeanDistance) &
           !(RegionalName %in% map_chr(1:6, ~paste0(lobes[[.x]], "_", lobes[[.x]]))) &
           !(RegionalName %in% map2_chr(c(1, 1, 2), c(1, 2, 2), ~paste0(hemispheres[[.x]], "_", hemispheres[[.y]])))) %>%
  write.csv(paste0(csvroot, 'DistancesRegionalSelect_ASNR2021.csv'))

distancedf = read.csv(paste0(csvroot, 'DistancesRegionalSelect_ASNR2021.csv')) %>% mutate(Consensus.Diagnosis = Consensus.Diagnosis %>% as.factor %>% fct_relevel('AD', after = Inf)) %>% 
  filter(!(SubjectID %in% c('ADCP_2066'))) # brain masking left eye tissue

#' # Statistical analysis
#' ## Linear models
mdls = distancedf %>%
  group_by(RegionalName) %>%
  group_map(~{
    mdl = lm(MeanDistance ~ Consensus.Diagnosis + Age + Sex,
             data = .x) %>%
      tidy %>%
      cbind(RegionalName = .y)
  }) %>%
  bind_rows(
    distancedf %>%
      group_by(RegionalName) %>%
      group_map(~{
        mdl = lm(MeanDistance ~ Age * Consensus.Diagnosis + Sex,
                 data = .x) %>%
          tidy %>%
          cbind(RegionalName = .y)
      }) %>%
      bind_rows %>% 
      filter(term %>% str_detect('Age:'))
  )

mdlcdsummary = mdls %>% 
  filter(term %in% c('Consensus.DiagnosisMCI', 'Consensus.DiagnosisAD')) %>%
  mutate(term = gsub("Consensus.Diagnosis", "", term),
         term = plyr::mapvalues(term, from = c('MCI',
                                               'AD'),
                                to = c('MCI-CU',
                                       'AD-CU')) %>%
           as.factor %>% 
           fct_relevel('AD-CU', after = Inf)) %>% 
  group_by(term) %>% 
  summarise(mest = mean(estimate), 
            sdest = sd(estimate))

orderInfoRegions = distancedf %>% 
  pivot_wider(names_from = 
                Consensus.Diagnosis, 
              values_from = MeanDistance) %>% 
  group_by(RegionalName) %>% 
  summarise(MmC = mean(MCI, na.rm = T) - mean(CU, na.rm = T), 
            AmC = mean(AD, na.rm = T) - mean(CU, na.rm = T), .groups = 'drop') %>% 
  inner_join(mdls %>% 
               filter(term %>% str_detect('Consensus'))) %>% select(-estimate, -std.error, -statistic) %>% 
  pivot_wider(names_from = term, 
              values_from = p.value) %>% 
  mutate(annotationcd = case_when(
    ((Consensus.DiagnosisMCI <= 0.05) & (Consensus.DiagnosisMCI <= 0.05)) ~ '**',
    (Consensus.DiagnosisMCI <= 0.05) ~ '*',
    (Consensus.DiagnosisAD <= 0.05) ~ '*',
    T ~ 'n.s.'
  ), 
  annotationagexcd = case_when(
    ((`Age:Consensus.DiagnosisMCI` <= 0.05) & (`Age:Consensus.DiagnosisMCI` <= 0.05)) ~ '**',
    (`Age:Consensus.DiagnosisMCI` <= 0.05) ~ '*',
    (`Age:Consensus.DiagnosisAD` <= 0.05) ~ '*',
    T ~ 'n.s.'
  )
  )

#' ## Consensus diagnosis (bar plot)
#+ fig.width=9.0, fig.height=5.75, warning=F
p = mdls %>%
  filter(term %in% c('Consensus.DiagnosisMCI', 'Consensus.DiagnosisAD')) %>%
  mutate(term = gsub("Consensus.Diagnosis", "", term),
         term = plyr::mapvalues(term, from = c('MCI',
                                               'AD'),
                                to = c('MCI-CU',
                                       'AD-CU')) %>%
           as.factor %>% 
           fct_relevel('AD-CU', after = Inf)) %>%
  ggplot(aes(x = reorder_within(RegionalName, 
                                estimate, term), 
             y = estimate)) +
  geom_point() +
  geom_text(aes(label = ifelse(p.value <= 0.05, "*", ""),
                y = estimate + std.error), vjust = 0.5,
            size = 8) +
  geom_errorbar(aes(ymin = estimate - std.error,
                    ymax = estimate + std.error),
                width = 0.5) +
  geom_hline(yintercept = 0,
             linetype = 2) +
  facet_rep_grid(. ~ term, scales = "free") +
  geom_hline(data = mdlcdsummary, aes(yintercept = mest)) +
  geom_hline(data = mdlcdsummary, aes(yintercept = mest - sdest), linetype = 'longdash') + 
  geom_hline(data = mdlcdsummary, aes(yintercept = mest + sdest), linetype = 'longdash') + 
  geom_rect(data = mdlcdsummary, 
            aes(xmin = -Inf, xmax = Inf, 
                ymin = mest - sdest, 
                ymax = mest + sdest), alpha = 0.1, 
            fill = 'gray64', 
            inherit.aes = F) +
  gtheme +
  labs(x = "",
       y = TeX("$\\Delta$ \\[RGPL\\] a.u.")) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_reordered() +
  guides(color = guide_legend(label.position = 'right')) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), size = 15),
        axis.title.y = element_text(margin = margin(0, -15, 0, 0), size = 20),
        panel.spacing.x = unit(-2.5, 'lines'),
        legend.title = element_blank(),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), size = 20),
        panel.grid.major.x = element_blank(),
        strip.text.x = element_text(size = 20, margin = margin(0, 0, -1.7, 0)))
p
ggsave(paste0(figroot, 'ConsensusDiagnosisEffects_ASNR2021', '.pdf'),
       width = 6.7,
       height = 5.3,
       p)

#' ## Consensus diagnosis (density)
#+ fig.width=12.25, fig.height=5.0, warning=F
p = orderInfoRegions %>% right_join(distancedf) %>%  mutate(RegionalName = RegionalName %>% as.factor) %>%
  ggplot(aes(x = MeanDistance,
             color = Consensus.Diagnosis)) +
  geom_density() +
  geom_text(aes(label = annotationcd,
                x = Inf,
                y = Inf), 
            hjust = 1,
            size = 16,
            vjust = 1,
            inherit.aes = F) +
  facet_rep_wrap(fct_reorder(RegionalName, AmC, .desc = T) ~ .,
                 ncol = 5) +
  scale_x_log10(breaks = c(0.5, 0.75, 1.0, 1.3, 1.75, 2.0, 3.0)) + annotation_logticks(sides = 'b', scaled = T) +
  labs(x = 'Relative geodesic path length (RGPL) a.u.',
       y = TeX('Sample density $a.u.^-^1$')) +
  geom_vline(data = distancedf %>%
               group_by(RegionalName,
                        Consensus.Diagnosis) %>%
               summarise(MD2 = mean(MeanDistance, na.rm = T)) %>% 
               ungroup() %>% left_join(orderInfoRegions),
             aes(xintercept = MD2,
                 color = Consensus.Diagnosis),
             linetype = 2,
             size = 1.0) +
  scale_color_manual(values = c('Blue', 'Orange', 'Red')) +
  gtheme +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 11,
                                   margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 12),
        panel.spacing.x = unit(-2.0, "lines"),
        panel.spacing.y = unit(-2.5, 'lines'),
        axis.ticks.length = unit(-0.25, "cm"),
        legend.margin = margin(0, 0, 0, 0), 
        legend.box.margin = margin(5, -10, -10, -10),
        axis.title.y = element_text(margin = margin(0, -15, 0, 0)),
        axis.title.x = element_text(margin = margin(-10, 0, 0, 0)),
        axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")))
p
ggsave(paste0(figroot, 'DistributionsOrdered_ConsensusDiagnosis_RGND_ASNR2021', '.pdf'),
       width = 12.25,
       height = 5.0,
       p)

#' ## Age x consensus diagnosis
#+ fig.width=8.0, fig.height=6.75, warning=F
p = orderInfoRegions %>% right_join(distancedf) %>% ggplot(aes(x = Age, y = MeanDistance, color = Consensus.Diagnosis)) +
  geom_point(size = 1, shape = 21, fill = NA) +
  geom_smooth(method = 'lm', alpha = 0.1, size = 0.5) +
  geom_text(aes(label = '',
                x = Inf,
                y = Inf), 
            size = 8,
            hjust = 1, 
            vjust = 1,
            inherit.aes = F) +
  facet_rep_wrap(fct_reorder(RegionalName, AmC, .desc = T) ~ .,
                 ncol = 5) +
  labs(x = 'Age y',
       y = 'Relative geodesic path length (RGPL) a.u.') +
  gtheme +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.margin = margin(0, 0, -5, 0), 
        legend.box.margin = margin(5, -10, -10, -10),
        strip.text.x = element_text(size = 10, margin = margin(0, 0, 0, 0)),
        panel.spacing.x = unit(-1.5, "lines"),
        panel.spacing.y = unit(-0.5, 'lines')) +
  scale_color_manual(values = c('Blue', 'Orange', 'Red')) 
p
ggsave(paste0(figroot, 'AgeXConsensusDiagnosis_RGND_ASNR2021', '.pdf'),
       width = 8.0,
       height = 6.75,
       p)

#' # Subject level visualization
# Setting up data and some visualization parameters
col_fun = colorRamp2(c(0, 1, 2), c('cyan', 'white', 'red'))
wideDistanceDF = distancedf %>%
  mutate(Sex = Sex %>% plyr::mapvalues(from = c('Male', 'Female'),
                                       to = c('M', 'F'))) %>% 
  mutate(CDSex = paste0(Consensus.Diagnosis, '-',
                        Sex) %>% as.factor %>% 
           fct_relevel('CU-F',
                       'CU-M',
                       'MCI-F',
                       'MCI-M',
                       'AD-F',
                       'AD-M')) %>% 
  select(SubjectID, Age, Sex, Consensus.Diagnosis, CDSex,
         RegionalName, MeanDistance) %>%
  pivot_wider(names_from = RegionalName, values_from = MeanDistance)
meanDistanceMatrix =  wideDistanceDF %>%
  select(-SubjectID, -Age, -Sex, -Consensus.Diagnosis, -CDSex) %>% 
  data.matrix
# rownames(meanDistanceMatrix) = wideDistanceDF$SubjectID %>% gsub('ADCP_', '', .)
rownames(meanDistanceMatrix) = meanDistanceMatrix %>% nrow %>% seq(1, .) %>% str_pad(3, pad = '0')
ageColFun = colorRamp2(c(min, median, max) %>% map_dbl(rlang::exec, wideDistanceDF$Age, na.rm = T), 
                       topo.colors(3)) 
ha = rowAnnotation(
  Age = anno_simple(wideDistanceDF$Age, col = ageColFun))
lgd_age = Legend(title = "Age y", 
                 legend_height = unit(10, "cm"),
                 col_fun = ageColFun, at = seq(55, 90, length.out = 10) %>% round(2), 
                 labels = seq(55, 90, length.out = 10) %>% round(2) %>% as.character,
                 border = 'black')

#+ fig.width=3.5, fig.height=11, warning=F
p = Heatmap(meanDistanceMatrix,
            rect_gp = gpar(col = "white", 
                           lwd = 0.5), 
            border = T, col = col_fun,
            column_names_gp = gpar(fontsize = 8), 
            row_names_gp = gpar(fontsize = 4),
            show_column_dend = F, show_row_dend = F, 
            cluster_rows = F, cluster_columns = F,
            row_split = wideDistanceDF$CDSex,
            row_order = order(wideDistanceDF$Age),
            column_order = order(orderInfoRegions$AmC, decreasing = T),
            right_annotation = ha,
            left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c(3,2)),
                                                             labels = wideDistanceDF %>% group_by(CDSex) %>% 
                                                               summarise(n = n()) %>% ungroup() %>% pull(n) %>% 
                                                               map_chr(~paste('n =', .x)), 
                                                             labels_gp = gpar(col = "white", fontsize = 12))),
            heatmap_legend_param = 
              list(legend_height = unit(10, "cm"), 
                   title = expression("RGPL a.u"), 
                   at = seq(0, 2, by = 0.2), 
                   border = "black", 
                   col_fun = col_fun))
draw(p, annotation_legend_list = list(lgd_age))
pdf(paste0(figroot, 'SubjectLevel_ASNR2021.pdf'),
    width = 3.5,
    height = 11)
draw(p, annotation_legend_list = list(lgd_age))
dev.off()