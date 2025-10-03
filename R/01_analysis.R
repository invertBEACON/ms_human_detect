## Code for: White et al. (2025) Seeing through the light: Colour contrasts drive human detection in the wild

## --- Clear out --- ##
  rm(list = ls())
  
## --- Libraries --- ##
  library(tidyverse)
  library(patchwork)
  library(lme4)
  library(easystats)
  library(DHARMa)
  library(pavo)
  library(MuMIn)
  library(broom.mixed)
  
## --- Data --- ##
  detect <- read.csv('../data/data_detections.csv')
  target_specs <- as.rspec(read.csv('../data/target_specs.csv'), lim = c(300, 700))
  irrad_specs <- getspec('../data/irrad/', ext = 'JazIrrad', subdir = TRUE, subdir.names = TRUE)

## --- Checks & processing--- ##
  
  # Smooth target specs & remove zeros
  target_specs <- procspec(target_specs, fixneg = 'zero', opt = 'smooth')
  irrad_specs <- irrad2flux(procspec(irrad_specs, fixneg = 'zero', opt = 'smooth'))
  
  # Set reference levels for stats convenience
  detect$colour <- relevel(factor(detect$colour), ref = "C")
  detect$treatment <- relevel(factor(detect$treatment), ref = "C")
  detect$lightness <- relevel(factor(detect$lightness), ref = "low")
  
  # Tidy names of irrad specs
  names(irrad_specs) <- str_remove(names(irrad_specs), "/OUTPUTFILE")
  
  # Tack on weather identifier to irrads
  names(irrad_specs)
  
  # NA outlier irrad_chromas
  detect <- 
    detect |> 
    mutate(irrad_a = replace(irrad_a, irrad_a == "i00_0032", NA),
           irrad_a = replace(irrad_a, irrad_a == "i00_0089", NA))
  
## --- Visual modelling --- #
  
## Segment analysis for irrad
  irrad_group <- 
    detect |> 
      filter(irrad_a %in% names(irrad_specs)[-1]) |> 
      group_by(part_id) |> 
      summarise(weather = unique(weather),
                irrad = unique(irrad_a)) |> 
    rename(irrad_a = irrad) |> 
    select(-part_id)
  
  # Subset specs to only those relevant
  sub_irrad <- subset(irrad_specs, irrad_group$irrad_a)

  # And create a grouping variable
  irrad_groupvar <- irrad_group$weather[match(names(sub_irrad)[-1], irrad_group$irrad_a)]
  
  seg_irrad <- 
    sub_irrad |> 
    vismodel(visual = "segment", achromatic = "all") |> 
    colspace(space = 'segment') |> 
    rownames_to_column(var = "irrad_a") |> 
    rename(irrad_chroma = C,
           irrad_bright = B) |> 
    select(irrad_a, irrad_chroma, irrad_bright)
  
  seg_irrad <- 
    sub_irrad |> 
    summary() |> 
    rownames_to_column(var = "irrad_a") |> 
    rename(irrad_chroma = S5,
           irrad_bright = B2) |> 
    select(irrad_a, irrad_chroma, irrad_bright)
  
  # Calculate means per weather type
  seg_irrad <- left_join(seg_irrad, irrad_group)
  seg_irrad_means <- 
    seg_irrad |> 
    group_by(weather) |> 
    summarise(irrad_chroma = mean(irrad_chroma),
              irrad_bright = mean(irrad_bright))
  
  # Fill missing data with means
  detect_incomplete <- 
    detect |> 
      filter(is.na(irrad_a)) |> 
      left_join(y = seg_irrad_means)
  
  # Add irradiant chroma, brightness to detection data
  detect_complete <- 
    detect |> 
    filter(!is.na(irrad_a)) |> 
    left_join(y = seg_irrad)
  
  # Full set
  detect <- rbind(detect_incomplete, detect_complete)

## Targets  
  # CIELCh model
  vmod_lch <- colspace(vismodel(target_specs, visual = 'cie10', 
                                illum = 'D65'), 
                       'cielch')
  vmod_lch_raw <- vmod_lch
  
  # CIE2000
  vmod_lchdist <- 
    vmod_lch |> 
    coldist(subset = 'bkg') |> 
    rename(treatment = patch1, cie2000 = dS) |> 
    select(treatment, cie2000)
  
  # Add treatment variable
  vmod_lch$treatment <- rownames(vmod_lch)
  
  # Calc diffs
  vmod_lch$L_diff <- vmod_lch$L - filter(vmod_lch, treatment == 'bkg')$L
  vmod_lch$C_diff <- vmod_lch$C - filter(vmod_lch, treatment == 'bkg')$C
  vmod_lch$h_diff <- vmod_lch$h - filter(vmod_lch, treatment == 'bkg')$h
  vmod_lch <- select(vmod_lch, -X, -Y, -Z)
  
  # Combine into single dataset
  detect <- left_join(detect, vmod_lch)
  detect <- left_join(detect, vmod_lchdist)
  
  ## --- Descriptive stats --- ##
  
  # Overall detection probability
  desc_overall <- 
    detect |> 
    summarise(
      mean_detect = mean(detected, na.rm = TRUE),
      sd_detect   = sd(detected, na.rm = TRUE),
      se_detect   = sd_detect / sqrt(n())
    )
  
  # By treatment
  desc_treatment <- 
    detect |> 
    group_by(treatment) |> 
    summarise(
      mean_detect = mean(detected, na.rm = TRUE),
      sd_detect   = sd(detected, na.rm = TRUE),
      se_detect   = sd_detect / sqrt(n())
    )
  
  # By treatment × lightness if useful
  desc_treat_light <- 
    detect |> 
    group_by(treatment, lightness) |> 
    summarise(
      mean_detect = mean(detected, na.rm = TRUE),
      sd_detect   = sd(detected, na.rm = TRUE),
      se_detect   = sd_detect / sqrt(n())
    )
  
  # Quick peek
  desc_overall
  desc_treatment
  desc_treat_light
  
## --- Plots --- #
  
  # Arrange detection frequencies
  detect_freq <- 
    detect |> 
      group_by(part_id, treatment) |> 
      summarise(n = n(),
                detections = sum(detected)) |> 
      mutate(freq = detections / n,
             treatment = factor(treatment, levels = c('G', 'B', 'Y', 'GL', 'BL', 'YL', 'C')))
  
  # Plots
  png('../figs/fig_stimuli.png', width = 28, height = 10, units = 'cm', res = 300)
  par(mfrow = c(1, 3))
    
  plot(target_specs, col = c('forestgreen', 'forestgreen', 
                             'blue', 'blue', 
                             'goldenrod', 'goldenrod', 'brown'), lwd = 1.5)
  
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  plot(vmod_lch_raw, col = c('forestgreen', 'forestgreen', 
                             'blue', 'blue', 
                             'goldenrod', 'goldenrod', 'brown', 'brown'))
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  plot(sub_irrad,
       ylab = c('Irradiance'))
  
  dev.off()
  
  # Treatments
  ggplot(detect_freq, aes(x = treatment, y = freq)) +
    geom_boxplot() +
    geom_jitter(alpha = 0.2, width = 0.2) +
    theme_classic() +
    ylab ('detection frequency')
  
## --- Stats --- ##
  
  # Model
  mod <- glmer(detected ~ scale(L_diff) + scale(C_diff) + scale(h_diff) + scale(distance_m)  + scale(irrad_chroma) + (1 | part_id), 
                family = binomial(),
                na.action = 'na.fail',
                data = detect)
  r2(mod)
  simulateResiduals(mod, plot = TRUE)
  summary(mod)
  
  # Plot
  (plot_hdiff <-   
    ggplot(detect, aes(x = h_diff, y = detected)) +
    geom_jitter(height = 0.03, alpha = 0.1) +
    geom_smooth(method = 'glm') +
    ylab('Detected') +
    xlab('Hue diff.') +
    theme_classic())
  
  (plot_cdiff <-   
    ggplot(detect, aes(x = C_diff, y = detected)) +
      geom_jitter(height = 0.03, alpha = 0.1) +
      geom_smooth(method = 'glm') +
      ylab('') +
      xlab('Chroma diff.') +
      theme_classic())
  
  (plot_ldiff <-   
    ggplot(detect, aes(x = L_diff, y = detected)) +
    geom_jitter(height = 0.03, alpha = 0.1) +
    #geom_smooth(method = 'glm') +
    ylab('') +
    xlab('Luminance diff.') +
    theme_classic())
  
  plot_hdiff + plot_cdiff + plot_ldiff

  # Save
  ggsave('../figs/fig_statmodel.tiff', height = 5, width = 12)
  ggsave('../figs/fig_statmodel.png', height = 5, width = 12)
  
  # Tidy and coefficient plot
  mod_res <- 
    mod |> 
    tidy(conf.int = TRUE) |> 
    filter(term != "(Intercept)", effect == 'fixed') |> 
    mutate(term = factor(term, levels = c('scale(distance_m)', 
                                          'scale(irrad_chroma)', 
                                          'scale(L_diff)', 
                                          'scale(C_diff)', 
                                          'scale(h_diff)')))
  
  (plot_coef <- 
    ggplot(mod_res, aes(x = term, y = estimate)) +
      geom_point() +
      coord_flip() +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
      ylab('beta (95% CI)') +
      ylim(-1.5, 1.5) +
      geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
      scale_x_discrete(labels = c('Distance (m)', 'Irrad. chroma', 'Lum. diff.', 'Chroma diff.', 'Hue diff.')) +
      theme_classic())
  
  # Save
  ggsave('../figs/fig_coefplot.tiff', height = 5, width = 6)
  ggsave('../figs/fig_coefplot.png', height = 5, width = 6)
  