# Cell 3 gene panel: pooled mean normalized RNA, with a zero baseline.
# Error bars: cell-level SEM; no biological-replicate inference.
make_my_list_barplots <- function(summary, output_dir = 'results/carm_dr18_two_panels',
                                 colors = c(CARM = '#2878A5', DR18 = '#CE5151'),
                                 prefix = 'my_list', title = 'My List') {
  stopifnot(requireNamespace('ggplot2', quietly = TRUE),
            !anyDuplicated(summary$gene), nrow(summary) > 0L)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!'category' %in% names(summary)) summary$category <- 'Exploration'
  categories <- unique(as.character(summary$category))
  plots <- setNames(vector('list', length(categories)), categories)
  values <- do.call(rbind, lapply(c('CARM', 'DR18'), function(group) {
    data.frame(gene = summary$gene, category = summary$category,
               available = summary$available, group = group,
               mean_normalized = summary[[paste0(group, '_mean_normalized')]],
               sem_normalized = summary[[paste0(group, '_sem_normalized')]],
               sd_normalized = summary[[paste0(group, '_sd_normalized')]],
               n_cells = summary[[paste0(group, '_n_cells')]])
  }))
  values$p_value <- summary$p_value[match(values$gene, summary$gene)]
  values$p_adj_BH <- summary$p_adj_BH[match(values$gene, summary$gene)]
  values$group <- factor(values$group, levels = c('CARM', 'DR18'))
  utils::write.csv(values, file.path(output_dir, paste0(prefix, '_barplot_values.csv')), row.names = FALSE)
  labels <- c(Mono_derived = 'Monocyte-derived', Chemotaxis = 'Chemotaxis',
              Phagocytic_Mac = 'Phagocytic macrophage', MHC_I = 'MHC I',
              MHC_II = 'MHC II', Ag_processing = 'Antigen processing', Co_stim = 'Co-stimulation')
  for (category in categories) {
    panel <- values[values$category == category, ]
    panel$gene <- factor(panel$gene, levels = summary$gene[summary$category == category])
    annotations <- do.call(rbind, lapply(levels(panel$gene), function(g) {
      z <- panel[panel$gene == g, ]
      top <- max(z$mean_normalized + z$sem_normalized, na.rm = TRUE)
      if (!is.finite(top)) top <- 0
      pv <- z$p_value[1]; qv <- z$p_adj_BH[1]
      label <- if (is.na(qv)) 'NA' else if (qv <= 0.0001) '****' else if (qv <= 0.001) '***' else if (qv <= 0.01) '**' else if (qv <= 0.05) '*' else 'ns'
      data.frame(gene = factor(g, levels = levels(panel$gene)), y = if (top > 0) top * 1.25 else 0.08, label = label)
    }))
    missing <- unique(panel[!panel$available, c('gene') , drop = FALSE])
    plot_title <- if (length(categories) == 1L) title else paste0(title, ': ', if (category %in% names(labels)) labels[[category]] else category)
    p <- ggplot2::ggplot(panel[panel$available, ], ggplot2::aes(group, mean_normalized, fill = group)) +
      ggplot2::geom_col(width = 0.65) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_normalized - sem_normalized,
                                        ymax = mean_normalized + sem_normalized),
                             width = 0.18, linewidth = 0.5) +
      ggplot2::geom_text(data = annotations, ggplot2::aes(x = 1.5, y = y, label = label),
                         inherit.aes = FALSE, size = 5, vjust = 1) +
      ggplot2::facet_wrap(~gene, ncol = 4, scales = 'free_y', drop = FALSE) +
      ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
      ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.12))) +
      ggplot2::labs(title = plot_title, x = NULL,
                    y = 'Mean normalized RNA expression',
                    caption = sprintf('Mean +/- cell-level SEM; CARM n=%d, DR18 n=%d cells.\nExploratory cell-level Student t-test; BH within panel. No biological-replicate inference.\nBH-adjusted: ns > 0.05; * <= 0.05; ** <= 0.01; *** <= 0.001; **** <= 0.0001; NA = undefined.\nY-axis range varies by gene.',
                      summary$CARM_n_cells[which(summary$available)[1]],
                      summary$DR18_n_cells[which(summary$available)[1]])) +
      ggplot2::theme_classic(base_size = 14, base_family = 'sans') +
      ggplot2::theme(legend.position = 'none', strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(face = 'italic', size = 16),
                     plot.title = ggplot2::element_text(size = 18),
                     plot.caption = ggplot2::element_text(size = 14, hjust = 0),
                     panel.spacing = grid::unit(1.2, 'lines'))
    if (nrow(missing)) p <- p + ggplot2::geom_text(data = missing,
      ggplot2::aes(x = 1.5, y = 0, label = 'Not available'), inherit.aes = FALSE, vjust = -1)
    plots[[category]] <- p
    height <- 2.8 * ceiling(nlevels(panel$gene) / 4) + 1.5
    for (extension in c('pdf', 'png')) ggplot2::ggsave(
      file.path(output_dir, paste0(prefix, '_', category, '_barplot.', extension)),
      p, width = 14, height = height, units = 'in', dpi = 200, bg = 'white')
  }
  grDevices::pdf(file.path(output_dir, paste0(prefix, '_all_genes_barplots.pdf')), width = 14, height = 2.8 * ceiling(max(table(summary$category)) / 4) + 1.5)
  tryCatch(for (p in plots) print(p), finally = grDevices::dev.off())
  plots
}

# All 26 genes in one faceted figure, matching the signature plot layout in Cell 9.
make_combined_my_list_barplot <- function(summary, output_dir, colors) {
  values <- do.call(rbind,lapply(c('CARM','DR18'),function(g) data.frame(
    gene=factor(summary$gene,levels=summary$gene),group=factor(g,levels=c('CARM','DR18')),
    mean=summary[[paste0(g,'_mean_normalized')]],sem=summary[[paste0(g,'_sem_normalized')]])))
  q <- summary$p_adj_BH
  stars <- ifelse(is.na(q),'NA',ifelse(q<=0.0001,'****',ifelse(q<=0.001,'***',
             ifelse(q<=0.01,'**',ifelse(q<=0.05,'*','ns')))))
  top <- pmax(summary$CARM_mean_normalized+summary$CARM_sem_normalized,
              summary$DR18_mean_normalized+summary$DR18_sem_normalized)
  annotations <- data.frame(gene=factor(summary$gene,levels=summary$gene),
                            y=ifelse(top>0,top*1.22,0.08),label=stars)
  p <- ggplot2::ggplot(values,ggplot2::aes(group,mean,fill=group))+
    ggplot2::geom_hline(yintercept=0,linewidth=0.3,color='grey60')+
    ggplot2::geom_col(width=0.65)+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-sem,ymax=mean+sem),width=0.18,linewidth=0.5)+
    ggplot2::geom_text(data=annotations,ggplot2::aes(x=1.5,y=y,label=label),
      inherit.aes=FALSE,size=5,vjust=1)+
    ggplot2::facet_wrap(~gene,ncol=4,scales='free_y')+
    ggplot2::scale_fill_manual(values=colors)+
    ggplot2::scale_y_continuous(limits=c(0,NA),expand=ggplot2::expansion(mult=c(0,0.15)))+
    ggplot2::labs(title='My List: CARM vs DR18',x=NULL,y='Mean normalized RNA expression +/- SEM',
      caption='CARM n=165; DR18 n=273 cells. Error bars: cell-level SEM.\nExploratory unpaired Student t-tests; BH adjusted within this panel.\nns > 0.05; * <= 0.05; ** <= 0.01; *** <= 0.001; **** <= 0.0001. No biological-replicate inference.')+
    ggplot2::theme_classic(base_size=14)+
    ggplot2::theme(legend.position='none',strip.background=ggplot2::element_blank(),
      strip.text=ggplot2::element_text(size=13,face='italic'),
      plot.caption=ggplot2::element_text(hjust=0,size=11),panel.spacing=grid::unit(1.3,'lines'))
  for(ext in c('png','pdf')) ggplot2::ggsave(file.path(output_dir,paste0('my_list_combined_barplot.',ext)),
    p,width=14,height=20.5,dpi=200,bg='white')
  p
}
