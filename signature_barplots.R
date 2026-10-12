# Exploratory cell-level Student tests of Seurat AddModuleScore.
signature_statistics <- function(data, score_columns, output_dir, prefix) {
  summaries <- tests <- list()
  for (signature in score_columns) {
    x <- data[[signature]][data$group == 'CARM']
    y <- data[[signature]][data$group == 'DR18']
    fit <- tryCatch(stats::t.test(x, y, paired = FALSE, var.equal = TRUE,
                                 alternative = 'two.sided'), error = function(e) NULL)
    pooled_sd <- sqrt(((length(x)-1)*sd(x)^2+(length(y)-1)*sd(y)^2)/(length(x)+length(y)-2))
    sed <- pooled_sd * sqrt(1/length(x)+1/length(y))
    tests[[signature]] <- data.frame(signature=signature, CARM_n=length(x), DR18_n=length(y),
      mean_difference_CARM_minus_DR18=mean(x)-mean(y), SED=sed,
      t_statistic=if(is.null(fit)) NA_real_ else unname(fit$statistic),
      df=if(is.null(fit)) NA_real_ else unname(fit$parameter),
      p_value=if(is.null(fit)) NA_real_ else fit$p.value,
      method='Exploratory cell-level two-sided unpaired Student t-test; equal variances')
    summaries[[signature]] <- do.call(rbind, lapply(c('CARM','DR18'), function(g) {
      scores <- data[[signature]][data$group==g]
      data.frame(signature=signature, group=g, n_cells=length(scores), mean_score=mean(scores),
                 sd_score=sd(scores), sem_score=sd(scores)/sqrt(length(scores)), SED=sed)
    }))
  }
  tests <- do.call(rbind,tests); rownames(tests) <- NULL
  tests$p_adj_BH <- p.adjust(tests$p_value,'BH')
  summary <- do.call(rbind,summaries); rownames(summary) <- NULL
  write.csv(tests,file.path(output_dir,paste0(prefix,'_student_t_tests.csv')),row.names=FALSE)
  write.csv(summary,file.path(output_dir,paste0(prefix,'_barplot_values.csv')),row.names=FALSE)
  list(summary=summary,tests=tests)
}

plot_signature_statistics <- function(result, output_dir, prefix, colors,
                                      title = 'IL-18-related signature scores', ncol = 2) {
  summary <- result$summary; tests <- result$tests
  summary$group <- factor(summary$group, levels=c('CARM','DR18'))
  order <- unique(summary$signature)
  summary$signature <- factor(summary$signature,levels=order)
  annotations <- do.call(rbind,lapply(order,function(s) {
    z <- summary[summary$signature==s,]
    lo <- min(0,z$mean_score-z$sem_score); hi <- max(0,z$mean_score+z$sem_score)
    pv <- tests$p_value[tests$signature==s]
    data.frame(signature=factor(s,levels=order), y=hi+0.22*(hi-lo),
               label=if(is.na(pv)) 'p = NA' else sprintf('p = %.3g',pv))
  }))
  labels <- c(WP_IL18_SIGNALING='WikiPathways IL-18 signaling',
    REACTOME_INTERLEUKIN_18_SIGNALING='Reactome IL-18 signaling',
    GOBP_RESPONSE_TO_INTERLEUKIN_18='Response to IL-18',
    GOBP_INTERLEUKIN_18_PRODUCTION='IL-18 production',
    REACTOME_THE_NLRP3_INFLAMMASOME='NLRP3 inflammasome',
    REACTOME_INFLAMMASOMES='Inflammasomes',
    REACTOME_INTERFERON_GAMMA_SIGNALING='Interferon-gamma signaling')
  p <- ggplot2::ggplot(summary,ggplot2::aes(group,mean_score,fill=group))+
    ggplot2::geom_hline(yintercept=0,linewidth=0.3,color='grey60')+
    ggplot2::geom_col(width=0.65)+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_score-sem_score,ymax=mean_score+sem_score),width=0.18,linewidth=0.5)+
    ggplot2::geom_text(data=annotations,ggplot2::aes(x=1.5,y=y,label=label),inherit.aes=FALSE,size=4,vjust=1)+
    ggplot2::facet_wrap(~signature,ncol=ncol,scales='free_y',labeller=ggplot2::as_labeller(labels))+
    ggplot2::scale_fill_manual(values=colors,drop=FALSE)+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0.08,0.15)))+
    ggplot2::labs(title=title,x=NULL,y='Mean AddModuleScore +/- SEM',
      caption=sprintf('CARM n=%d; DR18 n=%d cells. Error bars: cell-level SEM.\nRaw p: two-sided unpaired Student t-test, equal variances.\nExploratory cell-level analysis; no biological-replicate inference.',
        summary$n_cells[summary$group=='CARM'][1],summary$n_cells[summary$group=='DR18'][1]))+
    ggplot2::theme_classic(base_size=14)+
    ggplot2::theme(legend.position='none',strip.background=ggplot2::element_blank(),
      strip.text=ggplot2::element_text(size=13),plot.caption=ggplot2::element_text(hjust=0,size=11),
      panel.spacing=grid::unit(1.3,'lines'))
  width <- if(length(order)==1) 6.5 else 9.6
  height <- if(length(order)==1) 5.5 else 11.2
  for(ext in c('pdf','png')) ggplot2::ggsave(file.path(output_dir,paste0(prefix,'_barplot.',ext)),p,width=width,height=height,dpi=200,bg='white')
  p
}
