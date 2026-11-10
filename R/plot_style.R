review_colours <- c(ink='#111111',secondary='#777777',guide='#333333')
review_theme <- function(base_size=10) {
 ggplot2::theme_classic(base_size=base_size,base_family='serif') +
  ggplot2::theme(
   text=ggplot2::element_text(colour=review_colours['ink']),
   plot.title=ggplot2::element_text(face='plain',size=base_size,hjust=0,
                                  margin=ggplot2::margin(b=5)),
   strip.text=ggplot2::element_text(face='plain',size=base_size,
                                  margin=ggplot2::margin(t=3,b=5)),
   strip.background=ggplot2::element_blank(),
   axis.title=ggplot2::element_text(size=base_size),
   axis.title.x=ggplot2::element_text(margin=ggplot2::margin(t=5)),
   axis.title.y=ggplot2::element_text(margin=ggplot2::margin(r=5)),
   axis.text=ggplot2::element_text(size=base_size-1,colour=review_colours['ink']),
   axis.line=ggplot2::element_line(linewidth=.3,colour=review_colours['ink']),
   axis.ticks=ggplot2::element_line(linewidth=.3,colour=review_colours['ink']),
   axis.ticks.length=grid::unit(2,'pt'),
   panel.spacing=grid::unit(10,'pt'),
   plot.margin=ggplot2::margin(4,5,4,4),
   legend.position='none')
}
review_pdf_device <- function(filename,width,height,bg='white',...) {
 if(identical(Sys.info()[['sysname']],'Darwin') && capabilities('aqua'))
  grDevices::quartz(file=filename,type='pdf',width=width,height=height,
                    family='serif',bg=bg,...)
 else
  grDevices::pdf(filename,width=width,height=height,family='Times',
                 useDingbats=FALSE,bg=bg,...)
}
save_review_plot <- function(plot,stem,width,height) {
 ggplot2::ggsave(paste0(stem,'.png'),plot,width=width,height=height,dpi=600,bg='white')
 ggplot2::ggsave(paste0(stem,'.pdf'),plot,width=width,height=height,
                 device=review_pdf_device,bg='white')
}
save_review_panels <- function(plots,stem,width,height,widths=rep(1,length(plots))) {
 draw <- function() {
  grid::grid.newpage()
  grobs <- lapply(plots,ggplot2::ggplotGrob)
  stopifnot(length(unique(vapply(grobs,function(g)length(g$heights),0L)))==1L)
  heights <- do.call(grid::unit.pmax,lapply(grobs,function(g)g$heights))
  for(i in seq_along(grobs))grobs[[i]]$heights <- heights
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,length(grobs),widths=widths)))
  for(i in seq_along(grobs)) {
   grid::pushViewport(grid::viewport(layout.pos.row=1,layout.pos.col=i))
   grid::grid.draw(grobs[[i]])
   grid::popViewport()
  }
  grid::popViewport()
 }
 grDevices::png(paste0(stem,'.png'),width=width,height=height,units='in',res=600,bg='white')
 draw();grDevices::dev.off()
 review_pdf_device(paste0(stem,'.pdf'),width=width,height=height)
 draw();grDevices::dev.off()
}

