
review_colours <- c(blue='#3C6078',grey='#666666',ink='#222222')
review_theme <- function(base_size=11) {
 ggplot2::theme_classic(base_size=base_size,base_family='sans') +
  ggplot2::theme(
   plot.title=ggplot2::element_text(face='plain',size=11),
   strip.text=ggplot2::element_text(face='plain'),
   strip.background=ggplot2::element_blank(),
   axis.text=ggplot2::element_text(colour=review_colours['ink']))
}
finalize_review_pdf <- function(path) {
 if(nzchar(Sys.which('gs'))) {
  out <- tempfile(fileext='.pdf')
  grDevices::embedFonts(path,format='pdfwrite',outfile=out,
                       options='-dNoOutputFonts')
  stopifnot(file.exists(out),file.copy(out,path,overwrite=TRUE));unlink(out)
 }
 invisible(path)
}
save_review_plot <- function(plot,stem,width,height) {
 ggplot2::ggsave(paste0(stem,'.png'),plot,width=width,height=height,dpi=450,bg='white')
 ggplot2::ggsave(paste0(stem,'.pdf'),plot,width=width,height=height,
                  device=grDevices::pdf,useDingbats=FALSE,bg='white')
 finalize_review_pdf(paste0(stem,'.pdf'))
}
save_review_panels <- function(plots,stem,width,height,widths=rep(1,length(plots))) {
 draw <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,length(plots),widths=widths)))
  for(i in seq_along(plots))
   print(plots[[i]],vp=grid::viewport(layout.pos.row=1,layout.pos.col=i))
 }
 grDevices::png(paste0(stem,'.png'),width=width,height=height,units='in',res=450,bg='white')
 draw();grDevices::dev.off()
 grDevices::pdf(paste0(stem,'.pdf'),width=width,height=height,bg='white',
                family='Helvetica',useDingbats=FALSE)
 draw();grDevices::dev.off()
 finalize_review_pdf(paste0(stem,'.pdf'))
}
