# Figure 1 layout approved by the author on 8 September 2026.
# Eight DGM nodes and ten directed edges; no numerical simulation is run.
make_figure1 <- function(root) {
  library(grid)
  stopifnot(requireNamespace("ragg", quietly = TRUE),
            requireNamespace("jsonlite", quietly = TRUE))
  out <- file.path(root, "figures")
  qa_dir <- file.path(root, "verification", "figure1")
  dir.create(qa_dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  stem <- "Figure1_DGM_observation_boundary"
  W <- 1320; H <- 1320; width_in <- 6; height_in <- width_in * H / W
  ink <- "#344653"; blue <- "#40647B"; orange <- "#A56819"
  latent_fill <- "#F0F5F8"; recorded_fill <- "#FFF2DC"
  nodes <- data.frame(
    id = c("U", "X", "S", "Tc", "To", "C", "Obs", "Cens"),
    x = c(250, 650, 650, 250, 650, 1030, 650, 1080),
    y = c(240, 240, 600, 600, 880, 600, 1140, 1140),
    w = c(240, 190, 220, 260, 260, 260, 340, 240),
    h = c(130, 130, 110, 110, 110, 110, 120, 110),
    label = c("U\nhealth reserve\n(unmeasured)", "X\nage proxy\n(measured)",
              "S*\ntrue strategy", "Cancer-death\ntime", "Other-cause\ndeath time",
              "C\nsummary code", "Observed endpoint\nand follow-up", "Independent\ncensoring"),
    recorded = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE)
  edge <- function(id, from, to, xy, dashed = FALSE, accent = FALSE) {
    list(id = id, from = from, to = to, xy = matrix(xy, ncol = 2, byrow = TRUE),
         dashed = dashed, accent = accent)
  }
  edges <- list(
    edge("F1-E01", "U", "S", c(370,265,570,545)),
    edge("F1-E02", "X", "S", c(650,305,650,545)),
    edge("F1-E03", "U", "Tc", c(250,305,250,545)),
    edge("F1-E04", "U", "To", c(250,175,250,110,1230,110,1230,880,780,880)),
    edge("F1-E05", "S", "Tc", c(540,600,380,600)),
    edge("F1-E06", "S", "To", c(650,655,650,825), TRUE),
    edge("F1-E07", "S", "C", c(760,600,900,600), accent = TRUE),
    edge("F1-E08", "Tc", "Obs", c(250,655,250,1140,480,1140)),
    edge("F1-E09", "To", "Obs", c(650,935,650,1080)),
    edge("F1-E10", "Cens", "Obs", c(960,1140,820,1140)))
  expected <- c("U>S", "X>S", "U>Tc", "U>To", "S>Tc", "S>To", "S>C", "Tc>Obs", "To>Obs", "Cens>Obs")
  stopifnot(nrow(nodes) == 8L, length(edges) == 10L,
            identical(vapply(edges, function(e) paste(e$from, e$to, sep = ">"), ""), expected),
            identical(which(vapply(edges, function(e) e$dashed, logical(1))), 6L))
  on_boundary <- function(p, n) {
    x0 <- n$x - n$w/2; x1 <- n$x + n$w/2
    y0 <- n$y - n$h/2; y1 <- n$y + n$h/2; tol <- 1e-7
    (min(abs(p[1] - c(x0,x1))) < tol && p[2] >= y0 && p[2] <= y1) ||
      (min(abs(p[2] - c(y0,y1))) < tol && p[1] >= x0 && p[1] <= x1)
  }
  # Geometry supports both straight diagonals and orthogonal polylines.
  # Shrinking node interiors excludes intentional boundary contacts.
  hits_interior <- function(a,b,n) {
    lo <- c(n$x-n$w/2,n$y-n$h/2)+1e-6
    hi <- c(n$x+n$w/2,n$y+n$h/2)-1e-6
    tlo <- 0; thi <- 1; d <- b-a
    for (k in 1:2) {
      if (abs(d[k]) < 1e-10) {
        if (a[k] < lo[k] || a[k] > hi[k]) return(FALSE)
      } else {
        tt <- sort(c((lo[k]-a[k])/d[k],(hi[k]-a[k])/d[k]))
        tlo <- max(tlo,tt[1]); thi <- min(thi,tt[2])
      }
    }
    tlo <= thi
  }
  orientation <- function(a,b,c) (b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])
  intersects <- function(a,b,c,d) {
    for (k in 1:2) if(max(min(a[k],b[k]),min(c[k],d[k])) > min(max(a[k],b[k]),max(c[k],d[k]))+1e-8) return(FALSE)
    orientation(a,b,c)*orientation(a,b,d) <= 1e-8 &&
      orientation(c,d,a)*orientation(c,d,b) <= 1e-8
  }
  stopifnot(intersects(c(0,0),c(10,10),c(0,10),c(10,0)),
            !intersects(c(0,0),c(10,10),c(0,1),c(10,11)),
            intersects(c(0,0),c(10,0),c(5,0),c(15,0)),
            !intersects(c(0,0),c(4,0),c(5,0),c(15,0)))
  test_node <- data.frame(x=5,y=5,w=4,h=4)
  stopifnot(hits_interior(c(0,0),c(10,10),test_node),
            !hits_interior(c(0,3),c(10,3),test_node))
  segments <- list(); blockers <- character()
  for (e in edges) {
    stopifnot(on_boundary(e$xy[1,],nodes[nodes$id==e$from,]),
              on_boundary(e$xy[nrow(e$xy),],nodes[nodes$id==e$to,]))
    for (i in seq_len(nrow(e$xy)-1)) {
      a <- e$xy[i,]; b <- e$xy[i+1,]
      segments[[length(segments)+1L]] <- list(id=e$id,a=a,b=b)
      for (j in seq_len(nrow(nodes))) {
        if(hits_interior(a,b,nodes[j,])) blockers <- c(blockers,paste(e$id,"through",nodes$id[j]))
      }
    }
  }
  crossings <- character()
  for (i in seq_along(segments)) for (j in seq_along(segments)) {
    if(j <= i || segments[[i]]$id==segments[[j]]$id) next
    a <- segments[[i]]; b <- segments[[j]]
    if(intersects(a$a,a$b,b$a,b$b)) crossings <- c(crossings,paste(a$id,b$id,sep=" / "))
  }
  stopifnot(length(blockers)==0L,length(crossings)==0L)
  straight_count <- sum(vapply(edges,function(e)nrow(e$xy)==2L,logical(1)))
  bend_count <- sum(vapply(edges,function(e)nrow(e$xy)-2L,integer(1)))
  stopifnot(straight_count==8L,bend_count==4L)
  alignment_pass <- nodes$y[nodes$id=="U"]==nodes$y[nodes$id=="X"] &&
    nodes$y[nodes$id=="Tc"]==nodes$y[nodes$id=="S"]
  obs_node <- nodes[nodes$id=="Obs",]
  tc_path <- edges[[8]]$xy
  tc_entry_pass <- nrow(tc_path)==3L && tc_path[1,1]==tc_path[2,1] &&
    tc_path[2,2]==tc_path[3,2] && tc_path[3,1]==obs_node$x-obs_node$w/2 &&
    tc_path[3,2]==obs_node$y && tc_path[2,1]<tc_path[3,1]
  stopifnot(alignment_pass,tc_entry_pass)
  svg_parts <- character()
  xml <- function(s) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", s, fixed=TRUE), fixed=TRUE), fixed=TRUE)
  put <- function(s) svg_parts <<- c(svg_parts, s)
  rect <- function(x,y,w,h,fill,stroke,lwd=1.25) {
    grid.rect(x/W, 1-y/H, w/W, h/H, gp=gpar(fill=fill,col=stroke,lwd=lwd))
    put(sprintf('<rect x="%g" y="%g" width="%g" height="%g" fill="%s" stroke="%s" stroke-width="%g"/>',
                x-w/2,y-h/2,w,h,fill,stroke,lwd*0.75*W/(width_in*72)))
  }
  label <- function(s,x,y,pt=9,col=ink,bold=FALSE,anchor="middle") {
    grid.text(s,x/W,1-y/H,just=if(anchor=="start") "left" else "centre",
              gp=gpar(fontfamily="Arial",fontsize=pt,col=col,fontface=if(bold)"bold" else"plain"))
    put(sprintf('<text x="%g" y="%g" fill="%s" font-family="Arial, Helvetica, sans-serif" font-size="%g" font-weight="%s" text-anchor="%s" dominant-baseline="central">%s</text>',
                x,y,col,pt*W/(width_in*72),if(bold)"bold" else"normal",anchor,xml(s)))
  }
  draw <- function() {
    svg_parts <<- character(); grid.newpage()
    put(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%gin" height="%gin" viewBox="0 0 %g %g">',width_in,height_in,W,H))
    put('<title>Simulation assumptions and the summary-code observation boundary</title>')
    put('<desc>Eight nodes and ten directed edges. The only dashed edge is S to other-cause death, present in H2 only. No edge connects C to the observed endpoint.</desc>')
    put('<defs><marker id="arrow" markerWidth="7" markerHeight="7" refX="7" refY="3.5" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L7,3.5 L0,7 Z" fill="#344653"/></marker><marker id="accent" markerWidth="7" markerHeight="7" refX="7" refY="3.5" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L7,3.5 L0,7 Z" fill="#A56819"/></marker></defs>')
    put(sprintf('<rect width="%g" height="%g" fill="white"/>',W,H))
    rect(75,42,26,26,latent_fill,blue)
    label("Simulated processes",110,42,10,blue,TRUE,"start")
    rect(755,42,26,26,recorded_fill,orange)
    label("Recorded variables",790,42,10,orange,TRUE,"start")
    for (e in edges) {
      ec <- if(e$accent) orange else ink
      grid.lines(e$xy[,1]/W,1-e$xy[,2]/H,
                 arrow=arrow(length=unit(1.65,"mm"),type="closed"),
                 gp=gpar(col=ec,fill=ec,lwd=1.4,lty=if(e$dashed)"dashed" else"solid",linejoin="round",lineend="butt"))
      put(sprintf('<polyline id="%s" data-from="%s" data-to="%s" points="%s" fill="none" stroke="%s" stroke-width="3.20833" stroke-linejoin="round"%s marker-end="url(#%s)"/>',
                  e$id,e$from,e$to,paste(apply(e$xy,1,paste,collapse=","),collapse=" "),ec,
                  if(e$dashed)' stroke-dasharray="14 10"' else "",if(e$accent)"accent" else"arrow"))
    }
    for (i in seq_len(nrow(nodes))) {
      n <- nodes[i,]; fill <- if(n$recorded) recorded_fill else latent_fill
      rect(n$x,n$y,n$w,n$h,fill,if(n$recorded)orange else blue)
      lines <- strsplit(n$label,"\n",fixed=TRUE)[[1]]
      yy <- n$y + (seq_along(lines)-(length(lines)+1)/2)*35
      for (j in seq_along(lines)) label(lines[j],n$x,yy[j],9.2,
                                       bold=j==1 && n$id %in% c("U","X","S","C"))
    }
    label("eff",460,558,8.6)
    label("H2 only",755,745,8.6)
    label("Se / Sp",830,558,8.6,orange)
    label("Analysis exposure",1030,695,8.5,orange)
    label("Solid arrows: core DGM assumptions; dashed arrow: H2 assumption-violation scenarios.",65,1240,8.0,ink,FALSE,"start")
    label("These are simulation assumptions, not causal effects established in SEER.",65,1285,8.0,ink,FALSE,"start")
    put('</svg>')
  }
  ragg::agg_png(file.path(out,paste0(stem,".png")),width=width_in,height=height_in,units="in",res=300,background="white")
  draw(); dev.off()
  ragg::agg_tiff(file.path(out,paste0(stem,".tiff")),width=width_in,height=height_in,units="in",res=600,compression="lzw",background="white")
  draw(); dev.off()
  grDevices::cairo_pdf(file.path(out,paste0(stem,".pdf")),width=width_in,height=height_in,family="Arial",bg="white")
  draw(); dev.off()
  writeLines(svg_parts,file.path(out,paste0(stem,".svg")),useBytes=TRUE)
  write.csv(nodes,file.path(qa_dir,"figure1_nodes.csv"),row.names=FALSE)
  edge_table <- do.call(rbind,lapply(edges,function(e)data.frame(id=e$id,from=e$from,to=e$to,line=if(e$dashed)"dashed" else"solid",label=switch(e$id,"F1-E05"="eff","F1-E06"="H2 only","F1-E07"="Se / Sp",""))))
  write.csv(edge_table,file.path(qa_dir,"figure1_edges.csv"),row.names=FALSE)
  qa <- list(status="author_approved_2026-09-08",date="2026-09-08",node_count=nrow(nodes),edge_count=length(edges),
             exact_edge_set_pass=TRUE,border_anchoring_pass=TRUE,edge_node_intersections=blockers,
             connector_crossings=crossings,dashed_edge="S>To (H2 only)",C_to_Obs_absent=TRUE,
             straight_edges=straight_count,polyline_bends=bend_count,geometry_helper_tests_pass=TRUE,
             layout="aligned U/X and cancer-death/S rows; Tc to Obs enters left boundary",
             requested_alignment_pass=alignment_pass,
             Tc_to_Obs_single_bend=tc_entry_pass,dimensions_inches=c(width=width_in,height=height_in),png_dpi=300,
             note="Geometry checks only; author approved the v04 layout; export review is recorded separately. No statistical analysis run.")
  jsonlite::write_json(qa,file.path(qa_dir,"figure1_geometry_qa.json"),pretty=TRUE,auto_unbox=TRUE)
  writeLines(capture.output(sessionInfo()),file.path(qa_dir,"figure1_R_sessionInfo.txt"))
  cat("Created Figure 1 bundle from author-approved layout.\n")
}
