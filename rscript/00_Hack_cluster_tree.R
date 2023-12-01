#!/staging/biology/ls807terra/0_Programs/anaconda3/envs/RNAseq_quantTERRA/bin/Rscript

# This script refine tree on the heatmap
# -----------------------------------------------------------------------------------------------------------------------------------------#
draw_dendrogram <- function(hc, gaps, horizontal = T) {
  # Define equal-length branches
  hc$height <- cumsum(rep(1/length(hc$height), length(hc$height)))
  h = hc$height/max(hc$height)/1.05
  m = hc$merge
  o = hc$order
  n = length(o)
  m[m > 0] = n + m[m > 0]
  m[m < 0] = abs(m[m < 0])
  dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, 
                                                               c("x", "y")))
  dist[1:n, 1] = 1/n/2 + (1/n) * (match(1:n, o) - 1)
  for (i in 1:nrow(m)) {
    dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1])/2
    dist[n + i, 2] = h[i]
  }
  draw_connection = function(x1, x2, y1, y2, y) {
    res = list(x = c(x1, x1, x2, x2), y = c(y1, y, y, y2))
    return(res)
  }
  x = rep(NA, nrow(m) * 4)
  y = rep(NA, nrow(m) * 4)
  id = rep(1:nrow(m), rep(4, nrow(m)))
  for (i in 1:nrow(m)) {
    c = draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], 
                        dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
    k = (i - 1) * 4 + 1
    x[k:(k + 3)] = c$x
    y[k:(k + 3)] = c$y
  }
  x = pheatmap:::find_coordinates(n, gaps, x * n)$coord
  y = unit(y, "npc")
  if (!horizontal) {
    a = x
    x = unit(1, "npc") - y
    y = unit(1, "npc") - a
  }
  res = grid::polylineGrob(x = x, y = y, id = id)
  return(res)
}
# Replace the non-exported function `draw_dendrogram` in `pheatmap`:
assignInNamespace(x="draw_dendrogram", value=draw_dendrogram, ns="pheatmap")
# -----------------------------------------------------------------------------------------------------------------------------------------#
