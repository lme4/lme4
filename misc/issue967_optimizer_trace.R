## issue967_optimizer_trace.R
##
## GH #967: structured-covariance glmer (cs(1 + age | district)) fails from the
## default starting value but succeeds from start = c(0.1, 0, 0.1).
##
## This script does two things:
##   (1) reproduces Ben Bolker's deviance-surface plot over the theta grid
##       (par1, par2 = the two SD-scale params; par3 = the correlation), and
##   (2) the "to do" from the issue thread: runs *instrumented* optimizations
##       from the two starting points, recording the trajectory of every theta
##       value the optimizer evaluates, then overlays those paths on the surface
##       to show where each run goes (and why the default start gets stuck).
##
## The instrumentation wraps the deviance function and drives it with lme4's own
## optwrap()/bobyqa -- i.e. exactly the optimizer glmer uses for nAGQ = 0 -- so
## the recorded trajectory is the real optimization path, not a re-implementation.

suppressMessages(devtools::load_all(quiet = TRUE))   # or: library(lme4)
library(ggplot2)
theme_set(theme_bw())

data("Contraception", package = "mlmRev")

## ---- model + deviance function ------------------------------------------
## nAGQ = 0: devfun() is a function of theta only (length 3 here).
devfun <- glmer(use ~ cs(1 + age | district), Contraception, binomial,
                nAGQ = 0, devFunOnly = TRUE)
rho   <- environment(devfun)
lower <- rho$lower                              # c(0, 0, -1)
upper <- rho$upper %||% rep(Inf, length(lower)) # c(Inf, Inf, 1)

## ---- (1) deviance surface over the theta grid ---------------------------
## Coarser than Ben's 31x31x25 by default so the script runs in a few seconds;
## bump `n` back up for a publication-quality landscape.
## par1/par2 floor extended down to 10^-4: the optimizer drives the SDs toward
## the boundary 0, so a -2 floor left the surface undrawn exactly where the
## trajectories go (the white band that looked like a "vanishing gradient").
n <- c(30, 30, 9)
dd_par <- expand.grid(
  par1 = 10^seq(-4, 1, length.out = n[1]),
  par2 = 10^seq(-4, 1, length.out = n[2]),
  par3 = round(seq(-0.9, 0.9, length.out = n[3]), 3))

message(sprintf("evaluating deviance on %d grid points ...", nrow(dd_par)))
system.time(dvec <- apply(dd_par, 1, devfun))
dd_par$dev  <- dvec
dd_par$dev2 <- dvec - min(dvec) + 0.01           # shift so log scale is happy
par3_levels <- sort(unique(dd_par$par3))

## ---- (2) instrumented optimization --------------------------------------
## Wrap devfun so every evaluation appends (theta, dev) to a log, then optimize
## with the same machinery glmer uses (optwrap + bobyqa, adj = FALSE for nAGQ 0).
trace_optim <- function(devfun, start, lower, upper, label) {
  log <- new.env()
  log$rows <- list()
  traced <- function(par) {
    val <- devfun(par)
    log$rows[[length(log$rows) + 1L]] <-
      c(par1 = par[1], par2 = par[2], par3 = par[3], dev = val)
    val
  }
  opt <- lme4:::optwrap("bobyqa", traced, start,
                        lower = lower, upper = upper,
                        control = list(), adj = FALSE, verbose = 0L)
  traj <- as.data.frame(do.call(rbind, log$rows))
  names(traj) <- c("par1", "par2", "par3", "dev")
  traj$step  <- seq_len(nrow(traj))
  traj$start <- label
  list(traj = traj, opt = opt)
}

starts <- list(
  "default: (1, 1, 0)"   = c(1, 1, 0),
  "working: (0.1, 0, 0.1)" = c(0.1, 0, 0.1))

runs <- Map(function(s, lab) trace_optim(devfun, s, lower, upper, lab),
            starts, names(starts))

## report outcome of each run
for (r in runs) {
  o <- r$opt
  cat(sprintf("\nstart %-22s : %d evals, final dev = %.3f, conv = %s\n",
              r$traj$start[1], nrow(r$traj), o$fval, o$conv))
  cat("   final theta:", paste(signif(o$par, 4), collapse = ", "), "\n")
}

traj <- do.call(rbind, lapply(runs, `[[`, "traj"))

## snap each trajectory point to the nearest par3 facet level so it can be
## drawn on top of the faceted raster.
traj$par3_facet <- par3_levels[
  vapply(traj$par3, function(p) which.min(abs(par3_levels - p)), integer(1))]

## starting points (for emphasis) and final points
start_pts <- do.call(rbind, lapply(runs, function(r) r$traj[1, ]))
start_pts$par3_facet <- par3_levels[
  vapply(start_pts$par3, function(p) which.min(abs(par3_levels - p)), integer(1))]
final_pts <- do.call(rbind, lapply(runs, function(r) r$traj[nrow(r$traj), ]))
final_pts$par3_facet <- par3_levels[
  vapply(final_pts$par3, function(p) which.min(abs(par3_levels - p)), integer(1))]

## Clamp trajectory coords to the evaluated grid for plotting.  The optimizer
## pushes the SDs to the boundary par = 0, i.e. log10(par) -> -Inf, which has no
## location on a log axis (you can't draw -Inf).  So we pin any point at/below
## the grid floor to the floor value and flag it (offgrid = TRUE).
##
## NB: the floor itself (here log10 = -4) is NOT meaningful -- it is just the
## smallest par we evaluated the surface on (par = 10^seq(-4, 1, ...) above), i.e.
## the bottom edge of the canvas.  Widen the grid to 10^seq(-6, 1) and the same
## boundary points would sit at -6 instead.  A clamped y-coord is therefore a
## placeholder, NOT the value; the true par (often exactly 0) is kept in `traj`
## and the saved .rds.  The `offgrid` flag is what carries the real information
## ("this point is at/below the floor, heading to the singular boundary 0").
flo1 <- min(dd_par$par1); flo2 <- min(dd_par$par2)
add_plot_coords <- function(d) {
  d$par1_plot <- pmax(d$par1, flo1)
  d$par2_plot <- pmax(d$par2, flo2)
  d$offgrid   <- d$par1 < flo1 | d$par2 < flo2
  d
}
traj      <- add_plot_coords(traj)
start_pts <- add_plot_coords(start_pts)
final_pts <- add_plot_coords(final_pts)

## ======================= HOW TO READ THESE PLOTS =========================
## The deviance is a function of THREE params: theta = (par1, par2, par3), where
## par1, par2 are SD-scale params and par3 is the correlation.  A heatmap can
## only show two, so we plot deviance over (log10 par1, log10 par2) and handle
## par3 separately:
##   - faceted plot  : one panel per par3 level (a slice of the 3-D surface);
##                     each evaluated theta is drawn in the panel nearest its par3.
##   - slice plot    : a single par3 = 0 panel (the faceted plot already shows the
##                     par3 dimension), used here to compare the two runs head on.
##
## COLOUR = optimizer run in both plots: each run (default vs working start) has
## ONE colour shared by its start circle, its evaluated points, and its final
## star, so you can always tell which markers belong to the same run.
##
## Marker key (same in both plots):
##   white circle, run-coloured ring = that run's STARTING theta
##   star (asterisk), run-coloured   = that run's FINAL theta (where it stopped)
##   small dot / triangle, run-col.  = every other theta the optimizer evaluated
##   black x over a point            = "off grid": its par1 or par2 is at/below the
##                                     floor (heading to the singular boundary 0),
##                                     so its plotted height is a clamped
##                                     placeholder, NOT the real value (clamp note).
## Points are evaluation samples, NOT a connected descent path (see NB below).
## =========================================================================
## facet_wrap(~par3) facets on the column named `par3`, so for the trajectory
## layers we feed in copies whose `par3` IS the snapped facet level (otherwise
## each continuous trajectory par3 would spawn its own empty facet).
to_facet <- function(d) { d$par3 <- d$par3_facet; d }

## NB: no connecting lines/arrows.  These are the thetas bobyqa *evaluated* in
## call order -- the first several are its initial coordinate stencil, not a
## descent step -- so drawing a directed path would overstate the structure.
p <- ggplot(dd_par, aes(log10(par1), log10(par2))) +
  geom_raster(aes(fill = dev2)) +
  geom_contour(aes(z = log10(dev2)), colour = "grey40", alpha = 0.5) +
  facet_wrap(~par3, labeller = label_both) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(trans = "log10", name = "dev - min") +
  ## each evaluated theta, placed in the facet whose par3 is nearest its value
  geom_point(data = to_facet(traj),
             aes(log10(par1_plot), log10(par2_plot), colour = start,
                 shape = offgrid), size = 1, alpha = 0.85) +
  ## start = white-filled circle with a RUN-COLOURED ring (so each start matches
  ## its own evaluated points and final star), final = run-coloured star.
  geom_point(data = to_facet(start_pts),
             aes(log10(par1_plot), log10(par2_plot), colour = start),
             shape = 21, fill = "white", size = 3.5, stroke = 1.4) +
  geom_point(data = to_facet(final_pts),
             aes(log10(par1_plot), log10(par2_plot), colour = start),
             shape = 8, size = 3.5, stroke = 1.2) +
  scale_colour_manual(values = c("default: (1, 1, 0)" = "red",
                                 "working: (0.1, 0, 0.1)" = "cyan"),
                      name = "optimizer run") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                     labels = c(`FALSE` = "on grid", `TRUE` = "SD->0 (clamped)"),
                     name = NULL) +
  labs(title = "GH #967: glmer cs(1 + age | district) -- evaluated thetas",
       subtitle = paste("colour = optimizer run; circle = start, star = final, dot = evaluated theta",
                        "(incl. bobyqa's stencil), snapped to nearest par3 facet"))

ggsave(file.path("misc", "issue967_devsurf_trajectories.png"), p,
       height = 8, width = 10)

## ---- par3 = 0 slice; everything coloured by optimizer RUN ----------------
## Each run gets ONE colour shared by its start circle, its evaluated points, and
## its final star (and a matching point shape), so you can see at a glance which
## markers belong together -- the two runs really do take different routes.
## (par3 drift off this slice is shown in the faceted plot above; here the focus
## is run identity.)  No connecting lines: these are evaluation samples, not a walk.
slice0 <- subset(dd_par, par3 == par3_levels[which.min(abs(par3_levels))])
p2 <- ggplot(slice0, aes(log10(par1), log10(par2))) +
  geom_raster(aes(fill = dev2)) +
  geom_contour(aes(z = log10(dev2)), colour = "grey40", alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(trans = "log10", name = "dev - min") +
  geom_point(data = traj,
             aes(log10(par1_plot), log10(par2_plot), colour = start, shape = start),
             size = 1.8) +
  geom_point(data = subset(traj, offgrid),               # black x: clamped (SD->0)
             aes(log10(par1_plot), log10(par2_plot)),
             shape = 4, colour = "black", size = 1.3, stroke = 0.7) +
  geom_point(data = start_pts,                            # white-filled run-coloured ring
             aes(log10(par1_plot), log10(par2_plot), colour = start),
             shape = 21, fill = "white", size = 4, stroke = 1.5) +
  geom_point(data = final_pts,                            # run-coloured star
             aes(log10(par1_plot), log10(par2_plot), colour = start),
             shape = 8, size = 4, stroke = 1.3) +
  scale_colour_manual(values = c("default: (1, 1, 0)" = "red",
                                 "working: (0.1, 0, 0.1)" = "cyan"),
                      name = "optimizer run") +
  scale_shape_manual(values = c("default: (1, 1, 0)" = 16,
                                "working: (0.1, 0, 0.1)" = 17),
                     name = "optimizer run") +
  labs(title = "GH #967: optimizer-evaluated thetas projected onto par1/par2",
       subtitle = paste("background = deviance at par3 = 0; colour = optimizer run;",
                        "circle = start, star = final, x = clamped (SD->0, singular)"))

ggsave(file.path("misc", "issue967_devsurf_slice.png"), p2, height = 7, width = 10)

saveRDS(list(grid = dd_par, traj = traj, runs = runs),
        file.path("misc", "issue967_optimizer_trace.rds"))

message("done: wrote misc/issue967_devsurf_trajectories.png, ",
        "misc/issue967_devsurf_slice.png, misc/issue967_optimizer_trace.rds")
