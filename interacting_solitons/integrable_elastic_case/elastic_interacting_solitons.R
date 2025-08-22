# ============================================================
# Integrable NLS (focusing): two-soliton elastic collision
# PDE: i psi_t = - (1/2) psi_xx - |psi|^2 psi
# Outputs:
#   - nls_collision_space.mp4     (|psi(x,t)| vs x)
#   - nls_collision_spectrum.mp4  (|hat{psi}(k,t)|^2 vs k)
# Robust: guards against NA/empty, safer dt, even PNG sizes, yuv420p MP4
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(av)
})

# ---------- Helpers ----------
fftshift_idx <- function(n){
  if (n %% 2 != 0) stop("fftshift_idx: n must be even")
  c((n/2 + 1):n, 1:(n/2))
}
make_k_from_x <- function(x){
  n <- length(x); if (n %% 2 != 0) stop("make_k_from_x: even length required")
  dx <- mean(diff(x))
  L  <- (max(x) - min(x) + dx)/2     # grid is [-L, L)
  (pi / L) * c(0:(n/2 - 1), -n/2:-1)
}
safe_range <- function(v){
  r <- range(v, na.rm = TRUE)
  if (!all(is.finite(r))) r <- c(0, 1)
  if (diff(r) == 0) r <- r + c(-1, 1) * 1e-6
  r
}

# ---------- Domain ----------
L  <- 80                  # half-domain: x in [-L, L)
N  <- 1024                # even (FFT-friendly)
x  <- seq(-L, L, length.out = N + 1); x <- x[-(N+1)]
dx <- x[2] - x[1]
k  <- make_k_from_x(x)
k2 <- k^2

# ---------- Two bright solitons, well-separated & counter-propagating ----------
bright <- function(x, eta, v, x0, phi0 = 0){
  eta / cosh(eta * (x - x0)) * exp(1i * (v * x + phi0))
}
eta1 <- 1.0; v1 <- +1.0; x01 <- -20; phi1 <- 0
eta2 <- 0.8; v2 <- -0.8; x02 <- +20; phi2 <- 0
psi  <- bright(x, eta1, v1, x01, phi1) + bright(x, eta2, v2, x02, phi2)

# Guard: finite initial field?
if (!all(is.finite(Re(psi))) || !all(is.finite(Im(psi)))) stop("Initial psi not finite")

# ---------- Time stepping (Strang split-step Fourier) ----------
dt         <- 0.05   # safer than 0.001 if your FFT lib is touchy
Tend       <- 60
steps      <- ceiling(Tend / dt)
save_every <- 50       # ~300 frames total at smaller dt

# Linear propagator: ψ̂(t+dt)=ψ̂(t) * exp(-i * (k^2/2) dt)
Lfac <- exp(-1i * 0.5 * k2 * dt)

# ---------- Simulate & collect frames ----------
frames_space <- list()
frames_spec  <- list()
fid <- 1L

add_frame <- function(psi, fid){
  if (!all(is.finite(Re(psi))) || !all(is.finite(Im(psi)))) return(invisible(FALSE))
  y_space <- Mod(psi)
  Ahat <- fft(psi) * dx
  Pw   <- Mod(Ahat)^2
  if (!all(is.finite(y_space)) || !all(is.finite(Pw))) return(invisible(FALSE))
  idx <- fftshift_idx(length(Pw))
  frames_space[[fid]] <<- data.frame(frame = fid, x = x, y = y_space)
  frames_spec [[fid]] <<- data.frame(frame = fid, k = k[idx], y = Pw[idx])
  TRUE
}

ok <- add_frame(psi, fid); fid <- fid + 1L
if (isFALSE(ok)) stop("Failed to add initial frame; psi may be invalid")

for (n in 1:steps){
  # Nonlinear half-step
  amp2 <- Mod(psi)^2
  if (!all(is.finite(amp2))) break
  psi <- psi * exp(1i * amp2 * (dt/2))
  
  # Linear full-step in Fourier
  psi_hat <- fft(psi)
  psi_hat <- psi_hat * Lfac
  psi <- fft(psi_hat, inverse = TRUE) / N
  
  # Nonlinear half-step
  amp2 <- Mod(psi)^2
  if (!all(is.finite(amp2))) break
  psi <- psi * exp(1i * amp2 * (dt/2))
  
  # Save frame
  if (n %% save_every == 0L){
    add_frame(psi, fid); fid <- fid + 1L
  }
}

frames_space <- frames_space[!vapply(frames_space, is.null, logical(1))]
frames_spec  <- frames_spec [!vapply(frames_spec,  is.null, logical(1))]

if (length(frames_space) == 0L || length(frames_spec) == 0L) {
  stop("No frames captured — try reducing dt further or check your R FFT setup.")
}

df_space <- bind_rows(frames_space)
df_spec  <- bind_rows(frames_spec)

# Guards: data sanity
if (nrow(df_space) == 0L || nrow(df_spec) == 0L) {
  stop("Empty data after binding frames.")
}
if (all(!is.finite(df_space$y))) stop("SPACE data all non-finite")
if (all(!is.finite(df_spec$y)))  stop("SPECTRUM data all non-finite")

cat("Frames captured:", length(unique(df_space$frame)),
    " | rows(space):", nrow(df_space), " rows(spec):", nrow(df_spec), "\n")

# ---------- Write PNG frames (even pixel sizes) ----------
# width*DPI and height*DPI must be even: 10*150=1500, 8*150=1200 (both even)
w_in <- 10; h_in <- 8; dpi <- 150

dir_space <- "frames_nls_space"
dir_spec  <- "frames_nls_spec"
for (d in c(dir_space, dir_spec)) {
  if (dir.exists(d)) unlink(d, recursive = TRUE, force = TRUE)
  dir.create(d)
}

# Fixed y-limits across frames (with guards)
ylim_space <- safe_range(df_space$y)
ylim_spec  <- safe_range(df_spec$y)
use_log_spec <- FALSE  # set TRUE to plot log10 spectrum

# Save SPACE frames
uf <- sort(unique(df_space$frame))
for (f in uf){
  d <- df_space[df_space$frame == f, , drop = FALSE]
  p <- ggplot(d, aes(x = x, y = y)) +
    geom_line(linewidth = 0.9) +
    labs(title = sprintf("Integrable NLS: two-soliton collision — frame %d", f),
         x = "x", y = "|psi(x)|") +
    coord_cartesian(ylim = ylim_space) +
    theme_minimal(base_size = 14)
  ggsave(file.path(dir_space, sprintf("frame_%04d.png", f)),
         p, width = w_in, height = h_in, dpi = dpi, bg = "white")
}

# Save SPECTRUM frames
uf2 <- sort(unique(df_spec$frame))
if (use_log_spec) {
  df_spec$y_plot <- log10(df_spec$y + 1e-12)
  ylim_spec_plot <- safe_range(df_spec$y_plot)
} else {
  df_spec$y_plot <- df_spec$y
  ylim_spec_plot <- ylim_spec
}
for (f in uf2){
  d <- df_spec[df_spec$frame == f, , drop = FALSE]
  p <- ggplot(d, aes(x = k, y = y_plot)) +
    geom_line(linewidth = 0.9) +
    labs(title = sprintf("Spectrum |hat{psi}(k)|^2 — frame %d", f),
         x = "k", y = if (use_log_spec) "log10 power" else "power") +
    coord_cartesian(ylim = ylim_spec_plot) +
    theme_minimal(base_size = 14)
  ggsave(file.path(dir_spec, sprintf("frame_%04d.png", f)),
         p, width = w_in, height = h_in, dpi = dpi, bg = "white")
}

pngs_space <- list.files(dir_space, pattern = "^frame_\\d+\\.png$", full.names = TRUE)
pngs_space <- pngs_space[order(pngs_space)]
pngs_spec  <- list.files(dir_spec,  pattern = "^frame_\\d+\\.png$", full.names = TRUE)
pngs_spec  <- pngs_spec[order(pngs_spec)]

if (length(pngs_space) == 0L || length(pngs_spec) == 0L) {
  stop("No PNG frames found to encode — check write permissions / ggplot.")
}

cat("Encoding", length(pngs_space), "space frames; ",
    length(pngs_spec), "spectrum frames...\n")

# ---------- Encode MP4s (yuv420p -> widely compatible) ----------
av::av_encode_video(
  input     = pngs_space,
  output    = "nls_collision_space.mp4",
  framerate = 30,
  vfilter   = "format=yuv420p"
)
av::av_encode_video(
  input     = pngs_spec,
  output    = "nls_collision_spectrum.mp4",
  framerate = 30,
  vfilter   = "format=yuv420p"
)

cat("Done -> nls_collision_space.mp4 & nls_collision_spectrum.mp4\n")
