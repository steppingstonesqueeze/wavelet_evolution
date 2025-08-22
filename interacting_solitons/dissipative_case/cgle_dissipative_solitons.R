# ============================================================
# Dissipative solitons (laser-style CGLE/Haus model): interactions of two pulses
# PDE (dimensionless):
#   ∂_z A = (D + i*β) ∂_{tt} A + (g(E) - ℓ) A - (s + i*γ) |A|^2 A - ν |A|^4 A
#   with saturable gain g(E) = g0 / (1 + E/Esat),  E = ∫ |A|^2 dt
# Numerics: Strang split-step (linear in ω; nonlinear/gain via Heun RK2)
# Outputs:
#   - cgle_space.mp4        (|A(t,z)| vs t)
#   - cgle_spectrum.mp4     (|Â(ω,z)|^2 vs ω)
# Pipeline: even-sized PNG frames -> av::av_encode_video(yuv420p)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(av)
})

# ---------- Helpers ----------
fftshift_idx <- function(n){ if(n%%2!=0) stop("fftshift_idx: need even n"); c((n/2+1):n, 1:(n/2)) }
make_w_from_t <- function(t){
  n <- length(t); if(n%%2!=0) stop("make_w_from_t: need even n")
  dt <- mean(diff(t))
  T  <- (max(t) - min(t) + dt)/2    # grid intended ~[-T, T)
  (pi/T) * c(0:(n/2-1), -n/2:-1)    # angular frequency ω
}
energy <- function(A, dt){ sum(Mod(A)^2) * dt }

# Nonlinear/gain RHS at field A (pointwise vector), given energy-dependent gain
rhs_nl <- function(A, g_eff, s, gamma, nu){
  # (g_eff)A - (s + iγ)|A|^2 A - ν |A|^4 A
  mod2 <- Mod(A)^2
  (g_eff) * A - (s + 1i*gamma) * mod2 * A - nu * (mod2^2) * A
}

# One Heun (RK2) step for the nonlinear/gain part over dz_nl
step_nl_heun <- function(A, dz_nl, dt, g0, Esat, ell, s, gamma, nu){
  E    <- energy(A, dt)
  g1   <- g0 / (1 + E / Esat) - ell
  f1   <- rhs_nl(A, g1, s, gamma, nu)
  Atil <- A + dz_nl * f1
  E2   <- energy(Atil, dt)
  g2   <- g0 / (1 + E2 / Esat) - ell
  f2   <- rhs_nl(Atil, g2, s, gamma, nu)
  A + 0.5 * dz_nl * (f1 + f2)
}

# ---------- Grid (retarded time) ----------
Twin <- 250                        # half-window: t in [-Twin, Twin)
N     <- 1024                     # even (FFT-friendly)
t     <- seq(-Twin, Twin, length.out = N + 1); t <- t[-(N+1)]
dt    <- t[2] - t[1]
w     <- make_w_from_t(t)
w2    <- w^2

# ---------- Model parameters (choose a stable DS window) ----------
# Linear: spectral filtering D (>0), dispersion β (sign per convention: β>0 ~ anomalous)
D      <- 0.02
beta   <- 0.8

# Nonlinear: Kerr γ, saturable absorber s (>0), quintic ν (stabilizer; can be 0)
gamma  <- 1.0
s      <- 0.6
nu     <- 0.02

# Gain/loss: saturable gain g(E) and linear loss ℓ
g0     <- 1.0
Esat   <- 6.0
ell    <- 0.3

# ---------- Propagation (z) ----------
dz         <- 0.005
Zend       <- 200
steps      <- ceiling(Zend / dz)
save_every <- 80                 # ~1000 frames

# Linear propagator in ω (full step)
# ∂_z A = (D + iβ) ∂_{tt} A  ->  Â(z+dz) = Â(z) * exp( - (D + iβ) w^2 dz )
Lfac <- exp( - (D + 1i*beta) * w2 * dz )

# ---------- Initial condition: two pulses (sech) with separation/phase ----------
sech <- function(x) 1 / cosh(x)
A0_1 <- function(t, eta, t0, phi){ eta * sech(eta*(t - t0)) * exp(1i*phi) }
A0_2 <- function(t, eta, t0, phi){ eta * sech(eta*(t - t0)) * exp(1i*phi) }

eta1 <- 10.0; t1 <- -12; phi1 <- 0.0
eta2 <- 9; t2 <- +12; phi2 <- 0.0        # try phi2=pi for anti-bonding dynamics

A <- A0_1(t, eta1, t1, phi1) + A0_2(t, eta2, t2, phi2)

# ---------- Frame storage ----------
frames_space <- list()
frames_spec  <- list()
fid <- 1L

add_frame <- function(A, fid){
  # space magnitude
  frames_space[[fid]] <<- data.frame(frame=fid, t=t, y=Mod(A))
  # spectrum power (FFT * dt for scaling)
  Ahat <- fft(A) * dt
  Pw   <- Mod(Ahat)^2
  idx  <- fftshift_idx(length(Pw))
  frames_spec[[fid]] <<- data.frame(frame=fid, w=w[idx], y=Pw[idx])
}

add_frame(A, fid); fid <- fid + 1L

# ---------- March in z (Strang: NL half + LIN full + NL half) ----------
for (n in 1:steps){
  # Nonlinear + gain/loss half-step (Heun)
  A <- step_nl_heun(A, dz/2, dt, g0, Esat, ell, s, gamma, nu)
  
  # Linear full-step in Fourier
  Ahat <- fft(A)
  Ahat <- Ahat * Lfac
  A    <- fft(Ahat, inverse=TRUE) / N
  
  # Nonlinear + gain/loss half-step again
  A <- step_nl_heun(A, dz/2, dt, g0, Esat, ell, s, gamma, nu)
  
  # Save a frame
  if (n %% save_every == 0L){
    add_frame(A, fid); fid <- fid + 1L
  }
}

frames_space <- frames_space[!vapply(frames_space, is.null, logical(1))]
frames_spec  <- frames_spec [!vapply(frames_spec,  is.null, logical(1))]

df_space <- bind_rows(frames_space)
df_spec  <- bind_rows(frames_spec)

cat("Frames captured:", length(unique(df_space$frame)),
    " | rows(space):", nrow(df_space), " rows(spec):", nrow(df_spec), "\n")

# ---------- Write PNG frames (even pixel sizes) ----------
# width*DPI and height*DPI must be even -> 10*150=1500, 8*150=1200
w_in <- 10; h_in <- 8; dpi <- 150

dir_space <- "frames_cgle_space"
dir_spec  <- "frames_cgle_spec"
for (d in c(dir_space, dir_spec)) {
  if (dir.exists(d)) unlink(d, recursive=TRUE, force=TRUE)
  dir.create(d)
}

# Fixed y-limits across frames (avoid flicker)
ylim_space <- range(df_space$y); if (diff(ylim_space)==0) ylim_space <- ylim_space + c(-1,1)*1e-6
ylim_spec  <- range(df_spec$y);  if (diff(ylim_spec) ==0) ylim_spec  <- ylim_spec  + c(-1,1)*1e-6

# (Optional) log spectrum toggle
use_log_spec <- FALSE
if (use_log_spec) {
  df_spec$y <- log10(df_spec$y + 1e-12)
  ylim_spec <- range(df_spec$y)
}

# Save SPACE frames
uf <- sort(unique(df_space$frame))
for (f in uf){
  d <- df_space[df_space$frame==f, , drop=FALSE]
  p <- ggplot(d, aes(x=t, y=y)) +
    geom_line(linewidth=0.9) +
    labs(title=sprintf("Dissipative solitons (CGLE) — frame %d", f),
         x="t", y="|A(t)|") +
    coord_cartesian(ylim = ylim_space) +
    theme_minimal(base_size=14)
  ggsave(file.path(dir_space, sprintf("frame_%04d.png", f)),
         p, width=w_in, height=h_in, dpi=dpi, bg="white")
}

# Save SPECTRUM frames
uf2 <- sort(unique(df_spec$frame))
for (f in uf2){
  d <- df_spec[df_spec$frame==f, , drop=FALSE]
  p <- ggplot(d, aes(x=w, y=y)) +
    geom_line(linewidth=0.9) +
    labs(title=sprintf("Spectrum |Â(ω)|^2 — frame %d", f),
         x=expression(omega), y=if(use_log_spec) "log10 power" else "power") +
    coord_cartesian(ylim = ylim_spec) +
    theme_minimal(base_size=14)
  ggsave(file.path(dir_spec, sprintf("frame_%04d.png", f)),
         p, width=w_in, height=h_in, dpi=dpi, bg="white")
}

pngs_space <- list.files(dir_space, pattern="^frame_\\d+\\.png$", full.names=TRUE)
pngs_space <- pngs_space[order(pngs_space)]
pngs_spec  <- list.files(dir_spec,  pattern="^frame_\\d+\\.png$", full.names=TRUE)
pngs_spec  <- pngs_spec[order(pngs_spec)]

cat("Encoding", length(pngs_space), "space frames; ",
    length(pngs_spec), "spectrum frames...\n")

# ---------- Encode MP4s (yuv420p -> plays everywhere) ----------
av::av_encode_video(
  input     = pngs_space,
  output    = "cgle_space.mp4",
  framerate = 30,
  vfilter   = "format=yuv420p"
)
av::av_encode_video(
  input     = pngs_spec,
  output    = "cgle_spectrum.mp4",
  framerate = 30,
  vfilter   = "format=yuv420p"
)

cat("Done -> cgle_space.mp4 & cgle_spectrum.mp4\n")
