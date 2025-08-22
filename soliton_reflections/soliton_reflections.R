# Travelling bright soliton with reflecting barriers & per-bounce loss (alpha)
# Writes two MP4s (space & spectrum). Frame sizes are even → no ffmpeg yuv420p errors.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(av)
})

# ---------- Helpers ----------
fftshift_idx <- function(n){ if(n%%2) stop("need even n"); c((n/2 + 1):n, 1:(n/2)) }
make_k_from_x <- function(x){
  n <- length(x); if(n%%2) stop("need even length x")
  dx <- mean(diff(x)); L <- (max(x) - min(x) + dx)/2
  (pi/L)*c(0:(n/2 - 1), -n/2:-1)
}

# Analytic bright soliton on a segment (eta, v, xs, ts, phi_s)
soliton_segment <- function(x, t, eta, v, xs, ts, phi_s){
  tau  <- t - ts
  xc   <- xs + v * tau
  env  <- eta / cosh(eta * (x - xc))
  phase <- v * (x - xs - 0.5 * v * tau) - (eta^2 - 0.5 * v^2) * tau + phi_s
  env * exp(1i * phase)
}

# Phase alignment at reflection so psi is continuous at the wall center
align_phase <- function(psi_old, x, x_hit){
  idx <- which.min(abs(x - x_hit))
  z   <- psi_old[idx]
  if (Mod(z) == 0) 0 else Arg(z)
}

# Mirror center + velocity + loss if at wall
reflect_if_needed <- function(xc, v, eta, xL, xR, alpha){
  hit <- FALSE; wall <- NA; xc_new <- xc; v_new <- v; eta_new <- eta
  if (xc >= xR){
    d <- xc - xR; xc_new <- xR - d; v_new <- -abs(v); eta_new <- (1 - alpha)*eta; hit <- TRUE; wall <- "R"
  } else if (xc <= xL){
    d <- xL - xc; xc_new <- xL + d; v_new <-  abs(v); eta_new <- (1 - alpha)*eta; hit <- TRUE; wall <- "L"
  }
  list(hit=hit, wall=wall, xc=xc_new, v=v_new, eta=eta_new)
}

# ---------- Domain, barriers, soliton params ----------
L <- 50
N <- 1024                       # even
x <- seq(-L, L, length.out = N + 1); x <- x[-(N + 1)]
dx <- x[2] - x[1]
k  <- make_k_from_x(x)

xL <- -L + 6
xR <-  L - 6

eta0  <- 1.0
v0    <- 10.0                    # >0 means moving right
x0    <- -25
phi0  <- 0
alpha <- 0.20                   # fraction of L2 mass lost per bounce

# Timeline
dt         <- 0.01
Tend       <- 120
steps      <- ceiling(Tend/dt)
save_every <- 10

# Segment state (updates on reflections)
eta   <- eta0
v     <- v0
xs    <- x0
ts    <- 0.0
phi_s <- phi0

# ---------- Simulate (analytic between reflections) & collect frames ----------
frames_space <- list()
frames_spec  <- list()
fid <- 1L

add_frame <- function(psi, fid){
  frames_space[[fid]] <<- data.frame(frame=fid, x=x, y=Mod(psi))
  Fx <- fft(psi) * dx
  Pk <- Mod(Fx)^2
  idx <- fftshift_idx(length(Pk))
  frames_spec[[fid]] <<- data.frame(frame=fid, k=k[idx], y=Pk[idx])
}

t <- 0.0
psi <- soliton_segment(x, t, eta, v, xs, ts, phi_s)
add_frame(psi, fid); fid <- fid + 1L

for (n in 1:steps){
  t_next  <- t + dt
  xc_next <- xs + v * (t_next - ts)
  refl    <- reflect_if_needed(xc_next, v, eta, xL, xR, alpha)
  
  if (refl$hit){
    # exact hit time: xs + v*(thit - ts) = x_wall
    x_w  <- if (refl$wall == "R") xR else xL
    thit <- ts + (x_w - xs)/v
    
    # frame at the instant of hit
    psi <- soliton_segment(x, thit, eta, v, xs, ts, phi_s)
    if (n %% save_every == 0) { add_frame(psi, fid); fid <- fid + 1L }
    
    # start new segment at the wall with updated params
    psi_before <- psi
    xs    <- x_w
    ts    <- thit
    eta   <- refl$eta
    v     <- refl$v
    phi_s <- align_phase(psi_before, x, xs)  # continuity at wall center
    
    # continue to t_next under new segment
    psi <- soliton_segment(x, t_next, eta, v, xs, ts, phi_s)
  } else {
    psi <- soliton_segment(x, t_next, eta, v, xs, ts, phi_s)
  }
  
  t <- t_next
  if (n %% save_every == 0) { add_frame(psi, fid); fid <- fid + 1L }
  if (eta < 1e-3) break
}

df_space <- bind_rows(frames_space)
df_spec  <- bind_rows(frames_spec)

cat("Frames:", length(unique(df_space$frame)),
    "| rows(space):", nrow(df_space), "rows(spec):", nrow(df_spec), "\n")

# ---------- Write PNG frames (EVEN pixel sizes) ----------
# width*DPI and height*DPI must both be even. 10*150=1500 (even), 8*150=1200 (even).
w_in <- 10; h_in <- 8; dpi <- 150

dir_space <- "frames_soliton_space_even"
dir_spec  <- "frames_soliton_spec_even"
for (d in c(dir_space, dir_spec)) {
  if (dir.exists(d)) unlink(d, recursive=TRUE, force=TRUE)
  dir.create(d)
}

# fixed y-limits (avoid flicker)
ylim_space <- range(df_space$y); if (diff(ylim_space)==0) ylim_space <- ylim_space + c(-1,1)*1e-6
ylim_spec  <- range(df_spec$y);  if (diff(ylim_spec)==0)  ylim_spec  <- ylim_spec  + c(-1,1)*1e-6

uf <- sort(unique(df_space$frame))
for (f in uf){
  ds <- df_space[df_space$frame==f,]
  p  <- ggplot(ds, aes(x=x, y=y)) +
    geom_line(linewidth=0.9) +
    geom_vline(xintercept=c(xL,xR), linetype=2) +
    coord_cartesian(ylim=ylim_space) +
    labs(title=sprintf("Travelling soliton — frame %d", f), x="x", y="|psi(x)|") +
    theme_minimal(base_size=14)
  ggsave(file.path(dir_space, sprintf("frame_%04d.png", f)),
         p, width=w_in, height=h_in, dpi=dpi, bg="white")
}

uf2 <- sort(unique(df_spec$frame))
for (f in uf2){
  ds <- df_spec[df_spec$frame==f,]
  p  <- ggplot(ds, aes(x=k, y=y)) +
    geom_line(linewidth=0.9) +
    coord_cartesian(ylim=ylim_spec) +
    labs(title=sprintf("Spectrum |psi_hat(k)|^2 — frame %d", f), x="k", y="power") +
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

# ---------- Encode MP4s (yuv420p; even sizes so no error) ----------
av::av_encode_video(
  input     = pngs_space,
  output    = "soliton_space.mp4",
  framerate = 30,
  vfilter   = "format=yuv420p"
)
av::av_encode_video(
  input     = pngs_spec,
  output    = "soliton_spectrum.mp4",
  framerate = 30,
  vfilter   = "format=yuv420p"
)

cat("Done -> soliton_space.mp4 & soliton_spectrum.mp4\n")
