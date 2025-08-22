
library(ggplot2)
library(dplyr)
library(tidyr)
library(gganimate)
library(av)

# ---- Mothers ----
morlet <- function(x, w0 = 6)  pi^(-1/4) * exp(1i*w0*x) * exp(-x^2/2)
haar   <- function(x) ifelse(x>=0 & x<0.5, 1, ifelse(x>=0.5 & x<1, -1, 0))

# ---- Explicit heat solver (Neumann BC), returns frames for animation ----
diffuse_fd_movie <- function(u0, D, dx, tmax, nframes, cfl = 0.3) {
  s  <- cfl
  dt <- s * dx^2 / D
  steps_total     <- ceiling(tmax / dt)
  steps_per_frame <- max(1L, ceiling(steps_total / nframes))
  true_nframes    <- floor(steps_total / steps_per_frame)
  
  u <- u0
  N <- length(u)
  frames <- vector("list", true_nframes + 1L)
  
  frames[[1]] <- data.frame(xi = seq_len(N), Re = Re(u), Im = Im(u), t = 0)
  
  for (f in 1:true_nframes) {
    for (k in 1:steps_per_frame) {
      left  <- c(u[2], u[1:(N-1)])     # Neumann mirror
      right <- c(u[2:N], u[N-1])
      u <- u + s * (left - 2*u + right)
    }
    frames[[f+1]] <- data.frame(xi = seq_len(N), Re = Re(u), Im = Im(u), t = f * steps_per_frame * dt)
  }
  bind_rows(frames)
}

# ---- Grid & params ----
L <- 6; N <- 801
x <- seq(-L, L, length.out = N); dx <- x[2] - x[1]
D <- 0.15
tmax <- 0.35
nframes <- 90  # target; solver adapts to match step size

# ---- ICs ----
u0_m <- morlet(x, w0 = 6)
u0_h <- haar(x)

# ---- Run solver ----
df_m <- diffuse_fd_movie(u0_m, D, dx, tmax, nframes) %>% mutate(x = x[xi], kind = "Morlet") %>%
  select(-xi)
df_h <- diffuse_fd_movie(u0_h, D, dx, tmax, nframes) %>% mutate(x = x[xi], kind = "Haar") %>%
  select(-xi)

# ---- Tidy for plotting ----
df_long <- bind_rows(
  df_m %>% select(x, t, Re, Im, kind) %>% pivot_longer(c(Re, Im), names_to = "part", values_to = "val"),
  df_h %>% select(x, t, Re, kind) %>% mutate(Im = 0) %>% pivot_longer(c(Re), names_to = "part", values_to = "val")
)

# ---- Animation plot ----
p <- ggplot(df_long, aes(x = x, y = val, color = part)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ kind, ncol = 1, scales = "free_y") +
  labs(
    title = "Wavelet diffusion (heat eqn) — t = {sprintf('%.3f', frame_time)}",
    x = "x", y = "value", color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top") +
  transition_time(t) +
  ease_aes("linear")

# ---- Render with av (MP4) ----
# Requires ffmpeg provided by the 'av' package; gganimate will call av::av_encode_video under the hood.
anim_save(
  filename = "wavelet_diffusion.mp4",
  animation = p,
  renderer = av_renderer(),
  nframes = length(unique(df_long$t)),  # one frame per captured time
  fps = 24, width = 960, height = 720
)

# Optional: also save a lightweight GIF (slower, larger file)
# anim_save("wavelet_diffusion.gif", p, renderer = gifski_renderer(), nframes = 80, fps = 20, width = 800, height = 600)
