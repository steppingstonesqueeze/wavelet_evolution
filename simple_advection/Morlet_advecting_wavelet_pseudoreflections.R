# Travelling packet (Morlet-like), reflecting at two walls with per-bounce loss alpha.
# Pseudo reflections that cause discontinjuities in the wave function
# Two MP4s: (1) space amplitude |psi(x)|, (2) spectrum |psi_hat(k)|^2
# Same animation pattern as the earlier diffusion movie: explicit frames + transition_manual() + av.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(gganimate); library(av)
})

#---------- Helpers ----------
fftshift_idx <- function(n){ if(n%%2) stop("need even n"); c((n/2+1):n, 1:(n/2)) }
make_k_from_x <- function(x){
  n <- length(x); if(n%%2) stop("need even length x")
  dx <- mean(diff(x)); L <- (max(x)-min(x)+dx)/2
  (pi/L)*c(0:(n/2-1), -n/2:-1)
}

#---------- Grid & packet ----------
L <- 50; N <- 1024
x <- seq(-L, L, length.out=N+1); x <- x[-(N+1)]; dx <- x[2]-x[1]
k <- make_k_from_x(x)

# Morlet-like travelling packet (not solving NLS here; just clean transport + reflections)
w0 <- 6
psi0 <- pi^(-1/4) * exp(1i*w0*(x+25)/6) * exp(-((x+25)/6)^2/2)  # centered at x=-25, to the right

# transport speed (units of x per unit time)
c_speed <- 1.0
dt <- 0.01
Tend <- 60
steps <- ceiling(Tend/dt)
save_every <- 4

# walls at +/- (L-6), lose alpha of L2 mass per bounce
xL <- -L+6; xR <- L-6
alpha <- 0.20
bounce_margin <- 0.2

#---------- Advection by exact fractional shift (linear interp), then reflect with loss ----------
shift_linear <- function(u, s){
  # shift by s grid points (can be fractional). periodic for transport step
  n <- length(u)
  i0 <- floor(s); frac <- s - i0
  idxA <- ((1:n - i0 - 1) %% n) + 1
  idxB <- ((idxA) %% n) + 1
  (1-frac)*u[idxA] + frac*u[idxB]
}

reflect_with_loss <- function(u, x, xL, xR, alpha){
  # reflect any "overflow" beyond walls back inside and apply sqrt(1-alpha) loss.
  # We do it geometrically on the grid.
  scale <- sqrt(1 - alpha)
  # left overflow
  maskL <- x < xL
  if(any(maskL)){
    x_over <- x[maskL]
    x_ref  <- 2*xL - x_over       # mirror
    u_ref  <- rev(u[maskL]) * scale
    # add reflected into nearest grid points
    for(j in seq_along(x_over)){
      # find nearest grid index to x_ref[j]
      idx <- which.min(abs(x - x_ref[j]))
      u[idx] <- u[idx] + u_ref[j]
    }
    u[maskL] <- 0
  }
  # right overflow
  maskR <- x > xR
  if(any(maskR)){
    x_over <- x[maskR]
    x_ref  <- 2*xR - x_over
    u_ref  <- rev(u[maskR]) * scale
    for(j in seq_along(x_over)){
      idx <- which.min(abs(x - x_ref[j]))
      u[idx] <- u[idx] + u_ref[j]
    }
    u[maskR] <- 0
  }
  u
}

#---------- Simulate & collect explicit frames ----------
psi <- psi0
frames_space <- list(); frames_spec <- list()
fid <- 1L

add_frames <- function(psi, fid){
  # space
  frames_space[[fid]] <<- data.frame(frame=fid, x=x, y=Mod(psi))
  # spectrum
  Fx <- fft(psi)*dx
  Pk <- Mod(Fx)^2
  idx <- fftshift_idx(length(Pk))
  frames_spec[[fid]] <<- data.frame(frame=fid, k=k[idx], y=Pk[idx])
}

add_frames(psi, fid); fid <- fid + 1L

for(n in 1:steps){
  # exact transport by shift s = c*dt/dx grid cells
  s <- c_speed*dt/dx
  psi <- shift_linear(psi, s)
  
  # reflect with loss at walls
  psi <- reflect_with_loss(psi, x, xL, xR, alpha)
  
  if(n %% save_every == 0){
    add_frames(psi, fid); fid <- fid + 1L
  }
}

df_space <- bind_rows(frames_space)
df_spec  <- bind_rows(frames_spec)

cat("Frames:", length(unique(df_space$frame)), "  rows(space):", nrow(df_space),
    "  rows(spec):", nrow(df_spec), "\n")

#---------- VIDEO 1: |psi(x)| ----------
p_space <- ggplot(df_space, aes(x=x, y=y)) +
  geom_line(linewidth=0.9) +
  geom_vline(xintercept=c(xL,xR), linetype=2) +
  labs(title="Travelling wavelet with lossy reflections — frame {current_frame}",
       x="x", y="|psi(x)|") +
  theme_minimal(base_size=14)

anim_space <- p_space + transition_manual(frame)

anim_save(
  "travelling_wavelet_space.mp4",
  animation = anim_space,
  renderer  = av_renderer(),
  nframes   = length(unique(df_space$frame)),
  fps       = 30, width=1000, height=450
)

#---------- VIDEO 2: |psi_hat(k)|^2 ----------
p_spec <- ggplot(df_spec, aes(x=k, y=y)) +
  geom_line(linewidth=0.9) +
  labs(title="Spectrum |psi_hat(k)|^2 — frame {current_frame}",
       x="k", y="power") +
  theme_minimal(base_size=14)

anim_spec <- p_spec + transition_manual(frame)

anim_save(
  "travelling_wavelet_spectrum.mp4",
  animation = anim_spec,
  renderer  = av_renderer(),
  nframes   = length(unique(df_spec$frame)),
  fps       = 30, width=1000, height=450
)
