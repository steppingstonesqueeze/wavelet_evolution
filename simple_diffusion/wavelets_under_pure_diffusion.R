# Wavelets under diffusion using a closed form convolution kernel for computing
# wavelet amplitudes

# --- Diffusing Morlet & Haar wavelets with ggplot2 ---
library(ggplot2)
library(dplyr)
library(tidyr)

# ===== Wavelets (mothers) =====
morlet <- function(x, w0 = 6) {
  # Complex Morlet (unit L2): psi(x) = pi^(-1/4) * exp(i w0 x) * exp(-x^2/2)
  c0 <- pi^(-1/4)
  c0 * exp(1i * w0 * x) * exp(-x^2/2)
}

haar <- function(x) {
  # Haar on [0,1): +1 on [0,0.5), -1 on [0.5,1), 0 otherwise
  ifelse(x >= 0 & x < 0.5,  1,
         ifelse(x >= 0.5 & x < 1, -1, 0))
}

# ===== Heat kernel and exact (linear) convolution on a grid =====
heat_kernel <- function(x, D, t) {
  if (t <= 0) return(as.numeric(x == 0)) # Dirac on-grid for t=0 (handled separately anyway)
  1/sqrt(4*pi*D*t) * exp(-(x^2)/(4*D*t))
}

# Linear convolution (no wrap-around) with grid spacing dx.
# Returns (g * u)(x) where g is centered (symmetric) on 0, sampled on the same grid length.
linconv_same <- function(u, g_full, dx) {
  # 'u' length N, 'g_full' length (2N-1), centered at index N
  # open convolution length = 3N-2; take the middle slice N : (2N-1)
  y_open <- convolve(u, rev(g_full), type = "open") * dx
  N <- length(u)
  y_open[N:(2*N-1)]
}

# Build a centered kernel sampled on differences {-N+1,...,0,...,N-1}*dx
centered_kernel <- function(N, dx, kernel_fun) {
  # positions for the full kernel
  k <- seq(-(N-1), (N-1))
  xk <- k * dx
  g_full <- kernel_fun(xk)
  # normalize discretely so sum ~ 1
  g_full / sum(g_full * dx)
}

# ===== Simulation grid & parameters =====
L  <- 6            # half-width of domain
N  <- 2001         # odd length (so we have a center)
x  <- seq(-L, L, length.out = N)
dx <- x[2] - x[1]

D  <- 0.15         # diffusion coefficient
times <- c(0.00, 0.03, 0.10, 0.30)  # diffusion times

# ===== Initial conditions =====
psi_m0 <- morlet(x, w0 = 6)        # complex
psi_h0 <- haar(x)                  # real

# ===== Diffuse each IC for requested times =====
diffuse_series <- function(psi0, is_complex = FALSE) {
  out_list <- list()
  for (t in times) {
    if (t == 0) {
      # no change
      if (is_complex) {
        val_re <- Re(psi0); val_im <- Im(psi0)
      } else {
        val_re <- psi0;    val_im <- rep(0, length(psi0))
      }
    } else {
      # build centered heat kernel for this t
      g_full <- centered_kernel(N, dx, function(z) heat_kernel(z, D, t))
      if (is_complex) {
        val_re <- linconv_same(Re(psi0), g_full, dx)
        val_im <- linconv_same(Im(psi0), g_full, dx)
      } else {
        val_re <- linconv_same(psi0,     g_full, dx)
        val_im <- rep(0, length(psi0))
      }
    }
    out_list[[length(out_list)+1]] <- data.frame(
      x = x, t = t, Re = val_re, Im = val_im
    )
  }
  bind_rows(out_list)
}

df_m <- diffuse_series(psi_m0, is_complex = TRUE)  %>%
  mutate(kind = "Morlet")

df_h <- diffuse_series(psi_h0, is_complex = FALSE) %>%
  mutate(kind = "Haar")

# Stack & tidy for plotting
df_long <- bind_rows(df_m, df_h) %>%
  mutate(t = factor(t, levels = times, labels = paste0("t = ", times))) %>%
  pivot_longer(cols = c(Re, Im), names_to = "part", values_to = "value")

# ===== Plots =====

# Morlet: show real & imaginary parts over time
p_m <- df_long %>%
  filter(kind == "Morlet") %>%
  ggplot(aes(x = x, y = value, color = part)) +
  geom_line() +
  facet_wrap(~ t, ncol = 1, scales = "free_y") +
  labs(title = "Diffusion of Morlet wavelet (heat equation)",
       x = "x", y = "value") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# Haar: single real signal over time
p_h <- df_long %>%
  filter(kind == "Haar", part == "Re") %>%
  ggplot(aes(x = x, y = value)) +
  geom_line() +
  facet_wrap(~ t, ncol = 1, scales = "free_y") +
  labs(title = "Diffusion of Haar wavelet (heat equation)",
       x = "x", y = "value") +
  theme_minimal(base_size = 13)

print(p_m)
Sys.sleep(5)
print(p_h)

# --- (Optional) energy decay check under diffusion ---
energy <- function(v) sum(v^2) * dx
E_m <- df_m |>
  group_by(t) |>
  summarise(E = energy(Re) + energy(Im))
E_h <- df_h |>
  group_by(t) |>
  summarise(E = energy(Re))  # Im is zero

print(E_m)
print(E_h)
