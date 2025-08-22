library(ggplot2)
library(dplyr)
library(tidyr)

# --- Wavelet primitives ---
morlet_wavelet <- function(t, a = 1, b = 0, w0 = 6) {
  u <- (t - b) / a
  c0 <- pi^(-1/4)
  (1/sqrt(a)) * c0 * exp(1i * w0 * u) * exp(-u^2 / 2)
}

mexican_hat_wavelet <- function(t, a = 1, b = 0) {
  u <- (t - b) / a
  c0 <- 2 / (sqrt(3) * pi^(1/4))
  (1/sqrt(a)) * c0 * (1 - u^2) * exp(-u^2 / 2)
}

haar_wavelet <- function(t, a = 1, b = 0) {
  u <- (t - b) / a
  psi0 <- ifelse(u >= 0 & u < 0.5,  1,
                 ifelse(u >= 0.5 & u < 1, -1, 0))
  (1/sqrt(a)) * psi0
}

# --- Simulation grid ---
t <- seq(-4, 4, length.out = 4096)

# --- Build dataframe for ggplot ---
df <- bind_rows(
  data.frame(t = t, value = Re(morlet_wavelet(t, a=0.6, b=0.5)), part = "RealMor", wavelet = "Morlet"),
  data.frame(t = t, value = Im(morlet_wavelet(t, a=0.6, b=0.5)), part = "ImagMor", wavelet = "Morlet"),
  data.frame(t = t, value = mexican_hat_wavelet(t, a=1.2, b=-0.6), part = "Mexican", wavelet = "Mexican hat"),
  data.frame(t = t, value = haar_wavelet(t, a=0.8, b=0), part = "Haar", wavelet = "Haar")
)

# --- Plot ---
ggplot(df, aes(x = t, y = value, color = part)) +
  geom_line() +
  facet_wrap(~wavelet, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top") +
  labs(title = "Canonical Wavelets",
       x = "t",
       y = "ψ(t)")
