# play with various Morlet wavelets -- illustrates usee in geophysical processing

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# --- Wavelet primitives ---
morlet_wavelet <- function(t, a = 1, b = 0, w0 = 6) {
  u <- (t - b) / a
  c0 <- pi^(-1/4)
  (1/sqrt(a)) * c0 * exp(1i * w0 * u) * exp(-u^2 / 2)
}


# --- Simulation grid ---
TIME_STEPS <- 4096

t <- seq(-5, 5, length.out = TIME_STEPS)

len_t <- TIME_STEPS

# make some wavelets
mor1 <- morlet_wavelet(t, a=0.01, b=0, w0=6)
mor2 <- morlet_wavelet(t, a=0.1, b=0, w0=6)
mor3 <- morlet_wavelet(t, a=1, b=0, w0=6)
mor4 <- morlet_wavelet(t, a=10, b=0, w0=6)

df_mor1 <- data.frame(
  t = t,
  remor = Re(mor1),
  immor = Im(mor1),
  value = "a = 0.01, b = 0, w0 = 6"
)

df_mor2 <- data.frame(
  t = t,
  remor = Re(mor2),
  immor = Im(mor2),
  value = "a = 0.1, b = 0, w0 = 6"
)

df_mor3 <- data.frame(
  t = t,
  remor = Re(mor3),
  immor = Im(mor3),
  value = "a = 1, b = 0, w0 = 6"
)

df_mor4 <- data.frame(
  t = t,
  remor = Re(mor4),
  immor = Im(mor4),
  value = "a = 10, b = 0, w0 = 6"
)

df <- bind_rows(
  df_mor1,
  df_mor2,
  df_mor3,
  df_mor4
)
## Plot ##

ggplot(
  df, aes(x = t)
) + geom_line(aes(y = remor, colour = value), linetype = 1) + 
  geom_line(aes(y = immor, colour = value), linetype = 2) + 
  facet_wrap(~value)

