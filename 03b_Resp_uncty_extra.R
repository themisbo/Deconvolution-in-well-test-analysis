mean_resp <- as.numeric(summarise_all(df_z, list(mean)))

df_data3 <- tibble(tau = seq(log(t0), log(tN), length.out = N_Nodes), mean_resp = mean_resp)

g_formean <- as_tibble(t(exp(z)/t))
mean_g <- as.numeric(summarise_all(g_formean, list(mean)))
z_mean_g <- log(t*mean_g)
df_data4 <- tibble(tau = seq(log(t0), log(tN), length.out = N_Nodes), mean_zg = z_mean_g)

mparam <- as.numeric(summarise_all(dat, list(mean)))
df_data5 <- tibble(tau, mean_param = get_z(mparam))

ggplot(df_z2, aes(x = tau, y = value)) +
  geom_line(aes(group = id, col = "Uncty"), alpha = 0.1) +
  geom_line(data=df_data1, aes(x = seq(log(t0), log(tN), length.out = N_Nodes), y = MAP_z, col = "MAP")) +
  geom_line(data=df_data3, aes(x = seq(log(t0), log(tN), length.out = N_Nodes), y = mean_resp, linetype = "mean_z")) +
  geom_line(data=df_data4, aes(x = seq(log(t0), log(tN), length.out = N_Nodes), y = mean_zg, linetype = "mean_g")) +
  geom_line(data=df_data5, aes(x = seq(log(t0), log(tN), length.out = N_Nodes), y = mean_param, linetype = "mean_par")) +
  geom_line(data=df_data2, aes(x = tau, y = true_z, col = "TLS")) +
  scale_color_manual(values = c(Uncty = uncty_col, TLS = 'grey30', MAP = MAP_col)) +
  scale_linetype_manual(values = c(mean_par = "dashed", mean_g = 'dotted', mean_z = "dotdash"))
