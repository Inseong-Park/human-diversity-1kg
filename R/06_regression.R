# 06_regression.R — 아디스아바바로부터의 거리로 He를 설명하는 회귀
# 하는 일:
#   1. He(표 3)와 거리(표 4)를 합친다
#   2. 포함 집단을 바꿔 가며 직선을 맞추고 표 5에 저장한다
#   3. 산점도(그림 2)와 잔차 막대그래프(그림 2-1)를 저장한다
# 실행 순서: 04_diversity.R → 05_distance.R → 이 스크립트

source("R/00_shared.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)


# ── 1. 다양성과 거리 합치기 ──────────────────────────────────
diversity <- read_csv("results/tables/table3_diversity.csv", show_col_types = FALSE) |>
  select(POP, SUPER, admixed, He, Ho_mean)

distance <- read_csv("results/tables/table4_coordinates_distance.csv", show_col_types = FALSE) |>
  select(POP, distance_km)

d <- inner_join(diversity, distance, by = "POP")
stopifnot(nrow(d) == 26)


# ── 2. 표 5: 회귀 버전별 결과 ────────────────────────────────
# R2_africa, R2_outside: Day 6에서 아프리카 안의 거리와 He가 따로 노는 것을 보고 추가했다
r2_data <- filter(d, !admixed)

versions <- list(
  R1         = list(data = d,                                     y = "He",      label = "26개 전부"),
  R2         = list(data = r2_data,                               y = "He",      label = "혼혈 5개 제외, PEL 포함"),
  R3         = list(data = filter(r2_data, POP != "PEL"),         y = "He",      label = "혼혈 5개와 PEL 제외"),
  R4         = list(data = r2_data,                               y = "Ho_mean", label = "R2를 관측 이형접합도로 반복"),
  R2_africa  = list(data = filter(r2_data, SUPER == "AFR"),       y = "He",      label = "R2 중 아프리카 안"),
  R2_outside = list(data = filter(r2_data, SUPER != "AFR"),       y = "He",      label = "R2 중 아프리카 밖")
)

table5 <- imap(versions, \(v, name) {
  fit_distance(v$data, v$y) |>
    mutate(version = name, included = v$label, response = v$y, .before = 1)
}) |>
  list_rbind()

write_csv(table5, "results/tables/table5_regression.csv")
print(table5, width = Inf)


# ── 3. R2 직선 기준 잔차 ─────────────────────────────────────
# 잔차 = 실제 He - 직선이 예상한 He. 혼혈 5개는 직선을 맞출 때 쓰지 않았지만 같은 직선으로 잔차를 계산한다
r2_fit <- filter(table5, version == "R2")

residuals_r2 <- d |>
  mutate(
    He_pred = r2_fit$intercept + r2_fit$slope_per_1000km * distance_km / 1000,
    residual = He - He_pred,
    SUPER = factor(SUPER, levels = super_levels),
    group_type = if_else(admixed, "최근 혼혈 집단", "그 외 집단")
  ) |>
  arrange(distance_km)

write_csv(residuals_r2, "results/tables/table5b_residuals.csv")
print(select(residuals_r2, POP, SUPER, distance_km, He, He_pred, residual), n = Inf)

# Day 6 예측 3: CEU를 유타 좌표(05_distance.R 출력값 16,748km)로 넣었을 때 직선이 예상하는 He
ceu_he <- d$He[d$POP == "CEU"]
ceu_utah_pred <- r2_fit$intercept + r2_fit$slope_per_1000km * 16748 / 1000
cat("\nCEU 실제 He:", round(ceu_he, 4), " / 유타 좌표일 때 직선이 예상하는 He:", round(ceu_utah_pred, 4), "\n")


# ── 4. 그림 2: 거리와 He ─────────────────────────────────────
# 글자가 겹치지 않도록 대륙 안에서 거리순으로 위아래를 번갈아 놓는다
# 점이 붙어 있는 집단은 글자 위치를 직접 정한다
label_offsets <- tribble(
  ~POP,  ~dx,   ~dy,
  "ASW",    0,  0.0035,
  "ACB",  520,  0,
  "YRI",  500,  0,
  "ESN", -500,  0,
  "PUR", -560,  0,
  "CLM",  560,  0,
  "CHB",    0,  0.0035,
  "CDX",    0, -0.0035,
  "CHS",    0, -0.0035,
  "KHV",    0,  0.0035,
  "JPT",    0,  0.0035,
  "CEU", -500,  0,
  "FIN",    0, -0.0035,
  "GBR",  500,  0
)

plot_data <- residuals_r2 |>
  group_by(SUPER) |>
  mutate(default_dy = if_else(row_number() %% 2 == 1, 0.0035, -0.0035)) |>
  ungroup() |>
  left_join(label_offsets, by = "POP") |>
  mutate(label_x = distance_km + coalesce(dx, 0),
         label_y = He + coalesce(dy, default_dy))

p2 <- ggplot(plot_data, aes(distance_km, He)) +
  geom_abline(intercept = r2_fit$intercept, slope = r2_fit$slope_per_1000km / 1000,
              colour = "grey45", linewidth = 0.6) +
  geom_point(aes(colour = SUPER, shape = group_type), size = 2.8, stroke = 1) +
  geom_text(aes(x = label_x, y = label_y, label = POP), size = 2.7, colour = "grey20") +
  expand_limits(y = 0.222) +
  scale_colour_manual(values = super_colours, name = "대륙") +
  scale_shape_manual(values = c("그 외 집단" = 16, "최근 혼혈 집단" = 1), name = NULL) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "아디스아바바로부터 경유지를 거친 거리 (km)", y = "기대 이형접합도 (He)",
       title = "거리에 따른 유전적 다양성",
       subtitle = "회색 선: 혼혈 5개 집단을 뺀 21개 집단의 직선 (R2)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave("results/figures/fig2_distance_he.png", p2, width = 9, height = 5.5, dpi = 200, bg = "white")


# ── 5. 그림 2-1: 잔차 ────────────────────────────────────────
# 왼쪽일수록 가까운 집단이다. 막대가 위면 직선보다 He가 높다
p2b <- residuals_r2 |>
  mutate(POP = fct_inorder(POP)) |>
  ggplot(aes(POP, residual, fill = SUPER, alpha = group_type)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = super_colours, name = "대륙") +
  scale_alpha_manual(values = c("그 외 집단" = 0.95, "최근 혼혈 집단" = 0.4), name = NULL) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey30"))) +
  labs(x = "집단 (왼쪽일수록 가까움)", y = "잔차 (실제 He − 직선 예상 He)",
       title = "R2 직선에서 벗어난 정도") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave("results/figures/fig2b_residuals.png", p2b, width = 9, height = 4.5, dpi = 200, bg = "white")
