# 12_loo_he.R — 추가 분석 ②: 보지 않은 집단의 He 예측
# 하는 일:
#   1. 21개 집단(혼혈 5개 제외, PEL 포함)에서 집단 하나를 빼고 직선을 맞춰 빠진 집단의 He를 예측한다 (21번 반복)
#   2. 모델 A(거리만)와 모델 B(거리 + 대륙)의 예측 오차를 비교한다
#   3. 표 9와 그림 8을 저장한다
# 실행 순서: 04_diversity.R(He), 05_distance.R(거리) → 이 스크립트

source("R/00_shared.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

pops <- read_csv("results/tables/table3_diversity.csv", show_col_types = FALSE) |>
  filter(!admixed) |>
  select(SUPER, POP, He) |>
  left_join(read_csv("results/tables/table4_coordinates_distance.csv", show_col_types = FALSE) |>
              select(POP, distance_km), by = "POP") |>
  mutate(dist_1000 = distance_km / 1000,
         SUPER = factor(SUPER, levels = super_levels))

stopifnot(nrow(pops) == 21)


# ── 1. 전체 21개로 맞춘 두 모델 (기울기 비교용) ──────────────
fit_a <- lm(He ~ dist_1000, data = pops)
fit_b <- lm(He ~ dist_1000 + SUPER, data = droplevels(filter(pops, POP != "PEL")))
# 모델 B는 PEL 하나뿐인 아메리카 높이가 PEL에 딱 맞춰지므로, 기울기는 PEL을 뺀 20개로 맞춘다

slopes <- bind_rows(
  tibble(model = "A: 거리만 (21개)", term = "dist_1000",
         estimate = coef(fit_a)[["dist_1000"]],
         ci_low = confint(fit_a)["dist_1000", 1], ci_high = confint(fit_a)["dist_1000", 2],
         r2 = summary(fit_a)$r.squared),
  tibble(model = "B: 거리 + 대륙 (PEL 제외 20개)", term = "dist_1000",
         estimate = coef(fit_b)[["dist_1000"]],
         ci_low = confint(fit_b)["dist_1000", 1], ci_high = confint(fit_b)["dist_1000", 2],
         r2 = summary(fit_b)$r.squared)
)
cat("\n거리 기울기 (1,000km당 He 변화, 95% 신뢰구간)\n")
print(slopes, width = Inf)


# ── 2. 하나씩 빼고 예측 ──────────────────────────────────────
predict_left_out <- function(pop) {
  train <- filter(pops, POP != pop)
  test  <- filter(pops, POP == pop)
  pred_a <- predict(lm(He ~ dist_1000, data = train), test)
  # 빠진 집단의 대륙이 학습 데이터에 없으면(PEL) 모델 B로 예측할 수 없다
  pred_b <- if (test$SUPER %in% train$SUPER) {
    predict(lm(He ~ dist_1000 + SUPER, data = droplevels(train)), test)
  } else NA_real_
  tibble(POP = pop, pred_A = pred_a, pred_B = pred_b)
}

loo <- map(pops$POP, predict_left_out) |>
  list_rbind() |>
  left_join(pops, by = "POP") |>
  mutate(error_A = pred_A - He, error_B = pred_B - He) |>
  select(SUPER, POP, distance_km, He, pred_A, error_A, pred_B, error_B) |>
  arrange(SUPER, distance_km)

write_csv(loo, "results/tables/table9_loo_prediction.csv")
cat("\n표 9: 집단별 예측 (error = 예측 - 실제)\n")
print(loo, n = Inf, width = Inf)


# ── 3. 요약: 평균 오차 ───────────────────────────────────────
summary_table <- bind_rows(
  loo |> summarise(set = "21개 (모델 A만)", MAE_A = mean(abs(error_A)), MAE_B = NA_real_),
  loo |> filter(POP != "PEL") |>
    summarise(set = "PEL 제외 20개 (두 모델 비교)", MAE_A = mean(abs(error_A)), MAE_B = mean(abs(error_B)))
)
cat("\n평균 오차 (절댓값)\n")
print(summary_table)

cat("\n모델 A에서 오차가 가장 큰 집단 3개\n")
loo |> slice_max(abs(error_A), n = 3) |> select(POP, He, pred_A, error_A) |> print()

cat("\n대륙별 평균 오차 (절댓값)\n")
loo |> group_by(SUPER) |>
  summarise(n = n(), MAE_A = mean(abs(error_A)), MAE_B = mean(abs(error_B)), .groups = "drop") |>
  print()


# ── 4. 그림 8: 예측 vs 실제 ──────────────────────────────────
plot_data <- loo |>
  pivot_longer(c(pred_A, pred_B), names_to = "model", values_to = "pred") |>
  filter(!is.na(pred)) |>
  mutate(model = recode(model, pred_A = "모델 A: 거리만", pred_B = "모델 B: 거리 + 대륙"))

lims <- range(c(plot_data$He, plot_data$pred))

p <- ggplot(plot_data, aes(He, pred, colour = SUPER)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed") +
  geom_point(size = 2.2) +
  # 같은 대륙 점들이 겹쳐 이름이 읽히지 않으므로, 모델 A에서 오차가 가장 큰 2개 집단에만 이름을 붙인다
  geom_text(data = plot_data |> filter(model == "모델 A: 거리만") |> slice_max(abs(pred - He), n = 2),
            aes(label = POP), size = 3, hjust = 1.3, show.legend = FALSE) +
  facet_wrap(~model) +
  scale_colour_manual(values = super_colours, name = NULL) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(x = "실제 He", y = "예측 He (그 집단을 빼고 맞춘 직선)",
       title = "보지 않은 집단의 He 예측",
       subtitle = "점선 위면 정확히 맞힌 것. 이름은 모델 A에서 오차가 가장 큰 2개 집단. 모델 B는 PEL을 예측할 수 없어 20개") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave("results/figures/fig8_loo_prediction.png", p, width = 9, height = 5.5, dpi = 200, bg = "white")
