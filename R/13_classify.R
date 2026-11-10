# 13_classify.R — 추가 분석 ①: 사람 단위로 대륙과 집단 맞히기
# 하는 일:
#   1. 기준 B에서 혼혈 5개 집단을 뺀 21개 집단의 유전자형을 R로 읽는다
#   2. 집단마다 70%는 학습용, 30%는 시험용으로 나눈다
#   3. SNP 수를 바꿔 가며 (학습용으로만 계산한 PCA → 가까운 5명의 다수결) 대륙과 집단 정확도, 계산 시간을 잰다
#   4. 표 10과 그림 9를 저장한다
# 실행 순서: 02_qc.R(기준 B) → 이 스크립트

source("R/00_shared.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("work/cls", showWarnings = FALSE)


# ── 1. 유전자형 읽기 ─────────────────────────────────────────
pops <- read_tsv("work/pops.txt", col_types = cols(.default = "c")) |>
  rename(IID = `#IID`) |>
  filter(!POP %in% admixed_pops)

write_tsv(select(pops, `#IID` = IID), "work/cls/keep.txt")

# --export A : 사람 x SNP 표로 저장한다. 값은 ALT 대립유전자 수(0, 1, 2)
run_plink(c("--pfile", "work/B", "--keep", "work/cls/keep.txt",
            "--export", "A", "--out", "work/cls/geno"))

raw <- data.table::fread("work/cls/geno.raw")
ids <- as.character(raw$IID)
geno <- as.matrix(raw[, -(1:6)])        # 앞 6개 열은 FID, IID 등 사람 정보
rm(raw)

labels <- tibble(IID = ids) |> left_join(pops, by = "IID")
stopifnot(!anyNA(labels$POP), sum(is.na(geno)) == 0)
cat("사람", nrow(geno), "명, SNP", ncol(geno), "개\n")


# ── 2. 학습용 / 시험용 나누기 ────────────────────────────────
# 집단마다 70%를 학습용으로 뽑는다. 모든 조건에서 같은 나눔을 쓴다
set.seed(1)
train_ids <- labels |> group_by(POP) |> slice_sample(prop = 0.7) |> pull(IID)
is_train <- labels$IID %in% train_ids
cat("학습용", sum(is_train), "명, 시험용", sum(!is_train), "명\n")


# ── 3. 분류 한 번 ────────────────────────────────────────────
# PCA는 학습용으로만 계산하고, 시험용은 학습용의 평균과 축으로 옮긴다 (데이터 누수 방지)
# 학습용 안에서 변이가 없는 SNP는 PCA에 쓸 수 없어 뺀다
n_pc <- 10
k <- 5

classify <- function(snp_idx) {
  start <- Sys.time()
  x_train <- geno[is_train, snp_idx, drop = FALSE]
  keep <- apply(x_train, 2, var) > 0
  x_train <- x_train[, keep, drop = FALSE]
  x_test <- geno[!is_train, snp_idx, drop = FALSE][, keep, drop = FALSE]

  pca <- prcomp(x_train, center = TRUE, scale. = FALSE, rank. = n_pc)
  pc_train <- pca$x
  pc_test <- scale(x_test, center = pca$center, scale = FALSE) %*% pca$rotation

  pred_pop <- class::knn(pc_train, pc_test, labels$POP[is_train], k = k)
  pred_super <- class::knn(pc_train, pc_test, labels$SUPER[is_train], k = k)
  seconds <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  list(
    summary = tibble(
      acc_super = mean(pred_super == labels$SUPER[!is_train]),
      acc_pop = mean(pred_pop == labels$POP[!is_train]),
      seconds = seconds
    ),
    predictions = tibble(IID = labels$IID[!is_train], true_pop = labels$POP[!is_train],
                         pred_pop = as.character(pred_pop))
  )
}


# ── 4. SNP 수를 바꿔 가며 반복 ───────────────────────────────
snp_counts <- c(100, 500, 1000, 2000, ncol(geno))
n_rep <- 5

runs <- expand_grid(n_snp = snp_counts, rep = seq_len(n_rep)) |>
  # 전체 SNP를 쓰면 무작위로 고를 것이 없어 결과가 같으므로 한 번만 한다
  filter(n_snp < ncol(geno) | rep == 1) |>
  mutate(result = map2(n_snp, rep, \(n, r) {
    set.seed(100 + r)
    classify(sort(sample(ncol(geno), n)))
  }))

accuracy <- runs |>
  mutate(summary = map(result, "summary")) |>
  select(n_snp, rep, summary) |>
  unnest(summary)

table10 <- accuracy |>
  group_by(n_snp) |>
  summarise(n_runs = n(),
            acc_super_min = min(acc_super), acc_pop_min = min(acc_pop),   # 평균으로 덮어쓰기 전에 최솟값부터
            acc_super = mean(acc_super), acc_pop = mean(acc_pop),
            seconds = mean(seconds), .groups = "drop") |>
  select(n_snp, n_runs, acc_super, acc_super_min, acc_pop, acc_pop_min, seconds)

write_csv(table10, "results/tables/table10_classification.csv")
cat("\n표 10: SNP 수별 정확도 (5번 평균, min은 5번 중 가장 낮은 값)\n")
print(table10, width = Inf)


# ── 5. 전체 SNP일 때 무엇을 헷갈렸나 ─────────────────────────
pred_all <- runs |> filter(n_snp == ncol(geno)) |> pull(result) |> pluck(1, "predictions")

pop_accuracy <- pred_all |>
  left_join(distinct(pops, true_pop = POP, SUPER), by = "true_pop") |>
  group_by(SUPER, true_pop) |>
  summarise(n_test = n(), accuracy = mean(pred_pop == true_pop), .groups = "drop") |>
  arrange(accuracy)

confused <- pred_all |>
  filter(pred_pop != true_pop) |>
  count(true_pop, pred_pop, name = "n_wrong") |>
  arrange(desc(n_wrong))

write_csv(pop_accuracy, "results/tables/table10b_population_accuracy.csv")
write_csv(confused, "results/tables/table10c_confused_pairs.csv")
cat("\n전체 SNP: 집단별 정확도\n")
print(pop_accuracy, n = Inf)
cat("\n전체 SNP: 가장 많이 헷갈린 경우 (실제 → 예측) 상위 10개\n")
print(head(confused, 10))


# ── 6. 그림 9: SNP 수와 정확도 ───────────────────────────────
p9 <- accuracy |>
  pivot_longer(c(acc_super, acc_pop), names_to = "target", values_to = "acc") |>
  mutate(target = recode(target, acc_super = "대륙 (5개 중 하나)", acc_pop = "집단 (21개 중 하나)")) |>
  ggplot(aes(n_snp, acc * 100, colour = target)) +
  geom_point(alpha = 0.5, size = 1.8) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  scale_x_log10(breaks = snp_counts, labels = scales::comma) +
  scale_colour_manual(values = c("대륙 (5개 중 하나)" = "#2a78d6", "집단 (21개 중 하나)" = "#eb6834"), name = NULL) +
  labs(x = "사용한 SNP 수 (로그 눈금)", y = "시험용 정확도 (%)",
       title = "SNP 수에 따른 분류 정확도",
       subtitle = paste0("혼혈 제외 21개 집단, 학습 ", sum(is_train), "명 / 시험 ", sum(!is_train),
                         "명, PC ", n_pc, "개 + 가까운 ", k, "명 다수결\n점: 무작위 SNP 선택 5번, 선: 평균")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave("results/figures/fig9_classification_accuracy.png", p9, width = 8, height = 5, dpi = 200, bg = "white")


# ── 7. 그림 9b: 전체 SNP일 때 집단 혼동 행렬 ──────────────────
pop_order <- pops |> distinct(POP, SUPER) |>
  mutate(SUPER = factor(SUPER, levels = super_levels)) |> arrange(SUPER, POP) |> pull(POP)

p9b <- pred_all |>
  count(true_pop, pred_pop) |>
  group_by(true_pop) |> mutate(share = n / sum(n)) |> ungroup() |>
  mutate(true_pop = factor(true_pop, levels = rev(pop_order)),
         pred_pop = factor(pred_pop, levels = pop_order)) |>
  ggplot(aes(pred_pop, true_pop, fill = share)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "#eef4fc", high = "#123a6b", name = "비율", labels = scales::percent) +
  coord_fixed() +
  labs(x = "예측한 집단", y = "실제 집단",
       title = "집단 혼동 행렬 (SNP 8,551개)",
       subtitle = "줄마다 실제 집단의 시험용 사람이 어느 집단으로 예측됐는지의 비율") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid = element_blank())

ggsave("results/figures/fig9b_confusion.png", p9b, width = 7.5, height = 7, dpi = 200, bg = "white")
