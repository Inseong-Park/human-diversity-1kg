# 10_compare.R — 전처리 조건 비교 (B vs V1~V5)와 민감도 2×2
# 하는 일:
#   1. 6개 데이터(B, V1~V5)마다 집단별 He, 개인별 Ho, PCA를 계산한다
#   2. 거리-He 회귀(R2: 혼혈 5개 제외 21개 집단)를 데이터마다 다시 맞춘다
#   3. 표 6(전처리 비교), 표 7(민감도 2×2)과 비교 그림 3종을 저장한다
# 실행 순서: 02_qc.R(B, V1~V4), 03_kinship.R(친척 쌍), 04b_excess_het.R(V5), 05_distance.R(거리) → 이 스크립트

source("R/00_shared.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("work/cmp", showWarnings = FALSE)

datasets <- c("B", "V1", "V2", "V3", "V4", "V5")
dataset_info <- tibble(
  dataset = datasets,
  condition = c("기준", "친척 포함", "MAF 0.01", "pruning 생략", "MAF 0.01 + pruning 생략", "이상 SNP 제외")
)

pops <- read_tsv("work/pops.txt", col_types = cols(.default = "c")) |>
  rename(IID = `#IID`)

distance <- read_csv("results/tables/table4_coordinates_distance.csv", show_col_types = FALSE) |>
  select(POP, distance_km)


# ── 1. 집단별 He ─────────────────────────────────────────────
# --loop-cats POP : 집단마다 따로 빈도를 계산한다 (Day 5의 --keep 반복과 같은 결과)
calc_he <- function(dataset) {
  out <- paste0("work/cmp/freq_", dataset)
  run_plink(c("--pfile", paste0("work/", dataset), "--pheno", "work/pops.txt",
              "--loop-cats", "POP", "--freq", "--out", out))
  map(unique(pops$POP), \(p) {
    freq <- read_tsv(paste0(out, ".", p, ".afreq"),
                     col_types = cols(.default = "c", ALT_FREQS = "d"))
    tibble(dataset = dataset, POP = p, n_snp = nrow(freq),
           He = mean(2 * freq$ALT_FREQS * (1 - freq$ALT_FREQS)))
  }) |> list_rbind()
}

# 개인별 Ho의 집단 평균
calc_ho <- function(dataset) {
  out <- paste0("work/cmp/het_", dataset)
  run_plink(c("--pfile", paste0("work/", dataset), "--het", "--out", out))
  read_tsv(paste0(out, ".het"), col_types = cols(`#IID` = "c", .default = "d")) |>
    transmute(IID = `#IID`, Ho = (OBS_CT - `O(HOM)`) / OBS_CT) |>
    inner_join(pops, by = "IID") |>
    group_by(POP) |>
    summarise(n = n(), Ho_mean = mean(Ho), .groups = "drop") |>
    mutate(dataset = dataset)
}

diversity <- map(datasets, calc_he) |> list_rbind() |>
  left_join(map(datasets, calc_ho) |> list_rbind(), by = c("dataset", "POP")) |>
  left_join(distinct(pops, POP, SUPER), by = "POP") |>
  left_join(distance, by = "POP") |>
  mutate(admixed = POP %in% admixed_pops,
         SUPER = factor(SUPER, levels = super_levels))

stopifnot(nrow(diversity) == 26 * length(datasets))


# ── 2. 회귀 (R2 버전) ────────────────────────────────────────
regression <- diversity |>
  filter(!admixed) |>
  group_by(dataset) |>
  group_modify(\(d, key) fit_distance(d, "He")) |>
  ungroup()


# ── 3. PCA ───────────────────────────────────────────────────
# V3는 로딩 비교에 쓰므로 allele-wts도 저장한다
run_pca <- function(dataset) {
  out <- paste0("work/cmp/pca_", dataset)
  extra <- if (dataset %in% c("B", "V3")) "allele-wts" else character()
  run_plink(c("--pfile", paste0("work/", dataset), "--pca", extra, "10", "--out", out))
  read_tsv(paste0(out, ".eigenvec"), show_col_types = FALSE) |>
    rename(IID = `#IID`) |>
    mutate(IID = as.character(IID), dataset = dataset)
}

pca <- map(datasets, run_pca) |> list_rbind()

# B와 같은 사람끼리 PC 점수의 상관. PC의 부호는 임의이므로 절댓값을 쓴다
# 1에 가까우면 그 PC가 B와 같은 축이다
pca_b <- filter(pca, dataset == "B")
pca_similarity <- pca |>
  filter(dataset != "B") |>
  group_by(dataset) |>
  group_modify(\(d, key) {
    m <- inner_join(d, pca_b, by = "IID", suffix = c("", "_B"))
    tibble(n_common = nrow(m),
           PC1_abs_r = abs(cor(m$PC1, m$PC1_B)),
           PC2_abs_r = abs(cor(m$PC2, m$PC2_B)),
           PC3_abs_r = abs(cor(m$PC3, m$PC3_B)))
  }) |>
  ungroup()


# ── 4. 표 6: 전처리 비교 ─────────────────────────────────────
# 아프리카 격차: 비혼혈 아프리카 5개 집단 평균 He − 아프리카 밖 비혼혈 16개 집단(PEL 포함) 평균 He
he_summary <- diversity |>
  filter(!admixed) |>
  group_by(dataset) |>
  summarise(
    n_snp = first(n_snp),
    He_AFR = mean(He[SUPER == "AFR"]),
    He_outside = mean(He[SUPER != "AFR"]),
    afr_gap = He_AFR - He_outside,
    outside_over_AFR = He_outside / He_AFR,
    Ho_minus_He = mean(Ho_mean - He),
    .groups = "drop"
  )

# 26개 집단 He 순서가 B와 얼마나 같은가 (스피어만 순위 상관)
he_rank <- diversity |>
  select(dataset, POP, He) |>
  left_join(filter(diversity, dataset == "B") |> select(POP, He_B = He), by = "POP") |>
  group_by(dataset) |>
  summarise(He_rank_r_with_B = cor(He, He_B, method = "spearman"), .groups = "drop")

table6 <- dataset_info |>
  left_join(he_summary, by = "dataset") |>
  left_join(he_rank, by = "dataset") |>
  left_join(select(regression, dataset, slope_per_1000km, r2, p), by = "dataset") |>
  left_join(pca_similarity, by = "dataset")

write_csv(table6, "results/tables/table6_preprocessing_compare.csv")
write_csv(diversity |> arrange(dataset, SUPER, desc(He)),
          "results/tables/table6b_he_by_dataset.csv")
cat("\n표 6: 전처리 비교\n")
print(table6, width = Inf)


# ── 5. 표 7: 민감도 2×2 ──────────────────────────────────────
table7 <- tibble(dataset = c("B", "V2", "V3", "V4"),
                 MAF = c("0.05", "0.01", "0.05", "0.01"),
                 pruning = c("함", "함", "생략", "생략")) |>
  left_join(he_summary |> select(dataset, n_snp), by = "dataset") |>
  left_join(regression, by = "dataset") |>
  mutate(slope_sign = if_else(slope_per_1000km < 0, "음수(멀수록 감소)", "양수"))

write_csv(table7, "results/tables/table7_sensitivity_2x2.csv")
cat("\n표 7: 민감도 2×2\n")
print(table7, width = Inf)


# ── 6. 비교 1: 친척 쌍이 PCA에서 가깝게 보이는가 ─────────────
# Day 4에서 찾은 친족계수 0.177 이상 15쌍. 이 쌍들은 V1에만 둘 다 들어 있다
pairs <- read_tsv("work/king.kin0", col_types = cols(`#IID1` = "c", IID2 = "c", .default = "d")) |>
  rename(IID1 = `#IID1`) |>
  filter(KINSHIP >= 0.177) |>
  left_join(select(pops, IID1 = IID, POP1 = POP), by = "IID1") |>
  left_join(select(pops, IID2 = IID, POP2 = POP), by = "IID2") |>
  # 두 사람의 집단이 다르면(STU와 ITU) 두 집단을 합쳐 한 묶음으로 본다
  mutate(group = if_else(POP1 == POP2, POP1, paste(POP1, POP2, sep = "+")))

v1 <- filter(pca, dataset == "V1") |> inner_join(pops, by = "IID")
pc_cols <- paste0("PC", 1:10)

# 각 쌍의 PC 거리가 같은 집단 안의 모든 두 사람 거리 중 몇 % 위치인가 (작을수록 유난히 가까움)
pair_percentile <- function(id1, id2, group_pops, cols) {
  members <- filter(v1, POP %in% group_pops)
  x <- as.matrix(members[, cols])
  all_d <- as.vector(dist(x))
  d <- sqrt(sum((x[members$IID == id1, ] - x[members$IID == id2, ])^2))
  100 * mean(all_d <= d)
}

relatives_check <- pairs |>
  rowwise() |>
  mutate(pct_PC1_2 = pair_percentile(IID1, IID2, c(POP1, POP2), c("PC1", "PC2")),
         pct_PC1_10 = pair_percentile(IID1, IID2, c(POP1, POP2), pc_cols)) |>
  ungroup() |>
  select(IID1, IID2, POP1, POP2, KINSHIP, pct_PC1_2, pct_PC1_10) |>
  arrange(desc(KINSHIP))

cat("\n비교 1: 친척 쌍의 PC 거리 (같은 집단 안 모든 두 사람 거리 중 하위 몇 %)\n")
print(relatives_check, n = Inf)

# PC마다 가장 끝에 있는 두 사람이 친척인지 (친척 쌍이 PC 하나를 차지했는지)
related_ids <- unique(c(pairs$IID1, pairs$IID2))
cat("\nV1의 PC별 가장 끝에 있는 2명 중 친척 쌍에 속한 사람 수:\n")
map(pc_cols, \(pc) {
  top <- v1 |> slice_max(abs(.data[[pc]]), n = 2)
  tibble(PC = pc, top2 = paste(top$IID, top$POP, collapse = ", "),
         n_related = sum(top$IID %in% related_ids))
}) |> list_rbind() |> print()


# ── 7. 그림 6a: V1 PCA 위의 친척 쌍 ─────────────────────────
# 친척 쌍이 있는 집단만 확대해서 그린다 (집단이 다른 쌍은 두 집단을 한 칸에)
group_members <- pairs |>
  distinct(group, POP1, POP2) |>
  pivot_longer(c(POP1, POP2), values_to = "POP") |>
  distinct(group, POP) |>
  inner_join(v1, by = "POP", relationship = "many-to-many")

segments <- pairs |>
  left_join(select(v1, IID1 = IID, x = PC1, y = PC2), by = "IID1") |>
  left_join(select(v1, IID2 = IID, xend = PC1, yend = PC2), by = "IID2")

p6a <- ggplot(group_members, aes(PC1, PC2)) +
  geom_point(colour = "grey70", size = 0.9) +
  geom_segment(data = segments, aes(x = x, y = y, xend = xend, yend = yend),
               colour = "#eb6834", linewidth = 0.5) +
  geom_point(data = segments, aes(x, y), colour = "#eb6834", size = 1.4) +
  geom_point(data = segments, aes(xend, yend), colour = "#eb6834", size = 1.4) +
  facet_wrap(~group, scales = "free") +
  labs(title = "비교 1: 친척 쌍의 PC1-PC2 위치 (V1, 친척 포함)",
       subtitle = "회색: 같은 집단의 다른 사람, 주황 선: 친족계수 0.177 이상 15쌍. 집단마다 축 범위가 다름") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        axis.text = element_text(size = 7))

ggsave("results/figures/fig6a_relatives_pca.png", p6a, width = 9, height = 6, dpi = 200, bg = "white")


# ── 8. 그림 6b: 비교별 집단 He ───────────────────────────────
pop_order <- diversity |> filter(dataset == "B") |> arrange(SUPER, desc(He)) |> pull(POP)
comparisons <- tibble(dataset = c("V1", "V2", "V3", "V5"),
                      panel = c("비교 1: 친척 포함 (V1)", "비교 2: MAF 0.01 (V2)",
                                "비교 3: pruning 생략 (V3)", "비교 4: 이상 SNP 제외 (V5)"))

he_plot <- comparisons |>
  rowwise() |>
  reframe(panel = panel,
          filter(diversity, dataset %in% c("B", .env$dataset)) |>
            mutate(version = if_else(dataset == "B", "B (기준)", "바꾼 조건"))) |>
  mutate(POP = factor(POP, levels = pop_order))

p6b <- ggplot(he_plot, aes(POP, He, colour = version, group = version)) +
  geom_line(linewidth = 0.4, alpha = 0.6) +
  geom_point(size = 1.3) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("B (기준)" = "grey45", "바꾼 조건" = "#2a78d6"), name = NULL) +
  labs(x = "집단 (대륙 순서: AFR, EUR, SAS, EAS, AMR / 대륙 안은 B의 He 순)", y = "He",
       title = "조건을 바꿨을 때 집단별 He") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
        legend.position = "top")

ggsave("results/figures/fig6b_he_compare.png", p6b, width = 9, height = 9, dpi = 200, bg = "white")


# ── 9. 그림 6c: 로딩 B vs V3 ─────────────────────────────────
read_loadings <- function(dataset) {
  read_tsv(paste0("work/cmp/pca_", dataset, ".eigenvec.allele"),
           col_types = cols(.default = "c", PC1 = "d", PC2 = "d")) |>
    filter(A1 == ALT) |>
    transmute(dataset = dataset, pos = as.numeric(sub(".*:", "", ID)), PC1, PC2)
}

loadings <- bind_rows(read_loadings("B"), read_loadings("V3"))
stopifnot(sum(loadings$dataset == "B") == 8551, sum(loadings$dataset == "V3") == 97088)

# 로딩 크기가 상위 1%인 SNP가 염색체 1Mb 구간 중 몇 곳에 몰려 있는가
cat("\n로딩 절댓값 상위 1% SNP가 가장 많이 몰린 1Mb 구간 (상위 3개)\n")
loadings |>
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "weight") |>
  group_by(dataset, PC) |>
  filter(abs(weight) >= quantile(abs(weight), 0.99)) |>
  count(bin_Mb = floor(pos / 1e6)) |>
  mutate(share = n / sum(n)) |>
  slice_max(n, n = 3, with_ties = FALSE) |>
  print(n = Inf)

p6c <- loadings |>
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "weight") |>
  mutate(dataset = recode(dataset, B = "B (pruning 함, 8,551개)", V3 = "V3 (pruning 생략, 97,088개)")) |>
  ggplot(aes(pos / 1e6, weight)) +
  geom_point(size = 0.3, alpha = 0.4, colour = "#2a78d6") +
  facet_grid(PC ~ dataset, scales = "free_y") +
  labs(x = "22번 염색체 위치 (Mb)", y = "로딩",
       title = "비교 3: PC1, PC2 로딩의 위치 분포",
       subtitle = "PC의 부호는 임의라서 B와 V3의 위아래 방향은 비교하지 않는다") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))

ggsave("results/figures/fig6c_loadings_B_V3.png", p6c, width = 10, height = 6, dpi = 200, bg = "white")
