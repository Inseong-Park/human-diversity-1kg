# 04_diversity.R — 기준 데이터 B로 집단별 유전적 다양성 계산
# 하는 일:
#   1. 개인별 관측 이형접합도(Ho)를 계산한다
#   2. 집단별 대립유전자 빈도로 기대 이형접합도(He)를 계산한다
#   3. 집단별 표(표 3)와 boxplot(그림 1)을 저장한다

source("R/00_shared.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

pops <- read_tsv("work/pops.txt", col_types = cols(.default = "c")) |>
  rename(IID = `#IID`)


# ── 1. 개인별 관측 이형접합도 ─────────────────────────────────
# Ho = (관측된 유전자형 수 - 동형접합 수) / 관측된 유전자형 수
# 사람 한 명의 유전자형만으로 정해지는 값이다
# 같은 파일의 F는 Day 10(월룬드 효과)에서 쓴다
run_plink(c("--pfile", "work/B", "--het", "--out", "work/B_het"))

het <- read_tsv("work/B_het.het", show_col_types = FALSE) |>
  rename(IID = `#IID`, O_HOM = `O(HOM)`) |>
  mutate(IID = as.character(IID),
         Ho = (OBS_CT - O_HOM) / OBS_CT) |>
  inner_join(pops, by = "IID")

stopifnot(nrow(het) == 2491)   # 기준 B의 인원과 같아야 한다


# ── 2. 집단별 기대 이형접합도 ─────────────────────────────────
# 집단마다 대립유전자 빈도 p를 구하고, SNP마다 2pq를 계산해 평균을 낸다
# keep 파일에 친척 13명이 들어 있어도, work/B에 없는 사람이므로 계산에 쓰이지 않는다
calc_he <- function(pop) {
  out <- paste0("work/freq_B_", pop)
  run_plink(c("--pfile", "work/B",
              "--keep", paste0("work/keep_", pop, ".txt"),
              "--freq", "--out", out))

  # ID("22:16052080")를 R이 시각으로 잘못 읽지 않도록 모든 열을 글자로 읽고, 빈도만 숫자로 읽는다
  freq <- read_tsv(paste0(out, ".afreq"), col_types = cols(.default = "c", ALT_FREQS = "d"))
  stopifnot(nrow(freq) == 8551)   # 모든 집단이 같은 SNP 목록을 써야 비교할 수 있다

  tibble(POP = pop, He = mean(2 * freq$ALT_FREQS * (1 - freq$ALT_FREQS)))
}

he <- map(unique(het$POP), calc_he) |> list_rbind()


# ── 3. 표 3: 집단별 다양성 ───────────────────────────────────
table3 <- summarize_groups(het, Ho, c(SUPER, POP)) |>
  rename(Ho_mean = mean, Ho_sd = sd) |>
  left_join(he, by = "POP") |>
  mutate(
    admixed = POP %in% admixed_pops,
    SUPER = factor(SUPER, levels = super_levels)
  ) |>
  arrange(SUPER, desc(He)) |>
  select(SUPER, POP, admixed, n, He, Ho_mean, Ho_sd)

write_csv(table3, "results/tables/table3_diversity.csv")
print(table3, n = Inf)


# ── 4. 그림 1: 집단별 관측 이형접합도 ─────────────────────────
# 대륙 칸을 나누고, 칸 안에서는 Ho 중앙값이 높은 집단부터 놓는다
# 혼혈 집단은 색을 옅게 칠해 색을 구분하기 어려워도 알아볼 수 있게 한다
plot_data <- het |>
  mutate(
    SUPER = factor(SUPER, levels = super_levels),
    POP = fct_reorder(POP, Ho, .fun = median, .desc = TRUE),
    group_type = if_else(POP %in% admixed_pops, "최근 혼혈 집단", "그 외 집단")
  )

p <- ggplot(plot_data, aes(POP, Ho, fill = SUPER, alpha = group_type)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.5, outlier.alpha = 0.6, width = 0.7) +
  facet_grid(~SUPER, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = super_colours, name = "대륙") +
  scale_alpha_manual(values = c("그 외 집단" = 0.95, "최근 혼혈 집단" = 0.35), name = NULL) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey30"))) +
  labs(x = NULL, y = "관측 이형접합도 (Ho)",
       title = "집단별 개인 관측 이형접합도",
       subtitle = "기준 B: 2,491명, 22번 염색체 SNP 8,551개") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(4, "pt"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("results/figures/fig1_heterozygosity_box.png", p,
       width = 9, height = 5, dpi = 200, bg = "white")
