# 09_wahlund.R — 월룬드 효과: 전체 기준 F와 집단 기준 F 비교
# 하는 일:
#   1. 전체 기준 F: 2,491명 전체의 빈도로 계산한 근교계수
#   2. 집단 기준 F: 집단마다 자기 집단의 빈도로 계산한 근교계수 (B와 V5)
#   3. 집단별 표와 boxplot(그림 5)을 저장하고, 집단 기준 F가 높은 사람을 출력한다
# 실행 순서: 04b_excess_het.R(V5) → 이 스크립트

source("R/00_shared.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

pops <- read_tsv("work/pops.txt", col_types = cols(.default = "c")) |>
  rename(IID = `#IID`)
all_pops <- sort(unique(pops$POP))

# PLINK --het 결과 읽기. F = (관측 동형접합 - 기대 동형접합) / (전체 - 기대 동형접합)
read_het <- function(path) {
  read_tsv(path, col_types = cols(`#IID` = "c", .default = "d")) |>
    transmute(IID = `#IID`, F)
}


# ── 1. 전체 기준 F ───────────────────────────────────────────
# 빈도를 2,491명 전체에서 계산한다 (Day 5의 work/B_het.het와 같은 계산)
run_plink(c("--pfile", "work/B", "--het", "--out", "work/B_het"))
run_plink(c("--pfile", "work/V5", "--het", "--out", "work/V5_het"))

f_overall <- read_het("work/B_het.het") |> rename(F_overall = F)


# ── 2. 집단 기준 F ───────────────────────────────────────────
# --keep으로 한 집단만 남기면 PLINK가 그 집단 안에서 빈도를 다시 계산한다
het_within <- function(dataset, pop) {
  out <- paste0("work/het_within_", dataset, "_", pop)
  run_plink(c("--pfile", paste0("work/", dataset),
              "--keep", paste0("work/keep_", pop, ".txt"),
              "--het", "--out", out))
  read_het(paste0(out, ".het"))
}

f_within_B  <- map(all_pops, \(p) het_within("B", p))  |> list_rbind() |> rename(F_within = F)
f_within_V5 <- map(all_pops, \(p) het_within("V5", p)) |> list_rbind() |> rename(F_within_V5 = F)

f_all <- pops |>
  inner_join(f_overall, by = "IID") |>
  inner_join(f_within_B, by = "IID") |>
  inner_join(f_within_V5, by = "IID") |>
  mutate(SUPER = factor(SUPER, levels = super_levels),
         admixed = POP %in% admixed_pops)

stopifnot(nrow(f_all) == 2491)


# ── 3. 표: 집단별 평균 F ─────────────────────────────────────
f_table <- f_all |>
  group_by(SUPER, POP, admixed) |>
  summarise(
    n = n(),
    F_overall_mean   = mean(F_overall),
    F_within_mean    = mean(F_within),
    F_within_V5_mean = mean(F_within_V5),
    wahlund_gap      = F_overall_mean - F_within_mean,   # 집단을 섞어서 생긴 몫
    n_F_within_over_0.0625 = sum(F_within > 0.0625),     # 사촌 간 결혼의 자식 수준(1/16) 이상
    .groups = "drop"
  ) |>
  arrange(SUPER, desc(F_overall_mean))

write_csv(f_table, "results/tables/table_wahlund_f.csv")
print(f_table, n = Inf, width = Inf)

cat("\n대륙별 평균 (혼혈 5개 제외):\n")
f_table |>
  filter(!admixed) |>
  group_by(SUPER) |>
  summarise(F_overall = mean(F_overall_mean), F_within = mean(F_within_mean),
            F_within_V5 = mean(F_within_V5_mean), .groups = "drop") |>
  print()

cat("\n집단 기준 F가 가장 높은 10명:\n")
f_all |> slice_max(F_within, n = 10) |> select(IID, POP, SUPER, F_overall, F_within, F_within_V5) |> print()


# ── 4. 그림 5: 두 F를 나란히 ─────────────────────────────────
plot_data <- f_all |>
  select(IID, SUPER, POP, F_overall, F_within) |>
  pivot_longer(c(F_overall, F_within), names_to = "basis", values_to = "F") |>
  mutate(
    basis = recode(basis, F_overall = "전체 기준", F_within = "집단 기준"),
    POP = fct_reorder(POP, F, .fun = median, .desc = TRUE)
  )

p <- ggplot(plot_data, aes(POP, F, fill = basis)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.4, outlier.alpha = 0.5, width = 0.75,
               position = position_dodge(width = 0.8)) +
  facet_grid(~SUPER, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("전체 기준" = "#2a78d6", "집단 기준" = "#eb6834"), name = "빈도 기준") +
  labs(x = NULL, y = "근교계수 F",
       title = "전체 기준 F와 집단 기준 F",
       subtitle = "기준 B: 2,491명, 22번 염색체 SNP 8,551개. 같은 사람을 빈도 기준만 바꿔 계산") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.spacing.x = unit(4, "pt"), strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave("results/figures/fig5_wahlund.png", p, width = 10, height = 5.5, dpi = 200, bg = "white")
