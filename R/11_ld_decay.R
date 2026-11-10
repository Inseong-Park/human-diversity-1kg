# 11_ld_decay.R — 대륙별 대표 집단의 LD decay
# 하는 일:
#   1. 대표 집단 5개(YRI, CEU, GIH, CHB, PEL)와 유럽 5개 혼합에서 같은 인원을 무작위로 뽑는다
#   2. 집단마다 MAF 0.05를 다시 적용하고 SNP를 솎아 낸 뒤, 500kb 안 SNP 쌍의 r²를 계산한다
#   3. 거리 구간별 평균 r²와 반감 거리를 표 8, 그림 7로 저장한다
# LD pruning은 하지 않는다 (LD 자체를 재는 분석이므로)
# r²는 PLINK 2의 --r2-unphased로 계산한다. PLINK 1.9의 --r2 기본값과 같은 계산이다

source("R/00_shared.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("work/ld", showWarnings = FALSE)

# 기준 B에 남은 사람만 쓴다 (친척 13명 제외)
pops <- read_tsv("work/pops.txt", col_types = cols(.default = "c")) |>
  rename(IID = `#IID`) |>
  semi_join(read_tsv("work/B.psam", col_types = cols(.default = "c")) |> rename(IID = `#IID`),
            by = "IID")

ld_pops <- c("YRI", "CEU", "GIH", "CHB", "PEL")
eur_pops <- c("CEU", "GBR", "FIN", "IBS", "TSI")


# ── 1. 같은 인원 뽑기 ────────────────────────────────────────
# r²는 인원이 적을수록 우연히 커지므로 모든 집단을 가장 작은 집단 인원에 맞춘다
n_each <- pops |> filter(POP %in% ld_pops) |> count(POP) |> pull(n) |> min()
cat("집단당 인원:", n_each, "\n")

set.seed(1)
samples <- map(ld_pops, \(p) {
  pops |> filter(POP == p) |> slice_sample(n = n_each) |> mutate(group = p)
}) |> list_rbind()

# 유럽 5개 혼합: 같은 인원을 집단마다 똑같이 나눠 뽑는다 (85명이면 17명씩)
eur_mix <- map(eur_pops, \(p) {
  pops |> filter(POP == p) |> slice_sample(n = n_each %/% length(eur_pops))
}) |> list_rbind() |> mutate(group = "EUR_mix")
samples <- bind_rows(samples, eur_mix)

groups <- c(ld_pops, "EUR_mix")
walk(groups, \(g) {
  samples |> filter(group == g) |> select(`#IID` = IID) |>
    write_tsv(paste0("work/ld/keep_", g, ".txt"))
})


# ── 2. r² 계산 ───────────────────────────────────────────────
# 1단계 --maf 0.05 : 그 집단 안에서 MAF를 다시 계산해 적용한다
# 2단계 --thin-count 10000 : 계산량을 줄이려고 SNP를 무작위로 1만 개만 남긴다 (약 3.5kb에 하나)
#   한 명령에 같이 쓰면 PLINK가 솎아내기를 먼저 해서 MAF 뒤에 SNP가 너무 적게 남는다. 그래서 두 단계로 나눈다
# --ld-window-r2 0 : r²가 작은 쌍도 모두 기록한다 (기본값은 0.2 이상만 기록)
calc_ld <- function(g) {
  out <- paste0("work/ld/", g)
  run_plink(c("--pfile", "work/raw", "--keep", paste0("work/ld/keep_", g, ".txt"),
              "--maf", "0.05", "--make-pgen", "--out", paste0(out, "_maf")))
  run_plink(c("--pfile", paste0(out, "_maf"),
              "--thin-count", "10000", "--seed", "1",
              "--r2-unphased", "cols=pos",
              "--ld-window-kb", "500", "--ld-window", "999999", "--ld-window-r2", "0",
              "--out", out))
  data.table::fread(paste0(out, ".vcor"), col.names = c("POS_A", "POS_B", "r2")) |>
    as_tibble() |>
    mutate(group = g, dist_kb = (POS_B - POS_A) / 1000)
}

ld <- map(groups, calc_ld) |> list_rbind()


# ── 3. 구간별 평균 r² ────────────────────────────────────────
# width kb 간격으로 묶는다. 그림 전체는 5kb, 반감 거리는 1kb 간격으로 계산한다
bin_ld <- function(width) {
  ld |>
    filter(dist_kb < 500) |>    # 500kb에 딱 걸린 몇 쌍만 있는 마지막 구간을 뺀다
    mutate(bin =floor(dist_kb / width) * width + width / 2) |>   # 구간 가운데 값
    group_by(group, bin) |>
    summarise(mean_r2 = mean(r2), n_pairs = n(), .groups = "drop")
}
decay <- bin_ld(5)
decay_fine <- bin_ld(1)


# ── 4. 표 8: 반감 거리 ───────────────────────────────────────
# 시작값: 첫 1kb 구간(0~1kb)의 평균 r²
# 바닥값: 400~500kb의 평균 r². 인원이 유한해서 관계없는 SNP 사이에도 남는 r²다
# 반감 거리(half_kb): 평균 r²가 처음으로 '바닥값 + (시작값 - 바닥값) / 2' 아래로 떨어지는 거리
# half_kb_simple: 타임테이블의 정의('시작값의 절반'). 바닥값을 0으로 본 것이라 함께 기록한다
# 떨어지는 지점은 앞뒤 구간 사이를 직선으로 이어 읽는다
cross_point <- function(x, y, threshold) {
  i <- which(y < threshold)[1]
  if (is.na(i) || i == 1) return(NA_real_)
  x[i - 1] + (threshold - y[i - 1]) / (y[i] - y[i - 1]) * (x[i] - x[i - 1])
}

half_distance <- decay_fine |>
  group_by(group) |>
  arrange(bin, .by_group = TRUE) |>
  summarise(
    r2_start = first(mean_r2),
    r2_floor = mean(mean_r2[bin >= 400]),
    half_kb = cross_point(bin, mean_r2, r2_floor + (r2_start - r2_floor) / 2),
    half_kb_simple = cross_point(bin, mean_r2, r2_start / 2),
    n_pairs = sum(n_pairs),
    .groups = "drop"
  ) |>
  left_join(read_csv("results/tables/table3_diversity.csv", show_col_types = FALSE) |>
              select(group = POP, He), by = "group") |>
  mutate(group = factor(group, levels = groups)) |>
  arrange(group)

write_csv(half_distance, "results/tables/table8_ld_half_distance.csv")
cat("\n표 8: 반감 거리 (half_kb: 바닥값 보정, half_kb_simple: 시작값의 절반)\n")
print(half_distance, width = Inf)


# ── 5. 그림 7: LD decay ──────────────────────────────────────
group_colours <- c(YRI = super_colours[["AFR"]], CEU = super_colours[["EUR"]],
                   GIH = super_colours[["SAS"]], CHB = super_colours[["EAS"]],
                   PEL = super_colours[["AMR"]], EUR_mix = "grey35")
group_labels <- c(YRI = "YRI (아프리카)", CEU = "CEU (유럽)", GIH = "GIH (남아시아)",
                  CHB = "CHB (동아시아)", PEL = "PEL (아메리카)", EUR_mix = "유럽 5개 혼합")

p <- decay |>
  mutate(group = factor(group, levels = groups)) |>
  ggplot(aes(bin, mean_r2, colour = group, linetype = group == "EUR_mix")) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = group_colours, labels = group_labels, name = NULL) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
  labs(x = "두 SNP 사이 거리 (kb)", y = "평균 r²",
       title = "대표 집단별 LD decay",
       subtitle = paste0("22번 염색체, 집단당 ", n_each, "명, 집단 안 MAF 0.05, SNP 1만 개, 5kb 구간 평균")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

ggsave("results/figures/fig7_ld_decay.png", p, width = 8, height = 5, dpi = 200, bg = "white")

# 앞쪽 100kb만 확대
ggsave("results/figures/fig7b_ld_decay_100kb.png",
       p + coord_cartesian(xlim = c(0, 100)) + labs(title = "대표 집단별 LD decay (0~100kb 확대)"),
       width = 8, height = 5, dpi = 200, bg = "white")
