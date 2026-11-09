# 03_kinship.R — KING 친족계수로 2,504명 안에 남은 친척 찾기
# 하는 일:
#   1. 모든 사람 쌍의 친족계수를 계산한다 (PLINK 2 KING)
#   2. 친족계수로 관계를 나누고, 같은 집단/다른 집단/다른 대륙 쌍이 몇 개인지 표로 만든다
#   3. 친족계수 분포를 그림으로 저장한다
#   4. 친척 쌍마다 한 명씩 골라 work/related.txt에 쓴다

suppressPackageStartupMessages(library(tidyverse))

plink2 <- "bin/plink2.exe"

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# PLINK를 실행하고, 실패하면 로그 마지막 부분을 보여주며 멈추는 함수
run_plink <- function(args) {
  out <- system2(plink2, args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cat(tail(out, 20), sep = "\n")
    stop("PLINK 실행 실패 (위 메시지 확인)")
  }
  invisible(out)
}


# ── 1. KING 친족계수 계산 ────────────────────────────────────
# --maf 0.05            : 흔한 SNP만 사용한다
# --make-king-table     : 모든 사람 쌍의 친족계수를 계산한다
# --king-table-filter   : 친족계수가 이 값 이상인 쌍만 파일에 적는다 (3촌 기준)
run_plink(c(
  "--pfile", "work/raw",
  "--maf", "0.05",
  "--make-king-table",
  "--king-table-filter", "0.0442",
  "--out", "work/king"
))


# ── 2. 결과 읽기: 쌍마다 집단 정보와 관계 붙이기 ────────────────
pops <- read_tsv("work/pops.txt", col_types = cols(.default = "c")) |>
  rename(IID = `#IID`)

king <- read_tsv("work/king.kin0",
                 col_types = cols(`#IID1` = "c", IID2 = "c", .default = "d")) |>
  rename(IID1 = `#IID1`) |>
  left_join(rename(pops, IID1 = IID, POP1 = POP, SUPER1 = SUPER), by = "IID1") |>
  left_join(rename(pops, IID2 = IID, POP2 = POP, SUPER2 = SUPER), by = "IID2") |>
  mutate(
    degree = case_when(
      KINSHIP >= 0.354  ~ "0. 동일인/일란성 쌍둥이",
      KINSHIP >= 0.177  ~ "1. 1촌",
      KINSHIP >= 0.0884 ~ "2. 2촌",
      TRUE              ~ "3. 3촌"
    )
  )


# ── 3. 표 2: 관계별 쌍 수 ───────────────────────────────────
# 진짜 친척이라면 대부분 같은 집단 안에 있어야 한다.
table2 <- king |>
  group_by(degree) |>
  summarise(
    n_pairs      = n(),
    same_pop     = sum(POP1 == POP2),
    diff_pop     = sum(POP1 != POP2),
    diff_super   = sum(SUPER1 != SUPER2),
    .groups = "drop"
  )

write_csv(table2, "results/tables/table2_kinship_pairs.csv")
print(table2)


# ── 4. 친족계수 분포 그림 ────────────────────────────────────
# 쌍 수가 관계별로 크게 달라서 세로축을 로그로 그린다
p <- ggplot(king, aes(KINSHIP)) +
  geom_histogram(bins = 80) +
  geom_vline(xintercept = c(0.0884, 0.177, 0.354), linetype = "dashed", colour = "red") +
  scale_y_log10() +
  labs(x = "KING 친족계수 (22번 염색체)", y = "쌍 수 (로그)",
       title = "친족계수 0.0442 이상인 쌍의 분포",
       subtitle = "빨간 선: 2촌, 1촌, 동일인 기준")

ggsave("results/figures/fig_kinship_hist.png", p, width = 7, height = 4.5, dpi = 200)


# ── 5. 친척 쌍마다 한 명씩 고르기 ─────────────────────────────
# 여러 쌍에 걸친 사람을 먼저 뺀다. 그래야 가장 적은 인원으로 모든 쌍을 없앨 수 있다.
# 횟수가 같으면 ID 순서로 정해서, 다시 돌려도 항상 같은 결과가 나오게 한다.
pick_relatives <- function(cutoff) {
  pairs <- king |> filter(KINSHIP >= cutoff) |> select(IID1, IID2)
  removed <- character(0)
  while (nrow(pairs) > 0) {
    counts <- table(c(pairs$IID1, pairs$IID2))
    pick <- sort(names(counts)[counts == max(counts)])[1]
    removed <- c(removed, pick)
    pairs <- pairs |> filter(IID1 != pick, IID2 != pick)
  }
  sort(removed)
}

# 기준을 바꾸면 몇 명이 빠지는지 비교한다
cutoff_table <- tibble(cutoff = c(0.0442, 0.0884, 0.177, 0.354)) |>
  mutate(n_removed = map_int(cutoff, \(x) length(pick_relatives(x))))
print(cutoff_table)

# 제외 기준: 1촌 이상
# 22번 염색체 하나로는 약한 관계를 믿을 수 없어서 2촌, 3촌은 제외 기준으로 쓰지 않는다.
# 이유는 notes.md 비교 1 참고.
cutoff <- 0.177

related <- pick_relatives(cutoff)
write_tsv(tibble(`#IID` = related), "work/related.txt")
cat("\n제외 기준", cutoff, "→ work/related.txt에", length(related), "명 기록\n")
