# 08_fst_tree.R — 집단 쌍별 F_ST와 Neighbor-Joining 계통수
# 하는 일:
#   1. PLINK로 26개 집단의 모든 쌍에 대해 F_ST를 계산한다
#   2. 쌍별 표와 대륙 조합별 요약표를 저장한다
#   3. heatmap과 계통수 그림을 저장한다

source("R/00_shared.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)


# ── 1. F_ST 계산 ─────────────────────────────────────────────
# --pheno : 집단 이름이 든 파일 (POP 열)
# --fst POP method=hudson : POP의 모든 집단 쌍에 대해 Hudson 방식으로 계산한다
run_plink(c("--pfile", "work/B",
            "--pheno", "work/pops.txt",
            "--fst", "POP", "method=hudson",
            "--out", "work/B_fst"))

pops <- read_tsv("work/pops.txt", col_types = cols(.default = "c")) |>
  rename(IID = `#IID`) |>
  distinct(POP, SUPER)

fst <- read_tsv("work/B_fst.fst.summary", show_col_types = FALSE) |>
  rename(POP1 = `#POP1`, FST = HUDSON_FST) |>
  left_join(rename(pops, POP1 = POP, SUPER1 = SUPER), by = "POP1") |>
  left_join(rename(pops, POP2 = POP, SUPER2 = SUPER), by = "POP2")

stopifnot(nrow(fst) == 325)   # 26개 집단의 쌍 수


# ── 2. 쌍별 표와 대륙 조합별 요약 ────────────────────────────
# 혼혈 5개는 요약에서 빼고, PEL은 한 집단으로 따로 본다
group_of <- function(pop, super) if_else(pop == "PEL", "PEL", super)
group_levels <- c(setdiff(super_levels, "AMR"), "PEL")

pair_name <- function(a, b) {
  ia <- match(a, group_levels)
  ib <- match(b, group_levels)
  if_else(ia == ib, paste0(a, " 안"),
          paste(group_levels[pmin(ia, ib)], group_levels[pmax(ia, ib)], sep = "–"))
}

fst_summary <- fst |>
  filter(!POP1 %in% admixed_pops, !POP2 %in% admixed_pops) |>
  mutate(G1 = group_of(POP1, SUPER1), G2 = group_of(POP2, SUPER2),
         pair = pair_name(G1, G2)) |>
  group_by(pair) |>
  summarise(n_pairs = n(), mean_fst = mean(FST), min_fst = min(FST), max_fst = max(FST),
            .groups = "drop") |>
  arrange(desc(mean_fst))

write_csv(select(fst, POP1, POP2, SUPER1, SUPER2, FST), "results/tables/table_fst_pairs.csv")
write_csv(fst_summary, "results/tables/table_fst_summary.csv")
print(fst_summary, n = Inf)

cat("\nF_ST가 가장 큰 집단 쌍 5개:\n")
fst |> slice_max(FST, n = 5) |> select(POP1, POP2, FST) |> print()


# ── 3. 26 x 26 행렬 ─────────────────────────────────────────
# 대륙 순서로 정렬한다. 혼혈 집단은 대륙 안 뒤쪽에 둔다
pop_order <- pops |>
  mutate(SUPER = factor(SUPER, levels = super_levels),
         admixed = POP %in% admixed_pops) |>
  arrange(SUPER, admixed, POP) |>
  pull(POP)

fst_long <- bind_rows(
  select(fst, A = POP1, B = POP2, FST),
  select(fst, A = POP2, B = POP1, FST),
  tibble(A = pop_order, B = pop_order, FST = 0)
)

fst_matrix <- fst_long |>
  pivot_wider(names_from = B, values_from = FST) |>
  column_to_rownames("A") |>
  as.matrix()
fst_matrix <- fst_matrix[pop_order, pop_order]

# 가까운 집단 쌍에서는 추정값이 0보다 조금 작게 나올 수 있어 거리로 쓸 때 0으로 맞춘다
fst_dist <- pmax(fst_matrix, 0)


# ── 4. 그림 4: heatmap ───────────────────────────────────────
p <- fst_long |>
  mutate(A = factor(A, levels = pop_order), B = factor(B, levels = rev(pop_order))) |>
  ggplot(aes(A, B, fill = pmax(FST, 0))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "#eef4fc", high = "#123a6b", name = expression(F[ST])) +
  coord_fixed() +
  labs(x = NULL, y = NULL, title = "집단 쌍별 F_ST (기준 B)",
       subtitle = "대륙 순서: AFR, EUR, SAS, EAS, AMR. 대륙 안에서 혼혈 집단은 뒤쪽") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank())

ggsave("results/figures/fig4_fst_heatmap.png", p, width = 8, height = 7.5, dpi = 200, bg = "white")


# ── 5. 그림 4: NJ 계통수 ─────────────────────────────────────
# NJ는 가지마다 변화 속도가 같다고 가정하지 않는다
# NJ 결과에는 뿌리가 없다. 집단 이름을 읽기 쉽도록 아프리카 묶음과 나머지를 가르는 가지에 뿌리를 놓아 그린다
# 뿌리 위치는 그림을 읽기 위한 표시이고, 이 분석으로 뿌리를 검증한 것은 아니다
tree <- ape::nj(as.dist(fst_dist))
afr_tips <- pops$POP[pops$SUPER == "AFR"]
stopifnot(ape::is.monophyletic(ape::unroot(tree), afr_tips) |
            ape::is.monophyletic(ape::unroot(tree), setdiff(pops$POP, afr_tips)))
tree <- ape::root(tree, outgroup = afr_tips, resolve.root = TRUE)

# NJ는 가지 길이가 음수로 나올 수 있다. 나무 모양으로 설명되지 않는 집단(혼혈 등)에서 생긴다
# 음수 가지가 어디인지 기록하고, 그림에서는 0으로 그린다
negative <- which(tree$edge.length < 0)
tip_edges <- tree$edge[negative, 2] <= length(tree$tip.label)
cat("\n음수 가지:", length(negative), "개 / 집단에 바로 붙은 가지:",
    paste(tree$tip.label[tree$edge[negative[tip_edges], 2]], collapse = ", "), "\n")
tree$edge.length[negative] <- 0

tip_colour <- super_colours[pops$SUPER[match(tree$tip.label, pops$POP)]]

png("results/figures/fig4_nj_tree.png", width = 7.5, height = 7, units = "in", res = 200, bg = "white")
par(mar = c(1, 1, 2, 1))
ape::plot.phylo(tree, type = "phylogram", tip.color = tip_colour, cex = 0.8,
                edge.color = "grey55", no.margin = FALSE, main = "NJ tree (F_ST, baseline B)")
ape::add.scale.bar(cex = 0.8)
legend("topright", legend = super_levels, text.col = super_colours[super_levels], bty = "n", cex = 0.9)
dev.off()
invisible(NULL)
