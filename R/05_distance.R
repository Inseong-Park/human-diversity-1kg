# 05_distance.R — 집단 좌표와 아디스아바바로부터의 이동 거리
# 하는 일:
#   1. 26개 집단의 좌표와 선택 근거를 표로 만든다
#   2. 대륙별 경로의 경유지를 따라 거리를 더한다
#   3. 표 4를 저장하고, 예측을 확인할 집단의 거리를 출력한다
# 좌표 선택 기준은 notes.md Day 6 참고 (표본을 모은 곳이 아니라 조상이 살던 곳)

source("R/00_shared.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)


# ── 1. 출발점과 경유지 ───────────────────────────────────────
waypoints <- tribble(
  ~name,       ~lat,   ~lon,
  "Addis",      9.03,   38.74,
  "Cairo",     30.04,   31.24,
  "Istanbul",  41.01,   28.98,
  "Anadyr",    64.73,  177.50,
  "PrinceRupert", 54.31, -130.32
)

# 대륙별로 지나가는 경유지 (도착지 직전까지)
routes <- list(
  AFR = "Addis",
  EUR = c("Addis", "Cairo", "Istanbul"),
  SAS = c("Addis", "Cairo", "Istanbul"),
  EAS = c("Addis", "Cairo", "Istanbul"),
  AMR = c("Addis", "Cairo", "Istanbul", "Anadyr", "PrinceRupert")
)

# 경유지를 차례로 지나 (lat, lon)까지 가는 거리(km)
route_km <- function(super, lat, lon) {
  stops <- tibble(name = routes[[super]]) |> left_join(waypoints, by = "name")
  pts <- rbind(cbind(stops$lon, stops$lat), c(lon, lat))
  n <- nrow(pts)
  sum(geosphere::distHaversine(pts[-n, , drop = FALSE], pts[-1, , drop = FALSE])) / 1000
}


# ── 2. 집단 좌표표 ───────────────────────────────────────────
coords <- tribble(
  ~POP,  ~SUPER, ~place,          ~lat,    ~lon,    ~sampled_in,           ~reason,
  "YRI", "AFR",  "이바단",          7.38,    3.90,  "나이지리아 이바단",   "표본지",
  "ESN", "AFR",  "Esan 지역",       6.75,    6.13,  "나이지리아",          "Esan 거주 지역 근사",
  "GWD", "AFR",  "브리카마",       13.27,  -16.65,  "감비아 서부",         "Western Division 안의 도시로 근사",
  "MSL", "AFR",  "보",              7.96,  -11.74,  "시에라리온",          "Mende 거주 지역의 도시로 근사",
  "LWK", "AFR",  "웨부예",          0.62,   34.77,  "케냐 웨부예",         "표본지",
  "ASW", "AFR",  "이바단",          7.38,    3.90,  "미국 남서부",         "아프리카 조상 쪽 서아프리카 근사. 유럽 조상도 있음",
  "ACB", "AFR",  "이바단",          7.38,    3.90,  "바베이도스",          "아프리카 조상 쪽 서아프리카 근사. 유럽 조상도 있음",
  "CEU", "EUR",  "브뤼셀",         50.85,    4.35,  "미국 유타",           "북서유럽계 조상. 북서유럽 근사",
  "GBR", "EUR",  "리즈",           53.80,   -1.55,  "잉글랜드·스코틀랜드", "도시 미기재. 중간 지점 근사",
  "FIN", "EUR",  "헬싱키",         60.17,   24.94,  "핀란드",              "도시 미기재. 수도로 근사",
  "IBS", "EUR",  "마드리드",       40.42,   -3.70,  "스페인",              "도시 미기재. 중심부로 근사",
  "TSI", "EUR",  "피렌체",         43.77,   11.26,  "이탈리아 토스카나",   "토스카나의 중심 도시",
  "GIH", "SAS",  "아마다바드",     23.02,   72.57,  "미국 휴스턴",         "구자라트계 조상",
  "PJL", "SAS",  "라호르",         31.55,   74.34,  "파키스탄 라호르",     "표본지",
  "BEB", "SAS",  "다카",           23.81,   90.41,  "방글라데시",          "도시 미기재. 수도로 근사",
  "STU", "SAS",  "자프나",          9.66,   80.01,  "영국",                "스리랑카 타밀 거주 지역",
  "ITU", "SAS",  "하이데라바드",   17.39,   78.49,  "영국",                "텔루구어 사용 지역",
  "CHB", "EAS",  "베이징",         39.90,  116.41,  "중국 베이징",         "표본지",
  "JPT", "EAS",  "도쿄",           35.68,  139.69,  "일본 도쿄",           "표본지",
  "CHS", "EAS",  "창사",           28.23,  112.94,  "중국 남부",           "지역 미기재. 후난으로 근사 (미확인)",
  "CDX", "EAS",  "징훙",           22.01,  100.80,  "중국 시솽반나",       "표본지",
  "KHV", "EAS",  "호찌민",         10.82,  106.63,  "베트남 호찌민",       "표본지",
  "MXL", "AMR",  "멕시코시티",     19.43,  -99.13,  "미국 LA",             "멕시코계 조상",
  "PUR", "AMR",  "산후안",         18.47,  -66.11,  "푸에르토리코",        "도시 미기재. 수도로 근사",
  "CLM", "AMR",  "메데인",          6.24,  -75.57,  "콜롬비아 메데인",     "표본지",
  "PEL", "AMR",  "리마",          -12.05,  -77.04,  "페루 리마",           "표본지"
)

# 26개 집단이 빠짐없이, 대륙까지 맞게 들어갔는지 확인한다
pops <- read_tsv("work/pops.txt", col_types = cols(.default = "c")) |> distinct(POP, SUPER)
stopifnot(nrow(coords) == 26,
          setequal(paste(coords$POP, coords$SUPER), paste(pops$POP, pops$SUPER)))


# ── 3. 거리 계산과 표 4 ──────────────────────────────────────
table4 <- coords |>
  mutate(
    route = map_chr(SUPER, \(s) paste(routes[[s]], collapse = " → ")),
    distance_km = round(pmap_dbl(list(SUPER, lat, lon), route_km)),
    admixed = POP %in% admixed_pops
  ) |>
  arrange(distance_km)

write_csv(table4, "results/tables/table4_coordinates_distance.csv")
print(select(table4, POP, SUPER, place, distance_km), n = Inf)


# ── 4. 예측 확인 ─────────────────────────────────────────────
cat("\n예측 1, 2에 나온 집단의 거리:\n")
table4 |>
  filter(POP %in% c("GWD", "TSI", "GIH", "FIN", "PJL")) |>
  select(POP, place, distance_km) |>
  print()

# 예측 3: CEU를 표본지(유타, 솔트레이크시티)로 넣으면 아메리카 경로를 지난다
ceu_utah <- round(route_km("AMR", 40.76, -111.89))
cat("\nCEU 거리: 브뤼셀", table4$distance_km[table4$POP == "CEU"], "km, 유타", ceu_utah, "km\n")
