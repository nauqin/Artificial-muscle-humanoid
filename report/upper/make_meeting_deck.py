"""미팅 발표 PDF (16:9, matplotlib). 숫자는 report/upper/README.md 9·10절과 같은 값."""
import sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import font_manager
import matplotlib.image as mpimg

for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f:
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False
DATASETS = [
    ["UE_ADLdatabase (16명)", "무료(MIT), 팔꿈치·손목·손 마커 다 있음", "물체 무게 없음"],
    ["Kinematic+EMG (2025)", "무료(CC BY), 파이프 돌리기로 회내·회외", "마커 적음"],
    ["Yale ADL (12명)", "Herneth 2024 가 MoBL-ARMS 로 IK·ID 완료", "유료 (검산 기준으로만)"],
]
ROOT = "/Users/nyubin/Desktop/scHooL/의료로봇연구실/addb_opensim/"
OUT = sys.argv[1]
NAVY, ACC, GREY, LIGHT, GOOD, BAD = "#1f3a5f", "#d1622b", "#5b6573", "#eef2f7", "#2e7d4f", "#b3261e"
W, H = 13.333, 7.5
page_no = [0]


def slide(title, sub=None):
    fig = plt.figure(figsize=(W, H)); fig.patch.set_facecolor("white")
    page_no[0] += 1
    fig.add_artist(plt.Rectangle((0, 0.885), 1, 0.115, transform=fig.transFigure, color=NAVY))
    fig.text(0.04, 0.935, title, fontsize=24, color="white", weight="bold", va="center")
    if sub: fig.text(0.04, 0.845, sub, fontsize=13.5, color=GREY, va="center")
    fig.text(0.96, 0.025, f"{page_no[0]}", fontsize=10, color=GREY, ha="right")
    fig.text(0.04, 0.025, "세종대 의료로봇연구실 · 인공근육 휴머노이드 근육 배치 · 2026-11", fontsize=9, color=GREY)
    return fig


def table(fig, rect, header, rows, col_w=None, fs=12.5, hl_rows=(), hl_color=LIGHT, bold_first=True):
    ax = fig.add_axes(rect); ax.axis("off")
    t = ax.table(cellText=rows, colLabels=header, loc="upper center", cellLoc="center", colLoc="center", colWidths=col_w)
    t.auto_set_font_size(False); t.set_fontsize(fs); t.scale(1, 1.9)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("#c9d1db"); cell.set_linewidth(0.8)
        if r == 0: cell.set_facecolor(NAVY); cell.get_text().set_color("white"); cell.get_text().set_weight("bold")
        elif (r - 1) in hl_rows: cell.set_facecolor(hl_color)
        if r > 0 and c == 0 and bold_first: cell.get_text().set_weight("bold")
    return t


def bullets(fig, x, y, items, fs=15, dy=0.065, color="#222"):
    for it in items:
        if isinstance(it, tuple): txt, c = it
        else: txt, c = it, color
        fig.text(x, y, "•  " + txt, fontsize=fs, color=c, va="top", wrap=True); y -= dy
    return y


def stat(fig, x, y, big, small, color=NAVY):
    fig.text(x, y, big, fontsize=40, weight="bold", color=color, ha="center")
    fig.text(x, y - 0.065, small, fontsize=13, color=GREY, ha="center")


with PdfPages(OUT) as pdf:
    # 1 표지
    fig = plt.figure(figsize=(W, H)); fig.patch.set_facecolor(NAVY); page_no[0] += 1
    fig.text(0.07, 0.62, "인공근육 휴머노이드 근육 배치", fontsize=38, color="white", weight="bold")
    fig.text(0.07, 0.53, "키 160 cm · 40 kg 전신 설계와 검증", fontsize=26, color="#cfd8e3")
    fig.text(0.07, 0.40, "하체 8 · 몸통 4 · 어깨 6 · 팔꿈치 2 (한쪽 20개, 양쪽 36개)", fontsize=18, color="white")
    fig.text(0.07, 0.16, "세종대학교 의료로봇연구실  나유빈\n2026-11 미팅", fontsize=15, color="#cfd8e3", linespacing=1.6)
    pdf.savefig(fig); plt.close(fig)

    # 2 요약
    fig = slide("요약", "11/3 미팅 결정(160 cm 복구, 근육 길이 단축)을 반영해 전신을 다시 설계했고, 모든 부위가 검증을 통과했다")
    stat(fig, 0.17, 0.62, "36개", "케이블 근육 (양쪽)")
    stat(fig, 0.40, 0.62, "0개", "모듈 배정 불가")
    stat(fig, 0.63, 0.62, "≤ 2.6 %", "케이블이 못 낸 토크 (95퍼센타일)", GOOD)
    stat(fig, 0.86, 0.62, "약 5.0 kg", "모듈 무게 (park2019 기준, 냉각·구조 제외)", ACC)
    bullets(fig, 0.06, 0.42, [
        "근육 길이 단축 = 과제 계획서 모듈표의 힘·스트로크 상한으로 경로를 설계. 하체·몸통은 케이블마다 등급을 먼저 정함",
        "어깨는 경로 2단계를 신뢰 구간으로 제대로 돌리니 6개로 걷기 2·들기 2시행 통과. 회전 지렛대(하드웨어)는 불필요",
        "시너지 배분비와 비교: 같은 개수에서 직접 최적화가 필요한 총 토크 용량이 작거나 같다",
        ("제안: 팀원이 찾은 손목 데이터셋으로 손목 최적화 추가.  질문: 들기 하중 목표 (3.79 kg 이면 충분한지)", ACC),
    ], fs=15, dy=0.075)
    pdf.savefig(fig); plt.close(fig)

    # 3 지난 미팅 대비
    fig = slide("지난 미팅(11/3) 결정 대비", "결정 7·8 을 먼저 끝냈고, 결정 요청 2건 중 어깨 지렛대는 소프트웨어로 해결")
    table(fig, [0.04, 0.12, 0.92, 0.68], ["항목", "결정 / 요청", "결과"], [
        ["7", "모델 키 160 cm 복구, 40 kg 유지", "완료. 뼈대 길이비 0.988, 상체 ID 16건·하체 ID 3건 재계산"],
        ["8", "근육 길이 단축", "완료. 모듈표 힘·스트로크 상한으로 경로 설계, 케이블별 모듈 배정 (배정 불가 0)"],
        ["②", "어깨 회전용 지렛대", "불필요. 어깨 6개 + 신뢰 구간 2단계로 다른 사람 들기까지 0.6 %"],
        ["①", "손목 포함 여부", "팀원이 손목 데이터셋 조사 완료 → 손목만 최적화해 추가하는 것을 제안"],
        ["추가", "상체 근육 시너지 분석", "완료. 몸통(Gait2392 근육 6개)·팔(MoBL-ARMS 24개)로 직접 최적화와 비교"],
    ], col_w=[0.08, 0.30, 0.62], fs=13, hl_rows=(2,), hl_color="#e7f3ec")
    pdf.savefig(fig); plt.close(fig)

    # 4 방법
    fig = slide("설계 방법", "공통 과정은 같고, 모듈을 정하는 순서에 따라 방식 1·2·3 으로 갈린다")
    def box(x, y, w, h, title, body, fc=LIGHT, tc=NAVY):
        fig.add_artist(plt.Rectangle((x, y), w, h, transform=fig.transFigure, color=fc, ec=NAVY, lw=1.3))
        fig.text(x + w / 2, y + h * 0.68, title, fontsize=12.5, weight="bold", color=tc, ha="center", va="center")
        fig.text(x + w / 2, y + h * 0.30, body, fontsize=10, color="#333", ha="center", va="center", linespacing=1.3)
    def arrow(x0, y0, x1, y1):
        fig.add_artist(matplotlib.patches.FancyArrowPatch((x0, y0), (x1, y1), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=14, color=ACC, lw=1.6))
    # 공통 앞
    box(0.03, 0.485, 0.13, 0.17, "관절 토크", "걷기·스쿼트·들기\n역동역학(ID)")
    box(0.19, 0.485, 0.15, 0.17, "관절 배분비·\n토크 용량", "직접 최적화\n(선형계획법 + Powell)")
    arrow(0.16, 0.57, 0.19, 0.57)
    # 세 갈래
    lanes = [("방식 1  사후 배정", "경로 설계", "부위 전체에 상한 하나\n(하체 115·상체 80 mm)", "모듈 배정", "결과 보고 등급"),
             ("방식 2  모듈 먼저", "모듈 미리 정함", "W = 1.1·c·Δθ 로 등급", "경로 설계", "케이블별 등급 상한"),
             ("방식 3  보정 모듈 먼저", "모듈 미리 정함", "W × 실측 비율로 등급", "경로 설계", "케이블별 등급 상한")]
    ys = [0.66, 0.49, 0.32]
    for (name, a1, b1, a2, b2), y in zip(lanes, ys):
        fig.text(0.375, y + 0.135, name, fontsize=12, weight="bold", color=ACC)
        box(0.375, y, 0.15, 0.12, a1, b1); box(0.56, y, 0.15, 0.12, a2, b2)
        arrow(0.34, 0.57, 0.375, y + 0.06); arrow(0.525, y + 0.06, 0.56, y + 0.06); arrow(0.71, y + 0.06, 0.745, 0.57)
    box(0.745, 0.485, 0.14, 0.17, "검증", "MocoInverse 로 재현\n못 낸 토크 < 10 %")
    fig.text(0.815, 0.465, "경로 2단계: 기본 또는\n신뢰 구간(±1 cm 씩)", fontsize=10, color=GREY, ha="center", va="top")
    # 어디에 쓰였나
    table(fig, [0.04, 0.05, 0.92, 0.2], ["설계", "하체", "몸통", "어깨", "팔꿈치"], [
        ["이전 설계 (맥, 한쪽 22개)", "방식 1", "방식 1", "방식 1, 8개 (2단계 기본)", "방식 1"],
        ["최종 (한쪽 20개)", "방식 3", "방식 2", "방식 1, 6개 (2단계 신뢰 구간)", "방식 1"],
    ], col_w=[0.2, 0.2, 0.14, 0.3, 0.16], fs=12)
    pdf.savefig(fig); plt.close(fig)

    # 4-1 배분비와 모듈 연결
    fig = slide("배분비에서 모듈 등급까지", "배분비가 정해지면 경로를 찾기 전에 필요한 모듈 등급과 모멘트암 범위를 계산할 수 있다")
    chain = [("관절 배분비 v", "어느 관절을\n어떤 비율로"), ("토크 용량 c", "모든 프레임을 내는\n최소 세기 (N·m)"), ("움직이는 각도 Δθ", "동작 중 v 방향\n각도 범위 (rad)"),
             ("필요한 일 W", "W = 1.1 × c × Δθ\n(J)"), ("모듈 등급", "힘 상한 × 스트로크\n≥ W 인 가장 작은 등급"), ("모멘트암 범위", "1.1c ÷ 힘 상한 ~\n스트로크 상한 ÷ Δθ")]
    for i, (h, d) in enumerate(chain):
        x = 0.035 + i * 0.158
        fig.add_artist(plt.Rectangle((x, 0.52), 0.14, 0.23, transform=fig.transFigure, color=LIGHT if i not in (3,) else "#fbe9df", ec=NAVY, lw=1.3))
        fig.text(x + 0.07, 0.70, h, fontsize=14, weight="bold", color=NAVY, ha="center", va="center")
        fig.text(x + 0.07, 0.59, d, fontsize=11, color="#333", ha="center", va="center", linespacing=1.4)
        if i < 5: fig.text(x + 0.149, 0.635, "▶", fontsize=13, color=ACC, ha="center", va="center")
    bullets(fig, 0.05, 0.44, [
        "모멘트암을 키우면 힘은 줄고 길이 변화는 는다 → 힘 × 길이 변화 ≈ W 로 거의 일정. 그래서 경로 전에 등급을 고를 수 있다",
        "예 (설명용): c = 50 N·m, Δθ = 1 rad → W = 55 J → M2 (1,200 N × 80 mm = 96 J), 모멘트암 4.6~8 cm 안에서 부착점을 찾는다",
        ("한계: 모멘트암이 자세마다 변해 실제 힘 × 길이 변화는 W 의 1.1~3.4배 → 실측 비율을 곱한 '보정' 등급이 필요 (윈도우 실험)", GREY),
    ], fs=13.5, dy=0.08)
    pdf.savefig(fig); plt.close(fig)

    # 4-2 왜 이렇게 설정했나
    fig = slide("왜 이렇게 설정했나 — 두 가지를 검토", "① 배분비를 시너지로 정할까, 관절 토크에서 직접 정할까   ② 모듈을 경로 뒤에 붙일까, 먼저 정하고 상한을 걸까")
    fig.text(0.04, 0.78, "① 배분비: 시너지 vs 직접 최적화 (같은 개수, 같은 평가)", fontsize=15, weight="bold", color=NAVY)
    table(fig, [0.04, 0.50, 0.44, 0.25], ["부위", "직접 최적화", "시너지"], [
        ["하체 8", "237 N·m", "311 N·m"],
        ["몸통 4", "156 N·m", "304 N·m"],
        ["어깨 6", "81 N·m", "344 N·m (토크 부족)"],
    ], col_w=[0.3, 0.3, 0.4], fs=12.5)
    fig.text(0.04, 0.44, "→ 직접 최적화 채택. 시너지는 비교 기준", fontsize=13.5, weight="bold", color=GOOD)
    fig.text(0.52, 0.78, "② 모듈: 사후 배정 vs 먼저 정하기 (윈도우 실험)", fontsize=15, weight="bold", color=NAVY)
    table(fig, [0.52, 0.50, 0.45, 0.25], ["부위", "먼저 정하기 결과"], [
        ["몸통", "더 가벼움 (573 → 220 g), 모두 0 %"],
        ["하체", "보정하면 배정 불가 0, 모두 통과"],
        ["팔꿈치 / 어깨", "차이 없음 / 들기 미달"],
    ], col_w=[0.3, 0.7], fs=12.5)
    fig.text(0.52, 0.44, "→ 부위별로 선택 (몸통·하체에 유리)", fontsize=13.5, weight="bold", color=GOOD)
    bullets(fig, 0.04, 0.34, [
        "두 단계는 순서대로 붙는다: 배분비 (시너지 또는 직접) → 모듈 (사후 또는 먼저) → 부착점",
        "시너지 배분비는 필요한 토크 용량이 1.3~2배라, 모듈을 먼저 정해도 더 무거운 등급이 나올 것으로 예상 (시너지 + 먼저 정하기 조합은 아직 안 돌림)",
        "VAF 로 고른 시너지 수가 필요한 케이블 수와 꼭 같지는 않다 (몸통 4 = 4, 어깨는 시너지 6개로 토크를 다 못 냄)",
    ], fs=12.5, dy=0.075)
    pdf.savefig(fig); plt.close(fig)

    # 5 최종 결과
    fig = slide("최종 설계와 검증 (160 cm · 40 kg)", "검증 = 케이블이 못 내 보조 액추에이터가 떠맡은 토크의 95퍼센타일 ÷ 측정 토크 피크. 기준 10 %")
    table(fig, [0.04, 0.30, 0.92, 0.50], ["부위", "케이블 (한쪽)", "케이블 힘 (N)", "길이 변화 (mm)", "걷기", "들기 · 스쿼트"], [
        ["하체 (방식 3)", "8", "198~2,926", "36~119", "2.6 % · 1.4 %", "스쿼트 1.5 %"],
        ["몸통 (방식 2)", "4", "447~1,173", "54~66", "0 % · 0 %", "0 % · 0 %"],
        ["어깨 (방식 1, 신뢰 구간)", "6", "1,314~2,145", "65~83", "0.6 % · 0.0 %", "0.1 % · 0.6 %"],
        ["팔꿈치", "2", "1,835~2,327", "81~84", "0 % · 0 %", "0 % · 0 %"],
    ], col_w=[0.22, 0.11, 0.15, 0.14, 0.18, 0.2], fs=13)
    bullets(fig, 0.06, 0.40, [
        "걷기는 설계에 쓴 사람(subject3)과 다른 사람(subject9, 1.6 m/s), 들기는 subject1·subject6 찬 물병. 하체는 LabValidation 걷기·Moore 걷기·스쿼트",
        "어깨 힘은 ×1.5 (상한 2,700 N 가정), 팔꿈치 폄은 팔꿈치 뒤 지렛대(부착 자리) 전제",
    ], fs=13, dy=0.065)
    pdf.savefig(fig); plt.close(fig)

    # 6 모듈
    fig = slide("모듈 배정 (최종)", "하체·몸통은 등급을 먼저 정해 케이블마다 그 상한을 걸었고, 어깨·팔꿈치는 경로 뒤에 들어가는 가장 작은 등급")
    table(fig, [0.04, 0.38, 0.55, 0.42], ["등급 (힘 / 스트로크)", "하체", "몸통", "어깨", "팔꿈치", "양쪽 합"], [
        ["M1 (1,000~5,000 N / 40~120 mm)", "2", "0", "4", "2", "16"],
        ["M2 (400~1,200 N / 30~80 mm)", "4", "2", "0", "0", "10"],
        ["M3 (150~450 N / 25~60 mm)", "0", "2", "0", "0", "2"],
        ["M4 (40~200 N / 35~90 mm)", "1", "0", "0", "0", "2"],
        ["M6 (50~1,500 N / 20~80 mm)", "1", "0", "2", "0", "6"],
    ], col_w=[0.40, 0.10, 0.10, 0.10, 0.12, 0.16], fs=12)
    fig.text(0.63, 0.75, "모듈 무게 추정", fontsize=17, weight="bold", color=NAVY)
    bullets(fig, 0.63, 0.69, [
        "케이블 힘 합 약 53,700 N (양쪽)",
        "park2019 SMA 스프링 묶음 기준",
        "  130 N / 12 g ≈ 10.8 N/g → 약 5.0 kg",
        "(계획서 목표 15 N/g 이면 3.6 kg)",
        ("냉각·구조·케이블 무게는 제외", GREY),
    ], fs=14, dy=0.06)
    bullets(fig, 0.06, 0.28, [
        "하체 케이블 대부분은 전체 길이 80 cm 대지만 모듈이 실제로 움직이는 길이는 36~119 mm → 나머지는 늘어나지 않는 케이블",
        "이전 설계(사후 배정)보다 M1 이 23 → 16개로 줄고, 몸통에는 M1 이 없다 (M3 2 + M2 2)",
    ], fs=13, dy=0.065)
    pdf.savefig(fig); plt.close(fig)

    # 7 어깨
    fig = slide("어깨: 경로를 제대로 맞추니 6개로 된다", "처음 6개는 팔을 뒤로 젖힐 때(-22°) 앞으로 당기는 케이블이 어깨 중심 가까이 지나가 걷기를 못 냈다")
    try:
        img = mpimg.imread(ROOT + "report/upper/fig_shoulder_flex_flip_40kg_h160_s80_f15.png")
        ax = fig.add_axes([0.53, 0.30, 0.45, 0.53]); ax.imshow(img); ax.axis("off")
    except Exception:
        pass
    table(fig, [0.03, 0.42, 0.47, 0.40], ["어깨 설계", "걷기", "들기"], [
        ["6개, 2단계 기본 (11/4)", "87 %", "17 %"],
        ["8개, 2단계 기본 + 힘 재계산", "0 % · 0 %", "0.1 % · 2.1 %"],
        ["6개, 신뢰 구간 2단계 (최종)", "0.6 % · 0.0 %", "0.1 % · 0.6 %"],
        ["  (시드 6개 중 통과)", "4개", ""],
    ], col_w=[0.52, 0.24, 0.24], fs=12.5, hl_rows=(2,), hl_color="#e7f3ec")
    bullets(fig, 0.04, 0.31, [
        "기존 2단계는 1단계 답이 스트로크 상한 근처면 시작점보다\n    나쁜 답을 내서 1단계 답으로 되돌아갔다 (어깨에서 한 번도 안 돎)",
        "신뢰 구간: 지금 위치 ±1 cm 안에서 찾고 옮기기를 반복",
        "시드 여러 개 중 2단계 잔차가 가장 작은 설계를 고름\n    (시드 4: 8.9 %). 걷기 0.6·0.0 %, 들기 0.1·0.6 %",
    ], fs=12.5, dy=0.08)
    pdf.savefig(fig); plt.close(fig)

    # 8 시너지
    fig = slide("시너지 배분비와 비교", "상체 근육 시너지를 처음 뽑아 같은 입력·같은 평가로 비교 (근육 힘 분배 → NMF → 근육 힘을 시너지에 귀속)")
    table(fig, [0.04, 0.40, 0.92, 0.40], ["부위", "개수", "직접 최적화 (N·m)", "시너지 (N·m)", "근육 모델"], [
        ["하체", "8", "237.4", "310.6", "Rajagopal 근육 80개 (9월)"],
        ["몸통", "4", "155.5", "304.2", "Gait2392 허리 근육 6개"],
        ["어깨", "6", "80.8", "344.1 (15 N·m 못 냄)", "MoBL-ARMS 어깨·팔꿈치 24개"],
        ["어깨", "8 (일하는 수)", "79.4", "93.7 (후보 11개)", "같음"],
    ], col_w=[0.1, 0.14, 0.2, 0.26, 0.3], fs=13.5)
    bullets(fig, 0.06, 0.33, [
        "같은 토크를 내는 데 필요한 총 토크 용량이 직접 최적화가 작거나 같다 → 설계는 직접 최적화, 시너지는 비교 기준",
        "VAF 90 % 로 고른 시너지 수는 필요한 케이블 수와 다르다 (어깨: 시너지 6개로는 토크를 다 못 냄)",
        "사람 어깨(MoBL)는 팔 -22° 에서도 앞삼각근 +4.1 cm 굴곡 모멘트암이 남는다 → 부착점 시작값으로 쓸 여지",
    ], fs=13.5, dy=0.07)
    pdf.savefig(fig); plt.close(fig)

    # 9 전신
    fig = slide("전신 걷기 (최종, 한쪽 20개)", "부위별로 푼 걷기 해를 한 모델로 합침 (부위끼리 같은 관절을 당기는 케이블이 없어 한 번에 푸는 것과 같은 문제)")
    try:
        img = mpimg.imread(ROOT + "report/frame_upper_robot40h160_wincombo_walk_s3_T047_12.png")
        ax = fig.add_axes([0.04, 0.08, 0.92, 0.73]); ax.imshow(img); ax.axis("off")
    except Exception:
        pass
    fig.text(0.5, 0.055, "영상: report/video_upper_robot40h160_wincombo_walk_s3_T047_12.mp4", fontsize=11, color=GREY, ha="center")
    pdf.savefig(fig); plt.close(fig)

    # 9-1 맥 vs 윈도우
    fig = slide("이전 설계(맥)와 최종 비교", "같은 걷기(Moore subject3 1.2 m/s)를 두 설계로 재현. 비교 영상: report/video_compare_mac_vs_win_walk_s3_T047_12.mp4")
    for k, (fn, lab) in enumerate([("report/frame_upper_robot40h160_final_walk_s3_T047_12.png", "이전 설계 (한쪽 22개)"), ("report/frame_upper_robot40h160_wincombo_walk_s3_T047_12.png", "최종 (한쪽 20개)")]):
        try:
            img = mpimg.imread(ROOT + fn); ax = fig.add_axes([0.02 + k * 0.32, 0.36, 0.31, 0.44]); ax.imshow(img); ax.axis("off")
            fig.text(0.175 + k * 0.32, 0.81, lab, fontsize=13, weight="bold", color=NAVY, ha="center")
        except Exception:
            pass
    table(fig, [0.66, 0.36, 0.32, 0.45], ["", "이전 (맥)", "최종"], [
        ["케이블 (양쪽)", "40", "36"],
        ["M1 (양쪽)", "23", "16"],
        ["모듈 무게 (park)", "5.2 kg", "5.0 kg"],
        ["하체", "사후 배정", "보정 모듈 먼저"],
        ["몸통", "사후 배정", "모듈 먼저 (M1 0)"],
        ["어깨", "8개", "6개 (신뢰 구간)"],
    ], col_w=[0.36, 0.30, 0.34], fs=11.5)
    bullets(fig, 0.04, 0.29, [
        "두 설계 모두 걷기 2시행·들기(스쿼트) 2시행을 10 % 미만으로 통과",
        "최종이 케이블 4개·M1 7개 적다. 차이는 경로 2단계를 신뢰 구간으로 제대로 돌린 것과 몸통·하체 모듈을 먼저 정한 것",
        ("맥 어깨 8개는 경로 2단계가 실제로 돌지 않은 설계였다 (힘 재계산으로 통과)", GREY),
    ], fs=12.5, dy=0.065)
    pdf.savefig(fig); plt.close(fig)

    # 10 제안·질문·다음
    fig = slide("다음 단계")
    fig.text(0.04, 0.80, "제안: 손목 최적화 추가", fontsize=18, weight="bold", color=ACC)
    fig.text(0.04, 0.755, "팀원들이 손목(회내·회외, 굽힘) 데이터셋을 찾아 줘서, 같은 파이프라인으로 손목만 최적화해 더 해 보려 한다", fontsize=12.5, color="#333")
    table(fig, [0.04, 0.47, 0.56, 0.25], ["데이터셋", "좋은 점", "한계"], DATASETS, col_w=[0.27, 0.43, 0.30], fs=11.5)
    fig.text(0.04, 0.36, "모듈 사양·무게 (이렇게 진행)", fontsize=15, weight="bold", color=NAVY)
    bullets(fig, 0.04, 0.31, [
        "힘·스트로크는 계획서 모듈표 범위 안에서 설계 (M1 ≤ 5,000 N·120 mm 등)",
        "무게는 park2019 SMA 묶음 기준 (10.8 N/g)",
    ], fs=12.5, dy=0.055)
    fig.text(0.64, 0.80, "질문", fontsize=18, weight="bold", color=ACC)
    bullets(fig, 0.64, 0.74, [
        "들기 하중 목표는 어느 정도인지",
        "지금 기준 찬 물병 3.79 kg 이면 충분한지",
        "강병현 교수님 뼈대 모델은 언제 받을 수 있는지",
        ("  (상체 크기는 들기가 정함: 어깨 81·몸통 155·팔꿈치 34 N·m)", GREY),
    ], fs=13, dy=0.06)
    fig.text(0.64, 0.43, "다음 단계", fontsize=18, weight="bold", color=NAVY)
    bullets(fig, 0.64, 0.37, [
        "손목 최적화 (위 제안)",
        "MuJoCo 에 160 cm 파일 연결 → RL 걷기",
        "경로 개선: 신뢰 구간 2단계, 시드 여러 개",
        "하체 무릎·발목 감는 면, 왼쪽 검증",
    ], fs=12.5, dy=0.06)
    pdf.savefig(fig); plt.close(fig)
print("saved", OUT, page_no[0], "pages")
