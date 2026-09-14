"""현황 보고 PPT 생성. 9/8 48-trial 요약 + subject10 결과를 병합해 7명 42 trial로 정리한다.

    conda activate osim
    python report/make_ppt.py
"""
import csv, os, re, statistics as st
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "report")
KEEP = ("subject3", "subject4", "subject5", "subject7", "subject8", "subject10", "subject11")
EXCLUDED = {"subject2": "manualReview + AddB 잔차 5~13 %", "subject6": "AddB 잔차 13~18 %BW",
            "subject9": "AddB 잔차 16~30 %BW"}
FONT = "Apple SD Gothic Neo"
for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f or "AppleGothic" in f:
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False

NAVY, GREY, RED, GREEN = RGBColor(0x1F, 0x3A, 0x5F), RGBColor(0x55, 0x55, 0x55), RGBColor(0xC0, 0x39, 0x2B), RGBColor(0x2E, 0x7D, 0x32)


def num(v):
    try: return float(v)
    except (TypeError, ValueError): return None


def cond(t): return "walkingTS" if "TS" in t else "walking"


# ── 1. 병합 ────────────────────────────────────────────────────────────────
def load_rows():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "out/summary_2026-09-08_48trial.csv"))))
    s10 = os.path.join(ROOT, "out/real/summary.csv")
    if os.path.exists(s10):
        rows += [r for r in csv.DictReader(open(s10)) if r["subject"] == "subject10"]
    keep = [r for r in rows if r["subject"] in KEEP]
    fields = list(rows[0].keys())
    with open(os.path.join(ROOT, "out/summary_7subjects.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(keep)
    return rows, keep


# ── 2. AddB 자체 잔차 (inspect 로그) ───────────────────────────────────────
def addb_residuals():
    out, sub, mass = {}, None, None
    for line in open(os.path.join(ROOT, "out/inspect_uhlrich2023.txt")):
        m = re.match(r"=== (subject\d+) ===", line)
        if m: sub = m.group(1); continue
        m = re.search(r"질량 ([\d.]+) kg", line)
        if m: mass = float(m.group(1)); continue
        m = re.match(r"\s+(walking\S+)\s+\d+\s+\d+\s+평균([\d.]+) 최대([\d.]+)", line)
        if m and sub:
            out.setdefault(sub, []).append(float(m.group(2)) / (mass * 9.80665) * 100)
    return out


# ── 3. 그림 ────────────────────────────────────────────────────────────────
def fig_k_hist(keep):
    g = [int(num(r[c])) for r in keep for c in ("syn_kglobal_r", "syn_kglobal_l") if num(r[c])]
    cnt = Counter(g)
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ks = sorted(cnt); ax.bar(ks, [cnt[k] for k in ks], color="#1F3A5F", width=0.6)
    for k in ks: ax.text(k, cnt[k] + 1, str(cnt[k]), ha="center", fontsize=11)
    ax.set_xlabel("시너지 개수 k (global: 전체 VAF ≥ 90 %)"); ax.set_ylabel("분해 수 (trial × 좌우)")
    ax.set_xticks(ks); ax.set_title(f"{len(g)}개 분해 · 평균 {st.mean(g):.2f} · 중앙값 {st.median(g):.0f}", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
    p = os.path.join(OUT, "fig_k_hist.png"); fig.savefig(p, dpi=200); plt.close(fig); return p, g


def fig_condition(keep):
    subs = sorted({r["subject"] for r in keep}, key=lambda s: int(s[7:]))
    a, b = [], []
    for s in subs:
        for c, dst in (("walking", a), ("walkingTS", b)):
            v = [num(r[k]) for r in keep if r["subject"] == s and cond(r["trial"]) == c
                 for k in ("syn_kglobal_r", "syn_kglobal_l") if num(r[k])]
            dst.append(st.mean(v) if v else np.nan)
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    x = np.arange(len(subs))
    for i in range(len(subs)): ax.plot([0, 1], [a[i], b[i]], color="#999", lw=1, zorder=1)
    ax.scatter([0] * len(subs), a, s=60, color="#1F3A5F", zorder=2, label="walking")
    ax.scatter([1] * len(subs), b, s=60, color="#C0392B", zorder=2, label="walkingTS")
    ax.plot([0, 1], [np.nanmean(a), np.nanmean(b)], color="k", lw=2.5, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels([f"walking\n평균 {np.nanmean(a):.2f}", f"walkingTS\n평균 {np.nanmean(b):.2f}"])
    ax.set_ylabel("subject 평균 global k"); ax.set_xlim(-0.4, 1.4); ax.set_ylim(3.5, 7)
    ax.set_title(f"조건 차이 {np.nanmean(b)-np.nanmean(a):+.2f}  (n={len(subs)})", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
    p = os.path.join(OUT, "fig_condition.png"); fig.savefig(p, dpi=200); plt.close(fig); return p, a, b, subs


def fig_residual(res):
    subs = sorted(res, key=lambda s: int(s[7:]))
    means = [st.mean(res[s]) for s in subs]; maxs = [max(res[s]) for s in subs]
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    col = ["#C0392B" if s in EXCLUDED else ("#E67E22" if s == "subject10" else "#1F3A5F") for s in subs]
    ax.bar(range(len(subs)), means, color=col, width=0.6)
    ax.scatter(range(len(subs)), maxs, color="k", s=14, zorder=3, label="최악 trial")
    ax.axhline(2.0, color="#888", ls="--", lw=1); ax.text(len(subs) - 0.5, 2.3, "2 %BW", ha="right", fontsize=9, color="#555")
    ax.set_yscale("log"); ax.set_xticks(range(len(subs))); ax.set_xticklabels([s[7:] for s in subs])
    ax.set_xlabel("subject"); ax.set_ylabel("AddBiomechanics 자체 잔차\n(walking trial 평균, %BW, log)")
    ax.legend(loc="upper left", fontsize=9, frameon=False)
    ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
    p = os.path.join(OUT, "fig_residual.png"); fig.savefig(p, dpi=200); plt.close(fig); return p, dict(zip(subs, means))


def fig_synergy_w():
    """subject3 walking1 오른쪽, k=5 W 히트맵."""
    import sys; sys.path.insert(0, ROOT)
    import synergy as sy
    act = os.path.join(ROOT, "out/real/subject3/walking1_segment_0/so/walking1_segment_0_StaticOptimization_activation.sto")
    mus = sy.muscle_names(os.path.join(ROOT, "data/subject3/subject3_scaled.osim"))
    A, names, t, _ = sy.load_activations(act, "r", muscles=mus, log=lambda *a: None)
    C, W = sy.nmf(A, 5, seed=0); C, W, contrib = sy.sort_synergies(C, W)
    Wn = W / W.max(axis=1, keepdims=True)
    order = np.argsort(-Wn.argmax(axis=0) * 100 - Wn.max(axis=0))  # 소속 시너지별 정렬
    order = sorted(range(len(names)), key=lambda j: (int(np.argmax(Wn[:, j])), -Wn[:, j].max()))
    fig, ax = plt.subplots(figsize=(11, 2.6))
    ax.imshow(Wn[:, order], aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_yticks(range(5)); ax.set_yticklabels([f"S{i+1} ({contrib[i]*100:.0f}%)" for i in range(5)], fontsize=9)
    ax.set_xticks(range(len(names))); ax.set_xticklabels([names[j].rsplit("_", 1)[0] for j in order], rotation=90, fontsize=7)
    ax.set_title("subject3 walking1 (오른쪽) k=5 — 시너지 가중치 W (행별 최대 1로 정규화)", fontsize=10)
    fig.tight_layout(); p = os.path.join(OUT, "fig_synergy_w.png"); fig.savefig(p, dpi=200); plt.close(fig); return p


def synergy_torque_rows(sub="subject3", trial="walking1_segment_0", side="r"):
    """시너지별 관절 토크 서명 — 새 근육 i 가 각 관절에 내야 할 피크 Nm.
    근육 j 의 SO 힘을 C_i W_ij / Σ_k C_k W_kj 비율로 시너지에 귀속시켜 모멘트암을 곱한다."""
    import sys; sys.path.insert(0, ROOT)
    import synergy as sy, addb_io as io, run_muscle_analysis as rma
    base = os.path.join(ROOT, f"out/real/{sub}/{trial}")
    frc = io.read_mot(f"{base}/so/{trial}_StaticOptimization_force.sto")
    idm = io.read_mot(f"{base}/{trial}_id.sto")
    mus = sy.muscle_names(os.path.join(ROOT, f"data/{sub}/{sub}_scaled.osim"))
    coords = [f"hip_flexion_{side}", f"hip_adduction_{side}", f"knee_angle_{side}", f"ankle_angle_{side}"]
    ma = rma.load_moment_arms(f"{base}/ma", trial, coords)
    A, names, _, _ = sy.load_activations(f"{base}/so/{trial}_StaticOptimization_activation.sto", side, muscles=mus, log=lambda *a: None)
    C, W = sy.nmf(A, 5, seed=0); C, W, contrib = sy.sort_synergies(C, W)
    share = C[:, :, None] * W[None, :, :]; share /= np.maximum(share.sum(axis=1, keepdims=True), 1e-12)
    F = np.column_stack([frc.column(n) for n in names]); t = frc.time
    R = {c: np.column_stack([np.interp(t, ma[c].time, ma[c].column(n)) for n in names]) for c in coords}
    rows = [["시너지", "고관절 굴곡(+)", "고관절 내전(+)", "무릎 굴곡(+)", "발목 배굴(+)", "구성근 F_max 합"]]
    for i in range(5):
        r = [f"S{i+1} ({contrib[i]*100:.0f}%)"]
        for c in coords:
            T = (share[:, i, :] * F * R[c]).sum(axis=1); r.append(f"{T[np.argmax(np.abs(T))]:+.0f} Nm")
        r.append(f"{(W[i] / W[i].max() * F.max(axis=0)).sum():,.0f} N")
        rows.append(r)
    idpk = [np.abs(np.interp(t, idm.time, idm.column(c + "_moment"))).max() for c in coords]
    rows.append(["ID 모멘트 피크"] + [f"{v:.0f} Nm" for v in idpk] + ["(5개 합 = ID)"])
    return rows


# ── 4. PPT ─────────────────────────────────────────────────────────────────
def txt(slide, x, y, w, h, text, size=14, bold=False, color=GREY, align=PP_ALIGN.LEFT):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)); tf = tb.text_frame; tf.word_wrap = True
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run(); r.text = line; r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = color; r.font.name = FONT
        p.space_after = Pt(4)
    return tb


def bullets(slide, x, y, w, h, items, size=13):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)); tf = tb.text_frame; tf.word_wrap = True
    for i, it in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        lvl = 1 if it.startswith("  ") else 0
        r = p.add_run(); r.text = ("• " if lvl == 0 else "– ") + it.strip(); r.font.size = Pt(size - 2 * lvl)
        r.font.color.rgb = GREY; r.font.name = FONT; p.level = lvl; p.space_after = Pt(5)
    return tb


def table(slide, x, y, w, rows, col_w=None, size=10, header_fill=NAVY):
    nr, nc = len(rows), len(rows[0])
    shp = slide.shapes.add_table(nr, nc, Inches(x), Inches(y), Inches(w), Inches(0.28 * nr)); t = shp.table
    if col_w:
        for j, cw in enumerate(col_w): t.columns[j].width = Inches(cw)
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            c = t.cell(i, j); c.text = ""; p = c.text_frame.paragraphs[0]; r = p.add_run(); r.text = str(val)
            r.font.size = Pt(size); r.font.name = FONT; r.font.bold = i == 0
            r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF) if i == 0 else GREY
            c.margin_top = c.margin_bottom = Emu(20000)
            if i == 0: c.fill.solid(); c.fill.fore_color.rgb = header_fill
            elif i % 2 == 0: c.fill.solid(); c.fill.fore_color.rgb = RGBColor(0xF2, 0xF4, 0xF7)
            else: c.fill.solid(); c.fill.fore_color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
    return shp


def header(slide, title, n):
    bar = slide.shapes.add_shape(1, 0, 0, Inches(13.333), Inches(0.9)); bar.fill.solid(); bar.fill.fore_color.rgb = NAVY; bar.line.fill.background()
    txt(slide, 0.5, 0.15, 11.5, 0.6, title, size=24, bold=True, color=RGBColor(0xFF, 0xFF, 0xFF))
    txt(slide, 12.2, 0.25, 0.9, 0.5, f"{n} / 5", size=12, color=RGBColor(0xDD, 0xDD, 0xDD), align=PP_ALIGN.RIGHT)


def build(keep, g, k_hist, k_cond, cond_a, cond_b, subs, resid_png, resid_mean, w_png, torque_rows):
    prs = Presentation(); prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    blank = prs.slide_layouts[6]

    # 1 ──────────────────────────────────────────
    s = prs.slides.add_slide(blank); header(s, "보행 근육 시너지 추출 파이프라인 — 현황 (2026-09-14)", 1)
    txt(s, 0.5, 1.1, 12, 0.5, "목표: 마커·힘판 데이터에서 관절 모멘트 → 근육별 힘 → 근육 시너지를 자동 계산하고, 인공근육 휴머노이드의 액추에이터 배치 근거로 쓴다", size=15, bold=True, color=NAVY)
    steps = ["마커·힘판", "IK\n관절 각도", "ID\n관절 모멘트", "SO\n근육별 힘", "NMF\n근육 시너지"]
    for i, sname in enumerate(steps):
        x = 0.7 + i * 2.5
        b = s.shapes.add_shape(1, Inches(x), Inches(1.9), Inches(2.0), Inches(0.9)); b.fill.solid()
        b.fill.fore_color.rgb = NAVY if i else GREY; b.line.fill.background()
        tf = b.text_frame; tf.paragraphs[0].alignment = PP_ALIGN.CENTER
        for k, line in enumerate(sname.split("\n")):
            p = tf.paragraphs[0] if k == 0 else tf.add_paragraph(); p.alignment = PP_ALIGN.CENTER
            r = p.add_run(); r.text = line; r.font.size = Pt(13 if k == 0 else 11); r.font.bold = k == 0; r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF); r.font.name = FONT
        if i < len(steps) - 1: txt(s, x + 2.0, 2.1, 0.5, 0.5, "→", size=22, bold=True, color=NAVY, align=PP_ALIGN.CENTER)
    txt(s, 5.7, 2.9, 5, 0.4, "MuscleAnalysis(모멘트암) — ID·SO 검산에 사용", size=11, color=GREY)
    bullets(s, 0.5, 3.5, 6.2, 3.5, [
        "데이터: AddBiomechanics 공개 데이터셋 Uhlrich2023 (With_Arm)",
        "  무릎 부하를 줄이는 보행 수정 중재 — 같은 사람이 평상(walking)·수정(walkingTS) 둘 다",
        "  Rajagopal 모델 37 DOF, 근육 80개 → SO·시너지까지 가능",
        "  스케일된 모델과 IK가 함께 배포되므로 Scale/IK는 직접 하지 않음",
        "각 단계마다 물리 검산 (PASS / WARN / FAIL)",
        "  계산은 입력이 틀려도 숫자를 내놓는다 — 검산 없이는 조용히 틀린 채로 간다",
    ])
    bullets(s, 7.0, 3.5, 6, 3.5, [
        "환경: macOS(Apple Silicon) 한 대에서 변환 + 분석",
        "  예전 Windows + Linux 이원화 폐기 (nimblephysics에 Windows 휠이 없어서였음)",
        "  conda 환경 2개: osim (OpenSim 4.6) / b3d (nimblephysics)",
        "작업물 소실 후 README·NOTES만으로 재구축",
        "  드라이런(번들 Rajagopal 보행) 전 단계 재현 완료 → 2쪽",
        "  9/8 실제 데이터 48 trial 결과는 남아 있어 그대로 활용",
    ])

    # 2 ──────────────────────────────────────────
    s = prs.slides.add_slide(blank); header(s, "검산 — '조용히 틀리는' 지점을 어떻게 막았나", 2)
    txt(s, 0.5, 1.05, 12, 0.4, "에러 없이 그럴듯한 그래프가 나오는 함정이 네 곳 있다. 전부 이름표 대신 물리로 판정한다.", size=13, color=GREY)
    table(s, 0.5, 1.55, 7.6, [
        ["함정", "방법", "드라이런 결과"],
        ["GRF를 어느 발에 붙이나", "컬럼 이름 무시. IK 자세에서 COP–발 거리로 배정", "r→calcn_r  0.137 vs 0.483 m (여유 34 cm)"],
        ["어느 시간 구간을 쓰나", "IK∩GRF ∩ 지지구간(≥0.7 BW) + 필터 마진(1/6 Hz)", "잔차 93.6 % → 41.9 % → 7.8 %"],
        ["관절 모멘트 자릿수", "밴드: 발목 1.2~1.6, 무릎 0.4~0.8, 고관절 0.7~1.2 Nm/kg", "전 관절 자릿수 일치"],
        ["SO가 진짜 근육으로 풀었나", "보조 액추에이터 개입 비율 (reserve) < 10 %", "0.19 % — 근육 힘 그대로 해석 가능"],
        ["GRF 파일이 맞나", "힘 xyz + 작용점 xyz 둘 다 있어야 GRF", "자유모멘트(_mx) 오분류 차단"],
    ], col_w=[1.9, 3.3, 2.4], size=9.5)
    txt(s, 8.4, 1.5, 4.6, 0.4, "재현 확인 — Windows(예전) vs macOS(지금)", size=12, bold=True, color=NAVY)
    table(s, 8.4, 1.9, 4.6, [
        ["항목", "예전", "지금"],
        ["발목 r/l (Nm/kg)", "1.644 / 1.680", "1.643 / 1.681"],
        ["무릎 r/l", "0.670 / 0.700", "0.670 / 0.701"],
        ["고관절 r/l", "0.726 / 0.703", "0.726 / 0.702"],
        ["잔차", "7.82 %", "7.81 %"],
        ["분석 구간", "0.539~1.892 s", "0.539~1.892 s"],
        ["보조 액추에이터", "29개", "29개"],
        ["SO recon / reserve", "0.07 / 0.19 %", "0.072 / 0.188 %"],
    ], col_w=[1.8, 1.4, 1.4], size=9.5)
    bullets(s, 0.5, 4.7, 12.3, 2.5, [
        "Scale·IK를 새로 돌렸는데도 관절 모멘트가 소수 셋째 자리까지 같다 → 파이프라인 회귀 없음",
        "변환기 출력 형식(라디안, ground_force_<body>_v/p/m)으로 드라이런을 다시 써서 넣어도 summary 전 항목 동일 → 변환·매핑 검증",
        "재구축 중 새로 잡은 함정: 활성도 파일에 팔 토크 액추에이터(shoulder_flex 등)가 근육처럼 섞여 있었다 → 모델 근육 목록으로 걸러냄 (41 → 38개)",
    ], size=12)

    # 3 ──────────────────────────────────────────
    s = prs.slides.add_slide(blank); header(s, "데이터 품질 — 이름표가 아니라 물리로 subject를 고른다", 3)
    s.shapes.add_picture(resid_png, Inches(0.4), Inches(1.1), width=Inches(6.4))
    rows = [["subject", "AddB 자체 잔차\n(%BW, 평균)", "우리 잔차\n(% peak GRF)", "판정"]]
    ours = {}
    for sname in sorted(resid_mean, key=lambda x: int(x[7:])):
        v = [num(r["residual"]) for r in keep if r["subject"] == sname if num(r["residual"])]
        ours[sname] = f"{st.mean(v):.1f}" if v else "—"
    for sname in sorted(resid_mean, key=lambda x: int(x[7:])):
        verdict = EXCLUDED.get(sname, "사용")
        if sname == "subject10": verdict = "사용 (manualReview 복구)"
        rows.append([sname[7:], f"{resid_mean[sname]:.1f}", ours[sname], verdict])
    table(s, 7.0, 1.1, 6.0, rows, col_w=[0.8, 1.4, 1.3, 2.5], size=9.5)
    bullets(s, 0.4, 4.6, 12.5, 2.7, [
        "AddBiomechanics가 스스로 보고한 골반 잔차(--inspect)를 체중으로 정규화 — 양호 subject는 전부 1.3 %BW 미만, subject6·9는 13~30 %BW (20~40배)",
        "  subject6·9의 높은 잔차는 우리 계산 문제가 아니라 원 데이터의 동역학 적합 실패 → 제외. 우리 검산기가 이 둘만 골라낸 것이 정상 작동의 증거",
        "subject10은 manualReview(사람이 검토 안 함) 플래그로 자동 제외됐었지만 물리 품질은 subject3보다 좋다 → 플래그를 열고 우리 검산기로 판정해 복구",
        "  subject2는 같은 플래그인데 잔차도 실제로 나쁘다 → 제외. 두 사람을 같은 이유로 묶어 제외했던 것이 잘못",
        f"결과: 10명 → 7명 × 6 trial = {len(keep)} trial  (walking 3 + walkingTS 3 / 명)",
    ], size=12)

    # 4 ──────────────────────────────────────────
    s = prs.slides.add_slide(blank); header(s, "결과 — 시너지는 5개, 보행 수정으로 개수는 바뀌지 않는다", 4)
    s.shapes.add_picture(k_hist, Inches(0.4), Inches(1.05), width=Inches(5.0))
    s.shapes.add_picture(k_cond, Inches(5.5), Inches(1.05), width=Inches(5.0))
    cnt = Counter(g)
    bullets(s, 10.6, 1.1, 2.6, 3.3, [
        f"{len(g)}개 분해 (trial × 좌우)",
        f"k=5가 {cnt[5]}개 ({cnt[5]/len(g)*100:.0f} %)",
        f"범위 {min(g)}~{max(g)}, 평균 {st.mean(g):.2f}",
        "문헌 보행 시너지 4~6개와 일치",
        f"walking {np.nanmean(cond_a):.2f} vs TS {np.nanmean(cond_b):.2f}",
        "subject별 차이 −0.5~+0.3, 방향 뒤섞임",
    ], size=11)
    s.shapes.add_picture(w_png, Inches(0.4), Inches(4.4), width=Inches(8.3))
    table(s, 8.9, 4.45, 4.2, [
        ["시너지 (k=5, subject 공통)", "묶이는 근육"],
        ["고관절 외전·골반 안정", "glmed1-3, glmin1-3, tfl"],
        ["push-off + 고관절 굴곡", "soleus, psoas, iliacus, fhl, perlong, tibpost"],
        ["고관절 신전", "glmax1-3, bflh, semimem, semiten, addmagIsch"],
        ["배굴 (유각·초기접지)", "tibant, edl, ehl (± vasti)"],
        ["내전", "addlong, addbrev, addmagProx/Mid, grac, piri"],
    ], col_w=[1.7, 2.5], size=9)
    txt(s, 8.9, 6.4, 4.2, 0.8, "개수로는 조건 차이를 말할 수 없다 (n=7, k는 정수).\n바뀔 수 있는 것은 구성(W)과 타이밍(C) → 5쪽", size=11, bold=True, color=RED)

    # 5 ──────────────────────────────────────────
    s = prs.slides.add_slide(blank); header(s, "해석 방법과 다음 단계", 5)
    txt(s, 0.5, 1.05, 6.2, 0.4, "Q. 예전엔 8개가 안정적이라 했는데 5개?", size=14, bold=True, color=NAVY)
    table(s, 0.5, 1.5, 6.2, [
        ["기준", "조건", "실제 데이터", "뜻"],
        ["global", "전체 VAF ≥ 90 %", "4~7, 평균 5.2", "제어 구조의 차원"],
        ["dual", "+ 모든 근육 개별 VAF ≥ 75 %", "5~12, 평균 8.0", "SO 임의 배분까지 재현"],
    ], col_w=[0.9, 2.1, 1.4, 1.8], size=9.5)
    bullets(s, 0.5, 2.5, 6.2, 2.4, [
        "같은 데이터의 두 기준. '8'은 dual, '5'는 global",
        "추가 3개는 gaslat/gasmed, glmed1/2/3처럼 모멘트암이 겹치는 근육 사이에 SO가 임의로 나눈 부하를 재현하는 데 쓰인다 — 근육 하나짜리 성분이 나오는 것이 증거",
        "trial을 이어붙여 표본을 4배로 늘려도 dual은 안 내려간다 → 표본 부족이 아니라 SO의 성질",
        "본 결과는 global(5). dual은 '이 근육은 개별 해석 금지' 진단 지표로 남긴다",
        "하드웨어 설계에 8을 쓰면 SO의 동전 던지기를 기계에 새기게 된다",
    ], size=11.5)
    txt(s, 0.5, 4.75, 6.2, 0.4, "새 근육 5개의 스펙 — 시너지별 관절 토크 서명 (subject3 walking1, 오른쪽)", size=12, bold=True, color=NAVY)
    table(s, 0.5, 5.15, 6.2, torque_rows, col_w=[1.0, 1.05, 1.05, 1.0, 1.0, 1.1], size=9)
    txt(s, 0.5, 7.0, 6.2, 0.4, "S2·S4가 둘 다 발목 저굴 −66/−72 Nm — soleus/gasmed를 SO가 갈라놓은 흔적. 묶기(②) 없이는 하드웨어에 그대로 새겨진다", size=9.5, color=RED)
    txt(s, 7.0, 1.05, 6, 0.4, "다음 단계", size=14, bold=True, color=NAVY)
    bullets(s, 7.0, 1.5, 6.1, 5.8, [
        "① subject10 결과 병합 → 7명 42 trial 확정 (완료)",
        "② 중복 협응근 묶기를 synergy.py에 추가 (gasmed+gaslat, glmed1-3, glmax1-3, addmag*)",
        "  묶은 뒤 dual이 global 쪽으로 내려오면 → SO 임의성이 원인이라는 증명이자 W 비교의 전제",
        "③ walking vs walkingTS를 W(구성)·C(타이밍)로 비교 — 근육 집합을 교집합으로 맞춘 뒤",
        "④ 시너지당 새 근육 1개(다리당 5개)로 걷기 — 그대로는 안 된다",
        "  이론: 당김 전용 액추에이터가 6 DOF 토크를 양의 결합으로 덮으려면 다리당 ≥7 (Phase 0). 5개는 subtalar·hip_rotation 잠근 4 DOF에서만 n+1",
        "  실측: 5개 시너지 토크 방향벡터가 4 DOF를 양으로 덮는지 12개 분해로 검사 → 1/12만 성립. 걷기 방향만 있고 균형·외란 대응 방향은 못 냄",
        "  되는 구조: 하드웨어 근육 K≥7~8 (스패닝+스트로크, Phase 0~5 계획) + 제어 신호 5개 (SynergyController에 W). '8'은 근육 수, '5'는 제어 차원",
        "  걷기 재현만 목적이면 시뮬레이션은 가능: SynergyController(80근육·5신호) → 5개 lumped 액추에이터 + Moco 추적 (설계검증 흐름도 'CMC 추적 검증')",
        "  아래 표가 그 출발 스펙. S2·S4 저굴 중복(soleus/gasmed 갈림)은 묶기(②) 후 다시 뽑는다",
    ], size=11)

    p = os.path.join(OUT, "보행시너지_현황_2026-09-14.pptx"); prs.save(p); return p


if __name__ == "__main__":
    rows, keep = load_rows()
    print(f"병합: 전체 {len(rows)} → 7명 {len(keep)} trial")
    res = addb_residuals()
    k_hist, g = fig_k_hist(keep)
    k_cond, a, b, subs = fig_condition(keep)
    resid_png, resid_mean = fig_residual(res)
    w_png = fig_synergy_w()
    torque_rows = synergy_torque_rows()
    print("PPT:", build(keep, g, k_hist, k_cond, a, b, subs, resid_png, resid_mean, w_png, torque_rows))
