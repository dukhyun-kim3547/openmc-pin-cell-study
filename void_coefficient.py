# =============================================================================
# 보이드 계수 분석 (Void Coefficient Analysis)
# 목적: 냉각수 밀도 감소에 따른 k-infinity 변화로 PWR 부의 보이드 계수 확인
# 물리적 의미: 사고 시 냉각수 손실 → 자동 출력 억제 → 고유 안전성(Inherent Safety)
# =============================================================================

import openmc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =============================================================================
# 분석할 냉각수 밀도 목록 (g/cm3)
# 0.71: PWR 정상 운전 (325°C, 155 bar)
# 0.50: 부분 보이드
# 0.30: 심각한 보이드
# 0.10: 거의 증기 상태
# 0.001: 완전 보이드 (증기만 존재)
# =============================================================================
water_densities = [0.75, 0.71, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05, 0.001]
void_fractions  = [(0.75 - d) / 0.75 * 100 for d in water_densities]  # 보이드율 (%)

keff_values = []
keff_errors = []

# =============================================================================
# 밀도별 반복 계산
# =============================================================================
for density in water_densities:
    void_pct = (0.75 - density) / 0.75 * 100
    print(f"\n{'='*50}")
    print(f"  냉각수 밀도 {density:.3f} g/cc (보이드율 {void_pct:.1f}%) 계산 중...")
    print(f"{'='*50}")

    # --- 재료 정의 (농축도 3.1% 고정) ---
    uo2 = openmc.Material(name='UO2')
    uo2.add_nuclide('U235', 3.1,          percent_type='wo')
    uo2.add_nuclide('U238', 100-3.1-13.5, percent_type='wo')
    uo2.add_element('O',    13.5,         percent_type='wo')
    uo2.set_density('g/cm3', 10.97)

    zircaloy = openmc.Material(name='Zircaloy-4')
    zircaloy.add_element('Zr', 98.0, percent_type='wo')
    zircaloy.add_element('Sn',  1.5, percent_type='wo')
    zircaloy.add_element('Fe',  0.2, percent_type='wo')
    zircaloy.add_element('Cr',  0.1, percent_type='wo')
    zircaloy.set_density('g/cm3', 6.56)

    water = openmc.Material(name='Water')
    water.add_nuclide('H1',  2.0, percent_type='ao')
    water.add_nuclide('O16', 1.0, percent_type='ao')
    water.set_density('g/cm3', density)       # ← 밀도 변화
    if density > 0.01:                         # 거의 증기 상태에선 S(a,b) 불필요
        water.add_s_alpha_beta('c_H_in_H2O')

    materials = openmc.Materials([uo2, zircaloy, water])
    materials.export_to_xml()

    # --- 기하학 (동일) ---
    fuel_or = openmc.ZCylinder(r=0.4096)
    clad_or = openmc.ZCylinder(r=0.4750)
    pitch   = 1.26
    left   = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
    right  = openmc.XPlane(x0=+pitch/2, boundary_type='reflective')
    bottom = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
    top    = openmc.YPlane(y0=+pitch/2, boundary_type='reflective')

    fuel_cell  = openmc.Cell(fill=uo2,      region=-fuel_or)
    clad_cell  = openmc.Cell(fill=zircaloy, region=+fuel_or & -clad_or)
    water_cell = openmc.Cell(fill=water,    region=+clad_or & +left & -right & +bottom & -top)

    universe = openmc.Universe(cells=[fuel_cell, clad_cell, water_cell])
    geometry = openmc.Geometry(universe)
    geometry.export_to_xml()

    # --- 계산 설정 ---
    settings = openmc.Settings()
    settings.run_mode  = 'eigenvalue'
    settings.batches   = 150
    settings.inactive  = 50
    settings.particles = 5000

    bounds = [-0.4096, -0.4096, -1, 0.4096, 0.4096, 1]
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(bounds[:3], bounds[3:])
    )
    settings.export_to_xml()

    # --- 실행 ---
    openmc.run(output=False)

    # --- 결과 읽기 ---
    sp = openmc.StatePoint('statepoint.150.h5')
    keff = sp.keff
    keff_values.append(keff.nominal_value)
    keff_errors.append(keff.std_dev)
    sp.close()

    print(f"  k-infinity = {keff.nominal_value:.5f} ± {keff.std_dev:.5f}")

# =============================================================================
# 보이드 계수 계산
# α_void = Δk / Δ(void fraction)  [단위: pcm/%void = 10^-5 / %void]
# 부의 값 → 안전 (보이드 증가 시 반응도 감소)
# =============================================================================
print(f"\n{'='*55}")
print("  최종 결과 요약")
print(f"{'='*55}")
print(f"{'밀도(g/cc)':>10} | {'보이드율(%)':>10} | {'k-infinity':>12} | {'불확도':>8}")
print("-" * 55)
for d, v, k, e in zip(water_densities, void_fractions, keff_values, keff_errors):
    marker = " ← 정상운전" if d == 0.71 else ""
    print(f"{d:>10.3f} | {v:>10.1f} | {k:>12.5f} | ±{e:.5f}{marker}")

# 보이드 계수 (정상운전 → 완전보이드)
k_nominal = keff_values[1]   # 0.71 g/cc
k_void    = keff_values[-1]  # 0.001 g/cc
delta_k   = k_void - k_nominal
delta_void = void_fra
