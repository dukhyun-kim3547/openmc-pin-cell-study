# =============================================================================
# 농축도 민감도 분석 (Enrichment Sensitivity Study)
# 목적: U-235 농축도 변화에 따른 k-infinity 변화 계산
# =============================================================================

import openmc
import numpy as np
import matplotlib
matplotlib.use('Agg')  # GUI 없는 환경(WSL)에서 그래프 저장용
import matplotlib.pyplot as plt

# =============================================================================
# 분석할 농축도 목록 (wt%)
# 2%: 저농축 / 3.1%: 일반 PWR / 4.5%: 고연소도 / 5%: 법적 상한선
# 20%: HALEU (차세대 SMR용, 참고용)
# =============================================================================
enrichments = [1.0, 2.0, 3.1, 4.0, 5.0, 10.0, 19.75]

keff_values  = []  # k-effective 결과 저장
keff_errors  = []  # 불확도 저장

# =============================================================================
# 농축도별 반복 계산
# =============================================================================
for enr in enrichments:
    print(f"\n{'='*50}")
    print(f"  농축도 {enr:.2f} wt% 계산 중...")
    print(f"{'='*50}")

    # --- 재료 정의 ---
    uo2 = openmc.Material(name=f'UO2_{enr}pct')
    uo2.add_nuclide('U235', enr,         percent_type='wo')  # 농축도 변화
    uo2.add_nuclide('U238', 100-enr-13.5, percent_type='wo') # U238 (나머지)
    uo2.add_element('O',    13.5,         percent_type='wo') # 산소 고정
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
    water.set_density('g/cm3', 0.71)
    water.add_s_alpha_beta('c_H_in_H2O')

    materials = openmc.Materials([uo2, zircaloy, water])
    materials.export_to_xml()

    # --- 기하학 정의 (동일) ---
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
    uniform_dist   = openmc.stats.Box(bounds[:3], bounds[3:])
    settings.source = openmc.IndependentSource(space=uniform_dist)
    settings.export_to_xml()

    # --- 실행 ---
    openmc.run(output=False)  # output=False: 로그 숨김 (깔끔한 출력)

    # --- 결과 읽기 ---
    sp = openmc.StatePoint('statepoint.150.h5')
    keff = sp.keff
    keff_values.append(keff.nominal_value)
    keff_errors.append(keff.std_dev)
    sp.close()

    print(f"  ✅ k-infinity = {keff.nominal_value:.5f} ± {keff.std_dev:.5f}")

# =============================================================================
# 결과 출력 및 그래프
# =============================================================================
print(f"\n{'='*50}")
print("  최종 결과 요약")
print(f"{'='*50}")
print(f"{'농축도(%)':>10} | {'k-infinity':>12} | {'불확도':>10}")
print("-" * 40)
for enr, k, err in zip(enrichments, keff_values, keff_errors):
    marker = " ← 현재 PWR" if enr == 3.1 else ""
    marker = " ← SMR(HALEU)" if enr == 19.75 else marker
    print(f"{enr:>10.2f} | {k:>12.5f} | ±{err:.5f}{marker}")

# --- 그래프 저장 ---
plt.figure(figsize=(8, 5))
plt.errorbar(enrichments, keff_values, yerr=keff_errors,
             fmt='o-', color='steelblue', capsize=5, linewidth=2, markersize=8)
plt.axhline(y=1.0, color='red', linestyle='--', linewidth=1.5, label='임계 (k=1.0)')
plt.axvline(x=5.0, color='orange', linestyle=':', linewidth=1.5, label='민수용 농축 상한 (5%)')
plt.xlabel('U-235 농축도 (wt%)', fontsize=12)
plt.ylabel('k-infinity', fontsize=12)
plt.title('PWR 핀 셀: 농축도에 따른 k-infinity 변화', fontsize=13)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('enrichment_sensitivity.png', dpi=150)
print("\n📊 그래프 저장됨: enrichment_sensitivity.png")
