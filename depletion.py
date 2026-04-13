# =============================================================================
# 연소 계산 (Depletion Simulation)
# 목적: PWR 핀 셀에서 연료 연소에 따른 핵종 조성 변화 및 k-infinity 추적
# 물리: U-235 감소, Pu-239 생성, 핵분열 생성물 축적
# =============================================================================

import openmc
import openmc.deplete
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math

# =============================================================================
# PART 1: 재료 정의
# =============================================================================

pitch  = 1.26
fuel_r = 0.4096

uo2 = openmc.Material(name='UO2 Fuel')
uo2.add_nuclide('U235', 3.1,           percent_type='wo')
uo2.add_nuclide('U238', 100-3.1-13.5,  percent_type='wo')
uo2.add_element('O',    13.5,          percent_type='wo')
uo2.set_density('g/cm3', 10.97)
uo2.depletable = True
uo2.volume = math.pi * fuel_r**2 * 1.0
print(f"연료 부피: {uo2.volume:.4f} cm³")

zircaloy = openmc.Material(name='Zircaloy-4')
zircaloy.add_element('Zr', 98.0, percent_type='wo')
zircaloy.add_element('Sn',  1.5, percent_type='wo')
zircaloy.add_element('Fe',  0.2, percent_type='wo')
zircaloy.add_element('Cr',  0.1, percent_type='wo')
zircaloy.set_density('g/cm3', 6.56)

water = openmc.Material(name='Coolant Water')
water.add_nuclide('H1',  2.0, percent_type='ao')
water.add_nuclide('O16', 1.0, percent_type='ao')
water.set_density('g/cm3', 0.71)
water.add_s_alpha_beta('c_H_in_H2O')

materials = openmc.Materials([uo2, zircaloy, water])
materials.export_to_xml()
print("✅ 재료 정의 완료")

# =============================================================================
# PART 2: 핀 셀 기하학
# =============================================================================

fuel_or = openmc.ZCylinder(r=fuel_r)
clad_or = openmc.ZCylinder(r=0.4750)

left   = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
right  = openmc.XPlane(x0=+pitch/2, boundary_type='reflective')
bottom = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
top    = openmc.YPlane(y0=+pitch/2, boundary_type='reflective')

fuel_cell  = openmc.Cell(name='Fuel',      fill=uo2,      region=-fuel_or)
clad_cell  = openmc.Cell(name='Clad',      fill=zircaloy, region=+fuel_or & -clad_or)
water_cell = openmc.Cell(name='Moderator', fill=water,    region=+clad_or & +left & -right & +bottom & -top)

universe = openmc.Universe(cells=[fuel_cell, clad_cell, water_cell])
geometry = openmc.Geometry(universe)
geometry.export_to_xml()
print("✅ 기하학 정의 완료")

# =============================================================================
# PART 3: 계산 설정
# =============================================================================

settings = openmc.Settings()
settings.run_mode  = 'eigenvalue'
settings.batches   = 100
settings.inactive  = 30
settings.particles = 5000

bounds = [-fuel_r, -fuel_r, -1, fuel_r, fuel_r, 1]
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(bounds[:3], bounds[3:])
)
settings.export_to_xml()
print("✅ 계산 설정 완료")

# =============================================================================
# PART 4: 연소 계산 설정
# =============================================================================

chain_file = openmc.deplete.Chain.from_xml(
    '/home/kevin/miniconda3/envs/openmc/lib/python3.11/site-packages/openmc_data/depletion/chain_endf_b8.0_pwr.xml'
)

model = openmc.Model(
    geometry=geometry,
    materials=materials,
    settings=settings
)

operator = openmc.deplete.CoupledOperator(
    model,
    chain_file
)

power = 174.0  # W

timesteps_days = [1,30,180,360,720,1080]

timesteps_seconds = [d * 86400 for d in timesteps_days]

integrator = openmc.deplete.PredictorIntegrator(
    operator,
    timesteps_seconds,
    power=power,
    timestep_units='s'
)

print("\n🚀 연소 계산 시작...")
print(f"   총 {len(timesteps_days)}개 스텝, 최대 {timesteps_days[-1]}일")
print("   (약 30~40분 소요)")

integrator.integrate()

# =============================================================================
# PART 5: 결과 분석 및 시각화
# =============================================================================

print("\n📊 결과 분석 중...")

results = openmc.deplete.Results('depletion_results.h5')

time_k, k_data = results.get_keff()
time_days = time_k / 86400

keff_vals = [k.nominal_value for k in k_data]
keff_errs = [k.std_dev for k in k_data]

_, u235  = results.get_atoms('1', 'U235')
_, u238  = results.get_atoms('1', 'U238')
_, pu239 = results.get_atoms('1', 'Pu239')
_, xe135 = results.get_atoms('1', 'Xe135')

u235_ratio  = u235  / u235[0]  * 100
pu239_ratio = pu239 / u235[0]  * 100

print(f"\n{'='*55}")
print(f"  연소 계산 결과 요약")
print(f"{'='*55}")
print(f"{'시간(일)':>8} | {'k-eff':>8} | {'U235(%)':>8} | {'Pu239(%)':>9}")
print("-" * 45)
for i in range(0, len(time_days), max(1, len(time_days)//8)):
    print(f"{time_days[i]:>8.1f} | {keff_vals[i]:>8.5f} | "
          f"{u235_ratio[i]:>8.2f} | {pu239_ratio[i]:>9.4f}")

# =============================================================================
# PART 6: 그래프
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 8))

axes[0,0].errorbar(time_days, keff_vals, yerr=keff_errs,
                   fmt='o-', color='steelblue', capsize=3, linewidth=2, markersize=4)
axes[0,0].axhline(y=1.0, color='red', linestyle='--', linewidth=1.5, label='Critical')
axes[0,0].set_xlabel('Time (days)')
axes[0,0].set_ylabel('k-infinity')
axes[0,0].set_title('k-infinity vs Burnup')
axes[0,0].legend()
axes[0,0].grid(True, alpha=0.3)

axes[0,1].plot(time_days, u235_ratio, 'o-', color='darkorange', linewidth=2, markersize=4)
axes[0,1].set_xlabel('Time (days)')
axes[0,1].set_ylabel('U-235 remaining (%)')
axes[0,1].set_title('U-235 Depletion')
axes[0,1].grid(True, alpha=0.3)

axes[1,0].plot(time_days, pu239_ratio, 'o-', color='purple', linewidth=2, markersize=4)
axes[1,0].set_xlabel('Time (days)')
axes[1,0].set_ylabel('Pu-239 / initial U-235 (%)')
axes[1,0].set_title('Pu-239 Production')
axes[1,0].grid(True, alpha=0.3)

axes[1,1].plot(time_days, xe135, 'o-', color='darkred', linewidth=2, markersize=4)
axes[1,1].set_xlabel('Time (days)')
axes[1,1].set_ylabel('Xe-135 atoms')
axes[1,1].set_title('Xe-135 Buildup (Reactor Poison)')
axes[1,1].grid(True, alpha=0.3)

plt.suptitle('PWR Pin Cell Depletion Analysis', fontsize=13)
plt.tight_layout()
plt.savefig('depletion_results.png', dpi=150, bbox_inches='tight')
print("\n✅ 그래프 저장됨: depletion_results.png")
