# =============================================================================
# 중성자 플럭스 분포 시각화 (Neutron Flux Distribution)
# 목적: 17×17 집합체 내 열중성자/빠른중성자 flux 2D 분포 계산 및 시각화
# 핵심: Mesh Tally를 이용해 공간별 중성자 밀도 측정
# =============================================================================

import openmc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# =============================================================================
# PART 1~4: 집합체 모델 (assembly.py와 동일)
# =============================================================================

# --- 재료 ---
uo2 = openmc.Material(name='UO2 Fuel')
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

water = openmc.Material(name='Coolant Water')
water.add_nuclide('H1',  2.0, percent_type='ao')
water.add_nuclide('O16', 1.0, percent_type='ao')
water.set_density('g/cm3', 0.71)
water.add_s_alpha_beta('c_H_in_H2O')

guide_tube_mat = openmc.Material(name='Guide Tube Water')
guide_tube_mat.add_nuclide('H1',  2.0, percent_type='ao')
guide_tube_mat.add_nuclide('O16', 1.0, percent_type='ao')
guide_tube_mat.set_density('g/cm3', 0.71)
guide_tube_mat.add_s_alpha_beta('c_H_in_H2O')

materials = openmc.Materials([uo2, zircaloy, water, guide_tube_mat])
materials.export_to_xml()

# --- 핀 유니버스 ---
pitch = 1.26
fuel_or = openmc.ZCylinder(r=0.4096)
clad_or = openmc.ZCylinder(r=0.4750)
fuel_cell  = openmc.Cell(fill=uo2,      region=-fuel_or)
clad_cell  = openmc.Cell(fill=zircaloy, region=+fuel_or & -clad_or)
water_cell = openmc.Cell(fill=water,    region=+clad_or)
fuel_pin_universe = openmc.Universe(cells=[fuel_cell, clad_cell, water_cell])

gt_inner = openmc.ZCylinder(r=0.5610)
gt_outer = openmc.ZCylinder(r=0.6020)
gt_inner_cell = openmc.Cell(fill=guide_tube_mat, region=-gt_inner)
gt_clad_cell  = openmc.Cell(fill=zircaloy,       region=+gt_inner & -gt_outer)
gt_outer_cell = openmc.Cell(fill=water,          region=+gt_outer)
guide_tube_universe = openmc.Universe(cells=[gt_inner_cell, gt_clad_cell, gt_outer_cell])

# --- 17×17 격자 ---
layout = [
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0],
 [0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0],
 [0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
]
universes = []
for row in layout:
    uni_row = []
    for ct in row:
        uni_row.append(fuel_pin_universe if ct == 0 else guide_tube_universe)
    universes.append(uni_row)

assembly_lattice = openmc.RectLattice(name='17x17 Assembly')
assembly_lattice.pitch      = (pitch, pitch)
assembly_lattice.lower_left = (-17*pitch/2, -17*pitch/2)
assembly_lattice.universes  = universes

half  = 17 * pitch / 2
left   = openmc.XPlane(x0=-half, boundary_type='reflective')
right  = openmc.XPlane(x0=+half, boundary_type='reflective')
bottom = openmc.YPlane(y0=-half, boundary_type='reflective')
top    = openmc.YPlane(y0=+half, boundary_type='reflective')
bot_z  = openmc.ZPlane(z0=-20,  boundary_type='reflective')
top_z  = openmc.ZPlane(z0=+20,  boundary_type='reflective')

assembly_cell = openmc.Cell(fill=assembly_lattice,
                            region=+left & -right & +bottom & -top & +bot_z & -top_z)
root_universe = openmc.Universe(cells=[assembly_cell])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# =============================================================================
# PART 5: Mesh Tally 정의 ← 이번 핵심!
# Mesh Tally = 공간을 격자로 나눠서 각 격자 셀의 flux를 측정
# =============================================================================

# 170×170 메쉬 (핀 하나당 10×10 = 100개 격자)
# 해상도가 높을수록 정밀하지만 계산 시간 증가
mesh = openmc.RegularMesh(name='flux_mesh')
mesh.dimension    = [170, 170]          # x, y 방향 격자 수
mesh.lower_left   = [-half, -half]      # 메쉬 시작점
mesh.upper_right  = [+half, +half]      # 메쉬 끝점

mesh_filter = openmc.MeshFilter(mesh)

# 에너지 필터: 열중성자 / 빠른중성자 분리
# 열중성자: E < 0.625 eV (감속된 중성자, 핵분열 주체)
# 빠른중성자: E > 0.625 eV (핵분열에서 갓 태어난 중성자)
energy_filter = openmc.EnergyFilter([0, 0.625, 20.0e6])
#                                     ^ 열   ^빠른  ^ MeV 단위

# Tally 1: 전체 flux (에너지 구분 없음)
tally_total = openmc.Tally(name='total_flux')
tally_total.filters = [mesh_filter]
tally_total.scores  = ['flux']

# Tally 2: 열중성자 / 빠른중성자 분리 flux
tally_energy = openmc.Tally(name='energy_flux')
tally_energy.filters = [mesh_filter, energy_filter]
tally_energy.scores  = ['flux']

tallies = openmc.Tallies([tally_total, tally_energy])
tallies.export_to_xml()
print("✅ Mesh Tally 정의 완료 (170×170 메쉬)")

# =============================================================================
# PART 6: 계산 설정 (입자 수 증가 → 높은 통계 필요)
# =============================================================================

settings = openmc.Settings()
settings.run_mode  = 'eigenvalue'
settings.batches   = 200
settings.inactive  = 50
settings.particles = 20000   # 플럭스 분포는 더 많은 입자 필요

bounds = [-half, -half, -18, half, half, 18]
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(bounds[:3], bounds[3:])
)
settings.export_to_xml()
print(f"✅ 계산 설정 완료 (총 {settings.particles*(settings.batches-settings.inactive):,}개 중성자)")

# =============================================================================
# PART 7: 시뮬레이션 실행
# =============================================================================

print("\n🚀 플럭스 분포 계산 시작...")
print("   (약 15~30분 소요)")
openmc.run()

# =============================================================================
# PART 8: 결과 읽기 및 시각화
# =============================================================================

print("\n📊 결과 처리 및 그래프 생성 중...")

sp = openmc.StatePoint('statepoint.200.h5')

# 전체 flux 데이터 추출
tally_total  = sp.get_tally(name='total_flux')
tally_energy = sp.get_tally(name='energy_flux')

flux_total = tally_total.get_values(scores=['flux']).reshape(170, 170)

# 열/빠른 중성자 분리
flux_all    = tally_energy.get_values(scores=['flux']).reshape(170, 170, 2)
flux_thermal = flux_all[:, :, 0]   # 첫 번째 에너지 그룹 (열중성자)
flux_fast    = flux_all[:, :, 1]   # 두 번째 에너지 그룹 (빠른중성자)

sp.close()

# --- 그래프 생성 ---
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

plots = [
    (flux_total,   'Total Flux',          'hot'),
    (flux_thermal, 'Thermal Flux (E<0.625eV)', 'Blues'),
    (flux_fast,    'Fast Flux (E>0.625eV)',    'Reds'),
]

for ax, (data, title, cmap) in zip(axes, plots):
    extent = [-half, half, -half, half]
    im = ax.imshow(data, extent=extent, origin='lower',
                   cmap=cmap, norm=mcolors.LogNorm(
                       vmin=max(data.max()*1e-4, 1e-30),
                       vmax=data.max()
                   ))
    ax.set_title(title, fontsize=11)
    ax.set_xlabel('x (cm)', fontsize=10)
    ax.set_ylabel('y (cm)', fontsize=10)
    plt.colorbar(im, ax=ax, label='Flux (n/cm²/src)', shrink=0.8)

plt.suptitle('PWR 17×17 Assembly Neutron Flux Distribution', fontsize=13)
plt.tight_layout()
plt.savefig('flux_distribution.png', dpi=150, bbox_inches='tight')
print("\n✅ 그래프 저장됨: flux_distribution.png")
print("\n[플럭스 분포 통계]")
print(f"  전체 flux 최대/최소 비: {flux_total.max()/flux_total[flux_total>0].min():.1f}")
print(f"  열중성자 flux 최대:  {flux_thermal.max():.3e}")
print(f"  빠른중성자 flux 최대: {flux_fast.max():.3e}")
