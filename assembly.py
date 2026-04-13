# =============================================================================
# PWR 17×17 연료 집합체 시뮬레이션
# 목적: 실제 Westinghouse 17×17 집합체의 k-effective 계산
# Phase 2(핀 셀)와 비교: 누설, 안내관, 집합체 효과 확인
# =============================================================================

import openmc
import numpy as np

# =============================================================================
# PART 1: 재료 정의
# =============================================================================

# UO2 연료 (3.1% 농축)
uo2 = openmc.Material(name='UO2 Fuel')
uo2.add_nuclide('U235', 3.1,          percent_type='wo')
uo2.add_nuclide('U238', 100-3.1-13.5, percent_type='wo')
uo2.add_element('O',    13.5,         percent_type='wo')
uo2.set_density('g/cm3', 10.97)

# Zircaloy-4 피복재
zircaloy = openmc.Material(name='Zircaloy-4')
zircaloy.add_element('Zr', 98.0, percent_type='wo')
zircaloy.add_element('Sn',  1.5, percent_type='wo')
zircaloy.add_element('Fe',  0.2, percent_type='wo')
zircaloy.add_element('Cr',  0.1, percent_type='wo')
zircaloy.set_density('g/cm3', 6.56)

# 냉각수 (PWR 운전 조건)
water = openmc.Material(name='Coolant Water')
water.add_nuclide('H1',  2.0, percent_type='ao')
water.add_nuclide('O16', 1.0, percent_type='ao')
water.set_density('g/cm3', 0.71)
water.add_s_alpha_beta('c_H_in_H2O')

# 안내관 재료 (Zircaloy, 빈 관 = 물로 채워짐)
# 제어봉 삽입 시에는 B4C/Ag-In-Cd가 들어오지만
# 여기서는 제어봉 미삽입(ARO: All Rods Out) 상태로 가정
guide_tube_mat = openmc.Material(name='Guide Tube Water')
guide_tube_mat.add_nuclide('H1',  2.0, percent_type='ao')
guide_tube_mat.add_nuclide('O16', 1.0, percent_type='ao')
guide_tube_mat.set_density('g/cm3', 0.71)
guide_tube_mat.add_s_alpha_beta('c_H_in_H2O')

materials = openmc.Materials([uo2, zircaloy, water, guide_tube_mat])
materials.export_to_xml()
print("✅ 재료 정의 완료")

# =============================================================================
# PART 2: 핀 셀 유니버스 정의
# 두 종류의 핀: 연료봉 / 안내관
# =============================================================================

pitch = 1.26  # cm, 핀 피치

# --- 연료봉 핀 유니버스 ---
fuel_or = openmc.ZCylinder(r=0.4096)  # 연료 반지름
clad_or = openmc.ZCylinder(r=0.4750)  # 피복재 외경

fuel_cell  = openmc.Cell(fill=uo2,      region=-fuel_or)
clad_cell  = openmc.Cell(fill=zircaloy, region=+fuel_or & -clad_or)
water_cell = openmc.Cell(fill=water,    region=+clad_or)

fuel_pin_universe = openmc.Universe(cells=[fuel_cell, clad_cell, water_cell])

# --- 안내관 핀 유니버스 ---
# 안내관은 피복재만 있고 내부가 비어있음 (물로 채워짐)
# 내경 0.5610 cm, 외경 0.6020 cm (Westinghouse 17×17 기준)
gt_inner = openmc.ZCylinder(r=0.5610)  # 안내관 내경
gt_outer = openmc.ZCylinder(r=0.6020)  # 안내관 외경

gt_inner_cell = openmc.Cell(fill=guide_tube_mat, region=-gt_inner)
gt_clad_cell  = openmc.Cell(fill=zircaloy,       region=+gt_inner & -gt_outer)
gt_outer_cell = openmc.Cell(fill=water,          region=+gt_outer)

guide_tube_universe = openmc.Universe(cells=[gt_inner_cell, gt_clad_cell, gt_outer_cell])

print("✅ 핀 유니버스 정의 완료")

# =============================================================================
# PART 3: 17×17 격자(Lattice) 정의
# Westinghouse 17×17 집합체 배치도
# F = 연료봉, G = 안내관, I = 계측관(중앙, 안내관과 동일 구조)
# =============================================================================

# 17×17 배치 맵 (0 = 연료봉, 1 = 안내관/계측관)
# Westinghouse 표준 배치
layout = [
# 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 1
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 2
 [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0],  # 3
 [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],  # 4
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 5
 [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],  # 6
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 7
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 8
 [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],  # 9  ← 중앙(9,9)=계측관
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 10
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 11
 [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],  # 12
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 13
 [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],  # 14
 [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0],  # 15
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 16
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 17
]

# 배치 맵을 유니버스 배열로 변환
universes = []
for row in layout:
    uni_row = []
    for cell_type in row:
        if cell_type == 0:
            uni_row.append(fuel_pin_universe)
        else:
            uni_row.append(guide_tube_universe)
    universes.append(uni_row)

# RectLattice: 17×17 격자 생성
assembly_lattice = openmc.RectLattice(name='17x17 Assembly')
assembly_lattice.pitch      = (pitch, pitch)           # 핀 피치 (x, y)
assembly_lattice.lower_left = (-17*pitch/2, -17*pitch/2)  # 격자 좌하단 좌표
assembly_lattice.universes  = universes                # 유니버스 배열

print("✅ 17×17 격자 정의 완료")
print(f"   연료봉 수: {sum(row.count(0) for row in layout)}개")
print(f"   안내관 수: {sum(row.count(1) for row in layout)}개")

# =============================================================================
# PART 4: 집합체 외부 경계 정의
# =============================================================================

# 집합체 전체 크기: 17 × 1.26 cm = 21.42 cm
half = 17 * pitch / 2  # 10.71 cm

# 집합체 바깥 경계면 (반사 경계: 무한 격자 가정)
left   = openmc.XPlane(x0=-half, boundary_type='reflective')
right  = openmc.XPlane(x0=+half, boundary_type='reflective')
bottom = openmc.YPlane(y0=-half, boundary_type='reflective')
top    = openmc.YPlane(y0=+half, boundary_type='reflective')

# Z방향 경계 (축방향 무한 가정)
bot_z = openmc.ZPlane(z0=-20, boundary_type='reflective')
top_z = openmc.ZPlane(z0=+20, boundary_type='reflective')

# 집합체 셀
assembly_cell = openmc.Cell(name='Assembly')
assembly_cell.fill   = assembly_lattice
assembly_cell.region = +left & -right & +bottom & -top & +bot_z & -top_z

root_universe = openmc.Universe(cells=[assembly_cell])
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()
print("✅ 기하학 정의 완료")

# =============================================================================
# PART 5: 계산 설정
# 집합체는 핀 셀보다 크므로 입자 수 증가
# =============================================================================

settings = openmc.Settings()
settings.run_mode  = 'eigenvalue'
settings.batches   = 200           # 배치 수 증가 (집합체는 더 많은 통계 필요)
settings.inactive  = 50
settings.particles = 10000         # 배치당 입자 수 증가 (핀 셀의 2배)

# 초기 소스: 연료 영역 전체에 균일 분포
bounds = [-half, -half, -18, half, half, 18]
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(bounds[:3], bounds[3:])
)
settings.export_to_xml()
print("✅ 계산 설정 완료")
print(f"   총 중성자 수: {settings.particles * (settings.batches - settings.inactive):,}개")

# =============================================================================
# PART 6: 시뮬레이션 실행
# =============================================================================

print("\n🚀 17×17 집합체 시뮬레이션 시작...")
print("   (핀 셀보다 오래 걸려요, 약 5~15분)")
openmc.run()

# 결과 읽기
sp = openmc.StatePoint('statepoint.200.h5')
keff = sp.keff
print(f"\n{'='*50}")
print(f"  17×17 집합체 k-effective = {keff.nominal_value:.5f} ± {keff.std_dev:.5f}")
print(f"{'='*50}")
print(f"\n  Phase 2 핀 셀 k-infinity  = 1.34052")
print(f"  Phase 3 집합체 k-effective = {keff.nominal_value:.5f}")
print(f"  차이 (안내관 + 격자 효과)  = {keff.nominal_value - 1.34052:+.5f}")
sp.close()
