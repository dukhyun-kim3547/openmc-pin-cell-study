# =============================================================================
# 무한 핀 셀 시뮬레이션 (Infinite Pin Cell Simulation)
# 목적: 경수로(PWR) 연료봉의 중성자 증배계수 k-infinity 계산
# 코드: OpenMC 0.15.x (Monte Carlo 중성자 수송)
# =============================================================================

import openmc# =============================================================================
# PART 1: 재료(Materials) 정의
# 원자로에 들어가는 물질들의 핵적 성질을 정의한다.
# =============================================================================

# --- 재료 1: UO2 핵연료 ---
# 실제 PWR에서 사용하는 이산화우라늄 소결체
# 농축도 3.1 wt% = U-235가 전체 우라늄의 3.1%
uo2 = openmc.Material(name='UO2 Fuel')
uo2.add_nuclide('U235', 3.1,    percent_type='wo')  # 농축 우라늄-235 (핵분열 주체)
uo2.add_nuclide('U238', 96.9,   percent_type='wo')  # 열화 우라늄-238 (대부분 차지)
uo2.add_element('O', 13.5, percent_type='wo')  # UO2에서 산소 질량% = ~13.5%
uo2.set_density('g/cm3', 10.97)                     # UO2 이론밀도 ~10.97 g/cc

# --- 재료 2: Zircaloy-4 피복재 ---
# 중성자 흡수 단면적이 매우 작아 연료봉 피복재로 최적
# Zr이 주성분, 소량의 Sn/Fe/Cr 합금
zircaloy = openmc.Material(name='Zircaloy-4 Cladding')
zircaloy.add_element('Zr', 98.0, percent_type='wo') # 지르코늄 (주성분, 낮은 중성자 흡수)
zircaloy.add_element('Sn',  1.5, percent_type='wo') # 주석 (기계적 강도 향상)
zircaloy.add_element('Fe',  0.2, percent_type='wo') # 철 (내식성)
zircaloy.add_element('Cr',  0.1, percent_type='wo') # 크롬 (내식성)
zircaloy.set_density('g/cm3', 6.56)

# --- 재료 3: 냉각수(경수) ---
# 중성자 감속재 + 냉각재 역할
# PWR 운전 조건: 약 325°C, 155 bar → 밀도 0.71 g/cc
water = openmc.Material(name='Coolant Water')
water.add_nuclide('H1',  2.0, percent_type='ao')  # 수소-1 (최고의 감속재, 질량 ≈ 중성자)
water.add_nuclide('O16', 1.0, percent_type='ao')  # 산소-16
water.set_density('g/cm3', 0.71)                  # PWR 운전온도 밀도
water.add_s_alpha_beta('c_H_in_H2O')              # 열중성자 산란 보정 (S(α,β) 처리)
                                                   # 물 분자 결합으로 인한 열중성자 거동 보정

# 재료들을 하나의 컬렉션으로 묶기
materials = openmc.Materials([uo2, zircaloy, water])
materials.export_to_xml()
print("✅ 재료 정의 완료")

# =============================================================================
# PART 2: 기하학(Geometry) 정의
# 핀 셀의 물리적 구조를 CSG(Constructive Solid Geometry)로 정의한다.
# 동심원 실린더로 연료/피복재/냉각수 영역을 구분.
# =============================================================================

# --- 실린더 표면 정의 (무한 Z방향 원통) ---
# 실제 PWR 연료봉 치수 기준 (Westinghouse 17x17 집합체)
fuel_outer_radius  = openmc.ZCylinder(r=0.4096)  # 연료 펠렛 반지름 [cm]
clad_outer_radius  = openmc.ZCylinder(r=0.4750)  # 피복재 외경 반지름 [cm]
                                                   # 피복재 두께 = 0.4750 - 0.4096 = 0.0654 cm

# --- 경계 조건 정의 ---
# 핀 피치(pin pitch): 핀 셀 중심 간 거리 = 1.26 cm (17x17 격자)
# 반사 경계 = 무한 격자 모사 (주변이 동일한 핀 셀로 둘러싸인 것처럼 처리)
pitch = 1.26  # cm, PWR 표준 핀 피치
left   = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
right  = openmc.XPlane(x0=+pitch/2, boundary_type='reflective')
bottom = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
top    = openmc.YPlane(y0=+pitch/2, boundary_type='reflective')

# --- 셀(Cell) 정의: 각 영역에 재료 할당 ---
fuel_cell = openmc.Cell(name='Fuel')
fuel_cell.fill   = uo2
fuel_cell.region = -fuel_outer_radius  # 연료 반지름 내부

clad_cell = openmc.Cell(name='Cladding')
clad_cell.fill   = zircaloy
clad_cell.region = +fuel_outer_radius & -clad_outer_radius  # 연료 외부 & 피복재 내부

water_cell = openmc.Cell(name='Moderator')
water_cell.fill   = water
water_cell.region = +clad_outer_radius & +left & -right & +bottom & -top  # 피복재 외부 사각형

# --- 유니버스(Universe)와 지오메트리 완성 ---
universe  = openmc.Universe(cells=[fuel_cell, clad_cell, water_cell])
geometry  = openmc.Geometry(universe)
geometry.export_to_xml()
print("✅ 기하학 정의 완료")

# =============================================================================
# PART 3: 설정(Settings) 정의
# 몬테카를로 시뮬레이션의 계산 조건을 설정한다.
# =============================================================================

settings = openmc.Settings()
settings.run_mode    = 'eigenvalue'    # 임계 계산 모드 (k-eff 계산)
settings.batches     = 150            # 총 배치(batch) 수
settings.inactive    = 50             # 비활성 배치 수 (소스 분포 수렴용, 통계에서 제외)
settings.particles   = 5000          # 배치당 중성자 수
# 총 계산 중성자 = 5000 × (150-50) = 500,000개 → 통계 불확도 ~0.05%

# 초기 중성자 소스: 연료 영역 내 균일 분포
bounds = [-0.4096, -0.4096, -1, 0.4096, 0.4096, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
settings.source = openmc.IndependentSource(space=uniform_dist)

settings.export_to_xml()
print("✅ 계산 설정 완료")

# =============================================================================
# PART 4: 시뮬레이션 실행
# =============================================================================

print("\n🚀 시뮬레이션 시작...")
openmc.run()
print("\n✅ 시뮬레이션 완료!")

