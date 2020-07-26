#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <vector>

#include "cell_index.hpp"
#include "neighbor_info.hpp"
#include "vector2.hpp"
using namespace std;
/*
csvファイルを書き出すインターバル
*/
constexpr int INTERVAL = 10;
constexpr int ARRAY_SIZE = 30000;
constexpr double MYU = 0.01;
constexpr double PRESSURE_COEFF = 1;
constexpr double PARTICLE_DISTANCE = 0.005;
constexpr double SMOOTHING_LENGTH = 3.3 * PARTICLE_DISTANCE;
constexpr double SMOOTHING_LENGTH_SQ = SMOOTHING_LENGTH * SMOOTHING_LENGTH;
constexpr double SMOOTHING_LENGTH_SQ_INV = 1.0 / SMOOTHING_LENGTH_SQ;
constexpr double COLLISION_COEFF = 0.3;

constexpr int NEIGHBOR_SIZE = 100;
constexpr double PI = 3.1415926535897932384626433832;
constexpr double MASS =
    PI * PARTICLE_DISTANCE * PARTICLE_DISTANCE * 0.25 * 1000;
/*
cell linked listに用いる正方形のセルの一辺の長さ
*/
constexpr double cell_len = 1.01 * SMOOTHING_LENGTH;
/*
初期状態のセルの存在する範囲のx,y座標の最大/最小値
*/
constexpr double initial_min_x = -1 * PARTICLE_DISTANCE;
constexpr double initial_max_x = 260 * PARTICLE_DISTANCE;
constexpr double initial_min_y = -1 * PARTICLE_DISTANCE;
constexpr double initial_max_y = 260 * PARTICLE_DISTANCE;
constexpr double RHO = 1000.0;
constexpr double DELTA_T = 0.0005;
constexpr int cell_num =
    max((initial_max_x - initial_min_x), (initial_max_y - initial_min_y)) /
    cell_len;
int particle_num;
/*
壁と流体粒子の初期配置に使う
*/
constexpr int initial_wall_width = 4;
constexpr int initial_fluid_x = 100;
constexpr int initial_fluid_y = 150;
constexpr int initial_wall_y = 250;
constexpr int initial_wall_x = 250;
/*
セルの存在する範囲のx,y座標の最大/最小値(セルの存在範囲は移動する)
*/
double min_x, max_x, min_y, max_y;
/*
壁と流体の衝突処理を行うx,y座標
但し時間とともに変化する
*/
double floor_wall = (initial_wall_width - 1 + 0.01) * PARTICLE_DISTANCE;
double left_wall = (initial_wall_width - 1 + 0.01) * PARTICLE_DISTANCE;
double right_wall = (initial_wall_x + 1 - 0.01) * PARTICLE_DISTANCE;
double ceil_wall = (initial_wall_y + 1 - 0.01) * PARTICLE_DISTANCE;

vector2 Position[ARRAY_SIZE], Velocity[ARRAY_SIZE], Acceleration[ARRAY_SIZE];

double Pressure[ARRAY_SIZE];
/*
セルは粒子番号をリストで保持する
tmp_cellはセルの更新で一時的に用いる
*/
list<int> cell[cell_num][cell_num], tmp_cell[cell_num][cell_num];
/*
セルが存在しない範囲に位置する粒子の番号を保持するリスト
*/
list<int> out_of_region;
/*
近傍のセル番号を保持する
*/
cell_index neighbor_cell[9];
/*
セル内の粒子の更新に用いる
*/
cell_index tmp_cell_index;
/*
近傍の粒子情報を保持する
*/
neighbor_info neighbor[NEIGHBOR_SIZE];
double density[ARRAY_SIZE];
double density_0;
const vector2 g = vector2(0, -1.0);
/*
粒子の種類
*/
enum particle_type { fluid, wall };

particle_type ptype[ARRAY_SIZE];

double polynomial_kernel(const double r_sq, const double h_sq,
                         const double h_sq_inv) {
  if (r_sq > h_sq) {
    return 0;
  }
  const double a = 1 - r_sq * h_sq_inv;
  return 4.0 / PI * a * a * a * h_sq_inv;
}
vector2 grad_spiky_kernel(const vector2& r_vec, const double r,
                          const double h) {
  if (r > h) {
    return vector2(0, 0);
  }
  const double a = 1 - r / h;
  return r_vec * (-30 / PI * a * a / (h * h * h * r));
}
double laplacian_viscosity_kernel(const double r, const double h,
                                  const double h_sq_inv) {
  if (r > h) {
    return 0;
  }
  return 20 / (3 * PI) * (1 - r / h) * h_sq_inv * h_sq_inv;
}
vector2 calc_viscosity(const double r, const double r_sq, const vector2& r_ji,
                       const int i_index, const int j_index, const double h,
                       const double h_sq_inv) {
  const vector2 vij = Velocity[j_index] - Velocity[i_index];
  // density[i_index]で割るべきだが、割らないほうが安定したので割らずにパラメータを調整した
  return vij * (MYU * MASS * laplacian_viscosity_kernel(r, h, h_sq_inv) /
                (/*density[i_index] * */ density[j_index]));
}
vector2 calc_pressure_gradient(const double r, const double r_sq,
                               const vector2& rji, const int i_index,
                               const int j_index, const double h) {
  const double p_ij_average = (Pressure[i_index] + Pressure[j_index]) * 0.5;
  // density[i_index]で割るべきだが、割らないほうが安定した割らないほうが安定したので割らずにパラメータを調整した
  return grad_spiky_kernel(rji, r, h) *
         (-p_ij_average * MASS / (/*density[i_index] * */ density[j_index]));
}
/*
セルiがセルの存在範囲内にあるかどうか
*/
bool is_in_cell(const cell_index& i) {
  return 0 <= i.x && i.x < cell_num && 0 <= i.y && i.y < cell_num;
}
bool is_in_cell(const int x, const int y) {
  return 0 <= x && x < cell_num && 0 <= y && y < cell_num;
}
/*
近傍セルの探索
*/
int neighbor_cell_search(const int x, const int y) {
  int res = 0;
  for (int i = -1; i < 2; i++) {
    for (int j = -1; j < 2; j++) {
      const int nx = x + i;
      const int ny = y + j;
      if (is_in_cell(nx, ny) == true) {
        neighbor_cell[res].set(nx, ny);
        res++;
      }
    }
  }
  return res;
}
/*
近傍粒子の探索
*/
int neighbor_particle_search(const int particle_index,
                             const int neighbor_cell_num,
                             const double radius_sq) {
  int res = 0;
  const vector2& self_pos = Position[particle_index];
  for (int i = 0; i < neighbor_cell_num; i++) {
    for (const auto& itr : cell[neighbor_cell[i].x][neighbor_cell[i].y]) {
      const vector2 xij = Position[itr] - self_pos;
      const double d_sq = xij.square();
      if (d_sq < radius_sq) {
        neighbor[res].set(itr, d_sq, xij);
        res++;
      }
    }
  }
  return res;
}
/*
粘性項と圧力勾配項
*/
void calc_force() {
  for (int i = 0; i < cell_num; i++) {
    for (int j = 0; j < cell_num; j++) {
      const int cell_num = neighbor_cell_search(i, j);
      for (const auto& itr : cell[i][j]) {
        if (ptype[itr] == wall) {
          continue;
        }
        Acceleration[itr] = g;
        // if (density[itr] < 0.8 * density_0) {
        //   continue;
        // }
        const int neighbor_num =
            neighbor_particle_search(itr, cell_num, SMOOTHING_LENGTH_SQ);
        for (int k = 0; k < neighbor_num; k++) {
          if (neighbor[k].index == itr) {
            continue;
          }
          const double r_sq = neighbor[k].dist_sq;
          const double r = sqrt(r_sq);
          const vector2 rji = neighbor[k].v * -1.0;
          Acceleration[itr] +=
              calc_viscosity(r, r_sq, rji, itr, neighbor[k].index,
                             SMOOTHING_LENGTH, SMOOTHING_LENGTH_SQ_INV);
          Acceleration[itr] += calc_pressure_gradient(
              r, r_sq, rji, itr, neighbor[k].index, SMOOTHING_LENGTH);
        }
      }
    }
  }
}
double calc_density(const int neighbor_num) {
  double res = 0;
  for (int i = 0; i < neighbor_num; i++) {
    res += polynomial_kernel(neighbor[i].dist_sq, SMOOTHING_LENGTH_SQ,
                             SMOOTHING_LENGTH_SQ_INV);
  }
  return MASS * res;
}
void update_density() {
  for (int i = 0; i < cell_num; i++) {
    for (int j = 0; j < cell_num; j++) {
      const int neighbor_cell_num = neighbor_cell_search(i, j);
      for (const auto& itr : cell[i][j]) {
        const int neighbor_particle_num = neighbor_particle_search(
            itr, neighbor_cell_num, SMOOTHING_LENGTH_SQ);
        density[itr] = calc_density(neighbor_particle_num);
      }
    }
  }
}
void get_cell_index(const vector2& v) {
  const int x_index = (v.x - min_x) / cell_len;
  const int y_index = (v.y - min_y) / cell_len;
  tmp_cell_index.set(x_index, y_index);
}
/*
セルの初期化を行う
粒子をセルの格納する
*/
void initialize_cell() {
  for (int i = 0; i < particle_num; i++) {
    get_cell_index(Position[i]);
    cell[tmp_cell_index.x][tmp_cell_index.y].push_back(i);
  }
}
/*
セルの存在範囲外に存在する粒子を格納するout_of_regionを更新する
*/
void update_out_of_region() {
  for (auto itr = out_of_region.begin(); itr != out_of_region.end();) {
    get_cell_index(Position[*itr]);
    if (is_in_cell(tmp_cell_index) == false) {
      itr++;
      continue;
    }
    cell[tmp_cell_index.x][tmp_cell_index.y].push_back(*itr);
    itr = out_of_region.erase(itr);
  }
}
/*
全セルを更新する
*/
void update_cell() {
  for (int i = 0; i < cell_num; i++) {
    for (int j = 0; j < cell_num; j++) {
      for (auto itr = cell[i][j].begin(); itr != cell[i][j].end();) {
        get_cell_index(Position[*itr]);
        if (is_in_cell(tmp_cell_index) == false) {
          out_of_region.push_back(*itr);
          itr = cell[i][j].erase(itr);
          continue;
        }
        if (tmp_cell_index.is_same(i, j) == false) {
          tmp_cell[tmp_cell_index.x][tmp_cell_index.y].push_back(*itr);
          itr = cell[i][j].erase(itr);
          continue;
        }
        itr++;
      }
    }
  }
  for (int i = 0; i < cell_num; i++) {
    for (int j = 0; j < cell_num; j++) {
      if (tmp_cell[i][j].empty() == true) {
        continue;
      }
      cell[i][j].splice(cell[i][j].end(), tmp_cell[i][j]);
    }
  }
  update_out_of_region();
}
void update_position() {
  for (int i = 0; i < particle_num; i++) {
    Position[i] += Velocity[i] * DELTA_T;
  }
  update_cell();
}
void update_velocity() {
  for (int i = 0; i < particle_num; i++) {
    if (ptype[i] != wall) {
      Velocity[i] += Acceleration[i] * DELTA_T;
    }
  }
}
void print_csv(const int time_step) {
  ofstream ofs("sph" + std::to_string(time_step) + ".csv");
  if (!ofs) {
    cout << "file open failed" << std::endl;
    exit(true);
  }
  ofs << "x,y,Velocity_x,Velocity_y,Pressure,Particle_type" << std::endl;
  for (int i = 0; i < particle_num; i++) {
    ofs << Position[i].x << "," << Position[i].y << "," << Velocity[i].x << ","
        << Velocity[i].y << "," << Pressure[i] << "," << ptype[i] << endl;
  }
  ofs.close();
}
void calc_pressure() {
  for (int i = 0; i < particle_num; i++) {
    if (density[i] < density_0) {
      Pressure[i] = 0;
      continue;
    }
    const double a = density[i] / density_0;
    const double a_sq = a * a;
    Pressure[i] = PRESSURE_COEFF * (a * a_sq * a_sq * a_sq - 1.0);
  }
}
/*
流体粒子の初期配置に使う
*/
bool is_in_fluid(const int x, const int y) {
  return initial_wall_width <= x && x <= initial_fluid_x &&
         initial_wall_width <= y && y <= initial_fluid_y;
}
/*
壁粒子の初期配置に使う
*/
bool is_not_in_wall(const int x, const int y) {
  return initial_wall_width <= x && x <= initial_wall_x &&
         initial_wall_width <= y && y <= initial_wall_y;
}
/*
壁と流体粒子の衝突処理
*/
void collision_wall(const double floor_wall, const double left_wall,
                    const double right_wall, const double ceil_wall,
                    const double x_velocity) {
  for (int i = 0; i < particle_num; i++) {
    if (ptype[i] == wall) {
      continue;
    }
    const auto& pos = Position[i];
    if (pos.y < floor_wall) {
      Velocity[i].y *= -COLLISION_COEFF;
      Position[i].y = floor_wall;
    }
    if (pos.y > ceil_wall) {
      Velocity[i].y *= -COLLISION_COEFF;
      Position[i].y = ceil_wall;
    }
    if (pos.x < left_wall) {
      Velocity[i].x =
          x_velocity + COLLISION_COEFF * (x_velocity - Velocity[i].x);
      Position[i].x = left_wall;
    }
    if (pos.x > right_wall) {
      Velocity[i].x =
          x_velocity - COLLISION_COEFF * (Velocity[i].x - x_velocity);
      Position[i].x = right_wall;
    }
  }
}
/*
壁速度の設定
*/
void set_wall_velocity(const vector2& wall_velocity) {
  for (int i = 0; i < particle_num; i++) {
    if (ptype[i] == wall) {
      Velocity[i] = wall_velocity;
    }
  }
}
/*
シミュレーションを1ステップ更新する
*/
void update_simulation(const int time, const vector2& wall_vel) {
  set_wall_velocity(wall_vel);
  /*
  セルの存在範囲を更新
  */
  min_x += wall_vel.x * DELTA_T;
  max_x += wall_vel.x * DELTA_T;
  /*
  衝突処理を行う座標を更新
  */
  left_wall += wall_vel.x * DELTA_T;
  right_wall += wall_vel.x * DELTA_T;

  update_density();
  calc_pressure();
  calc_force();
  update_velocity();
  update_position();
  collision_wall(floor_wall, left_wall, right_wall, ceil_wall, wall_vel.x);
  if (time % INTERVAL == 0) {
    print_csv(time / INTERVAL);
  }
}
int main() {
  /*
  セルの存在するx,y座標の最大/最小値(時間とともに動く)
  */
  min_x = initial_min_x;
  max_x = initial_max_x;
  max_y = initial_max_y;
  min_y = initial_min_y;

  /*
  基準密度の計算
  */
  density_0 = 0;
  for (int i = -10; i < 10; i++) {
    for (int j = -10; j < 10; j++) {
      const double dist_sq =
          PARTICLE_DISTANCE * PARTICLE_DISTANCE * (i * i + j * j);
      const double dist = sqrt(dist_sq);
      if (dist < SMOOTHING_LENGTH) {
        density_0 += polynomial_kernel(dist_sq, SMOOTHING_LENGTH_SQ,
                                       SMOOTHING_LENGTH_SQ_INV);
      }
    }
  }
  density_0 *= MASS;
  cout << density_0 << endl;
  cout << cell_num << endl;
  /*
  粒子の初期配置をする
  */
  int res = 0;
  for (int i = 0; i < initial_wall_x + 4; i++) {
    for (int j = 0; j < initial_wall_y + 4; j++) {
      if (is_in_fluid(i, j) == true) {
        Position[res].set(i * PARTICLE_DISTANCE, j * PARTICLE_DISTANCE);
        ptype[res] = fluid;
      } else if (!is_not_in_wall(i, j)) {
        Position[res].set(i * PARTICLE_DISTANCE, j * PARTICLE_DISTANCE);
        ptype[res] = wall;
      } else {
        continue;
      }
      res++;
    }
  }
  /*
  粒子数の設定
  */
  particle_num = res;
  fill(Velocity, Velocity + ARRAY_SIZE, vector2(0, 0));
  /*
  セルの初期化
  */
  initialize_cell();

  cout << particle_num << endl;
  const vector2 right_vel = vector2(0.2, 0);
  const vector2 left_vel = right_vel * -1;
  const vector2 static_vel = vector2(0, 0);
  vector2 wall_vel;
  int cnt = 0;
  while (cnt < 12000) {
    cout << "time=" << cnt << endl;
    /*
    7200回タイムステップあとは静止
    */
    if (cnt < 7200) {
      /*
      1200回タイムステップごとに壁の速度を変える
      */
      const int wall_move = (cnt / 1200) % 3;
      if (wall_move == 0) {
        wall_vel = static_vel;
      } else if (wall_move == 1) {
        wall_vel = right_vel;
      } else {
        wall_vel = left_vel;
      }
    } else {
      wall_vel = static_vel;
    }
    update_simulation(cnt, wall_vel);
    cnt++;
  }

  return 0;
}
