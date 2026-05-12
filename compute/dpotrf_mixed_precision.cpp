/**
 *
 * @file dpotrf.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated d Fri Apr  1 11:02:54 2016
 *
 **/
#ifdef PLASMA_WITH_CUDA
#if PLASMA_WITH_HPCASIA24
#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>

#include "common.h"
#include "cuda_fp8.h"
#include "mixed_precision.h"
#include "plasma_d_mixed.h"
#include "mx_wrapper.h"

// Env toggle: when set to 1, build tiles directly on GPU to avoid host
// column-major -> tile conversion. Fallback to host path otherwise.
static bool use_gpu_tiling() {
  static int init = 0;
  static bool enabled = false;
  if (!init) {
    const char *env = getenv("MX_GPU_TILING");
    enabled = (env && env[0] == '1');
    init = 1;
  }
  return enabled;
}

static std::string forced_format() {
  static int init = 0;
  static std::string fmt;
  if (!init) {
    if (const char *env = getenv("MX_FORCE_FORMAT")) {
      fmt = env;
      std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    }
    init = 1;
  }
  return fmt;
}

static std::string bucket_format(const char *env_name) {
  if (const char *env = getenv(env_name)) {
    std::string fmt = env;
    std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return fmt;
  }
  return {};
}

static long double format_unit_roundoff(const std::string &fmt_in) {
  std::string fmt = fmt_in;
  std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (fmt == "fp32" || fmt == "mx_fp32" || fmt == "mx_f32") {
    return static_cast<long double>(std::numeric_limits<float>::epsilon());
  }
  if (fmt == "fp16" || fmt == "mx_fp16" || fmt == "mx_f16") {
    return static_cast<long double>(std::numeric_limits<Eigen::half>::epsilon());
  }
  if (fmt == "bf16") {
    return 1.0L / 128.0L;
  }
  if (fmt == "mx_e4m3" || fmt == "mx_fp8_e4m3" || fmt == "e4m3" ||
      fmt == "fp8_e4m3" || fmt == "fp8e4m3") {
    return 1.0L / 16.0L;
  }
  if (fmt == "mx_e5m2" || fmt == "mx_fp8_e5m2" || fmt == "e5m2" ||
      fmt == "fp8_e5m2" || fmt == "fp8e5m2" || fmt == "e3m2") {
    return 1.0L / 8.0L;
  }
  if (fmt == "e2m3") {
    return 1.0L / 16.0L;
  }
  if (fmt == "e2m1") {
    return 1.0L / 4.0L;
  }
  return 1.0L / 16.0L;
}

static bool format_uses_shared_scale(const std::string &fmt_in) {
  std::string fmt = fmt_in;
  std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return (fmt == "mx_fp16" || fmt == "mx_f16" || fmt == "mx_fp32" ||
          fmt == "mx_f32" || fmt == "mx_e4m3" || fmt == "mx_fp8_e4m3" ||
          fmt == "mx_e5m2" || fmt == "mx_fp8_e5m2" || fmt == "e3m2" ||
          fmt == "e2m3" || fmt == "e2m1");
}

static long double format_fmin(const std::string &fmt_in) {
  std::string fmt = fmt_in;
  std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (fmt == "fp32" || fmt == "mx_fp32" || fmt == "mx_f32") {
    return std::ldexp(1.0L, -126);
  }
  if (fmt == "fp16" || fmt == "mx_fp16" || fmt == "mx_f16") {
    return std::ldexp(1.0L, -14);
  }
  if (fmt == "bf16") {
    return std::ldexp(1.0L, -126);
  }
  if (fmt == "mx_e4m3" || fmt == "mx_fp8_e4m3" || fmt == "e4m3" ||
      fmt == "fp8_e4m3" || fmt == "fp8e4m3") {
    return std::ldexp(1.0L, -6);
  }
  if (fmt == "mx_e5m2" || fmt == "mx_fp8_e5m2" || fmt == "e5m2" ||
      fmt == "fp8_e5m2" || fmt == "fp8e5m2") {
    return std::ldexp(1.0L, -14);
  }
  if (fmt == "e3m2") {
    return std::ldexp(1.0L, -2);
  }
  if (fmt == "e2m3") {
    return std::ldexp(1.0L, 0);
  }
  if (fmt == "e2m1") {
    return std::ldexp(1.0L, 0);
  }
  return 0.0L;
}

static long double format_fmax(const std::string &fmt_in) {
  std::string fmt = fmt_in;
  std::transform(fmt.begin(), fmt.end(), fmt.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (fmt == "fp32" || fmt == "mx_fp32" || fmt == "mx_f32") {
    return (std::numeric_limits<float>::max)();
  }
  if (fmt == "fp16" || fmt == "mx_fp16" || fmt == "mx_f16") {
    return 65504.0L;
  }
  if (fmt == "bf16") {
    return std::ldexp(1.0L, 127) * (2.0L - std::ldexp(1.0L, -7));
  }
  if (fmt == "mx_e4m3" || fmt == "mx_fp8_e4m3" || fmt == "e4m3" ||
      fmt == "fp8_e4m3" || fmt == "fp8e4m3") {
    return 448.0L;
  }
  if (fmt == "mx_e5m2" || fmt == "mx_fp8_e5m2" || fmt == "e5m2" ||
      fmt == "fp8_e5m2" || fmt == "fp8e5m2") {
    return 57344.0L;
  }
  if (fmt == "e3m2") {
    return 28.0L;
  }
  if (fmt == "e2m3") {
    return 7.5L;
  }
  if (fmt == "e2m1") {
    return 6.0L;
  }
  return 0.0L;
}

static long double format_scaled_underflow_term(const std::string &fmt_in) {
  const long double fmin = format_fmin(fmt_in);
  const long double fmax = format_fmax(fmt_in);
  if (!(fmin > 0.0L) || !(fmax > 0.0L)) return 0.0L;

  static int mode_init = 0;
  static bool use_gu = false;
  if (!mode_init) {
    if (const char *env = getenv("MX_UNDERFLOW_MODE")) {
      std::string m = env;
      std::transform(m.begin(), m.end(), m.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      use_gu = (m == "gu" || m == "gradual");
    }
    mode_init = 1;
  }

  const long double u = format_unit_roundoff(fmt_in);
  const long double gmin = use_gu ? (u * fmin) : (0.5L * fmin);
  return gmin / fmax;
}

static long double format_gmin(const std::string &fmt_in) {
  const long double fmin = format_fmin(fmt_in);
  if (!(fmin > 0.0L)) return 0.0L;

  static int mode_init = 0;
  static bool use_gu = false;
  if (!mode_init) {
    if (const char *env = getenv("MX_UNDERFLOW_MODE")) {
      std::string m = env;
      std::transform(m.begin(), m.end(), m.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      use_gu = (m == "gu" || m == "gradual");
    }
    mode_init = 1;
  }

  const long double u = format_unit_roundoff(fmt_in);
  return use_gu ? (u * fmin) : (0.5L * fmin);
}

static bool use_explicit_tile_bound_selection() {
  static int init = 0;
  static bool enabled = false;
  if (!init) {
    if (const char *env = getenv("MX_SELECTION_CRITERIA")) {
      std::string v = env;
      std::transform(v.begin(), v.end(), v.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      enabled = (v == "bound" || v == "explicit_bound" || v == "nf");
    }
    init = 1;
  }
  return enabled;
}

enum class MxMode { Tile, Block, Vec1D };

static int8_t computeScaleFromMax(float max_val) {
  if (max_val <= 0.0f) return 0;
  int scale = static_cast<int>(std::floor(std::log2(max_val)));
  if (scale > 127) scale = 127;
  if (scale < -128) scale = -128;
  return static_cast<int8_t>(scale);
}

static int8_t computeScale(const float *data, size_t n) {
  float max_val = 0.0f;
  for (size_t i = 0; i < n; ++i) {
    max_val = std::fmax(max_val, std::fabs(data[i]));
  }
  if (max_val <= 0.0f) return 0;
  int scale = static_cast<int>(std::floor(std::log2(max_val)));
  if (scale > 127) scale = 127;
  if (scale < -128) scale = -128;
  return static_cast<int8_t>(scale);
}

static int computeScaleBits(const float *data, size_t n, int scale_bits) {
  if (scale_bits <= 0) {
    return static_cast<int>(computeScale(data, n));
  }
  float max_val = 0.0f;
  for (size_t i = 0; i < n; ++i) {
    max_val = std::fmax(max_val, std::fabs(data[i]));
  }
  if (max_val <= 0.0f) return 0;
  int scale = static_cast<int>(std::floor(std::log2(max_val)));
  const int max_scale = (1 << (scale_bits - 1)) - 1;
  const int min_scale = - (1 << (scale_bits - 1));
  if (scale > max_scale) scale = max_scale;
  if (scale < min_scale) scale = min_scale;
  return scale;
}

static int computeScaleBitsFromMax(float max_val, int scale_bits) {
  if (max_val <= 0.0f) return 0;
  int scale = static_cast<int>(std::floor(std::log2(max_val)));
  if (scale_bits <= 0) {
    if (scale > 127) scale = 127;
    if (scale < -128) scale = -128;
    return scale;
  }
  const int max_scale = (1 << (scale_bits - 1)) - 1;
  const int min_scale = -(1 << (scale_bits - 1));
  if (scale > max_scale) scale = max_scale;
  if (scale < min_scale) scale = min_scale;
  return scale;
}

static float pow2(int exp) { return std::ldexp(1.0f, exp); }

// Underflow policy for the host-side reference FP8 (E4M3 / E5M2) emulator.
// Default: flush-to-zero (matches legacy behaviour and the pre-OCP MX kernel).
// Set MX_FP8_SUBNORMAL=1 to switch to OCP-spec gradual underflow (subnormals).
static bool fp8_use_subnormals() {
  static int init = 0;
  static bool enabled = false;
  if (!init) {
    if (const char *env = getenv("MX_FP8_SUBNORMAL")) {
      enabled = (env[0] == '1' || env[0] == 'y' || env[0] == 'Y' ||
                 env[0] == 't' || env[0] == 'T');
    }
    init = 1;
  }
  return enabled;
}

static float quantizeFp(float x, int ebits, int mbits, float max_norm) {
  if (x == 0.0f) return 0.0f;
  const float sign = std::signbit(x) ? -1.0f : 1.0f;
  float ax = std::fabs(x);
  if (ax > max_norm) ax = max_norm;
  int exp = static_cast<int>(std::floor(std::log2(ax)));
  const int bias = (1 << (ebits - 1)) - 1;
  const int exp_enc = exp + bias;
  if (exp_enc <= 0) {
    if (!fp8_use_subnormals()) {
      return 0.0f;  // FTZ (default)
    }
    // OCP subnormal: encode as mant_int * 2^(1 - bias - mbits).
    // Smallest representable subnormal:
    //   E4M3 (bias=7, mbits=3) -> 2^-9 ≈ 1.95e-3
    //   E5M2 (bias=15, mbits=2) -> 2^-16 ≈ 1.53e-5
    // Round-to-nearest-even via std::lrint (with FE_TONEAREST default).
    // mant_int = 0 -> underflow to zero.
    // mant_int = 2^mbits -> "rounded up to smallest normal" — keep the natural
    // value (mi * subn_step = 2^(1-bias) = smallest normal). No cap needed; ax
    // is bounded above by 2^(1-bias) when exp_enc <= 0, so mi <= 2^mbits.
    const float subn_step = std::ldexp(1.0f, 1 - bias - mbits);
    int mi = static_cast<int>(std::lrint(ax / subn_step));
    if (mi <= 0) return 0.0f;
    return sign * static_cast<float>(mi) * subn_step;
  }
  const int exp_max = (1 << ebits) - 1;
  if (exp_enc >= exp_max) {
    return sign * max_norm;
  }
  const float base = pow2(exp);
  float mant = ax / base - 1.0f;
  const float step = 1.0f / static_cast<float>(1 << mbits);
  mant = std::round(mant / step) * step;
  float out = (1.0f + mant) * base;
  if (out > max_norm) out = max_norm;
  return sign * out;
}

static bool debug_tile_enabled() {
  static int init = 0;
  static bool enabled = false;
  if (!init) {
    const char *env = getenv("MX_DEBUG_TILE_DUMP");
    enabled = (env && env[0] == '1');
    init = 1;
  }
  return enabled;
}

static bool debug_tile_match(int row, int col) {
  static int init = 0;
  static int want_row = -1;
  static int want_col = -1;
  if (!init) {
    if (const char *env = getenv("MX_DEBUG_TILE_ROW")) {
      want_row = atoi(env);
    }
    if (const char *env = getenv("MX_DEBUG_TILE_COL")) {
      want_col = atoi(env);
    }
    init = 1;
  }
  if (want_row < 0 && want_col < 0) return true;
  if (want_row >= 0 && want_col >= 0) return (row == want_row && col == want_col);
  if (want_row >= 0) return row == want_row;
  return col == want_col;
}

static std::string fp32_bits_with_spaces(uint32_t bits) {
  std::string out;
  out.reserve(1 + 1 + 8 + 1 + 23);
  out.push_back(((bits >> 31) & 1) ? '1' : '0');
  out.push_back(' ');
  for (int i = 30; i >= 23; --i) {
    out.push_back(((bits >> i) & 1) ? '1' : '0');
  }
  out.push_back(' ');
  for (int i = 22; i >= 0; --i) {
    out.push_back(((bits >> i) & 1) ? '1' : '0');
  }
  return out;
}

static void debug_dump_sample(const char *label, const float *data, size_t n) {
  size_t limit = n < 8 ? n : 8;
  std::cout << "[DEBUG] " << label << " sample (first " << limit << ")" << std::endl;
  for (size_t i = 0; i < limit; ++i) {
    uint32_t bits = 0;
    std::memcpy(&bits, &data[i], sizeof(bits));
    std::cout << "  [" << i << "] " << data[i]
              << " bits=" << fp32_bits_with_spaces(bits) << std::endl;
  }
}

static std::string debug_dump_path() {
  static int init = 0;
  static std::string path;
  if (!init) {
    if (const char *env = getenv("MX_DEBUG_TILE_DUMP_FILE")) {
      path = env;
    }
    init = 1;
  }
  return path;
}

static bool debug_dump_fp8() {
  static int init = 0;
  static bool enabled = false;
  if (!init) {
    const char *env = getenv("MX_DEBUG_TILE_DUMP_FP8");
    enabled = (env && env[0] == '1');
    init = 1;
  }
  return enabled;
}

static void debug_dump_tile_to_file(const char *label,
                                    int row, int col,
                                    const float *pre,
                                    const float *post,
                                    size_t m, size_t n) {
  const std::string path = debug_dump_path();
  if (path.empty()) return;

  const bool include_fp8 = debug_dump_fp8();
  std::ofstream out(path, std::ios::app);
  if (!out) return;

  static bool header_written = false;
  if (!header_written) {
    out << "label,row,col,r,c,idx,pre,pre_bits,post,post_bits";
    if (include_fp8) {
      out << ",fp8_e4m3,fp8_e4m3_bits,fp8_e5m2,fp8_e5m2_bits";
    }
    out << "\n";
    header_written = true;
  }

  for (size_t r = 0; r < m; ++r) {
    for (size_t c = 0; c < n; ++c) {
      const size_t idx = r * n + c;
      uint32_t pre_bits = 0;
      uint32_t post_bits = 0;
      std::memcpy(&pre_bits, &pre[idx], sizeof(pre_bits));
      std::memcpy(&post_bits, &post[idx], sizeof(post_bits));
      out << label << "," << row << "," << col << "," << r << "," << c
          << "," << idx << "," << pre[idx] << ","
          << fp32_bits_with_spaces(pre_bits) << "," << post[idx] << ","
          << fp32_bits_with_spaces(post_bits);
      if (include_fp8) {
        const float fp8_e4m3 = quantizeFp(pre[idx], 4, 3, 448.0f);
        const float fp8_e5m2 = quantizeFp(pre[idx], 5, 2, 57344.0f);
        uint32_t fp8_e4m3_bits = 0;
        uint32_t fp8_e5m2_bits = 0;
        std::memcpy(&fp8_e4m3_bits, &fp8_e4m3, sizeof(fp8_e4m3_bits));
        std::memcpy(&fp8_e5m2_bits, &fp8_e5m2, sizeof(fp8_e5m2_bits));
        out << "," << fp8_e4m3 << "," << fp32_bits_with_spaces(fp8_e4m3_bits)
            << "," << fp8_e5m2 << "," << fp32_bits_with_spaces(fp8_e5m2_bits);
      }
      out << "\n";
    }
  }
}

static MxMode mx_mode() {
  static int init = 0;
  static MxMode mode = MxMode::Tile;
  if (!init) {
    if (const char *env = getenv("MX_MX_MODE")) {
      std::string v = env;
      std::transform(v.begin(), v.end(), v.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      if (v == "block" || v == "blocked") {
        mode = MxMode::Block;
      } else if (v == "vec1d" || v == "vector" || v == "vec") {
        mode = MxMode::Vec1D;
      }
    }
    init = 1;
  }
  return mode;
}

// Optional per-format override of MxMode for "wider" formats (mx_fp16 / mx_fp32).
// Set MX_MODE_FP16={tile|block|vec1d} to force a different mode just for those
// tiles, while leaving the smaller MX formats (e2m1/e4m3/e5m2/...) on the
// globally-selected MX_MX_MODE. Returns true and writes *out if an override
// applies for (ebits, mbits); otherwise returns false.
static bool mx_mode_override_for_format(int ebits, int mbits, MxMode *out) {
  static int init = 0;
  static bool have_override = false;
  static MxMode override_mode = MxMode::Tile;
  if (!init) {
    if (const char *env = getenv("MX_MODE_FP16")) {
      std::string v = env;
      std::transform(v.begin(), v.end(), v.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      if (v == "tile") { override_mode = MxMode::Tile; have_override = true; }
      else if (v == "block" || v == "blocked") { override_mode = MxMode::Block; have_override = true; }
      else if (v == "vec1d" || v == "vector" || v == "vec") { override_mode = MxMode::Vec1D; have_override = true; }
    }
    init = 1;
  }
  if (!have_override) return false;
  // Apply to mx_fp16 (ebits=5, mbits=10) and mx_fp32 (ebits=8, mbits=23).
  const bool is_fp16 = (ebits == 5 && mbits == 10);
  const bool is_fp32 = (ebits == 8 && mbits == 23);
  if (!(is_fp16 || is_fp32)) return false;
  if (out) *out = override_mode;
  return true;
}

static int mx_block_subtile() {
  static int init = 0;
  static int subtile = 0;
  if (!init) {
    if (const char *env = getenv("MX_BLOCK_SUBTILE")) {
      int v = atoi(env);
      subtile = v > 0 ? v : 0;
    }
    init = 1;
  }
  return subtile;
}

#ifdef __CUDACC__
// Simple SPD fill on GPU: A = alpha * I + beta * v v^T with v_i = 0.001*(i+1)
// Written in row-major layout.
__global__ void fill_spd_tile(float *dst, int ld, int m, int n, float alpha,
                              float beta) {
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  if (row < m && col < n) {
    float vi = 0.001f * (row + 1);
    float vj = 0.001f * (col + 1);
    float val = beta * vi * vj;
    if (row == col) val += alpha;
    dst[row * ld + col] = val;
  }
}
#endif

static bool is_device_or_managed(const void *p) {
  if (!p) return false;
  cudaPointerAttributes attr;
  if (cudaPointerGetAttributes(&attr, p) != cudaSuccess) return false;
  return attr.type == cudaMemoryTypeDevice || attr.type == cudaMemoryTypeManaged;
}

static const char *dtype_to_string(int dtype) {
  switch (dtype) {
    case CUDA_R_64F:
      return "fp64";
    case CUDA_R_32F:
      return "fp32";
    case CUDA_R_16F:
      return "fp16";
    case CUDA_R_16BF:
      return "bf16";
    case CUDA_R_8F_E4M3:
      return "fp8_e4m3";
    default:
      return "unknown";
  }
}

struct QuantWorkspace {
  float *d_in = nullptr;
  float *d_out = nullptr;
  float *d_max = nullptr;
  size_t capacity = 0;  // elements
  size_t max_capacity = 0;  // max-values elements
};

static QuantWorkspace g_quant_ws;
static size_t g_pinned_bytes = 0;
static size_t g_pinned_tiles = 0;
static size_t g_quant_ws_bytes = 0;

static void ensureQuantWorkspace(size_t elems, size_t max_elems) {
  if (g_quant_ws.capacity >= elems && g_quant_ws.max_capacity >= max_elems) {
    return;
  }
  if (g_quant_ws.d_in) {
    cudaFree(g_quant_ws.d_in);
    cudaFree(g_quant_ws.d_out);
    cudaFree(g_quant_ws.d_max);
  }
  g_quant_ws.capacity = elems;
  g_quant_ws.max_capacity = max_elems;
  cudaMalloc(reinterpret_cast<void **>(&g_quant_ws.d_in), elems * sizeof(float));
  cudaMalloc(reinterpret_cast<void **>(&g_quant_ws.d_out), elems * sizeof(float));
  cudaMalloc(reinterpret_cast<void **>(&g_quant_ws.d_max), max_elems * sizeof(float));
  g_quant_ws_bytes = elems * sizeof(float) * 2 + max_elems * sizeof(float);
}
#define FMULS_POTRF_ULL(__n) (__n * __n / 2 * (__n / 3 + 1) + __n / 3)
#define FADDS_POTRF_ULL(__n) (__n * ((__n * __n - 1) / 6))
#define FLOPS_DPOTRF_ULL(__n) (FMULS_POTRF_ULL((__n)) + FADDS_POTRF_ULL((__n)))

/***************************************************************************/
/**
 *
 * @ingroup double
 *
 *  PLASMA_dpotrf - Computes the Cholesky factorization of a symmetric positive
 *definite (or Hermitian positive definite in the complex case) matrix A. The
 *factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = PlasmaLower}^{U^H\times U, if uplo =
 *PlasmaUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric positive definite (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly
 *lower triangular part of A is not referenced. If UPLO = 'L', the leading
 *N-by-N lower triangular part of A contains the lower triangular part of the
 *matrix A, and the strictly upper triangular part of A is not referenced. On
 *exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U**T*U or A = L*L**T.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, the leading minor of order i of A is not positive
 *definite, so the factorization could not be completed, and the solution has
 *not been computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dpotrf_Tile
 * @sa PLASMA_dpotrf_Tile_Async
 * @sa PLASMA_cpotrf
 * @sa PLASMA_dpotrf
 * @sa PLASMA_spotrf
 * @sa PLASMA_dpotrs
 *
 ******************************************************************************/
void printMixedPrecisionTiledArray(MixedPrecisionTiledArray *array);

void makeMixedPrecisionTiledArray(MixedPrecisionTiledArray *array,
                                  const int uplo, const double *data,
                                  const size_t m, const size_t n,
                                  const size_t ld, const size_t nb) {

    // DEBUG CHECK: Confirms the function is being executed
    std::cout << "\n[DEBUG] --- ENTERING makeMixedPrecisionTiledArray ---" << std::endl; 
  g_pinned_tiles = 0;
  g_pinned_bytes = 0;
    
    array->nb = nb;
    array->mb = nb;
    array->m = array->n = n;
    array->nt = ((long)n - 1) / nb + 1;
  std::cout << "[DEBUG] nb used: " << nb << ", nt: " << array->nt << std::endl;
    array->tiles = new MixedPrecisionTile[array->nt * array->nt]();
    
    // Assuming Eigen headers are included:
    // If GPU tiling is requested, build tiles directly on GPU and skip
    // column-major mapping. Otherwise, keep the existing host path.
    if (use_gpu_tiling()) {
  #ifdef __CUDACC__
      array->uplo = uplo;
      if (uplo == PlasmaLower) {
        dim3 block(16, 16);
        for (int row = 0; row < array->nt; ++row) {
          for (int col = 0; col <= row; ++col) {
            auto &tile = array->tiles[col + row * array->nt];
            size_t nbRow = row == array->nt - 1 ? array->m - row * nb : nb;
            size_t nbCol = col == array->nt - 1 ? array->n - col * nb : nb;
            tile.m = nbRow;
            tile.n = nbCol;
            tile.ld = nbCol;
            tile.layout = CblasRowMajor;
            tile.dtype = CUDA_R_32F; // GPU path currently stores FP32 tiles
            cudaMalloc(reinterpret_cast<void **>(&tile.data),
                   sizeof(float) * tile.m * tile.n);
            dim3 grid((nbCol + block.x - 1) / block.x,
                  (nbRow + block.y - 1) / block.y);
            fill_spd_tile<<<grid, block>>>(static_cast<float *>(tile.data),
                             static_cast<int>(tile.ld),
                             static_cast<int>(nbRow),
                             static_cast<int>(nbCol),
                             1.0f, 0.01f);
          }
        }
        cudaDeviceSynchronize();
      } else {
        std::cout << "Only Lower matrix is supported so far" << std::endl;
      }
      return;
  #else
      std::cout << "[WARN] MX_GPU_TILING=1 requires CUDA compilation; falling back to host tiling." << std::endl;
  #endif
    }

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
           Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, 1>>
      mappedA{data,
          static_cast<Eigen::Index>(m),
          static_cast<Eigen::Index>(n),
          {static_cast<Eigen::Index>(ld), 1}};

    long double normA = mappedA.norm();
    array->uplo = uplo;

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties_v2(&deviceProp, 0);
    const bool hasFp8 = deviceProp.major >= 9;

    // Assuming CUDA headers are included:
    if (uplo == PlasmaLower) {
      for (int row = 0; row < array->nt; ++row) {
        for (int col = 0; col <= row; ++col) {
          auto &tile = array->tiles[col + row * array->nt];
          size_t nbRow = row == array->nt - 1 ? array->m - row * nb : nb;
          size_t nbCol = col == array->nt - 1 ? array->n - col * nb : nb;
          tile.m = nbRow;
          tile.n = nbCol;
          tile.ld = nbCol;
          tile.layout = CblasRowMajor;       
          auto mappedBlock = mappedA.block(row * nb, col * nb, nbRow, nbCol);
          long double normTile = mappedBlock.norm();
          auto epsilonRatio = array->nt * normTile / normA;

          const long double sourceEpsilon = []() {
            const char *env = getenv("MX_SOURCE_EPSILON");
            if (!env || !env[0]) {
              return static_cast<long double>(
                  std::numeric_limits<float>::epsilon());
            }
            char *end = nullptr;
            const long double v = std::strtold(env, &end);
            if (end == env || *end != '\0' || !(v > 0.0L)) {
              std::cout << "[WARN] Invalid MX_SOURCE_EPSILON='" << env
                        << "', using default float epsilon." << std::endl;
              return static_cast<long double>(
                  std::numeric_limits<float>::epsilon());
            }
            return v;
          }();
          static const int plain_band = []() {
            const char *env = getenv("MX_PLAIN_BAND");
            if (!env) return 1;  // default: one super-/sub-diagonal stays plain FP32
            int v = atoi(env);
            return v < 0 ? 0 : v;
          }();
          const std::string fmt = forced_format();
          const bool has_forced = !fmt.empty();
          bool force_fp32_all = false;
          std::string target_low;
          const std::string fp32_bucket = bucket_format("MX_BUCKET_FP32");
          const std::string fp16_bucket = bucket_format("MX_BUCKET_FP16");
          static const int fp32_scale_bits = []() {
            const char *env = getenv("MX_FP32_SCALE_BITS");
            if (!env) return 11;
            int v = atoi(env);
            return v > 0 ? v : 11;
          }();
          if (has_forced) {
            if (fmt == "fp32") {
              force_fp32_all = true;
            } else if (fmt == "fp16") {
              target_low = "fp16";
            } else if (fmt == "bf16") {
              target_low = "bf16";
            } else if (fmt == "mx_fp16" || fmt == "mx_f16") {
              target_low = "mx_fp16";
            } else if (fmt == "mx_fp32" || fmt == "mx_f32") {
              target_low = "mx_fp32";
            } else if (fmt == "mx_e4m3" || fmt == "mx_fp8_e4m3" || fmt == "e4m3") {
              target_low = "mx_e4m3";
            } else if (fmt == "mx_e5m2" || fmt == "mx_fp8_e5m2" || fmt == "e5m2") {
              target_low = "mx_e5m2";
            } else if (fmt == "fp8_e4m3" || fmt == "fp8e4m3") {
              target_low = "fp8_e4m3";
            } else if (fmt == "fp8_e5m2" || fmt == "fp8e5m2") {
              target_low = "fp8_e5m2";
            } else if (fmt == "e3m2" || fmt == "fp6_e3m2" || fmt == "fp6e3m2") {
              target_low = "e3m2";
            } else if (fmt == "e2m3" || fmt == "fp6_e2m3" || fmt == "fp6e2m3") {
              target_low = "e2m3";
            } else if (fmt == "e2m1" || fmt == "fp4_e2m1" || fmt == "fp4e2m1") {
              target_low = "e2m1";
            } else {
              std::cout << "[WARN] Unknown MX_FORCE_FORMAT='" << fmt
                        << "', falling back to heuristic." << std::endl;
            }
          }
          const bool allow_low = !target_low.empty() || hasFp8;
          const std::string selected_low_format =
              target_low.empty() ? std::string("mx_e4m3") : target_low;
          static const int scale_aware_epsilon = []() {
            const char *env = getenv("MX_SCALE_AWARE_EPSILON");
            if (!env) return 1;
            return (env[0] == '0') ? 0 : 1;
          }();
          long double lowUnitRoundoff = format_unit_roundoff(selected_low_format);
          if (scale_aware_epsilon && format_uses_shared_scale(selected_low_format)) {
            lowUnitRoundoff += format_scaled_underflow_term(selected_low_format);
          }
          if (!(lowUnitRoundoff > 0.0L)) {
            lowUnitRoundoff = 1.0L / 16.0L;
          }
          const long double lowCutoff = sourceEpsilon / lowUnitRoundoff;

          auto apply_mx_quant = [&](const char *label, int ebits, int mbits,
                                    float max_norm) {
            std::cout << "--- Tile (" << row << ", " << col
                      << ") selected for CUDA_R_32F (" << label
                      << " Quantization) ---" << std::endl;

            tile.dtype = CUDA_R_32F;
            tile.mx_ebits = ebits;
            tile.mx_mbits = mbits;
            tile.mx_max_norm = max_norm;
            tile.mx_scale_bits = 0;
            MxMode cur_mode = mx_mode();
            // Per-format override: e.g. MX_MODE_FP16=tile keeps mx_fp16/mx_fp32
            // tiles on a single shared scale even when the global mode is vec1d.
            mx_mode_override_for_format(ebits, mbits, &cur_mode);
            const bool block_mode  = (cur_mode == MxMode::Block);
            const bool vec1d_mode  = (cur_mode == MxMode::Vec1D);
            tile.mx_mode_block = block_mode ? 1 : (vec1d_mode ? 2 : 0);
            tile.mx_block_subtile = (block_mode || vec1d_mode) ? mx_block_subtile() : 0;

            size_t tile_size = tile.m * tile.n;
            float *h_output = nullptr;

            cudaMallocHost((void **)&tile.data, sizeof(uint32_t) * tile_size);
            h_output = static_cast<float *>(tile.data);
            g_pinned_tiles++;
            g_pinned_bytes += sizeof(uint32_t) * tile_size;

            std::vector<float> tmp(tile_size);
            Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic,
                         Eigen::RowMajor>>{
              tmp.data(),
              static_cast<Eigen::Index>(tile.m),
              static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<float>();

            const long total_elements = static_cast<long>(tile_size);
            float global_max = 0.0f;
            for (long i = 0; i < total_elements; ++i) {
              global_max = std::fmax(global_max, std::fabs(tmp[i]));
            }

            std::cout << "  Tile Size: " << total_elements
                      << ", MaxVals: " << global_max
                      << ", Mode: "
                      << (block_mode ? "block" : (vec1d_mode ? "vec1d" : "tile")) << std::endl;

            const bool mx_fp16 = (ebits == 5 && mbits == 10);
            const bool debug_tile = debug_tile_enabled() && debug_tile_match(row, col);
            std::vector<float> pre_tile;
            if (debug_tile) {
              pre_tile.assign(tmp.begin(), tmp.end());
              debug_dump_sample("pre-quant", tmp.data(), tile_size);
            }
            if (block_mode) {
              const int subtile = mx_block_subtile();
              if (subtile > 0) {
                const size_t tile_m = static_cast<size_t>(tile.m);
                const size_t tile_n = static_cast<size_t>(tile.n);
                for (size_t r0 = 0; r0 < tile_m; r0 += static_cast<size_t>(subtile)) {
                  const size_t r_max = (static_cast<size_t>(subtile) < (tile_m - r0))
                                         ? static_cast<size_t>(subtile)
                                         : (tile_m - r0);
                  for (size_t c0 = 0; c0 < tile_n; c0 += static_cast<size_t>(subtile)) {
                    const size_t c_max = (static_cast<size_t>(subtile) < (tile_n - c0))
                                           ? static_cast<size_t>(subtile)
                                           : (tile_n - c0);
                    float max_val = 0.0f;
                    for (size_t r = 0; r < r_max; ++r) {
                      const float *rp = tmp.data() + (r0 + r) * tile_n + c0;
                      for (size_t c = 0; c < c_max; ++c)
                        max_val = std::fmax(max_val, std::fabs(rp[c]));
                    }
                    const int8_t st_scale = computeScaleFromMax(max_val);
                    const float scale = pow2(static_cast<int>(st_scale));
                    const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
                    for (size_t r = 0; r < r_max; ++r) {
                      float *rp = tmp.data() + (r0 + r) * tile_n + c0;
                      for (size_t c = 0; c < c_max; ++c) {
                        const float x = rp[c] * inv_scale;
                        if (mx_fp16) {
                          rp[c] = static_cast<float>(static_cast<Eigen::half>(x)) * scale;
                        } else {
                          rp[c] = quantizeFp(x, ebits, mbits, max_norm) * scale;
                        }
                      }
                    }
                  }
                }
              } else {
                for (size_t r = 0; r < static_cast<size_t>(tile.m); ++r) {
                  float *rp = tmp.data() + r * static_cast<size_t>(tile.n);
                  const int8_t row_scale = computeScale(rp, tile.n);
                  const float scale = pow2(static_cast<int>(row_scale));
                  const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
                  for (size_t c = 0; c < static_cast<size_t>(tile.n); ++c) {
                    const float x = rp[c] * inv_scale;
                    if (mx_fp16) {
                      rp[c] = static_cast<float>(static_cast<Eigen::half>(x)) * scale;
                    } else {
                      rp[c] = quantizeFp(x, ebits, mbits, max_norm) * scale;
                    }
                  }
                }
              }
            } else if (vec1d_mode) {
              // 1D row-vector: each row split into vectors of `vec_sz` elements, each with own scale
              const size_t vec_sz = (mx_block_subtile() > 0)
                                      ? static_cast<size_t>(mx_block_subtile())
                                      : 32;
              const size_t tile_m = static_cast<size_t>(tile.m);
              const size_t tile_n = static_cast<size_t>(tile.n);
              for (size_t r = 0; r < tile_m; ++r) {
                float *rp = tmp.data() + r * tile_n;
                for (size_t c0 = 0; c0 < tile_n; c0 += vec_sz) {
                  const size_t c_end = (c0 + vec_sz < tile_n) ? c0 + vec_sz : tile_n;
                  float max_val = 0.0f;
                  for (size_t c = c0; c < c_end; ++c)
                    max_val = std::fmax(max_val, std::fabs(rp[c]));
                  const int8_t vec_scale = computeScaleFromMax(max_val);
                  const float scale = pow2(static_cast<int>(vec_scale));
                  const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
                  for (size_t c = c0; c < c_end; ++c) {
                    const float x = rp[c] * inv_scale;
                    if (mx_fp16)
                      rp[c] = static_cast<float>(static_cast<Eigen::half>(x)) * scale;
                    else
                      rp[c] = quantizeFp(x, ebits, mbits, max_norm) * scale;
                  }
                }
              }
            } else {
              const int8_t tile_scale = computeScale(tmp.data(), tile_size);
              const float scale = pow2(static_cast<int>(tile_scale));
              const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
              for (size_t i = 0; i < tile_size; ++i) {
                const float x = tmp[i] * inv_scale;
                if (mx_fp16) {
                  tmp[i] = static_cast<float>(static_cast<Eigen::half>(x)) * scale;
                } else {
                  tmp[i] = quantizeFp(x, ebits, mbits, max_norm) * scale;
                }
              }
            }
            for (size_t i = 0; i < tile_size; ++i)
              h_output[i] = tmp[i];
            if (debug_tile) {
              debug_dump_sample("post-quant", tmp.data(), tile_size);
              debug_dump_tile_to_file(label, row, col,
                                       pre_tile.data(), tmp.data(),
                                       static_cast<size_t>(tile.m),
                                       static_cast<size_t>(tile.n));
            }
          };
          auto apply_mx_quant_fp64 = [&](const char *label, int ebits, int mbits,
                                         float max_norm, int scale_bits) {
            std::cout << "--- Tile (" << row << ", " << col
                      << ") selected for CUDA_R_64F (" << label
                      << " Quantization) ---" << std::endl;

            tile.dtype = CUDA_R_64F;
            tile.mx_ebits = ebits;
            tile.mx_mbits = mbits;
            tile.mx_max_norm = max_norm;
            tile.mx_scale_bits = scale_bits;
            const bool block_mode_fp64 = (mx_mode() == MxMode::Block);
            const bool vec1d_mode_fp64 = (mx_mode() == MxMode::Vec1D);
            tile.mx_mode_block = block_mode_fp64 ? 1 : (vec1d_mode_fp64 ? 2 : 0);
            tile.mx_block_subtile = (block_mode_fp64 || vec1d_mode_fp64) ? mx_block_subtile() : 0;

            size_t tile_size = tile.m * tile.n;
            double *h_output = nullptr;

            cudaMallocHost((void **)&tile.data, sizeof(double) * tile_size);
            h_output = static_cast<double *>(tile.data);
            g_pinned_tiles++;
            g_pinned_bytes += sizeof(double) * tile_size;

            std::vector<float> tmp(tile_size);
            Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic,
                         Eigen::RowMajor>>{
              tmp.data(),
              static_cast<Eigen::Index>(tile.m),
              static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<float>();

            const long total_elements = static_cast<long>(tile_size);
            const bool block_mode = (mx_mode() == MxMode::Block);
            const bool vec1d_mode = (mx_mode() == MxMode::Vec1D);

            const bool debug_tile = debug_tile_enabled() && debug_tile_match(row, col);
            std::vector<float> pre_tile;
            if (debug_tile) {
              pre_tile.assign(tmp.begin(), tmp.end());
              debug_dump_sample("pre-quant", tmp.data(), tile_size);
            }

            if (block_mode) {
              const int subtile = mx_block_subtile();
              if (subtile > 0) {
                const size_t tile_m = static_cast<size_t>(tile.m);
                const size_t tile_n = static_cast<size_t>(tile.n);
                for (size_t r0 = 0; r0 < tile_m; r0 += static_cast<size_t>(subtile)) {
                  const size_t r_max = (static_cast<size_t>(subtile) < (tile_m - r0))
                                         ? static_cast<size_t>(subtile)
                                         : (tile_m - r0);
                  for (size_t c0 = 0; c0 < tile_n; c0 += static_cast<size_t>(subtile)) {
                    const size_t c_max = (static_cast<size_t>(subtile) < (tile_n - c0))
                                           ? static_cast<size_t>(subtile)
                                           : (tile_n - c0);
                    float max_val = 0.0f;
                    for (size_t r = 0; r < r_max; ++r) {
                      const float *row = tmp.data() + (r0 + r) * tile_n + c0;
                      for (size_t c = 0; c < c_max; ++c) {
                        max_val = std::fmax(max_val, std::fabs(row[c]));
                      }
                    }
                    const int st_scale = computeScaleBitsFromMax(max_val, scale_bits);
                    const float scale = pow2(st_scale);
                    const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
                    for (size_t r = 0; r < r_max; ++r) {
                      float *row = tmp.data() + (r0 + r) * tile_n + c0;
                      for (size_t c = 0; c < c_max; ++c) {
                        const float x = row[c] * inv_scale;
                        const float q = quantizeFp(x, ebits, mbits, max_norm);
                        row[c] = q * scale;
                      }
                    }
                  }
                }
              } else {
                for (size_t r = 0; r < static_cast<size_t>(tile.m); ++r) {
                  float *row = tmp.data() + r * static_cast<size_t>(tile.n);
                  const int row_scale = computeScaleBits(row, tile.n, scale_bits);
                  const float scale = pow2(row_scale);
                  const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
                  for (size_t c = 0; c < static_cast<size_t>(tile.n); ++c) {
                    const float x = row[c] * inv_scale;
                    const float q = quantizeFp(x, ebits, mbits, max_norm);
                    row[c] = q * scale;
                  }
                }
              }
            } else if (vec1d_mode) {
              // 1D row-vector mode for FP32-emulated MX: each row split into
              // groups of `vec_sz` elements, each with its own scale_bits-wide
              // shared scale (default vec_sz=32, matching canonical MX).
              const size_t vec_sz = (mx_block_subtile() > 0)
                                      ? static_cast<size_t>(mx_block_subtile())
                                      : 32;
              const size_t tile_m_ = static_cast<size_t>(tile.m);
              const size_t tile_n_ = static_cast<size_t>(tile.n);
              for (size_t r = 0; r < tile_m_; ++r) {
                float *rp = tmp.data() + r * tile_n_;
                for (size_t c0 = 0; c0 < tile_n_; c0 += vec_sz) {
                  const size_t c_end = (c0 + vec_sz < tile_n_) ? c0 + vec_sz : tile_n_;
                  float max_val = 0.0f;
                  for (size_t c = c0; c < c_end; ++c)
                    max_val = std::fmax(max_val, std::fabs(rp[c]));
                  const int vec_scale = computeScaleBitsFromMax(max_val, scale_bits);
                  const float scale = pow2(vec_scale);
                  const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
                  for (size_t c = c0; c < c_end; ++c) {
                    const float x = rp[c] * inv_scale;
                    const float q = quantizeFp(x, ebits, mbits, max_norm);
                    rp[c] = q * scale;
                  }
                }
              }
            } else {
              const int tile_scale = computeScaleBits(tmp.data(), tile_size, scale_bits);
              const float scale = pow2(tile_scale);
              const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
              for (size_t i = 0; i < tile_size; ++i) {
                const float x = tmp[i] * inv_scale;
                const float q = quantizeFp(x, ebits, mbits, max_norm);
                tmp[i] = q * scale;
              }
            }

            for (size_t i = 0; i < tile_size; ++i) {
              h_output[i] = static_cast<double>(tmp[i]);
            }

            if (debug_tile) {
              debug_dump_sample("post-quant", tmp.data(), tile_size);
              debug_dump_tile_to_file(label, row, col,
                                       pre_tile.data(), tmp.data(),
                                       static_cast<size_t>(tile.m),
                                       static_cast<size_t>(tile.n));
            }
          };
          auto apply_plain_fp_quant = [&](const char *label, int ebits, int mbits,
                                          float max_norm) {
            std::cout << "--- Tile (" << row << ", " << col
                      << ") selected for CUDA_R_32F (" << label
                      << " Quantization) ---" << std::endl;

            tile.dtype = CUDA_R_32F;
            // mx_scale_bits = -1 signals plain FP8 (no MX block scale) to requantizeTileHost
            tile.mx_ebits = ebits;
            tile.mx_mbits = mbits;
            tile.mx_max_norm = max_norm;
            tile.mx_scale_bits = -1;
            tile.mx_mode_block = 0;
            tile.mx_block_subtile = 0;

            size_t tile_size = tile.m * tile.n;
            float *h_output = nullptr;

            cudaMallocHost((void **)&tile.data, sizeof(uint32_t) * tile_size);
            h_output = static_cast<float *>(tile.data);
            g_pinned_tiles++;
            g_pinned_bytes += sizeof(uint32_t) * tile_size;

            Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic,
                         Eigen::RowMajor>>{
              h_output,
              static_cast<Eigen::Index>(tile.m),
              static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<float>();

            const long total_elements = static_cast<long>(tile_size);
            const bool debug_tile = debug_tile_enabled() && debug_tile_match(row, col);
            std::vector<float> pre_tile;
            if (debug_tile) {
              pre_tile.assign(h_output, h_output + tile_size);
            }
            if (debug_tile) {
              debug_dump_sample("pre-quant", h_output, tile_size);
            }
            for (long i = 0; i < total_elements; ++i) {
              h_output[i] = quantizeFp(h_output[i], ebits, mbits, max_norm);
            }
            if (debug_tile) {
              debug_dump_sample("post-quant", h_output, tile_size);
              debug_dump_tile_to_file(label, row, col,
                                       pre_tile.data(), h_output,
                                       static_cast<size_t>(tile.m),
                                       static_cast<size_t>(tile.n));
            }
          };
          auto log_tile_target = [&](const char *target) {
            std::cout << "[TILE_TARGET] (" << row << ", " << col << ") "
                      << target << std::endl;
          };

          auto handle_bucket_value = [&](const std::string &v) -> bool {
            if (v == "fp32") {
              log_tile_target("fp32");
              tile.dtype = CUDA_R_32F;
              cudaMallocHost(&tile.data, sizeof(uint32_t) * tile.m * tile.n);
              Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic,
                                      Eigen::RowMajor>>{
                  static_cast<float *>(tile.data),
                  static_cast<Eigen::Index>(tile.m),
                  static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<float>();
              return true;
            }
            if (v == "fp16") {
              log_tile_target("fp16");
              tile.dtype = CUDA_R_16F;
              cudaMallocHost(&tile.data, sizeof(uint16_t) * tile.m * tile.n);
              Eigen::Map<Eigen::Matrix<Eigen::half, Eigen::Dynamic, Eigen::Dynamic,
                                      Eigen::RowMajor>>{
                  static_cast<Eigen::half *>(tile.data),
                  static_cast<Eigen::Index>(tile.m),
                  static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<Eigen::half>();
              return true;
            }
            if (v == "bf16") {
              log_tile_target("bf16");
              tile.dtype = CUDA_R_16BF;
              cudaMallocHost(&tile.data, sizeof(uint16_t) * tile.m * tile.n);
              Eigen::Map<Eigen::Matrix<Eigen::bfloat16, Eigen::Dynamic, Eigen::Dynamic,
                                      Eigen::RowMajor>>{
                  static_cast<Eigen::bfloat16 *>(tile.data),
                  static_cast<Eigen::Index>(tile.m),
                  static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<Eigen::bfloat16>();
              return true;
            }
            if (v == "mx_fp16") {
              log_tile_target("mx_fp16");
              apply_mx_quant("MX_FP16", 5, 10, 65504.0f);
              return true;
            }
            if (v == "mx_fp32" || v == "mx_f32") {
              log_tile_target("mx_fp32");
              apply_mx_quant_fp64("MX_FP32", 8, 23,
                                  (std::numeric_limits<float>::max)(),
                                  fp32_scale_bits);
              return true;
            }
            if (v == "mx_e4m3" || v == "e4m3") {
              log_tile_target("mx_e4m3");
              apply_mx_quant("MX_E4M3", 4, 3, 448.0f);
              return true;
            }
            if (v == "mx_e5m2" || v == "e5m2") {
              log_tile_target("mx_e5m2");
              apply_mx_quant("MX_E5M2", 5, 2, 57344.0f);
              return true;
            }
            if (v == "fp8_e4m3" || v == "fp8e4m3") {
              log_tile_target("fp8_e4m3");
              apply_plain_fp_quant("FP8_E4M3", 4, 3, 448.0f);
              return true;
            }
            if (v == "fp8_e5m2" || v == "fp8e5m2") {
              log_tile_target("fp8_e5m2");
              apply_plain_fp_quant("FP8_E5M2", 5, 2, 57344.0f);
              return true;
            }
            if (v == "e3m2") {
              log_tile_target("e3m2");
              apply_mx_quant("MX_E3M2", 3, 2, 28.0f);
              return true;
            }
            if (v == "e2m3") {
              log_tile_target("e2m3");
              apply_mx_quant("MX_E2M3", 2, 3, 7.5f);
              return true;
            }
            if (v == "e2m1") {
              log_tile_target("e2m1");
              apply_mx_quant("MX_E2M1", 2, 1, 6.0f);
              return true;
            }
            return false;
          };

          auto apply_bucket_target = [&](const std::string &bucket,
                                         const char *fallback) {
            std::string v = bucket.empty() ? std::string(fallback) : bucket;
            if (!handle_bucket_value(v)) {
              std::cout << "[WARN] Unknown bucket format '" << v
                        << "', using " << fallback << std::endl;
              handle_bucket_value(fallback);
            }
          };

          auto canonical_format_name = [&](const std::string &raw,
                                           const char *fallback) {
            std::string v = raw.empty() ? std::string(fallback) : raw;
            std::transform(v.begin(), v.end(), v.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (v == "mx_f32") v = "mx_fp32";
            if (v == "mx_f16") v = "mx_fp16";
            if (v == "mx_fp8_e4m3") v = "mx_e4m3";
            if (v == "mx_fp8_e5m2") v = "mx_e5m2";
            if (v == "fp8e4m3") v = "fp8_e4m3";
            if (v == "fp8e5m2") v = "fp8_e5m2";
            if (v == "fp6_e3m2" || v == "fp6e3m2") v = "e3m2";
            if (v == "fp6_e2m3" || v == "fp6e2m3") v = "e2m3";
            if (v == "fp4_e2m1" || v == "fp4e2m1") v = "e2m1";

            if (v == "fp32" || v == "fp16" || v == "bf16" ||
                v == "mx_fp16" || v == "mx_fp32" ||
                v == "mx_e4m3" || v == "mx_e5m2" ||
                v == "fp8_e4m3" || v == "fp8_e5m2" ||
                v == "e3m2" || v == "e2m3" || v == "e2m1") {
              return v;
            }
            return std::string(fallback);
          };

          static const int bound_debug = []() {
            const char *env = getenv("MX_BOUND_DEBUG");
            if (!env) return 0;
            return (env[0] == '1') ? 1 : 0;
          }();
          static const int bound_cheap_first = []() {
            const char *env = getenv("MX_BOUND_RANKING");
            if (!env) return 0;  // default keeps legacy fp32->fp16->low behavior
            std::string v = env;
            std::transform(v.begin(), v.end(), v.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return (v == "cheap" || v == "cheap_first" || v == "ascending_cost") ? 1 : 0;
          }();
          enum class BoundLadderMode { Legacy, Full, IeeeOnly };

          static const BoundLadderMode bound_ladder_mode = []() {
            const char *env = getenv("MX_BOUND_LADDER");
            if (!env) {
              // default: enable full low->high ladder in bound mode
              return BoundLadderMode::Full;
            }
            std::string v = env;
            std::transform(v.begin(), v.end(), v.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (v == "full" || v == "all" || v == "1" || v == "on") {
              return BoundLadderMode::Full;
            }
            if (v == "ieee" || v == "ieee_only" || v == "plain_ieee" ||
                v == "ieee_fp" || v == "ieee_fp_only") {
              return BoundLadderMode::IeeeOnly;
            }
            return BoundLadderMode::Legacy;
          }();
          static const int bound_e5m2_first = []() {
            const char *env = getenv("MX_BOUND_LOW_ORDER");
            if (!env) return 0;  // default keeps mx_e4m3 before mx_e5m2
            std::string v = env;
            std::transform(v.begin(), v.end(), v.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return (v == "e5m2_first" || v == "mx_e5m2_first" || v == "reverse") ? 1 : 0;
          }();

          struct BoundEvalResult {
            bool fits = false;
            bool overflow = false;
            bool uses_scale = false;
            long double u = 0.0L;
            long double fmin = 0.0L;
            long double fmax = 0.0L;
            long double gmin = 0.0L;
            long double tol = 0.0L;
            long double nrmE = 0.0L;
            long double nrmF = 0.0L;
            long double min_scale = 1.0L;
            long double max_scale = 1.0L;
            long double avg_scale = 1.0L;
            long long scale_groups = 0;
            long long underflows = 0;
          };

          auto evaluate_format_tile_bound = [&](const std::string &fmt_name) {
            BoundEvalResult out;
            const long double u = format_unit_roundoff(fmt_name);
            const long double fmin = format_fmin(fmt_name);
            const long double fmax = format_fmax(fmt_name);
            const long double gmin = format_gmin(fmt_name);
            out.u = u;
            out.fmin = fmin;
            out.fmax = fmax;
            out.gmin = gmin;
            if (!(u > 0.0L)) return out;

            const long double tol_tile =
                sourceEpsilon * normA / static_cast<long double>(array->nt);
            const long double nrmE = u * normTile;
            long double nrmF2 = 0.0L;
            bool overflow = false;
            out.tol = tol_tile;
            out.nrmE = nrmE;

            const bool uses_scale = format_uses_shared_scale(fmt_name);
            const bool use_fp32_bits = (fmt_name == "mx_fp32");
            const bool block_mode = (mx_mode() == MxMode::Block);
            const int subtile = mx_block_subtile();
            const size_t tile_m = static_cast<size_t>(tile.m);
            const size_t tile_n = static_cast<size_t>(tile.n);
            out.uses_scale = uses_scale;

            long double sum_scale = 0.0L;
            long double min_scale = (std::numeric_limits<long double>::max)();
            long double max_scale = 0.0L;
            long long scale_groups = 0;
            long long underflows = 0;

            auto note_scale = [&](long double pre_scale) {
              sum_scale += pre_scale;
              if (pre_scale < min_scale) min_scale = pre_scale;
              if (pre_scale > max_scale) max_scale = pre_scale;
              ++scale_groups;
            };

            // pre_scale is the multiplier used before quantization (x_qin = pre_scale * x).
            // This matches the implementation path where x_qin = x * inv_scale.
            auto consume_val = [&](long double ax, long double pre_scale) {
              const long double scaled = pre_scale * ax;
              if (fmax > 0.0L && scaled > fmax) {
                overflow = true;
                return;
              }
              if (fmin > 0.0L && gmin > 0.0L && scaled < fmin) {
                long double eta = (pre_scale > 0.0L) ? (gmin / pre_scale) : gmin;
                if (eta > ax) eta = ax;
                nrmF2 += eta * eta;
                ++underflows;
              }
            };

            if (!uses_scale) {
              for (size_t r = 0; r < tile_m && !overflow; ++r) {
                for (size_t c = 0; c < tile_n; ++c) {
                  const long double ax = std::fabs(static_cast<long double>(
                      mappedBlock(static_cast<Eigen::Index>(r),
                                  static_cast<Eigen::Index>(c))));
                  consume_val(ax, 1.0L);
                  if (overflow) break;
                }
              }
            } else if (block_mode && subtile > 0) {
              const size_t bs = static_cast<size_t>(subtile);
              for (size_t r0 = 0; r0 < tile_m && !overflow; r0 += bs) {
                const size_t r_max = (bs < (tile_m - r0)) ? bs : (tile_m - r0);
                for (size_t c0 = 0; c0 < tile_n && !overflow; c0 += bs) {
                  const size_t c_max = (bs < (tile_n - c0)) ? bs : (tile_n - c0);
                  float max_val = 0.0f;
                  for (size_t r = 0; r < r_max; ++r) {
                    for (size_t c = 0; c < c_max; ++c) {
                      const float v = static_cast<float>(std::fabs(
                          mappedBlock(static_cast<Eigen::Index>(r0 + r),
                                      static_cast<Eigen::Index>(c0 + c))));
                      max_val = std::fmax(max_val, v);
                    }
                  }
                  const int st_scale = use_fp32_bits
                                           ? computeScaleBitsFromMax(max_val,
                                                                     fp32_scale_bits)
                                           : static_cast<int>(
                                                 computeScaleFromMax(max_val));
                  const long double pre_scale = std::ldexp(1.0L, -st_scale);
                  note_scale(pre_scale);
                  for (size_t r = 0; r < r_max; ++r) {
                    for (size_t c = 0; c < c_max; ++c) {
                      const long double ax = std::fabs(static_cast<long double>(
                          mappedBlock(static_cast<Eigen::Index>(r0 + r),
                                      static_cast<Eigen::Index>(c0 + c))));
                      consume_val(ax, pre_scale);
                      if (overflow) break;
                    }
                    if (overflow) break;
                  }
                }
              }
            } else if (block_mode) {
              for (size_t r = 0; r < tile_m && !overflow; ++r) {
                float max_val = 0.0f;
                for (size_t c = 0; c < tile_n; ++c) {
                  const float v = static_cast<float>(std::fabs(
                      mappedBlock(static_cast<Eigen::Index>(r),
                                  static_cast<Eigen::Index>(c))));
                  max_val = std::fmax(max_val, v);
                }
                const int row_scale = use_fp32_bits
                                          ? computeScaleBitsFromMax(max_val,
                                                                    fp32_scale_bits)
                                          : static_cast<int>(
                                                computeScaleFromMax(max_val));
                const long double pre_scale = std::ldexp(1.0L, -row_scale);
                note_scale(pre_scale);
                for (size_t c = 0; c < tile_n; ++c) {
                  const long double ax = std::fabs(static_cast<long double>(
                      mappedBlock(static_cast<Eigen::Index>(r),
                                  static_cast<Eigen::Index>(c))));
                  consume_val(ax, pre_scale);
                  if (overflow) break;
                }
              }
            } else {
              float max_val = 0.0f;
              for (size_t r = 0; r < tile_m; ++r) {
                for (size_t c = 0; c < tile_n; ++c) {
                  const float v = static_cast<float>(std::fabs(
                      mappedBlock(static_cast<Eigen::Index>(r),
                                  static_cast<Eigen::Index>(c))));
                  max_val = std::fmax(max_val, v);
                }
              }
              const int tile_scale = use_fp32_bits
                                         ? computeScaleBitsFromMax(max_val,
                                                                   fp32_scale_bits)
                                         : static_cast<int>(
                                               computeScaleFromMax(max_val));
              const long double pre_scale = std::ldexp(1.0L, -tile_scale);
              note_scale(pre_scale);
              for (size_t r = 0; r < tile_m && !overflow; ++r) {
                for (size_t c = 0; c < tile_n; ++c) {
                  const long double ax = std::fabs(static_cast<long double>(
                      mappedBlock(static_cast<Eigen::Index>(r),
                                  static_cast<Eigen::Index>(c))));
                  consume_val(ax, pre_scale);
                  if (overflow) break;
                }
              }
            }

            out.overflow = overflow;
            const long double nrmF = std::sqrt(nrmF2);
            out.nrmF = nrmF;
            out.underflows = underflows;
            out.scale_groups = scale_groups;
            if (scale_groups > 0) {
              out.min_scale = min_scale;
              out.max_scale = max_scale;
              out.avg_scale = sum_scale / static_cast<long double>(scale_groups);
            }
            out.fits = (!overflow) && ((nrmE + nrmF) <= tol_tile);
            return out;
          };

          if (row == col ) {
            std::cout << "\n[DEBUG] Epsilon Ratio for Tile (" << row << ", " << col << "): " << epsilonRatio << std::endl;
            std::cout << "[DEBUG] Source Epsilon / FP32 Epsilon: " << sourceEpsilon / std::numeric_limits<float>::epsilon() << std::endl;
            std::cout << "[DEBUG] Low format cutoff (" << selected_low_format
                      << "): " << lowCutoff
                      << " (scale-aware=" << scale_aware_epsilon << ")"
                      << std::endl;
            tile.dtype = CUDA_R_64F;
            cudaMallocHost(&tile.data, sizeof(uint64_t) * tile.m * tile.n);
            g_pinned_tiles++;
            g_pinned_bytes += sizeof(uint64_t) * tile.m * tile.n;
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                         Eigen::RowMajor>>{
              static_cast<double *>(tile.data),
              static_cast<Eigen::Index>(tile.m),
              static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<double>();
            log_tile_target("fp64");
            std::cout << "[TILE_DTYPE] (" << row << ", " << col << ") "
                      << dtype_to_string(tile.dtype) << std::endl;
            continue;
          }

          if (force_fp32_all) {
            apply_bucket_target("fp32", "fp32");
          } else if (use_explicit_tile_bound_selection()) {
            const std::string fp32_choice =
                canonical_format_name(fp32_bucket, "fp32");
            const std::string fp16_choice =
                canonical_format_name(fp16_bucket, "fp16");
            const std::string low_choice =
                canonical_format_name(selected_low_format, "mx_e4m3");

            if (bound_debug) {
              std::cout << "[BOUND_SELECT] tile(" << row << "," << col << ")"
                        << " fp32_bucket=" << fp32_choice
                        << " fp16_bucket=" << fp16_choice
                        << " low_bucket=" << low_choice
                        << " ladder="
                        << (bound_ladder_mode == BoundLadderMode::Full
                          ? "full"
                          : (bound_ladder_mode == BoundLadderMode::IeeeOnly
                                 ? "ieee_only"
                                 : "legacy"))
                        << " low_order=" << (bound_e5m2_first ? "e5m2_first" : "e4m3_first")
                        << " ranking=" << (bound_cheap_first ? "cheap_first" : "legacy")
                        << " mode=" << (mx_mode() == MxMode::Block ? "block" : "tile")
                        << " subtile=" << mx_block_subtile() << std::endl;
            }

            std::vector<std::pair<std::string, BoundEvalResult>> evals;

            auto eval_one = [&](const std::string &fmt) {
              const auto ev = evaluate_format_tile_bound(fmt);
              evals.push_back({fmt, ev});
              return ev;
            };

            auto print_eval = [&](const std::string &name,
                                  const BoundEvalResult &ev) {
              if (!bound_debug) return;
              std::cout << "[BOUND_EVAL] tile(" << row << "," << col << ")"
                        << " fmt=" << name
                        << " fits=" << (ev.fits ? 1 : 0)
                        << " overflow=" << (ev.overflow ? 1 : 0)
                        << " uses_scale=" << (ev.uses_scale ? 1 : 0)
                        << " u=" << ev.u
                        << " fmin=" << ev.fmin
                        << " fmax=" << ev.fmax
                        << " gmin=" << ev.gmin
                        << " tol=" << ev.tol
                        << " nrmE=" << ev.nrmE
                        << " nrmF=" << ev.nrmF
                        << " nrmE+nrmF=" << (ev.nrmE + ev.nrmF)
                        << " underflows=" << ev.underflows
                        << " scale_groups=" << ev.scale_groups
                        << " s_min=" << ev.min_scale
                        << " s_avg=" << ev.avg_scale
                        << " s_max=" << ev.max_scale
                        << std::endl;
            };

            auto get_eval = [&](const std::string &fmt) {
              for (const auto &p : evals) {
                if (p.first == fmt) return p.second;
              }
              return BoundEvalResult{};
            };

            std::string decision_fmt;
            std::string decision_reason;
            int decision_step = -1;

            if (bound_ladder_mode == BoundLadderMode::Full) {
              // Full ascending-cost ladder requested by user:
              // FP4/FP6/FP8 -> FP16/MXFP16 -> FP32/MXFP32 -> FP64 fallback.
              const char *ladder_raw_default[] = {
                  "e2m1", "mx_e4m3", "mx_fp16", "fp32"};
              const char *ladder_raw_e5m2_first[] = {
                  "e2m1", "mx_e4m3", "mx_fp16", "fp32"};
              const char **ladder_raw =
                bound_e5m2_first ? ladder_raw_e5m2_first : ladder_raw_default;
              const int ladder_len = 4;

              int step = 0;
              bool placed = false;
              for (int idx = 0; idx < ladder_len; ++idx) {
              const char *cand_raw = ladder_raw[idx];
                std::string cand = canonical_format_name(cand_raw, cand_raw);
                const auto ev = eval_one(cand);
                print_eval(cand, ev);
                if (ev.fits) {
                  apply_bucket_target(cand, cand.c_str());
                  decision_fmt = cand;
                  decision_reason = "ladder_accept";
                  decision_step = step;
                  placed = true;
                  break;
                }
                ++step;
              }
              if (!placed) {
                log_tile_target("fp64");
                tile.dtype = CUDA_R_64F;
                cudaMallocHost(&tile.data, sizeof(uint64_t) * tile.m * tile.n);
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                        Eigen::RowMajor>>{
                    static_cast<double *>(tile.data),
                    static_cast<Eigen::Index>(tile.m),
                    static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<double>();
                decision_fmt = "fp64";
                decision_reason = "ladder_reject_all";
              }
            } else if (bound_ladder_mode == BoundLadderMode::IeeeOnly) {
              // IEEE-only ladder (no shared-scale MX formats):
              // FP8_E4M3 -> FP16 -> FP32 -> FP64 fallback.
              const char *ladder_raw[] = {"fp8_e4m3", "fp16", "fp32"};
              const int ladder_len = 3;

              int step = 0;
              bool placed = false;
              for (int idx = 0; idx < ladder_len; ++idx) {
                const char *cand_raw = ladder_raw[idx];
                std::string cand = canonical_format_name(cand_raw, cand_raw);
                const auto ev = eval_one(cand);
                print_eval(cand, ev);
                if (ev.fits) {
                  apply_bucket_target(cand, cand.c_str());
                  decision_fmt = cand;
                  decision_reason = "ieee_ladder_accept";
                  decision_step = step;
                  placed = true;
                  break;
                }
                ++step;
              }
              if (!placed) {
                log_tile_target("fp64");
                tile.dtype = CUDA_R_64F;
                cudaMallocHost(&tile.data, sizeof(uint64_t) * tile.m * tile.n);
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                        Eigen::RowMajor>>{
                    static_cast<double *>(tile.data),
                    static_cast<Eigen::Index>(tile.m),
                    static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<double>();
                decision_fmt = "fp64";
                decision_reason = "ieee_ladder_reject_all";
              }
            } else {
              const auto fp32_eval = eval_one(fp32_choice);
              const auto fp16_eval = eval_one(fp16_choice);
              const auto low_eval = allow_low
                                        ? eval_one(low_choice)
                                        : BoundEvalResult{};
              print_eval(fp32_choice, fp32_eval);
              print_eval(fp16_choice, fp16_eval);
              if (allow_low) {
                print_eval(low_choice, low_eval);
              } else if (bound_debug) {
                std::cout << "[BOUND_EVAL] tile(" << row << "," << col
                          << ") fmt=" << low_choice
                          << " skipped=1 reason=allow_low_false" << std::endl;
              }

              if (bound_cheap_first) {
                if (allow_low && low_eval.fits) {
                  apply_bucket_target(low_choice, "mx_e4m3");
                  decision_fmt = low_choice;
                  decision_reason = "cheap_first_accept_low";
                } else if (fp16_eval.fits) {
                  apply_bucket_target(fp16_choice, "fp16");
                  decision_fmt = fp16_choice;
                  decision_reason = allow_low ? "cheap_first_accept_fp16_reject_low"
                                              : "cheap_first_accept_fp16_low_not_allowed";
                } else if (fp32_eval.fits) {
                  apply_bucket_target(fp32_choice, "fp32");
                  decision_fmt = fp32_choice;
                  decision_reason = "cheap_first_accept_fp32";
                } else {
                  log_tile_target("fp64");
                  tile.dtype = CUDA_R_64F;
                  cudaMallocHost(&tile.data, sizeof(uint64_t) * tile.m * tile.n);
                  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                          Eigen::RowMajor>>{
                      static_cast<double *>(tile.data),
                      static_cast<Eigen::Index>(tile.m),
                      static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<double>();
                  decision_fmt = "fp64";
                  decision_reason = "cheap_first_reject_all";
                }
              } else if (!fp32_eval.fits) {
                log_tile_target("fp64");
                tile.dtype = CUDA_R_64F;
                cudaMallocHost(&tile.data, sizeof(uint64_t) * tile.m * tile.n);
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                        Eigen::RowMajor>>{
                    static_cast<double *>(tile.data),
                    static_cast<Eigen::Index>(tile.m),
                    static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<double>();
                decision_fmt = "fp64";
                decision_reason = "reject_fp32_bound";
              } else if (!fp16_eval.fits) {
                apply_bucket_target(fp32_choice, "fp32");
                decision_fmt = fp32_choice;
                decision_reason = "accept_fp32_reject_fp16";
              } else if (!allow_low || !low_eval.fits) {
                apply_bucket_target(fp16_choice, "fp16");
                decision_fmt = fp16_choice;
                decision_reason = allow_low ? "accept_fp16_reject_low"
                                            : "accept_fp16_low_not_allowed";
              } else {
                apply_bucket_target(low_choice, "mx_e4m3");
                decision_fmt = low_choice;
                decision_reason = "accept_low";
              }
            }

            if (bound_debug) {
              const auto fp32_eval = get_eval("fp32");
              const auto fp16_eval = get_eval("fp16");
              const auto low_eval = get_eval("mx_e4m3");
              std::cout << "[BOUND_DECISION] tile(" << row << "," << col << ")"
                        << " selected=" << decision_fmt
                        << " reason=" << decision_reason
                        << " step=" << decision_step
                        << " fp32_fit=" << (fp32_eval.fits ? 1 : 0)
                        << " fp16_fit=" << (fp16_eval.fits ? 1 : 0)
                        << " low_fit=" << (low_eval.fits ? 1 : 0)
                        << " fp32_overflow=" << (fp32_eval.overflow ? 1 : 0)
                        << " fp16_overflow=" << (fp16_eval.overflow ? 1 : 0)
                        << " low_overflow=" << (low_eval.overflow ? 1 : 0)
                        << std::endl;
            }
          } else if (epsilonRatio > sourceEpsilon / std::numeric_limits<float>::epsilon()) {
            log_tile_target("fp64");
            tile.dtype = CUDA_R_64F;
            cudaMallocHost(&tile.data, sizeof(uint64_t) * tile.m * tile.n);
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                    Eigen::RowMajor>>{
                static_cast<double *>(tile.data),
                static_cast<Eigen::Index>(tile.m),
                static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<double>();
          } else if (epsilonRatio > sourceEpsilon / std::numeric_limits<Eigen::half>::epsilon()) {
            std::cout << "---------------- Tile (" << row << ", " << col
                        << ") selected for CUDA_R_32F -------------------"
                        << std::endl;
            apply_bucket_target(fp32_bucket, "fp32");
          } else if (!allow_low ||
                     epsilonRatio > lowCutoff) {
            apply_bucket_target(fp16_bucket, "fp16");
          } else {
            if (target_low == "fp16") {
              log_tile_target("fp16");
              tile.dtype = CUDA_R_16F;
              cudaMallocHost(&tile.data, sizeof(uint16_t) * tile.m * tile.n);
              Eigen::Map<Eigen::Matrix<Eigen::half, Eigen::Dynamic, Eigen::Dynamic,
                                      Eigen::RowMajor>>{
                  static_cast<Eigen::half *>(tile.data),
                  static_cast<Eigen::Index>(tile.m),
                  static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<Eigen::half>();
            } else if (target_low == "bf16") {
              log_tile_target("bf16");
              tile.dtype = CUDA_R_16BF;
              cudaMallocHost(&tile.data, sizeof(uint16_t) * tile.m * tile.n);
              Eigen::Map<Eigen::Matrix<Eigen::bfloat16, Eigen::Dynamic, Eigen::Dynamic,
                                      Eigen::RowMajor>>{
                  static_cast<Eigen::bfloat16 *>(tile.data),
                  static_cast<Eigen::Index>(tile.m),
                  static_cast<Eigen::Index>(tile.n)} = mappedBlock.cast<Eigen::bfloat16>();
            } else if (target_low == "mx_fp16") {
              log_tile_target("mx_fp16");
              apply_mx_quant("MX_FP16", 5, 10, 65504.0f);
            } else if (target_low == "mx_fp32") {
              log_tile_target("mx_fp32");
              apply_mx_quant_fp64("MX_FP32", 8, 23,
                                  (std::numeric_limits<float>::max)(),
                                  fp32_scale_bits);
            } else if (target_low == "mx_e5m2") {
              log_tile_target("mx_e5m2");
              apply_mx_quant("MX_E5M2", 5, 2, 57344.0f);
            } else if (target_low == "fp8_e4m3") {
              log_tile_target("fp8_e4m3");
              apply_plain_fp_quant("FP8_E4M3", 4, 3, 448.0f);
            } else if (target_low == "fp8_e5m2") {
              log_tile_target("fp8_e5m2");
              apply_plain_fp_quant("FP8_E5M2", 5, 2, 57344.0f);
            } else if (target_low == "e3m2") {
              log_tile_target("e3m2");
              apply_mx_quant("MX_E3M2", 3, 2, 28.0f);
            } else if (target_low == "e2m3") {
              log_tile_target("e2m3");
              apply_mx_quant("MX_E2M3", 2, 3, 7.5f);
            } else if (target_low == "e2m1") {
              log_tile_target("e2m1");
              apply_mx_quant("MX_E2M1", 2, 1, 6.0f);
            } else {
              log_tile_target("mx_e4m3");
              apply_mx_quant("MX_E4M3", 4, 3, 448.0f);
            }
          }

          std::cout << "[TILE_DTYPE] (" << row << ", " << col << ") "
                    << dtype_to_string(tile.dtype) << std::endl;
        }
      }
    } else {
      std::cout << "Only Lower matrix is supported so far" << std::endl;
    }
}

void uncompressMixedPrecisionTiledArray(const MixedPrecisionTiledArray *array,
                                        double *data, const size_t ld) {
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
             Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, 1>>
      mappedA{data,
              static_cast<Eigen::Index>(array->m),
              static_cast<Eigen::Index>(array->n),
              {static_cast<Eigen::Index>(ld), 1}};
  if (array->uplo == PlasmaLower) {
    for (int row = 0; row < array->nt; ++row) {
      for (int col = 0; col <= row; ++col) {
        const auto &tile = array->tiles[col + row * array->nt];
        auto mappedBlock =
            mappedA.block(row * array->nb, col * array->nb, tile.m, tile.n);
        // If tile data is on device/managed, stage to host for the map.
        bool deviceTile = is_device_or_managed(tile.data);
        switch (tile.dtype) {
          case CUDA_R_64F: {
            const double *src = static_cast<const double *>(tile.data);
            std::vector<double> hostBuf;
            const double *hostPtr = src;
            if (deviceTile) {
              hostBuf.resize(tile.m * tile.n);
              cudaMemcpy(hostBuf.data(), src, tile.m * tile.n * sizeof(double),
                         cudaMemcpyDeviceToHost);
              hostPtr = hostBuf.data();
            }
            mappedBlock = Eigen::Map<const Eigen::Matrix<
                double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>{
                hostPtr, static_cast<Eigen::Index>(tile.m),
                static_cast<Eigen::Index>(tile.n)};
            break;
          }
          case CUDA_R_32F: {
            const float *src = static_cast<const float *>(tile.data);
            std::vector<float> hostBuf;
            const float *hostPtr = src;
            if (deviceTile) {
              hostBuf.resize(tile.m * tile.n);
              cudaMemcpy(hostBuf.data(), src, tile.m * tile.n * sizeof(float),
                         cudaMemcpyDeviceToHost);
              hostPtr = hostBuf.data();
            }
            mappedBlock =
                Eigen::Map<const Eigen::Matrix<
                    float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>{
                    hostPtr, static_cast<Eigen::Index>(tile.m),
                    static_cast<Eigen::Index>(tile.n)}
                    .cast<double>();
            break;
          }
          case CUDA_R_16F: {
            const Eigen::half *src = static_cast<const Eigen::half *>(tile.data);
            std::vector<Eigen::half> hostBuf;
            const Eigen::half *hostPtr = src;
            if (deviceTile) {
              hostBuf.resize(tile.m * tile.n);
              cudaMemcpy(hostBuf.data(), src, tile.m * tile.n * sizeof(Eigen::half),
                         cudaMemcpyDeviceToHost);
              hostPtr = hostBuf.data();
            }
            mappedBlock =
                Eigen::Map<
                    const Eigen::Matrix<Eigen::half, Eigen::Dynamic,
                                        Eigen::Dynamic, Eigen::RowMajor>>{
                    hostPtr, static_cast<Eigen::Index>(tile.m),
                    static_cast<Eigen::Index>(tile.n)}
                    .cast<double>();
            break;
          }
          case CUDA_R_16BF: {
            const Eigen::bfloat16 *src =
                static_cast<const Eigen::bfloat16 *>(tile.data);
            std::vector<Eigen::bfloat16> hostBuf;
            const Eigen::bfloat16 *hostPtr = src;
            if (deviceTile) {
              hostBuf.resize(tile.m * tile.n);
              cudaMemcpy(hostBuf.data(), src,
                         tile.m * tile.n * sizeof(Eigen::bfloat16),
                         cudaMemcpyDeviceToHost);
              hostPtr = hostBuf.data();
            }
            mappedBlock =
                Eigen::Map<
                    const Eigen::Matrix<Eigen::bfloat16, Eigen::Dynamic,
                                        Eigen::Dynamic, Eigen::RowMajor>>{
                    hostPtr, static_cast<Eigen::Index>(tile.m),
                    static_cast<Eigen::Index>(tile.n)}
                    .cast<double>();
            break;
          }
          case CUDA_R_8F_E4M3: {
            const __nv_fp8_e4m3 *src =
                static_cast<const __nv_fp8_e4m3 *>(tile.data);
            std::vector<__nv_fp8_e4m3> hostBuf;
            const __nv_fp8_e4m3 *hostPtr = src;
            if (deviceTile) {
              hostBuf.resize(tile.m * tile.n);
              cudaMemcpy(hostBuf.data(), src,
                         tile.m * tile.n * sizeof(__nv_fp8_e4m3),
                         cudaMemcpyDeviceToHost);
              hostPtr = hostBuf.data();
            }
            mappedBlock =
                Eigen::Map<
                    const Eigen::Matrix<__nv_fp8_e4m3, Eigen::Dynamic,
                                        Eigen::Dynamic, Eigen::RowMajor>>{
                    hostPtr, static_cast<Eigen::Index>(tile.m),
                    static_cast<Eigen::Index>(tile.n)}
                    .cast<double>();
            break;
          }
          default:
            printf("Wrong datatype in %s:%d\n", __FILE__, __LINE__);
            break;
        }
      }
    }
  } else {
    printf("Only Lower matrix is supported so far\n");
  }
}

void printMixedPrecisionTiledArray(MixedPrecisionTiledArray *array) {
  if (array->uplo == PlasmaLower) {
    for (int row = 0; row < array->nt; ++row) {
      for (int col = 0; col <= row; ++col) {
        auto &tile = array->tiles[col + row * array->nt];
        switch (tile.dtype) {
          case CUDA_R_64F:
            printf("\033[31mFP64\033[0m    ");
            break;
          case CUDA_R_32F:
            printf("\033[33mFP32\033[0m    ");
            break;
          case CUDA_R_16F:
            printf("\033[34mFP16\033[0m    ");
            break;
          case CUDA_R_16BF:
            printf("\033[34mBF16\033[0m    ");
            break;
          case CUDA_R_8F_E4M3:
            printf("\033[32mFP8\033[0m     ");
            break;
          default:
            printf("Not supported data type\n");
        }
      }
      printf("\n");
    }
  } else {
    printf("Only Lower matrix is supported so far\n");
  }
}

void freeMixedPrecisionTiledArray(MixedPrecisionTiledArray *array) {
  if (array->uplo == PlasmaLower) {
    for (int row = 0; row < array->nt; ++row) {
      for (int col = 0; col <= row; ++col) {
        auto &tile = array->tiles[col + row * array->nt];
        if (!tile.data) continue;
        if (is_device_or_managed(tile.data)) {
          CHECK_CUDA(cudaFree(tile.data));
        } else {
          CHECK_CUDA(cudaFreeHost(tile.data));
        }
      }
    }
  } else {
    printf("Only Lower matrix is supported so far\n");
  }

  delete[] array->tiles;

  if (g_quant_ws.capacity) {
    cudaFree(g_quant_ws.d_in);
    cudaFree(g_quant_ws.d_out);
    cudaFree(g_quant_ws.d_max);
    g_quant_ws = {};
  }
}

void requantizeTileHost(const MixedPrecisionTile *tile) {
  if (tile->mx_ebits == 0) return;
  if (tile->dtype != CUDA_R_32F) return;

  float *data = static_cast<float *>(tile->data);
  const size_t tile_size = tile->m * tile->n;
  const int ebits = tile->mx_ebits;
  const int mbits = tile->mx_mbits;
  const float max_norm = tile->mx_max_norm;
  const int scale_bits = tile->mx_scale_bits;
  const bool mx_fp16 = (ebits == 5 && mbits == 10);

  // Plain FP8: no MX block scaling, just clamp/round each element
  if (scale_bits < 0) {
    for (size_t i = 0; i < tile_size; ++i)
      data[i] = quantizeFp(data[i], ebits, mbits, max_norm);
    return;
  }

  if (tile->mx_mode_block == 2) {
    // 1D row-vector mode: each row split into groups of mx_block_subtile elements
    const size_t vec_sz = (tile->mx_block_subtile > 0)
                            ? static_cast<size_t>(tile->mx_block_subtile)
                            : 32;
    for (size_t r = 0; r < tile->m; ++r) {
      float *rp = data + r * tile->n;
      for (size_t c0 = 0; c0 < tile->n; c0 += vec_sz) {
        const size_t c_end = (c0 + vec_sz < tile->n) ? c0 + vec_sz : tile->n;
        float max_val = 0.0f;
        for (size_t c = c0; c < c_end; ++c)
          max_val = std::fmax(max_val, std::fabs(rp[c]));
        const int sc = computeScaleBitsFromMax(max_val, scale_bits);
        const float scale = pow2(sc);
        const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
        for (size_t c = c0; c < c_end; ++c) {
          const float x = rp[c] * inv_scale;
          if (mx_fp16)
            rp[c] = static_cast<float>(static_cast<Eigen::half>(x)) * scale;
          else
            rp[c] = quantizeFp(x, ebits, mbits, max_norm) * scale;
        }
      }
    }
  } else if (tile->mx_mode_block == 1 && tile->mx_block_subtile > 0) {
    const size_t tile_m = tile->m;
    const size_t tile_n = tile->n;
    const size_t st = static_cast<size_t>(tile->mx_block_subtile);
    for (size_t r0 = 0; r0 < tile_m; r0 += st) {
      const size_t r_end = (r0 + st < tile_m) ? r0 + st : tile_m;
      for (size_t c0 = 0; c0 < tile_n; c0 += st) {
        const size_t c_end = (c0 + st < tile_n) ? c0 + st : tile_n;
        float max_val = 0.0f;
        for (size_t r = r0; r < r_end; ++r)
          for (size_t c = c0; c < c_end; ++c)
            max_val = std::fmax(max_val, std::fabs(data[r * tile_n + c]));
        const int sc = computeScaleBitsFromMax(max_val, scale_bits);
        const float scale = pow2(sc);
        const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
        for (size_t r = r0; r < r_end; ++r)
          for (size_t c = c0; c < c_end; ++c) {
            float &v = data[r * tile_n + c];
            const float x = v * inv_scale;
            if (mx_fp16)
              v = static_cast<float>(static_cast<Eigen::half>(x)) * scale;
            else
              v = quantizeFp(x, ebits, mbits, max_norm) * scale;
          }
      }
    }
  } else if (tile->mx_mode_block) {
    for (size_t r = 0; r < tile->m; ++r) {
      float *rp = data + r * tile->n;
      const int sc = computeScaleBits(rp, tile->n, scale_bits);
      const float scale = pow2(sc);
      const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
      for (size_t c = 0; c < tile->n; ++c) {
        const float x = rp[c] * inv_scale;
        if (mx_fp16)
          rp[c] = static_cast<float>(static_cast<Eigen::half>(x)) * scale;
        else
          rp[c] = quantizeFp(x, ebits, mbits, max_norm) * scale;
      }
    }
  } else {
    const int sc = computeScaleBits(data, tile_size, scale_bits);
    const float scale = pow2(sc);
    const float inv_scale = (scale == 0.0f) ? 1.0f : 1.0f / scale;
    for (size_t i = 0; i < tile_size; ++i) {
      const float x = data[i] * inv_scale;
      if (mx_fp16)
        data[i] = static_cast<float>(static_cast<Eigen::half>(x)) * scale;
      else
        data[i] = quantizeFp(x, ebits, mbits, max_norm) * scale;
    }
  }
}

size_t getSizeofTileElement(int dtype) {
  switch (dtype) {
    case CUDA_R_64F:
      return sizeof(uint64_t);
    case CUDA_R_32F:
      return sizeof(uint32_t);
    case CUDA_R_16F:
      return sizeof(uint16_t);
    case CUDA_R_16BF:
      return sizeof(uint16_t);
    case CUDA_R_8F_E4M3:
      return sizeof(uint8_t);
    default:
      printf("Not supported data type\n");
      return 0;
  }
}

int getComputeType(int dtype) {
  switch (dtype) {
    case CUDA_R_64F:
      return CUBLAS_COMPUTE_64F;
    case CUDA_R_32F:
      return CUBLAS_COMPUTE_32F;
    case CUDA_R_16F:
      return CUBLAS_COMPUTE_32F;
    case CUDA_R_16BF:
      return CUBLAS_COMPUTE_32F;
    case CUDA_R_8F_E4M3:
      return CUBLAS_COMPUTE_32F;
    default:
      printf("Not supported data type\n");
      return 0;
  }
}

// int PLASMA_dpotrf_gpu_async_copy(PLASMA_enum uplo, int N,
//                                  double *A, int LDA)
//{
//   int NB;
//   int status;
//   plasma_context_t *plasma;
//   PLASMA_sequence *sequence = NULL;
//   PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
//   PLASMA_desc descA;
//
//   plasma = plasma_context_self();
//   if (plasma == NULL) {
//     plasma_fatal_error("PLASMA_dpotrf", "PLASMA not initialized");
//     return PLASMA_ERR_NOT_INITIALIZED;
//   }
//   /* Check input arguments */
//   if (uplo != PlasmaUpper && uplo != PlasmaLower) {
//     plasma_error("PLASMA_dpotrf", "illegal value of uplo");
//     return -1;
//   }
//   if (N < 0) {
//     plasma_error("PLASMA_dpotrf", "illegal value of N");
//     return -2;
//   }
//   if (LDA < max(1, N)) {
//     plasma_error("PLASMA_dpotrf", "illegal value of LDA");
//     return -4;
//   }
//   /* Quick return */
//   if (max(N, 0) == 0)
//     return PLASMA_SUCCESS;
//
//   /* Tune NB depending on M, N & NRHS; Set NBNB */
//   status = plasma_tune(PLASMA_FUNC_DPOSV, N, N, 0);
//   if (status != PLASMA_SUCCESS) {
//     plasma_error("PLASMA_dpotrf", "plasma_tune() failed");
//     return status;
//   }
//
//   /* Set NT */
//   NB   = PLASMA_NB;
//
//   plasma_sequence_create(plasma, &sequence);
//
//   if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
//     plasma_dooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence,
//     &request,
//                         plasma_desc_mat_free(&(descA)) );
//   } else {
//     plasma_diplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
//                        sequence, &request);
//   }
//
//   /* Call the tile interface */
//   PLASMA_dpotrf_gpu_async_copy_Tile_Async(uplo, &descA, sequence, &request);
//
//   if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
//     plasma_dooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
//     plasma_dynamic_sync();
//     plasma_desc_mat_free(&descA);
//   } else {
//     plasma_diptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
//     plasma_dynamic_sync();
//   }
//
//   status = sequence->status;
//   plasma_sequence_destroy(plasma, sequence);
//
//   return status;
// }

int PLASMA_dpotrf_gpu_reuse_data_mixed_precision(PLASMA_enum uplo, int N,
                                                 double *A, int LDA) {
  int NB;
  int status;
  plasma_context_t *plasma;
  PLASMA_sequence *sequence = NULL;
  PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
  printf("Starting PLASMA_dpotrf_gpu_reuse_data_mixed_precision...\n");
  plasma = plasma_context_self();
  if (plasma == NULL) {
    plasma_fatal_error("PLASMA_dpotrf", "PLASMA not initialized");
    return PLASMA_ERR_NOT_INITIALIZED;
  }
  /* Check input arguments */
  if (uplo != PlasmaUpper && uplo != PlasmaLower) {
    plasma_error("PLASMA_dpotrf", "illegal value of uplo");
    return -1;
  }
  if (N < 0) {
    plasma_error("PLASMA_dpotrf", "illegal value of N");
    return -2;
  }
  if (LDA < max(1, N)) {
    plasma_error("PLASMA_dpotrf", "illegal value of LDA");
    return -4;
  }
  /* Quick return */
  if (max(N, 0) == 0) return PLASMA_SUCCESS;

  /* Tune NB depending on M, N & NRHS; Set NBNB */
  status = plasma_tune(PLASMA_FUNC_DPOSV, N, N, 0);
  if (status != PLASMA_SUCCESS) {
    plasma_error("PLASMA_dpotrf", "plasma_tune() failed");
    return status;
  }

  /* Set NT */
  NB = (plasma->autotuning_enabled ? PLASMA_NB : plasma->nb);
  std::cout << "[DEBUG] using NB=" << NB << " (autotune="
            << plasma->autotuning_enabled << ")" << std::endl;

  plasma_sequence_create(plasma, &sequence);

  MixedPrecisionTiledArray array;
  printf("Making mixed-precision tiled array...\n");
  auto t0 = std::chrono::high_resolution_clock::now();
  makeMixedPrecisionTiledArray(&array, uplo, A, N, N, LDA, NB);
  auto t1 = std::chrono::high_resolution_clock::now();

  //    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
  //        plasma_dooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence,
  //        &request,
  //                            plasma_desc_mat_free(&(descA)) );
  //    } else {
  //        plasma_diplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
  //                           sequence, &request);
  //    }

  auto start = std::chrono::high_resolution_clock::now();

  /* Call the tile interface */
  PLASMA_dpotrf_gpu_reuse_data_mixed_precision_Tile_Async(uplo, &array,
                                                          sequence, &request);

  auto duration = std::chrono::high_resolution_clock::now() - start;
  auto t2 = std::chrono::high_resolution_clock::now();
  const auto micros =
      std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
  const long double flops =
      (long double)(FLOPS_DPOTRF_ULL(static_cast<unsigned long long>(N)));
  const long double seconds = static_cast<long double>(micros) * 1e-6L;
  const long double tflops = flops / seconds / 1e12L;
  printf("FLOPS: %Lf\n", tflops);

  const auto tile_us =
      std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
  const auto factor_us =
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - start)
          .count();
  const auto uncompress_us =
      std::chrono::duration_cast<std::chrono::microseconds>(
          std::chrono::high_resolution_clock::now() - t2)
          .count();
  std::cout << "[TIMES] tile+quant(us): " << tile_us
            << ", factor(us): " << factor_us
            << ", uncompress(us): " << uncompress_us << std::endl;

  const long expected_tiles = (long)array.nt * (array.nt + 1) / 2;
  std::cout << "[STATS] nt: " << array.nt
            << ", expected lower tiles: " << expected_tiles
            << ", pinned tiles: " << g_pinned_tiles
            << ", pinned bytes: " << g_pinned_bytes
            << ", quant workspace bytes: " << g_quant_ws_bytes << std::endl;

  uncompressMixedPrecisionTiledArray(&array, A, LDA);
  freeMixedPrecisionTiledArray(&array);

  //    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
  //        plasma_dooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
  //        plasma_dynamic_sync();
  //        plasma_desc_mat_free(&descA);
  //    } else {
  //        plasma_diptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
  //        plasma_dynamic_sync();
  //    }

  status = sequence->status;
  plasma_sequence_destroy(plasma, sequence);

  return status;
}

int PLASMA_dpotrf_gpu_reuse_data_table_mixed_precision(PLASMA_enum uplo, int N,
                                                       double *A, int LDA) {
  int NB;
  int status;
  plasma_context_t *plasma;
  PLASMA_sequence *sequence = NULL;
  PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

  plasma = plasma_context_self();
  if (plasma == NULL) {
    plasma_fatal_error("PLASMA_dpotrf", "PLASMA not initialized");
    return PLASMA_ERR_NOT_INITIALIZED;
  }
  /* Check input arguments */
  if (uplo != PlasmaUpper && uplo != PlasmaLower) {
    plasma_error("PLASMA_dpotrf", "illegal value of uplo");
    return -1;
  }
  if (N < 0) {
    plasma_error("PLASMA_dpotrf", "illegal value of N");
    return -2;
  }
  if (LDA < max(1, N)) {
    plasma_error("PLASMA_dpotrf", "illegal value of LDA");
    return -4;
  }
  /* Quick return */
  if (max(N, 0) == 0) return PLASMA_SUCCESS;

  /* Tune NB depending on M, N & NRHS; Set NBNB */
  status = plasma_tune(PLASMA_FUNC_DPOSV, N, N, 0);
  if (status != PLASMA_SUCCESS) {
    plasma_error("PLASMA_dpotrf", "plasma_tune() failed");
    return status;
  }

  /* Set NT */
  NB = (plasma->autotuning_enabled ? PLASMA_NB : plasma->nb);
  std::cout << "[DEBUG] using NB=" << NB << " (autotune="
            << plasma->autotuning_enabled << ")" << std::endl;

  plasma_sequence_create(plasma, &sequence);

  MixedPrecisionTiledArray array;
  auto t0 = std::chrono::high_resolution_clock::now();
  makeMixedPrecisionTiledArray(&array, uplo, A, N, N, LDA, NB);
  auto t1 = std::chrono::high_resolution_clock::now();

  //    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
  //        plasma_dooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence,
  //        &request,
  //                            plasma_desc_mat_free(&(descA)) );
  //    } else {
  //        plasma_diplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
  //                           sequence, &request);
  //    }

  auto start = std::chrono::high_resolution_clock::now();

  /* Call the tile interface */
  PLASMA_dpotrf_gpu_reuse_data_table_mixed_precision_Tile_Async(
      uplo, &array, sequence, &request);

  auto duration = std::chrono::high_resolution_clock::now() - start;
  const auto micros =
      std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
  const long double flops =
      (long double)(FLOPS_DPOTRF_ULL(static_cast<unsigned long long>(N)));
  const long double seconds = static_cast<long double>(micros) * 1e-6L;
  const long double tflops = flops / seconds / 1e12L;
  printf("FLOPS: %Lf\n", tflops);

  auto t2 = std::chrono::high_resolution_clock::now();
  const auto tile_us =
      std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
  const auto factor_us =
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - start)
          .count();
  const auto uncompress_us =
      std::chrono::duration_cast<std::chrono::microseconds>(
          std::chrono::high_resolution_clock::now() - t2)
          .count();
  std::cout << "[TIMES] tile+quant(us): " << tile_us
            << ", factor(us): " << factor_us
            << ", uncompress(us): " << uncompress_us << std::endl;

  uncompressMixedPrecisionTiledArray(&array, A, LDA);
  freeMixedPrecisionTiledArray(&array);

  //    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
  //        plasma_dooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
  //        plasma_dynamic_sync();
  //        plasma_desc_mat_free(&descA);
  //    } else {
  //        plasma_diptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
  //        plasma_dynamic_sync();
  //    }

  status = sequence->status;
  plasma_sequence_destroy(plasma, sequence);

  return status;
}

int PLASMA_dpotrf_gpu_reuse_data_table_all_managed_mixed_precision(
    PLASMA_enum uplo, int N, double *A, int LDA) {
  int NB;
  int status;
  plasma_context_t *plasma;
  PLASMA_sequence *sequence = NULL;
  PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

  plasma = plasma_context_self();
  if (plasma == NULL) {
    plasma_fatal_error("PLASMA_dpotrf", "PLASMA not initialized");
    return PLASMA_ERR_NOT_INITIALIZED;
  }
  /* Check input arguments */
  if (uplo != PlasmaUpper && uplo != PlasmaLower) {
    plasma_error("PLASMA_dpotrf", "illegal value of uplo");
    return -1;
  }
  if (N < 0) {
    plasma_error("PLASMA_dpotrf", "illegal value of N");
    return -2;
  }
  if (LDA < max(1, N)) {
    plasma_error("PLASMA_dpotrf", "illegal value of LDA");
    return -4;
  }
  /* Quick return */
  if (max(N, 0) == 0) return PLASMA_SUCCESS;

  /* Tune NB depending on M, N & NRHS; Set NBNB */
  status = plasma_tune(PLASMA_FUNC_DPOSV, N, N, 0);
  if (status != PLASMA_SUCCESS) {
    plasma_error("PLASMA_dpotrf", "plasma_tune() failed");
    return status;
  }

  /* Set NT */
  NB = (plasma->autotuning_enabled ? PLASMA_NB : plasma->nb);
  std::cout << "[DEBUG] using NB=" << NB << " (autotune="
            << plasma->autotuning_enabled << ")" << std::endl;

  plasma_sequence_create(plasma, &sequence);

  MixedPrecisionTiledArray array;
  auto t0 = std::chrono::high_resolution_clock::now();
  makeMixedPrecisionTiledArray(&array, uplo, A, N, N, LDA, NB);
  auto t1 = std::chrono::high_resolution_clock::now();

  //    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
  //        plasma_dooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence,
  //        &request,
  //                            plasma_desc_mat_free(&(descA)) );
  //    } else {
  //        plasma_diplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
  //                           sequence, &request);
  //    }

  auto start = std::chrono::high_resolution_clock::now();

  /* Call the tile interface */
  PLASMA_dpotrf_gpu_reuse_data_table_all_managed_mixed_precision_Tile_Async(
      uplo, &array, sequence, &request);

  auto duration = std::chrono::high_resolution_clock::now() - start;
  const auto micros =
      std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
  const long double flops =
      (long double)(FLOPS_DPOTRF_ULL(static_cast<unsigned long long>(N)));
  const long double seconds = static_cast<long double>(micros) * 1e-6L;
  const long double tflops = flops / seconds / 1e12L;
  printf("FLOPS: %Lf\n", tflops);

  auto t2 = std::chrono::high_resolution_clock::now();
  const auto tile_us =
      std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
  const auto factor_us =
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - start)
          .count();
  const auto uncompress_us =
      std::chrono::duration_cast<std::chrono::microseconds>(
          std::chrono::high_resolution_clock::now() - t2)
          .count();
  std::cout << "[TIMES] tile+quant(us): " << tile_us
            << ", factor(us): " << factor_us
            << ", uncompress(us): " << uncompress_us << std::endl;

  uncompressMixedPrecisionTiledArray(&array, A, LDA);
  freeMixedPrecisionTiledArray(&array);

  //    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
  //        plasma_dooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
  //        plasma_dynamic_sync();
  //        plasma_desc_mat_free(&descA);
  //    } else {
  //        plasma_diptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
  //        plasma_dynamic_sync();
  //    }

  status = sequence->status;
  plasma_sequence_destroy(plasma, sequence);

  return status;
}

/***************************************************************************/
/**
 *
 * @ingroup double_Tile
 *
 *  PLASMA_dpotrf_Tile - Computes the Cholesky factorization of a symmetric
 *positive definite or Hermitian positive definite matrix. Tile equivalent of
 *PLASMA_dpotrf(). Operates on matrices stored by tiles. All matrices are passed
 *through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] A
 *          On entry, the symmetric positive definite (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly
 *lower triangular part of A is not referenced. If UPLO = 'L', the leading
 *N-by-N lower triangular part of A contains the lower triangular part of the
 *matrix A, and the strictly upper triangular part of A is not referenced. On
 *exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = U**T*U or A = L*L**T.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, the leading minor of order i of A is not positive
 *definite, so the factorization could not be completed, and the solution has
 *not been computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dpotrf
 * @sa PLASMA_dpotrf_Tile_Async
 * @sa PLASMA_cpotrf_Tile
 * @sa PLASMA_dpotrf_Tile
 * @sa PLASMA_spotrf_Tile
 * @sa PLASMA_dpotrs_Tile
 *
 ******************************************************************************/
// int PLASMA_dpotrf_gpu_reuse_data_Tile(PLASMA_enum uplo, PLASMA_desc *A)
//{
//     plasma_context_t *plasma;
//     PLASMA_sequence *sequence = NULL;
//     PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
//     int status;
//
//     plasma = plasma_context_self();
//     if (plasma == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile", "PLASMA not initialized");
//         return PLASMA_ERR_NOT_INITIALIZED;
//     }
//     plasma_sequence_create(plasma, &sequence);
//     PLASMA_dpotrf_gpu_reuse_data_Tile_Async(uplo, A, sequence, &request);
//     plasma_dynamic_sync();
//     status = sequence->status;
//     plasma_sequence_destroy(plasma, sequence);
//     return status;
// }
//
// int PLASMA_dpotrf_gpu_reuse_data_table_Tile(PLASMA_enum uplo, PLASMA_desc *A)
//{
//     plasma_context_t *plasma;
//     PLASMA_sequence *sequence = NULL;
//     PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
//     int status;
//
//     plasma = plasma_context_self();
//     if (plasma == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile", "PLASMA not initialized");
//         return PLASMA_ERR_NOT_INITIALIZED;
//     }
//     plasma_sequence_create(plasma, &sequence);
//     PLASMA_dpotrf_gpu_reuse_data_table_Tile_Async(uplo, A, sequence,
//     &request); plasma_dynamic_sync(); status = sequence->status;
//     plasma_sequence_destroy(plasma, sequence);
//     return status;
// }
//
// int PLASMA_dpotrf_gpu_reuse_data_table_all_managed_Tile(PLASMA_enum uplo,
// PLASMA_desc *A)
//{
//     plasma_context_t *plasma;
//     PLASMA_sequence *sequence = NULL;
//     PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
//     int status;
//
//     plasma = plasma_context_self();
//     if (plasma == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile", "PLASMA not initialized");
//         return PLASMA_ERR_NOT_INITIALIZED;
//     }
//     plasma_sequence_create(plasma, &sequence);
//     PLASMA_dpotrf_gpu_reuse_data_table_all_managed_Tile_Async(uplo, A,
//     sequence, &request); plasma_dynamic_sync(); status = sequence->status;
//     plasma_sequence_destroy(plasma, sequence);
//     return status;
// }
//
// int PLASMA_dpotrf_gpu_reuse_data_table_no_free_Tile(PLASMA_enum uplo,
// PLASMA_desc *A)
//{
//     plasma_context_t *plasma;
//     PLASMA_sequence *sequence = NULL;
//     PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
//     int status;
//
//     plasma = plasma_context_self();
//     if (plasma == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile", "PLASMA not initialized");
//         return PLASMA_ERR_NOT_INITIALIZED;
//     }
//     plasma_sequence_create(plasma, &sequence);
//     PLASMA_dpotrf_gpu_reuse_data_table_no_free_Tile_Async(uplo, A, sequence,
//     &request); plasma_dynamic_sync(); status = sequence->status;
//     plasma_sequence_destroy(plasma, sequence);
//     return status;
// }
//
// int PLASMA_dpotrf_gpu_reuse_data_static_table_Tile(PLASMA_enum uplo,
// PLASMA_desc *A)
//{
//     plasma_context_t *plasma;
//     PLASMA_sequence *sequence = NULL;
//     PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
//     int status;
//
//     plasma = plasma_context_self();
//     if (plasma == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile", "PLASMA not initialized");
//         return PLASMA_ERR_NOT_INITIALIZED;
//     }
//     plasma_sequence_create(plasma, &sequence);
//     PLASMA_dpotrf_gpu_reuse_data_static_table_Tile_Async(uplo, A, sequence,
//     &request); plasma_dynamic_sync(); status = sequence->status;
//     plasma_sequence_destroy(plasma, sequence);
//     return status;
// }

/***************************************************************************/
/**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dpotrf_Tile_Async - Computes the Cholesky factorization of a
 *symmetric positive definite or Hermitian positive definite matrix.
 *Non-blocking equivalent of PLASMA_dpotrf_Tile(). May return before the
 *computation is finished. Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa PLASMA_dpotrf
 * @sa PLASMA_dpotrf_Tile
 * @sa PLASMA_cpotrf_Tile_Async
 * @sa PLASMA_dpotrf_Tile_Async
 * @sa PLASMA_spotrf_Tile_Async
 * @sa PLASMA_dpotrs_Tile_Async
 *
 ******************************************************************************/
// int PLASMA_dpotrf_gpu_async_copy_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A,
//                                             PLASMA_sequence *sequence,
//                                             PLASMA_request *request)
//{
//     PLASMA_desc descA;
//     plasma_context_t *plasma;
//
//     plasma = plasma_context_self();
//     if (plasma == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "PLASMA not
//         initialized"); return PLASMA_ERR_NOT_INITIALIZED;
//     }
//     if (sequence == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL sequence");
//         return PLASMA_ERR_UNALLOCATED;
//     }
//     if (request == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL request");
//         return PLASMA_ERR_UNALLOCATED;
//     }
//     /* Check sequence status */
//     if (sequence->status == PLASMA_SUCCESS)
//         request->status = PLASMA_SUCCESS;
//     else
//         return plasma_request_fail(sequence, request,
//         PLASMA_ERR_SEQUENCE_FLUSHED);
//
//     /* Check descriptors for correctness */
//     if (plasma_desc_check(A) != PLASMA_SUCCESS) {
//         plasma_error("PLASMA_dpotrf_Tile_Async", "invalid descriptor");
//         return plasma_request_fail(sequence, request,
//         PLASMA_ERR_ILLEGAL_VALUE);
//     } else {
//         descA = *A;
//     }
//     /* Check input arguments */
//     if (descA.nb != descA.mb) {
//         plasma_error("PLASMA_dpotrf_Tile_Async", "only square tiles
//         supported"); return plasma_request_fail(sequence, request,
//         PLASMA_ERR_ILLEGAL_VALUE);
//     }
//     if (uplo != PlasmaUpper && uplo != PlasmaLower) {
//         plasma_error("PLASMA_dpotrf_Tile_Async", "illegal value of uplo");
//         return plasma_request_fail(sequence, request, -1);
//     }
//     /* Quick return */
///*
//    if (max(N, 0) == 0)
//        return PLASMA_SUCCESS;
//*/
//    plasma_parallel_call_4(plasma_pdpotrf_gpu_async_copy,
//        PLASMA_enum, uplo,
//        PLASMA_desc, descA,
//        PLASMA_sequence*, sequence,
//        PLASMA_request*, request);
//
//    return PLASMA_SUCCESS;
//}

int PLASMA_dpotrf_gpu_reuse_data_mixed_precision_Tile_Async(
    PLASMA_enum uplo, MixedPrecisionTiledArray *A, PLASMA_sequence *sequence,
    PLASMA_request *request) {
  MixedPrecisionTiledArray descA;
  plasma_context_t *plasma;

  plasma = plasma_context_self();
  if (plasma == NULL) {
    plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "PLASMA not initialized");
    return PLASMA_ERR_NOT_INITIALIZED;
  }
  if (sequence == NULL) {
    plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL sequence");
    return PLASMA_ERR_UNALLOCATED;
  }
  if (request == NULL) {
    plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL request");
    return PLASMA_ERR_UNALLOCATED;
  }
  /* Check sequence status */
  if (sequence->status == PLASMA_SUCCESS)
    request->status = PLASMA_SUCCESS;
  else
    return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

  /* Check descriptors for correctness */
  if (PLASMA_SUCCESS /*plasma_desc_check(A)*/ != PLASMA_SUCCESS) {
    plasma_error("PLASMA_dpotrf_Tile_Async", "invalid descriptor");
    return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
  } else {
    descA = *A;
  }
  /* Check input arguments */
  if (descA.nb != descA.mb) {
    plasma_error("PLASMA_dpotrf_Tile_Async", "only square tiles supported");
    return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
  }
  if (uplo != PlasmaUpper && uplo != PlasmaLower) {
    plasma_error("PLASMA_dpotrf_Tile_Async", "illegal value of uplo");
    return plasma_request_fail(sequence, request, -1);
  }
  /* Quick return */
  /*
      if (max(N, 0) == 0)
          return PLASMA_SUCCESS;
  */
  plasma_parallel_call_4(plasma_pdpotrf_gpu_reuse_data_mixed_precision,
                         PLASMA_enum, uplo, MixedPrecisionTiledArray, descA,
                         PLASMA_sequence *, sequence, PLASMA_request *,
                         request);

  return PLASMA_SUCCESS;
}

int PLASMA_dpotrf_gpu_reuse_data_table_mixed_precision_Tile_Async(
    PLASMA_enum uplo, MixedPrecisionTiledArray *A, PLASMA_sequence *sequence,
    PLASMA_request *request) {
  MixedPrecisionTiledArray descA;
  plasma_context_t *plasma;

  plasma = plasma_context_self();
  if (plasma == NULL) {
    plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "PLASMA not initialized");
    return PLASMA_ERR_NOT_INITIALIZED;
  }
  if (sequence == NULL) {
    plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL sequence");
    return PLASMA_ERR_UNALLOCATED;
  }
  if (request == NULL) {
    plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL request");
    return PLASMA_ERR_UNALLOCATED;
  }
  /* Check sequence status */
  if (sequence->status == PLASMA_SUCCESS)
    request->status = PLASMA_SUCCESS;
  else
    return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

  /* Check descriptors for correctness */
  if (PLASMA_SUCCESS /*plasma_desc_check(A)*/ != PLASMA_SUCCESS) {
    plasma_error("PLASMA_dpotrf_Tile_Async", "invalid descriptor");
    return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
  } else {
    descA = *A;
  }
  /* Check input arguments */
  if (descA.nb != descA.mb) {
    plasma_error("PLASMA_dpotrf_Tile_Async", "only square tiles supported");
    return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
  }
  if (uplo != PlasmaUpper && uplo != PlasmaLower) {
    plasma_error("PLASMA_dpotrf_Tile_Async", "illegal value of uplo");
    return plasma_request_fail(sequence, request, -1);
  }
  /* Quick return */
  /*
      if (max(N, 0) == 0)
          return PLASMA_SUCCESS;
  */
  plasma_parallel_call_4(plasma_pdpotrf_gpu_reuse_data_table_mixed_precision,
                         PLASMA_enum, uplo, MixedPrecisionTiledArray, descA,
                         PLASMA_sequence *, sequence, PLASMA_request *,
                         request);

  return PLASMA_SUCCESS;
}

int PLASMA_dpotrf_gpu_reuse_data_table_all_managed_mixed_precision_Tile_Async(
    PLASMA_enum uplo, MixedPrecisionTiledArray *A, PLASMA_sequence *sequence,
    PLASMA_request *request) {
  MixedPrecisionTiledArray descA;
  plasma_context_t *plasma;

  plasma = plasma_context_self();
  if (plasma == NULL) {
    plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "PLASMA not initialized");
    return PLASMA_ERR_NOT_INITIALIZED;
  }
  if (sequence == NULL) {
    plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL sequence");
    return PLASMA_ERR_UNALLOCATED;
  }
  if (request == NULL) {
    plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL request");
    return PLASMA_ERR_UNALLOCATED;
  }
  /* Check sequence status */
  if (sequence->status == PLASMA_SUCCESS)
    request->status = PLASMA_SUCCESS;
  else
    return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

  /* Check descriptors for correctness */
  if (PLASMA_SUCCESS /*plasma_desc_check(A)*/ != PLASMA_SUCCESS) {
    plasma_error("PLASMA_dpotrf_Tile_Async", "invalid descriptor");
    return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
  } else {
    descA = *A;
  }
  /* Check input arguments */
  if (descA.nb != descA.mb) {
    plasma_error("PLASMA_dpotrf_Tile_Async", "only square tiles supported");
    return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
  }
  if (uplo != PlasmaUpper && uplo != PlasmaLower) {
    plasma_error("PLASMA_dpotrf_Tile_Async", "illegal value of uplo");
    return plasma_request_fail(sequence, request, -1);
  }
  /* Quick return */
  /*
      if (max(N, 0) == 0)
          return PLASMA_SUCCESS;
  */
  plasma_parallel_call_4(
      plasma_pdpotrf_gpu_reuse_data_table_all_managed_mixed_precision,
      PLASMA_enum, uplo, MixedPrecisionTiledArray, descA, PLASMA_sequence *,
      sequence, PLASMA_request *, request);

  return PLASMA_SUCCESS;
}

// int PLASMA_dpotrf_gpu_reuse_data_table_no_free_Tile_Async(PLASMA_enum uplo,
// PLASMA_desc *A,
//                                                           PLASMA_sequence
//                                                           *sequence,
//                                                           PLASMA_request
//                                                           *request)
//{
//     PLASMA_desc descA;
//     plasma_context_t *plasma;
//
//     plasma = plasma_context_self();
//     if (plasma == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "PLASMA not
//         initialized"); return PLASMA_ERR_NOT_INITIALIZED;
//     }
//     if (sequence == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL sequence");
//         return PLASMA_ERR_UNALLOCATED;
//     }
//     if (request == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL request");
//         return PLASMA_ERR_UNALLOCATED;
//     }
//     /* Check sequence status */
//     if (sequence->status == PLASMA_SUCCESS)
//         request->status = PLASMA_SUCCESS;
//     else
//         return plasma_request_fail(sequence, request,
//         PLASMA_ERR_SEQUENCE_FLUSHED);
//
//     /* Check descriptors for correctness */
//     if (plasma_desc_check(A) != PLASMA_SUCCESS) {
//         plasma_error("PLASMA_dpotrf_Tile_Async", "invalid descriptor");
//         return plasma_request_fail(sequence, request,
//         PLASMA_ERR_ILLEGAL_VALUE);
//     } else {
//         descA = *A;
//     }
//     /* Check input arguments */
//     if (descA.nb != descA.mb) {
//         plasma_error("PLASMA_dpotrf_Tile_Async", "only square tiles
//         supported"); return plasma_request_fail(sequence, request,
//         PLASMA_ERR_ILLEGAL_VALUE);
//     }
//     if (uplo != PlasmaUpper && uplo != PlasmaLower) {
//         plasma_error("PLASMA_dpotrf_Tile_Async", "illegal value of uplo");
//         return plasma_request_fail(sequence, request, -1);
//     }
//     /* Quick return */
//     /*
//         if (max(N, 0) == 0)
//             return PLASMA_SUCCESS;
//     */
//     plasma_parallel_call_4(plasma_pdpotrf_gpu_reuse_data_table_no_free,
//                            PLASMA_enum, uplo,
//                            PLASMA_desc, descA,
//                            PLASMA_sequence*, sequence,
//                            PLASMA_request*, request);
//
//     return PLASMA_SUCCESS;
// }
//
// int PLASMA_dpotrf_gpu_reuse_data_static_table_Tile_Async(PLASMA_enum uplo,
// PLASMA_desc *A,
//                                                               PLASMA_sequence
//                                                               *sequence,
//                                                               PLASMA_request
//                                                               *request)
//{
//     PLASMA_desc descA;
//     plasma_context_t *plasma;
//
//     plasma = plasma_context_self();
//     if (plasma == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "PLASMA not
//         initialized"); return PLASMA_ERR_NOT_INITIALIZED;
//     }
//     if (sequence == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL sequence");
//         return PLASMA_ERR_UNALLOCATED;
//     }
//     if (request == NULL) {
//         plasma_fatal_error("PLASMA_dpotrf_Tile_Async", "NULL request");
//         return PLASMA_ERR_UNALLOCATED;
//     }
//     /* Check sequence status */
//     if (sequence->status == PLASMA_SUCCESS)
//         request->status = PLASMA_SUCCESS;
//     else
//         return plasma_request_fail(sequence, request,
//         PLASMA_ERR_SEQUENCE_FLUSHED);
//
//     /* Check descriptors for correctness */
//     if (plasma_desc_check(A) != PLASMA_SUCCESS) {
//         plasma_error("PLASMA_dpotrf_Tile_Async", "invalid descriptor");
//         return plasma_request_fail(sequence, request,
//         PLASMA_ERR_ILLEGAL_VALUE);
//     } else {
//         descA = *A;
//     }
//     /* Check input arguments */
//     if (descA.nb != descA.mb) {
//         plasma_error("PLASMA_dpotrf_Tile_Async", "only square tiles
//         supported"); return plasma_request_fail(sequence, request,
//         PLASMA_ERR_ILLEGAL_VALUE);
//     }
//     if (uplo != PlasmaUpper && uplo != PlasmaLower) {
//         plasma_error("PLASMA_dpotrf_Tile_Async", "illegal value of uplo");
//         return plasma_request_fail(sequence, request, -1);
//     }
//     /* Quick return */
//     /*
//         if (max(N, 0) == 0)
//             return PLASMA_SUCCESS;
//     */
//     plasma_parallel_call_4(plasma_pdpotrf_gpu_reuse_data_static_table,
//                            PLASMA_enum, uplo,
//                            PLASMA_desc, descA,
//                            PLASMA_sequence*, sequence,
//                            PLASMA_request*, request);
//
//     return PLASMA_SUCCESS;
// }
#endif
#endif