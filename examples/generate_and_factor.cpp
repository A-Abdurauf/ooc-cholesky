#include <api/ExaGeoStat.hpp>
#include <common/Definitions.hpp>
#include <plasma.h>
#include <plasma_d_mixed.h>
#include <chameleon.h>
#include <Eigen/Dense>
#include <vector>

using namespace exageostat::api;
using namespace exageostat::configurations;
using namespace exageostat::common;

void generate_and_factor() {
    const int N = 1024;                 // problem size (points), adjust as needed
    Configurations cfg;
    cfg.SetProblemSize(N);                 // number of points
    cfg.SetDimension(exageostat::common::Dimension2D);                // or Dimension3D
    cfg.SetKernelName("UnivariateMaternStationary"); // common climate kernel
    cfg.SetDenseTileSize(256);             // match PLASMA nb if you want
    cfg.SetSeed(1234);
    cfg.SetTimeSlot(1);                    // single time slot (default synthetic case)
    cfg.SetCoresNumber(2);                 // aligns with PLASMA_Init below
    // add bounds/initial_theta as needed

#if DEFAULT_RUNTIME
    auto hw = ExaGeoStatHardware(cfg.GetComputation(), cfg.GetCoresNumber(),
                                 cfg.GetGPUsNumbers(), cfg.GetPGrid(), cfg.GetQGrid());
#else
    auto hw = ExaGeoStatHardware(cfg);
#endif

    std::unique_ptr<ExaGeoStatData<double>> data;
    ExaGeoStat<double>::ExaGeoStatLoadData(cfg, data);   // fills covariance + Z

    const int n = cfg.GetProblemSize() * cfg.GetTimeSlot();

    // Export generated observations vector Z and build an SPD matrix: A = Z * Z^T + n*I
    auto *descZ = data->GetDescriptorData()
                      ->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z)
                      .chameleon_desc;
    std::vector<double> z(n);
    CHAMELEON_Desc2Lap(ChamUpperLower, descZ, z.data(), n);

    Eigen::Map<Eigen::VectorXd> zvec(z.data(), n);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> SPD = zvec * zvec.transpose();
    SPD.diagonal().array() += n; // shift to ensure positive definiteness

    std::vector<double> A(n * n);
    std::memcpy(A.data(), SPD.data(), sizeof(double) * A.size());

    // Factor with ooc-cholesky / PLASMA
    PLASMA_Init(cfg.GetCoresNumber());
    plasma_context_t *plasma = plasma_context_self();
    plasma->autotuning_enabled = 0;
    plasma->nb = cfg.GetDenseTileSize();   // keep tiles consistent

    const auto status = PLASMA_dpotrf_gpu_reuse_data_table_mixed_precision(
        PlasmaLower, n, A.data(), n);
    if (status != 0) { /* handle error */ }

    PLASMA_Finalize();
}

int main() {
    generate_and_factor();
    return 0;
}