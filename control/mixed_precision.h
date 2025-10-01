#pragma once
#include <cublas_v2.h>

#ifdef __cplusplus
extern "C" {
#endif

struct MixedPrecisionTile {
  int dtype;
  void *data;
  int layout;
  size_t m, n, ld;
};

struct MixedPrecisionTiledArray {
  MixedPrecisionTile *tiles;
  int nt;
  int uplo;
  int nb, mb;
  size_t m, n;
};

size_t getSizeofTileElement(int dtype);

int getComputeType(int dtype);

void *castToBuffer(int m, int n, const void *source, cudaDataType sourceType,
                   int ldSource, void *target, cudaDataType targetType,
                   int ldTarget, cudaStream_t stream);

#ifdef __cplusplus
}
#endif