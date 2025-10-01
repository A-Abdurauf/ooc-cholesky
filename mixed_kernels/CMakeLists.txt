add_library(mixed_kernels STATIC casting.cu gemm.cpp syrk.cpp)
target_link_libraries(mixed_kernels PRIVATE CUDA::cudart CUDA::cublas)
target_include_directories(mixed_kernels PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(mixed_kernels PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)