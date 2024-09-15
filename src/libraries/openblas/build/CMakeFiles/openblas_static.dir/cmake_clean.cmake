file(REMOVE_RECURSE
  "lib/libopenblas.a"
  "lib/libopenblas.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang ASM C Fortran)
  include(CMakeFiles/openblas_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
