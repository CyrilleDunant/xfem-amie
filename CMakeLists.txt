cmake_minimum_required (VERSION 2.6) 
project(AMIE)
include(CheckCXXCompilerFlag)
include (CheckCXXSourceCompiles)
find_package(OpenMP)
if( OpenMP_FOUND EQUAL 1 )
	include(FindOpenMP)
endif( OpenMP_FOUND EQUAL 1 )
find_package(GLUT)
if( GLUT_FOUND EQUAL 1 )
include(GLUT)
endif( GLUT_FOUND EQUAL 1 )

set(ALL_SRC 
elements/elements.cpp 
elements/integrable_entity.cpp 
features/feature_base.cpp
features/crack.cpp 
features/crackinitiation.cpp 
features/enrichmentbehaviour.cpp 
features/enrichmentInclusion3d.cpp 
features/enrichmentInclusion.cpp 
features/enrichmentRing.cpp 
features/expansiveRing.cpp 
features/expansiveZone3d.cpp 
features/expansiveZone.cpp 
features/features.cpp 
features/inclusion3d.cpp 
features/inclusion.cpp 
features/layeredinclusion.cpp 
features/pore3d.cpp 
features/pore.cpp 
features/sample3d.cpp 
features/sample.cpp 
features/vibratingcircularmembrane.cpp 
filters/voxelfilter.cpp 
filters/voxelporefilter.cpp 
geometry/geometry_2D.cpp 
geometry/geometry_3D.cpp 
geometry/geometry_base.cpp 
geometry/level_set.cpp
mesher/delaunay_3d.cpp 
mesher/delaunay.cpp 
mesher/structuredmesh.cpp 
physics/damagemodels/isotropiclineardamage.cpp 
physics/damagemodels/lineardamage.cpp 
physics/fracturecriteria/vonmises.cpp 
physics/diffusion.cpp 
physics/dual_behaviour.cpp 
physics/fracturecriteria/fracturecriterion.cpp 
physics/linearstiffnessgradient.cpp
physics/radialstiffnessgradient.cpp
physics/kelvinvoight.cpp
physics/fracturecriteria/maxstrain.cpp 
physics/fracturecriteria/mohrcoulomb.cpp 
physics/fracturecriteria/ruptureenergy.cpp 
physics/radial_distributed_stiffness.cpp
physics/spatially_distributed_stiffness.cpp
physics/stress_defined_stiffness.cpp
physics/homogeneised_behaviour.cpp 
physics/stiffness_and_fracture.cpp 
physics/stiffness.cpp 
physics/stiffness_with_imposed_deformation.cpp 
physics/stiffness_with_imposed_deformation_and_fracture.cpp 
physics/stiffness_with_variable_imposed_deformation_and_fracture.cpp 
physics/stiffness_with_variable_imposed_deformation.cpp 
physics/physics_base.cpp
physics/physics.cpp
physics/triple_behaviour.cpp 
physics/void_form.cpp 
physics/wave.cpp 
physics/weibull_distributed_stiffness.cpp 
physics/weibull_distributed_stiffness_with_variable_imposed_deformation_and_fracture.cpp 
polynomial/scalarfield.cpp 
polynomial/vm_base.cpp 
polynomial/vm_function_base.cpp 
polynomial/vm_function_matrix.cpp 
polynomial/vm_refcount_token.cpp 
polynomial/vm_token.cpp 
solvers/assembly.cpp 
solvers/biconjugategradientstabilized.cpp 
solvers/choleskidecomposed.cpp 
solvers/conjugategradient.cpp
solvers/eigenvalues.cpp
solvers/gausseidell.cpp
solvers/gaussseidellstep.cpp
solvers/incompletecholeskidecomposition.cpp
solvers/inversediagonal.cpp
solvers/polakribiereconjugategradient.cpp
solvers/preconditionners.cpp
solvers/solver.cpp
solvers/multigrid.cpp
solvers/tridiagonal.cpp
sparse/sparse_matrix.cpp
sparse/sparse_vector.cpp
utilities/configuration.cpp
utilities/grid.cpp
utilities/granulo.cpp
utilities/homogenizer.cpp
utilities/matrixops.cpp
utilities/parser.cpp
utilities/pid.cpp
utilities/placement.cpp
utilities/samplingcriterion.cpp
utilities/xml.cpp
)

find_package(GLUT)
find_package(OpenMP)

message (STATUS "checking for C++ TR1 extensions")
try_run (TR1_TEST_EXITCODE TR1_TEST_COMPILED
  ${CMAKE_BINARY_DIR}
  ${PROJECT_SOURCE_DIR}/cmake/tr1-test.cpp
)

if(TR1_TEST_COMPILED)
	if (TR1_TEST_EXITCODE EQUAL 0)
		message (STATUS "  C++ TR1 extensions found")
		set (TR1_FOUND 1)
		add_definitions(-DHAVE_TR1)
	else (TR1_TEST_EXITCODE EQUAL 0)
		message (STATUS "  C++ TR1 extensions not found")
		set (TR1_FOUND 0)
	endif (TR1_TEST_EXITCODE EQUAL 0)
else()
	message (STATUS "  C++ TR1 extensions not found")
        set (TR1_FOUND 0)
endif(TR1_TEST_COMPILED)

CHECK_CXX_COMPILER_FLAG(-ftree-vectorize COMPILER_HAS_VECTORISER )
CHECK_CXX_COMPILER_FLAG(-floop-strip-mine COMPILER_HAS_LOOP_INTERCHANGE )
CHECK_CXX_COMPILER_FLAG(-msse4a COMPILER_HAS_SSE4A )
CHECK_CXX_COMPILER_FLAG(-msse4.1 COMPILER_HAS_SSE41 )
CHECK_CXX_COMPILER_FLAG(-msse4.2 COMPILER_HAS_SSE42 )
CHECK_CXX_COMPILER_FLAG(-msse3 COMPILER_HAS_SSE3 )

 if(COMPILER_HAS_SSE3)
 	add_definitions("-msse3")
 	try_run (SSE3_TEST_EXITCODE SSE3_TEST_COMPILED
 		${CMAKE_BINARY_DIR}
 		${PROJECT_SOURCE_DIR}/cmake/sse3-test.cpp
 		COMPILE_DEFINITIONS "-msse3"
 	)
 
 	if (SSE3_TEST_EXITCODE EQUAL 0)
 		message (STATUS "  SSE3 found")
 		set (SSE3_FOUND 1)
 		add_definitions(-DHAVE_SSE3)
 	else (SSE3_TEST_EXITCODE EQUAL 1)
 		message (STATUS "  SSE3 not found")
 		set (SSE3_FOUND 0)
 	endif (SSE3_TEST_EXITCODE EQUAL 0)
 
	if(COMPILER_HAS_SSE4A)
		add_definitions("-msse4a")
		try_run (SSE4_TEST_EXITCODE SSE4_TEST_COMPILED
			${CMAKE_BINARY_DIR}
			${PROJECT_SOURCE_DIR}/cmake/sse4-test.cpp
			COMPILE_DEFINITIONS "-msse4a"
		)

		if (SSE4_TEST_EXITCODE EQUAL 0)
			message (STATUS "  SSE 4a found")
			set (SSE4_FOUND 1)
			add_definitions(-DHAVE_SSE4)
		else (SSE4_TEST_EXITCODE EQUAL 1)
			message (STATUS "  SSE 4a not found")
			set (SSE4_FOUND 0)
		endif (SSE4_TEST_EXITCODE EQUAL 0)

	endif(COMPILER_HAS_SSE4A)
	if(COMPILER_HAS_SSE41)
		add_definitions("-msse4.1")
		try_run (SSE4_TEST_EXITCODE SSE4_TEST_COMPILED
			${CMAKE_BINARY_DIR}
			${PROJECT_SOURCE_DIR}/cmake/sse4-test.cpp
			COMPILE_DEFINITIONS "-msse4.1"
		)

		if (SSE4_TEST_EXITCODE EQUAL 0)
			message (STATUS "  SSE 4.1 found")
			set (SSE4_FOUND 1)
			add_definitions(-DHAVE_SSE4)
		else (SSE4_TEST_EXITCODE EQUAL 1)
			message (STATUS "  SSE 4.1 not found")
			set (SSE4_FOUND 0)
		endif (SSE4_TEST_EXITCODE EQUAL 0)

	endif(COMPILER_HAS_SSE41)
	if(COMPILER_HAS_SSE42)
		add_definitions("-msse4.2")
		try_run (SSE4_TEST_EXITCODE SSE4_TEST_COMPILED
			${CMAKE_BINARY_DIR}
			${PROJECT_SOURCE_DIR}/cmake/sse4-test.cpp
			COMPILE_DEFINITIONS "-msse4.2"
		)

		if (SSE4_TEST_EXITCODE EQUAL 0)
			message (STATUS "  SSE 4.2 found")
			set (SSE4_FOUND 1)
			add_definitions(-DHAVE_SSE4)
		else (SSE4_TEST_EXITCODE EQUAL 1)
			message (STATUS "  SSE 4.2 not found")
			set (SSE4_FOUND 0)
		endif (SSE4_TEST_EXITCODE EQUAL 0)

	endif(COMPILER_HAS_SSE42)
 endif(COMPILER_HAS_SSE3)

if(COMPILER_HAS_VECTORISER)
	add_definitions(-ftree-vectorize -ftree-vectorizer-verbose=0)
endif(COMPILER_HAS_VECTORISER)

if(COMPILER_HAS_LOOP_INTERCHANGE)
	add_definitions(-floop-strip-mine -floop-block -floop-interchange)
endif(COMPILER_HAS_LOOP_INTERCHANGE)

# openMP is terrible with g++ (as of 4.4). Unfortunately, icc pretends it is g++
if(DISABLE_OPENMP EQUAL 1)
	set (OPENMP_FOUND 0) 
	message (STATUS "  Disabling OpenMP")
endif(DISABLE_OPENMP EQUAL 1)

if(OPENMP_FOUND)
	CHECK_CXX_COMPILER_FLAG(-pthread COMPILER_HAS_PTHREAD )
# 	if(COMPILER_HAS_PTHREAD AND CMAKE_COMPILER_IS_GNUCXX)
		add_definitions(${OpenMP_CXX_FLAGS} " -parallel -DHAVE_OPENMP")
# 		set(OMP_LINK_FLAGS "-pthread -lgomp")
# 	else()
		set(OMP_LINK_FLAGS ${OpenMP_CXX_FLAGS})
# 	endif()
endif(OPENMP_FOUND)

if(DEBUG_BUILD EQUAL 1)
	add_definitions(-g3 -fno-exceptions -O0 -DNDEBUG -march=native -funroll-loops)
	message (STATUS "  Debug Build")
else()
	add_definitions(-g3 -fno-exceptions -O3 -DNDEBUG -march=native -funroll-loops)
	message (STATUS "  Optimised Build")
endif(DEBUG_BUILD EQUAL 1)

#if we are under windows, we'll assume freeglut...
if(CMAKE_HOST_WIN32)
	set(GLUT_LIBRARY "-lfreeglut -lopengl32 -lglu32")
endif(CMAKE_HOST_WIN32)

add_library (Amie ${ALL_SRC})

add_executable(small examples/main_small.cpp )
target_link_libraries( small Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(3d  examples/main_3d.cpp  )
target_link_libraries( 3d Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(3d_diffusion  examples/main_3d_diffusion.cpp )
target_link_libraries( 3d_diffusion Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(cl examples/main_cl.cpp )
target_link_libraries( cl Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(example examples/main.cpp )
target_link_libraries( example Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(def examples/main_def.cpp )
target_link_libraries( def Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(kill examples/main_kill.cpp )
target_link_libraries( kill Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(multicrack examples/main_multicrack.cpp )
target_link_libraries( multicrack Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(alain examples/main_alain.cpp )
target_link_libraries( alain Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(rag examples/main_rag.cpp )
target_link_libraries( rag Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(rag_simple examples/main_rag_simple.cpp )
target_link_libraries( rag_simple Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(ellipse_test examples/main_ellipse_test.cpp )
target_link_libraries( ellipse_test Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(thermal_expansion examples/main_thermalexpansion.cpp )
target_link_libraries( thermal_expansion Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(ellipse examples/main_ellipse.cpp )
target_link_libraries( ellipse Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(benchmark examples/main_benchmark.cpp )
target_link_libraries( benchmark Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(rag_ellipse examples/main_ellipse_rag.cpp )
target_link_libraries( rag_ellipse Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(iterative examples/main_iterative.cpp )
target_link_libraries( iterative Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(simple examples/main_simple.cpp )
target_link_libraries( simple Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(tripoint examples/main_tripoint.cpp )
target_link_libraries( tripoint Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable(visco examples/main_visco.cpp )
target_link_libraries( visco Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable( wave examples/main_wave.cpp )
target_link_libraries( wave Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)

add_executable( cyrille_test examples/main_cyrille.cpp )
target_link_libraries( cyrille_test Amie ${GLUT_LIBRARY} -lglut -lGLU ${OMP_LINK_FLAGS} m)


