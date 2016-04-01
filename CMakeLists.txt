cmake_minimum_required (VERSION 2.6) 
project(AMIE)
include(CheckCXXCompilerFlag)
include (CheckCXXSourceCompiles)
find_package(OpenMP)
find_package(MPI)

if (CMAKE_CXX_COMPILER MATCHES ".*clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_COMPILER_IS_CLANGXX 1)
endif ()

if( OpenMP_FOUND EQUAL 1 )
	include(FindOpenMP)
endif( OpenMP_FOUND EQUAL 1 )
    
include_directories(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_PATH})

if (USE_MPI EQUAL 0)
    set (MPI_COMPILE_FLAGS "")
    set (MPI_LINK_FLAGS "")
else()
    set(CXX MPI_CXX_COMPILER)
    add_definitions(${MPI_COMPILE_FLAGS} ${Boost_COMPILE_FLAGS})
endif(USE_MPI EQUAL 0)

set(ALL_SRC 
utilities/matrixops.cpp
geometry/geometry_base.cpp 
elements/elements.cpp 
elements/integrable_entity.cpp 
elements/generalized_spacetime_viscoelastic_element_state.cpp 
features/feature_base.cpp
features/crack.cpp 
features/crackinitiation.cpp 
features/enrichmentmanagers/crackmanager.cpp 
features/enrichmentmanagers/gelmanager.cpp 
features/enrichmentbehaviour.cpp 
features/enrichmentInclusion3d.cpp 
features/enrichmentInclusion.cpp 
features/timeDependentEnrichmentInclusion.cpp 
features/enrichmentRing.cpp 
features/expansiveRing.cpp 
features/expansiveZone3d.cpp 
features/boundarycondition.cpp 
features/expansiveZone.cpp 
features/features.cpp 
features/growingExpansiveZone.cpp
features/inclusion3d.cpp 
features/inclusion.cpp 
features/layeredinclusion.cpp 
features/pore3d.cpp 
features/pore.cpp 
features/sample3d.cpp 
features/sample.cpp
features/polygonSample.cpp
features/polygonSample3d.cpp
features/loftedPolygonSample3d.cpp
features/vibratingcircularmembrane.cpp 
filters/voxelfilter.cpp
filters/voxelporefilter.cpp 
polynomial/vm_base.cpp 
polynomial/vm_function_base.cpp 
polynomial/vm_function_extra.cpp 
polynomial/vm_function_matrix.cpp 
polynomial/vm_token.cpp 
geometry/geometry_2D.cpp 
geometry/geometry_3D.cpp 
geometry/geometry_with_exclusion.cpp 
geometry/level_set.cpp
geometry/sampler/sampler.cpp
geometry/sampler/gradient_sampler.cpp
geometry/sampler/regular_sampler.cpp
mesher/delaunay_3d.cpp 
mesher/delaunay.cpp 
mesher/element_checker.cpp 
mesher/parallel_delaunay.cpp 
mesher/parallel_delaunay_3d.cpp 
mesher/structuredmesh.cpp 
physics/fracturecriteria/vonmises.cpp 
physics/fracturecriteria/fractionmcft.cpp 
physics/fracturecriteria/confinedmohrcoulomb.cpp
physics/fracturecriteria/confinedmohrcoulombwithstrain.cpp
physics/fracturecriteria/confinedvonmises.cpp
physics/fracturecriteria/nonlocalvonmises.cpp 
physics/fracturecriteria/mazars.cpp 
physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.cpp 
physics/fracturecriteria/spacetimemultisurfacefracturecriterion.cpp 
physics/fracturecriteria/spacetimelimitsurfacefracturecriterion.cpp 
physics/fracturecriteria/spacetimeflowrule.cpp 
physics/fracturecriteria/druckerprager.cpp 
physics/fracturecriteria/boundedvonmises.cpp 
physics/fracturecriteria/creeprupture.cpp 
physics/damagemodels/damagemodel.cpp
physics/damagemodels/isotropiclineardamage.cpp 
physics/damagemodels/nonlocalisotropiclineardamage.cpp 
physics/damagemodels/fiberbasedisotropiclineardamage.cpp 
physics/damagemodels/spacetimefiberbasedisotropiclineardamage.cpp 
physics/damagemodels/spacetimebadisotropiclineardamage.cpp 
physics/damagemodels/spacetimefiberbasedbilineardamage.cpp 
physics/damagemodels/spacetimeratebasedisotropiclineardamage.cpp 
physics/damagemodels/spacetimeisotropiclineardamage.cpp 
physics/damagemodels/anisotropicdamage.cpp
physics/damagemodels/fractiondamage.cpp
physics/damagemodels/strainbrokenisotropiclineardamage.cpp
physics/damagemodels/lineardamage.cpp 
physics/damagemodels/plasticstrain.cpp
physics/damagemodels/rotatingcrack.cpp
physics/material_laws/logcreep_accumulator.cpp
physics/material_laws/material_laws.cpp
physics/material_laws/humidity_material_laws.cpp
physics/material_laws/mechanical_material_laws.cpp
physics/material_laws/temperature_material_laws.cpp
#physics/mad.cpp 
physics/diffusion.cpp 
physics/mass.cpp 
physics/laplacian.cpp 
physics/logarithmic_creep.cpp 
physics/logarithmic_creep_with_external_parameters.cpp 
physics/dual_behaviour.cpp 
physics/viscoelasticity.cpp 
physics/viscoelasticity_and_fracture.cpp 
physics/viscoelasticity_and_imposed_deformation.cpp 
physics/finite_difference_viscoelasticity.cpp 
physics/fracturecriteria/fracturecriterion.cpp 
physics/fracturecriteria/limitstrains.cpp 
physics/linearstiffnessgradient.cpp
physics/radialstiffnessgradient.cpp
physics/fracturecriteria/maxstrain.cpp 
physics/fracturecriteria/mohrcoulomb.cpp
physics/fracturecriteria/mcft.cpp
physics/fracturecriteria/mazars.cpp
physics/non_linear_stiffness.cpp 
physics/fracturecriteria/ruptureenergy.cpp 
physics/radial_distributed_stiffness.cpp
physics/spatially_distributed_stiffness.cpp
physics/homogeneised_behaviour.cpp 
physics/stiffness_and_fracture.cpp 
physics/stiffness.cpp 
#physics/maxwell.cpp 
#physics/maxwell_and_fracture.cpp 
#physics/parallel_behaviour.cpp REVIRON
#physics/serial_behaviour.cpp 
physics/orthotropicstiffness.cpp 
physics/stiffness_and_indexed_fracture.cpp
physics/fraction_stiffness_and_fracture.cpp
physics/stiffness_with_imposed_deformation.cpp 
physics/stiffness_with_imposed_stress.cpp 
physics/stiffness_with_imposed_deformation_and_fracture.cpp 
physics/stiffness_with_imposed_stress_and_fracture.cpp 
physics/stiffness_with_kinetics_imposed_deformation.cpp
physics/stiffness_with_variable_imposed_deformation_and_fracture.cpp 
physics/stiffness_with_diffusion_deformation_and_fracture.cpp
physics/stiffness_with_variable_imposed_deformation.cpp 
physics/physics_base.cpp
physics/twod_cohesive_force.cpp
physics/triple_behaviour.cpp 
physics/void_form.cpp 
physics/wave.cpp 
physics/weibull_distributed_stiffness.cpp 
physics/weibull_distributed_stiffness_with_variable_imposed_deformation_and_fracture.cpp 
physics/homogenization/phase.cpp
physics/homogenization/composite.cpp
physics/materials/aggregate_behaviour.cpp
physics/materials/gel_behaviour.cpp
physics/materials/paste_behaviour.cpp
physics/materials/rebar_behaviour.cpp
physics/materials/steel_behaviour.cpp
physics/materials/concrete_behaviour.cpp
physics/materials/c3s_behaviour.cpp
physics/materials/csh_behaviour.cpp
physics/materials/ch_behaviour.cpp
polynomial/scalarfield.cpp 
solvers/assembly.cpp 
geometry/space_time_geometry_2D.cpp 
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
solvers/multigridstep.cpp
solvers/tridiagonal.cpp
solvers/ssor.cpp
sparse/sparse_matrix.cpp
sparse/sparse_vector.cpp
utilities/configuration.cpp
utilities/grid.cpp
utilities/optimizer.cpp
utilities/granulo.cpp
utilities/inclusion_family.cpp
utilities/tensor.cpp
utilities/parser.cpp
utilities/pid.cpp
utilities/itoa.cpp
utilities/placement.cpp
utilities/behaviour_converter.cpp
utilities/postprocessor.cpp
utilities/samplingcriterion.cpp
features/microstructuregenerator.cpp
utilities/writer/triangle_writer.cpp
utilities/writer/voxel_writer.cpp
utilities/writer/exodus_writer.cpp
utilities/object_translator.cpp
utilities/mineral.cpp
)


CHECK_CXX_COMPILER_FLAG(-ftree-vectorize COMPILER_HAS_VECTORISER )
CHECK_CXX_COMPILER_FLAG(-floop-strip-mine COMPILER_HAS_LOOP_INTERCHANGE )
CHECK_CXX_COMPILER_FLAG(-msse4a COMPILER_HAS_SSE4A )
CHECK_CXX_COMPILER_FLAG(-msse4.1 COMPILER_HAS_SSE41 )
CHECK_CXX_COMPILER_FLAG(-msse4.2 COMPILER_HAS_SSE42 )
CHECK_CXX_COMPILER_FLAG(-msse3 COMPILER_HAS_SSE3 )

if(BUILD_32 EQUAL 1)
	set(BIT32_FLAGS "-m32")
	message (STATUS "  32 bit build")
endif(BUILD_32 EQUAL 1)


if(DEBUG_BUILD EQUAL 0)
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
endif(DEBUG_BUILD EQUAL 0)

if(COMPILER_HAS_VECTORISER)
	add_definitions(-ftree-vectorize )
endif(COMPILER_HAS_VECTORISER)

if(COMPILER_HAS_LOOP_INTERCHANGE)
	add_definitions(-floop-strip-mine -floop-block -floop-interchange)
endif(COMPILER_HAS_LOOP_INTERCHANGE)

if(DISABLE_OPENMP EQUAL 1)
	set (OPENMP_FOUND 0) 
	message (STATUS "  Disabling OpenMP")
endif(DISABLE_OPENMP EQUAL 1)

if(OPENMP_FOUND)
	CHECK_CXX_COMPILER_FLAG(-pthread COMPILER_HAS_PTHREAD )
message (STATUS "  found OpenMP")
# 	if(COMPILER_HAS_PTHREAD AND CMAKE_COMPILER_IS_GNUCXX)
		add_definitions(${OpenMP_CXX_FLAGS} " -DHAVE_OPENMP ")
# 		set(OMP_LINK_FLAGS "-pthread -lgomp")
# 	else()
		set(OMP_LINK_FLAGS ${OpenMP_CXX_FLAGS})
# 	endif()
endif(OPENMP_FOUND)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANGXX)
	add_definitions ( "-std=c++11" )
endif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANGXX)

if(DEBUG_BUILD EQUAL 1)
#     set(CMAKE_BUILD_TYPE Debug)
	add_definitions(-g3 -Og -DNDEBUG -funroll-loops -Wall -Wpedantic )
	message (STATUS "  Debug Build")
else()
#     set(CMAKE_BUILD_TYPE Release)
	add_definitions(-g3 -O3 ${BIT32_FLAGS} -DNDEBUG -funroll-loops -Wall -Wpedantic )
	message (STATUS "  Optimised Build")
endif(DEBUG_BUILD EQUAL 1)

set(ALL_LINK_FLAGS  ${OMP_LINK_FLAGS} m ${BIT32_FLAGS} ${MPI_LINK_FLAGS} ${MPI_LIBRARIES} ${Boost_LINK_FLAGS} ${Boost_LIBRARIES})

add_library (Amie ${ALL_SRC})

add_executable(dam examples/main_dam.cpp )
target_link_libraries( dam Amie ${ALL_LINK_FLAGS} )

add_executable(back_homogenisation examples/main_back_homogenisation.cpp )
target_link_libraries( back_homogenisation Amie ${ALL_LINK_FLAGS} )

add_executable(visco_CSH_CASH examples/main_visco_CSH_CASH.cpp )
target_link_libraries( visco_CSH_CASH Amie ${ALL_LINK_FLAGS} )

add_executable(single_2d_inclusion examples/main_single_2d_inclusion.cpp )
target_link_libraries( single_2d_inclusion Amie ${ALL_LINK_FLAGS} )

add_executable(single_2d_viscoelasticity examples/main_single_2d_viscoelasticity.cpp )
target_link_libraries( single_2d_viscoelasticity Amie ${ALL_LINK_FLAGS} )

add_executable(3d_asr  examples/main_3d_asr.cpp  )
target_link_libraries( 3d_asr Amie  ${ALL_LINK_FLAGS})

add_executable(hyd_paste examples/main_hydrating_paste.cpp )
target_link_libraries( hyd_paste Amie ${ALL_LINK_FLAGS})

add_executable(def examples/main_def.cpp )
target_link_libraries( def Amie ${ALL_LINK_FLAGS})

add_executable(multicrack examples/main_multicrack.cpp )
target_link_libraries( multicrack Amie ${ALL_LINK_FLAGS})

add_executable(steel examples/main_steel.cpp )
target_link_libraries( steel Amie ${ALL_LINK_FLAGS})

add_executable(asr examples/main_asr.cpp )
target_link_libraries( asr Amie ${ALL_LINK_FLAGS})

add_executable(asr_simple examples/main_asr_simple.cpp )
target_link_libraries( asr_simple Amie ${ALL_LINK_FLAGS})

add_executable(wedge examples/main_wedge_splitting.cpp )
target_link_libraries( wedge Amie ${ALL_LINK_FLAGS})

add_executable(geometry examples/main_geometry.cpp )
target_link_libraries( geometry Amie ${ALL_LINK_FLAGS} )

add_executable(2d_composite examples/main_2d_composite.cpp )
target_link_libraries( 2d_composite Amie ${ALL_LINK_FLAGS} )

add_executable(asr_shape examples/main_asr_shape.cpp )
target_link_libraries( asr_shape Amie ${ALL_LINK_FLAGS})

add_executable(tripoint examples/main_tripoint.cpp )
target_link_libraries( tripoint Amie ${ALL_LINK_FLAGS})

add_executable(tripoint_simple examples/main_tripoint_simple.cpp )
target_link_libraries( tripoint_simple Amie ${ALL_LINK_FLAGS})

add_executable(irradiation examples/main_irradiation.cpp )
target_link_libraries( irradiation Amie ${ALL_LINK_FLAGS})

add_executable(asr_visco examples/main_asr_visco.cpp )
target_link_libraries( asr_visco Amie ${ALL_LINK_FLAGS} )

#add_executable(from_mic examples/main_from_mic.cpp )
#target_link_libraries( from_mic Amie ${ALL_LINK_FLAGS} )

add_executable(tension examples/main_tension.cpp )
target_link_libraries( tension Amie ${ALL_LINK_FLAGS})

add_executable( pore_pressure examples/main_crystal_pressure.cpp )
target_link_libraries( pore_pressure Amie ${ALL_LINK_FLAGS})

add_executable( visu examples/main_visu_behavior.cpp )
target_link_libraries( visu Amie ${ALL_LINK_FLAGS})

add_executable( visco examples/main_visco.cpp )
target_link_libraries( visco Amie ${ALL_LINK_FLAGS})

add_executable( tension_visco examples/main_tension_visco.cpp )
target_link_libraries( tension_visco Amie ${ALL_LINK_FLAGS})

add_executable( tension_visco_damage examples/main_tension_visco_damage.cpp )
target_link_libraries( tension_visco_damage Amie ${ALL_LINK_FLAGS})

# benchmarks

add_executable(3d_benchmark EXCLUDE_FROM_ALL examples/main_3d_benchmark.cpp )
target_link_libraries( 3d_benchmark Amie ${ALL_LINK_FLAGS} )

add_executable(tension_benchmark EXCLUDE_FROM_ALL examples/main_tension_benchmark.cpp )
target_link_libraries( tension_benchmark Amie ${ALL_LINK_FLAGS})

add_executable( performance EXCLUDE_FROM_ALL examples/main_performance.cpp )
target_link_libraries( performance Amie ${ALL_LINK_FLAGS})

ADD_CUSTOM_TARGET(benchmark)
ADD_DEPENDENCIES(benchmark 
3d_benchmark 
tension_benchmark 
performance
)


# playground

add_executable(alain EXCLUDE_FROM_ALL examples/playground/main_alain.cpp )
target_link_libraries( alain Amie ${ALL_LINK_FLAGS})

add_executable(michaela EXCLUDE_FROM_ALL examples/playground/main_michaela.cpp )
target_link_libraries( michaela Amie ${ALL_LINK_FLAGS})

add_executable(function_test EXCLUDE_FROM_ALL examples/playground/main_function_test.cpp )
target_link_libraries( function_test Amie ${ALL_LINK_FLAGS})

add_executable(get_fields_from_trg EXCLUDE_FROM_ALL examples/playground/main_get_fields_from_trg.cpp )
target_link_libraries( get_fields_from_trg Amie ${ALL_LINK_FLAGS} )

add_executable(homogenization EXCLUDE_FROM_ALL examples/playground/main_homogenization.cpp )
target_link_libraries( homogenization Amie ${ALL_LINK_FLAGS})

ADD_CUSTOM_TARGET(playground)
ADD_DEPENDENCIES(playground 
alain 
michaela 
function_test
get_fields_from_trg 
homogenization
)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test ${CMAKE_BINARY_DIR}/mesh)

# tests

add_executable(check_test EXCLUDE_FROM_ALL examples/main_check_test.cpp )
target_link_libraries( check_test Amie ${ALL_LINK_FLAGS} )

add_executable(check_behaviour EXCLUDE_FROM_ALL examples/main_check_behaviour.cpp )
target_link_libraries( check_behaviour Amie ${ALL_LINK_FLAGS} )

add_executable( test/cauchy_green EXCLUDE_FROM_ALL examples/test/main_test_cauchy_green.cpp )
target_link_libraries( test/cauchy_green Amie ${ALL_LINK_FLAGS})

add_executable( test/spacetime_bilateral_damage EXCLUDE_FROM_ALL examples/test/main_test_spacetime_bilateral_damage.cpp )
target_link_libraries( test/spacetime_bilateral_damage Amie ${ALL_LINK_FLAGS})

add_executable( test/stxfem EXCLUDE_FROM_ALL examples/test/main_test_stxfem.cpp )
target_link_libraries( test/stxfem Amie ${ALL_LINK_FLAGS})

add_executable( test/time_continuity EXCLUDE_FROM_ALL examples/test/main_test_time_continuity.cpp )
target_link_libraries( test/time_continuity Amie ${ALL_LINK_FLAGS})

add_executable( test/viscodamage_multilinearsoftening EXCLUDE_FROM_ALL examples/test/main_test_viscodamage_multilinearsoftening.cpp )
target_link_libraries( test/viscodamage_multilinearsoftening Amie ${ALL_LINK_FLAGS})

add_executable( test/viscodamage_multisurface_criterion EXCLUDE_FROM_ALL examples/test/main_test_viscodamage_multisurface_criterion.cpp )
target_link_libraries( test/viscodamage_multisurface_criterion Amie ${ALL_LINK_FLAGS})

add_executable( test/xfem EXCLUDE_FROM_ALL examples/test/main_test_xfem.cpp )
target_link_libraries( test/xfem Amie ${ALL_LINK_FLAGS})


ADD_CUSTOM_TARGET(unitTests)
ADD_DEPENDENCIES(unitTests 
check_test 
check_behaviour 
test/cauchy_green
test/spacetime_bilateral_damage
test/stxfem 
test/time_continuity
test/viscodamage_multilinearsoftening 
test/viscodamage_multisurface_criterion
test/xfem 
)

# mesh

add_executable(check_mesh EXCLUDE_FROM_ALL examples/main_check_mesh.cpp )
target_link_libraries( check_mesh Amie ${ALL_LINK_FLAGS} )

add_executable( mesh/circle EXCLUDE_FROM_ALL examples/mesh/main_mesh_circle.cpp )
target_link_libraries( mesh/circle Amie ${ALL_LINK_FLAGS})

add_executable( mesh/layered_circles EXCLUDE_FROM_ALL examples/mesh/main_mesh_layered_circles.cpp )
target_link_libraries( mesh/layered_circles Amie ${ALL_LINK_FLAGS})

add_executable( mesh/circle_circle_intersection EXCLUDE_FROM_ALL examples/mesh/main_mesh_circle_circle_intersection.cpp )
target_link_libraries( mesh/circle_circle_intersection Amie ${ALL_LINK_FLAGS})

add_executable( mesh/circle_circle_mask EXCLUDE_FROM_ALL examples/mesh/main_mesh_circle_circle_mask.cpp )
target_link_libraries( mesh/circle_circle_mask Amie ${ALL_LINK_FLAGS})

add_executable( mesh/circle_ellipse_polygon EXCLUDE_FROM_ALL examples/mesh/main_mesh_circle_ellipse_polygon.cpp )
target_link_libraries( mesh/circle_ellipse_polygon Amie ${ALL_LINK_FLAGS})

add_executable( mesh/circle_ellipse_intersection EXCLUDE_FROM_ALL examples/mesh/main_mesh_circle_ellipse_intersection.cpp )
target_link_libraries( mesh/circle_ellipse_intersection Amie ${ALL_LINK_FLAGS})

add_executable( mesh/circle_ellipse_mask EXCLUDE_FROM_ALL examples/mesh/main_mesh_circle_ellipse_mask.cpp )
target_link_libraries( mesh/circle_ellipse_mask Amie ${ALL_LINK_FLAGS})

add_executable( mesh/circle_sample_intersection EXCLUDE_FROM_ALL examples/mesh/main_mesh_circle_sample_intersection.cpp )
target_link_libraries( mesh/circle_sample_intersection Amie ${ALL_LINK_FLAGS})

add_executable( mesh/concrete_circle EXCLUDE_FROM_ALL examples/mesh/main_mesh_concrete_circle.cpp )
target_link_libraries( mesh/concrete_circle Amie ${ALL_LINK_FLAGS})

add_executable( mesh/concrete_crushed EXCLUDE_FROM_ALL examples/mesh/main_mesh_concrete_crushed.cpp )
target_link_libraries( mesh/concrete_crushed Amie ${ALL_LINK_FLAGS})

add_executable( mesh/concrete_ellipse EXCLUDE_FROM_ALL examples/mesh/main_mesh_concrete_ellipse.cpp )
target_link_libraries( mesh/concrete_ellipse Amie ${ALL_LINK_FLAGS})

add_executable( mesh/concrete_gravel EXCLUDE_FROM_ALL examples/mesh/main_mesh_concrete_gravel.cpp )
target_link_libraries( mesh/concrete_gravel Amie ${ALL_LINK_FLAGS})

add_executable( mesh/ellipse EXCLUDE_FROM_ALL examples/mesh/main_mesh_ellipse.cpp )
target_link_libraries( mesh/ellipse Amie ${ALL_LINK_FLAGS})

add_executable( mesh/ellipse_ellipse_intersection EXCLUDE_FROM_ALL examples/mesh/main_mesh_ellipse_ellipse_intersection.cpp )
target_link_libraries( mesh/ellipse_ellipse_intersection Amie ${ALL_LINK_FLAGS})

add_executable( mesh/ellipse_ellipse_mask EXCLUDE_FROM_ALL examples/mesh/main_mesh_ellipse_ellipse_mask.cpp )
target_link_libraries( mesh/ellipse_ellipse_mask Amie ${ALL_LINK_FLAGS})

add_executable( mesh/ellipse_sample_intersection EXCLUDE_FROM_ALL examples/mesh/main_mesh_ellipse_sample_intersection.cpp )
target_link_libraries( mesh/ellipse_sample_intersection Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polycrystal EXCLUDE_FROM_ALL examples/mesh/main_mesh_polycrystal.cpp )
target_link_libraries( mesh/polycrystal Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polycrystal_itz EXCLUDE_FROM_ALL examples/mesh/main_mesh_polycrystal_itz.cpp )
target_link_libraries( mesh/polycrystal_itz Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polygon_regular EXCLUDE_FROM_ALL examples/mesh/main_mesh_polygon_regular.cpp )
target_link_libraries( mesh/polygon_regular Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polygon_convex EXCLUDE_FROM_ALL examples/mesh/main_mesh_polygon_convex.cpp )
target_link_libraries( mesh/polygon_convex Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polygon_polygon_edge EXCLUDE_FROM_ALL examples/mesh/main_mesh_polygon_polygon_edge.cpp )
target_link_libraries( mesh/polygon_polygon_edge Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polygon_circle_intersection EXCLUDE_FROM_ALL examples/mesh/main_mesh_polygon_circle_intersection.cpp )
target_link_libraries( mesh/polygon_circle_intersection Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polygon_circle_mask EXCLUDE_FROM_ALL examples/mesh/main_mesh_polygon_circle_mask.cpp )
target_link_libraries( mesh/polygon_circle_mask Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polygon_ellipse_intersection EXCLUDE_FROM_ALL examples/mesh/main_mesh_polygon_ellipse_intersection.cpp )
target_link_libraries( mesh/polygon_ellipse_intersection Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polygon_ellipse_mask EXCLUDE_FROM_ALL examples/mesh/main_mesh_polygon_ellipse_mask.cpp )
target_link_libraries( mesh/polygon_ellipse_mask Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polygon_polygon_intersection EXCLUDE_FROM_ALL examples/mesh/main_mesh_polygon_polygon_intersection.cpp )
target_link_libraries( mesh/polygon_polygon_intersection Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polygon_polygon_mask EXCLUDE_FROM_ALL examples/mesh/main_mesh_polygon_polygon_mask.cpp )
target_link_libraries( mesh/polygon_polygon_mask Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polygon_sample_edge EXCLUDE_FROM_ALL examples/mesh/main_mesh_polygon_sample_edge.cpp )
target_link_libraries( mesh/polygon_sample_edge Amie ${ALL_LINK_FLAGS})

add_executable( mesh/polygon_sample_intersection EXCLUDE_FROM_ALL examples/mesh/main_mesh_polygon_sample_intersection.cpp )
target_link_libraries( mesh/polygon_sample_intersection Amie ${ALL_LINK_FLAGS})

add_executable( mesh/sample EXCLUDE_FROM_ALL examples/mesh/main_mesh_sample.cpp )
target_link_libraries( mesh/sample Amie ${ALL_LINK_FLAGS})

add_executable( mesh/sample_gradient EXCLUDE_FROM_ALL examples/mesh/main_mesh_sample_gradient.cpp )
target_link_libraries( mesh/sample_gradient Amie ${ALL_LINK_FLAGS})

add_executable( mesh/sample_weibull EXCLUDE_FROM_ALL examples/mesh/main_mesh_sample_weibull.cpp )
target_link_libraries( mesh/sample_weibull Amie ${ALL_LINK_FLAGS})

ADD_CUSTOM_TARGET(meshTests)
ADD_DEPENDENCIES(meshTests 
check_mesh 
mesh/circle 
mesh/circle_circle_intersection 
mesh/circle_circle_mask 
mesh/circle_ellipse_intersection 
mesh/circle_ellipse_mask
mesh/circle_ellipse_polygon
mesh/circle_sample_intersection 
mesh/concrete_circle
mesh/concrete_crushed
mesh/concrete_ellipse
mesh/concrete_gravel
mesh/ellipse 
mesh/ellipse_ellipse_intersection 
mesh/ellipse_ellipse_mask 
mesh/ellipse_sample_intersection
mesh/layered_circles 
mesh/polycrystal
mesh/polycrystal_itz
mesh/polygon_circle_intersection
mesh/polygon_circle_mask
mesh/polygon_convex
mesh/polygon_ellipse_intersection
mesh/polygon_ellipse_mask
mesh/polygon_polygon_edge
mesh/polygon_polygon_intersection
mesh/polygon_polygon_mask
mesh/polygon_regular
mesh/polygon_sample_edge
mesh/polygon_sample_intersection
mesh/sample
mesh/sample_gradient
mesh/sample_weibull
)

if(EXISTS ${CMAKE_SOURCE_DIR}/local.cmake )
    message( STATUS "found local.cmake" )
    include(local.cmake)
else()
    message( STATUS "file local.cmake not found, nothing to do" )
endif()

