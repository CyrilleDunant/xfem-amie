# Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C} 2005-2007
#
# Copyright: See COPYING file that comes with this distribution
#

CXX = g++

# CXXFLAGS = -g3 -Wall -O2 -DNDEBUG 
#  CXXFLAGS = -g3 -Wall -O0 -DNDEBUG 
## use -O0 for mem-checking with valgrind
## valgrind --num-callers=8 ./your_target
# CXXFLAGS = -g3 -Wall -fno-exceptions -O0 -DNDEBUG -fstrict-aliasing
CXXFLAGS = -g3 -Wall -ftree-vectorize -fno-exceptions -O3 -DNDEBUG 

LDFLAGS = -Wl -lm -lGL -lglut

SOURCE_PHYSICS = physics/physics.cpp physics/physics_base.cpp physics/stiffness.cpp physics/diffusion.cpp physics/stiffness_with_imposed_deformation.cpp physics/weibull_distributed_stiffness.cpp physics/stiffness_and_fracture.cpp  physics/void_form.cpp physics/fracturecriterion.cpp physics/mohrcoulomb.cpp physics/vonmises.cpp physics/maxstrain.cpp physics/dual_behaviour.cpp physics/radialstiffnessgradient.cpp physics/linearstiffnessgradient.cpp physics/ruptureenergy.cpp physics/kelvinvoight.cpp
OBJECTS_PHYSICS = physics/physics.o physics/physics_base.o physics/stiffness.o physics/stiffness_and_fracture.o  physics/void_form.o physics/weibull_distributed_stiffness.o physics/stiffness_with_imposed_deformation.o physics/fracturecriterion.o physics/mohrcoulomb.o physics/vonmises.o physics/maxstrain.o physics/diffusion.o physics/dual_behaviour.o physics/radialstiffnessgradient.o  physics/linearstiffnessgradient.o physics/ruptureenergy.o physics/kelvinvoight.o

SOURCE_FILTERS = filters/voxelfilter.cpp filters/voxelporefilter.cpp
OBJECTS_FILTERS = filters/voxelfilter.o filters/voxelporefilter.o

SOURCE_ELEMENTS = elements/elements.cpp elements/integrable_entity.cpp
OBJECTS_ELEMENTS = elements/elements.o elements/integrable_entity.o

SOURCE_SOLVERS = solvers/assembly.cpp solvers/choleskidecomposed.cpp solvers/conjugategradient.cpp solvers/gausseidell.cpp solvers/polakribiereconjugategradient.cpp solvers/solver.cpp solvers/preconditionners.cpp solvers/gaussseidellstep.cpp solvers/inversediagonal.cpp solvers/incompletecholeskidecomposition.cpp solvers/biconjugategradientstabilized.cpp
OBJECTS_SOLVERS = solvers/assembly.o solvers/choleskidecomposed.o solvers/conjugategradient.o solvers/gausseidell.o solvers/polakribiereconjugategradient.o solvers/solver.o solvers/preconditionners.o solvers/gaussseidellstep.o solvers/inversediagonal.o solvers/incompletecholeskidecomposition.o solvers/biconjugategradientstabilized.o

SOURCE_NEW_GEO = geometry/geometry_base.cpp geometry/geometry_2D.cpp geometry/geometry_3D.cpp
OBJECTS_NEW_GEO = geometry/geometry_base.o geometry/geometry_2D.o geometry/geometry_3D.o

SOURCE_FEATURE = features/features.cpp features/pore.cpp features/pore3d.cpp features/inclusion.cpp features/inclusion3d.cpp features/sample.cpp features/sample3d.cpp features/crack.cpp features/enrichmentInclusion.cpp features/expansiveZone.cpp features/enrichmentbehaviour.cpp features/crackinitiation.cpp features/layeredinclusion.cpp
OBJECTS_FEATURE = features/features.o features/pore.o features/pore3d.o features/inclusion.o features/inclusion3d.o features/sample.o features/sample3d.o features/crack.o features/enrichmentInclusion.o features/expansiveZone.o features/enrichmentbehaviour.o features/crackinitiation.o features/layeredinclusion.o

SOURCE_FUNCTIONS = polynomial/vm_function_matrix.cpp polynomial/vm_token.cpp polynomial/vm_refcount_token.cpp polynomial/vm_function_base.cpp polynomial/vm_base.cpp 
OBJECTS_FUNCTIONS = polynomial/vm_function_matrix.o polynomial/vm_token.o polynomial/vm_refcount_token.o polynomial/vm_function_base.o polynomial/vm_base.o 

SOURCE_SPARSE = sparse/sparse_matrix.cpp sparse/sparse_vector.cpp
OBJECTS_SPARSE = sparse/sparse_matrix.o sparse/sparse_vector.o

SOURCE_MAIN_2D =  ${SOURCE_FUNCTIONS} ${SOURCE_PHYSICS} ${SOURCE_ELEMENTS} samplingcriterion.cpp  main.cpp matrixops.cpp utilities/granulo.cpp utilities/placement.cpp ${SOURCE_NEW_GEO}  ${SOURCE_SOLVERS} configuration.cpp ${SOURCE_FEATURE} delaunay_3d.cpp delaunay.cpp ${SOURCE_SPARSE}
OBJECTS_MAIN_2D =  ${OBJECTS_FUNCTIONS} ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS}   samplingcriterion.o main.o matrixops.o utilities/granulo.o utilities/placement.o ${OBJECTS_NEW_GEO} ${OBJECTS_SOLVERS} configuration.o ${OBJECTS_FEATURE} delaunay_3d.o delaunay.o ${OBJECTS_SPARSE}
TARGET_MAIN_2D = tryit

# TARGET_MAIN = steph main_simple 
SOURCE_MAIN_SIMPLE =  ${SOURCE_FUNCTIONS} ${SOURCE_PHYSICS} ${SOURCE_ELEMENTS} samplingcriterion.cpp  main_simple.cpp matrixops.cpp ${SOURCE_NEW_GEO}  ${SOURCE_SOLVERS} configuration.cpp ${SOURCE_FEATURE} delaunay_3d.cpp delaunay.cpp ${SOURCE_SPARSE}
OBJECTS_MAIN_SIMPLE =  ${OBJECTS_FUNCTIONS} ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS}   samplingcriterion.o main_simple.o matrixops.o  ${OBJECTS_NEW_GEO} ${OBJECTS_SOLVERS} configuration.o ${OBJECTS_FEATURE} delaunay_3d.o delaunay.o ${OBJECTS_SPARSE}
TARGET_MAIN_SIMPLE = simple

# TARGET_MAIN = steph main_kill
SOURCE_MAIN_KILL =  ${SOURCE_FUNCTIONS} ${SOURCE_PHYSICS} ${SOURCE_ELEMENTS} samplingcriterion.cpp  main_kill.cpp matrixops.cpp ${SOURCE_NEW_GEO}  ${SOURCE_SOLVERS} configuration.cpp ${SOURCE_FEATURE} delaunay_3d.cpp delaunay.cpp ${SOURCE_SPARSE}
OBJECTS_MAIN_KILL =  ${OBJECTS_FUNCTIONS} ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS}   samplingcriterion.o main_kill.o matrixops.o  ${OBJECTS_NEW_GEO} ${OBJECTS_SOLVERS} configuration.o ${OBJECTS_FEATURE} delaunay_3d.o delaunay.o ${OBJECTS_SPARSE}
TARGET_MAIN_KILL = killelements


SOURCE_MAIN_DIFFUSION =  ${SOURCE_FUNCTIONS} ${SOURCE_PHYSICS} ${SOURCE_ELEMENTS} samplingcriterion.cpp  main_diffusion.cpp matrixops.cpp ${SOURCE_NEW_GEO}  ${SOURCE_SOLVERS} configuration.cpp ${SOURCE_FEATURE} delaunay_3d.cpp delaunay.cpp ${SOURCE_SPARSE}
OBJECTS_MAIN_DIFFUSION =  ${OBJECTS_FUNCTIONS} ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS}   samplingcriterion.o main_diffusion.o matrixops.o  ${OBJECTS_NEW_GEO} ${OBJECTS_SOLVERS} configuration.o ${OBJECTS_FEATURE} delaunay_3d.o delaunay.o ${OBJECTS_SPARSE}
TARGET_MAIN_DIFFUSION = tryit_diffusion

SOURCE_MAIN_JEROME =  ${SOURCE_FUNCTIONS} ${SOURCE_PHYSICS} ${SOURCE_ELEMENTS} utilities/granulo.cpp samplingcriterion.cpp  main_jerome.cpp matrixops.cpp ${SOURCE_NEW_GEO}  ${SOURCE_SOLVERS} configuration.cpp ${SOURCE_FEATURE} delaunay_3d.cpp delaunay.cpp ${SOURCE_SPARSE}
OBJECTS_MAIN_JEROME =  ${OBJECTS_FUNCTIONS} ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS} utilities/granulo.o samplingcriterion.o main_jerome.o matrixops.o  ${OBJECTS_NEW_GEO} ${OBJECTS_SOLVERS} configuration.o ${OBJECTS_FEATURE} delaunay_3d.o delaunay.o ${OBJECTS_SPARSE}
TARGET_MAIN_JEROME = statistique

SOURCE_MAIN =  ${SOURCE_FUNCTIONS} ${SOURCE_PHYSICS} ${SOURCE_ELEMENTS}  samplingcriterion.cpp   main_3d.cpp matrixops.cpp ${SOURCE_NEW_GEO}  ${SOURCE_SOLVERS} configuration.cpp ${SOURCE_FEATURE} delaunay_3d.cpp delaunay.cpp ${SOURCE_SPARSE} ${SOURCE_FILTERS} utilities/granulo.cpp utilities/placement.cpp
OBJECTS_MAIN =  ${OBJECTS_FUNCTIONS} ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS}   samplingcriterion.o main_3d.o matrixops.o  ${OBJECTS_NEW_GEO} ${OBJECTS_SOLVERS} configuration.o ${OBJECTS_FEATURE} delaunay_3d.o delaunay.o ${OBJECTS_SPARSE} ${OBJECTS_FILTERS} utilities/granulo.o utilities/placement.o
TARGET_MAIN = tryit_3d

SOURCE_MAIN_3D_DIFFUSION =  ${SOURCE_FUNCTIONS} ${SOURCE_PHYSICS} ${SOURCE_ELEMENTS}  samplingcriterion.cpp   main_3d_diffusion.cpp matrixops.cpp ${SOURCE_NEW_GEO}  ${SOURCE_SOLVERS} configuration.cpp ${SOURCE_FEATURE} delaunay_3d.cpp delaunay.cpp ${SOURCE_SPARSE} ${SOURCE_FILTERS}
OBJECTS_MAIN_3D_DIFFUSION =  ${OBJECTS_FUNCTIONS} ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS}   samplingcriterion.o main_3d_diffusion.o matrixops.o  ${OBJECTS_NEW_GEO} ${OBJECTS_SOLVERS} configuration.o ${OBJECTS_FEATURE} delaunay_3d.o delaunay.o ${OBJECTS_SPARSE} ${OBJECTS_FILTERS}
TARGET_MAIN_3D_DIFFUSION = tryit_3d_diffusion

SOURCE_NURB =  main_nurbs.cpp matrixops.cpp ${SOURCE_NEW_GEO} 
OBJECTS_NURB = main_nurbs.o matrixops.o ${OBJECTS_NEW_GEO} 
TARGET_NURB = tryit_nurbs
 
SOURCE_AMOR_3D =  ${SOURCE_PHYSICS} ${SOURCE_ELEMENTS} main_amor.cpp matrixops.cpp ${SOURCE_FUNCTIONS} ${SOURCE_NEW_GEO} ${SOURCE_SOLVERS} samplingcriterion.cpp configuration.cpp ${SOURCE_FEATURE} delaunay.cpp delaunay_3d.cpp parser.cpp ${SOURCE_SPARSE}
OBJECTS_AMOR_3D = ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS} main_amor.o matrixops.o ${OBJECTS_FUNCTIONS} ${OBJECTS_NEW_GEO}  ${OBJECTS_SOLVERS} samplingcriterion.o configuration.o ${OBJECTS_FEATURE} delaunay.o  delaunay_3d.o parser.o ${OBJECTS_SPARSE}
TARGET_AMOR_3D = tryit_amor_3d

SOURCE_elasticity_c3s =  ${SOURCE_PHYSICS} ${SOURCE_ELEMENTS} main_elasticity_c3s.cpp matrixops.cpp ${SOURCE_FUNCTIONS} ${SOURCE_NEW_GEO} ${SOURCE_SOLVERS} samplingcriterion.cpp configuration.cpp ${SOURCE_FEATURE} delaunay.cpp delaunay_3d.cpp parser.cpp ${SOURCE_SPARSE}
OBJECTS_elasticity_c3s = ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS} main_elasticity_c3s.o matrixops.o ${OBJECTS_FUNCTIONS} ${OBJECTS_NEW_GEO} ${OBJECTS_SOLVERS} samplingcriterion.o configuration.o ${OBJECTS_FEATURE} delaunay_3d.o delaunay.o  parser.o ${OBJECTS_SPARSE}
TARGET_elasticity_c3s = tryit_elasticity_c3s 

SOURCE_periodic_c3s =  ${SOURCE_PHYSICS} ${SOURCE_ELEMENTS} main_periodic_c3s.cpp matrixops.cpp ${SOURCE_FUNCTIONS} ${SOURCE_NEW_GEO} ${SOURCE_SOLVERS} samplingcriterion.cpp configuration.cpp ${SOURCE_FEATURE} delaunay.cpp delaunay_3d.cpp parser.cpp ${SOURCE_SPARSE}
OBJECTS_periodic_c3s = ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS} main_periodic_c3s.o matrixops.o ${OBJECTS_FUNCTIONS} ${OBJECTS_NEW_GEO}  ${OBJECTS_SOLVERS} samplingcriterion.o configuration.o ${OBJECTS_FEATURE} delaunay_3d.o delaunay.o  parser.o ${OBJECTS_SPARSE}
TARGET_periodic_c3s = tryit_periodic_c3s

SOURCE_viscoelasticity_c3s =  ${SOURCE_PHYSICS} ${SOURCE_ELEMENTS} main_viscoelasticity_c3s.cpp matrixops.cpp ${SOURCE_FUNCTIONS} ${SOURCE_NEW_GEO} ${SOURCE_SOLVERS} samplingcriterion.cpp configuration.cpp ${SOURCE_FEATURE} delaunay.cpp delaunay_3d.cpp parser.cpp ${SOURCE_SPARSE}
OBJECTS_viscoelasticity_c3s = ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS} main_viscoelasticity_c3s.o matrixops.o ${OBJECTS_FUNCTIONS} ${OBJECTS_NEW_GEO} ${OBJECTS_SOLVERS} samplingcriterion.o configuration.o ${OBJECTS_FEATURE} delaunay_3d.o delaunay.o  parser.o ${OBJECTS_SPARSE}
TARGET_viscoelasticity_c3s = tryit_viscoelasticity_c3s 


${TARGET_MAIN}: ${OBJECTS_MAIN}
	${CXX} ${LDFLAGS} -o $@ $+

${TARGET_MAIN_3D_DIFFUSION}: ${OBJECTS_MAIN_3D_DIFFUSION}
	${CXX} ${LDFLAGS} -o $@ $+

${TARGET_MAIN_DIFFUSION}: ${OBJECTS_MAIN_DIFFUSION}
	${CXX} ${LDFLAGS} -o $@ $+

${TARGET_MAIN_2D}: ${OBJECTS_MAIN_2D}
	${CXX} ${LDFLAGS} -o $@ $+ 

${TARGET_MAIN_SIMPLE}: ${OBJECTS_MAIN_SIMPLE}
	${CXX} ${LDFLAGS} -o $@ $+ 

${TARGET_MAIN_KILL}: ${OBJECTS_MAIN_KILL}
	${CXX} ${LDFLAGS} -o $@ $+ 

${TARGET_MAIN_JEROME}: ${OBJECTS_MAIN_JEROME}
	${CXX} ${LDFLAGS} -o $@ $+

${TARGET_NEW_GEO}: ${OBJECTS_NEW_GEO}
	${CXX} ${LDFLAGS} -o $@ $+

${TARGET_AMOR_3D}: ${OBJECTS_AMOR_3D}
	${CXX} ${LDFLAGS} -o $@ $+ 

${TARGET_elasticity_c3s}: ${OBJECTS_elasticity_c3s}
	${CXX} ${LDFLAGS} -o $@ $+

${TARGET_periodic_c3s}: ${OBJECTS_periodic_c3s}
	${CXX} ${LDFLAGS} -o $@ $+

${TARGET_viscoelasticity_c3s}: ${OBJECTS_viscoelasticity_c3s} 
	${CXX} ${LDFLAGS} -o $@ $+

${TARGET_NURB}: ${OBJECTS_NURB}
	${CXX} ${LDFLAGS} -o $@ $+


depend_main: Makefile
	${CXX} -MM ${CXXFLAGS} ${SOURCE_MAIN} > $@

depend_main_3d_diffusion: Makefile
	${CXX} -MM ${CXXFLAGS} ${SOURCE_MAIN_3D_DIFFUSION} > $@
    
depend_main_2d: Makefile 
	${CXX} -MM ${CXXFLAGS} ${SOURCE_MAIN_2D} > $@

depend_main_diffusion: Makefile
	${CXX} -MM ${CXXFLAGS} ${SOURCE_MAIN_2D} > $@

depend_main_jerome: Makefile
	${CXX} -MM ${CXXFLAGS} ${SOURCE_MAIN_JEROME} > $@ 

depend_amor: Makefile
	${CXX} -MM ${CXXFLAGS} ${SOURCE_AMOR_3D} > $@

depend_elasticity_c3s: Makefile
	${CXX} -MM ${CXXFLAGS} ${SOURCE_elasticity_c3s} > $@

depend_periodic_c3s: Makefile 
	${CXX} -MM ${CXXFLAGS} ${SOURCE_periodic_c3s} > $@

depend_viscoelasticity_c3s: Makefile
	${CXX} -MM ${CXXFLAGS} ${SOURCE_viscoelasticity_c3s} > $@

depend_nurb: Makefile
	${CXX} -MM ${CXXFLAGS} ${SOURCE_NURB} > $@ 

depend_main_simple: Makefile
	${CXX} -MM ${CXXFLAGS} ${SOURCE_MAIN_SIMPLE} > $@

depend_main_kill: Makefile
	${CXX} -MM ${CXXFLAGS} ${SOURCE_MAIN_KILL} > $@

clean:
	rm -f ${OBJECTS_MAIN} ${TARGET_MAIN} ${TARGET_MAIN_3D_DIFFUSION} ${OBJECTS_AMOR_3D} ${TARGET_AMOR_3D}  ${OBJECTS_elasticity_c3s} ${TARGET_elasticity_c3s} ${OBJECTS_periodic_c3s} ${TARGET_periodic_c3s} ${OBJECTS_viscoelasticity_c3s} ${TARGET_viscoelasticity_c3s} ${OBJECTS_NEW_GEO} ${OBJECTS_FEATURE} ${OBJECTS_PHYSICS} ${OBJECTS_ELEMENTS} ${OBJECTS_SOLVERS} ${OBJECTS_SPARSE} ${OBJECTS_FILTERS} utilities/granulo.o depend_main  depend_amor depend_elasticity_c3s depend_periodic_c3s depend_viscoelasticity_c3s  depend_main_jerome depend_main_simple depend_main_kill 

include depend_main
include depend_main_3d_diffusion
include depend_amor
include depend_elasticity_c3s
include depend_periodic_c3s
include depend_viscoelasticity_c3s
include depend_nurb
include depend_main_jerome
include depend_main_2d 
include depend_main_simple
include depend_main_kill
